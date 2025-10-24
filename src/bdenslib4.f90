!===================================================================
!
!  file Bdenslib4.f90
!  modified from bapplyhlib.f90 
!  computes two-body densities between two vectors
!
!  started 2/2019 by CWJ
!
!===================================================================


module applytwobodydens
contains

!
!  boss routine
!
! note: default is going from vecin to vecout but this can be reversed depending on hchar
!

subroutine boss_den2b_bundled
  use nodeinfo
  use flagger
  use flags3body
  use precisions
  use opbundles
  use fragments
  use interaction
  use basis
  use localvectors
  use bmpi_mod
  use lanczos_info
  implicit none

  integer iprocs, procstart,procstop
  character(1) :: vchar  = 'n'
  real(kind=8) :: sum
  integer :: tid
  integer(kind=basis_prec) :: i, j
  integer :: ierr

!........ OPTION TO SIMULATE MPI ON A SINGLE "CORE".....

  if(distributeMPI .and. nproc == 1)then
      procstart = 0
      procstop  = nprocs -1
  else
      procstart = iproc
      procstop  = iproc
  end if

  call clocker('hmu','sta')
  do iprocs = procstart,procstop
     call proc_clock(iprocs,'sta')
	
!........... PP .................

     if(noisy0) print *, "Starting applyHbundled_g PP"
     call clocker('ppo','sta')
     call procOP_clock(iprocs,'sta','PPO')
     call applyhPPbundled_den('f',opbundlestart(iprocs), opbundleend(iprocs))
     call applyhPPbundled_den('b',opbundlestart(iprocs), opbundleend(iprocs))
	 
     call procOP_clock(iprocs,'end','PPO')

     call clocker('ppo','end')

!............NN ..................
     call clocker('nno','sta')
     call procOP_clock(iprocs,'sta','NNO')
	 
     call applyhNNbundled_den('f',opbundlestart(iprocs), opbundleend(iprocs))
     call applyhNNbundled_den('b',opbundlestart(iprocs), opbundleend(iprocs))
	 
     call clocker('nno','end')
     call procOP_clock(iprocs,'end','NNO')

!........... PN....................................
     call clocker('one','sta')
	 
	     call clocker('pno','sta')
	     call procOP_clock(iprocs,'sta','PNO')
        call applyhPNbundled_den('f',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNbundled_den('h',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNbundled_den('b',opbundlestart(iprocs), opbundleend(iprocs))

		call clocker('pno','end')
        call procOP_clock(iprocs,'end','PNO')

    call proc_clock(iprocs,'end')
	call clocker('hmu','end')
	
!........... NOW REDUCE....................

   call clocker('blr','sta')
   
#ifdef _MPI
   if(size(dmatpp) > 0)call BMPI_REDUCE(dmatpp,size(dmatpp),MPI_SUM,0,MPI_COMM_WORLD,ierr)
   if(size(dmatpphc) > 0)call BMPI_REDUCE(dmatpphc,size(dmatpp),MPI_SUM,0,MPI_COMM_WORLD,ierr)

   if(size(dmatnn) > 0)call BMPI_REDUCE(dmatnn,size(dmatnn),MPI_SUM,0,MPI_COMM_WORLD,ierr)
   if(size(dmatnnhc) > 0)call BMPI_REDUCE(dmatnnhc,size(dmatnn),MPI_SUM,0,MPI_COMM_WORLD,ierr)

   if(size(dmatpn) > 0)call BMPI_REDUCE(dmatpn,size(dmatpn),MPI_SUM,0,MPI_COMM_WORLD,ierr)
   if(size(dmatpnhc) > 0)call BMPI_REDUCE(dmatpnhc,size(dmatpnhc),MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif
   call clocker('blr','end')
   

  end do  ! iprocs

  return

end subroutine boss_den2b_bundled
!==================================================

!  subroutine applyhPPbundled
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)
!   vchar = 'n' (normal), 'r' (reverse)
!      (fragments) of lanczos vectors stored in module localvectors
!        in vec1 and vec2; if normal  H vec1 = vec2
!                          if reverse H vec2 = vec1
!
! cleaned up in 7.7.9: removed some loops that have become disused
!
!===================================================================
subroutine applyhPPbundled_den (hchar,startbundle,endbundle )

  use nodeinfo
  use localvectors
  use system_parameters
  use jumpNbody
  use precisions
  use interaction
  use opbundles
  use fragments
  use basis
  use lanczos_info
  use flagger
  use bmpi_mod
  use butil_mod
  use contigpointervectors, only : vecin,vecout, p2b_1sd,p2b_2sd
  implicit none

  logical :: pinfo
  integer :: ibundle
  character(1) :: hchar
  integer :: startbundle,endbundle

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=4)   xme, prod
  integer(kind=8) :: statei, statef,nsd,psdi,psdf
  integer(kind=basis_prec) :: statefoff, stateistart, stateistop
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(kind=8) :: vs
  integer(4) :: mythread

!..............................................................
!
!  SET UP POINTERS

vecin  => vec1
vecout => vec2

! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
if( hchar == 'f')then
           p2b_1sd => p2b_isd
           p2b_2sd => p2b_fsd
else
           p2b_1sd => p2b_fsd
           p2b_2sd => p2b_isd
endif
		
do ibundle = startbundle,endbundle

   if(opbundle(ibundle)%optype /= 'PP')cycle
   if(opbundle(ibundle)%hchar /= hchar )cycle

   if(diagonalsectorsonly .and. opbundle(ibundle)%isector /= opbundle(ibundle)%fsector)cycle
	 
   call bundle_clock(ibundle,'sta')

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
  csdstart = opbundle(ibundle)%nxstart
  csdend   = opbundle(ibundle)%nxend
  xjmpstart = opbundle(ibundle)%pxstart
  xjmpend   = opbundle(ibundle)%pxend
  cstride   = opbundle(ibundle)%cstride

  ncstates = (csdend +cstride -csdstart)/cstride
  
  select case (hchar)
  
  case ('f')  
  
!--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
!          this makes for simple OpenMP threading
!$omp parallel private(vs,xjmp, Xoplabel, xme, num_threads, mythread,nsd, pinfo)    &
!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!$omp          private(statefoff, stateistart, stateistop, prod)  &
!$omp          firstprivate(cstride, ncstates, xjmpstart,xjmpend)  &
!$omp          shared(vecin, vecout)  &
!$omp          shared(vec2threadchunkm1) &
!$omp          shared(p2b_op, p2b_1sd, p2b_2sd, p2b_phase), reduction(+:dmatpp)

     num_threads =  omp_get_num_threads()
     mythread = omp_get_thread_num()

! KSM:  chunks are guarenteed not to overlap on statef, so we
! KSM:  don't have to worry about about collisions between threads.
! KSM:  each thread gets a different range of neutron SDs.
  pinfo = ncstates > 10
  chunk = (ncstates + num_threads - 1)/num_threads
  istart = mythread*chunk + 1
  iend = bmin((mythread + 1)*chunk,ncstates)
  csd_index = csdstart + (istart - 1)*cstride - cstride
  if(istart <= iend)then	
     do csd = istart, iend
        csd_index = csd_index + cstride
        nsd = nstart(csd_index)

!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
        do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
            Xoplabel = p2b_op(xjmp)
!---------- GET INITIAL, FINAL SDs and place in basis..............

            statei = p2b_1sd(xjmp)+ nsd !csd_index
            statef = p2b_2sd(xjmp)+nsd !csd_index
		    xme =vecout(statef)*vecin(statei)*p2b_phase(xjmp)
		    dmatpp(Xoplabel)=dmatpp(Xoplabel)+xme

        end do  ! xjmp
      end do  ! csd
   end if
!$omp end parallel

  case ('b')  
  
!--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
!          this makes for simple OpenMP threading
!$omp parallel private(vs,xjmp, Xoplabel, xme, num_threads, mythread,nsd, pinfo)    &
!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!$omp          private(statefoff, stateistart, stateistop, prod)  &
!$omp          firstprivate(cstride, ncstates, xjmpstart,xjmpend)  &
!$omp          shared(vecin, vecout)  &
!$omp          shared(vec2threadchunkm1) &
!$omp          shared(p2b_op, p2b_1sd, p2b_2sd, p2b_phase), reduction(+:dmatpphc)

  num_threads =  omp_get_num_threads()
  mythread = omp_get_thread_num()

! KSM:  chunks are guarenteed not to overlap on statef, so we
! KSM:  don't have to worry about about collisions between threads.
! KSM:  each thread gets a different range of neutron SDs.
  pinfo = ncstates > 10
  chunk = (ncstates + num_threads - 1)/num_threads
  istart = mythread*chunk + 1
  iend = bmin((mythread + 1)*chunk,ncstates)
  csd_index = csdstart + (istart - 1)*cstride - cstride
  if(istart <= iend)then
	
  do csd = istart, iend
     csd_index = csd_index + cstride
     nsd = nstart(csd_index)

!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
          Xoplabel = p2b_op(xjmp)
!---------- GET INITIAL, FINAL SDs and place in basis..............

          statei = p2b_1sd(xjmp)+ nsd !csd_index
          statef = p2b_2sd(xjmp)+nsd !csd_index
		  xme =vecout(statef)*vecin(statei)*p2b_phase(xjmp)
		  dmatpphc(Xoplabel)=dmatpphc(Xoplabel)+xme

      end do  ! xjmp
   end do  ! csd
end if
   
!$omp end parallel

	case default
	
	print*,hchar, 'bad hchar !'
	stop
	
    end select
	
!--------------OR DO HERMITIAN/BACKWARDS APPLICATION----------

  call bundle_clock(ibundle,'end')

  end do ! ibundle
  return
end subroutine applyhPPbundled_den



!===================================================================
!  subroutine applyhNNbundled
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., NN between two (sub) sectors, etc)
!   vchar = 'n' (normal), 'r' (reverse)
!      (fragments) of lanczos vectors stored in module localvectors
!        in vec1 and vec2; if normal  H vec1 = vec2
!                          if reverse H vec2 = vec1
!   
!
!===================================================================
subroutine applyhNNbundled_den (hchar,startbundle,endbundle )
   use localvectors
   use nodeinfo
   use system_parameters
   use jumpNbody
   use precisions
   use interaction
   use opbundles
   use fragments
   use basis
   use lanczos_info
   use flagger
   use bmpi_mod
   use butil_mod
   use contigpointervectors, only : vecin,vecout, n2b_1sd,n2b_2sd
   implicit none

   ! arguments
   character(1),intent(in) :: hchar
   integer,intent(in) :: startbundle,endbundle

! --------- NOTE: basestart, basestop stored in module fragments

!------------------------------------------------------------

   integer :: ibundle
   integer(kind=8) csdstart, csdend,csd, csd_index,cstride,ncstates, pstride
   integer(kind=8) xjmp,xjmpstart,xjmpend
   integer(kind=8):: Xoplabel
   real(kind=4)   xme
   integer(kind=8) :: statei, statef,psd,nsdi,nsdf,statef_start,statef_end
   integer(kind=8) :: istart, iend, chunk
   integer(kind=8) :: vs

!-------- OpenMP functions ---------------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: num_threads
   integer :: mythread

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)

         vecin  => vec1
         vecout => vec2

         if( hchar == 'f' )then
            n2b_1sd => n2b_isd
            n2b_2sd => n2b_fsd
         else
            n2b_1sd => n2b_fsd
            n2b_2sd => n2b_isd
         endif


      do ibundle = startbundle,endbundle
         if(opbundle(ibundle)%optype /= 'NN')cycle
         if(opbundle(ibundle)%hchar /= hchar )cycle
		 
         if(diagonalsectorsonly .and. opbundle(ibundle)%isector /= opbundle(ibundle)%fsector)cycle
	     call bundle_clock(ibundle,'sta')
		 
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
         csdstart = opbundle(ibundle)%pxstart
         csdend   = opbundle(ibundle)%pxend
         cstride  = opbundle(ibundle)%cstride   !
         xjmpstart = opbundle(ibundle)%nxstart
         xjmpend   = opbundle(ibundle)%nxend
         ncstates = (csdend +cstride -csdstart)/cstride
		 
		 select case(hchar)
		 
!--------- OUTER LOOP OVER CONJUGATE PROTON SDs---------
!          this makes for simple OpenMP threading
!       NOTE CSTRIDE OVER PROTON SDs 


case ('f')
! firstprivate gives each thread its own copy, but initializes it
!    better than shared for read-only vars
! private gives each thread its own copy
!$omp parallel private(vs, xjmp, Xoplabel, xme, num_threads, mythread,psd)         &
!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!$omp          firstprivate(cstride, ncstates, xjmpstart, xjmpend)  &
!$omp          shared(vecin, vecout)  &
!$omp          shared(vec2threadchunkm1) &
!$omp          shared(n2b_op, n2b_1sd, n2b_2sd, n2b_phase), reduction(+:dmatnn)
         num_threads =  omp_get_num_threads()
         mythread = omp_get_thread_num()

         chunk = (ncstates + num_threads - 1)/num_threads
         istart = mythread*chunk + 1
         iend = bmin((mythread + 1)*chunk,ncstates)
         csd_index = csdstart + (istart - 1)*cstride - cstride
         if(istart <= iend)then

!..... THIS FOLLOWING IS TO TRY TO FIND AN OPTIMAL ORDERING OF LOOPS
            pstride = pstridecut +1 
            if(iend-istart > 0)pstride   = pstart(csdstart+cstride)-pstart(csdstart)


                  do csd = istart, iend
                     csd_index = csd_index + cstride
                     psd = pstart(csd_index)
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
                     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
                        Xoplabel = n2b_op(xjmp)

!---------- GET INITIAL, FINAL SDs and place in basis..............
                        statei = n2b_1sd(xjmp)+psd ! csd_index
                        statef = n2b_2sd(xjmp)+psd  !csd_index
			  		  xme =vecout(statef)*vecin(statei)*n2b_phase(xjmp)
			  		  dmatnn(Xoplabel)=dmatnn(Xoplabel)+xme			
!					  if(Xoplabel==2)print*,' f ',statei,statef,xme,ibundle	  
                     end do  ! xjmp
                  end do  ! csd

			  end if

!$omp end parallel

case ('b')
!$omp parallel private(vs, xjmp, Xoplabel, xme, num_threads, mythread,psd)         &
!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!$omp          firstprivate(cstride, ncstates, xjmpstart, xjmpend)  &
!$omp          shared(vecin, vecout)  &
!$omp          shared(vec2threadchunkm1) &
!$omp          shared(n2b_op, n2b_1sd, n2b_2sd, n2b_phase), reduction(+:dmatnnhc)
         num_threads =  omp_get_num_threads()
         mythread = omp_get_thread_num()

         chunk = (ncstates + num_threads - 1)/num_threads
         istart = mythread*chunk + 1
         iend = bmin((mythread + 1)*chunk,ncstates)
         csd_index = csdstart + (istart - 1)*cstride - cstride
         if(istart <= iend)then

!..... THIS FOLLOWING IS TO TRY TO FIND AN OPTIMAL ORDERING OF LOOPS
            pstride = pstridecut +1 
            if(iend-istart > 0)pstride   = pstart(csdstart+cstride)-pstart(csdstart)

                  do csd = istart, iend
                     csd_index = csd_index + cstride
                     psd = pstart(csd_index)
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
                     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
                        Xoplabel = n2b_op(xjmp)

!---------- GET INITIAL, FINAL SDs and place in basis..............
                        statei = n2b_1sd(xjmp)+psd ! csd_index
                        statef = n2b_2sd(xjmp)+psd  !csd_index
			  		  xme =vecout(statef)*vecin(statei)*n2b_phase(xjmp)
			  		  dmatnnhc(Xoplabel)=dmatnnhc(Xoplabel)+xme	
!					  if(Xoplabel==2)print*,' b ',statei,statef,xme	,ibundle	  
					  				  
                     end do  ! xjmp
                  end do  ! csd

			  end if

!$omp end parallel


end select


call bundle_clock(ibundle,'end')

   end do ! ibundle
   return
end subroutine applyhNNbundled_den

!=================================================================
!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPNbundled_den (hchar,startbundle,endbundle )
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use interaction
  use opbundles
  use fragments
  use lanczos_info
  use flagger
  use bmpi_mod
  use contigpointervectors, only : vecin,vecout
  implicit none

  integer :: ibundle,startbundle,endbundle
  character(1) :: hchar

!------------------------------------------------------------
  integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer :: a,b,c,d
  integer(kind=8) :: coplabel,doplabel
  integer :: phasep,phasen
  real(kind=4) ::   xme
  integer(kind=basis_prec) :: statei, statef
  integer(kind=4) num_threads

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  integer(kind=basis_prec) :: vs

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)

        vecin  => vec1
        vecout => vec2

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
	 if(diagonalsectorsonly .and. opbundle(ibundle)%isector /= opbundle(ibundle)%fsector)cycle
     call bundle_clock(ibundle,'sta')

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend
  
!  print*,ibundle, pjmpstart,pjmpend,njmpstart,njmpend, 'jmps'
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads

! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( hchar /= 'b' )then

!$omp parallel do private(vs, mythread,startp_thread,npjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)   &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          shared(vec2threadchunkm1) &
!$omp          firstprivate(njmpstart, njmpend)       &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)              &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!$omp          shared(cpnpair,dpnpair,vecin,vecout), reduction(+:dmatpn)
     do mythread = 0,numpthreads -1

     startp_thread = opbundle(ibundle)%startp_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     npjmps_thread = opbundle(ibundle)%startp_thread(mythread+1) - startp_thread

! KSM:  start/stop set up so that each proton final state appears on only one thread
! KSM:  prevents collison over update of voutp(statef) below
!---------   Forward direction ------------------------------
     do pjmp = startp_thread + 1, startp_thread + npjmps_thread
        psdi = p1b_isd(pjmp)       ! initial proton slater determinant
        psdf = p1b_fsd(pjmp)       ! final proton SD
        phasep = p1b_phase(pjmp)   ! phase of proton jumps
        a = p1b_cop(pjmp)     ! KSM: Proton 1-body creation label
        c = p1b_dop(pjmp)     ! KSM: Proton 1-body destruction label
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
        do njmp = njmpstart,njmpend
!----------- FIND MATRIX ELEMTN --------------------------------------------
           b = n1b_cop(njmp)  ! KSM: Neutron 1-body creation label
           d = n1b_dop(njmp)  ! KSM: Neutron 1-body destruction label
           phasen = n1b_phase(njmp)
           coplabel = cpnpair(b,a)
           doplabel = dpnpair(d,c)
!           xme = hmatpn(coplabel + doplabel)   ! get matrix element
!           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
           nsdi = n1b_isd(njmp)
           nsdf = n1b_fsd(njmp)
           statei = nsdi + psdi                ! initial state in combined basis
           statef = nsdf + psdf                ! final state in combined basis
 		  xme =vecout(statef)*vecin(statei)*phasep*phasen
!		  if(abs(a)==1 .and. abs(b)==1 .and. abs(c)==1 .and. abs(d)==1)then
!		  print*,' f/h ',a,b,c,d,xme, coplabel+doplabel
!		  endif 
		  		  dmatpn(coplabel+doplabel)=dmatpn(coplabel+doplabel)+xme
        end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do
else
!---- Backward direction using hermiticity ------------------- 

!$omp parallel do private(vs, mythread,startn_thread,nnjmps_thread)           &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)   &
!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!$omp          firstprivate(pjmpstart, pjmpend)                     &
!$omp          shared(vec2threadchunkm1)                            &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)    &
!$omp          shared(ibundle,opbundle)    &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)    &
!$omp          shared(cpnpair,dpnpair,vecin,vecout),reduction(+:dmatpnhc)
  do mythread = 0, numnthreads-1

     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread

!...... OPTION TO SWITCH ORDER OF LOOPS WHEN NO OpenMP.......

     if(numnthreads > 1 .or. disableNoOMPloopswitch)then

     do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
        nsdi = n1b_isd(njmp)
        nsdf = n1b_fsd(njmp)
        b  = n1b_cop(njmp)
        d  = n1b_dop(njmp)
        phasen = n1b_phase(njmp)
        do pjmp = pjmpstart,pjmpend
           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
           psdf = p1b_fsd(pjmp)       ! final proton SD
           phasep = p1b_phase(pjmp)   ! phase of proton jumps
           a  = p1b_cop(pjmp) 
           c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
!----------- FIND MATRIX ELEMENT -------------------------------------------
           coplabel = cpnpair(d,c)
           doplabel = dpnpair(b,a)		   
!           xme = hmatpn(coplabel + doplabel)     ! get matrix element
!           xme = xme*phasep*phasen               ! multiply matrix element by jump phases
           statei = nsdi + psdi                  ! initial state in combined basis
           statef = nsdf + psdf                  ! final state in combined basis
!           voutp(statei) = voutp(statei) + xme*vecin(statef)
		   
  		  xme =vecout(statei)*vecin(statef)*phasep*phasen
!		  print*,xme, statei,statef,coplabel+doplabel,' test pn  b'
!if(abs(a)==1 .and. abs(b)==1 .and. abs(c)==1 .and. abs(d)==1)then
!print*,' b1 ',a,b,c,d,xme,coplabel+doplabel,cpnpair(d,c)+dpnpair(b,a)
!endif		  
  		  dmatpn(coplabel+doplabel)=dmatpn(coplabel+doplabel)+xme
        end do  ! pjmp
     end do  ! njmp

     else


	     startp_thread = opbundle(ibundle)%startp_thread(mythread)     !  starting position for proton 1-body jumps for this thread
	     npjmps_thread = opbundle(ibundle)%startp_thread(mythread+1) - startp_thread
	     do pjmp = startp_thread + 1, startp_thread + npjmps_thread
	        psdi = p1b_isd(pjmp)       ! initial proton slater determinant
	        psdf = p1b_fsd(pjmp)       ! final proton SD
	        phasep = p1b_phase(pjmp)   ! phase of proton jumps
	        a = p1b_cop(pjmp)     ! KSM: Proton 1-body creation label
	        c = p1b_dop(pjmp)     ! KSM: Proton 1-body destruction label
	!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
	        do njmp = njmpstart,njmpend
	!----------- FIND MATRIX ELEMTN --------------------------------------------
	           b = n1b_cop(njmp)  ! KSM: Neutron 1-body creation label
	           d = n1b_dop(njmp)  ! KSM: Neutron 1-body destruction label
	           phasen = n1b_phase(njmp)
	           coplabel = cpnpair(d,c)
	           doplabel = dpnpair(b,a)
!	           xme = hmatpn(coplabel + doplabel)   ! get matrix element
!	           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
	           nsdi = n1b_isd(njmp)
	           nsdf = n1b_fsd(njmp)
	           statei = nsdi + psdi                ! initial state in combined basis
	           statef = nsdf + psdf                ! final state in combined basis
!	           voutp(statei) = voutp(statei) + xme*vecin(statef)
			   
	  		  xme =vecout(statei)*vecin(statef)*phasep*phasen
	  		  dmatpn(coplabel+doplabel)=dmatpn(coplabel+doplabel)+xme
!			  if(abs(a)==1 .and. abs(b)==1 .and. abs(c)==1 .and. abs(d)==1)then
!			  print*,' b2 ',a,b,c,d,xme,coplabel+doplabel
!			  endif
	        end do  ! njmp
	     end do  ! pjmp        

     end if
  end do
!$omp end parallel do


  end if
  call bundle_clock(ibundle,'end')

  end do  ! ibundle

  return
end subroutine applyhPNbundled_den
!===============================================================


end module applytwobodydens
