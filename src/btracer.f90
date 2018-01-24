!===================================================================
!
!  file BTRACER.f90
!
!  routines to compute traces
!  started 5/2012 by CWJ
!
!===================================================================


subroutine tracemaster
   use tracy
   use opbundles
   use basis
   use system_parameters
   use menu_choices
   use io
   implicit none
   integer(kind=8) :: icol, jcol, ncol
   real(kind=8) :: tr1,tr2,tmp
   integer(kind=basis_prec) :: ibasis
   integer :: aerr

   allocate( metrace( dimbasis), stat=aerr)
   if(aerr /= 0) call memerror("tracemaster 1")

   do ibasis = 1,dimbasis
      metrace(ibasis)%ncol = 0

   end do

   tracecount = .true.

   call applyHtrace

   noperation = 0
   do ibasis = 1,dimbasis
      noperation = noperation + metrace(ibasis)%ncol
      allocate( metrace(ibasis)%v(  metrace(ibasis)%ncol), stat=aerr)
      if(aerr /= 0) call memerror("tracemaster 2")
      allocate( metrace(ibasis)%col(  metrace(ibasis)%ncol), stat=aerr)
      if(aerr /= 0) call memerror("tracemaster 3")

      metrace(ibasis)%ncol = 0

   end do
   tracecount = .false.
   call applyHtrace

!............................. NOW SIMPLIFY...........
   print*,'check :',noperation

   noperation = 0
   do ibasis = 1, dimbasis

     do icol = 1, metrace(ibasis)%ncol

        if( metrace(ibasis)%v(icol)==0.0)cycle
        do jcol = icol+1,metrace(ibasis)%ncol
           if( metrace(ibasis)%col(icol) == metrace(ibasis)%col(jcol))then
               metrace(ibasis)%v(icol)= metrace(ibasis)%v(icol) + metrace(ibasis)%v(jcol)
               metrace(ibasis)%v(jcol) = 0.0
               metrace(ibasis)%col(jcol) = 0
           end if
        end do
     end do
!............ NOW COMPACT
     do icol = 1,metrace(ibasis)%ncol
        if( metrace(ibasis)%v(icol) == 0.0) then
           do jcol = icol+1, metrace(ibasis)%ncol
              if(metrace(ibasis)%v(jcol)/= 0.0)then
                  metrace(ibasis)%v(icol) = metrace(ibasis)%v(jcol)
                  metrace(ibasis)%col(icol) = metrace(ibasis)%col(jcol)
                  metrace(ibasis)%v(jcol) = 0.0
                  metrace(ibasis)%col(jcol) = 0
                  exit
              end if
           end do
        end if

     end do
     ncol = 0
     do icol = 1,metrace(ibasis)%ncol
         if(metrace(ibasis)%v(icol) /=0.0)ncol = ncol +1
     end do
     metrace(ibasis)%ncol = ncol
     noperation = ncol + noperation
   end do
   print*,' simplified; there are ',noperation,' unique matrix elements '


  tr1 = 0.d0
  tr2 = 0.d0

  do ibasis = 1,dimbasis
    do icol = 1,metrace(ibasis)%ncol
      if(ibasis == metrace(ibasis)%col(icol))then
         tr1 = tr1 + real(metrace(ibasis)%v(icol),8)
         tr2 = tr2 + real(metrace(ibasis)%v(icol),8)**2

      else
         tr2 = tr2 + 2.d0*real(metrace(ibasis)%v(icol),8)**2


      end if

    end do

  end do
  tr1 = tr1/dfloat(dimbasis)
  tr2 = tr2/dfloat(dimbasis)-tr1*tr1

202  format(' for dimension ',i8,', centroid = ',f12.6,', width = ',f12.6)
     write(6,202)dimbasis,tr1,dsqrt(tr2)
     write(6,*)' (saved to file trace.bigstick )'
     if ( writeout ) write(resultfile,202)dimbasis,tr1,dsqrt(tr2)
  open(unit=89,file='trace.bigstick',status='unknown')
  write(89,'(3i5,i8,2f20.6)')jz,np(1),np(2),dimbasis,tr1,tr2
  write(89,*)' 2xM    Z    N    basis           centroid       width^2 '
  close(89)
   return

end subroutine tracemaster
!===================================================================


!
! subroutine applyHtrace
!
! note: default is going from vecin to vecout but this can be reversed depending on hchar
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)
!   vchar = 'n' (normal), 'r' (reverse)
!      (fragments) of lanczos vectors stored in module localvectors
!        in vec1 and vec2; if normal  H vec1 = vec2
!                          if reverse H vec2 = vec1

subroutine applyHtrace
  use tracy
  use nodeinfo
  use flags3body
  use precisions
  use opbundles
  use fragments
!  use shampoo
  use interaction
  use basis
  use localvectors

  implicit none

  integer iprocs, procstart,procstop
  character(1) :: vchar  

!........ OPTION TO SIMULATE MPI ON A SINGLE "CORE".....

  if(distributeMPI .and. nproc == 1)then
      procstart = 0
      procstop  = nprocs -1
  else
      procstart = iproc
      procstop  = iproc
  end if

  ioperation = 0
  do iprocs = procstart,procstop
!............SPE......................................
     call applySPEbundled_trace(opbundlestart(iprocs), opbundleend(iprocs))

     if(.not.threebody)then
!........... PP .................

     call applyhPPbundled_trace('f',opbundlestart(iprocs), opbundleend(iprocs))
     call applyhPPbundled_trace('b',opbundlestart(iprocs), opbundleend(iprocs))

!............NN ..................

     call applyhNNbundled_trace('f',opbundlestart(iprocs), opbundleend(iprocs))
     call applyhNNbundled_trace('b',opbundlestart(iprocs), opbundleend(iprocs))


!........... PN....................................


        call applyhPNbundled_trace('f',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNbundled_trace('h',opbundlestart(iprocs), opbundleend(iprocs))
        call applyhPNbundled_trace('b',opbundlestart(iprocs), opbundleend(iprocs))


     else
!............PPP.......................................
 

    end if  ! if threebody



  end do  ! iprocs
  return

end subroutine applyHtrace

!===================================================================
!  subroutine applyhPPbundled
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)
!   vchar = 'n' (normal), 'r' (reverse)
!      (fragments) of lanczos vectors stored in module localvectors
!        in vec1 and vec2; if normal  H vec1 = vec2
!                          if reverse H vec2 = vec1

!
!===================================================================
subroutine applyhPPbundled_trace  (hchar,startbundle,endbundle )
  use tracy

  use nodeinfo
  use localvectors
  use system_parameters
!  use sectors
  use jumpNbody
  use precisions
  use interaction
!  use lanczos_info
  use opbundles
  use fragments
  use basis
  use lanczos_info
  use bmpi_mod
  use butil_mod
  implicit none

  integer :: ibundle
  character(1) :: hchar,vchar
  integer :: startbundle,endbundle
! --------- NOTE: basestart, basestop stored in module fragments
  real (kind = lanc_prec), pointer :: vecin(: ) 
  real (kind = lanc_prec), pointer :: vecout(: ) 

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=4)   xme
  integer(kind=8) :: statei, statef,nsd
  integer(kind=basis_prec), pointer :: p2b_1sd(:), p2b_2sd(:)
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=basis_prec) :: istart, iend, chunk
  integer(4) :: mythread
  integer irepeat


!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!

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

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
  csdstart = opbundle(ibundle)%nxstart
  csdend   = opbundle(ibundle)%nxend
  xjmpstart = opbundle(ibundle)%pxstart
  xjmpend   = opbundle(ibundle)%pxend
  cstride   = opbundle(ibundle)%cstride

  ncstates = (csdend +cstride -csdstart)/cstride

!--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
!       REMOVE OPENMP
!          this makes for simple OpenMP threading
!!$omp parallel private(xjmp, Xoplabel, xme, num_threads, mythread,nsd)         &
!!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!!$omp          shared(vecin, cstride, vecout,  ncstates,xjmpstart,xjmpend)  &
!!$omp          shared(p2b_op, p2b_1sd, p2b_2sd, p2b_phase, hmatpp)  &
!!$omp          reduction(+: metrace(:)%ncol)
  num_threads =  1  ! omp_get_num_threads()
  mythread = 0  ! omp_get_thread_num()
  chunk = (ncstates + num_threads - 1)/num_threads
  istart = mythread*chunk + 1
  iend = bmin((mythread + 1)*chunk,ncstates)
  csd_index = csdstart + (istart - 1)*cstride - cstride
  if(istart <= iend)then

  do csd = istart, iend
     csd_index = csd_index + cstride
     nsd = nstart(csd_index)

!  do csd = csdstart,csdend,cstride
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
          Xoplabel = p2b_op(xjmp)
          xme = hmatpp(Xoplabel)
!--------- GET PHASE.........................................
          xme = xme*p2b_phase(xjmp)
!---------- GET INITIAL, FINAL SDs and place in basis..............

          statei = p2b_1sd(xjmp)+ nsd !csd_index
          statef = p2b_2sd(xjmp)+nsd !csd_index

          if(statei > statef)cycle
          metrace(statei)%ncol = metrace(statei)%ncol + 1
          if(.not.tracecount)then
             metrace(statei)%v( metrace(statei)%ncol ) = xme
             metrace(statei)%col( metrace(statei)%ncol ) = statef

          end if


      end do  ! xjmp
   end do  ! csd
   end if
!!$omp end parallel
!--------------OR DO HERMITIAN/BACKWARDS APPLICATION----------

  end do ! ibundle
  return
end subroutine applyhPPbundled_trace
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

subroutine applyhNNbundled_trace (hchar,startbundle,endbundle )
  use tracy

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
  use bmpi_mod
  use butil_mod
  implicit none

  integer :: irepeat
  integer :: ibundle
   character(1) :: hchar,vchar
  integer startbundle,endbundle
! --------- NOTE: basestart, basestop stored in module fragments
  real (kind = lanc_prec), pointer :: vecin(: ) 
  real (kind = lanc_prec), pointer :: vecout(: ) 

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend,csd, csd_index,cstride,ncstates
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=4)   xme
  integer(kind=8) :: statei, statef,psd
  integer(kind=basis_prec), pointer :: n2b_1sd(:), n2b_2sd(:)

  integer(kind=4) num_threads
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(4) :: mythread

!..............................................................
!..............................................................
!
!  SET UP POINTERS
!   IF vchar = 'n' then H vec1 = vec2 
!        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
!            = 'b' then multiply H_ji vec1_i = vec2_j)
!
!   if vchar = 'r' then H vec2 = vec1)
!
!         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
!            = 'b' then multiply H_ij vec2_j = vec1_i  )
!  HERE i and j imply jumps between (sub) sectors
!
!
!
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
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
  csdstart = opbundle(ibundle)%pxstart
  csdend   = opbundle(ibundle)%pxend
  cstride  = opbundle(ibundle)%cstride   !
  xjmpstart = opbundle(ibundle)%nxstart
  xjmpend   = opbundle(ibundle)%nxend
  ncstates = (csdend +cstride -csdstart)/cstride
!--------- OUTER LOOP OVER CONJUGATE PROTON SDs---------
!          this makes for simple OpenMP threading
!       NOTE CSTRIDE OVER PROTON SDs 

!---- ELIMINATE OPENMP
!!$omp parallel private(xjmp, Xoplabel, xme, num_threads, mythread,psd)         &
!!$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
!!$omp          shared(vecin, cstride, vecout,  ncstates,xjmpstart,xjmpend)  &
!!$omp          shared(n2b_op, n2b_1sd, n2b_2sd, n2b_phase, hmatnn)
  num_threads =  1 !omp_get_num_threads()
  mythread = 0 !  omp_get_thread_num()
  chunk = (ncstates + num_threads - 1)/num_threads
  istart = mythread*chunk + 1
  iend = bmin((mythread + 1)*chunk,ncstates)
  csd_index = csdstart + (istart - 1)*cstride - cstride
  if(istart <= iend)then
  do csd = istart, iend
     csd_index = csd_index + cstride
     psd = pstart(csd_index)
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
        Xoplabel = n2b_op(xjmp)
        xme = hmatnn(Xoplabel)
!--------- GET PHASE.........................................
        xme = xme*n2b_phase(xjmp)
!---------- GET INITIAL, FINAL SDs and place in basis..............
        statei = n2b_1sd(xjmp)+psd ! csd_index
        statef = n2b_2sd(xjmp)+psd  !csd_index


          if(statei > statef)cycle
          metrace(statei)%ncol = metrace(statei)%ncol + 1
          if(.not.tracecount)then
             metrace(statei)%v( metrace(statei)%ncol ) = xme
             metrace(statei)%col( metrace(statei)%ncol ) = statef

          end if

        end do  ! xjmp
   end do  ! csd
   end if
!!$omp end parallel

  end do ! ibundle
  return
end subroutine applyhNNbundled_trace
!=================================================================
!
! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPNbundled_trace (hchar,startbundle,endbundle )
  use tracy
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use interaction
  use opbundles
  use fragments
  use lanczos_info
  use bmpi_mod
  implicit none

  integer :: irepeat
  integer :: ibundle,startbundle,endbundle
   character(1) :: hchar,vchar
! --------- NOTE: basestart, basestop stored in module fragments
  real (kind = lanc_prec), pointer :: vecin(: ) 
  real (kind = lanc_prec), pointer :: vecout(: ) 

!------------------------------------------------------------
  integer(kind=8) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer a,b,c,d
  integer(kind=8) coplabel,doplabel
  integer :: phasep,phasen
  real(kind=4)   xme
  integer(kind=8) :: statei, statef
  integer(kind=4) num_threads

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  logical    :: launched(0:3)

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads

!  print*,numberthreads,' threads '
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( ( hchar /= 'b')  )then
!--- ELIMINATE OPENMP
!!$omp parallel do private(mythread,startp_thread,npjmps_thread)                  &
!!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)          &
!!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)              &
!!$omp          shared(ibundle,opbundle)    &
!!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)              &
!!$omp          shared(cpnpair,dpnpair,vecin,vecout)
!     mythread = omp_get_thread_num()
     do mythread = 0,numpthreads -1
     startp_thread = opbundle(ibundle)%startp_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     npjmps_thread = opbundle(ibundle)%startp_thread(mythread+1) - startp_thread

!---------   Forward direction ------------------------------
     do pjmp = startp_thread + 1, startp_thread + npjmps_thread
        psdi = p1b_isd(pjmp)       ! initial proton slater determinant
        psdf = p1b_fsd(pjmp)       ! final proton SD
        phasep = p1b_phase(pjmp)   ! phase of proton jumps
        a = p1b_cop(pjmp) 
        c = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
        do njmp = njmpstart,njmpend
!----------- FIND MATRIX ELEMTN --------------------------------------------
           b = n1b_cop(njmp)
           d = n1b_dop(njmp)
           phasen = n1b_phase(njmp)
           coplabel = cpnpair(b,a)
           doplabel = dpnpair(d,c)
           xme = hmatpn(coplabel + doplabel)   ! get matrix element
           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
           nsdi = n1b_isd(njmp)
           nsdf = n1b_fsd(njmp)
           statei = nsdi + psdi                ! initial state in combined basis
           statef = nsdf + psdf                ! final state in combined basis


          if(statei > statef)cycle
          metrace(statei)%ncol = metrace(statei)%ncol + 1
          if(.not.tracecount)then
             metrace(statei)%v( metrace(statei)%ncol ) = xme
             metrace(statei)%col( metrace(statei)%ncol ) = statef
          end if

        end do  ! njmp
     end do  ! pjmp        
  end do
!!$omp end parallel
else
!---- Backward direction using hermiticity ------------------- 

!!$omp parallel do private(mythread,startn_thread,nnjmps_thread)                  &
!!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)      &
!!$omp          private(coplabel,doplabel,xme,statei,statef)                   &
!!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)    &
!!$omp          shared(ibundle,opbundle)    &
!!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)    &
!!$omp          shared(cpnpair,dpnpair,vecin,vecout)
!     mythread = omp_get_thread_num()
     do mythread = 0, numnthreads-1
     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread
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
           coplabel = cpnpair(b,a)
           doplabel = dpnpair(d,c)
           xme = hmatpn(coplabel + doplabel)     ! get matrix element
           xme = xme*phasep*phasen               ! multiply matrix element by jump phases
           statei = nsdi + psdi                  ! initial state in combined basis
           statef = nsdf + psdf                  ! final state in combined basis

          if(statef > statei)cycle
          metrace(statef)%ncol = metrace(statef)%ncol + 1
          if(.not.tracecount)then
             metrace(statef)%v( metrace(statef)%ncol ) = xme
             metrace(statef)%col( metrace(statef)%ncol ) = statei

          end if

        end do  ! pjmp
     end do  ! njmp
  end do
!!$omp end parallel


  end if

  end do  ! ibundle

  return
end subroutine applyhPNbundled_trace

!==========================================================


!==================================================================
!
!  subroutine applyspes
!
!  applies single-particle energies -- purely diagonal
!
!====================================================================
subroutine applySPEbundled_trace (startbundle,endbundle )
  use tracy

  use basis
  use sectors
  use diagh
  use precisions
  use lanczos_info
  use localvectors
  use nodeinfo
  use system_parameters
  use opbundles
  use fragments
  use bmpi_mod
  implicit none

  integer :: irepeat
  integer :: ibundle,startbundle,endbundle
  character(1) :: vchar
! --------- NOTE: basestart, basestop stored in module fragments
  real (kind = lanc_prec), pointer :: vecin(: ) 
  real (kind = lanc_prec), pointer :: vecout(: ) 

!------------------------------------------------------------

  integer(kind=8) nsdstart,nsdend,psdstart,psdend
  integer(kind=8) xjmp,xjmpstart,xjmpend

  integer(kind=8) :: ibasis
  integer(kind=4) num_threads
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads

  real(4)    ::  pspe,nspe
  integer(kind=basis_prec)   :: ip,in



  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'SPE')cycle

  psdstart = opbundle(ibundle)%pxstart
  psdend   = opbundle(ibundle)%pxend
  nsdstart = opbundle(ibundle)%nxstart
  nsdend   = opbundle(ibundle)%nxend

!!omp parallel do private(ip,in,pspe,nspe,ibasis) & 
!!omp  shared(psdstart,psdend,nsdstart,nsdend,vecin,vecout) & 
!!omp  shared(pstart,nstart,pspe_h,nspe_h)

!................ LOOP OVER PROTON SDs in that sector..........
  do ip = psdstart,psdend
           pspe = pspe_h(ip)   ! the proton contribution to the single-particle energies
           ibasis = pstart(ip) + nstart(nsdstart)

!............... LOOP OVER NEUTRON SDS in sector jsc........................
           do in = nsdstart,nsdend
                 nspe = nspe_h(in)  ! neutron contributions to s.p.e.
!             ibasis = pstart(ip) + nstart(in)  ! find basis state
          metrace(ibasis)%ncol = metrace(ibasis)%ncol + 1
          if(.not.tracecount)then
             metrace(ibasis)%v( metrace(ibasis)%ncol ) = pspe+nspe
             metrace(ibasis)%col( metrace(ibasis)%ncol ) = ibasis
          end if
                 ibasis = ibasis+1
           end do  !in
  end do  !ip


  end do ! ibundle
  return
end subroutine applySPEbundled_trace

!===============================================================
! ADDED IN 7.6.3 for computing the sector trace (PN only)
! ADDED IN 7.6.6 also PP, NN
!=================================================================
!
!  master routine for compute pn traces over sectors
!  for use in generating pivots
!
!  CALLED BY
!
!  SUBROUTINES CALLED:
!    applyPNtraceonly
!
subroutine masterpnsectortrace

	use tracy
	use nodeinfo
	use opbundles
	use sectors
	use bmpi_mod
	
	implicit none
	integer :: ibundle
	integer :: isector
	integer(kind=8) :: localdim
	real, allocatable :: tracedummy(:)
	integer :: ierr
	
	allocate(sectortrace(nsectors(1)),tracedummy(nsectors(1)))
	
	sectortrace(:)=0.0
	
    call applyhPNtraceonly('f',opbundlestart(iproc), opbundleend(iproc))
    call applyhPNtraceonly('h',opbundlestart(iproc), opbundleend(iproc))
    call applyhPNtraceonly('b',opbundlestart(iproc), opbundleend(iproc))
	
!........... REDUCE FROM ALL MPI PROCESSES AND DIVIDE BY DIMENSIONS  .......	
	
	call BMPI_ALLREDUCE(sectortrace,tracedummy, nsectors(1),  MPI_SUM, icomm, ierr)
    do isector = 1,nsectors(1)
		localdim = 1+xsd(1)%sector(isector)%basisend-xsd(1)%sector(isector)%basisstart
		sectortrace(isector)=tracedummy(isector)/real( localdim)
		
    end do ! isector
	return

end subroutine masterpnsectortrace
!=================================================================
!
!  master routine for compute traces over sectors
!  for use in generating preconditioners
!
!  CALLED BY
!
!  SUBROUTINES CALLED:
!    applyPNtraceonly
!
subroutine mastersectortrace

	use tracy
	use nodeinfo
	use opbundles
	use sectors
	use bmpi_mod
	use basis
	
	implicit none
	integer :: ibundle
	integer :: isector
	integer(kind=8) :: localdim
	real, allocatable :: tracedummy(:)
	integer :: ierr
	real(kind=8) :: testtrace
	
	allocate(sectortrace(nsectors(1)),tracedummy(nsectors(1)))
	
	sectortrace(:)=0.0
    call applyhSPEtraceonly(opbundlestart(iproc), opbundleend(iproc))
	
    call applyhPNtraceonly('f',opbundlestart(iproc), opbundleend(iproc))
    call applyhPNtraceonly('h',opbundlestart(iproc), opbundleend(iproc))
    call applyhPNtraceonly('b',opbundlestart(iproc), opbundleend(iproc))
	
    call applyhPPtraceonly('f',opbundlestart(iproc), opbundleend(iproc))
    call applyhPPtraceonly('b',opbundlestart(iproc), opbundleend(iproc))
    call applyhNNtraceonly('f',opbundlestart(iproc), opbundleend(iproc))
    call applyhNNtraceonly('b',opbundlestart(iproc), opbundleend(iproc))
	
!........... REDUCE FROM ALL MPI PROCESSES AND DIVIDE BY DIMENSIONS  .......	
	
	call BMPI_ALLREDUCE(sectortrace,tracedummy, nsectors(1),  MPI_SUM, icomm, ierr)
	testtrace = 0.d0
    do isector = 1,nsectors(1)
		localdim = 1+xsd(1)%sector(isector)%basisend-xsd(1)%sector(isector)%basisstart
		testtrace = testtrace + tracedummy(isector)
		sectortrace(isector)=tracedummy(isector)/real( localdim)
    end do ! isector
	testtrace= testtrace / real(dimbasis,kind=8)
	if(iproc==0)print*,' as a test, centroid = ',testtrace
	return

end subroutine mastersectortrace
!=================================================================
subroutine applyhSPEtraceonly(startbundle,endbundle )
   use basis
   use sectors
   use diagh
   use precisions
   use lanczos_info
   use nodeinfo
   use system_parameters
   use opbundles
   use fragments
   use tracy
   implicit none

   integer :: ibundle,startbundle,endbundle
   character(1) :: vchar
!------------------------------------------------------------

   integer(kind=8) nsdstart,nsdend,psdstart,psdend
   integer(kind=8) xjmp,xjmpstart,xjmpend

   integer(kind=8) :: ibasis
   integer(kind=4) num_threads
!-------- OpenMP functions ---------------------------------
   integer(kind=4) :: omp_get_thread_num, omp_get_num_threads

   real(4)    ::  pspe,nspe
   integer(kind=8) :: ip,in
   
   real(kind=4) :: localsum
   integer      :: localsector
   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'SPE')cycle
      localsum = 0.0
 	  localsector = opbundle(ibundle)%isector

      psdstart = opbundle(ibundle)%pxstart
      psdend   = opbundle(ibundle)%pxend
      nsdstart = opbundle(ibundle)%nxstart
      nsdend   = opbundle(ibundle)%nxend

! note: first private says that each such private copy is initialized with 
! the value before the parallel pragma
!................ LOOP OVER PROTON SDs in that sector..........
!$omp parallel do private(ip,in,pspe,nspe,ibasis) & 
!$omp  firstprivate(nsdstart, nsdend) &
!$omp  shared(psdstart,psdend) & 
!$omp  shared(pstart,nstart,pspe_h,nspe_h) &
!$omp      reduction(+:localsum)

      do ip = psdstart,psdend
           pspe = pspe_h(ip)   ! the proton contribution to the single-particle energies
           ibasis = pstart(ip) + nstart(nsdstart)
!............... LOOP OVER NEUTRON SDS in sector jsc........................
           do in = nsdstart,nsdend
              nspe = nspe_h(in)  ! neutron contributions to s.p.e.
			  localsum = localsum + (pspe+nspe)
              ibasis = ibasis+1   ! neutron SDs are contiguous
           end do  !in
      end do  !ip
!$omp end parallel do
  sectortrace(localsector)=sectortrace(localsector)+localsum

   end do ! ibundle
   return
end subroutine applyhSPEtraceonly
!========================================================

! NOTE for OpenMP:  the 1-body jumps are sorted as follows:
!      protons on final states
!      neutrons on "initial" states
! 
!
subroutine applyhPNtraceonly(hchar,startbundle,endbundle )
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
  use tracy
  implicit none

  integer :: irepeat
  integer :: ibundle,startbundle,endbundle
   character(1) :: hchar,vchar

!------------------------------------------------------------
  integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer :: a,b,c,d
  integer(kind=8) :: coplabel,doplabel
  integer :: phasep,phasen
  real(kind=4) :: localsum
  integer      :: localsector
  real(kind=4) ::   xme
  integer(kind=basis_prec) :: statei, statef
  integer(kind=4) num_threads

!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(4) :: mythread,numpthreads,numnthreads
  integer(8) :: startp_thread, npjmps_thread
  integer(8) :: startn_thread, nnjmps_thread
  logical    :: launched(0:3)

  localsum = 0.0

  do ibundle = startbundle,endbundle
     if(opbundle(ibundle)%optype /= 'PN')cycle
     if(opbundle(ibundle)%hchar /= hchar )cycle
	 
	 if(opbundle(ibundle)%isector/=opbundle(ibundle)%fsector)cycle   
     localsum = 0.0
	 localsector = opbundle(ibundle)%isector
	 
!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend  
  numpthreads = opbundle(ibundle)%numpthreads
  numnthreads = opbundle(ibundle)%numnthreads

!  print*,numberthreads,' threads '
! NOTE: 
!   hchar = 'f' (forwards), 'b' (backwards)
!       This relates to v_i = H_ij v_j (forwards)
!       and its conjugate v_j = H_ji v_i
!
  if( hchar /= 'b' )then

!$omp parallel do private(mythread,startp_thread,npjmps_thread) &
!$omp      private(pjmp,njmp,psdi,psdf,nsdi,nsdf)               &
!$omp      private(phasep,phasen,a,b,c,d)                       &
!$omp      private(coplabel,doplabel,xme,statei,statef)         &
!$omp      firstprivate(njmpstart, njmpend)                     &
!$omp      shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)    &
!$omp      shared(ibundle,opbundle)                             &
!$omp      shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)    &
!$omp      shared(cpnpair,dpnpair)                              & 
!$omp      reduction(+:localsum)
     do mythread = 0,numpthreads -1
     startp_thread = opbundle(ibundle)%startp_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     npjmps_thread = opbundle(ibundle)%startp_thread(mythread+1) - startp_thread

! KSM:  start/stop set up so that each proton final state appears on only one thread
! KSM:  prevents collison over update of vecout(statef) below
!---------   Forward direction ------------------------------
     do pjmp = startp_thread + 1, startp_thread + npjmps_thread
        psdi = p1b_isd(pjmp)       ! initial proton slater determinant
        psdf = p1b_fsd(pjmp)       ! final proton SD
		if(psdf/=psdi)cycle
        phasep = p1b_phase(pjmp)   ! phase of proton jumps
        a = p1b_cop(pjmp)     ! KSM: Proton 1-body creation label
        c = p1b_dop(pjmp)     ! KSM: Proton 1-body destruction label
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
        do njmp = njmpstart,njmpend
           nsdi = n1b_isd(njmp)
           nsdf = n1b_fsd(njmp)
		   if(nsdi/=nsdf)cycle
!----------- FIND MATRIX ELEMTN --------------------------------------------
           b = n1b_cop(njmp)  ! KSM: Neutron 1-body creation label
           d = n1b_dop(njmp)  ! KSM: Neutron 1-body destruction label
           phasen = n1b_phase(njmp)
           coplabel = cpnpair(b,a)
           doplabel = dpnpair(d,c)
           xme = hmatpn(coplabel + doplabel)   ! get matrix element
           xme = xme*phasep*phasen             ! multiply matrix element by jump phases
!           statei = nsdi + psdi                ! initial state in combined basis
!           statef = nsdf + psdf                ! final state in combined basis
           localsum = localsum + xme
        end do  ! njmp
     end do  ! pjmp        
  end do
!$omp end parallel do
else
!---- Backward direction using hermiticity ------------------- 

!$omp parallel do private(mythread,startn_thread,nnjmps_thread)              &
!$omp          private(pjmp,njmp,psdi,psdf,nsdi,nsdf,phasep,phasen,a,b,c,d)  &
!$omp          private(coplabel,doplabel,xme,statei,statef)                  &
!$omp          firstprivate(pjmpstart, pjmpend)                              &
!$omp          shared(p1b_isd,p1b_fsd,p1b_phase,p1b_cop,p1b_dop)             &
!$omp          shared(ibundle,opbundle)                                      &
!$omp          shared(n1b_isd,n1b_fsd,n1b_phase,n1b_cop,n1b_dop)             &
!$omp          shared(cpnpair,dpnpair)                                       &
!$omp          reduction(+:localsum)

     do mythread = 0, numnthreads-1
     startn_thread = opbundle(ibundle)%startn_thread(mythread)     !  starting position for proton 1-body jumps for this thread
     nnjmps_thread = opbundle(ibundle)%startn_thread(mythread+1) - startn_thread

!...... OPTION TO SWITCH ORDER OF LOOPS WHEN NO OpenMP.......

     if(numnthreads > 1 .or. disableNoOMPloopswitch)then

     do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
        nsdi = n1b_isd(njmp)
        nsdf = n1b_fsd(njmp)
		if(nsdf/=nsdi)cycle
        b  = n1b_cop(njmp)
        d  = n1b_dop(njmp)
        phasen = n1b_phase(njmp)
        do pjmp = pjmpstart,pjmpend
           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
           psdf = p1b_fsd(pjmp)       ! final proton SD
		   if(psdi/=psdf)cycle
           phasep = p1b_phase(pjmp)   ! phase of proton jumps
           a  = p1b_cop(pjmp) 
           c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
!----------- FIND MATRIX ELEMENT -------------------------------------------
           coplabel = cpnpair(b,a)
           doplabel = dpnpair(d,c)
!		   if(coplabel+doplabel==0)then
!			   print*,iproc,'(b)',njmp,coplabel,doplabel
!			   print*,a,b,c,d,ibundle
!		   end if		   
           xme = hmatpn(coplabel + doplabel)     ! get matrix element
           xme = xme*phasep*phasen               ! multiply matrix element by jump phases
!           statei = nsdi + psdi                  ! initial state in combined basis
!           statef = nsdf + psdf                  ! final state in combined basis
           localsum = localsum + xme
        end do  ! pjmp
     end do  ! njmp

     else

     do pjmp = pjmpstart,pjmpend
           psdi = p1b_isd(pjmp)       ! initial proton slater determinant
           psdf = p1b_fsd(pjmp)       ! final proton SD
		   if(psdi/=psdf)cycle
           phasep = p1b_phase(pjmp)   ! phase of proton jumps
           a  = p1b_cop(pjmp) 
           c  = p1b_dop(pjmp)
!--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
           do njmp = startn_thread + 1, startn_thread + nnjmps_thread     
              nsdi = n1b_isd(njmp)
              nsdf = n1b_fsd(njmp)
	  		  if(nsdf/=nsdi)cycle
			  
              b  = n1b_cop(njmp)
              d  = n1b_dop(njmp)
              phasen = n1b_phase(njmp)

!----------- FIND MATRIX ELEMENT -------------------------------------------
              coplabel = cpnpair(b,a)
              doplabel = dpnpair(d,c)
              xme = hmatpn(coplabel + doplabel)     ! get matrix element
              xme = xme*phasep*phasen               ! multiply matrix element by jump phases
              localsum = localsum + xme
           end do
      end do  ! pjmp


     end if
  end do
!$omp end parallel do


  end if
  sectortrace(localsector)=sectortrace(localsector)+localsum
  end do  ! ibundle

  return
end subroutine applyhPNtraceonly

!==========================================================
subroutine applyhPPtraceonly (hchar,startbundle,endbundle )


    use nodeinfo
    use localvectors
    use system_parameters
  !  use sectors
    use jumpNbody
    use precisions
    use interaction
  !  use lanczos_info
    use opbundles
    use fragments
    use basis
    use lanczos_info
    use flagger
    use bmpi_mod
    use butil_mod
	use tracy
    use contigpointervectors, only : p2b_1sd,p2b_2sd
	
    implicit none

    logical :: pinfo
    integer :: ibundle
    character(1) :: hchar
    integer :: startbundle,endbundle
    real(kind=4) :: localsum
    integer      :: localsector

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
!    real(kind=lanc_prec), pointer :: voutp(:)


  !..............................................................
  !
  !  SET UP POINTERS
  !   IF vchar = 'n' then H vec1 = vec2 
  !        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
  !            = 'b' then multiply H_ji vec1_i = vec2_j)
  !
  !   if vchar = 'r' then H vec2 = vec1)
  !
  !         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
  !            = 'b' then multiply H_ij vec2_j = vec1_i  )
  !  HERE i and j imply jumps between (sub) sectors
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
	 
	 if(opbundle(ibundle)%isector/=opbundle(ibundle)%fsector)cycle   
     localsum = 0.0
	 localsector = opbundle(ibundle)%isector
	 

  !...... EXTRACT INFORMATION FROM OPBUNDLE ........
    csdstart = opbundle(ibundle)%nxstart
    csdend   = opbundle(ibundle)%nxend
    xjmpstart = opbundle(ibundle)%pxstart
    xjmpend   = opbundle(ibundle)%pxend
    cstride   = opbundle(ibundle)%cstride

    ncstates = (csdend +cstride -csdstart)/cstride

  !--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
  !          this makes for simple OpenMP threading
  !$omp parallel private(vs,xjmp, Xoplabel, xme, num_threads, mythread,nsd, pinfo)    &
  !$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
  !$omp          private(statefoff, stateistart, stateistop, prod)  &
  !$omp          firstprivate(cstride, ncstates, xjmpstart,xjmpend)  &
  !$omp          shared(vec2threadchunkm1) &
  !$omp          shared(p2b_op, p2b_1sd, p2b_2sd, p2b_phase, hmatpp) & 
  !$omp          reduction(+:localsum)
  
    num_threads =  omp_get_num_threads()
    mythread = omp_get_thread_num()
    ! thread local vec2, reduce at end
    if(useVec2Thread) then
      !! voutp(v2s:v2e) => vec2thread(:, mythread)
 !     vs = mythread * vec2threadchunk;
  !    voutp(v2s:v2e) => vec2threadflat(vs : vs + vec2threadchunkm1)
    else
   !   voutp(v2s:v2e) => vecout
    end if

  ! KSM:  chunks are guarenteed not to overlap on statef, so we
  ! KSM:  don't have to worry about about collisions between threads.
  ! KSM:  each thread gets a different range of neutron SDs.
    pinfo = ncstates > 10
    chunk = (ncstates + num_threads - 1)/num_threads
    istart = mythread*chunk + 1
    iend = bmin((mythread + 1)*chunk,ncstates)
    csd_index = csdstart + (istart - 1)*cstride - cstride
    if(istart <= iend)then

  !......... THERE ARE TWO VERSIONS.....
  !          1st way: store in jumps index to PP matrix element;
  !          this has faster setup
  !          2nd way: store PP matrix elements directly;
  !          slower set up, but on MPI nodes reduced memory load

  if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO PP MATRIX ELEMENTS

    if(num_threads > 1 .or. disableNoOMPloopswitch)then
  ! if(pinfo) print *, "plpx: thread=", mythread

  if(.true.) then
    do csd = istart, iend
       csd_index = csd_index + cstride
       nsd = nstart(csd_index)

  !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
       do xjmp = xjmpstart,xjmpend

!---------- GET INITIAL, FINAL SDs and place in basis..............

		    statei = p2b_1sd(xjmp)+ nsd !csd_index
		    statef = p2b_2sd(xjmp)+nsd !csd_index		   
			if(statei/=statef)cycle
  !--------- FETCH MATRIX ELEMENT...............................
            Xoplabel = p2b_op(xjmp)
            xme = hmatpp(Xoplabel)
  !--------- GET PHASE.........................................
            xme = xme*p2b_phase(xjmp)

			localsum = localsum + xme
   !         voutp(statef) = voutp(statef) + xme*vecin(statei)
  !write(54,'(4i4,f10.5)')ibundle,xjmp,statei,statef,xme
  !write(6,'(4i4,f10.5)')ibundle,xjmp,statei,statef,xme
        end do  ! xjmp
     end do  ! csd
  else
     ! Trial speedup
     ! result so for is mysteriously slower
  !    if(pinfo) print *, "plpt: thread=", mythread
       csd_index = nstart(csd_index+1)-1
       do xjmp = xjmpstart,xjmpend
  !--------- FETCH MATRIX ELEMENT...............................
            Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
            xme = hmatpp(Xoplabel)
  !--------- GET PHASE.........................................
            xme = xme*p2b_phase(xjmp)
            psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
            psdf = p2b_2sd(xjmp)   ! KSM: final P SD
			if(psdi/=psdf)cycle
       
            ! nsd = csd_index
  !---------- GET INITIAL, FINAL SDs and place in basis..............
            ! do csd = istart, iend
            !     nsd = nsd + cstride
            !     statei = psdi+ nsd !csd_index
            !     statef = psdf+nsd !csd_index
            !     voutp(statef) = voutp(statef) + xme*vecin(statei)
            ! end do

            statefoff = psdf - psdi;
            stateistart = csd_index + psdi + cstride;
            stateistop = stateistart + (iend - istart)*cstride
            ! really iterating over chunk of neutrons.   Input
            ! and output have the same state offset due to neutrons
            ! because they are unaffected by PP
            ! There is no way to generate overlap
  !          if(pinfo) then
  !            print *, "thread=", mythread, ", stateistart=", stateistart, ", stateistop=", stateistop
  !             call sleep(1)
  !             stop 20
  !          end if
            do statei=stateistart, stateistop, cstride
                statef = statei + statefoff;				
				localsum = localsum + xme
 !               prod = xme * vecin(statei)
!                voutp(statef) = voutp(statef) + prod
            end do
        end do  ! xjmp

  end if


     else    ! only 1 thread   ! ADDED IN VERSION 7.1.1 June 2013

  !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............

  ! KSM: get start point for neutron slater det.   This is the beginning
  ! KSM: of a diagonal block.
       csd_index = nstart(csd_index+1)-1
       cstride = 1     !   wouldn't be 1 for apply NNN
       do xjmp = xjmpstart,xjmpend
  !--------- FETCH MATRIX ELEMENT...............................
            Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
            xme = hmatpp(Xoplabel)
  !--------- GET PHASE.........................................
            xme = xme*p2b_phase(xjmp)
            psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
            psdf = p2b_2sd(xjmp)   ! KSM: final P SD
			if(psdi/=psdf)cycle
       
  !---------- GET INITIAL, FINAL SDs and place in basis..............
            ! do csd = istart, iend
            !     nsd = nsd + cstride
            !     statei = psdi+ nsd !csd_index
            !     statef = psdf+nsd !csd_index
            !     voutp(statef) = voutp(statef) + xme*vecin(statei)
            ! end do

            statefoff = psdf - psdi;
            stateistart = csd_index + psdi;
            stateistop = stateistart + (iend - istart)*cstride
            do statei=stateistart, stateistop, cstride
                statef = statei + statefoff;				
				localsum = localsum + xme
 !               voutp(statef) = voutp(statef) + xme*vecin(statei)
            end do
        end do  ! xjmp

     end if

  else    ! USE STORED PP MATRIX ELEMENTS DIRECTLY........
    if(num_threads > 1 .or. disableNoOMPloopswitch)then
    do csd = istart, iend
       csd_index = csd_index + cstride
       nsd = nstart(csd_index)

  !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
       do xjmp = xjmpstart,xjmpend
  !--------- FETCH MATRIX ELEMENT...............................
            xme = p2b_me(xjmp)
  !          Xoplabel = p2b_op(xjmp)
  !          xme = hmatpp(Xoplabel)
  !--------- GET PHASE.........................................
            xme = xme*p2b_phase(xjmp)
  !---------- GET INITIAL, FINAL SDs and place in basis..............
            statei = p2b_1sd(xjmp)+ nsd !csd_index
            statef = p2b_2sd(xjmp)+nsd !csd_index
			if(statei/=statef)cycle
			
			localsum = localsum + xme
			
!            voutp(statef) = voutp(statef) + xme*vecin(statei)
        end do  ! xjmp
     end do  ! csd

     else    ! only 1 thread   ! ADDED IN VERSION 7.1.1 June 2013
  !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
       csd_index = nstart(csd_index+1)-1
       cstride = 1
       do xjmp = xjmpstart,xjmpend
  !--------- FETCH MATRIX ELEMENT...............................
            xme = p2b_me(xjmp)
  !          Xoplabel = p2b_op(xjmp)
  !          xme = hmatpp(Xoplabel)
  !--------- GET PHASE.........................................
            xme = xme*p2b_phase(xjmp)
            psdi = p2b_1sd(xjmp)
            psdf = p2b_2sd(xjmp)
			if(psdi/=psdf)cycle
            nsd = csd_index
       
  !---------- GET INITIAL, FINAL SDs and place in basis..............
            do csd = istart, iend
                nsd = nsd + cstride
                statei = psdi+ nsd !csd_index
                statef = psdf+nsd !csd_index				
				localsum = localsum + xme
 !               voutp(statef) = voutp(statef) + xme*vecin(statei)
            end do
        end do  ! xjmp

     end if

  end if

     end if
  !$omp end parallel
  !--------------OR DO HERMITIAN/BACKWARDS APPLICATION----------
  sectortrace(localsector)=sectortrace(localsector)+localsum

    end do ! ibundle
    return
  end subroutine applyhPPtraceonly
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
subroutine applyhNNtraceonly(hchar,startbundle,endbundle )
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
	use tracy
    use contigpointervectors, only : n2b_1sd,n2b_2sd
    implicit none

    ! arguments
    character(1),intent(in) :: hchar
    integer,intent(in) :: startbundle,endbundle
    real(kind=4) :: localsum
    integer      :: localsector
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
!    real(kind=lanc_prec), pointer :: voutp(:)

 !..............................................................
 !..............................................................
 !
 !  SET UP POINTERS
 !   IF vchar = 'n' then H vec1 = vec2 
 !        (IF hchar = 'f' then multiply H_ij vec1_j = vec2_i
 !            = 'b' then multiply H_ji vec1_i = vec2_j)
 !
 !   if vchar = 'r' then H vec2 = vec1)
 !
 !         (IF hchar = 'f' then multiply H_ji vec2_i = vec1_j
 !            = 'b' then multiply H_ij vec2_j = vec1_i  )
 !  HERE i and j imply jumps between (sub) sectors
 !

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
		  
	 	 if(opbundle(ibundle)%isector/=opbundle(ibundle)%fsector)cycle   
	      localsum = 0.0
	 	 localsector = opbundle(ibundle)%isector
		 
 !...... EXTRACT INFORMATION FROM OPBUNDLE ........
          csdstart = opbundle(ibundle)%pxstart
          csdend   = opbundle(ibundle)%pxend
          cstride  = opbundle(ibundle)%cstride   !
          xjmpstart = opbundle(ibundle)%nxstart
          xjmpend   = opbundle(ibundle)%nxend
          ncstates = (csdend +cstride -csdstart)/cstride
 !--------- OUTER LOOP OVER CONJUGATE PROTON SDs---------
 !          this makes for simple OpenMP threading
 !       NOTE CSTRIDE OVER PROTON SDs 

 ! firstprivate gives each thread its own copy, but initializes it
 !    better than shared for read-only vars
 ! private gives each thread its own copy
 !$omp parallel private(vs, xjmp, Xoplabel, xme, num_threads, mythread,psd)         &
 !$omp          private(istart, iend, chunk, csd, csd_index, statef,statei)  &
 !$omp          firstprivate(cstride, ncstates, xjmpstart, xjmpend)  &
 !$omp          shared(vec2threadchunkm1) &
 !$omp          shared(n2b_op, n2b_1sd, n2b_2sd, n2b_phase, hmatnn) &
 !$omp          reduction(+:localsum)
 
          num_threads =  omp_get_num_threads()
          mythread = omp_get_thread_num()
          ! thread local vec2, reduce at end
          if(useVec2Thread) then
             !! voutp(v2s:v2e) => vec2thread(:, mythread)
             vs = mythread * vec2threadchunk
!             voutp(v2s:v2e) => vec2threadflat(vs: vs + vec2threadchunkm1)
          else
!             voutp(v2s:v2e) => vecout
          end if
          chunk = (ncstates + num_threads - 1)/num_threads
          istart = mythread*chunk + 1
          iend = bmin((mythread + 1)*chunk,ncstates)
          csd_index = csdstart + (istart - 1)*cstride - cstride
          if(istart <= iend)then

 !..... THIS FOLLOWING IS TO TRY TO FIND AN OPTIMAL ORDERING OF LOOPS
             pstride = pstridecut +1 
             if(iend-istart > 0)pstride   = pstart(csdstart+cstride)-pstart(csdstart)

 !............ THERE ARE TWO VERSIONS.....
 !             1st way: store in jumps index to NN matrix element;
 !             this has faster setup
 !             2nd way: store NN matrix elements directly;
 !             slower set up, but on MPI nodes reduced memory load
             if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO NN MATRIX ELEMENTS

                if(num_threads > 1 .or. pstride > pstridecut .or. disableNoOMPloopswitch)then
                   do csd = istart, iend
                      csd_index = csd_index + cstride
                      psd = pstart(csd_index)
 !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
                      do xjmp = xjmpstart,xjmpend
 !--------- FETCH MATRIX ELEMENT...............................
                         Xoplabel = n2b_op(xjmp)
 !		if(Xoplabel==0)print*,' ZERO LABEL ',iproc,ibundle,xjmp
                         xme = hmatnn(Xoplabel)
 !--------- GET PHASE.........................................
                         xme = xme*n2b_phase(xjmp)
 !---------- GET INITIAL, FINAL SDs and place in basis..............
                         statei = n2b_1sd(xjmp)+psd ! csd_index
                         statef = n2b_2sd(xjmp)+psd  !csd_index
						 if(statei/=statef)cycle
						 localsum = localsum + xme
   !                      voutp(statef) = voutp(statef) +  xme*vecin(statei)
                      end do  ! xjmp
                   end do  ! csd
                else  ! only 1 thread  ADDED V7.1.1 June 2013
 !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
                   csd_index = pstart(csdstart)
     
                   do xjmp = xjmpstart,xjmpend
 !--------- FETCH MATRIX ELEMENT...............................
                      Xoplabel = n2b_op(xjmp)
 !		if(Xoplabel==0)print*,' ZERO LABEL ',iproc,ibundle,xjmp
		
                      xme = hmatnn(Xoplabel)
 !--------- GET PHASE.........................................
                      xme = xme*n2b_phase(xjmp)
                      statei = n2b_1sd(xjmp) + csd_index
                      statef_start = n2b_2sd(xjmp) +csd_index
 !---------- GET INITIAL, FINAL SDs and place in basis..............
                      statef_end = statef_start+pstride*(iend-istart)
                      do statef = statef_start,statef_end,pstride
 						 if(statei==statef)localsum = localsum + xme
!                         voutp(statef) = voutp(statef) +  xme*vecin(statei)
                         statei = statei+pstride
                      end do  !
                   end do  ! xjmp
                end if
             else
                if(num_threads > 1 .or. pstride > pstridecut .or. disableNoOMPloopswitch)then
                do csd = istart, iend
                   csd_index = csd_index + cstride
                   psd = pstart(csd_index)
 !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
                   do xjmp = xjmpstart,xjmpend
 !--------- FETCH MATRIX ELEMENT...............................
                      xme = n2b_me(xjmp)
                      !  Xoplabel = n2b_op(xjmp)
                      !  xme = hmatnn(Xoplabel)
 !--------- GET PHASE.........................................
                      xme = xme*n2b_phase(xjmp)
                      !---------- GET INITIAL, FINAL SDs and place in basis..............
                      statei = n2b_1sd(xjmp)+psd ! csd_index
                      statef = n2b_2sd(xjmp)+psd  !csd_index
					  if(statei/=statef)cycle
					 localsum = localsum + xme	  
            !          voutp(statef) = voutp(statef) +  xme*vecin(statei)
                   end do  ! xjmp
                end do  ! csd
             else  ! only 1 thread  ADDED V7.1.1 June 2013
 !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
                csd_index = pstart(csdstart)
                do xjmp = xjmpstart,xjmpend
 !--------- FETCH MATRIX ELEMENT...............................
                   xme = n2b_me(xjmp)
 !--------- GET PHASE.........................................
                   xme = xme*n2b_phase(xjmp)
                   statei = n2b_1sd(xjmp) + csd_index
                   statef_start = n2b_2sd(xjmp) +csd_index
 !---------- GET INITIAL, FINAL SDs and place in basis..............
                   statef_end = statef_start+pstride*(iend-istart)
                   do statef = statef_start,statef_end,pstride
					 if(statei==statef)localsum = localsum + xme
!                      voutp(statef) = voutp(statef) +  xme*vecin(statei)
                      statei = statei+pstride
                   end do  !
                end do  ! xjmp
             end if
          end if 
       end if
 !$omp end parallel
   sectortrace(localsector)=sectortrace(localsector)+localsum
 
    end do ! ibundle
    return
end subroutine applyhNNtraceonly
