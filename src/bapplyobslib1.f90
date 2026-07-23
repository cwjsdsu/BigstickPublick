!========================================================
!
!  routines for computing expectation values and observables
!  such as J^2, T^2, in the "new" parallelization scheme
!
!=======================================================
module apply_obs_mod
contains
	
!
!============================================================

!   subroutine applicator_h
!       master subroutine to control application of a scalar Hamiltonian-like operator
!
! subroutines called:
!  setup_localvectors
!  wfn_read_nkeep
!  wfn_write_nkeep
!  wfn_readeigenvec
!  initialize_final
!  applyhbundled_g(
!
   subroutine applicator_h

   use precisions
   use basis
   use io
   use localvectors
   use fragments
   use nodeinfo
   use opmatrixelements, only: subtract_enabled
   use wfn_mod
   use bmpi_mod
   use bvectorlib_mod
!   use lanczos_util
   use apply_ham
   implicit none
!   include 'binterfaces.inc'

   integer nkeep
   integer ikeep
   integer i,j,n
   integer(kind=basis_prec) :: k
   real :: e,xj,xt2
   real (kind = obs_prec) :: dot  !, exph
   real (kind = 8) :: dot8, dot8a
   integer ierr
   integer :: aerr

   call setup_localvectors

   call wfn_read_nkeep(oldwfnfile, nkeep)  ! does BCAST inside
   if(iproc==0) print *, nkeep, ' states '
   call wfn_write_nkeep(nkeep) ! write number of vectors to wfn file
   call   procOP_clock(0,'set','all')
   do ikeep = 1,nkeep
      !MKGK edit
      if (ikeep==1 .and. iproc==0 .and. subtract_enabled) then
           write(*,*) '** Subtracting elastic peak - bdenslib2 **'
      end if

      i = ikeep
      call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,e,xj,xt2) 
      if(iproc==0) print*, i, e, xj,xt2
      call initialize_final('n')
      call applyhbundled_g('n')

      ! subtract elastic peak
      if(subtract_enabled)then
	      if(iproc == 0) then
		      print *, "CALVIN: applicator_h with subtrack_enabled"
		   end if
         dot = 0.0
         if(isfragroot) then
            ! each fragment has one selected root node where frag1==frag2
            do k=v1s,v1e
              dot = dot + vec1(k)*vec2(k)
            end do
            ! perform reduction if we are using fragments
            dot8 = dot;
         end if
		   ! next step has to run on all nodes
         if(nfragments > 1) then ! this nfragments reference is ok
            dot8a = dot8
#ifdef _MPI
            call BMPI_ALLREDUCE(dot8a, dot8, 1, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
		   end if
		   ! vec2 and vec1 have same format on root nodes
         if(isfragroot) vec2(:) = vec2(:) - real(dot8*vec1(:), kind=lanc_prec)
      end if
	   ! write only uses fragment root nodes
      if(writeout)call wfn_writeeigenvec(wfnfile,frag2, vec2, i,e,xj,xt2)
   end do ! ikeep

   return
   end subroutine applicator_h   

!===================================================================
!
! subroutine applyobsbundled
!
! note: default is going from vecin to vecout but this can be reversed depending on hchar
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)
!   ivec : 1 = vec1, 2 = vec2

subroutine applyobsbundled(ivec)
  use nodeinfo
  use precisions
  use opbundles
  use fragments
!  use shampoo
  use interaction
  use obs
  use bmpi_mod
  implicit none

!  include 'binterfaces.inc'

  integer :: ibundle,ivec
  integer :: ierr
  integer :: procstart,procstop,iprocs
  call clocker('obs','sta')
!........ OPTION TO SIMULATE MPI ON A SINGLE "CORE".....

  if(distributeMPI .and. nproc == 1)then
      procstart = 0
      procstop  = nprocs -1
  else
      procstart = iproc
      procstop  = iproc
  end if

  do iprocs = procstart,procstop

  do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
!	  if(opbundle(ibundle)%annexed)cycle

  select case ( opbundle(ibundle)%optype )

    case ('PP')
        if( opbundle(ibundle)%hchar == 'f' )then
            call applyobsPPbundled(ibundle,ivec)

        end if
    case ('NN')

        if( opbundle(ibundle)%hchar == 'f' )then

           call applyobsNNbundled(ibundle,ivec)

        endif
    case ('SPE')

        call applyobsSPEbundled(ibundle,ivec)

    case ('PN') 
        if( opbundle(ibundle)%hchar /= 'b'  )then
           if(useTR)then

           call applyobsPNbundledTR(ibundle,ivec)
           else

           call applyobsPNbundled(ibundle,ivec)

           end if
        end if
    case ('PPP')

    case ('PPN')

    case ('PNN')

    case ('NNN')


  end select

  end do  !ibundle
  end do
!--------- REDUCE-----
#ifdef _MPI
  if(nproc > 1)then
     call BMPI_Allreduce(xj2,xj2tmp,1,MPI_SUM,MPI_COMM_WORLD,ierr)
     call BMPI_Allreduce(xt2,xt2tmp,1,MPI_SUM,MPI_COMM_WORLD,ierr)
     xj2 = xj2tmp
     xt2 = xt2tmp
  end if
#endif
  call clocker('obs','end')

  return

end subroutine applyobsbundled

!===================================================================
!
!
! subroutine applyspoccbundled   ADDED IN 7.3.7
!    based upon applyobsbundled
! note: default is going from vecin to vecout but this can be reversed depending on hchar
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)
!   ivec : 1 = vec1, 2 = vec2

subroutine applyspoccbundled(ivec)
  use nodeinfo
  use precisions
  use opbundles
  use fragments
!  use shampoo
  use interaction
  use obs
  use bmpi_mod
  implicit none

!  include 'binterfaces.inc'

  integer :: ibundle,ivec
  integer :: ierr
  integer :: procstart,procstop,iprocs
  call clocker('obs','sta')
!........ OPTION TO SIMULATE MPI ON A SINGLE "CORE".....

  if(distributeMPI .and. nproc == 1)then
      procstart = 0
      procstop  = nprocs -1
  else
      procstart = iproc
      procstop  = iproc
  end if

  do iprocs = procstart,procstop

  do ibundle = opbundlestart(iprocs), opbundleend(iprocs)

  select case ( opbundle(ibundle)%optype )

    case ('SPE')

        call applyobsSPEbundled(ibundle,ivec)

  end select

  end do  !ibundle
  end do
!--------- REDUCE-----

#ifdef _MPI
  if(nproc > 1)then
     call BMPI_Allreduce(xj2,xj2tmp,1,MPI_SUM,MPI_COMM_WORLD,ierr)   ! note here xj2 and xt2 are just dummy variables
     call BMPI_Allreduce(xt2,xt2tmp,1,MPI_SUM,MPI_COMM_WORLD,ierr)
     xj2 = xj2tmp
     xt2 = xt2tmp
  end if
#endif
  call clocker('obs','end')

  return

end subroutine applyspoccbundled

!===================================================================
!  subroutine applyobsPPbundled
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)

!
!===================================================================
subroutine applyobsPPbundled ( ibundle,ivec )

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
  use obs
  use basis
  use bmpi_mod
  use contigpointervectors,only : vecin,vecout 
  implicit none

  integer :: ibundle
! --------- NOTE: basestart, basestop stored in module fragments

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,nsd
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=obs_prec)   xme, xtme
  real(kind=obs_prec) :: vmult,vmult1,vmult2
  integer hfactor  ! to take care of HERMITICITY
  integer(kind=basis_prec) :: state1,state2
  integer(kind=4) num_threads
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  real, pointer :: otherobs(:)
  integer ivec
!  real(kind=lanc_prec), pointer, contiguous :: veci(:), vecf(:)

  if(ivec ==1)then
     vecin => vec1
     vecout => vec2
  else
     vecin => vec2
     vecout => vec1
  endif
!...... EXTRACT INFORMATION FROM OPBUNDLE ........

  csdstart = opbundle(ibundle)%nxstart
  csdend   = opbundle(ibundle)%nxend
  xjmpstart = opbundle(ibundle)%pxstart
  xjmpend   = opbundle(ibundle)%pxend
!  cstride   = opbundle(ibundle)%cstride
  hfactor = 2   ! always for this
  if(twoobsflag)then
     otherobs => obsmatpp
  else
     otherobs => hmatpp
  end if
!--------- OUTER LOOP OVER CONJUGATE NEUTRON SDs---------
!          this makes for simple OpenMP threading

!$omp parallel do reduction(+:xj2,xt2) &
!$omp private(csd,xjmp,xoplabel,xme,xtme,state1,state2,vmult,nsd) &
!$omp shared(csdstart,csdend,hmatpp,otherobs,hfactor) & 
!$omp shared(p2b_op,p2b_phase,p2b_isd,p2b_fsd,vec1,vec2)
  do csd = csdstart,csdend
     nsd = nstart(csd)
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
        Xoplabel = p2b_op(xjmp)
        xme = real(hmatpp(Xoplabel),obs_prec)
        xtme= real(otherobs(Xoplabel),obs_prec)

!--------- GET PHASE.........................................
        xme = xme*p2b_phase(xjmp)
        xtme = xtme*p2b_phase(xjmp)

!---------- GET INITIAL, FINAL SDs and place in basis..............
        state1 = p2b_isd(xjmp)+nsd
        state2 = p2b_fsd(xjmp)+nsd

        vmult = real(vecin(state1),obs_prec)*real(vecout(state2),obs_prec)
!        print*,vmult,xme,hfactor,xj2
        xj2 = xj2 + xme*vmult*hfactor
        xt2 = xt2 + xtme*vmult*hfactor

        end do  ! xjmp
   end do  ! csd
!$omp end parallel do
!--------------OR DO HERMITIAN/BACKWARDS APPLICATION----------

  return
end subroutine applyobsPPbundled
!===================================================================
!  subroutine applyobsNNbundled
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., NN between two (sub) sectors, etc)
!   
!
!===================================================================

subroutine applyobsNNbundled ( ibundle,ivec )
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use interaction
  use opbundles
  use fragments
  use obs
  use basis
  use bmpi_mod
  use contigpointervectors,only : vecin,vecout 
 
  implicit none

  integer :: ibundle
!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,psd
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=obs_prec)   xme, xtme
  real(kind=obs_prec) :: vmult
  integer hfactor  ! to take care of HERMITICITY
  integer(kind=8), pointer :: statei, statef
  integer(kind=8), target :: state1,state2
  integer(kind=4) num_threads
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  real, pointer :: otherobs(:)
  integer ivec
!  real(kind=lanc_prec), pointer, contiguous :: veci(:), vecf(:)

  if(ivec ==1)then
     vecin => vec1
     vecout => vec2
  else
     vecin => vec2
     vecout => vec1
  endif

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
  csdstart = opbundle(ibundle)%pxstart
  csdend   = opbundle(ibundle)%pxend
!  cstride  = opbundle(ibundle)%cstride   !
  xjmpstart = opbundle(ibundle)%nxstart
  xjmpend   = opbundle(ibundle)%nxend

  hfactor =2
  if(twoobsflag)then
     otherobs => obsmatnn
  else
     otherobs => hmatnn
  end if
!--------- OUTER LOOP OVER CONJUGATE PROTON SDs---------
!          this makes for simple OpenMP threading
!       NOTE CSTRIDE OVER PROTON SDs 
!$omp parallel do reduction(+:xj2,xt2) &
!$omp private(csd,xjmp,xoplabel,xme,xtme,state1,state2,vmult,psd) &
!$omp shared(csdstart,csdend,hmatnn,otherobs,hfactor) & 
!$omp shared(n2b_op,n2b_phase,n2b_isd,n2b_fsd,vec1,vec2)
  do csd = csdstart,csdend
     psd = pstart(csd)
!--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
     do xjmp = xjmpstart,xjmpend
!--------- FETCH MATRIX ELEMENT...............................
        Xoplabel = n2b_op(xjmp)
        xme = real(hmatnn(Xoplabel),obs_prec)
        xtme= real(otherobs(Xoplabel),obs_prec)

!--------- GET PHASE.........................................
        xme = xme*n2b_phase(xjmp)
        xtme = xtme*n2b_phase(xjmp)

!---------- GET INITIAL, FINAL SDs and place in basis..............
        state1 = n2b_isd(xjmp)+psd
        state2 = n2b_fsd(xjmp)+psd

        vmult = real(vecin(state1),obs_prec)*real(vecout(state2),obs_prec)

        xj2 = xj2 + xme*vmult*hfactor
        xt2 = xt2 + xtme*vmult*hfactor
        end do  ! xjmp
   end do  ! csd


  return
end subroutine applyobsNNbundled
!==================================================================

subroutine applyobsPNbundled( ibundle,ivec )
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use interaction
  use opbundles
  use fragments
  use obs
  use bmpi_mod
  use contigpointervectors,only : vecin,vecout 
  
  implicit none

  integer :: ibundle


!------------------------------------------------------------
  integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer a,b,c,d
  integer(kind=8) coplabel,doplabel
  integer :: phasep,phasen
  real(kind=obs_prec)   xme, xtme
  real(kind=obs_prec) :: vmult
  integer hfactor  ! to take care of HERMITICITY
  integer(kind=basis_prec) :: statei, statef
  integer(kind=4) num_threads
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  real, pointer :: otherobs(:)
  integer ivec
!  real(kind=lanc_prec), pointer, contiguous :: veci(:), vecf(:)

  if(ivec ==1)then
     vecin => vec1
     vecout=> vec2
  else
     vecin => vec2
     vecout=> vec1
  endif

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

  if(opbundle(ibundle)%hchar == 'h')then
      hfactor = 1
  else
      hfactor = 2
  end if

  if(twoobsflag)then
     otherobs => obsmatpn
  else
     otherobs => hmatpn
  end if

!$omp parallel do reduction(+:xj2,xt2) &
!$omp private(pjmp,njmp,psdi,psdf,nsdi,nsdf,statei,statef,vmult) &
!$omp private(a,b,c,d,coplabel,doplabel,phasep,phasen,xme,xtme) &
!$omp shared(pjmpstart,pjmpend,hmatpn,otherobs,hfactor) & 
!$omp shared(n1b_cop,n1b_dop,n1b_phase,n1b_isd,n1b_fsd,cpnpair,dpnpair) & 
!$omp shared(p1b_cop,p1b_dop,p1b_phase,p1b_isd,p1b_fsd,vec1,vec2)
!--------- LOOP OVER PROTON JUMPS................................
  do pjmp = pjmpstart,pjmpend
     psdi = p1b_isd(pjmp)       ! initial proton slater determinant
     psdf = p1b_fsd(pjmp)       ! final proton SD
     phasep = p1b_phase(pjmp)   ! phase of proton jumps

     a    = p1b_cop(pjmp) 
     c    = p1b_dop(pjmp)
  
!--------- LOOP OVER NEUTRON JUMPS
        do njmp = njmpstart,njmpend
      
!----------- FIND MATRIX ELEMTN
           b    = n1b_cop(njmp)
           d    = n1b_dop(njmp)

           coplabel = cpnpair(b,a)
           doplabel = dpnpair(d,c)
           phasen = n1b_phase(njmp)

           xme = real(hmatpn(coplabel+doplabel),obs_prec)     ! get matrix element
           xtme= real(otherobs(coplabel+doplabel),obs_prec)
           xme = xme*phasep*phasen    ! multiply matrix element by jump phases
           xtme = xtme*phasep*phasen  ! multiply matrix element by jump phases

           nsdi = n1b_isd(njmp)
           nsdf = n1b_fsd(njmp)
    
           statei = nsdi+psdi     ! initial state in combined basis
           statef = nsdf+psdf     ! final state in combined basis
           vmult = real(vecin(statei),obs_prec)*real(vecout(statef),obs_prec)
           xj2 = xj2 + xme*vmult*hfactor
           xt2 = xt2 + xtme*vmult*hfactor

        end do  ! njmp
  end do  ! pjmp

  return
end subroutine applyobsPNbundled
!==================================================================
subroutine applyobsPNbundledTR( ibundle,ivec )
  use localvectors
  use nodeinfo
  use system_parameters
  use jumpNbody
  use precisions
  use interaction
  use opbundles
  use fragments
  use obs
  use bmpi_mod
  use contigpointervectors,only : vecin,vecout 
  
  implicit none

  integer :: ibundle


!------------------------------------------------------------
  integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

  integer(kind=8) pjmp,pjmpstart,pjmpend
  integer(kind=8) njmp,njmpstart,njmpend
  integer a,b,c,d
  integer(kind=8) :: coplabel,doplabel
  integer :: phasep,phasen
  real(kind=obs_prec)   xme, xtme
  real(kind=obs_prec) :: vmult
  integer hfactor  ! to take care of HERMITICITY
  integer(kind=basis_prec) :: statei, statef
  integer(kind=4) num_threads
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  real, pointer :: otherobs(:)
  integer ivec
!  real(kind=lanc_prec), pointer, contiguous :: veci(:), vecf(:)

  if(ivec ==1)then
     vecin => vec1
     vecout => vec2
  else
     vecin => vec2
     vecout => vec1
  endif

!...... EXTRACT INFORMATION FROM OPBUNDLE ........
!
  pjmpstart = opbundle(ibundle)%pxstart
  pjmpend   = opbundle(ibundle)%pxend
  njmpstart = opbundle(ibundle)%nxstart
  njmpend   = opbundle(ibundle)%nxend

  if(opbundle(ibundle)%hchar == 'h')then
      hfactor = 1
  else
      hfactor = 2
  end if

  if(twoobsflag)then
     otherobs => obsmatpn
  else
     otherobs => hmatpn
  end if

!$omp parallel do reduction(+:xj2,xt2) &
!$omp private(pjmp,njmp,psdi,psdf,nsdi,nsdf,statei,statef,vmult) &
!$omp private(a,b,c,d,coplabel,doplabel,phasep,phasen,xme,xtme) &
!$omp shared(pjmpstart,pjmpend,hmatpn,otherobs,hfactor) & 
!$omp shared(n1b_cop,n1b_dop,n1b_phase,n1b_isd,n1b_fsd,cpnpair,dpnpair) & 
!$omp shared(p1b_cop,p1b_dop,p1b_phase,p1b_isd,p1b_fsd,vec1,vec2)
!--------- LOOP OVER PROTON JUMPS................................
  do pjmp = pjmpstart,pjmpend
     psdi = p1b_isd(pjmp)       ! initial proton slater determinant
     psdf = p1b_fsd(pjmp)       ! final proton SD
     phasep = p1b_phase(pjmp)   ! phase of proton jumps

     a    = p1b_cop(pjmp) 
     c    = p1b_dop(pjmp)
  
!--------- LOOP OVER NEUTRON JUMPS
        do njmp = njmpstart,njmpend
      
!----------- FIND MATRIX ELEMTN
           b    = n1b_cop(njmp)
           d    = n1b_dop(njmp)

           coplabel = cpnpair(b,a)
           doplabel = dpnpair(d,c)
           phasen = n1b_phase(njmp)*phasepnpair(b,a)*phasepnpair(d,c)

           xme = real(hmatpn(coplabel+doplabel),obs_prec)     ! get matrix element
           xtme= real(otherobs(coplabel+doplabel),obs_prec)
           xme = xme*phasep*phasen    ! multiply matrix element by jump phases
           xtme = xtme*phasep*phasen  ! multiply matrix element by jump phases

           nsdi = n1b_isd(njmp)
           nsdf = n1b_fsd(njmp)
    
           statei = nsdi+psdi     ! initial state in combined basis
           statef = nsdf+psdf     ! final state in combined basis
           vmult = real(vecin(statei),obs_prec)*real(vecout(statef),obs_prec)
           xj2 = xj2 + xme*vmult*hfactor
           xt2 = xt2 + xtme*vmult*hfactor

        end do  ! njmp
  end do  ! pjmp

  return
end subroutine applyobsPNbundledTR
!==================================================================
!
!  subroutine applyspes
!
!  applies single-particle energies -- purely diagonal
!
!====================================================================
subroutine applyobsSPEbundled(ibundle,ivec)
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
  use obs
  use bmpi_mod
  use contigpointervectors,only : vecin
  
  implicit none

  integer :: ibundle

!------------------------------------------------------------

  integer(kind=8) nsdstart,nsdend,psdstart,psdend
  integer(kind=8) xjmp,xjmpstart,xjmpend

  integer(kind=4) num_threads
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  real(kind=obs_prec) :: v
  real(kind=obs_prec) pspe,nspe,pT2spe,nT2spe
  integer(kind=basis_prec)  :: ip,in, ibasis
  real(4), pointer :: otherpspe(:), othernspe(:)
  integer ivec
!  real(kind=lanc_prec), pointer, contiguous :: vec(:)

  if(ivec ==1)then
     vecin => vec1
  else
     vecin => vec2
  endif
  psdstart = opbundle(ibundle)%pxstart
  psdend   = opbundle(ibundle)%pxend
  nsdstart = opbundle(ibundle)%nxstart
  nsdend   = opbundle(ibundle)%nxend

  if(twoobsflag)then
    otherpspe =>  pspe_obs
    othernspe =>  nspe_obs
  else
    otherpspe =>  pspe_h
    othernspe =>  nspe_h
  end if

!omp parallel do reduction(+:xj2,xt2) & 
!omp private(ip,in,pspe,nspe,ibasis,pT2spe,nT2spe,v) & 
!omp  shared(psdstart,psdend,nsdstart,nsdend,vec1) & 
!omp  shared(pstart,nstart,pspe_h,nspe_h,otherpspe,othernspe)
!................ LOOP OVER PROTON SDs in that sector..........
  do ip = psdstart,psdend
           pspe = real(pspe_h(ip),obs_prec)   ! the proton contribution to the single-particle energies
           pT2spe = real(otherpspe(ip),obs_prec)   ! the proton contribution to the single-particle energies
!............... LOOP OVER NEUTRON SDS in sector jsc........................
           do in = nsdstart,nsdend
                 nspe = real(nspe_h(in),obs_prec)  ! neutron contributions to s.p.e.
                 nT2spe = real(othernspe(in),obs_prec)   ! the proton contribution to the single-particle energies
                 ibasis = pstart(ip) + nstart(in)  ! find basis state
                 v = real(vecin(ibasis),obs_prec)
                 xj2 = xj2 + v*v*(pspe+nspe)
                 xt2 = xt2 + v*v*(pT2spe+nT2spe)
           end do  !in
  end do  !ip

  return
end subroutine applyobsSPEbundled

end module apply_obs_mod
