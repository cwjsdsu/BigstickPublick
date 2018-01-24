!====================================================================
!  Particle occupation report support routines
!  Update by KSM, based on code by Calvin in blanczoslib1.f90
!
!  Original code was not totaling to all the nucleons
!
!  versions for 'new' parallelization scheme -- FALL 2011
!===========================================================================
module pocc_mod
   ! single particle occupation table
   real(kind=4), allocatable :: spoccx(:,:,:) ! spoccx(species, eigenstate, orbit)

   integer :: aerr
contains


subroutine pocc_cleanup()
   implicit none

   if(allocated(spoccx)) deallocate(spoccx)
end subroutine pocc_cleanup

subroutine pocc_init_spoccx()
   use sporbit
   use io
   use lanczos_info
   use flagger
   implicit none
   integer :: maxorb

   call pocc_cleanup()
   if(spoccflag .and. .not. densityflag) then
      maxorb = MAX(numorb(1), numorb(2))
      allocate(spoccx(2, nkeep, maxorb), stat=aerr)
      if(aerr /= 0) call memerror("pocc_init_spoccx 1")
   else
      allocate(spoccx(1,1,1), stat=aerr)  ! always init
      if(aerr /= 0) call memerror("pocc_init_spoccx 2")
   end if
   spoccx = 0.0
end subroutine pocc_init_spoccx

!=============================================================
! Print out an occupation vector.
!   fn is the file number to write to
!   spoccx(pn, vec, orbit) has the occupation data
!   pn selects proton-1, neutron-2
!   lbl is either "   p occ:  ", or "  n occ:  "
subroutine pocc_write_occvec(fn, spoccx, i, pn, lbl)
   use sporbit
   use nodeinfo
   implicit none
   integer, intent(in) :: fn
   real(kind=4), intent(in)  :: spoccx(:,:,:)  ! single-particle occupations
   integer, intent(in) :: i, pn
   character(len=*), intent(in) :: lbl
   integer :: iorb
   if(iproc/=0)return
   do iorb = 1, numorb(1)
      if(MOD(iorb-1, 10) == 0) then
         if(iorb == 1) then
            write(fn, "(X,A10)", advance='no') lbl
         else
            write(fn, "(X)") ! newline
            write(fn, "(11X)", advance='no')
         end if
      else if(iorb /= 1) then
!         write(fn, "(A)", advance='no') ", "
         write(fn, "(A)", advance='no') "  "    ! putting in comma makes it more difficult to read with fortran
      end if
      write(fn, "(f7.3)", advance='no') spoccx(pn, i, iorb)
   end do
   write(fn, "(X)") ! cause advance
end subroutine pocc_write_occvec

! write out orbits in a form suitable for loading into Mathematica
subroutine pocc_write_orbits_m(mfn)
   use sporbit
   implicit none
   integer :: mfn
   integer :: i
   character(len=2) :: endb = "} ", endbc="},", endv

   write(mfn, '(A)') 'densSetOrbits[{ (*  {#, n, j*2, L } *)'
   do i=1,numorb(1)
      endv = endbc
      if(i == numorb(1)) endv = endb
      write(mfn, '(A,I4,A,I4,A,I4,A,I4,A)') "   {", i, ",", orbqn(1,i)%nr, ",", orbqn(1,i)%j, ",", orbqn(1,i)%l, endv
   end do
   write(mfn, '(A)') '}];'
end subroutine pocc_write_orbits_m

!=============================================================
! write out orbit quantum numbers with label so that
! occupation report is understandable
subroutine pocc_write_orbits_sub(fn)
   use sporbit
   implicit none
   integer :: fn ! file to write to
   integer :: i
   write(fn,*) ' Single particle state quantum numbers'
   write(fn,"(A)", advance='NO') 'ORBIT : '
   do i=1,numorb(1)
      write(fn, 100, advance='NO') i
   end do
   write(fn, *) " "  ! newline at end
   write(fn, "(A)", advance='NO') '   N : '
   do i=1,numorb(1)
      write(fn, 100, advance='NO') orbqn(1,i)%nr
   end do
   write(fn, *) " "  ! newline at end
   write(fn, "(A)", advance='NO') '   J : '
   do i=1,numorb(1)
      write(fn, 100, advance='NO') orbqn(1,i)%j
   end do
   write(fn, *) " "  ! newline at end
   write(fn, "(A)", advance='NO') '   L : '
   do i=1,numorb(1)
      write(fn, 100, advance='NO') orbqn(1,i)%l
   end do
   write(fn, *) " "  ! newline at end
   write(fn, *) " "  ! blank line
100 format(i6)
end subroutine pocc_write_orbits_sub

!=============================================================
! write out orbit quantum numbers with label so that
! occupation report is understandable
! write to both the log file and the .res file
subroutine pocc_write_orbits()
   use io
   use nodeinfo
   use flagger
   implicit none
   if(iproc == 0) then
      call pocc_write_orbits_sub(6) ! to log file
      if(writeout)  call pocc_write_orbits_sub(resultfile)
   end if
end subroutine pocc_write_orbits

subroutine pocc_write_nkeep_msg_sub(fn)
   use nodeinfo
   use lanczos_info
   implicit none
   integer :: fn

   if(iproc /= 0) return
   write(fn, '(A,i8)') 'Number of kept vectors=', nkeep
   write(fn, '(A)') ' '
end subroutine pocc_write_nkeep_msg_sub

subroutine pocc_write_nkeep_msg()
   use io
   implicit none

   call pocc_write_nkeep_msg_sub(6)
   if(writeout) call pocc_write_nkeep_msg_sub(resultfile)
end subroutine pocc_write_nkeep_msg

subroutine pocc_write_table_header()
   use nodeinfo
   use io
   implicit none
   character (len=*), parameter :: msg1 = ' State      E        Ex         J       T '
   character (len=*), parameter :: msg2 = ' State      E        Ex         '     ! added in 7.6.8

   if(iproc /= 0) return
   if(get_JT)then    !added in 7.6.8
       write(6,*) msg1
       if(writeout) write(resultfile, *) msg1
   else
       write(6,*) msg2
       if(writeout) write(resultfile, *) msg2	   
   endif
   return
end subroutine pocc_write_table_header

! write out table entry i
!
subroutine pocc_write_ejt(i, e, xj, xt)
   use nodeinfo
   use io
   implicit none
   integer :: i
   real(kind=egv_prec) :: e(*)
   real    :: xj, xt
   if(iproc == 0) then
      write(6, 100) i, e(i), e(i)-e(1), xj, xt
      if(writeout) write(resultfile, 100) i, e(i), e(i)-e(1), xj, xt
   end if
100 format(i5,3x,2f10.5,2x,2f8.3)
end subroutine pocc_write_ejt
!============================================================
!
! CALLED BY:
!   lanczos_output
!   exactdiag_p
!   particle_occupation_p
!
! SUBROUTINES CALLED:
!   makespeme
!   applyspoccbundled
!   setup4obsspe
subroutine pocc_compute_spocc(i, restorejt)
   use sporbit
   use io
   use nodeinfo
   use flagger
   use obs
   use coupledmatrixelements !interaction
   use system_parameters
   use btbme_mod
   use diagh
   use apply_obs_mod  
   implicit none
   logical, intent(in) :: restorejt
   integer, intent(in) :: i  ! which eigenvector
   integer :: iorb
   integer :: it
   
   if(spoccflag .and. .not.densityflag)then
      do iorb = 1,MAX(numorb(1),numorb(2) )
         if(numorb(1) > 0 .and. iorb <= numorb(1))then
            ! set up interaction
            pspe = 0.0
            if(numorb(2)>1)nspe = 0.0
            pspe(iorb) = 1.0
            call makespeme(1,'J')
            call makespeme(2,'J')
         end if
         if(numorb(2) > 0 .and. iorb <= numorb(2))then
            ! set up interaction
            nspe = 0.0
            if(numorb(1)>0)pspe = 0.0
            nspe(iorb) = 1.0
            call makespeme(1,'T')

            call makespeme(2,'T')
         end if
         xj2 = 0.0e0_obs_prec
         xt2 = 0.0e0_obs_prec
         call applyspoccbundled(1)
         if(np(1) > 0 .and. iorb <= numorb(1)) spoccx(1,i,iorb)= xj2   ! NOTE here xj2,xt2 are just dummy variables
         if(np(2) > 0 .and. iorb <= numorb(2)) spoccx(2,i,iorb)= xt2   ! they are NOT J2 and T2
		 
!.......... ADDED IN 7.6.0...... FOR PARTICLE-HOLE CONJUGATION..............
         do it =1,2
            if(phconj(it))then
				spoccx(it,i,iorb)= orbqn(it,iorb)%j+1-spoccx(it,i,iorb)    ! subtract off from filled
			end if
		end do  ! it		 
      end do
      if(restorejt) then
   !------------- RESTORE SINGLE PARTICLE ENERGIES FOR J2, T2 OPERATORS--------
         call setup4obsspe('J')
         call makespeme(1,'J')
         call makespeme(2,'J')
		 call setup4obsmaster('T')
      end if
   end if
end subroutine pocc_compute_spocc
!..==============================================================
!
!  CALLED BY: main routine
! 
!  SUBROUTINES CALLED:
!
!   setup_localvectors
subroutine particle_occupation_p
   use io
   use localvectors
   use obs
   use sporbit
   use lanczos_info
   use bvectorlib_mod
   use wfn_mod
   use system_parameters
   use coupledmatrixelements !interaction
	implicit none
   integer :: i,ikeep  ! which eigenvector
   real (kind=egv_prec), allocatable :: evec(:)
   real :: e, xj, xt
   character (len=*), parameter :: hdr_msg = "Generating new occupation report from saved eigenvectors (*.wfn)"
   if(iproc == 0) then
      write(6,*) hdr_msg
      if(writeout) write(resultfile,*) hdr_msg
   end if
   call setup_localvectors

   if(iproc == 0) then
      print *, "nfragments=", nfragments
   end if
   ! print *, "iproc=", iproc, ", v1s=", v1s, ", v1e=", v1e, ", v2s=", v2s, ", v2e=", v2e
   call wfn_read_nkeep(oldwfnfile, nkeep)
   call pocc_write_nkeep_msg()
   call pocc_write_orbits()
   call pocc_write_table_header()
   allocate(evec(nkeep), stat=aerr)
   if(aerr /= 0) call memerror("particle_occupation_p 1")
   if(allocated(pspe))deallocate(pspe,nspe)
   allocate(pspe(numorb(1)), nspe(numorb(2)), stat=aerr)
   if(aerr /= 0) call memerror("particle_occupation_p 2")
   call pocc_init_spoccx()
   twoobsflag = .true.  ! compute both two observables at same time

   ! outchoice = 'd'
   ! if(np(1)>0 .and. np(2) > 0)outchoice='b'
   ! if(np(1)>0 .and. np(2) == 0)outchoice='p'
   ! if(np(1)==0 .and. np(2) > 0)outchoice='n'
   do ikeep=1,nkeep
      ! note: write and read actually store xt^2, have to convert
      ! new interface, we say which vec to read and it checks
      i = ikeep
      call wfn_readeigenvec(oldwfnfile, frag1, fcomm1, vec1, i, e, xj, xt) ! KSM: updated
      evec(ikeep) = e
      xt2 = xt
      xt = real(-0.5 + sqrt(xt2 + 0.25))
      call pocc_compute_spocc(i, .false.)
      if(iproc /= 0) cycle
      ! now we have to figure out what to print 
      call pocc_write_ejt(i, evec, xj, xt)
      if(npeff(1) > 0) call pocc_write_occvec(6, spoccx, i, 1, "   p_occ:")
      if(npeff(2) > 0) call pocc_write_occvec(6, spoccx, i, 2, "   n_occ:")
      if(writeout) then
         if(npeff(1) > 0) call pocc_write_occvec(resultfile, spoccx, i, 1, "   p_occ:")
         if(npeff(2) > 0) call pocc_write_occvec(resultfile, spoccx, i, 2, "   n_occ:")
      end if
   end do
   ! add blank line for formatting
   if(iproc == 0) then
      write(6,*) " "
      if(writeout) write(resultfile, *) " "
   end if
   deallocate(evec)
   deallocate(pspe)
   deallocate(nspe)
   call pocc_cleanup()
end subroutine particle_occupation_p

!=============================================================
!  subroutine  particle_occupation_p_orig
!  based upon subroutine expectator_p
!       master subroutine to control calculation of particle occupations
!
!  NOTE: There will have to be modifications when sections of the lanczos vector are broken into pieces
!
! routines called:
!  setup_localvectors
!  pocc_write_orbits
!   wfn_read_nkeep
!  pocc_write_nkeep_msg()
!   pocc_write_table_header
!

subroutine particle_occupation_p_orig

  use sporbit
  use precisions
  use basis
  use io
  use localvectors
  use obs
  use fragments
  use nodeinfo
  use coupledmatrixelements ! interaction
  use wfn_mod
  use bvectorlib_mod
  use diagh

  implicit none
!  include 'binterfaces.inc'
  
  integer nkeep
  integer ikeep
  integer i,j,n
  real :: e, xj, xt
  real (kind=egv_prec), allocatable :: evec(:)
  real (kind = 8) :: dnorm
  logical zeroflag
  integer ierr
  real (kind = lanc_prec), allocatable :: exph(:)
   integer it
   integer iorb
  character (len=*), parameter :: hdr_msg = "Generating occupation report from saved eigenvectors (*.wfn)"

  if(iproc == 0) then
	  write(6,*) hdr_msg
	  if(writeout) write(resultfile,*) hdr_msg
  end if

!........... TEMP...... WILL NEED TO FIX LATER....
  basestart(1) = 1
  basestop(1) = dimbasis
  frag1 = 1
  frag2 = 1
  call setup_localvectors

  call pocc_write_orbits()

  call wfn_read_nkeep(oldwfnfile, nkeep)

  allocate(evec(nkeep), stat=aerr)
  if(aerr /= 0) call memerror("particle_occupation_p_orig 1")

  call pocc_write_nkeep_msg()
  call pocc_write_table_header()

  twoobsflag = .false.  ! may change
  !print *, "iproc=", iproc, ", Reading eigenvectors, nkeep=", nkeep
  do ikeep = 1,nkeep
     ! new interface, we say which vec to read and it checks
     i = ikeep
     call wfn_readeigenvec(oldwfnfile, frag1, fcomm1, vec1, i, e, xj, xt) ! KSM: updated
     evec(ikeep) = e
     xt2 = xt
     xt = real(-0.5 + sqrt(xt2 + 0.25))
     ! print *, "iproc=", iproc, ", read ev, xj=", xj, ", xt=", xt
     call dnormvec_p('n','i',dnorm,zeroflag)
     if(allocated(pspe))deallocate(pspe,nspe)
     allocate(pspe(numorb(1)), nspe(numorb(2)), stat=aerr)
     if(aerr /= 0) call memerror("particle_occupation_p_orig 2")
     do it = 1,2

        allocate(exph( numorb(it)), stat=aerr)
        if(aerr /= 0) then
           call memerror("particle_occupation_p_orig 3")
           stop 5 ! for -Wuninitialized
        end if
        exph = 0.0
        do iorb = 1,numorb(it)
!--------------------- SET UP ARRAY
        pspe = 0.0
        nspe = 0.0
        if(it == 1)pspe(iorb) = 1.0
        if(it == 2)nspe(iorb) = 1.0
        call makespeme(1,'H')
        call makespeme(2,'H')

!---------------------- COMPUTE OCCUPATION
        call applyspes_expect(vec1,exph(iorb))
     enddo  ! iorb

!------------------ WRITE OUT
     if(iproc==0)then

       if(it == 1)then
          call pocc_write_ejt(i, evec, xj, xt)

          write(6,101)(exph(iorb),iorb = 1,numorb(it))
          if(writeout)write(resultfile,101)(exph(iorb),iorb = 1,numorb(it))
       else
          write(6,102)(exph(iorb),iorb = 1,numorb(it))
          if(writeout)write(resultfile,102)(exph(iorb),iorb = 1,numorb(it))
       endif

    end if
101 format(' p : ',10f7.3)
102 format(' n : ',10f7.3)
1001  format(i5,2x,f11.4,2x,f5.1,2x,f5.1)

     deallocate(exph)
   enddo  ! it
  end do ! ikeep
  deallocate(evec)
  
  return
end subroutine particle_occupation_p_orig

!==================================================================

!
!  subroutine applyspes_obs
!
!  applies single-particle energies -- purely diagonal
!

subroutine applyspes_expect(vecin,exph)

  use basis
  use sectors
  use diagh
  use precisions
  use nodeinfo
  use localvectors
  use fragments
  use bmpi_mod
  implicit none

  real (kind = lanc_prec) :: vecin(dimbasis)
  real (kind = lanc_prec) :: v
  real (kind = lanc_prec) :: exph
  integer :: is,isc,jsc
  integer(kind=basis_prec) :: ip, in, ibasis
  real pspe,nspe
  integer:: ierr
  
!---------- LOOP OVER PROTON SECTORS-----------

  do is = 1,nsectors(1)
!................ LOOP OVER PROTON SDs in that sector

     do ip = xsd(1)%sector(is)%xsdstart,xsd(1)%sector(is)%xsdend
        pspe = pspe_h(ip)   ! the proton contribution to the single-particle energies
!-------- LOOP OVER NEUTRON SECTORS conjugate to is
        do isc = 1,xsd(1)%sector(is)%ncsectors
           jsc = xsd(1)%sector(is)%csector(isc)
!............... LOOP OVER NEUTRON SDS in sector jsc
           do in = xsd(2)%sector(jsc)%xsdstart,xsd(2)%sector(jsc)%xsdend
              nspe = nspe_h(in)  ! neutron contributions to s.p.e.s
              ibasis = pstart(ip)+nstart(in)  ! find basis state
			  if(ibasis >= v1s .and. ibasis <= v1e) then
				  v = vecin(ibasis)
				  exph = exph + v*v*(pspe+nspe)
			  end if
           enddo
        enddo !isc
     enddo  !ip

  enddo  !is

  if(.not. isfragroot) exph = 0.0
  call BMPI_ALLREDUCE(exph, 1,  MPI_SUM, icomm, ierr) ! in place

  return
end subroutine applyspes_expect
!======================================================================
end module pocc_mod
