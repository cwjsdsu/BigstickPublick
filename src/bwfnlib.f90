! NOTE: reorganized in 7.6.6 by CWJ to group like subroutines together
!
!  modules with global data and subroutines to open/close/read/write wfn files
!

module wfn_mod
   use io
   use nodeinfo
   use fragments
   use bmpi_mod
   implicit none

   logical :: wfn_mpiio

   integer, parameter :: wfn_magic = 31415926
   integer, parameter :: wfn_vecmagic = 27182818
   integer, parameter :: wfn_version = 1
   integer, parameter :: wfn_after_header_word = 3
   ! hack - kind's do not have to be the number of 
   ! bytes.  They are symbolic labels for options
   ! could use bit_size to test.
   integer, parameter :: wfn_integer_kind = 4

   ! remember locations in file for later update
   
#ifdef _MPI   
   integer(kind=MPI_OFFSET_KIND) :: wfn_after_header_pos  ! location just after header
   type (MPI_file) :: mpifile_in, mpifile_out
#else
     integer :: wfn_after_header_pos  ! location just after header

#endif
	   
contains
	
! CONTAINS THE FOLLOWING SUBROUTINES
!
!
!
!	

!-------------------- WAVEFUNCTION OPEN/CLOSE/READ/WRITE SUBROUTINES---------------
!
! Check if we should be doing MPI IO
subroutine wfn_checkmpiio()
   use localvectors
   implicit none
   if(iproc == 0) return ! ok
   if(isfragroot) return ! ok
   print *, "iproc=", iproc, "Should not be doing MPI IO"
   stop 1 ! has msg
end subroutine wfn_checkmpiio

!============== MASTER ROUTINES FOR WRITING WFN FILE(S) =================
!-------------------------------------------------------------
!   subroutine wrfn_writeeigenvec (formerly: writeeigenvec_p)
!
!  writes eigenvector to disk
!
!  INPUT
!       filenumber = which file to write to
!       fg = fragment held in this processor, frag1 or frag2
!      v(:) either (v1s:v1e) or (v2s:v2e)
!       i  = index of eigenvector
!       e  = energy of eigenvector
!       xj = j (as real number) of eigenvector
!       xt2 = T(T+1) [as real number] of eigenvector
!
!   CALLS:
!     wfn_seekvecfragpos
!   CALLED BY:
!    exactdiag_p, exactdiag_MPI
!    lanczos_output
!    applicator1b
!    applicator_h
!   
subroutine wfn_writeeigenvec(filenumber, fg, v, i,e,xj,xt2)  ! formerly writeeigenvec_p
   use lanczos_info
   use basis
   use precisions
   use io
   use nodeinfo
   use localvectors
   use fragments
   implicit none

   character(1) :: vchar
   real(kind=lanc_prec), intent(in) :: v(:) ! either (v1s:v1e) or (v2s:v2e)
   real :: e,xj,xt2
   integer(4):: i,j
   integer :: filenumber
   integer :: fg
   integer:: ierr
   integer(kind=basis_prec) :: nw
   integer(kind=basis_prec) :: jl
#ifdef _MPI
   integer (kind=MPI_OFFSET_KIND) :: offset
   type(MPI_file) :: my_mpi_file
   type(MPI_status) :: status
#else
    integer :: offset

#endif   

   if(.not.writeout)return   
   if(wfn_mpiio) then
#ifdef _MPI	   
      if(.not. isfragroot) return
	  my_mpi_file = mpifile_out
!---------------- move to start of vec/fragment
      call wfn_seekvecfragpos(filenumber, i, fg)
      call BMPI_FILE_get_position(my_mpi_file, offset, ierr)
    	! start writing
      if(fg == 1) then
         ! first fragment writes the overhead
         call BMPI_FILE_WRITE(my_mpi_file, wfn_vecmagic, 1, status, ierr)
         call BMPI_FILE_WRITE(my_mpi_file, i, 1,status, ierr)
         call BMPI_FILE_WRITE(my_mpi_file, e, 1, status, ierr)
         call BMPI_FILE_WRITE(my_mpi_file, xj, 1, status, ierr)
         call BMPI_FILE_WRITE(my_mpi_file, xt2, 1,  status, ierr)
      end if
      nw = basestop(fg)- basestart(fg) + 1;
      call BMPI_FILE_WRITE(my_mpi_file, v, nw, status, ierr)
      call BMPI_FILE_GET_POSITION(my_mpi_file, offset, ierr)
#endif	  
   else 
      if(iproc/=0) return
      inquire(filenumber, pos=offset)
      write(filenumber) wfn_vecmagic,i,e,xj,xt2
      write(filenumber) (v(jl),jl=1,dimbasis)
      inquire(filenumber, pos=offset)
   end if
   return
end subroutine wfn_writeeigenvec  ! formerly:  writeeigenvec_p

!=========================================================
! ====  MASTER ROUTINE(S) FOR READING WFN FILE(S) ========
!=========================================================
!   subroutine wfn_readeigenvec
!
!  reads eigenvector from disk
!
!  always reads into a vector formatted like vec1.
!
!  INPUT:
!       filenumber = which unit to fetch vector from
!       fg = fragment held in this processor, frag1 or frag2
!       fc = communicator for update of other nodes belonging to fragment fg
!            either fcomm1 or fcomm2
!       i  = index of eigenvector - checks that we get the right one
!  OUTPUT
!    v  ! Either v1s:v1e, or v2s:v2e
!       e  = energy of eigenvector
!       xj = j (as real number) of eigenvector
!
!  SUBROUTINES CALLED:
!   wfn_seekvecfragpos
!
subroutine wfn_readeigenvec(filenumber, fg, fc_index, v, i, e, xj, xt2)
	use lanczos_info
	use basis
	use precisions
	use io
	use nodeinfo
	use localvectors
   use mod_reorthog
   use bmpi_mod
	implicit none
	! arguments
	integer, intent(in) :: filenumber
	integer, intent(in) :: fc_index ! communicator for vec1/2 slicing
	integer, intent(in) :: fg ! fragment this node works on
	integer, intent(in) :: i
	real, intent(out) :: e,xj,xt2
   ! we have an interface, so should pick up bounds of vector v
   real (kind=lanc_prec) :: v(:)  ! Either v1s:v1e, or v2s:v2e
	! locals
	real(kind=4) :: ejtbuf(3) ! read/bcast buffer for e, xj, xt2
	integer(4) :: j
	integer(kind=basis_prec) :: jl
	integer :: ierr
!	integer :: status(MPI_STATUS_SIZE)
   integer(4) :: vmagic, vnum
   integer :: rootproc
   character(len=*), parameter :: msg_wrongvec = "wfn_readeigenvec: read wrong vector"
#ifdef _MPI
      type(MPI_file) :: my_mpi_file
      type(MPI_status) :: status
	  type(MPI_comm) :: fc
		  	
	  my_mpi_file = mpifile_in	  
	  
	  select case(fc_index)
	  
	  case(1)
	  fc = fcomm1
	  case(2)
	  fc = fcomm2
	  case default
	  print*,' DISASTER! BAD FC_INDX ',fc_index
	  stop
      end select
	  
#endif   
	if(wfn_mpiio) then
      if(isfragroot) then
			! move to start of vec/fragment
			call wfn_seekvecfragpos(filenumber, i, fg)
			! start reading
         if(fg == 1) then
            ! fragment 1 owns the vector overhead
#ifdef _MPI
            call BMPI_FILE_READ(my_mpi_file, vmagic, 1, status, ierr)
            call BMPI_FILE_READ(my_mpi_file, vnum, 1, status, ierr)
            if(i /= vnum) then
				print*,i,vnum, '(1)'
				call BMPI_ERR(msg_wrongvec)
			end if 
            call BMPI_FILE_READ(my_mpi_file, ejtbuf, 3, status, ierr)
#endif
         end if
         ! read the vector
#ifdef _MPI
         call BMPI_FILE_READ(my_mpi_file, v, size(v,1,basis_prec), status, ierr)
#endif
		end if
      ! All procs have to use the fg=1, isfragroot process
      ! as rootproc for BCAST:   br_rank_of_frag1_root
      rootproc = br_rank_of_frag1_root
      ! note size with dflt kind here is ok, we are only sending fragment
#ifdef _MPI
      call BMPI_BCAST(v, size(v), rootproc, fc, ierr)
#endif
   else
		if(iproc==0)then
		    ! no seeking needed, we read the whole thing
         call wfn_seekvecfragpos(filenumber, i, fg)
		    ! read i, e, xj, xt2 from file
			read(filenumber) vmagic, vnum, (ejtbuf(j), j = 1,3)
!			print*,' Test of reading magic ',vmagic,vnum,ejtbuf(:)
         if(i /= vnum) then
			print*,i,vnum, '(2)'
			 
#ifdef _MPI
            if(nproc > 1) call BMPI_ERR(msg_wrongvec)
#endif
            print *, msg_wrongvec
            stop 1 ! has msg
         end if
			read(filenumber)(v(jl),jl=1,dimbasis)
		end if
      ! BMPI_BCAST does chunking
      rootproc = 0
#ifdef _MPI
      call BMPI_BCAST(v, size(v), rootproc, fc, ierr)
#endif
	end if
	! update all nodes with vec properties . This works  with or without fragments.
   if(iproc == rootproc) then
      if(vmagic /= wfn_vecmagic) then
         print *, "vector", i, " is missing magic number"
         stop 1 ! has msg
      end if
   end if
#ifdef _MPI
	call BMPI_BCAST(ejtbuf, size(ejtbuf), rootproc, MPI_COMM_WORLD, ierr)
#endif
	e = ejtbuf(1)
	xj = ejtbuf(2)
	xt2 = ejtbuf(3)
		
	return
end subroutine wfn_readeigenvec

!=====================================================

!   subroutine readpivot
!
!   modified in 7.6.6 by CWJ to work in MPI
!   
!   the file of eigenvectors must be previously opened
!   and the header reader;
!   reads the # of kept vectors (nkeep), 
!   prints out E J, T^2 of each vector
!   then asks which to use as a pivot
!
!   OUTPUT
!      v(n), n = 1,dimbasis = vector read from disk
!      dnorm = norm of pivot vector
! 
!  CALLED BY: initialize_lanczos_vector in blanczoslib.f90
!
!   CALLS
!     wnf_read_nkeep
!     wfn_readeigenvec
!     wfn_close_file
!
subroutine readpivot
   use localvectors
   use precisions
   use basis
   use io
   use nodeinfo
   use bvectorlib_mod
   implicit none

   real (kind = 8) :: dnorm
   integer nkeep
   integer ikeep
   integer(4):: i,j,n, ierr
   real :: e,xj,xt2
   logical zeroflag
   integer :: aerr

!   if(iproc==0)print*,' Need to setup local vectors ??? '
!   call setup_localvectors
   
!   if(useNewReorthog) then
!      if(iproc == 0) print *, "readpiviot not supported with new orthog"
!      stop 1 ! has msg
!   end if

   call wfn_read_nkeep(oldwfnfile, nkeep)  ! does BCAST
   if(iproc==0) print*,' There are ', nkeep,' wavefunctions '

   if(iproc==0)then
        print*,nkeep,' states '
        print*,' '
        print*,' STATE      E         J          <H > '
	endif
    do ikeep = 1,nkeep
      i = ikeep
      ! new interface, we say which vector to read - it checks
      call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,e,xj,xt2) ! KSM: updated
      if(iproc==0)then
          write(6,101)i,e,xj,xt2

 101      format(i5,2x,4(1x,f11.4))
      end if
   end do ! ikeep

   if (iproc==0) then
      print*,' Which do you want as pivot? '
      read(5,*) ikeep
      if(ikeep < 1 .or. ikeep > nkeep) then
         print *, "Selected vector is not in range 1:", nkeep
         stop 1 ! has msg
      end if
 !     nkeep = ikeep;
   end if
!   call wfn_rewind(oldwfnfile)
!   call read_wfn_header(oldwfnfile,.false.) 
!   call BMPI_Barrier(icomm,ierr)
!   call wfn_read_nkeep(oldwfnfile, ikeep)  !  dummy read
   ! broadcast the selected value by the user
#ifdef _MPI
   call BMPI_BCAST(ikeep,1,0,MPI_COMM_WORLD,ierr)
#endif
!   do ikeep = 1,nkeep
      call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1, ikeep, e, xj, xt2) ! KSM: updated
!   enddo

   call wfn_close_file(oldwfnfile)

   return
end subroutine readpivot
! ============================================================
!   HEADER WRITE/READ ROUTINES
! ============================================================
!
!  writes header to .wfn file with information about the basis
!
subroutine write_wfn_header(filenumber)
  use system_parameters
  use sporbit
  use spstate
  use w_info
  use basis
  implicit none
  integer :: ierr
  integer :: filenumber
  integer :: it,n, dummy, v
#ifdef _MPI
        integer (kind=MPI_OFFSET_KIND) :: offset, dimbasispos
        type(MPI_file) :: my_mpi_file
		my_mpi_file = mpifile_out
#else
integer  :: offset, dimbasispos
		
#endif
			
  if(.not.writeout) return
  dimbasispos = 0
#ifdef _MPI
  call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
! KSM:  should really use root process of frag1 with mpi, which may not be iproc == 
  if(iproc==0) then
     !
     ! Write file type magic number and
     ! version number
     ! print *, "writing from iproc=0"	 

     call wfn_write_int4(filenumber, wfn_magic)   ! ids file
     call wfn_write_int4(filenumber, wfn_version) ! version number
     ! save position of offset to after header
     call wfn_write_int4(filenumber, 0) ! write over later
     !
     ! Write some zeroed space at the end for offsets
     ! The first offset will be to after the header
     ! write offset past header
     dummy = 0
     do n = 1, 5
       call wfn_write_int4(filenumber, dummy)
     end do
 
     !------------- WRITE INFO ON VALENCE PARTICLES
     !              PAY ATTENTION TO P-H CONJUGATION
     !        
     if(phconj(1))np(1) = -np(1)
     if(phconj(2))np(2) = -np(2)
     ! write(filenumber)np(1),np(2)
     call wfn_write_int4(filenumber, np(1))
     call wfn_write_int4(filenumber, np(2))
     if(phconj(1))np(1) = -np(1)
     if(phconj(2))np(2) = -np(2)
 
     !------------ WRITE OUT INFORMATION ON S.P. ORBITS 
     !  write(filenumber)isoflag
     !  write(filenumber)numorb(1),numorb(2)
     call wfn_write_logical(filenumber, isoflag)
     call wfn_write_int4(filenumber, numorb(1))
     call wfn_write_int4(filenumber, numorb(2))
     do it = 1,2
       do n = 1,numorb(it)
       	 ! all fields are integers
         !  write(filenumber)orbqn(it,n)%nr,orbqn(it,n)%j,orbqn(it,n)%l,  & 
         !            orbqn(it,n)%par, orbqn(it,n)%w
         call wfn_write_int4(filenumber, orbqn(it,n)%nr)
         call wfn_write_int4(filenumber, orbqn(it,n)%j)
         call wfn_write_int4(filenumber, orbqn(it,n)%l)
         call wfn_write_int4(filenumber, orbqn(it,n)%par)
         call wfn_write_int4(filenumber, orbqn(it,n)%w)
       enddo
     enddo
 
     !-------------- WRITE INFO ON S.P. STATES   -- CWJ 11/09  
     !               fixes a subtle bug when reading in without w-cuts
     !  write(filenumber)nsps(1),nsps(2)
     call wfn_write_int4(filenumber, nsps(1))
     call wfn_write_int4(filenumber, nsps(2))
      do it = 1,2
        do n = 1,nsps(it)
          !  write(filenumber) spsqn(it,n)%nr, spsqn(it,n)%j, spsqn(it,n)%m, & 
          !             spsqn(it,n)%l, spsqn(it,n)%w, spsqn(it,n)%par, & 
          !             spsqn(it,n)%orb, spsqn(it,n)%group
          call wfn_write_int4(filenumber, spsqn(it,n)%nr)
          call wfn_write_int4(filenumber, spsqn(it,n)%j)
          call wfn_write_int4(filenumber, spsqn(it,n)%m)
          call wfn_write_int4(filenumber, spsqn(it,n)%l)
          call wfn_write_int4(filenumber, spsqn(it,n)%w)
          call wfn_write_int4(filenumber, spsqn(it,n)%par)
          call wfn_write_int4(filenumber, spsqn(it,n)%orb)
          call wfn_write_int4(filenumber, spsqn(it,n)%group)
        enddo
     enddo
 
     !------------- WRITE OUT INFORMATION ON JZ, PARITY, W
     !  write(filenumber)jz,cparity,maxWtot
     !  write(filenumber)allsameparity, allsameW,spinless
     call wfn_write_int4(filenumber, jz)
     v = ICHAR(cparity)
     call wfn_write_int4(filenumber, v)
     call wfn_write_int4(filenumber, maxWtot)
     call wfn_write_logical(filenumber, allsameparity)
     call wfn_write_logical(filenumber, allsameW)
     call wfn_write_logical(filenumber, spinless)
 
     !------------ WRITE OUT BASIS DIMENSION
     !  write(filenumber)dimbasis
     ! print *, "dimbasis=", dimbasis
     if(nfragments == 1) inquire(filenumber, POS=dimbasispos)
     call wfn_write_bkint(filenumber, dimbasis)
  end if
#ifdef _MPI
  call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
  ! All processes reach here.  Broadcast offset if we are using
  ! MPI
  if(wfn_mpiio) then
    ! get offset
    offset = 13

    if(iproc==0) then
#ifdef _MPI	
!      call BMPI_FILE_get_position(filenumber, offset, ierr)
       call BMPI_FILE_get_position(my_mpi_file, offset, ierr)
#endif	   
       offset = offset + 4  ! allocate space for nkeep between header and vectors
       dummy = int(offset,4)  ! use 4 bytes for offset.
       ! KSM - FIXME Something wrong here
       call wfn_write_int4_at(filenumber, wfn_after_header_word*4, dummy)
#ifdef _MPI	   
!      call BMPI_FILE_seek(filenumber, offset, MPI_SEEK_SET, ierr)
       call BMPI_FILE_seek(my_mpi_file, offset, MPI_SEEK_SET, ierr)
#endif
       call wfn_write_int4(filenumber, dummy)
    end if
	
    ! sync all processes
    ! if(iproc == 0) print *, "dummy=", dummy
#ifdef _MPI
    call BMPI_BCAST(dummy, 1, 0, MPI_COMM_WORLD, ierr)
#endif	
    offset = dummy
    ! save position after header
    wfn_after_header_pos = offset;
    ! print *, "iproc=", iproc, ", position after header=", offset
    ! Set position for writing of nkeep
#ifdef _MPI	
!   call BMPI_FILE_SEEK(filenumber, offset-4, MPI_SEEK_SET, ierr)
    call BMPI_FILE_SEEK(my_mpi_file, offset-4, MPI_SEEK_SET, ierr)
#endif	
  else if(iproc == 0) then
    ! no fragments, iproc==0 does all the work
    inquire(UNIT=filenumber, POS=offset)
    offset = offset + 4  ! allocate space for nkeep between header and vectors
    dummy = int(offset-1,4)  ! positions start at 1, MPI_positions start at 0
    v = wfn_after_header_word*4+1
    write(filenumber, pos=v) dummy
    ! print *, "Wrote vec offset at ", v, ", value=", dummy
    ! move back to end
	
    write(filenumber, pos=dimbasispos) dimbasis
	
    ! write(filenumber) 1234
    wfn_after_header_pos = offset;
    ! print *, "no-frag, position after header=", offset
  end if
  return
end subroutine write_wfn_header
!===========================================================================
!
!  writes header to .wfn file with information about the basis
!  human-readable ASCII version added 7.9.1 
! 
subroutine write_wfn_header_ASCII(filenumber)
  use system_parameters
  use sporbit
  use spstate
  use w_info
  use basis
  implicit none
  integer :: ierr
  integer :: filenumber
  integer :: it,n, dummy, v
#ifdef _MPI
  integer (kind=MPI_OFFSET_KIND) :: offset, dimbasispos
#else
   integer  :: offset, dimbasispos
#endif
  
  if(.not.writeout) return
  if(wfn_mpiio)then
	  print*,' Should not be using routine write_wfn_header_ASCII in MPI mode '
	  stop
  end if
  dimbasispos = 0
  ! KSM:  should really use root process of frag1 with mpi, which may not be iproc == 0
  if(iproc==0) then
     !
     ! Write file type magic number and
     ! version number
     ! print *, "writing from iproc=0"	 
     
	 write(filenumber,*)wfn_magic
	 write(filenumber,*)wfn_version
!	 write(filenumber,*)0
	 
!     call wfn_write_int4(filenumber, wfn_magic)   ! ids file
!     call wfn_write_int4(filenumber, wfn_version) ! version number
     ! save position of offset to after header
!     call wfn_write_int4(filenumber, 0) ! write over later
     !
     ! Write some zeroed space at the end for offsets
     ! The first offset will be to after the header
     ! write offset past header
!     dummy = 0
!     do n = 1, 5
!       call wfn_write_int4(filenumber, dummy)
!     end do
 
     !------------- WRITE INFO ON VALENCE PARTICLES
     !              PAY ATTENTION TO P-H CONJUGATION
     !        
     if(phconj(1))np(1) = -np(1)
     if(phconj(2))np(2) = -np(2)
     ! write(filenumber)np(1),np(2)
	 write(filenumber,*)np(1),np(2)
!     call wfn_write_int4(filenumber, np(1))
!     call wfn_write_int4(filenumber, np(2))
     if(phconj(1))np(1) = -np(1)
     if(phconj(2))np(2) = -np(2)
 
     !------------ WRITE OUT INFORMATION ON S.P. ORBITS 
     !  write(filenumber)isoflag
     !  write(filenumber)numorb(1),numorb(2)
	 write(filenumber,*)isoflag
	 write(filenumber,*)numorb(1),numorb(2)
!     call wfn_write_logical(filenumber, isoflag)
!     call wfn_write_int4(filenumber, numorb(1))
!     call wfn_write_int4(filenumber, numorb(2))
     do it = 1,2
       do n = 1,numorb(it)
       	 ! all fields are integers
           write(filenumber,'(i5,4i5)')orbqn(it,n)%nr,orbqn(it,n)%j,orbqn(it,n)%l,  & 
                    orbqn(it,n)%par, orbqn(it,n)%w
!         call wfn_write_int4(filenumber, orbqn(it,n)%nr)
!         call wfn_write_int4(filenumber, orbqn(it,n)%j)
!         call wfn_write_int4(filenumber, orbqn(it,n)%l)
!         call wfn_write_int4(filenumber, orbqn(it,n)%par)
!         call wfn_write_int4(filenumber, orbqn(it,n)%w)
       enddo
     enddo
 
     !-------------- WRITE INFO ON S.P. STATES   -- CWJ 11/09  
     !               fixes a subtle bug when reading in without w-cuts
     !  write(filenumber)nsps(1),nsps(2)
	 write(filenumber,*)nsps(1),nsps(2)
!     call wfn_write_int4(filenumber, nsps(1))
!     call wfn_write_int4(filenumber, nsps(2))
      do it = 1,2
        do n = 1,nsps(it)
            write(filenumber,'(i5,8i5)') spsqn(it,n)%nr, spsqn(it,n)%j, spsqn(it,n)%m, & 
                       spsqn(it,n)%l, spsqn(it,n)%w, spsqn(it,n)%par, & 
                       spsqn(it,n)%orb, spsqn(it,n)%group
!          call wfn_write_int4(filenumber, spsqn(it,n)%nr)
!          call wfn_write_int4(filenumber, spsqn(it,n)%j)
!          call wfn_write_int4(filenumber, spsqn(it,n)%m)
!          call wfn_write_int4(filenumber, spsqn(it,n)%l)
!          call wfn_write_int4(filenumber, spsqn(it,n)%w)
!          call wfn_write_int4(filenumber, spsqn(it,n)%par)
!          call wfn_write_int4(filenumber, spsqn(it,n)%orb)
!          call wfn_write_int4(filenumber, spsqn(it,n)%group)
        enddo
     enddo
 
     !------------- WRITE OUT INFORMATION ON JZ, PARITY, W
       write(filenumber,*)jz,cparity,maxWtot
       write(filenumber,*)allsameparity, allsameW,spinless
!     call wfn_write_int4(filenumber, jz)
!     v = ICHAR(cparity)
!     call wfn_write_int4(filenumber, v)
!     call wfn_write_int4(filenumber, maxWtot)
!     call wfn_write_logical(filenumber, allsameparity)
!     call wfn_write_logical(filenumber, allsameW)
!     call wfn_write_logical(filenumber, spinless)
 
     !------------ WRITE OUT BASIS DIMENSION
       write(filenumber,*)dimbasis
     ! print *, "dimbasis=", dimbasis
!     if(nfragments == 1) inquire(filenumber, POS=dimbasispos)
!     call wfn_write_bkint(filenumber, dimbasis)
  end if
  ! All processes reach here.  Broadcast offset if we are using
  ! MPI

  return
end subroutine write_wfn_header_ASCII
!===========================================================================
!
!
!  reads header from .wfn file with information about the basis
!
! CALLED BY :
!    menu
!
subroutine read_wfn_header(filenumber,verbose)
  use system_parameters
  use sporbit
  use spstate
  use w_info
  use basis
  use butil_mod
  use reporter
  implicit none
  integer:: filenumber
  logical:: verbose,dummyread
  integer:: it,n, v
  integer:: ierr
  integer :: fpos
  integer :: aerr


#ifdef _MPI
    call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif   
  if(iproc == 0) then
     ! print *, "reading wfn header, filenumber=", filenumber
     ! read initial 8 word table
     call wfn_read_int4(filenumber, v)
	 !print*,' testing magic number ',v,wfn_magic
     ! if(iproc == 0) print *, "read magic=", v
     if(v /= wfn_magic) then
        if(iproc == 0) print *,v,wfn_magic, "WFN file is missing magic number"
        stop 1 ! has msg
     end if
	 
     call wfn_read_int4(filenumber, v) ! version number
     ! print *, "read version=", v
     if(v /= wfn_version) then
        if(iproc == 0) print *, "WFN file is wrong version"
        stop 1 ! has msg
     end if
     call wfn_read_int4(filenumber, v) ! word not used yet
     call wfn_read_int4(filenumber, v) ! word 3 (starting with word 0)
     ! print *, "read vec pos=", v
     if(.not. wfn_mpiio) v = v + 1 ! NON MPI read/write starts at 1
     wfn_after_header_pos = v ! diff size
     do n=1, 4
       call wfn_read_int4(filenumber, v)
     end do
	 
  !------------- READ INFO ON VALENCE PARTICLES
     ! read(filenumber)np(1),np(2)
     call wfn_read_int4(filenumber, np(1))
     call wfn_read_int4(filenumber, np(2))
!	  do it = 1,2
!		 if(np(it)< 0)then
!			 print*,' NEED TO FIX READING IN p-h conjugation '

!		 end if
!	  end do

     !------------ READ OUT INFORMATION ON S.P. ORBITS 
     ! read(filenumber)isoflag
     call wfn_read_logical(filenumber, isoflag)
     ! read(filenumber)numorb(1),numorb(2)
     call wfn_read_int4(filenumber, numorb(1))
     call wfn_read_int4(filenumber, numorb(2))
  endif  
#ifdef _MPI  
  call BMPI_Barrier(MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(wfn_after_header_pos, 1, 0, MPI_COMM_WORLD, ierr)
  
  call BMPI_BCAST(np(1), 1, 0, MPI_COMM_WORLD, ierr)
  call BMPI_BCAST(np(2), 1, 0, MPI_COMM_WORLD, ierr)
  call BMPI_BCAST(numorb(1), 1, 0, MPI_COMM_WORLD, ierr)
  call BMPI_BCAST(numorb(2),1,0, MPI_COMM_WORLD, ierr)
  call BMPI_BCAST(isoflag, 1, 0, MPI_COMM_WORLD, ierr)
#endif
  numorbmax = bmax(numorb(1),numorb(2))
!------------------ALLOCATE MEMORY------------------------------------------
  if(.not.allocated(orbqn)) then
     allocate ( orbqn(2,numorbmax), stat=aerr )
     if(aerr /= 0) call memerror("read_wfn_header 1")
  end if
  do it = 1,2
      do n = 1,numorb(it)
        if(iproc==0) then
           ! read(filenumber)orbqn(it,n)%nr,orbqn(it,n)%j,orbqn(it,n)%l,  & 
           !           orbqn(it,n)%par, orbqn(it,n)%w 
    	     call wfn_read_int4(filenumber, orbqn(it, n)%nr)
    	     call wfn_read_int4(filenumber, orbqn(it, n)%j)
    	     call wfn_read_int4(filenumber, orbqn(it, n)%l)
    	     call wfn_read_int4(filenumber, orbqn(it, n)%par)
    	     call wfn_read_int4(filenumber, orbqn(it, n)%w)
        end if  
    	! update all other nodes
#ifdef _MPI
        call BMPI_BCAST(orbqn(it,n)%nr,1,0,MPI_COMM_WORLD,ierr)
        call BMPI_BCAST(orbqn(it,n)%j ,1,0,MPI_COMM_WORLD,ierr)
        call BMPI_BCAST(orbqn(it,n)%l ,1,0,MPI_COMM_WORLD,ierr)
        call BMPI_BCAST(orbqn(it,n)%par,1,0,MPI_COMM_WORLD,ierr)
        call BMPI_BCAST(orbqn(it,n)%w ,1,0,MPI_COMM_WORLD,ierr)
#endif
      enddo

  enddo
  
  if(iproc==0)then
     ! read(filenumber)nsps(1),nsps(2)
     call wfn_read_int4(filenumber, nsps(1))
     call wfn_read_int4(filenumber, nsps(2))
  end if
#ifdef _MPI
  call BMPI_BCAST(nsps,2,0,MPI_COMM_WORLD,ierr)
#endif
  if(.not.allocated(spsqn)) then
     allocate(spsqn(2,bmax(nsps(1),nsps(2))), stat=aerr)
     if(aerr /= 0) call memerror("read_wfn_header 2")
  end if
  ! print *, "Abount to read spsqn"
  do it = 1,2
     do n = 1,nsps(it)
       if(iproc==0) then
          ! read(filenumber) spsqn(it,n)%nr, spsqn(it,n)%j, spsqn(it,n)%m, & 
          !          spsqn(it,n)%l, spsqn(it,n)%w, spsqn(it,n)%par, & 
          !          spsqn(it,n)%orb, spsqn(it,n)%group
    	  call wfn_read_int4(filenumber, spsqn(it, n)%nr)
    	  call wfn_read_int4(filenumber, spsqn(it, n)%j)
    	  call wfn_read_int4(filenumber, spsqn(it, n)%m)
    	  call wfn_read_int4(filenumber, spsqn(it, n)%l)
    	  call wfn_read_int4(filenumber, spsqn(it, n)%w)
    	  call wfn_read_int4(filenumber, spsqn(it, n)%par)
    	  call wfn_read_int4(filenumber, spsqn(it, n)%orb)
    	  call wfn_read_int4(filenumber, spsqn(it, n)%group)
       end if
       ! update other nodes
#ifdef _MPI
       call BMPI_BCAST(spsqn(it,n)%nr,1,0,MPI_COMM_WORLD,ierr)
       call BMPI_BCAST(spsqn(it,n)%j ,1,0,MPI_COMM_WORLD,ierr)
       call BMPI_BCAST(spsqn(it,n)%m ,1,0,MPI_COMM_WORLD,ierr)
       call BMPI_BCAST(spsqn(it,n)%l ,1,0,MPI_COMM_WORLD,ierr)
       call BMPI_BCAST(spsqn(it,n)%w,1,0,MPI_COMM_WORLD,ierr)
       call BMPI_BCAST(spsqn(it,n)%par,1,0,MPI_COMM_WORLD,ierr)
       call BMPI_BCAST(spsqn(it,n)%orb,1,0,MPI_COMM_WORLD,ierr)
       call BMPI_BCAST(spsqn(it,n)%group,1,0,MPI_COMM_WORLD,ierr)
#endif
    enddo

  enddo
!------- ALLOCATE MEMORY  -- CWJ 11/09
!------------- READ IN INFORMATION ON JZ, PARITY, W
  if(iproc==0) then
     ! read(filenumber)jz,cparity,maxWtot
     call wfn_read_int4(filenumber, jz)
     call wfn_read_int4(filenumber, n)
     cparity = ACHAR(n) ! convert integer back to character
     call wfn_read_int4(filenumber, maxWtot)
  end if
#ifdef _MPI
  call BMPI_BCAST(jz,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(maxWtot,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(cparity,1,0,MPI_COMM_WORLD,ierr)
#endif
  select case (cparity)
    case ('+','0')
       iparity = 1
    case ('-')
       iparity = 2
    case default
       write(6,*)' Problem with parity ',cparity
       stop
  end select

  if(iproc==0) then 
     ! read(filenumber)allsameparity, allsameW,spinless
     call wfn_read_logical(filenumber, allsameparity)
     call wfn_read_logical(filenumber, allsameW)
     call wfn_read_logical(filenumber, spinless)
  end if
#ifdef _MPI
  call BMPI_BCAST(allsameparity,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(allsameW,1,0,MPI_COMM_WORLD,ierr)
  call BMPI_BCAST(spinless,1,0,MPI_COMM_WORLD,ierr)
#endif
!------------ WRITE OUT BASIS DIMENSION

  if(iproc==0)then
     ! read(filenumber)dimbasischeck
     call wfn_read_bkint(filenumber, dimbasischeck)
     print *, "dimbasischeck=", dimbasischeck
  end if
#ifdef _MPI
  call BMPI_BCAST(dimbasischeck,1,0,MPI_COMM_WORLD,ierr)
#endif
!-------------------------------
  if(verbose .and. iproc==0)then
     print*,' Valence Z, N = ',np(1),np(2)
     print*,' Single particle space : '
     if(isoflag)then
         print*,'   N    L  2xJ '
         do n = 1,numorb(1)
            write(6,'(3i5)')orbqn(1,n)%nr,orbqn(1,n)%l,orbqn(1,n)%j
         end do
         print*,' '
         print*,' Total # of orbits = ',numorb(1)
     else
         print*,' Oops not yet set up '
     endif
     if(allsameparity)then
        print*,' 2 x Jz = ',jz
     else
        if(iparity == 1)then
           print*,' 2 x Jz, parity = ',Jz,'+'
        else
           print*,' 2 x Jz, parity = ',Jz,'-'
        endif
     endif
     it = 1
     if(writeout)then
        if(isoflag)then
            write(resultfile,*) ' N       L       2xJ '
            do n = 1,numorb(1)
               write(resultfile,*) orbqn(it,n)%nr,orbqn(it,n)%l,orbqn(it,n)%j
            end do
        else
            print*,' Oops not yet set up '
        endif
     endif
  endif
  return
end subroutine read_wfn_header

!=====================================================================
!   subroutine wfn_seekvecfragpos(filenumber, i, fg)
!
!  seeks to start of fragment fg in vector i in the wfn file
!
!  note:  The first fragment carries the vector overhead.  This way
!         the file is independent of fragmentation choices
!
!  INPUT
!       filenumber = file id to operate on
!       i = vector number, starting with vector 1
!       fragment = piece of vector - see module fragments
! 
!  CALLED BY:
!   wfn_readeigenvec
!   wfn_writeeigenvec
!
!=====================================================================
! MPI UPDATE 7.10.9: this needs to be fixed, so that
! filenumber will trigger the correct MPI file
subroutine wfn_seekvecfragpos(filenumber, i, fg)
   use fragments
	use localvectors
	use io
	implicit none
	integer, intent(in) :: filenumber, i, fg
	integer :: fi
	integer:: ierr
   integer :: overhead
   integer(kind=4) :: dummy4
!   integer :: vs, pos, fragpos, veclen, nw, offset
!   integer(kind=MPI_OFFSET_KIND) :: vs, pos, fragpos, veclen, nw, offset

#ifdef _MPI
        type(MPI_file) :: my_mpi_file
   integer(kind=MPI_OFFSET_KIND) :: vs, pos, fragpos, veclen, nw, offset

		select case(filenumber)
		
			case(13)	  ! wfnfile 
	   		  my_mpi_file = mpifile_out
		  
			case (14)  !oldwfnfile
			my_mpi_file = mpifile_in
		  
			case default
		
			print*,' bad file number ',filenumber
			stop
		
	   	end select
#else
   integer(8) :: vs, pos, fragpos, veclen, nw, offset

#endif
   

   overhead = 5 * 4 ! (i, e, xj, xt2) + vec data 
   if(lanc_prec == 4) then
      vs = 4
   else
      vs = 8
   end if
	! should do this part once, could call from end of header read
   pos = 0
   fragpos = 0
   do fi = 1, nfragments
      if(fi == fg) fragpos = pos;
      nw = basestop(fi)- basestart(fi) + 1;
      pos = pos + nw*vs ! (i, e, xj, xt2) + vec data 
   end do
   veclen = pos + overhead  ! total vector size, all fragments
   ! compute offset to what we want to write
   offset = wfn_after_header_pos + (i-1)*veclen + fragpos
   if(fg /= 1) offset = offset + overhead
   if(wfn_mpiio) then
#ifdef _MPI	   
      call BMPI_FILE_seek(my_mpi_file, offset, MPI_SEEK_SET, ierr)
#endif
   else
      ! seek to start of vector

	  ! rewind necessary to clear state of filenumber before read when using ifort
	  ! otherwise, read will hang the second time around. - KSM
      rewind(filenumber)  
      read(filenumber, pos=offset-4) dummy4
   end if
	return
end subroutine wfn_seekvecfragpos

!========================================================
! Wrapper for rewinding to the front of the file that
! handles the difference between MPI-IO and normal IO
subroutine wfn_rewind(fn)
	use io
	use localvectors
   implicit none
   integer, intent(in) :: fn
   integer:: ierr
#ifdef _MPI
           type(MPI_file) :: my_mpi_file
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
		if(.not.isfragroot)return ! ADDED in 7.11.2 to fix errors in calling
#endif
   if(wfn_mpiio) then
      call wfn_checkmpiio()
      ! pos 0 is beginning for MPI_FILE_seek
#ifdef _MPI	  
      call BMPI_FILE_seek(my_mpi_file, int(0, MPI_OFFSET_KIND), MPI_SEEK_SET, ierr)
#endif
   else
      ! pos 1 is beginning for inquire/fseek/... on normal files
      ! all file I/O happens in iproc==0 for normal files
      if(iproc == 0) rewind(fn)
   end if
end subroutine wfn_rewind

!===========================================================
!  MISC READ/WRITE ROUTINES
!===========================================================
! sometimes the file is read manually and nkeep is
! discarded without overwriting the global
subroutine wfn_write_nkeep(nkeep)
   implicit none
   integer, intent(in):: nkeep
   integer :: ierr
   integer(kind=4) :: nkeep4
   if(writeout .and. write_wfn .and. iproc == 0) then
      nkeep4 = nkeep
      call wfn_write_int4(wfnfile, nkeep4)
   end if
end subroutine wfn_write_nkeep

subroutine wfn_read_nkeep(fn, nkeep)
   implicit none
   integer, intent(in) :: fn
   integer, intent(out) :: nkeep
   integer :: ierr
   integer(kind=4) :: nkeep4  
   if(iproc == 0) then
      call wfn_read_int4(fn, nkeep4)
      nkeep = nkeep4
   end if
#ifdef _MPI	    
   if(nproc > 1) call BMPI_BCAST(nkeep,1, 0, MPI_Comm_World, ierr)
#endif
end subroutine wfn_read_nkeep

!===========================================================================

!
! Write one integer to the wfn file. 
! 
subroutine wfn_write_int4(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer, intent(in) :: v
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
   integer(kind=4) :: data4

#ifdef _MPI
    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif

   data4 = v
   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI	  
      call BMPI_FILE_WRITE(my_mpi_file, data4, 1, status, ierr)
#endif
   else
      write(fn) data4
   end if 
end subroutine wfn_write_int4
! write one integer8 to the file
subroutine wfn_write_int8(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) :: v
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
!
#ifdef _MPI
    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif

   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI	  
      call BMPI_FILE_WRITE(my_mpi_file, v, 1, status, ierr)
#endif
   else
      write(fn) v
   end if 
end subroutine wfn_write_int8
! read one integer8 from the file
subroutine wfn_read_int8(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(out) :: v
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
!
#ifdef _MPI
    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif

   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI	  
      call BMPI_FILE_READ(my_mpi_file, v, 1, status, ierr)
#endif
   else
      read(fn, err=911, iostat=ierr) v
      return
911   print *, "read error: ", ierr
      call printstacktrace
      stop 1 ! has msg
   end if
end subroutine wfn_read_int8
!========================================================

subroutine wfn_write_real4(fn, v)
   implicit none
   integer, intent(in) :: fn
   real, intent(in) :: v
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
   real(kind=4) :: data4
#ifdef _MPI
    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif
   data4 = v
   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI	  
      call BMPI_FILE_WRITE(my_mpi_file, data4, 1, status, ierr)
#endif
   else
      write(fn) data4
   end if 
end subroutine wfn_write_real4


! read one real from the file
subroutine wfn_read_real4(fn, v)
   implicit none
   integer, intent(in) :: fn
   real, intent(out) :: v
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
   real(kind=4) :: data4
#ifdef _MPI
    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif
   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI	  
      call BMPI_FILE_READ(my_mpi_file, data4, 1, status, ierr)
#endif
      v = data4
   else
      read(fn, err=911, iostat=ierr) data4
      v = data4
      return
911   print *, "read error: ", ierr
      call printstacktrace
      stop 1 ! has msg
   end if
end subroutine wfn_read_real4
!========================================================


! This writes a "basis kind" integer
! We just convert to 8 byte
subroutine wfn_write_bkint(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer (kind=basis_prec) :: v
   integer (kind=8) :: v8

   v8 = v
   call wfn_write_int8(fn, v8)
end subroutine wfn_write_bkint
!========================================================
! read one integer from the file
! just read as 8 byte integer and assign
subroutine wfn_read_bkint(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer (kind = basis_prec), intent(out) :: v
   integer (kind=8) :: v8

!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr

   call wfn_read_int8(fn, v8)
   v = v8
end subroutine wfn_read_bkint
!========================================================
! write integer at specified position.   Used to update header of
! file
subroutine wfn_write_int4_at(fn, fpos, v)
   implicit none
   integer, intent(in) :: fn, fpos
   integer(kind=4), intent(in) :: v
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
#ifdef _MPI
    integer (kind=MPI_OFFSET_KIND) :: offset
    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif
   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI	  
      offset = fpos
      call BMPI_FILE_seek(my_mpi_file, offset, MPI_SEEK_SET, ierr)
      call BMPI_FILE_WRITE(my_mpi_file, v, 1, status, ierr)
#endif
   else
      write(fn, pos=fpos) v
   end if 
end subroutine wfn_write_int4_at
!========================================================
! read one integer from the file
subroutine wfn_read_int4(fn, v)
   implicit none
   integer, intent(in) :: fn
   integer, intent(out) :: v
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
   integer(kind=4) :: data4
#ifdef _MPI
    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
				
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif

   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI
      call BMPI_FILE_READ(my_mpi_file, data4, 1, status, ierr)
#endif
      v = data4
   else
      read(fn, err=911, iostat=ierr) data4
      v = data4
      return
911   print *, "read error: ", ierr
      call printstacktrace
      stop 1 ! has msg
   end if
end subroutine wfn_read_int4
!========================================================
! read one integer at a position
subroutine wfn_read_int4_at(fn, fpos, v)
   implicit none
   integer, intent(in) :: fn, fpos
   integer, intent(out) :: v
   integer(kind=4) :: data4
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
#ifdef _MPI
    integer (kind=MPI_OFFSET_KIND) :: offset

    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif

   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI	  
      offset = fpos;
      call BMPI_FILE_SEEK(my_mpi_file, offset, MPI_SEEK_SET, ierr)
      call BMPI_FILE_READ(my_mpi_file, data4, 1, status, ierr)
      v = data4
#endif
   else
      read(fn, pos=fpos, err=911, iostat=ierr) data4
      v = data4
      return
911   print *, "read error: ", ierr
      call printstacktrace
      stop 1 ! has msg
   end if
end subroutine wfn_read_int4_at
!========================================================
! write a boolean into a file.  Uses same space as integer
subroutine wfn_write_logical(fn, lg)
   implicit none
   logical, intent(in) :: lg
   integer, intent(in) :: fn
   integer(kind=4) :: v
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
#ifdef _MPI
    integer (kind=MPI_OFFSET_KIND) :: offset

    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif

   v = 0
   if(lg) v = 1
   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI	  
      call BMPI_FILE_WRITE(my_mpi_file, v, 1, status, ierr)
#endif
   else
      write(fn) v
   endif 
end subroutine wfn_write_logical
!========================================================
! read boolean from file
subroutine wfn_read_logical(fn, lg)
   implicit none
   integer, intent(in) :: fn
   logical, intent(out) :: lg
   integer(kind=4) :: v
!   integer :: status(MPI_STATUS_SIZE)
   integer :: ierr
#ifdef _MPI
    integer (kind=MPI_OFFSET_KIND) :: offset

    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif   
   if(wfn_mpiio) then
      call wfn_checkmpiio()
#ifdef _MPI	  
      call BMPI_FILE_READ(my_mpi_file, v, 1,  status, ierr)
#endif
      lg = v .ne. 0
   else
      read(fn, err=911, iostat=ierr) lg
      return
911   print *, "read error: ", ierr
      call printstacktrace
      stop 1 ! has msg
   end if
end subroutine wfn_read_logical

!----------------------------------------------------
! SUBROUTINES TO OPEN FILES
!========================================================
! fragment mode:  Used to open an MPI file for writing out the
! eigenvectors (wavefunctions) for the first N states
! Only iproc==0 has afilename correct at entry so we
! have to broadcast it.
! Then we open the file.
!
subroutine wfn_wopen_file_frag(fn, afilename)
   implicit none
   integer, intent(inout) :: fn
   integer:: ierr
   character(len=*), intent(in) :: afilename
   character(len=1024) :: wfnfilename

#ifdef _MPI
    integer(kind=MPI_OFFSET_KIND), parameter :: zsize = 0
    type(MPI_file) :: my_mpi_file
		
!		select case(fn)
		
!		case(13)	   != wfnfile
!   		  my_mpi_file = mpifile_out
		  
!		case (14) ! = oldwfnfile
!		my_mpi_file = mpifile_in
		  
!		case default
		
!		print*,' bad file number ',fn
!		stop
		
!   	    end select
#endif   
   ! distribute the file name
   wfnfilename = afilename   
!   print*,iproc,'(X1a)',wfnfilename,len(wfnfilename)
#ifdef _MPI
   call BMPI_BCAST(wfnfilename, LEN(wfnfilename) , 0, MPI_COMM_WORLD, ierr)
   select case (fn)
   case (13) 
   call BMPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(wfnfilename), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpifile_out, ierr)
   case (14)
   call BMPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(wfnfilename), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpifile_in, ierr)
   case default
   print*,' bad file number ',fn
   		stop
	end select
#endif   
   if(ierr .ne. 0) then
      if(iproc == 0) print *, "Can't wopen .wfn file: ", TRIM(wfnfilename)
      stop 1 ! has msg
   end if
   ! truncate the file
#ifdef _MPI   
   call BMPI_FILE_SET_SIZE(my_mpi_file, zsize, ierr)
#endif
   
end subroutine wfn_wopen_file_frag
!========================================================
! fragment mode:  Used to open the MPI file for the WFN file
subroutine wfn_ropen_file_frag(fn, afilename)
   implicit none
   integer, intent(inout) :: fn
   integer:: ierr
   character(len=*), intent(in) :: afilename
   character(len=1024) :: wfnfilename
#ifdef _MPI
    integer (kind=MPI_OFFSET_KIND) :: offset

    type(MPI_status) :: status
    type(MPI_file) :: my_mpi_file
		
#endif   

   ! distribute the file name
   wfnfilename = afilename
#ifdef _MPI   
   call BMPI_BCAST(wfnfilename, LEN(wfnfilename), 0, MPI_COMM_WORLD, ierr)
   select case(fn)
   
   case(13)
   call MPI_FILE_open(MPI_COMM_WORLD, wfnfilename, MPI_MODE_RDONLY, MPI_INFO_NULL, mpifile_out, ierr)
   case (14)
   call MPI_FILE_open(MPI_COMM_WORLD, wfnfilename, MPI_MODE_RDONLY, MPI_INFO_NULL, mpifile_in, ierr)
   case default
   print*,' bad file number ',fn
   		stop
	end select
#endif
   if(ierr .ne. 0) then
      if(iproc == 0) print *, "Can't ropen .wfn file: ", TRIM(wfnfilename)
      stop 1 ! has msg
   end if

end subroutine wfn_ropen_file_frag

!=========================================================================
!   Open the wfn file for reading only.
!
subroutine wfn_ropen_file(filenumber)
	use reporter
   implicit none
   integer, intent(inout) :: filenumber
   integer :: ilast
   logical :: found
   character (len=1024) :: wfnfilename
   integer :: ierr

   found = .false.
   if(iproc == 0) then
      do while(.not.found)
         print*,' Enter input name of .wfn file '
         if(auto_input)then
            read(autoinputfile,'(a)')outfile
         else
            read(5,'(a)')outfile
            write(autoinputfile,'(a)')outfile
         end if
         ilast = index(outfile,' ')-1
         wfnfilename = outfile(1:ilast)//'.wfn'
		 wfn_in_filename = outfile(1:ilast)//'.wfn'
!...... MODIFIED in 7.6.5 to INQUIRE statement		 
         inquire(file=wfnfilename, &
             exist=found)
	     if(.not.found)then
				 print*,outfile(1:ilast)//'.wfn does not exist '
		elseif(auto_input)then
		    print*,' reading from ',outfile(1:ilast)//'.wfn  '
					 
		 end if

      end do
   end if
   ! all procs get here
   ! KSM - open on all mpi nodes
   ! The open above is just a trial run.
   wfn_mpiio = .false.
   if(nproc == 1) then
       open(unit=filenumber, file=wfnfilename, &
   	      action='READ', access='stream', form='unformatted', &
           status = 'old')	   
   else
      wfn_mpiio = .true.
      ! We successfully opened it on the root.  Now we close,
      ! broadcast the filename, and open on all processes.
      ! with MPI IO
      call wfn_ropen_file_frag(filenumber, wfnfilename)
   end if
   return
end subroutine wfn_ropen_file

!=====================================================================
! Open wfn file for writing.
! If we have fragments, open with MPI call
! We still do trial open on iproc==0 to make
! sure file is writeable.
!
subroutine wfn_wopen_file(filenumber,basisfile)
   implicit none
   integer, intent(inout) :: filenumber
   logical, intent(in) :: basisfile   ! added in 7.7.2; this way we can reuse the same routines
   integer :: ilast
   character (len=1024) :: wfnfilename
   character (len=100) :: basewfnfilename
   logical :: found
   integer :: stat
   logical :: dolink

   ilast = index(outfile,' ')-1
   dolink = .false.
   if(basisfile)then
      basewfnfilename = outfile(1:ilast)//'.bas'
   else
      basewfnfilename = outfile(1:ilast)//'.wfn'
  end if


  
   if(scratch_dir == '.') then
      wfnfilename = basewfnfilename
   else
      wfnfilename = TRIM(scratch_dir) // '/' // TRIM(basewfnfilename)
      dolink = .true.
   end if
   found = .false.
   
   if ( iproc == 0 .and. (write_wfn .or. basisfile)) then
      open(unit=wfnfile, iostat=stat, file=wfnfilename, status='old')
	  
      if(stat == 0) close(wfnfile, status='delete')
	  
	  if(baseASCII)then                                   ! ADDED 7.9.1
		  if(.not.basisfile)then
			  print*,' Something wrong in opening basis file? '
			  stop
		  end if
	      open(unit=wfnfile, file=wfnfilename, &
	            action='WRITE',  status = 'unknown', err=101)		  
	  else	  
         open(unit=wfnfile, file=wfnfilename, &
           access='stream', form='unformatted', action='WRITE', & 
           status = 'unknown', err=101)
	   end if
	   
      found = .true.
      ! make link from working directory to file in scratch directory
      ! so the file is where we expect it
      if(dolink) then
         call system("ln -s " // TRIM(wfnfilename) // " .")
      end if
101   continue

      if(.not.found)then
    	 print *, "Can't open .wfn file for writing: ", wfnfilename
         print*," (in wfn_wopen_file) "
    	 stop 1 ! has msg
      end if
   end if
   
   ! all procs get here
   ! KSM - open on all mpi nodes
   wfn_mpiio = .false.
   if(nproc > 1) then
	   
      wfn_mpiio = .true.
	  
      ! We successfully opened it on the root.  Now we close
      ! and delete it, broadcast the filename, and open on all proceses.
      ! with MPI IO
      if(iproc == 0) close(unit = filenumber, status='delete')	  

      call wfn_wopen_file_frag(filenumber, wfnfilename)
	  
	  
   end if
   return
end subroutine wfn_wopen_file
!=====================================================================
!
! opens files for restarting
! started 7/2010 by CWJ @SDSU
!
! NOTE: not yet parallelized
!
  subroutine open_restart_files

  use program_info
  use io
  use lanczos_info
  use nodeinfo
  implicit none
  integer filenumber
  integer ilast
  logical found

  found = .false.
  if(iproc==0)then
     print*,' Note: this option will reuse existing files and overwrite them '
     print*,' IMPORTANT: If you want to save previous data, COPY AND RENAME '
     print*,' '

     do while(.not.found)
        print*,' Enter previous name of output file '
        read(5,'(a)')outfile
        ilast = index(outfile,' ')-1
!..................... wfn file
        found = .false.
        open(unit=   wfnfile,file=outfile(1:ilast)//'.wfn', &
    	   access='stream', form='unformatted', & 
           status = 'old',err=101)
        found = .true.
101  continue
        if(.not.found)then
          print*,' .wfn file does not exist '
           exit
        endif
!..................... lanczos vector file
!  NOTE for parallel operations this has to be modified
!
        found = .false.
        open(unit=   lvec_file,file=outfile(1:ilast)//'.lvec',form='unformatted', & 
           status = 'old',err=102)
        found = .true.
102  continue
        if(.not.found)then
            print*,' .lvec file does not exist '
            call wfn_close_file(wfnfile)
            exit
        endif
!..................... lanczos coeficient file
!
        found = .false.
        open(unit=   coef_file,file=outfile(1:ilast)//'.lcoef',form='formatted', & 
          status = 'old',err=103)
        found = .true.
103  continue
        if(.not.found)then
           print*,' .lcoef file does not exist '
            call wfn_close_file(wfnfile)
            close(lvec_file)
           exit
        endif

       end do
       open(unit=resultfile,file=outfile(1:ilast)//'.res',status = 'unknown')
       write(6,*)'  Writing results over old .res file '
       write(resultfile,*)' BIGSTICK Version ',version,lastmodified
       writeout = .true.

  end if
!........... OPEN RESULTS FILE...................

  writetodisk = .true.

  return
  end subroutine open_restart_files
!-------------------------------------
! SUBROUTINE TO CLOSE FILE(S)
!========================================================
subroutine wfn_close_file(fn)
   implicit none
   integer, intent(inout) :: fn
   integer :: ierr
#ifdef _MPI
    type(MPI_file) :: my_mpi_file
		
		select case(fn)
		
		case(13)	   != wfnfile
   		  my_mpi_file = mpifile_out
		  
		case (14) ! = oldwfnfile
		my_mpi_file = mpifile_in
		  
		case default
		
		print*,' bad file number ',fn
		stop
		
   	    end select
#endif   

   if(fn <= 0) then
      print *, "wfn_close_file: wfnfile is not open"
      stop 1 ! has msg
   end if
   if(wfn_mpiio) then
	   if(iproc==0)print*,' CLOSING WFN FILE ',fn
#ifdef _MPI	   
      call BMPI_FILE_close(my_mpi_file, ierr)
#endif
   else
      close(fn)
   end if
   fn = 0
end subroutine wfn_close_file

!----------------------------------------------------------
! MISC FILES
!=====================================================================
!
! check dimension of basis when reading in wavefunction
!
  subroutine checkbasisdim
  use basis
  use io
  implicit none
  if(dimbasis /= dimbasischeck)then
    print*,' hmm problem with basis dimension '
    print*,' should be ',dimbasischeck,' but I get ',dimbasis
    stop
  endif

  return
  end subroutine checkbasisdim
  
!=====================================================================

end module wfn_mod



