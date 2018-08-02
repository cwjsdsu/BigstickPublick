!===================================================================
!
! MPI wrapper module
!
! This module accomplishes several things
!  1) There is a single point where MPI is included
!  2) We implement type specific wrappers so that there are interfaces
!     for all MPI calls used in Bigstick that implement type checking.  We've
!     had a number of bugs related to integer sizes ...
!  3) There are reported bugs on some platforms for the size of reduce
!     and bcast transfers.  This interface will break them up if needed.
!  4) Type selection happens automatically vial overload, so it is no longer
!     possible to do bad things like have an integer(kind=4) buffer sent
!     with an MPI_INTEGER  type (which won't work on machines with 
!     8 byte default integers
!
! I decided to use overloading so that the polymorphic calls are
! automatically bound to the right wrapper.
!
! Initial version:  Ken McElvain   17Jan2015
!
module bmpi_mod
   implicit none
   ! import constants, which we still use directly
   include 'mpif.h'

   ! limit on size of BCAST/REDUCE/ALLREDUCE calls
   ! such calls are automatically chunked at the following byte size
   integer(kind=8), parameter :: bmpi_maxbcastsize = 50000000
   integer(kind=8), parameter :: bmpi_iosize = 50000000

   !============================================================================
   ! interfaces to perform vector/scalar overloading for read/write/bcast/reduce
   ! Make use of procedure overloading to select the right function
   ! based on the type of the arguments
   ! Also implements MPI_IN_PLACE when the input buffer is droppped
   ! Function names ending in IP implement MPI_IN_PLACE

   !============================================================================
   ! READ
   !
   interface BMPI_FILE_READ
      MODULE PROCEDURE BMPI_FILE_READ_INTEGER4_VEC_LC, BMPI_FILE_READ_INTEGER4_VEC_SC, &
         BMPI_FILE_READ_INTEGER4_SCALAR, &
         BMPI_FILE_READ_INTEGER8_VEC_LC, BMPI_FILE_READ_INTEGER8_VEC_SC, &
         BMPI_FILE_READ_INTEGER8_SCALAR, &
         BMPI_FILE_READ_REAL4_VEC_LC, BMPI_FILE_READ_REAL4_VEC_SC, &
         BMPI_FILE_READ_REAL4_SCALAR, &
         BMPI_FILE_READ_REAL8_VEC_LC, BMPI_FILE_READ_REAL8_VEC_SC, &
         BMPI_FILE_READ_REAL8_SCALAR
   end interface BMPI_FILE_READ

   !============================================================================
   ! WRITE
   !
   interface BMPI_FILE_WRITE
      MODULE PROCEDURE BMPI_FILE_WRITE_INTEGER4_VEC_LC, BMPI_FILE_WRITE_INTEGER4_VEC_SC, &
      BMPI_FILE_WRITE_INTEGER4_SCALAR, &
      BMPI_FILE_WRITE_INTEGER8_VEC_LC, BMPI_FILE_WRITE_INTEGER8_VEC_SC, &
      BMPI_FILE_WRITE_INTEGER8_SCALAR, &
      BMPI_FILE_WRITE_REAL4_VEC_LC, BMPI_FILE_WRITE_REAL4_VEC_SC, &
      BMPI_FILE_WRITE_REAL4_SCALAR, &
      BMPI_FILE_WRITE_REAL8_VEC_LC, BMPI_FILE_WRITE_REAL8_VEC_SC, &
      BMPI_FILE_WRITE_REAL8_SCALAR
   end interface BMPI_FILE_WRITE

   !============================================================================
   ! BCAST
   !
   interface BMPI_BCAST
      MODULE PROCEDURE BMPI_BCAST_CHARACTER, &
         BMPI_BCAST_LOGICAL_VEC, BMPI_BCAST_LOGICAL_SCALAR, &
         BMPI_BCAST_INTEGER4_VEC, BMPI_BCAST_INTEGER4_SCALAR, &
         BMPI_BCAST_INTEGER8_VEC, BMPI_BCAST_INTEGER8_SCALAR, &
         BMPI_BCAST_REAL4_VEC_SC, BMPI_BCAST_REAL4_VEC_LC, BMPI_BCAST_REAL4_SCALAR, &
         BMPI_BCAST_REAL8_VEC_SC, BMPI_BCAST_REAL8_VEC_LC, BMPI_BCAST_REAL8_SCALAR
   end interface BMPI_BCAST

   !============================================================================
   ! REDUCE
   !
   interface BMPI_REDUCE
      MODULE PROCEDURE BMPI_REDUCE_INTEGER4_VEC, BMPI_REDUCE_INTEGER4_SCALAR, &
         BMPI_REDUCE_INTEGER4_VECIP, BMPI_REDUCE_INTEGER4_SCALARIP, &
         BMPI_REDUCE_INTEGER8_VEC, BMPI_REDUCE_INTEGER8_SCALAR, &
         BMPI_REDUCE_INTEGER8_VECIP, BMPI_REDUCE_INTEGER8_SCALARIP, &
         BMPI_REDUCE_REAL4_VEC, BMPI_REDUCE_REAL4_SCALAR, &
         BMPI_REDUCE_REAL4_VECIP, BMPI_REDUCE_REAL4_SCALARIP, &
         BMPI_REDUCE_REAL8_VEC, BMPI_REDUCE_REAL8_SCALAR, &
         BMPI_REDUCE_REAL8_VECIP, BMPI_REDUCE_REAL8_SCALARIP
   end interface BMPI_REDUCE
   interface BMPI_ALLREDUCE
      MODULE PROCEDURE BMPI_ALLREDUCE_INTEGER4_VEC, BMPI_ALLREDUCE_INTEGER4_SCALAR, &
         BMPI_ALLREDUCE_INTEGER4_VECIP, BMPI_ALLREDUCE_INTEGER4_SCALARIP, &
         BMPI_ALLREDUCE_INTEGER8_VEC, BMPI_ALLREDUCE_INTEGER8_SCALAR, &
         BMPI_ALLREDUCE_INTEGER8_VECIP, BMPI_ALLREDUCE_INTEGER8_SCALARIP, &
         BMPI_ALLREDUCE_REAL4_VEC, BMPI_ALLREDUCE_REAL4_SCALAR, &
         BMPI_ALLREDUCE_REAL4_VECIP, BMPI_ALLREDUCE_REAL4_SCALARIP, &
         BMPI_ALLREDUCE_REAL8_VEC, BMPI_ALLREDUCE_REAL8_SCALAR, &
         BMPI_ALLREDUCE_REAL8_VECIP, BMPI_ALLREDUCE_REAL8_SCALARIP
   end interface BMPI_ALLREDUCE

   !============================================================================
   ! SEND/RECV
   !
   interface BMPI_SEND
      MODULE PROCEDURE BMPI_SEND_INTEGER4_VEC, BMPI_SEND_INTEGER4_SCALAR, &
         BMPI_SEND_INTEGER8_VEC, BMPI_SEND_INTEGER8_SCALAR, &
         BMPI_SEND_REAL4_VEC, BMPI_SEND_REAL4_SCALAR, &
         BMPI_SEND_REAL8_VEC, BMPI_SEND_REAL8_SCALAR
   end interface BMPI_SEND
   interface BMPI_RECV
      MODULE PROCEDURE BMPI_RECV_INTEGER4_VEC, BMPI_RECV_INTEGER4_SCALAR, &
         BMPI_RECV_INTEGER8_VEC, BMPI_RECV_INTEGER8_SCALAR, &
         BMPI_RECV_REAL4_VEC, BMPI_RECV_REAL4_SCALAR, &
         BMPI_RECV_REAL8_VEC, BMPI_RECV_REAL8_SCALAR
   end interface BMPI_RECV

   !============================================================================
   ! SCATTERV/GATHERV
   !
   interface BMPI_SCATTERV
      MODULE PROCEDURE BMPI_SCATTERV_INTEGER4_VEC, BMPI_GATHERV_INTEGER8_VEC, &
         BMPI_SCATTERV_REAL4_VEC, BMPI_GATHERV_REAL8_VEC
   end interface BMPI_SCATTERV
   interface BMPI_GATHERV
      MODULE PROCEDURE BMPI_GATHERV_INTEGER4_VEC, BMPI_GATHERV_INTEGER8_VEC, &
         BMPI_GATHERV_REAL4_VEC, BMPI_GATHERV_REAL8_VEC
   end interface BMPI_GATHERV
contains

DOUBLE PRECISION FUNCTION  BMPI_WTIME()
   implicit none
   interface
      DOUBLE PRECISION FUNCTION MPI_WTIME()
      end function MPI_WTIME
   end interface

   BMPI_WTIME = MPI_WTIME()
   RETURN
end FUNCTION BMPI_WTIME

subroutine BMPI_INIT(ierr)
   implicit none
   integer, intent(out) :: ierr
   call MPI_INIT(ierr)
end subroutine BMPI_INIT

subroutine BMPI_FINALIZE(ierr)
   implicit none
   integer, intent(out) :: ierr

   call MPI_FINALIZE(ierr)
end subroutine BMPI_FINALIZE

subroutine BMPI_COMM_RANK(comm, rank, ierr)
   implicit none
   integer, intent(in) :: comm
   integer, intent(out) :: rank
   integer, intent(out) :: ierr

   call MPI_COMM_RANK(comm, rank, ierr)
end subroutine BMPI_COMM_RANK

subroutine BMPI_TRACE0()
   use nodeinfo
   implicit none
   integer :: rank, ierr
   call MPI_COMM_RANK(icomm, rank, ierr)

   if(rank == 0) then
     ! call printstacktrace
      flush(6)
   endif
end subroutine BMPI_TRACE0

subroutine BMPI_COMM_SIZE(comm, csize, ierr)
   implicit none
   integer, intent(in) :: comm
   integer, intent(out) :: csize
   integer, intent(out) :: ierr

   call MPI_COMM_SIZE(comm, csize, ierr)
end subroutine BMPI_COMM_SIZE

subroutine BMPI_COMM_SPLIT(comm, color, key, newcomm, ierr)
   implicit none
   integer, intent(in) :: comm, color, key
   integer, intent(out) :: newcomm
   integer, intent(out) :: ierr

   call MPI_COMM_SPLIT(comm, color, key, newcomm, ierr)
end subroutine BMPI_COMM_SPLIT

! Collects physical name of processor for this MPI process
subroutine BMPI_GET_PROCESSOR_NAME(mpiprocname, plen, ierr)
   implicit none
   character (len=*), intent(out) :: mpiprocname
   integer, intent(out) :: plen ! returns size of string
   integer, intent(out) :: ierr

   call MPI_GET_PROCESSOR_NAME(mpiprocname, plen, ierr)
end subroutine BMPI_GET_PROCESSOR_NAME

subroutine BMPI_TYPE_SIZE(type, tsize, ierr)
   implicit none
   integer, intent(in) :: type
   integer, intent(out) :: tsize
   integer, intent(out) :: ierr

   call MPI_TYPE_SIZE(type, tsize, ierr)
end subroutine BMPI_TYPE_SIZE

subroutine BMPI_ABORT(comm, rc, ierr)
   implicit none
   integer, intent(in) :: comm, rc
   integer, intent(out) :: ierr

   call MPI_ABORT(comm, rc, ierr)
end subroutine BMPI_ABORT

subroutine BMPI_ERR(msg)
	implicit none
   character (len=*), intent(in) :: msg
   integer :: rank, ierr
   ! don't want to include other modules, get iproc manually
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	if(rank == 0) print *, msg
   ! identify location of error
 !  call printstacktrace
   call MPI_ABORT(MPI_COMM_WORLD, 53, ierr)
   stop 2 ! just in case
end subroutine BMPI_ERR

subroutine BMPI_BARRIER(comm, ierr)
   implicit none
   integer, intent(in) :: comm
   integer, intent(out) :: ierr

   call MPI_BARRIER(comm, ierr)
end subroutine BMPI_BARRIER

subroutine BMPI_FILE_OPEN(icomm, filename, mode, info, fn, ierr)
   implicit none
   integer, intent(in) :: icomm, mode, fn, info
   integer, intent(out) :: ierr
   character (len=*) :: filename

   call MPI_FILE_OPEN(icomm, filename, mode, info, fn, ierr)
end subroutine BMPI_FILE_OPEN


subroutine BMPI_FILE_CLOSE(fn, ierr)
   implicit none
   integer, intent(in) :: fn
   integer, intent(out) :: ierr

   call MPI_FILE_CLOSE(fn, ierr)
end subroutine BMPI_FILE_CLOSE

subroutine BMPI_FILE_SET_SIZE(fn, newsize, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind = MPI_OFFSET_KIND), intent(in) :: newsize
   integer, intent(out) :: ierr
   call MPI_FILE_SET_SIZE(fn, newsize, ierr)
end subroutine BMPI_FILE_SET_SIZE


subroutine BMPI_FILE_GET_POSITION(fn, offset, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind=MPI_OFFSET_KIND), intent(out) :: offset
   integer, intent(out) :: ierr

   call MPI_FILE_GET_POSITION(fn, offset, ierr)
end subroutine BMPI_FILE_GET_POSITION

subroutine BMPI_FILE_SEEK(fn, offset, whence, ierr)
   implicit none
   integer, intent(in) :: fn, whence
   integer(kind=MPI_OFFSET_KIND), intent(in) :: offset
   integer, intent(out) :: ierr

   call MPI_FILE_SEEK(fn, offset, whence, ierr)
end subroutine BMPI_FILE_SEEK

! Note:   We only want to use the I/O routines that
! use fixed sizes.  Otherwise our files might not be portable
! Still have to think about bit order
subroutine BMPI_FILE_WRITE_INTEGER4_VEC_LC(fn, buf, count8, status, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) ::  count8
   integer(kind=4), intent(in) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: ubnd, a, na, b, chunk
   integer :: csize ! for MPI call

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("BMPI_WRITE with count>size")
   chunk = bmpi_iosize / 4
   ubnd = UBOUND(buf,1, 8)
   a = LBOUND(buf,1, 8)
   do while(a <= ubnd)
      na = a + chunk
      b = na - 1
      if(b > ubnd) b = ubnd
      csize = int(b - a,4) + 1
      call MPI_FILE_WRITE(fn, buf(a:b), csize, MPI_INTEGER4, status, ierr)
      a = na
   end do
end subroutine BMPI_FILE_WRITE_INTEGER4_VEC_LC
subroutine BMPI_FILE_WRITE_INTEGER4_VEC_SC(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   integer(kind=4), intent(in) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: count8
   count8 = count
   call BMPI_FILE_WRITE_INTEGER4_VEC_LC(fn, buf, count8, status, ierr)
end subroutine BMPI_FILE_WRITE_INTEGER4_VEC_SC
subroutine BMPI_FILE_WRITE_INTEGER4_SCALAR(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   integer(kind=4), intent(in) :: buf
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_WRITE of scalar with count/=1")
   call MPI_FILE_WRITE(fn, buf, count, MPI_INTEGER4, status, ierr)
end subroutine BMPI_FILE_WRITE_INTEGER4_SCALAR
subroutine BMPI_FILE_WRITE_INTEGER8_VEC_LC(fn, buf, count8, status, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) :: count8
   integer(kind=8), intent(in) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: ubnd, a, na, b, chunk
   integer :: csize ! for MPI call

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("BMPI_WRITE with count>size")
   chunk = bmpi_iosize / 8
   ubnd = UBOUND(buf,1, 8)
   a = LBOUND(buf,1, 8)
   do while(a <= ubnd)
      na = a + chunk
      b = na - 1
      if(b > ubnd) b = ubnd
      csize = int(b - a,4) + 1
      call MPI_FILE_WRITE(fn, buf(a:b), csize, MPI_INTEGER8, status, ierr)
      a = na
   end do
end subroutine BMPI_FILE_WRITE_INTEGER8_VEC_LC
subroutine BMPI_FILE_WRITE_INTEGER8_VEC_SC(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   integer(kind=8), intent(in) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: count8

   count8 = count
   call BMPI_FILE_WRITE_INTEGER8_VEC_LC(fn, buf, count8, status, ierr)
end subroutine BMPI_FILE_WRITE_INTEGER8_VEC_SC
subroutine BMPI_FILE_WRITE_INTEGER8_SCALAR(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   integer(kind=8), intent(in) :: buf
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_WRITE of scalar with count/=1")
   call MPI_FILE_WRITE(fn, buf, count, MPI_INTEGER8, status, ierr)
end subroutine BMPI_FILE_WRITE_INTEGER8_SCALAR
subroutine BMPI_FILE_WRITE_REAL4_VEC_LC(fn, buf, count8, status, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) :: count8
   real(kind=4), intent(in) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: ubnd, a, na, b, chunk
   integer :: csize ! for MPI call

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("BMPI_WRITE with count>size")
   chunk = bmpi_iosize / 4
   ubnd = UBOUND(buf,1, 8)
   a = LBOUND(buf,1, 8)
   do while(a <= ubnd)
      na = a + chunk
      b = na - 1
      if(b > ubnd) b = ubnd
      csize = int(b - a,4) + 1
      call MPI_FILE_WRITE(fn, buf(a:b), csize, MPI_REAL4, status, ierr)
      a = na
   end do
end subroutine BMPI_FILE_WRITE_REAL4_VEC_LC
subroutine BMPI_FILE_WRITE_REAL4_VEC_SC(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   real(kind=4), intent(in) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: count8

   count8 = count
   call BMPI_FILE_WRITE_REAL4_VEC_LC(fn, buf, count8, status, ierr)
end subroutine BMPI_FILE_WRITE_REAL4_VEC_SC
subroutine BMPI_FILE_WRITE_REAL4_SCALAR(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   real(kind=4), intent(in) :: buf
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_WRITE of scalar with count/=1")
   call MPI_FILE_WRITE(fn, buf, count, MPI_REAL4, status, ierr)
end subroutine BMPI_FILE_WRITE_REAL4_SCALAR
subroutine BMPI_FILE_WRITE_REAL8_VEC_LC(fn, buf, count8, status, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) :: count8
   real(kind=8), intent(in) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: ubnd, a, na, b, chunk
   integer :: csize ! for MPI call

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("BMPI_WRITE with count>size")
   chunk = bmpi_iosize / 8
   ubnd = UBOUND(buf,1, 8)
   a = LBOUND(buf,1, 8)
   do while(a <= ubnd)
      na = a + chunk
      b = na - 1
      if(b > ubnd) b = ubnd
      csize = int(b - a,4) + 1
      call MPI_FILE_WRITE(fn, buf(a:b), csize, MPI_REAL8, status, ierr)
      a = na
   end do
end subroutine BMPI_FILE_WRITE_REAL8_VEC_LC
subroutine BMPI_FILE_WRITE_REAL8_VEC_SC(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn
   integer, intent(in) :: count
   real(kind=8), intent(in) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: count8
   count8 = count

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("BMPI_WRITE with count>size")
   call BMPI_FILE_WRITE_REAL8_VEC_LC(fn, buf, count8, status, ierr)
end subroutine BMPI_FILE_WRITE_REAL8_VEC_SC
subroutine BMPI_FILE_WRITE_REAL8_SCALAR(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   real(kind=8), intent(in) :: buf
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_WRITE of scalar with count/=1")
   call MPI_FILE_WRITE(fn, buf, count, MPI_REAL8, status, ierr)
end subroutine BMPI_FILE_WRITE_REAL8_SCALAR


!==================================================================
! READ routines
!
subroutine BMPI_FILE_READ_INTEGER4_VEC_LC(fn, buf, count8, status, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) :: count8
   integer(kind=4), intent(out) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: ubnd, a, na, b, chunk
   integer :: csize ! for MPI call

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("Read with count>size")
   chunk = bmpi_iosize / 4
   ubnd = UBOUND(buf,1, 8)
   a = LBOUND(buf,1, 8)
   do while(a <= ubnd)
      na = a + chunk
      b = na - 1
      if(b > ubnd) b = ubnd
      csize = int(b - a,4) + 1
      call MPI_FILE_READ(fn, buf(a:b), csize, MPI_INTEGER4, status, ierr)
      a = na
   end do
end subroutine BMPI_FILE_READ_INTEGER4_VEC_LC
subroutine BMPI_FILE_READ_INTEGER4_VEC_SC(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   integer(kind=4), intent(out) :: buf(:)
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: count8
   count8 = count
   call BMPI_FILE_READ_INTEGER4_VEC_LC(fn, buf, count8, status, ierr)
end subroutine BMPI_FILE_READ_INTEGER4_VEC_SC
subroutine BMPI_FILE_READ_INTEGER4_SCALAR(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   integer(kind=4), intent(out) :: buf
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("Read into scalar with count/=1")
   call MPI_FILE_READ(fn, buf, count, MPI_INTEGER4, status, ierr)
end subroutine BMPI_FILE_READ_INTEGER4_SCALAR
subroutine BMPI_FILE_READ_INTEGER8_VEC_LC(fn, buf, count8, status, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) :: count8
   integer(kind=8), intent(out) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: ubnd, a, na, b, chunk
   integer :: csize ! for MPI call

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("Read with count>size")
   chunk = bmpi_iosize / 8
   ubnd = UBOUND(buf,1, 8)
   a = LBOUND(buf,1, 8)
   do while(a <= ubnd)
      na = a + chunk
      b = na - 1
      if(b > ubnd) b = ubnd
      csize = int(b - a,4) + 1
      call MPI_FILE_READ(fn, buf(a:b), csize, MPI_INTEGER8, status, ierr)
      a = na
   end do
end subroutine BMPI_FILE_READ_INTEGER8_VEC_LC
subroutine BMPI_FILE_READ_INTEGER8_VEC_SC(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   integer(kind=8), intent(out) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: count8
   count8 = count
   call BMPI_FILE_READ_INTEGER8_VEC_LC(fn, buf, count8, status, ierr)
end subroutine BMPI_FILE_READ_INTEGER8_VEC_SC
subroutine BMPI_FILE_READ_INTEGER8_SCALAR(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   integer(kind=8), intent(out) :: buf
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("Read into scalar with count/=1")
   call MPI_FILE_READ(fn, buf, count, MPI_INTEGER8, status, ierr)
end subroutine BMPI_FILE_READ_INTEGER8_SCALAR
subroutine BMPI_FILE_READ_REAL4_VEC_LC(fn, buf, count8, status, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) :: count8
   real(kind=4), intent(out) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: ubnd, a, na, b, chunk
   integer :: csize ! for MPI call

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("Read with count>size")
   chunk = bmpi_iosize / 4
   ubnd = UBOUND(buf,1, 8)
   a = LBOUND(buf,1, 8)
   do while(a <= ubnd)
      na = a + chunk
      b = na - 1
      if(b > ubnd) b = ubnd
      csize = int(b - a,4) + 1
      call MPI_FILE_READ(fn, buf(a:b), csize, MPI_REAL4, status, ierr)
      a = na
   end do
end subroutine BMPI_FILE_READ_REAL4_VEC_LC
subroutine BMPI_FILE_READ_REAL4_VEC_SC(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   real(kind=4), intent(out) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: count8
   count8 = count

   call BMPI_FILE_READ_REAL4_VEC_LC(fn, buf, count8, status, ierr)
end subroutine BMPI_FILE_READ_REAL4_VEC_SC
subroutine BMPI_FILE_READ_REAL4_SCALAR(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   real(kind=4), intent(out) :: buf
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("Read into scalar with count/=1")
   call MPI_FILE_READ(fn, buf, count, MPI_REAL4, status, ierr)
end subroutine BMPI_FILE_READ_REAL4_SCALAR
subroutine BMPI_FILE_READ_REAL8_VEC_LC(fn, buf, count8, status, ierr)
   implicit none
   integer, intent(in) :: fn
   integer(kind=8), intent(in) :: count8
   real(kind=8), intent(out) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: ubnd, a, na, b, chunk
   integer :: csize ! for MPI call

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("Read with count>size")
   chunk = bmpi_iosize / 8
   ubnd = UBOUND(buf,1, 8)
   a = LBOUND(buf,1, 8)
   do while(a <= ubnd)
      na = a + chunk
      b = na - 1
      if(b > ubnd) b = ubnd
      csize = int(b - a,4) + 1
      call MPI_FILE_READ(fn, buf(a:b), csize, MPI_REAL8, status, ierr)
      a = na
   end do
end subroutine BMPI_FILE_READ_REAL8_VEC_LC
subroutine BMPI_FILE_READ_REAL8_VEC_SC(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   real(kind=8), intent(out) :: buf(:)
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr
   integer(kind=8) :: count8
   count8 = count

   call BMPI_FILE_READ_REAL8_VEC_LC(fn, buf, count8, status, ierr)
end subroutine BMPI_FILE_READ_REAL8_VEC_SC
subroutine BMPI_FILE_READ_REAL8_SCALAR(fn, buf, count, status, ierr)
   implicit none
   integer, intent(in) :: fn, count
   real(kind=8), intent(out) :: buf
   ! looks to be default integer type
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("Read into scalar with count/=1")
   call MPI_FILE_READ(fn, buf, count, MPI_REAL8, status, ierr)
end subroutine BMPI_FILE_READ_REAL8_SCALAR

!================================================================
! BCAST ROUTINES
!

! no need for vec/scalar for character type
subroutine BMPI_BCAST_CHARACTER(buf, count, root, comm, ierr)
   implicit none
   CHARACTER (LEN=*), intent(inout) :: buf
   integer, intent(in):: count, root, comm
   integer, intent(out):: ierr

! KSM: not sure if this is right
   if(count > LEN(buf)) call BMPI_ERR("BMPI_BCAST with count>len")
   call MPI_BCAST(buf, count, MPI_CHARACTER, root, comm, ierr)
end subroutine BMPI_BCAST_CHARACTER

subroutine BMPI_BCAST_LOGICAL_VEC(buf, count, root, comm, ierr)
   implicit none
   logical, intent(inout) :: buf(:)
   integer, intent(in):: count, root, comm
   integer, intent(out):: ierr

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_BCAST with count>size")
   call MPI_BCAST(buf, count, MPI_LOGICAL, root, comm, ierr)
end subroutine BMPI_BCAST_LOGICAL_VEC
subroutine BMPI_BCAST_LOGICAL_SCALAR(buf, count, root, comm, ierr)
   implicit none
   logical, intent(inout) :: buf
   integer, intent(in):: count, root, comm
   integer, intent(out):: ierr

   if(count /= 1) call BMPI_ERR("BMPI_BCAST of scalar with count/=1")
   call MPI_BCAST(buf, count, MPI_LOGICAL, root, comm, ierr)
end subroutine BMPI_BCAST_LOGICAL_SCALAR

subroutine BMPI_BCAST_INTEGER4_VEC(buf, count, root, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: buf(:)
   integer, intent(in):: count, root, comm
   integer, intent(out):: ierr
   ! locals
   integer(kind=8) :: s, e, ns, slicesize, lim
   integer :: castsize  ! must be plain integer to match MPI interface

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_BCAST with count>size")
   slicesize = bmpi_maxbcastsize / 4
   s = LBOUND(buf, 1, 8)
   lim = s + count - 1
   do while(s <= lim)
      ns = s + slicesize
      e = ns - 1
      if(e > lim) e = lim
      ! convert down to 4 byte count
      castsize = int(e - s,4) + 1
      call MPI_BCAST(buf(s:e), castsize, MPI_INTEGER4, root, comm, ierr)
      s = ns
   end do
end subroutine BMPI_BCAST_INTEGER4_VEC
subroutine BMPI_BCAST_INTEGER4_SCALAR(buf, count, root, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: buf
   integer, intent(in):: count, root, comm
   integer, intent(out):: ierr

   if(count /= 1) call BMPI_ERR("BMPI_BCAST of scalar with count/=1")
   call MPI_BCAST(buf, count, MPI_INTEGER4, root, comm, ierr)
end subroutine BMPI_BCAST_INTEGER4_SCALAR

subroutine BMPI_BCAST_INTEGER8_VEC(buf, count, root, comm, ierr)
   implicit none
   integer (kind=8), intent(inout) :: buf(:)
   integer, intent(in):: count, root, comm
   integer, intent(out):: ierr
   ! locals
   integer(kind=8) :: s, e, ns, slicesize, lim
   integer :: castsize  ! must be plain integer to match MPI interface

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_BCAST with count>size")
   slicesize = bmpi_maxbcastsize / 8
   s = LBOUND(buf, 1, 8)
   lim = s + count - 1
   do while(s <= lim)
      ns = s + slicesize
      e = ns - 1
      if(e > lim) e = lim
      ! convert down to 4 byte count
      castsize = int(e - s,4) + 1
      call MPI_BCAST(buf(s:e), castsize, MPI_INTEGER8, root, comm, ierr)
      s = ns
   end do
end subroutine BMPI_BCAST_INTEGER8_VEC
subroutine BMPI_BCAST_INTEGER8_SCALAR(buf, count, root, comm, ierr)
   implicit none
   integer (kind=8), intent(inout) :: buf
   integer, intent(in):: count, root, comm
   integer, intent(out):: ierr

   if(count /= 1) call BMPI_ERR("BMPI_BCAST of scalar with count/=1")
   call MPI_BCAST(buf, count, MPI_INTEGER8, root, comm, ierr)
end subroutine BMPI_BCAST_INTEGER8_SCALAR

subroutine BMPI_BCAST_REAL4_VEC_LC(buf, count8, root, comm, ierr)
   implicit none
   real(kind=4), intent(inout)  :: buf(:)
   integer, intent(in):: root, comm
   integer(kind=8), intent(in) :: count8
   integer, intent(out):: ierr
   integer(kind=8) :: s, e, ns, slicesize, lim
   integer :: castsize  ! must be plain integer to match MPI interface

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("BMPI_BCAST with count>size")
   slicesize = bmpi_maxbcastsize / 4
   s = LBOUND(buf, 1, 8)
   lim = s + count8 - 1
   do while(s <= lim)
      ns = s + slicesize
      e = ns - 1
      if(e > lim) e = lim
      ! cut down to 4 bytes for MPI interface
      castsize = int(e - s,4) + 1
      call MPI_BCAST(buf(s:e), castsize, MPI_REAL4, root, comm, ierr)
      s = ns
   end do
end subroutine BMPI_BCAST_REAL4_VEC_LC

subroutine BMPI_BCAST_REAL4_VEC_SC(buf, count, root, comm, ierr)
   implicit none
   real(kind=4), intent(inout)  :: buf(:)
   integer, intent(in):: root, comm
   integer, intent(in) :: count
   integer, intent(out):: ierr
   integer(kind=8) :: count8

   count8 = count
   call BMPI_BCAST_REAL4_VEC_LC(buf, count8, root, comm, ierr)
end subroutine BMPI_BCAST_REAL4_VEC_SC
! should add count8 version
subroutine BMPI_BCAST_REAL4_SCALAR(buf, count, root, comm, ierr)
   implicit none
   real (kind=4), intent(inout)  :: buf
   integer, intent(in):: count, root, comm
   integer, intent(out):: ierr

   if(count /= 1) call BMPI_ERR("BMPI_BCAST of scalar with count/=1")
   call MPI_BCAST(buf, count, MPI_REAL4, root, comm, ierr)
end subroutine BMPI_BCAST_REAL4_SCALAR

! Implemented chunking to avoid bug with large BCAST
subroutine BMPI_BCAST_REAL8_VEC_LC(buf, count8, root, comm, ierr)
   implicit none
   real (kind=8), intent(inout)  :: buf(:)
   integer, intent(in):: root, comm
   integer(kind=8),  intent(in) :: count8
   integer, intent(out):: ierr
   integer(kind=8) :: s, e, ns, slicesize, lim
   integer :: castsize  ! must be plain integer to match MPI interface

   if(count8 > size(buf, 1, 8)) call BMPI_ERR("BMPI_BCAST with count>size")
   slicesize = bmpi_maxbcastsize / 8
   s = LBOUND(buf, 1, 8)
   lim = s + count8 - 1
   do while(s <= lim)
      ns = s + slicesize
      e = ns - 1
      if(e > lim) e = lim
      ! convert down to 4 byte count
      castsize = int(e - s,4) + 1
      call MPI_BCAST(buf(s:e), castsize, MPI_REAL8, root, comm, ierr)
      s = ns
   end do
end subroutine BMPI_BCAST_REAL8_VEC_LC

subroutine BMPI_BCAST_REAL8_VEC_SC(buf, count, root, comm, ierr)
   implicit none
   real (kind=8), intent(inout)  :: buf(:)
   integer, intent(in):: root, comm
   integer,  intent(in) :: count
   integer, intent(out):: ierr
   integer(kind=8) :: count8
   count8 = count
   call BMPI_BCAST_REAL8_VEC_LC(buf, count8, root, comm, ierr)
end subroutine BMPI_BCAST_REAL8_VEC_SC
subroutine BMPI_BCAST_REAL8_SCALAR(buf, count, root, comm, ierr)
   implicit none
   real (kind=8), intent(inout)  :: buf
   integer, intent(in):: count, root, comm
   integer, intent(out):: ierr

   if(count /= 1) call BMPI_ERR("BMPI_BCAST of scalar with count/=1")
   call MPI_BCAST(buf, count, MPI_REAL8, root, comm, ierr)
end subroutine BMPI_BCAST_REAL8_SCALAR

!================================================================
! REDUCE ROUTINES

! INTEGER4 ROUTINES
subroutine BMPI_REDUCE_INTEGER4_VEC(sendbuf, recvbuf, count, op, root, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: sendbuf(:), recvbuf(:)
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer(kind=8) :: ss, se, rs, re, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_REDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   slicesize = bmpi_maxbcastsize / 4
   ss = LBOUND(sendbuf, 1, 8)
   rs = LBOUND(recvbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      re = rs + slicesize
      over = se-lim
      if(over > 0) then
         se = se - over
         re = re - over
      end if
      scnt = int(se - ss,4) + 1
      call MPI_REDUCE(sendbuf(ss:se), recvbuf(rs:re), scnt, MPI_INTEGER4, op, root, comm, ierr)
      ss = se + 1
      rs = re + 1
   end do
end subroutine BMPI_REDUCE_INTEGER4_VEC
subroutine BMPI_REDUCE_INTEGER4_SCALAR(sendbuf, recvbuf, count, op, root, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: sendbuf, recvbuf
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   call MPI_REDUCE(sendbuf, recvbuf, count, MPI_INTEGER4, op, root, comm, ierr)
end subroutine BMPI_REDUCE_INTEGER4_SCALAR
! In place versions
! Note that MPI_IN_PLACE must replace the sendbuf.  Then
! all processes use the recvbuf as the source data
subroutine BMPI_REDUCE_INTEGER4_VECIP(sendbuf,  count, op, root, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: sendbuf(:)
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, slicesize, over, lim
   integer :: scnt
   integer :: rank
   integer :: rierr
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_REDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_COMM_RANK(comm, rank, rierr)
   slicesize = bmpi_maxbcastsize / 4
   ss = LBOUND(sendbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      over = se-lim
      if(over > 0) se = se - over
      scnt = int(se - ss,4) + 1
      ! this is truely awful - can only use MPI_IN_PLACE in root
      if(rank == root) then
         call MPI_REDUCE(MPI_IN_PLACE, sendbuf(ss:se), scnt, MPI_INTEGER4, op, root, comm, ierr)
      else
         call MPI_REDUCE(sendbuf(ss:se), sendbuf(ss:se), scnt, MPI_INTEGER4, op, root, comm, ierr)
      endif
      ss = se + 1
   end do
end subroutine BMPI_REDUCE_INTEGER4_VECIP
subroutine BMPI_REDUCE_INTEGER4_SCALARIP(sendbuf, count, op, root, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: sendbuf
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: rank
   integer :: rierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_COMM_RANK(comm, rank, rierr)
   if(rank == root) then
      call MPI_REDUCE(MPI_IN_PLACE, sendbuf, count, MPI_INTEGER4, op, root, comm, ierr)
   else
      call MPI_REDUCE(sendbuf, sendbuf, count, MPI_INTEGER4, op, root, comm, ierr)
   end if
end subroutine BMPI_REDUCE_INTEGER4_SCALARIP

! INTEGER8 routines
subroutine BMPI_REDUCE_INTEGER8_VEC(sendbuf, recvbuf, count, op, root, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: sendbuf(:), recvbuf(:)
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, rs, re, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_REDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   slicesize = bmpi_maxbcastsize / 8
   ss = LBOUND(sendbuf, 1, 8)
   rs = LBOUND(recvbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      re = rs + slicesize
      over = se-lim
      if(over > 0) then
         se = se - over
         re = re - over
      end if
      scnt = int(se - ss,4) + 1
      call MPI_REDUCE(sendbuf(ss:se), recvbuf(rs:re), scnt, MPI_INTEGER8, op, root, comm, ierr)
      ss = se + 1
      rs = re + 1
   end do
end subroutine BMPI_REDUCE_INTEGER8_VEC
subroutine BMPI_REDUCE_INTEGER8_SCALAR(sendbuf, recvbuf, count, op, root, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: sendbuf, recvbuf
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   call MPI_REDUCE(sendbuf, recvbuf, count, MPI_INTEGER8, op, root, comm, ierr)
end subroutine BMPI_REDUCE_INTEGER8_SCALAR
! In place versions
subroutine BMPI_REDUCE_INTEGER8_VECIP(sendbuf,  count, op, root, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: sendbuf(:)
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, slicesize, over, lim
   integer :: scnt
   integer :: rank
   integer :: rierr
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_REDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_COMM_RANK(comm, rank, rierr)
   slicesize = bmpi_maxbcastsize / 8
   ss = LBOUND(sendbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      over = se-lim
      if(over > 0) se = se - over
      scnt = int(se - ss,4) + 1
      if(rank == root) then
         call MPI_REDUCE(MPI_IN_PLACE, sendbuf(ss:se), scnt, MPI_INTEGER8, op, root, comm, ierr)
      else
         call MPI_REDUCE(sendbuf(ss:se), sendbuf(ss:se), scnt, MPI_INTEGER8, op, root, comm, ierr)
      end if
      ss = se + 1
   end do
end subroutine BMPI_REDUCE_INTEGER8_VECIP
subroutine BMPI_REDUCE_INTEGER8_SCALARIP(sendbuf, count, op, root, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: sendbuf
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: rank
   integer :: rierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_COMM_RANK(comm, rank, rierr)
   if(rank == root) then
      call MPI_REDUCE(MPI_IN_PLACE, sendbuf, count, MPI_INTEGER8, op, root, comm, ierr)
   else
      call MPI_REDUCE(sendbuf, sendbuf, count, MPI_INTEGER8, op, root, comm, ierr)
   end if
end subroutine BMPI_REDUCE_INTEGER8_SCALARIP

! REAL4 routines
subroutine BMPI_REDUCE_REAL4_VEC(sendbuf, recvbuf, count, op, root, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: sendbuf(:), recvbuf(:)
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, rs, re, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_REDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   slicesize = bmpi_maxbcastsize / 4
   ss = LBOUND(sendbuf, 1, 8)
   rs = LBOUND(recvbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      re = rs + slicesize
      over = se-lim
      if(over > 0) then
         se = se - over
         re = re - over
      end if
      scnt = int(se - ss,4) + 1
      call MPI_REDUCE(sendbuf(ss:se), recvbuf(rs:re), scnt, MPI_REAL4, op, root, comm, ierr)
      ss = se + 1
      rs = re + 1
   end do
end subroutine BMPI_REDUCE_REAL4_VEC
subroutine BMPI_REDUCE_REAL4_SCALAR(sendbuf, recvbuf, count, op, root, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: sendbuf, recvbuf
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL4, op, root, comm, ierr)
end subroutine BMPI_REDUCE_REAL4_SCALAR
! In place versions
subroutine BMPI_REDUCE_REAL4_VECIP(sendbuf,  count, op, root, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: sendbuf(:)
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, slicesize, over, lim
   integer :: scnt
   integer :: rank
   integer :: rierr
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_REDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_COMM_RANK(comm, rank, rierr)
   slicesize = bmpi_maxbcastsize / 4
   ss = LBOUND(sendbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      over = se-lim
      if(over > 0) se = se - over
      scnt = int(se - ss,4) + 1
      if(rank == root) then
         call MPI_REDUCE(MPI_IN_PLACE, sendbuf(ss:se), scnt, MPI_REAL4, op, root, comm, ierr)
      else
         call MPI_REDUCE(sendbuf(ss:se), sendbuf(ss:se), scnt, MPI_REAL4, op, root, comm, ierr)
      end if
      ss = se + 1
   end do
end subroutine BMPI_REDUCE_REAL4_VECIP
subroutine BMPI_REDUCE_REAL4_SCALARIP(sendbuf, count, op, root, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: sendbuf
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: rank
   integer :: rierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_COMM_RANK(comm, rank, rierr)
   if(rank == root) then
      call MPI_REDUCE(MPI_IN_PLACE, sendbuf, count, MPI_REAL4, op, root, comm, ierr)
   else
      call MPI_REDUCE(sendbuf, sendbuf, count, MPI_REAL4, op, root, comm, ierr)
   end if
end subroutine BMPI_REDUCE_REAL4_SCALARIP

! REAL8 routines
subroutine BMPI_REDUCE_REAL8_VEC(sendbuf, recvbuf, count, op, root, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: sendbuf(:), recvbuf(:)
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, rs, re, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_REDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   slicesize = bmpi_maxbcastsize / 8
   ss = LBOUND(sendbuf, 1, 8)
   rs = LBOUND(recvbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      re = rs + slicesize
      over = se-lim
      if(over > 0) then
         se = se - over
         re = re - over
      end if
      scnt = int(se - ss,4) + 1
      call MPI_REDUCE(sendbuf(ss:se), recvbuf(rs:re), scnt, MPI_REAL8, op, root, comm, ierr)
      ss = se + 1
      rs = re + 1
   end do
end subroutine BMPI_REDUCE_REAL8_VEC
subroutine BMPI_REDUCE_REAL8_SCALAR(sendbuf, recvbuf, count, op, root, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: sendbuf, recvbuf
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   call MPI_REDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, root, comm, ierr)
end subroutine BMPI_REDUCE_REAL8_SCALAR
! In place versions
subroutine BMPI_REDUCE_REAL8_VECIP(sendbuf,  count, op, root, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: sendbuf(:)
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, slicesize, over, lim
   integer :: scnt
   integer :: rank
   integer :: rierr
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_REDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_COMM_RANK(comm, rank, rierr)
   
   slicesize = bmpi_maxbcastsize / 8
   ss = LBOUND(sendbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      over = se-lim
      if(over > 0) se = se - over
      scnt = int(se - ss,4) + 1
      if(rank == root) then
         call MPI_REDUCE(MPI_IN_PLACE, sendbuf(ss:se), scnt, MPI_REAL8, op, root, comm, ierr)
      else
         call MPI_REDUCE(sendbuf(ss:se), sendbuf(ss:se), scnt, MPI_REAL8, op, root, comm, ierr)
      end if
 !     call MPI_REDUCE(MPI_IN_PLACE, sendbuf(ss:se), scnt, MPI_REAL8, op, root, comm, ierr)
      ss = se + 1
   end do
end subroutine BMPI_REDUCE_REAL8_VECIP
subroutine BMPI_REDUCE_REAL8_SCALARIP(sendbuf, count, op, root, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: sendbuf
   integer, intent(in) :: op, root, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: rank
   integer :: rierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_COMM_RANK(comm, rank, rierr)
   if(rank == root) then
      call MPI_REDUCE(MPI_IN_PLACE, sendbuf, count, MPI_REAL8, op, root, comm, ierr)
   else
      call MPI_REDUCE(sendbuf, sendbuf, count, MPI_REAL8, op, root, comm, ierr)
   end if
end subroutine BMPI_REDUCE_REAL8_SCALARIP

!================================================================
! ALLREDUCE ROUTINES
! These are almost the same is REDUCE.  The difference is that
! they don't specify a root since the result is delivered to all ranks

! INTEGER4 ROUTINES
subroutine BMPI_ALLREDUCE_INTEGER4_VEC(sendbuf, recvbuf, count, op, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: sendbuf(:), recvbuf(:)
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer(kind=8) :: ss, se, rs, re, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_ALLREDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   slicesize = bmpi_maxbcastsize / 4
   ss = LBOUND(sendbuf, 1, 8)
   rs = LBOUND(recvbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      re = rs + slicesize
      over = se-lim
      if(over > 0) then
         se = se - over
         re = re - over
      end if
      scnt = int(se - ss,4) + 1
      call MPI_ALLREDUCE(sendbuf(ss:se), recvbuf(rs:re), scnt, MPI_INTEGER4, op, comm, ierr)
      ss = se + 1
      rs = re + 1
   end do
end subroutine BMPI_ALLREDUCE_INTEGER4_VEC
subroutine BMPI_ALLREDUCE_INTEGER4_SCALAR(sendbuf, recvbuf, count, op, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: sendbuf, recvbuf
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_INTEGER4, op, comm, ierr)
end subroutine BMPI_ALLREDUCE_INTEGER4_SCALAR
! In place versions
subroutine BMPI_ALLREDUCE_INTEGER4_VECIP(sendbuf,  count, op, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: sendbuf(:)
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_ALLREDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   slicesize = bmpi_maxbcastsize / 4
   ss = LBOUND(sendbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      over = se-lim
      if(over > 0) se = se - over
      scnt = int(se - ss,4) + 1
      call MPI_ALLREDUCE(MPI_IN_PLACE, sendbuf(ss:se), scnt, MPI_INTEGER4, op, comm, ierr)
      ss = se + 1
   end do
end subroutine BMPI_ALLREDUCE_INTEGER4_VECIP
subroutine BMPI_ALLREDUCE_INTEGER4_SCALARIP(sendbuf, count, op, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: sendbuf
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_ALLREDUCE(MPI_IN_PLACE, sendbuf, count, MPI_INTEGER4, op, comm, ierr)
end subroutine BMPI_ALLREDUCE_INTEGER4_SCALARIP

! INTEGER8 routines
subroutine BMPI_ALLREDUCE_INTEGER8_VEC(sendbuf, recvbuf, count, op, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: sendbuf(:), recvbuf(:)
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, rs, re, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_ALLREDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   slicesize = bmpi_maxbcastsize / 8
   ss = LBOUND(sendbuf, 1, 8)
   rs = LBOUND(recvbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      re = rs + slicesize
      over = se-lim
      if(over > 0) then
         se = se - over
         re = re - over
      end if
      scnt = int(se - ss,4) + 1
      call MPI_ALLREDUCE(sendbuf(ss:se), recvbuf(rs:re), scnt, MPI_INTEGER8, op, comm, ierr)
      ss = se + 1
      rs = re + 1
   end do
end subroutine BMPI_ALLREDUCE_INTEGER8_VEC
subroutine BMPI_ALLREDUCE_INTEGER8_SCALAR(sendbuf, recvbuf, count, op, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: sendbuf, recvbuf
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_INTEGER8, op, comm, ierr)
end subroutine BMPI_ALLREDUCE_INTEGER8_SCALAR
! In place versions
subroutine BMPI_ALLREDUCE_INTEGER8_VECIP(sendbuf,  count, op, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: sendbuf(:)
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_ALLREDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   slicesize = bmpi_maxbcastsize / 8
   ss = LBOUND(sendbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      over = se-lim
      if(over > 0) se = se - over
      scnt = int(se - ss,4) + 1
      call MPI_ALLREDUCE(MPI_IN_PLACE, sendbuf(ss:se), scnt, MPI_INTEGER8, op, comm, ierr)
      ss = se + 1
   end do
end subroutine BMPI_ALLREDUCE_INTEGER8_VECIP
subroutine BMPI_ALLREDUCE_INTEGER8_SCALARIP(sendbuf, count, op, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: sendbuf
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_ALLREDUCE(MPI_IN_PLACE, sendbuf, count, MPI_INTEGER8, op, comm, ierr)
end subroutine BMPI_ALLREDUCE_INTEGER8_SCALARIP

! REAL4 routines
subroutine BMPI_ALLREDUCE_REAL4_VEC(sendbuf, recvbuf, count, op, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: sendbuf(:), recvbuf(:)
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, rs, re, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_ALLREDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   slicesize = bmpi_maxbcastsize / 4
   ss = LBOUND(sendbuf, 1, 8)
   rs = LBOUND(recvbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      re = rs + slicesize
      over = se-lim
      if(over > 0) then
         se = se - over
         re = re - over
      end if
      scnt = int(se - ss,4) + 1
      call MPI_ALLREDUCE(sendbuf(ss:se), recvbuf(rs:re), scnt, MPI_REAL4, op, comm, ierr)
      ss = se + 1
      rs = re + 1
   end do
end subroutine BMPI_ALLREDUCE_REAL4_VEC
subroutine BMPI_ALLREDUCE_REAL4_SCALAR(sendbuf, recvbuf, count, op, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: sendbuf, recvbuf
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL4, op, comm, ierr)
end subroutine BMPI_ALLREDUCE_REAL4_SCALAR
! In place versions
subroutine BMPI_ALLREDUCE_REAL4_VECIP(sendbuf,  count, op, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: sendbuf(:)
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_ALLREDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   slicesize = bmpi_maxbcastsize / 4
   ss = LBOUND(sendbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      over = se-lim
      if(over > 0) se = se - over
      scnt = int(se - ss,4) + 1
      call MPI_ALLREDUCE(MPI_IN_PLACE, sendbuf(ss:se), scnt, MPI_REAL4, op, comm, ierr)
      ss = se + 1
   end do
end subroutine BMPI_ALLREDUCE_REAL4_VECIP
subroutine BMPI_ALLREDUCE_REAL4_SCALARIP(sendbuf, count, op, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: sendbuf
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_ALLREDUCE(MPI_IN_PLACE, sendbuf, count, MPI_REAL4, op, comm, ierr)
end subroutine BMPI_ALLREDUCE_REAL4_SCALARIP

! REAL8 routines
subroutine BMPI_ALLREDUCE_REAL8_VEC(sendbuf, recvbuf, count, op, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: sendbuf(:), recvbuf(:)
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, rs, re, slicesize, over, lim
   integer :: scnt
   integer :: csize

   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_ALLREDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   slicesize = bmpi_maxbcastsize / 8
   ss = LBOUND(sendbuf, 1, 8)
   rs = LBOUND(recvbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      re = rs + slicesize
      over = se-lim
      if(over > 0) then
         se = se - over
         re = re - over
      end if
      scnt = int(se - ss,4) + 1
      call MPI_ALLREDUCE(sendbuf(ss:se), recvbuf(rs:re), scnt, MPI_REAL8, op, comm, ierr)
      ss = se + 1
      rs = re + 1
   end do
end subroutine BMPI_ALLREDUCE_REAL8_VEC
subroutine BMPI_ALLREDUCE_REAL8_SCALAR(sendbuf, recvbuf, count, op, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: sendbuf, recvbuf
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      recvbuf = sendbuf
      return
   end if
   call MPI_ALLREDUCE(sendbuf, recvbuf, count, MPI_REAL8, op, comm, ierr)
end subroutine BMPI_ALLREDUCE_REAL8_SCALAR
! In place versions
subroutine BMPI_ALLREDUCE_REAL8_VECIP(sendbuf,  count, op, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: sendbuf(:)
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   ! locals
   integer(kind=8) :: ss, se, slicesize, over, lim
   integer :: scnt
   integer :: csize


   if(count > size(sendbuf, 1, 8)) call BMPI_ERR("BMPI_ALLREDUCE count>size")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   slicesize = bmpi_maxbcastsize / 8
   ss = LBOUND(sendbuf, 1, 8)
   lim = ss + count - 1
   do while(ss <= lim)
      se = ss + slicesize
      over = se-lim
      if(over > 0) se = se - over
      scnt = int(se - ss,4) + 1
      call MPI_ALLREDUCE(MPI_IN_PLACE, sendbuf(ss:se), scnt, MPI_REAL8, op, comm, ierr)
      ss = se + 1
   end do
end subroutine BMPI_ALLREDUCE_REAL8_VECIP
subroutine BMPI_ALLREDUCE_REAL8_SCALARIP(sendbuf, count, op, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: sendbuf
   integer, intent(in) :: op, comm
   integer, intent(in) :: count
   integer, intent(out) :: ierr
   integer :: csize

   if(count /= 1) call BMPI_ERR("REDUCE of scalar with count/=1")
   call MPI_COMM_SIZE(comm, csize, ierr)
   if(csize == 1) then
      return
   end if
   call MPI_ALLREDUCE(MPI_IN_PLACE, sendbuf, count, MPI_REAL8, op, comm, ierr)
end subroutine BMPI_ALLREDUCE_REAL8_SCALARIP

!==================================================================
! SEND

subroutine BMPI_SEND_INTEGER4_VEC(buf, count, dest, tag, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: buf(:)
   integer, intent(in) :: count
   integer, intent(in) ::  dest, tag, comm
   integer, intent(out) :: ierr

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_SEND  count > size")
   call MPI_SEND(buf, count, MPI_INTEGER4, dest, tag, comm, ierr)
end subroutine BMPI_SEND_INTEGER4_VEC
subroutine BMPI_SEND_INTEGER4_SCALAR(buf, count, dest, tag, comm, ierr)
   implicit none
   integer(kind=4), intent(inout) :: buf
   integer, intent(in) :: count
   integer, intent(in) ::  dest, tag, comm
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_SEND scalar  count /= 1")
   call MPI_SEND(buf, count, MPI_INTEGER4, dest, tag, comm, ierr)
end subroutine BMPI_SEND_INTEGER4_SCALAR
subroutine BMPI_SEND_INTEGER8_VEC(buf, count, dest, tag, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: buf(:)
   integer, intent(in) :: count
   integer, intent(in) ::  dest, tag, comm
   integer, intent(out) :: ierr

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_SEND  count > size")
   call MPI_SEND(buf, count, MPI_INTEGER8, dest, tag, comm, ierr)
end subroutine BMPI_SEND_INTEGER8_VEC
subroutine BMPI_SEND_INTEGER8_SCALAR(buf, count, dest, tag, comm, ierr)
   implicit none
   integer(kind=8), intent(inout) :: buf
   integer, intent(in) :: count
   integer, intent(in) ::  dest, tag, comm
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_SEND scalar  count /= 1")
   call MPI_SEND(buf, count, MPI_INTEGER8, dest, tag, comm, ierr)
end subroutine BMPI_SEND_INTEGER8_SCALAR
subroutine BMPI_SEND_REAL4_VEC(buf, count, dest, tag, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: buf(:)
   integer, intent(in) :: count
   integer, intent(in) ::  dest, tag, comm
   integer, intent(out) :: ierr

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_SEND  count > size")
   call MPI_SEND(buf, count, MPI_REAL4, dest, tag, comm, ierr)
end subroutine BMPI_SEND_REAL4_VEC
subroutine BMPI_SEND_REAL4_SCALAR(buf, count, dest, tag, comm, ierr)
   implicit none
   real(kind=4), intent(inout) :: buf
   integer, intent(in) :: count
   integer, intent(in) ::  dest, tag, comm
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_SEND scalar  count /= 1")
   call MPI_SEND(buf, count, MPI_REAL4, dest, tag, comm, ierr)
end subroutine BMPI_SEND_REAL4_SCALAR
subroutine BMPI_SEND_REAL8_VEC(buf, count, dest, tag, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: buf(:)
   integer, intent(in) :: count
   integer, intent(in) ::  dest, tag, comm
   integer, intent(out) :: ierr

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_SEND  count > size")
   call MPI_SEND(buf, count, MPI_REAL8, dest, tag, comm, ierr)
end subroutine BMPI_SEND_REAL8_VEC
subroutine BMPI_SEND_REAL8_SCALAR(buf, count, dest, tag, comm, ierr)
   implicit none
   real(kind=8), intent(inout) :: buf
   integer, intent(in) :: count
   integer, intent(in) ::  dest, tag, comm
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_SEND scalar  count /= 1")
   call MPI_SEND(buf, count, MPI_REAL8, dest, tag, comm, ierr)
end subroutine BMPI_SEND_REAL8_SCALAR

!==================================================================
! RECV

subroutine BMPI_RECV_INTEGER4_VEC(buf, count, source, tag, comm, status, ierr)
   implicit none
   integer(kind=4), intent(inout) :: buf(:)
   integer, intent(in) :: count
   integer, intent(in) ::  source, tag, comm
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_RECV  count > size")
   call MPI_RECV(buf, count, MPI_INTEGER4, source, tag, comm, status, ierr)
end subroutine BMPI_RECV_INTEGER4_VEC
subroutine BMPI_RECV_INTEGER4_SCALAR(buf, count, source, tag, comm, status, ierr)
   implicit none
   integer(kind=4), intent(inout) :: buf
   integer, intent(in) :: count
   integer, intent(in) ::  source, tag, comm
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_RECV scalar  count /= 1")
   call MPI_RECV(buf, count, MPI_INTEGER4, source, tag, comm, status, ierr)
end subroutine BMPI_RECV_INTEGER4_SCALAR
subroutine BMPI_RECV_INTEGER8_VEC(buf, count, source, tag, comm, status, ierr)
   implicit none
   integer(kind=8), intent(inout) :: buf(:)
   integer, intent(in) :: count
   integer, intent(in) ::  source, tag, comm
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_RECV  count > size")
   call MPI_RECV(buf, count, MPI_INTEGER8, source, tag, comm, status, ierr)
end subroutine BMPI_RECV_INTEGER8_VEC
subroutine BMPI_RECV_INTEGER8_SCALAR(buf, count, source, tag, comm, status, ierr)
   implicit none
   integer(kind=8), intent(inout) :: buf
   integer, intent(in) :: count
   integer, intent(in) ::  source, tag, comm
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_RECV scalar  count /= 1")
   call MPI_RECV(buf, count, MPI_INTEGER8, source, tag, comm, status, ierr)
end subroutine BMPI_RECV_INTEGER8_SCALAR
subroutine BMPI_RECV_REAL4_VEC(buf, count, source, tag, comm, status, ierr)
   implicit none
   real(kind=4), intent(inout) :: buf(:)
   integer, intent(in) :: count
   integer, intent(in) ::  source, tag, comm
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_RECV  count > size")
   call MPI_RECV(buf, count, MPI_REAL4, source, tag, comm, status, ierr)
end subroutine BMPI_RECV_REAL4_VEC
subroutine BMPI_RECV_REAL4_SCALAR(buf, count, source, tag, comm, status, ierr)
   implicit none
   real(kind=4), intent(inout) :: buf
   integer, intent(in) :: count
   integer, intent(in) ::  source, tag, comm
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_RECV scalar  count /= 1")
   call MPI_RECV(buf, count, MPI_REAL4, source, tag, comm, status, ierr)
end subroutine BMPI_RECV_REAL4_SCALAR
subroutine BMPI_RECV_REAL8_VEC(buf, count, source, tag, comm, status, ierr)
   implicit none
   real(kind=8), intent(inout) :: buf(:)
   integer, intent(in) :: count
   integer, intent(in) ::  source, tag, comm
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count > size(buf, 1, 8)) call BMPI_ERR("BMPI_RECV  count > size")
   call MPI_RECV(buf, count, MPI_REAL8, source, tag, comm, status, ierr)
end subroutine BMPI_RECV_REAL8_VEC
subroutine BMPI_RECV_REAL8_SCALAR(buf, count, source, tag, comm, status, ierr)
   implicit none
   real(kind=8), intent(inout) :: buf
   integer, intent(in) :: count
   integer, intent(in) ::  source, tag, comm
   integer, intent(out) :: status(MPI_STATUS_SIZE)
   integer, intent(out) :: ierr

   if(count /= 1) call BMPI_ERR("BMPI_RECV scalar  count /= 1")
   call MPI_RECV(buf, count, MPI_REAL8, source, tag, comm, status, ierr)
end subroutine BMPI_RECV_REAL8_SCALAR

!==================================================================
! GatherV

subroutine BMPI_GATHERV_INTEGER4_VEC(sendbuf, sendcount, recvbuf, recvcounts, displs, root, comm, ierr)
   implicit none
   integer(kind=4), intent(in) :: sendbuf(:)
   integer(kind=4), intent(out) :: recvbuf(:)
   integer, intent(in) :: sendcount, recvcounts(:), displs(:)
   integer, intent(in) :: root, comm
   integer, intent(out) :: ierr

   call MPI_GATHERV(sendbuf, sendcount, MPI_INTEGER4, recvbuf, recvcounts, displs, MPI_INTEGER4, root, comm, ierr)
end subroutine BMPI_GATHERV_INTEGER4_VEC
subroutine BMPI_GATHERV_INTEGER8_VEC(sendbuf, sendcount, recvbuf, recvcounts, displs, root, comm, ierr)
   implicit none
   integer(kind=8), intent(in) :: sendbuf(:)
   integer(kind=8), intent(out) :: recvbuf(:)
   integer, intent(in) :: sendcount, recvcounts(:), displs(:)
   integer, intent(in) :: root, comm
   integer, intent(out) :: ierr

   call MPI_GATHERV(sendbuf, sendcount, MPI_INTEGER8, recvbuf, recvcounts, displs, MPI_INTEGER8, root, comm, ierr)
end subroutine BMPI_GATHERV_INTEGER8_VEC
subroutine BMPI_GATHERV_REAL4_VEC(sendbuf, sendcount, recvbuf, recvcounts, displs, root, comm, ierr)
   implicit none
   real(kind=4), intent(in) :: sendbuf(:)
   real(kind=4), intent(out) :: recvbuf(:)
   integer, intent(in) :: sendcount, recvcounts(:), displs(:)
   integer, intent(in) :: root, comm
   integer, intent(out) :: ierr

   call MPI_GATHERV(sendbuf, sendcount, MPI_REAL4, recvbuf, recvcounts, displs, MPI_REAL4, root, comm, ierr)
end subroutine BMPI_GATHERV_REAL4_VEC
subroutine BMPI_GATHERV_REAL8_VEC(sendbuf, sendcount, recvbuf, recvcounts, displs, root, comm, ierr)
   implicit none
   real(kind=8), intent(in) :: sendbuf(:)
   real(kind=8), intent(out) :: recvbuf(:)
   integer, intent(in) :: sendcount, recvcounts(:), displs(:)
   integer, intent(in) :: root, comm
   integer, intent(out) :: ierr

   call MPI_GATHERV(sendbuf, sendcount, MPI_REAL8, recvbuf, recvcounts, displs, MPI_REAL8, root, comm, ierr)
end subroutine BMPI_GATHERV_REAL8_VEC

!==================================================================
! ScatterVV

subroutine BMPI_SCATTERV_INTEGER4_VEC(sendbuf, sendcounts, displs, recvbuf, recvcount, root, comm, ierr)
   implicit none
   integer(kind=4), intent(in) :: sendbuf(:)
   integer(kind=4), intent(out) :: recvbuf(:)
   integer, intent(in) :: sendcounts(:), recvcount, displs(:)
   integer, intent(in) :: root, comm
   integer, intent(out) :: ierr

   call MPI_SCATTERV(sendbuf, sendcounts, displs, MPI_INTEGER4, recvbuf, recvcount,  MPI_INTEGER4, root, comm, ierr)
end subroutine BMPI_SCATTERV_INTEGER4_VEC
subroutine BMPI_SCATTERV_INTEGER8_VEC(sendbuf, sendcounts, displs, recvbuf, recvcount, root, comm, ierr)
   implicit none
   integer(kind=8), intent(in) :: sendbuf(:)
   integer(kind=8), intent(out) :: recvbuf(:)
   integer, intent(in) :: sendcounts(:), recvcount, displs(:)
   integer, intent(in) :: root, comm
   integer, intent(out) :: ierr

   call MPI_SCATTERV(sendbuf, sendcounts, displs, MPI_INTEGER8, recvbuf, recvcount, MPI_INTEGER8, root, comm, ierr)
end subroutine BMPI_SCATTERV_INTEGER8_VEC
subroutine BMPI_SCATTERV_REAL4_VEC(sendbuf, sendcounts, displs, recvbuf, recvcount, root, comm, ierr)
   implicit none
   real(kind=4), intent(in) :: sendbuf(:)
   real(kind=4), intent(out) :: recvbuf(:)
   integer, intent(in) :: sendcounts(:), recvcount, displs(:)
   integer, intent(in) :: root, comm
   integer, intent(out) :: ierr

   call MPI_SCATTERV(sendbuf, sendcounts, displs, MPI_REAL4, recvbuf, recvcount, MPI_REAL4, root, comm, ierr)
end subroutine BMPI_SCATTERV_REAL4_VEC
subroutine BMPI_SCATTERV_REAL8_VEC(sendbuf, sendcounts, displs, recvbuf, recvcount,  root, comm, ierr)
   implicit none
   real(kind=8), intent(in) :: sendbuf(:)
   real(kind=8), intent(out) :: recvbuf(:)
   integer, intent(in) :: sendcounts(:), recvcount, displs(:)
   integer, intent(in) :: root, comm
   integer, intent(out) :: ierr

   call MPI_SCATTERV(sendbuf, sendcounts, displs, MPI_REAL8, recvbuf, recvcount, MPI_REAL8, root, comm, ierr)
end subroutine BMPI_SCATTERV_REAL8_VEC

end module bmpi_mod
