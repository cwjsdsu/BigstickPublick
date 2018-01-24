!---------------------------------------------------------------------
!     *******  The following subroutines and functions are     *******
!     *******  to be used in the sequential environment.       *******
!---------------------------------------------------------------------
subroutine MPI_INIT(ierr)
  return
end subroutine MPI_INIT
 
subroutine MPI_COMM_RANK(icomm,iproc,ierr)
  iproc=0
  return
end subroutine MPI_COMM_RANK
 
subroutine MPI_COMM_SIZE(icomm,nproc,ierr)
  nproc=1
  return
end subroutine MPI_COMM_SIZE
 
subroutine MPI_Barrier(icomm,ierr)
  return
end subroutine MPI_Barrier
 
subroutine MPI_Finalize(ierr)
  return
end subroutine MPI_Finalize

subroutine MPI_GET_PROCESSOR_NAME(mpiprocname, maxlen, ierr)
   implicit none
   character (len=*) :: mpiprocname
   integer :: maxlen, ierr
   return
end subroutine MPI_GET_PROCESSOR_NAME

subroutine MPI_Bcast()
  return
end subroutine MPI_Bcast
 
subroutine MPI_Reduce()
  return
end subroutine MPI_Reduce
 
subroutine MPI_Allreduce()
  return
end subroutine MPI_Allreduce
 
subroutine MPI_ALLGATHER()
  return
end subroutine MPI_ALLGATHER

subroutine MPI_GATHER()
  return
end subroutine MPI_GATHER

subroutine MPI_ALLGATHERV()
  return
end subroutine MPI_ALLGATHERV

subroutine MPI_GATHERV()
  return
end subroutine MPI_GATHERV

subroutine MPI_SCATTERV()
   return
end subroutine MPI_SCATTERV

subroutine MPI_Abort()
  stop
  return
end subroutine MPI_Abort

subroutine MPI_INFO_CREATE()
  return
end subroutine MPI_INFO_CREATE

subroutine MPI_INFO_SET()
  return
end subroutine MPI_INFO_SET

subroutine MPI_FILE_OPEN()
  return
end subroutine MPI_FILE_OPEN

subroutine MPI_FILE_CLOSE()
  return
end subroutine MPI_FILE_CLOSE

subroutine MPI_FILE_SET_VIEW()
  return
end subroutine MPI_FILE_SET_VIEW

subroutine MPI_FILE_READ_AT()
  return
end subroutine MPI_FILE_READ_AT

subroutine MPI_FILE_WRITE_AT()
  return
end subroutine MPI_FILE_WRITE_AT

subroutine MPI_COMM_CREATE()
  return
end subroutine MPI_COMM_CREATE

subroutine MPI_COMM_GROUP()
  return
end subroutine MPI_COMM_GROUP

subroutine MPI_GROUP_EXCL()
  return
end subroutine MPI_GROUP_EXCL

subroutine MPI_GROUP_INCL()
  return
end subroutine MPI_GROUP_INCL

subroutine MPI_GROUP_RANK()
  return
end subroutine MPI_GROUP_RANK

subroutine MPI_GROUP_SIZE()
  return
end subroutine MPI_GROUP_SIZE

subroutine MPI_COMM_SPLIT
  return
end subroutine MPI_COMM_SPLIT

subroutine MPI_Send()
  return
end subroutine MPI_Send

subroutine MPI_Recv()
  return
end subroutine MPI_Recv


double precision function MPI_Wtime()
!  integer(4) :: ierr

  real(8)           :: timenow
  call cpu_time(timenow) 
  MPI_Wtime = timenow
  return
end function MPI_Wtime

subroutine MPI_TYPE_SIZE(type, size, ierr)
  integer(kind=4), intent(in):: type
  integer(kind=4), intent(out):: size
  integer(kind=4), intent(out):: ierr
  size = -1
  ierr = -2
end subroutine MPI_TYPE_SIZE

subroutine MPI_FILE_GET_POSITION(FH, OFFSET, IERR)
implicit none

integer :: FH
integer(KIND=selected_int_kind(8)) :: OFFSET
integer :: ierr
end subroutine MPI_FILE_GET_POSITION

subroutine MPI_FILE_READ(FH, BUF, COUNT, DATATYPE, STATUS, IERR)
integer :: FH, COUNT, DATATYPE
integer (kind=1) :: STATUS
integer :: IERR
end subroutine MPI_FILE_READ

subroutine MPI_FILE_WRITE(FH, BUF, COUNT, DATATYPE, STATUS, IERR)
integer :: FH, COUNT, DATATYPE
integer (kind=1) :: STATUS
integer :: IERR
end subroutine MPI_FILE_WRITE

subroutine MPI_FILE_SEEK(FH, OFFSET, WHENCE, IERR)
integer :: FH, WHENCE, IERR
integer(KIND=selected_int_kind(8)) :: OFFSET
end subroutine MPI_FILE_SEEK

subroutine MPI_FILE_SET_SIZE(FH, SIZE, IERR)
implicit none
integer :: FH, IERR
integer (kind=selected_int_kind(8)) :: SIZE
end subroutine MPI_FILE_SET_SIZE

