!===================================================================
!
! Simple utilities for math/porting

! I decided to use overloading so that the polymorphic calls are
! automatically bound to the right wrapper.
!
! Initial version:  Ken McElvain   17Jan2015
!
module butil_mod
   implicit none

   interface BMIN
      MODULE PROCEDURE BMIN_I4I4, BMIN_I4I8, BMIN_I8I4, BMIN_I8I8, &
      BMIN_R4R4, BMIN_R4R8, BMIN_R8R4, BMIN_R8R8
   end interface BMIN

   interface BMAX
      MODULE PROCEDURE BMAX_I4I4, BMAX_I4I8, BMAX_I8I4, BMAX_I8I8, &
      BMAX_R4R4, BMAX_R4R8, BMAX_R8R4, BMAX_R8R8
   end interface BMAX

contains

! min for kind=4 integers
integer(kind=4) function BMIN_I4I4(a,b)
   implicit none
   integer(kind=4),intent(in) :: a, b
   if(a < b) then
      BMIN_I4I4 = a
   else
      BMIN_I4I4 = b
   end if
   RETURN
end function BMIN_I4I4

! min for mixed case
integer(kind=8) function BMIN_I4I8(a,b)
   implicit none
   integer(kind=4),intent(in) :: a
   integer(kind=8),intent(in) :: b
   integer(kind=8) :: a8
   a8 = a
   if(a8 < b) then
      BMIN_I4I8 = a8
   else
      BMIN_I4I8 = b
   end if
   RETURN
end function BMIN_I4I8

! min for other mixed case
integer(kind=8) function BMIN_I8I4(a,b)
   implicit none
   integer(kind=8),intent(in) :: a
   integer(kind=4),intent(in) :: b
   integer(kind=8) :: b8
   b8 = b
   if(a < b8) then
      BMIN_I8I4 = a
   else
      BMIN_I8I4 = b8
   end if
   RETURN
end function BMIN_I8I4

! min for kind=8 integers
integer(kind=8) function BMIN_I8I8(a,b)
   implicit none
   integer(kind=8),intent(in) :: a
   integer(kind=8),intent(in) :: b
   if(a < b) then
      BMIN_I8I8 = a
   else
      BMIN_I8I8 = b
   end if
   RETURN
end function BMIN_I8I8

! min for kind=4 reals
real(kind=4) function BMIN_R4R4(a,b)
   implicit none
   real(kind=4),intent(in) :: a, b
   if(a < b) then
      BMIN_R4R4 = a
   else
      BMIN_R4R4 = b
   end if
   RETURN
end function BMIN_R4R4

! min for mixed kind=4,8 reals
real(kind=8) function BMIN_R4R8(a,b)
   implicit none
   real(kind=4),intent(in) :: a
   real(kind=8),intent(in) :: b
   real(kind=8) :: a8
   a8 = a
   if(a8 < b) then
      BMIN_R4R8 = a8
   else
      BMIN_R4R8 = b
   end if
   RETURN
end function BMIN_R4R8

! min for mixed kind=8,4 reals
real(kind=8) function BMIN_R8R4(a,b)
   implicit none
   real(kind=8),intent(in) :: a
   real(kind=4),intent(in) :: b
   real(kind=8) :: b8
   b8 = b
   if(a < b8) then
      BMIN_R8R4 = a
   else
      BMIN_R8R4 = b8
   end if
   RETURN
end function BMIN_R8R4

! min for kind=8 reals
real(kind=8) function BMIN_R8R8(a,b)
   implicit none
   real(kind=8),intent(in) :: a
   real(kind=8),intent(in) :: b
   if(a < b) then
      BMIN_R8R8 = a
   else
      BMIN_R8R8 = b
   end if
   RETURN
end function BMIN_R8R8

! max for kind=4 integers
integer(kind=4) function BMAX_I4I4(a,b)
   implicit none
   integer(kind=4),intent(in) :: a, b
   if(a > b) then
      BMAX_I4I4 = a
   else
      BMAX_I4I4 = b
   end if
   RETURN
end function BMAX_I4I4

! max for mixed case
integer(kind=8) function BMAX_I4I8(a,b)
   implicit none
   integer(kind=4),intent(in) :: a
   integer(kind=8),intent(in) :: b
   integer(kind=8) :: a8
   a8 = a
   if(a8 > b) then
      BMAX_I4I8 = a8
   else
      BMAX_I4I8 = b
   end if
   RETURN
end function BMAX_I4I8

! max for other mixed case
integer(kind=8) function BMAX_I8I4(a,b)
   implicit none
   integer(kind=8),intent(in) :: a
   integer(kind=4),intent(in) :: b
   integer(kind=8) :: b8
   b8 = b
   if(a > b8) then
      BMAX_I8I4 = a
   else
      BMAX_I8I4 = b8
   end if
   RETURN
end function BMAX_I8I4

! max for kind=8 integers
integer(kind=8) function BMAX_I8I8(a,b)
   implicit none
   integer(kind=8),intent(in) :: a
   integer(kind=8),intent(in) :: b
   if(a > b) then
      BMAX_I8I8 = a
   else
      BMAX_I8I8 = b
   end if
   RETURN
end function BMAX_I8I8

! max for kind=4 reals
real(kind=4) function BMAX_R4R4(a,b)
   implicit none
   real(kind=4),intent(in) :: a, b
   if(a > b) then
      BMAX_R4R4 = a
   else
      BMAX_R4R4 = b
   end if
   RETURN
end function BMAX_R4R4

! max for mixed kind=4,8 reals
real(kind=8) function BMAX_R4R8(a,b)
   implicit none
   real(kind=4),intent(in) :: a
   real(kind=8),intent(in) :: b
   real(kind=8) :: a8
   a8 = a
   if(a8 > b) then
      BMAX_R4R8 = a8
   else
      BMAX_R4R8 = b
   end if
   RETURN
end function BMAX_R4R8

! max for mixed kind=8,4 reals
real(kind=8) function BMAX_R8R4(a,b)
   implicit none
   real(kind=8),intent(in) :: a
   real(kind=4),intent(in) :: b
   real(kind=8) :: b8
   b8 = b
   if(a > b8) then
      BMAX_R8R4 = a
   else
      BMAX_R8R4 = b8
   end if
   RETURN
end function BMAX_R8R4

! max for kind=8 reals
real(kind=8) function BMAX_R8R8(a,b)
   implicit none
   real(kind=8),intent(in) :: a
   real(kind=8),intent(in) :: b
   if(a > b) then
      BMAX_R8R8 = a
   else
      BMAX_R8R8 = b
   end if
   RETURN
end function BMAX_R8R8

subroutine overflowerr(msg)
   implicit none
   character(len=*), intent(in) :: msg
   print *, msg
   call printstacktrace  ! if enabled in build
   stop 1
end subroutine overflowerr

! utility to convert an 8 byte integer to a 4 byte integer
! with overflow detection
integer(kind=4) function bint8to4(i8, id)
   implicit none
   integer(kind=8), intent(in) :: i8
   integer :: id ! identification of caller
   integer(kind=4) :: i4
   i4 = int(i8, kind(i4))
   if(i4 /= i8) then
      print *, "Stuffing oversized number in integer(kind=4)", i8, ", at id=", id
      call printstacktrace  ! if enabled in build
      stop 1
   end if
   bint8to4 = i4
   return
end function bint8to4

! generic way to stop on error
! should we print iproc?
subroutine errstop(msg)
   implicit none
   character(len=*), intent(in) :: msg
   write(6,"(A)") msg
   flush(6)
   call printstacktrace  ! if enabled in build
   stop 1
end subroutine errstop

!==========================================================
      real function zeta(i,j)

      implicit none

      integer i,j

      zeta = 1.0
      if(i ==j)zeta = sqrt(2.)
      return
      end function zeta

end module butil_mod
