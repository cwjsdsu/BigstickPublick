!---------------------------------------------------------------------------
!  information about sections
!  sections are partitions of the single-particle-state space
!  partitioned by W
!---------------------------------------------------------------------------
module spsections
  implicit none

  integer :: nsections(-2:2)

  type sectionbase
     integer :: w
     integer :: start,fin
     integer :: nsps   
  end type sectionbase

  type sectionit
     type (sectionbase),pointer :: list(:)
  end type sectionit

  type (sectionit) :: section(-2:2)

contains
	

!===========================================================================
!  Calls routines to divide hsp spaces into "sections", by W
!
! CALLED BY: basismaster
!
! CALLS: section_hsp
!===========================================================================
!
! SUBROUTINES CALLED
!	section_hsp
!
subroutine master_sections

  implicit none
  
  call section_hsp(1)
  call section_hsp(-1)
  call section_hsp(2)
  call section_hsp(-2)

  return
end subroutine master_sections
      
!===========================================================================
!  Divides up hsp of species it by W into sections
!
!  CALLED BY: master_sections
!
!===========================================================================
subroutine section_hsp(it)

  use spstate
  use W_info
  use haiku_info
  implicit none
  
  integer :: it         ! species

  integer :: i
  integer :: n
  
!------------------COUNT UP # OF SECTIONS-----------------------------------
  nsections(it) = 1
  do i = 2, nhsps(it)
     if ( hspsqn(it,i)%w > hspsqn(it,i-1)%w ) then
        nsections(it) = nsections(it) + 1
     end if
  end do ! it
!      print*,' there are ',nsections(it),' sections ',nhsps(it)
!------------------ALLOCATE-------------------------------------------------

  allocate ( section(it)%list(nsections(it)) )

!------------------Find START, FIN of each section--------------------------
  n = 1
  section(it)%list(1)%w = hspsqn(it,1)%w
  section(it)%list(1)%start = 1
  do i = 2, nhsps(it)
     if ( hspsqn(it,i)%w > hspsqn(it,i-1)%w ) then
        n = n + 1
        section(it)%list(n)%w = hspsqn(it,i)%w
        section(it)%list(n)%start = i
        section(it)%list(n-1)%fin = i - 1
     end if
  end do
  section(it)%list(n)%fin = nhsps(it)

!------------------DETERMINE HOW MANY IN EACH SECTION-----------------------
  do n = 1, nsections(it)
     section(it)%list(n)%nsps = 1 + section(it)%list(n)%fin - &
          section(it)%list(n)%start
  end do

  return
end subroutine section_hsp

end module spsections
  
