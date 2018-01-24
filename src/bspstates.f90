!===========================================================================
!  BIGSTICK configuration-interaction shell-model code
!
!===========================================================================
!
!  File bspstates.f90
!
!  routines to compute limits on W
!  W = excitation Weighting; 
!   (equivalent to Nhw in REDSTICK, t in ANTOINE)
!  
!  Subroutines:
!   sort_spsqn     --> sorts the sp.states by W, par, Jz
!   make_hspsqn    --> creates left/right single-particle spaces and q#s
!        (must then call W information)
!   limit_hspspace 	:  limits length of hspspace and computes # of words
!   writeouthsps   : writes out hsps (optional) 
!   section_hsp		:
!   grouper		: master routine, calls group_hsp
!   group_hsp		: divides up hsp of species it by W, jz, par into groups
!   writeoutgroups
!===========================================================================
! 
!  NOTES: 
!  Main purpose of these routines is to prepare the single-particle space
!
!  GOAL 1: Sort the single particle states by (first) W, then parity, then Jz
!          single particle state quantum numbers found in structure "spsqn"
!
!  GOAL 2: split single-particle space into half spaces, "left" with jz < 0 
!          and "right" with jz >= 0
!    
!          hspsqn(it,istate) = structure with quantum numbers
!          it = species (1=proton, 2 neutron)
!          convention:  -1 = left proton, +1 = right proton etc
!
!  GOAL 2A: set up max num of bits needed/ # of words needed
!          nhsps(it) = # of (required) single-particle states 
!          nword(it) = # of words required to store the haiku of species it
!            (I assume similar sizes for left, right haikus)
!===========================================================================
!===================================================================== 
!  Takes single-particle orbit information and inflates into states 
!  (by generating all the m-states) 
!  Takes q#s in derived type orbqn from module sporbit 
!  and fills in q#s into derived type spsqn in module spstate 
!
!  CALLED BY: main program
!
!===================================================================== 
subroutine extract_sps 
 
  use sporbit 
  use spstate 
  use nodeinfo 
  use bmpi_mod 
  implicit none 
 
  integer(4) :: ierr 
 
  integer  :: i, j, k, m 
  integer  :: it 
 
!------------------SET UP S.P. STATE INFO----------------------------------- 
 
  do it = 1,2 
     nsps(it) = 0 
     do i = 1, numorb(it) 
        nsps(it) = nsps(it) + orbqn(it,i)%j + 1 
     end do 
  end do 
  nspsmax = MAX(nsps(1),nsps(2)) 
!  if ( isoflag ) nspsmax = nsps(1) 
!  if ( isoflag ) nsps(2) = nsps(1) 
  allocate ( spsqn(2,nspsmax) ) 
  do it = 1, 2 
     k = 0 
     do i = 1, numorb(it) 
        j = orbqn(it,i)%j 
        do m = -j, j, 2 
           k = k + 1 
           spsqn(it,k)%nr = orbqn(it,i)%nr 
           spsqn(it,k)%l = orbqn(it,i)%l 
           spsqn(it,k)%par = orbqn(it,i)%par 
           spsqn(it,k)%j = orbqn(it,i)%j 
           spsqn(it,k)%w = orbqn(it,i)%w 
           spsqn(it,k)%m = m 
           spsqn(it,k)%orb = i 
           spsqn(it,k)%group = 0   ! KSM: not initialized, set to 0 
        end do ! m 
     end do ! i 
  end do ! it 
 
  return 
   
end subroutine extract_sps 

!===========================================================================
!  sorts the single-particle quantum numbers by W; then by parity; 
!  then finally by Jz
!
!  CALLED BY
!   define_system (only after particle-hole conjugation)
!   getWlimits
!   reconstruct_W_limits

!===========================================================================
subroutine sort_spsqn(it)
  use spstate
  use sporbit
  implicit none
  integer :: it
!..................variables for sorting ...................................

  type (spst) :: spstmp
  integer     :: is
  integer     :: j, k
  integer     :: w, wtmp, par, partmp, m, mtmp
  integer     :: istart, iend

!--------- KSM - initialize for -Wuninitialized
! set to values that would cause problem if used
  iend = 100000000

!------------------SORT ON W------------------------------------------------
  if ( .not. allsameW) then
     do is = 1, nsps(it)
        k = is
        w = spsqn(it,is)%w
        spstmp = spsqn(it,is)

        do j = is + 1, nsps(it)    
           wtmp = spsqn(it,j)%w
           if ( wtmp < w ) then
              k = j
              w = wtmp
              spstmp = spsqn(it,j)
           end if
        end do  !j
        if ( k /= is ) then  !swap
           spsqn(it,k) = spsqn(it,is)
           spsqn(it,is) = spstmp
        end if
     end do  ! is
  end if

!----------------SORT ON PAR within each W----------------------------------

  if ( .not. allsameparity ) then

     istart = 1
     do while( istart <= nsps(it) )

!-----------search for end -------------------------------------------------

        do is = istart, nsps(it)
           if (spsqn(it,is)%w == spsqn(it,istart)%w ) then
              iend = is
           end if
        end do !is
        do is = istart, iend
           k = is
           par = spsqn(it,is)%par
           spstmp = spsqn(it,is)

           do j = is + 1, iend    
              partmp = spsqn(it,j)%par
              if ( partmp < par ) then
                 k = j
                 par = partmp
                 spstmp = spsqn(it,j)
              end if
           end do  !j
           if ( k/=is ) then  !swap
              spsqn(it,k) = spsqn(it,is)
              spsqn(it,is) = spstmp
           end if
        end do  ! is
        istart = iend + 1
     end do  ! while
  end if
  
!------------------SORT ON abs(JZ/m) within each W, par--------------------------

  istart = 1
  do while ( istart <= nsps(it) )
     
!-----------search for end--------------------------------------------------

     do is = istart, nsps(it)
        if ( spsqn(it,is)%w == spsqn(it,istart)%w .and.    &
             spsqn(it,is)%par == spsqn(it,istart)%par ) then
           iend = is
        end if
     end do !is
     do is = istart, iend
        k = is
        m = spsqn(it,is)%m
        spstmp = spsqn(it,is)

        do j = is + 1, iend    
           mtmp = spsqn(it,j)%m
           if ( abs(mtmp) < abs(m) ) then
              k = j
              m = mtmp
              spstmp = spsqn(it,j)
           end if
        end do  !j
        if ( k/=is ) then  !swap
           spsqn(it,k) = spsqn(it,is)
           spsqn(it,is) = spstmp
        end if
     end do  ! is
     istart = iend+1
  end do  ! while

  return
end subroutine sort_spsqn

!===========================================================================
!  Splits up single-particle states into "left" (jz < 0) and "right" 
!  (jz > 0) and then puts into structures hspsqn
!
!  Here  hspsqn(it,i) labels ith state for it;
!  if it < 0, then "left" states
!  else            "right" states
!
!  CALLED BY:
!      main routine
!      getWlimits
!      reconstruct_W_limits
!===========================================================================
subroutine make_hspsqn
  use sporbit
  use spstate
  implicit none

  integer :: it
  integer :: n, nc
  integer :: i
  integer :: aerr
  
!------------count up # of left, right single-particle states---------------

  if(allsameparity)spsqn(:,:)%par = 1
  nhspsmax = 0
  do it = 1, 2
     n = 0
     do i = 1, nsps(it)
        if ( spsqn(it,i)%m >= 0 ) n = n + 1
     end do ! i
     nhsps0(it) = n
     nhsps0(-it) = nsps(it) - n
     nhspsmax = MAX(n,nhspsmax)
  end do  ! it
!------------------ALLOCATE SPACE-------------------------------------------
  if(.not.allocated(hspsqn)) then
     allocate(hspsqn(-2:2,nhspsmax), stat=aerr)
     if(aerr /= 0) call memerror("make_hspsqn")
  end if

!---------------FILL UP-----------------------------------------------------

  do it = 1, 2
     n  = 0
     nc = 0
     
     do i = 1, nsps(it)
        if ( spsqn(it,i)%m >= 0 ) then
           n = n + 1
           hspsqn(it,n) = spsqn(it,i)
        else
           nc = nc + 1
           hspsqn(-it,nc) = spsqn(it,i)
        end if

     end do ! i
  end do  ! it

!----------------- FIND TIME REVERSE OF HSPSQN  added Aug 2010 by CWJ
!                  needed in order to exploit TR in storing matrix elements

   do it = 1,2
      hspsqn(it,:)%tr = -1
      hspsqn(-it,:)%tr = -1
      do i = 1,nhsps0(it)
          if( hspsqn(it,i)%m == 0)then
            hspsqn(it,i)%tr = i
          else
            do n = 1,nhsps0(-it)
              if( hspsqn(it,i)%nr == hspsqn(-it,n)%nr .and. & 
                  hspsqn(it,i)%j == hspsqn(-it,n)%j .and. & 
                  hspsqn(it,i)%l == hspsqn(-it,n)%l .and. & 
                  hspsqn(it,i)%m == -hspsqn(-it,n)%m)then
                     if(hspsqn(it,i)%tr /= -1 .or. hspsqn(-it,n)%tr /= -1)then
                          print*,' problem w time reverse states'
                          print*,it,i,n
                          stop
                     endif
                     hspsqn(it,i)%tr = n
                     hspsqn(-it,n)%tr = i
                     exit
               end if
            end do
          end if

      end do

   end do
  return
end subroutine make_hspsqn

!===========================================================================    

  subroutine writeouthsps(ith)

  use verbosity
  use spstate
  use haiku_info

  implicit none

  integer ith

  integer i
  if(.not.print_hsps)return

  write(hsps_file,*)' SPS Q#s ',ith
  do i = 1,nhsps(ith)
    write(hsps_file,1033)i,hspsqn(ith,i)%j,hspsqn(ith,i)%m, & 
                           hspsqn(ith,i)%par,  hspsqn(ith,i)%w
1033 format(i3,' : ',4i4)
  enddo  !i

  return

  end subroutine writeouthsps

!===========================================================================    
! Finds limits on hspstates based upon W truncation
!     -> this is hspsmax(it)
!     -> then finds nwords
!   must call the W-limit routines first 
!
!  CALLED BY:
!   master_w_limits
!===========================================================================
subroutine limit_hspspace(it)

  use system_parameters
  use sporbit
  use spstate
  use W_info
  use bitstuff
  use haiku_info

  implicit none
  integer :: it
  integer :: i
  if ( allsameW .or. phconj(it)) then      ! modified 7.6.0 to account for p-h conjugatoin
     nhsps(it) = nhsps0(it)
     nhsps(-it) = nhsps0(-it)
  else
     nhsps(it) = 0
     do i = 1, nhsps0(it)
        if ( hspsqn(it,i)%w <= wceiling(it) ) nhsps(it) = nhsps(it) + 1
     end do ! i
     nhsps(-it) = 0
     do i = 1, nhsps0(-it)
        if ( hspsqn(-it,i)%w<=wceiling(it) ) nhsps(-it) = nhsps(-it) + 1
     end do ! i
  end if
  Nword(it) = nhsps(it) / max_bit_word + 1
!  print*,it,nhsps(-it),nhsps(it),nhsps0(-it),nhsps0(it)
!  print*,wceiling(it)
!....... FIND MAX ORBIT LABEL USED
  if(it==1) maxorblabel=0
  do i = 1,nhsps(it)
     maxorblabel = MAX(maxorblabel,hspsqn(it,i)%orb)
  end do
  do i = 1,nhsps(-it)
     maxorblabel = MAX(maxorblabel,hspsqn(-it,i)%orb)
  end do
  return
end subroutine limit_hspspace

!===========================================================================
!  Calls routines to divide hsp spaces into "groups", by W,jz,par
!===========================================================================
!
! CALLED BY:
!  basismaster
! SUBROUTINES CALLED
!	group_hsp
!
subroutine grouper

  implicit none
  
  call group_hsp(1)
  call group_hsp(-1)
  call group_hsp(2)
  call group_hsp(-2)

  return
end subroutine grouper
      
!===========================================================================
!  Divides up hsp of species it by W, jz, par into groups
!
! called by: grouper
!
!===========================================================================
subroutine group_hsp(ith)

  use spstate
  use W_info
  use haiku_info

  implicit none
  
  integer :: ith         ! species

  integer :: i
  integer :: n
  
!------------------COUNT UP # OF SECTIONS-----------------------------------
  ngroups(ith) = 1
  do i = 2, nhsps(ith)
     if ( hspsqn(ith,i)%w /= hspsqn(ith,i-1)%w .or. hspsqn(ith,i)%m /= hspsqn(ith,i-1)%m & 
         .or. hspsqn(ith,i)%par /= hspsqn(ith,i-1)%par ) then
        ngroups(ith) = ngroups(ith) + 1
     end if
  end do ! it
!      print*,' there are ',nsections(it),' sections ',nhsps(it)
!------------------ALLOCATE-------------------------------------------------

  allocate ( group(ith)%jz(ngroups(ith)),  group(ith)%par(ngroups(ith)), group(ith)%w(ngroups(ith)) )
  allocate ( group(ith)%start(ngroups(ith)), group(ith)%fin(ngroups(ith)) )

!------------------Find START, FIN of each section--------------------------
  n = 1
  group(ith)%w(1) = hspsqn(ith,1)%w
  group(ith)%jz(1) = hspsqn(ith,1)%m
  group(ith)%par(1) = hspsqn(ith,1)%par
  hspsqn(ith,1)%group = n
  group(ith)%start(1) = 1
  do i = 2, nhsps(ith)
     if ( hspsqn(ith,i)%w > hspsqn(ith,i-1)%w .or.  hspsqn(ith,i)%m /= hspsqn(ith,i-1)%m & 
         .or. hspsqn(ith,i)%par /= hspsqn(ith,i-1)%par ) then
        n = n + 1
        group(ith)%w(n) = hspsqn(ith,i)%w
        group(ith)%par(n) = hspsqn(ith,i)%par
        group(ith)%jz(n) = hspsqn(ith,i)%m

        group(ith)%start(n) = i
        group(ith)%fin(n-1) = i - 1
     end if
     hspsqn(ith,i)%group = n

  end do
  group(ith)%fin(n) = nhsps(ith)

  return
end subroutine group_hsp
!===========================================================================
  subroutine writeoutgroups(ith)
  use spstate
  use haiku_info
  use verbosity
  implicit none
  integer ith
  integer g

  if (.not.print_groups)return
   write(group_file,*)' '
   write(group_file,*)' GROUPS ',ith

   do g = 1,ngroups(ith)
     write(group_file,2444)g,group(ith)%jz(g),group(ith)%par(g),group(ith)%w(g)
2444 format(i3,' : ',3i3)

   enddo ! g
  return
  end subroutine writeoutgroups


!===========================================================================


!End of subroutines in this library.
!===========================================================================
