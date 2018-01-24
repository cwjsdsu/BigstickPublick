!==============================================================
!
!  file B3BME.f90
!
!  routines to set up to handle and store 3-body matrix elements
!
!  initiated Sept 2009 CWJ @ SDSU
!
!================================================================
!
! subroutine threebodyquestion :  asks if 3-body forced used
!
! subroutine threebodysetup : sets up arrays for decoupled 3-body matrix elements
! CALLS:	 XXXsetup
!           master_cross_ntuple (in bjumplib4.f90)
!           XXYsetup
!
! subroutine XXXsetup -- PPP/NNN equivalent of prepare_to_uncoupleXXtbme
!  CALLS:  maxwtrip
!          count_create_triplets
!	     sort3s
!          count_XXX_3bmes
!          find_triplet_start
!
! subroutine count_create_triplets
!
! subroutine count_XXX_3bmes
!
! subroutine find_triplet_start
!
! subroutine maxWtrip  - finds max W for triplets
!
! subroutine sort3s   - heapsort of triplets
!
!               *** the next routines prepare to readin matrix elements * * *
!
! subroutine threebodymaster:  orchestrates setup
!
! called by main program (in bigstick_main.f90)
!
! CALLS:    fetch3bodyinput  (in b3bme_input.f90)
!           XXXmaster
!           XXYmaster   (in b3bmelib2.f90)
!
! subroutine XXXmaster:  allocates arrays for PPP/NNN matrix elements and zeroes out
!
! subroutine XXYmaster:  allocates arrays for PPN/PNN matrix elements and zeroes out

! module wrapper creates interfaces automatically
module b3bme_mod
contains

!
!================================================================
!
!  CALLED BY HOOPMASTER in bbasislib1.f90
!
subroutine threebodyquestion
  use menu_choices
  use flags3body
  use nodeinfo
  use io
  use bmpi_mod
  implicit none

  integer(4)   :: ierr
  character(1) :: ychar
  
  threebody = .false.
  
  if(.not.threebodycheck) return
  if( menu_char == 'g' .or. menu_char=='p')then
       threebody= .false.
       return
  end if
  if ( iproc == 0 ) then
     write(6,*)' Do you want 3-body forces (y/n) ?'

     if(auto_input)then
        read(autoinputfile,'(a)')ychar
     else
        read(5,'(a)')ychar
        write(autoinputfile,'(a)')ychar
     endif
  end if
 call BMPI_BCAST(ychar,1,0,icomm,ierr)
 call BMPI_BARRIER(icomm,ierr)

  if(ychar == 'y' .or. ychar == 'Y')then 
          threebody = .true.
  else
          threebody = .false.
  endif

  return

end subroutine threebodyquestion
!================================================================
!
! subroutine threebodysetup
!
! sets up arrays for decoupled 3-body matrix elements
!
! CALLED BY -hoopmaster (in bbasislib1.f90)
!
! CALLS:
!	XXXsetup
!     master_cross_ntuple
!     XXYsetup
!
  subroutine threebodysetup_mod
  use system_parameters
  use flags3body
  use ntuple_info
  use menu_choices
  use jump_mod
  implicit none
  integer it
!...............................................................
  if(.not.threebody)return
  do it = 1,2
    if(np(it) <3)cycle

    if(restrict_ntuples) call master_q_ntuple(3, it)
    if(menu_char /= 'm')then
        call XXXsetup(it)
    end if
  enddo

  if( np(1) > 1 .and. np(2) > 0) then
      call master_cross_ntuple(2,1)

      if(menu_char /= 'm')then
         call XXYsetup(1)
      endif
  end if

  if( np(1) > 0 .and. np(2) > 1)then
      call master_cross_ntuple(1,2)
      if(menu_char /= 'm')then
 
         call XXYsetup(2)
      end if
  endif

  return

  end subroutine threebodysetup_mod

!=================================================================
!
! subroutine XXXsetup
!
!  equivalent to prepare_to_uncoupleXXtbme
!
!  creates triplets, orders, and sets up for matrix elements
!
!  CALLS:  maxwtrip
!          count_create_triplets
!	     sort3s
!          count_XXX_3bmes
!          find_triplet_start
!
!  CALLED BY threebodysetup (above)
!
  subroutine XXXsetup(it)

  use tripdef
  use interactions3body
  use nodeinfo
  use btbme_mod

  implicit none
! include 'binterfaces.inc'

  integer it
  integer maxw3,maxw2
  integer n3s
  type (triplet_qn), pointer :: tripXXX(:) 
  integer n,n2,i

  call maxwtrip(it,maxw3)
  call maxWpair(it,maxw2)

  tripXXX => XXX3(it)%trip
  call count_create_triplets(it,maxw2,maxw3,.false.,n3s,tripXXX)
  if(iproc==0)print*,'For species ',it,' there are ',n3s,' triplets '
  nXXX(it) = n3s
  tripXXX => XXX3(it)%trip

  call count_create_triplets(it,maxw2,maxw3,.true.,n3s,tripXXX)
  call sort3s(n3s,tripXXX)

  call count_XXX_3bmes(it,n3s,tripXXX,n)
  call find_triplet_start(it,n3s,tripXXX,n2)
  if(n /= n2)then
     if(iproc==0)print*,' problem with # of 3BMES ',n,n2
     stop
  endif
!  deallocate(tripxxx)

  if(it == 1)then
     nmatppp = n
     if(iproc==0)print*,n,' PPP matrix elements '
  endif
  if(it == 2)then
     nmatNNN = n
     if(iproc==0)print*,n,' NNN matrix elements '

  endif

  return
  end subroutine XXXsetup

!=================================================================
!
! creates triplets of protons or neutrons
!
! INPUT:
!  it = species, 1 = proton, 2 = neutrons
!  maxw3 = max w that can be carried by this triplet
!  create = flag; if = .false. then just count, else fills arrays tripXXX
!
! OUTPUT:
!  n3s = # of triplets (counted up with create = .false.)
!  tripXXX(n) = information on triplets (when create = .true.)
!
! CALLED BY XXXsetup
!
   subroutine count_create_triplets(it,maxw2,maxw3,create,n3s,tripXXX)

   use spstate
   use haiku_info
   use system_parameters
   use interactions3body

   use tripdef
   use ntuple_info
   implicit none
   integer it
   integer maxw2,maxw3
   logical create
   integer n3s
   type (triplet_qn),pointer :: tripXXX(:) 

   integer isps,jsps,ksps
   integer i,j,k
   integer ith,jth,kth
   integer isgn,jsgn,ksgn

   integer indx
   integer i3s
   integer m3,w3,par3

   integer :: aerr

   if(create)then
      tripXXX(:)%m = -99  ! flag for error
      tripXXX(:)%indx = 0
   endif
   indx = 0 
  
   i3s = 0

   do isps = 3, nhsps(it)+nhsps(-it)
      if(isps > nhsps(-it))then
         ith = it
         i = isps - nhsps(-it)
         isgn = 1
      else
         ith = -it
         i = isps
         isgn = -1
      endif
      do jsps = 2, isps-1
         if(jsps > nhsps(-it))then
            jth = it
            j = jsps - nhsps(-it)
            jsgn = 1
         else
            jth = -it
            j = jsps
            jsgn = -1
         endif
         do ksps = 1,jsps-1
            if(ksps > nhsps(-it))then
               kth = it
               k = ksps - nhsps(-it)
               ksgn = 1
            else
               kth = -it
               k = ksps
               ksgn = -1
            endif
            if (hspsqn(ith,i)%w+hspsqn(jth,j)%w+hspsqn(kth,k)%w > maxw3)cycle
!................ also cut 2 at a time
            if (hspsqn(ith,i)%w+hspsqn(jth,j)%w > maxw2)cycle
            if (hspsqn(ith,i)%w+hspsqn(kth,k)%w > maxw2)cycle
            if (hspsqn(kth,k)%w+hspsqn(jth,j)%w > maxw2)cycle

!....................... RESTRICTION.......................................
            m3 =  hspsqn(ith,i)%m + hspsqn(jth,j)%m+hspsqn(kth,k)%m
            w3 = hspsqn(ith,i)%w+hspsqn(jth,j)%w+hspsqn(kth,k)%w
            par3 = hspsqn(ith,i)%par*hspsqn(jth,j)%par*hspsqn(kth,k)%par
            par3 = (3-  par3)/2

            if(restrict_ntuples)then
               if(w3 > tripletlim(it)%maxw .or. w3 < tripletlim(it)%minw)cycle
               if( par3 > tripletlim(it)%w(w3)%parmax) cycle
               if( par3 < tripletlim(it)%w(w3)%parmin) cycle
               if( m3 > tripletlim(it)%w(w3)%par(par3)%jzmax) cycle
               if( m3 < tripletlim(it)%w(w3)%par(par3)%jzmin) cycle
            end if
            i3s = i3s+1
            if(create)then
               if(tripXXX(i3s)%m /= -99)then
                 print*,' oops already occupied '
                 print*,isps,jsps,ksps,i3s  !,pairXX(ipair)%m
                 stop
               endif
               
               tripXXX(i3s)%m = m3 
               tripXXX(i3s)%w = w3
               tripXXX(i3s)%par = par3 
               tripXXX(i3s)%indx = (isps-1)*(isps-2)*(isps-3)/6 +  & 
                                   (jsps-1)*(jsps-2)/2 +ksps
               tripXXX(i3s)%ia = i*isgn
               tripXXX(i3s)%ib = j*jsgn
               tripXXX(i3s)%ic = k*ksgn

            endif
         enddo !  ksps
      enddo ! jsps
   enddo  ! isps
   if(.not.create)then
      n3s = i3s
      allocate(XXX3(it)%trip(n3s), stat=aerr )
      if(aerr /= 0) call memerror("count_create_triplets")
   endif
   return
   end subroutine count_create_triplets

!====================================================================
!
!  count up how many 3-body matrix elements of type XXX
!  NOTE: need to restrict on change on W later
!
! CALLED BY XXXsetup
!
   subroutine count_XXX_3bmes(it,n3s,triplet,nxxxmat)
   use tripdef
   use interactions3body
   use nodeinfo
   implicit none

   integer it
   integer n3s
   type (triplet_qn) :: triplet(n3s)
   integer nxxxmat

   integer m, mstart,mend
   integer par,parstart,parend
   integer istart,iend,ifinish
   integer i
   integer :: aerr

   nxxxmat = 0

   mstart = triplet(1)%m
   mend   = triplet(n3s)%m
   XXX3(it)%jzstart = mstart
   XXX3(it)%jzend = mend
   allocate( XXX3(it)%meref(mstart:mend,2), stat=aerr)
   if(aerr /= 0) call memerror("count_XXX_3bmes 1")
   allocate( XXX3(it)%mestart(mstart:mend,2), stat=aerr)
   if(aerr /= 0) call memerror("count_XXX_3bmes 2")

   istart  = 1

   do m = mstart,mend, 2   
!......given m, search for start,finish of parity
      ifinish = istart
      do i = istart,n3s
         ifinish = i
         if(triplet(i)%m /= m)then
            ifinish = i-1
            exit
         endif
      enddo
      parstart = triplet(istart)%par
      parend   = triplet(ifinish)%par

      do par = parstart,parend
!............ NOW FIND ACTUAL START, FINISH FOR THIS M, PARITY
         if(istart > ifinish)cycle
         iend = istart
         do i = istart,ifinish
            iend = i
            if(triplet(i)%par /= par)then
               iend = i-1
               exit
            endif
         enddo
!.......... COMPUTE # OF 3BMES
         nxxxmat = nxxxmat + (iend-istart+1)*(iend-istart+2)/2
!........... SET UP FOR NEXT ROUND
         istart = iend + 1
      end do  ! par

   end do  ! m
   if(it == 1)then
   if(iproc==0)print*, 'There are ',nxxxmat, ' PPP matrix elements '
   else
   if(iproc==0)print*, 'There are ',nxxxmat, ' NNN matrix elements '

   endif

   return
   end subroutine count_XXX_3bmes

!====================================================================

!
!  subroutine find_triplet_start
!
!  finds, in a list of triplets, where quantum numbers start, end
!
! CALLED BY: XXXsetup
!
   subroutine find_triplet_start(it,n3s,triplet,nxxxmat)
   use tripdef
   use interactions3body
   use haiku_info
   implicit none

   integer it
   integer n3s
   type (triplet_qn) :: triplet(n3s)
   integer nxxxmat

   integer m, mstart,mend
   integer par,parstart,parend
   integer istart,iend,ifinish
   integer i,j,iref,n
   integer, pointer :: map(:)
   integer :: aerr

   mstart = triplet(1)%m
   mend   = triplet(n3s)%m

   allocate( XXX3(it)%meref(mstart:mend,2), stat=aerr)
   if(aerr /= 0) call memerror("find_triplet_start 1")
   allocate( XXX3(it)%mestart(mstart:mend,2), stat=aerr)
   if(aerr /= 0) call memerror("find_triplet_start 2")

   nxxxmat = 1

   XXX3(it)%meref(mstart,:) = 0
   XXX3(it)%mestart(mstart,:) = 0
   iref = 1
   do i = 2,n3s
        m = triplet(i)%m
        par = triplet(i)%par
        if(m > triplet(i-1)%m .or. par /= triplet(i-1)%par )then
           iref = i
           XXX3(it)%meref(m,par) = iref -1
           XXX3(it)%mestart(m,par) = nxxxmat
           
        endif 
        do j = iref,i
            nxxxmat = nxxxmat+1
        enddo
    enddo

!-------- NOW FIND MAP

   if(it == 1)then
      n = nhsps(it) + nhsps(-it)
      allocate( mapPPP( n*(n-1)*(n-2)/6), stat=aerr)
      if(aerr /= 0) call memerror("find_triplet_start 3")
      map => mapPPP    
   else
      n = nhsps(it) + nhsps(-it)
      allocate( mapNNN( n*(n-1)*(n-2)/6), stat=aerr)
      if(aerr /= 0) call memerror("find_triplet_start 4")
      map => mapNNN 
   endif
   map(:) = 0
   do i = 1,n3s
     
     j = triplet(i)%indx
     if ( map(j) /=0 ) then
           print*,' problem, triplet already found '
           print*,i,j

        stop
     end if
     map(j) = i
   end do

   return
   end subroutine find_triplet_start


!==========================================================
!
!  CALLED BY: XXXsetup
!             XXYsetup
!
     subroutine sort3s(n3s,triplet)

!
!  sort triplets  first by
!    total M
!    parity
!    W
!
      use tripdef
      implicit none
      integer n3s
      type (triplet_qn) :: triplet(n3s)
      type (triplet_qn) :: triptmp
      integer i,j,k,l,n

      integer m,par,w

      integer istart,istop
      integer jstart,jstop


!--------- REVISED: HEAPSORT ON Jz------
      l = n3s/2+1
      n = n3s
      do while (n > 0)
          if(l > 1)then
             l = l- 1
             triptmp = triplet(l)
          else
             triptmp = triplet(n)
             triplet(n) = triplet(1)
             n = n -1
             if(n == 1)then
                 triplet(1) = triptmp
                 exit
              end if
          end if
          i = l
          j = l+l
          do while ( j <= n)
             if( j < n)then
                if(triplet(j)%m< triplet(j+1)%m)j = j+1
             endif
             if( triptmp%m < triplet(j)%m)then 
                 triplet(i) = triplet(j)
                 i = j
                 j = j + j
             else
                 j = n + 1
             end if
          end do      
          triplet(i) = triptmp
      end do
!------------- next sort on parity ----------

!---------- find start, stop for each M value

      m = triplet(1)%m
      istart = 1
      istop  = 1
      do while(m <= triplet(n3s)%m )
!................. ERROR TRAP...................
        if(istart > n3s)then
          print*,' problem in sorting triplets '
          print*,istart,n3s,istop
          print*,m,triplet(n3s)%m
          stop
        endif
!.............END TRAP ...........
        if(istart == n3s)then
          istop = istart
        else
          do while(triplet(istop+1)%m == m )
           istop = istop+1
           if(istop == n3s) goto 3
          enddo
3         continue
        endif
!--------- PERFORM ACTUAL SORT --------

        l = (istop-istart+1)/2+istart
        n = istop
        do while (n >= istart)
          if(l > istart)then
             l = l- 1
             triptmp = triplet(l)
          else
             triptmp = triplet(n)
             triplet(n) = triplet(istart)
             n = n -1
             if(n == istart)then
                 triplet(istart) = triptmp
                 exit
              end if
          end if
          i = l
          j = l+l-(istart-1)
          do while ( j <= n)
             if( j < n)then
                   if(triplet(j)%par< triplet(j+1)%par)j = j+1
             endif
             if( triptmp%par < triplet(j)%par)then 
                 triplet(i) = triplet(j)
                 i = j
                 j = j + j-(istart-1)
             else
                 j = n + 1
             end if
          end do      
          triplet(i) = triptmp
       end do

!----------- FINALLY, SORT ON W --------
!              NB: not clear that I need this any longer

        if(triplet(istart)%par == triplet(istop)%par)then

           jstart = istart
           jstop  = istop
        
           l = (jstop-jstart+1)/2+jstart
           n = jstop
           do while (n >= jstart)
              if(l > jstart)then
                 l = l- 1
                 triptmp = triplet(l)
              else
                 triptmp = triplet(n)
                 triplet(n) = triplet(jstart)
                 n = n -1
                 if(n == jstart)then
                     triplet(jstart) = triptmp
                     exit
                 end if
               end if
               i = l
               j = l+l-(jstart-1)
               do while ( j <= n)
                   if( j < n)then
                       if(triplet(j)%w< triplet(j+1)%w)j = j+1
                   endif
                   if( triptmp%w < triplet(j)%w)then 
                       triplet(i) = triplet(j)
                       i = j
                       j = j + j-(jstart-1)
                   else
                       j = n + 1
                   end if
               end do      
             triplet(i) = triptmp
          end do

        else

           jstart = istart
           jstop = istart
           do while(triplet(jstop+1)%par == triplet(jstart)%par .and. jstop < n3s)

              jstop = jstop+1
              if(jstop > n3s)then
                 print*,' error in jstop '
                 stop
              endif
           enddo

           l = (jstop-jstart+1)/2+jstart
           n = jstop
           do while (n >= jstart)
              if(l > jstart)then
                 l = l- 1
                 triptmp = triplet(l)
              else
                 triptmp = triplet(n)
                 triplet(n) = triplet(jstart)
                 n = n -1
                 if(n == jstart)then
                     triplet(jstart) = triptmp
                     exit
                 end if
               end if
               i = l
               j = l+l-(jstart-1)
               do while ( j <= n)
                   if( j < n)then
                       if(triplet(j)%w< triplet(j+1)%w)j = j+1
                   endif
                   if( triptmp%w < triplet(j)%w)then 
                       triplet(i) = triplet(j)
                       i = j
                       j = j + j-(jstart-1)
                   else
                       j = n + 1
                   end if
               end do      
             triplet(i) = triptmp
          end do


           jstart = jstop+1
           jstop = istop
           l = (jstop-jstart+1)/2+jstart
           n = jstop
           do while (n >= jstart)
              if(l > jstart)then
                 l = l- 1
                 triptmp = triplet(l)
              else
                 triptmp = triplet(n)
                 triplet(n) = triplet(jstart)
                 n = n -1
                 if(n == jstart)then
                     triplet(jstart) = triptmp
                     exit
                 end if
               end if
               i = l
               j = l+l-(jstart-1)
               do while ( j <= n)
                   if( j < n)then
                       if(triplet(j)%w< triplet(j+1)%w)j = j+1
                   endif
                   if( triptmp%w < triplet(j)%w)then 
                       triplet(i) = triplet(j)
                       i = j
                       j = j + j-(jstart-1)
                   else
                       j = n + 1
                   end if
               end do      
             triplet(i) = triptmp
          end do

        endif
!-----------------------------------
        m = m+2
        istart = istop+1
        istop = istart
      enddo  

      return
      end subroutine sort3s


!==========================================================
!  subroutine maxwtrip
!  finds the maximum W allowed for a triplet
!  by comparing max W for N, N-3 particles
!
!  CALLED BY: XXXsetup
!  
   subroutine maxwtrip(it,maxw3)

   use spstate
   use W_info
   use haiku_info
   use system_parameters

   implicit none
   integer it
   integer maxw3
   
   integer w
   integer i
   maxw3 = 0
   if(Np(it) < 3)return
   w = 0
   do i =1,Np(it)-3
        w = w + spsqn(it,i)%w
   enddo
   maxw3 = maxW(it) - w
   return
   end subroutine maxwtrip


!================================================================
!
! SUBROUTINE THREEBODYMASTER
!
! called by main program (in bigstick_main.f90)
!
! CALLS:    fetch3bodyinput
!           fetch3bodyinputcoupled
!           XXXmaster  -- allocates arrays and zeroes out
!           XXYmaster  -- allocates arrays and zeroes out
!
  subroutine threebodymaster
  use system_parameters
  use flags3body
  use ntuple_info
  use nodeinfo
  implicit none


  integer it
!................ SET UP FOR INPUT...............................
!  call fetch3bodyinput
  call fetch3bodyinputcoupled

!...............................................................
  if(.not.threebody)return
  do it = 1,2
    if(np(it) <3)cycle

    call XXXmaster(it)
  enddo

  if( np(1) > 1 .and. np(2) > 0) then

      call XXYmaster(1)

  end if

  if( np(1) > 0 .and. np(2) > 1)then
      call XXYmaster(2)
  endif

  return
  end subroutine threebodymaster 

!=================================================================
!
! subroutine XXXmaster
! 
! allocates arrays for PPP/NNN matrix elements and zeroes out
!
! called by threebodymaster
!
  subroutine XXXmaster(it)

  use tripdef
  use interactions3body

  implicit none


  integer :: it
  integer(8) :: i
  integer :: aerr


  if(it == 1)then
     if(nmatppp > 0) then
         allocate(hmatPPP(nmatppp), stat=aerr)
         if(aerr /= 0) call memerror("XXXmaster 1")
     end if
     do i = 1,nmatppp
         hmatPPP(i) = 0.0
     end do
  endif
  if(it == 2)then
     if(nmatNNN > 0) then
        allocate(hmatNNN(nmatNNN), stat=aerr)
        if(aerr /= 0) call memerror("XXXmaster 1")
     end if
     do i = 1,nmatnnn
         hmatNNN(i) = 0.0
     end do
  endif

  return
  end subroutine XXXmaster
!==============================================================
! master subroutine to allocate XXY matrix elements
!
  subroutine XXYmaster(itx)


  use spstate
  use system_parameters
  use haiku_info
  use interaction
  use tripdef
  use interactions3body
  use ntuple_info
  use nodeinfo
  use bmpi_mod
  implicit none

  integer :: itx, ity
  integer :: aerr

!----------- ALLOCATE ------
  if(itx == 1)then


     if(nmatppn >0)then
       allocate( hmatppn(nmatppn), stat=aerr)
       if(aerr /= 0) call memerror("XXYmaster 1")
       hmatppn(:) = 0.0

       if(iproc==0)print*,nmatppn,' PPN Matrix elements ! '

     else
       if(iproc==0)print*,' zero PPN matrix elements ! '
       stop
     endif 

  else

     if(nmatpnn >0)then
       allocate( hmatpnn(nmatpnn), stat=aerr)
       if(aerr /= 0) call memerror("XXYmaster 2")
       hmatpnn(:) = 0.0
       if(iproc==0)print*,nmatpnn,' PNN Matrix elements ! '
     else
       if(iproc==0)print*,' zero PNN matrix elements ! '
       stop
     endif 
  endif
  
  return
  end subroutine XXYmaster
!=================================================================

end module b3bme_mod

! This is here to break a module dependency loop
subroutine threebodysetup
   use b3bme_mod
   implicit none
   call threebodysetup_mod
end subroutine threebodysetup
