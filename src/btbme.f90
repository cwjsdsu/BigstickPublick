!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK
!
!  file BTBME.F90
!
module btbme_mod
contains
!==============================================================
!
!  ** prepare_to_uncoupleXXtbme -- routine to uncouple pp, nn matrix elements
!     CALLS:     maxwpair
!                count_create_pairs
!		     sortpair
!		     find_pair_start
!
!  ** count_create_pairs -- creates pairs of protons or neutrons
!!
!  ** sortpair -- sorts pairs by Jz, then parity, then W
!
!  ** find_pair_start -- finds where lists of pairs with same Jz, parity start
!                    and computes # of uncoupled tbmes
!
!  ** maxwpair -- computes max W a proton or neutron pair can have
!!
!  ** setup4obs  -- master routine to set up arrays for computing J^2, T^2
!       actually reuses Hamiltonian arrays for J^2, new array for T^2
!  CALLS:	setup_obs_tbmes
!	      uncoupleXXtbme
!	      count_uncouplepn
! 	      makespeme
!
!  ** setup_obs_tbmes -- computes COUPLED tbmes for J^2 or T^2
!
!  ** uncoupleXXtbme  -- actually uncouples pp, nn tbmes
!
!----------- GENERAL FUNCTIONS-------------
!
!  ** zeta -- = sqrt(1 + delta(i,j))
!
!=====================================================================
!  MOST RECENT EDIT: Nov 2011 by CWJ
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  subroutine prepare_to_uncoupleXXtbme
!
!  sets up elements needed to uncouple pp,nn matrix elements
!  specifically: XX2(it)
!      XX2(it)%npairs = number of pairs species it
!             %pair(n) = information on the pairs
!  		%meref(jz,par) = references point for matrix elements for pairs with Jz, par
!		%mestart(jz,par) = starting point in matrix elements
!
!  INPUT: it = 1 proton, =2 neutron
!
!  CALLED BY:
!   hoopmaster
!
!  CALLS: 
!       maxwpair
!	count_create_pairs
! 	sortpair
!	find_pair_start
!   findTRpairsXX
!
subroutine prepare_to_uncoupleXXtbme(it)

  use spstate
  use haiku_info
  use system_parameters
  use interaction
  use coupledmatrixelements
  use nodeinfo
  use verbosity
  use bmpi_mod
  use io
  use TRstuff
  use adv_tbme_info
  
  implicit none

  integer(4) :: ierr
  integer :: it
  integer(kind=8) :: isps,jsps
  integer(kind=8) :: ipair,jpair
  integer :: indx
  integer :: jzmax,jzmin
  integer :: m
  integer :: n
  integer :: i,j
  integer :: ith,jth
  integer :: isgn,jsgn
  integer :: maxw2
  integer(kind=8), pointer :: map(:)
  integer :: aerr

  if(Np(it) < 2  )return
  call maxwpair(it,maxw2)  
  call count_create_pairs(it,maxw2,.false.,npairXX(it),XX2(it)%pair)
  allocate( XX2(it)%pair( npairXX(it) ), stat=aerr)  
  if(aerr /= 0) call memerror("prepare_to_uncoupleXXtbme 1")
  call count_create_pairs(it,maxw2,.true.,npairXX(it),XX2(it)%pair)
    
!--------- NEXT, SORT THESE -------
  if(npairXX(it)==0)then
     nmatXX(it) = 0
     return
  endif  
  call sortpair(npairXX(it),XX2(it)%pair)
  if(sort_XXpairs_on_W)call sortXXpairsW(npairXX(it),XX2(it)%pair,.false.)

!  if(useTRphase)call findTRpairsXX(it)
  if(useTRphase)call findTRpairs_gen(it,npairXX(it),XX2(it)%pair)
  
!------------ MAP FROM ORIGINAL CODING TO SORTED 

  if(it == 1)then
     n = nhsps(it) + nhsps(-it)
     allocate( mappairPP( n*(n-1)/2), stat=aerr)
     if(aerr /= 0) call memerror("prepare_to_uncoupleXXtbme 2")
     map => mappairPP
  else
     n = nhsps(it) + nhsps(-it)
     allocate( mappairNN( n*(n-1)/2), stat=aerr)
     if(aerr /= 0) call memerror("prepare_to_uncoupleXXtbme 3")
     map => mappairNN
  endif
  map = 0

  do ipair = 1,npairXX(it)
     
     indx = XX2(it)%pair(ipair)%indx
     if ( map(indx) /=0 ) then
 !       if ( iproc == 0 ) then
           print*,' problem, pair already found XX ',it
           print*,ipair,indx
 !       end if
#ifdef _MPI	
        call BMPI_ABORT(MPI_COMM_WORLD,101,ierr)
#endif
        stop
     end if
     map(indx) = ipair
  end do

!-------------- SET UP START, STOP ARRAYS ------------



!jzmin = XX2(it)%pair(1)%m      !
!jzmax = XX2(it)%pair(npairXX(it))%m      !
! ------- SWAPPED IN 7.9.4 in preparation for time-reversal
  jzmax = XX2(it)%pair(1)%m      !
  jzmin = XX2(it)%pair(npairXX(it))%m      !

  if(jzmax < jzmin)then
     if ( iproc == 0 ) print*,it,npairXX(it),n
#ifdef _MPI	
     call BMPI_ABORT(MPI_COMM_WORLD,101,ierr)
#endif
     stop
  endif
  
  
  
  allocate(XX2(it)%meref(jzmin:jzmax,2),XX2(it)%mestart(jzmin:jzmax,2),XX2(it)%block(jzmin:jzmax,2), stat=aerr)
  if(aerr /= 0) call memerror("prepare_to_uncoupleXXtbme 4")

  call find_pair_start(npairXX(it),XX2(it)%pair,jzmin, jzmax, XX2(it)%meref,XX2(it)%mestart,XX2(it)%block,nmatXX(it))

  if ( iproc == 0 .and. verbose_uncouple ) then

     if(it == 1) then
   	  write(6,'(i14," uncoupled pp matrix elements (",f10.3," Gb)")')nmatxx(it),nmatxx(it)*4.e-9
   	  write(logfile,'(i14," uncoupled pp matrix elements (",f10.3," Gb)")')nmatxx(it),nmatxx(it)*4.e-9		 
!        print*,' There are ',nmatxx(it),' uncoupled pp matrix elements (',nmatxx(it)*4.e-9,' Gb )'
     else
      write(6,'(i14," uncoupled nn matrix elements (",f10.3," Gb)")')nmatxx(it),nmatxx(it)*4.e-9
      write(logfile,'(i14," uncoupled nn matrix elements (",f10.3," Gb)")')nmatxx(it),nmatxx(it)*4.e-9
		 
       ! print*,' There are ',nmatxx(it),' uncoupled nn matrix elements (',nmatxx(it)*4.e-9,'Gb )'
        
     endif
  end if
  return
end subroutine prepare_to_uncoupleXXtbme

!=============================================================
!
! subroutine count_create_pairs
!
! creates pairs of protons or neutrons
!
! INPUT:
!  it = species, 1 = proton, 2 = neutrons
!  maxw2 = max w that can be carried by this pair
!  create = flag; if = .false. then just count, else fills arrays pairXX
!
! OUTPUT:
!  npair = # of pairs (counted up with create = .false.)
!  pairXX(n) = information on pairs (when create = .true.)
!
!  called by:
!   prepare_to_uncoupleXXtbme
!
!  subroutines called:
!
   subroutine count_create_pairs(it,maxw2,create,npair,pairXX)

      use spstate
      use haiku_info
      use system_parameters
      use pairdef
      use ntuple_info
	  use interaction,only:useTRnew

      implicit none
      integer it
      integer maxw2
      logical create
      integer npair
      type (pair_qn), pointer :: pairXX(:) 

      integer asps,bsps
      integer ipair,jpair
      integer indx
      integer jzmax
      integer m,icount
      integer n
      integer a,b
      integer ath,bth
      integer asgn,bsgn

      integer mpair,parpair,wpair
      if(create)then
         pairXX(:)%m = -99  ! flag for error
      endif
      indx = 0 

      ipair = 0
      do asps = 2, nhsps(it)+nhsps(-it)
		  
         if(asps > nhsps(-it))then
            ath = it
            a = asps - nhsps(-it)
            asgn = 1
         else
            ath = -it
            a = asps
            asgn = -1
         endif

         do bsps = 1,asps-1
			 
            if(bsps > nhsps(-it))then
               bth = it
               b = bsps - nhsps(-it)
               bsgn = 1
            else
               bth = -it
               b = bsps
               bsgn = -1
            endif
            if (hspsqn(ath,a)%w+hspsqn(bth,b)%w > maxw2)cycle
!------------- RESTRICT CREATION OF PAIRS ON QUANTUM NUMBERS
!              added 9/09 by CWJ
!............................................................
            mpair = ( hspsqn(ath,a)%m + hspsqn(bth,b)%m)
			if(useTRnew .and. mpair < 0)cycle
            wpair = hspsqn(ath,a)%w+hspsqn(bth,b)%w
            parpair = hspsqn(ath,a)%par*hspsqn(bth,b)%par
            parpair = (3-  parpair)/2
            if(restrict_ntuples)then
               if(wpair > pairlim(it)%maxw .or. wpair < pairlim(it)%minw)cycle
               if( parpair > pairlim(it)%w(wpair)%parmax) cycle
               if( parpair < pairlim(it)%w(wpair)%parmin) cycle
               if( mpair > pairlim(it)%w(wpair)%par(parpair)%jzmax) cycle
               if( mpair < pairlim(it)%w(wpair)%par(parpair)%jzmin) cycle

            end if

!............................................................
            ipair = ipair + 1
            if(create)then
              if(pairXX(ipair)%m /= -99)then
                print*,' oops already occupied '
                print*,asps,bsps,ipair,pairXX(ipair)%m
                stop
              endif
              pairXX(ipair)%m = mpair/2 !( hspsqn(ith,i)%m + hspsqn(jth,j)%m)/2
              pairXX(ipair)%w = wpair ! hspsqn(ith,i)%w+hspsqn(jth,j)%w
              pairXX(ipair)%par = parpair !hspsqn(ith,i)%par*hspsqn(jth,j)%par
!              pairXX(ipair)%par = (3-  pairXX(ipair)%par)/2
              pairXX(ipair)%indx = (asps -1)*(asps-2)/2 + bsps  ! note isps > jsps 
              pairXX(ipair)%ia = a*asgn
              pairXX(ipair)%ib = b*bsgn
			  pairXX(ipair)%TRphase = (-1)**( (hspsqn(ath,a)%j+hspsqn(bth,b)%j )/2 )
			  pairXX(ipair)%TR = ipair
			  
!			  print*,'pair ',ipair,a,asps,b,bsps

            endif
         enddo  ! bsps
      enddo  ! asps
   if(.not.create)npair = ipair

   return
   end subroutine count_create_pairs

!=============================================================
!
!  subroutine count_create_pnpairs
!
!  create pairs p_a n_b
!
!  INPUT:
!   wpn:  max on change in w for a pn pair
!   create: logical flag : = .false. just count, = .true. fill array
!
!  OUTPUT:
!   npair = # of pn pairs
!   pair = derived type with pair info
!
   subroutine count_create_pnpairs(wpn,create,npair,pair)

      use spstate
      use haiku_info
      use system_parameters
      use pairdef
      use ntuple_info
      use interaction
	  use coupledmatrixelements,only:dens2bflag

      implicit none
      integer wpn
      logical create
      integer npair
      type (pair_qn), pointer :: pair(:) 

      integer asps,bsps
      integer ipair,jpair
      integer indx
      integer jzmax
      integer m,icount
!! interface --      integer parmult
      integer wpair,mpair,parpair
      integer n
      integer a,b
      integer ath,bth
      integer asgn,bsgn


      if(create)then
         pair(:)%m = -99  ! flag for error
      endif
      indx = 0 

      ipair = 0
      n = nhsps(2)+nhsps(-2)
      do asps = 1, nhsps(1)+nhsps(-1)
         if(asps > nhsps(-1))then
            a = asps - nhsps(-1)
            ath = 1
            asgn = 1
         else
            a = asps
            ath = -1
            asgn = -1
         endif
         do bsps = 1,nhsps(2)+nhsps(-2)
            if(bsps > nhsps(-2))then
               b = bsps - nhsps(-2)
               bth = 2
               bsgn = 1
            else
               b = bsps
               bth =-2
               bsgn = -1
            endif
            if (hspsqn(ath,a)%w+hspsqn(bth,b)%w > wpn)cycle
!------------- RESTRICT CREATION OF PAIRS ON QUANTUM NUMBERS
!              added 10/09 by CWJ
!............................................................
            mpair = ( hspsqn(ath,a)%m + hspsqn(bth,b)%m)
			if(useTRnew .and. mpair < 0 .and. .not. dens2bflag)cycle
            wpair = hspsqn(ath,a)%w+hspsqn(bth,b)%w
            parpair = hspsqn(ath,a)%par*hspsqn(bth,b)%par
            parpair = (3-  parpair)/2

            if(restrict_ntuples )then
               if(wpair > pnlim%maxw .or. wpair < pnlim%minw)cycle
               if( parpair > pnlim%w(wpair)%parmax) cycle
               if( parpair < pnlim%w(wpair)%parmin) cycle
               if( mpair > pnlim%w(wpair)%par(parpair)%jzmax) cycle
               if( mpair < pnlim%w(wpair)%par(parpair)%jzmin) cycle
            end if

!--------- OPTIONAL USE OF TIME-REVERSAL TO COMPRESS STORAGE--

            ipair = ipair + 1
            if(create)then

              if(pair(ipair)%m /= -99)then
                print*,' oops already occupied '
                print*,asps,bsps,ipair,pair(ipair)%m
                stop
              endif
              pair(ipair)%m = mpair/2 !( hspsqn(ith,i)%m + hspsqn(jth,j)%m)/2
              pair(ipair)%w = wpair ! hspsqn(ith,i)%w+hspsqn(jth,j)%w
              pair(ipair)%par = parpair ! hspsqn(ith,i)%par*hspsqn(jth,j)%par
!              pair(ipair)%par = (3-  pair(ipair)%par)/2
              pair(ipair)%indx = (asps-1)*n + bsps  ! 
              pair(ipair)%ia = a*asgn
              pair(ipair)%ib = b*bsgn
			  pair(ipair)%tr = ipair
			  pair(ipair)%TRphase=1
!              print*,isps,jsps,i*isgn,j*jsgn
            endif
         enddo  ! jsps
      enddo  ! isps
   if(.not.create)npair = ipair
   return
   end subroutine count_create_pnpairs

!=============================================================
!
! subroutine count_create_rhos
!
! creates density-like operators (a^+a) of protons or neutrons
!
! INPUT:
!  it = species, 1 = proton, 2 = neutrons
!  create = flag; if = .false. then just count, else fills arrays pairXX
!
! OUTPUT:
!  npair = # of pairs (counted up with create = .false.)
!  pairXX(n) = information on pairs (when create = .true.)
!
!  called by:
!   
!
!  subroutines called:
!
   subroutine count_create_rhos(it,create,nrho,rhoX)

      use spstate
      use haiku_info
      use system_parameters
      use pairdef
      use ntuple_info
	  use TRstuff
	  use qs_on_proc

      implicit none
      integer it
      integer maxw2
      logical create
      integer nrho
      type (pair_qn), pointer :: rhoX(:) 

      integer asps,bsps
      integer ipair,jpair
      integer indx
      integer jzmax
      integer m,icount
      integer n
      integer a,b
      integer ath,bth
      integer asgn,bsgn

      integer mpair,parpair,wpair
	  
!	  call survey_qs_X(it)
	  
!	  print*,' SURVEY 1'
      if(create)then
         rhoX(:)%m = -99  ! flag for error
      endif
      indx = 0 

      ipair = 0
      do asps = 1, nhsps(it)+nhsps(-it)
		  
         if(asps > nhsps(-it))then
            ath = it
            a = asps - nhsps(-it)
            asgn = 1
         else
            ath = -it
            a = asps
            asgn = -1
         endif

         do bsps = 1,nhsps(it)+nhsps(-it)
			 
            if(bsps > nhsps(-it))then
               bth = it
               b = bsps - nhsps(-it)
               bsgn = 1
            else
               bth = -it
               b = bsps
               bsgn = -1
            endif
!------------- RESTRICT CREATION OF PAIRS ON QUANTUM NUMBERS
!............................................................
            mpair = ( hspsqn(ath,a)%m - hspsqn(bth,b)%m)/2
            wpair = hspsqn(ath,a)%w-hspsqn(bth,b)%w
            parpair = hspsqn(ath,a)%par*hspsqn(bth,b)%par
            parpair = (3-  parpair)/2
			
			
!			if(create)print*,ipair+1,mpair,'testing '		,useTRphase	
			if(useTRphase .and. ((it ==1 .and. mpair < 0) .or. (it==2 .and. mpair > 0)))cycle
			
!			print*,check_allowed_qs_X(it,mpair,parpair,wpair),mpair
			if(.not.check_allowed_qs_X(it,mpair,parpair,wpair))cycle   ! check that quantum numbers are in list
!............................................................
            ipair = ipair + 1
            if(create)then
              if(rhoX(ipair)%m /= -99)then
                print*,' oops already occupied '
                print*,asps,bsps,ipair,rhoX(ipair)%m
                stop
              endif
              rhoX(ipair)%m = mpair
              rhoX(ipair)%w = wpair 
              rhoX(ipair)%par = parpair 
              rhoX(ipair)%indx =  (asps-1)*(nhsps(it)+nhsps(-it))+bsps !(asps -1)*(asps-2)/2 + bsps  ! note isps > jsps 
              rhoX(ipair)%ia = a*asgn
              rhoX(ipair)%ib = b*bsgn
			  rhoX(ipair)%TRphase = 1
			  rhoX(ipair)%TR = ipair
			  rhoX(ipair)%wb = hspsqn(bth,b)%w

            endif
         enddo  ! bsps
      enddo  ! asps
   if(.not.create)nrho = ipair
   
!   if(create)print*,' m values ',rhoX(:)%m
   
   return
   end subroutine count_create_rhos

!=============================================================

      subroutine sortpair(npair,pair)
!
!  sort pairs (ia,ib) first by
!    total M
!    parity
!    does NOT sort on W; not clear if needed
!
!  REVISED 10/2011 with modified HEAPSORT algorithm 
!  (adapted from Numerical Recipes) 
!
!  (optional) subroutine called: sortXXpairsW
!
      use pairdef
!	  use adv_tbme_info
	  use nodeinfo
      implicit none
      integer npair
      type (pair_qn) :: pair(npair)
      type (pair_qn) :: pairtmp
      integer i,j,k,l,n

      integer m,par,w

      integer istart,istop
      integer jstart,jstop
	  
	  logical :: sort_on_w = .true.
	  
	  logical :: check = .false.


!--------- REVISED: HEAPSORT ON Jz------
      l = npair/2+1
      n = npair
      do while (n > 0)
          if(l > 1)then
             l = l- 1
             pairtmp = pair(l)
          else
             pairtmp = pair(n)
             pair(n) = pair(1)
             n = n -1
             if(n == 1)then
                 pair(1) = pairtmp
                 exit
              end if
          end if
          i = l
          j = l+l
          do while ( j <= n)
             if( j < n)then
                if(pair(j)%m >  pair(j+1)%m)j = j+1  !ordering switched in 7.9.4 to prepare for time reversal
             endif
             if( pairtmp%m > pair(j)%m)then !ordering switched in 7.9.4 to prepare for time reversal
                 pair(i) = pair(j)
                 i = j
                 j = j + j
             else
                 j = n + 1
             end if
          end do      
          pair(i) = pairtmp
      end do
!------------- next sort on parity ----------

!---------- find start, stop for each M value
      m = pair(1)%m
      istart = 1
      istop  = 1
 !     do while(m >= pair(npair)%m )
      do while(istart <= npair)
        if(istart > npair)then
          print*,' problem in sorting pairs '
          print*,m, pair(npair)%m
          print*,pair(npair)%ia,pair(npair)%ib
          stop
        endif
		m = pair(istart)%m
        if(istart == npair)then
          istop = istart
        else
          do while(pair(istop+1)%m == m )
           istop = istop+1
           if(istop == npair) goto 3
          enddo
3         continue
        endif
		
!--------- PERFORM ACTUAL SORT --------

        l = (istop-istart+1)/2+istart
        n = istop
        do while (n >= istart)
          if(l > istart)then
             l = l- 1
             pairtmp = pair(l)
          else
             pairtmp = pair(n)
             pair(n) = pair(istart)
             n = n -1
             if(n == istart)then
                 pair(istart) = pairtmp
                 exit
              end if
          end if
          i = l
          j = l+l-(istart-1)
          do while ( j <= n)
             if( j < n)then
                   if(pair(j)%par< pair(j+1)%par)j = j+1
             endif
             if( pairtmp%par < pair(j)%par)then 
                 pair(i) = pair(j)
                 i = j
                 j = j + j-(istart-1)
             else
                 j = n + 1
             end if
          end do      
          pair(i) = pairtmp
      end do
	  
	  if(check)then
		  print*,istart,istop
		  do i=istart,istop-1
			  if(pair(i)%m/=pair(i+1)%m)then
				  print*,' not sorted on M '
				  stop
			  end if
			  if(pair(i)%par > pair(i+1)%par)then
				  print*,' not sorted on parity '
				  stop
			  end if
			  
		  end do
		  
		  
	  end if
!
! IF ONE WERE TO SORT ON W, IT WOULD BE HERE
!-----------------------------------
        m = m+1
        istart = istop+1
        istop = istart
      enddo  
!	  if(sort_XXpairs_on_W)call sortXXpairsW(npair,pair,.false.)
	  
	!  if(sort_XXpairs_on_W.and. nprocs > 1)call sortXXpairsW(npair,pair)
	
!	do i = 1,npair
!		print*,i,pair(i)%m,pair(i)%par!,pair(jpair)%w,pair(jpair)%wb

!	end do

      return
      end subroutine sortpair

!==========================================================
!
!  subroutine find_pair_start
!
!  REVISED 7.9.1 to be more rational
!
!  finds, in a list of pairs, where quantum numbers start, end
!  also computes the number of matrix elements
!
!  INPUT:
!     npair = # of pairs
!     pair(:) = derived type of pairs
!     jzmax   = max value for jz for these pairs
!  OUTPUT:
!     mestart(jz,par) = where in the list pairs with jz,par start
!     meref(jz,par) = reference mapping
!     block(jz,par) = # of pairs in a block of jz, par (ADDED 7.9.1)
!     nme = total # of matrix elements
!
! CALLED BY  prepare_to_uncoupleXXtbme
!
      subroutine find_pair_start(npair,pair,jzmin,jzmax, meref,mestart,block,nme)
 
      use pairdef
	  use coupledmatrixelements,only:dens2bflag
	  use TRstuff
      implicit none

      integer npair
      type (pair_qn) :: pair(npair)
 
      integer jzmin,jzmax
      integer(8) :: mestart(jzmin:jzmax,2)
      integer(8) :: meref(jzmin:jzmax,2)
      integer(8) :: block(jzmin:jzmax,2)

      integer(8) :: nme  !  ultimate # of matrix elements
      integer istart,iend
      integer m,i,j,iref
      integer par
	  
	  
	  integer :: jstop
	  

! --- COUNT UP # OF PAIRS --- ADDED IN 7.9.1
!     only really needed when doing 2-body densities
	  block(:,:)=0
	  do i = 1,npair
		  m= pair(i)%m
		  par=pair(i)%par
		  block(m,par)=block(m,par)+1
	  end do
	  nme=0
	  iref = 0
!      do m=jzmin,jzmax
	  do m=jzmax,jzmin, -1  ! order switched in 7.9.4 to prepare for using time-reversal
		  do par = 1,2
              if(m >= 0 .or. .not. useTRphase)then
			     meref(m,par)=iref
			     mestart(m,par)=nme
			     if(.not.dens2bflag)then
			        nme=nme+block(m,par)*(block(m,par)+1)/2
			     else

		            nme=nme+block(m,par)**2
			      end if
!			  print*,m,par,block(m,par),' block size '
			  
			      iref = iref+block(m,par)
			   else  ! useTRphase and m < 0 
				  meref(m,par)= meref(-m,par)
				  mestart(m,par)=mestart(-m,par)
			  end if
		  end do
	  end do	  

!	  print*,nme,' matrix elements found'

      return
      end subroutine find_pair_start

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  subroutine maxwpair
!  finds the maximum W allowed for a pair
!  by comparing max W for N, N-2 particles
!  
   subroutine maxwpair(it,maxw2)

   use spstate
   use sporbit
   use W_info
   use haiku_info
   use system_parameters

   implicit none
   integer it
   integer maxw2
   
   integer w
   integer i
   maxw2 = 0
   if(Np(it) < 2 )return

!   if( phconj(it)) then  ! then we must allow all pairs
!      w = 0
!      do i = 1,numorb(it)
!         w = bmax(orbqn(it,i)%w,w)
!      end do
!      maxw2 = 2*w
!      return
!   end if

   w = 0
   do i =1,Np(it)-2
        w = w + spsqn(it,i)%w
   enddo
   maxw2 = maxW(it) - w
   return
   end subroutine maxwpair
!===================================================
!  subroutine setup4obsmaster
!
!  sets up for either J^2 or T^2 observable
!  introduced in 7.3.7 
!
!  obschar = 'J' create J^2
!  obschar = 'T' create T^2
!  NEW in 7.11.3: obschar='K' create J^2, but do not uncouple
!
!  CALLS:
!	setup_obs_tbmes
!       pandamaster
!	uncoupleXXtbme
!	count_uncouplepn

   subroutine setup4obsmaster(obschar)
   use diagh
   use sp_potter
   use coupledmatrixelements
   implicit none
   integer i,j
   character(1) obschar
   
   call setup4obsspe(obschar)
   
!------------ SET UP TWO-BODY MATRIX ELEMENTS
   call setup_obs_tbmes(obschar)

   if(obschar=='T')then
     call pandamaster
   end if
   if(.not.call_spe)call subsume_sp_pot_into_2body
   
   if(obschar=='K')return

   call uncoupleXXtbme(1,obschar)
   call uncoupleXXtbme(2,obschar)
   call count_uncouplepn(obschar,.true.)

   if(call_spe)then
      call makespeme(1,obschar)
      call makespeme(2,obschar)
   end if
   end subroutine setup4obsmaster

!===================================================
!
!  subroutine setup4obsspe
!
!  sets up for either J^2 or T^2 observable
!  MODIFIED in V7.2.5 to change T^2 for p-h conjugation
!  MODIFIED in 7.3.7 to handle just spes
!  INPUT:
!   obschar = 'J' or 'T'
!
!  CALLS:

!	makespeme
!
   subroutine setup4obsspe(obschar)
   use system_parameters
   use coupledmatrixelements
   use sporbit
   use spstate
   use descendents
   implicit none
   integer i,j
   character(1) obschar
   integer :: tmax2  ! maximal 2T
   
!---------- NEED TO ZERO OUT ONE-BODY POTENTIAL
   
   nsppot = 0.0
   psppot = 0.0

!----------- NOW SET UP SINGLE-PARTICLE

   select case (obschar)

   case ('J','K')
   do i = 1,numorb(1)
     j = orbqn(1,i)%j
     pspe(i) = 0.25*j*(j+2)
   enddo

   do i = 1,numorb(2)
     j = orbqn(2,i)%j
     nspe(i) = 0.25*j*(j+2)
   enddo

   case ('T')
   if(.not.phconj(1) .or. phconj(2)) then
      do i = 1,numorb(1)
         pspe(i) = 0.75
      enddo
   else
      Tmax2 = nsps(1)
      do i = 1,numorb(1)
         pspe(i) = 0.75 ! 0.25*Tmax2*(Tmax2+2)/float(np(1))-0.75 -0.5*(Tmax2-1)
      enddo
   end if

   if(.not.phconj(2) .or. phconj(1))then 
      do i = 1,numorb(2)
         nspe(i) = 0.75
      enddo
   else
      Tmax2 = nsps(2)
      do i = 1,numorb(1)
         nspe(i) = 0.75 ! 0.25*Tmax2*(Tmax2+2)/float(np(2)) -0.75 -0.5*(Tmax2-1)
      enddo
   end if
   end select
   
   return
   end subroutine setup4obsspe

!===================================================

!
!  sets up arrays of COUPLED tbmes for observable either J or T
!  
!  CALLED BY:  setup4obsmaster
!
  subroutine setup_obs_tbmes(obschar)
  use system_parameters
  use coupledmatrixelements
  use interaction
  use sporbit

  implicit none
  integer it
  integer i
  integer ia,ib
  integer a,b
  integer pair1, pair2
  integer indx, altindx
  integer jmax,jmin,j    
  integer t
  logical same,okalt
  real sj2i  ! six-j called from libra.f
  real v, altv
  integer ja,jb
  real xja,xjb
  character(1) :: obschar  ! = J or T
  real  vtmp
  real factor
  integer K
  type (vjs), pointer :: tbme(:)
  logical jeven
  integer :: ipar,iref,istart
  integer :: aerr

  if(obschar == 'T')then
    if(.not.allocated(obsmatpn)) then
       allocate(obsmatpn(nmatpn), stat=aerr)
       if(aerr /= 0) call memerror("setup_obs_tbmes")
    end if
  endif

!--------------PP/NN INTERACTIONS -------------------
  
  do it = 1,2
  if(np(it) < 2 .and. .not. phconj(it))cycle
  if(it == 1)then
     tbme     => ppme
  else
     tbme     => nnme
  endif
  do i = 1,nv2bmedim(it)
    if(tbme(i)%jmin <= tbme(i)%jmax)tbme(i)%v(:) = 0.0
  enddo

!--------------PP/NN INTERACTIONS -------------------
  do pair1 = 1,ncouplesXX(it)
        a = XXcouples(it)%pairc(pair1)%ia
        b = XXcouples(it)%pairc(pair1)%ib
        ja = orbqn(it,a)%j
        xja = float(ja)*0.5
        jb = orbqn(it,b)%j
        xjb = float(jb)*0.5

        if(a == b )then
             same = .true.
        else
             same = .false.
        endif
        ipar = XXcouples(it)%pairc(pair1)%par
		
        iref = XXcouples(it)%meref(ipar)
        istart=XXcouples(it)%mestart(ipar)
        indx = (pair1-iref)*(pair1-iref-1)/2+pair1-iref+istart
        jmax = (ja+jb)/2
        jmin = abs(ja-jb)/2

        if(jmin < tbme(indx)%jmin .or. jmax > tbme(indx)%jmax)then
               print*,pair1,indx
               print*,jmin,tbme(indx)%jmin
               print*,jmax,tbme(indx)%jmax
               stop
         endif
         do j = jmin,jmax
             t = 1   ! required for pp , nn
             if(same)then
               if(.not.spinless .and. ((j+t)/2)*2 == j+t)cycle  ! j+t must be odd if a = b
               if(spinless .and. ((j+t)/2)*2 /= j+t)cycle  ! j+t must be even if a = b
             endif
             v = 0.0   ! prevent uninit msg
             select case (obschar)
             case ('J','K')
                v = float(j*(j+1)) - xja*(xja+1) - xjb*(xjb+1)
             case ('T')
                v = 0.5
             end select

             tbme(indx)%v(j) = v

         enddo

  enddo ! pair1

!......................................................
  enddo ! it
!------------------------------- PN-------------------
  tbme => pnme
  do i = 1,nv2bmedim(0)
    if(tbme(i)%jmin <= tbme(i)%jmax)tbme(i)%v(:) = 0.0
  enddo

  do pair1 = 1,ncouplesPN
        a = PNcouples%pairc(pair1)%ia
        b = PNcouples%pairc(pair1)%ib
        ja = orbqn(1,a)%j
        xja = float(ja)*0.5
        jb = orbqn(2,b)%j
        xjb = float(jb)*0.5


        if( ja == jb .and. orbqn(1,a)%l == orbqn(2,b)%l & 
                  .and. orbqn(1,a)%nr == orbqn(2,b)%nr)then
              factor = 1.0
        else
              factor = 0.5
        endif

        if(ja == jb .and. orbqn(1,a)%l == orbqn(2,b)%l .and. & 
              orbqn(1,a)%nr  == orbqn(2,b)%nr )then
             same = .true.
        else
             same = .false.
        endif
        ipar = PNcouples%pairc(pair1)%par
        iref = PNcouples%meref(ipar)
        istart=PNcouples%mestart(ipar)
        indx = (pair1-iref)*(pair1-iref-1)/2+pair1-iref+istart
           
        if(.not.same .and. obschar=='T')then
!............. EXCHANGE TERM FOR COMPUTING ISOSPIN; do not use if isospin violating
              pair2 = numorb(2)*(b-1)+a  ! actually, this is not quite right; need to map
              pair2 = PNcouplemap(pair2)
              if(pair2 == -1)then
                  okalt = .false.
              else
                  okalt = .true.
                 if(pair1 > pair2)then
                    altindx = (pair1-iref)*(pair1-iref-1)/2+pair2-iref+istart
                 else
                    altindx = (pair2-iref)*(pair2-iref-1)/2+pair1-iref+istart
                 endif
              endif
        else
              okalt = .false.
        endif
        jmax = (ja+jb)/2
        jmin = abs(ja-jb)/2

        if(indx > nv2bmedim(0))then
                 print*,indx,nv2bmedim(0)
                 stop
        endif
        if(.not. same .and. okalt .and. altindx > nv2bmedim(0))then
                 print*,altindx,nv2bmedim(0), ' alt '
                 print*,ia,ib,a,b,indx
                 print*,pair1,pair2
                 stop
        endif
        if( indx == 0) then
                 print*,indx
                 stop
        endif
        if(.not.same .and. okalt .and. altindx == 0)then
                 print*,altindx, ' ALT ',indx
                 print*,pair1,pair2
                 print*,ia,ib,a,b,same
                 stop
        endif

        if(jmin < tbme(indx)%jmin .or. jmax > tbme(indx)%jmax)then
               print*,pair1,indx,'(pn)'
               print*,ja,jb
               print*,a,b,ia,ib
               print*,jmin,tbme(indx)%jmin
               print*,jmax,tbme(indx)%jmax
               stop

        endif
        if(.not.same .and. okalt)then
              if((jmin < tbme(altindx)%jmin .or. jmax > tbme(altindx)%jmax))then
               print*,pair1,indx,'(pn) (ALT) '
               print*,ja,jb
               print*,a,b,ia,ib
               print*,jmin,tbme(indx)%jmin
               print*,jmax,tbme(indx)%jmax
               stop
              endif
        endif

           do j = jmin,jmax

             altv = 0.0 ! prevent uninitialized msg
             v = 0.0
             select case(obschar)
             case ('J','K')
               v = float(j*(j+1)) - xja*(xja+1) - xjb*(xjb+1)
             case ('T')
               if(same)then
                  if( (j /2 )*2 ==j)then
                    jeven = .true.
                  else
                    jeven = .false.
                  endif
                  if( (spinless .and..not. jeven) .or. (.not.spinless .and. jeven))then  !  
                     v = 0.5

                  else
                     v = -1.5

                  endif
               else
                  v = -0.5
                 altv = -(-1)**(J +(ja + jb)/2)
               endif

             end select

               tbme(indx)%v(j) =  v  !*factor
               if(.not.same  .and. obschar == 'T' .and. okalt)then
                 tbme(altindx)%v(j) = altv
               endif
!            enddo  ! t
           enddo  ! j

  enddo ! pair1

!........ IF PARTICLE-HOLE CONJUGATION, MUST CONVERT PN.............
!  if( (phconj(1) .and. .not.phconj(2)) .or. (phconj(2) .and. .not.phconj(1)) )then
!     if(obschar=='T')then
!        print*,' SWAPPING ISOSPIN PN ',obschar
!        call pandacross
!     end if
!  end if
  return
  end subroutine setup_obs_tbmes


!==========================================================
!
!      subroutine uncoupleXXtbme
!
!  uncouples the T =1 matrix elements for pp, nn
!  includes the factor zeta = sqrt(1+delta(a,b)) etc.
!
!  IMPORTANT: How pp/nn matrix elements are encoded; 
!  needed for setting up the operator arrays for 
!  two-body "jumps"
!
!  INPUT:
!   it = species, 1 = proton, 2 = neutrons
!   hchar = 'H' hamiltonian, use all matrix elements;
!         = 'J' or 'T' only a subset are nonzero, so faster loop
!
!  SUBROUTINES CALLED: readXXdeformedtbmes  (only for special case of nonspherial interactions)
!
      subroutine uncoupleXXtbme(it,hchar)

      use verbosity
      use spstate
      use haiku_info
      use system_parameters
      use interaction
      use coupledmatrixelements
      use interaction_deformed
      use butil_mod
	  use TRstuff

      implicit none
      integer it      ! which species
      character(1) hchar  ! signals what is operator
                          ! default is ordinary hamiltonian
                          ! 'J' = ang momentum
                          ! 'T' = isospin
      integer(kind=8) :: dpair,cpair
      integer Jtot,jmin,jmax
      integer T
      integer M
      integer par
      integer ik,il
      integer ii,ij

      integer i,j,k,l
      integer ith,jth,kth,lth
      integer isps,jsps,ksps,lsps

      integer nk,nl,ni,nj
      integer jk,jl,ji,jj
      integer mk,ml,mi,mj

      integer(kind=8) :: iref,istart
      integer(kind=8) :: itbme

      integer ia,ib,ic,id
      integer pair1,pair2
      integer pair11,pair22
      integer phasekl,phaseij

      integer indx
      integer itmp
      real vtmp

      !! real zeta   ! zeta(i,j) = sqrt(1+delta(i,j))
      real cleb
      real hfactor  ! FOR HERMICIITY
      real,pointer :: v(:)
      real, target :: v_dummy(1)

      type (vjs), pointer :: tbme(:)
      integer pcref, pcstart, ppar
      integer,pointer :: xxmap(:)
      logical :: foundit
      integer :: aerr

      if(Np(it) < 2)return
      v => v_dummy

!-----  SPECIAL CASE FOR UNCOUPLED/DEFORMED
!       NOT USUALLY CALLED

      if(hchar=='H' .and. deformed)then
          call readXXdeformedtbmes(it,foundit)
          if(foundit)return
      end if

!-------------- DEPENDING ON VALUE OF hchar SWITCH BETWEEN ARRAYS
      select case (hchar)
      case ('T')
      if(it == 1)then
        if(.not.allocated(obsmatpp)) then
          allocate(obsmatpp(nmatXX(it)), stat=aerr)
          if(aerr /= 0) call memerror("uncoupleXXtbme 1")
        end if
        v => obsmatpp
      else
        if(.not.allocated(obsmatnn)) then
          allocate(obsmatnn(nmatXX(it)), stat=aerr)
          if(aerr /= 0) call memerror("uncoupleXXtbme 2")
        end if
        v => obsmatnn
      endif
      case ('H', 'J')
!------------- I LET THE J^2 ARRAY REPLACE THE HAMILTONIAN MEs---------

      if(it ==1)then
        if(.not.allocated(hmatpp)) then
           allocate(hmatpp(nmatXX(it)), stat=aerr)
          if(aerr /= 0) call memerror("uncoupleXXtbme 3")
        end if
        v => hmatpp
      else
        if(.not.allocated(hmatnn)) then
           allocate(hmatnn(nmatXX(it)), stat=aerr)
          if(aerr /= 0) call memerror("uncoupleXXtbme 4")
        end if
        v => hmatnn

      endif

      end select

      if(it == 1)then
        tbme => ppme
        xxmap => ppcouplemap

      else
        tbme =>   nnme
        xxmap => nncouplemap

      endif

      v=0.0   ! zero out array

      T = 1        ! for pp, nn matrix elements
	  
	  
!--- NOTE THERE IS SOME BUG IN THE OPENMP WHICH I CANNOT TRACK DOWN--

!!$omp parallel do private (m,par,ppar,iref,istart,cpair,dpair,itbme, & 
!!$omp                     ia,ib,ic,id,i,j,k,l,mi,mj,mk,ml,ji,jj,jk,jl, &
!!$omp                      ith,jth,kth,lth,isps,jsps,ksps,lsps,phaseij,phasekl, &
!!$omp                      itmp,pcref,pcstart,indx,jmin,jmax,jtot, &
!!$omp                      pair1,pair2,pair11,pair22,hfactor,vtmp), &
!!$omp           shared(v,XX2,hspsqn,nmatXX,tbme,xxmap,hmatpp,hmatnn,npairXX)


      do dpair = 1,npairXX(it)        ! loop over destruction pair

        m = XX2(it)%pair(dpair)%m     ! find Jz, parity of pair
		if(useTRphase .and. m < 0)cycle  
		
        par = XX2(it)%pair(dpair)%par

        iref = XX2(it)%meref(m,par)      ! iref, istart used to encode index for 
        istart = XX2(it)%mestart(m,par)  ! the uncoupled matrix element

        if(iref > dpair)then
          print*,' problem with XX iref (d) ',it
          print*,dpair,iref
          print*,m,par
          print*,xx2(it)%mestart(:,1)
          stop
        endif

!---------------------- GET Q#s FOR PARTICLES IN DESTRUCTION PAIR --------
        k = XX2(it)%pair(dpair)%ia
        l = XX2(it)%pair(dpair)%ib
        if( k < 0)then
          ksps = -k
          kth = -it
        else
          ksps = k
          kth  = it
        endif
        ia = hspsqn(kth,ksps)%orb
        jk = hspsqn(kth,ksps)%j
        mk = hspsqn(kth,ksps)%m

        if( l < 0)then
          lsps = -l
          lth = -it
        else
          lsps = l
          lth  = it
        endif
        ib = hspsqn(lth,lsps)%orb
        jl = hspsqn(lth,lsps)%j
        ml = hspsqn(lth,lsps)%m

        phasekl = 1
!-------------- MUST HAVE ia > ib, else swap and pick up a phase ---

        if(ia < ib)then
             itmp = ia
             ia = ib
             ib = itmp

             itmp = jk
             jk = jl
             jl = itmp

             itmp = mk
             mk   = ml
             ml   = itmp

             phasekl = -1  ! the rest of the phase comes from the Clebsch Gordan

        endif

        do cpair = 1,dpair        ! loop over creation pair  ! note triangular storage
          if(cpair == dpair)then
              hfactor = 0.5
          else
              hfactor = 1.0
          endif
          if(XX2(it)%pair(cpair)%m /=m)cycle       ! enforce quantum numbers
          if(XX2(it)%pair(cpair)%par /= par)cycle

          if(iref > cpair)then
            print*,' problem with XX iref (c) '
            print*,cpair,iref
            stop
          endif

!--------------- HOW INDEX FOR UNCOUPLED (m-scheme) tbme IS COMPUTED---------

          itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref
		  
          if(itbme > nmatXX(it))then          ! error trap
             print*,' me label too large'
             print*,itbme,nmatXX(it)
             print*,cpair,dpair,iref,istart
             stop
          endif
          if(itbme <= 0)then
               print*,' problem itbme (xx) ',it,itbme
			   print*,istart,dpair,iref,cpair,iref
               stop
          endif

!--------------- GET Q#s FOR CREATION PAIRS ---------------

          i = XX2(it)%pair(cpair)%ia
          j = XX2(it)%pair(cpair)%ib

          if( i < 0)then
            isps = -i
            ith = -it
          else
            isps = i
            ith  = it
          endif
          ic = hspsqn(ith,isps)%orb
          ji = hspsqn(ith,isps)%j
          mi = hspsqn(ith,isps)%m

          if( j < 0)then
            jsps = -j
            jth = -it
          else
            jsps = j
            jth  = it
          endif
          id = hspsqn(jth,jsps)%orb
          jj = hspsqn(jth,jsps)%j
          mj = hspsqn(jth,jsps)%m		  
!------------- PUT INTO "STANDARD" ORDER; MAY GET PHASE -------

          phaseij = 1

          if(ic < id) then

             itmp = ic
             ic = id
             id = itmp

             itmp = ji
             ji = jj
             jj = itmp
             itmp = mi
             mi   = mj
             mj   = itmp

             phaseij = -1  ! the rest of the phase comes from the Clebsch Gordan

          endif
          pair1 = ia*(ia-1)/2 + ib  ! used to find index of coupled TBME

          pair2 = ic*(ic-1)/2 + id

!----------------- FOR J^2, T^2, must have pair1 = pair2

          if(hchar == 'J'.and. pair1/=pair2)cycle
          if(hchar == 'T'.and. pair1/=pair2)cycle

!-------------- ALSO MUST HAVE STANDARD ORDERING OF PAIRS -----------

          pair1 = xxmap(pair1)
          pair2 = xxmap(pair2)
          if(pair1==-1 .or. pair2 == -1)then
             print*,' hmm strange ... (A)',it
             print*,pair1,pair2 
             print*,ia,ib,ic,id
             print*,xxmap
             stop
          endif
          ppar = XXcouples(it)%pairc(pair1)%par
          if(ppar /= XXcouples(it)%pairc(pair2)%par)then
              print*,' sigh parity problem '
              print*,ia,ib,ic,id,ppar,par
              print*,pair1,pair2
              stop
          endif
          if(pair1 < pair2)then
            pair11 = pair2
            pair22 = pair1
          else
           pair11 = pair1
           pair22 = pair2
          endif
          pcref = XXcouples(it)%meref(ppar)
          pcstart = XXcouples(it)%mestart(ppar)
! this is the indx for coupled (J-scheme) matrix elements
          indx = (pair11-pcref)*(pair11-1-pcref)/2+pair22-pcref+pcstart
          if(indx <= 0)then
             print*,' problem indx ',indx
             stop
          endif

!-------------FIND JMIN,JMAX ALLOWED ----------------

          jmax = bmin( jk+jl,ji+jj) /2
          jmin = bmax( abs(jk-jl),abs(ji-jj))/2
          jmin = bmax(jmin,abs(m))

!-------------NOW EXTRACT FROM TBMEs ---------------
!-------      NOTE FACTOR ZETA = SQRT(1+DELTA(A,B)) INCLUDED
          vtmp = 0.
          if(v(itbme)/=0.)then
            print*,' already here (XX) ',it
            print*,itbme,nmatXX(it),v(itbme)
            print*,dpair,cpair,iref,istart
            stop
          endif

          do jtot = jmin,jmax

             vtmp = vtmp+tbme(indx)%v(jtot) *zeta(ia,ib)*zeta(ic,id)   * &
                    cleb(jk,mk,jl,ml,2*jtot,2*m)*cleb(ji,mi,jj,mj,2*jtot,2*m)
          enddo           
          v(itbme)  = vtmp*phasekl*phaseij*hfactor
		  
        enddo  ! loop over cpair

      enddo  ! loop over dpair
!!$omp end parallel do

      return
      end subroutine uncoupleXXtbme
!
!==========================================================

!  specialized subroutine--replaces uncoupled matrix elements
!  with random numbers
!

subroutine randomize_tbme
	use interaction
	implicit none
	integer ime
	real rv
	if(.not. random_mscheme)return   ! flag set in module interaction
	
	print*,'  '
	print*,' Attention! Randomizing uncoupled two-body matrix elements '
	print*,' '
		
	do ime = 1,nmatXX(1)
		call random_number(rv)
		hmatpp(ime) = 2.*(rv-0.5)   ! randomize pp
	end do !ime
	
	do ime = 1,nmatXX(2)
		call random_number(rv)
		hmatnn(ime) = 2.*(rv-0.5)   ! randomize nn
	end do !ime
	
	do ime = 1,nmatpn
		call random_number(rv)
		hmatpn(ime) = 2.*(rv-0.5)   ! randomize pn
	end do !ime
	
	return
	
end subroutine randomize_tbme

!=============================================================

end module btbme_mod
