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
  use onebodypot, only: meanie
  
  implicit none

  !! include 'binterfaces.inc'

  integer(4) :: ierr
  integer :: it
  integer(kind=8) :: isps,jsps
  integer(kind=8) :: ipair,jpair
  integer :: indx
  integer :: jzmax,jzmin
  integer :: m
! interface  integer :: parmult
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

!------------ MAP FROM ORIGINAL CODING TO SORTED 

  if(it == 1)then
     n = nhsps(it) + nhsps(-it)
!       print*,n
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
        if ( iproc == 0 ) then
           print*,' problem, pair already found '
           print*,ipair,indx
        end if
        call BMPI_ABORT(icomm,101,ierr)
        stop
     end if
     map(indx) = ipair
  end do

!-------------- SET UP START, STOP ARRAYS ------------

  jzmin = XX2(it)%pair(1)%m      !

  jzmax = XX2(it)%pair(npairXX(it))%m      !

  if(jzmax < jzmin)then
     if ( iproc == 0 ) print*,it,npairXX(it),n
     call BMPI_ABORT(icomm,101,ierr)
     stop
  endif
  
  
  allocate(XX2(it)%meref(jzmin:jzmax,2),XX2(it)%mestart(jzmin:jzmax,2), stat=aerr)
  if(aerr /= 0) call memerror("prepare_to_uncoupleXXtbme 4")

  call find_pair_start(npairXX(it),XX2(it)%pair,jzmin, jzmax, XX2(it)%meref,XX2(it)%mestart,nmatXX(it))
  if(meanie)return  ! doing self-consistent mean-field, do not need the rest

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
   subroutine count_create_pairs(it,maxw2,create,npair,pairXX)

      use spstate
      use haiku_info
      use system_parameters
!      use interaction
!      use coupledmatrixelements
      use pairdef
      use ntuple_info

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
! interface      integer parmult
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
!      print*,nhsps(-1),nhsps(1),nhsps(-2),nhsps(2)
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
            if(useTR .and. Jz >=0 .and. mpair < 0 )cycle
            if(useTR .and. Jz <0 .and. mpair > 0 )cycle

!............................................................


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
!              print*,isps,jsps,i*isgn,j*jsgn
            endif
         enddo  ! jsps
      enddo  ! isps
   if(.not.create)npair = ipair
   return
   end subroutine count_create_pnpairs

!===================================================

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
      use pairdef
      implicit none
      integer npair
      type (pair_qn) :: pair(npair)
      type (pair_qn) :: pairtmp
      integer i,j,k,l,n

      integer m,par,w

      integer istart,istop
      integer jstart,jstop


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
                if(pair(j)%m< pair(j+1)%m)j = j+1
             endif
             if( pairtmp%m < pair(j)%m)then 
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
      do while(m <= pair(npair)%m )
        if(istart > npair)then
          print*,' problem in sorting pairs '
          print*,m, pair(npair)%m
          print*,pair(npair)%ia,pair(npair)%ib
          stop
        endif
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
!
! IF ONE WERE TO SORT ON W, IT WOULD BE HERE
!
!-----------------------------------
        m = m+1
        istart = istop+1
        istop = istart
      enddo  

      return
      end subroutine sortpair

!==========================================================
!
!  subroutine find_pair_start
!
!  finds, in a list of pairs, where quantum numbers start, end
!
!  INPUT:
!     npair = # of pairs
!     pair(:) = derived type of pairs
!     jzmax   = max value for jz for these pairs
!  OUTPUT:
!     mestart(jz,par) = where in the list pairs with jz,par start
!     meref(jz,par) = reference mapping
!     nme = total # of matrix elements
!
      subroutine find_pair_start(npair,pair,jzmin,jzmax, meref,mestart,nme)
 
      use pairdef
      implicit none

      integer npair
      type (pair_qn) :: pair(npair)
 
      integer jzmin,jzmax
      integer(8) :: mestart(jzmin:jzmax,2)
      integer(8) :: meref(jzmin:jzmax,2)

      integer(8) :: nme  !  ultimate # of matrix elements

      integer istart,iend
      integer m,i,j,iref
      integer par

      nme = 1

      meref(jzmin,:) = 0
      mestart(jzmin,:) = 0
      iref = 1
      do i = 2,npair
        m = pair(i)%m
        par = pair(i)%par
        if(m > pair(i-1)%m .or. par /= pair(i-1)%par )then
           iref = i
           meref(m,par) = iref -1
           mestart(m,par) = nme
        endif 
        do j = iref,i
            nme = nme+1
        enddo
      enddo

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
!  CALLS:
!	setup_obs_tbmes
!       pandamaster
!	uncoupleXXtbme
!	count_uncouplepn

   subroutine setup4obsmaster(obschar)
   use diagh
   use sp_potter
!   use flagger,only:subsume_spe
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

   case ('J')
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
             case ('J')
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
             case ('J')
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
	  use onebodypot, only: meanie

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

      if(meanie)return   ! doing self-consistent mean-field
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

!--------------- HOW INDEX FOR UNCOUPLED tbme IS COMPUTED---------

          itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref
          if(itbme > nmatXX(it))then          ! error trap
             print*,' me label too large'
             print*,itbme,nmatXX(it)
             print*,cpair,dpair,iref,istart
             stop
          endif
          if(itbme <= 0)then
               print*,' problem itbme ',itbme
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

!==========================================================
!
! routine to fetch standalone uncoupled XX matrix element
! used when storing uncoupled matrix elements directly into jumps
!
! INPUT:
!    it:  species 1 = proton, 2 = neutron
!
! OUTPUT :
!    xxme : uncoupled 2-body PP/NN matrix element
!

subroutine fetchXXme(it,i,j,k,l,xxme)
   use interaction
   use coupledmatrixelements
   use haiku_info
   use spstate
   use butil_mod
   implicit none
!..... INPUT.....
   integer :: it    ! species
   integer :: i,j,k,l  ! s.p. state labels;  if < 0, then 'left-handed' else 'right-handed'
!..... OUTPUT.....
   real    :: xxme
!...... INTERMEDIATE.........
   integer :: isps,jsps,ksps,lsps
   integer :: ith,jth,kth,lth

   integer :: jtot,jmin,jmax   ! J-values
   integer :: indx             ! constructed index for coupled matrix elements

   integer :: pair1,pair2      ! pairs needed to construct indx
   integer :: ia,ib,ic,id      ! orbital labels

   real cleb
   real hfactor  ! FOR HERMICIITY

   type (vjs), pointer :: tbme(:)
   integer pcref, pcstart, ppar
   integer,pointer :: xxmap(:)
   integer t
   character(1) :: hchar
   integer jk,jl,ji,jj,m,mk,ml,mi,mj,phasekl,phaseij
   integer pair11, pair22
   integer itmp
   integer par
   !! real zeta

   print*,' WARNING I DO NOT THINK ROUTINE fetchXXme WORKS YET '
   print*,' PROBABLY A PROBLEM WITH THE PHASE '
   if(it == 1)then
        tbme => ppme
        xxmap => ppcouplemap

   else
        tbme =>   nnme
        xxmap => nncouplemap

   endif
   T = 1        ! for pp, nn matrix elements
!.........  INTERPRET THE S.P. STATE INFO............ 
  if(l < 0)then
            lsps = -l
            lth = -it
  else
            lsps = l 
            lth = it
  endif

  if(k < 0)then
            ksps = -k
            kth = -it
  else
            ksps = k 
            kth = it
  endif

  if(j < 0)then
            jsps = -j
            jth = -it
  else
            jsps = j
            jth = it
  endif

  if(i < 0)then
            isps = -i
            ith = -it
  else
            isps = i
            ith = it
  endif

  ia = hspsqn(kth,ksps)%orb
  ib = hspsqn(lth,lsps)%orb
  ic = hspsqn(ith,isps)%orb
  id = hspsqn(jth,jsps)%orb

  ji = hspsqn(ith,isps)%j
  mi = hspsqn(ith,isps)%m
  jj = hspsqn(jth,jsps)%j
  mj = hspsqn(jth,jsps)%m
  jk = hspsqn(kth,ksps)%j
  mk = hspsqn(kth,ksps)%m
  jl = hspsqn(lth,lsps)%j
  ml = hspsqn(lth,lsps)%m
  m = mi+mj
  if(m /= mk+ml)then
     print*,' MISMATCH IN MMMMssss '
     stop
  end if
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
!-------------- MUST HAVE ic > id, else swap and pick up a phase ---
  phaseij = 1
  if(ic < id)then
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

!          if(hchar == 'J'.and. pair1/=pair2)cycle
!          if(hchar == 'T'.and. pair1/=pair2)cycle

!-------------- ALSO MUST HAVE STANDARD ORDERING OF PAIRS -----------

    pair1 = xxmap(pair1)
    pair2 = xxmap(pair2)

    if(pair1==-1 .or. pair2 == -1)then
             print*,' hmm strange ... (B)',it
             print*,pair1,pair2 
             print*,ia,ib,ic,id
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

    indx = (pair11-pcref)*(pair11-1-pcref)/2+pair22-pcref+pcstart
    if(indx <= 0)then
             print*,' problem indx ',indx
             stop
    endif

!-------------FIND JMIN,JMAX ALLOWED ----------------

   jmax = bmin( jk+jl,ji+jj) /2
   jmin = bmax( abs(jk-jl),abs(ji-jj))/2
   jmin = bmax(jmin,abs(m/2))

   hfactor = 1.0
   if((isps==ksps .and. jsps == lsps ).or. (isps==lsps.and.jsps==ksps))then
              hfactor = 0.5
    endif
!-------------NOW EXTRACT FROM TBMEs ---------------
!-------      NOTE FACTOR ZETA = SQRT(1+DELTA(A,B)) INCLUDED
   xxme = 0.
   do jtot = jmin,jmax
             xxme = xxme+tbme(indx)%v(jtot) *zeta(ia,ib)*zeta(ic,id)   * &
                    cleb(jk,mk,jl,ml,2*jtot,m)*cleb(ji,mi,jj,mj,2*jtot,m)
   enddo           
   xxme  = xxme*phasekl*phaseij*hfactor

   return
end subroutine fetchXXme

!=============================================================

end module btbme_mod
