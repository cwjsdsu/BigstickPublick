!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK
!
!  file BTBMEPN.F90
!
!  uncouples pn part of the interaction
!
!==============================================================
!
!  ** prepare_to_uncouplePNtbme -- set ups arrays
!  CALLS:
!	maxWpn
!	count_create_pnpair
!	sortpair
!  	count_uncouplepn!
!
!  ** count_create_pnpair  -- creates proton-neutron pairs
!
!  ** count_uncouplepn  -- does actual decoupling
!
!  ** maxwpn	-- computes max W a pn pair can have
!
!

!==========================================================
!  LAST EDITED 2/17/09 by CWJ
!=============================================================!
!
!   subroutine prepare_to_uncouplePNtbme
!  routine to set up arrays for uncoupled pn tbmes
!
!  CALLED BY:
!   hoopmaster
!   onebodysetup
!
!  CALLS:
!	maxWpn
!	count_create_pnpair
!	sortpair
!  	count_uncouplepn
!
subroutine prepare_to_uncouplePNtbme

  use spstate
  use haiku_info
  use system_parameters
  use interaction
  use nodeinfo
  use verbosity
  use bmpi_mod
  use btbme_mod
  use io
  implicit none

  integer :: ierr

  integer :: it
  integer :: isps,jsps
  integer :: ipair,jpair
  integer :: indx
  integer :: jzmax
  integer :: m,icount
  integer :: n
  integer :: i,j
  integer :: ith,jth
  integer :: isgn,jsgn
  integer :: wmaxpn
  integer :: aerr
  
  if(Np(1)*Np(2) < 1)return
  call maxWpn(wmaxpn)
  call count_create_pnpairs(wmaxpn,.false.,npairpn,PN2%pair)
  allocate( PN2%pair( npairpn ), stat=aerr)
  if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 1")
  call count_create_pnpairs(wmaxpn,.true.,npairpn,PN2%pair)
  
!--------- NEXT, SORT THESE -------
  if(npairpn==0)then
     nmatpn = 0
     return
  endif
  call sortpair(npairpn,PN2%pair)

!------------ MAP FROM ORIGINAL CODING TO SORTED 

  n = (nhsps(1) + nhsps(-1))*(nhsps(2)+nhsps(-2))
  allocate( mappairpn( n), stat=aerr)
  if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 1")
  mappairpn = 0
  
  do ipair = 1,npairpn
     indx = PN2%pair(ipair)%indx
     if(indx == 0)then
        if ( iproc == 0 ) then
           print*,' indx is 0 ',ipair
           print*,PN2%pair(ipair)%ia,PN2%pair(ipair)%ib
        end if
        call BMPI_ABORT(icomm,101,ierr)
        stop
     end if
     if ( mappairpn(indx) /=0 ) then
        if ( iproc == 0 ) then
           print*,' problem, pair already found '
           print*,ipair,indx
        end if
        call BMPI_ABORT(icomm,101,ierr)
        stop
     end if
     mappairpn(indx) = ipair
  end do
!-------------- SET UP START, STOP ARRAYS ------------

  jzmax = PN2%pair(npairpn)%m

  allocate( cpnpair( -nhsps(-2):nhsps(2), -nhsps(-1):nhsps(1)), stat=aerr)
  if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 3")
  allocate( dpnpair( -nhsps(-2):nhsps(2),-nhsps(-1):nhsps(1)), stat=aerr)
  if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 4")


  if(useTR)then
      allocate( phasepnpair( -nhsps(-2):nhsps(2), -nhsps(-1):nhsps(1)), stat=aerr)
     if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 5")
  end if
  cpnpair = 0
  dpnpair = 0
  call count_uncouplepn('H',.false.)  ! just count up
  
  if ( iproc == 0 .and. verbose_uncouple )then
	  write(6,'(i14," uncoupled pn matrix elements (",f10.3," Gb)")')nmatpn,nmatpn*4.e-9
	  write(logfile,'(i14," uncoupled pn matrix elements (",f10.3," Gb)")')nmatpn,nmatpn*4.e-9	  
!	  print*,' There are ',nmatpn,' uncoupled pn matrix elements (',nmatpn*4.e-9,' Gb)'
  end if
  allocate(hmatpn(nmatpn), stat=aerr)
  if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 6")
  
  return
end subroutine prepare_to_uncouplePNtbme

!========================================================
!      subroutine count_uncouplepn
!  NEW version written 7/2013 by CWJ @ GSI 
!      to allow for OpenMP parallelization
!
!  INPUT:
!   hchar = default is "H" for hamiltonian
!     if = 'T' or 'J' then it simplifies
!   uncouple = flag to uncouple (true) or just count (false)
!
!  SUBROUTINES CALLED:
!    readPNdeformedtbmes  (only for non-spherical matrix elements)
!
      subroutine count_uncouplepn(hchar,uncouple)

      use verbosity
      use spstate
      use sporbit
      use haiku_info
      use system_parameters
	   use coupledmatrixelements
      use interaction
	   use interaction_deformed
      use nodeinfo
      use ntuple_info
      use butil_mod
	  use onebodypot, only: meanie

      implicit none
      character(1) hchar  ! signals what is operator
                          ! default is ordinary hamiltonian
                          ! 'J' = ang momentum
                          ! 'T' = isospin
      logical uncouple
      integer(8) :: nme

      integer dpair,cpair
      integer Jtot,jmin,jmax
      integer T
      integer M
      integer par
      integer ik,il
      integer ii,ij

      integer i,j,k,l
      integer ith,jth,kth,lth
      integer isps,jsps,ksps,lsps
      integer itr,jtr,ktr,ltr

      integer nk,nl,ni,nj
      integer jk,jl,ji,jj
      integer mk,ml,mi,mj

      integer iref,istart
      integer(8) :: itbme,itbmetr
      integer a,b,c,d
      integer ia,ib,ic,id
      integer pair1,pair2
      integer pair11,pair22
      integer ppar,pcref,pcstart
      integer phasekl,phaseij

      integer indx
      integer itmp
      real vtmp

!      real zeta   ! zeta(i,j) = sqrt(1+delta(i,j))
      real cleb
      logical first
      integer nsection
      integer ipair,pairstart
      integer cindx,dindx  ! creation/destruction indices
      integer gi,gj,gk,gl ! group indices
      integer n0,n0j,n0k,nsumk
      integer kmin,kmax

      integer countsection,isection
      integer, allocatable :: dpairstart(:), dpairstop(:)
	   integer(8), allocatable :: nmesum(:)

      integer :: aerr

      if(meanie)return  ! doing self-consistent mean-field
      n0 = 0
      n0j = 0
      n0k = 0
      nsumk = 0
      if(Np(1)*Np(2) < 1)return
      if(uncouple)then
         select case (hchar)
           case('H','J')
              vmatpn => hmatpn
           case('T')
              vmatpn => obsmatpn

         end select
         vmatpn=0.0
      endif

      nme = 0

      m = -9999
      par = 0
      nsection = 0
!---------------- THE FOLLOWING SETS UP ARRAYS TO MAKE POSSIBLE OpenMP PARALLELIZATION--
!-------- COUNT HOW MANY SECTIONS OF PAIRS by m, parity 
      isection = 0
      do dpair = 1,npairpn       ! loop over destruction pairs
!---------- FIND OUT IF THERE IS A CHANGE IN M OR PARITY
        if(par /= PN2%pair(dpair)%par .or. m /= PN2%pair(dpair)%m)then
          isection = isection +1
          m = PN2%pair(dpair)%m
          par = PN2%pair(dpair)%par
        end if
      end do
!---------- NOW ALLOCATE---------
      countsection = isection
      if(.not. allocated( dpairstart) )then
          allocate(dpairstart(countsection),dpairstop(countsection), stat=aerr )
          if(aerr /= 0) call memerror("count_uncouplepn 1")
          allocate(nmesum(countsection), stat=aerr )
          if(aerr /= 0) call memerror("count_uncouplepn 2")
      end if

      nme = 0

      m = -9999
      par = 0
      nsection = 0        
!-------- COUNT HOW MANY SECTIONS OF PAIRS by m, parity 
      isection = 0
      do dpair = 1,npairpn       ! loop over destruction pairs
!---------- FIND OUT IF THERE IS A CHANGE IN M OR PARITY
        if(par /= PN2%pair(dpair)%par .or. m /= PN2%pair(dpair)%m)then
          isection = isection +1

!------------- ADD TO TOTAL THOSE FROM LAST TIME
          nme = nme + (nsection)*(nsection)
          nmesum(isection) = nme

!---------IF SO, THEN COUNT UP HOW MANY HAVE THE SAME M, PARITY
          m = PN2%pair(dpair)%m
          par = PN2%pair(dpair)%par
          pairstart = dpair
          do ipair = pairstart,npairpn  ! loop upwards until M, PARITY change
                                      ! pairs have previously been sorted
             if(par == PN2%pair(ipair)%par .and. m == PN2%pair(ipair)%m)then
!------------ NSECTION IS HOW MANY PAIRS HAVE SAME M, PARITY
                nsection = ipair -pairstart+1
             else
                exit
             endif
          enddo
          dpairstart(isection) = dpair
          dpairstop(isection)  = dpair-1+nsection
          if(isection > 1)then
             if(dpairstart(isection)/=dpairstop(isection-1)+1)then
                print*,' Mismatch in dpairstart/stop '
                print*,isection, dpairstop(isection-1),dpairstart(isection)
                stop
             end if
          end if
        endif

      end do
!----------------- OKAY, NOW READY TO GO!!!

!$omp parallel do private(nme,dindx,cindx,first,nsection,pairstart,dpair,cpair, & 
!$omp          m, par,i,j,k,l,isps,jsps,ksps,lsps,ith,jth,kth,lth, & 
!$omp           ia,ib,ic,id,ji,jj,jk,jl,mi,mj,mk,ml,gi,gj,gk,gl, & 
!$omp          phaseij,phasekl,ktr,ltr,itbme,itbmetr, pair1,pair11,pair2,pair22, &
!$omp          indx,jmin,jmax,jtot, vtmp,ppar,pcref,pcstart), & 
!$omp  shared(cpnpair,dpnpair,phasepnpair,vmatpn)

      do isection = 1,countsection 
      nme = nmesum(isection)
      dindx = 1
      first = .true.
      nsection = dpairstop(isection)-dpairstart(isection)+1
      pairstart = dpairstart(isection)
      do dpair = dpairstart(isection),dpairstop(isection) ! loop over destruction pairs
          m = PN2%pair(dpair)%m
          par = PN2%pair(dpair)%par

        if(.not.uncouple)then  ! ALREADY COUNTED, CAN SKIP THE REST
          cycle
        endif
!---------------- FROM NOW ON, UNCOUPLING -------------------------

!-------------- EXTRACT QUANTUM NUMBERS OF DESTRUCTION PAIR--------
        k = PN2%pair(dpair)%ia
        l = PN2%pair(dpair)%ib
        if( k < 0)then
          ksps = -k
          kth = -1
        else
          ksps = k
          kth  = 1
        endif
        ia = hspsqn(kth,ksps)%orb
        jk = hspsqn(kth,ksps)%j
        mk = hspsqn(kth,ksps)%m
        gk = hspsqn(kth,ksps)%group

        if( l < 0)then
          lsps = -l
          lth = -2
        else
          lsps = l
          lth  = 2
        endif
        ib = hspsqn(lth,lsps)%orb
        jl = hspsqn(lth,lsps)%j
        ml = hspsqn(lth,lsps)%m
        gl = hspsqn(lth,lsps)%group
        phasekl = 1

!-------------- FIND INDEX "PAIR1" used to find coupled tbme

        pair1 = numorb(2)*(ia-1)+ib
        pair1 = PNcouplemap(pair1)
        if(pair1 == -1)then
            print*,' uh problem pn map boss 1 '
            stop
        endif

!----------- SET UP MAPPING ARRAYS----------------------
        dpnpair(l,k) = dindx + nme
!------------- TIME REVERSE -----------
        if(useTR .and. m== 0)then
          phasepnpair(l,k) = 1
        endif
        if(useTR .and. m/= 0)then
          phasepnpair(l,k) = 1
          ktr = hspsqn(kth,ksps)%tr
          if(kth >0 .and. mk /= 0)ktr = -ktr
          ltr = hspsqn(lth,lsps)%tr
          if(lth >0 .and. ml /= 0)ltr = -ltr
          phasepnpair(ltr,ktr) = int((-1)**( (jl+jk)/2),1)
          dpnpair(ltr,ktr) = dindx + nme
!          print*,l,k,dpnpair(l,k),ltr,ktr,dpnpair(ltr,ktr)
        endif
!----------------- LOOP OVER CREATION PAIRS--------------
        cindx = 0
        do cpair = pairstart,pairstart+nsection-1
     
           cindx = cindx + 1
           itbme = dindx + (cindx-1)*nsection + nme
           itbmetr = cindx + (dindx-1)*nsection +nme ! "time-reversed"

           if(itbme > nmatpn .or. itbmetr > nmatpn)then
             print*,' me label too large'
             print*,itbme,itbmetr, nmatpn
             print*,cindx,dindx,nme,m,par
             print*,cpair,dpair,iref,istart
             stop
           endif
           if(itbme <= 0)then
               print*,' problem itbme ',itbme
               stop
           endif
!-------------- EXTRACT QUANTUM NUMBERS OF CREATION PAIR--------

           i = PN2%pair(cpair)%ia
           j = PN2%pair(cpair)%ib

           if( i < 0)then
             isps = -i
             ith = -1
           else
             isps = i
             ith  = 1
           endif
           ic = hspsqn(ith,isps)%orb
           ji = hspsqn(ith,isps)%j
           mi = hspsqn(ith,isps)%m
           gi = hspsqn(ith,isps)%group
          if( j < 0)then
             jsps = -j
            jth = -2
          else
            jsps = j
            jth  = 2
          endif
          id = hspsqn(jth,jsps)%orb
          jj = hspsqn(jth,jsps)%j
          mj = hspsqn(jth,jsps)%m
          gj = hspsqn(jth,jsps)%group

!------------------ SET MAPPING ARRAY
          cpnpair(j,i) = (cindx-1)*nsection
          if(useTR .and. m/= 0)then
            itr = hspsqn(ith,isps)%tr
            if(ith >0 .and. mi/=0)itr = -itr
            jtr = hspsqn(jth,jsps)%tr
            if(jth >0 .and. mj /= 0)jtr = -jtr
            cpnpair(jtr,itr) = (cindx-1)*nsection
          endif

!---------------- BECAUSE I HAVE TIME REVERSAL, I CAN SKIP A LOT
          if(cpair > dpair)cycle

!------------- PUT INTO "STANDARD" ORDER; MAY GET PHASE -------

          phaseij = 1

          pair2 = numorb(2)*(ic-1)+id
          pair2 = PNcouplemap(pair2)
          if(pair2==-1)then
           print*,' oh boy should not have gotten that '
           stop
          endif

!-------------- ALSO MUST HAVE STANDARD ORDERING OF PAIRS -----------

          if(pair1 < pair2)then
            pair11 = pair2
            pair22 = pair1
          else
           pair11 = pair1
           pair22 = pair2
          endif

          ppar = PNcouples%pairc(pair1)%par
          pcref = PNcouples%meref(ppar)
          pcstart = PNcouples%mestart(ppar)

          indx = (pair11-pcref)*(pair11-pcref-1)/2+pair22-pcref +pcstart
          if(indx <= 0)then
             print*,' problem indx ',indx,pair11,pair22
             print*,ia,ib,ic,id
             print*,cpair,dpair
             stop
          endif

!-------------FIND JMIN,JMAX ALLOWED ----------------

          jmax = bmin( jk+jl,ji+jj) /2
          jmin = bmax( abs(jk-jl),abs(ji-jj))/2
          jmin = bmax(jmin,abs(m))

!-------------NOW EXTRACT FROM TBMEs ---------------
!-------      NOTE FACTOR ZETA = SQRT(1+DELTA(A,B)) INCLUDED
          vtmp = 0.
          if(vmatpn(itbme)/=0.)then
            print*,' already here (pn) '
            print*,itbme,nmatpn
            stop
          endif

          do jtot = jmin,jmax
                vtmp = vtmp+ pnme(indx)%v(jtot) & 
              * cleb(jk,mk,jl,ml,2*jtot,mk+ml)*cleb(ji,mi,jj,mj,2*jtot,mi+mj)
          enddo    ! jtot

          vmatpn(itbme)  = vtmp
          vmatpn(itbmetr)  = vtmp
          if(vtmp == 0.)then
             n0 = n0 +1
             if(itbme /= itbmetr)n0 = n0 +1
          endif
          if(print_matrix)then
          write(matrix_file,333)itbme,i,j,k,l,phasekl,phaseij,vmatpn(itbme)
333      format(7i5,f10.5)
          endif
        enddo  ! loop over cpair
          first = .false.
          dindx = dindx + 1
      enddo  ! loop over dpair
      end do ! loop over isection
!$omp end parallel do
!      print*,vmatpn(:)
      
      nme = nme + (nsection)*(nsection)  ! need to do for the last one
      if(iproc==0 .and. verbose_uncouple)print*,n0,' zeroes! '
      if(.not.uncouple)nmatpn = nme

      if(deformed .and. hchar=='H'.and. uncouple)call readPNdeformedtbmes
      return
      end subroutine count_uncouplepn

!===================================================

!  subroutine maxwpair
!  finds the maximum W allowed for a pn pairt
!  by comparing max W for N, N-2 particles
!  
   subroutine maxwpn(wpn)

   use spstate
   use sporbit
   use W_info
   use haiku_info
   use system_parameters

   implicit none
  
   integer wpn
   
   integer wp,wn
   integer i
   wpn = 0
   if(Np(1) < 1 .or. Np(2) < 1)return

   wp = 0
   wn = 0
   do i =1,Np(1)-1
        wp = wp + spsqn(1,i)%w
   enddo
   do i =1,Np(2)-1
        wn = wn + spsqn(2,i)%w
   enddo

   wpn = maxWtot - wp - wn
   return
   end subroutine maxwpn
   
!============================================================
!
!  added in 7.7.9
!

subroutine delayedPNmatrixelements
	
    use ntuple_info
	
	
!    call master_cross_ntuple(1,1)  MUST CALL EARLIER
    call clocker('mes','sta')
    call prepare_to_uncouplePNtbme
    call clocker('mes','end')
	
	return
	
end subroutine delayedPNmatrixelements
   
!============================================================
!
! ROUTINES IN DEVELOPMENT; NOT CURRENTLY USED:
!
!
!=============================================================
!
! NOT CURRENTLY USED
!
! subroutine extraPNpairsort
!      added August 2010 by CWJ @ SDSU
! additional sorting for PN pairs
! for fixed jz, parity, sorts into two groups:
! those with proton jz < 0, and those with jz >= 0
! this will ultimately enable using time-reversal symmetry in 1-body jumps
! 

subroutine extraPNpairsort(npair,pair)
  use pairdef

  implicit none
  integer npair
  type (pair_qn) :: pair(npair)
  type (pair_qn) :: pairtmp

  integer ipair, jpair
  integer m,par

  integer pairstart,pairstop
  integer nval 
  integer nme
  pairstart = 1

  nme = 0
  do while(pairstart < npair +1)
     
      par = pair(pairstart)%par
      m   = pair(pairstart)%m
!........... FIND END OF THIS SECTION WITH SAME PAR, M...............
      pairstop = pairstart 
      do while( pair(pairstop)%m == m .and. pair(pairstop)%par == par)
          pairstop = pairstop +1
          if(pairstop == npair+1)exit
      end do
      pairstop = pairstop -1 
!.................... NOW COUNT UP HOW MANY IN LOWER SECTION.........
!               NOTE: the lower section are all pairs with protons with m < 0
!                     upper section are all pairs with protons with m >= 0

      nval = 0
      do ipair = pairstart, pairstop
          if(pair(ipair)%ia < 0 )nval = nval + 1
      end do

!.................... GO THROUGH AND SWAP.............................

      if( nval > 0 .and. nval < pairstop - pairstart +1)then
      jpair = pairstop
      do ipair = pairstart, pairstart + nval -1
          if( pair(ipair)%ia > 0)then  ! need to swap
              if( pair(jpair)%ia < 0) then  ! swap this
                  pairtmp = pair(jpair)
                  pair(jpair) = pair(ipair)
                  pair(ipair) = pairtmp
              else
                  jpair = jpair - 1
                  if(jpair <= ipair)then
                     print*,' gone too far ',jpair,ipair
                     stop
                  endif
              endif

          end if
      end do
      endif
!.......... COUNT UP MATRIX ELEMENTS
      nme = nme + (pairstop - pairstart + 1 -nval )*( pairstop - pairstart + nval +2)/2
      
!............... NOW THIS SECTION IS SORTED; PREPARE TO GO ON TO NEXT SECTION..........
      pairstart = pairstop + 1
  end do
  print*, ' expect ',nme, ' pn matrix elements '

  return
end subroutine extraPNpairsort


