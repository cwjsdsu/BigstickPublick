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
!  EDITED 2/17/09 by CWJ
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
subroutine prepare_to_uncouplePNtbme(countonly)

  use spstate
  use haiku_info
  use system_parameters
  use interaction
  use nodeinfo
  use verbosity
  use bmpi_mod
  use btbme_mod
  use io
  use TRstuff
  use menu_choices
  implicit none

  logical, intent(in) :: countonly  ! added in 7.9.2 to aid modeling
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
!  if(useTRphase)call findTRpairs_gen(0,npairpn,PN2%pair)  
  
!------------ MAP FROM ORIGINAL CODING TO SORTED 

  n = (nhsps(1) + nhsps(-1))*(nhsps(2)+nhsps(-2))
  if(.not.allocated(mappairpn))then
     allocate( mappairpn( n), stat=aerr)
     if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 2")
  endif
  mappairpn = 0
  
  do ipair = 1,npairpn
     indx = PN2%pair(ipair)%indx
     if(indx == 0)then
        if ( iproc == 0 ) then
           print*,' indx is 0 ',ipair
           print*,PN2%pair(ipair)%ia,PN2%pair(ipair)%ib
        end if
#ifdef _MPI	
        call BMPI_ABORT(MPI_COMM_WORLD,101,ierr)
#endif
        stop
     end if
     if ( mappairpn(indx) /=0 ) then
        if ( iproc == 0 ) then
           print*,' problem, pair already found PN '
           print*,ipair,indx
        end if
#ifdef _MPI	
        call BMPI_ABORT(MPI_COMM_WORLD,101,ierr)
#endif
        stop
     end if
     mappairpn(indx) = ipair
  end do
!-------------- SET UP START, STOP ARRAYS ------------

  jzmax = PN2%pair(npairpn)%m
  if(.not.countonly)then
     allocate( cpnpair( -nhsps(-2):nhsps(2), -nhsps(-1):nhsps(1)), stat=aerr)
     if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 3")
     allocate( dpnpair( -nhsps(-2):nhsps(2),-nhsps(-1):nhsps(1)), stat=aerr)
     if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 4")
     cpnpair = 0
     dpnpair = 0
	 
     if(useTRnew)then
         allocate( phasepnpair( -nhsps(-2):nhsps(2), -nhsps(-1):nhsps(1)), stat=aerr)
        if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 5")
		phasepnpair = 1
     end if
  end if
  
  call count_uncouplepn('H',.false.)  ! just count up

  if(countonly)return
  
  if ( iproc == 0 .and. verbose_uncouple )then
	  write(6,'(i14," uncoupled pn matrix elements (",f10.3," Gb)")')nmatpn,nmatpn*4.e-9
	  write(logfile,'(i14," uncoupled pn matrix elements (",f10.3," Gb)")')nmatpn,nmatpn*4.e-9	  
!	  print*,' There are ',nmatpn,' uncoupled pn matrix elements (',nmatpn*4.e-9,' Gb)'
  end if
  
  if(menu_char=='2b')then
	  allocate(dmatpn(nmatpn),dmatpnhc(nmatpn), stat=aerr)
  else
     allocate(hmatpn(nmatpn), stat=aerr)
  end if
  if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 6")
  
  
 ! if(countonly)call uncouplePNtbme_boss_new(countonly)
  
  
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
!    readPNdeformedtbmes  (only for non-spherical matrix elements, rarely used)
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
	  use menu_choices
	  use TRstuff

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

      integer(8) :: iref,istart
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
      integer(8) :: nsection
      integer(8) ipair,pairstart
      integer(8) cindx,dindx  ! creation/destruction indices
      integer gi,gj,gk,gl ! group indices
      integer n0,n0j,n0k,nsumk
      integer kmin,kmax

      integer countsection,isection
      integer, allocatable :: dpairstart(:), dpairstop(:)
	   integer(8), allocatable :: nmesum(:)

      integer :: aerr

      n0 = 0
      n0j = 0
      n0k = 0
      nsumk = 0
      if(Np(1)*Np(2) < 1)return
      if(uncouple)then
         select case (hchar)
           case('H','J')
              vmatpn => hmatpn
	          vmatpn=0.0
			  
           case('T')
              vmatpn => obsmatpn
	          vmatpn=0.0
			  
		   case('2')
			  dmatpn = 0.d0

         end select
		 
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
			
!		  if(useTRphase .and. .not.dens2bflag .and. m < 0 .and. dpair > 1)cycle
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
!         if(useTRphase .and. .not.dens2bflag .and. PN2%pair(dpair)%m < 0)cycle

        if(par /= PN2%pair(dpair)%par .or. m /= PN2%pair(dpair)%m)then
			
			
          isection = isection +1

!------------- ADD TO TOTAL THOSE FROM LAST TIME
          nme = nme + (nsection)*(nsection)
		  
		  if(nme < 0)then  ! error trap added in 7.10.6
			  print*,' problem in counting up PN matrix elements '
			  print*,nme,nsection,nsection*nsection
			  stop
			  
		  end if
		  
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
!                   END OF SET UP FOR OPENMP PARALLELISM  -------------	  
!----------------- OKAY, NOW READY TO GO!!!  --------------------
! OPENMP TEMP TURNED OFF WHILE I DEBUG TR OPTION (I think the problem is elsewhere, but this is to make sure)

!!$omp parallel do private(nme,dindx,cindx,first,nsection,pairstart,dpair,cpair, & 
!!$omp          m, par,i,j,k,l,isps,jsps,ksps,lsps,ith,jth,kth,lth, & 
!!$omp           ia,ib,ic,id,ji,jj,jk,jl,mi,mj,mk,ml,gi,gj,gk,gl, & 
!!$omp          phaseij,phasekl,ktr,ltr,itbme,itbmetr, pair1,pair11,pair2,pair22, &
!!$omp          indx,jmin,jmax,jtot, vtmp,ppar,pcref,pcstart), & 
!!$omp  shared(cpnpair,dpnpair,phasepnpair,vmatpn)

      do isection = 1,countsection   ! this is a dividing up for OpenMP parallelism
      nme = nmesum(isection)
      dindx = 1
      first = .true.
      nsection = dpairstop(isection)-dpairstart(isection)+1
      pairstart = dpairstart(isection)
      do dpair = dpairstart(isection),dpairstop(isection) ! loop over destruction pairs
          m = PN2%pair(dpair)%m
          par = PN2%pair(dpair)%par
		  if(useTRnew .and. m < 0 .and. .not. dens2bflag)then
			  print*,' I do not think I should be here with m < 0'
			  stop
		  end if		  

        if(.not.uncouple)then  ! ALREADY COUNTED, CAN SKIP THE REST
          cycle
        endif
!---------------- FROM NOW ON, UNCOUPLING -------------------------

!-------------- EXTRACT QUANTUM NUMBERS OF DESTRUCTION PAIR--------
        k = PN2%pair(dpair)%ia   ! extracting the single-particle STATES from an uncoupled PAIR
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

        pair1 = numorb(2)*(ia-1)+ib   ! find the COUPLED pair from single-particle ORBITALS
        pair1 = PNcouplemap(pair1)
        if(pair1 == -1)then
            print*,' uh problem pn map boss 1 '
            stop
        endif

!----------- SET UP MAPPING ARRAYS----------------------
        dpnpair(l,k) = dindx + nme   ! provides an index from proton, neutron single-particle STATES
!------------- TIME REVERSE ----------- OLD VERSION---
!       THIS CONNECTS single-particle states with a given m to -m
!       IDEA: let l, k be labels for single-particle states
!       then let ltr, ktr be labels for their 'time reverse' (m -> -m)
!       dnpnpair(l,k) = dnpairpair(ltr,ktr) but also pick up a phase
!
!       CAN ONLY BE  USED FOR HAMILTONIAN OPERATION, NOT FOR TWO-BODY DENSITIES

        if(useTRnew .and. m== 0)then
          phasepnpair(l,k) = 1
        endif
        if(useTRnew .and. m/= 0 .and. .not. dens2bflag)then
           phasepnpair(l,k) = 1
           ktr = hspsqn(kth,ksps)%tr
           if(kth >0 .and. mk /= 0)ktr = -ktr
           ltr = hspsqn(lth,lsps)%tr
           if(lth >0 .and. ml /= 0)ltr = -ltr
           phasepnpair(ltr,ktr) = int((-1)**( (jl+jk)/2),1)
           dpnpair(ltr,ktr) = dindx + nme   
        endif
!----------------- LOOP OVER CREATION PAIRS--------------
        cindx = 0
        do cpair = pairstart,pairstart+nsection-1
     
           cindx = cindx + 1
           itbme = dindx + (cindx-1)*nsection + nme
           itbmetr = cindx + (dindx-1)*nsection +nme ! "time-reversed"

           if(itbme > nmatpn .or. itbmetr > nmatpn)then
             print*,' me label too large (pn)'
             print*,itbme,itbmetr, nmatpn
             print*,cindx,dindx,nme,m,par
             print*,cpair,dpair,iref,istart
             stop
           endif
           if(itbme <= 0)then
               print*,' problem itbme (pn)',itbme
			   print*,dindx,cindx,nsection,nme
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
          if(useTRnew .and. m/= 0 .and. .not.dens2bflag)then   ! TR used for diagonalization, not for two-body densities
            itr = hspsqn(ith,isps)%tr
            if(ith >0 .and. mi/=0)itr = -itr
            jtr = hspsqn(jth,jsps)%tr
            if(jth >0 .and. mj /= 0)jtr = -jtr
            cpnpair(jtr,itr) = (cindx-1)*nsection
          endif

!---------------- BECAUSE I HAVE TIME REVERSAL, I CAN SKIP A LOT
          if(cpair > dpair)cycle
		  
!------------ IF TWO-BODY DENSITIES, DO NOT UNCOUPLE MATRIX ELEMENTS

          if(menu_char=='2b')cycle		  

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
		  
		  if(hchar/='2')then
             if(vmatpn(itbme)/=0.)then
               print*,' already here (pn) '
               print*,itbme,nmatpn
               stop
             endif
		   end if

          do jtot = jmin,jmax
                vtmp = vtmp+ pnme(indx)%v(jtot) & 
              * cleb(jk,mk,jl,ml,2*jtot,mk+ml)*cleb(ji,mi,jj,mj,2*jtot,mi+mj)
          enddo    ! jtot
          if(hchar/='2')then
             vmatpn(itbme)  = vtmp
             vmatpn(itbmetr)  = vtmp
		  end if
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
!!$omp end parallel do
!      print*,vmatpn(:)
      
      nme = nme + (nsection)*(nsection)  ! need to do for the last one
      if(iproc==0 .and. verbose_uncouple)print*,n0,' zeroes! '
      if(.not.uncouple)then
		  nmatpn = nme
		  if(iproc==0)print*,nmatpn,' uncoupled PN matrix elements '
	  end if

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
    call prepare_to_uncouplePNtbme(.false.)
    call clocker('mes','end')
	
	return
	
end subroutine delayedPNmatrixelements
   
!============================================================
!
!  NEW (7.9.4) APPROACH TO UNCOUPLING PN MATRIX ELEMENTS
!
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
subroutine uncouplePNtbme_boss_new(countonly)

  use spstate
  use haiku_info
  use system_parameters
  use interaction
  use nodeinfo
  use verbosity
  use bmpi_mod
  use btbme_mod
  use io
  use TRstuff
  use qs_on_proc
  use adv_tbme_info
  implicit none

  logical, intent(in) :: countonly  
  integer :: ierr

!  integer :: it
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
  
  
  if(.not.sort_XXpairs_on_W)then
	  print*,' Need to enable sorting on W to get full power '
	  print*,' set sort_XXpairs_on_W = .true. '
	  stop
  end if
!  print*,"(Y1)"
  
  call maxWpn(wmaxpn)
!  print*,"(Y2)"
  
!  print*,' SURVEY 2'
!  call survey_qs_X(1)
!  call survey_qs_X(2)
!  print*,' DONE SURVEYING '
!print*,"(Y3)"
  
  call count_create_rhos(1,.false.,nrhophX(1),rhophX(1)%pair)
 ! call count_create_pnpairs(wmaxpn,.false.,npairpn,PN2%pair)
  allocate( rhophX(1)%pair( nrhophX(1) ), stat=aerr)
  if(aerr /= 0) call memerror("uncouplePNtbme_boss_new 1")
  call count_create_rhos(1,.true.,nrhophX(1),rhophX(1)%pair)
  
!  print*,nrhophX(1),rhophX(1)%pair(:)%m
    
!  call count_create_pnpairs(wmaxpn,.true.,npairpn,PN2%pair)
  
!--------- NEXT, SORT THESE -------
  if(nrhophX(1)==0)then
     nmatpn = 0
     return
  endif
  call sortpair(nrhophX(1),rhophX(1)%pair)
!  print*,"(Y5)",sort_XXpairs_on_W
  
  if(sort_XXpairs_on_W)call sortXXpairsW(nrhophX(1),rhophX(1)%pair,.true.)
!  print*,"(Y6)"
  
  if(useTRphase)call findTRpairs_gen(1,nrhophX(1),rhophX(1)%pair)
!  print*,"(Y7)"
  
  call count_create_rhos(2,.false.,nrhophX(2),rhophX(2)%pair)
 ! call count_create_pnpairs(wmaxpn,.false.,npairpn,PN2%pair)
  allocate( rhophX(2)%pair( nrhophX(2) ), stat=aerr)
  if(aerr /= 0) call memerror("uncouplePNtbme_boss_new 2")
  call count_create_rhos(2,.true.,nrhophX(2),rhophX(2)%pair)
  call sortpair(nrhophX(2),rhophX(2)%pair)
  if(sort_XXpairs_on_W)call sortXXpairsW(nrhophX(2),rhophX(2)%pair,.true.)

  
  if(useTRphase)call findTRpairs_gen(2,nrhophX(2),rhophX(2)%pair)  
  
!  call survey_XY_qs(0,nprocs-1)
  
!............ ALL THE REST NEEDS TO BE REWRITTEN....
  
!------------ MAP FROM ORIGINAL CODING TO SORTED 

!  n = (nhsps(1) + nhsps(-1))*(nhsps(2)+nhsps(-2))
!  if(.not.allocated(mappairpn))then
!     allocate( mappairpn( n), stat=aerr)
!     if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 2")
!  endif
!  mappairpn = 0
  
!  do ipair = 1,npairpn
!     indx = PN2%pair(ipair)%indx
!     if(indx == 0)then
!        if ( iproc == 0 ) then
!           print*,' indx is 0 ',ipair
!           print*,PN2%pair(ipair)%ia,PN2%pair(ipair)%ib
!        end if
!        call BMPI_ABORT(icomm,101,ierr)
!        stop
!     end if
!     if ( mappairpn(indx) /=0 ) then
!        if ( iproc == 0 ) then
!           print*,' problem, pair already found PN '
!           print*,ipair,indx
!        end if
!        call BMPI_ABORT(icomm,101,ierr)
!        stop
!     end if
!     mappairpn(indx) = ipair
!  end do
!-------------- SET UP START, STOP ARRAYS ------------

!  jzmax = PN2%pair(npairpn)%m
!  if(.not.countonly)then
!     allocate( cpnpair( -nhsps(-2):nhsps(2), -nhsps(-1):nhsps(1)), stat=aerr)
!     if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 3")
!     allocate( dpnpair( -nhsps(-2):nhsps(2),-nhsps(-1):nhsps(1)), stat=aerr)
!     if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 4")
  
!     cpnpair = 0
!     dpnpair = 0
!  end if

call start_stop_count_PNme(wmaxpn)
  
!  call count_uncouplepn('H',.false.)  ! just count up
  
  if(countonly)return
  
  if ( iproc == 0 .and. verbose_uncouple )then
	  write(6,'(i14," uncoupled pn matrix elements (",f10.3," Gb)")')nmatpn,nmatpn*4.e-9
	  write(logfile,'(i14," uncoupled pn matrix elements (",f10.3," Gb)")')nmatpn,nmatpn*4.e-9	  
!	  print*,' There are ',nmatpn,' uncoupled pn matrix elements (',nmatpn*4.e-9,' Gb)'
  end if
  
!  allocate(hmatpn(nmatpn), stat=aerr)
!  if(aerr /= 0) call memerror("prepare_to_uncouplePNtbme 6")
    
  
  return
end subroutine uncouplePNtbme_boss_new

!========================================================

!
!  set up arrays for advanced storage (which includes blocks on W)
!  to count up PN atrix elements in new storage scheme
!

subroutine start_stop_count_PNme(wmaxpn)
	use TRstuff
	use qs_on_proc
	use haiku_info,only:nhsps
	use adv_tbme_info
	use spstate,only:hspsqn
	use interaction
!	use system_parameters
!	use sporbits,only:parmult
	implicit none
	integer,intent(in) :: wmaxpn
	
	integer :: minmp,maxmp,minparp,maxparp,mindWp,maxdWp,minWb,maxWb,minWd,maxWd
	integer :: iq 
	integer :: b,d
	integer :: wb
	integer :: mypJz,myppar,mydWp,myWb
	integer :: mynJz,mynpar,mydWn,myWd
	
	integer(8) :: nmexy
	integer :: rhophp,rhophn
	integer :: pairblocksize_p,pairstart_p,pairstop_p
	integer :: pairblocksize_n,pairstart_n,pairstop_n
	logical :: notfound
	
	integer :: summy
	
	
	
	
	minmp = 1000
	maxmp = -1000
	mindwp = 10000
	maxdwp = -1000 
	minparp = 2
	maxparp = 0
	
	do iq = 1,dqs_p%ndqs
		minmp = min(minmp,dqs_p%dJz(iq))
		maxmp = min(maxmp,dqs_p%dJz(iq))
		minparp = min(minparp,dqs_p%dpar(iq))
		maxparp = min(maxparp,dqs_p%dpar(iq))
		mindWp = min(minmp,dqs_p%dWp(iq))
		maxdWp= min(maxmp,dqs_p%dWp(iq))
		
		
	end do !iq
	
!............... NOW FIND MIN, MAX on W for protons	
	minWb = 1000
	maxWb = -1000
	minWd = 1000
	maxWd = -1000
	
	do b = 1,nhsps(-1)
		minWb = min(minWb,hspsqn(-1,b)%w)
		minWb = min(minWb,hspsqn(-1,b)%w)
		minWb = min(minWb,hspsqn(-1,b)%w)
		maxWb = max(maxWb,hspsqn(-1,b)%w)
	end do
	do b = 1,nhsps(1)
		minWb = min(minWb,hspsqn(1,b)%w)
		minWb = min(minWb,hspsqn(1,b)%w)
		minWb = min(minWb,hspsqn(1,b)%w)
		maxWb = max(maxWb,hspsqn(1,b)%w)
	end do
	do d = 1,nhsps(-2)
		minWd = min(minWb,hspsqn(-2,d)%w)
		minWd = min(minWb,hspsqn(-2,d)%w)
		minWd = min(minWb,hspsqn(-2,d)%w)
		maxWd = max(maxWb,hspsqn(-2,d)%w)
	end do
	do d = 1,nhsps(2)
		minWd = min(minWd,hspsqn(2,d)%w)
		minWd = min(minWd,hspsqn(2,d)%w)
		minWd = min(minWd,hspsqn(2,d)%w)
		maxWd = max(maxWd,hspsqn(2,d)%w)
	end do
	
	
	allocate(PNmestart(minmp:maxmp,minparp:maxparp, mindWp:maxdWp,minWb:maxWb))
	allocate(PNmeref(minmp:maxmp,minparp:maxparp, mindWp:maxdWp,minWb:maxWb))
	allocate(PNmeconj(minmp:maxmp,minparp:maxparp, mindWp:maxdWp,minWb:maxWb))
		
	PNmestart = -999
	PNmeref = -999
	PNmeconj = -999
	nmexy = 0
	rhophp = 1
	
!........... LOOK TO LIST OF ALLOWED dqs IN PN BUNDLES.................

    nmexy = 0
!	print*,' INTERMEDIATE ',nrhophX(1),nrhophX(2)
	summy = 0
    do iq = 1,XY_dqs%ndqs
		
		mypJz = XY_dqs%dJz(iq)
		myppar = XY_dqs%dpar(iq)
		mydWp  = XY_dqs%dWp(iq)
		mydWn  = XY_dqs%dWn(iq)
		
!		print*,mypJz,' JZ in section ',myppar,mydWp,mydWn

!		do mywb = minWb,maxWb
        maxWb=wmaxpn-minwd
		do wb = minWb,maxWb	 
			 
			call find_pair_block_size(nrhophX(1),rhophX(1)%pair,mypJz,myppar,.true.,mydWp, & 
			   .true.,wb,.false.,pairblocksize_p,pairstart_p,pairstop_p,notfound)
!			   print*,' found? ',notfound		
			if(notfound)cycle   
!            print*,mypJz,myppar,mydWp,mydWn,pairblocksize_p		
			if(mydWn==0)summy=summy+pairblocksize_p	
			myNjz = - mypJz
			mynpar = myppar
			
			maxwd = wmaxpn-wb
			if(maxwd < minwd)cycle
			
			call find_pair_block_size(nrhophX(2),rhophX(2)%pair,mynJz,mynpar,.true.,mydWn, & 
				   .false.,maxwd,.true.,pairblocksize_n,pairstart_n,pairstop_n,notfound)			
!	               print*,mypJz,myppar,mydWp,mydWn,pairblocksize_p,notfound		

			if(notfound)cycle   
!................ THE BLOCK OF STORED MATRIX ELEMENTS IS pairblocksize_p  X pairblocksize_n	
!print*,' intermediate ', 	pairblocksize_p,pairblocksize_n		,mypJz,myppar,mydWp
			nmexy = nmexy + int(pairblocksize_p,8)*int(pairblocksize_n,8)				 	
							
			
		end do
		
		
	end do	
	print*,' sum of pairs ', summy
	print*,' alternate count of PN matrix elements in new storage = ',nmexy
	print*,' or using ',nmexy*4.e-6,' Mb '

	return
end subroutine start_stop_count_PNme

!========================================================
