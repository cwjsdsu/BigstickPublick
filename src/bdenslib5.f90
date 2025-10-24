!
!  file BDENSLIB5.f90
!
! utilities for data structures for two-body densities
!

module twobodydens_util
	use interaction
	use coupledmatrixelements
	
	implicit none
	
	integer :: dens2bfile = 47
	
	integer :: npnsection
    integer, allocatable :: pnpairstart(:), pnpairstop(:)
    integer(8), allocatable :: npnmesum(:)
	
contains
	
! ROUTINE FOR COUPLING TOGETHER pp/nn TWO-BODY DENSITIES	
! here it is important to distinguish between ITBME, which labels uncoupled (m-scheme) matrix elements
! which are in turn stored in hmatpp/hmatnn
! and INDEX, which labels coupled (J-scheme) matrix elements
! stored in pp2bden and nn2bden
!	
	subroutine couple_2bdensXX(it,JJi,JJf,zerome)
		
      use verbosity
      use spstate
      use haiku_info
      use system_parameters
      use interaction
      use coupledmatrixelements
      use butil_mod

      implicit none
      integer, intent(in)::  it      ! which species
	  integer, intent(in) :: JJi,JJf  ! 2 x initial, final J of states
      logical, intent(in) ::zerome  ! needed to reset
      integer(kind=8) :: dpair,cpair
	  integer Jabmin,Jabmax,Jcdmin,Jcdmax,Jab,Jcd
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
	  integer :: ccouple,dcouple
	  integer :: couple1,couple2

      integer phasekl,phaseij

      integer indx,indxhc
      integer itmp
      real vtmp

      real cleb
	  real cleb4reduce
      real hfactor  ! Accounts for ordering for identical orbits
      real(8),pointer :: v(:),vhc(:)

	  type (vjs2b), pointer :: x2bden(:)
      integer pcref, pcstart, ppar
      integer,pointer :: xxmap(:)
      logical :: foundit
      integer :: aerr
	  integer(8) :: coupleblocksize,pairblocksize
	  
	  logical :: absame,cdsame   ! if pairs have identical orbits, this restricts Js
	  logical :: hcflag   ! using hermitian conjugates--added 7.10.4
	  
      if(Np(it) < 2)return

      if(it ==1)then

        v => dmatpp
		vhc => dmatpphc
      else

        v => dmatnn
		vhc => dmatnnhc

      endif

      if(it == 1)then
        x2bden => pp2bden
        xxmap => ppcouplemap

      else
        x2bden =>   nn2bden
        xxmap => nncouplemap

      endif

!----------------------------- ZERO OUT ARRAY-----------------

      T = 1        ! for pp, nn matrix elements


      do dpair = 1,npairXX(it)        ! loop over destruction m-scheme pair

        m = XX2(it)%pair(dpair)%m     ! find Jz, parity of pair
        par = XX2(it)%pair(dpair)%par

        iref = XX2(it)%meref(m,par)      ! iref, istart used to encode index for 
        istart = XX2(it)%mestart(m,par)  ! the uncoupled matrix element
		pairblocksize= XX2(it)%block(m,par)

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
        ic = hspsqn(kth,ksps)%orb
        jk = hspsqn(kth,ksps)%j
        mk = hspsqn(kth,ksps)%m

        if( l < 0)then
          lsps = -l
          lth = -it
        else
          lsps = l
          lth  = it
        endif
        id = hspsqn(lth,lsps)%orb
        jl = hspsqn(lth,lsps)%j
        ml = hspsqn(lth,lsps)%m
		
		cdsame =.false.
		if(ic == id)cdsame = .true.

        phasekl = 1
!-------------- MUST HAVE ia > ib, else swap and pick up a phase ---

        if(ic < id)then
             itmp = ic
             ic = id
             id = itmp

             itmp = jk
             jk = jl
             jl = itmp

             itmp = mk
             mk   = ml
             ml   = itmp

             phasekl = -1  ! the rest of the phase comes from the Clebsch Gordan

        endif
 

        do cpair = 1,dpair !npairXX(it)       ! loop over creation m-scheme pair  

          if(XX2(it)%pair(cpair)%m /=m)cycle       ! enforce quantum numbers
          if(XX2(it)%pair(cpair)%par /= par)cycle

          if(iref > cpair)then
            print*,' problem with XX iref (c) '
            print*,cpair,iref
            stop 
          endif

!--------------- HOW INDEX FOR UNCOUPLED (m-scheme) tbme IS COMPUTED---------
!          itbme = istart + pairblocksize*(dpair-iref-1)+cpair-iref
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

		  if(v(itbme)==0.0 .and. vhc(itbme) ==0.0.and. .not.zerome)cycle

!--------------- GET Q#s FOR CREATION PAIRS ---------------

          i = XX2(it)%pair(cpair)%ia
          j = XX2(it)%pair(cpair)%ib

!          print*,' operators ',i,j,k,l
          if( i < 0)then
            isps = -i
            ith = -it
          else
            isps = i
            ith  = it
          endif
          ia = hspsqn(ith,isps)%orb
          ji = hspsqn(ith,isps)%j
          mi = hspsqn(ith,isps)%m

          if( j < 0)then
            jsps = -j
            jth = -it
          else
            jsps = j
            jth  = it
          endif
          ib = hspsqn(jth,jsps)%orb
          jj = hspsqn(jth,jsps)%j
          mj = hspsqn(jth,jsps)%m
		  
!-------------- account for double counting
          hfactor = 1.0
          !if(cpair == dpair)then
          !   hfactor = 0.5
          !else
          !   hfactor = 1.0
          !endif	
		  
		  if(ia==ib)hfactor = hfactor*2. ! this is because   a^+_i  a^+j is ordered
		  if(ic==id)hfactor = hfactor*2.
	  
!------------- PUT INTO "STANDARD" ORDER; MAY GET PHASE -------

          phaseij = 1
		  
		  absame = .false.
		  if(ia==ib)absame =.true.

          if(ia < ib) then

             itmp = ia
             ia = ib
             ib = itmp

             itmp = ji
             ji = jj
             jj = itmp
             itmp = mi
             mi   = mj
             mj   = itmp

             phaseij = -1  ! the rest of the phase comes from the Clebsch Gordan

          endif
          dcouple = ic*(ic-1)/2 + id  ! used to find index of coupled TBME
          ccouple = ia*(ia-1)/2 + ib

!-------------- ALSO MUST HAVE STANDARD ORDERING OF PAIRS -----------

          ccouple = xxmap(ccouple)
          dcouple = xxmap(dcouple)
          if(ccouple==-1 .or. dcouple == -1)then
             print*,' hmm strange ... (A)',it
             print*,ccouple,dcouple 
             print*,ia,ib,ic,id
             print*,xxmap
             stop  
          endif
          ppar = XXcouples(it)%pairc(ccouple)%par
          if(ppar /= XXcouples(it)%pairc(dcouple)%par)then
              print*,' sigh parity problem '
              print*,ia,ib,ic,id,ppar,par
              print*,ccouple,dcouple
              stop  
          endif

          pcref = XXcouples(it)%meref(ppar)
          pcstart = XXcouples(it)%mestart(ppar)
		  
          if(ccouple <= XXcouples(it)%meref(2))then
			coupleblocksize = XXcouples(it)%meref(2)
		  else
			coupleblocksize = ncouplesXX(it)-XXcouples(it)%meref(2)
			
		  end if	
!.......... GET INDEX OF COUPLED (J-scheme) MATRIX ELEMENT.....		  

           indx = (dcouple - pcref -1 )*coupleblocksize +  		 & 
		          ccouple - pcref + pcstart			  		  
				  
          if(indx <= 0)then
             print*,' problem indx ',indx
             stop  
          endif

!-------------FIND JMIN,JMAX ALLOWED ----------------
          Jabmax = (ji+jj)/2
		  Jabmin =abs(ji-jj)/2
		  Jcdmax = (jk+jl)/2
		  Jcdmin = abs(jk-jl)/2
		  
!.... CHECK THESE AGAINST PRE-DETERMINED LIMITS

          if(Jabmax /= x2bden(indx)%Jabmax)then
			  print*,' mismatch in Jabmax ',Jabmax, x2bden(indx)%Jabmax,' indx ',indx 
			  print*,ia,ib,ic,id
			  print*,dcouple,pcref,coupleblocksize,ccouple,pcstart

			  stop
		  end if
		  
          if(Jabmin /= x2bden(indx)%Jabmin)then
			  print*,' mismatch in Jabmin ',Jabmin, x2bden(indx)%Jabmin
			  print*,i,j,ji,jj,indx
			  print*,coupleblocksize
			  stop
		  end if
          if(Jcdmax /= x2bden(indx)%Jcdmax)then
			  print*,' mismatch in Jcdmax ',Jcdmax, x2bden(indx)%Jcdmax
			  stop
		  end if
          if(Jcdmin /= x2bden(indx)%Jcdmin)then
			  print*,' mismatch in Jcdmin ',Jcdmin, x2bden(indx)%Jcdmin
			  stop
		  end if

!-------------NOW EXTRACT FROM TBMEs ---------------
!-------      NOTE FACTOR ZETA = SQRT(1+DELTA(A,B)) INCLUDED
          vtmp = v(itbme)*phasekl*phaseij *hfactor
!		  if(zerome)print*,i,j,k,l,itbme,v(itbme),vhc(itbme),cpair,dpair
		  if(vtmp==0.0 .and. vhc(itbme)==0.0 .and. .not.zerome)cycle

!		  if(.not.zerome)print*,' here ',mk,mi,v(itbme),vhc(itbme),phasekl*phaseij *hfactor
do Jab = Jabmin,Jabmax
	if(absame .and. (mod(Jab,2)/=0))cycle
	do Jcd = Jcdmin,Jcdmax
		if(cdsame .and. (mod(Jcd,2)/=0))cycle
		
        Jmin = x2bden(indx)%Jmin
		Jmin = max(Jmin,abs(Jab-Jcd))
        Jmax = x2bden(indx)%Jmax
		Jmax = min(Jmax,Jab+Jcd)
		if(Jmax < Jmin)cycle
		if(zerome)then
			x2bden(indx)%v(:,:,:) = 0.0
			cycle
		end if
          do jtot = jmin,jmax
			  
			  cleb4reduce = cleb(jtot*2,0,JJi,Jz,JJf,Jz)
			  if(abs(cleb4reduce) < 0.00001)then
				  x2bden(indx)%v(Jab,Jcd,Jtot)= -999.
				  cycle
			  end if
			  x2bden(indx)%v(Jab,Jcd,Jtot)=x2bden(indx)%v(Jab,Jcd,Jtot) + & 
			  vtmp * & 
                    cleb(jk,mk,jl,ml,2*jcd,2*m)*cleb(ji,mi,jj,mj,2*jab,2*m)*cleb(2*jab,2*m,2*jcd,-2*m,2*Jtot,0)  & 
					* (-1)**(Jcd -m)   & ! PHASE FROM time-reversal of second pair
					* sqrt(JJf+1.0)*(-1)**(Jtot+ (JJf-JJi)/2)/ cleb4reduce &    ! reduce
					/sqrt(2.0*Jtot+1.0)    ! needed for proper definition of density matrix
!					print*,' huh ', jab,jcd,jtot,&
!                    cleb(jk,mk,jl,ml,2*jcd,2*m)*cleb(ji,mi,jj,mj,2*jab,2*m)*cleb(2*jab,2*m,2*jcd,-2*m,2*Jtot,0)  & 
!					* (-1)**(Jcd -m)   & ! PHASE FROM time-reversal of second pair
!					* sqrt(JJf+1.0)*(-1)**(Jtot+ (JJf-JJi)/2)/ cleb4reduce &    ! reduce
!					/sqrt(2.0*Jtot+1.0) 
					
!					print*,jab,jcd,jtot, x2bden(indx)%v(Jab,Jcd,Jtot)
          enddo      

    end do  ! Jcd
end do   ! Jab
!if(.not.zerome)print*,indx,x2bden(indx)%v(2,0,2), x2bden(1)%v(2,0,2)
!................. NOW DO HERMITIAN CONJUGATE... ADDED 7.10.4.....

          if(cpair==dpair)cycle
          if(dcouple <= XXcouples(it)%meref(2))then
			coupleblocksize = XXcouples(it)%meref(2)
		  else
			coupleblocksize = ncouplesXX(it)-XXcouples(it)%meref(2)
			
		  end if	
!.......... GET INDEX OF COUPLED (J-scheme) MATRIX ELEMENT.....		  


           indx = (ccouple - pcref -1 )*coupleblocksize +  		 & 
		          dcouple - pcref + pcstart			  		  
				  
          if(indx <= 0)then
             print*,' problem indx ',indx
             stop  
          endif

!-------------FIND JMIN,JMAX ALLOWED ----------------
          Jabmax = (ji+jj)/2
		  Jabmin =abs(ji-jj)/2
		  Jcdmax = (jk+jl)/2
		  Jcdmin = abs(jk-jl)/2
		  
!.... CHECK THESE AGAINST PRE-DETERMINED LIMITS

          if(Jabmax /= x2bden(indx)%Jcdmax)then
			  print*,' mismatch in Jabmax ',Jabmax, x2bden(indx)%Jabmax,' indx ',indx 
			  print*,ia,ib,ic,id
			  print*,dcouple,pcref,coupleblocksize,ccouple,pcstart

			  stop
		  end if
		  
          if(Jabmin /= x2bden(indx)%Jcdmin)then
			  print*,' mismatch in Jabmin ',Jabmin, x2bden(indx)%Jabmin
			  print*,i,j,ji,jj,indx
			  print*,coupleblocksize
			  stop
		  end if
          if(Jcdmax /= x2bden(indx)%Jabmax)then
			  print*,' mismatch in Jcdmax ',Jcdmax, x2bden(indx)%Jcdmax
			  stop
		  end if
          if(Jcdmin /= x2bden(indx)%Jabmin)then
			  print*,' mismatch in Jcdmin ',Jcdmin, x2bden(indx)%Jcdmin
			  stop
		  end if

!-------------NOW EXTRACT FROM TBMEs ---------------
!-------      NOTE FACTOR ZETA = SQRT(1+DELTA(A,B)) INCLUDED
          vtmp = vhc(itbme)*phasekl*phaseij *hfactor
		  
		  if(vtmp==0.0 .and. .not.zerome)cycle

do Jab = Jabmin,Jabmax
	if(absame .and. (mod(Jab,2)/=0))cycle
	do Jcd = Jcdmin,Jcdmax
		if(cdsame .and. (mod(Jcd,2)/=0))cycle
		
        Jmin = x2bden(indx)%Jmin
		Jmin = max(Jmin,abs(Jab-Jcd))
        Jmax = x2bden(indx)%Jmax
		Jmax = min(Jmax,Jab+Jcd)
		if(Jmax < Jmin)cycle
		if(zerome)then
			x2bden(indx)%v(:,:,:) = 0.0
			cycle
		end if
          do jtot = jmin,jmax
			  
			  cleb4reduce = cleb(jtot*2,0,JJi,Jz,JJf,Jz)
			  if(abs(cleb4reduce) < 0.00001)then
				  x2bden(indx)%v(Jcd,Jab,Jtot)= -999.
				  cycle
			  end if
			  x2bden(indx)%v(Jcd,Jab,Jtot)=x2bden(indx)%v(Jcd,Jab,Jtot) + & 
			  vtmp * & 
                    cleb(jk,mk,jl,ml,2*jcd,2*m)*cleb(ji,mi,jj,mj,2*jab,2*m)*cleb(2*jcd,-2*m,2*jab,2*m,2*Jtot,0)  & 
					* (-1)**(Jab -m)   & ! PHASE FROM time-reversal of second pair
					* sqrt(JJf+1.0)*(-1)**(Jtot+ (JJf-JJi)/2)/ cleb4reduce &    ! reduce
					/sqrt(2.0*Jtot+1.0)    ! needed for proper definition of density matrix

          enddo      

    end do  ! Jcd
end do   ! Jab
!if(.not.zerome)print*,' again ',indx,x2bden(indx)%v(2,0,2), x2bden(1)%v(2,0,2)

        enddo  ! loop over cpair

      enddo  ! loop over dpair

      return		
		
		
	end subroutine couple_2bdensXX

!
!==================================================================================
subroutine setup_PNarrays
	
	use interaction
	implicit none
	integer :: m, par
	integer :: nme
	integer :: isection,dimsection
	
	integer :: xpair,ipair,xstart
	
	integer :: aerr
	
      nme = 0

      m = -9999
      par = 0
!---------------- THE FOLLOWING SETS UP ARRAYS TO MAKE POSSIBLE OpenMP PARALLELIZATION--
!-------- COUNT HOW MANY SECTIONS OF PAIRS by m, parity 
      isection = 0
      do xpair = 1,npairpn       ! loop over destruction pairs
!---------- FIND OUT IF THERE IS A CHANGE IN M OR PARITY
        if(par /= PN2%pair(xpair)%par .or. m /= PN2%pair(xpair)%m)then
          isection = isection +1
          m = PN2%pair(xpair)%m
          par = PN2%pair(xpair)%par
        end if
      end do
!---------- NOW ALLOCATE---------
      npnsection = isection
      if(.not. allocated( pnpairstart) )then
          allocate(pnpairstart(npnsection),pnpairstop(npnsection), stat=aerr )
          if(aerr /= 0) call memerror("count_uncouplepn 1")
          allocate(npnmesum(npnsection), stat=aerr )
          if(aerr /= 0) call memerror("count_uncouplepn 2")
      end if

      nme = 0

      m = -9999
      par = 0
!-------- COUNT HOW MANY SECTIONS OF PAIRS by m, parity 
      dimsection = 0
	  isection = 0
      do xpair = 1,npairpn       ! loop over destruction pairs
!---------- FIND OUT IF THERE IS A CHANGE IN M OR PARITY
        if(par /= PN2%pair(xpair)%par .or. m /= PN2%pair(xpair)%m)then
          isection = isection +1

!------------- ADD TO TOTAL THOSE FROM LAST TIME
          nme = nme + (dimsection)*(dimsection)
          npnmesum(isection) = nme

!---------IF SO, THEN COUNT UP HOW MANY HAVE THE SAME M, PARITY
          m = PN2%pair(xpair)%m
          par = PN2%pair(xpair)%par
          xstart = xpair
          do ipair = xstart,npairpn  ! loop upwards until M, PARITY change
                                      ! pairs have previously been sorted
             if(par == PN2%pair(ipair)%par .and. m == PN2%pair(ipair)%m)then
!------------ NSECTION IS HOW MANY PAIRS HAVE SAME M, PARITY
                dimsection = ipair -xstart+1
             else
                exit
             endif
          enddo
          pnpairstart(isection) = xpair
          pnpairstop(isection)  = xpair-1+dimsection
          if(isection > 1)then
             if(pnpairstart(isection)/=pnpairstop(isection-1)+1)then
                print*,' Mismatch in pnpairstart/stop '
                print*,isection, pnpairstop(isection-1),pnpairstart(isection)
                stop
             end if
          end if
        endif

      end do
	
	return
	
end subroutine setup_PNarrays

!==================================================================================
! ROUTINE FOR COUPLING TOGETHER pn TWO-BODY DENSITIES	
	
	subroutine couple_2bdensPN(JJi,JJf,zerome)
        use verbosity
        use spstate
        use haiku_info
        use system_parameters
        use interaction
        use coupledmatrixelements
        use butil_mod
		use sporbit
				
		implicit none
        integer, intent(in) :: JJi,JJf  ! 2 x initial, final J of states
        logical, intent(in) ::zerome  ! needed to reset		
		
!...................... INTERNAL..........................
        real cleb

        integer :: dpair,cpair
	    integer :: m, par
	    integer :: ia,ib,ic,id
	    integer :: i,j,k,l
	  	integer :: isps,ith,jsps,jth,ksps,kth,lsps,lth
		integer :: ji,jj,jk,jl,mi,mj,mk,ml
		integer :: phaseij,phasekl
		integer :: pair1,pair2,pair11,pair22
		integer :: itbme,indx,indxtr		
		integer :: pcref,pcstart
		integer :: cindx,dindx
		integer :: iref, istart
		
		logical :: first
		integer :: ppar
		integer :: nsection, nme
		
		integer :: Jab,Jcd,Jtot
		integer :: Jabmin,Jabmax,Jcdmin,Jcdmax,Jmin,Jmax
		
		real(8) :: vtmp,vtmphc
		
		integer :: pairblocksize
		
		integer :: isection
		real :: cleb4reduce

        do isection = 1,npnsection
	    dindx = 1
        do dpair = pnpairstart(isection),pnpairstop(isection) ! loop over destruction pairs
            m = PN2%pair(dpair)%m
            par = PN2%pair(dpair)%par


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
           ic = hspsqn(kth,ksps)%orb
           jk = hspsqn(kth,ksps)%j
           mk = hspsqn(kth,ksps)%m
!           gk = hspsqn(kth,ksps)%group

           if( l < 0)then
             lsps = -l
             lth = -2
           else
             lsps = l
             lth  = 2
           endif
           id = hspsqn(lth,lsps)%orb
           jl = hspsqn(lth,lsps)%j
           ml = hspsqn(lth,lsps)%m
           phasekl = 1

!-------------- FIND INDEX "PAIR1" used to find coupled tbme

           pair1 = numorb(2)*(ic-1)+id
		   
           pair1 = PNcouplemap(pair1)
           if(pair1 == -1)then
              print*,' uh problem pn map boss 1 '
              stop
           endif
		   
		   Jcdmax = (jk+jl)/2
		   Jcdmin = abs(jk-jl)/2
		   Jcdmin = max(Jcdmin,abs(m))
		   
	  	   if(par ==1)then
	           pairblocksize = PNcouples%meref(2)
	  	   else
	           pairblocksize = ncouplesPN-PNcouples%meref(2)
		 
	  	   end if

!----------- SET UP MAPPING ARRAYS----------------------
!           dpnpair(l,k) = dindx + nme
!------------- TIME REVERSE -----------
!           if(useTR .and. m== 0)then
!             phasepnpair(l,k) = 1
!           endif
!          if(useTR .and. m/= 0)then
!            phasepnpair(l,k) = 1
!            ktr = hspsqn(kth,ksps)%tr
!            if(kth >0 .and. mk /= 0)ktr = -ktr
!            ltr = hspsqn(lth,lsps)%tr
!            if(lth >0 .and. ml /= 0)ltr = -ltr
!            phasepnpair(ltr,ktr) = int((-1)**( (jl+jk)/2),1)
!            dpnpair(ltr,ktr) = dindx + nme
  !          print*,l,k,dpnpair(l,k),ltr,ktr,dpnpair(ltr,ktr)
!          endif
!----------------- LOOP OVER CREATION PAIRS--------------
           cindx = 0
           do cpair = pnpairstart(isection),pnpairstop(isection)
			   if(m/= PN2%pair(cpair)%m)cycle
			   if(par/= PN2%pair(cpair)%par)cycle
     
              cindx = cindx + 1
			  			  
              itbme = dindx + (cindx-1)*(pnpairstop(isection)+1-pnpairstart(isection)) + npnmesum(isection)
!              itbmetr = cindx + (dindx-1)*nsection +nme ! "time-reversed"

              if(itbme > nmatpn)then ! .or. itbmetr > nmatpn)then
                 print*,' me label too large'
                print*,itbme, nmatpn
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
              ia = hspsqn(ith,isps)%orb
              ji = hspsqn(ith,isps)%j
              mi = hspsqn(ith,isps)%m
              if( j < 0)then
                jsps = -j
                jth = -2
              else
               jsps = j
               jth  = 2
              endif
              ib = hspsqn(jth,jsps)%orb
              jj = hspsqn(jth,jsps)%j
              mj = hspsqn(jth,jsps)%m

!------------------ SET MAPPING ARRAY
!              cpnpair(j,i) = (cindx-1)*nsection
!            if(useTR .and. m/= 0)then   ! NOTE FOR NOW, TR NOT USED
!              itr = hspsqn(ith,isps)%tr
!              if(ith >0 .and. mi/=0)itr = -itr
!              jtr = hspsqn(jth,jsps)%tr
!              if(jth >0 .and. mj /= 0)jtr = -jtr
!              cpnpair(jtr,itr) = (cindx-1)*nsection
!            endif

!---------------- BECAUSE I HAVE TIME REVERSAL, I CAN SKIP A LOT
!            if(cpair > dpair)cycle
		  
!------------ IF TWO-BODY DENSITIES, DO NOT UNCOUPLE MATRIX ELEMENTS

!            if(menu_char=='2b')cycle		  

!------------- PUT INTO "STANDARD" ORDER; MAY GET PHASE -------

             phaseij = 1

             pair2 = numorb(2)*(ia-1)+ib
             pair2 = PNcouplemap(pair2)
             if(pair2==-1)then
                print*,' oh boy should not have gotten that '
                stop
             endif

             pcref = PNcouples%meref(par)
             pcstart = PNcouples%mestart(par)
			 
             indx = pairblocksize*(pair1-pcref-1)+pair2-pcref +pcstart

!             if(pair1 < pair2)then
!               pair11 = pair2
!               pair22 = pair1
!             else
!               pair11 = pair1
!               pair22 = pair2
!             endif
!             indx = (pair11-pcref)*(pair11-pcref-1)/2+pair22-pcref +pcstart
			 			 
             if(indx <= 0 .or. indx > nv2bmedim(0))then
                print*,' problem indx ',indx,pair1,pair2,nv2bmedim(0)
                print*,ia,ib,ic,id
                stop
             endif
			 
			 if(indx > nv2bmedim(0))then
				 print*,'error ',pairblocksize,dpair,cpair,indx,nv2bmedim(0)
				 
				 print*,PNcouples%meref(2) , ncouplesPN - PNcouples%meref(2) 
				 				 
			 end if

!-------------FIND JMIN,JMAX ALLOWED ----------------

!-------------NOW EXTRACT FROM TBMEs ---------------
!-------      NOTE FACTOR ZETA = SQRT(1+DELTA(A,B)) INCLUDED
              vtmp = dmatpn(itbme) !*phasekl*phaseij !*hfactor
!              vtmphc = dmatpnhc(itbme) !*phasekl*phaseij !*hfactor

              if(vtmp==0.0 .and. .not.zerome)cycle
!			  if(ia==1 .and. ib==1 .and. ic==1 .and. id==1 .and.zerome)print*,i,j,k,l,itbme,vtmp,vtmphc,zerome
			 
  		     Jabmax = (ji+jj)/2
  		     Jabmin = abs(ji-jj)/2
			 
			 Jabmin = max(Jabmin,abs(m))
			 			 
			 if(jcdmax > pn2bden(indx)%Jcdmax)then
				 
				 print*,' whoops ',indx
				 print*,pair1,pair2,pairblocksize
				 print*,pair1,ic,id,Jcdmax,pair2,pairblocksize
				 stop
			 end if
			 do Jab = Jabmin,Jabmax
			 	do Jcd = Jcdmin,Jcdmax
		
			        Jmin = pn2bden(indx)%Jmin
			 		Jmin = max(Jmin,abs(Jab-Jcd))
			        Jmax = pn2bden(indx)%Jmax
			 		Jmax = min(Jmax,Jab+Jcd)
			 		if(Jmax < Jmin)cycle
			 		if(zerome)then
			 			pn2bden(indx)%v(:,:,:) = 0.0
			 			cycle
			 		end if
			        do jtot = jmin,jmax
			  
			 			  cleb4reduce = cleb(jtot*2,0,JJi,Jz,JJf,Jz)
			 			  if(abs(cleb4reduce) < 0.00001)then
			 				  pn2bden(indx)%v(Jab,Jcd,Jtot)= -999.
			 				  cycle
			 			  end if
			 			  pn2bden(indx)%v(Jab,Jcd,Jtot)=pn2bden(indx)%v(Jab,Jcd,Jtot) + & 
			 			  vtmp * & 
			                     cleb(jk,mk,jl,ml,2*jcd,2*m)*cleb(ji,mi,jj,mj,2*jab,2*m)*cleb(2*jab,2*m,2*jcd,-2*m,2*Jtot,0)  & 
			 					* (-1)**(Jcd -m)   & ! PHASE FROM time-reversal of second pair
			 					* sqrt(JJf+1.0)*(-1)**(Jtot+ (JJf-JJi)/2)/ cleb4reduce &    ! reduce
			 					/sqrt(2.0*Jtot+1.0)    ! needed for proper definition of density matrix
 
			        enddo      

			     end do  ! Jcd
			 end do   ! Jab
          enddo  ! loop over cpair
          dindx = dindx + 1
        enddo  ! loop over dpair
	end do ! isection

		return
	end subroutine couple_2bdensPN
!==================================================================================
	
!........ SET UP TWO-BODY DENSITY MATRIX......
!         because these files can be large, create a separate file

    subroutine setup_2bdens_output
		
		use nodeinfo
		use io
		use program_info
		use sporbit
		use densities,only:pndensities
		use system_parameters
		implicit none
		integer :: ilast
	    character(10) :: date,time,zone
	    integer :: datetimeval(8)
		integer :: iorb
		
		if(iproc/=0)return
		
!.......... OPEN FILE...................

        ilast = index(outfile,' ')-1
		
        open(unit=dens2bfile,file=outfile(1:ilast)//".den2b",status = 'unknown')
		
		print*,' '
		print*,' Two-body densities will be written to: ',outfile(1:ilast),".den2b "
		if(diag_den2b)print*,' (Expectation values of two-body operators written to ',outfile(1:ilast),".res  )"
		print*, ' '		


!...........WRITE HEADER WITH S.P. INFORMATION.....		
        write(dens2bfile,*)'!#  Two-body densities from BIGSTICK version ',version,lastmodified
	    call date_and_time(date,time,zone,datetimeval)
	    write(dens2bfile,*)'!#  Run date: ',date(1:4),'-',date(5:6),'-',date(7:8)
			write(dens2bfile,*)'!# Densities written in explicit proton-neutron formalism '
			write(dens2bfile,*)'!# Single-particle orbits information follows '
			write(dens2bfile,*)'!# Number of single-particle orbits (same for both protons, neutrons ) '
			write(dens2bfile,*)numorb(1)
		    write(dens2bfile,*)'!#      '
			
			write(dens2bfile,*)'!# Proton, neutron    N      L   2xJ '
			do iorb = 1,numorb(1)
				write(dens2bfile,'(5x,2i5,2x,3i7)')iorb,iorb+numorb(1),orbqn(1,iorb)%nr,orbqn(1,iorb)%l,orbqn(1,iorb)%j	
			end do
			
			if(diag_den2b .and. pndensities)then
			    write(resultfile,*)'!#    orbit labels '
			
				write(resultfile,*)'!# Proton, neutron    N      L   2xJ '
				do iorb = 1,numorb(1)
					write(resultfile,'(5x,2i5,2x,3i7)')iorb,iorb+numorb(1),orbqn(1,iorb)%nr,orbqn(1,iorb)%l,orbqn(1,iorb)%j	
				end do				
				
			end if
			if(diag_den2b .and. .not.pndensities)then			
				write(resultfile,*)'!#    Orbit     N       L    2xJ '
				do iorb = 1,numorb(1)
					write(resultfile,'(5x,i5,2x,3i7)')iorb,orbqn(1,iorb)%nr,orbqn(1,iorb)%l,orbqn(1,iorb)%j	
				end do				
				
			end if

		write(dens2bfile,*)'!#   Z   N '
		write(dens2bfile,'(3x,2i4)')np(1),np(2)
		
		
		return
		
		
	end subroutine setup_2bdens_output	

	
! ROUTINE FOR PRINTING OUT COUPLED TWO-BODY DENSITIES	
! must have already called couple_2bdnesXX,PN
	
	subroutine print_out_2bdens(istate,fstate,Ei,Ef,xJi,xJf,xTi,xTf)
		
		use interaction
		use coupledmatrixelements
		use sporbit
		use densities,only:pndensities
		use nodeinfo
		use system_parameters
		use io
		implicit none
		integer,intent(in) :: istate,fstate  ! label of initial,final many-body state
		real(4),intent(in) :: Ei,Ef          ! energies of initial,final state
		real(4),intent(in) :: xJi,xJf,xTi,xTf  ! J,T (if not in pn-format) for initial,final state
		integer :: it
        integer(kind=4) :: abcouple,cdcouple  ! "couples" of orbits
		integer(kind=4) :: a,b,c,d
		integer(kind=4) :: Jab,Jcd
        integer Jtot,jmin,jmax
		integer(8) :: indx, parblocksize
		
		type (vjs2b), pointer :: x2bden(:)
			
		integer(4) :: ipar
		integer(4) :: couplestart,coupleend	
		
		real(8) :: test,sumrule,zetafactor,trfactor,zetaab,zetacd
		real(8) :: sumrulepp,sumrulenn,sumrulepn  ! intermediate for testing
		integer :: i
		
		sumrule =0.d0 
		sumrulepn = 0.d0
		sumrulepp = 0.d0
		sumrulenn = 0.d0
		if(iproc/=0)return
		write(dens2bfile,*)'!#  '
		write(dens2bfile,*)'!#  '
		
		
		if(diag_den2b)then
  		   write(dens2bfile,*)'!#   State      Energy         Jstate     '
  		   write(dens2bfile,'(i10,6x,f10.5,5x,f5.1)')istate,Ei,xJi
 		   write(dens2bfile,*)'!# a   b   c   d   Jpair   rho_(Jpair,Jpair:J=0)   '
		   
  		   write(resultfile,*)'!#   State      Energy         Jstate     '
  		   write(resultfile,'(i10,6x,f10.5,5x,f5.1)')istate,Ei,xJi
 		   if(pndensities)then
			   write(resultfile,*)'!# a   b   c   d     J     X_J(ab,cd) '
		   else
			   write(resultfile,*)'!# a   b   c   d     J     T      X_JT(ab,cd) '
		   endif 
		   
	   else
		   write(dens2bfile,*)'!# Ini state      Energy         Ji     '
		   write(dens2bfile,'(i10,6x,f10.5,5x,f5.1)')istate,Ei,xJi
		   write(dens2bfile,*)'!# Fin state      Energy         Jf     '
		   write(dens2bfile,'(i10,6x,f10.5,5x,f5.1)')fstate,Ef,xJf		
		   write(dens2bfile,*)'!# a   b    Jab    c   d    Jcd  Jmin Jmax  rho_(Jab,Jcd:J), J=Jmin,Jmax)'
		   
		end if   
		   
		
		do it = 1,2
			if(np(it) < 2)cycle
			
			if(it==1)then
				x2bden=>pp2bden
				
			else
				x2bden=>nn2bden
		
			end if
					
		    do abcouple = 1,ncouplesXX(it)
				a = XXcouples(it)%pairc(abcouple)%ia
				b = XXcouples(it)%pairc(abcouple)%ib
								
				if(it==2)then
					a = a+numorb(1)
					b = b+numorb(1)
					
				end if
				
!--- GET PARITY			

                if(abcouple <= XXcouples(it)%meref(2))then
					ipar = 1
					parblocksize = XXcouples(it)%meref(2)
					couplestart = 1
					coupleend   = XXcouples(it)%meref(2)
				else
					ipar = 2
					parblocksize = ncouplesXX(it)-XXcouples(it)%meref(2)
					couplestart = XXcouples(it)%meref(2)+1
					coupleend   = ncouplesXX(it)
					
				end if	

				
				do cdcouple= couplestart,coupleend !1,ncouplesXX(it)
					
					
					c = XXcouples(it)%pairc(cdcouple)%ia
					d = XXcouples(it)%pairc(cdcouple)%ib
					
					zetafactor = 1.0
					if(a==b)then   
						zetafactor = zetafactor/sqrt(2.)
					end if
					if(c==d)then   
						zetafactor = zetafactor/sqrt(2.)
					end if
					if(it==2)then
						c = c+numorb(1)
						d = d+numorb(1)
					
					end if	
					
					if(diag_den2b)then  ! prevent writing out time-reversed states
						if(a < c)cycle
						if(a==c .and. b < d)cycle
					end if			
					
					if(a==c .and. b==d)then
						trfactor = 1   
					else
						trfactor = 2   ! used to account for Time reversal  (ab,cd) = (cd,ab)
					end if	
					

!--- COMPUTE INDEX ........

                    indx = (cdcouple - XXcouples(it)%meref(ipar) -1 )*parblocksize +  		 & 
					       abcouple - XXcouples(it)%meref(ipar) + XXcouples(it)%mestart(ipar)
						   										
!----- LOOP OVER Js					
                    do Jab = x2bden(indx)%Jabmin,x2bden(indx)%Jabmax
						
	                    do Jcd = x2bden(indx)%Jcdmin,x2bden(indx)%Jcdmax
							if(diag_den2b)then
								Jmin=0
								Jmax= 0
							else
							    Jmin = abs(Jab-Jcd)
							    Jmin = max(Jmin, int(abs(xJi-xJf)))
						 	    Jmax = (Jab + Jcd)
							    Jmax = min(Jmax,int(xJi+xJf))
							end if

							if(Jmin > Jmax)cycle
							if(Jmin < x2bden(indx)%Jmin .or.   Jmax > x2bden(indx)%Jmax )then  !ERROR TRAP
							   print*,' ALARUM!  Some problem with min, max Js '
							   print*,Jab,Jcd,Jmin,Jmax

							   stop
						    end if		
							test =  0.d0
							do Jtot= Jmin,Jmax
								if(x2bden(indx)%v(Jab,Jcd,Jtot) < -900)cycle
								test = test + abs(x2bden(indx)%v(Jab,Jcd,Jtot))
							end do
							if(test < 1.0e-5)then
								cycle					
							endif
							if(diag_den2b)then
								write(dens2bfile,'(x,4i4,i6,2x,f13.7)')a,b,c,d,Jab, & 
								   (x2bden(indx)%v(Jab,Jcd,0))*zetafactor
								   
   								if(pndensities)write(resultfile,'(x,4i4,i6,2x,f13.7)')a,b,c,d,Jab, & 
   								   x2bden(indx)%v(Jab,Jcd,0)*sqrt(2*Jab+1.0)*zetafactor*trfactor/sqrt(2*xJi+1.0)
								   if(a==c .and. b==d)then
									   sumrule= sumrule + real(2*sqrt(2*Jab+1.0)*zetafactor*x2bden(indx)%v(Jab,Jcd,0)/sqrt(2*xJi+1.0),8)
									   if(it==1)sumrulepp= sumrulepp + 2*sqrt(2*Jab+1.0)*zetafactor*x2bden(indx)%v(Jab,Jcd,0)/sqrt(2*xJi+1.0)
									   if(it==2)sumrulenn= sumrulenn + 2*sqrt(2*Jab+1.0)*zetafactor*x2bden(indx)%v(Jab,Jcd,0)/sqrt(2*xJi+1.0)
									   
								   end if
							else
								
							write(dens2bfile,'(x,2i4,i6,2x,2i4,i6,2x,2i4,2x,20f13.7)')a,b,Jab,c,d,Jcd,Jmin,Jmax, & 
							   (x2bden(indx)%v(Jab,Jcd,Jtot)*zetafactor,Jtot=Jmin,Jmax)
							   if(a==c .and. b==d .and. istate==fstate)then
								   sumrule= sumrule + real(2*sqrt(2*Jab+1.0)*zetafactor*x2bden(indx)%v(Jab,Jcd,0)/sqrt(2*xJi+1.0),8)
								   if(it==1)sumrulepp= sumrulepp + 2*sqrt(2*Jab+1.0)*zetafactor*x2bden(indx)%v(Jab,Jcd,0)/sqrt(2*xJi+1.0)
								   
							   end if
							   
						    end if											
						end do ! Jcd
						
					end do !  Jab
		
				end do

			end do    !ab couple
		end do   !it
!...................  NOW FOR PROTON-NEUTRON

        if(np(1)*np(2) > 0)then		
		    do abcouple = 1,ncouplesPN
				a = PNcouples%pairc(abcouple)%ia
				b = PNcouples%pairc(abcouple)%ib
				
				zetafactor = 1.0

				
			    b = b+numorb(1)				
!--- GET PARITY			

                if(abcouple <= PNcouples%meref(2))then
					ipar = 1
					parblocksize = PNcouples%meref(2)
					couplestart = 1
					coupleend   = PNcouples%meref(2)
				else
					ipar = 2
					parblocksize = ncouplesPN-PNcouples%meref(2)
					couplestart = PNcouples%meref(2)+1
					coupleend   = ncouplesPN
					
				end if	
				
				do cdcouple= couplestart,coupleend !1,ncouplesXX(it)
					
					
					c = PNcouples%pairc(cdcouple)%ia
					d = PNcouples%pairc(cdcouple)%ib

					d = d+numorb(1)
					
					if(diag_den2b)then  ! prevent writing out time-reversed states
						if(a < c)cycle
						if(a==c .and. b < d)cycle
					end if	
					if(a==c .and. b==d)then
						trfactor = 1   
					else
						trfactor = 2   ! used to account for Time reversal  (ab,cd) = (cd,ab)
					end if	
!--- COMPUTE INDEX ........

                    indx = (abcouple - PNcouples%meref(ipar) -1 )*parblocksize +  		 & 
					       cdcouple - PNcouples%meref(ipar) + PNcouples%mestart(ipar)
						   										
!----- LOOP OVER Js					
                    do Jab = pn2bden(indx)%Jabmin,pn2bden(indx)%Jabmax
						
	                    do Jcd = pn2bden(indx)%Jcdmin,pn2bden(indx)%Jcdmax
							if(diag_den2b)then
							   Jmin=0
							   Jmax= 0
							else
							   Jmin = abs(Jab-Jcd)
							   Jmin = max(Jmin, int(abs(xJi-xJf)))
							   Jmax = (Jab + Jcd)
							   Jmax = min(Jmax,int(xJi+xJf))
						    end if
							if(Jmin > Jmax)cycle
							if(Jmin < pn2bden(indx)%Jmin .or.   Jmax > pn2bden(indx)%Jmax )then  !ERROR TRAP
							   print*,' ALARUM!  Some problem with min, max Js in pn'
							   print*,Jab,Jcd,Jmin,Jmax
							   stop
						    end if							
							! check to see if nonzero densities
							test =  0.d0
							do Jtot= Jmin,Jmax
								if(pn2bden(indx)%v(Jab,Jcd,Jtot) < -900)cycle
								test = test + abs(pn2bden(indx)%v(Jab,Jcd,Jtot))
							end do
!							print*,test
							if(test < 1.0e-5)then
								cycle					
							endif
							if(diag_den2b)then
								write(dens2bfile,'(x,4i4,i6,2x,f13.7)')a,b,c,d,Jab, & 
								   (pn2bden(indx)%v(Jab,Jcd,0))
      							if(pndensities)write(resultfile,'(x,4i4,i6,2x,f13.7)')a,b,c,d,Jab, & 
      								   pn2bden(indx)%v(Jab,Jcd,0)*sqrt(2*Jab+1.0)/sqrt(2*xJi+1.0)  *trfactor !*zetafactor
								   if(a==c .and. b==d)then
									   sumrule= sumrule + real(2*sqrt(2.0*Jab+1.0)*pn2bden(indx)%v(Jab,Jcd,0)/sqrt(2*xJi+1.0),8)
									   sumrulepn= sumrulepn + 2*sqrt(2.0*Jab+1.0)*pn2bden(indx)%v(Jab,Jcd,0)/sqrt(2*xJi+1.0)
									   
								   end if
								   
							else							
							  write(dens2bfile,'(x,2i4,i6,2x,2i4,i6,2x,2i4,2x,20f13.7)')a,b,Jab,c,d,Jcd,Jmin,Jmax, & 
							   (pn2bden(indx)%v(Jab,Jcd,Jtot),Jtot=Jmin,Jmax)
							   if(a==c .and. b==d .and. istate==fstate)then
								   sumrule= sumrule + real(2*sqrt(2.0*Jab+1.0)*pn2bden(indx)%v(Jab,Jcd,0)/sqrt(2*xJi+1.0),8)
								   sumrulepn= sumrulepn + 2*sqrt(2.0*Jab+1.0)*pn2bden(indx)%v(Jab,Jcd,0)/sqrt(2*xJi+1.0)
								   
							   end if
						    endif											
						end do ! Jcd
						
					end do !  Jab
		
				end do

			end do    !ab couple		
		end if
		sumrule = sumrule 
!		if(diag_den2b)then
        if(istate==fstate)then
			print*,'!# Sum rule = ',sumrule	, sumrulepp,sumrulenn,sumrulepn
			write(dens2bfile,*)'!# Sum rule = ',sumrule			
		end if
		write(dens2bfile,*)'!#  '
		write(dens2bfile,*)'!# Two-body density matrices defined in  '
		write(dens2bfile,*)'!# BIGSTICK Manual, section 5.2, eqn (5.19) '
		write(dens2bfile,*)'!# rho(Jab,Jcd:J) = ' 
		write(dens2bfile,*)'!# - < Jf || [ [a^+ x b^+]_Jab x [ ~c x ~d ]_Jcd ]_J || Ji> '
		write(dens2bfile,*)'!#  ---------------------------------------------------------'
		write(dens2bfile,*)'!#     sqrt[ (2J +1) ( 1+ delta_ab)(1+delta_cd ) ] '

		if(diag_den2b .and. .not.pndensities)call isospin_expectations(xJi)
		return
	end subroutine print_out_2bdens
	
!=============================================================
!
!  added 7.9.5 to allow us to write out diagonal (expectation value) matrix elements in isospin format
!

subroutine isospin_expectations(xJi)
	use interaction
	use coupledmatrixelements
	use sporbit
	use densities,only:pndensities
	use nodeinfo
	use system_parameters
	use io
	
	real,intent(in) :: xJi
	integer :: ntbmes,ncpair
	integer :: dpair,cpair
	integer :: dparity
	integer :: iref, istart
	integer :: a,b,c,d
	integer :: jmin,jmax
	integer :: ja,jb,jc,jd,jj
	integer :: parblocksize
	integer :: couplestart,coupleend
	integer :: pncouplestart,pncoupleend
	integer :: pnparblocksize
	integer :: ipar
	real :: trfactor,zetafactor
	logical :: same
	integer :: phase
	
	integer :: abcouple,cdcouple, pnabcouple,pncdcouple,pnabcoupleR,pncdcoupleR
	
	real(4) :: isoscale(0:30),isovec(0:30)  ! tmp storage for matrix element
	
	integer :: indx,pnindx
	
	integer :: aerr
	
	logical  :: printnonzeroflags=.true.
	
	if(ncouplesXX(1) /= ncouplesXX(2))then
		print*,' Cannot do isospin format in this case '
		stop
	end if
    do abcouple = 1,ncouplesXX(1)
		a = XXcouples(1)%pairc(abcouple)%ia
		b = XXcouples(1)%pairc(abcouple)%ib
		
		ja = orbqn(1,a)%j
		jb = orbqn(1,b)%j	
		
        if(abcouple <= XXcouples(1)%meref(2))then
			ipar = 1
			parblocksize = XXcouples(1)%meref(2)
			couplestart = 1
			coupleend   = XXcouples(1)%meref(2)
		else
			ipar = 2
			parblocksize = ncouplesXX(1)-XXcouples(1)%meref(2)
			couplestart = XXcouples(1)%meref(2)+1
			coupleend   = ncouplesXX(1)
			
		end if	

		
		do cdcouple= couplestart,coupleend !1,ncouplesXX(it)
			isovec= 0.0
			isoscale = 0.0
			
			c = XXcouples(1)%pairc(cdcouple)%ia
			d = XXcouples(1)%pairc(cdcouple)%ib
			
			if(a==b .or. c==d)then
				same = .true.
			else
				same = .false.
			end if
			
			jc = orbqn(1,c)%j
			jd = orbqn(1,d)%j
			
			jmin = max(abs(ja-jb),abs(jc-jd))/2
			jmax = min(abs(ja+jb),abs(jc+jd))/2
			if(jmax < jmin )cycle
			
			zetafactor = 1.0
			if(a==b)then   
				zetafactor = zetafactor/sqrt(2.)
			end if
			if(c==d)then   
				zetafactor = zetafactor/sqrt(2.)
			end if

			if(a < c)cycle
			if(a==c .and. b < d)cycle
			
			if(a==c .and. b==d)then
				trfactor = 1   
			else
				trfactor = 2   ! used to account for Time reversal  (ab,cd) = (cd,ab)
			end if	
			isoscale = 0.0
			
            indx = (cdcouple - XXcouples(1)%meref(ipar) -1 )*parblocksize +  		 & 
			       abcouple - XXcouples(1)%meref(ipar) + XXcouples(1)%mestart(ipar)
				   
			do jj = jmin,jmax
				if(same .and. mod(jj,2)==1)cycle  ! no T=1 matrix elements

				isovec(jj)= (pp2bden(indx)%v(jj,jj,0) + nn2bden(indx)%v(jj,jj,0)) & 
				          *sqrt(2*Jj+1.0)*zetafactor*trfactor/sqrt(2*xJi+1.0) 
				
			end do	   			
!..................... NOW NEED TO FETCH PROTON-NEUTRON INFORMATION.............
            if(ipar==1)then
	           pnparblocksize = PNcouples%meref(2)
	           pncouplestart = 1
	           pncoupleend   = PNcouples%meref(2)
            else
	           pnparblocksize = ncouplesPN-PNcouples%meref(2)
	           pncouplestart = PNcouples%meref(2)+1
	           pncoupleend   = ncouplesPN
	
            end if	
			do pnabcouple = pncouplestart,pncoupleend
				if(a==PNcouples%pairc(pnabcouple)%ia .and. b==PNcouples%pairc(pnabcouple)%ib)goto 1
			end do ! pnabcouples
			print*,' (1) Could not find couple ',a,b
			stop

			1 continue
			do pncdcouple = pncouplestart,pncoupleend
				if(c==PNcouples%pairc(pncdcouple)%ia .and. d==PNcouples%pairc(pncdcouple)%ib)goto 2
			end do ! pncdcouples
			print*,' (2) Could not find couple ',c,d
			stop 
			
			2 continue
            pnindx = (pnabcouple - PNcouples%meref(ipar) -1 )*pnparblocksize +  		 & 
			       pncdcouple - PNcouples%meref(ipar) + PNcouples%mestart(ipar)
				   
	   		do jj = jmin,jmax
	   			if(same .and. mod(jj,2)==1)cycle  ! no T=1 matrix elements	  
!				if(jj==0)print*,' TEST Y1 ', pn2bden(pnindx)%v(jj,jj,0) *sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)	

				isovec(jj) = isovec(jj) + (1./2.)*pn2bden(pnindx)%v(Jj,Jj,0)*sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)  *trfactor *zetafactor


			end do
	   		do jj = jmin,jmax
	   			if(same .and. mod(jj,2)==0)cycle  ! no T=0 matrix elements	   
				isoscale(jj) = isoscale(jj) + (1./2.)*pn2bden(pnindx)%v(Jj,Jj,0)*sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)  *trfactor *zetafactor
			end do
!............... SWAP   c -- > d
            phase = (-1)**( (jc+jd)/2)
			do pncdcoupleR = pncouplestart,pncoupleend
				if(d==PNcouples%pairc(pncdcoupleR)%ia .and. c==PNcouples%pairc(pncdcoupleR)%ib)goto 3
			end do ! pncdcouples
			print*,' (3) Could not find couple ',d,c
			stop   			
			
			3 continue
            pnindx = (pnabcouple - PNcouples%meref(ipar) -1 )*pnparblocksize +  		 & 
			       pncdcoupleR - PNcouples%meref(ipar) + PNcouples%mestart(ipar)
				   
	   		do jj = jmin,jmax
	   			if(same .and. mod(jj,2)==1)cycle  ! no T=1 matrix elements	   
!				if(jj==0)print*,' TEST Y2 ', pn2bden(pnindx)%v(jj,jj,0) , sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0),trfactor ,zetafactor, & 
!				(-1)**JJ*phase*0.5				
				isovec(jj) = isovec(jj) - (1./2.)*pn2bden(pnindx)%v(Jj,Jj,0)*sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)*trfactor *zetafactor*(-1)**JJ*phase
			end do
	   		do jj = jmin,jmax
	   			if(same .and. mod(jj,2)==0)cycle  ! no T=0 matrix elements	   
				isoscale(jj)=isoscale(jj)+(1./2.)*pn2bden(pnindx)%v(Jj,Jj,0)*sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)*trfactor *zetafactor*(-1)**JJ*phase
			end do		   				   
!............... SWAP   a -- > b, c <---> d
            phase = (-1)**( (ja+jb+jc+jd)/2)
			
			do pnabcoupleR = pncouplestart,pncoupleend
				if(b==PNcouples%pairc(pnabcoupleR)%ia .and. a==PNcouples%pairc(pnabcoupleR)%ib)goto 4
			end do ! pnabcouples
			print*,' (4) Could not find couple ',b,a
			stop  

			4 continue
            pnindx = (pnabcoupleR - PNcouples%meref(ipar) -1 )*pnparblocksize +  		 & 
			       pncdcoupleR - PNcouples%meref(ipar) + PNcouples%mestart(ipar)
				   
	   		do jj = jmin,jmax
	   			if(same .and. mod(jj,2)==1)cycle  ! no T=1 matrix elements	   
!				if(jj==0)print*,' TEST Y3 ', pn2bden(pnindx)%v(jj,jj,0) *sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)*trfactor*zetafactor*(-1)**JJ*phase*0.5				
				isovec(jj) = isovec(jj) + (1./2.)*pn2bden(pnindx)%v(Jj,Jj,0)*sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)  *trfactor *zetafactor*phase
			end do
	   		do jj = jmin,jmax
	   			if(same .and. mod(jj,2)==0)cycle  ! no T=0 matrix elements	   
				isoscale(jj)=isoscale(jj)+(1./2.)*pn2bden(pnindx)%v(Jj,Jj,0)*sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)  *trfactor *zetafactor*phase
			end do		
!............... SWAP   a -- > b  
			phase = (-1)**( (ja+jb)/2)
            pnindx = (pnabcoupleR - PNcouples%meref(ipar) -1 )*pnparblocksize +  		 & 
			       pncdcouple - PNcouples%meref(ipar) + PNcouples%mestart(ipar)
				   
	   		do jj = jmin,jmax
	   			if(same .and. mod(jj,2)==1)cycle  ! no T=1 matrix elements	   
!				if(jj==0)print*,' TEST Y4 ', pn2bden(pnindx)%v(jj,jj,0) *sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0) *trfactor*zetafactor*(-1)**JJ*phase*0.5
				
				isovec(jj) = isovec(jj) - (1./2.)*pn2bden(pnindx)%v(Jj,Jj,0)*sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)*trfactor*zetafactor*(-1)**JJ*phase
			end do
	   		do jj = jmin,jmax
	   			if(same .and. mod(jj,2)==0)cycle  ! no T=0 matrix elements	   
				isoscale(jj)=isoscale(jj)+(1./2.)*pn2bden(pnindx)%v(Jj,Jj,0)*sqrt(2*Jj+1.0)/sqrt(2*xJi+1.0)*trfactor *zetafactor*(-1)**JJ*phase
			end do		   		
!............... finally, write out			
			do jj = jmin,jmax
		   		if(.not. same .or. mod(jj,2) /= 1)then  ! no T=0 matrix elements	   
				   if((abs(isovec(jj)) >  1.e-5).or.printnonzeroflags)write(resultfile,'(x,4i4,2i6,2x,f13.7)')a,b,c,d,jj,1,isovec(jj)
			    end if
	   			if(.not. same .or. mod(jj,2) /= 0)then  ! no T=0 matrix elements	   
				if((abs(isoscale(jj)) > 1.e-5).or.printnonzeroflags)write(resultfile,'(x,4i4,2i6,2x,f13.7)')a,b,c,d,jj,0,isoscale(jj)
			    end if
			end do
				
			
				   											
		end do  ! cd couple
	end do  ! ab couple

	
	
end subroutine isospin_expectations	
	

end module twobodydens_util
	
	
	
