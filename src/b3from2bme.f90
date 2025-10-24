!========================================================================
!
! B3FROM2BME.f90
!
! routines to convert 2-body matrix elements into 3
!
!===========================================================================
!
!  subroutine to convert uncoupled TBMEs to 3-body MEs
!  
!  first version....only PPP and NNN
!  initiated 9/09 by CWJ @ SDSU
!
!  most of this is cribbed from uncoupleXXtbme in BTBME.F90
!
!  second version...includes PPN and PPN
!  started 3/11 by CWJ @ SDSU
!
  subroutine Masterconvert2to3bme
  use system_parameters

  if( np(1) > 2)call convert2to3bmeXXX(1)
  if( np(2) > 2)call convert2to3bmeXXX(2)

  if( np(1)*(np(2)-1) > 0)then
      call convertXXtoXXY(2)
      call convertPNtoXXY2(2)
  endif
  if( np(2)*(np(1)-1) > 0)then
      call convertXXtoXXY(1)
      call convertPNtoXXY2(1)

  endif

  end subroutine Masterconvert2to3bme

  subroutine convert2to3bmeXXX(it)

  use spstate
  use haiku_info
  use system_parameters
  use interaction
  use interactions3body
  use ntuple_info
  use jump_mod
  use welder

  implicit none
  integer it      ! which species

  integer(kind=8) :: dpair,cpair
  integer(kind=8) :: cindx,dindx
  integer T
  integer M,m3
  integer par,par3
  integer w3
  integer ik,il
  integer ii,ij

  integer i,j,k,l,x,ix
  integer ith,jth,kth,lth,xth
  integer isps,jsps,ksps,lsps,Xsps
  integer itmp,jtmp,ktmp,ltmp,xtmp
  integer cphase,dphase
  integer(kind=8) :: iref,istart,iref3,istart3
  integer(kind=8) :: itbme,i3bme
  integer(kind=8) :: nmat
  real fact
  real hfactor    ! HERMITICITY IN 2-BODY

  real,pointer :: v(:),v3(:)
  integer, pointer :: map(:)

  if(Np(it) < 3)return

  fact = 1.0/(np(1)+np(2)-2.0)

  if(it ==1)then
        v => hmatpp
        v3=> hmatppp
        nmat = nmatppp
        map => mapPPP
  else
        v => hmatnn
        v3=> hmatnnn
        nmat = nmatnnn
        map => mapNNN
  endif
  v3(:) = 0.0
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
          k = -k
        else
          ksps = k
          kth  = it
          k = k+nhsps(-it)
        endif
        if( l < 0)then
          lsps = -l
          lth = -it
          l = -l
        else
          lsps = l
          lth  = it
          l = l+nhsps(-it)
        endif

        do cpair = 1,dpair        ! loop over creation pair
!------------ hfactor NO LONGER NEEDED -- CWJ 3/2011
!          if(cpair == dpair)then
!              hfactor = 1.0 ! 2.0      ! this undoes a necessary hermiticity division by 2 back in subroutine uncoupleXXtbme in btbme.f90
!          else
              hfactor = 1.0
!          endif
          if(XX2(it)%pair(cpair)%m /=m)cycle       ! enforce quantum numbers
          if(XX2(it)%pair(cpair)%par /= par)cycle

          if(iref > cpair)then
            print*,' problem with XX iref (c) '
            print*,cpair,iref
            stop
          endif

!...........???? MAP????
!--------------- HOW INDEX FOR UNCOUPLED tbme IS COMPUTED---------

          itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref

          if(itbme > nmatXX(it))then          ! error trap
             print*,' me label too large (XX-> XXX)'
             print*,itbme,nmatXX(it)
             print*,cpair,dpair,iref,istart
             print*,m,par
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
            i = -i
          else
            isps = i
            ith  = it
            i = i + nhsps(-it)
          endif

          if( j < 0)then
            jsps = -j
            jth = -it
            j = -j
          else
            jsps = j
            jth  = it
            j = j+nhsps(-it)
          endif
!............... NOW LOOP OVER EXTRA PARTICLE.....................

          do ix = -nhsps(-it),nhsps(it)
             if(ix == 0)cycle
             if(ix < 0)then
                 xsps = -ix
                 xth = -it
                 x = -ix
             else
                 xsps = ix
                 xth  = it
                 x = ix+nhsps(-it)
             endif
!................ CHECK FOR DUPLICATION..........
             if(x == i .or. x == j ) cycle
             if(x == k .or. x == l ) cycle
!.................ORDER AND PICK UP PHASES...........
             itmp = i
             jtmp = j
             xtmp = x
             cphase = 1
             call swap(jtmp,xtmp,cphase)
             call swap(itmp,jtmp,cphase)
             call swap(jtmp,xtmp,cphase)
             cindx = (itmp-3)*(itmp-2)*(itmp-1)/6 + (jtmp-2)*(jtmp-1)/2 + xtmp
             cindx = map(cindx)

             if(cindx == 0)cycle

             itmp = k
             jtmp = l
             xtmp = x
             dphase = 1
             call swap(jtmp,xtmp,dphase)
             call swap(itmp,jtmp,dphase)
             call swap(jtmp,xtmp,dphase)
             dindx = (itmp-3)*(itmp-2)*(itmp-1)/6 + (jtmp-2)*(jtmp-1)/2 + xtmp
             dindx = map(dindx)

             if(dindx ==0)cycle

             w3 = hspsqn(ith,isps)%w + hspsqn(jth,jsps)%w + hspsqn(xth,xsps)%w
             m3 = hspsqn(ith,isps)%m + hspsqn(jth,jsps)%m + hspsqn(xth,xsps)%m
             par3 = hspsqn(ith,isps)%par* hspsqn(jth,jsps)%par* hspsqn(xth,xsps)%par
             par3 = (3-par3)/2
!................. CHECK TO SEE IF ALLOWED.....................
            if(restrict_ntuples)then
!               if(w3 > tripletlim(it)%maxw .or. w3 < tripletlim(it)%minw)cycle
               if( par3 > tripletlim(it)%w(w3)%parmax) cycle
               if( par3 < tripletlim(it)%w(w3)%parmin) cycle
               if( m3 > tripletlim(it)%w(w3)%par(par3)%jzmax) cycle
               if( m3 < tripletlim(it)%w(w3)%par(par3)%jzmin) cycle
            end if

             iref3 = XXX3(it)%meref(m3,par3)
             istart3 = XXX3(it)%mestart(m3,par3)


             if(cindx <= dindx) then
                 i3bme = istart3 + (dindx-iref3)*(dindx-iref3-1)/2+cindx-iref3
             else
                 i3bme = istart3 + (cindx-iref3)*(cindx-iref3-1)/2+dindx-iref3
             endif
 
             if(cindx < iref3 .or. dindx < iref3)then
                print*,' Problem with indices, refs '
                print*,cindx,dindx,iref3
                print*,m3,par3
                print*,i,j,x
                print*,k,l,x
                print*,ix
                print*,xx2(it)%pair(dpair)%ia, xx2(it)%pair(dpair)%ib
             itmp = i
             jtmp = j
             xtmp = x
             cphase = 1
             call swap(jtmp,xtmp,cphase)
             call swap(itmp,jtmp,cphase)
             call swap(jtmp,xtmp,cphase)
             cindx = (itmp-3)*(itmp-2)*(itmp-1)/6 + (jtmp-2)*(jtmp-1)/2 + xtmp
!             print*,itmp,jtmp,xtmp,cindx
!                print*,map
                stop
             endif

             if(i3bme > nmat)then
                print*,' index is too large ',i3bme,nmat,m3,par3
                print*,istart3,iref3
                print*,dindx,cindx
                stop
             endif
             if(i3bme < 1)then
                 print*,it
                 print*,i3bme,' i3bme ',dindx, cindx, iref3,istart3
                 print*,m3,par3,w3
                 print*,isps,jsps,xsps
                 print*,ksps,lsps,xsps
                 print*,itmp,jtmp,xtmp
             dindx = (itmp-3)*(itmp-2)*(itmp-1)/6 + (jtmp-2)*(jtmp-1)/2 + xtmp
                 print*,dindx,map(dindx)
!                 print*,map
                 stop
             endif 
             v3(i3bme) = v3(i3bme) +v(itbme)*cphase*dphase*fact*hfactor

          enddo  ! n

        enddo  ! loop over cpair
  enddo  ! loop over dpair

  return
  end subroutine convert2to3bmeXXX

!==========================================================================

  subroutine convertXXtoXXY(it)

  use spstate
  use haiku_info
  use system_parameters
  use interaction
  use interactions3body
  use ntuple_info
!  use hermit
  use nodeinfo
  use bmpi_mod
  implicit none

  integer it      ! which species
  integer itc
  integer dpair,cpair,cindx,dindx
  integer T
  integer M,m3
  integer par,par3
  integer w3
  integer ik,il
  integer ii,ij
  integer(kind=8) :: coplabel, doplabel
  integer,pointer :: cXXYtriplet(:,:), dXXYtriplet(:,:)
  integer i,j,k,l,x,ix
  integer ith,jth,kth,lth,xth
  integer isps,jsps,ksps,lsps,Xsps
  integer itmp,jtmp,ktmp,ltmp,xtmp
  integer cphase,dphase
  integer(kind=8) :: iref,istart,iref3,istart3
  integer(kind=8) :: itbme,i3bme, i3bmetr
  integer(kind=8) :: nmat,ncount
  real fact
  real hfactor    ! HERMITICITY IN 2-BODY

  real,pointer :: v(:),v3(:)
!  integer, pointer :: map(:,:)

  itc = 3 -it
  if( (np(itc)*(np(it)-1) < 1)) return

  fact = 1.0/(np(1)+np(2)-2.0)
  ncount = 0
  if(it ==1)then
        v => hmatpp
        v3=> hmatppn
        nmat = nmatppn
!        map => mapPPN
        cXXYtriplet => cppntriplet
        dXXYtriplet => dppntriplet
  else
        v => hmatnn
        v3=> hmatpnn
        nmat = nmatpnn
!        map => mapPNN
        cXXYtriplet => cpnntriplet
        dXXYtriplet => dpnntriplet
  endif
  v3(:) = 0.0
  do dpair = 1,npairXX(it)        ! loop over destruction pair

        m = XX2(it)%pair(dpair)%m     ! find Jz, parity of pair
!        if(useTR .and. m < 0 .and. jz >=0 )cycle
!        if(useTR .and. m > 0 .and. jz < 0 )cycle

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
          k = -k
        else
          ksps = k
          kth  = it
          k = k+nhsps(-it)
        endif
        if( l < 0)then
          lsps = -l
          lth = -it
          l = -l
        else
          lsps = l
          lth  = it
          l = l+nhsps(-it)
        endif
        dphase = 1
        do cpair = 1,dpair        ! loop over creation pair
!------------ hfactor NO LONGER NEEDED -- CWJ 3/2011
!          if(cpair == dpair .and. hermitian)then
!              hfactor = 1.0      ! this undoes a necessary hermiticity division by 2 back in subroutine uncoupleXXtbme in btbme.f90
!          else
              hfactor = 1.0
!          endif
          cphase = 1
          if(XX2(it)%pair(cpair)%m /=m)cycle       ! enforce quantum numbers
          if(XX2(it)%pair(cpair)%par /= par)cycle

          if(iref > cpair)then
            print*,' problem with XX iref (c) '
            print*,cpair,iref
            stop
          endif

!...........???? MAP????
!--------------- HOW INDEX FOR UNCOUPLED tbme IS COMPUTED---------

          itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref

          if(itbme > nmatXX(it))then          ! error trap
             print*,' me label too large (XX-> XXY)',it
             print*,itbme,nmatXX(it)
             print*,cpair,dpair,iref,istart
             print*,m,par
             stop
          endif
          if(itbme <= 0)then
               print*,' problem itbme (XX-> XXY) ',itbme
             print*,cpair,dpair,iref,istart
             print*,m,par
               stop
          endif

!--------------- GET Q#s FOR CREATION PAIRS ---------------

          i = XX2(it)%pair(cpair)%ia
          j = XX2(it)%pair(cpair)%ib

          if( i < 0)then
            isps = -i
            ith = -it
            i = -i
          else
            isps = i
            ith  = it
            i = i + nhsps(-it)
          endif

          if( j < 0)then
            jsps = -j
            jth = -it
            j = -j
          else
            jsps = j
            jth  = it
            j = j+nhsps(-it)
          endif
!............... NOW LOOP OVER EXTRA PARTICLE.....................

          do ix = -nhsps(-itc),nhsps(itc)
             if(ix == 0)cycle
             if(ix < 0)then
                 xsps = -ix
                 xth = -it
                 x = -ix
             else
                 xsps = ix
                 xth  = it
                 x = ix+nhsps(-itc)
             endif
!................ COMPRESSION VIA TR .............................
             if (useTR) then
                m3 = 2*m + hspsqn(xth,xsps)%m
                if( m3 < 0 .and. jz >= 0)cycle
                if( m3 > 0 .and. jz < 0 ) cycle
             end if
             coplabel = cXXYtriplet(cpair,ix)
             doplabel = dXXYtriplet(dpair,ix)
             if(doplabel /= -1 .and. coplabel /=-1)then
		!		 print*,' Umm ',cpair,ix,coplabel, dpair,ix,doplabel
				 
             i3bme = coplabel + doplabel 
             ncount = ncount + 1
             if(i3bme > nmat)then
                print*,' index is too large ',i3bme,nmat,m3,par3
                print*,istart3,iref3
                print*,dindx,cindx
                stop
             endif
             if(i3bme < 1)then
                 print*,it
                 print*,i3bme,' i3bme ',coplabel, doplabel
                 print*,cpair,dpair,ix
                 print*,m,par

                 stop
             endif 

             v3(i3bme) = v3(i3bme) +v(itbme)*cphase*dphase*fact*hfactor
		 else
                i3bme = -1
             endif
!-------------------- TIME REVERSED
             coplabel = cXXYtriplet(dpair,ix)
             doplabel = dXXYtriplet(cpair,ix)

             if(doplabel /= -1 .and. coplabel /=-1)then
             i3bmetr = coplabel + doplabel 

             ncount = ncount + 1
             if(i3bmetr > nmat)then
                print*,' index is too large ',i3bme,nmat,m3,par3
                print*,istart3,iref3
                print*,dindx,cindx
                stop
             endif
             if(i3bmetr < 1)then
                 print*,it
                 print*,i3bmetr,' i3bme ',coplabel, doplabel
                 print*,cpair,dpair,ix
                 print*,m,par

                 stop
             endif 
             if(i3bmetr /= i3bme )then
             v3(i3bmetr) = v3(i3bmetr) +v(itbme)*cphase*dphase*fact*hfactor
             endif
             endif
          enddo  ! n

        enddo  ! loop over cpair
  enddo  ! loop over dpair

  if(iproc == 0)then
  if(it == 1)then
     print*,ncount,' PPN matrix elements converted from PP '
  else
     print*,ncount,' PNN matrix elements converted from NN '
  endif
  endif
!.......... ZERO OUT OLD 2BMEs............

  return
  end subroutine convertXXtoXXY

!==========================================================================

!   based on   subroutine count_uncouplepn
!
!
      subroutine convertPNtoXXY2(it)

      use verbosity
!      use hermit
      use spstate
      use sporbit
      use haiku_info
      use system_parameters
      use interaction
      use interactions3body
      use ntuple_info
      use btbme_mod


      implicit none
      integer it,itc

      integer nme

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

      integer x,xsps, xth, mx,gx,jx,nx,ix,m3
      integer x1sps,x2sps
      integer phasexxy, cysps,dysps
      integer(kind=8) :: cxx, dxx
      integer i3bme   ! direct matrix element
      integer i3bmetr ! time-reverse
      integer nk,nl,ni,nj
      integer jk,jl,ji,jj
      integer mk,ml,mi,mj

      integer(kind=8) :: iref,istart
      integer(kind=8) :: itbme,itbmetr
      integer a,b,c,d
      integer ia,ib,ic,id
      integer ap,bn,cp,dn
      integer xtmp,atmp,btmp,ctmp,dtmp
      integer atr,btr,ctr,dtr
      integer itr,jtr,ktr,ltr
      integer pair1,pair2
      integer pair11,pair22
      integer ppar,pcref,pcstart
      integer phasekl,phaseij
      integer indx
      integer itmp,jtmp,ktmp,ltmp
      real vtmp

      ! in module now: real zeta   ! zeta(i,j) = sqrt(1+delta(i,j))
      real cleb
      logical first
      integer nsection
      integer ipair,pairstart
      integer cindx,dindx  ! creation/destruction indices
      integer gi,gj,gk,gl ! group indices
      integer n0
  integer(kind=8) :: coplabel, doplabel
  integer,pointer :: cXXYtriplet(:,:), dXXYtriplet(:,:)
  integer(kind=8), pointer :: xxmap(:)

  real fact
  real hfactor    ! HERMITICITY IN 2-BODY
  real,pointer :: v(:),v3(:)
  integer :: phasexy
  integer itest

  logical TRflip           ! when using TR symmetry to compactify matrix elements
                           ! we restrict the M of both 3-particle and 2-particle
                           ! sometimes these are in conflict (e.g., M2 of pair > 0 
                           ! but M3 of triplet is < 0; the flag TRflip alerts us to 
                           ! this case
  itc = 3 -it
  if( (np(itc)*(np(it)-1) < 1)) return

  fact = 1.0  /(np(1)+np(2)-2.0)
!print*,' factor is ',fact
  if(it ==1)then
        v3=> hmatppn
!        map => mapPPN
        cXXYtriplet => cppntriplet
        dXXYtriplet => dppntriplet
        xxmap => mappairPP
  else
        v3=> hmatpnn
!        map => mapPNN
        cXXYtriplet => cpnntriplet
        dXXYtriplet => dpnntriplet
        xxmap => mappairNN

  endif

  n0 = 0
  nme = 0
  m = -9999
  par = 0
  nsection = 0
  pairstart = -1000000

  do dpair = 1,npairpn       ! loop over PN destruction pairs
!---------- FIND OUT IF THERE IS A CHANGE IN M OR PARITY
        if(par /= PN2%pair(dpair)%par .or. m /= PN2%pair(dpair)%m)then

!------------- ADD TO TOTAL THOSE FROM LAST TIME
          nme = nme + (nsection)*(nsection)

!---------IF SO, THEN COUNT UP HOW MANY HAVE THE SAME M, PARITY
          m = PN2%pair(dpair)%m
          par = PN2%pair(dpair)%par
          pairstart = dpair
!------------------ FIND HOW FAR TO GO ------------------------------
          do ipair = pairstart,npairpn  ! loop upwards until M, PARITY change
                                      ! pairs have previously been sorted
             if(par == PN2%pair(ipair)%par .and. m == PN2%pair(ipair)%m)then
!------------ NSECTION IS HOW MANY PAIRS HAVE SAME M, PARITY
                nsection = ipair -pairstart+1
             else
                exit
             endif
          enddo
          first = .true.
          dindx = 1
        else
          first = .false.
          dindx = dindx + 1
        endif

!-------------- EXTRACT QUANTUM NUMBERS OF DESTRUCTION PAIR--------
        k = PN2%pair(dpair)%ia           ! proton 
        l = PN2%pair(dpair)%ib           ! neutron
        if( k < 0)then
          ksps = -k
          kth = -1
          cp = -k 

        else
          ksps = k
          kth  = 1
          cp = k+nhsps(-kth)

        endif
        ia = hspsqn(kth,ksps)%orb
        jk = hspsqn(kth,ksps)%j
        mk = hspsqn(kth,ksps)%m
        gk = hspsqn(kth,ksps)%group

        if( l < 0)then
          lsps = -l
          lth = -2
          dn = -l

        else
          lsps = l
          lth  = 2
          dn = l+nhsps(-lth)
        endif
        ib = hspsqn(lth,lsps)%orb
        jl = hspsqn(lth,lsps)%j
        ml = hspsqn(lth,lsps)%m
        gl = hspsqn(lth,lsps)%group
        phasekl = 1


!----------------- LOOP OVER CREATION PAIRS--------------
        cindx = 0
        do cpair = pairstart,pairstart+nsection-1
           hfactor = 1.
!           if(.not.hermitian .and. cpair==dpair)hfactor = 0.5
           if(cpair==dpair)hfactor = 0.5

!-------------- EXTRACT QUANTUM NUMBERS OF CREATION PAIR--------

           i = PN2%pair(cpair)%ia     ! proton
           j = PN2%pair(cpair)%ib     ! neutron

           if( i < 0)then
             isps = -i
             ith = -1
             ap = -i

           else
             isps = i
             ith  = 1
             ap = i+nhsps(-ith)

           endif

           ic = hspsqn(ith,isps)%orb
           ji = hspsqn(ith,isps)%j
           mi = hspsqn(ith,isps)%m
           gi = hspsqn(ith,isps)%group
          if( j < 0)then
             jsps = -j
            jth = -2
            bn = -j
          else
            jsps = j
            jth  = 2
            bn = j+nhsps(-jth)

        endif

        id = hspsqn(jth,jsps)%orb
        jj = hspsqn(jth,jsps)%j
        mj = hspsqn(jth,jsps)%m
        gj = hspsqn(jth,jsps)%group

!------------------ SET MAPPING ARRAY
!---------------- BECAUSE I HAVE TIME REVERSAL, I CAN SKIP A LOT
        if(cpair > dpair)cycle

!--------- ALTERNATE RECOVERY OF TWO-BODY MATRIX ELEMENT----

        itbme = cpnpair(j,i ) + dpnpair(l,k)
        phasexy = 1
        if(useTR)phasexy = phasepnpair(j,i)*phasepnpair(l,k)

        vtmp = hmatpn(itbme)*phasexy
        
!----------------- LOOP OVER ADDITIONAL PARTICLE ---------------------------
        do x = 1,nhsps(-it)+nhsps(it)

!................... FIND MATRIX ELEMENTS...............................
          if(it == 1 .and. (ap == x .or. cp == x))cycle
          if(it == 2 .and. (bn == x .or. dn == x))cycle

          phasexxy = 1
          if(it == 1)then
             call findXXYindex(it,ap,x,bn,cp,x,dn, i3bme,i3bmetr, phasexxy)
          else
             call findXXYindex(it,bn,x,ap,dn,x,cp, i3bme,i3bmetr, phasexxy)

          end if

          if(phasexxy ==0)cycle
           if(i3bme /= 0)then
!          if( cxxytriplet(cxx,cysps) /= -1 .and. dxxytriplet(dxx,dysps) /= -1)then


          v3(i3bme) = v3(i3bme)+ vtmp*phasexxy*hfactor*fact
          endif 


!
!  NOTE: ONLY DO HERMITIAN CONJUGATE 
!  IF a /= b for protons or c / =d for neutron

!          i3bmetr = dxxytriplet(cxx,cysps) + cxxytriplet(dxx,dysps)

          if(i3bme /= i3bmetr .and. .not.(it ==1 .and. ap == cp) .and.  & 
                                    .not.( it==2 .and. bn == dn) .and. & 
             i3bmetr /= 0)then


            v3(i3bmetr) = v3(i3bmetr)+ vtmp*phasexxy*hfactor*fact
          endif

        end do ! loop over x
!------------------
        enddo  ! loop over cpair

      enddo  ! loop over dpair
      return
      end subroutine convertPNtoXXY2
!=============================================================
!
! separate routine to find index of V(abc,def) including a phase
!
! INPUT: itx = species of X in XXY
! ax,bx = x-species labels for creation operators
! cy = y-species label for creation operator
! dx,ex = x-species labels for annihilation operators
! fy = y-species label for annihilation operators
!NOTE: all labels (ax...fy) are positive:
! they run from 1 through nhsps(it)+ nhsps(-it)
!     those from 1 to nhsps(-it) refer to those with m < 0
!     those from nhsps(-it) + 1 to nhsps(-it) + nhsps(it) have m>=0
! 
!  OUTPUT:
! i3bme = index of 3-body 
! i3bmehc = index of hermitian conjugate 3-body, i.e. V(def,abc)
! phase = picked up from any reordering
!       = 0 if there is no such index
! 
subroutine findXXYindex(itx,ax,bx,cy,dx,ex,fy,i3bme,i3bmehc,phase)
   use verbosity
!   use hermit
   use spstate
   use sporbit
   use haiku_info
   use system_parameters
   use interaction
   use interactions3body
   use ntuple_info
   implicit none
   integer itx
   integer ax,bx,cy,dx,ex,fy
   integer i3bme,i3bmehc
   integer phase

   integer(kind=8) :: cxx,dxx
   integer,pointer :: cXXYtriplet(:,:), dXXYtriplet(:,:)
   integer(kind=8),pointer :: xxmap(:)
   integer cysps, fysps
   integer(1), pointer :: phaseXXY(:,:)
   integer mxx
   integer asps,bsps,csps,dsps,esps,fsps
   integer ath,bth,cth,dth,eth,fth
   integer a1,b1,c2,d1,e1,f2
   logical flip
   integer ity

   ity = 3 -itx
!........... ERROR TRAP..........
   if(ax < 0 .or. bx < 0 .or. cy < 0 .or. &
      dx < 0 .or. ex < 0 .or. fy < 0)then
     print*,' All s.p. indices in findXXYindex must be > 0 '
     print*,ax,bx,cy,dx,ex,fy
     stop
   end if

   if(itx ==1)then

        cXXYtriplet => cppntriplet
        dXXYtriplet => dppntriplet
        xxmap => mappairPP
        phaseXXY => phaseppn

  else

        cXXYtriplet => cpnntriplet
        dXXYtriplet => dpnntriplet
        xxmap => mappairNN
        phaseXXY => phasepnn

  endif
  phase = 1
  a1 = ax
  b1 = bx
  c2 = cy
  d1 = dx
  e1 = ex
  f2 = fy

!....... CHECK FOR TIME-REVERSAL.........
  flip = .false.

  if(useTR)then
!........ CHECK M-value ..............
     if(ax > nhsps(-itx) )then
        ath = itx
        asps = ax -nhsps(-itx)
     else
        ath = -itx
        asps = ax
     end if
     mxx = hspsqn(ath,asps)%m

     if(bx > nhsps(-itx) )then
        bth = itx
        bsps = bx -nhsps(-itx)
     else
        bth = -itx
        bsps = bx
     end if
     mxx = mxx + hspsqn(bth,bsps)%m
!...........NEED TO FLIP??...........

     if( (mxx < 0 .and. jz >= 0 ) .or.  & 
          (mxx > 0 .and. jz < 0) )then
         flip = .true.
         phase = phase *( -1)**( (hspsqn(ath,asps)%j + hspsqn(bth,bsps)%j)/2)

     end if
     if(flip)then     ! flip everything else
         a1 = hspsqn(ath,asps)%tr
         if(ath < 0)a1 = a1 + nhsps(ath)
         b1 = hspsqn(bth,bsps)%tr
         if(bth < 0)b1 = b1 + nhsps(bth)

         if(dx > nhsps(-itx) )then
            dth = itx
            dsps = dx -nhsps(-itx)
         else
            dth = -itx
            dsps = dx
         end if
         d1 = hspsqn(dth,dsps)%tr
         if(dth < 0)d1 = d1 + nhsps(dth)

         if(ex > nhsps(-itx) )then
            eth = itx
            esps = ex -nhsps(-itx)
         else
            eth = -itx
            esps = ex
         end if
         e1 = hspsqn(eth,esps)%tr
         if(eth < 0)e1 = e1 + nhsps(eth)

         phase = phase*(-1)**( ( hspsqn(dth,dsps)%j + hspsqn(eth,esps)%j)/2)
         if(cy > nhsps(-ity) )then
            cth = ity
            csps = cy -nhsps(-ity)
         else
            cth = -ity
            csps = cy
         end if
         c2 = hspsqn(cth,csps)%tr
         if(cth < 0)c2 = c2 + nhsps(cth)

         if(fy > nhsps(-ity) )then
            fth = ity
            fsps = fy -nhsps(-ity)
         else
            fth = -ity
            fsps = fy
         end if
         f2 = hspsqn(fth,fsps)%tr
         if(fth < 0)f2 = f2 + nhsps(fth)

         phase = phase*(-1)**( ( hspsqn(cth,csps)%j + hspsqn(fth,fsps)%j)/2)

     end if
  end if

  if(a1 > b1)then
     cxx = (a1-1)*(a1-2)/2+b1
  else
     cxx = (b1-1)*(b1-2)/2+a1
     phase = -phase
  end if
  cxx = xxmap(cxx)
  if(cxx == 0)then
    phase = 0
    return
  endif


  if(d1 > e1)then
     dxx = (d1-1)*(d1-2)/2+e1
  else
     dxx = (e1-1)*(e1-2)/2+d1
     phase = -phase
  end if

  dxx = xxmap(dxx)
  if(dxx == 0)then
    phase = 0
    return
  endif

!...... MUST CONVERT cy, fy

  if(c2 <= nhsps(-ity))then
     cysps = -c2
  else
     cysps = c2 - nhsps(-ity)
  end if

  if(f2 <= nhsps(-ity))then
     fysps = -f2
  else
     fysps = f2 - nhsps(-ity)
  end if

  if(cxxytriplet(cxx,cysps) == -1 .or. dxxytriplet(dxx,fysps) == -1)then
      i3bme = 0
  else
      i3bme = cxxytriplet(cxx,cysps) + dxxytriplet(dxx,fysps)
  endif
  if(dxxytriplet(cxx,cysps) == -1 .or. cxxytriplet(dxx,fysps) == -1)then
      i3bmehc = 0
  else
      i3bmehc = dxxytriplet(cxx,cysps) + cxxytriplet(dxx,fysps)
  end if
  
  if(useTR)then
    phase = phase*phasexxy(cxx,cysps)*phasexxy(dxx,fysps)

  end if
  
  return
  end subroutine findXXYindex
!=============================================================
