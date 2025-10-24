!  BDEFORMEDLIB2.f90
!
!  SPECIAL SUBROUTINES for used with deformed/uncoupled 2-body matrix elements
!=============================================================
!
!  CALLED BY: uncoupleXXtbme
!
!  
!
   subroutine readXXdeformedtbmes(it,foundit)

   use verbosity
   use spstate
   use haiku_info
   use system_parameters
   use interaction

   implicit none
   integer it      ! which species   
   logical foundit
   integer ttmp, mtmp,jtmp,ntmp,ltmp
   integer(kind=8) :: mtot,partot
   integer(kind=8), allocatable :: mapsps(:),mapm(:),mappar(:)
   integer nspslocal
   integer(kind=8) :: i
   integer(kind=8) :: ntbme  ! # of uncoupled TBMEs in file
   integer ::  a,b,c,d
   integer(kind=8) :: isps,jsps,ksps,lsps
   integer(kind=8) :: cpair,dpair
   real hfactor
   real vxx
   logical foundsps
   integer(kind=8) :: iref, istart
   integer cphase,dphase
   integer(kind=8) :: itbme
   real,pointer :: v(:)
   integer(kind=8),pointer :: map(:)
   integer jzmin,jzmax
   integer :: aerr

   foundit = .false.
   open(unit=23,file='deformed.int',status='old',err=101)
   foundit = .true.
!...... FIRST READ PAST AND GET MAPPING OF SINGLE PARTICLE STATES.....
   
   if(it==1)then
      map => mappairPP
      if(.not.allocated(hmatpp))then
           allocate(hmatpp(nmatXX(it)), stat=aerr)
           if(aerr /= 0) call memerror("readXXdeformedtbmes 1");
           hmatpp=0.0
      end if
      v => hmatpp
   else
      map => mappairNN
      if(.not.allocated(hmatnn))then
            allocate(hmatnn(nmatXX(it)), stat=aerr)
            if(aerr /= 0) call memerror("readXXdeformedtbmes 2");
            hmatnn=0.0
      end if
      v => hmatnn
   end if
   jzmin = XX2(it)%pair(1)%m      !
   jzmax = XX2(it)%pair(npairXX(it))%m 

   read(23,*)nspslocal
   allocate(mapsps(nspslocal),mapm(nspslocal),mappar(nspslocal) , stat=aerr)
   if(aerr /= 0) then
      call memerror("readXXdeformedtbmes 3");
      stop 5  ! -Wuninitialized doesn't know that memerror won't return
   end if
   mapsps(:) = -1
   do i = 1,nspslocal
      read(23,*)ttmp, ntmp, ltmp,jtmp,mtmp
      if(ttmp /= it)then
         mapsps(i) = -1
         cycle
      end if
      mapm(i) = mtmp
      mappar(i) = (-1)**ltmp
      foundsps = .false.
!............ SEARCH FOR SINGLE-PARTICLE STATES.....
      if(mtmp > 0)then
        do isps = 1,nhsps(it)
          if( mtmp == hspsqn(it,isps)%m .and. jtmp == hspsqn(it,isps)%j  & 
          .and. ntmp == hspsqn(it,isps)%nr .and. ltmp == hspsqn(it,isps)%l )then
              foundsps = .true.
              mapsps(i) = isps+nhsps(-it)
              exit
          end if

        end do

      else
        do isps = 1,nhsps(-it)
          if( mtmp == hspsqn(-it,isps)%m .and. jtmp == hspsqn(-it,isps)%j  & 
          .and. ntmp == hspsqn(-it,isps)%nr .and. ltmp == hspsqn(-it,isps)%l )then
              foundsps = .true.
              mapsps(i) = isps
              exit
          end if

        end do

      end if
      if(.not.foundsps)then   ! ERROR TRAP  
         print*,' DID NOT FIND CORRESPONDING STATE '
         print*,i,ntmp,ltmp,jtmp,mtmp
         stop


      end if

   end do
!................. READ PAST SINGLE-PARTICLE ENERGIES  (added in V7.2.8)....

   do a = 1,nspslocal
      read(23,*)vxx
   end do

!................ READ THROUGH LIST OF MATRIX ELEMENTS.........

   read(23,*)ntbme
   print*,ntbme
   do i = 1,ntbme
      read(23,*)a,b,c,d,vxx
!      print*,a,b,c,d,vxx
      isps = mapsps(a)
      jsps = mapsps(b)
      ksps = mapsps(c)
      lsps = mapsps(d)
 
!....... CHECK ISOSPIN........
      if( isps == -1)cycle
      if( jsps == -1)cycle
      if( ksps == -1)cycle
      if( lsps == -1)cycle

!........MAP TO PAIRS........
      if(  ksps < lsps )then
            dpair = ksps + (lsps-1)*(lsps-2)/2
            dphase = -1
      else
            dpair = lsps + (ksps-1)*(ksps-2)/2
            dphase = 1
      endif
      if(  isps < jsps )then
            cpair = isps + (jsps-1)*(jsps-2)/2
            cphase = -1
      else
            cpair = jsps + (isps-1)*(isps-2)/2
            cphase = 1
      endif

      dpair = map(dpair)
      cpair = map(cpair)
      if(cpair == dpair)then
              hfactor = 0.5
      else
              hfactor = 1.0
      endif
      mtot = (mapm(a)+mapm(b))/2

      if(mtot> jzmax .or. mtot < jzmin)cycle

      if( mtot /= (mapm(c)+mapm(d))/2)then  ! ERROR TRAP
          print*,' MISMATCH IN M-VALUES '
          print*,a,b,c,d
          print*,mtot
          stop
      end if
      partot = mappar(a)*mappar(b)
      if( partot /= mappar(c)*mappar(d))then  ! ERROR TRAP
          print*,' MISMATCH IN PAR-VALUES '
          stop
      end if
      partot = (3-partot)/2

!....... FIND INDEX .............      
      iref = XX2(it)%meref(mtot,partot)
      istart = XX2(it)%mestart(mtot,partot)
      if(cpair <= dpair) then
             itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref
      else
             itbme = istart + (cpair-iref)*(cpair-iref-1)/2+dpair-iref

      endif
      v(itbme) = vxx*dphase*cphase*hfactor
   end do  ! i
   
   close(23)
   deallocate(mapsps,mapm,mappar )
!   print*,hmatnn
   return
101 continue
   print*,' file deformed.int does not exist '
   foundit = .false.
   return

   end subroutine readXXdeformedtbmes
!==========================================================
!
!  CALLED BY:
!

   subroutine readPNdeformedtbmes

   use verbosity
   use spstate
   use haiku_info
   use system_parameters
   use interaction

   implicit none
   integer it      ! which species   
   logical foundit
   integer ttmp, mtmp,jtmp,ntmp,ltmp
   integer mtot,partot
   integer, allocatable,target :: mapspsP(:),mapmP(:),mapparP(:)
   integer, allocatable,target :: mapspsN(:),mapmN(:),mapparN(:)
   integer, pointer :: mapsps(:),mapm(:),mappar(:)

   integer :: nspslocal
   integer :: i
   integer :: ntbme  ! # of uncoupled TBMEs in file
   integer :: a,b,c,d
   integer :: isps,jsps,ksps,lsps
   integer :: cpair,dpair
   real :: hfactor
   real :: vxx
   logical :: foundsps
   integer :: iref, istart
   integer :: cphase,dphase
   integer(8) :: itbme,itbmetr
   real,pointer :: v(:)
   integer,pointer :: map(:)
   integer :: jzmin,jzmax
   integer :: aerr

   foundit = .false.
   open(unit=23,file='deformed.int',status='old',err=101)
   foundit = .true.

!...... FIRST READ PAST AND GET MAPPING OF SINGLE PARTICLE STATES.....
   
   if(.not.allocated(hmatpn)) then
      allocate(hmatpn(nmatpn), stat=aerr)
      if(aerr /= 0) call memerror("readPNdeformedtbmes 1")
   end if
   v => hmatpn

   jzmin = PN2%pair(1)%m      !
   jzmax = PN2%pair(npairpn)%m 

   read(23,*)nspslocal
   allocate(mapspsP(nspslocal),mapmP(nspslocal),mapparP(nspslocal), stat=aerr )
   if(aerr /= 0) then
      call memerror("readPNdeformedtbmes 2")
      stop 5 ! -Wuninitialized doesn't know that memerror won't return
   end if
   allocate(mapspsN(nspslocal),mapmN(nspslocal),mapparN(nspslocal), stat=aerr )
   if(aerr /= 0) then
      call memerror("readPNdeformedtbmes 3")
      stop 5 ! -Wuninitialized doesn't know that memerror won't return
   end if

   mapspsP(:) = -1
   mapspsN(:) = -1
!---- FIND PROTON/NEUTRON STATES
   do it = 1,2
   if(it==1)then
      mapsps=> mapspsP
      mapm  => mapmP
      mappar=> mapparP
   else
      rewind(23)
      read(23,*)nspslocal
      mapsps=> mapspsN
      mapm  => mapmN
      mappar=> mapparN
   endif
   do i = 1,nspslocal
      read(23,*)ttmp, ntmp, ltmp,jtmp,mtmp
      if(ttmp /= it)then
         mapsps(i) = 0
         cycle
      end if
      mapm(i) = mtmp
      mappar(i) = (-1)**ltmp
      foundsps = .false.
!............ SEARCH FOR SINGLE-PARTICLE STATES.....
      if(mtmp > 0)then
        do isps = 1,nhsps(it)
          if( mtmp == hspsqn(it,isps)%m .and. jtmp == hspsqn(it,isps)%j  & 
          .and. ntmp == hspsqn(it,isps)%nr .and. ltmp == hspsqn(it,isps)%l )then
              foundsps = .true.
              mapsps(i) = isps
              exit
          end if

        end do

      else
        do isps = 1,nhsps(-it)
          if( mtmp == hspsqn(-it,isps)%m .and. jtmp == hspsqn(-it,isps)%j  & 
          .and. ntmp == hspsqn(-it,isps)%nr .and. ltmp == hspsqn(-it,isps)%l )then
              foundsps = .true.
              mapsps(i) = -isps
              exit
          end if

        end do

      end if
      if(.not.foundsps)then   ! ERROR TRAP  
         print*,' DID NOT FIND CORRESPONDING STATE '
         print*,i,ntmp,ltmp,jtmp,mtmp
         stop


      end if

   end do
   end do  ! it
!................. READ PAST SINGLE-PARTICLE ENERGIES  (added in V7.2.8)....

   do a = 1,nspslocal
      read(23,*)vxx
   end do

!................ READ THROUGH LIST OF MATRIX ELEMENTS.........

   read(23,*)ntbme
   do i = 1,ntbme
      read(23,*)a,b,c,d,vxx
      isps = mapspsP(a)
      jsps = mapspsN(b)
      ksps = mapspsP(c)
      lsps = mapspsN(d)
!      print*,a,b,c,d,vxx
!....... CHECK ISOSPIN; ASSUME ORDERING pnpn
      if( isps == 0)cycle
      if( jsps == 0)cycle
      if( ksps == 0)cycle
      if( lsps == 0)cycle

!....... ERROR TRAPS.....................

      if( mapmP(a) + mapmN(b) /= mapmP(c)+mapmN(d))then
         print*,' PN problem with m values '
         print*,a,b,c,d
         print*,mapmP(a),mapmN(b), mapmP(c), mapmN(d)
         stop
      end if

      if( mapparP(a) * mapparN(b) /= mapparP(c)*mapparN(d))then
         print*,' PN problem with par values '
         stop
      end if

!........MAP TO PAIRS........


      itbme = cpnpair(jsps,isps)+dpnpair(lsps,ksps)
      itbmetr = dpnpair(jsps,isps)+cpnpair(lsps,ksps)

      if(itbme/=0)v(itbme) = vxx
      if(itbmetr/=itbme .and. itbmetr/=0)v(itbmetr) = vxx
   end do  ! i
   
   close(23)
   deallocate(mapspsP,mapspsN,mapmP,mapmN,mapparP,mapparN )
   print*,' PN DONE '
   return
101 continue
   print*,' file deformed.int does not exist '
   foundit = .false.
   return

   end subroutine readPNdeformedtbmes
!==========================================================



