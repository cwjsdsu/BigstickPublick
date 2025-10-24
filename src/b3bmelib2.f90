

!==================================================================
!
!  file B3BMELIB2.f90
!
!  routines to set up to handle and store 3-body matrix elements
!  focus on XXY (ppn and pnn) 
!
!  initiated Feb 26 2010 CWJ @ SDSU
!
! ---------------------------------------------------------------
!  Goals of these routines:
!	figure out how many XXY (ppn/pnn) matrix elements there are
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! subroutine XXYsetup - master subroutine to create XXY (ppn/pnn) triplets
! CALLS: count_create_pairs
!         master_q_ntuple
!         count_create_XXY
!         count_XYZ_3bmes
!         assign_XXY_maps
!
! subroutine count_create_XXY
!
! subroutine count_XYZ_3bmes
!
! subroutine assign_XXY_maps
!  
!
module b3bmelib2_mod
contains
!================================================================
! subroutine XXYsetup
! master subroutine to create XXY (ppn/pnn) triplets
!
! first creates XX pairs and then adds in a Y particle
!
!  calls: count_create_pairs
!         master_q_ntuple
!         count_create_XXY
!         sort3s
!         count_XYZ_3bmes
!         assign_XXY_maps
!
! CALLED BY threebodysetup
!
  subroutine XXYsetupOrig(itx)

  use spstate
  use system_parameters
  use haiku_info
  use interaction
  use tripdef
  use interactions3body
  use ntuple_info
  use nodeinfo
  use verbosity
  use bmpi_mod
  use jump_mod
  use b3bme_mod
  use btbme_mod
  use ntuple_info
  use adv_tbme_info
  implicit none

  integer itx, ity
  integer maxw2
  integer n3
  integer(kind=8), pointer :: nxxymat
  integer(kind=8) :: ntmp
  integer,pointer :: mapXXY(:,:)
  integer :: aerr
  integer :: i

  ity = 3-itx
  if(itx == 1)then
      nxxymat => nmatppn
  else
      nxxymat => nmatpnn
  endif
  
  if( (Np(itx)-1) *Np(ity) < 1)return

!---------- GET XX PAIRS
  call master_q_ntuple(2, itx)
  maxw2 = pairlim(itx)%maxw

  call count_create_pairs(itx,maxw2,.false.,npairXX(itx),XX2(itx)%pair)
  allocate( XX2(itx)%pair( npairXX(itx) ), stat=aerr)
  if(aerr /= 0) call memerror("count_create_XXY 1")

  call count_create_pairs(itx,maxw2,.true.,npairXX(itx),XX2(itx)%pair)
!----------- SORT XX PAIRS

  if(npairXX(itx)==0)then
     if( itx == 1) then
         nmatppn = 0
     else
         nmatpnn = 0
     endif
     
     return
  endif
  call sortpair(npairXX(itx),XX2(itx)%pair)
  if(sort_XXpairs_on_W)call sortXXpairsW(npairXX(itx),XX2(itx)%pair,.false.)
  
!----------- ADD EXTRA Y and find limits
  n3 = 1
  allocate( XXY(itx)%trip(n3), stat=aerr)
  if(aerr /= 0) call memerror("count_create_XXY 2")
  call count_create_XXY(itx,.false.,npairXX(itx),XX2(itx)%pair,n3, XXY(itx)%trip)
  if(itx==1)then
     if(iproc==0)print*,' There are ',n3,' PPN triplets '
  else
     if(iproc==0)print*,' There are ',n3,' PNN triplets '
  endif
  deallocate(XXY(itx)%trip)
  allocate( XXY(itx)%trip(n3) , stat=aerr)
  if(aerr /= 0) call memerror("count_create_XXY 3")
  nXXY(itx) = n3
  call count_create_XXY(itx,.true.,npairXX(itx),XX2(itx)%pair,n3, XXY(itx)%trip)

!---------- MAP and SORT------
  call sort3s(n3,XXY(itx)%trip)

!----------- COUNT UP UNIQUE MATRIX ELEMENTS ------
  call count_XYZ_3bmes(itx,n3,XXY(itx),nxxymat)
  if(itx == 1)then
!     allocate( mapPPN( npairXX(itx),-nhsps(-ity) : nhsps(ity)  ), stat=aerr)
!     if(aerr /= 0) call memerror("count_create_XXY 5")
!     mapXXY => mapPPN
     allocate( cPPNtriplet( npairXX(itx), -nhsps(-ity): nhsps(ity)) , stat=aerr)
     if(aerr /= 0) call memerror("count_create_XXY 6")
     allocate( dPPNtriplet( npairXX(itx), -nhsps(-ity): nhsps(ity)) , stat=aerr)
     if(aerr /= 0) call memerror("count_create_XXY 7")
     if(useTR) then
        allocate(phaseppn( npairXX(itx), -nhsps(-ity): nhsps(ity)), stat=aerr )  
        if(aerr /= 0) call memerror("count_create_XXY 8")
     end if

!............... DO REVERSE TO MINIMIZE CACHE
!     allocate( cPPNtripletC(  -nhsps(-ity): nhsps(ity), npairXX(itx)) , stat=aerr)
!     if(aerr /= 0) call memerror("count_create_XXY 20")
!     allocate( dPPNtripletC(  -nhsps(-ity): nhsps(ity), npairXX(itx) ) , stat=aerr)
!     if(aerr /= 0) call memerror("count_create_XXY 21")

     if(verbose_uncouple .and. iproc == 0)then
        print*,' '
        print*,' dimensions of arrays cPPNtriplet, dPPNtriplet is '
        print*,npairXX(itx),' by ',nhsps(ity)+nhsps(-ity)+1
        print*,' '
     end if
  else
!     allocate( mapPNN( npairXX(itx), -nhsps(-ity): nhsps(ity) ) , stat=aerr)
!     if(aerr /= 0) call memerror("count_create_XXY 22")
!     mapXXY => mapPNN

     allocate( cPNNtriplet( npairXX(itx), -nhsps(-ity): nhsps(ity)) , stat=aerr)
     if(aerr /= 0) call memerror("count_create_XXY 30")
     allocate( dPNNtriplet( npairXX(itx), -nhsps(-ity): nhsps(ity)) , stat=aerr)
     if(aerr /= 0) call memerror("count_create_XXY 31")
!............... DO REVERSE TO MINIMIZE CACHE
!     allocate( cPNNtripletC(  -nhsps(-ity): nhsps(ity), npairXX(itx)) , stat=aerr)
!     if(aerr /= 0) call memerror("count_create_XXY 40")
!     allocate( dPNNtripletC(  -nhsps(-ity): nhsps(ity), npairXX(itx) ) , stat=aerr)
!     if(aerr /= 0) call memerror("count_create_XXY 41")

     if(useTR) then
        allocate(phasepnn( npairXX(itx), -nhsps(-ity): nhsps(ity)), stat=aerr )  
        if(aerr /= 0) call memerror("count_create_XXY 31")
     end if
     if(verbose_uncouple .and. iproc == 0)then
        print*,' '
        print*,' dimensions of arrays cPNNtriplet, dPNNtriplet is '
        print*,npairXX(itx),' by ',nhsps(ity)+nhsps(-ity)+1
        print*,' '
     end if
  endif
  call assign_XXY_maps(itx,n3,XXY(itx),nxxymat)
  return
  end subroutine XXYsetupOrig

!===================================================================
!
!  Count up and fill XXY triplets by reading in XX pairs and adding Y odd particles
!
! CALLED BY : XXYsetup
!
  subroutine count_create_XXY(itx,fill,nXXpair,pairlist,n3,trips)

  use haiku_info

  use spstate
  use pairdef
  use tripdef
  use ntuple_info
  use system_parameters
  use interactions3body
  use interaction

  implicit none
  integer itx, ity
  logical fill
  integer n3
  integer nXXpair
  type (pair_qn), pointer :: pairlist(:)
  type (triplet_qn) :: trips(n3) 
  type (tuple), pointer :: xxylim
  integer ipair
  integer wpair, parpair,mpair
  integer jsps,j,jth,jsgn
  integer w3,par3,m3
  integer parmult
  integer maxwxy,maxwxxy
  integer wy
  integer isps,i,ith,isgn,wx

  if( .not. restrict_ntuples)then
    print*,' I cannot set XXY without restricting ntuples '
    stop
  endif
  m3 = 0

  if(restrict_ntuples)then
       if(itx ==1) xxylim => ppnlim
       if(itx ==2) xxylim => pnnlim
  endif

  ity = 3 -itx
  n3 = 0
  call findmaxwxy(maxwxy)
  call findmaxwxxy(itx,maxwxxy)

!----- LOOP OVER PAIRS ----

  if(fill)trips(:)%indx = 0   ! set for error trap

  do ipair = 1,nXXpair
     wpair = pairlist(ipair)%w
     mpair = pairlist(ipair)%m*2
     parpair=pairlist(ipair)%par
     parpair = 3-parpair*2

!-------- LOOP OVER Y PARTICLES

    do jsps = 1, nhsps(ity)+nhsps(-ity)
         if(jsps > nhsps(-ity))then
            jth = ity
            j = jsps - nhsps(-ity)
            jsgn = 1
         else
            jth = -ity
            j = jsps
            jsgn = -1
         endif       
       w3 = wpair + hspsqn(jth,j)%w 
       if(w3 > maxwxxy)cycle
       wy =         hspsqn(jth,j)%w
!------------ POSSIBLY USE RESTRICTIONS HERE -------------
       
       if(restrict_ntuples)then
               if(w3 > xxylim%maxw .or. w3 < xxylim%minw)cycle
               par3 =  parpair*hspsqn(jth,j)%par
               par3 = (3-par3)/2

               if( par3 > xxylim%w(w3)%parmax) cycle
               if( par3 < xxylim%w(w3)%parmin) cycle
               m3 = mpair +hspsqn(jth,j)%m
               if( m3 > xxylim%w(w3)%par(par3)%jzmax) cycle
               if( m3 < xxylim%w(w3)%par(par3)%jzmin) cycle

       endif
!------------- further cuts-------------------------------

       i = pairlist(ipair)%ia
       if(i < 0)then
          i = -i
          ith = -itx
       else
          ith = itx
       endif
       wx = hspsqn(ith,i)%w
       if(wx+wy > maxwxy)cycle
       i = pairlist(ipair)%ib
       if(i < 0)then
          i = -i
          ith = -itx
       else
          ith = itx
       endif
       wx = hspsqn(ith,i)%w
       if(wx+wy > maxwxy)cycle

!--------- OPTIONAL USE OF TIME-REVERSAL TO COMPRESS STORAGE--
       if(useTR .and. Jz >=0 .and. m3 < 0 )cycle
       if(useTR .and. Jz <0 .and. m3 > 0 )cycle

       n3 = n3 + 1

       if(fill)then
            trips(n3)%m = mpair + hspsqn(jth,j)%m
            trips(n3)%w = wpair + hspsqn(jth,j)%w

            trips(n3)%par = parpair*hspsqn(jth,j)%par
            trips(n3)%par = (3 - trips(n3)%par)/2


!----------- FILL TRIPLETS---------------------------

            trips(n3)%ia = ipair 
            trips(n3)%ib = j*jsgn
            if(trips(n3)%indx > 0)then  ! error trap
		print*,' Whoops, some problem with triplets '
                stop
            endif
            trips(n3)%indx = jsps + (ipair -1)*( nhsps(ity)+nhsps(-ity))
       end if 

    end do ! jsps
  
  end do ! ipair

  return
  end subroutine count_create_XXY

!====================================================================

!
!  count up how many 3-body matrix elements of type XYZ
!  NOTE: need to restrict on change on W later
!
! used by PPN/PNN
!
! CALLED BY XXYsetup
!
   subroutine count_XYZ_3bmes(it,n3s,XYZ,nxyzmat)
   use tripdef
   use interactions3body
   use nodeinfo
   implicit none

   integer it
   integer n3s
   type (tripinfo) :: XYZ
   integer(kind=8):: nxyzmat

   integer m, mstart,mend
   integer par,parstart,parend
   integer istart,iend,ifinish
   integer i
   nxyzmat = 0
   ifinish = -1
   iend = -1

   mstart = XYZ%trip(1)%m
   mend   = XYZ%trip(n3s)%m

   istart  = 1
   do m = mstart,mend, 2   
      if(XYZ%trip(istart)%m /= m)then
           print*,' mismatch ',m,XYZ%trip(istart)%m
           stop
      endif
!......given m, search for start,finish of parity
      if(istart == n3s)then
         ifinish = n3s
      else
         do i = istart+1,n3s
            ifinish = i
            if(XYZ%trip(i)%m /= m)then
	       ifinish = i-1
               exit
            endif
         enddo
      endif
      parstart = XYZ%trip(istart)%par
      parend   = XYZ%trip(ifinish)%par
      
      do par = parstart,parend
!............ NOW FIND ACTUAL START, FINISH FOR THIS M, PARITY
         if(istart > ifinish)cycle
         if(istart == ifinish)then
	    iend = ifinish
          else
            do i = istart+1,ifinish
               iend = i
               if(XYZ%trip(i)%par /= par)then
		  iend = i-1
                  exit
               endif
            enddo
         endif
!.......... COMPUTE # OF 3BMES
!  NOTE: EVENTUALLY WOULD LIKE TRIANGLE STORAGE but this will require careful sorting
!         nxyzmat = nxyzmat + (iend-istart+1)*(iend-istart+2)/2  ! TRIANGLE STORAGE
         nxyzmat = nxyzmat + (iend-istart+1)*(iend-istart+1)  ! RECTANGLE STORAGE

!........... SET UP FOR NEXT ROUND
         istart = iend + 1
      end do  ! par

   end do  ! m
   if(iproc==0)print*, 'There are ',nxyzmat, ' 3-body matrix elements '

   return
   end subroutine count_XYZ_3bmes
!==================================================================

!
!  find assignments for mapping of XXY 3-body MEs
!  in particular assignment of CPPNTRIPLET, DPPNTRIPLET etc.
!
!  modeled on subroutine COUNT_XYZ_3BMES in file b3bme.f90
!
!  IMPORTANT: This will only work if triplets are sorted to reflect ordering of jumps.
!
! used by PPN/PNN
!
!  called by:
!    XXYsetupOrig
!
!  subroutines called:
!    findTRphaseXX
!
   subroutine assign_XXY_maps(itx,n3s,XYZ,nxyzmat)
   use tripdef
   use interactions3body
   use interaction
   use spstate
   use TRstuff
   implicit none

   integer itx,ity
   integer n3s
   type (tripinfo) :: XYZ
   integer(kind=8):: nxyzmat,nmat0
   integer, pointer :: mapXXY(:,:)
   integer,pointer :: cXXYtriplet(:,:), dXXYtriplet(:,:)
!....... TO MINIMIZE CACHE CALLS
   integer,pointer :: cXXYtripletC(:,:), dXXYtripletC(:,:)

   integer(1),pointer :: phaseXXY(:,:)
   integer(1), target :: phaseXXY_dummy(1,1)
   integer m, mstart,mend
   integer par,parstart,parend
   integer istart,iend,ifinish
   integer i
   integer cXXpair,dXXpair
   integer cYsps, dYsps
   integer jy,y,yth,ysgn,my
   integer cYspsTR, cXXpairTR  
   integer ctrip,dtrip
   integer :: aerr

   phaseXXY => phaseXXY_dummy
   ity = 3-itx
   nxyzmat = 0

   mstart = XYZ%trip(1)%m
   mend   = XYZ%trip(n3s)%m
   if(.not.associated(XYZ%meref)) then
      allocate( XYZ%meref(mstart:mend,2), stat=aerr)
      if(aerr /= 0) call memerror("assign_XXY_maps 1")
   end if
   if(.not.associated(XYZ%mestart)) then
      allocate( XYZ%mestart(mstart:mend,2), stat=aerr)
      if(aerr /= 0) call memerror("assign_XXY_maps 2")
   end if
   XYZ%meref(mstart,:) = 0
   XYZ%mestart(mstart,:) = 0

   if(itx ==1)then
        cXXYtriplet => cppntriplet
        dXXYtriplet => dppntriplet
!        cXXYtripletC => cppntripletC
!        dXXYtripletC => dppntripletC
        if(useTR)phaseXXY => phaseppn
   else
        cXXYtriplet => cpnntriplet
        dXXYtriplet => dpnntriplet
!        cXXYtripletC => cpnntripletC
!        dXXYtripletC => dpnntripletC
        if(useTR)phaseXXY => phasepnn

   end if

   iend = -1
   istart  = 1
   ifinish = -1
   cXXYtriplet=-1
   dXXYtriplet=-1
   if(useTR)phaseXXY = 1
   if(useTR)call findTRpairsXX(itx)
   do m = mstart,mend, 2   
      if(XYZ%trip(istart)%m /= m)then
           print*,' mismatch ',m,XYZ%trip(istart)%m
           stop
      endif
!......given m, search for start,finish of parity
      if(istart == n3s)then
         ifinish = n3s
      else
         do i = istart+1,n3s
            ifinish = i
            if(XYZ%trip(i)%m /= m)then
               ifinish = i-1
               exit
            endif
         enddo
      endif
      parstart = XYZ%trip(istart)%par
      parend   = XYZ%trip(ifinish)%par
      
      do par = parstart,parend
         XYZ%mestart(m,par) = nxyzmat
         XYZ%meref(m,par) = istart-1
!............ NOW FIND ACTUAL START, FINISH FOR THIS M, PARITY
         if(istart > ifinish)cycle
         if(istart == ifinish)then
            iend = ifinish
          else
            do i = istart+1,ifinish
               iend = i
               if(XYZ%trip(i)%par /= par)then
                  iend = i-1
                  exit
               endif
            enddo
         endif
         do ctrip = istart,iend   ! not sure if I want to start with ctrip or dtrip
!................. EXTRACT S.P. STATE INFO..........
             cYsps = XYZ%trip(ctrip)%ib

!.................. GET XXPAIRS...................
             cXXpair = XYZ%trip(ctrip)%ia
!             mapXXY(cXXpair,cYsps) = ctrip
             cXXYtriplet(cXXpair,cYsps) = nxyzmat
             dXXYtriplet(cXXpair,cYsps) = ctrip-istart+1
			 
!........... TO MINIMIZE CACHE CALLS LATER ON
!             cXXYtripletC(cYsps,cXXpair) = nxyzmat
!             dXXYtripletC(cYsps,cXXpair) = ctrip-istart+1
!................. USING TIME-REVERSAL...........
             if(useTR)then
                 if(cYsps < 0)then
                    yth = - ity
                    y   = -cYsps
                    ysgn= -1
                 else
                    yth =  ity
                    y   =  cYsps
                    ysgn= +1
                 endif
                 jy = hspsqn(yth,y)%j
                 my = hspsqn(yth,y)%m
                 if(my == 0)then
                   cYspsTR = cYsps
                 else
                   cYspsTR = -ysgn*hspsqn(yth,y)%tr
                 endif
                 cXXpairTR = XX2(itx)%pair(cXXpair)%tr
                 if(cXXpairTR > 0 )then
                    cXXYtriplet(cXXpairTR,cYspsTR) = nxyzmat
                    dXXYtriplet(cXXpairTR,cYspsTR) = ctrip-istart+1
!........... MINIMIZE CACHE????
                    phaseXXY(cXXpairTR,cYspsTR) = &
                        int( XX2(itx)%pair(cXXpair)%trphase*(-1)**((jy+1)/2), 1)
                 end if
             end if
             nxyzmat = nxyzmat + (iend-istart)+1


         end do ! ctrip
!........... SET UP FOR NEXT ROUND
         istart = iend + 1
      end do  ! par

   end do  ! m

   return
   end subroutine assign_XXY_maps


!=================================================================
!  subroutine maxwtrip
!  finds the maximum W allowed for a triplet
!  by comparing max W for N, N-3 particles
!
!  
   subroutine findmaxWxy(maxwxy)

   use spstate
   use W_info
   use haiku_info
   use system_parameters

   implicit none
   integer it
   integer maxwxy
   
   integer w
   integer i
   maxwxy = 0
   if(Np(1) < 1 .or. Np(2) < 1)return
   w = 0
   do it = 1,2
   if(np(it) ==1)cycle
   do i =1,Np(it)-1
        w = w + spsqn(it,i)%w
   enddo
   end do
   maxwxy = maxWtot - w
   return
   end subroutine findmaxWxy

!=================================================================
!  subroutine maxwtrip
!  finds the maximum W allowed for a triplet
!  by comparing max W for N, N-3 particles
!
!  
   subroutine findmaxWxxy(itx,maxwxxy)

   use spstate
   use W_info
   use haiku_info
   use system_parameters

   implicit none
   integer itx,ity
   integer maxwxxy
   
   integer w
   integer i
   maxwxxy = 0
   ity = 3 -itx
   if(Np(itx) < 2 .or. Np(ity) < 1)return
   w = 0
   do i =1,Np(itx)-2
        w = w + spsqn(itx,i)%w
   enddo
   do i =1,Np(ity)-1
        w = w + spsqn(ity,i)%w
   enddo
   maxwxxy = maxWtot - w
   return
   end subroutine findmaxWxxy


!================================================================
end module b3bmelib2_mod

! This routine is here to break a module dependency loop
subroutine XXYsetup(itx)
   use b3bmelib2_mod
   implicit none
   integer :: itx

   call XXYsetupOrig(itx)
end subroutine XXYsetup
