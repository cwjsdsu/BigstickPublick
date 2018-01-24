!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK CI shell-model code
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
module welder

contains	
!
!
!====================================================================
!  subroutine weld
!
! joins chains into jumps
!
!  INPUT:
!    it: = species (1= protons, 2 = neutrons)
!    is,fs: initial, final sector
!    nbody = 1 or 2 body
!    isj = index of sector jump
!    hflag : Hermiticity flag; if true then can set other flags
!    fill  : if true, then fill jump arrays; else just count
!    xjump  : data structure for jump control
!       .. ADDED IN 7.5.3 FOR COMPACT JUMPS STORAGE....
!    jumpstart, jumpend : start/stop for storing jumps
!    jumpshift          : shift in index for storing jumps
!  OUTPUT:
!    njumps  ! updated
!

!====================================================================
   subroutine weld(it,is,fs,nbody,isj,hflag,fill,xjump,jumpstart,jumpend,jumpshift, njumps)

   use verbosity
   use haiku_info
   use spstate
   use blocks
   use sectors
   use chains
   use jumpNbody
   use flags3body
   use jump3body
   use jumpdef
   use interaction
   use pairdef
   use precisions
   use basis
   use sporbit
   use flags3body
   use nodeinfo
   use jumplimits
   use flagger
   use butil_mod
   use btbme_mod

   implicit none
   integer it  ! which species
   integer is,fs
   integer nbody
   integer isj   ! which sector jump
   logical hflag
   logical fill   ! flag to just count up
   type (jumper) :: xjump
   integer(8) :: jumpstart,jumpend,jumpshift

!.........................
   integer(8):: njumps

  type (handchain),pointer :: mrchain
  type (chainbase), pointer :: lcurch,rcurch
  logical rdiag, ldiag
  integer igen
  logical usediag

  integer :: iblock,fblock
  integer(8) :: isdstart,fsdstart
  integer(8) :: isd,fsd
  integer(8) :: iradd,fradd,iladd,fladd
  integer(8) :: niradd,nfradd

  integer :: irblock,ilblock,frblock,flblock
  integer :: lphase,rphase

  integer :: opleft(maxgen),opright(maxgen),opsum(maxgen)

!--------------- NEEDED FOR FINDING INDICES OF OPERATORS ----------

  logical first
  integer :: m, par
  integer :: i,j,k,l
  integer(kind=8) :: isps,jsps,ksps,lsps
  integer :: ith,jth,kth,lth
  integer(kind=8) :: iref,istart
  integer(kind=8) :: cpair,dpair
  integer(kind=8) :: itbme
  integer :: cphase,dphase
  integer(kind=1) :: genphase
  integer(8) :: nzero
  real, pointer :: v(:)
  integer(8), pointer :: map(:)
  integer(8), target :: map_dummy(1)

  integer(kind=basis_prec),pointer :: xisd(:),xfsd(:)
  integer(kind=basis_prec), target :: xisd_dummy(1), xfsd_dummy(1)
  integer(kind=basis_prec),pointer :: xisd0(:),xfsd0(:)
  integer(kind=basis_prec), target :: xisd0_dummy(1), xfsd0_dummy(1)
  integer(kind = 2),pointer :: cop(:),dop(:)
  integer(kind = 2), target :: cop_dummy(1), dop_dummy(1)
  integer(4), pointer :: cop2b(:),dop2b(:)
  integer(4), target :: dop2b_dummy(1), cop2b_dummy(1)
  integer, pointer :: op(:)
  integer, target :: op_dummy(1)
  real(4),pointer :: me(:)
  real(4), target :: me_dummy(1)
  real :: xxme
  integer(kind = 1),pointer :: phase(:)
  integer(kind = 1), target :: phase_dummy(1)
  integer(kind=basis_prec), pointer :: xxstart(:)
  integer(kind=basis_prec), target :: xxstart_dummy(1)

!-------------- FOR SKIPPING UNUSED TWO-BODY JUMPS
!              because of failing the triangle rule

  integer jabmin,jabmax

  logical :: diagnose = .false.

  integer(8):: minjump,maxjump

  integer(4) :: iprocs

!
! init and use of par,m both conditioned by same condition
  par = -100000   ! prevent uninitialized msg
  m = -200000     ! prevent uninitialized msg
  itbme = -400000 ! prevent uninitialized msg
  ! failure to initialize pointers leads to many lines of msgs
  ! Should add some dynamic checks that pointers are associated.
  op => op_dummy
  me => me_dummy
  dop2b => dop2b_dummy
  cop2b => cop2b_dummy
  map => map_dummy
  xisd => xisd_dummy
  xfsd => xfsd_dummy
  xxstart => xxstart_dummy
  xisd0 => xisd0_dummy
  xfsd0 => xfsd0_dummy
  phase => phase_dummy
  cop => cop_dummy
  dop => dop_dummy
  minjump = -1
  maxjump = -1
  dphase = 0  ! bad value intentially, should be overwritten
  cphase = 0  ! bad value intentially, should be overwritten

!.......... SET UP LIMITS FOR RESTRICTING JUMPS....
!           MAY BE DUMMY LIMITS

   iprocs = iproc

  if((restrictjumps .and. nproc > 1 .and. fill) )then
	  if(compactjumpstorage)then
          if(nbody == 1 .and. it ==1)then
             if(.not.makep1bjumps(iprocs))return
			 minjump = jumpstart
			 maxjump = jumpend 
			 if(minjump-jumpshift <  startp1bjumps(iprocs)) then
				 print*,iproc,' Error in starting point for P jumps '
				 print*,minjump,maxjump,jumpshift, startp1bjumps(iprocs)
				 stop
			 end if
			 if(maxjump-jumpshift > stopp1bjumps(iprocs)) then
				 print*,iproc,' Error in stopping point for P jumps '
				 print*,minjump,maxjump,jumpshift, stopp1bjumps(iprocs)
				 stop
			 end if			 
          end if
          if(nbody == 1 .and. it ==2)then 
             if(.not.maken1bjumps(iprocs))return
 			 minjump = jumpstart
 			 maxjump = jumpend 
 			 if(minjump-jumpshift <  startn1bjumps(iprocs)) then
 				 print*,iproc,' Error in starting point for N jumps '
 				 print*,minjump,jumpshift, startn1bjumps(iprocs)
 				 stop
 			 end if
 			 if(maxjump-jumpshift > stopn1bjumps(iprocs)) then
 				 print*,iproc,' Error in stopping point for N jumps '
 				 print*,maxjump,jumpshift, stopn1bjumps(iprocs)
 				 stop
 			 end if					  
          end if
          if(nbody == 2 .and. it ==1)then
             if(.not.makeppjumps(iprocs))return
			 minjump = jumpstart
 			 maxjump = jumpend 
 			 if(minjump-jumpshift <  startppjumps(iprocs)) then
 				 print*,iproc,' Error in starting point for PP jumps '
 				 print*,minjump,jumpshift, startppjumps(iprocs)
 				 stop
 			 end if
 			 if(maxjump-jumpshift > stopppjumps(iprocs)) then
 				 print*,iproc,' Error in stopping point for PP jumps '
 				 print*,maxjump,jumpshift, stopppjumps(iprocs)
 				 stop
 			 end if		
          end if
          if(nbody == 2 .and. it ==2)then
             if(.not.makennjumps(iprocs))return
 			 minjump = jumpstart
 			 maxjump = jumpend 
 			 if(minjump-jumpshift <  startnnjumps(iprocs)) then
 				 print*,iproc,' Error in starting point for NN jumps '
 				 print*,minjump,jumpshift, startnnjumps(iprocs)
 				 stop
 			 end if
 			 if(maxjump-jumpshift > stopnnjumps(iprocs)) then
 				 print*,iproc,' Error in stopping point for NN jumps '
 				 print*,maxjump,jumpshift, stopnnjumps(iprocs)
 				 stop
 			 end if		
          end if
	  else

         if(nbody == 1 .and. it ==1)then
            if(.not.makep1bjumps(iprocs))return
            minjump = startp1bjumps(iprocs)
            maxjump = stopp1bjumps(iprocs)
         end if
         if(nbody == 1 .and. it ==2)then
            if(.not.maken1bjumps(iprocs))return
            minjump = startn1bjumps(iprocs)
            maxjump = stopn1bjumps(iprocs)
         end if
         if(nbody == 2 .and. it ==1)then
             if(.not.makeppjumps(iprocs))return
             minjump = startppjumps(iprocs)
             maxjump = stopppjumps(iprocs)
         end if
         if(nbody == 2 .and. it ==2)then
             if(.not.makennjumps(iprocs))return
             minjump = startnnjumps(iprocs)
             maxjump = stopnnjumps(iprocs)
         end if
         if(njumps+xjump%njumps < minjump)return
		 if(jumpshift/=0)then             ! ERROR TRAP
			 print*,iproc,' problem with jumps shift when not expecting ',jumpshift
			 stop
		 end if                           ! END ERROR TRAP
	  end if
  else
     minjump = 1
     if(nbody == 1 .and. it ==1)then
         maxjump = totp1bjumps
     end if
     if(nbody == 1 .and. it ==2)then
         maxjump = totn1bjumps
     end if
     if(nbody == 2 .and. it ==1)then
         maxjump = totp2bjumps
     end if
     if(nbody == 2 .and. it ==2)then
         maxjump = totn2bjumps
     end if
  end if
!  if(it==1 .and. fill)print*,iproc,' weld limits ',maxjump,stopppjumps(iproc)

  nzero = 0              ! count up how many zero matrix elements for 2-body
                        ! this speeds up application of 2-body matrix elements

  if((fill ) .and. xjump%njumps == 0)return
  if(nbody == 2)then
     if(it ==1 )then
        map => mappairPP
        v   => hmatpp
     else
        map => mappairNN
        v   => hmatnn
     endif
  endif

  if(njumps < 0)then
      print*,' oops neg jumsp ',njumps
      stop
  endif
  if(nbody == 1)then
  if(it == 1)then

    xisd => p1b_isd
    xfsd => p1b_fsd
    xisd0 => p1b_isd0    ! used for printing out jumps
    xfsd0 => p1b_fsd0
    cop => p1b_cop
    dop => p1b_dop
    phase=>p1b_phase
    xxstart=> pstart
  else

    xisd => n1b_isd
    xfsd => n1b_fsd
    xisd0 => n1b_isd0    ! used for printing out jumps
    xfsd0 => n1b_fsd0
    cop => n1b_cop
    dop => n1b_dop
    phase=>n1b_phase
    xxstart => nstart
  endif
  endif

  if(nbody == 2)then
  if(it == 1)then

    xisd => p2b_isd
    xfsd => p2b_fsd
    xisd0 => p2b_isd0    ! used for printing out jumps
    xfsd0 => p2b_fsd0
    phase=>p2b_phase
    xxstart=> pstart
    if(threebody)then
    dop2b => p2b_dop
    cop2b => p2b_cop
    else
       if(storeXXmesjumps)then
          me => p2b_me
       else
          op => p2b_op
       end if
    endif

  else
    xisd => n2b_isd
    xfsd => n2b_fsd
    xisd0 => n2b_isd0    ! used for printing out jumps
    xfsd0 => n2b_fsd0
!    me => n2b_me

    phase=>n2b_phase
    xxstart => nstart
    if(threebody)then
    dop2b => n2b_dop
    cop2b => n2b_cop
    else
       if(storeXXmesjumps)then
          me => n2b_me
       else
          op => n2b_op
       end if
    endif
  endif
  endif
  mrchain => drchain

     first = .true.
     iblock = mrchain%iblock
     fblock = mrchain%fblock
     genphase = mrchain%genphase
     ilblock = xsd(it)%sector(is)%lhblock(iblock)
     flblock = xsd(it)%sector(fs)%lhblock(fblock)
     irblock = xsd(it)%sector(is)%rhblock(iblock)
     frblock = xsd(it)%sector(fs)%rhblock(fblock)

     if(ilblock == flblock .and. irblock == frblock .and. hflag)then
! IMPORTANT FOR HERMITICITY
        usediag = .true.           
     else
        usediag = .false.
      endif

    lcurch => mrchain%lchain
!-----------------------------------      
    niradd = hblock(it)%list(irblock)%nhsd
    nfradd = hblock(it)%list(frblock)%nhsd

    isdstart =  xsd(it)%sector(is)%blockstart(iblock)-1 
    fsdstart =  xsd(it)%sector(fs)%blockstart(fblock)-1 

    do while(lcurch%chain%iadd /=0 )
      if(lcurch%chain%sterile)then
         lcurch => lcurch%next
         cycle
      endif

      iladd = lcurch%chain%iadd
      fladd = lcurch%chain%fadd
      lphase= lcurch%chain%phase
      opleft = lcurch%chain%oplist
      rcurch => mrchain%rchain

      do while(rcurch%chain%iadd /=0 )

         if(rcurch%chain%sterile)then    
              rcurch => rcurch%next
              cycle
         endif   

!--------------------- CONSTRUCT A JUMP -----------------

         iradd = rcurch%chain%iadd
         fradd = rcurch%chain%fadd
         rphase= rcurch%chain%phase

         isd = isdstart+ (iladd-1)*niradd + iradd
         fsd = fsdstart + (fladd-1)*nfradd+ fradd
   if(diagnose .and. nbody ==2)write(91,*)' (C) ',isd,fsd

         if(usediag .and. isd < fsd)then  ! could be the other way around
					! Important for HERMITICITY
              rcurch => rcurch%next
              cycle  
         endif
!-------------- ERROR TRAP ---------


!---------------- FIGURE OUT PHASE
!        IMPORTANT: RIGHT HAIKUS PICK UP AN EXTRA PHASE  FROM OPERATORS
!        COMMUTING PAST THE LEFT HAIKUS
!        this was set up in routine "shiva" in bhoplib.f90
!
!-----------------FIGURE OUT OPERATORS
         opsum = -opleft + rcurch%chain%oplist
         select case (nbody)
		 
!-------------ONE-BODY  --- a^+i a_j
          case  (1)

          j = opsum(1)
          i = opsum(2)
          if(j < 0)then
            jsps = -j
          else
            jsps = j + nhsps(-it)
          endif

          if(i < 0)then
            isps = -i
          else
            isps = i + nhsps(-it)
          endif

           njumps = njumps + 1
         if(fill .and. njumps-xjump%nstart > xjump%njumps)then
           print*,' error, too many jumps ',nbody,it
           print*,xjump%njumps
           stop
         endif
!	   if(fill .and. iproc==2)print*,iproc,' test filling ',it,njumps-jumpshift,minjump,maxjump

           if( fill .and. njumps >= minjump .and. njumps <=maxjump ) then
!			   if(iproc==1 .and. it==1)write(87,*),iproc,' 1body test filling ',njumps-jumpshift,i,j
              cop(njumps-jumpshift) = int(i, kind(cop(1)))
              dop(njumps-jumpshift) = int(j, kind(dop(1)))
              phase(njumps-jumpshift)= int(lphase*rphase*genphase, kind=1)
           endif
		   
!------------ TWO-BODY ---  a^+_i a^+_j a_l a_k
          case  (2)

          k = opsum(1)
          l = opsum(2)
          j = opsum(3)
          i = opsum(4)
          
!------- OPTION TO SKIP DIAGONALS ------------ ADDED in 7.6.0 ---
          if( sumdiagonalXXmesjumps .and. .not.threebody .and. & 
!		   ( (i == k .and. j==l) .or. (i == l .and. j==k)  ) )then
!		   ( i == k .or. j==l .or. i == l .or. j==k  ) )then
                       isd==fsd)then
               rcurch => rcurch%next
		       cycle
     	   end if
          if(l < 0)then
            lsps = -l
            lth = -it
          else
            lsps = l + nhsps(-it)
            lth = it
          endif

          if(k < 0)then
            ksps = -k
            kth = -it
          else
            ksps = k + nhsps(-it)
            kth = it
          endif

          if(j < 0)then
            jsps = -j
            jth = -it
          else
            jsps = j + nhsps(-it)
            jth = it
          endif

          if(i < 0)then
            isps = -i
            ith = -it
          else
            isps = i + nhsps(-it)
            ith = it
          endif


          if(fill )then   ! compute index to matrix element
!----------- THIS PHASES MAY BE BACKWARDS
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
          if(dpair <= 0)then
              print*,' dpair ',dpair,k,l
              print*,opsum
              stop
          endif
          if(cpair <= 0)then
              print*,'cpair ',cpair,i,j
              print*,opsum

              stop
          endif
          dpair = map(dpair)
          cpair = map(cpair)
!------------------ FIND M, PAR OF OPERATORS

          if(first)then

              if( i < 0 )then
                m = hspsqn(-it,-i)%m
                par=hspsqn(-it,-i)%par
              else
                m = hspsqn(it,i)%m
                par=hspsqn(it,i)%par
              endif
              if(j < 0) then
                 m = m+hspsqn(-it,-j)%m
                par=par*hspsqn(-it,-j)%par
              else
                m = m+hspsqn(it,j)%m
                par=par*hspsqn(it,j)%par
              endif
              m = m/2
              par = (3-par)/2
              first = .false.
           endif
            
           iref = XX2(it)%meref(m,par)
           istart = XX2(it)%mestart(m,par)
           if(cpair <= dpair) then
             itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref
           else
             itbme = istart + (cpair-iref)*(cpair-iref-1)/2+dpair-iref

           endif
           end if ! IF FILL
!------------- CHECK IF TRIANGLE RULE IS SATISFIED---------------
!              BUT NOT IF DOING 3-BODY

           if(skipzeros .and. .not. threebody )then
              jabmax = hspsqn(ith,abs(i))%j + hspsqn(jth,abs(j))%j
              jabmin = abs(hspsqn(ith,abs(i))%j - hspsqn(jth,abs(j))%j)
              if( hspsqn(kth,abs(k))%j + hspsqn(lth,abs(l))%j  < jabmin .or.  & 
              abs(hspsqn(kth,abs(k))%j - hspsqn(lth,abs(l))%j) > jabmax)then
                 rcurch => rcurch%next
                 nzero = nzero+1
                 cycle
              endif
           endif

!------------- CHECK FOR CLEBSCHS THAT VANISH
!              BUT NOT IF DOING 3-BODY
!  10/09   CWJ
!
!  rule: if  a = b, then the state | a  b J T > can only exist if  2ja + J +T is even.
!  for spinful fermions (j = half-integer) this means J + T is odd; for pp/nn this means J is even
!  for spinless fermions (j = integer) this means J+T is even and for pp/nn J is odd.
!
!  This scenario can only happen if i == j or k ==l and there is a unique value of J allowed
!
!  NB--there still seem to be zero matrix elements that I cannot predict.  I leave in these conditionals
!  but for now one must ultimately use the matrix elements themselves to cut on.
!

           if( skipzeros .and. .not. threebody .and. &
                  ( hspsqn(ith,abs(i))%orb == hspsqn(jth,abs(j))%orb .or. hspsqn(kth,abs(k))%orb == hspsqn(lth,abs(l))%orb) ) then

                jabmax = hspsqn(ith,abs(i))%j + hspsqn(jth,abs(j))%j
                jabmax = bmin(jabmax, hspsqn(kth,abs(k))%j + hspsqn(lth,abs(l))%j )

                jabmin = abs(hspsqn(ith,abs(i))%m + hspsqn(jth,abs(j))%m )
                jabmin = bmax(jabmin, abs( hspsqn(ith,abs(i))%j - hspsqn(jth,abs(j))%j ))
                jabmin = bmax(jabmin, abs( hspsqn(kth,abs(k))%j - hspsqn(lth,abs(l))%j ))

                if(jabmax < jabmin)then    ! triangle rule violated
                   rcurch => rcurch%next
                   nzero = nzero+1
                   cycle
                endif

                if(jabmax == jabmin)then
                    jabmax = jabmax/2
                    if( ( spinless .and. (jabmax/2)*2 == jabmax) .or. & 
                        ( .not.spinless .and. (jabmax/2)*2 /= jabmax ) ) then
                           rcurch => rcurch%next
                           nzero = nzero+1
                           cycle
                    
                    endif
                endif

!......CGCs of the form (j m, j m | JM ) = 0 if 2j+J is odd  

                if(  ( .not.spinless )  .and. & 
                   ( (hspsqn(ith,abs(i))%j == hspsqn(jth,abs(j))%j  .and. hspsqn(ith,abs(i))%m == hspsqn(jth,abs(j))%m ) .or. & 
                   (hspsqn(kth,abs(k))%j == hspsqn(lth,abs(l))%j  .and. hspsqn(kth,abs(k))%m == hspsqn(lth,abs(l))%m ) )) then

                      rcurch => rcurch%next
                      nzero = nzero+1
                      cycle
                 endif

           endif
!......CGCs of the form (j m, j m | JM ) = 0 if 2j+J is odd  -- ANOTHER CASE

           if(  ( skipzeros .and. .not.spinless .and. .not. threebody)  .and. & 
                   ( (hspsqn(ith,abs(i))%j == hspsqn(jth,abs(j))%j  .and. hspsqn(ith,abs(i))%m == hspsqn(jth,abs(j))%m ) .or. & 
                   (hspsqn(kth,abs(k))%j == hspsqn(lth,abs(l))%j  .and. hspsqn(kth,abs(k))%m == hspsqn(lth,abs(l))%m ) )) then

                jabmax = hspsqn(ith,abs(i))%j + hspsqn(jth,abs(j))%j
                jabmax = bmin(jabmax, hspsqn(kth,abs(k))%j + hspsqn(lth,abs(l))%j )

                jabmin = abs(hspsqn(ith,abs(i))%m + hspsqn(jth,abs(j))%m )
                jabmin = bmax(jabmin, abs( hspsqn(ith,abs(i))%j - hspsqn(jth,abs(j))%j ))
                jabmin = bmax(jabmin, abs( hspsqn(kth,abs(k))%j - hspsqn(lth,abs(l))%j ))

                if(jabmin == jabmax .and. (jabmin/4)*4 == jabmin)then
                      rcurch => rcurch%next
                      nzero = nzero+1
                      cycle
                endif
           endif
!------------------ OPTION TO SKIP OVER ZERO MATRIX ELEMENT JUMPS ----------

           njumps = njumps + 1
           if(fill .and. njumps-xjump%nstart > xjump%njumps)then
              print*,' error, too many jumps ',nbody,it
              print*,xjump%njumps
              stop
            endif
            if(fill .and. njumps >= minjump .and. njumps <=maxjump )then

                 if(threebody)then
                    dop2b(njumps-jumpshift) = int(dpair, kind(dop2b(1)))  ! placeholder
                    cop2b(njumps-jumpshift) = int(cpair, kind(cop2b(1))) ! placeholder
                 else
                    if(storeXXmesjumps)then   ! must decouple matrix element
                       call fetchXXme(it,i,j,k,l,xxme)
                       me(njumps-jumpshift) = xxme
                    else
                       op(njumps-jumpshift)= bint8to4(itbme, 123) ! downsize with test (loc=123)
                    end if
                 endif
                 phase(njumps-jumpshift)= int(lphase*rphase*cphase*dphase*genphase, kind=1)
            endif

         end select
         if( fill .and. njumps >= minjump .and. njumps <=maxjump )then
!			 if(it==1 )print*,iproc,njumps-jumpshift,stopppjumps(iproc)
            xisd(njumps-jumpshift) = xxstart(isd)     ! nzero is just used for 2-body
            xfsd(njumps-jumpshift) = xxstart(fsd)
            if(print_jumps)then
                xisd0(njumps-jumpshift) = isd
                xfsd0(njumps-jumpshift) = fsd
            endif
         endif
         if(fill .and. njumps > maxjump)return

!-----------------------------------------------------
!1963     continue
         rcurch => rcurch%next
      enddo  
      lcurch => lcurch%next
    enddo 

!------------------------------
    mrchain => mrchain%next

   return
   end subroutine weld
!====================================================================


!================================================================
!  subroutine weld3b
!
! joins chains into 3b-jumps
!
   subroutine weld3b(it,is,fs,nbody,isj,hflag,fill,xjump,jumpstart,jumpend,jumpshift,njumps)

   use haiku_info
   use spstate
   use blocks
   use sectors
   use chains
   use jump3body
   use jumpNbody
   use jumpdef
   use interaction
   use interactions3body
   use pairdef
   use precisions
   use basis
   use jumplimits
   use flagger
   use nodeinfo
   implicit none
   integer it  ! which species
   integer is,fs
   integer nbody
   integer isj   ! which sector jump
   logical hflag,fill
   integer(8) :: jumpstart,jumpend,jumpshift
   type (jumper) :: xjump

!.........................
   integer(8):: njumps

  type (handchain),pointer :: mrchain
  type (chainbase), pointer :: lcurch,rcurch
  logical rdiag, ldiag
  integer igen
  logical usediag

  integer :: iblock,fblock
  integer(8) :: isdstart,fsdstart
  integer(8) :: isd,fsd
  integer :: iradd,fradd,iladd,fladd
  integer :: niradd,nfradd

  integer :: irblock,ilblock,frblock,flblock
  integer :: lphase,rphase

  integer :: opleft(maxgen),opright(maxgen),opsum(maxgen)

!--------------- NEEDED FOR FINDING INDICES OF OPERATORS ----------

  logical first
  integer :: m, par
  integer :: i,j,k,l
  integer :: isps,jsps,ksps,lsps

  integer :: ic,jc,kc,id,jd,kd  ! creation/destruction labels
  integer :: icsps,jcsps,kcsps,idsps,jdsps,kdsps

  integer :: iref,istart
  integer :: cpair,dpair

  integer :: cindx,dindx

  integer :: itbme
  integer :: cphase,dphase
  integer(kind=1) :: genphase
  integer(8) :: nzero
  real, pointer :: v(:)
  integer nphx
  integer, pointer :: map(:)
  type (pair_qn), pointer :: phx(:)

  integer(kind=basis_prec),pointer :: xisd(:),xfsd(:)
  integer(kind = 2),pointer :: cop(:),dop(:)
  integer(kind=4), pointer :: op(:)
  integer(kind = 1),pointer :: phase(:)
  integer(kind=basis_prec), pointer :: xxstart(:)
  integer iprocs 

  integer(8) :: maxnmat
  integer(8):: minjump,maxjump


  if(nbody/=3)then
      print*,' Problem: supposed to be welding 3-body ',nbody
      stop
  end if
  minjump = 1
  maxjump = 1
!.......... SET UP LIMITS FOR RESTRICTING JUMPS....
!           MAY BE DUMMY LIMITS
  iprocs = iproc
  if(restrictjumps .and. nproc > 1 .and. fill)then
	  if(compactjumpstorage)then
          if(it ==1)then
             if(.not.makepppjumps(iprocs))return
			 minjump = jumpstart
			 maxjump = jumpend 
			 if(minjump-jumpshift <  startpppjumps(iproc)) then
				 print*,iproc,' Error in starting point for jumps '
				 print*,minjump,maxjump,jumpshift, startpppjumps(iprocs)
				 stop
			 end if
			 if(maxjump-jumpshift > stoppppjumps(iprocs)) then
				 print*,iproc,' Error in stopping point for jumps '
				 print*,minjump,maxjump,jumpshift, stopp1bjumps(iprocs)
				 stop
			 end if		
		 else
			 
             if(.not.makennnjumps(iprocs))return
			 minjump = jumpstart
			 maxjump = jumpend 
			 if(minjump-jumpshift <  startNNNjumps(iprocs)) then
				 print*,iproc,' Error in starting point for jumps '
				 print*,minjump,maxjump,jumpshift, startNNNjumps(iprocs)
				 stop
			 end if
			 if(maxjump-jumpshift > stopNNNjumps(iprocs)) then
				 print*,iproc,' Error in stopping point for jumps '
				 print*,minjump,maxjump,jumpshift, stopp1bjumps(iprocs)
				 stop
			 end if					 
			 
			 	 
          end if
	  else

         if(it ==1)then
            if(.not.makepppjumps(iproc))return
            minjump = startpppjumps(iproc)
            maxjump = stoppppjumps(iproc)
        end if
        if( it ==2)then
            if(.not.makennnjumps(iproc))return
            minjump = startnnnjumps(iproc)
            maxjump = stopnnnjumps(iproc)
        end if
        if(njumps+xjump%njumps < minjump)return
    end if
  else
     minjump = 1
     if( it ==1)then
         maxjump = totp3bjumps
     end if
     if( it ==2)then
         maxjump = totn3bjumps
     end if
  end if

  nzero = 0              ! count up how many zero matrix elements for 2-body
                        ! this speeds up application of 2-body matrix elements

  if(fill .and. xjump%njumps == 0)return
  if(nbody == 3)then
     if(it ==1 )then
        map => mapPPP
     else
        map => mapNNN
     endif
     if(it == 1)then

        xisd => p3b_isd
        xfsd => p3b_fsd
        op => p3b_op
        phase=>p3b_phase
        xxstart=> pstart

        maxnmat = nmatppp
!       print*,' There are ',nmatppp,' matrix elements here '
     else
        xisd => n3b_isd
        xfsd => n3b_fsd
        op => n3b_op
        phase=>n3b_phase
        xxstart => nstart

        maxnmat = nmatnnn
     endif
  endif

  mrchain => drchain
     first = .true.
     iblock = mrchain%iblock
     fblock = mrchain%fblock
     genphase = mrchain%genphase
     ilblock = xsd(it)%sector(is)%lhblock(iblock)
     flblock = xsd(it)%sector(fs)%lhblock(fblock)
     irblock = xsd(it)%sector(is)%rhblock(iblock)
     frblock = xsd(it)%sector(fs)%rhblock(fblock)

     if(ilblock == flblock .and. irblock == frblock .and. hflag)then
        usediag = .true.           
     else
        usediag = .false.
      endif

    lcurch => mrchain%lchain
!-----------------------------------      
    niradd = hblock(it)%list(irblock)%nhsd
    nfradd = hblock(it)%list(frblock)%nhsd

    isdstart =  xsd(it)%sector(is)%blockstart(iblock)-1 
    fsdstart =  xsd(it)%sector(fs)%blockstart(fblock)-1 

    do while(lcurch%chain%iadd /=0 )
      if(lcurch%chain%sterile)then
         lcurch => lcurch%next
         cycle
      endif

      iladd = lcurch%chain%iadd
      fladd = lcurch%chain%fadd
      lphase= lcurch%chain%phase
      opleft = lcurch%chain%oplist
      rcurch => mrchain%rchain

      do while(rcurch%chain%iadd /=0 )

         if(rcurch%chain%sterile)then    
              rcurch => rcurch%next
              cycle
         endif   

!--------------------- CONSTRUCT A JUMP -----------------

         iradd = rcurch%chain%iadd
         fradd = rcurch%chain%fadd
         rphase= rcurch%chain%phase

         isd = isdstart+ (iladd-1)*niradd + iradd
         fsd = fsdstart + (fladd-1)*nfradd+ fradd

         if(usediag .and. isd < fsd)then  ! could be the other way around
              rcurch => rcurch%next
              cycle  
         endif
         njumps = njumps + 1
!-------------- ERROR TRAP ---------

         if(fill .and. njumps-xjump%nstart > xjump%njumps)then
           print*,' error, too many jumps for 3 body '
           print*,xjump%njumps, njumps-xjump%nstart
           stop
         endif
         if(fill.and. njumps >=minjump .and. njumps <=maxjump)then
           xisd(njumps-jumpshift) = xxstart(isd)     ! nzero is just used for 2-body
           xfsd(njumps-jumpshift) = xxstart(fsd)
         endif
!---------------- FIGURE OUT PHASE
!        IMPORTANT: RIGHT HAIKUS PICK UP AN EXTRA PHASE  FROM OPERATORS
!        COMMUTING PAST THE LEFT HAIKUS
!        this was set up in routine "shiva" in bhoplib.f90
!
!-----------------FIGURE OUT OPERATORS
         opsum = -opleft + rcurch%chain%oplist
         if(fill .and. njumps >=minjump .and. njumps <=maxjump)then

!------------  a^+_ic a^+_jc  a^+ kc a_kd a_jd a_id
          id = opsum(1)
          jd = opsum(2)
          kd = opsum(3)
          kc = opsum(4)
          jc = opsum(5)
          ic = opsum(6)

          if(id < 0)then
            idsps = -id
          else
            idsps = id + nhsps(-it)
          endif
          if(ic < 0)then
            icsps = -ic
          else
            icsps = ic + nhsps(-it)
          endif
          if(jd < 0)then
            jdsps = -jd
          else
            jdsps = jd + nhsps(-it)
          endif
          if(jc < 0)then
            jcsps = -jc
          else
            jcsps = jc + nhsps(-it)
          endif
          if(kd < 0)then
            kdsps = -kd
          else
            kdsps = kd + nhsps(-it)
          endif
          if(kc < 0)then
            kcsps = -kc
          else
            kcsps = kc + nhsps(-it)
          endif
!............... DETERMINE DESTRUCTION/CREATION INDICES.........
!................ASSUME  i > j > k
          dphase = 1
          i = idsps
          j = jdsps
          k = kdsps

          call swap(j,k,dphase)
          call swap(i,j,dphase)
          call swap(j,k,dphase)

          dindx = (i-3)*(i-2)*(i-1)/6 + (j-2)*(j-1)/2 + k

          cphase = 1
          i = icsps
          j = jcsps
          k = kcsps

          call swap(j,k,cphase)
          call swap(i,j,cphase)
          call swap(j,k,cphase)

          cindx = (i-3)*(i-2)*(i-1)/6 + (j-2)*(j-1)/2 + k
          dindx = map(dindx)
          cindx = map(cindx)
!------------------ FIND M, PAR OF OPERATORS

          if(first)then
              if( ic < 0 )then
                m = hspsqn(-it,-ic)%m
                par=hspsqn(-it,-ic)%par
              else
                m = hspsqn(it,ic)%m
                par=hspsqn(it,ic)%par
              endif
              if(jc < 0) then
                 m = m+hspsqn(-it,-jc)%m
                par=par*hspsqn(-it,-jc)%par
              else
                m = m+hspsqn(it,jc)%m
                par=par*hspsqn(it,jc)%par
              endif
              if(kc < 0) then
                 m = m+hspsqn(-it,-kc)%m
                par=par*hspsqn(-it,-kc)%par
              else
                m = m+hspsqn(it,kc)%m
                par=par*hspsqn(it,kc)%par
              endif


              par = (3-par)/2
              first = .false.
           endif
            
           if(fill)then
           iref = XXX3(it)%meref(m,par)
           istart = XXX3(it)%mestart(m,par)
           if(cindx <= dindx) then
             itbme = istart + (dindx-iref)*(dindx-iref-1)/2+cindx-iref
           else
             itbme = istart + (cindx-iref)*(cindx-iref-1)/2+dindx-iref

           endif

           if(fill .and. itbme > maxnmat)then
             print*,' index is too large ',itbme,maxnmat,m,par
             print*,istart,iref
             print*,dindx,cindx
             print*,icsps,jcsps,kcsps
             print*,idsps,jdsps,kdsps
             print*,i,j,k
             stop
           endif
!           if(fill)then
              op(njumps-jumpshift)= itbme
!           write(33,*)njumps,itbme,m,par
          
              phase(njumps-jumpshift)= int(lphase*rphase*cphase*dphase*genphase,kind=1)
           end if

         end if
!-----------------------------------------------------
         rcurch => rcurch%next
      enddo  
      lcurch => lcurch%next
    enddo 

!------------------------------
    mrchain => mrchain%next


   return
   end subroutine weld3b
!====================================================================

   subroutine swap(i,j,phase)
   implicit none
   integer i,j,k
   integer phase

   if(i < j)then
     k = i
     i = j
     j = k
     phase = -phase
   endif
   return
   end subroutine swap

!====================================================================

end module welder

