!
! initiated in 7.9.3 Oct 2019 by CWJ @ SDSU
!
! routines aimed towards advanced storage
!
! The idea is to store uncoupled TBMEs restricted by:
! change (if any) in Jzp
! change (if any) in parp
! charge in W
!
! Note that change in Jzp, parp is easily handled, because they are blocks
! dW = Wf - Wi, so fixing requires several subblocks
!
!  DESIGN STEPS FOR ADVANCED STORAGE
!   * for each MPI rank/proc, determine the change(s) in quantum numbers
!       this doen by routine XXXXXXXXXX
!   * from this list determine what combinations of pairs are allowed
!
!=========================================================
!
!  added in 7.9.3 towards advanced storage
!  keep track of what quantum numbers are possible on MPI procs
!

module qs_on_proc
	
	implicit none
	
!	integer, allocatable :: ndqs(:)   ! for each MPI proc
	
!    type partner_list
		
		
!	end type partner_list
	type dq_list
		integer :: ndqs
		integer, allocatable :: dJz(:), dpar(:), dWp(:),dWn(:)
	end type dq_list
	
	type (dq_list), allocatable :: dqs_on_proc(:)  
	type (dq_list),target :: dqs_p,dqs_n	
	type (dq_list) :: XY_dqs	
		
contains

!
!  ROUTINE TO SURVEY the quantum numbers on MPI procs
!	
	subroutine survey_qs_on_MPIprocs
		
		use nodeinfo
		use opbundles
		use sectors
		
		implicit none
		integer iprocs
		integer ibundle,jbundle
		
		integer :: ipsector,insector,fpsector,fnsector
!		integer :: myminJz,mymaxJz, myminpar,mymaxpar,mymindWp,mymaxdWp,mymindWn,mymaxdWn
		integer :: bundleJz(nopbundles),bundlepar(nopbundles),bundledWp(nopbundles),bundledWn(nopbundles)
		integer :: ndqs
		logical :: duplicate
		logical :: myverbose = .false.
		
		allocate(dqs_on_proc(0:nprocs-1))
		
! HARVEST QUANTUM # INFO FROM BUNDLES-
		do ibundle = 1,nopbundles
			ipsector = opbundle(ibundle)%isector  ! FIRST FIND SECTOR INFO
			fpsector = opbundle(ibundle)%fsector
			insector = opbundle(ibundle)%insector
			fnsector = opbundle(ibundle)%fnsector
			bundleJz(ibundle)= Xsd(1)%sector(ipsector)%jzX
			bundlepar(ibundle)= Xsd(1)%sector(ipsector)%parX
			bundledWp(ibundle)= Xsd(1)%sector(fpsector)%wX -Xsd(1)%sector(ipsector)%wX
			bundledWn(ibundle)= Xsd(2)%sector(fnsector)%wX -Xsd(2)%sector(insector)%wX

		end do
		
		do iprocs = 0, nprocs-1
			
!-------------------- COUNT UP # OF UNIQUE Q#s ON AN MPI PROC
            ndqs = 1
            do ibundle = opbundlestart(iprocs)+1,opbundleend(iprocs)
				duplicate = .false.
				do jbundle= opbundlestart(iprocs),ibundle-1
					if( bundleJz(jbundle)==bundleJz(ibundle) .and. & 
					bundlepar(jbundle)==bundlepar(ibundle) .and. & 
					bundledWp(jbundle)==bundledWp(ibundle) .and. & 
					bundledWn(jbundle)==bundledWn(ibundle) ) then
					   duplicate= .true.
					   exit
				   end if
					
					
				end do ! jbundle
				if(.not.duplicate)ndqs = ndqs+1
				
			end do ! ibundle			
			dqs_on_proc(iprocs)%ndqs = ndqs
			
			if(iproc==0 .and. myverbose)then
				print*,' On MPI proc ',iprocs,' there are ',ndqs,' unique quantum numbers '
			end if
			
			allocate ( dqs_on_proc(iprocs)%dJz(ndqs))
			allocate ( dqs_on_proc(iprocs)%dpar(ndqs))
			allocate ( dqs_on_proc(iprocs)%dWp(ndqs))
			allocate ( dqs_on_proc(iprocs)%dWn(ndqs))
!-------------------- GO THROUGH AGAIN AND FILL
			
            ndqs = 1
			
			ibundle = opbundlestart(iprocs)
			dqs_on_proc(iprocs)%dJz(ndqs)= bundleJz(ibundle)
			dqs_on_proc(iprocs)%dpar(ndqs)= bundlepar(ibundle)
			dqs_on_proc(iprocs)%dWp(ndqs)= bundledWp(ibundle)
			dqs_on_proc(iprocs)%dWn(ndqs)= bundledWn(ibundle)
			
            do ibundle = opbundlestart(iprocs)+1,opbundleend(iprocs)
				duplicate = .false.
				do jbundle= opbundlestart(iprocs),ibundle-1
					if( bundleJz(jbundle)==bundleJz(ibundle) .and. & 
					bundlepar(jbundle)==bundlepar(ibundle) .and. & 
					bundledWp(jbundle)==bundledWp(ibundle) .and. & 
					bundledWn(jbundle)==bundledWn(ibundle) ) then
					   duplicate= .true.
					   exit
				   end if
					
					
				end do ! jbundle
				if(.not.duplicate)then
					ndqs = ndqs+1
					dqs_on_proc(iprocs)%dJz(ndqs)= bundleJz(ibundle)
					dqs_on_proc(iprocs)%dpar(ndqs)= bundlepar(ibundle)
					dqs_on_proc(iprocs)%dWp(ndqs)= bundledWp(ibundle)
					dqs_on_proc(iprocs)%dWn(ndqs)= bundledWn(ibundle)
					
					
				end if
				
			end do ! ibundle			
			
			
		end do ! iprocs
		
		
		return
	end subroutine survey_qs_on_MPIprocs
!===========================================	

!  ROUTINE TO SURVEY the quantum numbers on MPI procs
!  modified to focus on PN for a given range of processors	
!
	subroutine survey_XY_qs(procstart,procstop)
		
		use nodeinfo
		use opbundles
		use sectors
		
		implicit none
		integer,intent(in) :: procstart,procstop
		integer iprocs
		integer ibundle,jbundle
		integer :: bundlestart,bundlestop
		integer :: ipsector,insector,fpsector,fnsector
		integer,allocatable :: bundleJz(:),bundlepar(:),bundledWp(:),bundledWn(:)
		logical,allocatable :: okbundle(:)
		integer :: ndqs
		logical :: duplicate,foundit
		logical :: myverbose = .false.
				
		bundlestart = nopbundles
		bundlestop = 0
		do iprocs= procstart,procstop
			bundlestart = min(bundlestart,opbundlestart(iprocs))
			bundlestop = max(bundlestop,opbundleend(iprocs))
			
		end do		
!............... NOW GO THROUGH AND MAKE SURE ONE HAS PN BUNDLES.....

        do ibundle=bundlestart,bundlestop
			if(opbundle(ibundle)%optype=='PN')then
				foundit = .true.
				bundlestart = ibundle
				exit
			end if
			
		end do
		if(.not.foundit)then
			print*,' No PN bundles found on the range of processors ',procstart,procstop
			return
		end if		
		
		allocate(bundleJz(bundlestart:bundlestop))
		allocate(bundlepar(bundlestart:bundlestop))
		allocate(bundledWp(bundlestart:bundlestop))
		allocate(bundledWn(bundlestart:bundlestop))
		allocate(okbundle(bundlestart:bundlestop))
		okbundle=.false.
		
! HARVEST QUANTUM # INFO FROM BUNDLES-

        do iprocs = procstart,procstop
			do ibundle=opbundlestart(iprocs),opbundleend(iprocs)
				if(opbundle(ibundle)%optype=='PN')okbundle(ibundle)=.true.
				
			end do
		end do
		do ibundle = bundlestart,bundlestop
			if(.not.okbundle(ibundle))cycle
			
			ipsector = opbundle(ibundle)%isector  ! FIRST FIND SECTOR INFO
			fpsector = opbundle(ibundle)%fsector
			insector = opbundle(ibundle)%insector
			fnsector = opbundle(ibundle)%fnsector
			bundleJz(ibundle)= Xsd(1)%sector(ipsector)%jzX
			bundlepar(ibundle)= Xsd(1)%sector(ipsector)%parX
			bundledWp(ibundle)= Xsd(1)%sector(fpsector)%wX -Xsd(1)%sector(ipsector)%wX
			bundledWn(ibundle)= Xsd(2)%sector(fnsector)%wX -Xsd(2)%sector(insector)%wX

		end do
			
!-------------------- COUNT UP # OF UNIQUE Q#s ON AN MPI PROC
            ndqs = 1
            do ibundle = bundlestart+1,bundlestop
				if(.not.okbundle(ibundle))cycle
				
				duplicate = .false.
				do jbundle= bundlestart,ibundle-1
					if(.not.okbundle(jbundle))cycle
					
					if( bundleJz(jbundle)==bundleJz(ibundle) .and. & 
					bundlepar(jbundle)==bundlepar(ibundle) .and. & 
					bundledWp(jbundle)==bundledWp(ibundle) .and. & 
					bundledWn(jbundle)==bundledWn(ibundle) ) then
					   duplicate= .true.
					   exit
				   end if
					
					
				end do ! jbundle
				if(.not.duplicate)ndqs = ndqs+1
				
			end do ! ibundle			
			XY_dqs%ndqs = ndqs
						
			if(iproc==0 .and. myverbose)then
				print*,' On MPI procs ',procstart,procstop,' there are ',ndqs,' unique quantum numbers '
			end if
			
			allocate ( XY_dqs%dJz(ndqs))
			allocate (XY_dqs%dpar(ndqs))
			allocate ( XY_dqs%dWp(ndqs))
			allocate ( XY_dqs%dWn(ndqs))
!-------------------- GO THROUGH AGAIN AND FILL
			
            ndqs = 1
			
			ibundle = bundlestart
			XY_dqs%dJz(ndqs)= bundleJz(ibundle)
			XY_dqs%dpar(ndqs)= bundlepar(ibundle)
			XY_dqs%dWp(ndqs)= bundledWp(ibundle)
			XY_dqs%dWn(ndqs)= bundledWn(ibundle)
			
            do ibundle = bundlestart+1,bundlestop
				if(.not.okbundle(ibundle))cycle
				duplicate = .false.
				do jbundle= bundlestart,ibundle-1
					if(.not.okbundle(jbundle))cycle
					if( bundleJz(jbundle)==bundleJz(ibundle) .and. & 
					bundlepar(jbundle)==bundlepar(ibundle) .and. & 
					bundledWp(jbundle)==bundledWp(ibundle) .and. & 
					bundledWn(jbundle)==bundledWn(ibundle) ) then
					   duplicate= .true.
					   exit
				   end if
					
					
				end do ! jbundle
				if(.not.duplicate)then
					ndqs = ndqs+1
					XY_dqs%dJz(ndqs)= bundleJz(ibundle)
					XY_dqs%dpar(ndqs)= bundlepar(ibundle)
					XY_dqs%dWp(ndqs)= bundledWp(ibundle)
					XY_dqs%dWn(ndqs)= bundledWn(ibundle)
					
					
				end if
				
			end do ! ibundle			
			
		return
	end subroutine survey_XY_qs
!===========================================	
!
!  find all dq#s possible (on all MPI processes)
!  It does this by taking the difference of quantum numbers between sectors
!

	subroutine survey_qs_X(it)
		use nodeinfo
		use sectors
		use sporbit,only:parmult
		
		implicit none
		
		integer :: it
		integer :: is,js
		integer :: n,iq
		
		integer :: nqs
		logical :: duplicate
		type (dq_list),pointer :: xdq
	    integer, pointer :: dwX(:)
		
		integer, allocatable :: dJz(:),dpar(:),dW(:)
		
		
		n = nsectors(it)*(nsectors(it)+1)/2
		
		allocate( dJz(n),dpar(n),dW(n))
		
!		print*,it,' (A)'

!.......... FIND ALL DQs from sectors	
		iq = 0
		do is = 1,nsectors(it)
			do js = 1, is !,nsectors(it)
				iq = iq +1
				dJz(iq)= Xsd(it)%sector(js)%jzX - Xsd(it)%sector(is)%Jzx
				dW(iq)= Xsd(it)%sector(js)%wX - Xsd(it)%sector(is)%wX
				dpar(iq)= parmult (Xsd(it)%sector(js)%parx,Xsd(it)%sector(is)%parx )
				
				
			end do
		end do
		
!............ COUNT UP # OF UNIQUE dqs		
		nqs = 1
		do is = 2,n
			duplicate = .false.
			do js = 1,is-1
				if( dJz(is)==dJz(js) .and. dpar(is)==dpar(js) .and. dW(is)==dW(js))then
					duplicate = .true.
					exit
				end if
				
			end do   !
			if(.not.duplicate)nqs = nqs+1
			 
			
		end do ! is
		
		if(it==1)then
			allocate(dqs_p%djz(nqs))
			allocate(dqs_p%dpar(nqs))
			allocate(dqs_p%dWp(nqs))
			xdq => dqs_p
			dWx => dqs_p%dWp
			
		else
			allocate(dqs_n%djz(nqs))
			allocate(dqs_n%dpar(nqs))
			allocate(dqs_n%dWn(nqs))
			xdq => dqs_n
			dWx => dqs_n%dWn
			
		end if
		xdq%ndqs = nqs
		
		if(iproc==0)print*,it,nqs,' distinct quantum numbers '
		
		nqs = 1
		xdq%djz(1)= 0
		xdq%dpar(1)= 1
		dWx(1)= 0
		
		do is = 2,n
			duplicate = .false.
			do js = 1,is-1
				if( dJz(is)==dJz(js) .and. dpar(is)==dpar(js) .and. dW(is)==dW(js))then
					duplicate = .true.
					exit
				end if
				
			end do   !
			if(.not.duplicate)then
				nqs = nqs+1
				xdq%djz(nqs)= dJz(is)
				xdq%dpar(nqs)= dpar(is)
				dWx(nqs)= dW(is)	
				
			end if
			
		end do ! is
		
!		print*,' a check ',it,xdq%djz(:)
		
		deallocate(djz,dpar,dw)
		
		return

		
	end subroutine survey_qs_X
!----------------------------------------------------------------
!  checks that quantum numbers allowed
!	
	logical function check_allowed_qs_X(it,dMx,dparx,dWx)
		
		implicit none
		integer :: it
		integer :: dMx,dparx,dwx
!		logical,intent(out) :: ok
		
		integer:: iq
		type (dq_list),pointer :: xdq
	    integer, pointer :: mydwX(:)
				
		if(it==1)then
			xdq => dqs_p
			mydWx => dqs_p%dWp			
		else
			xdq => dqs_n
			mydWx => dqs_n%dWn	
			
		end if
		
		check_allowed_qs_X = .false.
		
		do iq = 1,xdq%ndqs
			if(abs(dMx) == abs(xdq%djz(iq)) .and. dparx==xdq%dpar(iq) .and. dWx == mydWx(iq))then
				check_allowed_qs_X = .true.
				return
			end if
		end do
		
		return
		
	end function check_allowed_qs_X

	
	
end module qs_on_proc
!=========================================================
!=========================================================
!=========================================================
module adv_tbme_info
	implicit none
	
	logical :: sort_XXpairs_on_W = .false.   !  NOTE THIS CAUSES PROBLEMS,
	
	integer(8), allocatable :: PNmestart(:,:,:,:)  ! starting point with Jz, par, dWp,dWn,Wb
	integer, allocatable :: PNmeref(:,:,:,:)  ! relative #me  with Jz, par, dWp,dWn,Wb
	integer, allocatable :: PNmeconj(:,:,:,:)  ! start conjugate p-n rho  with Jz, par, dWp,dWn,Wb
	
	
contains

! THIS ROUTINE SORTS PAIRS ON W
! assumes already sorted on Jz, par
!
!  INPUT:
!     npair = # of pairs
!     pair(:) = derived type of pairs


	subroutine sortXXpairsW(npair,pair,sortWb)
		use pairdef
		implicit none
!........ INPUT..............
!integer,intent(IN) :: it  ! species
		integer,intent(IN) :: npair
		logical,intent(IN) :: sortWb !  added to sort on the W of the 2nd particle only
        type (pair_qn) :: pair(npair)
			
		type (pair_qn) :: pairtmp
		
		integer :: jz,par
		integer :: pairstart,pairstop,ipair
		integer :: n,l,i,j
		integer :: wpairstart,wpairstop
		integer :: w, iw
		integer :: jpair
		
		logical :: check=.true.
		logical :: sortWbb=.false.
			
		if(.not.sort_XXpairs_on_W)return		
		
!		do jpair = 1,npair
!			print*,jpair,pair(jpair)%m,pair(jpair)%par!,pair(jpair)%w,pair(jpair)%wb
	
!		end do
		
		pairstop = 0
		pairstart = 1
!		print*,'w ',pair%w

!		print*,'before ',pair%indx
		do while(pairstart <= npair)
			
!............... FIND LIMITS...........................			
			jz = pair(pairstart)%m
			par= pair(pairstart)%par
			pairstop = npair
			do ipair=pairstart+1,npair
				if(pair(ipair)%m /= jz .or. pair(ipair)%par/=par)then
					pairstop = ipair-1
					exit
				end if
				
			end do
			
			if(check)then
				do ipair =pairstart,pairstop
					if(pair(ipair)%m /= pair(pairstart)%m .or. & 
					 pair(ipair)%par /= pair(pairstart)%par  )then
					 print*,' not sorted properly on m, par '
					 print*,pair(ipair)%m , pair(pairstart)%m
					 stop
				 end if
				end do
				
				
			end if
!			if(pairstop == -1)then  ! didn't find end
!				pairstop = npair
!			end if
			
!.............. NOW SORT WITHIN THESE LIMITS................
!     HEAPSORT ON W


           l = (pairstop-pairstart+1)/2+pairstart
           n = pairstop
           do while (n >= pairstart)
             if(l > pairstart)then
                l = l- 1
                pairtmp = pair(l)
             else
                pairtmp = pair(n)
                pair(n) = pair(pairstart)
                n = n -1
                if(n == pairstart)then
                    pair(pairstart) = pairtmp
                    exit
                 end if
             end if
             i = l
             j = l+l-(pairstart-1)
             do while ( j <= n)
                if( j < n)then
                      if(pair(j)%w< pair(j+1)%w)j = j+1
                endif
                if( pairtmp%w < pair(j)%w)then 
                    pair(i) = pair(j)
                    i = j
                    j = j + j-(pairstart-1)
                else
                    j = n + 1
                end if
             end do      
             pair(i) = pairtmp
         end do
!.............. IF ADDITIONAL SORTING ON wb DO IT HERE..................		
         if(sortWb)then 
			 
!................ SEARCH FOR A BLOCK OF W TO SORT.........			 
			 wpairstop = -1
			 wpairstart = pairstart
			 do while(wpairstart <= pairstop)   
			    w= pair(wpairstart)%w
                wpairstop = pairstop
				do ipair=wpairstart+1,pairstop
					if(pair(ipair)%w/= w )then
						wpairstop = ipair-1
						exit
					end if
				
				end do				
				if(wpairstop == -1)then  ! didn't find end
!					print*,w,'w  did not find end '
					wpairstop = pairstop
				end if				
				if(wpairstop < wpairstart)then
					print*,' some problem '
					stop
				end if
				
				if(check)then
					do ipair = wpairstart,wpairstop
						if(pair(ipair)%m/= pair(wpairstart)%m .or. & 
						pair(ipair)%par/= pair(wpairstart)%par .or. & 
						pair(ipair)%w/= pair(wpairstart)%w)then 
						print*,' not sorted '
						stop
						
					end if
				end do
					
					
				end if
				
!......................... NOW SORT THAT BLOCK			
                l = (wpairstop-wpairstart+1)/2+wpairstart
                n = wpairstop
                do while (n >= wpairstart)
                    if(l > wpairstart)then
                       l = l- 1
                       pairtmp = pair(l)
                    else
                       pairtmp = pair(n)
                       pair(n) = pair(wpairstart)
                       n = n -1
                       if(n == wpairstart)then
                           pair(wpairstart) = pairtmp
                           exit
                       end if
                    end if
                   i = l
                   j = l+l-(wpairstart-1)
                   do while ( j <= n)
                       if( j < n)then
                           if(pair(j)%wb< pair(j+1)%wb)j = j+1
                        endif
                        if( pairtmp%wb < pair(j)%wb)then 
                           pair(i) = pair(j)
                           i = j
                           j = j + j-(wpairstart-1)
                        else
                           j = n + 1
                        end if
                    end do      
                    pair(i) = pairtmp
                 end do  ! while
				 wpairstart=wpairstop+1
				 if(wpairstart>pairstop)exit
 			end do  !while 
			if(check)then
				do ipair = wpairstart,wpairstop
					if(pair(ipair)%m/= pair(wpairstart)%m .or. & 
					pair(ipair)%par/= pair(wpairstart)%par .or. & 
					pair(ipair)%w/= pair(wpairstart)%w)then 
					print*,' not sorted b '
					stop
					
				end if
			end do
			do ipair = wpairstart+1,wpairstop
				if(pair(ipair)%wb < pair(ipair-1)%wb)then
					print*,' wb not sorted properly '
					stop
				end if
				
			end do
				
				
			end if
		 end if
!......... END OF SORTING ON wb.................................
!.................................................................. 
		   pairstart = pairstop +1
		   if(pairstop==npair)exit	
			
		end do ! while (pairstop < npair)
!		print*,'after ',pair%indx

!do jpair = 1,npair
!	print*,jpair,pair(jpair)%m,pair(jpair)%par,pair(jpair)%w,pair(jpair)%wb
	
!end do
		
	    return
	
	end subroutine sortXXpairsW
!================================================================	
!
! finds how large a block of (sorted) pairs is  
!  ADDED 7.9.4 as part of campaign for advanced storage design for uncoupled matrix elements
!
! INPUT:
!    npair, # of "pairs" 
!   pair = array of pairs (in m-scheme) labeled a,b
!   mx, parx = Jz and parity of a pair
!   findW : if true, also restrict on W
!   Wx = W of pair
!   findmaxWb: if true, find the max W on the b particle (only used in particle-hole "pairs")
!   maxWb
!   restrictmaxWb: if true, use maxWb to restrict
!
!  OUTPUT:
!    pairblocksize -- number of pairs with matching quantum numbers, assumes sorted
!    pairstart, pairstop
!    notfound = flag that no such block found
!
subroutine find_pair_block_size(npair,pair,mx,parx,fixw,wx,fixwb,maxwb,restrictmaxWb,pairblocksize,pairstart,pairstop,notfound)
	use pairdef
	implicit none
!	integer,intent(in) :: it
	integer,intent(in) :: npair
	type (pair_qn),intent(in) :: pair(npair)
	integer,intent(in) :: mx
	integer,intent(in) :: parx
	integer,intent(in) :: Wx
	integer,intent(in) :: maxwb
	logical,intent(in) :: fixW,fixWb,restrictmaxWb  ! control flags
	integer,intent(out) :: pairblocksize
	logical,intent(out) :: notfound
	
	integer :: ipair,pairstart,pairstop
	integer :: jpair
	
	if(fixWb .and. restrictmaxWb)then
		print*,' In find_pair_block_size, cannot have both flags true '
		print*,' (findmaxwb and restrictmaxwb)'
		stop
	end if
	
	notfound = .true.
!	if(findmaxwb)maxwb= 0
	do ipair = 1,npair
!		print*,ipair,pair(ipair)%m,mx, ' cf M'
		if(pair(ipair)%m /= mx)cycle
!		if(pair(ipair)%w /= wx)cycle
		if(pair(ipair)%par /= parx)cycle
		
		if(fixW .and. pair(ipair)%w/=wx)cycle
		
		if(fixWb .and. pair(ipair)%wb /= maxwb)cycle
		if(restrictmaxwb .and. pair(ipair)%wb> maxwb)cycle
		
!		if(findmaxwb)maxwb = max(maxwb,pair(ipair)%wb)
!		if(findWb .and. pair(ipair)%wb/=wb)cycle
		
		pairstart = ipair
		notfound = .false.
		exit
		
	end do ! ipair
	
!    if(.not.notfound)
	
	if(notfound)return

	do ipair = npair,1,-1
		if(pair(ipair)%m /= mx)cycle
		if(fixW .and. pair(ipair)%w /= wx)cycle
		if(pair(ipair)%par /= parx)cycle
		if(fixWb .and. pair(ipair)%wb /= maxwb)cycle
		if(restrictmaxwb .and. pair(ipair)%wb> maxwb)cycle
		
!		if(findW .and. pair(ipair)%w/=wx)cycle
!		if(findWb .and. pair(ipair)%wb/=wb)cycle
		
		pairstop = ipair
		exit
		
	end do ! ipair
	
	pairblocksize = pairstop - pairstart+1
	
	do ipair = pairstart, pairstop
!		print*, pair(ipair)%m , pair(ipair)%par, pair(ipair)%w
		if(pair(ipair)%m /= mx .or. pair(ipair)%w /= wx .or. pair(ipair)%par/=parx)then
			print*,' mismatch ',ipair, pairstart,pairstop
!			print*, pair(pairstop)%m , pair(pairstop)%par, pair(pairstop)%w
do jpair = pairstart, pairstop
	print*,jpair,pair(jpair)%m,pair(jpair)%par,pair(jpair)%w,pair(jpair)%wb
	
end do

			stop
		end if
		
	end do
	
	return
end subroutine find_pair_block_size


!================================================================
	
end module adv_tbme_info
!=========================================================
!=========================================================
!=========================================================
module TRstuff
	implicit none
	
	logical :: useTRphase_def = .false.
	logical :: useTRphase 
	
contains
!===================================================================
! routine to find TR pairs (switching m -> -m) and induced phase
! Jan 2012 CWJ
! only used when storing 3-body with TR symmetry
!
! modified in 7.9.4
!

subroutine findTRpairsXX(it)
   use interaction
   use haiku_info
   use spstate

   implicit none

   integer it   ! which species

   integer ipair,jpair
   integer a,b,atr,btr,ma,mb,ja,jb
   integer asps,bsps,asgn,bsgn
   integer ath, bth
   integer c,d

   do ipair = 1,npairXX(it)
       XX2(it)%pair(ipair)%tr = ipair
       XX2(it)%pair(ipair)%trphase = 1
   end do
   do ipair = 1,npairXX(it)
	   if(XX2(it)%pair(ipair)%m >= 0)cycle
       a = XX2(it)%pair(ipair)%ia
       b = XX2(it)%pair(ipair)%ib
!......... UNPACK THE INDICES........
       if(a < 0)then
           ath = -it
           asps = -a
           asgn = -1
       else
           ath  = it
           asps = a
           asgn = +1
       endif
       ma = hspsqn(ath,asps)%m
       ja = hspsqn(ath,asps)%j
       if(ma == 0)then
         atr = a
       else
         atr = -asgn*hspsqn(ath,asps)%tr
       endif
       if(b < 0)then
           bth = -it
           bsps = -b
           bsgn = -1
       else
           bth  = it
           bsps = b
           bsgn =+1
       endif
       mb = hspsqn(bth,bsps)%m
       jb = hspsqn(bth,bsps)%j
       if(mb == 0)then
          btr = b
       else
         btr = -bsgn*hspsqn(bth,bsps)%tr
       endif
!........ CHECK IF SELF-CONJUGATE........
       if(ma == 0 .and. mb == 0)then
          XX2(it)%pair(ipair)%tr = ipair
          XX2(it)%pair(ipair)%trphase = int( (-1)**( (ja + jb)/2), kind=1)
          cycle
       endif

       do jpair = 1,npairXX(it)
           c = XX2(it)%pair(jpair)%ia
           d = XX2(it)%pair(jpair)%ib
!......... UNPACK THE INDICES........
           if( atr == c .and. btr == d)then
               XX2(it)%pair(ipair)%tr  = jpair
               XX2(it)%pair(ipair)%trphase = int( (-1)**( (ja + jb)/2), 1)
               XX2(it)%pair(jpair)%tr  = ipair
               XX2(it)%pair(jpair)%trphase = int( (-1)**( (ja + jb)/2), 1)
               go to 1111
           end if
           if( atr == d .and. btr == c)then  ! addition phase for swapping
               XX2(it)%pair(ipair)%tr  = jpair
               XX2(it)%pair(ipair)%trphase = int( -(-1)**( (ja + jb)/2), 1)
               XX2(it)%pair(jpair)%tr  = ipair
               XX2(it)%pair(jpair)%trphase = int( -(-1)**( (ja + jb)/2), 1)
               go to 1111
           end if

       end do ! jpair
!       print*,' Never found TR conjugate ',it,ipair
!       print*,ja,ma,jb,mb
!       stop
1111   continue

!print*,' testing ',XX2(it)%pair(ipair)%m ,XX2(it)%pair(jpair)%m,ipair,jpair
   end do  ! ipair
   return
end subroutine findTRpairsXX

!=========================================================

subroutine findTRpairs_gen(it,npair,pair)
   use interaction
   use haiku_info
   use spstate

   implicit none

   integer it   ! which species
   integer :: npair
   type (pair_qn) :: pair(npair)

   integer ipair,jpair
   integer a,b,atr,btr,ma,mb,ja,jb
   integer asps,bsps,asgn,bsgn
   integer ath, bth
   integer c,d

   do ipair = 1,npair
       pair(ipair)%tr = ipair
       pair(ipair)%trphase = 1
   end do
   do ipair = 1,npair
	   if(pair(ipair)%m >= 0)cycle
       a = pair(ipair)%ia
       b = pair(ipair)%ib
!......... UNPACK THE INDICES........
       if(a < 0)then
		   if(it > 0)then
               ath = -it
		   else   ! pn
			   ath = -1
		   end if
           asps = -a
           asgn = -1
       else
		   if(it > 0)then
              ath  = it
		   else
			  ath = 1
		   end if
           asps = a
           asgn = +1
       endif
       ma = hspsqn(ath,asps)%m
       ja = hspsqn(ath,asps)%j
       if(ma == 0)then
         atr = a
       else
         atr = -asgn*hspsqn(ath,asps)%tr
       endif
       if(b < 0)then
		   if(it > 0)then
              bth = -it
		   else
			  bth = -2
		   end if
           bsps = -b
           bsgn = -1
       else
		   if(it > 0)then
              bth  = it
		   else
			   bth = 2
		   end if
           bsps = b
           bsgn =+1
       endif
       mb = hspsqn(bth,bsps)%m
       jb = hspsqn(bth,bsps)%j
       if(mb == 0)then
          btr = b
       else
         btr = -bsgn*hspsqn(bth,bsps)%tr
       endif
!........ CHECK IF SELF-CONJUGATE........
       if(ma == 0 .and. mb == 0)then
          pair(ipair)%tr = ipair
          pair(ipair)%trphase = int( (-1)**( (ja + jb)/2), kind=1)
          cycle
       endif

       do jpair = 1,npair
           c = pair(jpair)%ia
           d = pair(jpair)%ib
!......... UNPACK THE INDICES........
           if( atr == c .and. btr == d)then
               pair(ipair)%tr  = jpair
               pair(ipair)%trphase = int( (-1)**( (ja + jb)/2), 1)
               pair(jpair)%tr  = ipair
               pair(jpair)%trphase = int( (-1)**( (ja + jb)/2), 1)
               go to 1111
           end if
           if( atr == d .and. btr == c)then  ! addition phase for swapping
               pair(ipair)%tr  = jpair
               pair(ipair)%trphase = int( -(-1)**( (ja + jb)/2), 1)
               pair(jpair)%tr  = ipair
               pair(jpair)%trphase = int( -(-1)**( (ja + jb)/2), 1)
               go to 1111
           end if

       end do ! jpair
!       print*,' Never found TR conjugate ',it,ipair
!       print*,ja,ma,jb,mb
!       stop
1111   continue

!print*,' testing ',XX2(it)%pair(ipair)%m ,XX2(it)%pair(jpair)%m,ipair,jpair
   end do  ! ipair
   return
end subroutine findTRpairs_gen


end module TRstuff
!=========================================================
!=========================================================
!=========================================================
