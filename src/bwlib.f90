!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK configuration-interaction shell-model code
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  file BWLIB.F
!
!  routines to compute limits on W
!  W = excitation Weighting; 
!   (equivalent to Nhw in REDSTICK, t in ANTOINE)
!  
!  subroutines:
!    getWlimits: asks for and controls limits on W excitatations
!       [CALLS THE FOLLOWING SUBROUTINES ]
!    minmax_W  : finds min,maxW possible for given species
!    setWzero  : if no truncations, sets all Ws = 0 for simplicity
!    master_w_limits : calls routines to find limits on nh (# of particles in a haiku), 
!                   limits on W for fixed nh
!    Nhlimits     :  finds limits on nh, W for fixed nh
!    find_Wceiling: finds "ceiling" on W
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   
!    NOTES:  Main goal of these routines is to prepare for cuts on W
!
!    GOAL 1: Find min, max W for proton, neutron states:
!        found in minW(it), maxW(it)
!    GOAL 2: set up for creating haikus: find min, max nh (= # particles in a haiku)
!            and find min, max W for each nh
!        found in: minNh(it), maxNh(it)
!                  minWh(it,nh), maxWh(it,nh)
!    GOAL 3: find "ceiling" of W accessible in single-particle states
!        output: wceiling(it)
!        this is used to determine max number of single-particle states actually used

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine getWlimits
!
!  finds out limits on W, if any
!
!  CALLED BY: 
!   define_system
!
!  SUBROUTINES CALLED:
!    minmax_W
!    setWzero
!    sort_spsqn
!    make_hspsqn
!
subroutine getWlimits

  use system_parameters
  use flagger
  use sporbit
  use W_info
  use io
  use nodeinfo
  use bmpi_mod
  use jumpstart,only:readjumpstart
  implicit none

  integer(kind=4) :: ierr
  integer(kind=4) :: tempW,tempWp,tempWn
  integer(kind=4) :: it
!----------- CONTROL CHARACTER ----------

  character(len=1) :: ychar
      
  if ( auto_readin .or. readjumpstart) then
     if ( iproc == 0 ) then
        print*,' Oops, should not be here (getWlimits)',auto_readin
     end if
	 ! KSM: fixed wrong number of arguments to MPI call
#ifdef _MPI
     call BMPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
     stop
  end if
!............ (Check if W truncation possible )
  maxWtot = 0
  minWtot = 0
  do it = 1,2
     
     call minmax_W(it,minW(it),maxW(it))
     maxWtot = maxWtot + maxW(it) - minW(it)
     minWtot = minWtot + minW(it)
  enddo
!............ (if possible, ask for truncation )
  if ( maxWtot > 0 ) then
     if ( iproc == 0 ) then
        if(auto_input)then
            read(autoinputfile,'(a)')ychar
        else
           print*,' Would you like to truncate ? (y/n/?=more information)'
           read(5,'(a)')ychar
           if(ychar=='?')then
              print*,' '
              print*,' If "n", then many-body space will not be truncated beyond '
              print*,' the specified quantum numbers. '
              print*,' If "y", then the many-body states will be truncated by ' 
              print*,' BIGSTICKs flexible truncation scheme,  '
              print*,' entering in the next step Max excitations.'
              print*,' For no-core shell-model calculations,  Max excite = Nmax = Nhw '
              print*,' For atoms, Max excite = 1-> single excitations, 2 -> doubles, etc.'
              print*,' '
              print*,' For separate proton/neutron limits, enter "p" '
              print*,' This will allow you to enter proton and neutron truncations '
              print*,' above and beyond the summed limit '
              print*,' '
              print*,' For more details consult the Inside Guide. '  
              print*,' '
              print*,' Would you like to truncate ? (y/n)'
              read(5,'(a)')ychar
           end if
           write(autoinputfile,'(a,"     ! truncate? ")')ychar
        endif
     end if
#ifdef _MPI
     call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
     ! call OLDMPI_BCAST(ychar,1,MPI_CHARACTER,0,icomm,ierr)
     call BMPI_BCAST(ychar,1,0,MPI_COMM_WORLD,ierr)
#endif	 
	  select case (ychar)
	 
	    case ('y','Y')
        if ( iproc == 0 ) then

           if(auto_input)then
              read(autoinputfile,*)tempW
           else
              print*,' Max excite (must be less than or equal to  ',maxWtot,' )'
              read*,tempW
              write(autoinputfile,*)tempW,'     ! truncation '
           endif
		   if(tempW > maxWtot .or. tempW < 0)tempW = maxWtot
        end if
#ifdef _MPI
        call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
        ! KSM:  FIXED MPI TYPE ERROR
        call BMPI_BCAST(tempW,1,0,MPI_COMM_WORLD,ierr)
#endif
        if ( tempW .le. maxWtot ) then
           maxWtot = tempW + minW(1)+minW(2)
           maxW(1) = tempW + minW(1)
           maxW(2) = tempW + minW(2)
        end if
		
		case ('p','P')       ! added in 7.6.1

        if ( iproc == 0 ) then

           if(auto_input)then
              read(autoinputfile,*)tempW,tempWp,tempWn
           else
              print*,' Max excite for sum, protons, neutrons? '
			  print*,' (must be less than or equal to  ',maxWtot,maxW(1) - minW(1),maxW(2) - minW(2),', respectively )'
              read*,tempW,tempWp,tempWn
              write(autoinputfile,*)tempW,tempWp,tempWn,'     ! truncation: total, p, n '
           endif
        end if
#ifdef _MPI
        call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
        ! KSM:  FIXED MPI TYPE ERROR
        call BMPI_BCAST(tempW,1,0,MPI_COMM_WORLD,ierr)
        call BMPI_BCAST(tempWp,1,0,MPI_COMM_WORLD,ierr)
        call BMPI_BCAST(tempWn,1,0,MPI_COMM_WORLD,ierr)
#endif		
        if ( tempW .le. maxWtot ) then
           maxWtot = tempW + minW(1)+minW(2)
           maxW(1) = min(tempW + minW(1),tempWp+minW(1))
           maxW(2) = min(tempW + minW(2),tempWn+minW(2))
        end if
		
		case default		
        allsameW = .true.
        maxWtot = 0
        minWtot = 0
        maxW(:) = 0
        minW(:) = 0
        call setWzero
!---------------- NEED TO REMAKE HSPSQN
        call sort_spsqn(1)     ! sort single particle states by W, parity, Jz
        call sort_spsqn(2)
        call make_hspsqn
		
	end select
	 
  else
     maxWtot = maxW(1)+maxW(2)
  end if
  
  return
end subroutine getWlimits
!================================================================
!      subroutine reconstruct_W_limits
!
!  finds out limits on W, if any, when "auto-readin" from prior file
! 
!  CALLED BY:
!   define_system
!
!  SUBROUTINES CALLED:
!    minmax_W
!    setWzero
!    sort_spsqn
!    make_hspsqn
!
subroutine reconstruct_W_limits

  use system_parameters
  use sporbit
  use W_info
  use menu_choices
  use io
  use nodeinfo
  use bmpi_mod
  use jumpstart,only:readjumpstart
  implicit none
  integer(4) :: ierr
  integer(4) :: tempW
  integer(4) :: it
!----------- CONTROL CHARACTER ----------

  character(len=1) :: ychar

      
  if ( .not. auto_readin .and. .not. readjumpstart) then
     if ( iproc == 0 ) print*,' Oops, should not be here in reconstruct ',auto_readin,menu_char
#ifdef _MPI
	 call BMPI_ERR(".not. auto_readin")
#endif
     stop
  end if

  if ( allsameW ) then
     maxWtot = 0
     minWtot = 0
     maxW(:) = 0
     minW(:) = 0
     call setWzero
!---------------- NEED TO REMAKE HSPSQN
     call sort_spsqn(1)     ! sort single particle states by W, parity, Jz
     call sort_spsqn(2)
     call make_hspsqn
     if ( iproc == 0 ) then
        print*,' No truncations on W '
        if(writeout)write(resultfile,*)' No truncations on W '
     end if
     return
  end if

  minWtot = 0
  do it = 1,2
     call minmax_W(it,minW(it),maxW(it))
     minWtot = minWtot + minW(it)
  enddo
  tempW = maxWtot-minWtot  ! # of excitations

  if(allsamew .and. tempW /= 0)then 
     print*,maxWtot,minWtot
     print*,minW
     stop
  endif
  if ( iproc == 0 ) then
     print*,' # of excitations = ',tempW,allsameW
     if(writeout)write(resultfile,*)' # of excitations = ',tempW
  end if
  do it = 1,2
     maxW(it) = tempW+minW(it)
  end do
  
  return
end subroutine reconstruct_W_limits

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine minmax_W
!
!  finds min, max W for either protons or neutrons
!
!  CALLED BY: getWlimits
!
!  INPUT:
!     it: species
!   OUTPUT: 
!   Wmin: min "excitation" possible; base excitation energy
!       from filling all the s.p. states with lowest w first
!   Wmax: max excitation possible
!       from filling all the s.p. states with highest w first
!
subroutine minmax_W(it,Wmin,Wmax)

  use system_parameters
  use spstate
  use sporbit ! temp
  use nodeinfo
  use bmpi_mod
  implicit none

  integer(4) :: ierr

  integer(4) :: N  
  integer(4) :: Wmin,Wmax
  
  integer(4) :: i,j,k,m,it
  integer(4) :: countmin,countmax

!----------------  ERROR TRAP
  if ( np(it) > nsps(it) ) then
     if ( iproc == 0 ) then
        print*,' too many particles for space ',np(it),nsps(it)
     end if
#ifdef _MPI
	 call BMPI_ERR("abort:  np(it) > nsps(it)")
#endif
     stop
  end if
!------------- END OF ERROR TRAP

!...........the single particles states are already ordered by W
!.......... now simply count up min, max W, by adding up
!  ..........either from high end (to get max) or from low end
!............to get min

  Wmin = 0
  Wmax = 0
!.......... OLD ALGORITHM... DOES NOT WORK FOR STERILE ORBITALS  
!  do i =1,Np(it)
!     Wmin = Wmin + spsqn(it,i)%w
!     Wmax = Wmax + spsqn(it,nsps(it)+1-i)%w
!  enddo

countmin = 0 ! keep track of how many particles counted
countmax = 0
do i = 1,nsps(it)
	if(spsqn(it,i)%w /= 99 .and. countmin < np(it))then
		countmin = countmin +1
		Wmin = Wmin + spsqn(it,i)%w
	end if
	if(spsqn(it,nsps(it)+1-i)%w /= 99 .and. countmax < np(it))then
		countmax = countmax +1
		Wmax = Wmax + spsqn(it,nsps(it)+1-i)%w
	end if	
end do
  return
end subroutine minmax_W

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine setWzero
!
!  used to set all Ws = 0 if no truncations
!
! CALLED BY:
!   check_W_orbits
!   reconstruct_W_limits
! 
!
subroutine setWzero

  use sporbit
  use spstate
  use nodeinfo
  use bmpi_mod
  implicit none
  integer(4) :: ierr
  integer(4) :: i,it
!  if ( iproc == 0 ) print*,' setting to zero '
  do it = 1,2
     do i = 1,numorb(it)
        orbqn(it,i)%w = 0
     end do ! i
     do i = 1,nsps(it)
        spsqn(it,i)%w = 0
     end do ! i
  end do ! it
  return
end subroutine setWzero
!============================================================
!  The following routine compute W limits for haikus
!  and are not broadly needed for other uses
!  (e.g., for export to other codes)
!============================================================
!
!      subroutine master_w_limits
!
!
!  master subroutine to determine:
!  min, max nh (= # of particles in a haiku) for protons, neutrons
!  min, max W as a function of nh 
!  CONVENTION:  nh < 0, = "left" haiku ( jz < 0); nh > 0 = "right" haiku (jz >= 0)
!  also finds the min and max nh 
!
!  CALLED BY: main program
!
!  SUBROUTINES CALLED:
!  nhlimits  -- finds limits on nh 
!  find_Wceiling -- finds "ceiling" on W
!  limit_hspspace (in bspstate.f) : cuts off hsp states, finds Nwords in haiku
!
subroutine master_w_limits
  
  use system_parameters
  use sporbit
  use spstate
  use haiku_info
  use nodeinfo
  use bmpi_mod
  use butil_mod
  implicit none
  integer(4) :: ierr
  integer(4) :: maxsize
  integer :: aerr
  
  maxsize = bmax(np(1),np(2))
  maxsize = bmin(maxsize,nhspsmax)
  
  allocate( minWh(-2:2,0:maxsize), stat=aerr )
  if(aerr /= 0) call memerror("master_w_limits 1");
  allocate( maxWh(-2:2,0:maxsize), stat=aerr )
  if(aerr /= 0) call memerror("master_w_limits 2");

  call Nhlimits(1)
  call Nhlimits(2)
  call find_Wceiling(1)

  call find_Wceiling(2)
  call limit_hspspace(1)
  call limit_hspspace(2)

  return
end subroutine master_w_limits

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine Nhlimits(it)
!
!  finds limits on # of particles based upon W
!
  use system_parameters
  use sporbit
  use spstate
  use W_info
  use haiku_info
  use butil_mod

  implicit none
  integer it  ! flavor of species
  integer nh  !  # of particles in a haiku
  integer nhc !  # of conjugate particles in a haiku
  integer :: nhcount
  if(allsameW)then  ! just simply divide up
     maxNh(it) = bmin(np(it),nhsps0(it))
     maxNh(-it)= bmin(np(it),nhsps0(-it))
     minNh(it) = np(it) - maxNh(it)
     minNh(-it) = np(it) - maxNh(-it)
     minWh = 0
     maxWh = 0
     return
  endif

  minWh(it,:) = 0
  maxWh(it,:) = 0
  minWh(-it,:) =0
  maxWh(-it,:)=0
!----------- LOOP OVER THE S.P. STATES adding the min . max  W -------
!  This has been rewritten 7.10.8 to address problems with 'sterile' W = 99 orbitals


!.......... add up the minimal Ws .....
nhcount = 0
do nh = 1,bmin(nhsps_eff(it),np(it))
	if(hspsqn(it,nh)%w==99)cycle
	nhcount = nhcount +1
	minWh(it,nhcount) = minWh(it,nhcount-1)+hspsqn(it,nh)%w
end do


!.... add up the maximal Ws
nhcount = 0
do nh = 1,bmin(nhsps_eff(it),np(it))
	if(hspsqn(it,nhsps_eff(it)-nh+1)%w==99)cycle
	nhcount = nhcount +1
	maxWh(it,nhcount) = maxWh(it,nhcount-1)+hspsqn(it,nhsps_eff(it)-nh+1)%w
end do

nhcount = 0
do nh = 1,bmin(nhsps_eff(-it),np(it))
	if(hspsqn(-it,nh)%w==99)cycle
	nhcount = nhcount +1
	minWh(-it,nhcount) = minWh(-it,nhcount-1)+hspsqn(-it,nh)%w
end do


!.... add up the maximal Ws
nhcount = 0
do nh = 1,bmin(nhsps_eff(-it),np(it))
	if(hspsqn(it,nhsps_eff(-it)-nh+1)%w==99)cycle
	nhcount = nhcount +1
	maxWh(-it,nhcount) = maxWh(-it,nhcount-1)+hspsqn(-it,nhsps_eff(-it)-nh+1)%w
end do

! old code.  causes problems with sterile orbits
!  do nh = 1,bmin(nhsps0(it),np(it))
!     if(hspsqn(it,nh)%w /=99)minWh(it,nh) = minWh(it,nh-1)+hspsqn(it,nh)%w
!     if(hspsqn(it,nhsps0(it)-nh+1)%w/=99)maxWh(it,nh) = maxWh(it,nh-1)+hspsqn(it,nhsps0(it)-nh+1)%w
!	 maxWh(it,nh)= bmin(maxWh(it,nh),maxW(it))   ! added in 7.10.8 to deal with 'sterile' orbits
!             print*,' limits ',nh,minWh(it,nh),maxWh(it,nh)
!  enddo ! nh

!  do nh = 1,  bmin(nhsps0(-it),np(it))
!     if(hspsqn(-it,nh)%w /=99)minWh(-it,nh) = minWh(-it,nh-1)+hspsqn(-it,nh)%w
!     if(hspsqn(-it,nhsps0(it)-nh+1)%w/=99)maxWh(-it,nh) = maxWh(-it,nh-1)+hspsqn(-it,nhsps0(-it)-nh+1)%w
!	 maxWh(-it,nh)= bmin(maxWh(-it,nh),maxW(it))   ! added in 7.10.8 to deal with 'sterile' orbits
	 
!  enddo ! nh

!---------- NOW REDUCE THE MAX 
!           is either the maximum you can fill, or subtracted off from max tot

  do nh = 0,bmin(nhsps_eff(it),np(it))
     
     nhc = np(it) - nh 
     if(nhc > bmin( nhsps_eff(-it),np(it) ) ) cycle
     maxWh(-it,nhc) = bmin(maxWh(-it,nhc), maxW(it)-minWh(it,nh))
     maxWh(-it,nhc) = bmax(maxWh(-it,nhc),0)

  enddo  !nh

  do nh = 0, bmin(nhsps_eff(-it),np(it))
     nhc = np(it) - nh 
     if(nhc > bmin( nhsps_eff(it),np(it) )) cycle
  
     maxWh(it,nhc) = bmin(maxWh(it,nhc), maxW(it)-minWh(-it,nh))
     maxWh(it,nhc) = bmax(maxWh(it,nhc),0)

  enddo  !nh

!C----------- NOW LOOK THROUGH AND ELIMINATE -------

  do nh = 0,bmin(nhsps_eff(it),np(it))
     if(minWh(it,nh)<= maxW(it) .and. minWh(it,nh) <= maxWh(it,nh)) maxNh(it) = nh
  enddo ! nh

  do nh = 0, bmin(nhsps_eff(-it),np(it))
     if(minWh(-it,nh)<= maxW(it).and. minWh(-it,nh) <= maxWh(-it,nh)) maxNh(-it) = nh
  enddo ! nh

!      print*,-it,maxWh(-it,:)
  minNh(it) = np(it) - maxNh(it) 
  minNh(-it) = np(it) - maxNh(-it) 
  
  return
end subroutine nhlimits

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine find_Wceiling
!
!  finds "ceiling", the maximal W on any one particle of species it
!
! CALLED BY: master_W_limits
!

subroutine find_Wceiling(it)

  use system_parameters
  use sporbit
  use spstate
  use W_info
  use haiku_info
  use butil_mod
  implicit none
  integer it
  integer nh
  integer wbase
  integer wtmp
  integer i

  wceiling(it) = 0
  if(allsameW)return
!C------ do for right haikus
  do nh = 1,maxNh(it)  ! loop over all the allowed number of particle
!C-------------- set up base W; lowest of nh - 1 particle
     wbase = 0
     do i = 1,nh-1
		 if(hspsqn(it,i)%w==99)cycle
        wbase = wbase + hspsqn(it,i)%w
     enddo ! i
!C----------------- now move last particle to see what max W is allowed 
     wtmp = 0
     do i = nh,nhsps0(it)
		 if(hspsqn(it,i)%w==99)cycle
		 
        if( wbase + hspsqn(it,i)%w <= maxWh(it,nh))then
           wtmp = hspsqn(it,i)%w
        endif
     enddo
     wceiling(it) = bmax(wtmp,wceiling(it))

  enddo ! nh
!C--------------- now repeat for left haikus -----------------
  do nh = 1,maxNh(-it)  ! loop over all the allowed number of particle
!C-------------- set up base W; lowest of nh - 1 particle
     wbase = 0
     do i = 1,nh-1
		 if(hspsqn(-it,i)%w==99)cycle
		 
        wbase = wbase + hspsqn(-it,i)%w
     enddo ! i
!C----------------- now move last particle to see what max W is allowed 
     wtmp = 0
     do i = nh,nhsps0(-it)
		 if(hspsqn(-it,i)%w==99)cycle
		 
        if( wbase + hspsqn(-it,i)%w <= maxWh(-it,nh))then
           wtmp = hspsqn(-it,i)%w
        endif
     enddo
     
     wceiling(it) = bmax(wtmp,wceiling(it))

  enddo ! nh

  return
end subroutine find_Wceiling

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
