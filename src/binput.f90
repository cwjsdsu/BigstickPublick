!===========================================================================
!  Program BIGSTICK
!
!  Based upon code REDSTICK by W. E. Ormand
!
!  File BINPUT.f90       started -- CWJ @ SDSU -- June 2008
!
!===================================================================== 
!
!  summary of subroutines found in this file
!
!  GET_ORBIT_INFO: fetches s.p. orbit quantum #s 
!     calls: autofillsps 
!  AUTOFILLSPS : skips reading s.p. file, creates single particle  
!   states in "standard" order; for use with MFDn input files 
!  DEFINE_SYSTEM: set Z, N, Jz, (parity), (W-truncation) 
!  EXTRACT_SPS   : inflates s.p. orbit q#s into sp state q#s (adding m's) 
!  CHECK_PARITY_ORBITS: checks if all orbits have same parity 
!  RESETPARITY : if parity is ignored, set all to "positive" 
!  CHECK_W_ORBITS: checks if all orbits have same W  
! 
!===================================================================== 
 
 
!===================================================================== 
!  Reads in s.p. space information; 
!  Fills in q#s into derived types: orbqn in module sporbit 
! 
!  Reads in either .spo or.sps files 
!  It first looks for .spo file; if that fails,  
!  automatically looks for .sps file of same name 
!
!  ADDED 7.8.1:
!  If it fails to find either of those, it looks for a .sp file (for NuShell)  
! 
!===================================================================== 
! 
! base routine to read in .sps file 
! 
! SUBROUTINES CALLED: 
!   autofillsps  
! 
subroutine get_orbit_info 
 
  use sporbit 
  use io 
  use nodeinfo 
  use reporter
  use bmpi_mod 
  use system_parameters
  implicit none 
 
  integer(4)         :: ierr 
 
!------------------FILE CONTROL--------------------------------------------- 
 
  character (len=15) :: filename 
  character (len=1)  :: achar 
        
  character (len=70) :: title 
  integer            :: ilast 
  logical            :: success,nushellsp 
 
 
!------------------TEMP----------------------------------------------------- 
 
  real               :: xn,xl,xj,xw,xwp,xwn    ! orbital q#s 
!--------------------------------------------------------------------------- 
!                  dummy counters 
!--------------------------------------------------------------------------- 
 
  integer            :: i,j,it 
 
  integer :: aerr
 
! KSM - initialize for -Wuninitialized 
  ilast = 0 
  Acore = 0
 
!------------------OPEN A FILE---------------------------------------------- 
  success = .false. 
 
  do while ( .not. success ) 
     if ( iproc == 0 ) then 
 
        if(auto_input)then 
            read(autoinputfile,'(a)')filename 
        else 
            write (6,*) ' Enter file with s.p. orbit information (.sps or .sp)' 
            write (6,*) ' (Enter "auto" to autofill s.p. orbit info for no-core shell model ) ' 
            write (6,*) ' (Enter "?" for more information ) ' 
            read (5,'(a)') filename 
            if( filename(1:1) /='?')then  
                ilast = index(filename,' ')-1 
                write(autoinputfile,'(a,"    !  name of .sps file ")')filename(1:ilast) 
            else 
                write(6,*)' ' 
                write(6,*)' BIGSTICK requires information on the single-particle space. ' 
                write(6,*)' For general model spaces, or if using OXBASH/NuShell-format interaction files, ' 
                write(6,*)' you need to use a file specifying the single particle states  ' 
                write(6,*)' (e.g., "sd" for the sd valence space indicating the sd.sps file) ' 
                write(6,*)' Can read either .sps (preferred) or .sp files. See manual for more information.  ' 
                write(6,*)' ' 
                write(6,*)' For no-core shell-model calculations, especially if using MFD-format files ' 
                write(6,*)' you can use "auto" to generate the single-particle space ' 
                write(6,*)' by specifying maximal N for the single-particle orbits ' 
                write(6,*)' ' 
                write(6,*)' You can set a path to a standard repository of .sps/.sp files ' 
                write(6,*)' by using the environmental variable BIG_SPS_DIR. ' 
                write(6,*)' Just do : ' 
                write(6,*)' export BIG_SPS_DIR = (directory name) ' 
                write(6,*)' export BIG_SPS_DIR=/Users/myname/sps_repo ' 
              
                if(sps_path=='/')then 
                  print*,' Currently BIG_SPS_DIR is not set ' 
                else 
                  print*,' Currently BIG_SPS_DIR = ',sps_path(:) 
                  print*,' However the current working directory is always checked first ' 
 
                end if 
                print*,' ' 
 
                cycle 
            endif 
        endif 
         
 
     end if 
     call BMPI_BARRIER(icomm,ierr) 
     call BMPI_BCAST(filename,LEN(filename),0,icomm,ierr) 
     ilast = index(filename,' ') - 1 
     if ( filename(1:ilast)=='auto' .or. filename(1:ilast) =='AUTO' ) then 
		 spfilename='auto'
        call autofillsps 
        call BMPI_BARRIER(icomm,ierr) 
        return 
 	 else
		spfilename=filename(1:ilast)
     end if 
 
!..................ATTEMPT TO OPEN .sps FILE.......................... 
!                 ADDED IN 7.3.3 BY WEO:  
!                 LOOKS FOR ENVIRONMENTAL VARIABLE "BIG_SPS_PATH" 
 
     if(iproc == 0)then 
        inquire(file=filename(1:ilast)//'.sps',exist=success) 
        if(success)then 
           open(unit=1,file=filename(1:ilast)//'.sps',status='old') 
        else 
           inquire(file=sps_path(1:length_sps_path)//filename(1:ilast)//'.sps',exist=success) 
           if(success)open(unit=1,file=sps_path(1:length_sps_path)//filename(1:ilast)//'.sps',status='old') 
        end if 
        if(.not.success)write(6,*)filename(1:ilast),'.spo/.sps file does not exist' 
     end if 
	 
!.............. ATTEMPT TO OPEN .sp (NuShell style) FILE...........................
!              ADDED IN 7.8.1 by CWJ  

     nushellsp = .false.
     if(iproc==0 .and. .not. success )then
		 print*,' ... Searching for NuShell-style file ',filename(1:ilast)//'.sp'
         inquire(file=filename(1:ilast)//'.sp',exist=success) 
         if(success)then 
            open(unit=1,file=filename(1:ilast)//'.sp',status='old') 
         else 
            inquire(file=sps_path(1:length_sps_path)//filename(1:ilast)//'.sp',exist=success) 
            if(success)open(unit=1,file=sps_path(1:length_sps_path)//filename(1:ilast)//'.sp',status='old') 
         end if 
         if(.not.success)write(6,*)filename(1:ilast),'.sp file does not exist either ' 
		 if(success)nushellsp = .true.		 
	 end if 	 
	 
     call BMPI_BCAST(success,1,0,icomm,ierr) 
     call BMPI_BCAST(nushellsp,1,0,icomm,ierr) 
     if(nushellsp)then
   	    call get_nushell_sp_info(success)
	  
     end if
	 
  end do 
!----------------  
 
  if ( writeout .and. iproc == 0 )write(resultfile,*)' single-particle file = ',filename(1:ilast) 
  if(nushellsp)return

!------------------READ PAST TITLE CARDS----------------------------- 
 
  success = .false. 
 
  if ( iproc == 0 ) then 
     do while ( .not. success ) 
        read (1,'(a)') achar 
        backspace(1) 
        if ( achar /= '#' .or. achar /= '!' ) then 
           success = .true. 
        else 
           read (1,'(a)') title 
           write (6,*) title 
        end if 
     end do 
  end if 
!------------------READ PAST POSSIBLE LABEL OF ISO/PN----------------------- 
  if ( iproc == 0 ) read (1,'(a)') achar 
  call BMPI_BARRIER(icomm,ierr) 
  call BMPI_BCAST(achar,1,0,icomm,ierr) 
  select case(achar) 
  case ('p','P') 
     isoflag  = .false. 
     !     write (6,*) ' .sps file in pn formalism, cannot handle '  
     !     stop 
  case ('i','I') 
     isoflag = .true.  
     pnwtflag = .false. 
     !              backspace(1) 
  case ('w','W')  ! isospin formalism but different p,n weights 
     isoflag = .true. 
     pnwtflag = .true. 
 
  case default 
     if ( iproc == 0 ) then 
        write (6,*) ' Cannot read type ',achar 
     end if 
     call BMPI_ABORT(icomm,101,ierr) 
     stop 
  end select 
!  call BMPI_BARRIER(icomm,ierr) 
!  call BMPI_BCAST(isoflag,1,0,icomm,ierr) 
!------------------READ # OF ORBITS----------------------------------------- 
  if ( isoflag .and. .not.pnwtflag) then   ! ISOSPIN FORMALISM 
     if ( iproc == 0 ) then 
        read (1,*) numorb(1) 
     end if 
     call BMPI_BARRIER(icomm,ierr) 
     call BMPI_BCAST(numorb(1),1,0,icomm,ierr) 
     numorbmax = numorb(1) 
	 if(numorb(1)<1)then
		 if(iproc==0)print*,' WARNING misread of .sps file -- zero states? '
		 call BMPI_ABORT(icomm,101,ierr)
	 end if
!------------------ALLOCATE MEMORY------------------------------------------ 
     allocate ( orbqn(2,numorbmax) , stat=aerr) 
     if(aerr /= 0) call memerror("get_orbit_info") 
!------------------READ IN-------------------------------------------------- 
     do i = 1, numorb(1) 
        if ( iproc == 0 ) then 
           read (1,*,end=2001) xn, xl, xj, xw 
           orbqn(1,i)%nr = int(xn) 
           orbqn(1,i)%l = int(xl) 
           orbqn(1,i)%par = (-1)**(orbqn(1,i)%l) 
           orbqn(1,i)%j = int(2*xj) 
           orbqn(1,i)%w = int(xw) 
        end if 
        call BMPI_BCAST(orbqn(1,i)%nr,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(1,i)%l,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(1,i)%par,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(1,i)%j,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(1,i)%w,1,0,icomm,ierr) 
     end do 
     numorb(2) = numorb(1) 
     orbqn(2,:) = orbqn(1,:) 
 
  else            ! pn-formalism 
     if ( iproc == 0 ) then 
        read (1,*) numorb(1),numorb(2) 
     end if 
     call BMPI_BARRIER(icomm,ierr) 
     call BMPI_BCAST(numorb(1),1,0,icomm,ierr) 
     call BMPI_BCAST(numorb(2),1,0,icomm,ierr) 
     numorbmax = MAX(numorb(1),numorb(2)) 
!------------------ALLOCATE MEMORY------------------------------------------ 
     allocate ( orbqn(2,numorbmax) ) 
!------------------READ IN-------------------------------------------------- 
     do i = 1, numorb(1) 
        if ( iproc == 0 ) then 
           read (1,*,end=2001) xn, xl, xj, xw 
           orbqn(1,i)%nr = int(xn) 
           orbqn(1,i)%l = int(xl) 
           orbqn(1,i)%par = (-1)**(orbqn(1,i)%l) 
           orbqn(1,i)%j = int(2*xj) 
           orbqn(1,i)%w = int(xw) 
        end if 
        call BMPI_BCAST(orbqn(1,i)%nr,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(1,i)%l,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(1,i)%par,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(1,i)%j,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(1,i)%w,1,0,icomm,ierr) 
     end do 
 
     do i = 1, numorb(2) 
        if ( iproc == 0 ) then 
           read (1,*,end=2001) xn, xl, xj, xw 
           orbqn(2,i)%nr = int(xn) 
           orbqn(2,i)%l = int(xl) 
           orbqn(2,i)%par = (-1)**(orbqn(1,i)%l) 
           orbqn(2,i)%j = int(2*xj) 
           orbqn(2,i)%w = int(xw) 
        end if 
        call BMPI_BCAST(orbqn(2,i)%nr,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(2,i)%l,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(2,i)%par,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(2,i)%j,1,0,icomm,ierr) 
        call BMPI_BCAST(orbqn(2,i)%w,1,0,icomm,ierr) 
     end do 
 
  end if 
  if ( iproc == 0 ) close(unit=1) 
!------------- CHECK FOR NON UNIQUE QUANTUM NUMBERS 
!              ADDED 6/2010 by CWJ 
! 
  do it = 1,2 
      if( numorb(it) < 2)cycle 
      do i = 2,numorb(it) 
          do j = 1, i-1 
             if( orbqn(it,i)%nr == orbqn(it,j)%nr .and. &  
                 orbqn(it,i)%j == orbqn(it,j)%j  .and. & 
                 orbqn(it,i)%l == orbqn(it,j)%l ) then 
                 print*,' Sorry, need unique N L J in .sps file ' 
                 stop 
              endif 
          end do 
 
      end do 
  end do 
 
  return 
!------------------ERROR TRAP FOR END OF FILE------------------------- 
2001 continue 
  if ( iproc == 0 ) then 
     write (6,*) ' sudden end of file in ',filename(1:ilast) 
  end if 
  call BMPI_ABORT(icomm,101,ierr) 
  stop   
end subroutine get_orbit_info 

!===================================================================== 
!
! added in 7.8.1 to read NuShell style .sp files to define single-particle orbits
!
! OUTPUT: 
!   success: logical flag, if F then did not succeed in reading in information
!
subroutine get_nushell_sp_info(success)
    use sporbit 
    use io 
    use nodeinfo 
    use reporter
    use bmpi_mod 
	use system_parameters
    implicit none 	
	logical :: success, test
	character :: testchar
	
	integer :: nraworbits,nspecies
	integer iorb,indx
	integer :: ierr
	integer :: it
	
	success = .false.
!....... READ PAST HEADER............	

if(iproc==0)then
	test = .false.
	do while(.not.test)
		read(1,'(a)')testchar
		if(testchar/='!' .and. testchar/='#')then
			test = .true.
			backspace(1)
		end if
	end do ! while not test
	
	read(1,'(a)')testchar
	
	select case (testchar)
	   case('t')
	       isoflag = .true.
	
	   case('p')
	   
	      isoflag = .false.
	
	   case default
	   
	       print*,' .sp file is not properly formatted, check t/pn '	   
	
    end select
	
	read(1,*)Acore
	print*,' Mass of core = ',Acore
	write(logfile,*)' Mass of core = ',Acore
	read(1,*)nraworbits
	if(isoflag)then
		read(1,*)nspecies,numorb(1)
		numorb(2)=numorb(1)
		success =.true.
	else
		read(1,*)nspecies,numorb(1),numorb(2)
		if(numorb(1)/=numorb(2))then
			print*,' For now, must have same number of proton and neutron orbits, sorry'
			go to 1001
		end if
		success = .true.
	end if
end if

   
1001 continue
    call BMPI_BCAST(success,1,0,icomm,ierr) 
	if(success)then
        call BMPI_BARRIER(icomm,ierr) 
        call BMPI_BCAST(numorb(1),1,0,icomm,ierr) 	
		call BMPI_BCAST(nspecies,1,0,icomm,ierr)
		if(nspecies==2)then
	        call BMPI_BCAST(numorb(2),1,0,icomm,ierr) 
			numorbmax = max(numorb(1),numorb(2))	
		else
			numorb(2)=numorb(1)
			numorbmax = numorb(1)
		end if
		allocate(orbqn(2,numorbmax))
		

		if(iproc==0)then
			do iorb = 1,numorb(1)
				read(1,*)indx,orbqn(1,iorb)%nr,orbqn(1,iorb)%l,orbqn(1,iorb)%j
				orbqn(1,iorb)%w = 0
			end do
			if(.not.isoflag)then
				do iorb = 1,numorb(2)
					read(1,*)indx,orbqn(2,iorb)%nr,orbqn(2,iorb)%l,orbqn(2,iorb)%j
					orbqn(2,iorb)%w = 0
				end do
			else
				do iorb = 1,numorb(1)
					orbqn(2,iorb)%nr = orbqn(1,iorb)%nr
					orbqn(2,iorb)%l  = orbqn(1,iorb)%l
					orbqn(2,iorb)%j  = orbqn(1,iorb)%j
					orbqn(2,iorb)%w  = orbqn(1,iorb)%w
				end do
				
			end if
			close(unit=1)
			
		
		end if
        call BMPI_BARRIER(icomm,ierr) 
		do it = 1,2
		  do iorb = 1,numorbmax
			  call BMPI_BCAST(orbqn(it,iorb)%nr,1,0,icomm,ierr)
			  call BMPI_BCAST(orbqn(it,iorb)%l,1,0,icomm,ierr)
			  call BMPI_BCAST(orbqn(it,iorb)%j,1,0,icomm,ierr)
			  call BMPI_BCAST(orbqn(it,iorb)%w,1,0,icomm,ierr)
			  orbqn(it,iorb)%par = (-1)**orbqn(it,iorb)%l
			  
		  end do
			
		end do
		
		
		
	end if

	
	return

end subroutine get_nushell_sp_info
 
!===================================================================== 
!
! automatically fills out single-particle orbits in "standard" 
! nuclear no-core shell nodel order; useful for interaction input files
! which do not specify the orbits
!
!  CALLED BY:
!    get_orbit_info
! 
subroutine autofillsps 
 
  use io 
  use sporbit 
  use nodeinfo
  use bmpi_mod 
  implicit none 
 
  integer :: ierr 
  integer :: nprinc 
  integer :: i,j,n,l,lparity 
  integer :: it 
  integer :: aerr 
 
  isoflag = .true. 
  if ( iproc == 0 ) then 
     if(auto_input)then 
        read(autoinputfile,*)nprinc 
     else 
        print*,' Enter maximum principle quantum number N ' 
        print*,' (starting with 0s = 0, 0p = 1, 1s0d = 2, etc. ) ' 
        read*,nprinc 
        write(autoinputfile,*)nprinc, '    ! max principal number for autofill of orbits ' 
     endif 
     if ( writeout .and. iproc == 0 )write(resultfile,*)' AUTO s.p. space with max Nprinc = ',Nprinc
	 
  end if 
  call BMPI_BARRIER(icomm,ierr) 
  call BMPI_BCAST(nprinc,1,0,icomm,ierr) 
 
  numorb(1) = (nprinc+1)*(nprinc+2)/2 
  numorb(2) = numorb(1) 
  allocate( orbqn(2,numorb(1) ), stat=aerr) 
  if(aerr /= 0) call memerror("autofillsps") 
  n = 0 
  do i = 0,nprinc 
     lparity = mod(i,2) 
 
     do j = 0,i 
        n = n + 1 
        do it = 1,2 
           orbqn(it,n)%j  = 2*j+1 
           l = j + 1 
           if(mod(l,2)/=lparity)l = l-1 
           orbqn(it,n)%l = l 
           orbqn(it,n)%par = (-1)**l 
           orbqn(it,n)%w  = i 
           orbqn(it,n)%nr = (i - l)/2 
        end do  ! it 
 
     end do   ! j  
  end do  ! i 
 
  return 
end subroutine autofillsps 
!================================================================== 
!  define system parameters: 
!  Z,N ( = np(1), np(2) ), jz, parity,  
!  any truncations by W 
! 
!  SUBROUTINES CALLED 
!      sort_spsqn
!      CHECK_PARITY_ORBITS: checks if all orbits have same parity 
!      resetparity:
!      CHECK_W_ORBITS : checks if all orbits have same W  
!      reconstruct_W_limits: used when autoinput
!      GETWLIMITS     : defines truncation by W 
!================================================================== 
subroutine define_system 
  use system_parameters 
  use io 
  use sporbit 
  use menu_choices 
  use nodeinfo 
  use spstate 
  use bmpi_mod 
  use flagger
  implicit none 
 
  integer(4)       :: ierr 
 
!------------------DUMMY VARIABLES----------------------------- 
  integer(4)       :: i,j 
  integer ::       it 
  integer ::       localwmax(2) 
 
!------------------# OF PROTONS AND NEUTRONS------------------- 
 
  do it = 1,2 
 
     localwmax(it) = 0 
     do i = 1,numorb(it) 
           npmax(it) = npmax(it) + orbqn(it,i)%j+1 
           localwmax(it) = MAX(localwmax(it),orbqn(it,i)%w) 
     end do  ! i 
  end do 
  if ( .not. auto_readin ) then 
 
     if (auto_input)then 
        read(autoinputfile,*)np(1),np(2) 
     else if ( iproc == 0 ) then 
        write(6,'(a34,i3,a18,i3,a1)')' Enter # of valence protons (max ',npmax(1),'), neutrons (max ',npmax(2),')' 
        read*,np(1),np(2) 
        write(autoinputfile,*)np(1),np(2),' ! # of valence protons, neutrons ' 
     end if 
     call BMPI_BARRIER(icomm,ierr) 
     call BMPI_BCAST(np(:),2,0,icomm,ierr)  
  end if 
  if ( writeout .and. iproc == 0 ) write(resultfile,*)np(1),np(2) 
 
!----------- CONVERT PARTICLES TO HOLE OPTION ----------- 
!            CAN BE DELETED WHEN EXPORTING
!
  do it = 1,2 
 
     phconj(it) = .false. 
	 npeff(it)  = np(it)
     if(np(it) < 0)then 
!....... GET MAXIMAL NUMBER OF PARTICLES 
        phconj(it) = .true. 
        if(iproc == 0)then 
           if(it == 1)then
              print*,-np(it),' proton holes = ',npmax(it)+np(it),' protons ' 
              write(logfile,*)-np(it),' proton holes = ',npmax(it)+np(it),' protons ' 
           end if
           if(it == 2)then
              print*,-np(it),' neutron holes = ',npmax(it)+np(it),' neutrons ' 
              write(logfile,*)-np(it),' neutron holes = ',npmax(it)+np(it),' neutrons ' 
           end if
 
           if( abs(np(it)) > npmax(it)/2)print*,' (Are you sure you want to do this? Less than half filled) ' 
        end if  
		npeff(it)=npmax(it)+np(it)
        np(it) = - np(it) 
        if(.not.auto_readin)then 
!................. SWAP W FACTORS................................ 
        do i = 1,numorb(it) 
           orbqn(it,i)%w = localwmax(it) - orbqn(it,i)%w 
        end do 
        do i = 1,nsps(it) 
           spsqn(it,i)%w = localwmax(it) - spsqn(it,i)%w 
        end do  
       call sort_spsqn(it)       ! need to resort 
     elseif( np(it) > npmax(it)/2 .and. iproc == 0 .and. np(it)<npmax(it))then 
         if(it==1)print*,'( You have more than half-filled of protons ' 
         if(it==2)print*,' (You have more than half-filled of neutrons ' 
         print*,' You can induce conversion to particle-hole formalism by entering ',-npmax(it)+np(it),' holes ' 
         print*,' This efficiency is only really helpful for large cases )' 
         end if 
     end if 
!........... OPTION TO USE PARTICLE-HOLE TRUNCATION ON FILLED CASES ........	 
	 if(iffulluseph .and. np(it)==npmax(it))then
		 phconj(it)=.true.
		 np(it) = 0
		 if(iproc==0)then
			 if(it==1)then
				 write(6,*)' Using particle-hole conjugation for filled proton space '
				 write(logfile,*)' Using particle-hole conjugation for filled proton space '
			 else
				 write(6,*)' Using particle-hole conjugation for filled neutron space '
				 write(logfile,*)' Using particle-hole conjugation for filled neutron space '				 
			 end if
		 end if
	 end if
 
  end do   ! END OF PARTICLE-HOLE CONJUGATION
  
!------------------READ IN JZ-------------------------------------- 
 
  if ( .not. auto_readin ) then 
 
     i = 1 
     do while ( (i/2) * 2  /= i ) 
 
        if(auto_input)then 
           read(autoinputfile,*)Jz 
        else if ( iproc == 0 ) then 
           print*,' Enter 2 x Jz of system ' 
           read*,Jz 
           write(autoinputfile,*)Jz,'       ! 2 x Jz of systems ' 
        end if 
        call BMPI_BARRIER(icomm,ierr) 
        call BMPI_BCAST(Jz,1,0,icomm,ierr) 
!..................error trap.................................. 
 
        j = orbqn(1,1)%j    ! check for spinful/spinless fermions 
        if ( (j/2)*2 == j ) then 
           spinless = .true. 
           i = jz            ! must be even independent of particle number 
        else  
           spinless = .false. 
           i = np(1) + np(2) + jz  ! must be even 
        end if 
        if ( (i/2) * 2 == i ) then 
           cycle 
        else 
           if ( iproc == 0 ) print*,' Incorrect value for 2 x Jz ' 
        end if 
     end do 
 
!------------------CHECK FOR PARITY--------------------------- 
     call check_parity_orbits 
     if ( .not. allsameparity ) then 
         
        iparity = 0 
      
        do while ( iparity == 0 ) 
           if(auto_input)then 
               read(autoinputfile,'(a)')cparity 
           else if ( iproc == 0 ) then 
              write (6,'(''Enter parity +/- :'')') 
              write (6,'(''(enter "0" if both parities wanted)'')') 
              read (5,'(a)')cparity 
              write(autoinputfile,'(a,"    ! parity of system ")')cparity 
           end if 
           call BMPI_BARRIER(icomm,ierr) 
           call BMPI_BCAST(cparity,1,0,icomm,ierr) 
!......... REVISED 8/2011 by CWJ TO ALLOW BOTH PARITIES........ 
            
           select case (cparity) 
            
           case ('+') 
             iparity = 1 
           case ('-')  
             iparity = 2 
 
           case ('0') 
             iparity = 1 
             allsameparity = .true. 
             call resetparity 
 
           case default 
               if ( iproc == 0 ) write (6,*) 'Incorrect parity type. Reenter' 
 
           end select 
 
        end do 
     else 
        iparity = 1 
        cparity ='+' 
     end if 
      
  else 
     if ( iparity == 1 ) cparity = '+' 
     if ( iparity == 2 ) cparity = '-' 
  endif 
   
  if ( writeout .and. iproc == 0 ) write(resultfile,*)jz,cparity 
 
!------------------W EXCITATION WEIGHTING------------------------- 
 
  call check_W_orbits   
  if ( .not. allsameW ) then 
     if ( auto_readin ) then 
        call reconstruct_W_limits 
     else 
        call getWlimits 
     end if 
  end if 
 
  return 
end subroutine define_system 
 
!===================================================================== 
!  Determines if all orbits have the same parity 
!  If so, set flag allsameparity in module sporbit = .true. 
!  CALLED BY:
!    define_system
!===================================================================== 
subroutine check_parity_orbits 
 
  use sporbit 
  implicit none 
 
  integer ::  it, i 
 
  allsameparity = .true. 
 
  do it = 1, 2 
     do i = 1, numorb(it) 
        if ( orbqn(it,i)%par /= orbqn(it,1)%par ) allsameparity = .false. 
     end do ! i 
  end do   ! it 
  if ( allsameparity ) orbqn(:,:)%par = 1 
  return 
end subroutine check_parity_orbits 
 
!===================================================================== 
!  Forces all parities to be the same 
!  CALLED BY:
!    define_system
!===================================================================== 
subroutine resetparity 
 
  use sporbit 
  implicit none 
 
  integer ::  it, i 
 
 
  do it = 1, 2 
     do i = 1, numorb(it) 
        orbqn(it,i)%par = 1 
     end do ! i 
  end do   ! it 
  return 
end subroutine resetparity 
 
!===================================================================== 
!  Determines if all orbits have the same W 
!  If so, set flag allsamew in module sporbit = .true. 

!  CALLED BY:
!    define_system
!
!  SUBROUTINES CALLED:
!    setWzero
!
!===================================================================== 
subroutine check_W_orbits 
 
  use sporbit 
  implicit none 
 
  integer ::  it, i 
 
  allsameW = .true. 
 
  do it = 1, 2 
     do i = 1, numorb(it) 
        if (orbqn(it,i)%w /= orbqn(it,1)%w ) allsamew = .false. 
     end do ! i 
  end do   ! it 
  if ( allsameW ) call setWzero 
   
  return 
end subroutine check_W_orbits 
 
!=============================================== 
