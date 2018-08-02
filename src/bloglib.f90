!
!  added in 7.6.0
! 
!  routines for log files, end-of-run report, etc.
!
!
!

!===========================================================================
!
! print out end-of-run information not already included in .res file
! added in 7.4.1 by CWJ
!
!  prints out the following information:
!  version, date and time
!  platform if known (this should go to another bigstick file)
!  run mode
!  .sps file, N, Z, etc
!  interaction files and parameters
!  MPI and OpenMP parameters
!  Lanczos run parameters 
!  optional comments (if not there skip over)
!
subroutine report
   use nodeinfo
   use reporter
   use io
   use program_info
   use menu_choices
   use system_parameters
   use W_info
   use sporbit
   use fragments
   use lanczos_info
   implicit none
   character(10) :: date,time,zone
   integer :: datetimeval(8)
   integer :: omp_get_num_threads 
   integer :: num_threads
   
   if(iproc/=0 .or. .not.writeout)return

   write(resultfile,*)' '
   write(resultfile,*)' -- END OF RUN REPORT -- END OF RUN REPORT -- END OF RUN REPORT -- '
   write(resultfile,*)' '
   write(resultfile,*)' BIGSTICK Version ',version,lastmodified
   call date_and_time(date,time,zone,datetimeval)
   write(resultfile,*)' Run date: ',date(1:4),'-',date(5:6),'-',date(7:8)
   write(resultfile,*)' Menu choice: ',menu_char
   write(resultfile,*) 
 !$omp parallel shared(num_threads,iproc)  
    num_threads = omp_get_num_threads()
 !$omp end parallel
 write(resultfile,*)'Number of MPI processors =',nproc, ', with # OpenMP threads = ', num_threads
   
!....... MODEL SPACE INFORMATION.......

   if(auto_readin)then
	   write(resultfile,*)' Single particle space read in from originating file ',wfn_in_filename
   else	
       if(spfilename(1:4)=='auto')then
	      write(resultfile,*)' Single particle space h.o. with max Nprincipal = ',nint((sqrt(8.*numorb(1)+1.) -3.)/2.)
       else
	      write(resultfile,*)' Single particle space found in ',spfilename
       endif
   end if
   write(resultfile,*)'Valence Z = ',np(1),', valence N = ',np(2)
   if((jz/2)*2==jz)then  ! even
	   write(resultfile,*)' M = ',jz/2,', parity = ',cparity
   else
	   write(resultfile,*)' M = ',jz,'/2, parity = ',cparity
   end if
   if(.not.allsameW)write(resultfile,*)' Max W excitations = ',maxWtot-minWtot 
   
!....... INFORMATION ON FRAGMENTATION....
   write(resultfile,*)' Number of fragments = ',nfragments
!....... INFORMAITON ON INTERACTION(S).....

!........ INFORMATION ON LANCZOS RUN.......
   write(resultfile,*)' Lanczos option ',lanczchar
   select case (lanczchar)
      case ('ld')
	     write(resultfile,*)nkeep,' vectors kept from ',niter,' max iterations '

	  case ('lf')
   	     write(resultfile,*)nkeep,' vectors kept from ',niter,' iterations (fixed)'

	  case ('td')
   end select
	
   write(resultfile,*)' '
   write(resultfile,*)' -- END OF REPORT -- GOOD-BYE! -- END OF REPORT -- '
   write(resultfile,*)' '
   return
end subroutine report
!====================================================================

!  added in 7.6.0: .log, .bin files
!
!  .log file is human-readable file of run
!  .bin file contains information on, for example, opbundles and # of jumps
!  this will allow for faster subsequent runs

subroutine openlogfiles
	
	use io
	use nodeinfo
	
	implicit none
	integer ilast
	
	if(iproc/=0)return
    ilast = index(outfile,' ')-1  
    if(outfile(1:ilast)=='none')then
		open(unit=logfile,file="logfile.bigstick",status='unknown')
	else
		open(unit=logfile,file=outfile(1:ilast)//".log",status='unknown')
	end if
	
	return
	
end subroutine openlogfiles
!====================================================================
!
!  routine to write various information to log file
!
! INPUT: writedirective is a 3-letter directive as to what data to write
!

subroutine writelogfile(writedirective)
    use system_parameters
	use io
	use nodeinfo
	use menu_choices
	use program_info
    use W_info
    use sporbit
	use jumplimits
	use fragments
	use basis
	use sectors
	use haikus
	use reporter
	use flagger
    use butil_mod
    use localvectors
	use timing
!	use annexations
	
	implicit none
	character*3 :: writedirective
    character(10) :: date,time,zone
    integer :: datetimeval(8)
    integer :: omp_get_num_threads 
    integer :: num_threads
	integer :: iorb,korb
	
	integer :: ifrag
	integer(8) :: minfrag,maxfrag
	real(8)  :: avgfrag,avg2frag
	
	if(iproc/=0)return
	select case (writedirective)
	
	case('ini')    ! initial information ------------------------------------------------------
	
	write(logfile,*)' = = = = = = BIGSTICK run log file = = = = = = '
	write(logfile,*)' Code and run information '
    write(logfile,*)' BIGSTICK Version ',version,lastmodified
    call date_and_time(date,time,zone,datetimeval)
    write(logfile,*)' Run date: ',date(1:4),'-',date(5:6),'-',date(7:8)
    write(logfile,*)' Menu choice: ',menu_char
    write(logfile,*) 
    write(logfile,*) "Running on NERSC_HOST: ", TRIM(nersc_host)
    write(logfile,*)"scratch_dir (*.wfn,...): ", TRIM(scratch_dir)
  !$omp parallel shared(num_threads,iproc)  
     num_threads = omp_get_num_threads()
  !$omp end parallel
    write(logfile,*)'Number of MPI processors =',nproc, ', with # OpenMP threads = ', num_threads
	if(num_threads > 1)then
		write(logfile,*)' Flag wantUseVec2Thread = ',wantUseVec2Thread,' (set in module localvectors )'
	    write(logfile,*)' Flag useHZSomp = ',useHZSomp,' (set in module localvectors )'
		write(6,*)' Flag wantUseVec2Thread = ',wantUseVec2Thread,' (set in module localvectors )'
	    write(6,*)' Flag useHZSomp = ',useHZSomp,' (set in module localvectors )'
	end if
!	if(.not.allsamew)write(logfile,*)' Flag annexation_is_go = ',annexation_is_go,' (set in module annexations)'

!..................... INFORMATION ON PRESET FLAGS...........	

    write(logfile,*)' '
    write(logfile,*)' Information on defaults/flags (found in module flagger in file bmodules_flags.f90):'
!..................... CHECK FOR DIFFERENCES FROM DEFAULT .........................
    if(.not.dosortjumps)write(logfile,*)' Not sorting jumps, dosortjumps = ',dosortjumps	
    if(.not.restrictjumps)write(logfile,*)' Not restricting jump storage, restrictjumps = ',dosortjumps	
    if(.not.setjumpceiling)write(logfile,*)' No ceiling on jump storage, setjumpceiling = ',setjumpceiling
	
	if(nproc > 1 .and. .not.enablelanczosincoreMPI )write(logfile,*)' Storing Lanczos vectors in core not enabled, '& 
	  , ' enablelanczosincoreMPI =',enablelanczosincoreMPI 
  	if(nproc == 1 .and. .not.enablelanczosincore1 )write(logfile,*)' Storing Lanczos vectors in core not enabled, '& 
  	  , ' enablelanczosincore1 =',enablelanczosincore1	  
	
    write(logfile,*)' Limits on storage (found in module flagger in file bmodules_flags.f90):'
	write(logfile,*)' maxjumpmemory_default = ',maxjumpmemory_default,' Gb '
	write(logfile,*)' maxlanczosstorage1 = ',maxlanczosstorage1 ,' Gb '
	
	if(nproc > 1 .and. .not.compactjumpstorage_enabled)write(logfile,*)' Intron deletion not enabled, '  & 
	, ' compactjumpstorage_enabled = ',compactjumpstorage_enabled
	
	if(nproc >1 .and. .not. break_vectors_enabled )write(logfile,*)' Fragments not enabled, ' & 
	, 'break_vectors_enabled = ',break_vectors_enabled 
	write(logfile,*)' flag iffulluseph (if full, use particle-hole transformation ) = ',iffulluseph
	
	write(logfile,*)' '
   flush(logfile)
	
	case('sys')    ! system information -------------------------------------------------
	write(logfile,*)' '
	write(logfile,*)' Many-body system information '
	
	if(auto_readin)then
		write(logfile,*)' Single particle space read in from file ',wfn_in_filename
	else
       if(spfilename(1:4)=='auto')then
 	      write(logfile,*)' Single particle space h.o. with max Nprincipal = ',nint((sqrt(8.*numorb(1)+1.) -3.)/2.)
       else
 	      write(logfile,*)' Single particle space found in ',spfilename
       endif
   end if
!......... WRITE OUT SINGLE PARTICLE ORBITS IN DETAIL.................
	
    if(isoflag .and. .not.pnwtflag)write(logfile,*)' Isospin formalism for single-particle space'
    if(isoflag .and. pnwtflag)write(logfile,*)' Isospin formalism for single-particle space (with different weights W for p,n) '
	if(.not.isoflag)write(logfile,*)'P/N formalism for single-particle space '
	
   ! KSM: FIXME  Temp turn off because korb+iorb is going out of bounds
   ! test case b10nmax4
	if(.false. .and. isoflag)then
		iorb = 0
		do while(iorb <= numorb(1))
		   ! upper case MIN for when we know we don't need bmin for mixed sizes
			write(logfile,'("Label: ",15i4)')(korb,korb = 1,MIN(iorb+15,numorb(1)))
			write(logfile,'("   N : ",15i4)')(orbqn(1,korb+iorb)%nr,korb = 1,MIN(iorb+15,numorb(1)))
			write(logfile,'("   L : ",15i4)')(orbqn(1,korb+iorb)%l,korb = 1,MIN(iorb+15,numorb(1)))
			write(logfile,'(" 2xJ : ",15i4)')(orbqn(1,korb+iorb)%j,korb = 1,MIN(iorb+15,numorb(1)))
			if(pnwtflag)then
		    write(logfile,'(" p W : ",15i4)')(orbqn(1,korb+iorb)%w,korb = 1,MIN(iorb+15,numorb(1)))
		    write(logfile,'(" n W : ",15i4)')(orbqn(1,korb+iorb)%w,korb = 1,MIN(iorb+15,numorb(1)))
			else
			write(logfile,'("   W : ",15i4)')(orbqn(1,korb+iorb)%w,korb = 1,MIN(iorb+15,numorb(1)))				
			end if
			iorb = iorb+15
			if(iorb/=numorb(1))write(logfile,*)' '
		end do 
	else
		write(logfile,*)' Proton-neutron formalism is not yet fully implemented...is it?'
	end if
    write(logfile,*)'Valence Z = ',np(1),', valence N = ',np(2)
    if((jz/2)*2==jz)then  ! even
 	   write(logfile,*)' M = ',jz/2,', parity = ',cparity
    else
 	   write(logfile,*)' M = ',jz,'/2, parity = ',cparity
    end if
    if(.not.allsameW)write(logfile,*)' Max W excitations = ',maxWtot-minWtot 
	write(logfile,*)' '
   flush(logfile)
	
	case ('bas')     !--------- BASIS INFORMATION--------------

	write(logfile,*)' Total basis dimension: ',dimbasis
	write(logfile,*)' # Slater determinants: proton ',nxsd(1),', neutron ',nxsd(2)	
	write(logfile,*)' # sectors: protons ',nsectors(1),', neutrons ',nsectors(2)
	!NEED TO WRITE ROUTINES TO COUNT UP # OF HAKUS!
!	write(logfile,*)' Number of haikus: proton L ',haiku(-1),', proton R ',haiku(1) & 
!	, ', neutron L  ',haiku(-2),', neutron : ',haiku(2)
	
   flush(logfile)
	case ('fra')     !--------- FRAGMENTATION INFORMATION -------
	
	if(nproc==1 .and. break_vectors_enabled )write(logfile,*)' Only 1 MPI proc, cannot break into fragments '
	if(nproc > 1)then
		if(.not.break_vectors_enabled)then
			write(logfile,*)' Fragments not enabled '
			go to 1001                                 ! APOLOGIES FOR THE CLUNKINESS
		end if
		write(logfile,*)' # fragments = ',nfragments
!............ FIND FRAGMENT MIN,MAX....................
        minfrag = 1000000000
        maxfrag=0
        avgfrag=0.d0
        avg2frag=0.d0
        do ifrag =1,nfragments
			print*,ifrag,fragmentlist(ifrag)%localdim
			minfrag = bmin(minfrag,fragmentlist(ifrag)%localdim)	
			maxfrag = bmax(maxfrag,fragmentlist(ifrag)%localdim)	
			avgfrag = avgfrag+real(fragmentlist(ifrag)%localdim,8)
			avg2frag = avg2frag+real(fragmentlist(ifrag)%localdim,8)**2
		end do
		avgfrag = avgfrag/real(nfragments,8)
		avg2frag = avg2frag/real(nfragments,8)
		avg2frag = dsqrt(avg2frag - avgfrag*avgfrag)
		write(logfile,*)' Fragment size min = ',minfrag,', max = ',maxfrag
		write(logfile,*)' Average fragment size = ',int(avgfrag,8),' +/- ',int(avg2frag,8)
        write(logfile,*)' Note: after fragmentation there are ',nsectors(1),' proton sectors (neutrons do not fragment)'			
		
	end if
1001 continue

    case('tim')
   	    write(logfile,*)' '
	    write(logfile,*)' **** Timing summary **** '
        write(logfile,*)' Total time to run : ',endall-startall
		write(logfile,*)' Time to compute lanczos : ',endlanczos-startlanczos
		write(logfile,*)' Time total in H mat-vec multiply  : ',time_hmult
		write(logfile,*)' Time in reorthogonalization : ',timereorthog
   	    write(logfile,*)' '
!   	    write(logfile,*)' *** END OF LOG *** END OF LOG *** END OF LOG *** '
		
	    flush(logfile)
		print*,' '
!		print*,' CLOSING LOG FILE '
!		print*,' '
!		close(logfile)
!
case('end')

    write(logfile,*)' '
    write(logfile,*)' *** END OF LOG *** END OF LOG *** END OF LOG *** '
	
    flush(logfile)
	print*,' '
	print*,' CLOSING LOG FILE '
	print*,' '
	close(logfile)
		
	case default    !---------- ERROR TRAP ------------------
	
	write(logfile,*)' When calling routine writelogfile, directive ',writedirective,' not yet defined '
	
    end select
	
	return

end subroutine writelogfile
!====================================================================


