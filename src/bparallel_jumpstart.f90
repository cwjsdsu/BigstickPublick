!
!  NOTE 7.9.12: I may deprecate these routines
!
!  JUMPSTART data and routines
!  ADDED 7.9.6 by CWJ
!  
!  'jumpstarting':  the code does an initial run and gets timing information, then stops
!  on second round of running, the code reads in fine-mesh timing, and uses that to compute
!  the distribution
!
!  file is OUTPUTNAME.jumpstart.bigstick
!

module jumpstart
	
	implicit none
	
	logical :: writejumpstart, readjumpstart	
	integer :: jumpstartfile = 76
	
	logical :: blockflag
	integer :: blocksize
	
contains
!====================================================================================================	
!	
! Routine to open the binary jumpstart file, with naming convention OUTPUTNAME.jumpstart.bigstick
!
	subroutine openjumpstartfile2write
		use io
		use nodeinfo
		
		implicit none
		integer :: ilast
		
		if(iproc /= 0)return
		if(nproc ==1 )return
		
		print*,' '
		print*,' Writing information to jumpstart (binary) file '
		ilast = index(outfile,' ')-1
		open(unit=jumpstartfile,file=outfile(1:ilast)//".jumpstart.bigstick",action='WRITE', access='stream', & 
		    form='unformatted',status='unknown',err=3033)
		
		return
		
3033    continue
        print*,' Hmm, some problem opening jumpstart file, do not know why '
		print*,' Giving up '
		stop 	
		
	end subroutine openjumpstartfile2write
!
!====================================================================================================	
!	
! Routine to open the binary jumpstart file, with naming convention OUTPUTNAME.jumpstart.bigstick
!
	subroutine openjumpstartfile2read
		use io
		use nodeinfo
		use bmpi_mod
		
		implicit none
		integer :: ilast
		logical :: found
		integer :: ierr
		
		if(nproc ==1 )return
		if(iproc == 0)then
		
		ilast = index(outfile,' ')-1
		
		inquire(file=outfile(1:ilast)//".jumpstart.bigstick",exist=found)
		
		if(found)then
			open(unit=jumpstartfile,file=outfile(1:ilast)//".jumpstart.bigstick",action='READ', access='stream', & 
			    form='unformatted',status='old',err=3034)
				write(6,*)' ........................... '
				write(6,*)' Jumpstart file found and opened!     '
				write(6,*)' ........................... '
				write(logfile,*)' ........................... '
				write(logfile,*)' Jumpstart file found  and opened!    '
				write(logfile,*)' ........................... '	
				readjumpstart = .true.
				goto 3000				
3034            continue
				
				write(6,*)' ....................... '
				write(6,*)' Jumpstart file found, but error in opening '
				write(6,*)' File is ',outfile(1:ilast)//".jumpstart.bigstick"
				write(6,*)' You must delete this file, it is confusing me '
				write(6,*)' Aborting run  '			
				write(6,*)' ....................... '
				write(logfile,*)' ....................... '
				write(logfile,*)' Jumpstart file found, but error in opening '
				write(logfile,*)' File is ',outfile(1:ilast)//".jumpstart.bigstick"
				write(logfile,*)' You must delete this file, it is confusing me '
				write(logfile,*)' Aborting run  '			
				write(logfile,*)' ....................... '
#ifdef _MPI	
				call BMPI_Abort(MPI_COMM_WORLD,101,ierr)
#endif				
				stop  
				
		else
			readjumpstart = .false.
			write(6,*)' ........................... '
			write(6,*)' No jumpstart file found     '
			write(6,*)' ........................... '
			write(logfile,*)' ........................... '
			write(logfile,*)' No jumpstart file found     '
			write(logfile,*)' ........................... '			
		end if
	    end if
3000    continue	

#ifdef _MPI	
		call BMPI_Bcast(readjumpstart,1,0,MPI_COMM_WORLD,ierr)
#endif		
		return
		
	end subroutine openjumpstartfile2read	
!====================================================================================================	
subroutine writeallinfo2jumpstart
    use system_parameters
    use sporbit
    use spstate
    use w_info
    use basis
	use fragments
	use nodeinfo
	use program_info
	use flagger
	use opbundles
	implicit none
	integer :: it
	integer :: n
	integer(4) :: v
	integer :: nob_draft,iob
	if(iproc==0)then
		write(jumpstartfile)version
		write(jumpstartfile)nproc
		write(jumpstartfile)num_threads_global
 
        !------------- WRITE INFO ON VALENCE PARTICLES
        !              PAY ATTENTION TO P-H CONJUGATION
        !        
        if(phconj(1))np(1) = -np(1)
        if(phconj(2))np(2) = -np(2)
        ! write(filenumber)np(1),np(2)
		write(jumpstartfile)np(1),np(2)
        if(phconj(1))np(1) = -np(1)
        if(phconj(2))np(2) = -np(2)
 
        !------------ WRITE OUT INFORMATION ON S.P. ORBITS 

		write(jumpstartfile)isoflag
		write(jumpstartfile)numorb(1),numorb(2)

        do it = 1,2
          do n = 1,numorb(it)
          	 ! all fields are integers
            !  write(filenumber)orbqn(it,n)%nr,orbqn(it,n)%j,orbqn(it,n)%l,  & 
            !            orbqn(it,n)%par, orbqn(it,n)%w
			write(jumpstartfile)orbqn(it,n)%nr,orbqn(it,n)%j,orbqn(it,n)%l,orbqn(it,n)%par,orbqn(it,n)%w
          enddo
        enddo
 
        !-------------- WRITE INFO ON S.P. STATES   -- CWJ 11/09  
        !               fixes a subtle bug when reading in without w-cuts
          write(jumpstartfile)nsps(1),nsps(2)
         do it = 1,2
           do n = 1,nsps(it)
             !  write(filenumber) spsqn(it,n)%nr, spsqn(it,n)%j, spsqn(it,n)%m, & 
             !             spsqn(it,n)%l, spsqn(it,n)%w, spsqn(it,n)%par, & 
             !             spsqn(it,n)%orb, spsqn(it,n)%group
			 write(jumpstartfile)spsqn(it,n)%nr,spsqn(it,n)%j,spsqn(it,n)%m,spsqn(it,n)%m, & 
			 spsqn(it,n)%l, spsqn(it,n)%w,spsqn(it,n)%par,spsqn(it,n)%orb,spsqn(it,n)%group
           enddo
        enddo
 
        !------------- WRITE OUT INFORMATION ON JZ, PARITY, W
        !  write(filenumber)jz,cparity,maxWtot
        !  write(filenumber)allsameparity, allsameW,spinless
		write(jumpstartfile)jz
        v = ICHAR(cparity)
		write(jumpstartfile)v
		write(jumpstartfile)maxWtot,allsameparity,allsameW,spinless
		write(jumpstartfile)dimbasis		
		
!.............. INFORMATION BEYOND THAT IN THE HEADER....
        write(jumpstartfile)maxfragmentsize
		write(jumpstartfile)blockflag
		write(jumpstartfile)blocksize
!............... NOW WRITE ALL THE OTHER INFORMATION........
        nob_draft = 0
        do iob =1, nopbundles
	       nob_draft = max(nob_draft,opbundle(iob)%origin)
        end do
		write(jumpstartfile)nob_draft

        do iob = 1,nob_draft
	        write(jumpstartfile)draft_opbundle(iob)%optype
			write(jumpstartfile)draft_opbundle(iob)%opwt
        end do
				
		
	end if	
	
	
	return
	
	
end subroutine writeallinfo2jumpstart


!====================================================================================================	
subroutine readbaseinfo4jumpstart
    use system_parameters
    use sporbit
    use spstate
    use w_info
    use basis
	use fragments
	use nodeinfo
	use program_info
	use flagger
	use opbundles
	use io
	use bmpi_mod
	use butil_mod
	implicit none
	integer :: it
	integer :: n
	integer(4) :: v
	integer :: nob_draft,iob
	integer :: testnproc,testnumthreads
	logical :: problem
	integer :: ierr,aerr
	character(6) :: xversion
	character(3) :: optype
	
	if(iproc==0)then
		read(jumpstartfile)xversion	
		if(xversion/=version)then
		   write(6,*)' Wrong version of Bigstick ',xversion
		   write(6,*)' Must match ',version
		   write(logfile,*)' Wrong version of Bigstick ',xversion
		   write(logfile,*)' Must match ',version		
		   print*,' Giving up!  '
		   problem = .true.
		   goto 1011
	    end if
		read(jumpstartfile)testnproc
		read(jumpstartfile)testnumthreads
		
 	   if(testnproc/=nproc )then
 		   write(6,*)' Mismatch in # MPI processes ',testnproc,' should be ',nproc
		   
 		   write(logfile,*)' Mismatch in # MPI processes ',testnproc,' should be ',nproc
 	       print*,' Giving up!  '
		   
 		   problem = .true.
		   goto 1011
 	   end if
 	   if( testnumthreads/=num_threads_global)then
 		   write(6,*)' Mismatch in # OpenMP threads ',testnumthreads,' should be ',num_threads_global
 		   write(logfile,*)' Mismatch in # OpenMP threads ',testnumthreads,' should be ',num_threads_global
 	       print*,' Giving up!  '
		   
 		   problem = .true.
		   goto 1011
 	   end if
 
        !------------- WRITE INFO ON VALENCE PARTICLES
        !              PAY ATTENTION TO P-H CONJUGATION
        !        

		read(jumpstartfile)np(1),np(2)

 
        !------------ WRITE OUT INFORMATION ON S.P. ORBITS 

		read(jumpstartfile)isoflag
		read(jumpstartfile)numorb(1),numorb(2)
		
	end if
#ifdef _MPI	
    call BMPI_BCAST(np(1), 1, 0, MPI_COMM_WORLD, ierr)
    call BMPI_BCAST(np(2), 1, 0, MPI_COMM_WORLD, ierr)
    call BMPI_BCAST(numorb(1), 1, 0, MPI_COMM_WORLD, ierr)
    call BMPI_BCAST(numorb(2),1,0, MPI_COMM_WORLD, ierr)
    call BMPI_BCAST(isoflag, 1, 0, MPI_COMM_WORLD, ierr)
#endif
    numorbmax = bmax(numorb(1),numorb(2))
	
    if(.not.allocated(orbqn)) then
       allocate ( orbqn(2,numorbmax), stat=aerr )
       if(aerr /= 0) call memerror("read_wfn_header 1")
    end if
    do it = 1,2
        do n = 1,numorb(it)
	
	       if(iproc==0)then

          	 ! all fields are integers
            !  write(filenumber)orbqn(it,n)%nr,orbqn(it,n)%j,orbqn(it,n)%l,  & 
            !            orbqn(it,n)%par, orbqn(it,n)%w
			read(jumpstartfile)orbqn(it,n)%nr,orbqn(it,n)%j,orbqn(it,n)%l,orbqn(it,n)%par,orbqn(it,n)%w
		   end if
#ifdef _MPI	
           call BMPI_BCAST(orbqn(it,n)%nr,1,0,MPI_COMM_WORLD,ierr)
           call BMPI_BCAST(orbqn(it,n)%j ,1,0,MPI_COMM_WORLD,ierr)
           call BMPI_BCAST(orbqn(it,n)%l ,1,0,MPI_COMM_WORLD,ierr)
           call BMPI_BCAST(orbqn(it,n)%par,1,0,MPI_COMM_WORLD,ierr)
           call BMPI_BCAST(orbqn(it,n)%w ,1,0,MPI_COMM_WORLD,ierr)
#endif		
			
        enddo
    enddo
 
        !-------------- WRITE INFO ON S.P. STATES   -- CWJ 11/09  
        !               fixes a subtle bug when reading in without w-cuts
   if(iproc==0)then
      read(jumpstartfile)nsps(1),nsps(2)
   end if
#ifdef _MPI	
   call BMPI_BCAST(nsps,2,0,MPI_COMM_WORLD,ierr)
#endif
   if(.not.allocated(spsqn)) then
      allocate(spsqn(2,bmax(nsps(1),nsps(2))), stat=aerr)
      if(aerr /= 0) call memerror("read_wfn_header 2")
   end if
   ! print *, "Abount to read spsqn"
   do it = 1,2
      do n = 1,nsps(it)
        if(iproc==0) then

             !  write(filenumber) spsqn(it,n)%nr, spsqn(it,n)%j, spsqn(it,n)%m, & 
             !             spsqn(it,n)%l, spsqn(it,n)%w, spsqn(it,n)%par, & 
             !             spsqn(it,n)%orb, spsqn(it,n)%group
			 read(jumpstartfile)spsqn(it,n)%nr,spsqn(it,n)%j,spsqn(it,n)%m,spsqn(it,n)%m, & 
			 spsqn(it,n)%l, spsqn(it,n)%w,spsqn(it,n)%par,spsqn(it,n)%orb,spsqn(it,n)%group
		 end if
#ifdef _MPI	
         call BMPI_BCAST(spsqn(it,n)%nr,1,0,MPI_COMM_WORLD,ierr)
         call BMPI_BCAST(spsqn(it,n)%j ,1,0,MPI_COMM_WORLD,ierr)
         call BMPI_BCAST(spsqn(it,n)%m ,1,0,MPI_COMM_WORLD,ierr)
         call BMPI_BCAST(spsqn(it,n)%l ,1,0,MPI_COMM_WORLD,ierr)
         call BMPI_BCAST(spsqn(it,n)%w,1,0,MPI_COMM_WORLD,ierr)
         call BMPI_BCAST(spsqn(it,n)%par,1,0,MPI_COMM_WORLD,ierr)
         call BMPI_BCAST(spsqn(it,n)%orb,1,0,MPI_COMM_WORLD,ierr)
         call BMPI_BCAST(spsqn(it,n)%group,1,0,MPI_COMM_WORLD,ierr)
#endif
	 end do 
   end do
   
   if(iproc==0)then
 
        !------------- WRITE OUT INFORMATION ON JZ, PARITY, W
        !  write(filenumber)jz,cparity,maxWtot
        !  write(filenumber)allsameparity, allsameW,spinless
		read(jumpstartfile)jz
		read(jumpstartfile)v
		read(jumpstartfile)maxWtot,allsameparity,allsameW,spinless
		read(jumpstartfile)dimbasischeck		
		
!.............. INFORMATION BEYOND THAT IN THE HEADER....
        read(jumpstartfile)maxfragmentsize
		read(jumpstartfile)blockflag
		read(jumpstartfile)blocksize				
		
	end if	
#ifdef _MPI	
    call BMPI_BCAST(jz,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(maxWtot,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(v,1,0,MPI_COMM_WORLD,ierr)
#endif
!	if(iproc==0)print*,' my parity ',cparity
    cparity = char(v)
    select case (cparity)
      case ('+','0')
         iparity = 1
      case ('-')
         iparity = 2
      case default
         if(iproc==0)write(6,*)' Problem with parity ',cparity
         stop  
    end select
#ifdef _MPI		
    call BMPI_BCAST(allsameparity,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(allsameW,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(spinless,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(dimbasischeck,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(maxfragmentsize,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(blockflag,1,0,MPI_COMM_WORLD,ierr)
    call BMPI_BCAST(blocksize,1,0,MPI_COMM_WORLD,ierr)
#endif		
	return

1011 continue
#ifdef _MPI	
  call BMPI_Abort(MPI_COMM_WORLD,321,ierr)	
#endif
	
end subroutine readbaseinfo4jumpstart
!====================================================================================================
subroutine readtimeinfo4jumpstart(nob_draft)
	use opbundles
	use io
	use bmpi_mod
	use butil_mod
	use nodeinfo
	implicit none
	integer :: nob_draft,iob,nob_test
	logical :: problem
	integer :: ierr,aerr
	character(3) :: optype
	real(4) :: opwt
	
	
	
	if(iproc==0)then
		
		print*,' ABOUT TO READ IN OP WTS'
		read(jumpstartfile)nob_test
	end if
#ifdef _MPI	
    call BMPI_BCAST(nob_test,1,0,MPI_COMM_WORLD,ierr)
#endif
	if(nob_draft /= nob_test)then
		if(iproc==0)print*,' Mismatch in # draft opbundles ',nob_draft,nob_test
		goto 1011
	end if
	

    do iob = 1,nob_draft
		problem = .false.
		if(iproc==0)then
	        read(jumpstartfile)optype
		    if(optype /= draft_opbundle(iob)%optype)then
				write(6,*)' Mismatch in draft optype ',iob,optype, ' expect ',draft_opbundle(iob)%optype
				write(logfile,*)' Mismatch in draft optype ',iob,optype, ' expect ',draft_opbundle(iob)%optype
				
				problem = .true.
				
			end if
			read(jumpstartfile)opwt
		end if
#ifdef _MPI	
	    call BMPI_BCAST(problem,1,0,MPI_COMM_WORLD,ierr)
#endif
		if(problem)goto 1011
#ifdef _MPI			
	    call BMPI_BCAST(opwt,1,0,MPI_COMM_WORLD,ierr)
#endif
		draft_opbundle(iob)%opwt =opwt
			
			
    end do
			
	return

1011 continue
#ifdef _MPI	
	call BMPI_Abort(MPI_COMM_WORLD,322,ierr)		
#endif

end	subroutine readtimeinfo4jumpstart

!====================================================================================================	

!
!  subroutine which writes information on the basis, and on jumps, to the jumpstartfile
!  NOTE: Must use specialized routines to write in MPI when using fragments
!
subroutine write_info2jumpstart1
	
	use fragments
	use nodeinfo
	use program_info
	use flagger
	use wfn_mod
!	use sectorjumps
	
    implicit none
	integer numthreads
	
	if(iproc /=0)return
	
!....... FIRST WRITE OUT VERSION OF BIGSTICK AS A CHECK........

!    write(jumpstartfile)version	

!	write(jumpstartfile)nproc
	call wfn_write_int4(jumpstartfile, nproc)

!	write(jumpstartfile)num_threads_global
	call wfn_write_int4(jumpstartfile, num_threads_global)
	call wfn_write_int8(jumpstartfile, maxfragmentsize)
 
!    write(jumpstartfile)maxfragmentsize
	
    return
end	subroutine write_info2jumpstart1
!====================================================================================================	
!
!  subroutine which writes information on the basis to the jumpstartfile
!
subroutine read_info4jumpstart1
	
	use fragments
	use nodeinfo
	use program_info
	use bmpi_mod
	use flagger
	use io
	use wfn_mod
!	use sectorjumps
	
    implicit none
!	character(6) :: xversion
	integer :: testnproc,testnumthreads
	logical :: problem
	integer :: ierr
	
    problem = .false.
    if (iproc==0)then
!        read(jumpstartfile)xversion	
!	    if(xversion/=version)then
!		   write(6,*)' Wrong version of Bigstick ',xversion
!		   write(6,*)' Must match ',version
!	 	   write(logfile,*)' Wrong version of Bigstick ',xversion
!		   write(logfile,*)' Must match ',version		
!	       print*,' Giving up!  '
!		   problem = .true.
!	   end if
	   
   	   call wfn_read_int4(jumpstartfile, testnproc)

   !	write(jumpstartfile)num_threads_global
   	   call wfn_read_int4(jumpstartfile,testnumthreads)
	 !  read(jumpstartfile)testnproc
	 !  read(jumpstartfile)testnumthreads
	   
	   if(testnproc/=nproc )then
		   write(6,*)' Mismatch in # MPI processes ',testnproc,' should be ',nproc
		   
		   write(logfile,*)' Mismatch in # MPI processes ',testnproc,' should be ',nproc
	       print*,' Giving up!  '
		   
		   problem = .true.
	   end if
	   if( testnumthreads/=num_threads_global)then
		   write(6,*)' Mismatch in # OpenMP threads ',testnumthreads,' should be ',num_threads_global
		   write(logfile,*)' Mismatch in # OpenMP threads ',testnumthreads,' should be ',num_threads_global
	       print*,' Giving up!  '
		   
		   problem = .true.
	   end if
	   
   end if    
#ifdef _MPI	
   call BMPI_Bcast(problem,1,0,MPI_COMM_WORLD,ierr)
#endif
   if(problem)then
#ifdef _MPI		   
		call BMPI_Abort(MPI_COMM_WORLD,133,ierr)	
#endif 
		stop
	end if
	
	if(iproc==0)call wfn_read_int8(jumpstartfile, maxfragmentsize)
#ifdef _MPI	
    call BMPI_Bcast(maxfragmentsize,1,0,MPI_COMM_WORLD,ierr)
#endif	
		
	return
	
end	subroutine read_info4jumpstart1

!====================================================================================================
!
! routine to write out the data on timing
!
subroutine write_opwt2jumpstart
	
	use opbundles
	use nodeinfo
	use wfn_mod
	implicit none
	
	integer :: iob,myob
	integer(4) :: opint
	integer :: nob_draft
	
!....... FIRST CHECK # OF DRAFT OPBUNDLES	
    if(iproc /=0 )return
    nob_draft = 0
	do iob =1, nopbundles
		nob_draft = max(nob_draft,opbundle(iob)%origin)
	end do
	
	do iob = 1,nob_draft
		call convert_optype2int(draft_opbundle(iob)%optype,opint)
		call wfn_write_int4(jumpstartfile, opint)
!............. ADDED 7.9.10 to try to improve jumpstarting........		
!		call wfn_write_int4(jumpstartfile,splinter(iob)%nsplits)		
!		if(splinter(iob)%nsplits>0)then
!			call wfn_write_int8(jumpstartfile,splinter(iob)%pxstart)
!			call wfn_write_int8(jumpstartfile,splinter(iob)%pxend)
!			call wfn_write_int8(jumpstartfile,splinter(iob)%nxstart)
!			call wfn_write_int8(jumpstartfile,splinter(iob)%nxstart)
!			do myob= 1,splinter(iob)%nsplits				
!				call wfn_write_int8(jumpstartfile,splinter(iob)%locsplit(myob))
!				call wfn_write_real4(jumpstartfile,splinter(iob)%opwtsplit(myob))				
!			end do
!			myob =splinter(iob)%nsplits+1
!			call wfn_write_real4(jumpstartfile,splinter(iob)%opwtsplit(myob))					
!		else    !.... ORIGINAL JUMPSTART WRITE OUT
		   call wfn_write_real4(jumpstartfile, draft_opbundle(iob)%opwt)		
!	    end if
	end do
	return

end subroutine write_opwt2jumpstart

!====================================================================================================
!
! routine to read in the data on timing
!
subroutine read_opwt4jumpstart(nob_draft)
	
	use opbundles
	use nodeinfo
	use io
	use bmpi_mod
	use wfn_mod
	implicit none
	
	integer :: iob
	
	integer :: nob_draft
	character(3) :: optype
	real(4) :: opwt
	real, allocatable :: opwtlist(:)
	integer(4) :: opint
	logical problem
	integer :: ierr
	
!....... FIRST CHECK # OF DRAFT OPBUNDLES	


    allocate(opwtlist(nob_draft))
    problem = .false.
    if(iproc==0)then
		print*,' READING IN JUMP START OP WEIGHTS'
	    do iob = 1,nob_draft
			call wfn_read_int4(jumpstartfile, opint)
			call convert_opint2type(opint,optype)
		    if(optype /= draft_opbundle(iob)%optype)then
				write(6,*)' Mismatch in draft optype ',iob,optype, ' expect ',draft_opbundle(iob)%optype
				write(logfile,*)' Mismatch in draft optype ',iob,optype, ' expect ',draft_opbundle(iob)%optype
				
				problem = .true.
				exit
		    end if
			call wfn_read_int4(jumpstartfile,splinter(iob)%nsplits )
			
			
			call wfn_read_real4(jumpstartfile, opwt)
			
			opwtlist(iob)=opwt
!			draft_opbundle(iob)%opwt = opwt
			
  	    end do
	end if
#ifdef _MPI	
    call BMPI_Bcast(problem,1,0,MPI_COMM_WORLD,ierr)
#endif
    if(problem)then
#ifdef _MPI	
 		call BMPI_Abort(MPI_COMM_WORLD,134,ierr)	
#endif
 		stop
 	end if	
#ifdef _MPI		
	call BMPI_Bcast(opwtlist,size(opwtlist),0,MPI_COMM_WORLD,ierr)
#endif
    do iob = 1,nob_draft

		draft_opbundle(iob)%opwt = opwtlist(iob)
		
    end do
	
	deallocate(opwtlist)
	return

end subroutine read_opwt4jumpstart

!===============================================================
!
! routine to read in the data on timing
!
subroutine read_opwt4jumpstart_OLD(nob_draft)
	
	use opbundles
	use nodeinfo
	use io
	use bmpi_mod
	use wfn_mod
	implicit none
	
	integer :: iob
	
	integer :: nob_draft
	character(3) :: optype
	real(4) :: opwt
	real, allocatable :: opwtlist(:)
	integer(4) :: opint
	logical problem
	integer :: ierr
	
!....... FIRST CHECK # OF DRAFT OPBUNDLES	


    allocate(opwtlist(nob_draft))
    problem = .false.
    if(iproc==0)then
		print*,' READING IN JUMP START OP WEIGHTS'
	    do iob = 1,nob_draft
			call wfn_read_int4(jumpstartfile, opint)
			call convert_opint2type(opint,optype)
		    if(optype /= draft_opbundle(iob)%optype)then
				write(6,*)' Mismatch in draft optype ',iob,optype, ' expect ',draft_opbundle(iob)%optype
				write(logfile,*)' Mismatch in draft optype ',iob,optype, ' expect ',draft_opbundle(iob)%optype
				
				problem = .true.
				exit
		    end if
			call wfn_read_real4(jumpstartfile, opwt)
			
			opwtlist(iob)=opwt
!			draft_opbundle(iob)%opwt = opwt
			
  	    end do
	end if
#ifdef _MPI	
    call BMPI_Bcast(problem,1,0,MPI_COMM_WORLD,ierr)
#endif
    if(problem)then
#ifdef _MPI	
 		call BMPI_Abort(MPI_COMM_WORLD,134,ierr)	
#endif
 		stop
 	end if	
#ifdef _MPI	
	call BMPI_Bcast(opwtlist,size(opwtlist),0,MPI_COMM_WORLD,ierr)
#endif
    do iob = 1,nob_draft

		draft_opbundle(iob)%opwt = opwtlist(iob)
		
    end do
	
	deallocate(opwtlist)
	return

end subroutine read_opwt4jumpstart_OLD

!===============================================================
subroutine convert_optype2int(optype,opint)
	implicit none
	character(3),intent(IN) :: optype
	integer(4),intent(OUT) :: opint
	
	select case (optype)
	
	
	case ('SPE')
	opint=0
	
	case('PP','PP0')
	opint=1
	
	case('NN','NN0')
	opint=2
	
	case('PN','PN0')
	opint=3
	
	
	case('PPP')
	opint=4
	
	case('NNN')
	opint=5
	
	case('PPN')
	opint=6	
	
	case('PNN')
	opint=7
	
	case default
	
	print*,' OOPS WRONG OPTYPE ',optype
	stop
	
    end select
    return
end subroutine convert_optype2int
!===============================================================
subroutine convert_opint2type(opint,optype)
	implicit none
	character(3),intent(OUT) :: optype
	integer(4),intent(IN) :: opint
	
	select case (opint)
	
	
	case (0)
	optype='SPE'
	
	case(1)
	optype='PP'
	
	case(2)
	optype='NN'
	
	case(3)
	optype='PN'
	
	
	case(4)
	optype='PPP'
	
	case(5)
	optype='NNN'
	
	case(6)
	optype='PPN'
	
	case(7)
	optype='PNN'
	
	case default
	
	print*,' OOPS WRONG OPINT ',opint
	stop
	
    end select
    return
end subroutine convert_opint2type



!====================================================================================================
end module jumpstart
