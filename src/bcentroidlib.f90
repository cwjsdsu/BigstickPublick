!
!  routines for computing and applying 3-body centroids
!  started in 7.9.7
!  many routines prototyped in tracer.x code, V8
!
!  TO APPLY TO ANY BASIS STATE
!
!  subroutine mapxSDxconfig, proton/neutron SDs to a configuration
!  which are stored in pmapSD2config and nmapSD2config in module configtools
!  and this is done using pcfindx
!

module threebodycentroids
	use sporbit
	use flags3body
	
	implicit none
	
	logical :: use3body
	real, allocatable,target :: W3ppp(:,:,:),  W3ppn(:,:,:), W3pnn(:,:,:),W3nnn(:,:,:)
	

!......... HOW TO FIND A CENTROID:
!     for state pstate, nstate,
!  centroid( pcfindx(pmapSD2config(pstate)+ nmaxSD2config(nstate))	
	real, allocatable :: centroid(:)
	logical :: enable_centroids = .false.
	
contains

!
!
!=======================================================	
!
!  CALLED BY:
!
subroutine read_and_fill_3bodycentroidpots
   	    use nodeinfo
		use bmpi_mod
		use io
		use reporter
		implicit none
		character ychar
		character*80 :: filename
		integer :: ilast
		logical :: success
		
		integer(8) :: ime
		integer :: a,b,c !,r,s,u
		integer :: mymaxorb
		real :: c1,c2,c3,c4,c5,c6
		integer :: d1,d2,d3,d4,d5,d6 ! not needed but included
		integer :: i
		
		real, allocatable :: tmpvec(:)
		
!		real :: xnorm
!		logical :: sameab
		
!		integer :: Jabmin,Jabmax,Jmin,Jmax
		integer :: ierr
		
		character(40) :: title
		
				
        if(iproc==0)then
		   print*,' '
		   success = .false.
		   do while (.not.success)
		      print*,' Enter full name of file containing three-body monopole potentials '
			  
			  if(auto_input)then
				  read(autoinputfile,'(a)')filename
			  else
			     read(5,'(a)')filename
				 write(autoinputfile,'(a)')filename
			  end if
			  
			  ilast = index(filename,' ')-1
			  if(filename(1:ilast)=='none')then
				  applycentroids=.false.
				  use3body =.false.
				  exit
			  end if
			  open(unit=11,file=filename(1:ilast),status='old',err=101)
			  success = .true.
			  print*,' Opened successfully ',filename(1:ilast)
			  
			  write(logfile,*)' Reading 3-body monopole potential from ',filename(1:ilast)
			  exit
			 
101           continue
              print*,filename(1:ilast),' does not appear to exist'
		
   	       end do ! while not success
		   use3body=.true.

	    end if
#ifdef _MPI
		call BMPI_BCAST(use3body,1,0,MPI_COMM_WORLD,ierr)
#endif		
		if(use3body)then
			applycentroids=.true.
		else
			applycentroids=.false.
			return
		end if
		
		numorbmax=max(numorb(1),numorb(2))
		
		allocate(W3ppp(numorb(1),numorb(1),numorb(1) ))
		allocate(W3ppn(numorb(1),numorb(1),numorb(2) ))
		allocate(W3pnn(numorb(1),numorb(2),numorb(2) ))
		allocate(W3nnn(numorb(2),numorb(2),numorb(2) ))		

		W3ppp = 0.000
		W3ppn = 0.000
		W3pnn = 0.000
		W3nnn = 0.000			
		
		
		if(iproc==0)then
		do ime = 1,10
			read(11,'(a)')title
			if(title(1:1)=='#')then
				write(6,'(a)')title
				write(logfile,'(a)')title
			else
				backspace(11)
				exit
			end if
		end do
		
		mymaxorb = 0
		do ime = 1,1000000000
			! note  actually 2Jab 2Jrs 2J Tab Trs   2T
			read(11,*,end=202)a,b,c,c1,c2,c3,c4,c5,c6,d1,d2,d3,d4,d5,d6
			
			mymaxorb=max(mymaxorb,a)
			mymaxorb=max(mymaxorb,b)
			mymaxorb=max(mymaxorb,c)
			
			if(a> numorbmax .or. b > numorbmax .or. c > numorbmax)cycle
			
			W3ppp(a,b,c)=c1
			W3nnn(a,b,c)=c1
			W3ppn(a,b,c)=c2
			W3pnn(c,a,b)=c2
			W3ppn(a,c,b)=c3
			W3ppn(c,a,b)=c3
			W3pnn(b,a,c)=c3
			W3ppn(b,c,a)=c4
			W3pnn(a,b,c)=c4

		end do
		
		print*,' Oops, did not reach end of file '
		stop
		
202     continue

        print*,ime,' 3-body matrix centroid pots read in'
		print*,' largest orbit read = ',mymaxorb,numorbmax
		close (11)
	    end if
        if(nproc==1)return
#ifdef _MPI		
		call BMPI_barrier(MPI_COMM_WORLD,ierr)
#endif		

		
!... NOW BROADCAST
        mymaxorb = min(mymaxorb,numorbmax)
#ifdef _MPI		
		call BMPI_Bcast(mymaxorb,1,0,MPI_COMM_WORLD,ierr)
		call BMPI_Bcast(numorbmax,1,0,MPI_COMM_WORLD,ierr)
#endif		
		
		allocate(tmpvec(mymaxorb*(mymaxorb+1)*(mymaxorb+2)/6))
		
		if(iproc==0)then
		i = 0
		do a = 1,mymaxorb
			do b = 1,a
				do c = 1,b
					i = i + 1
					tmpvec(i)=W3ppp(a,b,c)
				end do
			end do
		end do
     	end if
		
#ifdef _MPI
        call BMPI_BCAST(tmpvec,size(tmpvec),0,MPI_COMM_WORLD,ierr)
#endif
        if(iproc/=0)then
		   i = 0
		   do a = 1,mymaxorb
			  do b = 1,a
				 do c = 1,b
					i = i + 1
					W3ppp(a,b,c)=tmpvec(i)
				  end do
			  end do
		   end do	
	    end if	
#ifdef _MPI
		call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif		
		tmpvec = 0.0
		i = 0
		if(iproc==0)then
		do a = 1,mymaxorb
			do b = 1,a
				do c = 1,b
					i = i + 1
					tmpvec(i)=W3nnn(a,b,c)
				end do
			end do
		end do
     	end if
#ifdef _MPI		
        call BMPI_BCAST(tmpvec,size(tmpvec),0,MPI_COMM_WORLD,ierr)
#endif
        if(iproc/=0)then
		   i = 0
		   do a = 1,mymaxorb
			  do b = 1,a
				 do c = 1,b
					i = i + 1
					W3nnn(a,b,c)=tmpvec(i)
				  end do
			  end do
		   end do	
	    end if	
#ifdef _MPI
		call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif		
		tmpvec = 0.0	
		deallocate(tmpvec)
		allocate(tmpvec(mymaxorb**3))
		
		
		i = 0
		do a = 1,mymaxorb
			do b = 1,mymaxorb
				do c = 1,mymaxorb
					i = i + 1
					tmpvec(i)=W3ppn(a,b,c)
				end do
			end do
		end do
#ifdef _MPI		
        call BMPI_BCAST(tmpvec,size(tmpvec),0,MPI_COMM_WORLD,ierr)
#endif
        if(iproc/=0)then
		   i = 0
		   do a = 1,mymaxorb
			  do b = 1,mymaxorb
				 do c = 1,mymaxorb
					i = i + 1
					W3ppn(a,b,c)=tmpvec(i)
				  end do
			  end do
		   end do	
	    end if	
#ifdef _MPI
		call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
		tmpvec = 0.0	
		
		i = 0
		do a = 1,mymaxorb
			do b = 1,mymaxorb
				do c = 1,mymaxorb
					i = i + 1
					tmpvec(i)=W3pnn(a,b,c)
				end do
			end do
		end do
		
#ifdef _MPI
        call BMPI_BCAST(tmpvec,size(tmpvec),0,MPI_COMM_WORLD,ierr)
#endif
        if(iproc/=0)then
		   i = 0
		   do a = 1,mymaxorb
			  do b = 1,mymaxorb
				 do c = 1,mymaxorb
					i = i + 1
					W3pnn(a,b,c)=tmpvec(i)
				  end do
			  end do
		   end do	
	    end if	
#ifdef _MPI
		call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
		tmpvec = 0.0
		

        deallocate(tmpvec)
		return
		
	end subroutine read_and_fill_3bodycentroidpots
		
!=======================================================	
	subroutine e0config3b(it,orbarray,e03b)
		
		use system_parameters
		use sporbit
!		use configurations
		implicit none

		integer it
		integer orbarray( numorb(it) )
		real e03b
		integer iorb,jorb,korb
		real fijk
		real, pointer :: W3xxx(:,:,:)
		integer,pointer :: D3xxx(:,:,:)
		
		if(np(it) == 0)then
		   return
		endif
		if(it ==1)then
		    W3xxx => W3ppp
		else
		    W3xxx => W3nnn
		endif

		do iorb = 1,numorb(it)
		   if( orbarray(iorb) == 0)cycle
		   do jorb = iorb,numorb(it)
		      if(orbarray(jorb) == 0)cycle
			  
			  do korb = jorb,numorb(it)
			      if(orbarray(korb) == 0)cycle
				  
			  
				  fijk = float(orbarray(korb)*orbarray(jorb)*orbarray(iorb))
                  if(iorb==jorb .and. iorb ==korb)then
					  fijk = orbarray(iorb)*(orbarray(jorb)-1.)*(orbarray(korb)-2.)/6.
				  else
					  if(iorb==jorb)then
						  fijk = orbarray(korb)*orbarray(iorb)*(orbarray(jorb)-1.)/2.
					  end if
					  if(korb==jorb)then
						  fijk = orbarray(iorb)*orbarray(jorb)*(orbarray(korb)-1.)/2.
					  end if					  
				  end if
				  					  
		          e03b = e03b+ W3xxx(korb,jorb,iorb)*fijk
			  end do  !korb
				  
				  
		   end do  !jorb
		end do  !iorb
		
		
		return
	end subroutine e0config3b
	
	
	!----------------------------------------------------------------------
	!
	!  ADDS the protron-neutron 3-body contribution to the centroid
	!
	!  INPUT: porbarray( numorb(1) ): array of proton occupations
	!         norbarray( numorb(2) ): array of proton occupations
	!
	subroutine e0pn3b(porbarray,norbarray,e03b)

	use system_parameters
	use sporbit
	use interaction
!	use configurations
	implicit none

	integer porbarray( numorb(1) )
	integer norbarray( numorb(2) )
	
	real e03b
	integer iorb,jorb,korb
	real fab

    if(np(1)*np(2)==0)return

!........... FIRST DO PPN

    if(np(1) < 2)goto 1
	do iorb = 1,numorb(1)
	   if( porbarray(iorb) == 0)cycle
	   do jorb = 1,iorb
		   if( porbarray(jorb) == 0)cycle
		   
		   if(iorb==jorb)then
			   fab = porbarray(iorb)*(porbarray(jorb)-1)/2.
		   else
			   fab= float(porbarray(iorb)*porbarray(jorb))
		   end if
		   do korb = 1,numorb(2)
			   
			   e03b = e03b + W3ppn(iorb,jorb,korb)*fab*norbarray(korb)
			   
		   end do  ! korb
	   end do   !jorb
   end do ! iorb
	
1 continue	
!........... NEXT DO PNN

    if(np(2) < 2)return
	do iorb = 1,numorb(2)
	   if( norbarray(iorb) == 0)cycle
	   do jorb = 1,iorb
		   if( norbarray(jorb) == 0)cycle
		   
		   if(iorb==jorb)then
			   fab = norbarray(iorb)*(norbarray(jorb)-1)/2.
		   else
			   fab= float(norbarray(iorb)*norbarray(jorb))
		   end if
		   do korb = 1,numorb(1)
			  			   
			   e03b = e03b + W3pnn(korb,iorb,jorb)*fab*porbarray(korb)
			   
		   end do  ! korb
	   end do   !jorb
   end do ! iorb

	
	return

	end subroutine e0pn3b
!================================================================	
	
end module	threebodycentroids
	
	
	
