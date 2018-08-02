!
!  added in 7.8.3
!  routines to survey and 'annex' bundles
!  Begin by finding 'neighbors'--opbundles which have some jumps in common and others which are 'next door'
!

module annexations
	use opbundles
	implicit none
	
	logical :: annexneutronneighborsfirst = .true.
	logical :: annexation_is_go=.false.  
	
contains
	
!==========================================================
!
! THIS ROUTINE INITIALIZES ALL OPBUNDLES AS 'UNANNEXED'
!	
	subroutine set_unannexed
		use nodeinfo
		implicit none
		integer :: iob
		
		do iob = opbundlestart(iproc),opbundleend(iproc)
			opbundle(iob)%annexed=.false.
		end do
		
		return
	end subroutine set_unannexed


!=============================================================

subroutine annexleader
	use nodeinfo
	use sporbit
	use bmpi_mod
	implicit none
	integer, allocatable :: annexn(:),annexp(:)
	integer :: jproc
	integer :: ierr
	
	call set_unannexed
	
	if(.not.annexation_is_go)return
	if(allsameW)return
	
	if(iproc==0)then
		open(unit=35,file='annexinfo.bigstick',status='unknown')
	end if
	
	allocate(annexn(0:nproc-1),annexp(0:nproc-1))
	annexn = 0
	annexp = 0
	
	if(annexneutronneighborsfirst )then
		call annex_those_bundles('N',annexn(iproc))
		call annex_those_bundles('P',annexp(iproc))
	else
		call annex_those_bundles('P',annexp(iproc))		
		call annex_those_bundles('N',annexn(iproc))		
	end if
	
	call BMPI_ALLREDUCE(annexn,size(annexn),MPI_SUM,icomm,ierr)
	call BMPI_ALLREDUCE(annexp,size(annexn),MPI_SUM,icomm,ierr)
	
	if(iproc==0)then
     	write(35,*)' Proc     neutron annex      proton annex      total opbundles '
		do jproc= 0,nproc-1
			write(35,'(i6,3i18)')jproc,annexn(jproc),annexp(jproc),-opbundlestart(jproc)+1+opbundleend(jproc)
		end do
		close(35)
			
	end if
		
	deallocate(annexn,annexp)
	return
	
end subroutine annexleader
!=============================================================


subroutine annex_those_bundles(dir,nn)
	use nodeinfo
	implicit none
	character,intent(in) :: dir  ! = 'N' or 'P'
	integer :: iob,job
	integer,intent(out) :: nn
	
	nn = 0
	
	do iob = opbundlestart(iproc),opbundleend(iproc)-1
		if(opbundle(iob)%annexed)cycle
		
		do job = iob+1,opbundleend(iproc)
			if(opbundle(job)%annexed)cycle
			if(opbundle(iob)%optype /=opbundle(job)%optype)cycle
			if(opbundle(iob)%hchar  /=opbundle(job)%hchar)cycle
			
			if(dir=='N')then
			    if(opbundle(iob)%pxstart/=opbundle(job)%pxstart)cycle
			    if(opbundle(iob)%pxend  /=opbundle(job)%pxend)cycle				
				if(opbundle(iob)%nxend+1 == opbundle(job)%nxstart )then  ! annex !
					opbundle(iob)%nxend = opbundle(job)%nxend  ! I think this is the only information that 
																! needs to be modified
					opbundle(job)%annexed=.true.
					nn=nn+1
				end if
				if(opbundle(job)%nxend+1 == opbundle(iob)%nxstart )then  ! annex !
					opbundle(iob)%nxstart = opbundle(job)%nxstart  ! I think this is the only information that 
																! needs to be modified
					opbundle(job)%annexed=.true.
					nn=nn+1
				end if
				
			else
			    if(opbundle(iob)%nxstart/=opbundle(job)%nxstart)cycle
			    if(opbundle(iob)%nxend  /=opbundle(job)%nxend)cycle				
				if(opbundle(iob)%pxend+1 == opbundle(job)%pxstart )then  ! annex !
					opbundle(iob)%pxend = opbundle(job)%pxend  ! I think this is the only information that 
																! needs to be modified
					opbundle(job)%annexed=.true.
					nn=nn+1
				end if
				if(opbundle(job)%pxend+1 == opbundle(iob)%pxstart )then  ! annex !
					opbundle(iob)%pxstart = opbundle(job)%pxstart  ! I think this is the only information that 
																! needs to be modified
					opbundle(job)%annexed=.true.
					nn=nn+1
				end if				
				
			end if
			
		end do ! job
	
	end do ! iob
	
!	print*,' proc ',iproc,', ',nn,' annexations of type ',dir
	
	return
	
end subroutine annex_those_bundles



end module annexations
