!
!
!  subroutine output_TRDENS
!
!  routine to print out wfns in format useful for P. Navratil's TRDENS code
!  called after routine LANCZOS has finished
!
!  MODIFICATION 3/2013: write out proton and neutron states with the same number of 
!  s.p. states
!  MODIFICATION 5/2023 v 7.10.7: option to read in prior wave function
!.........................................................................

 subroutine output_TRDENS(from_old)

 use system_parameters
 use spstate
 use sporbit
 use W_info
 use haiku_info
 use io
 use basis
 use precisions
 use sectors
 use blocks
 use coupledmatrixelements
 use nodeinfo
 use obs
 use localvectors
 use lanczos_info
 use wfn_mod
 use bmpi_mod
 use butil_mod
 use basis
 use mod_reorthog
 use haikus
 use bvectorlib_mod
 implicit none

 logical :: from_old
 real(kind=lanc_prec), allocatable :: vamp(:)
 real(kind=lanc_prec), allocatable :: vec(:,:)
 integer i,j,nstate,n
 integer ith
 integer, parameter :: filenumber = 56 
 real e,xj,xt,xtt

 integer ps,cs,ns
 integer pblock, nblock
 integer prblock,plblock
 integer(8) :: psdstart
 integer nrblock,nlblock
 integer(8) :: nsdstart
 integer npladd,npradd,pladd,pradd
 integer nnladd,nnradd,nladd,nradd
 integer(8) :: ip,in
 integer, allocatable :: pocc(:), nocc(:)
 integer(kind=basis_prec) :: ibasis
 integer,allocatable :: pocctmp(:),nocctmp(:)
 integer,pointer  :: hsd(:)

 integer :: Nmshell,wspsmin,wspsmax,it,w
 integer ilast
 integer w0, A0
 integer :: itmax   ! to make uniform the number of proton and neutron s.p. states
  integer tag
  integer(4) :: ierr
  integer stat !(MPI_STATUS_SIZE)
  integer :: fromproc
  integer :: aerr
  
  logical :: write_pnwfn

  write_pnwfn = .true.

if(iproc==0)then
  print*,' '
  print*,' Writing out .trwfn file....'
  print*,' '
  
  if(write_pnwfn)then
	  print*,' Also writing out .pnwfn '
	  print*,' '
  else
	  print*,' Note: SKIPPING .pnwfn output'
	  print*,' '
  endif
end if

!................. ADDED 7.10.7.........

if(from_old)then
	if(nproc==1)then
		storelanczosincore1 = .true.
	else
		storelanczosincoreMPI = .true.
	end if
    call overlaptribution
    call wfn_read_nkeep(oldwfnfile, nkeep) ! does BCAST
	niter = nkeep
    if(iproc==0) print*,' There are ', nkeep,' wavefunctions '
    call setup_localvectors
	
end if

!.................. ALLOCATE OCCUPATION ............

 allocate(pocc(np(1)), nocc(np(2)), stat=aerr)
 if(aerr /= 0) call memerror("output_TRDENS 1")
 allocate(pocctmp(nhsps(1)),nocctmp(nhsps(2)), stat=aerr)
 if(aerr /= 0) call memerror("output_TRDENS 2")

 if(storelanczosincore1 .or. storelanczosincoreMPI)then	 
    allocate(vamp(nkeep), stat=aerr)
    if(aerr /= 0) call memerror("output_TRDENS 3")

 else
    allocate(vec(v1s:v1e,nkeep), stat=aerr)  ! allocate to match vec1
    if(aerr /= 0) then
       call memerror("output_TRDENS 10")
       stop 5
    end if
 end if
 
!..................BASIC INFORMATION ....................

 if(np(2) > np(1))then
   itmax = 2
 else
   itmax = 1
 end if
!..............COMPUTE # OF MAJOR SHELLS
!.............. NOTE THIS ASSUMES CONTIGUOUS........
 Nmshell = 0
 wspsmin = 10000
 wspsmax = 0
 do it = 1,2
 do i = 1,nhsps(it)
    w = hspsqn(it,i)%w
!    print*,w,wceiling(it),hspsqn(i,it)%nr, hspsqn(i,it)%l
    if( wceiling(it) < w)cycle
    wspsmin = bmin(wspsmin,w)
    wspsmax = bmax(wspsmax,w)
 enddo
 enddo
 nmshell = wspsmax-wspsmin +1
!...........................COMPUTE SIZE OF CORE, IF ANY.....
 w0 = 100
 do i = 1,numorb(1)
    w0 = bmin( w0,2*orbqn(1,i)%nr + orbqn(1,i)%l)
 end do
 A0 = 0
 if(w0 > 0)then
   do i = 0,w0-1
       A0 = A0 + (i+1)*(i+2)/2
   end do
   A0 = A0*4
 end if

!..............
 if(iproc==0)then
    ilast = index(outfile,' ') -1
    open(unit=filenumber,file=outfile(1:ilast)//".trwfn",status='unknown')
    if(write_pnwfn)open(unit=filenumber+1,file=outfile(1:ilast)//".pnwfn",status='unknown')
	
    ilast = index(intfilename,' ')-1
    write(filenumber,*)np(1)
    write(filenumber,*)np(2)
    write(filenumber,*)intfilename(1:ilast), ' ! INTERACTION FILE ' 
    if(hw == 0.0)then
       write(filenumber,*)41./float(A0+np(1)+np(2))**(0.3),' ! HW (approx) ',A0+np(1)+np(2)
    else
       write(filenumber,*)hw,'   !  HW  '
    end if
    write(filenumber,*)Nmshell,' ! # of majors shells '

    if(w0 > 0)then  ! include size of core
        write(filenumber,'(2i5  ," ! total p+n s.p.s, # shells core " ) ' ) & 
          nhsps(itmax)+nhsps(-itmax)+nhsps(itmax)+nhsps(-itmax),w0
    else
         write(filenumber,*)nhsps(itmax)+nhsps(-itmax)+nhsps(itmax)+nhsps(-itmax),' ! total number of p,n s.p. states '
    end if
    write(filenumber,*)MaxWtot-minWtot,' ! Nmax (excitations) '
    write(filenumber,*)dimbasis,' ! # of many-body configurations '
	if(write_pnwfn)write(filenumber+1,*)dimbasis,nxsd(1),nxsd(2)
    write(filenumber,*)iparity,' ! parity '
    write(filenumber,*)jz, ' ! 2 x Jz '
    write(filenumber,*)nkeep,' ! # of eigenstates '
	if(write_pnwfn)write(filenumber+1,*)nkeep
 end if
!------------------- PRINT OUT EIGENENERGIES AND J, T
if(from_old)then  ! added 7.10.7; for reading in previous file
    if(storelanczosincore1.or.storelanczosincoreMPI)then		
        call wfn_rewind(oldwfnfile)		
        call read_wfn_header(oldwfnfile,.false.)		
   	    call wfn_read_nkeep(oldwfnfile, i)  ! dummy to throw away
        ! pull all vectors into memory
		call br_set_histpos(0)
		call br_reset_histbuf()
        do i = 1,nkeep
   	    ! frag1 says that the slicing of the vector follows vec1
          ! new interface, we say which vec to read and it checks
           call wfn_readeigenvec(oldwfnfile, frag1, fcomm1_index, vec1,i,e,xj,xtt)
!           xt =(-0.5 + sqrt(xtt+0.25))
           xt = 0.0
		   
           if(iproc==0)write(filenumber,*)e,xj,xt
           if(iproc==0)print*,e,xj,xt
           call br_grab_vec1()! push vec1 -> br_reg		   
           call br_add2hist(i)		   
        enddo		
    else
        if(iproc ==0)then
           do i = 1,nkeep
!              xt =(-0.5 + sqrt(xtlist(i)+0.25))  
              write(filenumber,*)energy(i),xjlist(i),xt
           end do ! i
        end if
    end if
else
	
 if(.not. storelanczosincore1.and. .not.storelanczosincoreMPI)then
     call wfn_rewind(wfnfile)
     call read_wfn_header(wfnfile,.false.)
	 call wfn_read_nkeep(wfnfile, i)  ! dummy to throw away

     ! pull all vectors into memory
     do i = 1,nkeep
	    ! frag1 says that the slicing of the vector follows vec1
       ! new interface, we say which vec to read and it checks
        call wfn_readeigenvec(wfnfile, frag1, fcomm1_index, vec(:,i),i,e,xj,xtt)
        xt =(-0.5 + sqrt(xtt+0.25))
        if(iproc==0)write(filenumber,*)e,xj,xt
     enddo
 else
     if(iproc ==0)then
        do i = 1,nkeep
           xt =(-0.5 + sqrt(xtlist(i)+0.25))  
           write(filenumber,*)energy(i),xjlist(i),xt
        end do ! i
     end if
 end if
 
end if
 if(iproc==0)then
!------------------ WRITE OUT SINGLE-PARTICLE STATES --------------
 nstate = 0

!------------------- PROTON STATES -------------------------
 ith = -itmax
 do i = 1,nhsps(ith)  ! left states
   nstate = nstate + 1
   write(filenumber,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m, 1
   Nmshell = bmax(Nmshell, hspsqn(ith,i)%w)
 enddo  ! i
 
 ith = itmax
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(filenumber,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m, 1
   Nmshell = bmax(Nmshell, hspsqn(ith,i)%w)

 enddo  ! i
!------------------- NEUTRON STATES -------------------------
 ith = -itmax
 do i = 1,nhsps(ith)  ! left states
   nstate = nstate + 1
   write(filenumber,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m, -1
   Nmshell = bmax(Nmshell, hspsqn(ith,i)%w)

 enddo  ! i
 
 ith = itmax
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(filenumber,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m, -1
   Nmshell = bmax(Nmshell, hspsqn(ith,i)%w)

 enddo  ! i

 end if   
  
!--------------- MANY-BODY CONFIGURATIONS AND THE WAVEFUNCTIONS
!......... CF routine writeoutbasis in file bbasislib6.f90 for a model
!

 do ps = 1,nsectors(1)  ! loop over proton sectors
    do cs = 1,xsd(1)%sector(ps)%ncsectors  ! loop over conjugate neutron sectors
       ns = xsd(1)%sector(ps)%csector(cs)  ! index of conjugate neutron sector

!--------------- construct proton part
       do pblock = 1,xsd(1)%sector(ps)%nhblocks
          prblock = xsd(1)%sector(ps)%rhblock(pblock)
          plblock = xsd(1)%sector(ps)%lhblock(pblock)
          psdstart= xsd(1)%sector(ps)%blockstart(pblock)-1


!---------------- always make left address the outermost loop (a convention)
          npradd = hblock(1)%list(prblock)%nhsd
          npladd = hblock(-1)%list(plblock)%nhsd
          do pladd = 1,npladd
             do pradd = 1,npradd
                ip = psdstart + pradd+(pladd-1)*npradd
!--------------------- now construct proton occupation -- left haiku
                hsd => hblock(-1)%list(plblock)%hsd(:,pladd)
                call convert_haiku_array(-1,hsd,pocctmp)
                n = 0
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i
                  endif
                end do  ! i
!--------------------- now construct proton occupation -- right haiku

                hsd => hblock(1)%list(prblock)%hsd(:,pradd)
                call convert_haiku_array(1,hsd,pocctmp)
                
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i+nhsps(-itmax)
                  endif
                end do  ! i

!----------------construct neutron part

                do nblock = 1,xsd(2)%sector(ns)%nhblocks
                   nrblock = xsd(2)%sector(ns)%rhblock(nblock)
                   nlblock = xsd(2)%sector(ns)%lhblock(nblock)
                   nsdstart= xsd(2)%sector(ns)%blockstart(nblock)-1
!---------------- always make left address the outermost loop (a convention)
                   nnradd = hblock(2)%list(nrblock)%nhsd
                   nnladd = hblock(-2)%list(nlblock)%nhsd
                   do nladd = 1,nnladd
                      do nradd = 1,nnradd
                         in = nsdstart + nradd+(nladd-1)*nnradd
                         ibasis = nstart(in) + pstart(ip)
						 if(ibasis < 0)then
							 if(iproc == 0)then
								 print*,' Error in basis ',in,nstart(in),ip,pstart(ip)
								 
							 end if
							 stop 909
							 
						 end if

!--------------------- now construct neutron occupation -- left haiku
                         hsd => hblock(-2)%list(nlblock)%hsd(:,nladd)
                         call convert_haiku_array(-2,hsd,nocctmp)
                         n = 0
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                                n = n+1
                                nocc(n) = i+nhsps(-itmax)+nhsps(itmax)
                            endif
                         end do  ! i
!--------------------- now construct neutron occupation -- right haiku

                         hsd => hblock(2)%list(nrblock)%hsd(:,nradd)
                         call convert_haiku_array(2,hsd,nocctmp)
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                               n = n+1
                               nocc(n) = i+nhsps(-itmax)+nhsps(itmax)+nhsps(-itmax)
                            endif
                         end do  ! i
						 
!---------------------------- WRITE OUT OCCUPIED STATES -----------------------

                         if(iproc==0) write(filenumber,2222)(pocc(i),i=1,np(1)),(nocc(i),i=1,np(2))
						 
2222                     format(10i5)
!------------------------------- LOOP OVER WFNS -------------------------
!                          call BMPI_BARRIER(icomm,ierr)
 !                         call br_fetch_coef(ibasis,nkeep,vec(ibasis,:)) ! NOTE THIS MAY FAIL IN MPI
						  
                          if(storelanczosincore1 .or. storelanczosincoreMPI)then
#ifdef _MPI
	                          call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
	                          call br_fetch_coef(ibasis,nkeep,vamp) ! NOTE THIS MAY FAIL IN MPI
							  
                             if(iproc==0)then
								 write(filenumber,3333)(vamp(i),i=1,nkeep)
	                             if(write_pnwfn)write(filenumber+1,*)ip,in
								 if(write_pnwfn)write(filenumber+1,3333)(vamp(i),i=1,nkeep)	 
							 end if
                          end if
						  
						   
                           if(.not.storelanczosincore1 .and..not.storelanczosincoreMPI)then
                             if(iproc==0)write(filenumber,3333)(vec(ibasis,i),i=1,nkeep)
                           end if
 
						   
3333                     format(5f14.10)
3334                     format(2i8,5f14.10)

                      enddo ! nradd
                   enddo ! nladd
        
               enddo  ! nblock
             enddo ! pradd
          enddo ! pladd
       enddo ! pblocks
    enddo  ! cs
 enddo ! ps

 if(iproc==0)close(filenumber)

 if(iproc==0)print*,' finished writing .trwfn file '
 return

 end subroutine output_TRDENS
!------------------------------------------------------------------------------------------
! ADDED 7.10.7  -- option 'tx' reads in from .wfn and writes out as .trwfn file
! OBSOLETE-- REMOVE
!

subroutine output_TRDENS_from_oldwfn_boss1
    use localvectors
	use io
	use nodeinfo
    use bvectorlib_mod
	use lanczos_info
	use wfn_mod
	implicit none
	integer :: i
	real :: xe, xj,xt2
	real, allocatable :: energies(:),js(:),t2s(:)
	
    call overlaptribution
    call setup_localvectors

! --- SET UP HISTORY ----

    call wfn_read_nkeep(oldwfnfile, nkeep) ! does BCAST

    if(iproc==0) print*,' There are ', nkeep,' wavefunctions '
	allocate(energies(nkeep),js(nkeep),t2s(nkeep))
    do i = 1,nkeep
       ! new interface - we say which vec to read, it checks
       ! KSM:  This will be very slow, only need to read head part of each vector
       call wfn_readeigenvec(oldwfnfile,frag1, fcomm1_index, vec1,i,xe,xj,xt2)  
       if(iproc==0)print*,i,xe,xj,xt2   
!..... SAVE......
       energies(i)=xe
	   js(i)=xj
	   t2s(i)=xt2

!...... PUSH TO HISTORY......
	   
    enddo
	
!	call output_TRDENS
	return
end subroutine output_TRDENS_from_oldwfn_boss1
	

!------------------------------------------------------------------------------------------
!
!  subroutine WRITE_OUT_BASIS
!
!  added 4/2011 by CWJ @ SDSU
!
!  more user friendly output of basis; based upon subroutine output_TRDENS
!  NOT for post-processing; see routine basis_writeout4post
!
 subroutine write_out_basis

 use system_parameters
 use verbosity
 use spstate
 use W_info
 use haiku_info
 use haikus
 use io
 use basis
 use precisions
 use sectors
 use blocks
 use coupledmatrixelements
 use butil_mod
 use basis
 
 implicit none

 integer i,j,nstate,n
 integer ith
 integer, parameter :: filenumber = 56 
 real e,xj,xt,xtt

 integer ps,cs,ns
 integer pblock, nblock
 integer prblock,plblock
 integer nrblock,nlblock
 integer(kind=basis_prec) :: psdstart, nsdstart
 integer(kind=basis_prec) :: npladd,npradd,pladd,pradd
 integer(kind=basis_prec) :: nnladd,nnradd,nladd,nradd
 integer(kind=basis_prec) :: in,ip
 integer, allocatable :: pocc(:), nocc(:)
 integer(kind=basis_prec) :: ibasis
 integer,allocatable :: pocctmp(:),nocctmp(:)
 integer,pointer  :: hsd(:)

 integer :: Nmshell,wspsmin,wspsmax,it,w
 integer :: ilast
 character(1) :: paritychar
 integer :: aerr

 if(.not.print_basis)return

!................. OPEN OUTPUT FILE..................
 if(writeout)then

    ilast = index(outfile,' ') -1
    open(unit=basis_file,file=outfile(1:ilast)//".basis",status='unknown')
    open(unit=xsd_file,file=outfile(1:ilast)//".xsd",status='unknown')

 else
    open(unit=basis_file,file="basis.bigstick",status='unknown')
    open(unit=  xsd_file,file="xsd.bigstick"  ,status='unknown')

 endif
!.................. ALLOCATE OCCUPATION ............

 allocate(pocc(np(1)), nocc(np(2)), stat=aerr)
 if(aerr /= 0) call memerror("write_out_basis 1")
 allocate(pocctmp(nhsps(1)),nocctmp(nhsps(2)), stat=aerr)
 if(aerr /= 0) call memerror("write_out_basis 2")

!..................BASIC INFORMATION ....................

!..............COMPUTE # OF MAJOR SHELLS
!.............. NOTE THIS ASSUMES CONTIGUOUS........
 Nmshell = 0
 wspsmin = 10000
 wspsmax = 0
 do it = 1,2
 do i = 1,nhsps(it)
    w = hspsqn(it,i)%w
!    print*,w,wceiling(it),hspsqn(i,it)%nr, hspsqn(i,it)%l
    if( wceiling(it) < w)cycle
    wspsmin = bmin(wspsmin,w)
    wspsmax = bmax(wspsmax,w)
 enddo
 enddo
 nmshell = wspsmax-wspsmin +1
!...........................



 write(basis_file,*)np(1), ' ! valence Z '
 write(basis_file,*)np(2), ' ! valence N '
 write(xsd_file,*)np(1), ' ! valence Z '
 write(xsd_file,*)np(2), ' ! valence N '
 write(basis_file,*)Nmshell,' ! # of majors shells '
 write(basis_file,*)MaxWtot-minWtot,' ! Nmax (excitations) '
 write(basis_file,*)dimbasis,' ! # of many-body configurations '
 write(basis_file,*)iparity,' ! parity '
 write(basis_file,*)jz, ' ! 2 x Jz '


!------------------ WRITE OUT SINGLE-PARTICLE STATES --------------
 nstate = 0

!------------------- PROTON STATES -------------------------
  write(basis_file,*)' proton single-particle states '   
  write(basis_file,*)'       label          n          l          2j           2m     '
  write(  xsd_file,*)' proton single-particle states '   
  write(  xsd_file,*)'       label          n          l          2j           2m     '
 ith = -1
 do i = 1,nhsps(ith)  ! left states
   nstate = nstate + 1
   write(basis_file,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m
   write(  xsd_file,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m
   Nmshell = bmax(Nmshell, hspsqn(ith,i)%w)
 enddo  ! i
 
 ith = 1
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(basis_file,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m
   write(  xsd_file,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m
   Nmshell = bmax(Nmshell, hspsqn(ith,i)%w)

 enddo  ! i
!------------------- NEUTRON STATES -------------------------
  write(basis_file,*)' neutron single-particle states '   
  write(basis_file,*)'       label          n          l          2j           2m     '
  write(  xsd_file,*)' neutron single-particle states '   
  write(  xsd_file,*)'       label          n          l          2j           2m     '
 ith = -2
 do i = 1,nhsps(ith)  ! left states
   nstate = nstate + 1
   write(basis_file,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m
   write(  xsd_file,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m
   Nmshell = bmax(Nmshell, hspsqn(ith,i)%w)

 enddo  ! i
 
 ith = 2
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(basis_file,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m
   write(  xsd_file,*)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m
   Nmshell = bmax(Nmshell, hspsqn(ith,i)%w)

 enddo  ! i

!--------------- MANY-BODY CONFIGURATIONS AND THE WAVEFUNCTIONS
!......... CF routine writeoutbasis in file bbasislib6.f90 for a model
!


 do ps = 1,nsectors(1)  ! loop over proton sectors
    if( xsd(1)%sector(ps)%parX == 1)then
        paritychar = '+'
    else
        paritychar = '-'
    endif
    write(basis_file,*)ps,' ! proton sector : 2Jz = ', & 
xsd(1)%sector(ps)%jzX,', parity ',paritychar,', W = ', xsd(1)%sector(ps)%wX
   
    do cs = 1,xsd(1)%sector(ps)%ncsectors  ! loop over conjugate neutron sectors
       ns = xsd(1)%sector(ps)%csector(cs)  ! index of conjugate neutron sector

!--------------- construct proton part
       do pblock = 1,xsd(1)%sector(ps)%nhblocks
          prblock = xsd(1)%sector(ps)%rhblock(pblock)
          plblock = xsd(1)%sector(ps)%lhblock(pblock)
          psdstart= xsd(1)%sector(ps)%blockstart(pblock)-1


!---------------- always make left address the outermost loop (a convention)
          npradd = hblock(1)%list(prblock)%nhsd
          npladd = hblock(-1)%list(plblock)%nhsd
          do pladd = 1,npladd
             do pradd = 1,npradd
                ip = psdstart + pradd+(pladd-1)*npradd
!--------------------- now construct proton occupation -- left haiku
                hsd => hblock(-1)%list(plblock)%hsd(:,pladd)
                call convert_haiku_array(-1,hsd,pocctmp)
                n = 0
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i
                  endif
                end do  ! i
!--------------------- now construct proton occupation -- right haiku

                hsd => hblock(1)%list(prblock)%hsd(:,pradd)
                call convert_haiku_array(1,hsd,pocctmp)
                
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i+nhsps(-1)
                  endif
                end do  ! i

!----------------construct neutron part

                do nblock = 1,xsd(2)%sector(ns)%nhblocks
                   nrblock = xsd(2)%sector(ns)%rhblock(nblock)
                   nlblock = xsd(2)%sector(ns)%lhblock(nblock)
                   nsdstart= xsd(2)%sector(ns)%blockstart(nblock)-1
!---------------- always make left address the outermost loop (a convention)
                   nnradd = hblock(2)%list(nrblock)%nhsd
                   nnladd = hblock(-2)%list(nlblock)%nhsd
                   do nladd = 1,nnladd
                      do nradd = 1,nnradd
                         in = nsdstart + nradd+(nladd-1)*nnradd
                         ibasis = nstart(in) + pstart(ip)
!--------------------- now construct neutron occupation -- left haiku
                         hsd => hblock(-2)%list(nlblock)%hsd(:,nladd)
                         call convert_haiku_array(-2,hsd,nocctmp)
                         n = 0
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                                n = n+1
                                nocc(n) = i+nhsps(-1)+nhsps(1)
                            endif
                         end do  ! i

!--------------------- now construct neutron occupation -- right haiku

                         hsd => hblock(2)%list(nrblock)%hsd(:,nradd)
                         call convert_haiku_array(2,hsd,nocctmp)   ! BUG FIXED 8/2011 CWJ
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                               n = n+1
                               nocc(n) = i+nhsps(-1)+nhsps(1)+nhsps(-2)
                            endif
                         end do  ! i

!---------------------------- WRITE OUT OCCUPIED STATES -----------------------
                          write(basis_file,2222)ibasis,ip,in,(pocc(i),i=1,np(1)),(nocc(i),i=1,np(2))
2222  format('State:',i8,', pSD:',i6,', nSD:',i6,' Occupied:',10i4)


                      enddo ! nradd
                   enddo ! nladd
        
               enddo  ! nblock
             enddo ! pradd
          enddo ! pladd
       enddo ! pblocks
    enddo  ! cs
 enddo ! ps

 close(basis_file)

 write(xsd_file,*)' '
 write(xsd_file,*)' Proton slater determinants '

!.................................. NOW WRITE OUT PROTON, NEUTRON SLATER DETERMINANTS.....
 do ps = 1,nsectors(1)  ! loop over proton sectors
    write(xsd_file,*)ps,' ! sector '

!--------------- construct proton part
       do pblock = 1,xsd(1)%sector(ps)%nhblocks
          prblock = xsd(1)%sector(ps)%rhblock(pblock)
          plblock = xsd(1)%sector(ps)%lhblock(pblock)
          psdstart= xsd(1)%sector(ps)%blockstart(pblock)-1


!---------------- always make left address the outermost loop (a convention)
          npradd = hblock(1)%list(prblock)%nhsd
          npladd = hblock(-1)%list(plblock)%nhsd
          do pladd = 1,npladd
             do pradd = 1,npradd
                ip = psdstart + pradd+(pladd-1)*npradd
!--------------------- now construct proton occupation -- left haiku
                hsd => hblock(-1)%list(plblock)%hsd(:,pladd)
                call convert_haiku_array(-1,hsd,pocctmp)
                n = 0
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i
                  endif
                end do  ! i
!--------------------- now construct proton occupation -- right haiku

                hsd => hblock(1)%list(prblock)%hsd(:,pradd)
                call convert_haiku_array(1,hsd,pocctmp)
                
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i+nhsps(-1)
                  endif
                end do  ! i


!---------------------------- WRITE OUT OCCUPIED STATES -----------------------
                          write(xsd_file,2223)ip, pstart(ip),(pocc(i),i=1,np(1))
2223  format(' SD:',i8,', index:',i8,', Occupied:',10i4)


        
             enddo ! pradd
          enddo ! pladd
       enddo ! pblocks
 enddo ! ps
 write(xsd_file,*)' '
 write(xsd_file,*)' Neutron slater determinants '
 do ns = 1,nsectors(2)  ! loop over neutron sectors
    write(xsd_file,*)ns,' ! sector '


!----------------construct neutron part

                do nblock = 1,xsd(2)%sector(ns)%nhblocks
                   nrblock = xsd(2)%sector(ns)%rhblock(nblock)
                   nlblock = xsd(2)%sector(ns)%lhblock(nblock)
                   nsdstart= xsd(2)%sector(ns)%blockstart(nblock)-1
!---------------- always make left address the outermost loop (a convention)
                   nnradd = hblock(2)%list(nrblock)%nhsd
                   nnladd = hblock(-2)%list(nlblock)%nhsd
                   do nladd = 1,nnladd
                      do nradd = 1,nnradd
                         in = nsdstart + nradd+(nladd-1)*nnradd
                         ibasis = nstart(in) + pstart(ip)

!--------------------- now construct neutron occupation -- left haiku
                         hsd => hblock(-2)%list(nlblock)%hsd(:,nladd)
                         call convert_haiku_array(-2,hsd,nocctmp)
                         n = 0
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                                n = n+1
                                nocc(n) = i+nhsps(-1)+nhsps(1)
                            endif
                         end do  ! i
!--------------------- now construct neutron occupation -- right haiku

                         hsd => hblock(2)%list(nrblock)%hsd(:,nradd)
                         call convert_haiku_array(2,hsd,nocctmp)
                
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                               n = n+1
                               nocc(n) = i+nhsps(-1)+nhsps(1)+nhsps(-2)
                            endif
                         end do  ! i
!---------------------------- WRITE OUT OCCUPIED STATES -----------------------
                          write(xsd_file,2223)in, nstart(in),(nocc(i),i=1,np(2))



                      enddo ! nradd
                   enddo ! nladd
        
               enddo  ! nblock

    enddo  ! ns


 return

end subroutine write_out_basis
!=========================================================================================
!
!  subroutine BASIS_OUT4POSTPROCESSING
!
!  added 7.7.2  by CWJ @ SDSU
!
!  write out details of basis for postprocessing to a .bas file option 'b'
!
!  modified 7.9.1 by CWJ write out in ASCII -- option 'ba'
!
 subroutine basis_out4postprocessing

 use system_parameters
 use verbosity
 use spstate
 use W_info
 use haiku_info
 use haikus
 use io
 use basis
 use precisions
 use sectors
 use blocks
 use coupledmatrixelements
 use butil_mod
 use basis
 use nodeinfo
 use wfn_mod
 
 implicit none

 integer i,j,nstate,n
 integer ith
 real e,xj,xt,xtt

 integer ps,cs,ns
 integer pblock, nblock
 integer prblock,plblock
 integer nrblock,nlblock
 integer(kind=basis_prec) :: psdstart, nsdstart
 integer(kind=basis_prec) :: npladd,npradd,pladd,pradd
 integer(kind=basis_prec) :: nnladd,nnradd,nladd,nradd
 integer(kind=basis_prec) :: in,ip
 integer, allocatable :: pocc(:), nocc(:)
 integer(kind=basis_prec) :: ibasis
 integer,allocatable :: pocctmp(:),nocctmp(:)
 integer,pointer  :: hsd(:)

 integer :: Nmshell,wspsmin,wspsmax,it,w
 integer :: ilast
 character(1) :: paritychar
 integer :: aerr
  
 
 if(iproc/=0)return
 write(6,*)' This will write information on the basis useful for postprocessing '
 write(6,*)' (Written as a binary file, much like a .wfn file)'
 
 write(6,*)' Enter name of file (without extension )'
 read(5,'(a)')outfile
 
 nfragments=1
 print*,' assuming only a single fragment '
!................. OPEN OUTPUT FILE..................

    call wfn_wopen_file(wfnfile,.true.)
	if(baseASCII)then
	    call write_wfn_header_ASCII(wfnfile)
		
	else
        call write_wfn_header(wfnfile)
	end if
!.................. ALLOCATE OCCUPATION ............

 allocate(pocc(np(1)), nocc(np(2)), stat=aerr)
 if(aerr /= 0) call memerror("write_out_basis 1")
 allocate(pocctmp(nhsps(1)),nocctmp(nhsps(2)), stat=aerr)
 if(aerr /= 0) call memerror("write_out_basis 2")
 
 ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 if(.not.baseASCII)then   ! write out in binary  -- ADDED in 7.9.1 -- option 'b'

 nstate = 0

!------------------- PROTON STATES -------------------------
!  write(basis_file,*)' proton single-particle states '   
!  write(basis_file,*)'       label          n          l          2j           2m     '
!  write(  xsd_file,*)' proton single-particle states '   
!  write(  xsd_file,*)'       label          n          l          2j           2m     '
 ith = -1


 write(wfnfile)nhsps    ! MODIFIED IN 7.8.4---write all out at once
 do i = 1,nhsps(ith)  ! left states
   nstate = nstate + 1
   write(wfnfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
   hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par

 enddo  ! i
 
 ith = 1
! write(wfnfile)nhsps(ith)
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(wfnfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
 enddo  ! i
!------------------- NEUTRON STATES -------------------------
!  write(basis_file,*)' neutron single-particle states '   
!  write(basis_file,*)'       label          n          l          2j           2m     '
!  write(  xsd_file,*)' neutron single-particle states '   
!  write(  xsd_file,*)'       label          n          l          2j           2m     '
 ith = -2
! write(wfnfile)nhsps(ith)
 do i = 1,nhsps(ith)  ! left states
   nstate = nstate + 1
   write(wfnfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
 enddo  ! i
 
 ith = 2
! write(wfnfile)nhsps(ith)
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(wfnfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
 enddo  ! i

!.................................. NOW WRITE OUT PROTON, NEUTRON SLATER DETERMINANTS.....
write(wfnfile)nsectors(1),nxSD(1)
 do ps = 1,nsectors(1)  ! loop over proton sectors
    write(wfnfile)ps,xsd(1)%sector(ps)%jzX,xsd(1)%sector(ps)%parX,xsd(1)%sector(ps)%Wx
	write(wfnfile)xsd(1)%sector(ps)%xsdstart,xsd(1)%sector(ps)%xsdend,xsd(1)%sector(ps)%nxsd
	write(wfnfile)xsd(1)%sector(ps)%ncsectors
	write(wfnfile)(xsd(1)%sector(ps)%csector(i),i=1,xsd(1)%sector(ps)%ncsectors)
	write(wfnfile)xsd(1)%sector(ps)%basisstart,xsd(1)%sector(ps)%basisend

!--------------- construct proton part
       do pblock = 1,xsd(1)%sector(ps)%nhblocks
          prblock = xsd(1)%sector(ps)%rhblock(pblock)
          plblock = xsd(1)%sector(ps)%lhblock(pblock)
          psdstart= xsd(1)%sector(ps)%blockstart(pblock)-1


!---------------- always make left address the outermost loop (a convention)
          npradd = hblock(1)%list(prblock)%nhsd
          npladd = hblock(-1)%list(plblock)%nhsd
          do pladd = 1,npladd
             do pradd = 1,npradd
                ip = psdstart + pradd+(pladd-1)*npradd
!--------------------- now construct proton occupation -- left haiku
                hsd => hblock(-1)%list(plblock)%hsd(:,pladd)
                call convert_haiku_array(-1,hsd,pocctmp)
                n = 0
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i
                  endif
                end do  ! i
!--------------------- now construct proton occupation -- right haiku

                hsd => hblock(1)%list(prblock)%hsd(:,pradd)
                call convert_haiku_array(1,hsd,pocctmp)
                
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i+nhsps(-1)
                  endif
                end do  ! i


!---------------------------- WRITE OUT OCCUPIED STATES -----------------------
                          write(wfnfile)ip, pstart(ip),(pocc(i),i=1,np(1))
!2223  format(' SD:',i8,', index:',i8,', Occupied:',10i4)      
             enddo ! pradd
          enddo ! pladd
       enddo ! pblocks
 enddo ! ps
 print*,' finished writing proton SDs '
 write(wfnfile)nsectors(2),nXsd(2)
 do ns = 1,nsectors(2)  ! loop over neutron sectors
    write(wfnfile)ns,xsd(2)%sector(ns)%jzX,xsd(2)%sector(ns)%parX,xsd(2)%sector(ns)%Wx
	write(wfnfile)xsd(2)%sector(ns)%xsdstart,xsd(2)%sector(ns)%xsdend,xsd(2)%sector(ns)%nxsd
	write(wfnfile)xsd(2)%sector(ns)%ncsectors
	write(wfnfile)(xsd(2)%sector(ns)%csector(i),i=1,xsd(2)%sector(ns)%ncsectors)
	write(wfnfile)xsd(2)%sector(ns)%basisstart,xsd(2)%sector(ns)%basisend

!----------------construct neutron part

                do nblock = 1,xsd(2)%sector(ns)%nhblocks
                   nrblock = xsd(2)%sector(ns)%rhblock(nblock)
                   nlblock = xsd(2)%sector(ns)%lhblock(nblock)
                   nsdstart= xsd(2)%sector(ns)%blockstart(nblock)-1
!---------------- always make left address the outermost loop (a convention)
                   nnradd = hblock(2)%list(nrblock)%nhsd
                   nnladd = hblock(-2)%list(nlblock)%nhsd
                   do nladd = 1,nnladd
                      do nradd = 1,nnradd
                         in = nsdstart + nradd+(nladd-1)*nnradd
                         ibasis = nstart(in) + pstart(ip)

!--------------------- now construct neutron occupation -- left haiku
                         hsd => hblock(-2)%list(nlblock)%hsd(:,nladd)
                         call convert_haiku_array(-2,hsd,nocctmp)
                         n = 0
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                                n = n+1
                                nocc(n) = i+nhsps(-1)+nhsps(1)
                            endif
                         end do  ! i
!--------------------- now construct neutron occupation -- right haiku

                         hsd => hblock(2)%list(nrblock)%hsd(:,nradd)
                         call convert_haiku_array(2,hsd,nocctmp)
                
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                               n = n+1
                               nocc(n) = i+nhsps(-1)+nhsps(1)+nhsps(-2)
                            endif
                         end do  ! i
!---------------------------- WRITE OUT OCCUPIED STATES -----------------------
                          write(wfnfile)in, nstart(in),(nocc(i),i=1,np(2))
!						  print*,in,nstart(in),nocc



                      enddo ! nradd
                   enddo ! nladd
        
               enddo  ! nblock

    enddo  ! ns
    print*,' finished writing neutron SDs '
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

else       ! Write out in ASCII  -- option 'ba'
	
 nstate = 0

!------------------- PROTON STATES -------------------------
!  write(basis_file,*)' proton single-particle states '   
!  write(basis_file,*)'       label          n          l          2j           2m     '
!  write(  xsd_file,*)' proton single-particle states '   
!  write(  xsd_file,*)'       label          n          l          2j           2m     '
 ith = -1


 write(wfnfile,*)nhsps    ! MODIFIED IN 7.8.4---write all out at once
 do i = 1,nhsps(ith)  ! left states
   nstate = nstate + 1
   write(wfnfile,'(i6,6i5)')nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
   hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par

 enddo  ! i
 
 ith = 1
! write(wfnfile)nhsps(ith)
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(wfnfile,'(i6,6i5)')nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
 enddo  ! i
!------------------- NEUTRON STATES -------------------------
!  write(basis_file,*)' neutron single-particle states '   
!  write(basis_file,*)'       label          n          l          2j           2m     '
!  write(  xsd_file,*)' neutron single-particle states '   
!  write(  xsd_file,*)'       label          n          l          2j           2m     '
 ith = -2
! write(wfnfile)nhsps(ith)
 do i = 1,nhsps(ith)  ! left states
   nstate = nstate + 1
   write(wfnfile,'(i6,6i5)')nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
 enddo  ! i
 
 ith = 2
! write(wfnfile)nhsps(ith)
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(wfnfile,'(i6,6i5)')nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m,hspsqn(ith,i)%w,hspsqn(ith,i)%par
 enddo  ! i



!.................................. NOW WRITE OUT PROTON, NEUTRON SLATER DETERMINANTS.....
write(wfnfile,*)nsectors(1),nxSD(1)
 do ps = 1,nsectors(1)  ! loop over proton sectors
    write(wfnfile,*)ps,xsd(1)%sector(ps)%jzX,xsd(1)%sector(ps)%parX,xsd(1)%sector(ps)%Wx
	write(wfnfile,*)xsd(1)%sector(ps)%xsdstart,xsd(1)%sector(ps)%xsdend,xsd(1)%sector(ps)%nxsd
	write(wfnfile,*)xsd(1)%sector(ps)%ncsectors
	write(wfnfile,*)(xsd(1)%sector(ps)%csector(i),i=1,xsd(1)%sector(ps)%ncsectors)
	write(wfnfile,*)xsd(1)%sector(ps)%basisstart,xsd(1)%sector(ps)%basisend

!--------------- construct proton part
       do pblock = 1,xsd(1)%sector(ps)%nhblocks
          prblock = xsd(1)%sector(ps)%rhblock(pblock)
          plblock = xsd(1)%sector(ps)%lhblock(pblock)
          psdstart= xsd(1)%sector(ps)%blockstart(pblock)-1


!---------------- always make left address the outermost loop (a convention)
          npradd = hblock(1)%list(prblock)%nhsd
          npladd = hblock(-1)%list(plblock)%nhsd
          do pladd = 1,npladd
             do pradd = 1,npradd
                ip = psdstart + pradd+(pladd-1)*npradd
!--------------------- now construct proton occupation -- left haiku
                hsd => hblock(-1)%list(plblock)%hsd(:,pladd)
                call convert_haiku_array(-1,hsd,pocctmp)
                n = 0
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i
                  endif
                end do  ! i
!--------------------- now construct proton occupation -- right haiku

                hsd => hblock(1)%list(prblock)%hsd(:,pradd)
                call convert_haiku_array(1,hsd,pocctmp)
                
                do i = 1,nhsps(1)
                  if(pocctmp(i) == 1)then
                     n = n+1
                     pocc(n) = i+nhsps(-1)
                  endif
                end do  ! i


!---------------------------- WRITE OUT OCCUPIED STATES -----------------------
                          write(wfnfile,*)ip, pstart(ip),(pocc(i),i=1,np(1))
!2223  format(' SD:',i8,', index:',i8,', Occupied:',10i4)


        
             enddo ! pradd
          enddo ! pladd
       enddo ! pblocks
 enddo ! ps
 print*,' finished writing proton SDs '
 write(wfnfile,*)nsectors(2),nXsd(2)
 do ns = 1,nsectors(2)  ! loop over neutron sectors
    write(wfnfile,*)ns,xsd(2)%sector(ns)%jzX,xsd(2)%sector(ns)%parX,xsd(2)%sector(ns)%Wx
	write(wfnfile,*)xsd(2)%sector(ns)%xsdstart,xsd(2)%sector(ns)%xsdend,xsd(2)%sector(ns)%nxsd
	write(wfnfile,*)xsd(2)%sector(ns)%ncsectors
	write(wfnfile,*)(xsd(2)%sector(ns)%csector(i),i=1,xsd(2)%sector(ns)%ncsectors)
	write(wfnfile,*)xsd(2)%sector(ns)%basisstart,xsd(2)%sector(ns)%basisend

!----------------construct neutron part

                do nblock = 1,xsd(2)%sector(ns)%nhblocks
                   nrblock = xsd(2)%sector(ns)%rhblock(nblock)
                   nlblock = xsd(2)%sector(ns)%lhblock(nblock)
                   nsdstart= xsd(2)%sector(ns)%blockstart(nblock)-1
!---------------- always make left address the outermost loop (a convention)
                   nnradd = hblock(2)%list(nrblock)%nhsd
                   nnladd = hblock(-2)%list(nlblock)%nhsd
                   do nladd = 1,nnladd
                      do nradd = 1,nnradd
                         in = nsdstart + nradd+(nladd-1)*nnradd
                         ibasis = nstart(in) + pstart(ip)

!--------------------- now construct neutron occupation -- left haiku
                         hsd => hblock(-2)%list(nlblock)%hsd(:,nladd)
                         call convert_haiku_array(-2,hsd,nocctmp)
                         n = 0
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                                n = n+1
                                nocc(n) = i+nhsps(-1)+nhsps(1)
                            endif
                         end do  ! i
!--------------------- now construct neutron occupation -- right haiku

                         hsd => hblock(2)%list(nrblock)%hsd(:,nradd)
                         call convert_haiku_array(2,hsd,nocctmp)
                
                         do i = 1,nhsps(2)
                            if(nocctmp(i) == 1)then
                               n = n+1
                               nocc(n) = i+nhsps(-1)+nhsps(1)+nhsps(-2)
                            endif
                         end do  ! i
!---------------------------- WRITE OUT OCCUPIED STATES -----------------------
                          write(wfnfile,*)in, nstart(in),(nocc(i),i=1,np(2))
!						  print*,in,nstart(in),nocc



                      enddo ! nradd
                   enddo ! nladd
        
               enddo  ! nblock

    enddo  ! ns
    print*,' finished writing neutron SDs '	
	
end if
	
	
	close(wfnfile)

 return

end subroutine basis_out4postprocessing
!=========================================================================================
!
!
!  subroutine write_out_jumps
!
!  subroutine to write out selected jumps
!  initiated 4/2011 by CWJ @ SDSU
!
!

subroutine write_out_jumps(it,norder,fchar)

use jumpdef
use jumpNbody
use jump3body
use verbosity
use io
use nodeinfo

implicit none
integer norder   ! order of jumps
integer it       ! species
character(1) :: fchar

type(jumpsect), pointer :: xjump
integer sj      ! sector jump
integer(8) :: xjmp
integer(kind=basis_prec),pointer :: isd(:) => NULL(), fsd(:) => NULL()
integer,pointer ::  op(:) => NULL()
integer ilast

if(.not.print_jumps)return

print*,' PRINTING JUMPS '
if( nproc > 1)then
	if(iproc==0)print*,' Cannot write jumps to a file when in MPI modes'
	return
end if

if(fchar=='o')then   ! open a file
    if(writeout)then
        ilast = index(outfile,' ') -1
        open(unit=jump_file,file=outfile(1:ilast)//".jumps",status='unknown')
    else
        open(unit=jump_file,file="jumps.bigstick",status='unknown')
    endif
endif 

select case (norder)

    case (1)
        xjump => x1bjump(it)
        if(xjump%nsectjumps == 0)return
        if(it == 1)then 
		write(jump_file,*)' proton 1-body jumps '
                isd => p1b_isd0
                fsd => p1b_fsd0
        endif
        if(it == 2)then 
                 write(jump_file,*)' neutron 1-body jumps '
                isd => n1b_isd0
                fsd => n1b_fsd0
        endif
    case (2)
        xjump => x2bjump(it)
        if(xjump%nsectjumps == 0)return
        if(it == 1)then
             write(jump_file,*)' proton 2-body jumps '
                isd => p2b_isd0
                fsd => p2b_fsd0
        endif
        if(it == 2)then
              write(jump_file,*)' neutron 2-body jumps '
                isd => n2b_isd0
                fsd => n2b_fsd0
        endif
    case default
       print*,' Not yet selected ',norder
       return

end select

write(jump_file,*)' There are ',xjump%nsectjumps, ' sector jumps '
write(jump_file,*)xjump%sjmp(:)%nstart
write(jump_file,*)xjump%sjmp(:)%njumps

do sj = 1,xjump%nsectjumps
  write(jump_file,*)sj,': jumps from sector ',xjump%isector(sj),' to ',xjump%fsector(sj)
  write(jump_file,*)'  jump initial SD    final SD   '
  do xjmp = xjump%sjmp(sj)%nstart+1,xjump%sjmp(sj)%nstart+xjump%sjmp(sj)%njumps
       write(jump_file,'(i6,a,3i9)')xjmp,':', isd(xjmp),fsd(xjmp),n2b_op(xjmp)
  end do   ! xjmp

end do   ! sj

write(jump_file,*)' ' 
write(jump_file,*)' # # # '
write(jump_file,*)' '

if(fchar == 'c')close(jump_file)

return

end subroutine write_out_jumps

!----- ADDED 7.10.7 to facilitate option 'tw' -------

subroutine open_trwfn
	use io
	use nodeinfo
	implicit none
	character(80) :: filename
	integer :: ilast
	
	if(iproc/=0)return
	
	print*,' '
1   continue	
	print*,' Enter name of old .trwfn file (do not include extension)'
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	open(unit=oldtrwfnfile,file=filename(1:ilast)//".trwfn",status='old',err=331)
	print*,' Old file opened successfully '
	print*,' '
	print*,' NOTE: You will need to enter in the information for the basis;'
	print*,' it will not be read automatically.'
	print*,' '
	print*,' Output name MUST be different from name of .trwfn file '
	print*,' '
	return

331 continue
    print*,' that file does not exist '
	go to 1	
	
end subroutine open_trwfn
!===============================================================
! FIXED in 7.11.3 for nonsequential ordering of basis in .trwfn
!
subroutine readin_trwfn
    use system_parameters
    use spstate
    use sporbit
    use W_info
    use haiku_info
    use io
    use basis
    use precisions
    use sectors
    use blocks
    use coupledmatrixelements
    use nodeinfo
    use obs
    use localvectors
    use lanczos_info
    use wfn_mod
    use bmpi_mod
    use butil_mod
    use basis
    use mod_reorthog
    use haikus
    use bvectorlib_mod
	implicit none
	
	integer :: nx,ny
	integer(kind=8) :: ndx
	real    :: xx,yy,zz
	character(80) :: title
	logical :: mismatch
	integer :: ierr
	real, allocatable :: vamp(:),elist(:),jlist(:),tlist(:)
	integer(kind=8) :: ibasis,icount,mycount
	integer :: i
	integer, allocatable :: partlist(:)
	integer :: itmax
	integer(kind=basis_prec) :: map_size,map_start
	integer(kind=basis_prec),allocatable :: map(:)
	
	map_size = 500000000 ! I limit the size of the mapping so it does not get too big.
	
	allocate(map(map_size))
	
	mismatch = .false.

    if(np(2) > np(1))then
      itmax = 2
    else
      itmax = 1
    end if
	
	if(iproc==0)then
	   read(oldtrwfnfile,*)nx
	   if(nx /= np(1))then
		   print*,' Mismatch in Z ',np(1),nx
		   mismatch = .true.
		   goto 321
	   end if
	   read(oldtrwfnfile,*)nx
	   if(nx /= np(2))then
		   print*,' Mismatch in N ',np(2),nx
		   mismatch = .true.
		   goto 321
	   end if
	   read(oldtrwfnfile,'(a)')title
	   print*,title
	   read(oldtrwfnfile,*)xx
	   print*,' hw = ',xx
	   read(oldtrwfnfile,*)nx
	   print*,' # of major shells = ',nx
	   read(oldtrwfnfile,*)nx
	   if(nx /= nhsps(itmax)+nhsps(-itmax)+nhsps(itmax)+nhsps(-itmax))then
		   print*,' Mismatch in # p+n s.p. states ',nhsps(itmax)+nhsps(-itmax)+nhsps(itmax)+nhsps(-itmax),nx
		   mismatch = .true.
		   goto 321
	   end if   
	   read(oldtrwfnfile,*)nx
	   print*,' Nmax = ',nx
	   read(oldtrwfnfile,*)ndx
	   if(ndx /= dimbasis)then
		   print*,' Mismatch in basis dim ',dimbasis,ndx
		   mismatch = .true.
		   goto 321
	   end if
	   print*,' Basis dimension agreees: ',ndx
	   read(oldtrwfnfile,*)nx
	   if(nx /= iparity)then
		   print*,' Mismatch in parity ',iparity,nx
		   mismatch = .true.
		   goto 321
	   end if	   
	   read(oldtrwfnfile,*)nx
	   if(nx /= jz)then
		   print*,' Mismatch in 2x Jz ',jz,nx
		   mismatch = .true.
		   goto 321
	   end if	 	   
	   read(oldtrwfnfile,*)nkeep
	   print*,nkeep,' eigenvectors '

   end if
#ifdef _MPI
   call BMPI_Bcast(nkeep,1,0,MPI_COMM_WORLD,ierr)
#endif
   
   niter = nkeep
   call wfn_write_nkeep(nkeep)
   allocate(vamp(nkeep),elist(nkeep),jlist(nkeep),tlist(nkeep))
   allocate(partlist(np(1)+np(2)))
	   
   call overlaptribution
   call setup_localvectors
   	   
   if(iproc==0)then
!---- read in energies, Js -------------
      do i = 1,nkeep
		  read(oldtrwfnfile,*)xx,yy,zz
		  print*,xx,yy,zz
		  elist(i)=xx
		  jlist(i)=yy
		  tlist(i)=zz
	  end do ! i

!----- read in and compare .sps q#s ------
      do i = 1,nhsps(itmax)+nhsps(-itmax)+nhsps(itmax)+nhsps(-itmax)
		  read(oldtrwfnfile,*)nx
	  end do
   end if
!----- read in coefficients -------------	 
  map_start = 0  
  call trwfn_basis_map(map_size,map_start,map)
  do icount = 1,dimbasis
	  if(icount > map_start+map_size)then
		  map_start = map_start+map_size
		  call trwfn_basis_map(map_size,map_start,map)
	  end if
	  ibasis = map(icount-map_start)
!	  print*,' IBASIS ',icount,ibasis
	  if(ibasis == 0)stop
	  if(iproc==0)then
		  read(oldtrwfnfile,*)(partlist(i),i=1,np(1)+np(2))
		  read(oldtrwfnfile,*)(vamp(i),i=1,nkeep)
	  end if	
	  call br_push_coef(ibasis,nkeep,vamp)
  end do
!----- FINALLY WRITE OUT TO WFN ---------------
do i = 1,nkeep
	call br_set_histpos(i)

    call br_retrieve_hist(i)
    call br_restore_vec1()
	call wfn_writeeigenvec(wfnfile, frag1, vec1, i,elist(i),jlist(i),tlist(i)*(tlist(i)+1))
	
end do
	
	return
	
321 continue
    stop 2023
	
	
end subroutine readin_trwfn
!-------------------------------------------------------------
!  added 7.11.3 by CWJ @SDSU to fix a bug
!  the order of basis states in .trwfn is NOT the same as in .wfn
!  this maps from .trwfn -- .wfn
!  to save both time and space, the map is found in chunks

subroutine trwfn_basis_map(map_size,map_start,map)
    use system_parameters
    use spstate
    use sporbit
    use W_info
    use haiku_info
    use io
    use basis
    use precisions
    use sectors
    use blocks
!    use coupledmatrixelements
    use nodeinfo
!    use obs
!    use localvectors
!    use lanczos_info
!    use wfn_mod
    use bmpi_mod
    use butil_mod
!    use mod_reorthog
    use haikus
	implicit none
	integer(kind=basis_prec) :: map_size,map_start
	integer(kind=basis_prec) :: map(map_size)
    integer ps,cs,ns
    integer pblock, nblock
    integer prblock,plblock
    integer(8) :: psdstart
    integer nrblock,nlblock
    integer(8) :: nsdstart
    integer npladd,npradd,pladd,pradd
    integer nnladd,nnradd,nladd,nradd
    integer(8) :: ip,in
    integer(kind=basis_prec) :: ibasis,icount


	!--------------- MANY-BODY CONFIGURATIONS AND THE WAVEFUNCTIONS
	!......... CF routine writeoutbasis in file bbasislib6.f90 for a model
	!

	 icount = 0
	 do ps = 1,nsectors(1)  ! loop over proton sectors
	    do cs = 1,xsd(1)%sector(ps)%ncsectors  ! loop over conjugate neutron sectors
	       ns = xsd(1)%sector(ps)%csector(cs)  ! index of conjugate neutron sector

	!--------------- construct proton part
	       do pblock = 1,xsd(1)%sector(ps)%nhblocks
	          prblock = xsd(1)%sector(ps)%rhblock(pblock)
	          plblock = xsd(1)%sector(ps)%lhblock(pblock)
	          psdstart= xsd(1)%sector(ps)%blockstart(pblock)-1


	!---------------- always make left address the outermost loop (a convention)
	          npradd = hblock(1)%list(prblock)%nhsd
	          npladd = hblock(-1)%list(plblock)%nhsd
	          do pladd = 1,npladd
	             do pradd = 1,npradd
	                ip = psdstart + pradd+(pladd-1)*npradd
	


	!----------------construct neutron part

	                do nblock = 1,xsd(2)%sector(ns)%nhblocks
	                   nrblock = xsd(2)%sector(ns)%rhblock(nblock)
	                   nlblock = xsd(2)%sector(ns)%lhblock(nblock)
	                   nsdstart= xsd(2)%sector(ns)%blockstart(nblock)-1
	!---------------- always make left address the outermost loop (a convention)
	                   nnradd = hblock(2)%list(nrblock)%nhsd
	                   nnladd = hblock(-2)%list(nlblock)%nhsd
	                   do nladd = 1,nnladd
	                      do nradd = 1,nnradd
	                         in = nsdstart + nradd+(nladd-1)*nnradd
	                         ibasis = nstart(in) + pstart(ip)
							 icount = icount + 1
							 if(ibasis < 0)then
								 if(iproc == 0)then
									 print*,' Error in basis ',in,nstart(in),ip,pstart(ip)
								 
								 end if
								 stop 909
							 
							 end if
							 
							 if(icount > map_start+map_size)return
							 
							 if(icount > map_start)map(icount-map_start)=ibasis
							 
!							 print*,' counting ',icount, map_start,icount-map_start,ibasis

	                      enddo ! nradd
	                   enddo ! nladd
        
	               enddo  ! nblock
	             enddo ! pradd
	          enddo ! pladd
	       enddo ! pblocks
	    enddo  ! cs
	 enddo ! ps	
	 if(icount < dimbasis)then
		 print*,' Some problem, missing basis states ',icount,dimbasis
		 stop 9189
	 endif
	
	return
end subroutine trwfn_basis_map

!-------------------------------------------------------------
!  added 7.10.8 by CWJ @LLNL by request of OCG
!
!  routine to find the parity of a state when both parities are allowed
!  must allow for situations where an MPI process has zero amplitudes on it
!  a bit kludgey:
!     for each fragment, search for the first amplitude which is clearly nonzero
!     combine on the root; find the first
!     then find parity of that state
!
!
subroutine find_parity(myparity,whichvec)
	use sporbit
	use basis
	use precisions
	use localvectors
	use sectors
    use blocks
	use nodeinfo
    use haiku_info
	use haikus
	use spstate
	use bmpi_mod
	use butil_mod
	implicit none
	integer :: myparity
	integer :: whichvec
	integer(kind=basis_prec) :: ibasis,kbasis
	real(kind=8) :: thresh
    integer ps,cs,ns
    integer pblock, nblock
    integer prblock,plblock
    integer(8) :: psdstart
    integer nrblock,nlblock
    integer(8) :: nsdstart
    integer npladd,npradd,pladd,pradd
    integer nnladd,nnradd,nladd,nradd
    integer(8) :: ip,in
    integer, allocatable :: pocc(:), nocc(:)
    integer,allocatable :: pocctmp(:),nocctmp(:)
    integer,pointer  :: hsd(:)
	integer :: i,itmax,n,aerr
	integer :: tmppar
	integer, allocatable :: parlist(:)
	integer :: iip
	integer :: tag, ierr,kproc
#ifdef _MPI
    type(MPI_Status) :: stat	
!	integer :: stat(MPI_STATUS_SIZE)
#endif
	
	if(.not.allsameparity)then
		ibasis =1
		go to 123
	end if
		
    do i = 1,nhsps(-1)
		hspsqn(-1,i)%par = (-1)**hspsqn(-1,i)%l
	end do
    do i = 1,nhsps(1)
		hspsqn(1,i)%par = (-1)**hspsqn(1,i)%l
	end do
    do i = 1,nhsps(-2)
		hspsqn(-2,i)%par = (-1)**hspsqn(-2,i)%l
	end do
    do i = 1,nhsps(2)
		hspsqn(2,i)%par = (-1)**hspsqn(2,i)%l
	end do
!........ FIND AN AMPLITUDE WHICH IS NONZERO.....	
	thresh=1.d0/sqrt(float(dimbasis))
	ibasis = -1
	
	if(nproc==1 .or. isfragroot)then
		if(whichvec == 1)then
	    do kbasis = v1s,v1e
		   if(abs(vec1(kbasis)) > thresh)then
		     	ibasis = kbasis
			    exit
		   end if
	    end do
	    else
		    do kbasis = v2s,v2e
			   if(abs(vec2(kbasis)) > thresh)then
			     	ibasis = kbasis
				    exit
			   end if
		    end do			
		end if
   end if
!..... RETURN INFORMATION TO ROOT PROCESS......
    allocate(parlist(0:nproc-1))
	parlist = -1
	
	if(nproc > 1)then
		parlist(0)=ibasis
		if(iproc > 0  )then
			tag = iproc
#ifdef _MPI
		    call BMPI_SEND(ibasis,1,0,tag,MPI_COMM_WORLD,ierr)
#endif
		else
			do kproc = 1,nproc-1
				tag = kproc
#ifdef _MPI
				call BMPI_Recv(kbasis,1,kproc,tag,MPI_COMM_WORLD,stat,ierr)
#endif
				parlist(kproc)=kbasis
			end do
		end if
#ifdef _MPI		
		call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
		ibasis = parlist(0)
		do kproc = 1,nproc-1
			if(parlist(kproc)>0)ibasis = min(ibasis,parlist(kproc))
			
		end do
		
	end if
    
!......FIND proton, neutron decomposition - search
123 continue ! from beginning of routine

if(ibasis == 0 .or. iproc > 0)return
allocate(pocctmp(nhsps(1)),nocctmp(nhsps(2)), stat=aerr)

 do ps = 1,nsectors(1)  ! loop over proton sectors
    do cs = 1,xsd(1)%sector(ps)%ncsectors  ! loop over conjugate neutron sectors
       ns = xsd(1)%sector(ps)%csector(cs)  ! index of conjugate neutron sector

!--------------- construct proton part
       do pblock = 1,xsd(1)%sector(ps)%nhblocks
          prblock = xsd(1)%sector(ps)%rhblock(pblock)
          plblock = xsd(1)%sector(ps)%lhblock(pblock)
          psdstart= xsd(1)%sector(ps)%blockstart(pblock)-1


!---------------- always make left address the outermost loop (a convention)
          npradd = hblock(1)%list(prblock)%nhsd
          npladd = hblock(-1)%list(plblock)%nhsd
          do pladd = 1,npladd
             do pradd = 1,npradd
                ip = psdstart + pradd+(pladd-1)*npradd
				
                do nblock = 1,xsd(2)%sector(ns)%nhblocks
                   nrblock = xsd(2)%sector(ns)%rhblock(nblock)
                   nlblock = xsd(2)%sector(ns)%lhblock(nblock)
                   nsdstart= xsd(2)%sector(ns)%blockstart(nblock)-1
!---------------- always make left address the outermost loop (a convention)
                   nnradd = hblock(2)%list(nrblock)%nhsd
                   nnladd = hblock(-2)%list(nlblock)%nhsd
                   do nladd = 1,nnladd
                      do nradd = 1,nnradd
                         in = nsdstart + nradd+(nladd-1)*nnradd
                         kbasis = nstart(in) + pstart(ip)
						 if(ibasis == kbasis)then
							 myparity = 1
!
!--------------------- now construct proton occupation -- left haiku
                             hsd => hblock(-1)%list(plblock)%hsd(:,pladd)
                             call convert_haiku_array(-1,hsd,pocctmp)
                             do i = 1,nhsps(-1)
                                if(pocctmp(i) == 1)then
								   myparity = myparity*hspsqn(-1,i)%par
                                endif
                             end do  ! i
!--------------------- now construct proton occupation -- right haiku

                             hsd => hblock(1)%list(prblock)%hsd(:,pradd)
                             call convert_haiku_array(1,hsd,pocctmp)
                             do i = 1,nhsps(1)
                                if(pocctmp(i) == 1)then
								   myparity = myparity*hspsqn(1,i)%par
                                endif
                             end do  ! i						 
!
!--------------------- now construct neutron occupation -- left haiku
                             hsd => hblock(-2)%list(nlblock)%hsd(:,nladd)
                             call convert_haiku_array(-2,hsd,nocctmp)
                             do i = 1,nhsps(-2)
                                if(nocctmp(i) == 1)then
								   myparity = myparity*hspsqn(-2,i)%par
                                endif
                             end do  ! i

!--------------------- now construct neutron occupation -- right haiku

                             hsd => hblock(2)%list(nrblock)%hsd(:,nradd)
                             call convert_haiku_array(2,hsd,nocctmp)
                             do i = 1,nhsps(2)
                                if(nocctmp(i) == 1)then
								   myparity = myparity*hspsqn(2,i)%par
                                endif
                             end do  ! i	
							 					 
							 goto 999
						 end if
					  end do ! nradd
				   end do    ! nladd
			    end do  ! nblock
				
			end do  ! pradd
		  end do    ! pladd
	   end do ! pblock

   end do ! cs
end do ! ps
999 continue

	return
	
end subroutine find_parity
!-------------------------------------------------------------
