!
!
!  subroutine output_TRDENS
!
!  routine to print out wfns in format useful for P. Navratil's TRDENS code
!  called after routine LANCZOS has finished
!
!  MODIFICATION 3/2013: write out proton and neutron states with the same number of 
!  s.p. states
!.........................................................................

 subroutine output_TRDENS

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
! use interaction
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
 implicit none
! include 'binterfaces.inc'

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
  integer stat(MPI_STATUS_SIZE)
  integer :: fromproc
  integer :: aerr

!................. OPEN OUTPUT FILE..................


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
    open(unit=filenumber+1,file=outfile(1:ilast)//".pnwfn",status='unknown')
	
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
	write(filenumber+1,*)dimbasis,nxsd(1),nxsd(2)
    write(filenumber,*)iparity,' ! parity '
    write(filenumber,*)jz, ' ! 2 x Jz '
    write(filenumber,*)nkeep,' ! # of eigenstates '
	write(filenumber+1,*)nkeep
 end if
!------------------- PRINT OUT EIGENENERGIES AND J, T
 if(.not. storelanczosincore1.and. .not.storelanczosincoreMPI)then
     call wfn_rewind(wfnfile)
     call read_wfn_header(wfnfile,.false.)
	 call wfn_read_nkeep(wfnfile, i)  ! dummy to throw away

     ! pull all vectors into memory
     do i = 1,nkeep
	    ! frag1 says that the slicing of the vector follows vec1
       ! new interface, we say which vec to read and it checks
        call wfn_readeigenvec(wfnfile, frag1, fcomm1, vec(:,i),i,e,xj,xtt)
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
                         call convert_haiku_array(1,hsd,nocctmp)
                
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
!                         do i = 1,nkeep
!                            call br_retrieve_hist(i)
!                            call br_restore_vec1()
!						end do

!print*,ibasis,' here i B',iproc
!call BMPI_BARRIER(icomm,ierr)
                          call br_fetch_coef(ibasis,nkeep,vamp)
!						  print*,ibasis,' here i C',iproc
						  
                          if(storelanczosincore1)then
                             if(iproc==0)then
								 write(filenumber,3333)(vamp(i),i=1,nkeep)
	                             write(filenumber+1,*)ip,in
								 write(filenumber+1,3333)(vamp(i),i=1,nkeep)
								 
								 
							 end if
                          end if
                          if(storelanczosincoreMPI)then
!---- EXTRACT THE VECTOR ELEMENTS------------------  
!    find out which node the stored elements IBASIS reside on
!                              tag = 777
!                              if(ibasis >= Lstart .and. ibasis < Lstart +Ldim-1)then
!                                  do i = 1,nkeep
!		                              call br_retrieve_hist(i)
!		                              call br_restore_vec1()
 !                                     vamp(i) = Lvec(ibasis-Lstart+1,i)
 !                                 end do
!---- SEND TO ROOT NODE
 !                                  if(iproc /= 0)call BMPI_SEND(vamp,nkeep,0,tag,icomm,ierr) 

 !                              end if
 !                              if(iproc==0 .and. ibasis > Ldim)then
       
 !                                 do fromproc = 0,nproc-1
 !                                    if(ibasis >= Lstart_r(fromproc).and. ibasis <Lstart_r(fromproc+1))exit
 !                                 end do
 !                                 if(fromproc>0)call BMPI_Recv(vamp,nkeep,fromproc,tag,icomm,stat,  ierr) 
 !                              end if
                               if(iproc==0)write(filenumber,3333)(vamp(i),i=1,nkeep)
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


 return

 end subroutine output_TRDENS

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
!  write out details of basis for postprocessing to a .bas file
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
    call write_wfn_header(wfnfile)
!.................. ALLOCATE OCCUPATION ............

 allocate(pocc(np(1)), nocc(np(2)), stat=aerr)
 if(aerr /= 0) call memerror("write_out_basis 1")
 allocate(pocctmp(nhsps(1)),nocctmp(nhsps(2)), stat=aerr)
 if(aerr /= 0) call memerror("write_out_basis 2")


 nstate = 0

!------------------- PROTON STATES -------------------------
!  write(basis_file,*)' proton single-particle states '   
!  write(basis_file,*)'       label          n          l          2j           2m     '
!  write(  xsd_file,*)' proton single-particle states '   
!  write(  xsd_file,*)'       label          n          l          2j           2m     '
 ith = -1
! print*,' testing output ',nhsps
 write(wfnfile)nhsps    ! MODIFIED IN 7.8.4---write all out at once
 do i = 1,nhsps(ith)  ! left states
   nstate = nstate + 1
   write(wfnfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
   hspsqn(ith,i)%m,hspsqn(ith,i)%w

 enddo  ! i
 
 ith = 1
! write(wfnfile)nhsps(ith)
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(wfnfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m,hspsqn(ith,i)%w
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
                             hspsqn(ith,i)%m,hspsqn(ith,i)%w
 enddo  ! i
 
 ith = 2
! write(wfnfile)nhsps(ith)
 do i = 1,nhsps(ith)  ! right states
   nstate = nstate + 1
   write(wfnfile)nstate,hspsqn(ith,i)%nr,hspsqn(ith,i)%l,hspsqn(ith,i)%j, & 
                             hspsqn(ith,i)%m,hspsqn(ith,i)%w
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
       write(jump_file,'(i6,a,2i9)')xjmp,':', isd(xjmp),fsd(xjmp)
  end do   ! xjmp

end do   ! sj

write(jump_file,*)' ' 
write(jump_file,*)' # # # '
write(jump_file,*)' '

if(fchar == 'c')close(jump_file)

return

end subroutine write_out_jumps
