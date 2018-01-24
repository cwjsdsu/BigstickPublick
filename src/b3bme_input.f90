!===============================================================
!
!  file B3BME_INPUT.f90
!
!  routines to generate/read in 3-body matrix elements
!
!  initiated 9/09 by CWJ @ SDSU
!
!  some routines to read in 3-body matrix elements 
!  adapted from routines by E. Jurgenson, LLNL (Fall 2010)
!
!==================================================================
!
! THE FOLLOWING ADAPTED 11/2010 FROM E. JURGENSON, LLNL August 2010
! which in turn was adapted from routines by P. Navratil
!
module harmosc_TMP
   !!! general parameters of the space
   integer :: hbo
   integer :: nhom,nhom1sp,nhom2sp
!
! NOTE:  for principal quantum number N_a, 
!     any N_a <=  nhom1sp
!     any N_a + N_b <= nhom2sp
!     any N_a + N_b + N_c <= nhom
! or else the matrix element is skipped
!
   integer :: prot3eff,neut3eff,nucl3eff
   integer :: mjtot2,mttot2,parity,mtotlim=3
   integer(8) :: tot_num_of_3bmatel
   !!! input file vars
   character(len=240) :: input_file_name,output_file_name,basis_file_name
   logical :: input_file_name_exist
   integer :: input_unit=3,basis_unit=2,output_unit=4,str_pos,istat
end module harmosc_TMP

!!! module to keep track of basis state information
module spbasis_TMP
	integer :: isp1ntot,isp3ntot
	integer,allocatable :: isp1n(:,:),isp3n(:,:)
        integer, allocatable,target :: map2psps(:), map2nsps(:)   ! maps to "standard" bigstick s.p. states
                                                           ! protons and neutrons
        integer, allocatable,target :: cmap2psps(:), cmap2nsps(:)  ! maps to s.p. states with opposite m

        integer, allocatable :: mapXXX(:), phaseXXX(:)
        integer, allocatable :: cmapXXX(:), cphaseXXX(:)   ! maps to 3-body states with all M's reversed
        integer, allocatable :: tphaseXXX(:)  ! phase that arises from time-reversal of all M's
                                              ! (-1)**(j1+j2+j3 + 1)
end module spbasis_TMP
!------------------------------------------------------------------
subroutine fetch3bodyinput

   use flags3body
   use nodeinfo
   use harmosc_TMP
   use io
   use interactions3body
   use coupledmatrixelements
   use bmpi_mod
   implicit none

   integer :: ierr
   integer :: ifile

   if(.not.threebodycheck .or. .not.threebody) return
    
   if(iproc == 0)then
    write(6,*)' Enter filename for 3-body force '
    write(6,*)' (if none then enter "none"; still will convert 2-body to 3-body )'
    write(6,*)' (If filename is long, you can store in file FILENAME3BODY.TXT '
    write(6,*)'  and it will be auto read from there; enter "auto") '
    if(auto_input)then
        read(autoinputfile,'(a)')input_file_name
        print*,input_file_name
    else
       read(5,'(a)')input_file_name
        print*,input_file_name

   endif
   end if
  
     call BMPI_BARRIER(icomm,ierr)

     call BMPI_BCAST(input_file_name,12,0,icomm,ierr)
!    print*,iproc,' input ',input_file_name(1:5)
    select case (input_file_name(1:4))

      case ('none','NONE')   ! no 3-body
       if(iproc==0)then 
         print*,' no 3-body input '
         if(.not.auto_input)write(autoinputfile,'(a)')input_file_name
       endif
       return

      case ('auto','AUTO')
       threebody = .true.
       open(unit=4,file='FILENAME3BODY.TXT',status='old',err=111)
       read(4,'(a)')input_file_name

       goto 112

111 continue
       if(iproc==0)print*,' You have to create FILENAME3BODY.TXT '
       stop

      case default
       threebody = .true.
   
    end select
112 continue
    if(.not.auto_input .and. iproc==0)write(autoinputfile,'(a)')input_file_name

    if(iproc == 0)then
       if(auto_input)then
           read(autoinputfile,*)scale3body
       else
          write(6,*)' Enter scaling strength for 3-body '
          read(5,*)scale3body
          write(autoinputfile,*)scale3body
       endif
    end if 

    call get_file_pars
!
!  here's where one would autocount
!
   if(iproc==0)then
    open(input_unit,file=TRIM(ADJUSTL(input_file_name)),status='old',form='unformatted', action='read')

    end if
    emptyHam = .false.   ! successfully opened a file
   
  return
end subroutine fetch3bodyinput

!------------------------------------------------------------------

subroutine master_read_3bmes
   use spbasis_TMP
   use harmosc_TMP
   use nodeinfo
   implicit none
!------------------------ OPEN FILE!-----------------------------

	! Get the par file from the command line
!        print*,' Save that LOONNNGGG 3-body file name in fort.68 '
!        read(68,'(a)')input_file_name
!	call GET_COMMAND_ARGUMENT(1,input_file_name) ! parameter file name

	!!! parse the input file name for some info
        if(iproc==0)print*,'3 body file ',TRIM(ADJUSTL(input_file_name))

        if(input_file_name(1:4) == 'none' .or. input_file_name(1:4)=='NONE')return

        if(TRIM(ADJUSTL(input_file_name)) == 'none' .or. TRIM(ADJUSTL(input_file_name))=='NONE')return
	call get_file_pars

	!!! open the files
   open(input_unit,file=TRIM(ADJUSTL(input_file_name)),status='old',form='unformatted', action='read')
!----------------------- SET UP MAP BETWEEN BIGSTICK s.p. states and assumed ordering
        if(iproc==0)print*,' Max Nhw for s.p. states ',nhom1sp
 
!        call autocount_3bmesMFD(input_unit)

        call sp1nbas(nhom1sp)
        call mastermake3body

   return
end subroutine master_read_3bmes

!
! master subroutine, written by CWJ @ SDSU Nov 2010
! to control E. Jurgenson routines 
!
subroutine read_LLNL_3bmes
   use spbasis_TMP
   use harmosc_TMP
   implicit none
!------------------------ OPEN FILE!-----------------------------

	! Get the par file from the command line
        print*,' Save that LOONNNGGG 3-body file name in fort.68 '
        read(68,'(a)')input_file_name
!	call GET_COMMAND_ARGUMENT(1,input_file_name) ! parameter file name

	!!! parse the input file name for some info
	call get_file_pars

	!!! open the files
   open(input_unit,file=TRIM(ADJUSTL(input_file_name)),status='old',form='unformatted', action='read')
!----------------------- SET UP MAP BETWEEN BIGSTICK s.p. states and assumed ordering
        print*,' Max Nhw for s.p. states ',nhom1sp
 
        call autocount_3bmesMFD(input_unit)

        call sp1nbas(nhom1sp)
        call mastermake3body

   return
end subroutine read_LLNL_3bmes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! A routine to parse the file name for a couple of variables
!!! will stop if it detects an error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_file_pars()
   use harmosc_TMP
   use nodeinfo
	use bmpi_mod
	implicit none
   integer(4) :: ierr

   integer :: nmax
	logical :: input_file_exist
	real :: temp_input
   if(iproc == 0)then
      write(*,*) 'echoing input filename:'
      write(*,*) input_file_name
   end if
	!!! check the status of the input file   
   inquire(file=TRIM(ADJUSTL(input_file_name)), exist=input_file_exist)
   if (.not.input_file_exist) then
        if(iproc==0)then
           write(*,*) ' file not found: ',input_file_name
           write(*,*) ' Please make sure you have the filename right.'
        endif
		call BMPI_BARRIER(icomm,ierr)
		call BMPI_FINALIZE(ierr)
      stop
   end if
   
   open(input_unit,file=TRIM(ADJUSTL(input_file_name)),status='old',form='unformatted', action='read',iostat=istat)
   if (istat/=0) then 
   	write(*,*) ' *** Error: something wrong with initial input file form/name.'
      write(*,*) '             Input must have original name for parsing.'
   	stop
   end if
   
   read(input_unit,iostat=istat) temp_input
   if(istat/=0) then
   	write(*,*) ' *** Error: something wrong with initial input file form.'
   	write(*,*) '           Input file must be fortran binary matrix elements'
   	stop
   end if   

	close(input_unit)
   
   !!! assign the ouput file names
   str_pos = index(input_file_name,'.int')
   if(str_pos/=0) then 
		output_file_name = input_file_name(1:str_pos)//'form'
      basis_file_name = input_file_name(1:str_pos)//'basis'
   else
   	write(*,*) ' something wrong with input file ==> no .int'
      stop
   end if


	!!! get some parameters from the file name
   str_pos = index(input_file_name,'hw')+2
   read(input_file_name(str_pos:str_pos+1),fmt='(I2)',iostat=istat) hbo 
   if(str_pos/=2 .and. istat==0) then 
   	write(*,*) 'hw: ',input_file_name(str_pos:str_pos+1)
   else
      write(*,*) ' something wrong with input file ==> no hw'
      write(*,*) '   Input must have original name for parsing.'
      if (istat/=0) write(*,*) ' iostat = ',istat
      stop
   end if
   
   str_pos = index(input_file_name,'_nmax')+5
   read(input_file_name(str_pos:str_pos+1),fmt='(I2)',iostat=istat) nmax
   if(str_pos/=5 .and. istat==0) then 
   	write(*,*) 'nmax: ',input_file_name(str_pos:str_pos+1)
   else
   	write(*,*) ' something wrong with input file ==> no nmax'
      write(*,*) '   Input must have original name for parsing.'
      if (istat/=0) write(*,*) ' iostat = ',istat
   	stop
   end if
   
   !str_pos = index(input_file_name,'_CC')+3
   !if(str_pos/=3) then 
 	if(index(input_file_name,'_CC')/=0) then 
 		prot3eff=100
      neut3eff=100
   else if(index(input_file_name,'_li6_')/=0) then 
 		prot3eff=3
      neut3eff=3
   else if(index(input_file_name,'_pshell_')/=0) then 
 		prot3eff=3
      neut3eff=4
   else if(index(input_file_name,'_Z')/=0) then 
 	   str_pos = index(input_file_name,'_Z')+2
      read(input_file_name(str_pos:str_pos),fmt='(I1)',iostat=istat) prot3eff
      if(str_pos/=2 .and. istat==0) then 
   	   write(*,*) 'neutrons: ',input_file_name(str_pos:str_pos)
      else
   	   write(*,*) ' something wrong with input file ==> no Z'
         write(*,*) '   Input must have original name for parsing.'
         if (istat/=0) write(*,*) ' iostat = ',istat
   	   stop
      end if

      str_pos = index(input_file_name,'_Z3N')+4
      read(input_file_name(str_pos:str_pos),fmt='(I1)',iostat=istat) neut3eff
      if(str_pos/=4 .and. istat==0) then 
   	   write(*,*) 'protons: ',input_file_name(str_pos:str_pos)
      else
   	   write(*,*) ' something wrong with input file ==> no N'
         write(*,*) '   Input must have original name for parsing.'
         if (istat/=0) write(*,*) ' iostat = ',istat
   	   stop
      end if
   else
   	write(*,*) ' something wrong with input file ==> no N'
      write(*,*) '   Input must have original name for parsing.'
      if (istat/=0) write(*,*) ' iostat = ',istat
      stop  
   end if   
   nucl3eff=prot3eff+neut3eff
	write(*,*) 'nucl3eff: ',nucl3eff
	
	!!! choose the truncation scheme due to Pauli blocking
   if (nucl3eff<=5) then
   	! all extra nucleons are in the bottom shell
      nhom = nmax+0
   	nhom1sp = nhom
      nhom2sp = nhom
   else if (nucl3eff==200) then
   	!CC needs square truncations (all active nucleons can excite by nmax)
      nhom = nmax*3  !! all three up at nmax above the ground
   	nhom1sp = nmax  !! one at nmax above the ground
      nhom2sp = nmax*2  !! two at nmax above the ground
   else if (nucl3eff==6) then
   	! Li6 is first p-shell nuc - has only one extra nucleon at N=1
      nhom = nmax+2
   	nhom1sp = nhom-1
      nhom2sp = nhom
   else if (nucl3eff>6 .and. nucl3eff<=12) then
   	! rest of p-shell nuclei are satisfied by this for 3-body forces
      nhom = nmax+3
   	nhom1sp = nhom-2
      nhom2sp = nhom-1
   else
   	write(*,*) ' something wrong with input file ==> bad value for nucleons, A = ',nucl3eff
      stop
   end if
	
   ! check the decisions
   write(*,*) 'nhom: ',nhom,'  nhom2sp: ',nhom2sp,'  nhom1sp: ',nhom1sp

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! A routine to parse the file name for a couple of variables
!!! will stop if it detects an error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_file_parsOLD()
	use harmosc_TMP
	implicit none
   
   integer :: nmax
	logical :: input_file_exist
	real :: temp_input

   write(*,*) 'echoing input filename:'
   write(*,*) input_file_name

	!!! check the status of the input file   
   inquire(file=TRIM(ADJUSTL(input_file_name)), exist=input_file_exist)
   if (.not.input_file_exist) then
   	write(*,*) ' file not found: ',input_file_name
      write(*,*) ' Please make sure you have the filename right.'
   	stop
   end if
   
   open(input_unit,file=TRIM(ADJUSTL(input_file_name)),status='old',form='unformatted', action='read',iostat=istat)
   if (istat/=0) then 
   	write(*,*) ' *** Error: something wrong with initial input file form/name.'
      write(*,*) '             Input must have original name for parsing.'
   	stop
   end if
   
   read(input_unit,iostat=istat) temp_input
   if(istat/=0) then
   	write(*,*) ' *** Error: something wrong with initial input file form.'
   	write(*,*) '           Input file must be fortran binary matrix elements'
   	stop
   end if   

	close(input_unit)
   
   !!! assign the ouput file names
   str_pos = index(input_file_name,'.int')
   if(str_pos/=0) then 
		output_file_name = input_file_name(1:str_pos)//'form'
      basis_file_name = input_file_name(1:str_pos)//'basis'
   else
   	write(*,*) ' something wrong with input file ==> no .int'
      stop
   end if


	!!! get some parameters from the file name
   str_pos = index(input_file_name,'hw')+2
   read(input_file_name(str_pos:str_pos+1),fmt='(I2)',iostat=istat) hbo 
   if(str_pos/=2 .and. istat==0) then 
   	write(*,*) 'hw: ',input_file_name(str_pos:str_pos+1)
   else
   	write(*,*) ' something wrong with input file ==> no hw'
      write(*,*) '   Input must have original name for parsing.'
      if (istat/=0) write(*,*) ' iostat = ',istat
   	stop
   end if
   
   str_pos = index(input_file_name,'_nmax')+5
   read(input_file_name(str_pos:str_pos+1),fmt='(I2)',iostat=istat) nmax
   if(str_pos/=5 .and. istat==0) then 
   	write(*,*) 'nmax: ',input_file_name(str_pos:str_pos+1)
   else
   	write(*,*) ' something wrong with input file ==> no nmax'
      write(*,*) '   Input must have original name for parsing.'
      if (istat/=0) write(*,*) ' iostat = ',istat
   	stop
   end if
   
   str_pos = index(input_file_name,'_Z')+2
   read(input_file_name(str_pos:str_pos),fmt='(I1)',iostat=istat) prot3eff
   print*,prot3eff,str_pos,istat
   if(str_pos/=2 .and. istat==0) then 
   	write(*,*) 'neutrons: ',input_file_name(str_pos:str_pos)
   else
   	write(*,*) ' something wrong with input file ==> no Z'
      write(*,*) '   Input must have original name for parsing.'
      if (istat/=0) write(*,*) ' iostat = ',istat
!   	stop
   end if
	
   str_pos = index(input_file_name,'_Z3N')+4
   read(input_file_name(str_pos:str_pos),fmt='(I1)',iostat=istat) neut3eff
   if(str_pos/=4 .and. istat==0) then 
   	write(*,*) 'protons: ',input_file_name(str_pos:str_pos)
   else
   	write(*,*) ' something wrong with input file ==> no N'
      write(*,*) '   Input must have original name for parsing.'
      if (istat/=0) write(*,*) ' iostat = ',istat
!   	stop
   end if
   nucl3eff=prot3eff+neut3eff
	write(*,*) 'nucl3eff: ',nucl3eff
	
	!!! choose the truncation scheme due to Pauli blocking
   if (nucl3eff<=5) then
   	! all extra nucleons are in the bottom shell
      nhom = nmax+0
   	nhom1sp = nhom
      nhom2sp = nhom
   else if (nucl3eff==6) then
   	! Li6 is first p-shell nuc - has only one extra nucleon at N=1
      nhom = nmax+2
   	nhom1sp = nhom-1
      nhom2sp = nhom
   else if (nucl3eff>6 .and. nucl3eff<=12) then
   	! rest of p-shell nuclei are satisfied by this for 3-body forces
      nhom = nmax+3
   	nhom1sp = nhom-2
      nhom2sp = nhom-1
   else
   	write(*,*) ' something wrong with input file ==> bad value for nucleons, A = ',nucl3eff
      stop
   end if
	
   ! check the decisions
   write(*,*) 'nhom: ',nhom,'  nhom2sp: ',nhom2sp,'  nhom1sp: ',nhom1sp

end subroutine get_file_parsOLD
!===================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! routine that builds the A=1 sp basis
! + MAP TO INTRINSIC BIGSTICK STATES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sp1nbas(nhom)
   use spbasis_TMP
   use spstate
   use haiku_info
   implicit none
   integer,intent(IN) :: nhom
   integer :: n,l,j2,mj2,mt2
   integer :: ii,ntot
   integer :: isps
   integer :: aerr
   ii=0
   do ntot=0,nhom
      do l=mod(ntot,2),ntot,2
         n=(ntot-l)/2
         do j2=iabs(2*l-1),2*l+1,2
            do mj2=-j2,j2,2
               ii=ii+1
            end do
         end do
      end do
   end do
   isp1ntot=ii
   if(allocated(isp1n))deallocate(isp1n)
   allocate(isp1n(5,2*isp1ntot), stat=aerr)
   if(aerr /= 0) call memerror("sp1nbas 1")

!   write(iunitout,*) 'Number of sp states:',isp1ntot
!   write(iunitout,*) ' #(t_z=1/2)   n  l  j  j_z  #(t_z=-1/2)'
   ii=0

!....... ALLOCATE MAPS TO STANDARD BIGSTICK S.P. STATES
   if(allocated(map2psps))deallocate(map2psps,cmap2psps) 
   if(allocated(map2nsps))deallocate(map2nsps,cmap2nsps) 

   allocate( map2psps( isp1ntot ),map2nsps( isp1ntot+1:2*isp1ntot ) , stat=aerr)
   if(aerr /= 0) call memerror("sp1nbas 2")
   allocate( cmap2psps( isp1ntot ),cmap2nsps( isp1ntot+1:2*isp1ntot ) , stat=aerr)
   if(aerr /= 0) call memerror("sp1nbas 3")

   map2psps(:) = -1
   map2nsps(:) = -1
   cmap2psps(:)= -1
   cmap2nsps(:)= -1

   do ntot=0,nhom
      do l=mod(ntot,2),ntot,2
         n=(ntot-l)/2
         do j2=iabs(2*l-1),2*l+1,2
            do mj2=-j2,j2,2
               ii=ii+1
               isp1n(1,ii)=n
               isp1n(2,ii)=l
               isp1n(3,ii)=j2
               isp1n(4,ii)=mj2
               isp1n(5,ii)=1
               isp1n(1:4,ii+isp1ntot)=isp1n(1:4,ii)
               isp1n(5,ii+isp1ntot)=-1
!..................... FIND WHICH STANDARD STATE THIS IS...............

               if( mj2 >= 0 )then  
               do isps = 1,nhsps(1)  ! search over protons
                  if( n == hspsqn(1,isps)%nr  .and.  & 
                      l == hspsqn(1,isps)%l   .and.  & 
                      j2 == hspsqn(1,isps)%j   .and.  & 
                      mj2 == hspsqn(1,isps)%m )then
                           if(map2psps(ii) /= -1)then   ! check this hasn't been filled already
                                print*,' oops proton already filled ',ii,isps
                                stop
                           endif
                           map2psps(ii) = isps+nhsps(-1)  ! need this offset
                  endif
               end do
               do isps = 1,nhsps(-1)  ! search over protons

                  if( n == hspsqn(-1,isps)%nr  .and.  & 
                      l == hspsqn(-1,isps)%l   .and.  & 
                      j2 == hspsqn(-1,isps)%j   .and.  & 
                      mj2 == -hspsqn(-1,isps)%m )then
                           if(cmap2psps(ii) /= -1)then   ! check this hasn't been filled already
                                print*,' oops proton already filled ',ii,isps
                                stop
                           endif
                           cmap2psps(ii) = isps  
                  endif

               end do
               do isps = 1,nhsps(2)  ! search over neutrons
                  if( n == hspsqn(2,isps)%nr  .and.  & 
                      l == hspsqn(2,isps)%l   .and.  & 
                      j2 == hspsqn(2,isps)%j   .and.  & 
                      mj2 == hspsqn(2,isps)%m )then
                           if(map2nsps(ii+isp1ntot) /= -1)then   ! check this hasn't been filled already
                                print*,' oops neutron already filled ',ii,isps
                                stop
                           endif
                           map2nsps(ii+isp1ntot) = isps+nhsps(-2)
                  endif
               end do
               do isps = 1,nhsps(-2)  ! search over neutrons
                  if( n == hspsqn(-2,isps)%nr  .and.  & 
                      l == hspsqn(-2,isps)%l   .and.  & 
                      j2 == hspsqn(-2,isps)%j   .and.  & 
                      mj2 == -hspsqn(-2,isps)%m )then
                           if(cmap2nsps(ii+isp1ntot) /= -1)then   ! check this hasn't been filled already
                                print*,' oops neutron already filled ',ii,isps
                                stop
                           endif
                           cmap2nsps(ii+isp1ntot) = isps
                  endif

               end do



               else
               do isps = 1,nhsps(-1)  ! search over protons
                  if( n == hspsqn(-1,isps)%nr  .and.  & 
                      l == hspsqn(-1,isps)%l   .and.  & 
                      j2 == hspsqn(-1,isps)%j   .and.  & 
                      mj2 == hspsqn(-1,isps)%m )then
                           if(map2psps(ii) /= -1)then   ! check this hasn't been filled already
                                print*,' oops proton already filled ',ii,isps
                                stop
                           endif
                           map2psps(ii) = isps  
                  endif
               end do
               do isps = 1,nhsps(1)  ! search over protons

                  if( n == hspsqn(1,isps)%nr  .and.  & 
                      l == hspsqn(1,isps)%l   .and.  & 
                      j2 == hspsqn(1,isps)%j   .and.  & 
                      mj2 == -hspsqn(1,isps)%m )then
                           if(cmap2psps(ii) /= -1)then   ! check this hasn't been filled already
                                print*,' oops proton already filled ',ii,isps
                                stop
                           endif
                           cmap2psps(ii) = isps  + nhsps(-1)
                  endif

               end do

               do isps = 1,nhsps(-2)  ! search over neutrons
                  if( n == hspsqn(-2,isps)%nr  .and.  & 
                      l == hspsqn(-2,isps)%l   .and.  & 
                      j2 == hspsqn(-2,isps)%j   .and.  & 
                      mj2 == hspsqn(-2,isps)%m )then
                           if(map2nsps(ii+isp1ntot) /= -1)then   ! check this hasn't been filled already
                                print*,' oops neutron already filled ',ii,isps
                                stop
                           endif
                           map2nsps(ii+isp1ntot) = isps
                  endif
               end do
               do isps = 1,nhsps(2)  ! search over neutrons
                  if( n == hspsqn(2,isps)%nr  .and.  & 
                      l == hspsqn(2,isps)%l   .and.  & 
                      j2 == hspsqn(2,isps)%j   .and.  & 
                      mj2 == -hspsqn(2,isps)%m )then
                           if(cmap2nsps(ii+isp1ntot) /= -1)then   ! check this hasn't been filled already
                                print*,' oops neutron already filled ',ii,isps
                                stop
                           endif
                           cmap2nsps(ii+isp1ntot) = isps+nhsps(-2)
                  endif

               end do
               endif


!               write(iunitout,'(i9,i4,i3,i3,i4,i5)') ii,n,l,j2,mj2,ii+isp1ntot
            end do
         end do
      end do
   end do
end subroutine sp1nbas

!=======================================================================
!
!  given s.p. state indices a,b,c find the PPP or NNN index and a phase
!
! INPUT:
!  it = 1,2  species for triplet XXX
!  a,b,c = native BIGSTICK label of s.p.s. (but must convert to handed labels)
!
! OUTPUT:
!  indxtrip: index of triplet, to be used for input into fetching matrix element index
!  iphase  : phase arising from reordering into native BIGSTICK order
!  tphase  : phase arising from sending all m's -> -m
!
!
subroutine findXXXlabel(it,a,b,c,indxtrip,iphase,tphase)

   use interactions3body
   use spstate
   use haiku_info
   use jump_mod
   use welder
   implicit none

   integer :: it   ! species
   integer :: a, b, c
   integer :: atmp, btmp, ctmp
   integer :: indxtrip
   integer :: iphase
   integer :: tphase
   integer, pointer :: map(:)
   integer n
   integer ja,jb,jc

   n =nhsps(it) + nhsps(-it)
   n = n*(n-1)*(n-2)/6

  if(it ==1)then
        map => mapPPP
  else
        map => mapNNN
  endif
   iphase = 1
!...............GET TIME-REVERSED PHASE FOR A, B, C
   if( a > nhsps(-it))then
      ja = hspsqn(it,a -nhsps(-it))%j
   else
      ja = hspsqn(-it,a)%j
   end if
   if( b > nhsps(-it))then
      jb = hspsqn(it,b -nhsps(-it))%j
   else
      jb = hspsqn(-it,b)%j
   end if
   if( c > nhsps(-it))then
      jc = hspsqn(it,c -nhsps(-it))%j
   else
      jc = hspsqn(-it,c)%j
   end if

   tphase = (-1)**( (ja+jb+jc+1)/2)
   atmp = a
   btmp = b
   ctmp = c
!........ ORDER THESE SO THAT atmp > btmp > ctmp
   call swap(btmp,ctmp,iphase)
   call swap(atmp,btmp,iphase)
   call swap(btmp,ctmp,iphase)
   indxtrip = (atmp-3)*(atmp-2)*(atmp-1)/6 + (btmp-2)*(btmp-1)/2 + ctmp 
   if(indxtrip > n)then
       print*,' some problem ',n
       print*,atmp,btmp,ctmp,indxtrip    
       stop
   endif
   indxtrip = map(indxtrip)
   if(indxtrip == 0)indxtrip = -1

   return

end subroutine findXXXlabel
!==========================================================================


subroutine mastermake3body

   use system_parameters
   use spbasis_TMP
   use harmosc_TMP
   use interactions3body
   use ntuple_info
   use verbosity

   implicit none
   integer :: ipspi,ipspf,i
   integer c3state,d3state
   integer cc3state,cd3state

   integer it,par3
   integer :: iref3, istart3,i3bme ,i3bmetr
   integer :: ciref3, cistart3
   integer :: cphase, dphase
   integer :: ccphase, cdphase,tphasec,tphased

   integer :: c3statei,c3statef,d3statei,d3statef
   integer :: phasei,phasef,tphasei,tphasef

   integer :: mmax, mmin  ! limits for getting mref etc

   real(4),allocatable :: temp(:)
   real, pointer :: hmatXXX(:),hmatXXY(:)
   real, target :: hmatXXX_dummy(1), hmatXXY_dummy(1)
   logical :: semidiagflag
   integer nmex
   integer :: aerr
   real :: diagfactor    !  NOTE -- Diagonals are actually applied TWICE. (This is a relatively small inefficiency)
                         !  but this means need a factor of 1/2 for diagonal matrix elements

   cistart3 = -1000000 ! trigger bounds error if used
   ciref3 = -2000000
   iref3 = -10000000
   istart3 = -3000000
   hmatXXX => hmatXXX_dummy
   hmatXXY => hmatXXY_dummy

	!!! counter for the total number of A=3 states
   tot_num_of_3bmatel=0
   

   !!! loop over A=3 channels
    if(verbose_3body_readin)then
	write(*,*) 'looping over 3-body channels,  total number: ',4*(((nhom*2+3)/2)+1)*2
    endif
   nmex = 0
   do mttot2=-mtotlim,mtotlim,2

   !!! BREAK UP INTO PPP/NNN and PPN/PNN CASES

      if(abs(mttot2) == 3 )then  ! PPP/NNN

      do mjtot2=1,2*nhom+3,2   ! only loops over M > 0; need to also account for M < 0
         do parity=1,-1,-2
         
				!!! build the list of A=3 states
            call spXYZbas(parity,nhom,nhom2sp,mjtot2,mttot2,basis_unit,'a')
            if (isp3ntot==0) cycle
            call spXYZbas(parity,nhom,nhom2sp,mjtot2,mttot2,basis_unit,'f')
            
				!!! tally the total number of matrix elements
            tot_num_of_3bmatel = tot_num_of_3bmatel + isp3ntot*(isp3ntot+1)/2
!............... RETRIEVE SHIFTS.............      
             if( mttot2 == -3)then
                     it = 2
                     hmatXXX => hmatNNN
             endif
             if( mttot2 == 3 )then
                     it = 1
                     hmatXXX => hmatPPP
             endif
             mmin = XXX3(it)%jzstart
             mmax = XXX3(it)%jzend


             par3 = (3-parity)/2
             if(mjtot2 >= mmin .and. mjtot2 <= mmax)then
               iref3 = XXX3(it)%meref(mjtot2,par3)
               istart3 = XXX3(it)%mestart(mjtot2,par3)
             endif
             if(-mjtot2 >= mmin .and. -mjtot2 <= mmax)then

               ciref3 = XXX3(it)%meref(-mjtot2,par3)
               cistart3 = XXX3(it)%mestart(-mjtot2,par3)
             endif
            !!! read in unformatted version and write out to formatted file

            do ipspi=1,isp3ntot
               !!! read in v3trans row from current record
               allocate(temp(ipspi:isp3ntot), stat=aerr)
               if(aerr /= 0) call memerror("mastermake3body 10")
               read(input_unit,iostat=istat) temp
              
               if( np(it) < 3) then
                   deallocate(temp)
                   cycle
                endif
               !!! write out each element in a new row, along with single particle states
               !!! or take it in to your personal code 
               do ipspf=ipspi,isp3ntot

                  call getXXXindicia(it,isp3n(:,ipspi),.false.,c3state,cphase,tphasec)
                  call getXXXindicia(it,isp3n(:,ipspf),.false.,d3state,dphase,tphased)


                 if( c3state > 0 .and. d3state > 0)then
                    if(c3state <= d3state)then
                       i3bme = istart3 + (d3state-iref3)*(d3state-iref3-1)/2+c3state-iref3
                    else
                      i3bme = istart3 + (c3state-iref3)*(c3state-iref3-1)/2+d3state-iref3
                    endif
                    if(i3bme ==0)then
                        print*,d3state,c3state,mjtot2,parity
                        stop
                    endif
!................... BECAUSE HERMITICITY EITHER WAY, DIAGONALS APPLIED TWICE SO MULTIPLY BY 1/2 .....
                    if( ipspf == ipspi)then
                      diagfactor = 0.5
                    else
                      diagfactor = 1.0
                    end if
!....................... NOW SEND TO MY STORED 3-BODY MATRIX ELEMENTS...................  
                    hmatXXX(i3bme) = hmatXXX(i3bme) + temp(ipspf)*cphase*dphase*diagfactor*scale3body

                 endif

!---------------------- NOW DO TIME-REVERSE 
                 call getXXXindicia(it,isp3n(:,ipspi),.true.,c3state,cphase,tphasec)
                  call getXXXindicia(it,isp3n(:,ipspf),.true.,d3state,dphase,tphased)

                 if( c3state > 0 .and. d3state > 0)then

                    if(c3state <= d3state)then
                       i3bme = cistart3 + (d3state-ciref3)*(d3state-ciref3-1)/2+c3state-ciref3
                    else
                      i3bme = cistart3 + (c3state-ciref3)*(c3state-ciref3-1)/2+d3state-ciref3
                    endif
                    if(i3bme ==0)then
                        print*,d3state,c3state,mjtot2,parity
                        stop
                    endif


!................... BECAUSE HERMITICITY EITHER WAY, DIAGONALS APPLIED TWICE SO MULTIPLY BY 1/2 .....
                    if( ipspf == ipspi)then
                      diagfactor = 0.5
                    else
                      diagfactor = 1.0
                    end if
!....................... NOW SEND TO MY STORED 3-BODY MATRIX ELEMENTS...................  
                    hmatXXX(i3bme) = hmatXXX(i3bme) +temp(ipspf)*cphase*dphase*tphasec*(tphased)*diagfactor*scale3body

                 endif

               end do

               deallocate(temp)

            end do ! loop ipspi
            
         end do  ! loop parity
      end do   ! loop on mjtot2

   else    
!................ DO PPN/PNN................
      if(mttot2 == -1)then
          it = 2
          hmatXXY => hmatPNN
      else
          it = 1
          hmatXXY => hmatPPN

      endif

      do mjtot2=1,2*nhom+3,2   ! only loops over M > 0; need to also account for M < 0
         do parity=1,-1,-2
         
				!!! build the list of A=3 states
            call spXYZbas(parity,nhom,nhom2sp,mjtot2,mttot2,basis_unit,'a')
            if (isp3ntot==0) cycle
            call spXYZbas(parity,nhom,nhom2sp,mjtot2,mttot2,basis_unit,'f')
            
				!!! tally the total number of matrix elements
            tot_num_of_3bmatel = tot_num_of_3bmatel + isp3ntot*(isp3ntot+1)/2

             par3 = (3-parity)/2

            !!! read in unformatted version and write out to formatted file

            do ipspi=1,isp3ntot

               !!! read in v3trans row from current record
               allocate(temp(ipspi:isp3ntot), stat=aerr)
               if(aerr /= 0) call memerror("mastermake3body 30")
               read(input_unit,iostat=istat) temp

      
               if( np(it) < 2 .or. np(3-it) < 1)then
                   deallocate(temp)
                   cycle
               endif

               do ipspf = ipspi, isp3ntot
                 
                  call getXXYindicia(it,isp3n(:,ipspi),.false.,c3statei,d3statei,phasei,tphasei)
                  call getXXYindicia(it,isp3n(:,ipspf),.false.,c3statef,d3statef,phasef,tphasef)
                  call check_semidiag(it,isp3n(:,ipspi),isp3n(:,ipspf),semidiagflag)


!................... BECAUSE HERMITICITY EITHER WAY, DIAGONALS APPLIED TWICE SO MULTIPLY BY 1/2 .....
                    if( ipspf == ipspi .or. semidiagflag)then
!                    if( ipspf == ipspi )then
                      diagfactor = 0.5
                    else
                      diagfactor = 1.0
                    end if
!....................... NOW SEND TO MY STORED 3-BODY MATRIX ELEMENTS................ ...  

                    if(c3statei > -1 .and. d3statef > -1)then
                       i3bme = c3statei+d3statef
!                       if(hmatXXY(i3bme) /= 0.0)then
!                          print*,' what the heck ... ',i3bme,hmatXXY(i3bme)
!                          stop
!                       endif
                       hmatXXY(i3bme) = hmatXXY(i3bme) + temp(ipspf)*phasei*phasef*diagfactor*scale3body
                    else
                       i3bme = -1
                    endif

                    if(c3statef > -1 .and. d3statei > -1)then

                       i3bmetr = c3statef+d3statei
                       if(i3bmetr /= i3bme)then

!                       if(hmatXXY(i3bmetr) /= 0.0)then
!                          print*,' what the heck (tr) ... ',i3bme,hmatXXY(i3bmetr)
!                          stop
!                       endif
                           hmatXXY(i3bmetr) = hmatXXY(i3bmetr) + temp(ipspf)*phasef*phasei*diagfactor*scale3body
                       endif
                    endif
!......................... NOW FOR TIME-REVERSED STATES...SEND M TO -M...
                  call getXXYindicia(it,isp3n(:,ipspi),.true.,c3statei,d3statei,phasei,tphasei)
                  call getXXYindicia(it,isp3n(:,ipspf),.true.,c3statef,d3statef,phasef,tphasef)
                  call check_semidiag(it,isp3n(:,ipspi),isp3n(:,ipspf),semidiagflag)


!................... BECAUSE HERMITICITY EITHER WAY, DIAGONALS APPLIED TWICE SO MULTIPLY BY 1/2 .....
                    if( ipspf == ipspi .or. semidiagflag)then
!                    if( ipspf == ipspi )then
                      diagfactor = 0.5
                    else
                      diagfactor = 1.0
                    end if
!....................... NOW SEND TO MY STORED 3-BODY MATRIX ELEMENTS................ ...  

                    if(c3statei > -1 .and. d3statef > -1)then
                       i3bme = c3statei+d3statef
!                       if(hmatXXY(i3bme) /= 0.0)then
!                          print*,' what the heck ... ',i3bme,hmatXXY(i3bme)
!                          stop
!                       endif
                       hmatXXY(i3bme) = hmatXXY(i3bme) + temp(ipspf)*phasei*phasef*diagfactor*tphasei*tphasef*scale3body
                    else
                       i3bme = -1
                    endif

                    if(c3statef > -1 .and. d3statei > -1)then

                       i3bmetr = c3statef+d3statei
                       if(i3bmetr /= i3bme)then

!                       if(hmatXXY(i3bmetr) /= 0.0)then
!                          print*,' what the heck (tr) ... ',i3bme,hmatXXY(i3bmetr)
!                          stop
!                       endif
                           hmatXXY(i3bmetr) = hmatXXY(i3bmetr) + temp(ipspf)*phasef*phasei*diagfactor*tphasei*tphasef*scale3body
                       endif
                    endif


               end do
 
               deallocate(temp)

            end do ! loop ipspi
            
         end do  ! loop parity
      end do   ! loop on mjtot2

   endif

   end do  ! loop on mttot2
	write(*,*) 'total number of A=3 matrix elements processed: ',tot_num_of_3bmatel
end subroutine mastermake3body

!==========================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! routine that builds the A=3 sp basis
!  NOTE ultimately a bit inefficient for use in BIGSTICK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spXYZbas(parity,nhom,nhom2sp,mjtot2,mttot2,iunitout,fillchar)
   use spbasis_TMP
   use system_parameters
   use verbosity
   implicit none
   integer,intent(IN) :: parity,nhom,nhom2sp,mjtot2,mttot2,iunitout
   integer :: ntot_a,n_a,l_a,j2_a,mj2_a,mt2_a,isp1n_a
   integer :: ntot_b,n_b,l_b,j2_b,mj2_b,mt2_b,isp1n_b
   integer :: ntot_c,n_c,l_c,j2_c,mj2_c,mt2_c,isp1n_c
   integer :: ii
   character*1 fillchar

!.... choices for fillchar:    'a' means count and allocate; 'f' means fill; any other just counts

   integer :: asps, bsps, csps
   integer :: aerr

	!!! count the states
   ii=0
   do isp1n_a=1,2*isp1ntot
      n_a=isp1n(1,isp1n_a)
      l_a=isp1n(2,isp1n_a)
      j2_a=isp1n(3,isp1n_a)
      mj2_a=isp1n(4,isp1n_a)
      mt2_a=isp1n(5,isp1n_a)
      ntot_a=2*n_a+l_a
      if (ntot_a>nhom) cycle
      do isp1n_b=isp1n_a+1,2*isp1ntot
         n_b=isp1n(1,isp1n_b)
         l_b=isp1n(2,isp1n_b)
         j2_b=isp1n(3,isp1n_b)
         mj2_b=isp1n(4,isp1n_b)
         mt2_b=isp1n(5,isp1n_b)
         ntot_b=2*n_b+l_b
         if (ntot_a+ntot_b>nhom2sp) cycle
         do isp1n_c=isp1n_b+1,2*isp1ntot
            n_c=isp1n(1,isp1n_c)
            l_c=isp1n(2,isp1n_c)
            j2_c=isp1n(3,isp1n_c)
            mj2_c=isp1n(4,isp1n_c)
            mt2_c=isp1n(5,isp1n_c)
            ntot_c=2*n_c+l_c
            if (ntot_a+ntot_c>nhom2sp) cycle
            if (ntot_b+ntot_c>nhom2sp) cycle
            if (ntot_a+ntot_b+ntot_c>nhom) cycle
            if (mj2_a+mj2_b+mj2_c/=mjtot2) cycle
            if (mt2_a+mt2_b+mt2_c/=mttot2) cycle
            if ((-1)**(l_a+l_b+l_c)/=parity) cycle
            ii=ii+1 
            if(fillchar=='f')then
            isp3n(1,ii)=isp1n_a
            isp3n(2,ii)=isp1n_b
            isp3n(3,ii)=isp1n_c

            end if
         end do
      end do
   end do
   isp3ntot=ii
!   write(iunitout,1000) isp3ntot,mjtot2,parity,mttot2
   if(verbose_3body_readin)then
      write(*,1000) isp3ntot,mjtot2,parity,mttot2
   endif
   1000 format(' Number of 3N sp states:',i8, ' for mjtot=',i4,'/2   parity=',i2,'   mttot=',i4,'/2')

	!!! allocate and fill them in
   if(fillchar=='a')then
   if (allocated(isp3n)) deallocate(isp3n)
   allocate(isp3n(3,isp3ntot), stat=aerr)
   if(aerr /= 0) call memerror("spXYZbas")
   endif 


   return
        
end subroutine spXYZbas
!=====================================================================

subroutine getXXXindicia(it,ispXXXn,conflag,indx,iphase,tphase)

  use spbasis_TMP

  implicit none
  integer :: it  ! species of XXX
  integer :: ispXXXn(3)  ! 3 labels
  logical conflag   ! get conjugate states?
  integer :: indx
  integer iphase   ! from reordering into native bigstick order
  integer tphase   ! from sending all m -> -m
  integer asps, bsps, csps 
  integer, pointer :: map(:)


  indx = -1
  iphase = -999
  tphase = -999

  if(it == 1)then
       if(ispXXXn(1) > isp1ntot .or. ispXXXn(2) > isp1ntot .or. ispXXXn(3) > isp1ntot )then
           print*,' problem proton ', ispXXXn(:)
           stop
       endif
  else
       if(ispXXXn(1) <= isp1ntot .or. ispXXXn(2) <= isp1ntot .or.  ispXXXn(3) <= isp1ntot )then
           print*,' problem neutron ', ispXXXn(:)
           stop
       endif
  endif

  if(it ==1)then
    if(conflag)then
        map => cmap2psps
    else
        map   => map2psps
    endif
  else
    if(conflag)then
       map => cmap2nsps
    else
       map   => map2nsps
    endif
  endif

  asps = map(ispXXXn(1))
  bsps = map(ispXXXn(2))
  csps = map(ispXXXn(3))


  if(asps < 1 .or. bsps < 1 .or. csps < 1)then

      return
  endif

  call findXXXlabel(it,asps,bsps,csps,indx,iphase,tphase)

  return
end subroutine getXXXindicia

!=====================================================================
!
!
!
subroutine getXXYindicia(itx,ispXYZ,conflag,cindx,dindx,iphase,tphase)

  use spbasis_TMP

  implicit none
  integer :: itx  ! species of XX
  integer :: ity  ! species of Y
  integer :: ispXYZ(3)  ! 3 labels
  logical conflag   ! get conjugate states?
  integer :: cindx,dindx
  integer iphase   ! from reordering into native bigstick order
  integer tphase   ! from sending all m -> -m
  integer asps, bsps, csps 
  integer, pointer :: mapX(:),mapY(:)

  ity = 3-itx

  cindx = -1
  dindx = -1
  iphase = -999
  tphase = -999

!....... WE KNOW ispXYZ(1) < ispXYZ(2) < ispXYZ(3) ......USE AS ERROR TRAP......

  if(itx == 1)then
       if(ispXYZ(1) > isp1ntot .or. ispXYZ(2) > isp1ntot .or. ispXYZ(3) <= isp1ntot )then
           print*,' problem PPN ', ispXYZ(:)
           stop
       endif
  else
       if(ispXYZ(1) > isp1ntot .or. ispXYZ(2) <= isp1ntot .or.  ispXYZ(3) <= isp1ntot )then
           print*,' problem neutron PNN ', ispXYZ(:)
           stop
       endif
  endif

  if(itx ==1)then
    if(conflag)then
        mapX => cmap2psps
        mapY => cmap2nsps
    else
        mapX   => map2psps
        mapY   => map2nsps

    endif
  else
    if(conflag)then
       mapX => cmap2nsps
       mapY => cmap2psps
    else
       mapX   => map2nsps
       mapY   => map2psps

    endif
  endif
  bsps = mapX(ispXYZ(2) )
  if(itx == 1)then
  asps = mapX(ispXYZ(1))
  csps = mapY(ispXYZ(3))

  else
  asps = mapY(ispXYZ(1))
  csps = mapX(ispXYZ(3))

  endif

!  if(ispXYZ(1) == 2 .and. ispXYZ(2) == 8 .and. ispXYZ(3) == 41)then
!      write(82,*)ispXYZ(:),asps,bsps,csps
!  endif

  if(asps < 1 .or. bsps < 1 .or. csps < 1)then

      return
  endif

  call findXXYlabel(itx,asps,bsps,csps,cindx, dindx,iphase,tphase)

  return
end subroutine getXXYindicia
!=================================================
!
!  given s.p. state indices a,b,c find the PPP or NNN index and a phase
!
! INPUT:
!  itx = 1,2  species for pair XX
!  a,b,c = native BIGSTICK label of s.p.s. (but must convert to handed labels)
!
! OUTPUT:
!  cindxtrip, dindxtrip: index of creation/destruction triplets, for fetching matrix element index
!  iphase  : phase arising from reordering into native BIGSTICK order
!  tphase  : phase arising from sending all m's -> -m
!
!
subroutine findXXYlabel(itx,a,b,c,cindxtrip,dindxtrip,iphase,tphase)

  use interactions3body
  use spstate
  use haiku_info
  use interaction
   implicit none
   integer :: itx,ity   ! species
   integer :: a, b, c
   integer :: cindxtrip,dindxtrip
   integer :: iphase
   integer :: tphase
   integer, pointer :: cxxytrip(:,:), dxxytrip(:,:)
   integer ja,jb,jc
   integer ysps
   integer(kind=8) :: xxpair


  if(itx ==1)then
         cxxytrip => cPPNtriplet
         dxxytrip => dPPNtriplet
  else
         cxxytrip => cPNNtriplet
         dxxytrip => dPNNtriplet
  endif
   iphase = 1
!...............GET TIME-REVERSED PHASE FOR A, B, C
!....... WE KNOW isp3n(1) < isp3n(2) < isp3n(3) ......
!.......... SO a is always proton, c is always neutron, and b depends on itx

   if( a > nhsps(-1))then
      ja = hspsqn(1,a -nhsps(-1))%j
   else
      ja = hspsqn(-1,a)%j
   end if
   if( b > nhsps(-itx))then
      jb = hspsqn(itx,b -nhsps(-itx))%j
   else
      jb = hspsqn(-itx,b)%j
   end if
   if( c > nhsps(-2))then
      jc = hspsqn(2,c -nhsps(-2))%j
   else
      jc = hspsqn(-2,c)%j
   end if

   tphase = (-1)**( (ja+jb+jc+1)/2)

!.............. COMPUTE PAIRS

  if(itx== 1)then
     if(a > b)then
       xxpair = (a-1)*(a-2)/2+b
     else
       xxpair = (b-1)*(b-2)/2+a
       iphase = -iphase
     endif 
     xxpair = mappairPP(xxpair)
     if(c > nhsps(-2))then
        ysps = c-nhsps(-2)
     else
        ysps = -c
     endif
  else
     if(b > c)then
       xxpair = (b-1)*(b-2)/2+c
     else
       xxpair = (c-1)*(c-2)/2+b
       iphase = -iphase
     endif 
     xxpair = mappairNN(xxpair)
     if(a > nhsps(-1))then
        ysps = a-nhsps(-1)
     else
        ysps = -a
     endif
  endif
  if(xxpair == 0)return

   cindxtrip = cxxytrip(xxpair,ysps)
   dindxtrip = dxxytrip(xxpair,ysps)


   return

end subroutine findXXYlabel
!==========================================================================
!
!  an important constraint on matrix elements because of hermiticity
!  
!  if we have a three-body XXY matrix element, V(ax,bx,cy, dx, ex,fy)
!  and ax, bx = dx,ex  but cy /= fy, then the way jumps are set up can double count;
!  therefore must divide by two
!
subroutine check_semidiag(itx,isp3i,isp3f,semidiagflag)

   implicit none
   integer itx   ! majority species
   integer isp3i(3), isp3f(3)
   logical semidiagflag
   integer yspsi, yspsf  
   integer ity

   semidiagflag = .true.

   ity = 3-itx
   if(itx == 1)then
      if(isp3i(1) /= isp3f(1) .or. isp3i(2) /= isp3f(2) )then
          semidiagflag = .false.
          return
      endif
      if(isp3i(3) == isp3f(3))then
          semidiagflag = .false.
          return
      endif 
      yspsi = isp3i(3)
      yspsf = isp3f(3)
  

   else
      if(isp3i(3) /= isp3f(3) .or. isp3i(2) /= isp3f(2) )then
          semidiagflag = .false.
          return
      endif
      if(isp3i(1) == isp3f(1))then
          semidiagflag = .false.
          return
      endif 
      yspsi = isp3i(1)
      yspsf = isp3f(1)
   endif

!........ NOW I HAVE TO CHECK A FINAL THING.......
!  do the different minority singlets belong to the same group?  (same W)?  
!  this might also depends on the species
!  NB: not entirely sure I have to check this.

   return

end subroutine check_semidiag

!==========================================================================
!
! subroutine counts up # of 3-body matrix elements in MFD-format convention file
! and extracts relevant 
!
subroutine autocount_3bmesMFD(ifile)

implicit none

integer ifile
integer nhom, nhom1sp,nhom2sp

integer(8) nme,ime
real v3me
logical finished
nme = 0

do ime = 1,1000000
   read(ifile,end=101)v3me
   nme = nme+1
end do

101 continue
print*,nme,' 3-body matrix elements '


!......... LOOP OVER POSSIBLE #s OF MATRIX ELEMENTS......

   do nhom1sp = 1,30
   call sp1nbas(nhom1sp)

        do nhom2sp = nhom1sp,2*nhom1sp
           do nhom = nhom2sp, 3*nhom1sp
           call estimate3bodymes(nhom,nhom2sp, ime)
           if(ime == nme)then
              goto 1111
           endif
           if(ime > nme)exit
        end do
   end do
end do

1111 continue

print*,' VALUES : nhom ',nhom, ', nhom1sp: ', nhom1sp,', nhom2sp: ',nhom2sp

rewind(ifile)
return
end subroutine autocount_3bmesMFD

!==========================================================================
subroutine estimate3bodymes(nhom,nhom2sp,nmes)

use spbasis_TMP
implicit none

integer nhom,nhom2sp
integer(8) nmes
integer :: mtotlim =3
integer mttot2, mjtot2, ipspi, parity

   nmes = 0
   do mttot2=-mtotlim,mtotlim,2
      do mjtot2=1,2*nhom+3,2   ! only loops over M > 0;
         do parity=1,-1,-2
         
				!!! build the list of A=3 states
            call spXYZbas(parity,nhom,nhom2sp,mjtot2,mttot2,33,'c')
            if (isp3ntot==0) cycle

            do ipspi=1,isp3ntot
               nmes = nmes + 1 !isp3ntot-ipspi+1
            end do  ! ipspi
          end do ! parity
       end do  ! mjtot2
   end do  ! mttot2

return

end subroutine estimate3bodymes

!==========================================================================
!  
!  subroutine random_3bme
!
!  purely for testing purposes, fills uncoupled 3-body matrix elements
!  with random numbers
!
   subroutine random_3bme(it)

   use system_parameters
   use interactions3body
   implicit none

   integer it
   integer(8) :: i
   real rv
   if(it == 1 .and. nmatPPP > 0)then
      do i= 1,nmatPPP
        call random_number(rv)
        hmatPPP(i) = rv-0.5
      enddo
!     print*,nmatppp,hmatppp
   endif
   if(it == 1 .and. nmatppn > 0)then
      do i= 1,nmatppn
        call random_number(rv)
        hmatppn(i) = rv-0.5
      enddo

   end if
   if(it == 2 .and. nmatNNN > 0)then
      do i= 1,nmatNNN
        call random_number(rv)
        hmatNNN(i) = rv-0.5
      enddo
   endif

   if(it == 1 .and. nmatpnn > 0)then
      do i= 1,nmatpnn
        call random_number(rv)
        hmatpnn(i) = rv-0.5
      enddo

   end if
   return

   end subroutine random_3bme

!=============================================================

