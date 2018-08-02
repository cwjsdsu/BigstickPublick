!===========================================================================



!=====================================================================
!  Reads in s.p. space information;
!  Fills in q#s into derived types: orbqn in module sporbit
!
!  Reads in either .spo or.sps files
!  It first looks for .spo file; if that fails, 
!  automatically looks for .sps file of same name
!
!  AT THIS TIME ASSUMES ONLY TWO SPECIES
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
  implicit none

  integer(4)         :: ierr

!------------------FILE CONTROL---------------------------------------------

  character (len=15) :: filename
  character (len=1)  :: achar
       
  character (len=70) :: title
  integer            :: ilast
  logical            :: success


!------------------TEMP-----------------------------------------------------

  real               :: xn,xl,xj,xw    ! orbital q#s
!---------------------------------------------------------------------------
!                  dummy counters
!---------------------------------------------------------------------------

  integer            :: i,j,it

!------------------OPEN A FILE----------------------------------------------
  success = .false.

  do while ( .not. success )


     write (6,*) ' Enter file with s.p. orbit information (.sps)'
     write (6,*) ' (Enter "auto" to autofill s.p. orbit info ) '
     read (5,'(a)') filename

     ilast = index(filename,' ') - 1
     if ( filename(1:ilast)=='auto' .or. filename(1:ilast) =='AUTO' ) then
        call autofillsps
        return
     end if

!..................ATTEMPT TO OPEN .sps FILE..........................
     open (unit=1,file=filename(1:ilast)//'.sps',status='old',err=102)
     success = .true.
     cycle
102  continue
     write(6,*)filename(1:ilast),'.spo/.sps file does not exist '
  end do
!---------------- 

  write(6,*)' single-particle file = ',filename(1:ilast)
!------------------READ PAST TITLE CARDS-----------------------------

  success = .false.

     do while ( .not. success )
        read (1,'(a)') achar
        backspace(1)
        if ( achar /= '#' ) then
           success = .true.
        else
           read (1,'(a)') title
           write (6,*) title
        end if
     end do
!------------------READ PAST POSSIBLE LABEL OF ISO/PN-----------------------
  read (1,'(a)') achar

  select case(achar)
  case ('p','P')
     isoflag  = .false.
     !     write (6,*) ' .sps file in pn formalism, cannot handle ' 
     !     stop
  case ('i','I')
     isoflag = .true. 
     !              backspace(1)
  case default
     write (6,*) ' Cannot read type ',achar
     stop
  end select

!------------------READ # OF ORBITS-----------------------------------------
  if ( isoflag ) then   ! ISOSPIN FORMALISM
     read (1,*) numorb(1)

     numorbmax = numorb(1)
!------------------ALLOCATE MEMORY------------------------------------------
     allocate ( orbqn(2,numorbmax) )
!------------------READ IN--------------------------------------------------
     do i = 1, numorb(1)
           read (1,*,end=2001) xn, xl, xj, xw
           orbqn(1,i)%nr = int(xn)
           orbqn(1,i)%l = int(xl)
           orbqn(1,i)%par = (-1)**(orbqn(1,i)%l)
           orbqn(1,i)%j = int(2*xj)
           orbqn(1,i)%w = int(xw)

     end do
     numorb(2) = numorb(1)
     orbqn(2,:) = orbqn(1,:)

  else            ! pn-formalism
        read (1,*) numorb(1),numorb(2)

     numorbmax = max(numorb(1),numorb(2))
!------------------ALLOCATE MEMORY------------------------------------------
     allocate ( orbqn(2,numorbmax) )
!------------------READ IN--------------------------------------------------
     do i = 1, numorb(1)
           read (1,*,end=2001) xn, xl, xj, xw
           orbqn(1,i)%nr = int(xn)
           orbqn(1,i)%l = int(xl)
           orbqn(1,i)%par = (-1)**(orbqn(1,i)%l)
           orbqn(1,i)%j = int(2*xj)
           orbqn(1,i)%w = int(xw)
     end do

     do i = 1, numorb(2)
           read (1,*,end=2001) xn, xl, xj, xw
           orbqn(2,i)%nr = int(xn)
           orbqn(2,i)%l = int(xl)
           orbqn(2,i)%par = (-1)**(orbqn(1,i)%l)
           orbqn(2,i)%j = int(2*xj)
           orbqn(2,i)%w = int(xw)
     end do

  end if
  close(unit=1)
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
     write (6,*) ' sudden end of file in ',filename(1:ilast)
  stop  
end subroutine get_orbit_info

!=====================================================================

subroutine autofillsps

  use sporbit
  implicit none
  integer :: ierr
  integer :: nprinc
  integer :: i,j,n,l,lparity
  integer :: it


        print*,' Enter maximum principle quantum number N '
        print*,' (starting with 0s = 0, 0p = 1, 1s0d = 2, etc. ) '
        read*,nprinc

  numorb(1) = (nprinc+1)*(nprinc+2)/2
  numorb(2) = numorb(1)
  allocate( orbqn(2,numorb(1) ))
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

!  do n = 1,numorb(1)
!    print*,n, orbqn(1,n)%nr, orbqn(1,n)%l,orbqn(1,n)%j,orbqn(1,n)%w
!  enddo
  return
end subroutine autofillsps

!====================================================

  subroutine menu
  use phonon1b
  use menu_choices
  implicit none
  

  logical okay
  okay = .false.
  pnformal = .false.
  xpn =.false.
  do while(.not. okay)
      print*,' Choose one of the following '
      print*,' '
!      print*,' (0) isoscalar, scalar r^2 '
      print*,' (1) isoscalar/isovector E0 ( scalar r^2) '
      print*,' (2) isoscalar/isovector/p-n quadrupole (E2) '
      print*,' (3) proton or neutron spin '
      print*,' (4) isoscalar/vector spin-flip (sigma, not S) '
      print*,' (5) charge-changing Gamow-Teller'
      print*,' (6) isoscalar/isovector/p-n dipole (E1) '
      print*,' (7) isoscalar/isovector/p-n magnetic dipole (M1) '

      print*,' (8) ang mom vector J      '
      print*,' (9) Fermi operator tau_+/- (isospin raising/lowering) '
      
      read*,ichoice
      okay =.true.
      select case (ichoice)

!      case (0)
!        Jtrans = 0
!        Ttrans = 0
!        call r2scalarop
      case (1)
        Jtrans = 0
        print*,' Please enter T (0 for isoscalar, 1 for isovector) '
        read*,Ttrans
        call r2scalarop(ttrans)
      case (2)
        Jtrans = 2
        print*,' Please enter T (0 for isoscalar, 1 for isovector, 2 for proton-neutron) '
        read*,Ttrans		
        call quadop(ttrans)
      case (3)
        Jtrans = 1
        print*,' Please enter T ( 1 for proton only, 2 for neutron only) '
        read*,Ttrans
		Ttrans = ttrans + 1
        call spinop(ttrans)
      case (4)
        Jtrans = 1
        print*,' Please enter T (0 for isoscalar, 1 for isovector) '
        read*,Ttrans
        call spinop(ttrans)
      case (5)
        Jtrans = 1
        Ttrans = 1
        print*,' ATTENTION! I am using sigma tau_+/-'
        print*,' with the charge-changing, isospin lowering/raising operator '
        call gamowteller

      case (6)
        Jtrans = 1
        print*,' Please enter T (0 for isoscalar, 1 for isovector, 2 for proton-neutron) '
        read*,Ttrans
        call dipole(ttrans)

      case (7)
        Jtrans = 1
        print*,' Please enter T (0 for isoscalar, 1 for isovector, 2 for proton-neutron) '
        read*,Ttrans
         call m1(ttrans)
		 print*,' NOTE: If computing M1 *moments*, you must multiply by sqrt(4pi/3)=2.0466534'


      case (8)
        Jtrans = 1
        Ttrans = 0
        call jvec
      case (9)
        Jtrans = 0
        Ttrans = 1
        call Fermi

      case default
         print*,' That choice not available '
         okay = .false.
      end select

  end do ! while not okay
  return
  end subroutine menu
