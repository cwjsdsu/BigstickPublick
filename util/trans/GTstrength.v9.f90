!=================================================================
!
! code STRENGTH
!
! reads in one-body density matrices from BIGSTICK in file XXX.res
! and combines with reduced matrix elements for operator in file XXX.opme
! and computes the transition strength
!
! output is in the form B(O) where O is the operator
!
! B(O) = (average over initial)(sum over final) | < f | O | i > |^2
!
! started June 2012 by CWJ @ SDSU
!
! modified July 2012 by CWJ @ SDSU/ECT to read in multiple files
!        --> this is necessary because for certain M values the Clebsch-Gordan vanishes
!            and one needs to recalculate at a different M
!
! version 4
! modified August 2012 by CWJ @ LANL to allow for shifts due to T(T+1)
!
! version 5
! modified October 2012 by CWJ @ UW recapture "correct" ordering of states after T(T+1)
!                                     by use of reference parent, daughter information
! version 7: modified Feb 2015 to signal if a density matrix is empty
!
! version 8: modified May 2018 to allow for overall scaling
!================================================================

module stateinfo
   logical :: evenA
   real, allocatable,target :: energy(:)
   integer, allocatable,target :: Jx2state(:)
   integer, allocatable,target :: Tx2state(:)
   integer Mzi, Mzf
   integer Mzden    ! Mz of density file    
!......... THE FOLLOWING NEEDED WHEN I READ IN MULTIPLE FILES..........
   logical firstread  
   integer, allocatable :: mapstate(:)
   integer, allocatable :: map2parent(:), map2daughter(:)
   integer nlocalstates,maxlocalstates
   integer ncut
   real tshift   ! shift due to isospin
   logical reference
   integer :: nparent_ref, ndaughter_ref
   real, allocatable,target :: energy_parent(:), energy_daughter(:)
   integer, allocatable,target :: Jx2state_parent(:), Jx2state_daughter(:)
   integer, allocatable,target :: Tx2state_parent(:), Tx2state_daughter(:)

end module stateinfo

!....................................................................
! INFORMATION ON TRANSITION OPERATOR
module op_info
   implicit none

   integer nsporb   ! # of single-particle orbits
   integer Jt, Tt   ! spin, isospin of transition operator
   integer, allocatable :: jsporb(:)
   real, allocatable :: op1bme(:,:)  ! REDUCED MATRIX ELEMENTS

!........ ADDED VERSION 9 March 2017 .....
   real,allocatable :: p1bme(:,:),n1bme(:,:)   
   logical ::  pndens

   type onebdenmat
        logical good
        real, allocatable :: rho(:,:)
		real, allocatable :: rhop(:,:),rhon(:,:)  ! added in version 9
              
   end type onebdenmat

   type (onebdenmat), allocatable :: densitymats(:,:)
	   
   real :: gtscale

end module op_info

!================================================================
!
!  MAIN PROGRAM
!
    use stateinfo
    use op_info
    implicit none

    interface
       subroutine openresults(resfile,ctrlchar)
       use stateinfo
       implicit none
       integer resfile
       character(22):: filename
       character(1) :: ctrlchar
       end subroutine openresults

    end interface

    integer, parameter :: resfile = 33
    integer,parameter  :: outfile = 42
    integer tmin
    integer i,f
    real strength,sumstrength,ewsum
    integer n,z
    character ychar*1
    integer nparent,ndaughter
	
	print*,' '
	print*,' Gamow-Teller strength, version 9 '
	print*,' July 2018'

    call read_opme
	
!................

    print*,' Enter scaling (e.g., value of axial coupling g_A, quenching, etc.'
	read*,gtscale	

    firstread = .true.

!.......... READ IN REFERENCE STATES.......
    print*,' '
    print*,' For this version must first read in the energies '
    print*,' of the parents and daughter states '
    print*,' (For these, densities matrices not needed )'
    print*,' '
    print*,' PARENT *reference states*'
    call openresults(resfile,'p')
    call readheaderv2(resfile,'c')
    print*,' There are ',nlocalstates,' states '
    print*,' Enter cutoff in max number of parent states '
    print*,' (enter 0 to take all )'
    read*,nparent_ref
    if(nparent_ref < 1)nparent_ref = nlocalstates
    allocate( energy_parent(nparent_ref), Jx2state_parent(nparent_ref), Tx2state_parent(nparent_ref))
    call readheaderv2(resfile,'p')

    close(resfile)
    print*,' '
    print*,' DAUGHTER *reference states*'
    call openresults(resfile,'d')
    call readheaderv2(resfile,'c')
    print*,' There are ',nlocalstates,' states '
    print*,' Enter cutoff in max number of daughter states '
    print*,' (enter 0 to take all )'
    read*,ndaughter_ref
    if(ndaughter_ref < 1)ndaughter_ref = nlocalstates
    allocate( energy_daughter(ndaughter_ref), Jx2state_daughter(ndaughter_ref), Tx2state_daughter(ndaughter_ref))
    call readheaderv2(resfile,'d')
    close(resfile)
    call setupdensities
!..................................
    print*,'  '
    print*,' Now enter densities matrices. '
    print*,' IMPORTANT: Use largest file (most states) first '
    print*,' '
    ychar = 'y'
    do while(ychar == 'y' .or. ychar == 'Y')
       call openresults(resfile,'g')
       call readheaderv2(resfile,'c')
       if(firstread)then
         maxlocalstates = nlocalstates
         allocate( energy(nlocalstates), Jx2state(nlocalstates), Tx2state(nlocalstates))
         allocate(map2parent(nlocalstates),map2daughter(nlocalstates) )
       end if
       call readheaderv2(resfile,'f')
       call map2states
	   nparent=0
	   ndaughter=0

       do i = 1,nlocalstates
          if(map2parent(i) > 0)nparent = nparent + 1
          if(map2daughter(i) > 0)ndaughter = ndaughter + 1
		  
       end do

       print*,' This file contains ',nparent,' parents and '
       print*,' There are ',ndaughter,' daughters. '
	   
	   if(nparent == 0 .or. ndaughter==0)then
		   print*,' WARNING! To be useful a file needs BOTH parent and daughter states! '
		   if(nparent==0)   print*,' This file of densities matrices is missing parents '
		   if(ndaughter==0)   print*,' This file of densities matrices is missing daughters '

	   else
           call readalldensities(resfile)
       endif
       close(resfile)
       print*,' '
       print*,' Do you want to read another file (y/n) ? '
       read(5,'(a)')ychar
       firstread = .false.

    end do

    call prepoutput(outfile)

    print*,' Enter Z,N for initial states '
    read*,Z,N
    write(outfile,1)Z,N
1   format(2i3,'   ! Valence Z   N for parent ')
    Mzi = Z-N
    print*,' Enter Z,N for final states '
    read*,Z,N
    write(outfile,2)Z,N
2   format(2i3,'   ! Valence Z   N for daughter ')
    Mzf = Z-N

!.......COUNT UP ALLOWED PARENTS..

    nparent = 0
    ndaughter = 0
    do i = 1,nparent_ref
       if(tx2state_parent(i) >= abs(Mzi))nparent = nparent + 1
    end do
    do i = 1,ndaughter_ref
       if(tx2state_daughter(i) >= abs(Mzf))ndaughter = ndaughter + 1
    end do  

    print*,' There are max ',nparent,' parents and '
    print*,' There are max ',ndaughter,' daughters. '

!.......COUNT UP # OF NONZERO PARENTS.....
    nparent = 0
!    print*,nstates,ncut
    do i = 1,nparent_ref
       if(tx2state_parent(i) < abs(Mzi))cycle
! .... COUNT UP # OF DAUGHTERS.......
       ndaughter = 0
       do f = 1,ndaughter_ref
           if(tx2state_daughter(f) < abs(Mzf))cycle
           if(.not.densitymats(i,f)%good)cycle
           call tran_strength(resfile,i,f,strength,.false.)
           if(strength /= 0.0)then
              ndaughter = ndaughter + 1
           end if
       end do
       if(ndaughter > 0)then
             nparent=nparent +1
       else
            print*,' parent state ',i,' has no daughters '
            print*,'2 x J parent = ',Jx2state_parent(i)
       end if
       if(nparent == nparent_ref)exit
    end do
    print*,' actual # of parents is ',nparent
    write(outfile,99)nparent
99 format(i5,'    ! # parents ')
    nparent = 0  ! reset
    do i = 1,nparent_ref
       if(tx2state_parent(i) < abs(Mzi))cycle
       sumstrength = 0.
	   ewsum = 0.
! .... COUNT UP # OF DAUGHTERS.......
       ndaughter = 0
       do f = 1,ndaughter_ref
           
           if(tx2state_daughter(f) < abs(Mzf))cycle
           if(.not.densitymats(i,f)%good)cycle
           call tran_strength(resfile,i,f,strength,.false.)
           if(strength /= 0.0)then
              ndaughter = ndaughter + 1
           end if
           if(ndaughter > ndaughter_ref)exit
       end do
       if(ndaughter == 0)cycle
       nparent = nparent +1
       if(nparent > nparent_ref)exit
       write(outfile,100)i,energy_parent(i),&       
                 0.5*jx2state_parent(i),0.5*tx2state_parent(i)
100 format(i5,f10.3,f5.1,f5.1,'    ! parent  Energy  J     T')
       write(outfile,102)ndaughter
102 format(i5,'    ! # daughters ')

       do f = 1,ndaughter_ref
           if(tx2state_daughter(f) < abs(Mzf))cycle
           if(.not.densitymats(i,f)%good)cycle
           call tran_strength(resfile,i,f,strength,.false.)
		   strength = strength*gtscale**2
           if(strength /= 0.0)then
              write(outfile,101) f,energy_daughter(f),  & 
                         0.5*jx2state_daughter(f),0.5*tx2state_daughter(f),strength
101 format(i5,f10.3,f5.1,f5.1,f12.6, '  ! daughter  Energy  J     T   strength ')
              sumstrength = sumstrength+strength
			  ewsum = ewsum + strength*(energy_daughter(f)-energy_parent(i))
			  
           else
            print*,' zero strength '
            print*,i,energy_parent(i),' PARENT ',jx2state_parent(i),Tx2state_parent(i)
            print*,f, energy_daughter(f),' DAUGHTER ',jx2state_daughter(f),Tx2state_daughter(f)
!            print*,densitymats(i,f)%rho
           call tran_strength(resfile,i,f,strength,.true.)

           end if

       end do
       if(sumstrength> 0.0)write(outfile,'(i5,": Sum = ",f10.4,", EWSR = ",f10.4,", centroid = ",f10.4)') & 
	             i,sumstrength,ewsum,ewsum/sumstrength
    end do

    end

!================================================================

!............. READ IN TRANSITION MATRIX ELEMENTS.........

subroutine read_opme
   use op_info
   implicit none

   character*15 filename
   integer ilast
   character*15 temp

   integer i,j,k,l,m
   integer a,b
   real xj
   logical opened


   opened = .false.

   do while (.not.opened)

      print*,' Enter name of transition matrix file (.opme) '
      read(5,'(a)')filename
      ilast = index(filename,' ')-1
      open(unit=2,file=filename(1:ilast)//'.opme',status='old',err=1)
      opened = .true.
      exit
1     continue
      print*,filename(1:ilast),'.opme does not seem to exist '
   end do
   read(2,'(a)')temp
   read(2,*)nsporb
   print*,' There are ',nsporb,' orbits '
   allocate(jsporb(nsporb) )
   do i = 1,nsporb  
     read(2,*)j,k,l,xj !  ! read past orbit information
     jsporb(i) = nint(2*xj)
   end do
   read(2,*)jt,tt
   print*,' Transition matrix element has J = ',jt,', T = ',tt
   allocate( op1bme(nsporb,nsporb))
   op1bme(:,:) = 0.0
   do i = 1,nsporb*nsporb
       read(2,*,end=11)a,b,op1bme(a,b)
   end do
11 continue
   close(2)


   return
end subroutine read_opme

!----------------------------------------------------------
!
! routine to read in header to .res file
!
!  ctrlchar = 'c'  count up # of states
!  ctrlchar = 'f'  fill up info on states
!  ctrlchar = 'p'  fill up parent reference states
!  ctrlchar = 'd'  fill up daughter reference states

subroutine readheaderv2(resfile,ctrlchar)
   use stateinfo
   use op_info,only:pndens
   
   implicit none
   integer resfile  
   character(1) :: ctrlchar
   character(23) :: tmpline
   integer i,j,n,k
   real e,ex,xj,xt
   real etol
   real, pointer :: ee(:)
   integer, pointer :: jx2(:),tx2(:)
   integer nmax
   integer np,zp
   integer closest2J   ! function to convert poorly converged values

   print*,' reading file ',ctrlchar
   tshift = 0.0
   select case (ctrlchar)
     case ('f')
        ee => energy
        Jx2 => Jx2state
        Tx2 => Tx2state
        nmax = nlocalstates
        print*,' Enter any shift for isospin ( x T(T+1) for all states )'
        print*,' (Typical value = 0, no shift )'
        read*,tshift
!        tshift = 0.0
     case ('p')
        ee => energy_parent
        Jx2 => Jx2state_parent
        Tx2 => Tx2state_parent
        nmax = nparent_ref
     case ('d')
        ee => energy_daughter
        Jx2 => Jx2state_daughter
        Tx2 => Tx2state_daughter
        nmax = ndaughter_ref

     case default

   end select

   etol = 1.0e-3

   rewind(resfile)

!............ check whether EVEN or ODD..............
   if(ctrlchar=='p')then
      read(resfile,'(a)')tmpline
      read(resfile,'(a)')tmpline
      read(resfile,*)zp,np
	  Mzden = zp-np
      if( mod(np+zp,2) == 1)then
         evenA = .false.
      else
         evenA = .true.
      end if
      rewind(resfile)
  else
      read(resfile,'(a)')tmpline
      read(resfile,'(a)')tmpline
      read(resfile,*)zp,np
	  Mzden = zp-np	  
   end if
   do i = 1,20
      read(resfile,'(a)')tmpline
      if(tmpline(3:7)=='State')then

         select case (ctrlchar)
         case ('c')

            nlocalstates = 0
            do j = 1,50000
               read(resfile,*,err=3)n
               if(n==j)then
                   nlocalstates = nlocalstates +1
               else
                   exit
               end if
            end do

3           continue

          case ('f','p','d')
                print*,ctrlchar,nmax
              if(ctrlchar=='f' .and. nmax > maxlocalstates )then
                print*,' too many states in this file, read it first '
                stop

              end if
              do j = 1,nmax
                  read(resfile,*)n,e,ex,xj,xt
                  if(n/=j)then
                     print*,' error mismatch ',j,n
                     stop
                  end if

                  Jx2(j) = closest2J(evenA,xj)
                  Tx2(j) = closest2J(evenA,xt)
                     
!                  Jx2(j) = nint(2.*xj)
!                  Tx2(j) = nint(2.*xt)
                  ee(j) = e + tshift* 0.25* tx2(j)*(tx2(j)+2)

!.............. ERROR TRAP.....
		  if(pndens )then



                 if( abs(Jx2(j) -2*xj) > 0.1 )then
                        print*,j,'th state has not converged: j: ',xj
			  end if
		  else

			  if( abs(Jx2(j) -2*xj) > 0.1 .or. abs(Tx2(j) -2*xt) > 0.1)then
                    print*,j,'th state has not converged: j, t: ',xj,xt
!                    stop

             	  end if
		  end if
              end do

          end select
          return
      endif
   end do
   print*,' Did not find header '
   stop

end subroutine readheaderv2
!----------------------------------------------------------
subroutine map2states
   use stateinfo
   use op_info
   implicit none
   integer i,j
   real etol
   etol = 2.0e-3

!....... INITIALIZE......

   map2parent(:) = -1
   map2daughter(:) = -1

   do i = 1, nlocalstates
      do j = 1,nparent_ref
         if(abs(energy(i) - energy_parent(j))< etol .and. & 
               Jx2state(i) == Jx2state_parent(j) .and. Tx2state(i) ==Tx2state_parent(j) )then
            map2parent(i) = j
            exit
         end if

      end do ! j
      do j = 1,ndaughter_ref
         if(abs(energy(i) - energy_daughter(j))< etol .and. & 
               Jx2state(i) == Jx2state_daughter(j) .and. Tx2state(i) ==Tx2state_daughter(j) )then
            map2daughter(i) = j
            exit
         end if

      end do ! j

   end do  ! i

   return

end subroutine map2states
!----------------------------------------------------------
subroutine setupdensities
   
   use stateinfo
   use op_info
   implicit none
   integer istate,fstate
   integer jmin,jmax
   integer tmin,tmax
   allocate( densitymats(nparent_ref,ndaughter_ref))
   do istate = 1,nparent_ref
      do fstate = 1,ndaughter_ref
          jmin = nint( 0.5 * abs( jx2state_parent(istate) - jx2state_daughter(fstate) ) )
          jmax = nint( 0.5 * abs( jx2state_parent(istate) + jx2state_daughter(fstate) ) )
          tmin = nint( 0.5 * abs( tx2state_parent(istate) - tx2state_daughter(fstate) ) )
          tmax = nint( 0.5 * abs( tx2state_parent(istate) + tx2state_daughter(fstate) ) )

          if(jmax >= Jt .and. jmin <= Jt  .and. tmax >=Tt .and. tmin <=Tt)then
             densitymats(istate,fstate)%good = .true.
             allocate(densitymats(istate,fstate)%rho( nsporb,nsporb) )
             densitymats(istate,fstate)%rho(:,:) = 0.0
          else
             densitymats(istate,fstate)%good = .false.

          end if
      end do ! fstate
   end do ! istate

end subroutine setupdensities

!----------------------------------------------------------
subroutine openresults(resfile,ctrlchar)
   use stateinfo
   implicit none
   integer resfile

   character(22):: filename
   character(1) :: ctrlchar
   integer ilast

   logical success

   success = .false.
   print*,' '
   do while(.not.success)
       select case (ctrlchar)
       case('g')
       if(firstread)then
          print*,' Enter name of first results file (.res) '
       else
          print*,' Enter name of results file (.res) '
       end if
       case ('p')
          print*,' Enter name of parent results file (.res) '
       case ('d')
          print*,' Enter name of daughter results file (.res) ' 
       case default
           print*,' That case not chosen for openresults ',ctrlchar
           stop
       end select
       read(5,'(a)')filename
       ilast = index(filename,' ')-1
       open(unit=resfile,file=filename(1:ilast)//'.res',status='old',err=2)
       success = .true.
       return
2      continue
       print*,filename(1:ilast),'.res does not exist '

   end do

   return
end subroutine openresults
!----------------------------------------------------------
subroutine readalldensities(resfile)
   use op_info
   use stateinfo
   implicit none
   integer resfile
   integer istate,fstate
   integer a,b,i,j
   real ops,opv
   logical foundi,foundf,foundjt
   logical endoffile,endoflist
   logical finished,success
   logical nodensities  ! flag to make sure densities are not empty
   endoffile = .false.
   nodensities=.true.
   
!....... ADDED 3/2018 BY CWJ ....
!  CHECK THAT BOTH A PARENT AND A DAUGHTER EXIST
   

   do while(.not.endoffile)
      call read2state(resfile,'i',istate,foundi,finished)
      if(finished)exit
      if(.not.foundi)then
           endoffile = .true.
           exit
      end if

      call read2state(resfile,'f',fstate,foundf,finished)
 
      if(.not.foundf)then
           endoffile = .true.
           exit
      end if
      if(finished)exit
      endoflist = .false.
      do while(.not.endoflist)
          
          call read2Jtrans(resfile,j,foundjt)

          if(.not.foundjt)then
                endoflist = .true.
                backspace(resfile)
                exit
          end if
          if(map2parent(istate) > 0 .and. map2daughter(fstate) > 0)then
      
          call readdensity(resfile,map2parent(istate),map2daughter(fstate),j, success)
          if(success)nodensities=.false.
          end if
      end do ! endoflist

   end do  ! endoffile

   if(nodensities)then
          print*,' Wait! That density file held no densities ! '
   end if
   return
end subroutine readalldensities
!----------------------------------------------------------

subroutine prepoutput(outfile)
   implicit none
   integer outfile
   character(15) :: filename
   integer ilast


   print*,' Enter name for output (.str) file '
   read(5,'(a)')filename
   ilast = index(filename,' ')-1
   open(unit=outfile,file=filename(1:ilast)//'.str',status='unknown')
   return

end subroutine prepoutput

!----------------------------------------------------------

subroutine read2state(resfile,locchar,n,found,finished)
   implicit none
   integer resfile
   character(1) :: locchar
   integer n
   logical found,finished
   character(4) :: myloc
   character*40 :: dummyline
   integer m
   found = .false.
   select case (locchar)
     case ('i')
     do while(.not.found)
        read(resfile,1,end=111)myloc
1 format(a4)
        if(myloc(2:4)=='Ini')then
        backspace(resfile)
        read(resfile,11)myloc,n

11 format(a4,12x,i5)

           found = .true.
           exit
        endif
      end do

      case('f')

        read(resfile,11)myloc
        if(myloc(2:4)=='Fin')then
                found=.true.
                backspace(resfile)
                read(resfile,11,err=333)myloc,n
        end if


   end select
   finished = .false.
   return
111 continue
   finished = .true.
   return

333 continue
    backspace(resfile)
    read(resfile,'(a)')dummyline
    print*,dummyline
    stop

end subroutine read2state

!----------------------------------------------------------
subroutine read2Jtrans(resfile,jt,found)

   use op_info,only:pndens
   implicit none
   integer resfile
   integer jt
   logical found
   character(3) :: tmpchar
   integer j

   read(resfile,'(a3)',end=111)tmpchar
   if(tmpchar(2:3) == ' ' .or. tmpchar(2:3)=='In' .or. tmpchar(1:2)=='++')then
      found = .false.
      return
   endif
   if(tmpchar(2:3) == 'Jt')then
     backspace(resfile)
     read(resfile,'(5x,i4)')jt
   end if
   found = .true.
!..... ADDED IN VERSION 9 MAR 2017.... CHECK IF ISOSPIN OR PN FORMALISM...
   backspace(resfile)   
   read(resfile,'(11x,a3)')tmpchar
   if(tmpchar=='pro')then
	   pndens=.true.
   else
	   pndens=.false.
   end if

   return
111 continue
   found=.false.
   return   
end subroutine read2Jtrans
!----------------------------------------------------------
subroutine readdensity(resfile,istate,fstate,j,success)
   use op_info
   use stateinfo
   implicit none
   integer resfile
   integer istate,fstate
   integer a,b,i,j
   real ops,opv,rhop,rhon
   real fact0t,fact1t  ! isospin factors
   
   logical :: success
   real cleb !       ! function from LIBRA.f
   
   success=.false.
   fact0t = cleb(Tx2state_parent(istate),Mzden,0,0,Tx2state_daughter(fstate),Mzden)*sqrt(2.)/sqrt(Tx2state_daughter(fstate)+1.) 
   fact1t = cleb(Tx2state_parent(istate),Mzden,2,0,Tx2state_daughter(fstate),Mzden)*sqrt(6.)/sqrt(Tx2state_daughter(fstate)+1.) 
   success=.false.
   do i = 1,nsporb*nsporb

	   if(pndens)then
           read(resfile,*,err=1)a,b,rhop,rhon
	   else

          read(resfile,*,err=1)a,b,ops,opv
	  end if
          if(j==jt .and. istate >0 .and. fstate > 0)then
             if(densitymats(istate,fstate)%good)then
				 
				 if(pndens)then
					 ops = (rhop+rhon)/fact0t
					 opv = (rhop-rhon)/fact1t
				 end if
					 
                    if(tt == 0)then
                       if(ops /= 0.0)then
                          densitymats(istate,fstate)%rho(a,b)= ops
                          success=.true.
                       end if
                    else
                       if(opv /= 0.0)then
                           densitymats(istate,fstate)%rho(a,b) = opv
                           success=.true.
                       end if
                    end if
             end if
          end if
   end do
   return

1  continue
   backspace(resfile)
   return
end subroutine readdensity
!----------------------------------------------------------

subroutine tran_strength(resfile,i,f,strength,verbose)
    use stateinfo
    use op_info
    implicit none
    integer resfile
    integer i,f
    real strength
    logical foundi,foundf
    integer j
    integer a,b
    real cleb,tj2i !       ! function from LIBRA.f
    logical verbose

    strength = 0.0

    if(jx2state_parent(i) + jx2state_daughter(f) < 2*Jt)return
    if(abs(jx2state_parent(i) - jx2state_daughter(f)) > 2*Jt)return

    if(tx2state_parent(i) + tx2state_daughter(f) < 2*Tt)return
    if(abs(Tx2state_parent(i) - Tx2state_daughter(f)) > 2*Tt)return
    do a = 1,nsporb
       do b = 1,nsporb
                  strength = strength + densitymats(i,f)%rho(a,b)*op1bme(a,b)
!  if reversed order
!  use: < f || [a^+ x b ] || i > = 
!   (-1)**(ji-jf +ja -jb) < i || [ b^+ x a ] || f >
!

       end do
    end do
         
    strength = strength*strength/ (Jx2state_parent(i) +1.)
    if(verbose)print*,' intermediate strength ',strength
!...... CORRECT FOR ISOSPIN by undoing Wigner-Eckart...................
    strength = strength/(Tx2state_daughter(f)+1.) * ( cleb(2*Tt, Mzf-Mzi,Tx2state_parent(i),Mzi,Tx2state_daughter(f),Mzf))**2
!    strength = strength * ( TJ2I(Tx2state_parent(i),2,Tx2state_daughter(f),-Mzf,Mzf-Mzi,Mzi))**2
!     print*,f,tx2state
!   print*,f,i,Tx2state_parent(f),Tx2state_daughter(i),Mzf,Mzi, TJ2I(Tx2state_parent(f),2,Tx2state_daughter(i),-Mzf,Mzf-Mzi,Mzi)
   
    return

end subroutine tran_strength

!================================================
!
!  function to force conversion of unconverged xJ to integer J
!  that is, odd 2 x J for odd A, and even for even A
!
  function closest2J(evenA,xj)

  implicit none
  integer closest2J
  real xj
  logical evenA

  if(evenA)then
     closest2J = 2*nint(xj)
     if(closest2J < 0)closest2J = 0
  else
     closest2J = 2*nint(xj-0.5)+1
     if(closest2J < 1)closest2J = 1
  end if

  return
  end function closest2J
!================================================


