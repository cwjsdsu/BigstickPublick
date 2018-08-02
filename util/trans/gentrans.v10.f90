!================================================================2
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
! version 8: Nov 2015: backwritten as transition strength for non-charge changing operator
!   TURN OFF ISOSPIN SHIFTS
!   AUTO-READIN Z,N from file
!
! version 9:  March 2017: allows for reading in p-n formalism .opme files, plus p-n one-body densities
!================================================================

module stateinfo
   logical :: evenA
   real, allocatable,target :: energy(:)
   integer, allocatable,target :: Jx2state(:)
   integer, allocatable,target :: Tx2state(:)
   integer Mzi, Mzf    
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
   integer :: Z0,N0
   logical :: readinZN     ! flag to denote we've read in Z0 and N0
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
   logical :: pnops   ! flag for pnformalism in .opme
   real,allocatable :: p1bme(:,:),n1bme(:,:)   
   logical ::  pndens
   
   type onebdenmat
        logical good
        real, allocatable :: rho(:,:)
		real, allocatable :: rhop(:,:),rhon(:,:)  ! added in version 9
              
   end type onebdenmat

   type (onebdenmat), allocatable :: densitymats(:,:)
	   
   real :: opstrength	   

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
    print*,' General transition strengths '
    print*,' (non-charge-changing) VERSION 10 May 2018'
	readinZN = .false.
    call read_opme
	
	print*,' Enter scaling of transition operator '
	print*, '(can be = 1, but can include, e.g. osc. param. b)'
	read*,opstrength

    firstread = .true.

!.......... READ IN REFERENCE STATES.......
    print*,' '
    print*,' General transition strengths '
    print*,' (non-charge-changing) '
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
    Mzi = Z0-N0

    Mzf = Z0-N0
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
!		   print*, (tx2state_daughter(f) < abs(Mzf)),densitymats(i,f)%good
           if(tx2state_daughter(f) < abs(Mzf))cycle
           if(.not.densitymats(i,f)%good)cycle
           call tran_strength(resfile,i,f,strength,.false.)
		   strength = strength*opstrength**2
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
99 format(i4,'     ! # parents ')
    nparent = 0  ! reset
    do i = 1,nparent_ref
       if(tx2state_parent(i) < abs(Mzi))cycle
       sumstrength = 0.
	   ewsum       =0.
! .... COUNT UP # OF DAUGHTERS.......
       ndaughter = 0
       do f = 1,ndaughter_ref
           
           if(tx2state_daughter(f) < abs(Mzf))cycle
           if(.not.densitymats(i,f)%good)cycle
           call tran_strength(resfile,i,f,strength,.false.)
		   strength = strength*opstrength**2
		   
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
102 format(i4,'     ! # daughters ')

       do f = 1,ndaughter_ref
           if(tx2state_daughter(f) < abs(Mzf))cycle
           if(.not.densitymats(i,f)%good)cycle
           call tran_strength(resfile,i,f,strength,.false.)
		   strength = strength*opstrength**2
		   
           if(strength /= 0.0)then
              write(outfile,101) f,energy_daughter(f),  & 
                         0.5*jx2state_daughter(f),0.5*tx2state_daughter(f),strength
101 format(i5,f10.3,f5.1,f5.1,f10.4, '    ! daughter  Energy  J     T   strength ')
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
   logical xpn
   real cleb
   real xval


   opened = .false.
   pnops = .false.
   xpn = .false.
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
   
   select case (temp(1:3))
   
   case('iso')
   
   pnops = .false.
   print*,' Read in operator in isospin formalism '
   
   case ('pns')
   pnops = .true.
   print*,' Reading in operator in two-column proton-neutron formalism '
   
   case ('xpn')
   pnops = .true.
   xpn = .true.
   print*,' Reading in operator in single-column, explict proton-neutron formalism '   
   print*,' (Note, I will ignore charge-changing matrix elements)'
   
   case default
   
   print*,' I do not recognize the option ',temp(1:3)
   print*,' check your file format '
   stop
   
   end select
!   if(temp(1:3)=='pns')then
!	   pnops = .true.
!	   print*,' Reading in operator in proton-neutron formalism '
!  else
!	   pnops = .false.
!	   print*,' Read in operator in isospin formalism '
!   end if
   read(2,*)nsporb
   print*,' There are ',nsporb,' orbits '
   allocate(jsporb(nsporb) )
   if(xpn)then
	   
	   do i = 1,nsporb  
	     read(2,*)j,m,k,l,xj !  ! read past orbit information
		 if(m /= j+nsporb)then
			 print*,' I cannot seem to correctly interpret the s.p. labels in xpn format '
			 print*,' Please check the file '
			 stop
		 end if
			 
	     jsporb(i) = nint(2*xj)
	   end do
   else
   
   do i = 1,nsporb  
     read(2,*)j,k,l,xj !  ! read past orbit information
     jsporb(i) = nint(2*xj)
   end do
   
   end if
   if(pnops)then
      read(2,*)jt
      print*,' Transition matrix element has J = ',jt
	  tt = 2
   
   else
       read(2,*)jt,tt
       print*,' Transition matrix element has J = ',jt,', T = ',tt
! -- WARNING ABOUT SUM RULE 
       if(tt==1)then
 	     print*,' NOTE BENE! If you want to check isovector sum rules '
 	     print*,' You must undo the isospin Wigner-Eckart '
 	     print*,' (see comments at end of subroutine tran_strength) '
       end if	
   end if
   if(.not.pnops)then
      allocate( op1bme(nsporb,nsporb))
      allocate( p1bme(nsporb,nsporb),n1bme(nsporb,nsporb))  ! set these up anyway
      op1bme(:,:) = 0.0
      do i = 1,nsporb*nsporb
         read(2,*,end=11)a,b,op1bme(a,b)
		 
!.... NEW ! CONVERT!...............
         p1bme(a,b)=op1bme(a,b)* cleb(1,1,1,-1,2*tt,0)	/sqrt(2.*tt+1)	 
         n1bme(a,b)=op1bme(a,b)* cleb(1,-1,1,1,2*tt,0)	/sqrt(2.*tt+1)	 
		 
      end do
11 continue

   else
       allocate( p1bme(nsporb,nsporb),n1bme(nsporb,nsporb))
       p1bme(:,:) = 0.0
	   n1bme(:,:) = 0.0
	   
	   if(xpn)then
	       do i = 1,4*nsporb*nsporb
	          read(2,*,end=111)a,b,xval
			  if(a <=nsporb .and. b <= nsporb)p1bme(a,b)=xval
			  if(a > nsporb .and. b > nsporb)n1bme(a-nsporb,b-nsporb)=xval
			  
	       end do		   
		   
	   else
	   
       do i = 1,nsporb*nsporb
          read(2,*,end=111)a,b,p1bme(a,b),n1bme(a,b)
       end do
	   
      end if
111 continue	 
   end if
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
   use op_info,only:pndens,pnops
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
   logical askfortshift
   
   askfortshift=.false.

   print*,' reading file ',ctrlchar
   tshift = 0.0
   select case (ctrlchar)
     case ('f')
        ee => energy
        Jx2 => Jx2state
        Tx2 => Tx2state
        nmax = nlocalstates
		tshift = 0.0
		if(askfortshift)then
        print*,' Enter any shift for isospin ( x T(T+1) for all states )'
        print*,' (Typical value = 0, no shift )'
        read*,tshift
   	    end if
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
      print*,zp,np
	  if(.not.readinZN)then
		  Z0 = zp
		  N0 = np
	  else
		  if(Z0/=zp .or. N0/=np)then
			  print*,' Mismatch in valence particles '
			  print*,' Expect Z,N = ',Z0,N0
			  print*,' Found ',zp,np
			  stop
		  end if
	  end if
      if( mod(np+zp,2) == 1)then
         evenA = .false.
      else
         evenA = .true.
      end if
      rewind(resfile)
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
				  if(pndens .or. pnops)then
	
	
	
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
          if(pnops)then
              if(jmax >= Jt .and. jmin <= Jt  )then
                 densitymats(istate,fstate)%good = .true.
                 allocate(densitymats(istate,fstate)%rho( nsporb,nsporb) )
                 densitymats(istate,fstate)%rho(:,:) = 0.0
 			    allocate(densitymats(istate,fstate)%rhop(nsporb,nsporb))
    		   	    allocate(densitymats(istate,fstate)%rhon(nsporb,nsporb))
                 densitymats(istate,fstate)%rhop(:,:) = 0.0
                 densitymats(istate,fstate)%rhon(:,:) = 0.0

              else
                 densitymats(istate,fstate)%good = .false.
              end if			  
			  
		  else
             if(jmax >= Jt .and. jmin <= Jt  .and. tmax >=Tt .and. tmin <=Tt)then
                densitymats(istate,fstate)%good = .true.
                allocate(densitymats(istate,fstate)%rho( nsporb,nsporb) )
                densitymats(istate,fstate)%rho(:,:) = 0.0
			    allocate(densitymats(istate,fstate)%rhop(nsporb,nsporb))
   		   	    allocate(densitymats(istate,fstate)%rhon(nsporb,nsporb))
                densitymats(istate,fstate)%rhop(:,:) = 0.0
                densitymats(istate,fstate)%rhon(:,:) = 0.0

             else
                densitymats(istate,fstate)%good = .false.

             end if
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
!			  print*,' states ',istate,map2parent(istate),fstate,map2daughter(fstate)
      
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

11 format(a4,13x,i4)

           found = .true.
           exit
        endif
      end do

      case('f')

        read(resfile,11)myloc,n
        if(myloc(2:4)=='Fin')found=.true.


   end select
   finished = .false.
   return
111 continue
   finished = .true.
   return

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
	use stateinfo
   use op_info
   implicit none
   integer resfile
   integer istate,fstate
   integer a,b,i,j
   real ops,opv
   real fact0t,fact1t  ! isospin factors
   logical :: success
   real cleb !       ! function from LIBRA.f
   
   success=.false.
   fact0t = cleb(Tx2state_parent(istate),Mzi,0,0,Tx2state_daughter(fstate),Mzf)*sqrt(2.)/sqrt(Tx2state_daughter(fstate)+1.) 
   fact1t = cleb(Tx2state_parent(istate),Mzi,2,0,Tx2state_daughter(fstate),Mzf)*sqrt(6.)/sqrt(Tx2state_daughter(fstate)+1.) 

   do i = 1,nsporb*nsporb

          read(resfile,*,err=1,end=1)a,b,ops,opv
          if(j==jt .and. istate >0 .and. fstate > 0)then
             if(densitymats(istate,fstate)%good)then
				 if(pndens)then
					   if(ops/=0.0)then
						   densitymats(istate,fstate)%rhop(a,b)= ops
	                       success=.true.
					   end if
					   if(opv/=0.0)then
					      densitymats(istate,fstate)%rhon(a,b)= opv
                          success=.true.
					    end if
				 else
				   if(.not.pnops)then
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
!.......... CONVERT DENSITIES FROM ISOSPIN TO PROTON-NEUTRON......................
                       if(success)then
					  	densitymats(istate,fstate)%rhop(a,b) = 0.5*(fact0t*ops+fact1t*opv)
					 	densitymats(istate,fstate)%rhon(a,b) = 0.5*(fact0t*ops-fact1t*opv)
					  end if				
				  else  ! MUST CONVERT
					  if(ops/=0.0 .or. opv /=0.0)then
						  success=.true.
  					  	densitymats(istate,fstate)%rhop(a,b) = 0.5*(fact0t*ops+fact1t*opv)
  					 	densitymats(istate,fstate)%rhon(a,b) = 0.5*(fact0t*ops-fact1t*opv)						  
					  end if
					  
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

    if(.not. pndens .and. .not. pnops)then
		
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
! -- NOTE-- to check sum rules must remove isospin Wigner-Eckart
!    strength = strength/(Tx2state_parent(i)+1.) 
    
    return
else
    do a = 1,nsporb
       do b = 1,nsporb
                  strength = strength + densitymats(i,f)%rhop(a,b)*p1bme(a,b)
                  strength = strength + densitymats(i,f)%rhon(a,b)*n1bme(a,b)
!				  if(i==1)print*,i,f,a,b,densitymats(i,f)%rhon(a,b),n1bme(a,b),strength
				  
!  if reversed order
!  use: < f || [a^+ x b ] || i > = 
!   (-1)**(ji-jf +ja -jb) < i || [ b^+ x a ] || f >
!

       end do
    end do
         
    strength = strength*strength/ (Jx2state_parent(i) +1.)
    if(verbose)print*,' intermediate strength ',strength	
	return
end if

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


