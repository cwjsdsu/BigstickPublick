!  ROUTINES TO READ IN TWO-BODY HAMILTONIAN MATRIX ELEMENTS
!
!  SUBROUTINES IN THIS FILE
!
!  master_readin_tbmes
!  readv2bme_oxbash_iso
!
! revised in 7.7.7 to allow for more flexible proton-neutron input
!
!  BIGSTICK can read in isospin conserving or isospin breaking matrix elements
!  Internally, matrix elements are always store separately for pp, nn, or pn
!  Matrix elements are stored in the module "interaction"
!  Multiple files can be read in.  Matrix elements can be read in any order,
!  with orbital labels a,b,c,d in any order. 
!
!  STORAGE CONVENTIONS
!  PP, NN stored in derived type arrays ppme(indx),nnme(indx) 
! (each of type vjs defined in module interaction in bmodules_main.f90)
!  for matrix element V_J(a,b,c,d) = < ab; J | V | cd ; J >
!  where | ab; J > is normalized two-body wavefunction (see e.g. Brussaard and Glaudemans)
!  assume a >= b, c >= d, and the pair (a,b) >= (c,d)  
!  [If not read in this order, indices swapped with phases included ]
!  TO FIND indx:  pair1 =  a*(a-1)/2 + b, pair2 = c*(c-1)/2+d
!  [ if truncations, these get mapped: pair1 = PPcouplemap(pair1) ]
!  [ where PPcouplemap stored in module interaction 
!   and created in setup_tbme_XXcouples in this file] 
!  then indx = pair1*(pair1-1)+pair2  assuming pair1 >= pair2 (swap if needed)
!  [ more specifically indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart 
!  where iref and istart are shifts depending upon parity of the pairs
!  iref = XXcouples(it)%meref(ipar),    istart = XXcouples(it)%mestart(ipar)
!   (where it = 1, 2 meaning protons or neutrons, respectively, 
!   and ipar = parity of the pair) ]
!
!  PN stored in derived type array pnme(indx). 
!   for V_J^PN(a,b,c,d)  
!   assume a,c are proton labels, b,d neutron labels
!      pair1 = numorb(2)*(a-1) + b,   pair2 = numorb(2)*(c-1) + d
!   where numorb(2) = number of neutron orbitals 
!  [ if truncations, these get mapped: pair1 = PNcouplemap(pair1) ]
!  [ where PNcouplemap stored in module interaction 
!   and created in setup_tbme_PNcouples in this file] 
!  then indx = pair1*(pair1-1)+pair2  assuming pair1 >= pair2 (swap if needed)
!  [ more specifically indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart 
!  where iref and istart are shifts depending upon parity of the pairs
!  iref = PNcouples%meref(ipar),    istart = PNcouples%mestart(ipar)
!   and ipar = parity of the pair) ]
!
!  when extracting PN from isospin formalism
!  V_PN = sqrt((1+delta_ab)*(1+delta_cd))/2  * (V_T=0 + V_T=1)
!
!  for all the derived type arrays xyme(indx)%jmin, xyme(indx)%jmax = limits of J
!  store xyme(indx)%v(j) = V_J(a,b,c,d)  with unique and compact storage
!   
!

!
!=================================================================== 
subroutine master_readin_tbmes
  use io
  use nodeinfo
  use sporbit
  use coupledmatrixelements
  use menu_choices
  use bmpi_mod
  use interaction_deformed
!  use flagger,only: subsume_spe
!  use sp_potter
  implicit none

  integer(4) :: ierr
  logical    :: finished
  integer(4) :: intfiletype
  logical(4) :: formattedfile
  integer(4) :: bcastdim = 500000

! MKGK variables
  integer :: k, indmin, indmax, jminx,jmaxx,j
  integer :: aerr
  logical :: firstread

!---------------------------------------

  firstread = .true.
  allocate( pspe(numorb(1)), nspe(numorb(2)) , stat=aerr)
  if(aerr /= 0) call memerror("master_readin_tbmes 1")
  pspe(1:numorb(1)) = 0.0
  nspe(1:numorb(2)) = 0.0
  
  allocate( psppot(numorb(1),numorb(1)), nsppot(numorb(2),numorb(2)) , stat=aerr)
  if(aerr /= 0) call memerror("master_readin_tbmes 1b")
  
  psppot(:,:)=0.0
  nsppot(:,:)= 0.0
  call setup_tbme_XXcouples(1, .false.)
  call setup_tbme_XXcouples(1, .true.)
  call set_tbme_XXarray(1)

  call setup_tbme_XXcouples(2, .false.)
  call setup_tbme_XXcouples(2, .true.)
  call set_tbme_XXarray(2)

  call setup_tbme_PNcouples( .false.)
  call setup_tbme_PNcouples( .true.)
  call set_tbme_PNarray

  if( menu_char == 'g')return

  if ( iproc == 0 ) then
	 write(logfile,*)' '
	 write(logfile,*)' Interactions read in: '
   finished = .false.  
   do while(.not.finished) 
        call tbme_menu(firstread,finished,intfiletype,formattedfile)   ! NEW ROUTINE ADDED 7.7.7
        if ( finished ) then 
           exit 
        end if 
		print*,' file type = ',intfiletype

        select case (intfiletype) 
        case (0) 
 
           call readv2bme_oxbash_iso(firstread)   ! new and improved 
        case (1,2,3) 
           call readv2bme_mfd(firstread,formattedfile,intfiletype) 
		   
		   case (4,5)
		   call readv2bme_xpn(firstread,formattedfile,intfiletype) 
        end select 
        if(.not.firstread)then
		   print*,' '
		   print*,' Interaction file successfully read '
		   print*,' '
    	end if
 
     end do ! while .not. finished 
  end if 
  call BMPI_Bcast(emptyHam,1,0,icomm,ierr) 
 
  if(deformed)call deformedSPEs 
!  if(subsume_spe)call subsume_sp_pot_into_2body
  call broadcastTBMES 
 
  return 
end subroutine master_readin_tbmes 
!============================================================ 
!
! added 7.7.7 to give more flexibility in reading in files
! revision of open_tbme_file
!
!  do until either finished or file is opened:
!     read in filestring
!     if filestring == def or iso or mfd or xpn or upn
!        set intfiletype and cycle
!     else
!        attempt to open file
!        
!
!
!  OUTPUT:
!     finished: if true, then no more interaction files
!     intfiletype: type of file
!     formattedfile: if formatted or unformatted
!
!
subroutine tbme_menu(firstread,finished,intfiletype,formattedfile) 
    use io 
    use sporbit 
    use coupledmatrixelements 
    use nodeinfo 
    use bmpi_mod 
    use reporter
    implicit none 
 
    logical      :: firstread
    logical      :: finished 
    logical      :: formattedfile 
    integer      :: intfiletype 
!------------ INTERNAL -----------------------	 
    integer(4)   :: ierr 
	logical      :: success
    character*55 :: tmpfilename 
    character*3  :: formatchar
	integer      :: ilast
    logical      ::  gotit
    character   breakchar 
	integer      :: nme
    integer      :: idot 
	
if(iproc /= 0)return
formatchar='def'
success = .false.
do while(.not.success) 
    if(auto_input)then 
       read(autoinputfile,'(a)')tmpfilename 
    else 
		if(.not.firstread)print*,' Would you like to read in another interaction file? '
		select case (formatchar)
	
	     case ('mfd')
              print*,' Enter name of two-body interaction file in MFD format ' 
	
 	     case ('upn','xpn')
             print*,' Enter name of two-body interaction file in explicit proton-neutron format ' 
	
		 case('def','iso')
              print*,' Enter two-body interaction file name OR file format code (e.g., XPN) ' 
              print*,' (Enter "end" to finish; "opt" for file format options; "?" for  general info ) ' 
	     
		 case default
		    print*,' variable formatchar not properly set ',formatchar
			print*,' STOPPING RUN '
		    call BMPI_Abort(icomm,101,ierr)
			
		    stop
		
	    end select
        read(5,'(a)')tmpfilename 
    end  if 
	if(tmpfilename(1:3)=='opt')then
		print*,' '
		print*,' BIGSTICK TBME FILE FORMAT OPTIONS '
		print*,' (for more details see the Inside Guide with your distribution )'
		print*,' (iso) Isospin formalism, similar to OXBASH/NuShell format. DEFAULT '
		print*,' (mfd) MFD format with further options for isospin breaking. DEFAULT '
		print*,' (xpn) Explict p-n formalism (with normalized p-n matrix elements )'
		print*,' (upn) Like xpn, but with so-called unnormalized p-n matrix elements'
		print*,'  '
		print*,' DEFAULTS '
		print*,' If no file format specified,  BIGSTICK will assume iso/mfd format '
		print*,' Suppose you enter XXX as the file name.  '
		print*,' BIGSTICK will FIRST look for file XXX.int; if found, iso format assumed '
		print*,' If BIGSTICK cannot find XXX.int, it will look for file XXX with mfd format '
		print*,' '
		print*,' You can, however, use an MFD format file of name XXX.int: '
		print*,' Enter in full XXX.int as the file name.'
		print*,' '
		print*,' You may enter in explicit proton/neutron matrix elements in iso-like format '
		print*,' If the files have the form XXX.pp.int, XXX.nn.int, XXX.pn.int  '
		print*,' OR XXX.isos.int, XXX.isov.int, or XXX.isot.int, '
		print*,' Then specify XXX.pp etc when requested. '
		print*,' '
		print*,' For xpn/upn formats, you MUST specify the format. '
		print*,' In this format proton and neutron orbits are sequential and do not overlap, '
		print*,' E.g., proton orbits are 1,2,3 and neutron orbits are 4,5,6. '
		print*,' FOR NOW despite the distinct numbering the proton and neutron orbits '
		print*,' must encompass the same space. '
		print*,' NOTE:  upn format is typical for TBME files distributed with NuShell; '
		print*,' xpn/upn files must have the name XXX.int, but enter XXX when requested. '
		print*,' '			
		print*,' (Again, for more details see the Inside Guide with your distribution )'
		cycle
		
	end if
    if(tmpfilename(1:1) == '?')then 
       print*,' ' 
       print*,' BIGSTICK supports several interaction file conventions ' 
       print*,' and tries to automatically distinguish between them. ' 
       print*,' If using OXBASH-format input, leave off the .int suffix ' 
       print*,' If using MFD-format input, enter the entire filename (including suffix if any) ' 
       print*,' ' 
       write(6,*)' ' 
       write(6,*)' You can set a path to a standard repository of interaction files ' 
       write(6,*)' by using the environmental variable BIG_INT_DIR. ' 
       write(6,*)' Just do : ' 
       write(6,*)' export BIG_SPS_DIR = (directory name) ' 
       write(6,*)' export BIG_SPS_DIR=/Users/myname/int_repo ' 
          
       if(INT_path=='/')then 
              print*,' Currently BIG_INT_DIR is not set ' 
       else 
              print*,' Currently BIG_INT_DIR = ',int_path(:) 
              print*,' However the current working directory is always checked first ' 
       end if 

       cycle 
    endif
    if ( tmpfilename(1:3)== 'END' .or. tmpfilename(1:3) == 'end' ) then 
          finished = .true.
          if(.not.auto_input)write(autoinputfile,'(a)')'end'
		  
          return 
     end if 
	 
	 select case (tmpfilename(1:3))
	    case ('def','DEF')
		  formatchar = 'def'
		  if(.not.auto_input)write(autoinputfile,'(a)')'def'
	      cycle
  	    case ('iso','ISO')
  		  formatchar = 'iso'
		  if(.not.auto_input)write(autoinputfile,'(a)')'iso'
  	      cycle	 
  	    case ('xpn','XPN')
  		  formatchar = 'xpn'
		  if(.not.auto_input)write(autoinputfile,'(a)')'xpn'
  	      cycle	 
  	    case ('upn','UPN')
  		  formatchar = 'upn'
		  if(.not.auto_input)write(autoinputfile,'(a)')'upn'
  	      cycle
  	    case ('mfd','MFD')
  		  formatchar = 'mfd'
		  if(.not.auto_input)write(autoinputfile,'(a)')'mfd'
		  cycle
  	      		  		  
     end select

!	if(tmpfilename(1:3)=='def' .or. tmpfilename(1:3)=='iso' .or. tmpfilename(1:3)=='xpn' .or. tmpfilename(1:3)=='upn')cycle
	
    write(6,*)' Reading input file ',tmpfilename ,formatchar
    intfilename = tmpfilename 
    ilast = index(intfilename,' ')-1 		
!...........ATTEMPT TO OPEN .int FILE.................. 

     inquire(file=intfilename(1:ilast)//'.int',exist=success) 
     if(success)then 
        open(unit=1,file=intfilename(1:ilast)//'.int',status='old') 
		
		select case (formatchar)
		   case ('def','iso')
  	          write(6,*)' opening OXBASH/NuSHELL-style isospin formatted file ',intfilename(1:ilast)//'.int'
	  	      write(logfile,*)' opening OXBASH/NuSHELL-style isospin formatted file ',intfilename(1:ilast)//'.int'
			  formatchar = 'iso'
			  
		   case ('xpn')
  	          write(logfile,*)' opening explicit proton-neutron file ',intfilename(1:ilast)//'.int'
		      write(logfile,*)' (with normalized pn matrix elements )'
  	          write(6,*)' opening explicit proton-neutron file ',intfilename(1:ilast)//'.int'
		      write(6,*)' (with normalized pn matrix elements )'
			  
   		   case ('upn')
     	      write(logfile,*)' opening explicit proton-neutron file ',intfilename(1:ilast)//'.int'
   		      write(logfile,*)' (with unnormalized pn matrix elements )'
 	          write(6,*)' opening explicit proton-neutron file ',intfilename(1:ilast)//'.int'
	          write(6,*)' (with unnormalized pn matrix elements )'
			  
			  case default 
			  print*, ' should not have gotten to this point, bad file format choice ',formatchar
  			  print*,' STOPPING RUN '
  		      call BMPI_Abort(icomm,101,ierr)
			  stop
			  
  	     end select
      else 
         inquire(file=int_path(1:length_int_path)//intfilename(1:ilast)//'.int',exist=success) 
	     if(success)then 
            open(unit=1,file=int_path(1:length_int_path)//intfilename(1:ilast)//'.int',status='old') 
   	   		select case (formatchar)
   	   		   case ('def','iso')
   		          write(logfile,*)' opening OXBASH/NuSHELL-style file ',int_path(1:length_int_path)//intfilename(1:ilast)//'.int'
			  
	   		   case ('xpn')
	     	      write(logfile,*)' opening explicit proton-neutron file ',int_path(1:length_int_path)//intfilename(1:ilast)//'.int'
	   		      write(logfile,*)' (with normalized pn matrix elements )'
	      	   case ('upn')
	        	   write(logfile,*)' opening explicit proton-neutron file ',int_path(1:length_int_path)//intfilename(1:ilast)//'.int'
	      		   write(logfile,*)' (with unnormalized pn matrix elements )'
			  
	   		    case default 
	   			  print*, ' should not have gotten to this point, bad file format choice ',formatchar
	   			  stop
			  
	     	  end select
			   
			   
	      end if  ! success
      end if
	  if(success)then
		  if(.not.auto_input)write(autoinputfile,'(a)')intfilename(1:ilast)

		  emptyham =.false.
		  if(success)then
		      select case (formatchar) 
		      case ('def','iso')
		         intfiletype = 0
		         pp_int = .false. 
		         nn_int = .false. 
		         pn_int = .false. 
		         iso_int = .true.            !   Default interaction type 
		         isov_int = .false. 
		         isot_int = .false. 
		 	     coul_int = .false. 
		 	     mixed_T = .false. 
		         idot=0 
		         idot = index(intfilename,'.') 
		         if( idot > 0)then 
		            if(intfilename(idot+1:ilast) == 'pp')then          !  proton-proton interation 
		               pp_int = .true. 
		                iso_int = .false. 
		 			   write(logfile,*)' Reading only proton-proton matrix elements '
		            elseif(intfilename(idot+1:ilast) == 'nn')then      !  peutron-neutron interation 
		               nn_int = .true. 
		               iso_int = .false. 
		 		      write(logfile,*)' Reading only neutron-neutron matrix elements '
		            elseif(intfilename(idot+1:ilast) == 'pn')then      !  proton-neutron interation 
		               pn_int = .true. 
		               iso_int = .false. 
		 		      write(logfile,*)' Reading only proton-neutron matrix elements '
		            elseif(intfilename(idot+1:ilast) == 'iso')then     !  isoscalar interation 
		               iso_int = .true. 
		            elseif(intfilename(idot+1:ilast) == 'isov')then    !  isovector interation 
		               isov_int = .true. 
		               iso_int = .false. 
		            elseif(intfilename(idot+1:ilast) == 'isot')then    !  isotensor interation 
		               isot_int = .true. 
		               iso_int = .false. 
		            elseif(intfilename(idot+1:ilast) == 'coul')then    !  Coulomb scaled by hbw interation 
		               coul_int = .true. 
		               iso_int = .false. 
		            elseif(intfilename(idot+1:ilast) == 'mixed')then   !  Broken isospin <T=0|V|T=1> interation 
		               mixed_T = .true. 
		               iso_int = .false. 
		            end if 
		         end if 				 
	
		  	case ('xpn')
	
		  	intfiletype =  4
	
		  	case ('upn')
	
		  	intfiletype = 5
	
		  	case ('mfd')
	
		  	print*,' Expecting an MFD format file, should not have gotten here '
			print*,' STOPPING RUN '
		    call BMPI_Abort(icomm,101,ierr)
		  	stop
	
		  	case default
	
		  	print*,formatchar,' is not a recognized response'
		  	cycle
	
		      end select	
	
	
		  end if
	  end if
	  if(success)exit
	  if(.not.success .and.(formatchar /='mfd' .and. formatchar /='def') )then
		  write(6,*)' file ',intfilename(1:ilast)//'.int does not exist '
		  cycle
	  end if
	  if(.not.success)then  ! try to open a MFD file 	
	      inquire(file=intfilename(1:ilast),exist=success) 
	      if(success)then 
	         open(unit=1,file=intfilename(1:ilast),status='old') 
	 		 write(6,*)' opening MFD-style file ',intfilename(1:ilast)

	 		 write(logfile,*)' opening MFD-style file ',intfilename(1:ilast)
		
	      else 
	         inquire(file=int_path(1:length_int_path)//intfilename(1:ilast),exist=success) 
	 	      if(success)then 
	            open(unit=1,file=int_path(1:length_int_path)//intfilename(1:ilast),status='old') 
	    		write(6,*)' opening MFD-style file ',int_path(1:length_int_path)//intfilename(1:ilast)

	    		write(logfile,*)' opening MFD-style file ',int_path(1:length_int_path)//intfilename(1:ilast)
	 	      end if 
	      end if 		  
	  end if   ! not success
	  if(success)then   ! successfully opened mfd-style format
		  if(.not.auto_input)write(autoinputfile,'(a)')intfilename(1:ilast)
		  
		  formatchar='mfd'
		  emptyham=.false.
!........... CHECK WHETHER FORMATTED OR UNFORMATTED..... 
!  a bit kludgy... could be cleaned up 
! 
          formattedfile = .true.  ! default assumption
          gotit = .false. 
          do while(.not.gotit) 
             read(1,'(a)')breakchar 
             print*,breakchar 
             if(breakchar/='!')then 
                backspace(1) 
                gotit = .true. 
             endif 
          enddo 
          read(1,*,err = 202) nme 
          rewind(1) 
          goto 203 
202    continue 
          formattedfile = .false. 
          print*,' Unformatted MFD interaction file' 
	      write(logfile,*) ' Unformatted MFD interaction file' 
          close(1) 
          open(unit=1,file=intfilename(1:ilast),status='old', form='unformatted') 
203    continue		  
	      gotit = .false. 
	      intfiletype = -1 
	      do while (.not. gotit) 
 
	         print*,' For MFD-formatted input choose one of the following :' 
	         print*,' (I) No isospin breaking ' 
	         print*,' (P) Explicit proton-neutron formalism ' 
	         print*,' (C) Isospin breaking only through adding Coulomb ' 
	         if(auto_input)then 
	           read(autoinputfile,'(a)')breakchar 
	            if(breakchar == 'i')breakchar = 'I' 
	            if(breakchar == 'p')breakchar = 'P' 
	            if(breakchar == 'c')breakchar = 'C' 
	         else 
	            read(5,'(a)')breakchar 
	            if(breakchar == 'i')breakchar = 'I' 
	            if(breakchar == 'p')breakchar = 'P' 
	            if(breakchar == 'c')breakchar = 'C' 
 
	           write(autoinputfile,'(a)')breakchar 
	         endif 
	         if(breakchar =='I' ) intfiletype = 1 
	         if(breakchar == 'P' ) intfiletype = 3 
	         if(breakchar =='C') intfiletype = 2 
	         if(breakchar == 'P' .or. breakchar == 'C' .or. breakchar=='I')gotit = .true. 
          
	 !......... FOR OUTPUT MAY FLIP ISOFLAG TO FALSE..... 
	         if(breakchar == 'P' .or. breakchar == 'C')isoflag = .false. 
	      end do 
	      if(breakchar =='I' )write(logfile,*)' (I) No isospin breaking ' 
	      if(breakchar =='P' )write(logfile,*)' (P) Explicit proton-neutron formalism ' 
	      if(breakchar =='C' )write(logfile,*)' (C) Isospin breaking only through adding Coulomb ' 
		  
	  else
		  write(6,*)' file ',intfilename(1:ilast),' does not exist '
		  

	  end if
	
end do  ! while not success ------------------------------

!  AT THIS POINT EITHER A FILE IS OPENED, OR WE HAVE ENDED ENTERING FILES.....

if(.not.success .and. .not.finished)then
	print*,' some confusion on logic '
	stop
end if
	
	return
end subroutine tbme_menu 


!============================================================ 
 
subroutine open_tbme_file(finished,intfiletype,formattedfile) 
 
  use io 
  use sporbit 
  use coupledmatrixelements 
  use nodeinfo 
  use bmpi_mod 
  use reporter
  implicit none 
 
  integer(4)   :: ierr 
 
  logical      :: finished 
  logical      :: formattedfile 
  integer      :: intfiletype 
!----------------- INTERNAL VARIABLES --------------- 
!------ FILE CONTROL --------------------------- 
 
  integer      :: ilast,nme 
  logical      :: success, gotit
  character   breakchar 
  character*55 :: tmpfilename 
  integer      :: idot 
  character*3  :: formatchar
 
  if(iproc == 0)then 
  success = .false. 
  formattedfile = .true.  !default assumption 
  formatchar   = 'def'  
  
  
  do while(.not.success) 
 
     if(auto_input)then 
        read(autoinputfile,'(a)')tmpfilename 
     else 
 		select case (formatchar)
		
		case ('mfd')
           print*,' Enter name of two-body interaction file in MFD format ' 
		
		case ('upn','xpn')
           print*,' Enter name of two-body interaction file in explicit proton-neutron format ' 
		

        case default
            print*,' Enter two-body interaction file name ' 
            print*,' (Enter "end" to stop; "opt" for file format options; "?" for  general info ) ' 
			
		end select
        read(5,'(a)')tmpfilename 
		
		if(tmpfilename(1:3)=='opt')then
			print*,' '
			print*,' BIGSTICK TBME FILE FORMAT OPTIONS '
			print*,' (for more details see the Inside Guide with your distribution )'
			print*,' (iso) Isospin formalism, similar to OXBASH/NuShell format. DEFAULT '
			print*,' (mfd) MFD format with further options for isospin breaking. DEFAULT '
			print*,' (xpn) Explict p-n formalism (with normalized p-n matrix elements )'
			print*,' (upn) Like xpn, but with so-called unnormalized p-n matrix elements'
			print*,'  '
			print*,' DEFAULTS '
			print*,' If no file format specified,  BIGSTICK will assume iso/mfd format '
			print*,' Suppose you enter XXX as the file name.  '
			print*,' BIGSTICK will FIRST look for file XXX.int; if found, iso format assumed '
			print*,' If BIGSTICK cannot find XXX.int, it will look for file XXX with mfd format '
			print*,' '
			print*,' You can, however, use an MFD format file of name XXX.int: '
			print*,' Enter in full XXX.int as the file name.'
			print*,' '
			print*,' You may enter in explicit proton/neutron matrix elements in iso-like format '
			print*,' If the files have the form XXX.pp.int, XXX.nn.int, XXX.pn.int  '
			print*,' OR XXX.isos.int, XXX.isov.int, or XXX.isot.int, '
			print*,' Then specify XXX.pp etc when requested. '
			print*,' '
			print*,' For xpn/upn formats, you MUST specify the format. '
			print*,' In this format proton and neutron orbits are sequential and do not overlap, '
			print*,' E.g., proton orbits are 1,2,3 and neutron orbits are 4,5,6. '
			print*,' FOR NOW despite the distinct numbering the proton and neutron orbits '
			print*,' must encompass the same space. '
			print*,' NOTE:  upn format is typical for TBME files distributed with NuShell; '
			print*,' xpn/upn files must have the name XXX.int, but enter XXX when requested. '
			print*,' '			
			print*,' (Again, for more details see the Inside Guide with your distribution )'
			cycle
			
		end if
        if(tmpfilename(1:1) == '?')then 
           print*,' ' 
           print*,' BIGSTICK supports several interaction file conventions ' 
           print*,' and tries to automatically distinguish between them. ' 
           print*,' If using OXBASH-format input, leave off the .int suffix ' 
           print*,' If using MFD-format input, enter the entire filename (including suffix if any) ' 
           print*,' ' 
           write(6,*)' ' 
           write(6,*)' You can set a path to a standard repository of interaction files ' 
           write(6,*)' by using the environmental variable BIG_INT_DIR. ' 
           write(6,*)' Just do : ' 
           write(6,*)' export BIG_SPS_DIR = (directory name) ' 
           write(6,*)' export BIG_SPS_DIR=/Users/myname/int_repo ' 
              
           if(INT_path=='/')then 
                  print*,' Currently BIG_INT_DIR is not set ' 
           else 
                  print*,' Currently BIG_INT_DIR = ',int_path(:) 
                  print*,' However the current working directory is always checked first ' 
           end if 
 
           cycle 
        else 
           if(tmpfilename(1:3)/='end')then 
               write(autoinputfile,'(a,"    ! name of interaction file ")')tmpfilename 
           else 
               write(autoinputfile,'(a)')'end       ! end of interaction files ' 
           end if 
        endif  
     endif 
	 

     if ( tmpfilename(1:3)== 'END' .or. tmpfilename(1:3) == 'end' ) then 
 
           finished = .true. 
           return 
		   
	       if( tmpfilename(1:3)=='mfd' .or.tmpfilename(1:3)=='iso' .or. tmpfilename(1:3)=='xpn' .or. tmpfilename(1:3)=='upn' )then
			   write(6,*)' Format choice ',tmpfilename(1:3)
			   
			   
		   else
        write(6,*)' Reading input file ',tmpfilename 
	endif
      end if 
     intfilename = tmpfilename 
     ilast = index(intfilename,' ')-1 
 
!...........ATTEMPT TO OPEN .int FILE.................. 

         inquire(file=intfilename(1:ilast)//'.int',exist=success) 
         if(success)then 
           open(unit=1,file=intfilename(1:ilast)//'.int',status='old') 
		   select case (formatchar)
		   case ('def','iso')
	  	      write(logfile,*)' opening OXBASH/NuSHELL-style file ',intfilename(1:ilast)//'.int'
			  
		   case ('xpn')
  	          write(logfile,*)' opening explicit proton-neutron file ',intfilename(1:ilast)//'.int'
		      write(logfile,*)' (with normalized pn matrix elements )'
   		   case ('upn')
     	          write(logfile,*)' opening explicit proton-neutron file ',intfilename(1:ilast)//'.int'
   		      write(logfile,*)' (with unnormalized pn matrix elements )'
			  
			  case default 
			  print*, ' should not have gotten to this point, bad file format choice ',formatchar
			  stop
			  
  	       end select
        else 
            inquire(file=int_path(1:length_int_path)//intfilename(1:ilast)//'.int',exist=success) 
	        if(success)then 
               open(unit=1,file=int_path(1:length_int_path)//intfilename(1:ilast)//'.int',status='old') 
   		       write(logfile,*)' opening OXBASH/NuSHELL-style file ',int_path(1:length_int_path)//intfilename(1:ilast)//'.int'
	         end if 
         end if 
     if(.not.success) goto 102	   
     emptyHam = .false.   ! successfully opened a file 

     select case (formatchar) 
	    case ('def','iso')
	 
        intfiletype = 0
		
		case ('xpn')
		
		intfiletype =  4
		
		case ('upn')
		
		intfiletype = 5
		
		case ('mfd')
		
		print*,' Expecting an MFD format file, should not have gotten here '
		stop
		
		case default
		
		print*,formatchar,' is not a recognized response'
		cycle
		
     end select
	 print*,intfiletype

 
!...........Check if extensions control interaction type while in 
!...........isospin formalism 
     if(isoflag) then 
        pp_int = .false. 
        nn_int = .false. 
        pn_int = .false. 
        iso_int = .true.            !   Default interaction type 
        isov_int = .false. 
        isot_int = .false. 
	coul_int = .false. 
	mixed_T = .false. 
        idot=0 
        idot = index(intfilename,'.') 
        if( idot > 0)then 
           if(intfilename(idot+1:ilast) == 'pp')then          !  proton-proton interation 
              pp_int = .true. 
               iso_int = .false. 
			   write(logfile,*)' Reading only proton-proton matrix elements '
           elseif(intfilename(idot+1:ilast) == 'nn')then      !  peutron-neutron interation 
              nn_int = .true. 
              iso_int = .false. 
		      write(logfile,*)' Reading only neutron-neutron matrix elements '
           elseif(intfilename(idot+1:ilast) == 'pn')then      !  proton-neutron interation 
              pn_int = .true. 
              iso_int = .false. 
		      write(logfile,*)' Reading only proton-neutron matrix elements '
           elseif(intfilename(idot+1:ilast) == 'iso')then     !  isoscalar interation 
              iso_int = .true. 
           elseif(intfilename(idot+1:ilast) == 'isov')then    !  isovector interation 
              isov_int = .true. 
              iso_int = .false. 
           elseif(intfilename(idot+1:ilast) == 'isot')then    !  isotensor interation 
              isot_int = .true. 
              iso_int = .false. 
           elseif(intfilename(idot+1:ilast) == 'coul')then    !  Coulomb scaled by hbw interation 
              coul_int = .true. 
              iso_int = .false. 
           elseif(intfilename(idot+1:ilast) == 'mixed')then   !  Broken isospin <T=0|V|T=1> interation 
              mixed_T = .true. 
              iso_int = .false. 
           end if 
        end if 
     end if  
! end if
     cycle 
102  continue     


!if(formatchar == 'def' .or. formatchar == 'mfd')then
     inquire(file=intfilename(1:ilast),exist=success) 
     if(success)then 
        open(unit=1,file=intfilename(1:ilast),status='old') 
		write(logfile,*)' opening MFD-style file ',intfilename(1:ilast)
		
     else 
     inquire(file=int_path(1:length_int_path)//intfilename(1:ilast),exist=success) 
	 if(success)then 
           open(unit=1,file=int_path(1:length_int_path)//intfilename(1:ilast),status='old') 
   		   write(logfile,*)' opening MFD-style file ',int_path(1:length_int_path)//intfilename(1:ilast)
	 end if 
     end if 
! end if

     if(.not.success) goto 103	   
     emptyHam = .false.   ! successfully opened a file 
 
!........... CHECK WHETHER FORMATTED OR UNFORMATTED..... 
!  a bit kludgy... could be cleaned up 
! 
     gotit = .false. 
     do while(.not.gotit) 
        read(1,'(a)')breakchar 
        print*,breakchar 
        if(breakchar/='!')then 
           backspace(1) 
           gotit = .true. 
        endif 
     enddo 
     read(1,*,err = 202) nme 
     rewind(1) 
     goto 3 
202    continue 
     formattedfile = .false. 
     print*,' Unformatted MFD interaction file' 
	  write(logfile,*) ' Unformatted MFD interaction file' 
     close(1) 
     open(unit=1,file=intfilename(1:ilast),status='old', form='unformatted') 
 
!............. HAVE NOW OPENED AN MFD-formatted file 
3    continue 
     gotit = .false. 
     intfiletype = -1 
     do while (.not. gotit) 
 
        print*,' For MFD-formatted input choose one of the following :' 
        print*,' (I) No isospin breaking ' 
        print*,' (P) Explicit proton-neutron formalism ' 
        print*,' (C) Isospin breaking only through adding Coulomb ' 
        if(auto_input)then 
          read(autoinputfile,'(a)')breakchar 
           if(breakchar == 'i')breakchar = 'I' 
           if(breakchar == 'p')breakchar = 'P' 
           if(breakchar == 'c')breakchar = 'C' 
        else 
           read(5,'(a)')breakchar 
           if(breakchar == 'i')breakchar = 'I' 
           if(breakchar == 'p')breakchar = 'P' 
           if(breakchar == 'c')breakchar = 'C' 
 
          write(autoinputfile,'(a)')breakchar 
        endif 
        if(breakchar =='I' ) intfiletype = 1 
        if(breakchar == 'P' ) intfiletype = 3 
        if(breakchar =='C') intfiletype = 2 
        if(breakchar == 'P' .or. breakchar == 'C' .or. breakchar=='I')gotit = .true. 
          
!......... FOR OUTPUT MAY FLIP ISOFLAG TO FALSE..... 
        if(breakchar == 'P' .or. breakchar == 'C')isoflag = .false. 
     end do 
     if(breakchar =='I' )write(logfile,*)' (I) No isospin breaking ' 
     if(breakchar =='P' )write(logfile,*)' (P) Explicit proton-neutron formalism ' 
     if(breakchar =='C' )write(logfile,*)' (C) Isospin breaking only through adding Coulomb ' 
     cycle 
	 
103  continue 
     print*,intfilename(1:ilast),' TBME file does not exist ' 
  end do  ! loop over success 
  end if  !  End check on iproc == 0 to only read in int files 
 
  return 
end subroutine open_tbme_file 
!============================================================ 
subroutine readv2bme_oxbash_iso(firstread)
! 
!   READS IN single-particle energies, TBMEs 
!   CAN READ IN MULTIPLE FILES 
!   STORES IN DERIVED TYPE vtbme 
!   THE TBMEs V_JT(ab,cd) are ordered in a particular way: 
!   a >= b, c >= d;  a >= c 
!     VTBME(indx): 
!      indx =  pair1*(pair1-1)/2+pair2 
!    WHERE     pair1 = ia*(ia-1)/2 + ib 
!              pair2 = ic*(ic-1)/2 + id 
! 
  use system_parameters 
  use io 
  use sporbit 
  use W_info 
  use coupledmatrixelements 
  use nodeinfo 
  use flagger
  use bmpi_mod 
  use butil_mod
  
  implicit none 
  logical firstread
 
  integer(4) :: ierr 
  integer :: aerr 
 
!------------------------------ 
  character*1 ychar 
   
  character*70 title 
  character*3 :: formattest
  integer ia,ib,ic,id 
  integer a,b,c,d 
  integer na,nb,nc,nd 
  integer pair1,pair2,indx 
  integer j, t, t_ab, t_cd 
  real Vv 
  real, allocatable :: spetmp(:) 
  real spscale,ax,bx,x,vscale  ! for scaling interactions 
  real, allocatable :: ratio(:) 
  real hbw0, hbw 
 
  integer phase 
  integer nme 
  real factor 
  integer ipar,iref,istart 
  integer it 
 
  real(kind=4) :: factor_v,sc2
  integer :: maxsp    ! used for error trap  
!----------------------------------------------- 
!      dummy counters 
!-------------------------------------------------- 
  integer i,tmp,m,mmax 
     
  integer L 
  logical success 
  logical smint   !  a successor to .int  
  logical finished ,foundit ,formattedfile
  logical autoscale  ! added 7.8.1 for reading in NuShell files
  integer :: numval
  
!--------- KSM - initialize for -Wuninitialized 
! set to values that would cause problems if used 
  it = -10000000 
  factor_v = -1e20 
!---------------BEGIN ------------------------- 
!  if ( .not. isoflag ) then 
!     print*,' Whoa, not set up yet for pn input, friend! ' 
!     stop 
!  end if 
 
!-------------- READ PAST TITLE CARDS --------------------------- 
  rewind(1)
  success = .false. 
  do while(.not.success)
     read(1,'(a)')ychar 
     backspace(1) 
     if(ychar /= '#' .and. ychar /= '!')then 
        success = .true. 
     else 
        read(1,'(a)')title 
        write(6,*)title 
		write(logfile,*)title
		call fileformatcheck(title,foundit,formattest)
		if(foundit)then ! check format
		   if(( formattest /='def' .and. formattest/='iso')  )then
			   print*,' '
			   print*,'  WARNING !   '
			  print*,' File appears to be in wrong format '
			  print*,' Header denotes format as ',formattest
		      print*,'  SKIPPING FILE !   '
  		      print*,' '
			  write(logfile,*)' File appears to be in wrong format '
			  write(logfile,*)' Header denotes format as ',formattest
		      write(logfile,*)'  SKIPPING FILE !   '
			  return
		  end if
	    end if
     end if 
  end do 
!--------------- IF THIS IS A .SMINT FILE, COMPARE S.P. STATES  
!       not yet implemented 
 
!---------- CHECK FOR AUTOSCALING ---------

   read(1,*)nme
   if(nme < 0)then
	   autoscale = .true.
   else
	   autoscale = .false.
   end if
   backspace(1) 
!-------------- ENTER SCALING ----------------------------- 
  sc2 = 1.0 
if(.not.autoscale)then
  if(auto_input)then 
     if(coul_int)then 
        print*,' Not set up for auto-input with coulomb '
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr) 
        stop 
     end if 
     read(autoinputfile,*)spscale,ax,bx,x 
     if(.not.iso_int)read(autoinputfile,*)sc2 
     if ( bx == 0. .or. x == 0.0 ) then 
           vscale = ax 
     else 
           vscale = (ax/bx)**x 
     end if 
     vscale = vscale*sc2 
     if(iproc==0)then
		 print*,' Scaling ',spscale,ax,bx,x,vscale      
         write(logfile,*)   ' Scaling: spe ',spscale,' TBMEs by (',ax,'/',bx,')^',x,' = ',vscale      
	 end if
  else 
     if(.not. coul_int)then 
        print*,' Enter scaling for spes, A,B,X ( (A/B)^X ) for TBMEs ',iproc 
        print*,' (If B or X = 0, then scale by A ) ' 
        read*,spscale,ax,bx,x 
        if(.not.iso_int)then
			print*,' Enter overall scaling (enter 1.0 if none)'
			read*,sc2
		end if 
        write(autoinputfile,'(4f10.4,"        ! scaling: xspe, A0,A,x ")' )spscale,ax,bx,x 
        if(.not.iso_int)write(autoinputfile,*)sc2 
        if ( bx == 0. .or. x == 0.0 ) then 
           vscale = ax 
        else 
           vscale = (ax/bx)**x 
        end if 
        vscale = vscale*sc2 
        if(iproc==0)print*,' Scaling ',spscale,ax,bx,x,vscale 
        write(logfile,*)   ' Scaling: spe ',spscale,' TBMEs by (',ax,'/',bx,')^',x,' = ',vscale      
 
     else 
          write(6,*)'Enter overall strength and hbw0 and hbw - scaling = sqrt(hbw/hbw0)' 
          read(5,*) sc2, hbw0, hbw 
          vscale = sc2*sqrt(hbw/hbw0) 
          write(6,*)'Enter scaling ratios for each orbit' 
          allocate(ratio(numorb(1)), stat=aerr) 
          if(aerr /= 0) call memerror("readv2bme_oxbash_iso 1") 
          read(5,*)(ratio(i), i = 1, numorb(1)) 
     end if 
  endif 
end if
!-------------- READ IN SPEs --------------
  print*,' * * '
  print*,' * * NOTICE: I expect single-particles space with ',numorb(1),' orbits '
  print*,' * * ' 
  if(autoscale)then
	  spscale = 1.
	  vscale = 1.
	  numval = numorb(1)+3
  else
	  numval = numorb(1)
  end if
	  
   allocate(spetmp(numval), stat=aerr) 

  if(aerr /= 0) call memerror("readv2bme_oxbash_iso 2") 
  if(.not. mixed_T)then 
	  if(autoscale)then
         read(1,*,err=1899,end=1899)nme,(spetmp(i),i= 1,MIN(13,numval)) 
	 else
         read(1,*,err=1899,end=1899)nme,(spetmp(i),i= 1,MIN(10,numval)) 
		 
	 end if
  else 
     read(1,*)nme 
  end if 
  if(nme < 0)nme = -nme
 
  if(numorb(1) > 10 .and. .not. autoscale)then 
     do m = 10,numval-1,10 
        mmax = MIN(10+m,numval) 
        read(1,*,err=1899,end=1899)(spetmp(i),i=1+m,mmax) 
     end do 
  endif 
  if(pp_int .or. iso_int)pspe(:) = pspe(:) + spscale*spetmp(:) 
  if(nn_int .or. iso_int)nspe(:) = nspe(:) + spscale*spetmp(:) 
 
  if(isov_int)then 
     pspe(:) = pspe(:) + spscale*spetmp(:) 
     nspe(:) = nspe(:) - spscale*spetmp(:) 
  end if 
 
  if(coul_int)then 
     pspe(:) = pspe(:) + vscale*spetmp(:)*ratio(:) 
  end if     
 
  if(autoscale)then
	  bx = spetmp(numorb(1)+1)+ npeff(1)+npeff(2)
	  ax = spetmp(numorb(1)+2)
	  x = spetmp(numorb(1)+3)
      vscale = (ax/bx)**x 
      if(iproc==0)then
 		 print*,' Autoscaling ',ax,bx,x,vscale     
		 print*,' I think I have a core of ', spetmp(numorb(1)+1)
          write(logfile,*)   ' Autoscaling  TBMEs by (',ax,'/',bx,')^',x,' = ',vscale      
 	 end if	  
  end if
  deallocate( spetmp ) 
 
!----------------------------------- ERROR TRAP/REMINDER------ ADDED in 7.5.1-----

if(pn_int .and. pnsymorb)then
	print*,'         + X + X + X + X + X + X + X + '
	print*,'  Attention: you have flag "pnsymorb" in module flagger set TRUE '
	print*,'  This assumes, for proton-neutron part of interaction, '
	print*,'  the proton and neutron orbits are the same.'
	print*,'  This means, as an example, that V_pn(ab,cd) = V_pn(ba,dc) (up to a phase) '
	print*,'  WARNING this can be tricky, so this setting is not recommended '
	print*,'  However this works if you copy over an isospin conserving force '
	print*,' and relabel it as XXXX.pn.int. Use at your risk. '
	print*,'         + X + X + X + X + X + X + X + '
end if
if(pn_int .and. .not. pnsymorb)then
	print*,'         + X + X + X + X + X + X + X + '
	print*,'  Attention: you have flag "pnsymorb" in module flagger set FALSE '
	print*,'  This does NOT assume, for proton-neutron part of interaction, '
	print*,'  the proton and neutron orbits are the same.'
	print*,'  This means, as an example, you need to supply '
	print*,'  V_pn(ab,cd) and V_pn(ba,dc) separately.'
	print*,'         + X + X + X + X + X + X + X + '
end if

!-------------- READ IN TBMEs ---------------- 
  write(logfile,*)nme,' two-body matrix elements in this file '
  maxsp = 0 
  do i = 1,nme 
!---------- ERROR TRAP ADDED July 2011 CWJ ------------------' 
     
     if(.not. mixed_T)then 
        read(1,*,err=1899,end=1899)ia,ib,ic,id,j,t,vv 
 
!---------- ERROR TRAP ADDED July 2011 CWJ ------------------' 
 
        if( t /= 0 .and. t /= 1)then 
            print*,' bad TBME matrix element; T value is bad ',t 
            print*,' matrix element # ',i 
            print*,ia,ib,ic,id,j,t,vv 
            print*,' Check that single particle space matches hamiltonian ' 
            print*,' # of single-particle orbits = ',numorb(1) 
			print*,' STOPPING RUN '
		    call BMPI_Abort(icomm,101,ierr)
            stop 
         end if 
     else 
        read(1,*,err=1899,end=1899)ia,ib,ic,id,j,t_ab,t_cd,vv 
!---------- ERROR TRAP ADDED July 2011 CWJ ------------------' 
 
        if(.not. ( (t_ab == 0 .and. t_cd == 1) .or. (t_ab == 1 .and. t_cd == 0)) ) then 
            print*,' bad TBME matrix element; T values are bad ',t_ab, t_cd 
            print*,' matrix element # ',i 
            print*,ia,ib,ic,id,j,t_ab,t_cd,vv 
            print*,' Check that single particle space matches hamiltonian ' 
            print*,' # of single-particle orbits = ',numorb(1) 
			print*,' STOPPING RUN '
		    call BMPI_Abort(icomm,101,ierr)
            stop 
        end if 
     end if 
     maxsp = MAX(maxsp,ia)
     maxsp = MAX(maxsp,ib)
     maxsp = MAX(maxsp,ic)
     maxsp = MAX(maxsp,id)
 
!----------- PUT INTO CORRECT ORDER; PICK UP PHASES ------------- 
!           "CORRECT" ORDER: a >= b, c >= d 
 
!-------------------- PP interaction ------------------------- 
     if(T ==0 )goto 1002 
     if(nn_int)goto 1001           !  NN-interaction only 
     if(pn_int)goto 1002           !  PN-interaction only 
     if(mixed_T) goto 1002 
     if( np(1) < 2 .and. .not.phconj(1))goto 1001       ! order of this statement is important - 7/2011 CWJ 
     if(pp_int .or. iso_int) factor_v = 1.0 
     if(isov_int) factor_v = 0.5 
     if(isot_int) factor_v = 1./6. 
 
     if(coul_int) factor_v = ratio(ia)**0.25*ratio(ib)**0.25*                & 
                             ratio(ic)**0.25*ratio(id)**0.25 
 
     it = 1 
 
     a = ia 
     b = ib 
     c = ic  
     d = id 
!----------- CHECK PARITY ------------------- 
     L = orbqn(it,a)%l+ orbqn(it,b)%l+orbqn(it,c)%l+orbqn(it,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' Oops! Error in parity in (pp) matrix element # ',i
		print*,' Most likely this is a mismatch between single-particle orbits ' 
		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
        write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
		     orbqn(it,a)%j,orbqn(it,a)%l,  orbqn(it,b)%j,orbqn(it,b)%l,  orbqn(it,c)%j,orbqn(it,c)%l,  orbqn(it,d)%j,orbqn(it,d)%l
			 
	    write(logfile,*)' Oops! Error in parity in (pp) matrix element # ',i
	 	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
	    write(logfile,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
	 		     orbqn(it,a)%j,orbqn(it,a)%l,  orbqn(it,b)%j,orbqn(it,b)%l,  orbqn(it,c)%j,orbqn(it,c)%l,  orbqn(it,d)%j,orbqn(it,d)%l	
		close(logfile)		 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
     phase = 1 
 
 
     if(a < b)then 
        na = b 
        b = a 
        a = na 
        phase = (-1)**( J+T+(orbqn(it,a)%j+orbqn(it,b)%j)/2)  ! check 
     endif 
 
     if(c < d)then 
        nc = d 
        d = c 
        c = nc 
        phase = phase*(-1)**( J+T+(orbqn(it,c)%j+orbqn(it,d)%j)/2)  ! check 
     endif 
 
     if(a < c .or. (a==c .and. b < d))then 
        na = c 
        nb = d 
        c = a 
        d = b 
        a = na 
        b = nb 
     endif 
 
!---------- CONVERT ------------------------- 
 
     pair1 = a*(a-1)/2 + b 
     pair2 = c*(c-1)/2 + d 
     if( PPcouplemap(pair1) == -1) goto 1001 
     if( PPcouplemap(pair2) == -1) goto 1001 
     pair1 = PPcouplemap(pair1) 
     pair2 = PPcouplemap(pair2) 
 
     ipar = XXcouples(it)%pairc(pair1)%par 
     if(ipar /= XXcouples(it)%pairc(pair2)%par)then 
        print*,' problem with parity, boss  (PP) ',a,b,c,d 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     iref = XXcouples(it)%meref(ipar) 
     istart = XXcouples(it)%mestart(ipar) 
     if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif 
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!--------------- ERROR TRAP --------------- 
 
     if(j > ppme(indx)%jmax .or. j < ppme(indx)%jmin)then 
		 
        print*,' Oops! Error in adding Js in (pp) matrix element # ',i
 		print*,' Most likely this is a mismatch between single-particle orbits ' 
 		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
        write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
 		     orbqn(it,a)%j, orbqn(it,b)%j,  orbqn(it,c)%j, orbqn(it,d)%j
		write(6,*)' See log file for additional information '	 
	    write(logfile,*)' Oops! Error in adding Js in (pp) matrix element # ',i
	  	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
	    write(logfile,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
	  		     orbqn(it,a)%j, orbqn(it,b)%j,  orbqn(it,c)%j, orbqn(it,d)%j

		write(logfile,*)' Additional information ' 
        write(logfile,*)' J = ',j,' but min/max for these pairs are ',ppme(indx)%jmin,ppme(indx)%jmax 
		
!        print*,' error in Js (pp) ',pair1,pair2,indx 
!        print*,XXcouples(it)%pairc(pair1)%ia,XXcouples(it)%pairc(pair1)%ib 
!        print*,' orbits ',a,b,c,d,' J T ',j,t,' V ',vv
!        print*,a*(a-1)/2 + b, c*(c-1)/2 + d 
!        print*,PPcouplemap(a*(a-1)/2 + b),PPcouplemap( c*(c-1)/2 + d) 
!        print*,orbqn(it,a)%j,orbqn(it,b)%j,orbqn(it,c)%j,orbqn(it,d)%j 
 		close(logfile)		
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
         
        stop 
     endif 
     ppme(indx)%v(j)=ppme(indx)%v(j)+vv*vscale*phase*factor_v  
1001 continue 
     if(pp_int .or. coul_int)cycle 
 
!.......................................................................	 
!---------------------------------- NN INTERACTION ------------- 
!.......................................................................	 
	 

     if(nn_int .or. iso_int) factor_v = 1.0 
     if(np(2) < 2 .and. .not.phconj(2))goto 1002 
     if(isov_int) factor_v = -0.5 
     if(isot_int) factor_v = 1./6. 
 
     it = 2 
     a = ia 
     b = ib 
     c = ic  
     d = id 
!----------- CHECK PARITY ------------------- 
 
     L = orbqn(it,a)%l+ orbqn(it,b)%l+orbqn(it,c)%l+orbqn(it,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' Oops! Error in parity in matrix element # ',i
 		print*,' Most likely this is a mismatch between single-particle orbits ' 
 		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
         write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
 		     orbqn(it,a)%j,orbqn(it,a)%l,  orbqn(it,b)%j,orbqn(it,b)%l,  orbqn(it,c)%j,orbqn(it,c)%l,  orbqn(it,d)%j,orbqn(it,d)%l
			 
 	    write(logfile,*)' Oops! Error in parity in matrix element # ',i
 	 	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
 	 		     orbqn(it,a)%j,orbqn(it,a)%l,  orbqn(it,b)%j,orbqn(it,b)%l,  orbqn(it,c)%j,orbqn(it,c)%l,  orbqn(it,d)%j,orbqn(it,d)%l	
 		close(logfile)
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
     phase = 1 
      
     if(a < b)then 
        na = b 
        b = a 
        a = na 
        phase = (-1)**( J+T+(orbqn(it,a)%j+orbqn(it,b)%j)/2)  ! check 
     endif 
 
     if(c < d)then 
        nc = d 
        d = c 
        c = nc 
        phase = phase*(-1)**( J+T+(orbqn(it,c)%j+orbqn(it,d)%j)/2)  ! check 
     endif 
 
     if(a < c .or. (a==c .and. b < d))then 
        na = c 
        nb = d 
        c = a 
        d = b 
        a = na 
        b = nb 
     endif 
 
     pair1 = a*(a-1)/2 + b 
     pair2 = c*(c-1)/2 + d 
     if( NNcouplemap(pair1) == -1) goto 1002 
     if( NNcouplemap(pair2) == -1) goto 1002 
     pair1 = NNcouplemap(pair1) 
     pair2 = NNcouplemap(pair2) 
 
     ipar = XXcouples(it)%pairc(pair1)%par 
     if(ipar /= XXcouples(it)%pairc(pair2)%par)then 
        print*,' problem with parity, boss (NN oxb) ',a,b,c,d
	print*,ia,ib,ic,id,i
	print*,' STOPPING RUN '
    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     iref = XXcouples(it)%meref(ipar) 
     istart = XXcouples(it)%mestart(ipar) 
      if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif  
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!--------------- ERROR TRAP --------------- 
 
     if(j > nnme(indx)%jmax .or. j < nnme(indx)%jmin)then  
        print*,' Oops! Error in adding Js in (nn) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
        write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
  		     orbqn(it,a)%j, orbqn(it,b)%j,  orbqn(it,c)%j, orbqn(it,d)%j
 		write(6,*)' See log file for additional information '	 
 	    write(logfile,*)' Oops! Error in adding Js in (pp) matrix element # ',i
 	  	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
 	  		     orbqn(it,a)%j, orbqn(it,b)%j,  orbqn(it,c)%j, orbqn(it,d)%j

 		write(logfile,*)' Additional information ' 
         write(logfile,*)' J = ',j,' but min/max for these pairs are ',nnme(indx)%jmin,nnme(indx)%jmax 
		 
!        print*,' error in Js (nn) ',pair1,pair2,indx 
!        print*,orbqn(it,a)%j,orbqn(it,b)%j,orbqn(it,c)%j,orbqn(it,d)%j 
!        print*,nnme(indx)%jmax,nnme(indx)%jmin 
        close(logfile)
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     nnme(indx)%v(j)=nnme(indx)%v(j)+vv*vscale*phase*factor_v 
 
!............................................................... 
!     if ( nproc > 1 ) then 
!        n_nn = n_nn + 1 
!        nn_indx_J(1,n_nn) = indx 
!        nn_indx_J(2,n_nn) = j 
!        nn_bcast(n_nn) = nnme(indx)%v(j) 
!     end if 
!............................................................... 
1002 continue 
     if(nn_int .or. isov_int)cycle 
!.......................................................................	 
!--- PN --- PN --- PN  --- PN --- PN ---- PN ----- PN ---- PN ------
!.......................................................................	 
!
!  PN is more complicated, especially if one is unfolding from isospin formalism
!  or if the flag pnsymorb (found in module flagger) which signals 
!  proton and neutron orbits are identity
!
!  | a_p b_n > = sqrt((1+delta_ab))/sqrt(2) ( |a b, T=1> + |a,b, T = 0> )
!
!  V(aa,bb) store just this
!
!  V(ab,cc) or V(aa,bc)  : if isospin OR if pnsymorb  store also V(ba,cc)
!
!  V(ab,ab)              : if isospin OR if pnsymorb  store also V(ba,ba)
!
!  V(ab,ba)              : if isospin (unlikely) store also V(ab,ab)
!                        
!  V(ab,cd)              : if isospin store also V(ba,cd), V(ab,dc)
!                        : if isospin or pnsymorb store also V(ba,dc)
!
!.......................................................................	 
!--- PN --- PN --- PN  --- PN --- PN ---- PN ----- PN ---- PN ------
!.......................................................................	 

     if(pp_int .or. nn_int .or. isov_int .or. coul_int)cycle 
     if( np(1)*np(2) ==0 .and. .not.phconj(1) .and. .not.phconj(2))cycle 
     if(pn_int .or. iso_int .or. mixed_T) factor_v = 1.0 
!     if(pn_int .or. iso_int) factor_v = 1.0 
     if(isot_int) then 
         if(T == 1)then 
            factor_v = -1./3. 
         else 
            cycle 
         end if 
     end if 
 
!...........PN PART 1 ......................................
     a = ia 
     b = ib 
     c = ic  
     d = id 
     ! factor due to identical label 
     factor = 0.5 
     if(a == b)factor = factor*sqrt(2.) 
     if(c == d) factor = factor*sqrt(2.) 
 
     if(a < c .or. (a == c .and. b < d))then  ! swap 
        na = a 
        a = c 
        c = na 
        na = b 
        b = d 
        d = na 
     endif 
!----------- CHECK PARITY ------------------- 
 
     L = orbqn(1,a)%l+ orbqn(2,b)%l+orbqn(1,c)%l+orbqn(2,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' Oops! Error in parity in (pn 1) matrix element # ',i
 		print*,' Most likely this is a mismatch between single-particle orbits ' 
 		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
         write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
 		     orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
			 
 	    write(logfile,*)' Oops! Error in parity in matrix element # ',i
 	 	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
           orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
 		close(logfile)		
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
!-------- NO PHASE, ALREADY IN CORRECT ORDER 
     phase = 1 
 
     pair1 = numorb(2)*(a-1) + b 
     pair2 = numorb(2)*(c-1) + d 
     if( PNcouplemap(pair1) == -1) goto 1003 
     if( PNcouplemap(pair2) == -1) goto 1003 
     pair1 = PNcouplemap(pair1) 
     pair2 = PNcouplemap(pair2) 
 
     ipar = PNcouples%pairc(pair1)%par 
     if(ipar /= PNcouples%pairc(pair2)%par)then 
        print*,' problem with parity, boss (pn 1)  '
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif 
     iref = PNcouples%meref(ipar) 
     istart = PNcouples%mestart(ipar) 
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
 
!--------------- ERROR TRAP --------------- 
 
     if(j > pnme(indx)%jmax .or. j < pnme(indx)%jmin)then 
		 
         print*,' Oops! Error in adding Js in (pn1) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
         write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
    	     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j
 		write(6,*)' See log file for additional information '	 
 	    write(logfile,*)' Oops! Error in adding Js in (pp) matrix element # ',i
 	  	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
 	  		     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j

 		 write(logfile,*)' Additional information ' 
         write(logfile,*)' J = ',j,' but min/max for these pairs are ',pnme(indx)%jmin,pnme(indx)%jmax 		
		 close(logfile) 
!        print*,i,' error in Js (pn 1) ',pair1,pair2,indx 
!        print*,a,b,c,d,j,t,vv 
!        print*,orbqn(1,a)%j,orbqn(2,b)%j,orbqn(1,c)%j,orbqn(2,d)%j 
!        print*,j,pnme(indx)%jmin,pnme(indx)%jmax 
!        print*,PNcouples%pairc(pair1)%ia,PNcouples%pairc(pair1)%ib 
!        print*,PNcouples%pairc(pair2)%ia,PNcouples%pairc(pair2)%ib 
         print*,' STOPPING RUN '
        call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
!     pnme(indx)%v(j)=pnme(indx)%v(j)+vv*vscale*phase*factor*factor_v 
     if(.not. mixed_T)then 
        pnme(indx)%v(j)=pnme(indx)%v(j)+vv*vscale*phase*factor*factor_v 
     else 
        pnme(indx)%v(j)=pnme(indx)%v(j) - vv*vscale*phase*factor*factor_v 
     end if 

1003 continue 
     if( pn_int .and. .not.pnsymorb) cycle      ! none of them will work

 
!------------ second pn CASE ---swap a and b 
     b = ia 
     a = ib 
     c = ic  
     d = id 
!-------------- don't bother if they work out the same 
     if(a==b)goto 1004 
	 
!-------- CHECK CASE -----------

!	 if( pn_int .and. a==d .and. b==c .and. a/=b)then
!		 print*,a,b,c,d,j,t,vv
!		  goto 1004   !unlikely
!	  end if
     if(a < c .or. (a == c .and. b < d))then  ! swap 
        na = a 
        a = c 
        c = na 
        na = b 
        b = d 
        d = na 
     endif 
!----------- CHECK PARITY ------------------- 
 
     L = orbqn(2,a)%l+ orbqn(1,b)%l+orbqn(1,c)%l+orbqn(2,d)%l 
 
     if( (-1)**(L) ==-1)then 
         print*,' Oops! Error in parity in (pn 2) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
          write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
  		     orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
			 
  	    write(logfile,*)' Oops! Error in parity in (pn 2) matrix element # ',i
  	 	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
  	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
  	    write(logfile,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
            orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
  		close(logfile)		
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
!     phase = (-1)**( J+T+(orbqn(2,ia)%j+orbqn(1,ib)%j)/2)  ! check 
     if(.not. mixed_T)then 
        phase = (-1)**( J+T+(orbqn(2,ia)%j+orbqn(1,ib)%j)/2)  ! check 
     else 
        phase = (-1)**( J + t_ab + (orbqn(2,ia)%j+orbqn(1,ib)%j)/2)  ! check 
     end if 
 
     pair1 = numorb(2)*(a-1) + b 
     pair2 = numorb(2)*(c-1) + d 
     if( PNcouplemap(pair1) == -1) goto 1004 
     if( PNcouplemap(pair2) == -1) goto 1004 
     pair1 = PNcouplemap(pair1) 
     pair2 = PNcouplemap(pair2) 
 
     ipar = PNcouples%pairc(pair1)%par 
     if(ipar /= PNcouples%pairc(pair2)%par)then 
        print*,' problem with parity, boss (pn 2) ' 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif 
     iref = PNcouples%meref(ipar) 
     istart = PNcouples%mestart(ipar) 
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!--------------- ERROR TRAP --------------- 
      
     if(j > pnme(indx)%jmax .or. j < pnme(indx)%jmin)then 
         print*,' Oops! Error in adding Js in (pn2) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
         write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
    	     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j
 		write(6,*)' See log file for additional information '	 
 	    write(logfile,*)' Oops! Error in adding Js in (pn2) matrix element # ',i
 	  	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
 	  		     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j

 		 write(logfile,*)' Additional information ' 
         write(logfile,*)' J = ',j,' but min/max for these pairs are ',pnme(indx)%jmin,pnme(indx)%jmax 		
		 close(logfile) 
!        print*,' error in Js (pn 2) ',pair1,pair2,indx 
!        print*, a, b, c, d 
!        print*,orbqn(2,a)%j,orbqn(1,b)%j,orbqn(1,c)%j,orbqn(2,d)%j 
!        print*,pnme(indx)%jmax,pnme(indx)%jmin 
        print*,' STOPPING RUN '
        call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     pnme(indx)%v(j)=pnme(indx)%v(j)+vv*vscale*phase*factor*factor_v 

1004 continue 
 
!------------ third pn CASE ---swap c and d 
     a = ia 
     b = ib 
     d = ic  
     c = id 
!-------------- don't bother if they work out the same 
     if(c == d) goto 1005 
     if(c == b .and. d == a )goto 1005 
      
 
     if(a < c .or. (a == c .and. b < d))then  ! swap 
        na = a 
        a = c 
        c = na 
        na = b 
        b = d 
        d = na 
     endif 
!----------- CHECK PARITY ------------------- 
 
     L = orbqn(1,a)%l+ orbqn(2,b)%l+orbqn(2,c)%l+orbqn(1,d)%l 
      
     if( (-1)**(L) ==-1)then 
         print*,' Oops! Error in parity in (pn 3) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
          write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
  		     orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
			 
  	    write(logfile,*)' Oops! Error in parity in (pn 3) matrix element # ',i
  	 	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
  	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
  	    write(logfile,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
            orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
  		close(logfile)	
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
!     phase = (-1)**( J+T+(orbqn(2,ic)%j+orbqn(1,id)%j)/2)  ! check 
     if(.not. mixed_T)then 
     phase = (-1)**( J+T+(orbqn(2,ic)%j+orbqn(1,id)%j)/2)  ! check 
     else 
        phase = (-1)**( J + t_cd + (orbqn(2,ic)%j+orbqn(1,id)%j)/2)  ! check 
     end if 
 
     pair1 = numorb(2)*(a-1) + b 
     pair2 = numorb(2)*(c-1) + d 
     if( PNcouplemap(pair1) == -1) goto 1005 
     if( PNcouplemap(pair2) == -1) goto 1005 
     pair1 = PNcouplemap(pair1) 
     pair2 = PNcouplemap(pair2) 
     if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif 
     ipar = PNcouples%pairc(pair1)%par 
     if(ipar /= PNcouples%pairc(pair2)%par)then 
        print*,' problem with parity, boss (pn 3) ' 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     iref = PNcouples%meref(ipar) 
     istart = PNcouples%mestart(ipar) 
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!--------------- ERROR TRAP --------------- 
 
     if(j > pnme(indx)%jmax .or. j < pnme(indx)%jmin)then 
         print*,' Oops! Error in adding Js in (pn 3) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
         write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
    	     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j
 		write(6,*)' See log file for additional information '	 
 	    write(logfile,*)' Oops! Error in adding Js in (pn 3) matrix element # ',i
 	  	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
 	  		     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j

 		 write(logfile,*)' Additional information ' 
         write(logfile,*)' J = ',j,' but min/max for these pairs are ',pnme(indx)%jmin,pnme(indx)%jmax 		
		 close(logfile) 
		 
!        print*,' error in Js (pn 3) ',pair1,pair2,indx 
!        print*,orbqn(it,a)%j,orbqn(it,b)%j,orbqn(it,c)%j,orbqn(it,d)%j 
!        print*,j,pnme(indx)%jmin,pnme(indx)%jmax 
        print*,' STOPPING RUN '
        call BMPI_Abort(icomm,101,ierr)
             
        stop 
     endif 
     pnme(indx)%v(j)=pnme(indx)%v(j)+vv*vscale*phase*factor*factor_v 

1005 continue 
 
!------------ fourth pn CASE ---swap a and b and c and d 
     b = ia 
     a = ib 
     d = ic  
     c = id 
 
!-------------- don't bother if they work out the same 
     if(a == b)goto 1006 
     if(c==d)goto 1006 
      
     if(a < c .or. (a == c .and. b < d))then  ! swap 
        na = a 
        a = c 
        c = na 
        na = b 
        b = d 
        d = na 
     endif 
!----------- CHECK PARITY ------------------- 
 
     L = orbqn(2,a)%l+ orbqn(1,b)%l+orbqn(2,c)%l+orbqn(1,d)%l 
      
     if( (-1)**(L) ==-1)then 
         print*,' Oops! Error in parity in (pn 4) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
          write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
  		     orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
			 
  	    write(logfile,*)' Oops! Error in parity in matrix element # ',i
  	 	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
  	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
  	    write(logfile,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
            orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
  		close(logfile)		
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
     phase = (-1)**( (orbqn(2,a)%j+orbqn(1,b)%j+ orbqn(2,c)%j+orbqn(1,d)%j)/2)  ! check 
 
     pair1 = numorb(2)*(a-1) + b 
     pair2 = numorb(2)*(c-1) + d 
     if( PNcouplemap(pair1) == -1) goto 1006 
     if( PNcouplemap(pair2) == -1) goto 1006 
     pair1 = PNcouplemap(pair1) 
     pair2 = PNcouplemap(pair2) 
 
     ipar = PNcouples%pairc(pair1)%par 
     if(ipar /= PNcouples%pairc(pair2)%par)then 
        print*,' problem with parity, boss (pn 4) ' 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     iref = PNcouples%meref(ipar) 
     istart = PNcouples%mestart(ipar) 
     if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
 
!--------------- ERROR TRAP --------------- 
 
     if(j > pnme(indx)%jmax .or. j < pnme(indx)%jmin)then 
         print*,' Oops! Error in adding Js in (pn 4) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
         write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
    	     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j
 		write(6,*)' See log file for additional information '	 
 	    write(logfile,*)' Oops! Error in adding Js in (pp) matrix element # ',i
 	  	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
 	  		     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j

 		 write(logfile,*)' Additional information ' 
         write(logfile,*)' J = ',j,' but min/max for these pairs are ',pnme(indx)%jmin,pnme(indx)%jmax 		
		 close(logfile) 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
           
        stop 
     endif 
     pnme(indx)%v(j)=pnme(indx)%v(j)+vv*vscale*phase*factor*factor_v 

1006 continue 
 
  enddo  !nme 
  if(allocated(ratio))deallocate(ratio) 
!........... WARNING MESSAGES.... 
  if(maxsp < maxorblabel)then
     if(iproc==0)then
      print*,' '
      print*,' WARNING '
      print*,' This is missing matrix elements within the active single-particle space '
      print*,' May cause problems in runs ! '
      print*,' '
	  
      write(logfile,*)' WARNING '
      write(logfile,*)' This is missing matrix elements within the active single-particle space '
      write(logfile,*)' May cause problems in runs ! '
      end if
  end if
  if(maxsp > MAX(numorb(1),numorb(2)) )then
     if(iproc==0)then
      print*,' '
      print*,' WARNING '
      print*,' This file contains matrix elements beyond the active single-particle space '
      print*,' i.e., missing matrix elements coupling single particle orbits '
      print*,' from ',MAX(numorb(1),numorb(2)), ' to ',maxorblabel
      print*,' May cause problems in runs ! '
      print*,' END WARNING '
      print*,' '
	  
      write(logfile,*)' WARNING '
      write(logfile,*)' This file contains matrix elements beyond the active single-particle space '
      write(logfile,*)' i.e., missing matrix elements coupling single particle orbits '
      write(logfile,*)' from ',MAX(numorb(1),numorb(2)), ' to ',maxorblabel
      write(logfile,*)' May cause problems in runs ! '
      write(logfile,*)' END WARNING '
      end if
  end if
  firstread= .false.
!......... ADDED IN 7.8.2.... GO ON TO LOOK FOR OFF-DIAGONAL SINGLE-PARTICLE  
  formattedfile = .true.
  call readinsppot('iso',1,spscale,spscale,formattedfile,foundit)
  if(foundit)then
	  print*,' Finished reading in off-diagonal one-body potential '
	  write(logfile,*)' Finished reading in off-diagonal one-body potential '
  else
	  print*,' No one-body potential beyond single-particle energies '
	  write(logfile,*)' No one-body potential beyond single-particle energies '

  end if
  close(1)
  
  
  return

1899 continue
  if(iproc==0)then
    print*,' Error in reading file; problem may be incommensurate size of single-particle space '
    print*,' I expect single-particles space with ',numorb(1),' orbits '
    print*,nme,i
    backspace(1)
    read(1,'(a)')title
    write(6,*)title
  end if
    call BMPI_Abort(icomm,101,ierr)
   
end subroutine readv2bme_oxbash_iso 
 
!===================================================================== 
! 
!    a b c d J T  Trel Hcm  Vcoul  Vpn  Vpp Vnn 
! 
!    added May 2009 by CWJ 
!    modified July 2010 by CWJ to improve storage of coupled TBMEs 
! 
!   THE TBMEs V_JT(ab,cd) are ordered in a particular way: 
!   a >= b, c >= d;  a >= c 
!     VTBME(indx): 
!      indx =  pair1*(pair1-1)/2+pair2 
!    WHERE     pair1 = ia*(ia-1)/2 + ib 
!              pair2 = ic*(ic-1)/2 + id 
! 
!===================================================================== 
subroutine readv2bme_mfd(firstread,formattedfile,intcase) 
   
  use io 
  use system_parameters 
  use sporbit 
  use W_info 
  use coupledmatrixelements 
  use nodeinfo 
  use bmpi_mod 
  implicit none 
  
  logical :: firstread,formattedfile 
  integer :: intcase 
   
  integer :: a,b,c,d 
  integer :: ia,ib,ic,id 
  integer :: na,nb,nc,nd 
  integer :: pair1,pair2,indx 
  integer :: j,t 
  real :: Vv 
  real, allocatable :: spetmp(:) 
  real :: spscale,ax,bx,x,vscale  ! for scaling interactions 
  integer :: phase 
  integer :: nme 
  real :: factor 
 
  integer :: pcpar, pcref,pcstart 
  integer :: dw 
  integer :: it 
  integer :: ierr
   
  real :: Trel, Hrel,Hcm, Vcoul, Vpn, Vpp, Vnn, Vx 
  integer  :: maxsp   ! largest s.p. index, used for cross-checking
  character(80) :: dummyline
 
!--------------------------------------------------------------------- 
!      dummy counters 
!--------------------------------------------------------------------- 
  integer :: i,m,mmax 
     
  integer :: L 
 
  character :: bang 
  character(40) :: title
  character :: breakchar 
  logical :: finished ,foundit
  character(3) :: formattest
 
  real :: betacm,cm1,hcmspe 
  real :: hw_in_file      ! parameter found in file 
  integer ::  Nprinc_in_file,N2_in_file       ! parameters found in file 
  integer :: Nprinc_def   ! original definition
  integer ::  maxorbit
 
!------------- READ PAST HEADER--------------------------------------- 
 
  if ( formattedfile ) then 
     finished = .false. 
     do while(.not.finished) 
        read(1,'(a)')bang 
        if(bang/='!' .and. bang/='#')then 
           backspace(1) 
           finished = .true. 
        endif 
		if(.not.finished)then
			read(1,'(a)')title
			call fileformatcheck(title,foundit,formattest)
			if(foundit)then ! check format
		   	 	if(formattest /= 'mfd')then
	 			   print*,' '
	 			   print*,'  WARNING !   '
	 			  print*,' File appears to be in wrong format '
	 			  print*,' Header denotes format as ',formattest
	 		      print*,'  SKIPPING FILE !   '
	   		      print*,' '
				  write(logfile,*)' File appears to be in wrong format '
				  write(logfile,*)' Header denotes format as ',formattest
			      write(logfile,*)'  SKIPPING FILE !   '
	 			  return
			  	  return
		  		end if
	    	end if
		end if
     enddo 
  endif 
    
!--------------- GET # OF MATRIX ELEMENTS ------------------------ 
!                ALSO CHECK IF FORMATTED
!              IF POSSIBLE EXTRACT OTHER PARAMETERS (added 7.4.1)
 
 
  if ( formattedfile ) then 
!            MODiFIED IN 7.6.6 as some files do not include Nprinc_in_file, etc in first line	  
	  read(1,'(a)')dummyline        ! ADDED IN 7.6.6 to avoid errors in reading first line
	  print*,' First line of interaction file = ',dummyline
	  rewind(1)
     read(1,*,err=331) nme  !,Nprinc_in_file,N2_in_file,hw_in_file
!........... NOW SEARCH THROUGH FILE......................

     maxorbit = 0
     N2_in_file = 0
	 Nprinc_in_file = 0
!     read(1,*) nme 
     do i = 1,nme
       read(1,*)na,nb,nc,nd
	   maxorbit=max(maxorbit,na)
	   maxorbit=max(maxorbit,nb)
	   maxorbit=max(maxorbit,nc)
	   maxorbit=max(maxorbit,nd)
	   na = int((sqrt(8.*na-7.) -1.)/2.)
	   nb = int((sqrt(8.*nb-7.) -1.)/2.)
	   nc = int((sqrt(8.*nc-7.) -1.)/2.)
	   nd = int((sqrt(8.*nd-7.) -1.)/2.)
	   N2_in_file = max(N2_in_file,na+nb)
	   N2_in_file = max(N2_in_file,nc+nd)
	   Nprinc_in_file = max(Nprinc_in_file,na)
	   Nprinc_in_file = max(Nprinc_in_file,nb)
	   Nprinc_in_file = max(Nprinc_in_file,nc)
	   Nprinc_in_file = max(Nprinc_in_file,nd)
	   
     end do
	 
!........ NOW FIND MAX NPRINC IN ORIGINAL DEFINTION OF ORBITS....
     Nprinc_def = 0
     do i = 1,numorb(2)
		 Nprinc_def = max(Nprinc_def,orbqn(2,i)%l+2*orbqn(2,i)%nr)		 
	 end do
	 
	 if(Nprinc_def < Nprinc_in_file)then
		 print*,' '
		 print*,' Attention! You defined your single particle space '
		 print*,' to have a max principle quantum number N = ',Nprinc_def
		 print*,' BUT the file has a max N of ',Nprinc_in_file
		 print*,' This could lead to errors in interpreting the file '
		 print*,' You can solve this by increasing the defined max principle quantum number '
		 print*,' '
		 write(logfile,*)' Attention! S.P. space defined by max principal N = ',Nprinc_def
		 write(logfile,*)' But interaction file has max N of ',Nprinc_in_file
		 write(logfile,*)' This could lead to errors in interpreting the file '
		 write(logfile,*)' You can solve this by increasing the defined max principle quantum number '
		 
	 end if
	 
	 rewind(1)
	 read(1,*)nme

!............END OF SEARCH........................	 
     print*,' Highest orbital label in file = ',maxorbit
     print*,' Highest single-partcle shell has principal quantum number N = ',Nprinc_in_file
     print*,' ->',(Nprinc_in_file+1)*(Nprinc_in_file+2)/2,' orbits '
     print*,' You used ',maxorblabel,' orbits or a highest N of ',nint((sqrt(8.*maxorblabel+1.) -3.)/2.)
     write(logfile,*)' Highest orbital label in file = ',maxorbit	 
     write(logfile,*)' Highest single-partcle shell has principal quantum number N = ',Nprinc_in_file
     write(logfile,*)' ->',(Nprinc_in_file+1)*(Nprinc_in_file+2)/2,' orbits '
     write(logfile,*)' You used ',maxorblabel,' orbits or a highest N of ',nint((sqrt(8.*maxorblabel+1.) -3.)/2.)
     if(maxorblabel > (Nprinc_in_file+1)*(Nprinc_in_file+2)/2)then
        print*,' WARNING: This means there will be active single particle orbits with no matrix elements '
        write(logfile,*)' WARNING: This means there will be active single particle orbits with no matrix elements '
     end if
     print*,' Highest Na+Nb (sum of 2 particle Ns) = ',N2_in_file
!     print*, ' File has basis oscillator frequency = ',hw_in_file
     write(logfile,*)' Highest Na+Nb (sum of 2 particle Ns) = ',N2_in_file
!     write(logfile,*)' File has basis oscillator frequency = ',hw_in_file
     go to 333
331  continue  
     backspace (1)
     print*,' Could not extract basis oscillator frequency, max # of orbits '
     write(logfile,*)' Could not extract basis oscillator frequency, max # of orbits ' 
  else 
	  print*,' Be aware that this needs to rewritten to correspond to formatted file input'
     read(1,err=332) nme,Nprinc_in_file,N2_in_file,hw_in_file
     print*,' Highest single-partcle shell has principal quantum number N = ',Nprinc_in_file
     print*,' ->',(Nprinc_in_file+1)*(Nprinc_in_file+2)/2,' orbits '
     print*,' You used ',numorb(1),' orbits or a highest N of ',nint((sqrt(8.*numorb(1)+1.) -3.)/2.)
     if(maxorblabel > (Nprinc_in_file+1)*(Nprinc_in_file+2)/2)then
        print*,' WARNING: This means there will be active single particle orbits with no matrix elements '
     end if
     print*,' Highest Na+Nb (sum of 2 particle Ns) = ',N2_in_file
     print*, ' File has basis oscillator frequency = ',hw_in_file
     goto 333
332  continue
     backspace (1)
     print*,' Could not extract basis oscillator frequency, max # of orbits '
     read(1) nme 
  endif 

333 continue 

!------------------- GET PARAMETERS---------------------------------- 
 
  if(auto_input)then 
     read(autoinputfile,*)hw,betacm 
  else 
     print*,' Enter oscillator frequency (in MeV) and center-of-mass strength ' 
     read*,hw,betacm 
     write(autoinputfile,*)hw,betacm,'   ! hw,  beta for center of mass ' 
  end if 
  write(logfile,*)' hw = ',hw,' MeV; beta_cm = ',betacm
  cm1 = betacm*hw/float(np(1)+np(2)-1) 
!  cm2 = betacm*2*hw/float(np(1)+np(2)) 
 
!--------------- SINGLE-PARTICLE ENERGIES ARE ZERO-------------------- 
!  pspe(:) = 0.0 
!  nspe(:) = 0.0 
!..................................................................... 
  write(logfile,*)nme,' two-body matrix elements in this file '
  maxsp = 0
  do i = 1,nme 
 
     vpp = 0.0 
     vpn = 0.0 
     vnn = 0.0 
     select case(intcase) 
 
     case (1) 
        if(formattedfile)then 
           read(1,*)a,b,c,d,j,t,Trel,Hrel,Vcoul,  Vx 
        else 
           read(1)a,b,c,d,j,t,Trel,Hrel, Vcoul, Vx 
        endif 
        if(t == 1)then 
           vpp = 2*hw*Trel/float(np(1)+np(2) ) + vx - 2*hw*Hrel*betacm/float(np(1)+np(2)) 
           vnn = 2*hw*Trel/float(np(1)+np(2) ) +vx - 2*hw*Hrel*betacm/float(np(1)+np(2)) 
 
        endif 
        vpn = 2*hw*Trel/float(np(1)+np(2) ) +vx - 2*hw*Hrel*betacm / float(np(1)+np(2)) 
     case (2) 
        if(formattedfile)then 
           read(1,*)a,b,c,d,j,t,Trel, Hrel,Vcoul, Vx 
        else 
           read(1)a,b,c,d,j,t,Trel, Hrel,Vcoul, Vx 
        endif 
        Vcoul = sqrt(938.915*hw)/197.327*Vcoul 
        if(t == 1)then 
           vpp =2*hw*Trel/float(np(1)+np(2) )+ vx - 2*hw*Hrel*betacm/float(np(1)+np(2)) + Vcoul 
           vnn = 2*hw*Trel/float(np(1)+np(2) )+vx - 2*hw*Hrel*betacm/float(np(1)+np(2)) 
 
        endif 
        vpn = 2*hw*Trel/float(np(1)+np(2) )+ vx - 2*hw*Hrel*betacm / float(np(1)+np(2)) 
 
     case (3) 
        if(formattedfile)then 
           read(1,*)a,b,c,d,j,t,Trel,Hrel, Vcoul, Vpn, Vpp, Vnn 
        else 
           read(1)a,b,c,d,j,t,Trel,Hrel, Vcoul, Vpn, Vpp, Vnn 
        endif 
        if(t == 1)then 
           vpp = 2*hw*Trel/float(np(1)+np(2) )+ vpp - 2*hw*Hrel*betacm/float(np(1)+np(2)) 
           vnn = 2*hw*Trel/float(np(1)+np(2) )+ vnn - 2*hw*Hrel*betacm/float(np(1)+np(2)) 
 
        endif 
        vpn = 2*hw*Trel/float(np(1)+np(2) ) + vpn - 2*hw*Hrel*betacm / float(np(1)+np(2)) 
 
     case default      ! error trap 
        print*,' should not be here ',intcase 
        stop 
     end select 
	 
	 if(a > maxorblabel .or. b > maxorblabel .or. c > maxorblabel .or. d >maxorblabel)cycle ! this can avoid bad definitions

     maxsp = MAX(maxsp,a)
     maxsp = MAX(maxsp,b)
     maxsp = MAX(maxsp,c)
     maxsp = MAX(maxsp,d)
 
!------ ADDITIONAL CENTER OF MASS OPERATOR TO PUSH UP SPURIOUS STATES- 
 
     if((a == c .and. b == d) .or. (a == d .and. b == c))then 
        hcmspe = cm1*( orbqn(1,a)%w + orbqn(1,b)%w + 3.0 - 3.0/float( np(1)+np(2)) ) 
     else 
        hcmspe = 0.0 
     endif 
 
     vpp = vpp+hcmspe 
     vnn = vnn + hcmspe 
     vpn = vpn + hcmspe 
 
!----------- PUT INTO CORRECT ORDER; PICK UP PHASES ----------------- 
!           "CORRECT" ORDER: a >= b, c >= d 
 
!-------------------- PP interaction -------------------------------- 
     if ( T ==0 ) goto 1002 
     it = 1 
     ia = a 
     ib = b 
     ic = c 
     id = d 
!----------- CHECK PARITY ----------------------------------------- 
 
     L = orbqn(it,a)%l+ orbqn(it,b)%l+orbqn(it,c)%l+orbqn(it,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' error in parity ' 
        print*,a,b,c,d,j,t,vpp 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
     phase = 1 
 
     if(ia < ib)then 
        na = ib 
        ib = ia 
        ia = na 
        phase = (-1)**( J+T+(orbqn(it,a)%j+orbqn(it,b)%j)/2)  ! check 
     endif 
 
     if(ic < id)then 
        nc = id 
        id = ic 
        ic = nc 
        phase = phase*(-1)**( J+T+(orbqn(it,c)%j+orbqn(it,d)%j)/2)  ! check 
     endif 
 
     if(ia < ic .or. (ia==ic .and. ib < id))then 
        na = ic 
        nb = id 
        ic = ia 
        id = ib 
        ia = na 
        ib = nb 
     endif 
 
!---------- CONVERT -------------------------------------------- 
 
     pair1 = ia*(ia-1)/2 + ib 
     pair2 = ic*(ic-1)/2 + id 
     pair1 = PPcouplemap(pair1) 
     pair2 = PPcouplemap(pair2) 
     if(pair1 == -1 .or. pair2 == -1)goto 1001 
     pcpar = XXcouples(it)%pairc(pair1)%par 
     pcref = XXcouples(it)%meref(pcpar) 
     pcstart = XXcouples(it)%mestart(pcpar) 
 
     if(pair1 >= pair2)then 
       indx = (pair1-pcref)*(pair1-pcref-1)/2+pair2-pcref+pcstart 
     else 
       indx = (pair2-pcref)*(pair2-pcref-1)/2+pair1-pcref+pcstart 
     endif 
!--------------- ERROR TRAP --------------------------------------- 
 
     if(j > ppme(indx)%jmax .or. j < ppme(indx)%jmin)then 
        print*,' error in Js (pp) ',pair1,pair2,indx 
        print*,' orbits ',a,b,c,d 
        print*,ia,ib,ic,id,J,T,vpp 
        print*,orbqn(it,a)%j,orbqn(it,b)%j,orbqn(it,c)%j,orbqn(it,d)%j 
        print*,j,ppme(indx)%jmin,ppme(indx)%jmax 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
             
        stop 
     endif 
     ppme(indx)%v(j)=ppme(indx)%v(j)+vpp*phase 
 
1001 continue 
 
!---------------------------------- NN INTERACTION ------------------- 
 
     it = 2 
     ia = a 
     ib = b 
     ic = c 
     id = d 
!----------- CHECK PARITY -------------------------------------------- 
 
     L = orbqn(it,a)%l+ orbqn(it,b)%l+orbqn(it,c)%l+orbqn(it,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' error in parity ' 
        print*,a,b,c,d,j,t,vv
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
      
     phase = 1 
 
     if(ia < ib)then 
        na = ib 
        ib = ia 
        ia = na 
        phase = (-1)**( J+T+(orbqn(it,a)%j+orbqn(it,b)%j)/2)  ! check 
     endif 
 
     if(ic < id)then 
        nc = id 
        id = ic 
        ic = nc 
        phase = phase*(-1)**( J+T+(orbqn(it,c)%j+orbqn(it,d)%j)/2)  ! check 
     endif 
 
     if(ia < ic .or. (ia==ic .and. ib < id))then 
        na = ic 
        nb = id 
        ic = ia 
        id = ib 
        ia = na 
        ib = nb 
     endif 
 
!---------- CONVERT ------------------------------------------------- 
 
     pair1 = ia*(ia-1)/2 + ib 
     pair2 = ic*(ic-1)/2 + id 
     pair1 = NNcouplemap(pair1) 
     pair2 = NNcouplemap(pair2) 
     if(pair1 == -1 .or. pair2 == -1)goto 1002 
     pcpar = XXcouples(it)%pairc(pair1)%par 
     pcref = XXcouples(it)%meref(pcpar) 
     pcstart = XXcouples(it)%mestart(pcpar) 
 
     if(pair1 >= pair2)then 
       indx = (pair1-pcref)*(pair1-pcref-1)/2+pair2-pcref+pcstart 
     else 
       indx = (pair2-pcref)*(pair2-pcref-1)/2+pair1-pcref+pcstart 
     endif 
!--------------- ERROR TRAP ---------------------------------------- 
 
     if(j > nnme(indx)%jmax .or. j < nnme(indx)%jmin)then 
        print*,' error in Js (nn) ',pair1,pair2,indx 
        print*,ia,ib,ic,id,J,T,vv 
        print*,orbqn(it,a)%j,orbqn(it,b)%j,orbqn(it,c)%j,orbqn(it,d)%j 
        print*,nnme(indx)%jmax,nnme(indx)%jmin 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
         
        stop 
     endif 
     nnme(indx)%v(j)=nnme(indx)%v(j)+vnn*phase 
     
1002 continue 
 
!---------------------------------------------- PN-------------------- 
     ! factor due to identical label 
 
     ia = a 
     ib = b 
     ic = c 
     id = d 
 
     factor = 0.5 
     if(a == b)factor = factor*sqrt(2.) 
     if(c == d) factor = factor*sqrt(2.) 
 
     if(ia < ic .or. (ia == ic .and. ib < id))then  ! swap 
        na = ia 
        ia = ic 
        ic = na 
        na = ib 
        ib = id 
        id = na 
     endif 
!----------- CHECK PARITY ------------------------------------------- 
 
     L = orbqn(1,a)%l+ orbqn(2,b)%l+orbqn(1,c)%l+orbqn(2,d)%l 
      
     if( (-1)**(L) ==-1)then 
        print*,' error in parity ' 
        print*,a,b,c,d,j,t,vv 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
     phase = 1 
 
     pair1 = numorb(2)*(ia-1) + ib 
     pair2 = numorb(2)*(ic-1) + id 
     pair1 = PNcouplemap(pair1) 
     pair2 = PNcouplemap(pair2) 
     if(pair1 == -1 .or. pair2 == -1)goto 1003 
     pcpar = PNcouples%pairc(pair1)%par 
     pcref = PNcouples%meref(pcpar) 
     pcstart = PNcouples%mestart(pcpar) 
 
     if(pair1 >= pair2)then 
       indx = (pair1-pcref)*(pair1-pcref-1)/2+pair2-pcref+pcstart 
     else 
       indx = (pair2-pcref)*(pair2-pcref-1)/2+pair1-pcref+pcstart 
     endif 
!--------------- ERROR TRAP --------------- 
      
     if(j > pnme(indx)%jmax .or. j < pnme(indx)%jmin)then 
        print*,' error in Js (pn 1) ',pair1,pair2,indx 
        print*,a,b,c,d 
        print*,ia,ib,ic,id,J,T,vv 
        print*,orbqn(1,ia)%j,orbqn(2,ib)%j,orbqn(1,ic)%j,orbqn(2,id)%j 
        print*,j,pnme(indx)%jmin,pnme(indx)%jmax 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif      
     pnme(indx)%v(j)=pnme(indx)%v(j)+vpn*phase*factor 
 
1003 continue 
 
!------------ second pn CASE ---swap a and b 
 
     ia = b 
     ib = a 
     ic = c 
     id = d 
 
!-------------- don't bother if they work out the same--------------- 
     if(ia == a .and. ib == b)goto 1004 
 
     if(ia < ic .or. (ia == ic .and. ib < id))then  ! swap 
        na = ia 
        ia = ic 
        ic = na 
        na = ib 
        ib = id 
        id = na 
     endif 
!----------- CHECK PARITY -----open--------------------------------------------- 
 
     L = orbqn(2,a)%l+ orbqn(1,b)%l+orbqn(1,c)%l+orbqn(2,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' error in parity ' 
        print*,a,b,c,d,j,t,vv 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
      
     phase = (-1)**( J+T+(orbqn(2,a)%j+orbqn(1,b)%j)/2)  ! check 
 
     pair1 = numorb(2)*(ia-1) + ib 
     pair2 = numorb(2)*(ic-1) + id 
     pair1 = PNcouplemap(pair1) 
     pair2 = PNcouplemap(pair2) 
     if(pair1 == -1 .or. pair2 == -1)goto 1004 
     pcpar = PNcouples%pairc(pair1)%par 
     pcref = PNcouples%meref(pcpar) 
     pcstart = PNcouples%mestart(pcpar) 
 
     if(pair1 >= pair2)then 
       indx = (pair1-pcref)*(pair1-pcref-1)/2+pair2-pcref+pcstart 
     else 
       indx = (pair2-pcref)*(pair2-pcref-1)/2+pair1-pcref+pcstart 
     endif 
!--------------- ERROR TRAP ----------------------------------------- 
 
     if(j > pnme(indx)%jmax .or. j < pnme(indx)%jmin)then 
        print*,' error in Js (pn 2) ',pair1,pair2,indx 
        print*, a, b, c, d 
        print*,ia,ib,ic,id,J,T,vv 
        print*,orbqn(2,a)%j,orbqn(1,b)%j,orbqn(1,c)%j,orbqn(2,d)%j 
        print*,pnme(indx)%jmax,pnme(indx)%jmin 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
             
        stop 
     endif 
     pnme(indx)%v(j)=pnme(indx)%v(j)+vpn*phase*factor 
     
1004 continue 
 
!------------ third pn CASE ---swap c and d---------------------- 
     ia = a  
     ib = b  
     ic = d  
     id = c  
 
!-------------- don't bother if they work out the same--------------- 
     if(ic == c  .and. id == d )goto 1005 
     if(ic == b  .and. id == a )goto 1005    ! this prevents double counting!
      
 
     if(ia < ic .or. (ia == ic .and. ib < id))then  ! swap 
        na = ia 
        ia = ic 
        ic = na 
        na = ib 
        ib = id 
        id = na 
     endif 
!----------- CHECK PARITY -------------------------------------------- 
 
     L = orbqn(1,a)%l+ orbqn(2,b)%l+orbqn(2,c)%l+orbqn(1,d)%l 
      
     if( (-1)**(L) ==-1)then 
        print*,' error in parity ' 
        print*,a,b,c,d,j,t,vv 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
      
     phase = (-1)**( J+T+(orbqn(2,c)%j+orbqn(1,d)%j)/2)  ! check 
 
     pair1 = numorb(2)*(ia-1) + ib 
     pair2 = numorb(2)*(ic-1) + id 
     pair1 = PNcouplemap(pair1) 
     pair2 = PNcouplemap(pair2) 
     if(pair1 == -1 .or. pair2 == -1)goto 1005 
     pcpar = PNcouples%pairc(pair1)%par 
     pcref = PNcouples%meref(pcpar) 
     pcstart = PNcouples%mestart(pcpar) 
 
     if(pair1 >= pair2)then 
       indx = (pair1-pcref)*(pair1-pcref-1)/2+pair2-pcref+pcstart 
     else 
       indx = (pair2-pcref)*(pair2-pcref-1)/2+pair1-pcref+pcstart 
     endif 
 
!--------------- ERROR TRAP ------------------------------------- 
 
     if(j > pnme(indx)%jmax .or. j < pnme(indx)%jmin)then 
        print*,' error in Js (pn 3) ',pair1,pair2,indx 
        print*,ia,ib,ic,id,J,T,vv 
        print*,orbqn(it,a)%j,orbqn(it,b)%j,orbqn(it,c)%j,orbqn(it,d)%j 
        print*,j,pnme(indx)%jmin,pnme(indx)%jmax 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     pnme(indx)%v(j)=pnme(indx)%v(j)+vpn*phase*factor 
    
1005 continue 
 
!------------ fourth pn CASE ---swap a and b and c and d-------------- 
     ia = b  
     ib = a  
     ic = d  
     id = c  
 
!-------------- don't bother if they work out the same----------- 
     if(ia == a  .and. ib == b )goto 1006 
     if(ic == c  .and. id == d )goto 1006 
      
     if(ia < ic .or. (ia == ic .and. ib < id))then  ! swap 
        na = ia 
        ia = ic 
        ic = na 
        na = ib 
        ib = id 
        id = na 
     endif 
!----------- CHECK PARITY -------------------------------------------------- 
 
     L = orbqn(2,a)%l+ orbqn(1,b)%l+orbqn(2,c)%l+orbqn(1,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' error in parity ' 
        print*,a,b,c,d,j,t,vv
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
     phase = (-1)**( (orbqn(2,a)%j+orbqn(1,b)%j+ orbqn(2,c)%j+orbqn(1,d)%j)/2)  ! check 
 
 
 
     pair1 = numorb(2)*(ia-1) + ib 
     pair2 = numorb(2)*(ic-1) + id 
     pair1 = PNcouplemap(pair1) 
     pair2 = PNcouplemap(pair2) 
     if(pair1 == -1 .or. pair2 == -1)goto 1006 
     pcpar = PNcouples%pairc(pair1)%par 
     pcref = PNcouples%meref(pcpar) 
     pcstart = PNcouples%mestart(pcpar) 
 
     if(pair1 >= pair2)then 
       indx = (pair1-pcref)*(pair1-pcref-1)/2+pair2-pcref+pcstart 
     else 
       indx = (pair2-pcref)*(pair2-pcref-1)/2+pair1-pcref+pcstart 
     endif 
 
!--------------- ERROR TRAP -------------------------------------- 
     if(j > pnme(indx)%jmax .or. j < pnme(indx)%jmin)then 
        print*,' error in Js (pn 4) ',pair1,pair2,indx 
        print*,a,b,c,d 
        print*,ia,ib,ic,id,J,T,vv 
        print*,orbqn(it,a)%j,orbqn(it,b)%j,orbqn(it,c)%j,orbqn(it,d)%j 
        print*,j,pnme(indx)%jmin,pnme(indx)%jmax
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     pnme(indx)%v(j)=pnme(indx)%v(j)+vpn*phase*factor 
 
1006 continue 
 
  enddo  !i  
  close(1)

!........... WARNING MESSAGES.... 
  if(maxsp < maxorblabel)then
     if(iproc==0)then
      print*,' '
      print*,' WARNING '
      print*,' This is missing matrix elements within the active single-particle space '
      print*,' i.e., missing matrix elements to single particle orbits from ',maxsp+1, ' to ',maxorblabel
      print*,' May cause problems in runs ! '
      print*,' '
      write(logfile,*)' WARNING '
      write(logfile,*)' This is missing matrix elements within the active single-particle space '
      write(logfile,*)' i.e., missing matrix elements to single particle orbits from ',maxsp+1, ' to ',maxorblabel
      write(logfile,*)' May cause problems in runs ! '
      end if
  end if
  firstread= .false.
 
  return 
end subroutine readv2bme_mfd 
 
!============================================================ 
!
! reads two-body matrix elements in explicit proton-neutron formalism
!
! added 7.7.7 (April 2017)
! 
! 
!   THE TBMEs V_JT(ab,cd) are ordered in a particular way: 
!   a >= b, c >= d;  a >= c 
!     VTBME(indx): 
!      indx =  pair1*(pair1-1)/2+pair2 
!    WHERE     pair1 = ia*(ia-1)/2 + ib 
!              pair2 = ic*(ic-1)/2 + id 
! 
!===================================================================== 
subroutine readv2bme_xpn(firstread,formattedfile,intcase) 
   
  use io 
  use system_parameters 
  use sporbit 
  use W_info 
  use coupledmatrixelements 
  use nodeinfo 
  use bmpi_mod 
  implicit none 
 
  logical :: firstread,formattedfile 
  integer :: intcase 
   
  integer :: a,b,c,d 
  integer :: ia,ib,ic,id 
  integer :: na,nb,nc,nd 
  integer :: pair1,pair2,indx 
  real :: Vv 
  real, allocatable :: spetmp(:) 
  logical :: normalizedpn
  integer :: neutron_offset   ! offset for neutron indices
  

  integer(4) :: ierr 
  integer :: aerr 
  integer j, t, t_ab, t_cd 

  real :: ppscale,nnscale,pnscale,pspescale,nspescale
  real:: ax,bx,x,vscale,spscale
  logical :: usescale = .true. 
  integer phase 
  integer nme 
  real factor 
  integer ipar,iref,istart 
  integer it 
 
  real(kind=4) :: factor_v,sc2
  integer :: maxsp,minsp    ! used for error trap  
  logical autoscale  ! added 7.8.1 for reading in NuShell files
  integer :: numval
  
!----------------------------------------------- 
  character*1 ychar 
   
  character*70 title 
!-----------------------------------------------  
!      dummy counters 
!-------------------------------------------------- 
  integer i,tmp,m,mmax 
     
  integer L 
  logical success 
  logical smint   !  a successor to .int  
  logical finished,foundit
  character(3) :: formattest  
!--------- KSM - initialize for -Wuninitialized 
! set to values that would cause problems if used 
  it = -10000000 
  factor_v = -1e20 
  
 

  select case(intcase)
  case (4)
  
  if(iproc==0)print*,' Assuming normalized p-n matrix elements '
  normalizedpn = .true. 
  case (5)
  
  if(iproc==0)print*,' Assuming unnormalized p-n matrix elements '
  normalizedpn = .false.
  
  case default
  
  print*,' I do not recognize that value of intcase ',intcase
  stop 
     
  end select
  
!.......... PRINT OUT MAPPING................  


  neutron_offset = numorb(1)
  
  if(iproc==0)then
	  write(6,*)' Here is my understanding of the mapping of states '
	  if(numorb(1)==numorb(2))then
		  write(6,*)'    N    L    2J    proton label  neutron label  '
		  do i = 1,numorb(1)
			  write(6,'(2x,3i5,10x,i3,10x,i3)')orbqn(1,i)%nr,orbqn(1,i)%l, orbqn(1,i)%j, i,i+neutron_offset
			  
			  
		  end do
		  
	  else
		  print*,' Not yet implemented for unequal proton and neutron spaces '
		  
	  end if
		  
  end if

!-------------- READ PAST TITLE CARDS --------------------------- 
  success = .false. 
  do while(.not.success) 
     read(1,'(a)')ychar 
     backspace(1) 
     if(ychar /= '#' .and. ychar /= '!')then 
        success = .true. 
     else 
        read(1,'(a)')title 
        write(6,'(a)')title 
		write(logfile,'(a)')title
		call fileformatcheck(title,foundit,formattest)
		if(foundit)then ! check format
		   if((normalizedpn .and. formattest /= 'xpn') .or. & 
		      (.not.normalizedpn .and. formattest /='upn')  )then
		   print*,' '
		   print*,'  WARNING !   '
		   print*,' File appears to be in wrong format '
		   print*,' Header denotes format as ',formattest
		   print*,' Hint: If you want to read in xpn/upn format, '
		   print*,' you must first enter xpn/upn before entering file name. '
		   print*,' See manual for more details.'
	       print*,'  SKIPPING FILE !   '
 		      print*,' '
			  write(logfile,*)' File appears to be in wrong format '
			  write(logfile,*)' Header denotes format as ',formattest
		      write(logfile,*)'  SKIPPING FILE !   '
			  return
		  end if
	    end if
     end if 
  end do 

!---------- CHECK FOR AUTOSCALING ---------

   read(1,*)nme
   if(nme < 0)then
	   autoscale = .true.
   else
	   autoscale = .false.
   end if
   backspace(1)  
!-------------- ENTER SCALING ----------------------------- 
  ppscale = 1.0
  nnscale = 1.0
  pnscale = 1.0
  pspescale=1.0
  nspescale=1.0
  if(.not.autoscale)then
  if(auto_input)then 
	  
	  read(autoinputfile,*)spscale,ax,bx,x
      if ( bx == 0. .or. x == 0.0 ) then 
         vscale = ax 
      else 
         vscale = (ax/bx)**x 
      end if 
      print*,' Scaling ',spscale,ax,bx,x,' TBME: ',vscale 
      write(logfile,*)   ' Global scaling: spe ',spscale,' TBMEs by (',ax,'/',bx,')^',x,' = ',vscale     
	  read(autoinputfile,*)pspescale,nspescale,ppscale,nnscale,pnscale
	  write(logfile,*)' Additional scaling: proton spe: ',pspescale,', neutron spe: ',nspescale, & 
	  ', pp TBME: ',ppscale,', nn TBME: ',nnscale,', pn TBME: ',pnscale
	  write(6,*)' Additional scaling: proton spe: ',pspescale,', neutron spe: ',nspescale, & 
	  ', pp TBME: ',ppscale,', nn TBME: ',nnscale,', pn TBME: ',pnscale

  else 
	  if(usescale)then
          print*,' Enter global scaling for spes, A,B,X ( (A/B)^X ) for TBMEs '
          print*,' (If B or X = 0, then scale by A ) ' 
          read*,spscale,ax,bx,x 
          if ( bx == 0. .or. x == 0.0 ) then 
             vscale = ax 
          else 
             vscale = (ax/bx)**x 
          end if 

          write(autoinputfile,'(4f10.4,"        ! scaling: xspe, A0,A,x ")' )spscale,ax,bx,x 
          if(iproc==0)print*,' Scaling ',spscale,ax,bx,x,vscale 
          write(logfile,*)   ' Global scaling: spe ',spscale,' TBMEs by (',ax,'/',bx,')^',x,' = ',vscale     
		  print*,' Enter individual scaling for: proton spes, neutron spes, pp TBMEs, nn TBMEs, pn TBMES '
		  print*,' (If not sure, just enter 1 1 1 1 1 )'
		  read(5,*)pspescale,nspescale,ppscale,nnscale,pnscale
		  write(logfile,*)' Additional scaling: proton spe: ',pspescale,', neutron spe: ',nspescale, & 
		  ', pp TBME: ',ppscale,', nn TBME: ',nnscale,', pn TBME: ',pnscale
		  write(autoinputfile,'(5f10.4,"    ! scaling: pspe, nspe, pp TBME, nn TBME, pn TBME ")')pspescale,nspescale,&  
                  ppscale,nnscale,pnscale
 
		  
	  else
		  print*,' scaling turned off for xpn/upn formats '
	  end if
    

  endif 
end if
!-------------- READ IN SPEs --------------
  print*,' * * '
  print*,' * * NOTICE: I expect single-particles space with ',numorb(1)+numorb(2),' orbits '
  print*,' * * ' 
  
  if(autoscale)then
	  numval = numorb(1)+numorb(2)+3
  else
	  numval = numorb(1)+numorb(2)
	  
  endif
  allocate(spetmp(numval), stat=aerr) 
  if(aerr /= 0) call memerror("readv2bme_xpn 2") 
  
  if(autoscale)then
	  read(1,*,err=1899,end=1899)nme,(spetmp(i),i= 1,numval)
	  nme = -nme
	  spscale = 1
	  pspescale = 1
	  nspescale = 1
		  bx = spetmp(numorb(1)+numorb(2)+1)+ npeff(1)+npeff(2)
		  ax = spetmp(numorb(1)+numorb(2)+2)
		  x = spetmp(numorb(1)+numorb(2)+3)
	      vscale = (ax/bx)**x 
	      if(iproc==0)then
	 		 print*,' Autoscaling ',ax,bx,x,vscale     
			 print*,' I think I have a core of ', spetmp(numorb(1)+numorb(2)+1)
	          write(logfile,*)   ' Autoscaling  TBMEs by (',ax,'/',bx,')^',x,' = ',vscale      
	 	 end if	  
		 ppscale = 1
		 nnscale = 1
		 pnscale = 1
	     do i = 1,numorb(1)
	   	  pspe(i)=pspe(i) +spetmp(i)*pspescale*spscale
	     end do
  
	     do i = 1,numorb(2)
	   	  nspe(i)=nspe(i)+spetmp(i+neutron_offset)*nspescale*spscale
	     end do
	  
  else
	  
	  
  read(1,*,err=1899,end=1899)nme,(spetmp(i),i= 1,MIN(10,numorb(1)+numorb(2))) 

 
  if(numorb(1)+numorb(2) > 10)then 
     do m = 10,numorb(1)+numorb(2)-1,10 
        mmax = MIN(10+m,numorb(1)+numorb(2)) 
        read(1,*,err=1899,end=1899)(spetmp(i),i=1+m,mmax) 
     end do 
  endif 
  do i = 1,numorb(1)
	  pspe(i)=pspe(i) +spetmp(i)*pspescale*spscale
  end do
  
  do i = 1,numorb(2)
	  nspe(i)=nspe(i)+spetmp(i+neutron_offset)*nspescale*spscale
  end do
end if
 
  deallocate( spetmp ) 

!-------------- READ IN TBMEs ---------------- 
  write(logfile,*)nme,' two-body matrix elements in this file '
  write(6,*)nme,' two-body matrix elements in this file '
  
  maxsp = 0 
  

  do i = 1,nme 
!---------- ERROR TRAP ADDED July 2011 CWJ ------------------' 
     
      read(1,*,err=1899,end=1899)ia,ib,ic,id,j,t,vv 
 
!---------- ERROR TRAP ADDED July 2011 CWJ ------------------' 
 
      if( t /= 0 .and. t /= 1)then 
            print*,' bad TBME matrix element; T value is bad ',t 
            print*,' matrix element # ',i 
            print*,ia,ib,ic,id,j,t,vv 
            print*,' Check that single particle space matches hamiltonian ' 
            print*,' # of single-particle orbits = ',numorb(1) +numorb(2)
            stop 
       end if 

       maxsp = MAX(maxsp,ia)
       maxsp = MAX(maxsp,ib)
       maxsp = MAX(maxsp,ic)
       maxsp = MAX(maxsp,id)
 
!----------- PUT INTO CORRECT ORDER; PICK UP PHASES ------------- 
!           "CORRECT" ORDER: a >= b, c >= d 
 
!-------------------- PP interaction ------------------------- 
     if( ia > neutron_offset .or. ib > neutron_offset .or. ic > neutron_offset .or. id > neutron_offset)goto 1001
     if( np(1) < 2 .and. .not.phconj(1))cycle       ! order of this statement is important - 7/2011 CWJ 

 
     it = 1 
 
     a = ia 
     b = ib 
     c = ic  
     d = id 
!----------- CHECK PARITY ------------------- 
     L = orbqn(it,a)%l+ orbqn(it,b)%l+orbqn(it,c)%l+orbqn(it,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' Oops! Error in parity in (pp) matrix element # ',i
		print*,' Most likely this is a mismatch between single-particle orbits ' 
		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
        write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
		     orbqn(it,a)%j,orbqn(it,a)%l,  orbqn(it,b)%j,orbqn(it,b)%l,  orbqn(it,c)%j,orbqn(it,c)%l,  orbqn(it,d)%j,orbqn(it,d)%l
			 
	    write(logfile,*)' Oops! Error in parity in (pp) matrix element # ',i
	 	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
	    write(logfile,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
	 		     orbqn(it,a)%j,orbqn(it,a)%l,  orbqn(it,b)%j,orbqn(it,b)%l,  orbqn(it,c)%j,orbqn(it,c)%l,  orbqn(it,d)%j,orbqn(it,d)%l	
		close(logfile)		 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
     phase = 1 
 
 
     if(a < b)then 
        na = b 
        b = a 
        a = na 
        phase = (-1)**( J+T+(orbqn(it,a)%j+orbqn(it,b)%j)/2)  ! check 
     endif 
 
     if(c < d)then 
        nc = d 
        d = c 
        c = nc 
        phase = phase*(-1)**( J+T+(orbqn(it,c)%j+orbqn(it,d)%j)/2)  ! check 
     endif 
 
     if(a < c .or. (a==c .and. b < d))then 
        na = c 
        nb = d 
        c = a 
        d = b 
        a = na 
        b = nb 
     endif 
 
!---------- CONVERT ------------------------- 
 
     pair1 = a*(a-1)/2 + b 
     pair2 = c*(c-1)/2 + d 
     if( PPcouplemap(pair1) == -1) cycle !goto 1001 
     if( PPcouplemap(pair2) == -1) cycle !goto 1001 
     pair1 = PPcouplemap(pair1) 
     pair2 = PPcouplemap(pair2) 
 
     ipar = XXcouples(it)%pairc(pair1)%par 
     if(ipar /= XXcouples(it)%pairc(pair2)%par)then 
        print*,' problem with parity, boss  (PP) ',a,b,c,d
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     iref = XXcouples(it)%meref(ipar) 
     istart = XXcouples(it)%mestart(ipar) 
     if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif 
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!--------------- ERROR TRAP --------------- 
 
     if(j > ppme(indx)%jmax .or. j < ppme(indx)%jmin)then 
		 
        print*,' Oops! Error in adding Js in (pp) matrix element # ',i
 		print*,' Most likely this is a mismatch between single-particle orbits ' 
 		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
        write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
 		     orbqn(it,a)%j, orbqn(it,b)%j,  orbqn(it,c)%j, orbqn(it,d)%j
		write(6,*)' See log file for additional information '	 
	    write(logfile,*)' Oops! Error in adding Js in (pp) matrix element # ',i
	  	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
	    write(logfile,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
	  		     orbqn(it,a)%j, orbqn(it,b)%j,  orbqn(it,c)%j, orbqn(it,d)%j

		write(logfile,*)' Additional information ' 
        write(logfile,*)' J = ',j,' but min/max for these pairs are ',ppme(indx)%jmin,ppme(indx)%jmax 

 		close(logfile)		
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     ppme(indx)%v(j)=ppme(indx)%v(j)+vv*ppscale*phase  *vscale
	 cycle    ! here we can only have pp, OR nn, OR pn
1001 continue 
 
!.......................................................................	 
!---------------------------------- NN INTERACTION ------------- 
!.......................................................................	 
    if( ia <= neutron_offset .or.  ic <= neutron_offset )goto 1002
! print*,i,' reading NN '

    if(ic <= neutron_offset .or. id <= neutron_offset)then
		print*,' Some error should not have reached here '
		print*,' NN ?'
		print*,ia,ib,ic,id
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
		stop
	end if

     if(np(2) < 2 .and. .not.phconj(2))cycle
 
     it = 2 
     a = ia -neutron_offset
     b = ib -neutron_offset
     c = ic  -neutron_offset
     d = id -neutron_offset
!----------- CHECK PARITY ------------------- 
 
     L = orbqn(it,a)%l+ orbqn(it,b)%l+orbqn(it,c)%l+orbqn(it,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' Oops! Error in parity in matrix element # ',i
 		print*,' Most likely this is a mismatch between single-particle orbits ' 
 		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
         write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
 		     orbqn(it,a)%j,orbqn(it,a)%l,  orbqn(it,b)%j,orbqn(it,b)%l,  orbqn(it,c)%j,orbqn(it,c)%l,  orbqn(it,d)%j,orbqn(it,d)%l
			 
 	    write(logfile,*)' Oops! Error in parity in matrix element # ',i
 	 	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
 	 		     orbqn(it,a)%j,orbqn(it,a)%l,  orbqn(it,b)%j,orbqn(it,b)%l,  orbqn(it,c)%j,orbqn(it,c)%l,  orbqn(it,d)%j,orbqn(it,d)%l	
 		close(logfile)
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
     phase = 1 
      
     if(a < b)then 
        na = b 
        b = a 
        a = na 
        phase = (-1)**( J+T+(orbqn(it,a)%j+orbqn(it,b)%j)/2)  ! check 
     endif 
 
     if(c < d)then 
        nc = d 
        d = c 
        c = nc 
        phase = phase*(-1)**( J+T+(orbqn(it,c)%j+orbqn(it,d)%j)/2)  ! check 
     endif 
 
     if(a < c .or. (a==c .and. b < d))then 
        na = c 
        nb = d 
        c = a 
        d = b 
        a = na 
        b = nb 
     endif 
 
     pair1 = a*(a-1)/2 + b 
     pair2 = c*(c-1)/2 + d 
     if( NNcouplemap(pair1) == -1) cycle !goto 1002 
     if( NNcouplemap(pair2) == -1) cycle !goto 1002 
     pair1 = NNcouplemap(pair1) 
     pair2 = NNcouplemap(pair2) 
 
     ipar = XXcouples(it)%pairc(pair1)%par 
     if(ipar /= XXcouples(it)%pairc(pair2)%par)then 
        print*,' problem with parity, boss (NN in xpn) ',a,b,c,d 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     iref = XXcouples(it)%meref(ipar) 
     istart = XXcouples(it)%mestart(ipar) 
      if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif  
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
!--------------- ERROR TRAP --------------- 
 
     if(j > nnme(indx)%jmax .or. j < nnme(indx)%jmin)then  
        print*,' Oops! Error in adding Js in (nn) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
        write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
  		     orbqn(it,a)%j, orbqn(it,b)%j,  orbqn(it,c)%j, orbqn(it,d)%j
 		write(6,*)' See log file for additional information '	 
 	    write(logfile,*)' Oops! Error in adding Js in (pp) matrix element # ',i
 	  	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
 	  		     orbqn(it,a)%j, orbqn(it,b)%j,  orbqn(it,c)%j, orbqn(it,d)%j

 		write(logfile,*)' Additional information ' 
         write(logfile,*)' J = ',j,' but min/max for these pairs are ',nnme(indx)%jmin,nnme(indx)%jmax 
		 
!        print*,' error in Js (nn) ',pair1,pair2,indx 
!        print*,orbqn(it,a)%j,orbqn(it,b)%j,orbqn(it,c)%j,orbqn(it,d)%j 
!        print*,nnme(indx)%jmax,nnme(indx)%jmin 
        close(logfile)
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     nnme(indx)%v(j)=nnme(indx)%v(j)+vv*nnscale*phase *vscale
 
!............................................................... 
!     if ( nproc > 1 ) then 
!        n_nn = n_nn + 1 
!        nn_indx_J(1,n_nn) = indx 
!        nn_indx_J(2,n_nn) = j 
!        nn_bcast(n_nn) = nnme(indx)%v(j) 
!     end if 
!............................................................... 
     cycle
1002 continue 
!.......................................................................	 
!--- PN --- PN --- PN  --- PN --- PN ---- PN ----- PN ---- PN ------
!.......................................................................	 
 
!print*,i,' reading PN '

     if( np(1)*np(2) ==0 .and. .not.phconj(1) .and. .not.phconj(2))cycle 

!...........PN  ......................................
!
!  OBLIGATORY ORDERING:  ia,ic proton, ib,id neutron
     a = ia 
     b = ib - neutron_offset
     c = ic  
     d = id -neutron_offset
	 
	 if(a > numorb(1) .or. c > numorb(1))then
		 print*,' wrong proton labels in element ',i
		 print*,ia,ib,ic,id,j,t,vv
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
		 stop
	 end if
	 if(b > numorb(2) .or. d > numorb(2) .or. b < 1 .or. d < 1)then
		 print*,' wrong neutron labels in element ',i
		 print*,ia,ib,ic,id,j,t,vv
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
		 stop
	 end if
	 
     ! factor due to identical label 
	 ! factor should ALREADY be accounted for, if "normalized"
	 ! if "unnormalized" must take into account
	 if(normalizedpn)then
   	    factor=1.0
	else
     factor = 0.5 
     if(a == b)factor = factor*sqrt(2.) 
     if(c == d) factor = factor*sqrt(2.) 
	 
	 if(a /=b .and. a ==d .and. b == c)factor=0.25
    end if	 
!	 print*,i,a,b,c,d,factor
 
!----------- CHECK PARITY ------------------- 
 
     L = orbqn(1,a)%l+ orbqn(2,b)%l+orbqn(1,c)%l+orbqn(2,d)%l 
 
     if( (-1)**(L) ==-1)then 
        print*,' Oops! Error in parity in (pn 1) matrix element # ',i
 		print*,' Most likely this is a mismatch between single-particle orbits ' 
 		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
         write(6,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
 		     orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
			 
 	    write(logfile,*)' Oops! Error in parity in matrix element # ',i
 	 	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ,Ls are (respectively) ",4(i4,2x,i4,","))')  & 
           orbqn(1,a)%j,orbqn(1,a)%l,  orbqn(2,b)%j,orbqn(2,b)%l,  orbqn(1,c)%j,orbqn(1,c)%l,  orbqn(2,d)%j,orbqn(2,d)%l
 		close(logfile)		
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
 
!-------- NO PHASE, ALREADY IN CORRECT ORDER 
     phase = 1 
 
     pair1 = numorb(2)*(a-1) + b 
     pair2 = numorb(2)*(c-1) + d 
     if( PNcouplemap(pair1) == -1) cycle
     if( PNcouplemap(pair2) == -1) cycle
     pair1 = PNcouplemap(pair1) 
     pair2 = PNcouplemap(pair2) 
 
     ipar = PNcouples%pairc(pair1)%par 
     if(ipar /= PNcouples%pairc(pair2)%par)then 
        print*,' problem with parity, boss (pn 1)  ' 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
        stop 
     endif 
     if(pair1 < pair2)then 
         tmp = pair1 
         pair1 = pair2 
         pair2 = tmp 
     endif 
     iref = PNcouples%meref(ipar) 
     istart = PNcouples%mestart(ipar) 
      
     indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart    
 
!--------------- ERROR TRAP --------------- 
 
     if(j > pnme(indx)%jmax .or. j < pnme(indx)%jmin)then 
		 
         print*,' Oops! Error in adding Js in (pn1) matrix element # ',i
  		print*,' Most likely this is a mismatch between single-particle orbits ' 
  		write(6,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
         write(6,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
    	     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j
 		write(6,*)' See log file for additional information '	 
 	    write(logfile,*)' Oops! Error in adding Js in (pp) matrix element # ',i
 	  	write(logfile,*)' Most likely this is a mismatch between single-particle orbits ' 
 	    write(logfile,'(" a b c d = ",4i4," J T = ",2i4," V_JT(ab,cd) = ",f10.5)')a,b,c,d,j,t,vv 
 	    write(logfile,'(" I think orbital 2xJ are (respectively) ",4(i4,","))')  & 
 	  		     orbqn(1,a)%j, orbqn(2,b)%j,  orbqn(1,c)%j, orbqn(2,d)%j

 		 write(logfile,*)' Additional information ' 
         write(logfile,*)' J = ',j,' but min/max for these pairs are ',pnme(indx)%jmin,pnme(indx)%jmax 		
		 close(logfile) 
		print*,' STOPPING RUN '
	    call BMPI_Abort(icomm,101,ierr)
 
        stop 
     endif 
 
     pnme(indx)%v(j)=pnme(indx)%v(j)+vv*pnscale*phase*factor *vscale

  enddo  !nme 
!........... WARNING MESSAGES.... 
  if(maxsp < maxorblabel)then
     if(iproc==0)then
      print*,' '
      print*,' WARNING '
      print*,' This is missing matrix elements within the active single-particle space '
      print*,' May cause problems in runs ! '
      print*,' '
	  
      write(logfile,*)' WARNING '
      write(logfile,*)' This is missing matrix elements within the active single-particle space '
      write(logfile,*)' May cause problems in runs ! '
      end if
  end if
  if(maxsp > (numorb(1)+numorb(2)) )then
     if(iproc==0)then
      print*,' '
      print*,' WARNING '
      print*,' This file contains matrix elements beyond the active single-particle space '
      print*,' i.e., missing matrix elements coupling single particle orbits '
      print*,' from ',(numorb(1)+numorb(2)), ' to ',maxorblabel
      print*,' May cause problems in runs ! '
      print*,' END WARNING '
      print*,' '
	  
      write(logfile,*)' WARNING '
      write(logfile,*)' This file contains matrix elements beyond the active single-particle space '
      write(logfile,*)' i.e., missing matrix elements coupling single particle orbits '
      write(logfile,*)' from ',(numorb(1)+numorb(2)), ' to ',maxorblabel
      write(logfile,*)' May cause problems in runs ! '
      write(logfile,*)' END WARNING '
      end if
  end if

!......... ADDED IN 7.8.2.... GO ON TO LOOK FOR OFF-DIAGONAL SINGLE-PARTICLE  
  call readinsppot('xpn',1,pspescale*spscale,nspescale*spscale,formattedfile,foundit)
  if(foundit)then
	  print*,' Finished reading in off-diagonal one-body potential '
	  write(logfile,*)' Finished reading in off-diagonal one-body potential '
  else
	  print*,' No one-body potential beyond single-particle energies '
	  write(logfile,*)' No one-body potential beyond single-particle energies '

  end if
  
  close(1)
  firstread=.false.
  return

1899 continue
  if(iproc==0)then
    print*,' Error in reading file; problem may be incommensurate size of single-particle space '
    print*,' I expect single-particles space with ',numorb(1)+numorb(2),' orbits '
    print*,nme,i
    backspace(1)
    read(1,'(a)')title
    write(6,*)title
  end if
    call BMPI_Abort(icomm,101,ierr)
  	
	
	return
end subroutine readv2bme_xpn
 
! =======================================================================
!  subroutine to compute minimal set of coupled pairs for XX (pp or nn) TBMEs 
!  uses information from restricted quantum numbers to determine max W of pairs 
!  added 7/2010 by CWJ @ SDSU 
! 
! INPUT:  
!  it = species (1 = proton, 2 = neutron) 
!  create : if F then just count; if T then fill in 
! 
 
subroutine setup_tbme_XXcouples(it, create) 
 
  use sporbit 
  use system_parameters 
  use ntuple_info 
  use W_info 
  use coupledmatrixelements 
  use pairdef 
  implicit none 
  integer :: it ! = species of XX, either pp or nn 
  logical :: create  ! if count or just create 
   
  integer :: nrawpairs 
  integer :: a,b 
  integer :: ncount 
  integer :: W 
  integer, pointer :: cmap(:) 
  integer :: indx 
  integer :: par 
  integer :: neven,ieven 
  integer :: i,j,itmp 
  type (coupair_qn) :: cptemp 
  integer :: nv 
  integer :: aerr 
 
  nrawpairs = numorb(it)*(numorb(it) +1)/2 
 
  if(create)then 
     if( ncouplesXX(it) == 0)return 
     allocate( XXcouples(it)%pairc( ncouplesXX(it) ) , stat=aerr) 
     if(aerr /= 0) call memerror("setup_tbme_XXcouples 1") 
     if ( it == 1 ) then 
        allocate( PPcouplemap(nrawpairs) , stat=aerr) 
        if(aerr /= 0) call memerror("setup_tbme_XXcouples 2") 
        cmap => PPcouplemap 
     else 
        allocate( NNcouplemap(nrawpairs) , stat=aerr) 
        if(aerr /= 0) call memerror("setup_tbme_XXcouples 3") 
        cmap => NNcouplemap 
     end if 
     cmap(:) = -1 
 
  else 
     if ( it == 1 ) then 
        cmap => PPcouplemap 
     else 
        cmap => NNcouplemap 
     end if 
  endif 
 
  ncount = 0 
  do a = 1,numorb(it) 
     if( orbqn(it,a)%w > wceiling(it) .and. .not. phconj(it) .and. .not. phconj(3-it) ) cycle 
     do b = 1,a 
         if( orbqn(it,b)%w > wceiling(it) .and. .not. phconj(it) .and. .not. phconj(3-it) ) cycle 
         W =  orbqn(it,a)%w +  orbqn(it,b)%w 
         if(restrict_ntuples .and. W > pairlim(it)%maxW .and. .not. phconj(it) .and. .not. phconj(3-it) )cycle  ! KEY STEP 
                ! THROW OUT COUPLED PAIRS THAT EXCEED LIMIT ON W 
         ncount = ncount + 1 
         if(create)then 
            indx = a*(a-1)/2 + b 
            cmap(indx) = ncount 
            XXcouples(it)%pairc(ncount)%ia = a 
            XXcouples(it)%pairc(ncount)%ib = b 
            XXcouples(it)%pairc(ncount)%indx = indx 
            par = orbqn(it,a)%par*orbqn(it,b)%par 
            par = (3-par)/2 
            XXcouples(it)%pairc(ncount)%par = par 
         end if 
 
     end do ! b 
 
  end do ! a 
 
  if(.not. create)then 
    ncouplesXX(it) = ncount 
!    print*,ncouplesXX(it),' coupled pairs kept ',it 
  end if 
 
!---------------- SORT ! on parity only ------------------- 
  if(create)then 
     if(allsameparity)then 
         XXcouples(it)%mestart(1) = 0 
         XXcouples(it)%meref(1) = 0 
         XXcouples(it)%mestart(2) =-1 
         XXcouples(it)%meref(2) =ncouplesXX(it) 
     else 
!............ CLEVERSORT on parity............ 
!             FIRST count up # of even parities 
         neven = 0 
         do i = 1,ncouplesXX(it) 
            if( XXcouples(it)%pairc(i)%par ==1)neven = neven+1 
         end do 
         if(neven == 0 .or. neven == ncouplesXX(it))then 
           if(neven == ncouplesXX(it))then 
              XXcouples(it)%mestart(1) = 0 
              XXcouples(it)%meref(1) = 0 
              XXcouples(it)%mestart(2) =-1 
              XXcouples(it)%meref(2) =ncouplesXX(it) 
           else 
              XXcouples(it)%mestart(1) =-1 
              XXcouples(it)%meref(1) = 0 
              XXcouples(it)%mestart(2) =0 
              XXcouples(it)%meref(2) =0 
           endif 
         else 
            ieven = 0 
            i = 1 
            j = ncouplesXX(it) 
            do while( ieven /= neven) 
               if( XXcouples(it)%pairc(i)%par == 2)then 
                   do while( XXcouples(it)%pairc(j)%par ==2 ) 
                       j = j -1 
                       if(j <= i)then 
                           print*,' whoops missed ' 
                           stop 
                       endif 
                   end do 
!----------------------- SWAP 
                   cptemp = XXcouples(it)%pairc(i) 
                   XXcouples(it)%pairc(i) = XXcouples(it)%pairc(j) 
                   XXcouples(it)%pairc(j) = cptemp 
                   cmap(  XXcouples(it)%pairc(i)%indx )  =i 
                   cmap(  XXcouples(it)%pairc(j)%indx )  =j 
                    
               endif 
               ieven = ieven + 1 
               i = i + 1 
            end do 
!..................... FIND START, FINISH................ 
            XXcouples(it)%mestart(1) = 0 
            XXcouples(it)%meref(1) = 0    
            XXcouples(it)%mestart(2) = neven*(neven+1)/2 
            XXcouples(it)%meref(2) = neven            
 
         endif 
     endif 
 
  endif 
 
 
  return 
end subroutine setup_tbme_XXcouples 
 
!============================================================ 
! 
! set up storage of coupled PP,NN TBMEs 
!  
!  started 7/2010 by CWJ @SDSU 
! 
 
subroutine set_tbme_XXarray(it) 
 
  use sporbit 
  use coupledmatrixelements 
  use butil_mod
  implicit none 
 
  integer it  ! species, 1 = P, 2 =N 
 
  integer cpair,dpair,ddpair,ccpair 
  integer pair1,pair2,tmp 
  integer dparity 
  integer iref,istart 
  integer itbme 
  integer a,b,c,d 
  integer ja,jb,jc,jd 
  integer jmin,jmax 
  integer nv 
  integer ncpair 
  type (vjs), pointer :: xxme(:) 
  integer, pointer :: cmap(:) 
  integer :: aerr 
 
!.................... COMPUTE # OF MATRIX ELEMENTS 
  nv = 0 
  ncpair = XXcouples(it)%meref(2) 
  nv = ncpair*(ncpair+1)/2 
  ncpair = ncouplesXX(it) - XXcouples(it)%meref(2) 
  nv = nv + ncpair*(ncpair+1)/2 
  nv2bmedim(it) = nv 
  if(nv == 0)return 
  if(it == 1)then 
     allocate( ppme( nv2bmedim(it) ), stat=aerr) 
     if(aerr /= 0) call memerror("set_tbme_XXarray 1") 
     xxme => ppme 
     cmap => PPcouplemap 
  else 
     allocate( nnme( nv2bmedim(it) ), stat=aerr) 
     if(aerr /= 0) call memerror("set_tbme_XXarray 2") 
     xxme => nnme 
     cmap => NNcouplemap 
  endif 
 
  do dpair = 1,ncouplesXX(it) 
 
     dparity = XXcouples(it)%pairc(dpair)%par 
     iref = XXcouples(it)%meref(dparity) 
     istart = XXcouples(it)%mestart(dparity) 
     c = XXcouples(it)%pairc(dpair)%ia 
     d = XXcouples(it)%pairc(dpair)%ib 
     jc = orbqn(it,c)%j 
     jd = orbqn(it,d)%j 
 
     do cpair = 1, dpair 
 
          if( XXcouples(it)%pairc(cpair)%par /= dparity )cycle 
          a = XXcouples(it)%pairc(cpair)%ia 
          b = XXcouples(it)%pairc(cpair)%ib 
 
          ja = orbqn(it,a)%j 
          jb = orbqn(it,b)%j 
          jmin = MAX( abs(ja-jb), abs(jc-jd) )/2 
          jmax = MIN(ja+jb,jc+jd)/2 
           
          itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref 
          xxme(itbme)%jmin = jmin 
          xxme(itbme)%jmax = jmax 
          if(jmax >= jmin) then 
 
              allocate( xxme(itbme)%v(jmin:jmax), stat=aerr )  
              if(aerr /= 0) call memerror("set_tbme_XXarray 3") 
              xxme(itbme)%v(:) = 0.0 
 
          endif       
     end do  ! cpair 
 
  end do ! dpair 
 
 
  return 
end subroutine set_tbme_XXarray 
 
!============================================================ 
! 
!  subroutine to compute minimal set of coupled pairs for PN TBMEs 
!  uses information from restricted quantum numbers to determine max W of pairs 
!  added 7/2010 by CWJ @ SDSU 
!
!  NOTE: If we have particle-hole conjugation allow for all possible pairs
! 
! INPUT:  
!  create : if F then just count; if T then fill in 
! 
 
subroutine setup_tbme_PNcouples(create) 
 
  use system_parameters 
  use sporbit 
  use ntuple_info 
  use W_info 
  use coupledmatrixelements 
  use pairdefcoupled 
  implicit none 
  integer it ! = species of XX, either pp or nn 
  logical create  ! if count or just create 
   
  integer nrawpairs 
  integer a,b 
  integer ncount 
  integer W 
  integer indx 
  integer par 
  integer neven,ieven 
  integer i,j,itmp 
  type (coupair_qn) :: cptemp 
  integer :: aerr 
   
 
  nrawpairs = numorb(1)*numorb(2) 
  if(nrawpairs == 0)return 
 
!...... NEED TO COMPUTE NORB_ALLOW 
 
  if(.not.create)then 
     norb_allow(:) = 0 
     do it = 1,2 
        do a = 1,numorb(it) 
           if(orbqn(it,a)%w <= wceiling(it))norb_allow(it) = norb_allow(it)+1 
        end do 
     end do 
  end if  
 
!................. 
 
  if(create)then 
     if( ncouplesPN == 0)return 
     allocate( PNcouples%pairc( ncouplesPN ) , stat=aerr) 
     if(aerr /= 0) call memerror("setup_tbme_PNcouples 1") 
     allocate( PNcouplemap(nrawpairs) , stat=aerr) 
     if(aerr /= 0) call memerror("setup_tbme_PNcouples 2") 
     PNcouplemap(:) = -1 
  endif 
 
  ncount = 0 
  
  do a = 1,numorb(1) 
     if( orbqn(1,a)%w > wceiling(1) .and. .not.phconj(1) .and. .not.phconj(2)) cycle 
     do b = 1,numorb(2) 
         if( orbqn(2,b)%w > wceiling(2) .and. .not.phconj(1) .and. .not.phconj(2)) cycle 
         W =  orbqn(1,a)%w +  orbqn(2,b)%w 
         if(restrict_ntuples .and. W > pnlim%maxW .and. .not.phconj(1) .and. .not.phconj(2) )cycle  ! KEY STEP 
                ! THROW OUT COUPLED PAIRS THAT EXCEED LIMIT ON W 
         ncount = ncount + 1 
         if(create)then 
			
            indx = (a-1)*numorb(2) + b 
            PNcouplemap(indx) = ncount 
            PNcouples%pairc(ncount)%ia = a 
            PNcouples%pairc(ncount)%ib = b 
            PNcouples%pairc(ncount)%indx = indx 
            par = orbqn(1,a)%par*orbqn(2,b)%par 
            par = (3-par)/2 
            PNcouples%pairc(ncount)%par = par 
         end if 
 
     end do ! b 
 
  end do ! a 
 
  if(.not. create)then 
    ncouplesPN = ncount 
!    print*,ncouplesPN,' coupled PN pairs kept ',it 
  end if 
 
!---------------- SORT ! on parity only ------------------- 
  if(create)then 
     if(allsameparity)then 
         PNcouples%mestart(1) = 0 
         PNcouples%meref(1) = 0 
         PNcouples%mestart(2) =-1 
         PNcouples%meref(2) =ncouplesPN 
     else 
!............ CLEVERSORT on parity............ 
!             FIRST count up # of even parities 
         neven = 0 
         do i = 1,ncouplesPN 
            if( PNcouples%pairc(i)%par ==1)neven = neven+1 
         end do 
         if(neven == 0 .or. neven == ncouplesPN )then 
           if(neven == ncouplesPN )then 
              PNcouples%mestart(1) = 0 
              PNcouples%meref(1) = 0 
              PNcouples%mestart(2) =-1 
              PNcouples%meref(2) =ncouplesPN 
           else 
              PNcouples%mestart(1) =-1 
              PNcouples%meref(1) = 0 
              PNcouples%mestart(2) =0 
              PNcouples%meref(2) =0 
           endif 
         else 
            ieven = 0 
            i = 1 
            j = ncouplesPN 
            do while( ieven /= neven) 
               if( PNcouples%pairc(i)%par == 2)then 
                   do while( PNcouples%pairc(j)%par ==2 ) 
                       j = j -1 
                       if(j <= i)then 
                           print*,' whoops missed in pn ' 
                           stop 
                       endif 
                   end do 
!----------------------- SWAP 
                   cptemp = PNcouples%pairc(i) 
                   PNcouples%pairc(i) = PNcouples%pairc(j) 
                   PNcouples%pairc(j) = cptemp 
                   PNcouplemap(  PNcouples%pairc(i)%indx )  =i 
                   PNcouplemap(  PNcouples%pairc(j)%indx )  =j 
                    
               endif 
               ieven = ieven + 1 
               i = i + 1 
            end do 
!..................... FIND START, FINISH................ 
            PNcouples%mestart(1) = 0 
            PNcouples%meref(1) = 0    
            PNcouples%mestart(2) = neven*(neven+1)/2 
            PNcouples%meref(2) = neven            
 
         endif 
     endif 
 
  endif 
 
 
  return 
end subroutine setup_tbme_PNcouples 
 
!============================================================ 
! 
! set up storage of coupled PN TBMEs 
!  
!  started 7/2010 by CWJ @SDSU 
! 
 
subroutine set_tbme_PNarray 
 
  use sporbit 
  use coupledmatrixelements 
  use pairdefcoupled 
  use butil_mod
  implicit none 
 
  integer cpair,dpair,ccpair,ddpair,pair1,pair2 
  integer dparity 
  integer iref,istart 
  integer itbme 
  integer a,b,c,d 
  integer ja,jb,jc,jd 
  integer jmin,jmax 
  integer nv 
  integer ncpair 
  integer :: aerr 
 
!.................... COMPUTE # OF MATRIX ELEMENTS 
  nv = 0 
  ncpair = PNcouples%meref(2) 
  nv = ncpair*(ncpair+1)/2 
  ncpair = ncouplesPN - PNcouples%meref(2) 
  nv = nv + ncpair*(ncpair+1)/2 
  nv2bmedim(0) = nv 
  if(nv2bmedim(0) == 0)then 
     return 
  else 
     allocate( pnme( nv2bmedim(0) ) , stat=aerr) 
     if(aerr /= 0) call memerror("set_tbme_PNarray 1") 
  endif 
 
  do dpair = 1,ncouplesPN 
     dparity = PNcouples%pairc(dpair)%par 
     iref = PNcouples%meref(dparity) 
     istart = PNcouples%mestart(dparity) 
     c = PNcouples%pairc(dpair)%ia 
     d = PNcouples%pairc(dpair)%ib 
     jc = orbqn(1,c)%j 
     jd = orbqn(2,d)%j 
 
     do cpair = 1, dpair 
          if( PNcouples%pairc(cpair)%par /= dparity )cycle 
          a = PNcouples%pairc(cpair)%ia 
          b = PNcouples%pairc(cpair)%ib 
          ja = orbqn(1,a)%j 
          jb = orbqn(2,b)%j 
          jmin = MAX( abs(ja-jb), abs(jc-jd) )/2 
          jmax = MIN(ja+jb,jc+jd)/2 
          itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref 
          pnme(itbme)%jmin = jmin 
          pnme(itbme)%jmax = jmax 
          if( jmax >= jmin)then 
             allocate( pnme(itbme)%v(jmin:jmax) , stat=aerr) 
             if(aerr /= 0) call memerror("set_tbme_PNarray 2") 
             pnme(itbme)%v(:) = 0.0 
          endif 
 
     end do  ! cpair 
 
  end do ! dpair 
 
  return 
end subroutine set_tbme_PNarray 
 
!============================================================ 
 
 
! 
! provides approximate count of elements stored for coupled TBMEs 
! note: because stored as derived type, this dramatically underestimates 
! 
subroutine vtbmememory(ncount) 
 
use coupledmatrixelements 
implicit none 
integer ncount 
 
integer it,i 
type (vjs), pointer :: vme(:) 
 
ncount = 0 
do it = 0,2 
   select case (it) 
 
   case(0) 
       vme => pnme     
   case(1) 
       vme => ppme     
   case(2) 
       vme => nnme     
   end select 
   do i = 1,nv2bmedim(it) 
       if( vme(i)%jmax >= vme(i)%jmin)then 
           ncount = ncount+ vme(i)%jmax - vme(I)%jmin +3 
       else 
           ncount = ncount + 2 
       endif 
 
 
   end do ! i 
 
end do  ! it 
 
return 
 
end subroutine vtbmememory 
 
!============================================================ 

subroutine broadcastTBMEs 
   
  use nodeinfo 
  use sporbit 
  use coupledmatrixelements 
  use bmpi_mod 
  implicit none 
 
  integer(4) :: ierr 
  integer :: i, j, ii, num 
  integer :: num_me 
  integer, allocatable :: indx1(:), indx2(:) 
  real, allocatable  :: xxme(:) 
  integer :: aerr 
 
!------------ DISTRIBUTE TWO-BODY MATRIX ELEMENTS -------------------------- 
     if ( nproc > 1 ) then 
! Broadcast single particle energies........................................ 
        if ( .not. allocated( pspe ) ) then 
           allocate( pspe(numorb(1)) , stat=aerr) 
           if(aerr /= 0) call memerror("broadcastTBMEs 1") 
!		   allocate( psppot(numorb(1),numorb(1)),stat=aerr)
!           if(aerr /= 0) call memerror("broadcastTBMEs 1b") 
		   
        end if 
        if ( .not. allocated( nspe ) ) then 
           allocate( nspe(numorb(2)) , stat=aerr) 
           if(aerr /= 0) call memerror("broadcastTBMEs 2") 
!		   allocate( nsppot(numorb(2),numorb(2)),stat=aerr)
!           if(aerr /= 0) call memerror("broadcastTBMEs 2b") 
        end if 
!        call BMPI_BARRIER(icomm,ierr) 
 
        call BMPI_BCAST(pspe(:),size(pspe),0,icomm,ierr) 
        call BMPI_BCAST(nspe(:),size(nspe),0,icomm,ierr) 
!        call BMPI_BCAST(psppot(:,:),size(psppot),0,icomm,ierr) 
!        call BMPI_BCAST(nsppot(:,:),size(nsppot),0,icomm,ierr) 
		
! Broadcast 2-body matrix elements.......................................... 
! PP........................................................................ 
        num_me = 0 
        if(iproc == 0)then 
           do i = 1,  nv2bmedim(1) 
              do j = ppme(i)%jmin, ppme(i)%jmax 
                 num_me = num_me + 1 
              end do 
           end do               
        end if 
 
        call BMPI_Bcast(num_me, 1,  0, icomm, ierr) 
 
        allocate(indx1(num_me), indx2(num_me), xxme(num_me), stat=aerr) 
        if(aerr /= 0) call memerror("broadcastTBMEs 3") 
 
        if(iproc == 0) then 
           ii = 0 
           do i = 1, nv2bmedim(1) 
              do j = ppme(i)%jmin, ppme(i)%jmax 
                 ii = ii + 1 
                 indx1(ii) = i 
                 indx2(ii) = j 
                 xxme(ii) = ppme(i)%v(j) 
              end do 
           end do               
        end if 
 
        call BMPI_Bcast(indx1, num_me,  0, icomm, ierr) 
        call BMPI_Bcast(indx2, num_me,  0, icomm, ierr) 
        call BMPI_Bcast(xxme, num_me,  0, icomm, ierr) 
 
        if(iproc > 0)then 
           do ii = 1, num_me 
              i = indx1(ii) 
              j = indx2(ii) 
              ppme(i)%v(j) = xxme(ii) 
           end do 
        end if 
 
        deallocate(indx1, indx2, xxme) 
 
!        do i = 1, nv2bmedim(1) 
!           num = ppme(i)%jmax - ppme(i)%jmin + 1 
!           if(num >= 1)call BMPI_BCAST(ppme(i)%v(ppme(i)%jmin),num,0,icomm,ierr) 
!        end do 
 
 
!        call BMPI_BCAST(n_pp,1,0,icomm,ierr) 
!        if(n_pp > 0)then 
!        call block_bcast_Int(n_pp,pp_indx_J(1,1:n_pp)) 
!        call block_bcast_Int(n_pp,pp_indx_J(2,1:n_pp)) 
!        call block_bcast_Re(n_pp,pp_bcast(1:n_pp)) 
!        call BMPI_BARRIER(icomm,ierr) 
!        if ( iproc /= 0 ) then 
!           do i = 1, n_pp 
!              ppme(pp_indx_J(1,i))%v(pp_indx_J(2,i)) = pp_bcast(i) 
!           end do 
!        end if 
!        end if 
! NN........................................................................ 
        num_me = 0 
        if(iproc == 0)then 
           do i = 1,  nv2bmedim(2) 
              do j = nnme(i)%jmin, nnme(i)%jmax 
                 num_me = num_me + 1 
              end do 
           end do               
        end if 
 
        call BMPI_Bcast(num_me, 1, 0, icomm, ierr) 
 
        allocate(indx1(num_me), indx2(num_me), xxme(num_me), stat=aerr) 
        if(aerr /= 0) call memerror("broadcastTBMEs 4") 
 
        if(iproc == 0) then 
           ii = 0 
           do i = 1, nv2bmedim(2) 
              do j = nnme(i)%jmin, nnme(i)%jmax 
                 ii = ii + 1 
                 indx1(ii) = i 
                 indx2(ii) = j 
                 xxme(ii) = nnme(i)%v(j) 
              end do 
           end do               
        end if 
 
        call BMPI_Bcast(indx1, num_me,  0, icomm, ierr) 
        call BMPI_Bcast(indx2, num_me,  0, icomm, ierr) 
        call BMPI_Bcast(xxme, num_me,  0, icomm, ierr) 
 
        if(iproc > 0)then 
           do ii = 1, num_me 
              i = indx1(ii) 
              j = indx2(ii) 
              nnme(i)%v(j) = xxme(ii) 
           end do 
        end if 
 
        deallocate(indx1, indx2, xxme) 
 
!PN............................................................. 
        num_me = 0 
        if(iproc == 0)then 
           do i = 1,  nv2bmedim(0) 
              do j = pnme(i)%jmin, pnme(i)%jmax 
                 num_me = num_me + 1 
              end do 
           end do               
        end if 
 
        call BMPI_Bcast(num_me, 1,  0, icomm, ierr) 
 
        allocate(indx1(num_me), indx2(num_me), xxme(num_me), stat=aerr) 
        if(aerr /= 0) call memerror("broadcastTBMEs 5") 
 
        if(iproc == 0) then 
           ii = 0 
           do i = 1, nv2bmedim(0) 
              do j = pnme(i)%jmin, pnme(i)%jmax 
                 ii = ii + 1 
                 indx1(ii) = i 
                 indx2(ii) = j 
                 xxme(ii) = pnme(i)%v(j) 
              end do 
           end do               
        end if 
 
        call BMPI_Bcast(indx1, num_me, 0, icomm, ierr) 
        call BMPI_Bcast(indx2, num_me, 0, icomm, ierr) 
        call BMPI_Bcast(xxme, num_me, 0, icomm, ierr) 
 
        if(iproc > 0)then 
           do ii = 1, num_me 
              i = indx1(ii) 
              j = indx2(ii) 
              pnme(i)%v(j) = xxme(ii) 
           end do 
        end if 
 
        deallocate(indx1, indx2, xxme) 
!....... BROADCAST sppot.....	
        num_me = (numorb(1))**2+numorb(2)**2

       allocate(xxme(num_me), stat=aerr) 
       if(aerr /= 0) call memerror("broadcastTBMEs 6") 
       if(iproc == 0) then 
          ii = 0 
          do i = 1, numorb(1)
             do j = 1, numorb(1)
                ii = ii + 1 

                xxme(ii) = psppot(i,j)
             end do 
          end do               
          do i = 1, numorb(2)
             do j = 1, numorb(2)
                ii = ii + 1 

                xxme(ii) = nsppot(i,j)
             end do 
          end do   		  
       end if 
       call BMPI_Bcast(xxme, num_me, 0, icomm, ierr) 
       if(iproc > 0) then 
          ii = 0 
          do i = 1, numorb(1)
             do j = 1, numorb(1)
                ii = ii + 1 

                psppot(i,j)=xxme(ii) 
             end do 
          end do               
          do i = 1, numorb(2)
             do j = 1, numorb(2)
                ii = ii + 1 

                nsppot(i,j)=xxme(ii) 
             end do 
          end do   		  
       end if 	   
       deallocate(xxme) 
 
     end if 
 
  return 
end subroutine broadcastTBMEs 
!=============================================== 
! ADDED IN 7.8.2  Jan 2018
! looks for information in header on format
! 
!  INPUT
!    dummyline(:)  : array of characters from header
!
!  OUTPUT
!    foundit: if TRUE, then found a formattype
!    formattype = 'def','mfd','xpn', or 'upn'

subroutine fileformatcheck(dummyline,foundit,formattype)
	implicit none
	character(4),intent(in) :: dummyline
	logical,intent(out)  :: foundit
	character(3),intent(out) :: formattype
	
	foundit=.false.
	
	select case (dummyline(1:4))
	
	case ('!DEF','#DEF','!def','#def')
	
	foundit=.true.
	formattype='def'
	
	case ('!XPN','#XPN','!xpn','#xpn')
	
	foundit=.true.
	formattype='xpn'
	
	case ('!UPN','#UPN','!upn','#upn')
	
	foundit=.true.
	formattype='upn'
	
	case ('!MFD','#MFD','!mfd','#mfd')
	
	foundit=.true.
	formattype='mfd'
	
	case default
	
	foundit = .false.
	
    end select
	
	return
end subroutine fileformatcheck

!=============================================== 
! ADDED IN 7.8.2  Jan 2018
! CAPABILITY TO READ IN OFF DIAGONAL ONE-BODY SCALAR POTENTIAL
!
!  INPUT:
!    formatchar = 'iso' to read in proton, neutron equally
!               = 'pns' to read in proton, neutron separately in two-column format
!    filenumber = unit number of the file to be read in
!    spscalep,spscalen = scaling for proton, neutron
!    isformatted = logical, if true read in as formatted file, else unformatted
!
!  OUTPUT  
!    foundit = .true. if s.p. potential read,
!              .false. if none

subroutine readinsppot(formatchar,filenumber,spscalep,spscalen,isformatted,foundit)
	use nodeinfo
    use sporbit 
    use coupledmatrixelements 
	implicit none
	character(3),intent(in) :: formatchar
	integer,intent(in) :: filenumber
	real(4) :: spscalep,spscalen
	logical,intent(in)  :: isformatted
	logical,intent(out) ::    foundit
	character   :: xchar
	logical   :: success
	
	integer :: nme
	integer :: i,a,b
	real(4) ::xpot,ypot
	integer :: neutron_offset
		
	foundit =.false.
	
	if(iproc/=0)return
	
	if(.not.isformatted)then
		print*,' Cannot read unformatted single particle potential at this time '
		return
	end if
	
    neutron_offset = numorb(1)
	
	
!
!  READ PAST ANY HEADER LINES	

    success = .false.
	do while (.not. success)
       read(filenumber,'(a)',end=111,err=111)xchar
	   if(xchar =='!' .or. xchar == '#')cycle
	   success = .true.
	   backspace(filenumber)
   end do
   read(filenumber,*,err=111,end=111)nme
   print*,nme,' nonzero single particle scalar potential elements '
   
   if(nme==0)return   

   do i = 1,nme
	   select case (formatchar)
	   
	   case ('iso')
	       read(filenumber,*,err=112,end=112)a,b,xpot
		   psppot(a,b)=psppot(a,b)+xpot*spscalep
		   nsppot(a,b)=nsppot(a,b)+xpot*spscalen
		   if(a/=b)then
			   psppot(b,a)=psppot(b,a)+xpot*spscalep
			   nsppot(b,a)=nsppot(b,a)+xpot*spscalen			   
		   end if
		   
	   case ('xpn')
		       read(filenumber,*,err=112,end=112)a,b,xpot
			   if(a <= neutron_offset .and. b <= neutron_offset)then
			      psppot(a,b)=psppot(a,b)+xpot*spscalep
   			      if(a/=b)then
   				   psppot(b,a)=psppot(b,a)+xpot*spscalep
			      end if
   			   end if
			   if(a > neutron_offset .and. b > neutron_offset)then
			      nsppot(a-neutron_offset,b-neutron_offset)=nsppot(a-neutron_offset,b-neutron_offset)+xpot*spscalen
   			      if(a/=b)then
   				   nsppot(b-neutron_offset,a-neutron_offset)=nsppot(b-neutron_offset,a-neutron_offset)+xpot*spscalen
			      end if
   			   end if
			   
			   if(a <= neutron_offset .and. b > neutron_offset .or.  & 
			    b <= neutron_offset .and. a > neutron_offset)then
				print*,' Charge changing? check the s.p. potential at end of file ',a,b
			end if


		case('pns')
		   read(filenumber,*,err=112,end=112)a,b,xpot,ypot
		   psppot(a,b)=psppot(a,b)+xpot*spscalep
		   nsppot(a,b)=nsppot(a,b)+ypot*spscalen
		   if(a/=b)then
			   psppot(b,a)=psppot(b,a)+xpot*spscalep
			   nsppot(b,a)=nsppot(b,a)+ypot*spscalen			   
		   end if
	   end select
	   
   end do
   foundit = .true.
     
    return
	
111 continue
    return
	
112 continue
    print*,' WARNING, abrupt end of file'	
		

end subroutine readinsppot


