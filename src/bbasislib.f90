!...........................................................................
!
!  basis dimension
!  and information on how to find basis label from proton, neutron labels
!
!  if proton Slater determinant is ip and neutron Slater determinant is in
!  basis = pstart(ip)+nstart(in)
!
!  IMPORTANT:  ip and in must come from "conjugate" sectors
!  also note: neutrons form inner loop
!
module basis
  use precisions

  implicit none
  integer(kind=basis_prec) :: dimbasis 		    ! dimension of combined basis
  integer(kind=basis_prec) :: dimbasis_proc         ! local dimension of combined basis
  integer(kind=basis_prec) :: nXsd(2)         ! dimensions of proton, neutrons SDs
  integer(kind=basis_prec),allocatable,target :: pstart(:),nstart(:)   ! used for finding basi
  
contains

!----
!===========================================================================
!  BIGSTICK  configuration-interaction shell-model code
!
!===========================================================================
!
!  BBASISLIB1.f
!
!===========================================================================
!  SUBROUTINES
!     basismaster: master routine
!
!  FUNCTIONS
!      w_of_loc:  sums W for array loc (locations of particles in s.p. sections)
!      w_of_occ:  sums W for array occ (# of particles occupying each s.p. section)
!===========================================================================
!   NOTES: Main goal of these routines is to construct templates
!   structure: template(ith)%of_nh(nh)%list(i)%occ(is)
!      ith = species/hand;  nh = # of particles in a haiku
!      i = which template    occ(is) = template (occupation of is'th section)
!   Templates guide construction of haikus
!  Construction of templates is complicated. There are two paths to creating 
!  templates:
!  PATH 1: For each nh (# of particles in a haiku) there is a maximum W-cut
!  dictated by the truncation of the basis. These are generated via 
!   "temporary templates" which are not occupations but via ordered locations of
!   each particle: the ip'th particle is assigned a location loc(ip)
!   and one enforces loc(ip) >= loc(ip+1), while also enforcing a max summed W
!   and a max occupation in each s.p. section. Up to the max W, all possible 
!   templates are generated.
!   PATH 2: Because of the "hop" algorithm we use, we want to make sure that for 
!   every haiku, we can destroy one particle and also find that resulting haiku
!   in our list. Towards this end we generate additional templates by 
!   deleting particles from template occupation and keeping them even if they 
!   violate the W-cut (if they satisfy the W-cut they've already been created)
!
!   The templates are, for each nh,ordered by W, which will help with ordering 
!   of haikus later on.
!
!   Finally for each template a "parent" template is identified, by deleting 
!   a particle. Later on, in haiku generation, the daughter is created by 
!   adding one particle to the parent.  The parent is not necessarily unique, 
!   but does not need to be. 
!
!   ONCE THE PARENTS ARE IDENTIFIED, one no longer needs the explicit templates!
!
!=====================================================================
!
!  MASTER ROUTINE FOR GENERATING BASIS
!  
!  CALLED BY: bigstick_main program
!
!  SUBROUTINES CALLED:
!	master_sections 	: sets up sections in single-particle space (like quantum numbers)
!	writeouthsps 		: (optional)
!       grouper			: organizes hsps by "groups" with similar quantum numbers
!	writeoutgroups		: (optional)
!  	templatemaster		: templates are used when truncating on W
!	template_heritage
!	haiku_master		: creates haikus
!	blockmaster		: organizes haikus into blocks by quantum number
!	printblockhaikus	: (optional)
!  	mastersdlimits		: finds limits on proton, neutron Slater determinants
!	mastercombinedbasis	: combines proton, neutron states into final basis
!
subroutine basismaster
  use verbosity
  use menu_choices
  use nodeinfo
  use templates
  use spsections
  use haikus
  use blocks
  use slaterlimits

  implicit none

!----------- reorder hspsqn------------------------------------
  if(iproc == 0)then
     print*,' '
     print*,' .... Building basis ... '
     print*,' '
     print*,' Information about basis: '
  end if
  call master_sections
  if(print_hsps)call writeouthsps(1)
  if(print_hsps)call writeouthsps(-1)
!-----------set up groups of s.p. states----------------------------
  call grouper
  if(print_groups)call writeoutgroups(1)
  if(print_groups)call writeoutgroups(-1)
!------------------TEMPLATES-----------------------------------------
! templates are used in construction of the basis

  call templatemaster(1)         ! for right protons haikus
  call template_heritage(1)
  call templatemaster(-1)        ! for left protons haiksu
  call template_heritage(-1)
  call templatemaster(2)         ! for right neutron haikus
  call template_heritage(2)
  call templatemaster(-2)        !     left neutron haikus
  call template_heritage(-2)
!---------------- FIND JZ LIMITS OF TEMPLATES ---------------------
!  call master_of_jzsectplates()    
!  call assign_jz_to_templates()
!  call reconcile_templates_for_Jz()
!  call reconcile_templates_for_Jz()   ! CALL TWICE TO MAKE SURE
  
!------------------HAIKUS---------------------------------------------
 
  call haiku_master(1)           !   right proton haikus
  call haiku_master(-1)          !   left proton haikus
  call haiku_master(2)           !    right neutron haikus
  call haiku_master(-2)          !   left neutron haikus

!------------------HAIKU BLOCKS---------------------------------
  call blockmaster(1)            ! right proton haiku blocks
  call blockmaster(-1)           ! left proton haiku blocks
  call blockmaster(2)           ! right neutron haiku blocks
  call blockmaster(-2)           ! left neutron haiku blocks
  if ( print_haikus ) then
     print*,' Printing out haikus...'
     call printblockhaikus(1)
     call printblockhaikus(-1)
     call printblockhaikus(2)
     call printblockhaikus(-2)
  end if

!------------------COMBINED BASIS-----------------------------
  call mastersdlimits        ! find limits on SD quantum numbers
  call mastercombinedbasis

    if(iproc == 0)then
     print*,' '
     print*,' .... Basis built ... '
     print*,' '
  end if
  return

end subroutine basismaster
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK configuration-interaction shell-model code
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  bbasislib6.f90
!
!  sixth set of routines to generate the basis:
!  sets up for sectors of SDs
!
!===============================================================

!  In the previous routines we deduced the limits of q#s for quantum numbers
!  for proton/neutron SDs. Next we need to further refine these limits
!  to get the total basis 
!

!##############################################################################
! subroutine mastercombinedbasis
!
! SUBROUTINES CALLED:
!   count_create_sectors
!   count_fill_sectors_w_blocks
!   match_conjugate_sectors
!   combine_basis
!
subroutine mastercombinedbasis

  use sectors
  use verbosity
  use blocks
  implicit none
  integer it,is
  do it = 1,2
     call count_create_sectors(it,.false.)
     call count_create_sectors(it,.true.)
     do is = 1,nsectors(it)
        call count_fill_sectors_w_blocks(it,is,.false.)
        call count_fill_sectors_w_blocks(it,is,.true.)
     enddo ! is
  enddo  ! it
  if(verbose_haikus)call countusedhaikus
  do it = 1,2
     do is = 1,nsectors(it)
        call match_conjugate_sectors(it,is,.false.)
        call match_conjugate_sectors(it,is,.true.)
     enddo ! is
  enddo  ! it
  call combinebasis
!call checkbasis
  return
end subroutine mastercombinedbasis

!##############################################################################
! subroutine count_create_sections
!
! note, proton and neutrons sections have opposite order for jz
!
! CALLED BY: mastercombinedbasis
!

subroutine count_create_sectors(it,create)

  use system_parameters
  use sectors
  use slaterlimits
  use nodeinfo
  use bmpi_mod
  use io
  use sporbit
  implicit none

  integer(4) :: ierr
  integer :: it
  logical :: create
  integer :: icount
  integer :: jzx
  integer :: jzstart,jzend,jzstep
  integer :: parstart,parend,parstep,parx
  integer :: iw
  integer :: aerr

  if(it < 1)then  ! ERROR TRAP
     print*,' Doh! variable it must be = 1 or 2 ',it
     stop
  endif

  icount = 0

  if(it== 1)then  ! proton parity + first, then -
     parstart = limxsd(1)%parmin
     parend   = limxsd(1)%parmax
     parstep  = 1
  else
     if(iparity ==1)then   ! have neutron parity + first, then -
        parstart = limxsd(2)%parmin
        parend   = limxsd(2)%parmax
        parstep  = 1
     else                    ! neutron parities in retrograde order
        parstart = limxsd(2)%parmax
        parend   = limxsd(2)%parmin
        parstep  = -1

     endif
  endif

  do parx = parstart,parend,parstep

     if(it == 1)then
        jzstart = limxsd(1)%par(parx)%jzmin
        jzend =   limxsd(1)%par(parx)%jzmax
        jzstep = 2
     else
        jzstart = limxsd(2)%par(parx)%jzmax
        jzend =   limxsd(2)%par(parx)%jzmin
        jzstep = -2
     endif
     do jzx = jzstart,jzend,jzstep
        if(limxsd(it)%par(parx)%jz(jzx)%nw==0)cycle
        do iw = 1,limxsd(it)%par(parx)%jz(jzx)%nw

           icount = icount +1

           if(create)then
              xsd(it)%sector(icount)%jzX = jzx
!--------- fill in other quantum numbers
              xsd(it)%sector(icount)%parx = parx
              xsd(it)%sector(icount)%wx   = limxsd(it)%par(parx)%jz(jzx)%wlist(iw)
           endif
        enddo ! iw
     enddo  ! jzx
  enddo  ! par

  if( .not. create ) then
     nsectors(it) = icount
     if ( iproc == 0 ) print*,'there are ',icount,' sectors ',it
	  if(icount ==0)then
		  if(iproc==0)then
			  print*,' WARNING: No sectors for species ',it
			  write(logfile,*)' WARNING: No sectors for species ',it
			  if(allsameW .and. allsameparity)then
				  print*,' Most likely Jz is set too large'
				  write(logfile,*)' Most likely Jz is set too large'
			  end if
           if(.not.allsameW .and. .not. allsameparity .and. abs(jz) < 2)then
				  print*,' Likely possibility: bad choice of parity with truncation '			  
				  write(logfile,*)' Likely possibility: bad choice of parity with truncation '			  
			  end if

		  end if
		  stop
	  end if
     if( associated( xsd(it)%sector) ) nullify(xsd(it)%sector)
     allocate(xsd(it)%sector(icount), stat=aerr)
     if(aerr /= 0) call memerror("count_create_sectors")
  endif

  return
end subroutine count_create_sectors

!##############################################################################
! subroutine count_fill_sectors_w_blocks
!
!  finds which haiku blocks can make up a sector
!
! CALLED BY: mastercombinedbasis
!

subroutine count_fill_sectors_w_blocks(it,is,fill)

use system_parameters
use blocks
use sectors
use slaterlimits
use sporbit

implicit none
integer it  ! label of species
integer is  ! label of sector
logical fill
integer N
integer nh,nhc
integer jzx,jzh,jzhc
integer parx, parh,parhc
integer parhstart,parhend
integer parhcstart,parhcend
integer icount
integer rhblock,lhblock
integer blocksize
integer(8) :: xsdstart
integer wx
integer iwh,iwhc,wh,whc
logical okay
integer :: aerr

N = np(it)
jzx = xsd(it)%sector(is)%jzx
parx = xsd(it)%sector(is)%parx
wx   = xsd(it)%sector(is)%wx

icount = 0
if(is ==1)then
   xsdstart = 1
else
   xsdstart = xsd(it)%sector(is-1)%xsdstart+xsd(it)%sector(is-1)%nxsd
endif
if(fill)xsd(it)%sector(is)%blockstart(1) = xsdstart
do nh = blockmap(it)%nhmin,blockmap(it)%nhmax
   nhc = N -nh
   if(nhc < blockmap(-it)%nhmin .or. nhc > blockmap(-it)%nhmax)cycle
   do parh = blockmap(it)%nh(nh)%parmin,blockmap(it)%nh(nh)%parmax
     parhc = parmult(parx,parh)
     if( parhc > blockmap(-it)%nh(nhc)%parmax .or. parhc < blockmap(-it)%nh(nhc)%parmin)cycle

     do jzh = blockmap(it)%nh(nh)%par(parh)%jzmin, blockmap(it)%nh(nh)%par(parh)%jzmax,2
        jzhc = jzx - jzh
        if( jzhc > blockmap(-it)%nh(nhc)%par(parhc)%jzmax .or. jzhc < blockmap(-it)%nh(nhc)%par(parhc)%jzmin)cycle

        do iwh = 1,blockmap(it)%nh(nh)%par(parh)%jz(jzh)%nw
          wh = blockmap(it)%nh(nh)%par(parh)%jz(jzh)%wlist(iwh)
          iwhc = 1
          okay = .false.
          whc = blockmap(-it)%nh(nhc)%par(parhc)%jz(jzhc)%wlist(iwhc)
          if(wh+whc == wx)okay = .true.
          do while(.not.okay .and. iwhc < blockmap(-it)%nh(nhc)%par(parhc)%jz(jzhc)%nw)
            iwhc = iwhc +1
            whc = blockmap(-it)%nh(nhc)%par(parhc)%jz(jzhc)%wlist(iwhc)
            if(wh+whc == wx)okay = .true.

          enddo
          if(.not.okay)cycle
          icount = icount+1 
          if(fill)then
              rhblock = blockmap(it)%nh(nh)%par(parh)%jz(jzh)%hblock(iwh) 
              lhblock = blockmap(-it)%nh(nhc)%par(parhc)%jz(jzhc)%hblock(iwhc) 

              if(lhblock == 0 .or. rhblock==0)then
                print*,rhblock,lhblock
                print*,iwh,iwhc
                print*,blockmap(-it)%nh(nhc)%par(parhc)%jz(jzhc)%hblock
                print*,blockmap(-it)%nh(nhc)%par(parhc)%jz(jzhc)%nw

                stop
              endif
              xsd(it)%sector(is)%rhblock(icount) = rhblock
              xsd(it)%sector(is)%lhblock(icount) = lhblock

              blocksize = hblock(it)%list(rhblock)%nhsd * hblock(-it)%list(lhblock)%nhsd
              if(icount == 1)then
                 xsd(it)%sector(is)%blockstart(1) = xsdstart
              else
                 xsd(it)%sector(is)%blockstart(icount)= xsd(it)%sector(is)%blockend(icount-1)+1
              endif

              xsd(it)%sector(is)%blockend(icount) = xsd(it)%sector(is)%blockstart(icount)+blocksize - 1
			  
!...... CHECK FLAG THAT THESE BLOCKS ARE BEING USED.............
              hblock(it)%list(rhblock)%used=.true.
              hblock(-it)%list(lhblock)%used=.true.
			  			  
            endif
        enddo  ! iwh
     enddo ! jzh
  enddo ! parh
enddo  ! nh
if(.not.fill)then
  xsd(it)%sector(is)%nhblocks = icount
  if (associated( xsd(it)%sector(is)%rhblock ) ) nullify( xsd(it)%sector(is)%rhblock) 
  if (associated( xsd(it)%sector(is)%lhblock ) ) nullify( xsd(it)%sector(is)%lhblock) 
  allocate(xsd(it)%sector(is)%rhblock(icount),xsd(it)%sector(is)%lhblock(icount), stat=aerr)
  if(aerr /= 0) call memerror("count_fill_sectors_w_blocks 1")

  if (associated( xsd(it)%sector(is)%blockstart))nullify( xsd(it)%sector(is)%blockstart) 
  if (associated( xsd(it)%sector(is)%blockend))nullify( xsd(it)%sector(is)%blockend) 
  allocate(xsd(it)%sector(is)%blockstart(icount),xsd(it)%sector(is)%blockend(icount), stat=aerr)
  if(aerr /= 0) call memerror("count_fill_sectors_w_blocks 2")
endif
if(fill)then
   xsd(it)%sector(is)%xsdstart = xsdstart
   xsd(it)%sector(is)%nxsd = xsd(it)%sector(is)%blockend(icount) - xsdstart+1
   xsd(it)%sector(is)%xsdend = xsd(it)%sector(is)%blockend(icount) 
endif
return
end subroutine count_fill_sectors_w_blocks

!##############################################################################
!
!subroutine match_conjugate_sectors
!
!  The basis has fixed quantum numbers: Jz, parity, and max W. It is 
!  subdivided into "sectors" based upon proton and neutron quantum numbers.
!  For each proton sector (some Jzp, parp, Wp) there are a set of conjugate 
!  neutron sectors that can be coupled to make the basis, and vice versa. 
!  This routine, for a given species (it) and sector (is) finds the 
!  matching conjugate sectors.
!
!  INPUT:
!    it = species 
!    is = label of sector
!    fill = logical flag; if .true. then fill, else just count up
!
! CALLED BY: mastercombinedbasis
!
subroutine match_conjugate_sectors(it,is,fill)

use system_parameters
use W_info
use sectors
use sporbit

implicit none
integer it
integer is
logical fill

integer itc
integer isc
integer jzx,jzc
integer parx,parc
integer wx,wxc
integer icount
integer :: aerr

itc = 3 - it

jzx = xsd(it)%sector(is)%jzx
jzc = jz - jzx

parx = xsd(it)%sector(is)%parx
parc = parmult(parx, iparity)

wx =   xsd(it)%sector(is)%wx
wxc = maxWtot -wx 
icount = 0
do isc = 1,nsectors(itc)
   if(jzc == xsd(itc)%sector(isc)%jzx .and. parc == xsd(itc)%sector(isc)%parx .and. wxc >= xsd(itc)%sector(isc)%wx )then
      icount = icount + 1
      if(fill)then
        xsd(it)%sector(is)%csector(icount) = isc
      endif
   endif

enddo !isc
if(.not.fill)then
   if(icount ==0)then
     print*,' problem, no conjugates '
     print*,it,is
     print*,jzx,parx,wx,maxWtot
     print*,jzc,wxc
     stop
   endif
   xsd(it)%sector(is)%ncsectors = icount
   if( associated( xsd(it)%sector(is)%csector)) nullify( xsd(it)%sector(is)%csector)
   allocate(xsd(it)%sector(is)%csector(icount), stat=aerr)
   if(aerr /= 0) call memerror("match_conjugate_sectors")
endif
return 
end subroutine match_conjugate_sectors

!##############################################################################
!subroutine combinebasis
! key routine to generate the combined basis from proton and neutron SDs
!
! KEY RESULT is the construction of the arrays pstart, nstart,
!            where it = species and isd = which (proton/neutron) slater determinant
!            Given  a proton SD labeled by ip  and a neutron SD labeled by in
!            then pstart(ip) + nstart(in) = ibasis = label in combined basis
!
! The basis (that is, the xstart arrays) is generated by looping over 
!  the proton sectors and then over the conjugate neutron sectors
!  
!  CALLED BY: mastercombinedbasis
!
subroutine combinebasis

  use sectors
  use nodeinfo
  use bmpi_mod
  implicit none

  integer(4) :: ierr

  integer it
  integer ns
  integer nb,ib
  integer is, isc,jsc
  integer(8) :: ncount,pcount
  !integer pstart
  integer(8) :: ndim,pdim
  integer(8) :: ip,in
  integer :: aerr
!---------- # of SDs in each species
  nXsd = 0

! testing

!do it = 1,2
!  print*,' * * * '
!  do is = 1,nsectors(it)
!    print*,' sectors ',is,xsd(it)%sector(is)%nhblocks
!    do ib = 1,xsd(it)%sector(is)%nhblocks
!       print*,xsd(it)%sector(is)%blockstart(ib),xsd(it)%sector(is)%blockend(ib)
!    enddo
!  enddo
!enddo

  do it = 1,2
     ns = nsectors(it)
     nb = xsd(it)%sector(ns)%nhblocks
     nXsd(it) = xsd(it)%sector(ns)%blockend(nb)
     if ( iproc == 0 ) print*,nXsd(it),' SDs for species ',it
  end do ! it

! ---------- NOW COUPLE TOGETHER

  if( allocated( pstart )) deallocate( pstart) 
  if( allocated( nstart )) deallocate( nstart) 

  if(nxsd(1) >0) then
     allocate(pstart(nxsd(1)), stat=aerr)
     if(aerr /= 0) call memerror("combinebasis 1")
  end if
  if(nxsd(2) >0) then
     allocate(nstart(nxsd(2)), stat=aerr)
     if(aerr /= 0) call memerror("combinebasis 2")
  end if
  pstart = -1
  nstart = -1
  dimbasis = 1
  pcount = 0
  ncount = 0

  do is = 1,nsectors(1)  ! loop over proton sectors

     ndim = 0
     do isc = 1,xsd(1)%sector(is)%ncsectors  
        jsc = xsd(1)%sector(is)%csector(isc)
        ndim = ndim + xsd(2)%sector(jsc)%nxsd
     end do ! isc

!  ALTHOUGH LOOP OVER PROTONS SECTORS, FILL NEUTRON SECTORS FIRST
!  IF AND ONLY IF THE FIRST CONJUGATE NEUTRON SECTOR HAS THIS PROTON SECTOR 
!  AS FIRST CONJUGATE SECTOR 

     jsc = xsd(1)%sector(is)%csector(1)
     if( is == xsd(2)%sector(jsc)%csector(1))then  ! FILL NEUTRON

!  set up values in xstart(2,...) by looping over conjugate neutron sectors
!                            and then looping over neutron SDs in each sectors
        do in = 1,ndim
           ncount = ncount + 1
           if ( ncount > nxsd(2) ) then
              if ( iproc == 0 ) then
                 print*,' overstepped (2) ',ncount,nxsd(2)
                 print*,is,nsectors(1)
              end if
              call BMPI_ABORT(icomm,101,ierr)
              stop
           end if
           nstart(ncount) = in
        end do  !in
     end if

     pdim = xsd(1)%sector(is)%nxsd  ! get dimension of this proton sector
 ! set values in xstart(1,...) by looping over these proton SDs

     do ip = 1,pdim      
        pcount = pcount + 1

        if( pcount > nxsd(1) ) then
           if ( iproc == 0 ) print*,' overstepped (1) ',pcount,nxsd(1)
           call BMPI_ABORT(icomm,101,ierr)
           stop
        end if

        pstart(pcount) = dimbasis - 1
        dimbasis = dimbasis + ndim

     end do ! ip

  
  end do !is 
  dimbasis = dimbasis - 1
  if ( iproc == 0 ) print*,' Total basis = ',dimbasis

!--------- ERROR TRAPS

  do ip = 1,nxsd(1)
     if ( pstart(ip) < 0 ) then
        if ( iproc == 0 )print*,' error in proton basis ',ip,pstart(ip)
        call BMPI_ABORT(icomm,101,ierr)
        stop
     end if
  end do !ip
  
  do in = 1,nxsd(2)
     if ( nstart(in) < 0 ) then
        if ( iproc == 0 ) print*,' error in neutron basis ',in,nstart(in)
        call BMPI_ABORT(icomm,101,ierr)
        stop
     end if
  end do !ip
  return
end subroutine combinebasis

end module basis
!##############################################################################

!##############################################################################

!##############################################################################
!subroutine sector_stats
!
!  CALLED BY: main routine
!  (after basis is created) 
!
subroutine sector_stats(it)

use verbosity
use sectors
implicit none
integer it
integer is

integer nblocks,maxblocks

if(.not.verbose_sectors)return
nblocks = 0
maxblocks = 0
do is = 1,nsectors(it)
  nblocks = nblocks + xsd(it)%sector(is)%nhblocks
  maxblocks = max(maxblocks,nblocks + xsd(it)%sector(is)%nhblocks)
enddo  ! is

print*,' avg # of blocks per sector is ',nblocks/nsectors(it)
print*,' max # of blocks in a sector is ',maxblocks

return
end subroutine sector_stats
!##############################################################################
! subroutine checkbasis
! 
! optional subroutine used strictly to make sure basis arrays are correct
!  
!  called by main routine when reading in old wavefunction file
!
subroutine checkbasis

use basis
use sectors
implicit none

logical, allocatable :: check(:)

integer is,isc,jsc
integer(8) :: ip,in
integer :: aerr
integer(kind=basis_prec) :: ibasis

 allocate(check(dimbasis), stat=aerr)
 if(aerr /= 0) then
    call memerror("checkbasis 1")
    stop 5
 end if

 check(:) = .false.
 do is = 1,nsectors(1)
   do ip = xsd(1)%sector(is)%xsdstart,xsd(1)%sector(is)%xsdend

     do isc = 1,xsd(1)%sector(is)%ncsectors
        jsc = xsd(1)%sector(is)%csector(isc)

        do in = xsd(2)%sector(jsc)%xsdstart,xsd(2)%sector(jsc)%xsdend
          ibasis = pstart(ip)+nstart(in)
!----------------- ERROR TRAP --------------------
          if(check(ibasis))then
           print*,' state already occupied ! ',ibasis
           print*,ip,in
           print*,pstart(ip),nstart(in)
           stop
          endif

          check(ibasis) = .true.
        enddo
     enddo !isc
   enddo  !ip

 enddo  !is

 do ibasis = 1,dimbasis
   if(.not.check(ibasis))then
       print*,' missing basis state ',ibasis
       stop
   endif
 enddo
 print*,' All states present and accounted for '
 return
end subroutine checkbasis
!##############################################################################
! subroutine writeoutbasis
!
! optional subroutine to write out basis states to a file
!
! illustrates how the basis is created and implicitly stored
!
subroutine writeoutbasis

use spstate
use haiku_info
use blocks
use basis
use sectors
use verbosity
use precisions

implicit none

  character (len=15) :: filename
  character (len=1)  :: ychar
  integer ilast
  integer iunit
  integer it
  integer(8) :: i
  integer is  ! sector
  integer iblock, lblock,rblock
  integer(8) :: xsdstart, ix
  integer ladd,radd
  integer nradd,nladd
  logical, allocatable :: checkbasis(:)
  integer (kind=basis_prec), pointer :: start(:)
  integer :: aerr

  iunit = basis_file
  print*,' '
  print*,' Do you want to write out the basis (y/n)?'
  read(5,'(a)')ychar

  if(ychar=='n' .or. ychar=='N')return


  print*,' Print out basis (f) factorized (c) combined '
  do while( ychar /= 'f' .and. ychar /= 'c')
    read(5,'(a)')ychar
    if(ychar == 'F')ychar='f'
    if(ychar == 'C')ychar='c'
    if( ychar /= 'f' .and. ychar /= 'c')then
      print*,' (f) factorized means proton, neutron slaters printed out '
      print*,'     + start arrays '
      print*,' (c) combined means proton and neutron slaters side-by-side '
    endif
  enddo

!------------------- OPEN A FILE -------------------------------

  select case (ychar)

  case ('f')
    print*,' Enter output file (.basf) name '
    read(5,'(a)')filename
    ilast =  index(filename,' ') - 1
    open(unit=iunit,file=filename(1:ilast)//'.basf',status='unknown')
  case ('c')
    print*,' Enter output file (.basc) name '
    read(5,'(a)')filename
    ilast =  index(filename,' ') - 1
    open(unit=iunit,file=filename(1:ilast)//'.basc',status='unknown')

  end select

!------------------ WRITE OUT SINGLE PARTICLE STATES------------------------
!                   FOR NOW, ASSUME PROTON AND NEUTRON BASIS THE SAME 
!------------------- LEFT BASIS ----------------------------------
it = 1
do i = 1,nhsps(it)  
   write(iunit,101)hspsqn(it,i)%nr, hspsqn(it,i)%j,  hspsqn(it,i)%l, hspsqn(it,i)%m,  hspsqn(it,i)%w
101 format(5i4)
enddo
do i = 1,nhsps(-it)  
   write(iunit,101)hspsqn(-it,i)%nr, hspsqn(-it,i)%j,  hspsqn(-it,i)%l,  hspsqn(-it,i)%m,  hspsqn(-it,i)%w
enddo

!------------------- FACTORED BASIS ----------------------------------------
select case (ychar)

case ('f')

do it = 1,2
  allocate(checkbasis( nxsd(it)), stat=aerr)
  if(aerr /= 0) then
     call memerror("writeoutbasis 1")
     stop 5
  end if
  checkbasis = .false.
  if(it==1)then
    start=>pstart
  else
    start=>nstart
  endif
  do is = 1,nsectors(it) 
     do iblock = 1,xsd(it)%sector(is)%nhblocks
        rblock = xsd(it)%sector(is)%rhblock(iblock)
        lblock = xsd(it)%sector(is)%lhblock(iblock)
        xsdstart = xsd(it)%sector(is)%blockstart(iblock)-1

!---------------- always make left address the outermost loop (a convention)
        nradd = hblock(it)%list(rblock)%nhsd
        nladd = hblock(-it)%list(lblock)%nhsd

        do ladd = 1,nladd
            do radd = 1,nradd
                ix = xsdstart + radd+(ladd-1)*nradd
                if(checkbasis(ix))then
                   print*,' that state already exists '
                   print*,ix
                   stop
                endif
                checkbasis(ix) = .true.
!------------------- WRITE OUT -----------------------
write(iunit,102)ix,start(ix),(hblock(-it)%list(lblock)%hsd(i,ladd),i=1,nword(it)), & 
                            (hblock(it)%list(rblock)%hsd(i,radd),i=1,nword(it))
102 format(2i10,4x,4i7)
            enddo  ! radd
        enddo  ! ladd

     enddo ! iblock

  enddo ! is
  do i = 1,nxsd(it)
     if(.not.checkbasis(i))then
       print*,' missing a basis state '
     endif
  enddo
  deallocate(checkbasis)
enddo ! it
!--------------------COMBINED BASIS ----------------------------------------
case ('c')

end select
close(iunit)
return
end subroutine writeoutbasis

!##############################################################################
!subroutine writeoutsectors
!  optional subroutine to write out sector information
!
 subroutine writeoutsectors(it)

 use sectors
 use verbosity
 use io
 implicit none
 integer it
 integer ilast

 integer is,ib

 if(.not.print_sectors)return

 if(writeout .and. it == 1)then
    ilast = index(outfile,' ')-1
    open(unit=sector_file,file=outfile(1:ilast)//'.sectors',status='unknown')
 else
    open(unit=sector_file,file='sectors.bigstick',status='unknown')
 endif

 write(sector_file,*)' SECTORS ',it
 write(sector_file,*)nsectors(it)
 do is = 1,nsectors(it)
    write(sector_file,101)is, xsd(it)%sector(is)%jzX,xsd(it)%sector(is)%parX, xsd(it)%sector(is)%wX, &
	    xsd(it)%sector(is)%basisstart,xsd(it)%sector(is)%basisend
101 format(i4,' jz = ',i3,' par = ',i2,' W = ',i3,', basis from ',i12,' to ',i12)
       do ib = 1,xsd(it)%sector(is)%nhblocks
         write(sector_file,201)ib, xsd(it)%sector(is)%lhblock(ib), xsd(it)%sector(is)%rhblock(ib)
201 format(5x,' Block ',i3,' : ',2i4)
       enddo

 enddo ! is

 return

 end subroutine writeoutsectors

!##############################################################################
!
! added in 7.4.9 by CWJ
!
! finds start and stop in basis for each proton sector
!

subroutine sectorboundaries
	
	use sectors
	use basis
	implicit none
	
	integer :: is,ics,cs
	integer(8) :: dimsofar,npsd,nnsd
	
	dimsofar = 0   ! summary dimension
	
	do is = 1,nsectors(1)  ! loop over proton sectors
		xsd(1)%sector(is)%basisstart = dimsofar+1
		npsd = xsd(1)%sector(is)%nxsd  ! # of proton SDs in this proton sector
		nnsd = 0              
		do ics = 1,xsd(1)%sector(is)%ncsectors  ! loop over conjugate neutron sectors
			cs = xsd(1)%sector(is)%csector(ics) ! fetch index for conjugate neutron sector
			nnsd = nnsd + xsd(2)%sector(cs)%nxsd  ! count up # of neutron SDs in conjugate sectors
		end do  ! ics
		dimsofar = dimsofar + npsd*nnsd
		xsd(1)%sector(is)%basisend = dimsofar
	end do  ! is
! -----  CHECK IF DIMENSIONS AGREE ---

    if( dimsofar /= dimbasis)then
	    print*,' mismatch in dimensions ',dimsofar,dimbasis
		stop
    end if
	return
end subroutine sectorboundaries

!---
