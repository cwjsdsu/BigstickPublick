!===========================================================================
!  BIGSTICK  configuration-interaction shell-model code
!
!===========================================================================
!
!  BBASIS_TEMPLATES.f
!
!  first set of routines to generate the basis:
!  TEMPLATES
!
!===========================================================================
!  SUBROUTINES
!     templatemaster : coordinates construction of templates
!        calls: count_fill_templates
!               find_other_templates
!               fill_tmps
!     tmp_templates:  initial fill of templates based on W-cuts
!        calls: initial_template
!               count_fill_tmp
!     initial_template: fills in starting temporary templates for given nh
!     count_fill_tmp: constructs temporary templates by moving particles recursively
!        calls: w_of_loc
!     fill_tmps: fills in templates from temporary templates
!        calls: w_of_occ
!     find_other_templates: finds templates not from W-cuts, but from 
!                           deleting particles from template for nh+1
!        calls: w_of_occ
!     sort_templates: sorts by W; 
!     template_heritage: finds a "parent" for each template, to guide 
!                       recursive construction of haikus 
!     (write_templates: only used to print out templates for debugging)
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
 

module templates

  integer, allocatable :: ntemplate(:,:)

  type templatebase
     integer :: w
     integer :: parent         ! id of "parent" template
     integer :: add            ! which section in which to add
     integer, pointer :: occ(:)
     integer :: hstart,hfin    ! where to start,end counting haikus
	 integer :: jzmax,jzmin
  end type templatebase

  type templateofn
     type (templatebase), pointer :: list(:)
	 integer :: wmin,wmax
	 integer, allocatable :: jzmin_of_w(:),jzmax_of_w(:)
	 
  end type templateofn

  type templateit
     type (templateofn), pointer :: of_nh(:)
  end type templateit

  type (templateit) :: template(-2:2)
	  
!......... JZPLATES..................................................	  
!          jzplates (added in 7.7.2) adds the dimension of Jz to templates
!   Needed for large W cuts; otherwise original templates overproduce haikus
!   
!   Jzplates find min/max Jz for a given section (subset of single particle states with fixed W)

type jzplateofsect
 	                                                        ! for a given s.p. section
	integer, allocatable :: jzmax2_of_n(:), jzmin2_of_n(:)   ! 2 x max, min jz as a function of n
end type jzplateofsect

type jzplateit
	type ( jzplateofsect ), allocatable :: of_section(:)
end type jzplateit

type (jzplateit)  :: jzsectplate(-2:2)
	
!....... INTERMEDIATE USEFUL ARRAY OF JZ values in any given section

integer, allocatable :: jzarray(:)	

  contains

!===========================================================================
!
! subroutine templatemaster
!
!  INPUT:
!   ith = species/handness
!
!  SUBROUTINES CALLED:
!    tmp_templates
!     fill_tmps
!     find_other_templates
!     sort_templates
!
subroutine templatemaster(ith)
  use verbosity 
  use sporbit
  use haiku_info
  use builder
  implicit none
!  include 'binterfaces.inc'
!------------ INPUT
  integer ith        ! species, including left/right
!------------- INTERMEDIATE
  integer nh        ! # of particles in a haiku


  type (blist),pointer :: tmptemp(:)
  integer,pointer      :: nt(:)
  integer ntemp
  integer nsize
  integer nstart
  integer ih
  integer :: aerr
!------------ LOOP OVER NH DOWNWARDS 

  allocate(tmptemp(0:maxNh(ith)), stat=aerr) 
  if(aerr /= 0) call memerror("templatemaster 0")
  allocate(nt(0:maxNh(ith)), stat=aerr)
  if(aerr /= 0) call memerror("templatemaster 1")
  allocate(template(ith)%of_nh(0:maxNh(ith)), stat=aerr)
  if(aerr /= 0) call memerror("templatemaster 2")
  if ( .not. allocated( ntemplate ) ) then
     nsize = maxNh(1)
     nsize = MAX(nsize,maxNh(-1))
     nsize = MAX(nsize,maxNh(2))
     nsize = MAX(nsize,maxNh(-2))

     allocate( ntemplate(-2:2,0:nsize) , stat=aerr)
     if(aerr /= 0) call memerror("templatemaster 3")
  end if
!----------------- IF ALL SAME W THEN TRIVIAL ---------
  if ( allsameW ) then
     do nh  =  0, maxNh(ith)
        ntemplate(ith,nh) = 1
        allocate( template(ith)%of_nh(nh)%list(1) , stat=aerr)
        if(aerr /= 0) call memerror("templatemaster 10")
        allocate( template(ith)%of_nh(nh)%list(1)%occ(1) , stat=aerr)
        if(aerr /= 0) call memerror("templatemaster 11")
        template(ith)%of_nh(nh)%list(1)%occ(1) = nh
        template(ith)%of_nh(nh)%list(1)%W  =  0
     end do
     return
  end if

  if ( verbose_templates ) then
     print*,' For species ', ith, maxNh(ith)
  end if
!-----------------------------------------------------
  do nh = maxNh(ith), 1,-1
!------------ COUNT # TEMPLATES UP TO MAX W
     ntemp = 0
     nt = 0
!     if ( nh >= minNh(ith) ) then  ! otherwise there won't be any
     call tmp_templates(ith,nh,ntemp,nt,tmptemp)

!     end if
!------------ COUNT ADDITIONAL TEMPLATES DERIVED FROM NH+1
     if ( nh < maxNh(ith) ) then
        call find_other_templates(.false.,ith,nh,ntemp,0)
     end if

!------------ ALLOCATE FOR TEMPLATES ------------------
     ntemplate(ith,nh) = ntemp
     allocate( template(ith)%of_nh(nh)%list(ntemp) , stat=aerr)
     if(aerr /= 0) call memerror("templatemaster 20")
!------------ FILL IN TEMPLATES UP TO MAX W 

!     if ( nh >= minNh(ith) ) then
        call fill_tmps(ith,nh,nt,tmptemp)
!     end if

!------------- FILL IN TEMPLATES DERIVED FROM NH+1
     if ( nh < maxNh(ith) ) then
        ntemp = 0
        do ih = 0,nh
           ntemp = ntemp + nt(ih)
        end do
        nstart = ntemp + 1
        call find_other_templates(.true.,ith,nh,ntemp,nstart)
        ntemplate(ith,nh) = ntemp
     end if

     call sort_templates(ith,nh)
!--------------- OPTION TO PRINT OUT
     if ( verbose_templates ) then
        print*,  nh,' particles ', ntemp, ' templates ' 
     end if
  end do  ! nh
  deallocate( nt,tmptemp )
  return
end subroutine templatemaster

!===========================================================================
!
!      subroutine tmp_templates
!
!  recursively counts up (and fills) temporary templates for nh particles, up to wmax
!
!  intially, the templates are given by the location (array loc(1:nh)) of the 
!  particles, with loc(i) <= loc(i-1).  This allows one to generate a unique set
!  of templates. Then, recursively, the particles are moved upwards.
!
!  INPUT:
!    ith = species, handness (proton.neutron, left/right)
!    nh  = # of particles 
!
!  OUTPUT:
!   ntemp = # of templates for this nh, so far
!   nt = # tmplates distributed according to a shift
!   tmptemp = temporary array of templates
!
!  SUBROUTINES CALLED:
!     intial_template
!     count_fill_tmp
!
subroutine tmp_templates(ith,nh,ntemp,nt,tmptemp)

  use spsections
  use builder
  implicit none

!-------------- INPUT
  integer :: ith  ! species + handness
  integer :: nh   ! # of particles
!--------------OUTPUT
  integer :: ntemp ! # of templates for this nh, so far
  integer, pointer :: nt(:)
  type (blist),pointer :: tmptemp(:)
!------------ INTERMEDIATE --------
  integer ih,i
  integer :: aerr

!------------- ALLOCATE AND FILL INITIAL LOCATIONS
  allocate ( tmptemp(0)%list(1) , stat=aerr)
  if(aerr /= 0) call memerror("count_fill_tmp 1")
  allocate ( tmptemp(0)%list(1)%loc(nh) , stat=aerr)
  if(aerr /= 0) call memerror("count_fill_tmp 2")
  nt(0) = 1
  ntemp = 1
  call initial_template(ith,nh,tmptemp(0)%list(1)%loc)
!----------  LOOP OVER ----------------------
  do ih = 1, nh
     call count_fill_tmp(.false.,ith,nh,ih,nt,tmptemp)
     ntemp = ntemp + nt(ih)
     allocate ( tmptemp(ih)%list(nt(ih)) , stat=aerr)
     if(aerr /= 0) call memerror("count_fill_tmp 10")
     call count_fill_tmp(.true.,ith,nh,ih,nt,tmptemp)
  end do ! ih

!------------------ TESTING

!      print*,nh,' particles there are ',ntemp,' templates starting '
!      do ih = 0,nh
!         do i = 1,nt(ih)
!            print*,tmptemp(ih)%list(i)%loc
!         enddo
!      enddo
  return
end subroutine tmp_templates

!===========================================================================
!      subroutine initial_template
!
!  sets up initial template for nh particles of ith species/handness
!
!  INPUT:
!    ith: = species/handness
!    nh:    # of particles
!  OUTPUT:
!    loc(i) = location (section) of ith particle
!       requires  loc(i+1) <= loc(i)
!       also can only put a fixed number of particles into each section
! 
subroutine initial_template(ith,nh,loc)

  use spsections
  implicit none
!------------ INPUT
  integer ::  ith
  integer :: nh
!----------- OUTPUT
  integer, pointer :: loc(:)
!----------- INTERNAL

  integer :: i
  integer :: ncnt
  integer :: isection

  ncnt = 0
  isection = 1
  do i = nh, 1, -1
     if ( ncnt == section(ith)%list(isection)%nsps ) then
        isection = isection +1
        ncnt = 0
     end if
     ncnt = ncnt + 1
     loc(i) = isection
  end do ! nh

  return
end subroutine initial_template
!===========================================================================
!
!     subroutine count_fill_tmp
!
!  counts and fills temporary array
!  moves the ih'th particle (out of nh)
!
!  INPUT:
!    fill :logical flag; if .false. just count, else fill
!    ith = species/handness
!    nh = # of particles
!    ih = which particle being moved
!    nt(ih), ih = 0,nh = # of temporary templates for ih'th particle moved
!   INPUT/OUTPUT
!   tmp = structure that contains temporary templates
!
!  FUNCTIONS CALLED:
!    w_of_loc
!
subroutine count_fill_tmp(fill,ith,nh,ih,nt,tmp)
  
  use haiku_info
  use builder
  use spsections
  implicit none

  logical               :: fill
  integer               :: ith
  integer               :: nh
  integer               :: ih
  integer,pointer       :: nt(:)
  type (blist), pointer :: tmp(:)

  integer               :: i
  integer               :: jh, kh
  integer               :: locstart, locfin
  integer               :: occcnt
  
  integer               :: wstart
  logical               :: checking
  integer               :: icnt
  integer :: aerr
  
  icnt = 0
  do kh = 0,ih -1
     do i = 1,nt(kh)
!------------- FIND RANGE OVER WHICH I CAN MOVE THE PARTICLE ------

        locstart = tmp(kh)%list(i)%loc(ih)+1

        if(ih == 1)then
           locfin = nsections(ith)
        else
           locfin = tmp(kh)%list(i)%loc(ih-1)
!----------- CHECK NOT TOO MANY OCCUPYING FINAL STATE; IF SO, NOTCH BACK
           occcnt = 1
           do jh = 1,ih-1
             if(tmp(kh)%list(i)%loc(jh)==locfin)occcnt = occcnt+1
          enddo ! jh
          if(occcnt > section(ith)%list(locfin)%nsps)locfin=locfin-1
          
       endif

       if(locfin < locstart)cycle

!------------- FIND THE STARTING W --------

       wstart = w_of_loc(ith,nh,tmp(kh)%list(i)%loc)
       wstart = wstart - section(ith)%list(locstart-1)%w  ! subtract off moving

!------------CHECK THAT W's DON'T ADD TOO HIGH; IF SO, NOTCH BACK
       checking = .true.
       do while(checking)
          if(wstart + section(ith)%list(locfin)%w <= maxWh(ith,nh))then
             checking = .false.
          else
             locfin = locfin -1
             if(locfin < locstart)exit
          endif
       enddo
       if(locfin < locstart)cycle
       
!------------ CALCULATE NUMBER OF NEW TMPTEMPLATES GENERATED 
       if(.not.fill)then
          icnt = icnt + locfin-locstart+1
       else  ! fill ! 
          do jh = locstart,locfin
             icnt = icnt + 1
             allocate(tmp(ih)%list(icnt)%loc(nh), stat=aerr)
             if(aerr /= 0) call memerror("count_fill_tmp 1")
             tmp(ih)%list(icnt)%loc = tmp(kh)%list(i)%loc
             tmp(ih)%list(icnt)%loc(ih) = jh
          enddo ! jh
       endif

    enddo  ! i
 enddo ! kh
 if(.not.fill)nt(ih) = icnt
 return
end subroutine count_fill_tmp
!===========================================================================
!   function w_of_loc
!
!   computes the summed W for loc(:) 
!
!  INPUT:
!    ith = species/handness
!    nh =  # of particles
!    loc(1:nh) = location of each particle in the s.p. sections
!
!
integer function w_of_loc(ith,nh,loc)
  use spsections
  implicit none
  integer :: ith
  integer :: nh
  integer, pointer :: loc(:)

  integer ih
  w_of_loc = 0

  do ih = 1,nh
     w_of_loc = w_of_loc + section(ith)%list(loc(ih))%w
  enddo
  return
end function w_of_loc
!===========================================================================
!    subroutine fill_tmps
!
!    fills templates from temporary templates
!
!  INPUT:
!    ith = species, handness (proton.neutron, left/right)
!    nh  = # of particles 
!   nt = # tmplates distributed according to a shift
!   tmptemp = temporary array of templates
!
!  FUNCTIONS CALLED:
!    w_of_occ: sums W for a template occupation
!
!
subroutine fill_tmps(ith,nh,nt,tmptemp)
  
  use spsections
  use builder
  implicit none

!-------------- INPUT
  integer               :: ith  ! species + handness
  integer               :: nh   ! # of particles
  integer, pointer      :: nt(:)
  type (blist), pointer :: tmptemp(:)
!------------ INTERMEDIATE --------

  integer :: ih
  integer :: i
  integer :: ip
  integer :: icnt
  integer :: loc
  integer :: ncnt
  integer :: aerr

  icnt = 0
  do ih = 0,nh
     do i = 1,nt(ih)
        icnt = icnt + 1
        allocate(template(ith)%of_nh(nh)%list(icnt)%occ(nsections(ith)), stat=aerr)
        if(aerr /= 0) call memerror("fill_tmps")
        template(ith)%of_nh(nh)%list(icnt)%occ = 0

        do ip = 1,nh
           loc = tmptemp(ih)%list(i)%loc(ip)
           template(ith)%of_nh(nh)%list(icnt)%occ(loc) = &
                template(ith)%of_nh(nh)%list(icnt)%occ(loc)+1
!---------- ERROR TRAP
           if(loc > nsections(ith) .or. loc < 1)then
              print*,' problem with location of particle '
              print*,nh,ih,ip,loc,nsections(ith)
              stop
           endif
        enddo
        template(ith)%of_nh(nh)%list(icnt)%w =  &
             w_of_occ(ith,template(ith)%of_nh(nh)%list(icnt)%occ)
!        print*,nh,template(ith)%of_nh(nh)%list(icnt)%occ
     enddo !i
  enddo ! ih

  return
end subroutine fill_tmps
!===========================================================================
!
!    subroutine find_other_templates
!
!   finds additional templates  not generated by usual procedure
!   algorithm: loops through templates with nh+1, deletes a particle,
!   and then computes W; if greater than maxWh(nh) then add to list
!
!  INPUT:
!    fill: logical flag to fill
!    ith : species/hand
!    nh  : # of particles
!    ntemp: number of templates
!    nstart:  where to start when looking backwards
!
!  FUNCTIONS CALLED:
!    w_of_occ: sums W for a template occupation
!
subroutine find_other_templates(fill,ith,nh,ntemp,nstart)
  
  use haiku_info
  use spsections
  implicit none

  logical :: fill
  integer :: ith  ! species/hand
  integer :: nh   ! # of particles
  integer :: ntemp
  integer :: nstart
!--------------------
  integer :: i,j
  integer :: is,js
!      type (templatebase), pointer :: list(:)
  integer, pointer :: occ(:)
  integer ::wocc
!      integer w_of_occ
  logical :: unique,checker
  integer :: aerr

  allocate(occ(nsections(ith)), stat=aerr)
  if(aerr /= 0) call memerror("find_other_templates 1")
  do i = 1,ntemplate(ith,nh+1)
     occ = template(ith)%of_nh(nh+1)%list(i)%occ
     do is = 1,nsections(ith)
        occ = template(ith)%of_nh(nh+1)%list(i)%occ

        if(occ(is) > 0)then
           occ(is) = occ(is) -1
           wocc = w_of_occ(ith,occ)
           if(wocc <= maxWh(ith,nh))cycle
 
           ntemp = ntemp + 1
           if(fill)then  ! go back and check
              unique = .true.
              do j = nstart,ntemp-1
                 if(template(ith)%of_nh(nh)%list(j)%w == wocc)then  ! compare
                    checker = .true.
                    do js = 1,nsections(ith)
                       if(occ(js) /= &
                            template(ith)%of_nh(nh)%list(j)%occ(js))checker=.false.
                    enddo ! js
                    if(checker)unique = .false.
                 endif
              enddo ! j
              if(.not.unique)then
                 ntemp = ntemp -1
              else
                 allocate(template(ith)%of_nh(nh)%list(ntemp)% occ(nsections(ith)), stat=aerr)
                 if(aerr /= 0) call memerror("find_other_templates 1")
                 template(ith)%of_nh(nh)%list(ntemp)%occ=occ
                 template(ith)%of_nh(nh)%list(ntemp)%w =wocc
!                   print*,' Filled extra ',ntemp,
!     & template(ith)%of_nh(nh)%list(ntemp)%occ

              endif
           endif
        endif
     enddo !is
  enddo  ! i
  deallocate(occ)
  return
end subroutine find_other_templates
!===========================================================================
!
!  function w_of_occ
!  
!  computes summed W for a template occupation
!
!  INPUT:
!    ith:  species/hand
!    occ(:) : occupation 
!
integer function w_of_occ(ith,occ)
  
  use spsections
  implicit none
  integer          :: ith
  integer, pointer :: occ(:)
  integer          :: i

  w_of_occ = 0
  do i = 1,nsections(ith) 
     w_of_occ = w_of_occ + section(ith)%list(i)%w*occ(i)
  enddo! i

  return
end function w_of_occ

!===========================================================================
!
!  subroutine sort_templates
!
!  simple bubble sort on templates by W
!
!  INPUT:
!   ith = species/handness
!   nh = # of particles in a haiku
!
subroutine sort_templates(ith,nh)
  
  implicit none
  integer :: ith  ! species/handness
  integer :: nh
!------------------------------

  type (templatebase) :: tmp
  integer             :: i,j,k
  integer             :: w,wtmp

  do i = 1,ntemplate(ith,nh)
     k = i
     w = template(ith)%of_nh(nh)%list(i)%w
     tmp = template(ith)%of_nh(nh)%list(i)

     do j = i+1,ntemplate(ith,nh)
        wtmp = template(ith)%of_nh(nh)%list(j)%w
        if(wtmp < w)then
           k = j
           w = wtmp
           tmp = template(ith)%of_nh(nh)%list(j)
        endif
     enddo  !j
     if(k/=i)then  !swap
        template(ith)%of_nh(nh)%list(k) =  & 
             template(ith)%of_nh(nh)%list(i)
        template(ith)%of_nh(nh)%list(i) = tmp

     endif
  enddo  ! i

  return
end subroutine sort_templates

!===========================================================================
!      subroutine template_heritage(ith)
!
!  finds a "parent" of each template
!
!
subroutine template_heritage(ith)
  
  use haiku_info
  use spsections
  implicit none
  integer          :: ith
  integer          :: nh
  integer          :: i, j
  integer          :: is,js
  logical          :: match
  integer, pointer :: occ(:)
  integer :: aerr

  allocate(occ(nsections(ith)), stat=aerr)
  if(aerr /= 0) call memerror("template_heritage")
  do nh = 1,maxNh(ith)
     do i = 1,ntemplate(ith,nh)
        occ = template(ith)%of_nh(nh)%list(i)%occ
        do is = 1,nsections(ith)
           if(occ(is) > 0)then
              if(nh == 1)then
                 template(ith)%of_nh(nh)%list(i)%parent=1
                 template(ith)%of_nh(nh)%list(i)%add = is
              else
!------------------------- REMOVE A PARTICLE FROM THE IS'th SECTION
                 occ(is) = occ(is) -1
!--------------------------SEARCH FOR PARENT ----------------
                 do j = 1,ntemplate(ith,nh-1)
                    match = .true.
                    do js = 1,nsections(ith)
                       if(occ(js) /= template(ith)%of_nh(nh-1)%list(j)%occ(js))then
                          match = .false.
                          exit
                       endif
                    enddo
                    if ( match ) then
                       template(ith)%of_nh(nh)%list(i)%parent=j
                       template(ith)%of_nh(nh)%list(i)%add = is
                       exit
                    endif
                 enddo  ! j
                 occ(is) = occ(is) +1
              endif
           endif
        enddo !is
     enddo !i
  enddo ! nh
  deallocate(occ)
  return
       
end subroutine template_heritage


!===========================================================================
!
!  subroutine write_templates
!
!  writes out templates -- purely for testing
!
subroutine write_templates(ith)

  use haiku_info
  implicit none
  integer :: ith
  integer :: nh
  integer :: i

  do nh = 1,maxNh(ith)
     print*,nh,' particles ',maxWh(ith,nh)
     do i = 1,ntemplate(ith,nh)
        write(6,101)template(ith)%of_nh(nh)%list(i)%w,       &
             template(ith)%of_nh(nh)%list(i)%parent,         &
             template(ith)%of_nh(nh)%list(i)%add,            &
             template(ith)%of_nh(nh)%list(i)%occ(:)
101     format(i3,'( ',2i4,' ) ',10i3)
     enddo ! i
  enddo  ! nh
  
end subroutine write_templates

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  subroutine get_template_ws(ith,nh,nw,wlist,whmin,whmax,wmap)
!
!  extracts information about list of Ws within the templates
!
!  INPUT:
!    ith = species
!    nh = # of particles
!
!  OUTPUT:
!    nw = # of values of w allowed for this nh
!    wlist(:) = list of allowed ws
!    whmin,whmax : = min, max values of wlist
!    wmap(w) = where in wlist can be found (may not be needed)
!
      subroutine get_template_ws(ith,nh,nw,wlist,whmin,whmax,wmap)

      implicit none
      integer ith
      integer nh
      integer nw
      integer, pointer :: wlist(:),wmap(:)
      integer whmin,whmax
      integer w

      integer i
      integer :: aerr

!-------------- INITIALIZE
      nw = 1
      whmin = template(ith)%of_nh(nh)%list(1)%w
      whmax = template(ith)%of_nh(nh)%list(1)%w
!--------------- COUNT UP DISTINCT VALUES OF W
      do i = 2,ntemplate(ith,nh)
        whmax = MAX(whmax,template(ith)%of_nh(nh)%list(i)%w)
        if(whmax > template(ith)%of_nh(nh)%list(i-1)%w)then
          nw = nw+1
        endif
      enddo  !i

      if(associated(wmap))nullify(wmap)
      if(associated(wlist))nullify(wlist)

      allocate(wlist(nw),wmap(whmin:whmax), stat=aerr)
      if(aerr /= 0) call memerror("get_template_ws 1")
      wmap(:) = -1
      nw = 1
      wlist(1) = whmin
      wmap(whmin) = 1
      w = whmin
      do i = 2,ntemplate(ith,nh)
        w = template(ith)%of_nh(nh)%list(i)%w
        if(w > template(ith)%of_nh(nh)%list(i-1)%w)then
          nw = nw+1
          wlist(nw) = w
          wmap(w) = nw
        endif
      enddo  !i
      return
      end subroutine get_template_ws
!=---------------------------------------------------------------------
!
!  added 7.7.2: fill out : jzsectplates, which contain information
!                          on min/max jz in a 'section' of s.p. states


subroutine master_of_jzsectplates
	use spsections
	implicit none
	integer :: ith
	integer :: aerr
	integer :: isection
	integer :: nmysps
	
!................ ALLOCATE JZARRAY.........

    nmysps = 0
    do ith = -2,2
		if(ith==0)cycle
		if(nsections(ith)<1)cycle

		
		do isection = 1,nsections(ith)
			nmysps = max(nmysps,section(ith)%list(isection)%nsps)
		end do ! isection
		
	end do ! ith	
	
	allocate(jzarray(nmysps),stat = aerr)
	if(aerr/=0)call memerror('master of jzsectplates 0 ')
	
	
	do ith = -2,2
		if(ith==0)cycle
		if(nsections(ith)<1)cycle
		
		allocate(jzsectplate(ith)%of_section(nsections(ith)),stat=aerr)
		if(aerr/=0)call memerror('master of jzsectplates 1 ')
		
		do isection = 1,nsections(ith)
			
		    if(section(ith)%list(isection)%nsps < 1)cycle
			nmysps = section(ith)%list(isection)%nsps
			allocate( jzsectplate(ith)%of_section(isection)%jzmax2_of_n(0:nmysps), & 
			jzsectplate(ith)%of_section(isection)%jzmin2_of_n(0:nmysps), stat = aerr )
			if(aerr/=0)then
				print*,ith,isection,nmysps,nsections(ith)

				call memerror('master of jzsectplates 2 ')
			end if
			call fill_jzsectplates(ith,isection)
	    end do    ! isection
		
	end do ! ith
	return
	
end subroutine master_of_jzsectplates
!=---------------------------------------------------------------------

subroutine fill_jzsectplates(ith,isection)
	use spsections
	use spstate
	
	implicit none
	integer,intent(in) :: ith
	integer,intent(in) :: isection
	
	integer :: nmysps,imysps,jmysps,k,m,mtmp
	integer :: inh
	
	
!........ FILL JZARRAY.....

     jzarray(:) = 0
	 nmysps = 	section(ith)%list(isection)%nsps
	 if(nmysps < 1)return
	 do imysps = section(ith)%list(isection)%start, section(ith)%list(isection)%fin
		 jzarray(imysps +1 -section(ith)%list(isection)%start )=hspsqn(ith,imysps)%m
	 end do
!........... NOW DO SIMPLE BUBBLE SORT...................................	 
     do imysps = 1,nmysps
		 k = imysps
		 m = jzarray(imysps)
		 do jmysps = imysps+1,nmysps
			 mtmp = jzarray(jmysps)
			 if(mtmp < m)then
				 k = jmysps
				 m = mtmp
			 end if
		 end do ! jmysps
		 if(k /= imysps)then   ! swap
			 mtmp = jzarray(imysps)
			 jzarray(imysps)=jzarray(k)
			 jzarray(k)=mtmp
		 end if
		 
	 end do ! imysps
!	 print*,jzarray(1:nmysps)

!............ LOOP OVER # OF ALLOWED OCCUPATIONS AND FIND MIN,MAX JZ.....	
    jzsectplate(ith)%of_section(isection)%jzmax2_of_n(0)=0
    jzsectplate(ith)%of_section(isection)%jzmin2_of_n(0)=0
	
	do inh = 1,nmysps
		m = 0
		mtmp = 0
		do imysps= 1,inh
			m = m+jzarray(imysps)
			mtmp = mtmp + jzarray(nmysps+1-imysps)
		end do
	    jzsectplate(ith)%of_section(isection)%jzmin2_of_n(inh)=m
	    jzsectplate(ith)%of_section(isection)%jzmax2_of_n(inh)=mtmp
	end do	
	
	return
	
end subroutine fill_jzsectplates

!=---------------------------------------------------------------------
!
! using the jzsectplates, compute the max/min Jz values for templates
!
! also finds list of W values for templates and finds min/max for them
!
subroutine assign_jz_to_templates
	use haiku_info
	use spsections
	implicit none
	integer :: ith,nh,itemplate
	integer :: mmax,mmin,iocc,isection
	integer :: wlow,whi,whowmany,ww
	integer :: aerr
	
	
	do ith = -2,2
		if(ith==0)cycle
!		template(ith)%of_nh(0)%list(1)%jzmin = 0
!		template(ith)%of_nh(0)%list(1)%jzmax = 0
		template(ith)%of_nh(0)%wmin = 0
		template(ith)%of_nh(0)%wmax = 0
		
		allocate( template(ith)%of_nh(0)%jzmin_of_w(0:0),stat=aerr)
		if(aerr/=0)call memerror('assign jz to templates 1 ')
		
		allocate( template(ith)%of_nh(0)%jzmax_of_w(0:0),stat=aerr)
		if(aerr/=0)call memerror('assign jz to templates 2 ')
		template(ith)%of_nh(0)%jzmin_of_w(0) = 0
		template(ith)%of_nh(0)%jzmax_of_w(0) = 0

		do nh = 1,maxNh(ith)
			if(ntemplate(ith,nh)==0)cycle
			
			wlow = 10000
			whi  = -1
			do itemplate = 1,ntemplate(ith,nh)
				mmin = 0
				mmax = 0
				wlow = min(wlow, template(ith)%of_nh(nh)%list(itemplate)%W)
				whi =  max(whi, template(ith)%of_nh(nh)%list(itemplate)%W)

				do isection = 1,nsections(ith)
					iocc = template(ith)%of_nh(nh)%list(itemplate)%occ(isection)
				    mmin = mmin + jzsectplate(ith)%of_section(isection)%jzmin2_of_n(iocc)
				    mmax = mmax + jzsectplate(ith)%of_section(isection)%jzmax2_of_n(iocc)
				
				end do
				template(ith)%of_nh(nh)%list(itemplate)%jzmin = mmin
				template(ith)%of_nh(nh)%list(itemplate)%jzmax = mmax
					
			end do ! itemplate
			template(ith)%of_nh(nh)%wmin = wlow
			template(ith)%of_nh(nh)%wmax = whi		
			allocate( template(ith)%of_nh(nh)%jzmin_of_w(wlow:whi),stat=aerr)
			if(aerr/=0)call memerror('assign jz to templates 3 ')			
		
			allocate( template(ith)%of_nh(nh)%jzmax_of_w(wlow:whi),stat=aerr)
			if(aerr/=0)call memerror('assign jz to templates 4 ')	
			template(ith)%of_nh(nh)%jzmin_of_w(:) =  1000
			template(ith)%of_nh(nh)%jzmax_of_w(:) =  -1000
			
			do itemplate = 1,ntemplate(ith,nh)
				ww = template(ith)%of_nh(nh)%list(itemplate)%w
				mmin = template(ith)%of_nh(nh)%list(itemplate)%jzmin
				template(ith)%of_nh(nh)%jzmin_of_w(ww) = min(mmin,template(ith)%of_nh(nh)%jzmin_of_w(ww) )
				mmax = template(ith)%of_nh(nh)%list(itemplate)%jzmax
				template(ith)%of_nh(nh)%jzmax_of_w(ww) = max(mmax,template(ith)%of_nh(nh)%jzmax_of_w(ww) )
			end do ! itemplate			
			
!			do ww = template(ith)%of_nh(nh)%wmin,template(ith)%of_nh(nh)%wmax
!				print*,ith,' template ',nh,ww,template(ith)%of_nh(nh)%jzmin_of_w(ww),template(ith)%of_nh(nh)%jzmax_of_w(ww)
!			end do
						
		end do ! nh
		
	end do ! ith
	return
	
end subroutine assign_jz_to_templates
!=---------------------------------------------------------------------
!
!  last step: match up templates to find, for each ith, nh, and wh, min and max Jz
!
!  first, for the highest nh for each ith, for each wh find min and max Jz
!

subroutine reconcile_templates_for_Jz
	use system_parameters
	use W_info
	use haiku_info
	
	implicit none
	
	integer :: itx,ity
	integer :: nhx,nhxc,nhy,nhyc
	integer :: whx, whxc, why, whyc
	integer :: minjzh,maxjzh,tmpjzh
	
!	do itx = -2,2
!		if(itx == 0)cycle
!		do nhx = minNh(itx),maxNh(itx)
!			do whx = template(itx)%of_nh(nhx)%wmin, template(itx)%of_nh(nhx)%wmax
!		print*,' limits ',itx,nhx,whx,':',template(itx)%of_nh(nhx)%jzmin_of_w(whx),template(itx)%of_nh(nhx)%jzmax_of_w(whx)
!			end do
!		end do		
!	end do
	
	do itx = -2,2
		if(itx == 0)cycle
		
		ity = 3-abs(itx)
		do nhx = maxNh(itx),minNh(itx),-1
		nhxc = np(abs(itx))-nhx
		do whx = template(itx)%of_nh(nhx)%wmax, template(itx)%of_nh(nhx)%wmin,-1

	            minjzh = 10000
				maxjzh = -10000			
! LOOP OVER ALL OTHER POSSIBILITIES.....................................
 
                do whxc = template(-itx)%of_nh(nhxc)%wmin, template(-itx)%of_nh(nhxc)%wmax
					if(whx + whxc > maxW(abs(itx)) ) cycle
					
					do nhy = minNh(ity),maxNh(ity)
						nhyc = np(ity)-nhy
						do why = template(ity)%of_nh(nhy)%wmin, template(ity)%of_nh(nhy)%wmax
							if(why + whx +whxc > maxWtot)cycle
							
							do whyc = template(-ity)%of_nh(nhyc)%wmin, template(-ity)%of_nh(nhyc)%wmax
								if(why + whyc > maxW(ity))cycle
								if(why + whyc + whx +whxc > maxWtot)cycle
								
								tmpjzh = template(-itx)%of_nh(nhxc)%jzmax_of_w(whxc)
							
							    tmpjzh = tmpjzh + template(ity)%of_nh(nhy)%jzmax_of_w(why)
							    tmpjzh = tmpjzh + template(-ity)%of_nh(nhyc)%jzmax_of_w(whyc)
								
							!	print*,tmpjzh
								maxjzh = max(maxjzh,tmpjzh)
								
								tmpjzh = template(-itx)%of_nh(nhxc)%jzmin_of_w(whxc)
							
							    tmpjzh = tmpjzh + template(ity)%of_nh(nhy)%jzmin_of_w(why)
							    tmpjzh = tmpjzh + template(-ity)%of_nh(nhyc)%jzmin_of_w(whyc)
							!	print*,tmpjzh
								minjzh = min(minjzh,tmpjzh)

								
								
							end do  ! whyc
						end do ! why
					end do ! nhy
					
				end do !whxc
				
				if(itx < 0)then
					if(maxjzh==-10000)cycle
!					print*,' limits ',itx,nhx,whx,':',template(itx)%of_nh(nhx)%jzmin_of_w(whx),template(itx)%of_nh(nhx)%jzmax_of_w(whx),jz-maxjzh
					template(itx)%of_nh(nhx)%jzmin_of_w(whx)= max(jz-maxjzh,template(itx)%of_nh(nhx)%jzmin_of_w(whx))
					
				else
					if(minjzh==10000)cycle
					
!					print*,' limits ',itx,nhx,whx,':',template(itx)%of_nh(nhx)%jzmin_of_w(whx),template(itx)%of_nh(nhx)%jzmax_of_w(whx),jz-minjzh
					template(itx)%of_nh(nhx)%jzmax_of_w(whx)=min(jz-minjzh,template(itx)%of_nh(nhx)%jzmax_of_w(whx))
					
				end if 
			end do  ! whx
		end do ! nhx		
	end do ! itx

	
	return
	
	
end subroutine reconcile_templates_for_Jz


!======================================================================	  
  end module templates

