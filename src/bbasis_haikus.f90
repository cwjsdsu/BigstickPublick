!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK configuration interaction shell-model code
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  BBASIS_HAIKUS.f
!
!  third set of routines to generate the basis:
!  HAIKUS -- sorting 
!
! SUBROUTINES:
! sort_haikus -- sorts haikus by parity, Jz
! test_haiku_qns
!	calls: convert_haiku, jz_of_occ, par_of_occ
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
! NOTE: it is not clear that sorting is really needed, as I later 
! reorganize the haikus into blocks. This only slows down the code
! for very large cases, so it is acceptable for now
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!---------------------------------------------------------------------------
!  list of haikus
!---------------------------------------------------------------------------
module haikus
  implicit none

  type haikulist
     integer          :: nhaiku
     integer, pointer :: par(:),jz(:)
     integer, pointer :: hsd(:,:)     
  end type haikulist

  type haikuWlist
     integer          :: nw        ! # of ws allowed
     integer, pointer :: wlist(:)  ! list of allowed ws
     integer, pointer :: wmap(:)   ! map of w-value onto wlist
     type (haikulist), pointer :: of_w(:)
  end type haikuWlist

  type haikuNlist
     type ( haikuWlist ),pointer :: of_nh(:)
  end type haikuNlist

  type ( haikuNlist   ) :: haiku(-2:2)

contains
!=====
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK configuration-interaction shell-model code
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  BBASISLIB2.f
!
!  first set of routines to generate the basis:
!  HAIKUS
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  NOTES:  The goal of these routines is to generate the "haikus."
!  Haikus are half-Slater-determinants: left-haikus are constructed from 
!  s.p. states with jz < 0, and right-haikus with jz >= 0. Combining 
!  them, one constructs the proton and neutron Slater determinants. 
!  There are usually far fewer haikus than proton and neutron SDs.
!
!  Haikus are constructed recursively using the "parent" information for
!  templates, which are already organized by nh (particle number) and W. 
!  (After this point the templates can be discarded.) 
!  Once constructed, they are sorted by parity and jz and then organized 
!  into "blocks" which are labeled by the quantum numbers:  nh, wh, parh, and jzh.
!C  The location of an individual haiku within a block is its "address."
!C  Later, "hops" will encode the application of creation/destruction operators
!C  on haikus, listed by initial and final block/address. After this point, the 
!C  representation of haikus as binary words can be discarded.  Because we 
!C  only use the haikus briefly, we don't have to worry about fetching them from
!C  from memory efficiently and can point to a list out-of-order.
!C
!C  Step 1: For each nh, generate all haikus for each allowed W, using templates
!C  Step 2: For each nh, wh, sort by parity and then by jz
!C  Step 3: Compute haiku blocks (lists of nh-wh-parh-jzh and # of haikus)
!C  Step 4: create block_lists 
!C
!C  In order to know how to effectively store haikus, one must think ahead to how 
!c  the proton / neutron SDs will be constructed from haikus. 
!C
!C  The proton / neutron SDs are organized by sectors, labeled by Jz, parity, W.
!C  Each sector is defined by a list of paired haiku blocks. 
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  SUBROUTINES :
!C   ** haiku master **:  master routine to call others
!C     calls subroutines get_templates_ws
!C                    :  count_create_haikus:
!C        		haiku_qns : 
!C  (** test_haikus **   for testing only)
!   ** get_templates_ws **  extracts info on Ws from templates
!   ** count_create_haikus ** recursively creates haikus
!C     calls subroutines  convert_haiku  
!                        add_haiku_bit  : 
!     calls functions    bit_occ
!   **  add_haiku_bit **  adds a bit to a haiku
!     call subroutine  convert_loc_word
!   **  convert_loc_word ** converts the loc of a state to a word and specific bit
!   **  haiku_qns   **   computes q#s (jz, parity) of a set of haikus
!     calls subroutine convert_haiku
!     call functions  jzh_of_occ
!                     par_of_occ
!   ** convert_haiku **  converts a haiku in binary form to an occupation array
!     call function  bit_occ
!
!  (** convert_haiku_array **  converts a haiku in binary form to an occupation array
!                        does not assume a fix # of particles; seldom used)

!  FUNCTIONS:
!   ** bit_occ **        sees if a bit (state) in a haiku is occupied
!   ** jzh_of_occ **   computes jzh of an occ
!   ** par_of_occ0 **   computes parity of an occ
!   ** parmult  ** multiples two parities of the form 1,2
!
!  (OBSOLETE OR SELDOM USED FUNCTIONS )
!   ** nh_of_occ **     computes nh of an occ
!   ** jzh_of_occ0 **   computes jzh of an occ, checking ALL bits
!   ** wh_of_occ0 **   computes jzh of an occ, checking ALL bits
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!  subroutine haiku_master
!
!  coordinates creation of haikus
!
!  INPUT: ith = species/handness
!
!  SUBROUTINES CALLED:
!    get_templates_ws : extracs information on Ws allowed for fixed nh
!    count_create_haikus: creates the haikus recursively
!    haiku_qns : computes haiku quantum #s
!
!
subroutine haiku_master(ith)

use bitstuff
use haiku_info
use templates
use verbosity
implicit none

integer ith
integer it
integer nh
integer nw,whmin,whmax
integer, pointer :: wlist(:),wmap(:)
integer iw,wh
integer nhaiku,ntemp
integer, allocatable :: occ(:)
integer :: aerr

it = abs(ith)

allocate(occ(nhsps(abs(ith))), stat=aerr)
if(aerr /= 0) call memerror("haiku_master 1")
allocate(haiku(ith)%of_nh(0:maxNh(ith)), stat=aerr)
if(aerr /= 0) call memerror("haiku_master 2")

!------------- DERIVE THE HAIKU for zero particles

haiku(ith)%of_nh(0)%nw = 1
allocate(haiku(ith)%of_nh(0)%wlist(1), stat=aerr)
if(aerr /= 0) call memerror("haiku_master 3")
allocate(haiku(ith)%of_nh(0)%wmap(0:0), stat=aerr)
if(aerr /= 0) call memerror("haiku_master 4")
allocate(haiku(ith)%of_nh(0)%of_w(1), stat=aerr)
if(aerr /= 0) call memerror("haiku_master 5")
haiku(ith)%of_nh(0)%wlist(1) = 0
haiku(ith)%of_nh(0)%wmap(0) = 1
haiku(ith)%of_nh(0)%of_w(1)%nhaiku = 1

allocate(haiku(ith)%of_nh(0)%of_w(1)%hsd(nword(it),1 ), stat=aerr)
if(aerr /= 0) call memerror("haiku_master 6")
haiku(ith)%of_nh(0)%of_w(1)%hsd = 0
call haiku_qns(ith,0,1,occ)
!------------- loop over nh
do nh = 1,maxNh(ith)
!------------- FIGURE OUT W INFORMATION --------------

   call get_template_ws(ith,nh,nw,wlist,whmin,whmax,wmap)

!------------- ALLOCATE ------------------------
   haiku(ith)%of_nh(nh)%nw = nw
   allocate(haiku(ith)%of_nh(nh)%wlist(nw), stat=aerr)
   if(aerr /= 0) call memerror("haiku_master 10")
   allocate(haiku(ith)%of_nh(nh)%wmap(whmin:whmax), stat=aerr)
   if(aerr /= 0) call memerror("haiku_master 11")
   do iw = 1,nw
      haiku(ith)%of_nh(nh)%wlist(iw) = wlist(iw)
   enddo !iw
   do iw = whmin,whmax
      haiku(ith)%of_nh(nh)%wmap(iw) = wmap(iw)
   enddo  !iw

   allocate(haiku(ith)%of_nh(nh)%of_w(nw), stat=aerr)
   if(aerr /= 0) call memerror("haiku_master 12")

   do iw = 1,nw
      wh = wlist(iw)
!-------------- COUNT UP # OF HAIKUS
      call count_create_haikus(.false.,ith,nh,wh,iw,nhaiku)
!------------- ALLOCATE --------------------
      haiku(ith)%of_nh(nh)%of_w(iw)%nhaiku = nhaiku

      allocate(haiku(ith)%of_nh(nh)%of_w(iw)%hsd(nword(it),nhaiku), stat=aerr)
      if(aerr /= 0) call memerror("haiku_master 13")

!--------------- STORE THE HAIKUS ----------------------
      call count_create_haikus(.true.,ith,nh,wh,iw,nhaiku)
!--------------- COMPUTE QUANTUM NUMBERS --------------------
      call haiku_qns(ith,nh,iw,occ)
!----------------- SORT -------------------------------
   enddo  !iw
enddo ! nh
do nh = 1,maxNh(ith)
    nw = haiku(ith)%of_nh(nh)%nw
enddo  ! nh
deallocate(occ)
!---------------------- WRITE OUT IF VERBOSE --------------

if(verbose_haikus)then
   nhaiku = 0
   print*,' '
   print*,' species ',ith
   do nh = 0,maxNh(ith)
          ntemp = 0
         nw = haiku(ith)%of_nh(nh)%nw 
         do iw = 1,nw
           ntemp = ntemp + haiku(ith)%of_nh(nh)%of_w(iw)%nhaiku
         enddo  !iw
         print*,' For nh ',nh,' there are ',ntemp,' haikus '
         nhaiku = nhaiku + ntemp
   enddo
   print*,' There are a total of ',nhaiku,' haikus '
   print*,' '
endif
return
end subroutine haiku_master
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  for testing purposes only
!
!  CALLS: 
!   convert_haiku_array
!   nh_of_occ
!   wh_of_occ
!
      subroutine test_haikus(ith)

      use haiku_info
      implicit none
!      include 'binterfaces.inc'

      integer ith
      integer nh
      integer wh,iwh,nw
      integer nhaiku,ih
      
      integer,allocatable :: occ(:)
      integer,pointer  :: hsd(:)
      integer :: aerr

!      integer :: nh_of_occ,wh_of_occ

      allocate(occ(nhsps(abs(ith))), stat=aerr)
      if(aerr /= 0) call memerror("test_haikus 1")

      do nh = 0,maxNh(ith)
        nw = haiku(ith)%of_nh(nh)%nw
        do iwh = 1,nw
          wh = haiku(ith)%of_nh(nh)%wlist(iwh)
          nhaiku = haiku(ith)%of_nh(nh)%of_w(iwh)%nhaiku
          do ih = 1,nhaiku
           hsd => haiku(ith)%of_nh(nh)%of_w(iwh)%hsd(:,ih)
           call convert_haiku_array(ith,hsd,occ)
103      format(60i1)
           if(nh /= nh_of_occ(ith,occ))then
              print*,' nh problem ',nh,nh_of_occ(ith,occ),iwh,ih

           endif
           if(wh /= wh_of_occ(ith,nh,occ))then
              print*,' wh problem ',wh,ih
           endif

          enddo

        enddo ! iwh
      enddo ! nh

      deallocate(occ)
      return
      end subroutine test_haikus

! subroutine get_template_ws MOVED in 7.6.9 to bbasis_template.f90 containing the module templates
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  subroutine count_create_haikus
!
!  counts up and creates haikus recursively, using "parents" in templates
!
!  INPUT:
!   create : logical flag to create (.true.) or just count up (.false.)
!   ith = species/handness
!   nh = # of particles
!   wh = w-value of haikus
!   iwh = label of which wh for haikus
!
! OUTPUT (if create = .false.)
!  nhaiku = # of haikus for this nh,wh
!
! SUBROUTINES CALLED:
!  convert_haiku  : converts a haiku (binary word) to occupation array (occ)
!  add_haiku_bit  : directly adds a bit to a haiku
!
! FUNCTIONS CALLED:
!  bit_occ = .true. if a single particle state is occupied
!

      subroutine count_create_haikus(create,ith,nh,wh,iwh,nhaiku)

      use spsections
      use templates
      use haiku_info
      implicit none

!.....................................................
      logical create
      integer ith
      integer nh
      integer wh,iwh
      integer nhaiku
!......................
      integer it
      integer i
      integer parent,wparent,nparent,iwparent
      integer  ih
      integer padd
      integer pstart,pfin
      integer p
      integer iword, ibit
      integer hstart,hfin
      integer,pointer :: hsd(:),hsd2(:)
!      integer addbit
!      integer nh_of_occ
      integer occ(nhsps(abs(ith)))
      it = abs(ith)
      nhaiku = 0

!      if(create)print*,' creating ',ith,nh,wh,iwh,ntemplate(ith,nh)
      do i = 1,ntemplate(ith,nh)
        if(template(ith)%of_nh(nh)%list(i)%w < wh)cycle
        if(template(ith)%of_nh(nh)%list(i)%w > wh)exit

       if(.not.create)template(ith)%of_nh(nh)%list(i)%hstart = nhaiku+1

        padd    = template(ith)%of_nh(nh)%list(i)%add
        parent = template(ith)%of_nh(nh)%list(i)%parent

        if(nh == 1)then
           wparent = 0
           iwparent = 1
           nparent = 1
           hstart = 1
           hfin = 1
        else
           wparent = template(ith)%of_nh(nh-1)%list(parent)%w
           iwparent = haiku(ith)%of_nh(nh-1)%wmap(wparent)
           hfin = template(ith)%of_nh(nh-1)%list(parent)%hfin
           hstart = template(ith)%of_nh(nh-1)%list(parent)%hstart
        endif

        pstart = section(ith)%list(padd)%start
        pfin = section(ith)%list(padd)%fin
!		if(create)print*,nh,' p fin start ',ith,padd,pstart,pfin
        do ih = hstart,hfin

          hsd => haiku(ith)%of_nh(nh-1)%of_w(iwparent)%hsd(:,ih)

          do p = pfin,pstart,-1  ! loop over and find open space
             if(.not.bit_occ(hsd,p))then
                nhaiku = nhaiku+1
                if(create)then
                   haiku(ith)%of_nh(nh)%of_w(iwh)%hsd(:,nhaiku)=      haiku(ith)%of_nh(nh-1)%of_w(iwparent)%hsd(:,ih)
       call add_haiku_bit(haiku(ith)%of_nh(nh)%of_w(iwh)%hsd(:,nhaiku),nword(it),p)

!------------------------ TESTING ------------------------
          hsd2 => haiku(ith)%of_nh(nh)%of_w(iwh)%hsd(:,nhaiku)
           call convert_haiku(ith,nh,hsd2,occ)

                endif
              else
                exit  
             endif
          enddo  ! p
        enddo  ! ih
        if(.not.create)then  ! fill in # of haikus per template
          template(ith)%of_nh(nh)%list(i)%hfin = nhaiku
        endif
101     format('so far ',6i3)
      enddo ! i

      return
      end subroutine count_create_haikus

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  function bit_occ
!
!  determines if the isps state is occupied
!  NB: this function is called a LOT; some tweaking to optimize has been done
!  but possibly could be improved (however not currently a significant bottleneck)
!
!  INPUT:
!   hsd(:) = pointer of binary words that make up the haiku
!   isps = which s.p. state
!
!  OUTPUT:
!   bit_occ = .true. if there is a bit at location isps
!
      function bit_occ(hsd,isps)
      use bitstuff
      implicit none
      logical bit_occ
      integer,pointer :: hsd(:)
      integer isps

      integer iword,i
      integer imask
      
      iword = (isps-1)/(max_bit_word) + 1
      i = isps - (iword -1)*(max_bit_word)
      imask = ishft(1,i-1)
	  bit_occ = .not.(iand(imask,hsd(iword)) == 0)   ! should be faster than an if statement
!      if(iand(imask,hsd(iword)) == 0)then
!        bit_occ = .false.
!      else
!        bit_occ = .true.
!      endif
      return
      end function bit_occ

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  subroutine add_haiku_bit
!    adds a bit to a haiku to represent an occupied state
!
! INPUT:  
!   hsd(:) = array of binary words = haiku
!   n = nwords in haiku
!   isps = location of state
!
!  subroutines called: convert_loc_word
!
      subroutine add_haiku_bit(hsd,n,isps)
      implicit none
      integer :: hsd(n)
      integer isps
      integer n
      integer i
      integer iword
      integer imask
      call convert_loc_word(isps,iword,i)
      imask = ishft(1,i-1)
      hsd(iword) = hsd(iword) + imask
      return
      end subroutine add_haiku_bit

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  function addbit
!
!  adds a bit in the isps'th location to a binary word (haiku)
!
!  INPUT: isps  = which s.p. state
!
!  OUTPUT: addbit = binary word to be added directly
!
!  SUBROUTINES CALLED:
!   convert_loc_word: computes where the word, bit of the state is located
!
      integer function addbit(isps)
      implicit none
      integer isps

      integer iword
      integer i
      call convert_loc_word(isps,iword,i)
      addbit = ishft(1,i-1)
      return
      end function addbit

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  subroutine convert_loc_word
!
!  converts a location in abstract single-particle space to 
!  an absolute location in terms of words and bits
!
!  each haiku is stored in binary over several words
!
!  the location of a single-particle state is which bit on which word
!  We have to map an abstract index for the state (e.g., 29, 44, 83)
!
!  for example, suppose we have 31 bits per "word"
!  If isps = 29, then iword = 1,i = 29
!  If isps = 44, then iword = 2,i = 13
!  If sps = 83,  then iword = 3,i = 11
!
!  INPUT: isps  = abstract location
!
!  OUTPUT: 
!         iword = which word
!         i     = which bit of iword
!

      subroutine convert_loc_word(isps,iword,i)
      use bitstuff
      implicit none
      integer isps,iword,i

      iword = (isps-1)/(max_bit_word) + 1
      i = isps - (iword -1)*(max_bit_word)
      return

      end subroutine convert_loc_word
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine haiku_qns
!
!   computes the quantum numbers jz, wh, parh for a list of haikus
!   assumes the haikus have already been segregated by nh, w
!
!  INPUT:
!    ith  = species/handness
!    nh = # of particles, iwh = which W value
!    occ(:)  an array used repeatedly 
!
!  SUBROUTINES CALLED:
!    convert_haiku : converts a haiku (hsd) to occupation array (occ)
!
!  FUNCTIONS CALLED:
!    jzh_of_occ, par_of_occ
!
      subroutine haiku_qns(ith,nh,iwh,occ)

      use haiku_info
      implicit none

      integer ith
      integer nh
      integer iwh
      integer :: occ(:)
!............................................
      integer nhaiku
      integer ih
      integer,pointer :: hsd(:)
      integer :: aerr
!      integer jzh_of_occ,par_of_occ

      nhaiku = haiku(ith)%of_nh(nh)%of_w(iwh)%nhaiku
      allocate(haiku(ith)%of_nh(nh)%of_w(iwh)%par(nhaiku), stat=aerr)
      if(aerr /= 0) call memerror("haiku_qns 1")
      allocate(haiku(ith)%of_nh(nh)%of_w(iwh)%jz(nhaiku), stat=aerr)
      if(aerr /= 0) call memerror("haiku_qns 2")
!$omp parallel do private(ih,hsd,occ)
      do ih = 1,nhaiku
         hsd => haiku(ith)%of_nh(nh)%of_w(iwh)%hsd(:,ih)
         call convert_haiku(ith,nh,hsd,occ)
         haiku(ith)%of_nh(nh)%of_w(iwh)%jz(ih) = jzh_of_occ(ith,nh,occ)
         haiku(ith)%of_nh(nh)%of_w(iwh)%par(ih) =  par_of_occ(ith,nh,occ)
      enddo
!$omp end parallel do
      return
      end subroutine haiku_qns

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  subroutine convert_haiku_array
!
!  converts a haiku  (binary word) to an array 
!
!  INPUT:
!   ith = species/handness
!   hsd(:) = pointer to binary words that make up the haiku
!
!  OUTPUT:
!    occ(:)  array = 0. or 1 designating occupied s.p. states
!
!  FUNCTIONS CALLED:
!  bit_occ: computes whether or not a specific bit is occupied
!
      subroutine convert_haiku_array(ith,hsd,occ)

      use haiku_info
      implicit none

      integer :: ith
      integer,pointer :: hsd(:)
      integer :: occ(nhsps(abs(ith)))
      integer isps

      occ = 0
      do isps = 1,nhsps(abs(ith))
        if(bit_occ(hsd,isps))occ(isps) =1
      enddo
      return
      end subroutine convert_haiku_array
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  subroutine convert_haiku
!  converts a haiku  (binary word) to an array ASSUMING fixed nh
!  faster because it exits 
!  note this has a different definition of occ from convert_haiku_array
!
!  INPUT:
!   ith = species/handness
!   nh = # of paricles
!   hsd(:) = pointer to binary words that make up the haiku
!
!  OUTPUT:
!    occ(:)  array = 0. or 1 designating occupied s.p. states
!
!
      subroutine convert_haiku_OLD(ith,nh,hsd,occ)
      use bitstuff
      use haiku_info
      implicit none


      integer :: ith
      integer :: nh
      integer,pointer :: hsd(:)
      integer :: occ(:)

      integer isps
      integer n
      integer i, iword, imask
      occ = 0
      n = 0
      do isps = 1,nhsps(abs(ith))
         iword = (isps-1)/(max_bit_word) + 1
         i = isps - (iword -1)*(max_bit_word)
         imask = ishft(1,i-1)
         if(iand(imask,hsd(iword)) /= 0)then
            n = n + 1
            occ(n) = isps
            if(n==nh)return
         endif
      enddo
      if(n < nh)then
        print*,' problem converting ',n,nh
        stop
      endif
      return
      end subroutine convert_haiku_OLD

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
!  subroutine convert_haiku NEW
!  converts a haiku  (binary word) to an array ASSUMING fixed nh
!  NEW VERSION March 2013 to speed up usage for certain cases
!  reduces division which is slow
!  note this has a different definition of occ from convert_haiku_array
!
!  INPUT:
!   ith = species/handness
!   nh = # of paricles
!   hsd(:) = pointer to binary words that make up the haiku
!
!  OUTPUT:
!    occ(:)  array = 0. or 1 designating occupied s.p. states
!
!
      subroutine convert_haiku(ith,nh,hsd,occ)
      use bitstuff
      use haiku_info
      implicit none


      integer :: ith
      integer :: nh
      integer,pointer :: hsd(:)
      integer :: occ(:)

      integer isps
      integer n
      integer i, iword, imask,ilast,jhsd
      occ = 0
      n = 0

      isps = 0
      do iword = 1,Nword( abs(ith) )
         jhsd=hsd(iword)
!............. SET UP SIZE OF INNER DO LOOP .....
         if(iword == Nword( abs(ith) ))then
            ilast = nhsps(abs(ith))-(Nword(abs(ith)) -1)*(max_bit_word)
         else
            ilast = max_bit_word
         end if
         imask = 1
         do i = 1,ilast
            isps = isps + 1
            if(iand(imask,jhsd) /= 0)then
               n = n + 1
               occ(n) = isps
            endif
            imask = ishft(imask,1)  ! shift by 1
         end do
         if(n==nh)return

      end do ! iword

      if(n < nh)then
        print*,' problem converting ',n,nh
        stop
      endif
      return
      end subroutine convert_haiku

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  function nh_of_occ0  (mostly obsolete)
!
!  computes the nh of a haiku, using an occupation array
!  only used for testing, to verify that an occupation array has nh particles
!
!  input: ith = species/handness
!      occ(:) = occupation of s.p. states
!
!
      function nh_of_occ(ith,occ)
      use spsections
      use haiku_info

      implicit none
      integer :: ith
      integer :: occ(nhsps(abs(ith)))
      integer nh_of_occ
      integer isps

      nh_of_occ = 0
      do isps = 1,nhsps(abs(ith))
         nh_of_occ = nh_of_occ + occ(isps)
      enddo
      return
      end function nh_of_occ
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  function wh_of_occ0  (mostly obsolete)
!
!  computes the w of a haiku, using an occupation array
!  the difference from wh_of_occ is does not use nh to stop counting
!
!  input: ith = species/handness
!      occ(:) = occupation of s.p. states
!
!
      function wh_of_occ0(ith,occ)
      use spstate
      use haiku_info

      implicit none
      integer :: ith
      integer :: occ(nhsps(abs(ith)))
      integer wh_of_occ0
      integer isps

      wh_of_occ0 = 0
      do isps = 1,nhsps(abs(ith))
         wh_of_occ0 = wh_of_occ0 + occ(isps)*hspsqn(ith,isps)%w
      enddo
      return
      end function wh_of_occ0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  function jzh_of_occ0  (mostly obsolete)
!
!  computes the jz of a haiku, using an occupation array
!  the difference from jzh_of_occ is does not use nh to stop counting
!
!  input: ith = species/handness
!      occ(:) = occupation of s.p. states
!
!
      function jzh_of_occ0(ith,occ)
      use spstate
      use haiku_info

      implicit none
      integer :: ith
      integer :: occ(nhsps(abs(ith)))
      integer jzh_of_occ0
      integer isps

      jzh_of_occ0 = 0
      do isps = 1,nhsps(abs(ith))
         jzh_of_occ0 = jzh_of_occ0 + occ(isps)*hspsqn(ith,isps)%m
      enddo
      return
      end function jzh_of_occ0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  function wh_of_occ
!
!  computes the wh of a haiku, using an occupation array
!
!  input: ith = species/handness
!         nh = # of particles
!      occ(:) = occupation of s.p. states
!
!
      function wh_of_occ(ith,nh,occ)
      use spstate
      use haiku_info

      implicit none
      integer :: ith
      integer :: nh
      integer :: occ(nhsps(abs(ith)))
      integer wh_of_occ
      integer isps

      wh_of_occ = 0
      do isps = 1,nh
         wh_of_occ = wh_of_occ + hspsqn(ith,occ(isps))%w
      enddo
      return
      end function wh_of_occ

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  function jzh_of_occ
!
!  computes the jz of a haiku, using an occupation array
!
!  input: ith = species/handness
!         nh = # of particles
!      occ(:) = occupation of s.p. states
!
!
      function jzh_of_occ(ith,nh,occ)
      use spstate
      use haiku_info

      implicit none
      integer :: ith
      integer :: nh
      integer :: occ(nhsps(abs(ith)))
      integer jzh_of_occ
      integer isps

      jzh_of_occ = 0
      do isps = 1,nh
         jzh_of_occ = jzh_of_occ + hspsqn(ith,occ(isps))%m
      enddo
      return
      end function jzh_of_occ

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  function par_of_occ
!
!  computes the parity (1 + +, 2 = -) of a haiku, using an occupation array
!
!  input: ith = species/handness
!         nh = # of particles
!      occ(:) = occupation of s.p. states
!
!
      function par_of_occ(ith,nh,occ)

      
      use sporbit
      use spstate
      use haiku_info

      implicit none
      integer :: ith
      integer :: nh
      integer :: occ(nhsps(abs(ith)))
      integer par_of_occ
      integer isps

      par_of_occ = 1
      if(allsameparity)return
      if(nh == 0)return

      do isps = 1,nh
         par_of_occ = par_of_occ*hspsqn(ith,occ(isps))%par
      enddo
      par_of_occ = (3-par_of_occ)/2
      if(par_of_occ < 1 .or. par_of_occ > 2)then
        print*,' problem with parity ',par_of_occ,nh
        stop
      endif
      return
  end       function par_of_occ

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      subroutine test_haiku_qns
!
!   test the quantum numbers jz, wh, parh for a list of haikus
!   assumes the haikus have already been segregated by nh, w
!
!  INPUT:
!    ith  = species/handness
!    nh = # of particles, iwh = which W value
!    occ(:)  an array used repeatedly 
!
!  SUBROUTINES CALLED:
!    convert_haiku : converts a haiku (hsd) to occupation array (occ)
!
!  FUNCTIONS CALLED:
!    jzh_of_occ, par_of_occ
!
      subroutine test_haiku_qns(ith,nh,iwh,occ)

      use haiku_info
      implicit none

      integer ith
      integer nh
      integer iwh
      integer :: occ(:)
!............................................
      integer nhaiku
      integer ih
      integer,pointer :: hsd(:)
!      integer jzh_of_occ,par_of_occ

      nhaiku = haiku(ith)%of_nh(nh)%of_w(iwh)%nhaiku

      do ih = 1,nhaiku
         hsd => haiku(ith)%of_nh(nh)%of_w(iwh)%hsd(:,ih)
         call convert_haiku(ith,nh,hsd,occ)
         if(haiku(ith)%of_nh(nh)%of_w(iwh)%jz(ih) /=  jzh_of_occ(ith,nh,occ))then

             print*,' problem with sorted jz '
             print*,ith,nh,iwh
             print*,ih
             print*,haiku(ith)%of_nh(nh)%of_w(iwh)%jz(ih)
             print*,jzh_of_occ(ith,nh,occ)
             print*,' '
             print*,nhaiku
             print*,haiku(ith)%of_nh(nh)%of_w(iwh)%jz(:)

            stop
         endif
         if(haiku(ith)%of_nh(nh)%of_w(iwh)%par(ih) /= par_of_occ(ith,nh,occ))then

            print*,' problem with sorted par '
             print*,ith,nh,iwh
             print*,ih
            stop
         endif
      enddo
      return
      end subroutine test_haiku_qns

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

end module haikus