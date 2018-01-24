!========================================================================
!
! BJUMPLIB_NTUPLES.f90 
!
!
! initiated Sept 2009 by CWJ @ SDSU
!========================================================================

!
!  information on Q#s for n-tuples
!  e.g., for a given W, parity, what is the min,max M allowed?
!
!  AUTHOR: CWJ     INITIAL DATE: 9/09
!
! ntuples are computed in bjumplib_ntuples.f90
!

  module ntuple_info
  implicit none

  logical :: restrict_ntuples = .true.  ! flag whether or not to use restrictions

  type partup
    integer jzmin,jzmax
  end type partup

  type wtup
    integer :: parmin,parmax
    type (partup):: par(2)
  end type wtup

  type tuple
     integer :: maxw,minw
     type (wtup), allocatable :: w(:)
  end type tuple

  type (tuple), target :: pairlim(2),tripletlim(2)  
  type (tuple), target :: pnlim, ppnlim, pnnlim

!................................................................
!............. USED WITH MPI TO REDUCE STORAGE...................
!         OR WITH "ALTERNATE" STORAGE OF PN UNCOUPLED TBMES
!         (E.G. WITH NO TRUNCATION IN W )
!................................................................


  type locYtuple
     integer :: maxdWy,mindWy   ! this is all we need to know
  end type locYtuple

  type locXpartup
    integer jzmin,jzmax
    type (locYtuple), allocatable :: Mx(:)
  end type locXpartup

  type locwXtup
    integer :: parmin,parmax
    type (locXpartup):: par(0:1)
  end type locwXtup

  type locXtuple
     integer :: maxdWx,mindWx
     integer :: ndWx
     integer, allocatable :: dWlist(:)
     type (locWXtup), allocatable :: Wx(:)
  end type locXtuple

  type (locXtuple), target, allocatable :: localpnlim(:), localppnlim(:), localpnnlim(:)

!.......... DIFFERENT KIND OF STORAGE...INTRODUCED in 7.3.6...........
!.... store group pairs/triplets/etc organized first by M, then by parity, then by W

!........ADDED in 7.5.9: USE LINKED LISTS TO KEEP TRACK OF GROUPLES
  
!............... LINKED LIST TO KEEP TRACK OF SECTOR JUMPS ASSOCIATED WITH EACH GROUPLE...  
  type sectorjump4grouples
     type (sectorjump4grouples), pointer :: next
     integer :: sjmp
  end type sectorjump4grouples

!................ LINKED LIST TO KEEP TRACK OF UNIQUE GROUPLES..........  
  type grouple
	  type (grouple),pointer :: next
	  integer, allocatable :: cgrops(:)  ! list of creation group operators
	  integer, allocatable :: dgrops(:)  ! list of destruction group operators  
	  integer :: Wdops        !  W values for destruction operator(s) 	  
	  integer :: myindex             ! identifying index
	  integer :: nsjmps                  ! # of sector jumps associated with this grouple
	  type (sectorjump4grouples) :: sj  ! linked list of sector jumps
  end type grouple
   
  type wgtuple
     integer :: ngops   ! # of unique group operators
     type (grouple) :: groop		 
     integer :: minWconj,maxWconj    ! min, max of W of conjugate sectors
     integer :: minWdops,maxWdops
     integer :: maxWdopsum              ! max sum of Wdop for p,n
		 
  end type wgtuple

  type pargtuple
          integer :: maxdWx,mindWx
     integer :: ndWx
     integer, allocatable :: dWlist(:)
     type (wgtuple), allocatable :: Wx(:)
  end type pargtuple

  type mgtuple
    integer :: parmin,parmax
    type (pargtuple):: par(0:1)
  end type mgtuple

  type gtuple
    integer jzmin,jzmax
    type (mgtuple), allocatable :: Mx(:)
  end type gtuple

  type (gtuple), target :: XXgrouples(2), XXXXgrouples(2)

  contains
!
! includes: routines to compute q# limits on 2-body pairs, 3-body triplets
!
!  SUBROUTINES:
!      master_q_ntuple: computes the allowed Q#s for n-tuples (pairs, triplets)
!      ntuple_wlim: compute min, max w for a given descent
!      ntuple_otherlimes: computes limits on jz, w for given descent
!
!      master_cross_ntuple:  master subroutine to compute allowed Q#s
!         for cross-species n-tuples (pairs, triplets)
!         e.g., for pn, ppn, pnn combinations
!      cross_ntuple_wlim
!      cross_ntuple_otherlims
!=====================================================================

! master subroutine to compute the allowed Q#s for n-tuples (pairs, triplets)
!
! uses "descendents" technology to compute this
!
! INPUT:
!     ntple = n of n-tuple (e.g. 2 for pairs, 3 for triplets...
!     it    = species
!
! SUBROUTINES CALLED:
!   sectordescent: creates descendents for each sector
!   ntuple_wlim: compute min, max w for a given descent
!   ntuple_otherlims: computes limits on jz, w for given descent
!---------------------------------------------------------------------
  subroutine master_q_ntuple(ntple, it)
  use system_parameters
  use verbosity
  use sectors
  use descendents

  implicit none
!  include 'binterfaces.inc'

  integer :: ntple
  integer :: it
  type (groupjump) :: descent
  integer :: is
  type (tuple), pointer :: xtuple
  integer :: maxw,minw
  integer :: aerr

  if(.not.restrict_ntuples)return
  if(ntple > np(it))return

  select case (ntple)
  case (2)
    xtuple =>  pairlim(it)
  case (3)
    xtuple => tripletlim(it)
  case default
    print*,ntple, ' -tuple ; error '
    stop
  end select

  descent%nsectors = nsectors(it)
  allocate(descent%sector(nsectors(it)), stat=aerr)
  if(aerr /= 0) call memerror("master_q_tuple 1")
  do is = 1,nsectors(it)
     call sectordescent(it,is,ntple,descent%sector(is))
  enddo ! is

!.............. GET LIMITS ON CHANGE IN W..................
  call ntuple_wlim(it,descent,maxw,minw)

  xtuple%maxw = maxw
  xtuple%minw = minw
!  print*,' For ',ntple,'-body, min/max w = ',minw,maxw
  if(allocated(xtuple%w))deallocate(xtuple%w)
  allocate(xtuple%w(minw:maxw), stat=aerr)
  if(aerr /= 0) call memerror("master_q_tuple 2")
  call ntuple_otherlims(it,descent,xtuple)

!............. CLEAN UP....................
  deallocate(descent%sector)

  return
  end subroutine master_q_ntuple

!=====================================================================
!
! for a given descent of a species it, computes min, max w
!
! INPUT: 
!       it:  species, 1 = proton, 2 = neutron
!       descent: from a given sector
! 
! OUTPUT:
!       maxw, minw
!
  subroutine ntuple_wlim(it,descent,maxw,minw)
  use sporbit
  use sectors
  use blocks
  use descendents
  use butil_mod
  implicit none
  integer it
  type (groupjump) :: descent
  integer is
  integer blocki
  integer irblock,ilblock,frblock,flblock
  integer(8) :: idescend
  integer w,wi,wf
  integer maxw,minw

  minw = 10000
  maxw = 0
  
  if(allsamew)then
     minw = 0
     return
  endif

  do is = 1,nsectors(it)   ! loop over all sectors
     do blocki = 1,descent%sector(is)%nblocks   ! loop over blocks
        irblock = descent%sector(is)%block(blocki)%irblock
        ilblock = descent%sector(is)%block(blocki)%ilblock
        wi = hblock(it)%list(irblock)%wh + hblock(-it)%list(ilblock)%wh
!............ COMPUTE W of initial blocks.................
        do idescend = 1, descent%sector(is)%block(blocki)%ndescendents  ! loop over descendents
           frblock = descent%sector(is)%block(blocki)%frblock(idescend)
           flblock = descent%sector(is)%block(blocki)%flblock(idescend)
!.............COMPUTE W OF FINAL BLOCKS.............
           wf = hblock(it)%list(frblock)%wh + hblock(-it)%list(flblock)%wh
           w = wi - wf
           if( w < 0)then  ! error trap
              print*,' problem with change in w ',wi,wf
              stop
           endif
           minw = bmin(minw,w)
           maxw = bmax(maxw,w)
        enddo  ! idescend
     end do  ! blocki
  enddo  ! is

  return
  end subroutine ntuple_wlim
!=====================================================================
!
! for a given descent of a species it, computes limits on jz, parity,w
!
! 
  subroutine ntuple_otherlims(it,descent,xtuple)
  use sporbit
  use sectors
  use blocks
  use descendents
  use butil_mod
  use basis
  implicit none
  integer it
  type (groupjump) :: descent
  type (tuple), pointer :: xtuple
  integer is
  integer blocki
  integer irblock,ilblock,frblock,flblock
  integer(8) :: idescend
  integer w,wi,wf
  integer dpar,parl,parr,pari,parf
! interface --  integer parmult
  integer djz,jzi,jzf

!---------- SET DEFAULTS FOR MIN,MAX
  do w = xtuple%minw,xtuple%maxw
     xtuple%w(w)%parmin = 3
     xtuple%w(w)%parmax = 0
     xtuple%w(w)%par(1)%jzmin = 10000
     xtuple%w(w)%par(1)%jzmax = -10000
     xtuple%w(w)%par(2)%jzmin = 10000
     xtuple%w(w)%par(2)%jzmax = -10000
  enddo  !
!------------- NOW LOOP AND FIND MIN, MAX..............

  do is = 1,nsectors(it)   ! loop over all sectors
     do blocki = 1,descent%sector(is)%nblocks   ! loop over blocks
        irblock = descent%sector(is)%block(blocki)%irblock
        ilblock = descent%sector(is)%block(blocki)%ilblock
        wi = hblock(it)%list(irblock)%wh + hblock(-it)%list(ilblock)%wh
        parr = hblock(it)%list(irblock)%parh
        parl = hblock(-it)%list(ilblock)%parh
        pari = parmult(parr,parl)
        jzi = hblock(it)%list(irblock)%jzh + hblock(-it)%list(ilblock)%jzh

!............ COMPUTE W of initial blocks.................
        do idescend = 1, descent%sector(is)%block(blocki)%ndescendents  ! loop over descendents
           frblock = descent%sector(is)%block(blocki)%frblock(idescend)
           flblock = descent%sector(is)%block(blocki)%flblock(idescend)
!.............COMPUTE W OF FINAL BLOCKS.............
           wf = hblock(it)%list(frblock)%wh + hblock(-it)%list(flblock)%wh
           parr = hblock(it)%list(frblock)%parh
           parl = hblock(-it)%list(flblock)%parh
           parf = parmult(parr,parl)
           dpar = parmult(parf,pari)
           jzf = hblock(it)%list(frblock)%jzh + hblock(-it)%list(flblock)%jzh
           djz = jzi- jzf
           w = wi - wf

           xtuple%w(w)%parmin = bmin( dpar, xtuple%w(w)%parmin)
           xtuple%w(w)%parmax = bmax( dpar, xtuple%w(w)%parmax)

           xtuple%w(w)%par(dpar)%jzmin = bmin(djz, xtuple%w(w)%par(dpar)%jzmin)
           xtuple%w(w)%par(dpar)%jzmax = bmax(djz, xtuple%w(w)%par(dpar)%jzmax)

        enddo  ! idescend
     end do  ! blocki
  enddo  ! is

  return
  end subroutine ntuple_otherlims
!=====================================================================
!
! subroutine master_cross_ntuple
!
! master subroutine to compute the allowed Q#s for cross-species n-tuples (pairs, triplets)
! e.g., for pn, ppn, pnn combinations
!
! uses "descendents" technology to compute this
!
! IMPORTANT: must call BEFORE reading in files because certain important arrays are created
!
! INPUT:
!     p_tple = n of proton-tuple (e.g. 2 for pairs, 3 for triplets...
!     n_tple = n of neutron-tuple 
!
! SUBROUTINES CALLED:
!     sectordescent:
!     cross_ntuple_wlim
!     cross_ntuple_otherlims
!--------------------------------------------------------------------
  subroutine master_cross_ntuple(p_tple, n_tple)
  use system_parameters
  use verbosity
  use sectors
  use descendents

  implicit none

  integer :: p_tple, n_tple
 
  type (groupjump) :: p_descent,n_descent  ! proton, neutron descents
  integer :: is
  type (tuple), pointer :: xtuple
  integer :: maxw,minw
  logical :: selectok
  integer :: aerr

!  print*,' In CROSS CHECK ',p_tple,n_tple
  if(.not.restrict_ntuples)return
  if(p_tple > np(1))return
  if(n_tple > np(2))return

  selectok = .false.
  if( p_tple == 1 .and. n_tple ==1)then
      xtuple => pnlim
      selectok = .true.
  endif


  if( p_tple == 2 .and. n_tple ==1)then
      xtuple => ppnlim
      selectok = .true.
  endif


  if( p_tple == 1 .and. n_tple ==2)then
      xtuple => pnnlim
      selectok = .true.
  endif

  if(.not.selectok)then
     print*,' Didnt select right values ',p_tple,n_tple
    stop
  endif

  
  p_descent%nsectors = nsectors(1)
  allocate(p_descent%sector(nsectors(1)), stat=aerr)
  if(aerr /= 0) call memerror("master_cross_ntuple 1")
  do is = 1,nsectors(1)
     call sectordescent(1,is,p_tple,p_descent%sector(is))
  enddo ! is

  n_descent%nsectors = nsectors(2)
  allocate(n_descent%sector(nsectors(2)), stat=aerr)
  if(aerr /= 0) call memerror("master_cross_ntuple 2")
  do is = 1,nsectors(2)
     call sectordescent(2,is,n_tple,n_descent%sector(is))
  enddo ! is

!.............. GET LIMITS ON CHANGE IN W..................
  call cross_ntuple_wlim(p_descent,n_descent,maxw,minw)

  xtuple%maxw = maxw
  xtuple%minw = minw
!  print*,' For pn-body, min/max w = ',minw,maxw
  allocate(xtuple%w(minw:maxw), stat=aerr)
  if(aerr /= 0) call memerror("master_cross_ntuple 3")
  call cross_ntuple_otherlims(p_descent,n_descent,xtuple)
!............. CLEAN UP....................
  deallocate(p_descent%sector,n_descent%sector)

  return
  end subroutine master_cross_ntuple
!=====================================================================

  subroutine cross_ntuple_wlim(p_descent,n_descent,maxw,minw)
  use sporbit
  use sectors
  use blocks
  use descendents
  use butil_mod
  implicit none
  type (groupjump) :: p_descent,n_descent
  integer isp,csn, isn
  integer pblocki
  integer irpblock,ilpblock,frpblock,flpblock
  integer(8) :: ipdescend
  integer nblocki
  integer irnblock,ilnblock,frnblock,flnblock
  integer(8) :: indescend
  integer w,wi,wf,wpi,wni,wpf,wnf
  integer maxw,minw

  minw = 10000
  maxw = 0
  
  if(allsamew)then
     minw = 0
     return
  endif

  do isp = 1,nsectors(1)   ! loop over all proton sectors

     do  csn = 1,xsd(1)%sector(isp)%ncsectors                     ! loop over conjugate neutron sectors
         isn = xsd(1)%sector(isp)%csector(csn)
         do pblocki = 1,p_descent%sector(isp)%nblocks   ! loop over proton blocks
            irpblock = p_descent%sector(isp)%block(pblocki)%irblock
            ilpblock = p_descent%sector(isp)%block(pblocki)%ilblock

            wpi = hblock(1)%list(irpblock)%wh + hblock(-1)%list(ilpblock)%wh

            do nblocki = 1,n_descent%sector(isn)%nblocks     ! loop over neutron blocks
               irnblock = n_descent%sector(isn)%block(nblocki)%irblock
               ilnblock = n_descent%sector(isn)%block(nblocki)%ilblock

!............ COMPUTE W of initial blocks.................
               wni = hblock(2)%list(irnblock)%wh + hblock(-2)%list(ilnblock)%wh


               do ipdescend = 1, p_descent%sector(isp)%block(pblocki)%ndescendents  ! loop over proton descendents
                  frpblock = p_descent%sector(isp)%block(pblocki)%frblock(ipdescend)
                  flpblock = p_descent%sector(isp)%block(pblocki)%flblock(ipdescend)
!.............COMPUTE W OF FINAL BLOCKS.............
                  wpf = hblock(1)%list(frpblock)%wh + hblock(-1)%list(flpblock)%wh

                  do indescend = 1, n_descent%sector(isn)%block(nblocki)%ndescendents  ! loop over proton descendents
                     frnblock = n_descent%sector(isn)%block(nblocki)%frblock(indescend)
                     flnblock = n_descent%sector(isn)%block(nblocki)%flblock(indescend)
!.............COMPUTE W OF FINAL BLOCKS.............         
                     wnf = hblock(2)%list(frnblock)%wh + hblock(-2)%list(flnblock)%wh

  
                     w = wpi+wni - wpf-wnf
                     if( w < 0)then  ! error trap
                         print*,' problem with change in w ',wi,wf
                         stop
                      endif
                      minw = bmin(minw,w)
                      maxw = bmax(maxw,w)
                   enddo ! indescend
                enddo  ! ipdescend
       enddo  ! nblocki
     end do  ! pblocki
     enddo ! csn
  enddo  ! isp

  return
  end subroutine cross_ntuple_wlim


!=====================================================================

  subroutine cross_ntuple_otherlims(p_descent,n_descent,xtuple)
  use sporbit
  use sectors
  use blocks
  use descendents
  use butil_mod
  use basis

  implicit none
  type (groupjump) :: p_descent,n_descent
  type (tuple), pointer :: xtuple
  integer isp,csn, isn
  integer pblocki
  integer irpblock,ilpblock,frpblock,flpblock
  integer(8) :: ipdescend
  integer nblocki
  integer irnblock,ilnblock,frnblock,flnblock
  integer(8) :: indescend
  integer w,wi,wf,wpi,wni,wpf,wnf
  integer dpar,parl,parr,pari,parf
  integer parpi,parpf,parni,parnf
!! interface --  integer parmult
  integer djz,jzi,jzf
  integer jzpi,jzpf,jzni,jznf

!---------- SET DEFAULTS FOR MIN,MAX
  do w = xtuple%minw,xtuple%maxw
     xtuple%w(w)%parmin = 3
     xtuple%w(w)%parmax = 0
     xtuple%w(w)%par(1)%jzmin = 10000
     xtuple%w(w)%par(1)%jzmax = -10000
     xtuple%w(w)%par(2)%jzmin = 10000
     xtuple%w(w)%par(2)%jzmax = -10000
  enddo  !

  do isp = 1,nsectors(1)   ! loop over all proton sectors

     do  csn = 1,xsd(1)%sector(isp)%ncsectors                     ! loop over conjugate neutron sectors
         isn = xsd(1)%sector(isp)%csector(csn)
         do pblocki = 1,p_descent%sector(isp)%nblocks   ! loop over proton blocks
            irpblock = p_descent%sector(isp)%block(pblocki)%irblock
            ilpblock = p_descent%sector(isp)%block(pblocki)%ilblock

            wpi = hblock(1)%list(irpblock)%wh + hblock(-1)%list(ilpblock)%wh
            parr = hblock(1)%list(irpblock)%parh
            parl = hblock(-1)%list(ilpblock)%parh
            parpi = parmult(parr,parl)
            jzpi = hblock(1)%list(irpblock)%jzh + hblock(-1)%list(ilpblock)%jzh
            do nblocki = 1,n_descent%sector(isn)%nblocks     ! loop over neutron blocks
               irnblock = n_descent%sector(isn)%block(nblocki)%irblock
               ilnblock = n_descent%sector(isn)%block(nblocki)%ilblock

!............ COMPUTE W of initial blocks.................
               wni = hblock(2)%list(irnblock)%wh + hblock(-2)%list(ilnblock)%wh
               parr = hblock(2)%list(irnblock)%parh
               parl = hblock(-2)%list(ilnblock)%parh
               parni = parmult(parr,parl)
               jzni = hblock(2)%list(irnblock)%jzh + hblock(-2)%list(ilnblock)%jzh

               do ipdescend = 1, p_descent%sector(isp)%block(pblocki)%ndescendents  ! loop over proton descendents
                  frpblock = p_descent%sector(isp)%block(pblocki)%frblock(ipdescend)
                  flpblock = p_descent%sector(isp)%block(pblocki)%flblock(ipdescend)
!.............COMPUTE W OF FINAL BLOCKS.............
                  wpf = hblock(1)%list(frpblock)%wh + hblock(-1)%list(flpblock)%wh
                  parr = hblock(1)%list(frpblock)%parh
                  parl = hblock(-1)%list(flpblock)%parh
                  parpf = parmult(parr,parl)
                  jzpf = hblock(1)%list(frpblock)%jzh + hblock(-1)%list(flpblock)%jzh

                  do indescend = 1, n_descent%sector(isn)%block(nblocki)%ndescendents  ! loop over proton descendents
                     frnblock = n_descent%sector(isn)%block(nblocki)%frblock(indescend)
                     flnblock = n_descent%sector(isn)%block(nblocki)%flblock(indescend)
!.............COMPUTE W OF FINAL BLOCKS.............         
                     wnf = hblock(2)%list(frnblock)%wh + hblock(-2)%list(flnblock)%wh
                     parr = hblock(2)%list(frnblock)%parh
                     parl = hblock(-2)%list(flnblock)%parh
                     parnf = parmult(parr,parl)
                     jznf = hblock(2)%list(frnblock)%jzh + hblock(-2)%list(flnblock)%jzh
  
                     w = wpi+wni - wpf-wnf
                     djz = jzni + jzpi - jznf -jzpf
                     dpar = parmult(parpi,parpf)
                     dpar = parmult(dpar,parni)
                     dpar = parmult(dpar,parnf)

                     xtuple%w(w)%parmin = bmin( dpar, xtuple%w(w)%parmin)
                     xtuple%w(w)%parmax = bmax( dpar, xtuple%w(w)%parmax)

                     xtuple%w(w)%par(dpar)%jzmin = bmin(djz, xtuple%w(w)%par(dpar)%jzmin)
                     xtuple%w(w)%par(dpar)%jzmax = bmax(djz, xtuple%w(w)%par(dpar)%jzmax)

                   enddo ! indescend
                enddo  ! ipdescend
       enddo  ! nblocki
     end do  ! pblocki
     enddo ! csn
  enddo  ! isp
!  print*,' cross limits '
!  print*,xtuple%minw,xtuple%maxw
!  w = xtuple%minw
!  print*,xtuple%w(w)%par(1)%jzmin, xtuple%w(w)%par(1)%jzmax
!  print*,' '

  return
  end subroutine cross_ntuple_otherlims
  
end module ntuple_info
  
