!=================================================================================
!
!  bfragments.f90
!
!  routines to break up Lanczos vectors into "fragments"  
!
!============================================================
!
!  The next module is used for 
!  breaking up the Lanczos vectors 
!
!=========================================================
!
!  FRAGMENTS
!
!   contiguous (sub) sectors of Lanczos vectors grouped together
!   
! started 8/2011 by CWJ
!

module fragments
   use precisions
   use operation_stats
   implicit none

   logical :: fragment_vectors   ! flag to break up vectors into fragments
   integer :: nfragments         ! number of fragments
   logical :: useNewReorthog1 = .true.  ! Use new reorthog even for single fragment
   logical :: useNewReorthog  ! set to (nfragments > 1) .or. useNewReorthog1

   type frag
        integer :: nsectors  ! # of subsectors
        integer(8) :: ssectorstart,ssectorend  ! start, end of sectors; these are defined by "proton" sectors
		  integer(8) :: csectorstart,csectorend  ! start, end of conjugate sectors, defined by "neutron" sectors, must check contiguity
        integer(8):: basisstart,basisend  ! where the basis starts and stops for this sector
        integer(8):: localdim      ! local dimension = basisend-basisstart+1
        real(8)    :: totops    ! tot# of operations to/from this fragment

   end type frag

   type (frag), allocatable :: fragmentlist(:)

!   integer, allocatable :: sectormap2frag(:)  ! maps a sector to a fragment  SCHEDULED FOR OBSOLESCENCE

!--------- THIS ALLOWS US TO STORE JUST FRAGMENTS OF A VECTOR --

   integer (kind=basis_prec), allocatable :: basestart(:), basestop(:)
   
   logical :: test_fragments = .false.  ! added in 7.4.6, used to test new fragment scheme on a single processor
   
!------------- JUMP STORAGE ACROSS FRAGMENTS (f2f matvec blocks)-------------- added in 7.6.2 -----------------

   type f2fjump
	   integer(8) :: nX1bjump(2),nXX2bjump(2),nXXX3bjump(2)
	   real(8)    :: totjumpstorage
	   integer(4) :: minprocs     ! minimal assigned MPI procs/nodes
	   integer(4) :: nprocs       ! # of actual assigned procs in the end
   end type f2fjump   
   
   type (f2fjump), allocatable,target :: f2fjumpstat(:,:)

!------- INFORMATION ON WHAT FRAGMENTS ON WHAT MPI PROCESSES
!        previously in module tribution, moved in 7.6.7
	   
	type whatsonanode
	      integer :: ifragment,ffragment   ! which fragments are here
	      logical :: ifirst ! first node in group servicing ifragment
		  
!	      integer :: sfragment  ! which fragment is stored here; later might want to store more than one
!	      real :: frac_ops   ! what fraction of possible ops are stored here      
	end type whatsonanode

	type (whatsonanode), allocatable :: nodal(:)	   

contains
!======================================================================
!
!  SUBROUTINES CONTAINED IN THIS FILE
!    fragmenter:  master routine to break up Lanczos vectors into fragments
!        which calls:
!    vectorstorageceiling
!    defaultfragments:  routine to set fragment information when only 1 fragment
!    count_create_fragments_blocks: updated fragmentation routine which breaks fragments along haiku blocks
!    reconstruct_sector: when sectors are subdivided, reconstruct proton sectors
!    fragmentvectorsize:  sets base start/stop for a given fragment
!       also in this file (called by distro_opbundles_over_fragments )
!    divide_procs_over_fragments: routine to divide up available nodes over fragments
!
!=======================================================================
!
!  Notes 7.4.6 (March 2015) regarding "new" fragmentation based upon proton haiku blocks
!  To test/debug fragmentation, run on a single processor;
!  in order to do this, we must modify routines which divide up work over MPI processes
!   
!  key routines include: 
!          setnodaltribution  (in bparallel_lib1.f90)  assigns initial, final fragments to MPI procs/compute nodes
!          setup_localvectors (in blanczoslib1.f90)  allocates memory for (fragmented) lanczos vectors
!  If one simulates a fragmented case on a single processor, the above have to be modified
!


!======================================================================
!  subroutine set_nfragments
! Sets nfragments and side effects of assignment
subroutine set_nfragments(nf)
 !  use fragments
   use nodeinfo
   implicit none
   integer :: nf

   nfragments = nf
   useNewReorthog = useNewReorthog1 .or. nfragments > 1
   if(iproc == 0) print *, "useNewReorthog=", useNewReorthog
end subroutine set_nfragments

!======================================================================
!  subroutine fragmenter
!
!  master routine to break up Lanczos vectors into "fragments"
!  when necessary it breaks up sectors into "subsectors"
!  by breaking along the proton basis Slater Determinants
!
!  The "non-optimized" version simply breaks lanczos vectors by sectors,
!  because that is easiest to implement.
!
!  The "block" version breaks up sectors by proton haiku blocks, which is a natural division
!  and preserves most of the downstream routines intact; the divided sectors must be then redefined
!
!  first count up the number of break points
!
!  CALLED BY:
!    main routine
!
!  SUBROUTINES CALLED:
!    vectorstorageceiling   checks to see if vectors get broken up, and how
!    count_create_fragments_noopt
!        [this creates non-optimized fragments by simply bundling together sectors;
!         later will use optimized fragments with subroutine count_create_fragments ]
!       OR
!    count_create_fragments_blocks
!    reconstruct_sector
!    [map_sectors2noopt_fragments] 
!    defaultfragments      only one fragment
!    fragmentvectorsize    sets start and stop for basis in fragments
!    fragmentconjugatesectors  checks which neutron sectors belong (added 7.6.2)
!
subroutine fragmenter
   use flagger
!   use fragments
   use nodeinfo
   use sectors
   use menu_choices

   implicit none

   call vectorstorageceiling  !  check to see if vectors get broken up into fragments
   if(break_vectors .or. test_fragments)then
!	   call test_sector_breaks(83)  ! only for testing
       call count_create_fragments_blocks(.false.)
       call count_create_fragments_blocks(.true.)
       call reconstruct_sectors   
!   call test_sector_breaks(84)  ! only for testing
   else
      call defaultfragments
   endif
   call fragmentvectorsize
   call fragmentconjugatesectors
   return
end subroutine fragmenter

!===============================================================

!
!  VECTORSTORAGECEILING
!
!  computes a ceiling on storing pieces of vectors, if any
!
!  initiated 8/2011 by CWJ @ SDSU
!
!  CALLED BY:
!     fragmenter
!
!  CALLS:
!    fragmentsurveysays

subroutine vectorstorageceiling
   use nodeinfo
   use system_parameters
   use basis
   use menu_choices
   use flagger
   use nodeinfo
   use bmpi_mod
!   use fragments
   use io
   implicit none
   integer(8) :: minmaxfragsize  ! largest block found by splitting along haikus
   integer :: ierr
   character :: dummchar
   logical   :: interrflag  ! error in reading an integer

!...... parameters ask2break_vectors and maxfragmentsize_default found in module flagger..
!... MODIFIED IN 7.6.5 to ask for dummy fragmentsize even if not in MPI...

   if(ask2break_vectors)then  ! should always default to this...
	   if(nproc > 1 .or. menu_char=='m' .or. test_fragments)then
         if(iproc==0)then
		  call fragmentsurveysays(minmaxfragsize)
          print*,' '
          print*,' Enter desired limit on fragment size for breaking Lanczos vectors'
		  print*,' (Largest un-splittable block = ',minmaxfragsize,')'
		  if(minmaxfragsize < maxfragmentsize_default)then
             print*,'( Default = ',maxfragmentsize_default,')'
 		  else
			 print*,'( Default = ',minmaxfragsize,')'
!... PUT IN HINT IF Z /= N.....
             if( np(1)/= np(2))print*,' You might be able to reduce this by swapping (valence) protons and neutrons '
  		  end if
          print*,'( Enter 0 to use default )'
          read(5,*,err=1111)maxfragmentsize
          if(maxfragmentsize < 1 .and. minmaxfragsize <= maxfragmentsize_default)maxfragmentsize = maxfragmentsize_default
          if(maxfragmentsize < 1 .and. minmaxfragsize > maxfragmentsize_default)maxfragmentsize = minmaxfragsize
		  
          print*,' Using target Lanczos fragment dimension of ',maxfragmentsize
          print*,' '
		  write(autoinputfile,'(i12,"    !  FRAGMENT SIZE (0 = use default)")')maxfragmentsize
		  
         end if
         call BMPI_Bcast(maxfragmentsize,1,0,icomm,ierr)
		 
       else
!..... ADDED 7.6.5......NEED TO LEAVE ROOM FOR FRAGMENTSIZE IN AUTOINPUT
         if(auto_input)then
	         read(autoinputfile,'(a)')dummchar
		 else
			 write(autoinputfile,'("0     !  LANCZOS FRAGMENT SIZE (0 = use default)")')
		 end if
         maxfragmentsize=maxfragmentsize_default
	   
       end if
 
   else
      maxfragmentsize=maxfragmentsize_default

   end if

   if(break_vectors_enabled .and. dimbasis > maxfragmentsize)then
       if(menu_char/='m' .and. nproc ==1  .and. .not.test_fragments)then
          break_vectors=.false.
!		  if(.not.test_fragments)then
             print*,' Only one process; cannot break vector '
             print*,' You may run out of memory '
!		 end if
       else
          break_vectors = .true.
          if(iproc==0)then
             print*,' '
             print*,' Basis dimension greater than limit of ',maxfragmentsize
             print*,' Breaking up into fragments '
             print*,' (NOTE: You can set parameter maxfragmentsize in file bmodule_flags.f90 ) '
             print*,' '
          end if
        end if
   else
       break_vectors = .false.
   end if
 
   return
1111 continue
   if(iproc==0)then
	   print*,' ERROR you did not set the fragmentsize ERROR'
	   write(resultfile,*)' ERROR you did not set the fragmentsize ERROR'
	   write(logfile,*)' ERROR you did not set the fragmentsize ERROR'
	   close(resultfile)
	   close(logfile)
   endif
   call BMPI_ABORT(icomm,101,ierr)   
   stop
end subroutine vectorstorageceiling
!===============================================================
!
! routine to survery fragment sizes
! returns the largest haiku-limited fragment
!  added 7.5.0
!
! called by vectorstorageceiling
!
subroutine fragmentsurveysays(maxfragsize)
    use sectors
!    use fragments
    use flagger
    use nodeinfo
    use basis	
    use butil_mod

	implicit none
    integer(8) :: maxfragsize
    integer isector
    integer(8) :: npsd, nnsd
    integer(8) :: localpsd  ! how many proton SDs are described by a (right, left) pair of haiku blocks
    integer ic,ics
    integer i
    integer :: nphblocks  ! # of proton haiku blocks	
	
	maxfragsize=0
    do isector = 1,nsectors(1)  ! loop over proton sectors
       nphblocks = xsd(1)%sector(isector)%nhblocks  ! # of proton haiku blocks in this sector
 !.............. COUNT UP HOW MANY CONJUGATE NEUTRON SDs THERE ARE    
       nnsd = 0
       do ic = 1,xsd(1)%sector(isector)%ncsectors   ! LOOP OVER CONJUGATE NEUTRON SECTORS
          ics = xsd(1)%sector(isector)%csector(ic)
          nnsd = nnsd + xsd(2)%sector(ics)%nxsd
       end do   ! ic 
 !.............. COUNT UP ALONG PROTON HAIKU BLOCKS..................

       do i = 1,nphblocks
          localpsd = xsd(1)%sector(isector)%blockend(i)-xsd(1)%sector(isector)%blockstart(i)+1
		  maxfragsize=bmax(maxfragsize,localpsd*nnsd)
	  end do ! i
  end do ! isector
  return
	
end subroutine fragmentsurveysays

!===============================================================
!
!    defaultfragments
! routine to set fragment information when only 1 fragment
!
! CALLED BY
!    fragmenter
!
subroutine defaultfragments

   use sectors
!   use fragments
   use basis
   use nodeinfo
   implicit none
   integer :: aerr

   if(iproc==0)print*,' defaulting to 1 Lanczos vector fragment (no break-up) '
   call set_nfragments(1)
   if(allocated(fragmentlist))deallocate(fragmentlist)
   allocate(fragmentlist(1), stat=aerr)
   if(aerr /= 0) call memerror("defaultfragments")
   fragmentlist(1)%nsectors = nsectors(1)
   fragmentlist(1)%ssectorstart = 1
   fragmentlist(1)%ssectorend  = nsectors(1)
   fragmentlist(1)%basisstart = 1
   fragmentlist(1)%basisend   = dimbasis
   fragmentlist(1)%localdim   = dimbasis
   
   return
end subroutine defaultfragments

!======================================================================

!
!  this version optimizes the distribution by counting up  proton haiku blocks
!  started 7.4.1
!
!  This approach will minimize the changes needed in other routines downstream
!  Note: any sector is constructed by pairs of "left" and "right" "blocks" of "haikus"
!    -- a left haiku is a slater determinant constructed only by s.p. states with jz < 0
!    -- a right haiku is a slater determinant constructed only by s.p. states with jz >= 0
!    -- a block is a set of haikus with the same abelian quantum numbers
!  
!  Why this is useful, briefly (see detailed routines in bjumplib1,2,3....)
!  Jumps are created from geneologies, which are created from descendents, which descend from blocks
!  By splitting along blocks, all the jump-creation technology will continue to be used unchanged
!
!  INPUT:
!    create: logical flag to either count (create=F) or create (T)
!
!  CALLED BY: 
!     fragmenter (above)
!
subroutine count_create_fragments_blocks(create)
   use sectors
!   use fragments
   use flagger
   use nodeinfo
   use basis
   use butil_mod
   use basis

   implicit none
   logical create
   integer isector
   integer(8) :: localbasiscount
   integer(8) :: nnsd
   integer(8) :: localpsd  ! how many proton SDs are described by a (right, left) pair of haiku blocks
   integer(8) :: ic,ics
   integer(8) :: psd_splitsize
   integer :: nfrag        ! # of fragments
   integer :: nss          ! # of subsectors
   integer(8) :: minfragsize,maxfragsize
   integer(8) :: sumbasis
   integer :: nsplit       ! temporary # of splits in a sector
   integer :: i, nphblocks  ! # of proton haiku blocks
   integer :: aerr
   integer :: ii
   integer :: udbg

   logical, parameter :: writefrag = .true.

   udbg = 0
   if(writefrag .and. iproc==0 .and. create) udbg = 124;
   if(udbg .ne. 0) then
      open(unit=udbg,file='fraginfo.bigstick',status='unknown')
      if(create)write(udbg,*)' # of fragments ',nfragments
      if(create)write(udbg,*)' fragment label     /       dim '
   end if

   nfrag = 1
   nss = 0   
   localbasiscount = 0
   maxfragsize = 0
   sumbasis = 0
   minfragsize = maxfragmentsize
   
   if(iproc==0 .and. .not.create)print*,' Improved fragments (splitting via proton haiku blocks) '
   
   do isector = 1,nsectors(1)  ! loop over proton sectors
      nphblocks = xsd(1)%sector(isector)%nhblocks  ! # of proton haiku blocks in this sector
      if(create)then
         nsplit =  xsd(1)%sector(isector)%ncuts
         allocate(xsd(1)%sector(isector)%cut(nsplit), stat=aerr)
         if(aerr /= 0) call memerror("count_create_fragments_blocks 1")
      endif
      nsplit = 0 
      nss = nss+1   ! start a new sector
!.............. COUNT UP HOW MANY CONJUGATE NEUTRON SDs THERE ARE    
      nnsd = 0
      do ic = 1,xsd(1)%sector(isector)%ncsectors   ! LOOP OVER CONJUGATE NEUTRON SECTORS
         ics = xsd(1)%sector(isector)%csector(ic)
         nnsd = nnsd + xsd(2)%sector(ics)%nxsd
      end do   ! ic 
!.............. COUNT UP ALONG PROTON HAIKU BLOCKS..................

      do i = 1,nphblocks
         localpsd = xsd(1)%sector(isector)%blockend(i)-xsd(1)%sector(isector)%blockstart(i)+1
         if(localbasiscount + nnsd*localpsd > maxfragmentsize .and. i < nphblocks)then  ! reset
                                                            ! create a new fragment, split off a new sector
            if(i==1) then   ! close out fragment at prior sector 
               if(localbasiscount==0)then
                  print*,' SOME PROBLEM WITH FRAGMENTING '
                  print*,' FAILING IN FIRST SECTOR '
                  print*,localbasiscount,nnsd*localpsd,maxfragmentsize
                  print*,isector,nfrag
                  stop
  				   end if
               if(create)then  ! fill current fragment
					   if(nfrag==1)then
                     fragmentlist(nfrag)%basisstart = 1
                     fragmentlist(nfrag)%basisend = localbasiscount
                     fragmentlist(nfrag)%ssectorstart = 1
                  else
                     fragmentlist(nfrag)%basisstart = fragmentlist(nfrag-1)%basisend+1
                     fragmentlist(nfrag)%basisend = localbasiscount+fragmentlist(nfrag-1)%basisend
                     fragmentlist(nfrag)%ssectorstart = 1+fragmentlist(nfrag-1)%ssectorend
                  end if               
                  fragmentlist(nfrag)%ssectorend   = nss - 1 ! prev sector in this case  bug  fix KSM 7.6.2
               end if				 
            else 												
               nsplit = nsplit+1
               if(create)then  ! fill current fragment
                  fragmentlist(nfrag)%localdim = localbasiscount
                  xsd(1)%sector(isector)%cut(nsplit)= i-1
                  if(nfrag==1)then
                      fragmentlist(nfrag)%basisstart = 1
                      fragmentlist(nfrag)%basisend = localbasiscount
                      fragmentlist(nfrag)%ssectorstart = 1
                  else
                      fragmentlist(nfrag)%basisstart = fragmentlist(nfrag-1)%basisend+1
                      fragmentlist(nfrag)%basisend = localbasiscount+fragmentlist(nfrag-1)%basisend
                      fragmentlist(nfrag)%ssectorstart = 1+fragmentlist(nfrag-1)%ssectorend
                  end if               
                  fragmentlist(nfrag)%ssectorend   = nss
               end if
               nss = nss+1      ! start new subsector
            end if
            maxfragsize = bmax(maxfragsize, localbasiscount)
            minfragsize = bmin(minfragsize, localbasiscount)
            nfrag = nfrag +1  ! start new fragment
            sumbasis = sumbasis+localbasiscount    ! update basis so far
            if(udbg /= 0)then
               write(udbg,*)nfrag-1,localbasiscount
            end if
            localbasiscount = 0
         end if
         localbasiscount = localbasiscount+nnsd*localpsd
         xsd(1)%sector(isector)%ncuts = nsplit
      end do

      ! If we just finished the last sector
      if(isector==nsectors(1))then
         maxfragsize = bmax(maxfragsize, localbasiscount)
         minfragsize = bmin(minfragsize, localbasiscount)
         sumbasis = sumbasis+localbasiscount
         if(create)then
            if(nfrag==1)then
               fragmentlist(nfrag)%basisstart = 1
               fragmentlist(nfrag)%basisend = localbasiscount
               fragmentlist(nfrag)%ssectorstart = 1
            else
               fragmentlist(nfrag)%basisstart = fragmentlist(nfrag-1)%basisend+1
               fragmentlist(nfrag)%basisend = localbasiscount+fragmentlist(nfrag-1)%basisend
               fragmentlist(nfrag)%ssectorstart = 1+fragmentlist(nfrag-1)%ssectorend
            end if
            fragmentlist(nfrag)%ssectorend   = nss
         end if
      end if
   end do ! isector
   if(create)then
      do ii=1,nfrag
         fragmentlist(ii)%localdim = fragmentlist(ii)%basisend - fragmentlist(ii)%basisstart + 1
      end do
   end if
   if(iproc==0 .and. .not. create)then
      print*,' There are a total of ',nfrag,' Lanczos fragments'
      print*,' from a new total of ',nss,' sectors '
      print*,' min/max fragment dimension = ',minfragsize,maxfragsize
   end if
   if(udbg /= 0) then ! KSM
      write(udbg, *) "DEBUG - fragment dump1"
      do ii = 1, nfragments
         write(udbg, "(I3,A,I12,A,I12,A)") ii, ":(", fragmentlist(ii)%basisstart, ",", fragmentlist(ii)%basisend, ")"
         write(udbg, "(A,I0,A,I0)") "   sectorstart=", fragmentlist(ii)%ssectorstart, ", sectorend=", fragmentlist(ii)%ssectorend
         write(udbg, "(A,I0)")      "   nsectors=", fragmentlist(ii)%nsectors
      end do
   end if
   if(udbg /= 0) then
      close(udbg) 
      udbg = 0
   end if
   if(sumbasis/=dimbasis)then
      print*,' Mismatch in basis dimension '
      stop
   end if

   if(.not.create)then
      call set_nfragments(nfrag)
      if(.not.allocated(fragmentlist)) then
         allocate(fragmentlist(nfragments) , stat=aerr)
         if(aerr /= 0) call memerror("count_create_fragments_blocks 10")
      end if
   end if
   return
end subroutine count_create_fragments_blocks

!======================================================================
!
! when sectors are subdivided, reconstruct proton sectors
! added 7.4.1; used after count_create_fragments_blocks
!
!  CALLED BY:
!    fragmenter

subroutine reconstruct_sectors
   use sectors
   use basis
   use nodeinfo
   implicit none
   type(mastersect) :: draftpsd
   integer nnewsectors
   integer :: isector,ics,jsector,jcs,ncs
   integer :: icut
   integer :: nblocks,iblock,jblock
   integer :: aerr
   integer(8) :: npsd

!......... COUNT UP NEW SECTORS
   nnewsectors = 0
   do isector = 1,nsectors(1)
      nnewsectors = nnewsectors + xsd(1)%sector(isector)%ncuts + 1
   end do
   if(nnewsectors == nsectors(1)) return   ! nothing to do

!........ CREATE TEMPORARY DRAFT OF SECTORS......
   allocate( draftpsd%sector(nnewsectors) , stat=aerr)
   if(aerr /= 0) call memerror("reconstruct_sectors 1")
   
!............ ISECTOR labels old sectors, while NNEWSECTORS labels new sectors......
   
   nnewsectors = 0
   do isector = 1,nsectors(1)
      nnewsectors = nnewsectors+1   ! INCREMENT # OF SECTORS

      draftpsd%sector(nnewsectors)%xsdstart = xsd(1)%sector(isector)%xsdstart   ! SET BASIC QUANTUM NUMBERS
      draftpsd%sector(nnewsectors)%jzX      = xsd(1)%sector(isector)%jzX
      draftpsd%sector(nnewsectors)%parX     = xsd(1)%sector(isector)%parX
      draftpsd%sector(nnewsectors)%WX       = xsd(1)%sector(isector)%wX
!      draftpsd%sector(nnewsectors)%parentsector = isector
!................ COPY OVER CONJUGATE SECTORS...................
      draftpsd%sector(nnewsectors)%ncsectors= xsd(1)%sector(isector)%ncsectors
      allocate ( draftpsd%sector(nnewsectors)%csector( draftpsd%sector(nnewsectors)%ncsectors ), stat=aerr)
      if(aerr /= 0) call memerror("reconstruct_sectors 10")
      do ics = 1, xsd(1)%sector(isector)%ncsectors
         draftpsd%sector(nnewsectors)%csector(ics) = xsd(1)%sector(isector)%csector(ics) 
      end do
      
      if( xsd(1)%sector(isector)%ncuts > 0)then
         do icut = 1, xsd(1)%sector(isector)%ncuts
!............. BEFORE STARTING NEW "SECTOR" FINISH DATA FOR CURRENT ONE........
            nblocks  = xsd(1)%sector(isector)%cut(icut) 
            if(icut/=1)then
               nblocks = nblocks - xsd(1)%sector(isector)%cut(icut-1)
            end if
            draftpsd%sector(nnewsectors)%nhblocks  = nblocks
            draftpsd%sector(nnewsectors)%xsdend   = xsd(1)%sector(isector)%blockend( xsd(1)%sector(isector)%cut(icut)  )
            draftpsd%sector(nnewsectors)%nxsd = draftpsd%sector(nnewsectors)%xsdend - draftpsd%sector(nnewsectors)%xsdstart+1
            allocate (draftpsd%sector(nnewsectors)%rhblock(nblocks), draftpsd%sector(nnewsectors)%lhblock(nblocks), stat=aerr)
            if(aerr /= 0) call memerror("reconstruct_sectors 20")
            allocate (draftpsd%sector(nnewsectors)%blockstart(nblocks), draftpsd%sector(nnewsectors)%blockend(nblocks), stat=aerr)
            if(aerr /= 0) call memerror("reconstruct_sectors 21")

            do iblock = 1,nblocks
               jblock = iblock 
               if(icut > 1)jblock = jblock + xsd(1)%sector(isector)%cut(icut-1)
               draftpsd%sector(nnewsectors)%rhblock(iblock) = xsd(1)%sector(isector)%rhblock(jblock)
               draftpsd%sector(nnewsectors)%lhblock(iblock) = xsd(1)%sector(isector)%lhblock(jblock)
               draftpsd%sector(nnewsectors)%blockstart(iblock) = xsd(1)%sector(isector)%blockstart(jblock)
               draftpsd%sector(nnewsectors)%blockend(iblock) = xsd(1)%sector(isector)%blockend(jblock)
            end do
            draftpsd%sector(nnewsectors)%nxsd = draftpsd%sector(nnewsectors)%blockend(nblocks)- & 
                                                draftpsd%sector(nnewsectors)%blockstart(1)+1
            draftpsd%sector(nnewsectors)%xsdend = draftpsd%sector(nnewsectors)%xsdstart + draftpsd%sector(nnewsectors)%nxsd -1
!................... START NEW SECTOR..........................
            nnewsectors = nnewsectors+1
			draftpsd%sector(nnewsectors)%xsdstart= draftpsd%sector(nnewsectors-1)%xsdend +1
            draftpsd%sector(nnewsectors)%jzX      = xsd(1)%sector(isector)%jzX    ! SET BASIC QUANTUM NUMBERS
            draftpsd%sector(nnewsectors)%parX     = xsd(1)%sector(isector)%parX
            draftpsd%sector(nnewsectors)%WX       = xsd(1)%sector(isector)%wX
!           draftpsd%sector(nnewsectors)%parentsector = isector
!................ COPY OVER CONJUGATE SECTORS...................
            draftpsd%sector(nnewsectors)%ncsectors= xsd(1)%sector(isector)%ncsectors
            allocate ( draftpsd%sector(nnewsectors)%csector( draftpsd%sector(nnewsectors)%ncsectors ), stat=aerr)
            if(aerr /= 0) call memerror("reconstruct_sectors 30")
            do ics = 1, xsd(1)%sector(isector)%ncsectors
                draftpsd%sector(nnewsectors)%csector(ics) = xsd(1)%sector(isector)%csector(ics) 
            end do   
         end do
      end if
!............. FINISH DATA FOR CURRENT "SECTOR"...............

      draftpsd%sector(nnewsectors)%xsdend        = xsd(1)%sector(isector)%xsdend
      draftpsd%sector(nnewsectors)%nxsd = draftpsd%sector(nnewsectors)%xsdend - draftpsd%sector(nnewsectors)%xsdstart+1
      nblocks = xsd(1)%sector(isector)%nhblocks
      if(xsd(1)%sector(isector)%ncuts /=0)then
         nblocks = nblocks - xsd(1)%sector(isector)%cut(xsd(1)%sector(isector)%ncuts)
      end if
      draftpsd%sector(nnewsectors)%nhblocks  = nblocks
      allocate (draftpsd%sector(nnewsectors)%rhblock(nblocks), draftpsd%sector(nnewsectors)%lhblock(nblocks), stat=aerr)
      if(aerr /= 0) call memerror("reconstruct_sectors 40")
      allocate (draftpsd%sector(nnewsectors)%blockstart(nblocks), draftpsd%sector(nnewsectors)%blockend(nblocks), stat=aerr)
      if(aerr /= 0) call memerror("reconstruct_sectors 41")
      do iblock = 1,nblocks
               jblock = iblock 
               if(xsd(1)%sector(isector)%ncuts > 0)jblock = jblock + xsd(1)%sector(isector)%cut(xsd(1)%sector(isector)%ncuts)
               draftpsd%sector(nnewsectors)%rhblock(iblock) = xsd(1)%sector(isector)%rhblock(jblock)
               draftpsd%sector(nnewsectors)%lhblock(iblock) = xsd(1)%sector(isector)%lhblock(jblock)
               draftpsd%sector(nnewsectors)%blockstart(iblock) = xsd(1)%sector(isector)%blockstart(jblock)
               draftpsd%sector(nnewsectors)%blockend(iblock) = xsd(1)%sector(isector)%blockend(jblock)
      end do

   end do  ! isector
!...... TEST THE SDs ARE ALIGNED............
   npsd = 0
   do isector = 1,nnewsectors
      npsd = npsd+ draftpsd%sector(isector)%xsdend - draftpsd%sector(isector)%xsdstart+1
	  if(draftpsd%sector(isector)%xsdend - draftpsd%sector(isector)%xsdstart+1 /= draftpsd%sector(isector)%nxsd )then
	     if(iproc==0)print*,' mismatch in new sectors ', isector
		 stop
	  end if
	  if(isector > 1)then
		  if( draftpsd%sector(isector)%xsdstart /= draftpsd%sector(isector-1)%xsdend+1)then
			  if(iproc==0)then
				  print*,' Mismatch in start/stop of proton sds ',isector
				  print*,draftpsd%sector(isector)%xsdstart, draftpsd%sector(isector-1)%xsdend+1
			  end if
			  stop
		  end if
	  end if
	  
   end do
   if(npsd /= nxsd(1))then
	   if(iproc==0)print*,' Mismatch in new # of proton SDs ',npsd,nxsd(1)
	   stop
   end if
!........ CLEAR OUT OLD SECTORS 
   do isector = 1,nsectors(1)
      nullify(xsd(1)%sector(isector)%csector)
      nullify(xsd(1)%sector(isector)%rhblock,xsd(1)%sector(isector)%lhblock )
      nullify(xsd(1)%sector(isector)%blockstart,xsd(1)%sector(isector)%blockend )
      nullify(xsd(1)%sector(isector)%cut)
   end do
   nullify(xsd(1)%sector)

!........ COPY DRAFT BACK INTO SECTORS
   nsectors(1) = nnewsectors
   allocate( xsd(1)%sector(nnewsectors) , stat=aerr)
   if(aerr /= 0) call memerror("reconstruct_sectors 50")
   do isector = 1,nsectors(1)
      xsd(1)%sector(isector)%xsdstart = draftpsd%sector(isector)%xsdstart
      xsd(1)%sector(isector)%xsdend   = draftpsd%sector(isector)%xsdend
      xsd(1)%sector(isector)%nxsd     = draftpsd%sector(isector)%nxsd
      xsd(1)%sector(isector)%jzX      = draftpsd%sector(isector)%jzX
      xsd(1)%sector(isector)%parX     = draftpsd%sector(isector)%parX
      xsd(1)%sector(isector)%WX       = draftpsd%sector(isector)%wX
      xsd(1)%sector(isector)%ncsectors= draftpsd%sector(isector)%ncsectors
      xsd(1)%sector(isector)%nhblocks = draftpsd%sector(isector)%nhblocks
!      xsd(1)%sector(isector)%parentsector = draftpsd%sector(isector)%parentsector
      allocate (xsd(1)%sector(isector)%csector(  draftpsd%sector(isector)%ncsectors ) , stat=aerr)
      if(aerr /= 0) call memerror("reconstruct_sectors 60")
      allocate (xsd(1)%sector(isector)%lhblock(  draftpsd%sector(isector)%nhblocks ) , stat=aerr)
      if(aerr /= 0) call memerror("reconstruct_sectors 61")
      allocate (xsd(1)%sector(isector)%rhblock(  draftpsd%sector(isector)%nhblocks ) , stat=aerr)
      if(aerr /= 0) call memerror("reconstruct_sectors 62")
      allocate (xsd(1)%sector(isector)%blockstart(  draftpsd%sector(isector)%nhblocks ) , stat=aerr)
      if(aerr /= 0) call memerror("reconstruct_sectors 63")
      allocate (xsd(1)%sector(isector)%blockend(  draftpsd%sector(isector)%nhblocks ) , stat=aerr)
      if(aerr /= 0) call memerror("reconstruct_sectors 64")
      do ics = 1, draftpsd%sector(isector)%ncsectors 
          xsd(1)%sector(isector)%csector(ics) = draftpsd%sector(isector)%csector(ics)
      end do
      do ics = 1, draftpsd%sector(isector)%nhblocks
          xsd(1)%sector(isector)%rhblock(ics) = draftpsd%sector(isector)%rhblock(ics)
          xsd(1)%sector(isector)%lhblock(ics) = draftpsd%sector(isector)%lhblock(ics)
          xsd(1)%sector(isector)%blockstart(ics) = draftpsd%sector(isector)%blockstart(ics)
          xsd(1)%sector(isector)%blockend(ics) = draftpsd%sector(isector)%blockend(ics)
      end do  ! ics
   end do  ! isector

!........ ASSIGN NEW CONJUGATE SECTORS TO NEUTRONS

   do jsector = 1,nsectors(2)
      ncs = 0
      do isector = 1,nsectors(1)
         do ics = 1,xsd(1)%sector(isector)%ncsectors
            if(xsd(1)%sector(isector)%csector(ics) == jsector)then 
               ncs = ncs+1
               exit
            end if
         end do ! ics
      end do
      xsd(2)%sector(jsector)%ncsectors = ncs
      if(ncs ==0)cycle
      nullify( xsd(2)%sector(jsector)%csector )
      allocate ( xsd(2)%sector(jsector)%csector(ncs) , stat=aerr)
      if(aerr /= 0) call memerror("reconstruct_sectors 70")
      ncs = 0
      do isector = 1,nsectors(1)
         do ics = 1,xsd(1)%sector(isector)%ncsectors
            if(xsd(1)%sector(isector)%csector(ics) == jsector)then 
                ncs = ncs+1
                xsd(2)%sector(jsector)%csector(ncs) = isector
                exit
            end if
         end do ! ics
      end do

   end do  ! jsector
   return
end subroutine reconstruct_sectors
!======================================================================
!
!  check new (7.4.1) break up of sectors for fragments is okay
!  writes to file fort.unitout; run two different ways to compare
!
subroutine test_sector_breaks(unitout)
   use sectors
   use basis

   implicit none
   integer unitout
   integer(8) :: basiscount
   integer    :: ips,ins   
   integer    :: pblock,nblock
   integer    :: ncs   ! # of conjugate sectors
   integer    :: ics
   integer(8) :: pblocksize, nblocksize,nlocal
   basiscount = 0
   if(unitout ==83)then
       open(unit=83,file='sector.in',status='unknown')
   else
       open(unit=84,file='sector.out',status='unknown')

   endif
   do ips = 1,nsectors(1)  ! loop over proton sectors
      ncs = xsd(1)%sector(ips)%ncsectors   ! fetch the number of conjugate neutron sectors
	  nlocal = 0
      do pblock = 1,xsd(1)%sector(ips)%nhblocks  ! within this sector, loop over the pairs of proton haiku blocks
          pblocksize = int(xsd(1)%sector(ips)%blockend(pblock) - xsd(1)%sector(ips)%blockstart(pblock)+1,8)  
                 ! fetch the proton contribution 
		  nlocal = nlocal+pblocksize		
		  write(unitout,*)pblock,pblocksize
          do ics = 1,ncs
              ins = xsd(1)%sector(ips)%csector(ics)
              do nblock = 1,xsd(2)%sector(ins)%nhblocks
                 nblocksize = int(xsd(2)%sector(ins)%blockend(nblock) - xsd(2)%sector(ins)%blockstart(nblock)+1,8)  
                 basiscount = basiscount + pblocksize*nblocksize
!				 nlocal = nlocal + pblocksize*nblocksize
!                 write(unitout,*)basiscount


              end do ! nblock

          end do !ics

      end do ! pblock
	  write(unitout,*)ips,' proton sector SDs ',nlocal, xsd(1)%sector(ips)%nxsd


   end do
   write(unitout,*)basiscount,dimbasis
   if(basiscount /= dimbasis)then  ! error flag
      print*, '  v v v v v v v v v v v v v v v v v v v v v '
      print*,' Oh noes! Mismatch in basis ! '
      print*,dimbasis,basiscount
      print*, '  ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ '
      stop
   end if
   close(unitout)
   return

end subroutine test_sector_breaks
!======================================================================

!
!  fragmentvectorsize
!     sets base start/stop for a given fragment
!
! CALLED BY:
!  fragmenter
!
subroutine fragmentvectorsize
!   use fragments
   implicit none
   integer :: ifrag
   integer :: aerr

   allocate(basestart(nfragments), basestop(nfragments) , stat=aerr)
   if(aerr /= 0) call memerror("fragmentvectorsize")

   do ifrag = 1,nfragments
       basestart(ifrag) = fragmentlist(ifrag)%basisstart
       basestop(ifrag)  = fragmentlist(ifrag)%basisend

   end do
   return

end subroutine fragmentvectorsize

!======================================================================

!
!  printopfragstat
!     sets base start/stop for a given fragment
!
! CALLED BY:
!  fragmenter
!
subroutine printopfragstat
   use nodeinfo
   use operation_stats
!   use fragments
   implicit none
   integer :: ifrag, jfrag

   if(iproc==0)then
	   do ifrag = 1,nfragments
		   do jfrag=1,nfragments
			   print*,ifrag,jfrag,' nodes: ', opfragstat(ifrag,jfrag)%nnodes
		   end do
	   end do
   end if   
end subroutine printopfragstat
!===============================================================
!
!  initiated July 2014 @ LLNL by CWJ
!
!  routines to guide distributing work across and between 
!  Lanczos vectors "fragments."
!
!==========================================================================================
!   divide_procs_over_fragments
!
! routine to divide up available nodes over fragments
!
! main goal is to set up opfragstat(if,ff)%nnodes 
! which assign nodes to initial, final fragments if,ff
!
! INPUT:
!    nnodesx  : sample # of MPI processes 
!                This might differ the actual # of nodes available, 
!                if one needs to reserve nodes
!    assignnodes: flag to actually assign nodes
!
! CALLED BY:
!   distro_opbundles_over_fragments
!
subroutine divide_procs_over_fragments(nnodesx,assignnodes)
!   use fragments
   use operation_stats
   use nodeinfo
!   use tribution
   use verbosity
   use butil_mod
   implicit none
!.......INPUT........
   integer :: nnodesx  !  # of nodes made available
   logical :: assignnodes  ! flag to assign nodes

!....... INTERNAL............
   real(8) :: opfrac      ! frac of operations
   integer :: ifrag, jfrag,i,j
   real(8) :: optotall,opavg,optotsingle
   integer :: nodesum, nodediff,nodesingle,nodeadd
   real(8) :: diffcf
   ! KSM  - just use ifrag, jfrag for position:   integer :: inodecf
   integer :: si, sj  ! used for search for biggest discrepency
!.... ESTIMATE EFFICIENCY.....
   real(8) :: maxnodeops, avgnodeops

!...... KEEP TRACK OF WHICH FRAGMENT BLOCKS HAVE BIGGEST DISCREPANCIES....
   ! KSM:  real(8) :: opdisc(nfragments*nfragments) 
   real(8) :: opdisca(nfragments, nfragments)
   logical :: localverbose  = .false.   ! used for debugging

   real :: assignratio = 0.9
      
!........... DEFAULT IF NO FRAGMENTS.........
   
   opdisca = 0.0  ! initialize
 
!............. COMPUTE TOTAL # OF OPS...............
   optotall = 0.0d0

   do ifrag = 1,nfragments
      do jfrag = 1,nfragments
         ! shan: using new weight op
         if (iproc == 0 .and. assignnodes) write(*, 111) "FRAGOP ", ifrag, jfrag, " stateop ", opfragstat(ifrag,jfrag)%optot, &
              " fragop ", frag2fragops(ifrag, jfrag)
         
         opfragstat(ifrag,jfrag)%optot = frag2fragops(ifrag, jfrag)
          optotall = optotall + opfragstat(ifrag,jfrag)%optot
   end do
   end do

111 format(a7, i6, i6, a9, f16.2, a8, f16.2)
   
   opavg = optotall/real(nnodesx,8)
!!............ IF ALREADY ASSIGNED.................   

   if(allocated(nodal))then
!	   if(iproc==0)print*,' nodes already assigned '
!	   if(iproc==0)print*,opfragstat(:,:)%nnodes
	   return
	   
   end if   
!   print*,' operation totals for TESTING ',optotall,opavg
!........... COUNT THE NUMBER OF SINGLE-NODE FRAGMENT BLOCKS.......
   nodesingle= 0
   optotsingle = 0.0d0
   do ifrag = 1,nfragments
      do jfrag = 1,nfragments
          if(opfragstat(ifrag,jfrag)%optot < opavg .and. opfragstat(ifrag,jfrag)%optot> 0.0d0)then
              opfragstat(ifrag,jfrag)%nnodes=1
              ! KSM - already set to 0:   opdisc(nfragments*(ifrag-1)+jfrag) = 0.0d0
              nodesingle = nodesingle+1
              optotsingle=optotsingle + opfragstat(ifrag,jfrag)%optot
          else
            opfragstat(ifrag,jfrag)%nnodes=0
          end if
      end do  
   end do

   if (nnodesx * assignratio -nodesingle < 0) then
      assignratio = 1.0
   endif
   
   if(iproc==0 .and. verbose_fragments)then
        print*,nodesingle,' fragment blocks on single MPI nodes '
        print*,' fraction of operations on single MPI nodes = ',optotsingle, optotall, optotsingle/optotall, assignratio
   end if

!........... FIRST ATTEMPT AT A DISTRIBUTION........
   opavg = (optotall-optotsingle)/real(nnodesx-nodesingle,8)
   nodesum = 0
   do ifrag = 1,nfragments
      do jfrag = 1,nfragments
          if(opfragstat(ifrag,jfrag)%nnodes==1)then
              nodesum = nodesum+1
              opfragstat(ifrag,jfrag)%nnodes=1
              ! KSM:  already 0:  opdisc(nfragments*(ifrag-1)+jfrag) = 0.0d0
              cycle
          endif
          if(opfragstat(ifrag,jfrag)%optot== 0.0d0)then
              opfragstat(ifrag,jfrag)%nnodes=0
              ! KSM:  already 0:  opdisc(nfragments*(ifrag-1)+jfrag) = 0.0d0
              cycle
          end if
          opfrac = opfragstat(ifrag,jfrag)%optot/(optotall-optotsingle)
          opfragstat(ifrag,jfrag)%nnodes = nint(opfrac*real(nnodesx * assignratio -nodesingle,8)) 

          if(opfragstat(ifrag,jfrag)%nnodes==0 .and. opfragstat(ifrag,jfrag)%optot>0)then
              opfragstat(ifrag,jfrag)%nnodes=1
              ! KSM:  already 0:  opdisc(nfragments*(ifrag-1)+jfrag) = 0.0d0
          else
              ! opdisc(nfragments*(ifrag-1)+jfrag) = opfragstat(ifrag,jfrag)%optot/opavg-opfragstat(ifrag,jfrag)%nnodes
             !opdisca(ifrag,jfrag) = opfragstat(ifrag,jfrag)%optot/opavg-opfragstat(ifrag,jfrag)%nnodes
             opdisca(ifrag,jfrag) = (opfragstat(ifrag,jfrag)%optot/opavg- & 
			                        opfragstat(ifrag,jfrag)%nnodes)/opfragstat(ifrag,jfrag)%nnodes
          end if
          nodesum = nodesum + opfragstat(ifrag,jfrag)%nnodes
      end do  
   end do
   
!........ TOO FEW/TOO MANY NODES PASSED OUT........
   nodediff = nnodesx - nodesum
! -- CHECK DISTRIBUTION ---- ADDED in 7.4.9.....
   if(localverbose) then
      call printopfragstat
      if(iproc == 0) print*,' The difference is ',nodediff
   end if
!......... NOW ADD/SUBTRACT NODES FROM FRAGMENT BLOCKS.........................
!
!  because we won't have exactly used up all the nodes, find out how many extra
!  (or missing) nodes there are; and add or subtract from fragment blocks 
!
! ADDED 7.4.9: cannot subtract nodes if only one node assigned

   if (iproc == 0) print *, "NODEDIFF ", nnodesx, nodesum, assignratio

   if(nodediff == 0)return
   if(nodediff < 0)then
		   nodeadd = -1
   else
		   nodeadd = +1
   end if

   do i= 1,abs(nodediff)
	!....... SEARCH FOR LARGEST DIFFERENCE........................
		  diffcf = 0.00
		  ifrag = 0
		  jfrag = 0
		  do si = 1, nfragments
			 do sj = 1, nfragments
				if( opfragstat(si,sj)%nnodes  < 2 .and. nodeadd < 0)cycle  !don't subtract if only 1 nodes assigned
				if(opdisca(si, sj)*nodeadd > diffcf) then
				   diffcf = opdisca(si, sj) * nodeadd
				   ! save position
				   ifrag = si
				   jfrag = sj
				end if
			 end do
		  end do
		  if(jfrag .eq. 0) then
			  print *, "iproc=", iproc, ", Can't find suitable fragment pair to  add/remove nodes"
			  print*,nfragments,' fragments '
			  print*,nodediff,' nodes to be reassigned out of ',nprocs,' with ',nprocs-nfragments**2,' available '
			  stop
		  end if
		  opfragstat(ifrag,jfrag)%nnodes= opfragstat(ifrag,jfrag)%nnodes+nodeadd
		  ! KSM:  opdisc(nfragments*(ifrag-1)+jfrag) = opfragstat(ifrag,jfrag)%optot/opavg-opfragstat(ifrag,jfrag)%nnodes
		  ! opdisca(ifrag, jfrag) = opfragstat(ifrag,jfrag)%optot/opavg-opfragstat(ifrag,jfrag)%nnodes
                  opdisca(ifrag, jfrag) = (opfragstat(ifrag,jfrag)%optot/opavg- & 
				  opfragstat(ifrag,jfrag)%nnodes)/opfragstat(ifrag,jfrag)%nnodes
   end do
   if(iproc==0)write(24,*)'     init frag / final frag  / # nodes assigned'
	   do ifrag = 1,nfragments
              do jfrag = 1,nfragments
        if(iproc==0)write(24,*)ifrag,jfrag,opfragstat(ifrag,jfrag)%nnodes
              end do
	   end do
	   if(.not.assignnodes)return
	!....... NEXT STEP (AT LEAST WHEN MODELING) ASSIGN FRAGMENTS TO NODES.....

	   i = -1
   if(iproc==0)write(24,*)'       node /    init frag  / final frag  / operations load '
	   maxnodeops = 0.0
	   avgnodeops = 0.0
	   do ifrag = 1,nfragments
		  do jfrag = 1,nfragments
			  avgnodeops = avgnodeops + opfragstat(ifrag,jfrag)%optot
			  do j = 1, opfragstat(ifrag,jfrag)%nnodes
				   i = i+1
				   maxnodeops = bmax(maxnodeops, opfragstat(ifrag,jfrag)%optot/real(opfragstat(ifrag,jfrag)%nnodes,8))
				   
               if(iproc==0)write(24,*)i,ifrag,jfrag,opfragstat(ifrag,jfrag)%optot/opfragstat(ifrag,jfrag)%nnodes
			  end do
		  end do ! jfrag
	   end do  ! ifrag
   avgnodeops = avgnodeops/float(nprocs)
!   if(iproc==0)print*,' estimated efficiency = ',avgnodeops/maxnodeops
   close(24)
   return
end subroutine divide_procs_over_fragments
!==========================================================================================
!
!  sets variable "ibelong" in xsd(it)%sector to denote if a particular sector
!  belongs to fragments residing on a MPI process
!
!  added in 7.5.9
!
subroutine setsectorpossession(iprocs)
	use nodeinfo
	use sectors
!	use fragments
	implicit none
	integer :: iprocs   ! dummy label for MPI processes
	integer :: it   ! label of species
	integer(8) :: is,cs,ics  ! label of sector
	integer :: ifrag,ffrag  ! fragment labels
	
	do it =1,2
		if(nfragments==1)then    ! ALL SECTORS BELONG
			do is = 1,nsectors(it)
				xsd(it)%sector(is)%ibelong=.true.
			end do ! is
			cycle
		end if
		
		do is = 1,nsectors(it)        ! DEFAULT
			xsd(it)%sector(is)%ibelong=.false.
        end do		

	end do ! it
!-------------- "POSSESS" proton sectors
	
	ifrag = nodal(iprocs)%ifragment
	ffrag = nodal(iprocs)%ffragment
	do is = 1,fragmentlist(ifrag)%ssectorstart,fragmentlist(ifrag)%ssectorend
		xsd(1)%sector(is)%ibelong=.true.
	end do ! is
	do is = 1,fragmentlist(ffrag)%ssectorstart,fragmentlist(ffrag)%ssectorend
		xsd(1)%sector(is)%ibelong=.true.
	end do ! is

!......... NOW POSSESS NEUTRON FRAGMENTS-----------
    do is = 1,nsectors(1)
		if(xsd(1)%sector(is)%ibelong)then
			do ics = 1,xsd(1)%sector(is)%ncsectors
				cs = xsd(1)%sector(is)%csector(ics)
				xsd(2)%sector(cs)%ibelong=.true.
			end do ! ics
		end if
	end do	
	return
end subroutine setsectorpossession

!======================================================================
!
! added to 7.6.2
!
! determines which conjugate ("neutron") sectors belongs to a particular fragment
! and confirms contiguity
!
!  CALLED BY fragmenter

subroutine fragmentconjugatesectors
	use sectors
!	use fragments
    implicit none
	
    integer ifrag
	logical, allocatable :: sectornocc(:)   ! logical flag to denote occupied neutron sectors
	integer(8) :: is,cs,ics,firstsector,lastsector
	
	allocate(sectornocc( nsectors(2)))

	do ifrag = 1,nfragments
		sectornocc = .false.
		do is = 1,fragmentlist(ifrag)%ssectorstart,fragmentlist(ifrag)%ssectorend
			do ics = 1,xsd(1)%sector(is)%ncsectors
				cs = xsd(1)%sector(is)%csector(ics)
                sectornocc(cs)=.true.
			end do ! ics		
		end do ! is
!----------- FIND LIMITS		
        firstsector = nsectors(2)
		lastsector  = 1
		do is = 1,nsectors(2)
			if(sectornocc(is))then
				firstsector = min(firstsector,is)
				lastsector  = max(lastsector, is)
			end if
		end do
		
!--------- CHECK FOR CONTIGUITY--------------------
	 	
        do is = firstsector,lastsector
	        if(.not.sectornocc(is))then
				print*,' OOPS! Non continuous neutron sectors in fragment ',ifrag
				print*,' first, last sectors = ',firstsector,lastsector
				print*,' but neutrons sector ',is,' not included '
				stop
			end if
        end do
!----------- SET LIMITS ---------------
        fragmentlist(ifrag)%csectorstart = firstsector		
        fragmentlist(ifrag)%csectorend =   lastsector		
		
	end do  ! ifrag	

    deallocate(sectornocc)
    return
end	subroutine fragmentconjugatesectors


!===============================================
! debug print of tribution data structure.
! should print on all processes
subroutine printtribution
!   use fragments
   use nodeinfo
   use io
   integer :: i

   if(iproc /= 0) return
   write(logfile,*) "nprocs = ", nprocs
   write(logfile,*)"nodeal(proc) = {ifragment, ffragment }"
   do i = 0, nprocs-1
      write(logfile,*) "nodeal(", i, ") =", "{", nodal(i)%ifragment, ", ", nodal(i)%ffragment, "}"
   end do
   return
end subroutine printtribution

!===============================================================
!   setnodal!
!  assigns MPI processes / compute nodes to initial and final fragments
!  sets up "nodal" and fills data
!
!  IMPORTANT: opfragstat(ifrag,ffrag)%nnodes must be previous set
!            done in distro_opbundles_over_fragments and 
!             and in divide_procs_over_fragments
!
!  CALLED BY:
!    master_para_distribute
!    density1b_from_oldwfn
!    overlap
!
subroutine setnodaltribution
!   use fragments
   use nodeinfo
   use io
   use reporter
   implicit none
   integer ifrag,jfrag
   integer i,j,ii
   integer :: aerr
      
   allocate(nodal(0:nprocs-1), stat=aerr )
   if(aerr /= 0) call memerror("setnodaltribution")
   if(iproc==0 .and. nfragments > 1) then
      open(unit=fragfile,file='fraginfo.bigstick',status='unknown')
      write(fragfile,*)' # of fragments ',nfragments
      write(fragfile,*)' fragment label     /       dim '
      write(logfile,*)' # of fragments ',nfragments
      write(logfile,*)' fragment label     /       dim '
	  do ii = 1,nfragments
		  write(fragfile,*)ii,fragmentlist(ii)%localdim
		  write(logfile,*)ii,fragmentlist(ii)%localdim
		  
	  end do
      do ii = 1, nfragments
         write(fragfile, "(I3,A,I12,A,I12,A)") ii, ":(", fragmentlist(ii)%basisstart, ",", fragmentlist(ii)%basisend, ")"
         write(fragfile, "(A,I0,A,I0)") "   sectorstart=", fragmentlist(ii)%ssectorstart, ", sectorend=", & 
		                                                          fragmentlist(ii)%ssectorend
         write(fragfile, "(A,I0)")      "   nsectors=", fragmentlist(ii)%nsectors
      end do
   end if

   i = -1
   if(nprocs==1)then  ! KSM - change from nproc.  Was clearing MPI data during modeling
      opfragstat(:,:)%nnodes = 1
   end if
   do ifrag = 1,nfragments
      do jfrag = 1,nfragments
          if (iproc == 0) then
            write(logfile,*) "opfragstat(", ifrag, ", ", jfrag, ")%nnodes = ", opfragstat(ifrag,jfrag)%nnodes
          end if
		  if(iproc==0 .and.opfragstat(ifrag,jfrag)%nnodes >0)then
			  write(fragfile,*)' frags ',ifrag,'->',jfrag,' on procs ',i+1,i+opfragstat(ifrag,jfrag)%nnodes
			  write(logfile,*)' frags ',ifrag,'->',jfrag,' on procs ',i+1,i+opfragstat(ifrag,jfrag)%nnodes
			  
		  end if
          do j = 1, opfragstat(ifrag,jfrag)%nnodes
            i = i+1
            if(i >= nprocs) then    ! nprocs works even for modeling
               if(iproc == 0) then
                  print *, "setnodaltribution: nprocs=", nprocs, ", is too small"
                  write(logfile,*) "setnodaltribution: nprocs=", nprocs, ", is too small"
                  call printstacktrace
               end if
               stop 1
            end if
            nodal(i)%ifirst = j .eq. 1 .and. jfrag .eq. 1 ! marks one node which reads ifrag
            nodal(i)%ifragment=ifrag
            nodal(i)%ffragment=jfrag
!            nodal(i)%sfragment=0  ! KSM: init for debug checking
!            nodal(i)%frac_ops=0.0  ! KSM: init for debug checking
            ! don't overwrite ifirst here, look up 5 lines!
          end do
       end do
   end do

   call checknodaltribution
   if(nfragments > 1)close(fragfile)
   return
end subroutine setnodaltribution

!===============================================================
!
!  checks assignment of MPI processes/nodes to fragments against assignment of bundles
!  added 7.6.8
!
subroutine checknodaltribution
    use nodeinfo
    use io
	use opbundles
	use bmpi_mod
    implicit none
    integer ifrag,ffrag
	integer :: jprocs
    integer ibundle
	logical :: muststop
	
	integer :: ierr 
	
    muststop = .false.
	do jprocs = 0,nprocs-1
		ifrag = nodal(jprocs)%ifragment
		ffrag = nodal(jprocs)%ffragment
        if(opbundlestart(jprocs) < 1 .or. opbundleend(jprocs) < 1)cycle
		do ibundle = opbundlestart(jprocs),opbundleend(jprocs)
			if(opbundle(ibundle)%ifragment/=ifrag .or. opbundle(ibundle)%ffragment/=ffrag)then
				muststop = .true.
				if(iproc==0)then
					print*,' Mismatch of fragments on MPI process/node ',jprocs
					print*,' Variable NODAL expects initial, final ',ifrag,ffrag
					print*,' but assigned bundle ',ibundle,' expects ', opbundle(ibundle)%ifragment,opbundle(ibundle)%ffragment
					
				end if
			end if
			
			
		end do
		
	end do
	if(muststop)then
		call BMPI_ABORT(icomm,101,ierr)
		stop
	end if
	return
	
end subroutine checknodaltribution
!===============================================================
!
! alternate routine for assigning MPI processes
! based upon OPBUNDLESTART/END 
!
subroutine setnodaltribution_alt
	
!    use fragments
    use nodeinfo
    use io
	use opbundles
    implicit none
    integer ifrag,ffrag
	integer :: bstart,bend
	integer :: jprocs,kprocs
    integer ibundle
	integer :: aerr
	
    allocate(nodal(0:nprocs-1), stat=aerr )
    if(aerr /= 0) call memerror("setnodaltribution_alt")

    if(nprocs==1)then  ! KSM - change from nproc.  Was clearing MPI data during modeling
       opfragstat(:,:)%nnodes = 1
    end if
	
    if(iproc==0)print*,' Alternate nodal distribution routine '
	do jprocs = 0,nprocs-1
		bstart = opbundlestart(jprocs)
		bend   = opbundleend(jprocs)
		
		if(bstart==0 .or. bend ==0)then
			nodal(jprocs)%ifragment = nfragments  ! assign to last fragment
			nodal(jprocs)%ffragment = nfragments  ! assign to last fragment
			cycle
		end if
!............ CHECK ALL ASSIGNED OPBUNDLES HAVE SAME FRAGMENTS .............		
		
		ifrag = opbundle(bstart)%ifragment
		ffrag = opbundle(bend)%ffragment
		do ibundle = bstart,bend
			if(opbundle(ibundle)%ifragment /= ifrag .or. opbundle(ibundle)%ffragment/=ffrag)then
				if(iproc==0)then
					print*,' Bad choice of fragments on opbundles for process ',jprocs
					print*,' opbundle start/end = ',bstart,bend
					print*,' Expect initial/final fragments = ',ifrag,ffrag
					print*,' but bundle ',ibundle,' has ',opbundle(ibundle)%ifragment,opbundle(ibundle)%ffragment
				end if
				stop
			end if
			
		end do
		nodal(jprocs)%ifragment =ifrag
		nodal(jprocs)%ffragment =ffrag
		
!......... FIND 'FIRST' NODE WITH THESE FRAGMENTS ........
        nodal(jprocs)%ifirst = .true.
		
		do 	kprocs = 0,jprocs-1
			if(nodal(kprocs)%ifragment ==ifrag .and. nodal(kprocs)%ffragment==ffrag)then
		        nodal(jprocs)%ifirst = .false.
				exit
			end if
		end do	

		
	end do
	return
		
end subroutine setnodaltribution_alt
!======================================================================
!
! Debugging routine to dump the vec1/vec2 bounds for a given rank
!
!  CALLED BY: 
!     check_opbundle_bounds
!
subroutine print_node_bounds(fn, jnode)
   use precisions
   use nodeinfo
!   use fragments
   use basis
   use butil_mod
   implicit none
   integer, intent(in) :: fn  ! file number
   integer, intent(in) :: jnode
   integer :: frag1, frag2 ! fragments for vec1 and vec2
   integer(kind=basis_prec) :: f1s, f1e, f2s, f2e  ! fragment bounds

   frag1 = nodal(jnode)%ifragment
   frag2 = nodal(jnode)%ffragment
   ! get bounds for vec1 and vec2
   f1s = basestart(frag1)
   f1e = basestop(frag1)
   f2s = basestart(frag2)
   f2e = basestop(frag2)
   write(fn, "(A,I0,A,I0,A,I0,A,I0,A,I0,A)") "node=", jnode, & 
      ", vec1 bounds=(", f1s, ":", f1e, "), vec2 bounds=(", f2s, ":", f2e, ")"
   flush(fn)
end subroutine print_node_bounds
!===============================================================

end module fragments
