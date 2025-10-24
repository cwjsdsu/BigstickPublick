!========================================================
!
! file BPARALLEL_MAIN.f90
!
! information for modeling BIGSTICK data
! and distributing over MPI (distributed memory) nodes
!
! started 8/2011 by CWJ @ SDSU
!
!=======================================================
!
!  SOME GENERAL NOTES ON MPI PARALLELIZATION  (July 2014)
!
!  BIGSTICK has both MPI (distributed memory) and OpenMP (shared memory)
!  parallelization.  Most of the work to calculate and carry out the
!  distribution of work across MPI compute nodes are contained in 
!  files bparallel_libX.f90. 
!
!  One wants to distribute both work (load-balancing) and memory (memory caps)
!
!  3 kinds of information has to be distributed across MPI nodes:
!  -- the matrix-vector multiply operations  (load-balancing)
!  -- memory for jump arrays encoding the mat-vec operations (memory cap)
!    ( ++ also in here one needs to worry about the uncoupled matrix elements)
!  -- memory for the Lanczos vectors themselves (memory cap)
!
!  If the dimension of the many-body space is very large, then one must
!  break up the Lanczos vectors into "fragments." Information about fragments
!  are 
!
!  The first are the "jump" arrays which organize the mat-vec multiplications.
!  The jumps arrays are contained in module jumpNbody in bmodules_main.f90
!  They are generated in bjumplibX.f90 files.  
!
!  But the jumps are combined into mat-vec operations.  
!  (You should think of jumps as forming the x- and y-edges of a rectangle;
!   the area of the rectangle is the number of operations, while the perimeter
!   is the number of jumps needed to combine into those operations.)
!  The structures organizing and combining jumps are "opbundles," and 
!  become the basic structure for controlling mat-vec operations
!  and for load-balancing. 
!
!  DATA FOR LANCZOS VECTOR FRAGMENTS
!    information about fragments in fragmentlist found in module fragments 
!    in bmodule_parallel.f90
!
!    actual allocation for current Lanczos fragments (even if only one "fragment")
!    found in module localvectors in bmodules_parallel.f90
!
!  DATA FOR JUMPS, OPBUNDLES
!    jump arrays found in module jumpNbody in bmodules_main.f90
!    opbundles found in module opbundles in bmodules_parallel.f90
!
!  DATA FOR OPERATIONS
!    the structure opstat contains information on operations between sectors
!    the structure opfragstat contains information on operations between fragments
!    both are found in module operation_stats in bmodules_parallel.f90
!  
!  OVERVIEW OF WORK relevant to parallelization (note: does not cover all of BIGSTICK)
!    
!  !!!!   NOTE SOME OF THE FOLLOWING MAY HAVE CHANGED!!!!
!    
!          + +  preparatory work: + +
!
!     subroutine hoopmaster counts up the number of jumps 
!         CALLED in main program in bigstick_main.f90
!         FOUND in bjumplib1.f90
!     subroutine master_op_stat counts up the operations between sectors
!         CALLED in main program in bigstick_main.f90
!         FOUND in bparallel_lib1.f90
!
!         -->> CALCULATION DISTRIBUTION OF OPERATION, MEMORY  <<--
!
!     subroutine master_para_distribute calculates the distribution across MPI nodes
!     This is the main subroutine for parallelization.
!         CALLED in main program in bigstick_main.f90
!         FOUND in bparallel_lib1.f90
!
!         calls vectorstorageceiling (below) : checks to see if one should/can create fragments
!         if break_vectors then
!           calls defaultsectorbreakup (below)
!              (targeted for obsolescense; equates "subsectors" as "sectors"
!           calls fragmenter           ! ==> in bparallel_lib2.f90 ; defines fragments
!           calls getfragmentstats     !     in bparallel_lib2.f90 
!                  counts operations between fragments 
!                  (collates information from opstat into opfragstat)
!           calls master_distro_over_fragments   in bparallel_lib2.f90  NOT FINISHED
!                 which
!                 calls divide_nodes_over_fragments   in bparallel_lib2.f90
!                 calls distro_opbundles_1frag_g      in bparallel_lib3.f9,  
!                        (will get replaced by a general routine)
!                 calls howmuchstorage                in bparallel_lib2.f90
!
!        else
!          call defaultsectorbreakup  (below)
!          call defaultfragments       (below) : defines a single fragment
!        endif
!        call fragmentvectorsize      (below)
!        call subsectoropstats        (below; should go away)
!        call checksubopsum          (below; should go away)
!
!        calls proc_clock  starts instrument on timing all the MPI processes
!        if distributeMPI then
!            calls distro_opbundles_1frag_g   a key routine, found in bparallel_lib3.f90
!                which
!               calls count_create_draft_opbundles_g    creates draft versions of opbundles
!               calls analyze_draft_opbundle            analyzes how much work each opbundles does
!               calls split_draft_opbundles_g           splits the opbundles 
!               calls setMPIjumplimits                 
!               calls count_jumps_on_node               checks to see if too many jumps stored on
!                                                       any MPI compute node
!               calls redistributor                     redistribute jumps if needed
!
!        else
!            calls defaultopbundles
!       endif
!       calls setMPIjumplimits     
!
!          + + +  final work:  + + +
!
!     (read in and uncouple the two- and three-body matrix elements)
!
!    subroutine jumpmaster generates the jumps
!    Only creates the required jumps on a given MPI compute node
!         CALLED in main program in bigstick_main.f90
!         FOUND in bjumplib1.f90
!
!    (after this point the code carries out Lanczos)
!

!=====================================================================
!
!
!  MAIN ROUTINES FOR PARALLEL DISTRIBUTION
!
module para_main_mod
use para_bundles_mod
use para_util_mod	
use para_report_mod

implicit none


contains
	
!=====================================================================

!  
!  master routine for controlling distribution of work across MPI compute processes
!  for 'normal' runs
!
! CALLED BY
!     main routine in bigstick_main.f90
!
! SUBROUTINES CALLED:
!   print_memory_limits  (in this file)
!   master_op_stats     counts operations between sectors, fragments
!   proc_clock   (in file btimelib.f90)
!   distro_opbundles_over_fragments  (in this file)
!   setnodeltribution  (in file bfragments.f90)
!   setMPIjumpslimits  (in bparallel_util.f90)
!          THESE ARE FOR NEW, ADVANCED STORAGE
!   survey_qs_on_MPIprocs  ! added 7.9.3; found in btbme_adv.f90
!   uncouplePNtbme_boss_new   in btbmepn.f90
!             ----------------
!   print_all_distro  (in bparallel_report.f90)
!
subroutine master_para_distribute
   use nodeinfo
   use flagger
   use basis
   use menu_choices
   use jumplimits
   use fragments
   use operation_stats
   use qs_on_proc
   implicit none
   logical filljumpdirect   ! flag whether or not to set jump redirects
   
 !................................................
 
   call print_memory_limits  ! added in 7.9.0
   call master_op_stat    ! calculates # of operations in bundles; SLATED FOR OBSOLESCENCE
!.......... PRINT OUT TOTAL JUMPS HERE........................   
   call proc_clock(nprocs,'set')
   if(nprocs > 1 .and. compactjumpstorage_enabled)then
       compactjumpstorage = .true.
   else
       compactjumpstorage = .false.
   end if

   call distro_opbundles_over_fragments  !  main routine for distribution of workload
   
   call setnodaltribution
   if(menu_char=='m')then
	   filljumpdirect=.false.
   else
	   filljumpdirect=.true.
   end if
   call setMPIjumplimits(filljumpdirect)   ! final call for this routine
   
!------- IN PREPARATION FOR ADVANCED STORAGE---------   
!   call survey_qs_on_MPIprocs
!   call uncouplePNtbme_boss_new(.true.)
!------------------------------------------
   call print_all_distro
   return
end subroutine master_para_distribute

!=================================================================
!
!  master routine to distribute "opbundles" (and thus mat-vec operations) over
!  fragmented (distributed) Lanczos vectors.  
!
!  CALLED BY:
!     master_para_distribute
!
!  SUBROUTINES CALLED:
!     draft_opbundles  : routine to generate first draft of opbundles
!     divide_procs_over_fragments 
!     split_draft_opbundles_frag  
!     rebalance_MPI_distro
!

subroutine distro_opbundles_over_fragments
   use opbundles
   use fragments
   use verbosity
   use nodeinfo
   use operation_stats
   use flagger
   use menu_choices
   use bmpi_mod
   use jumplimits
   use para_bundles_mod
   implicit none

   integer ifrag,ffrag
   logical successlocal,successall
   integer :: nprocs_avail, nprocs_reserved_tot
   integer, allocatable :: nprocs_reserved(:,:)  
   integer :: distro_loop,max_distro_loops
   integer :: nob_draft,nob,nob_reserved
   integer :: Nmyprocs, myprocstart
   integer :: iprocs
   integer :: aerr

   integer ierr   
   logical fragments_assigned   !  a flag to indicate if fragments have already been assigned
   
   integer:: numob,numobtot
   integer(8) :: numops,numopssum,numopstot
   
   logical :: verbose_number_opbundles =.false.
   
   fragments_assigned = allocated(nodal)

!......... CREATE "DRAFT" OPBUNDLES......

   call draft_opbundles(nob_draft)   ! set up initial op_bundles

   nob_reserved = int(nob_draft*2)      ! allocated twice as many additional opbundles

!..................... ALLOCATE ARRAYS..........  
   allocate(opbundle(nob_reserved), stat=aerr )
   if(aerr /= 0) then
      call memerror("distro_opbundles_over_fragments 1")
      stop 5
   end if

   if(.not.allocated(opbundlestart)) then
      allocate( opbundlestart(0:nprocs-1), stat=aerr)
      if(aerr /= 0) then
         call memerror("distro_opbundles_over_fragments 2")
         stop 5
      end if
   end if
   if(.not.allocated(opbundleend)) then
      allocate( opbundleend(0:nprocs-1), stat=aerr)
      if(aerr /= 0) then
         call memerror("distro_opbundles_over_fragments 3")
         stop 5
      end if
   end if

!....... DISTRIBUTE OPBUNDLES OVER NODES.....
!  --   draft assignment of fragments over nodes
!  --   split the draft opbundles
!  --   check for storage of jumps
!  --   --  if surpass ceiling for storage of jumps, reserve nodes and go back
!
!  We want to distribute work over nprocs 
! (nprocs either = nproc, the physical # of MPI processes, 
!  OR, when "modeling" (menu option "m"), the # of MPI processes to be modeled )
!
!  The routine divide_procs_over_fragments distributes MPI processes
!  for work between fragments. 
!  The routine split_and_rebalance_distro will, for a given initial and final fragment,
!  split up that work, by splitting (if necessary) and assigning opbundles.
!  There is an additional distribution requirement: the number of jumps stored on a 
!  given MPI process may be capped.  If this cap is exceeded, then additional MPI procs
!  must be assigned and the jumps distributed over them.
!  Therefore each block of work between initial and final fragments, ifrag/ffrag,
!  has assigned a certain number of MPIprocess, a start of the processes, and 
!  some number reserved for overbooking of jumps
!
!  IF, however, fragments have already been assigned to various MPI processes,(fragments_assigned=.true.)
!  we can only work within the assigned MPI processes
!

   if(.not. allocated(nprocs_reserved) )allocate (nprocs_reserved(nfragments, nfragments), stat=aerr )
   if(aerr /= 0) then
      call memerror("distro_opbundles_over_fragments 10")
      stop 5
   end if
   
   max_distro_loops = 10
   nprocs_avail = nprocs 
   nprocs_reserved_tot = 0
   nprocs_reserved(:,:) = 0
   
   do distro_loop = 1,max_distro_loops
      myprocstart = 0
      if(iproc==0 .and. verbose_distro)print*,'Attempt ',distro_loop,' at distribution ',nprocs_reserved_tot
      successall = .true.
! .... ASSIGN MPI PROCESSES TO FRAGMENTS........................	  
      nprocs_avail = nprocs - nprocs_reserved_tot
	  if(iproc==0 .and. verbose_distro)print*,' Assigning using only ',nprocs_avail  
      call divide_procs_over_fragments(nprocs_avail,.true.)  ! NOTE: If we are computing J and T afterwards, we do not want to
	                                                         ! change the division of fragments
	   if(iproc==0 .and. verbose_distro)print*,' Nodes assigned between fragments '
															 
      nob = 0
!........................ DO INITIAL DISTRIBUTION.....................
      do ifrag = 1,nfragments
         do ffrag = 1,nfragments             
             Nmyprocs = opfragstat(ifrag,ffrag)%nnodes  ! value computed in divide_procs_over_fragments
             Nmyprocs = Nmyprocs + nprocs_reserved(ifrag,ffrag)
             call split_draft_opbundles_frag(.false.,nob_draft,ifrag,ffrag,  & 
                        Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag),nob)
             myprocstart = myprocstart+Nmyprocs
         end do  ! ffrag
      end do ! ifrag
#ifdef _MPI	
	 if( verbose_distro)call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
!................... SEE IF NEED TO RESERVE MORE OPBUNDLES..............j	  
      if(nob > nob_reserved)then  ! need to reserve more opbundles
        if(.not. allocated(opbundle)) then
            print *, "opbundle deallocated without prior allocation"
         end if
         deallocate(opbundle)
         nob_reserved = int(2*nob)
         allocate(opbundle(nob_reserved), stat=aerr )
         if(aerr /= 0) then
            call memerror("distro_opbundles_over_fragments 20")
            stop 5
         end if
      end if

!...................... NOW CHECK FOR OVERBOOKING OF STORAGE FOR JUMPS................
      nob = 0
      myprocstart = 0
      nprocs_reserved_tot = 0
#ifdef _MPI	
	 if( verbose_distro)call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif	  
      do ifrag = 1,nfragments
         do ffrag = 1,nfragments
!			 if(iproc==0 .and. verbose_distro)print*,' distro on fragments ',ifrag,ffrag
             
             Nmyprocs = opfragstat(ifrag,ffrag)%nnodes  ! value computed in divide_procs_over_fragments
             Nmyprocs = Nmyprocs + nprocs_reserved(ifrag,ffrag)
             opfragstat(ifrag,ffrag)%nnodes= opfragstat(ifrag,ffrag)%nnodes+ nprocs_reserved(ifrag,ffrag)
             call split_draft_opbundles_frag(.true.,nob_draft,ifrag,ffrag,  & 
                        Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag),nob)
			 if(iproc==0 .and. verbose_distro)print*,ifrag,ffrag,nob,' opbundles so far '
             if(setjumpceiling)then
!				     if(iproc==0 .and. verbose_distro)print*,' attempting to rebalance...'
                call rebalance_MPI_distro(ifrag,ffrag, successlocal,Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag))
             end if
			 
             if(.not.successlocal)then
				 successall =.false.
				 if(iproc==0)print*,' failure to rebalance ',ifrag,ffrag,', need to reserve ',& 
				 nprocs_reserved(ifrag,ffrag)
			 else
!				 if(iproc==0 .and. verbose_distro)print*,' rebalanced ',ifrag,ffrag				 
			 end if
             nprocs_reserved_tot = nprocs_reserved_tot + nprocs_reserved(ifrag,ffrag)
             myprocstart = myprocstart+Nmyprocs
#ifdef _MPI	
			 if( verbose_distro)call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif			 
			 if( iproc==0 .and. verbose_distro)print*,' finished distro on fragments ',ifrag,ffrag
			 
         end do  ! ffrag
      end do ! ifrag
	 if(iproc==0 .and. verbose_distro)print*,' Finished looping over frag2frag ',successall
	  
      if(successall)then   ! one final round with fill the opbundles
		 if(nob < nob_draft)then  ! some error has occurred
			 if(iproc==0)print*,' Error in splitting opbundles , ended with only ',nob,' out of ',nob_draft
			 stop
		 end if 
		  
         nopbundles = nob
         if(iproc==0)then
            print*,' '
            print*,nopbundles,' final opbundles '
            print*,' '
         end if
!..............   ALLOCATE FINAL OPBUNDLES................
        if(nob > nob_reserved)then  ! need to reserve more opbundles
           if(.not. allocated(opbundle)) then
              print *, "opbundle deallocated without prior allocation"
           end if
           deallocate(opbundle)
           nob_reserved = int(1.2*nob)
           allocate(opbundle(nob_reserved), stat=aerr )
           if(aerr /= 0) call memerror("distro_opbundles_over_fragments 21")
         end if
         myprocstart = 0
         nprocs_avail = nprocs - nprocs_reserved_tot
         call divide_procs_over_fragments(nprocs_avail,.true.)
         nob = 0
         do ifrag = 1,nfragments
            do ffrag = 1,nfragments
               opfragstat(ifrag,ffrag)%nnodes= opfragstat(ifrag,ffrag)%nnodes+ nprocs_reserved(ifrag,ffrag)
            
               Nmyprocs = opfragstat(ifrag,ffrag)%nnodes  ! value computed in divide_procs_over_fragments
               if(Nmyprocs==0)cycle
               call split_draft_opbundles_frag(.true.,nob_draft,ifrag,ffrag,Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag),nob)
               if(setjumpceiling)then
                  call rebalance_MPI_distro(ifrag,ffrag, successlocal,Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag))
               end if
               myprocstart = myprocstart+Nmyprocs
            end do  ! ffrag
         end do ! ifrag
!................ CHECK TO SEE IF WE LOST ANY...........
         call check_ops('SPE',nob_draft)
         call check_ops('PP ',nob_draft)
         call check_ops('NN ',nob_draft)
         call check_ops('PN ',nob_draft)
         call check_ops('PPP',nob_draft)
         call check_ops('PPP',nob_draft)
         call check_ops('PPN',nob_draft)
         call check_ops('PNN',nob_draft)
         call check_ops('NNN',nob_draft)

         exit

      else
		 if(iproc==0 .and. verbose_distro)print*,' Failed, need to try again'
		  

      end if
    end do
	if(iproc==0 .and. verbose_distro_report)print*,' # distro loops = ',distro_loop
    if(.not.successall)then
       if(iproc==0)then   ! TRY TO END GRACEFULLY
          print*,' FAILED to distribute '
       end if

       if(menu_char/= 'm')then
		   if(iproc==0)print*,' Ending run '
#ifdef _MPI	
		  call BMPI_Barrier(MPI_COMM_WORLD,ierr)		   
          call BMPI_Abort(MPI_COMM_WORLD, 106, ierr)
#endif
       end if
    end if
!........ PRINTOUT INTRONS.......

    if(compactjumpstorage) call intron_master(.false.,.true.)   ! PRINT OUT INTRONS
	
!........ CHECK THE NUMBER OF OPBUNDLES OF EACH TYPE...ADDED 7.10.2 to diagnose TBDM
if (verbose_number_opbundles )then
   numobtot= 0
   numopssum = 0
   numopstot = 0
   call count_opbundle_types('PP ','f',numob,numops)	
   print*,' Number of PP f ',numob
   numobtot=numobtot+numob
   numopstot = numopstot+numops
   numopssum = numopssum+numops
   
   call count_opbundle_types('PP ','b',numob,numops)	
   print*,' Number of PP b ',numob
   numobtot=numobtot+numob
   numopstot = numopstot+numops
   numopssum = numopssum+numops
   
   call count_opbundle_types('PP ','h',numob,numops)	
   print*,' Number of PP h ',numob
   numobtot=numobtot+numob
   numopstot = numopstot+numops
   numopssum = numopssum+numops
   print*,' PP tot ',numopssum
   numopssum =0
   
   call count_opbundle_types('NN ','f',numob,numops)	
   print*,' Number of NN f ',numob
   numobtot=numobtot+numob
   numopstot = numopstot+numops
   numopssum = numopssum+numops
   
   call count_opbundle_types('NN ','b',numob,numops)	
   print*,' Number of NN b ',numob
   numobtot=numobtot+numob
   numopstot = numopstot+numops
   numopssum = numopssum+numops
   
   call count_opbundle_types('NN ','h',numob,numops)	
   print*,' Number of NN h ',numob
   numobtot=numobtot+numob
   numopstot = numopstot+numops
   numopssum = numopssum+numops
   print*,' NN tot ',numopssum
   numopssum =0
   call count_opbundle_types('PN ','f',numob,numops)	
   print*,' Number of PN f ',numob
   numobtot=numobtot+numob
   numopstot = numopstot+numops
   numopssum = numopssum+numops
   
   call count_opbundle_types('PN ','b',numob,numops)	
   print*,' Number of PN b ',numob
   numobtot=numobtot+numob
   numopstot = numopstot+numops
   numopssum = numopssum+numops
   
   call count_opbundle_types('PN ','h',numob,numops)	
   print*,' Number of PN h ',numob
   numobtot=numobtot+numob
   numopstot = numopstot+numops
   numopssum = numopssum+numops
   print*,' PN tot ',numopssum
   numopssum =0
   print*,' total = ',numobtot,nopbundles,numopstot
end if
!...... ALL FINISHED; PUBLISH THE DISTRIBUTION....

   return

end subroutine distro_opbundles_over_fragments

!=================================================================
!
!  master routine to distribute "opbundles" (and thus mat-vec operations) over
!  fragmented (distributed) Lanczos vectors.  
!
!  NEW routine in 7.9.6 to account for 'greedy' opbundles
!
!  CALLED BY:
!     master_para_distribute
!
!  SUBROUTINES CALLED:
!     draft_opbundles  : routine to generate first draft of opbundles
!     divide_procs_over_fragments 
!     split_draft_opbundles_frag  
!     rebalance_MPI_distro
!

subroutine distro_opbundles_over_fragments_new
   use opbundles
   use fragments
   use verbosity
   use nodeinfo
   use operation_stats
   use flagger
   use menu_choices
   use bmpi_mod
   use jumplimits
   use para_bundles_mod
   implicit none

   integer ifrag,ffrag
   logical successlocal,successall
   integer :: nprocs_avail, nprocs_reserved_tot
   integer, allocatable :: nprocs_reserved(:,:)  
   integer :: distro_loop,max_distro_loops
   integer :: nob_draft,nob,nob_reserved
   integer :: Nmyprocs, myprocstart
   integer :: iprocs
   integer :: aerr

   integer ierr   
   logical fragments_assigned   !  a flag to indicate if fragments have already been assigned
   
   fragments_assigned = allocated(nodal)

!......... CREATE "DRAFT" OPBUNDLES......

   call draft_opbundles(nob_draft)   ! set up initial op_bundles

   nob_reserved = int(nob_draft*2)      ! allocated twice as many additional opbundles

!..................... ALLOCATE ARRAYS..........  
   allocate(opbundle(nob_reserved), stat=aerr )
   if(aerr /= 0) then
      call memerror("distro_opbundles_over_fragments 1")
      stop 5
   end if

   if(.not.allocated(opbundlestart)) then
      allocate( opbundlestart(0:nprocs-1), stat=aerr)
      if(aerr /= 0) then
         call memerror("distro_opbundles_over_fragments 2")
         stop 5
      end if
   end if
   if(.not.allocated(opbundleend)) then
      allocate( opbundleend(0:nprocs-1), stat=aerr)
      if(aerr /= 0) then
         call memerror("distro_opbundles_over_fragments 3")
         stop 5
      end if
   end if

!....... DISTRIBUTE OPBUNDLES OVER NODES.....
!  --   draft assignment of fragments over nodes
!  --   split the draft opbundles
!  --   check for storage of jumps
!  --   --  if surpass ceiling for storage of jumps, reserve nodes and go back
!
!  We want to distribute work over nprocs 
! (nprocs either = nproc, the physical # of MPI processes, 
!  OR, when "modeling" (menu option "m"), the # of MPI processes to be modeled )
!
!  The routine divide_procs_over_fragments distributes MPI processes
!  for work between fragments. 
!  The routine split_and_rebalance_distro will, for a given initial and final fragment,
!  split up that work, by splitting (if necessary) and assigning opbundles.
!  There is an additional distribution requirement: the number of jumps stored on a 
!  given MPI process may be capped.  If this cap is exceeded, then additional MPI procs
!  must be assigned and the jumps distributed over them.
!  Therefore each block of work between initial and final fragments, ifrag/ffrag,
!  has assigned a certain number of MPIprocess, a start of the processes, and 
!  some number reserved for overbooking of jumps
!
!  IF, however, fragments have already been assigned to various MPI processes,(fragments_assigned=.true.)
!  we can only work within the assigned MPI processes
!

   if(.not. allocated(nprocs_reserved) )allocate (nprocs_reserved(nfragments, nfragments), stat=aerr )
   if(aerr /= 0) then
      call memerror("distro_opbundles_over_fragments 10")
      stop 5
   end if
   
   max_distro_loops = 10
   nprocs_avail = nprocs 
   nprocs_reserved_tot = 0
   nprocs_reserved(:,:) = 0
   
   do distro_loop = 1,max_distro_loops
      myprocstart = 0
      if(iproc==0 .and. verbose_distro)print*,'Attempt ',distro_loop,' at distribution ',nprocs_reserved_tot
      successall = .true.
! .... ASSIGN MPI PROCESSES TO FRAGMENTS........................	  
      nprocs_avail = nprocs - nprocs_reserved_tot
	  if(iproc==0 .and. verbose_distro)print*,' Assigning using only ',nprocs_avail  
      call divide_procs_over_fragments(nprocs_avail,.true.)  ! NOTE: If we are computing J and T afterwards, we do not want to
	                                                         ! change the division of fragments
	   if(iproc==0 .and. verbose_distro)print*,' Nodes assigned between fragments '
															 
      nob = 0
!........................ DO INITIAL DISTRIBUTION.....................
      do ifrag = 1,nfragments
         do ffrag = 1,nfragments             
             Nmyprocs = opfragstat(ifrag,ffrag)%nnodes  ! value computed in divide_procs_over_fragments
             Nmyprocs = Nmyprocs + nprocs_reserved(ifrag,ffrag)
             call split_draft_opbundles_frag(.false.,nob_draft,ifrag,ffrag,  & 
                        Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag),nob)
             myprocstart = myprocstart+Nmyprocs
         end do  ! ffrag
      end do ! ifrag
#ifdef _MPI	
	 if( verbose_distro)call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
!................... SEE IF NEED TO RESERVE MORE OPBUNDLES..............j	  
      if(nob > nob_reserved)then  ! need to reserve more opbundles
        if(.not. allocated(opbundle)) then
            print *, "opbundle deallocated without prior allocation"
         end if
         deallocate(opbundle)
         nob_reserved = int(2*nob)
         allocate(opbundle(nob_reserved), stat=aerr )
         if(aerr /= 0) then
            call memerror("distro_opbundles_over_fragments 20")
            stop 5
         end if
      end if

!...................... NOW CHECK FOR OVERBOOKING OF STORAGE FOR JUMPS................
      nob = 0
      myprocstart = 0
      nprocs_reserved_tot = 0
#ifdef _MPI	
	 if( verbose_distro)call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif	  
      do ifrag = 1,nfragments
         do ffrag = 1,nfragments
!			 if(iproc==0 .and. verbose_distro)print*,' distro on fragments ',ifrag,ffrag
             
             Nmyprocs = opfragstat(ifrag,ffrag)%nnodes  ! value computed in divide_procs_over_fragments
             Nmyprocs = Nmyprocs + nprocs_reserved(ifrag,ffrag)
             opfragstat(ifrag,ffrag)%nnodes= opfragstat(ifrag,ffrag)%nnodes+ nprocs_reserved(ifrag,ffrag)
             call split_draft_opbundles_frag(.true.,nob_draft,ifrag,ffrag,  & 
                        Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag),nob)
			 if(iproc==0 .and. verbose_distro)print*,ifrag,ffrag,nob,' opbundles so far '
             if(setjumpceiling)then
!				     if(iproc==0 .and. verbose_distro)print*,' attempting to rebalance...'
                call rebalance_MPI_distro(ifrag,ffrag, successlocal,Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag))
             end if
			 
             if(.not.successlocal)then
				 successall =.false.
				 if(iproc==0)print*,' failure to rebalance ',ifrag,ffrag,', need to reserve ',& 
				 nprocs_reserved(ifrag,ffrag)
			 else
!				 if(iproc==0 .and. verbose_distro)print*,' rebalanced ',ifrag,ffrag				 
			 end if
             nprocs_reserved_tot = nprocs_reserved_tot + nprocs_reserved(ifrag,ffrag)
             myprocstart = myprocstart+Nmyprocs
#ifdef _MPI	
			 if( verbose_distro)call BMPI_Barrier(MPI_COMM_WORLD,ierr)
#endif			 
			 if( iproc==0 .and. verbose_distro)print*,' finished distro on fragments ',ifrag,ffrag
			 
         end do  ! ffrag
      end do ! ifrag
	 if(iproc==0 .and. verbose_distro)print*,' Finished looping over frag2frag ',successall
	  
      if(successall)then   ! one final round with fill the opbundles
		 if(nob < nob_draft)then  ! some error has occurred
			 if(iproc==0)print*,' Error in splitting opbundles , ended with only ',nob,' out of ',nob_draft
			 stop
		 end if 
		  
         nopbundles = nob
         if(iproc==0)then
            print*,' '
            print*,nopbundles,' final opbundles (A)'
            print*,' '
         end if
!..............   ALLOCATE FINAL OPBUNDLES................
        if(nob > nob_reserved)then  ! need to reserve more opbundles
           if(.not. allocated(opbundle)) then
              print *, "opbundle deallocated without prior allocation"
           end if
           deallocate(opbundle)
           nob_reserved = int(1.2*nob)
           allocate(opbundle(nob_reserved), stat=aerr )
           if(aerr /= 0) call memerror("distro_opbundles_over_fragments 21")
         end if
         myprocstart = 0
         nprocs_avail = nprocs - nprocs_reserved_tot
         call divide_procs_over_fragments(nprocs_avail,.true.)
         nob = 0
         do ifrag = 1,nfragments
            do ffrag = 1,nfragments
               opfragstat(ifrag,ffrag)%nnodes= opfragstat(ifrag,ffrag)%nnodes+ nprocs_reserved(ifrag,ffrag)
            
               Nmyprocs = opfragstat(ifrag,ffrag)%nnodes  ! value computed in divide_procs_over_fragments
               if(Nmyprocs==0)cycle
               call split_draft_opbundles_frag(.true.,nob_draft,ifrag,ffrag,Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag),nob)
               if(setjumpceiling)then
                  call rebalance_MPI_distro(ifrag,ffrag, successlocal,Nmyprocs,myprocstart,nprocs_reserved(ifrag,ffrag))
               end if
               myprocstart = myprocstart+Nmyprocs
            end do  ! ffrag
         end do ! ifrag
!................ CHECK TO SEE IF WE LOST ANY...........
         call check_ops('SPE',nob_draft)
         call check_ops('PP ',nob_draft)
         call check_ops('NN ',nob_draft)
         call check_ops('PN ',nob_draft)
         call check_ops('PPP',nob_draft)
         call check_ops('PPP',nob_draft)
         call check_ops('PPN',nob_draft)
         call check_ops('PNN',nob_draft)
         call check_ops('NNN',nob_draft)

         exit

      else
		 if(iproc==0 .and. verbose_distro)print*,' Failed, need to try again'
		  

      end if
    end do
	
    if(.not.successall)then
       if(iproc==0)then   ! TRY TO END GRACEFULLY
          print*,' FAILED to distribute '
       end if

       if(menu_char/= 'm')then
		   if(iproc==0)print*,' Ending run '
#ifdef _MPI	
		  call BMPI_Barrier(MPI_COMM_WORLD,ierr)
		   
          call BMPI_Abort(MPI_COMM_WORLD, 106, ierr)
#endif
       end if
    end if
!........ PRINTOUT INTRONS.......

    if(compactjumpstorage) call intron_master(.false.,.true.)   ! PRINT OUT INTRONS

!...... ALL FINISHED; PUBLISH THE DISTRIBUTION....

   return

end subroutine distro_opbundles_over_fragments_new

!=================================================================
! subroutine draft_opbundles MOVED TO bparallel_opbundles.f90

!===============================================================
!
! routine to split of opbundles within a block of fragments
! and then rebalance if necessary;
! checks if additional nodes needed for distributing jumps
!
! INPUT
!   ifrag,ffrag :: initial, final fragments
!
! OUTPUT:
!   logical :: success
!   nnodesneeded :: how many additional nodes needed if not success
!
! CALLED BY:
!   distro_opbundles_over_fragments
!
! SUBROUTINES CALLED:
!   setMPIjumplimits
!   count_jumps_on_node
!   redistributor_frag
!
subroutine rebalance_MPI_distro(ifrag,ffrag, success,Nmyprocs,myprocstart,myreserveprocs)

   use opbundles
   use fragments
   use verbosity
   use nodeinfo
   use operation_stats
   use flagger
   use bmpi_mod
   use jumpNbody
   implicit none

!..... INPUT........
   integer ifrag,ffrag
   integer :: Nmyprocs,myreserveprocs,myprocstart

!...... OUTPUT........
   logical :: success
!........ INTERNAL......
   integer nob,nob_draft
   real(4) :: localwt
   integer(8) :: localjumps,maxlocaljumps
   real(8) :: localjumpstorage
   integer iprocs
   integer :: nprocsoverbooked, nprocsexcess

   if(opfragstat(ifrag,ffrag)%nnodes == 0)return
   if(nprocs == 1)then   ! CAN'T REBALANCE, GIVE UP AND RETURN
	   success = .true.
	   return
   end if

   call setMPIjumplimits(.false.)  ! this includes accounting for "introns" or unused jumps
   nprocsoverbooked = 0
   nprocsexcess     = 0
   maxlocaljumps = int( maxjumpmemory_default/ (bytesper2Bjump*1.e-9),8 )
   do iprocs = myprocstart,myprocstart+Nmyprocs-1
         call count_jumps_on_node(iprocs,0,localjumps,localjumpstorage)  ! if intron cuts are enabled, this accounts for them
         if( localjumps > maxlocaljumps) then
			 
	         if(iproc==0)print*,' rebalancing, found MPI rank ',iprocs,' overruns jumps by a factor of ',& 
		          real(localjumps,8)/real(maxlocaljumps,8)
			 
             nprocsoverbooked = nprocsoverbooked +1
             nprocsexcess  = nprocsexcess + int(localjumps/maxlocaljumps+1, 4)
         end if
   end do
   success = .true.
   if( nprocsexcess - nprocsoverbooked > 0)then  ! need to redistribute
      if ( nprocsexcess - nprocsoverbooked <= myreserveprocs)then
           call redistributor_frag(Nmyprocs,myprocstart,ifrag,ffrag)
      else
           myreserveprocs = nprocsexcess - nprocsoverbooked
           success = .false.
      end if
   end if
   return

end subroutine rebalance_MPI_distro

!===============================================================
!
! this routine reads through the opbundles and splits them 
! when a max number of ops has been reached
!
! splits for a particle initial and final fragment
!
! uses data gathers in routine analyze_draft_opbundle
! splits along "nsetops"
!
! INPUT:
!   fill: flag to actually populate the final opbundle
!   nob_draft:  # of draft
!   ifrag,ffrag : initial/final fragment labels
!   myprocs   : MPI processes (nodes) assigned 
!   myreserveprocs : reserved MPI processes/nodes
!
! INPUT/OUTPUT: nob = (updated) # of opbundles
!
! CALLED BY:
!          distro_opbundles_over_fragments
!
! SUBROUTINES CALLED: 
!         search_draft_opbundles4split in baparallel_opbundles.f90
!

subroutine split_draft_opbundles_frag(fill,nob_draft,ifrag,ffrag,myprocs,myprocstart,myreserveprocs,nob)
   use opbundles
   use nodeinfo
   use operation_stats
   use io
   use bmpi_mod
   use butil_mod
   implicit none

!............. INPUT
   logical, intent(in) :: fill
   integer nob_draft   ! appears unused
   integer, intent(in) :: ifrag,ffrag
   integer myprocs,myprocstart,myprocstop,myreserveprocs    ! # of MPI processes assigned to this set of fragments
!..............OUTPUT.......
   integer(4) :: nob,dnob   ! change in the # of opbundles
!........... INTERNAL.....
   integer :: nodestartob(0:nprocs), nodestopob(0:nprocs)    ! which opbundles to start/stop on
   integer(8) :: nodestartz(0:nprocs), nodestopz(0:nprocs)  !which operations to start/stop on
   integer(8) :: localops,localjumps
   integer(4):: iob
   integer jnode
!   integer lastob
   integer(8) :: startz,stopz,laststart
   integer(8) :: maxlocalops
   integer iunit
   integer ierr
   logical,parameter :: localverbose = .false.

   integer(4) :: my_nob_start, my_nob_end
   
   if(localverbose)then
	   print*,' # '
	   print*,ifrag,ffrag,opfragstat(ifrag,ffrag)%nnodes
	   print*,' # '
   end if
   if(opfragstat(ifrag,ffrag)%nnodes ==0)return  ! no operations here
   dnob = nob   ! will be used in second half
!..........................

   if(myprocs <= myreserveprocs)then
     if(iproc==0)then
         print*,' trying to reserve additional MPI procs which do not exist ',myprocs, myreserveprocs
         print*,' on fragments ',ifrag,' --> ',ffrag
         print*,' Probably may (?) be not enough nodes to store jumps '
     end if
     stop
   end if
   my_nob_start = frag_draftopbundlestart(ifrag,ffrag)   ! start of opbundles assigned to this frag2frag
   my_nob_end   = frag_draftopbundleend(ifrag,ffrag)     ! end of opbundles assigned to this frag2frag
   if(iproc==0 .and. localverbose)then
     print*,' from process ',myprocstart, ' to ',myprocs+myprocstart-1
     print*,' draft opbundles from ',my_nob_start,' to ',my_nob_end
   end if
   if(my_nob_end < my_nob_start)then
      print*,' problem with fragments ',ifrag,ffrag
      print*, my_nob_start,my_nob_end
      stop
   end if
   nopsnode = int(frag2fragops(ifrag,ffrag)/real(myprocs-myreserveprocs,8),8)
   if(iproc==0 .and. localverbose )then
          print*,' On fragment ',ifrag,'--> ',ffrag,' approx ',nopsnode,' operations per MPI process '
   end if
   do jnode = myprocstart,myprocstart +myprocs-1
      opbundlestart(jnode) = my_nob_start
      opbundleend(jnode) = -1
   end do 
   nodestartob(myprocstart) = my_nob_start
   nodestartz(myprocstart)  = 1

   maxlocalops = 0

!....... LOOP OVER # OF PROCESSORS;
!        FOR EACH PROCESSOR COUNT UP OPERATIONS
!        UNTIL REACHED TARGET, AND THEN SPLIT

   myprocstop = myprocstart+ myprocs-myreserveprocs-1 
   do jnode = myprocstart,myprocstop  ! loop over nodes
	   if(nodestartob(jnode)> my_nob_end)then
		   print*,' some problem with splitting on proc',jnode, myprocstart,myprocstop
		   print*,' between fragments ',ifrag,ffrag
		   print*,my_nob_end,nodestartob(jnode)
		   stop
	   end if
	   
      call search_draft_opbundles4split_g(my_nob_end,nodestartob(jnode),  &   ! search for a suitable place to split an opbundle
             nodestopob(jnode) ,nodestartz(jnode),nodestopz(jnode))
      if( nodestopz(jnode) == draft_opbundle( nodestopob(jnode))%nsetops)then
           nodestartob(jnode+1) = nodestopob(jnode)+1
           nodestartz(jnode+1) = 1
      else
           nodestartob(jnode+1) = nodestopob(jnode)
           nodestartz(jnode+1) = nodestopz(jnode) + 1
      end if

!......... MAKE SURE ALL OPBUNDLES INCLUDED.....
      if(jnode == myprocstart+myprocs-myreserveprocs -1 .and. nodestopob(jnode)/= my_nob_end)then
          nodestopob(jnode) = my_nob_end
          nodestopz(jnode) = draft_opbundle(my_nob_end)%nsetops
      end if

!.......... ADD TO NUMBER OF "NEW" OPBUNDLES ........
      if(iproc==0 .and. localverbose)print*, ' process ',jnode, ' assigned opbundle range: ', nodestartob(jnode),nodestopob(jnode)

      nob = nob + nodestopob(jnode) - nodestartob(jnode)+1
      if(nodestopob(jnode)> my_nob_end)then
          if(iproc==0)print*,jnode,nodestopob(jnode),' too many !! '
          stop
      end if
	  
      if(nodestartob(jnode+1) == my_nob_end+1 .and. jnode -myprocstart < myprocs-myreserveprocs-1)then   ! have run out bundles
		  if(iproc==0)then
			  print*,myprocs+myprocstart-myreserveprocs-1-jnode,' procs left over between fragments ',ifrag,ffrag
			  
		  end if
		  myprocs = myprocs - (jnode -myprocstart - myprocs+myreserveprocs+1 )   ! reduce the number of processors allocated
		  exit
		  
	  end if
	  
      if(nodestartob(jnode+1) > my_nob_end+1 .and. jnode -myprocstart < myprocs-myreserveprocs-1)then  ! FIXED A BUG HERE 7.9.3
		  ! what happens here: trying to assign bundles beyond what is allowed
          print*,' gone toooooo far '
          print*,jnode+1,nodestartob(jnode+1), my_nob_end
		  print*,jnode,myprocstart, myprocs-myreserveprocs-1
		  print*,nodestopz(jnode), draft_opbundle( nodestopob(jnode))%nsetops,nodestopob(jnode)
          stop
      end if
!......... check not missing any
      if(jnode == myprocstart+myprocs-myreserveprocs-1 .and. nodestopz(jnode) /= draft_opbundle(my_nob_end)%nsetops)then
         if( real(draft_opbundle(my_nob_end)%nsetops,kind=8) -real(nodestopz(jnode),kind=8 )< .5*nopsnode)then
            nodestopz(jnode) = draft_opbundle(my_nob_end)%nsetops
         else
            if(iproc==0)then
				print*,jnode,' missing ops from ', nodestopz(jnode),' to ', draft_opbundle(my_nob_end)%nsetops
				print*,my_nob_end, draft_opbundle(my_nob_end)%optype
			end if
#ifdef _MPI	
			call BMPI_Barrier(MPI_COMM_WORLD,ierr)
            call BMPI_ABORT(MPI_COMM_WORLD,101,ierr)
#endif
            stop
         end if
      end if 

   end do    ! jnode

   if(iproc==0 .and. localverbose) print*,' progress in splitting ',my_nob_start,my_nob_end,dnob,nob
   if(.not.fill)return

!..................... SET UP FINAL OPBUNDLES ................
!                      ONLY DONE IF NOT MODELING
!                      ONLY DONE FOR NODES THAT NEED IT
   
   do jnode = myprocstart,myprocstart+myprocs-myreserveprocs -1 
      opbundlestart(jnode) = dnob+1
      localops = 0
      do iob = nodestartob(jnode),nodestopob(jnode)
         dnob = dnob + 1
         opbundle(dnob)%node = jnode
         opbundle(dnob)%optype = draft_opbundle(iob)%optype
         opbundle(dnob)%hchar  = draft_opbundle(iob)%hchar
         opbundle(dnob)%cstride= draft_opbundle(iob)%cstride
         opbundle(dnob)%ifragment  = draft_opbundle(iob)%ifragment
         opbundle(dnob)%ffragment  = draft_opbundle(iob)%ffragment
         opbundle(dnob)%isector  = draft_opbundle(iob)%isector
         opbundle(dnob)%fsector  = draft_opbundle(iob)%fsector
         opbundle(dnob)%insector  = draft_opbundle(iob)%insector
         opbundle(dnob)%fnsector  = draft_opbundle(iob)%fnsector
         opbundle(dnob)%numpthreads  = draft_opbundle(iob)%numpthreads  !actually  = 0
         opbundle(dnob)%numnthreads  = draft_opbundle(iob)%numnthreads
         opbundle(dnob)%min_nop = draft_opbundle(iob)%min_nop
		 opbundle(dnob)%origin  = iob    ! added 7.9.6 so we can determine where this opbundle came from
		 opbundle(dnob)%opwt    = draft_opbundle(iob)%opwt
		 opbundle(dnob)%greedy  = draft_opbundle(iob)%greedy

         if(iob == nodestartob(jnode))then
            startz = nodestartz(jnode)
         else
            startz = 1
         endif
         if(iob == nodestopob(jnode))then
            stopz = nodestopz(jnode)
            if(stopz < 0)then
               print*,' HEY WAIT A MINUTE ',jnode
            end if
         else
            stopz = draft_opbundle(iob)%nsetops
            if(stopz < 0)then
               print*,' HEY WAIT A SECOND ',iob
            end if
         endif

         select case ( draft_opbundle(iob)%optype )
            case('PP ', 'PPP','PPN', 'SPE')
               opbundle(dnob)%pxstart = draft_opbundle(iob)%pxstart+startz-1
               opbundle(dnob)%pxend   = draft_opbundle(iob)%pxstart+stopz-1
               opbundle(dnob)%nxstart = draft_opbundle(iob)%nxstart
               opbundle(dnob)%nxend   = draft_opbundle(iob)%nxend
               opbundle(dnob)%nsortstart = draft_opbundle(iob)%nsortstart
               opbundle(dnob)%nsortend   = draft_opbundle(iob)%nsortend

            case('NN ','PN ','NNN','PNN','CEN')
!		   if(dnob==22)then
!			   print*,' testing copying '
!			   print*,startz,stopz, draft_opbundle(iob)%nxend
!		   end if
               opbundle(dnob)%nxstart = draft_opbundle(iob)%nxstart+startz-1
               opbundle(dnob)%nxend   = draft_opbundle(iob)%nxstart+stopz-1
               opbundle(dnob)%pxstart = draft_opbundle(iob)%pxstart
               opbundle(dnob)%pxend   = draft_opbundle(iob)%pxend
               opbundle(dnob)%nsortstart = draft_opbundle(iob)%nsortstart
               opbundle(dnob)%nsortend   = draft_opbundle(iob)%nsortend
         end select
         localops = localops + int(stopz-startz+1,8)*int(draft_opbundle(iob)%min_nop,8)
      end do
      maxlocalops = bmax(maxlocalops, localops)
      opbundleend(jnode) =dnob

      if(iproc==0)then
        do iob = opbundlestart(jnode),opbundleend(jnode)
          select case ( opbundle(iob)%optype )
             case('SPE')
                localjumps = 0
             case('PP','PPP')
                localjumps = int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)
             case('NN','NNN')
                localjumps = int( opbundle(iob)%nsortend -opbundle(iob)%nsortstart+1,8)
             case('PN','PPN','PNN','CEN')
                localjumps = int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)    & 
               + int( opbundle(iob)%nsortend -opbundle(iob)%nsortstart+1,8)

          end select
        end do
      endif 
   end do   ! jnode
   if(iproc==0 .and. localverbose) print*,' comparing numbers ',nob,dnob

   return
end subroutine split_draft_opbundles_frag

!===============================================================
!
!  when the number of jumps predicted to be stored on a node go over a limit,
!  (maxjumpmemory_default, set in module flagger in bmodules_flags.f90 )
!  this routine redistributes excess jumps onto additional nodes. 
!  Simplest version, assumes no further breakup of opbundles is needed
!  (Oct 2013)
!  Modified Aug 2014 to include fragmented Lanczos vectors
!  CALLED BY
!    rebalance_MPI_distro
!
!  CALLS
!    count_jumps_on_node
!    get_jumpcount
!
subroutine redistributor_frag(Nmyprocs,myprocstart,ifrag,ffrag)

   use opbundles
   use fragments
   use verbosity
   use nodeinfo
   use operation_stats
   use flagger
!   use tribution
   use bmpi_mod
   use jumpNbody
   implicit none

   integer :: Nmyprocs   ! # of MPI processes assigned
   integer :: myprocstart ! where MPI processes start
   integer :: ifrag,ffrag

   integer(8) :: localjumps,maxlocaljumps,targetlocaljumps,jumpsofar,oldlocaljumps
   integer(8),allocatable :: originaljumps(:)  ! for debugging purposes, store original setting of jumps
   integer iprocs,jprocs,kprocs
   integer :: nnodesoverbooked, nnodesexcess,nreservetest
   integer :: iob
   integer(8) :: jumpstartX,jumpendX,jumpstartY,jumpendY
   character(3) :: current_optype
   logical :: origproc    ! flag to denote processor where jumps originally stored
   integer :: bundlestart,bundlestop
   integer :: ierr
   integer :: aerr
   real(8) :: opjumpstorage,localjumpstorage
   
   
   maxlocaljumps = int( maxjumpmemory_default/ (bytesper2Bjump*1.0e-9) )   ! computes # of jumps to store on a node

   nreservetest = 0

   allocate (originaljumps(myprocstart: myprocstart+Nmyprocs-1), stat=aerr)
   if(aerr /= 0) then
      call memerror("redistributor_frag")
      stop 5
   end if

   do iprocs = myprocstart,myprocstart+Nmyprocs-1         ! count up how many reserve nodes there are
                                    ! reserve nodes do not yet have jumps/opbundles assigned to them
      call count_jumps_on_node(iprocs,0,localjumps,localjumpstorage)
      originaljumps(iprocs) = localjumps
      if(localjumps ==0)nreservetest = nreservetest+1  ! if no localjumps yets, then it is a reserve node
   end do
   if(iproc == 0 .and. verbose_distro_report) print *, "Number of reserved(unused) nodes=", nreservetest

   jprocs=nprocs -nreservetest-1

   do iprocs = myprocstart,myprocstart+Nmyprocs -1-nreservetest   ! loop over all nodes with opbundles/jumps assigned

      call count_jumps_on_node(iprocs,0,localjumps,localjumpstorage)

      if( localjumps > maxlocaljumps) then     ! too many jumps on this node, redistribute
             nnodesexcess  = int(localjumps/maxlocaljumps, kind=4)
             targetlocaljumps = localjumps/(nnodesexcess+1)  ! roughly how many jumps to store locally
!....................... SPLIT UP.............................................

             current_optype = opbundle( opbundlestart(iprocs) )%optype
             origproc = .true.  
!............... LOOP OVER AND COUNT UP JUMPS.....
             bundlestart = opbundlestart(iprocs)
             bundlestop  = opbundleend  (iprocs)
!................ INITIALIZE JUMP COUNT ..........
             oldlocaljumps = localjumps
             call get_jumpcount(bundlestart,jumpstartX,jumpendX,jumpstartY,jumpendY,.true., localjumps)   
! KSM:   change comparison from localjumps > targetlocaljumps to localjumps > maxlocaljumps
! I was getting errors like    10 > 16 GB, problem cannot divide ..
! Seems like we want to know if we are under the max for this chunk

             if(localjumps > maxlocaljumps)then   ! ERROR TRAP -- SIMPLE DIVISION WILL NOT WORK
                 if(iproc==0)then
                     print*,' problem cannot divide up simply over ',Nmyprocs,' MPI procs '
                     print*,' storing jumps on proc ',iprocs,' requires ', & 
                        localjumps*bytesper2Bjump*1.0e-9,' > allowed ',maxjumpmemory_default,' Gb '
                     print*,' If possible change maxjumpmemory_default '
					 print*,' Bundle ',bundlestart,' is the problem '
					 print*,opbundle(bundlestart)%optype,' P: ',opbundle(bundlestart)%pxend-opbundle(bundlestart)%pxstart+1, & 
					 ', N: ',opbundle(bundlestart)%nsortend-opbundle(bundlestart)%nsortstart+1
                  end if
#ifdef _MPI	
				  call BMPI_Barrier(MPI_COMM_WORLD,ierr)
                  call BMPI_Abort(MPI_COMM_WORLD, 102, ierr)
#endif
                  stop
             end if  

             jumpsofar = 0           
             do iob = bundlestart,bundlestop
                
                if( current_optype /= opbundle(iob)%optype)then  ! reset jumpstarts
                    jumpsofar = jumpsofar+localjumps
                    call get_jumpcount(iob,jumpstartX,jumpendX,jumpstartY,jumpendY,.true., localjumps)                
                    current_optype = opbundle(iob)%optype
                end if

                call get_jumpcount(iob,jumpstartX,jumpendX,jumpstartY,jumpendY,.false.,localjumps)
                if(jumpsofar+localjumps > maxlocaljumps)then    ! hit reset
                    jprocs = jprocs + 1
                    if (jprocs > nprocs -1)then   ! error trap
                        if(iproc==0)then
                            print *, ' oops, too many nodes needed ',current_optype
							print*,  " fragments ",ifrag,'-> ',ffrag
                            print *, "nprocs=", nprocs
							print *, "procs reserved ",Nmyprocs,myprocstart
                            print *, "jprocs=", jprocs
                            print *, "iob=", iob
                            print *, "jumpssofar=", jumpsofar
                            print *, "localjumps=", localjumps
                            print *, "maxlocaljumps", maxlocaljumps
							call analyze_opbundle(iob,.false.,opjumpstorage)
							print*, " operations = ", & 
						(opbundle(iob)%pxend-opbundle(iob)%pxstart+1)*(opbundle(iob)%nxend-opbundle(iob)%nxstart+1)
						    print*,' jump storage = ',opjumpstorage
                        end if
                        STOP
                    end if
                    if(origproc)then
                        opbundleend(iprocs)   = iob-1     ! end of storage on original 
                        opbundlestart(jprocs) = iob     ! start of storage on next available
                        origproc = .false.
                     else
                        opbundleend(jprocs-1) = iob-1    !end of storage on previous
                        opbundlestart(jprocs) = iob     ! start of storage on next available
                     end if
                     jumpsofar = 0
                     call get_jumpcount(iob,jumpstartX,jumpendX,jumpstartY,jumpendY,.true.,localjumps)
                     jumpsofar = 0
                end if
             end do
             opbundleend(jprocs) = bundlestop   !  finish up
      end if
      
   end do
   deallocate(originaljumps)
   return
end subroutine redistributor_frag

!===============================================================
!
! CALLED BY master_para_distro
!
subroutine print_memory_limits
	use flagger
	use nodeinfo
	use jumpNbody
	use menu_choices
	implicit none
	
	if(nproc==1 .and. menu_char/='m')return
	if(iproc==0)then
		print*,' - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
		print*,' Memory limits per (MPI) processor: '
		print*,' Storage of all Lanczos vectors: ', maxlanczosstorage1,' Gb '
		print*,'(see  maxlanczosstorage1 in module flagger )'
		print*,' Storage of jumps: ',maxjumpmemory_default ,' Gb '
		print*,' or ~ ',maxjumpmemory_default /real( bytesper2Bjump,8),' billion jumps '
		print*,' (see maxjumpmemory_default in module flagger )'
		print*,' - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
		print*,' '
	end if
	
	
end subroutine print_memory_limits

!===============================================================

end module para_main_mod


   

