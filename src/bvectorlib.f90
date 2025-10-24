!========================================================================
!
!  BVECTORLIB.f90
!
!  vector routines for BIGSTICK
!
!  routines to carry out vector algebra and storage
!  including for parallel implementations  
!  initiated June 2013 @ SDSU
!
module bvectorlib_mod
contains
!=======================================================================
!  SUBROUTINES IN THIS FILE
!
!

!
!  ALLOCATE "LOCAL" FRAGMENTS OF VECTORS
!
!  here basestart(f), basestop(f) are start and stop
!  of the basis states for fragment f
!
!  Default for 1 processor is frag1 = frag2 = 1 (only one "fragment")
!
! Also builds fragment level communicators
!
! CALLED BY
!    setup_for_lanczos
!    exactdiag_p
!    density1b_from_oldwfn
!    overlap
!    applicator_h   ( apply scalar 2-body opeator)
!    applicator1b   ( apply nonscalar 1-body operator)
!    particle_occupation_p
!    particle_occupation_p_orig
!
!  SUBROUTINES CALLED:
!   setup_mpi_hist
!   memreport
!
subroutine setup_localvectors
   use flagger
   use localvectors
   use localblocks
   use io
   use nodeinfo
   use fragments
   use lanczos_info
   use basis
!   use tribution
   use mod_reorthog
   use bmpi_mod
   use butil_mod
   implicit none

   integer :: key, ierr
   integer :: aerr
   integer :: ii
   integer :: num_threads
   integer(kind=basis_prec) :: vi
   integer :: tid

   ! OpenMP functions
   integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
         		 
   frag1 = nodal(iproc)%ifragment
   frag2 = nodal(iproc)%ffragment
   
   if(frag1 < 1 .or. frag1 > nfragments) then ! This nfragments reference is ok
      print *, "setup_localvectors: ", iproc, " out of range, frag1=", frag1
      write(logfile,*) "setup_localvectors: ", iproc, " out of range, frag1=", frag1
!      flush(logfile)
      stop 1
   end if
   if(frag2 < 1 .or. frag2 > nfragments) then ! This nfragments reference is ok
      print *, "setup_localvectors: ", iproc, " out of range, frag2=", frag2
      write(logfile,*) "setup_localvectors: ", iproc, " out of range, frag2=", frag2
!      flush(logfile)
      stop 1
   end if
   
   ! use to access slice regardless of how vec1/vec2 are allocated
   v1s = basestart(frag1)
   v1e = basestop(frag1)
   v2s = basestart(frag2)
   v2e = basestop(frag2)
!      write(logfile,*), "setuplocalvectors: ", iproc, ", frag1=", frag1, ", frag2=", frag2
!      write(logfile,*), "setuplocalvectors: ", iproc, ", v1s=", v1s, ", v1e=", v1e, ", v2s=", v2s, ", v2e=", v2e 

   if(v1e < v1s .or. v2e < v2s) then
      print *, "setuplocalvectors: ", iproc, ", range error ve < vs, stopping"
      stop 1
   end if
   if(.not. useNewReorthog)then
      allocate( vec1 ( dimbasis ), stat=aerr )
      if(aerr /= 0) call memerror("setup_localvectors 1")
      allocate( vec2 ( dimbasis ), stat=aerr )
      if(aerr /= 0) call memerror("setup_localvectors 2")
   else
 	  ! this crashes when using option TW in MPI
      allocate( vec1(v1s:v1e), stat=aerr )       ! holds previous vector  
      if(aerr /= 0) call memerror("setup_localvectors 3")
      allocate( vec2(v2s:v2e), stat=aerr )       ! holds new vector	  
      if(aerr /= 0) call memerror("setup_localvectors 4")
   end if

   ! allocate memory for each thread to write to in OpenMP mode
!$omp parallel
   num_threads = omp_get_num_threads();
!$omp end parallel
   ompNumThreads = num_threads
   useVec2Thread = useNewReorthog .and. (ompNumThreads > 1) .and. wantUseVec2Thread
      
   if(useVec2Thread) then
      if(iproc == 0) print *, "Allocating vec2thread, threadcount=", ompNumThreads
      !! allocate(vec2thread(v2s:v2e, 0:ompNumThreads-1), stat=aerr)
      vec2threadchunk = v2e - v2s + 1          ! size of chunk
      vec2threadchunkm1 = vec2threadchunk - 1  ! to find ent of slice
      vec2threadend = ompNumThreads * (v2e - v2s + 1) - 1  ! end of buffer
      allocate(vec2threadflat(0 : vec2threadend), stat=aerr)
      if(aerr /= 0) call memerror("setup_localvectors vec2thread - set wantUseVec2Thread=.false. to save memory");
      !! vec2thread(:, :) = 0.0
      vec2threadflat(:) = 0.0
   end if

   if(storevectorsenabled .and. nproc>=niter)then
      writetodisk = .false.
   else
      writetodisk = .true.
   end if

   ! Seems like a good place to set up history mechanism
!   if(iproc == 0) print *, "KSM:  Calling setup_mpi_hist(", niter, ")"

   if(block_flag==0)then
      call setup_mpi_hist(niter)
   else
	   call setup_mpi_hist(niter+dimblock)
   end if
   
   if(noisy0) print *, "Returned from setup_mpi_hist"
   call memreport     ! generates a report on storage ADDED 7.3.8
   if(noisy0) print *, "Generated memory report"

   ! KSM - 15Aug2014
   ! Build communicators for allreduce and for broadcast from root node of each fragment
   ! will work fine even if only one fragment
   key = 1
   if(isfragroot) key = 0  ! root has to get rank 0
   ! These communicators are split by the fragment they serve
#ifdef _MPI			     
   call BMPI_COMM_SPLIT(MPI_COMM_WORLD, frag2, key, fcomm2, ierr)
   call BMPI_COMM_SPLIT(MPI_COMM_WORLD, frag1, key, fcomm1, ierr)
   if(noisy0) print *, "Done making fcomm1, and fcomm2"
#endif   
!....these are indices to tell MPI which fcomm to use
   fcomm1_index = 1
   fcomm2_index = 2

   return
end subroutine setup_localvectors
!
!=======================================================================
!  double-precision normalization
!
!  -- for 'new' parallelization
!
! vchar :  information on vector (normal or reverse ordering
! ipoint:  point to 'i' (initial) or 'f' (final) vector
!
! dnorm : double-precision norm of dvec
! smallflag: logical flag = .true. if norm is smaller than dtol
!
! CALLED BY:
!    initialize_lanczos_vector
!    lanczos_p
!    random_restart_p
!    expectator_p
!==========================================================
subroutine dnormvec_p(vchar,ipoint,dnorm,smallflag)
  use precisions
  use nodeinfo
  use timing
  use localvectors
  use fragments
  use bmpi_mod
  implicit none

  integer(4) :: ierr
  character(1) :: vchar,ipoint
  real(kind=lanc_prec), pointer :: dvec(:)
  real(kind=8) :: dnorm
  logical smallflag
!------------------ INTERMEDIATE -------------------------------------------
  real(kind=8) :: d
  real(kind=8) :: dtol
  real(kind=8) :: tdnorm

  integer(4)   :: i

  integer(kind=basis_prec) :: il,vstart,vstop

  call clocker('dot','sta')
  dtol = 1.d-8
  dnorm = 0.d0
  smallflag = .false.

!..... choose which vector to normalize
!      complicated because during lanczos 
!      we switch between vec1 and vec2 (stored in module localvectors
!      also can normalize initial or final vector

  if( (vchar == 'n' .and. ipoint == 'i' ) .or. & 
        (vchar == 'r' .and. ipoint == 'f') )then
        dvec => vec1
        vstart = v1s
        vstop  = v1e
  else
        dvec => vec2
        vstart = v2s
        vstop  = v2e
  end if
  
!$omp parallel do private(il,d), shared(dvec, vstart,vstop), reduction(+ : dnorm) 
   do il= vstart,vstop
        d = real(dvec(il),kind=8)
        dnorm = dnorm+d*d
   end do
!$omp end parallel do

!..... IF USING MPI with fragments, MUST COMBINE HERE
! KSM:  We only want to provide a result from the nodes
!       labled isfragroot.  This selects one representative
!       for each fragment.  Otherwise we would be double counting.
! This routine should only be used sparingly because it wastes
! many nodes.   On the other hand, the cost of distributing the
! vector slices using breorthog.f90 can be high.
!
  if(.not. isfragroot) dnorm = 0.d0
#ifdef _MPI			      
  call BMPI_ALLREDUCE(dnorm, tdnorm, 1,  MPI_SUM, MPI_COMM_WORLD, ierr)
#else
tdnorm = dnorm
#endif
  dnorm = tdnorm

  if ( dnorm < dtol ) then
     smallflag = .true.
     if ( iproc == 0 )write(6,*)' zero vector ',dnorm
     return
  end if

  dnorm = dsqrt(dnorm)
  d = 1.d0 / dnorm

!$omp parallel do private(il), shared(dvec, d, vstart,vstop)
     do il = vstart,vstop
        dvec(il) = real(dvec(il)*d,kind=lanc_prec)
     end do
!$omp end parallel do

  call clocker('dot','end')
  return
end subroutine dnormvec_p

!=================================================================
!  double-precision projection
!
!
! n: dimension of vector
! dvec1: double-precision vector
! dvec2: double-preciscion vector
!
! dsclrprod  = dvec1*dvec2
!
! dvec1 -> dvec1 - dvec2*dsclrprod
! returns \alpha - the overlap with the previous vector
!
! CALLED BY: reorthogonalize_a
!====================================================================
subroutine dvecproj_p(vchar,dsclrprod)
  use precisions
  use nodeinfo
  use timing
  use lanczos_info
  use localvectors
  use fragments
  use bmpi_mod
  implicit none

  character(1) :: vchar
  integer(4)               :: ierr
  real(kind=8)             :: dsclrprod
  real(kind=lanc_prec), pointer :: dvec1(:),dvec2(:)

!------------------ INTERMEDIATE -------------------------------------------
  real(kind=8)             :: d1,d2
  integer(4)               :: i

  integer(kind=basis_prec) :: il,vstart,vstop

  if(nfragments > 1) then   ! This nfragments reference is ok, not about new reorthog
     if(iproc == 0) print *, "dvecproj_p: not supported with new reorthog"
     stop 1 ! has msg
  end if

  if(.not.storelanczosincoreMPI .and. iproc > 0)return
!...............................

  call clocker('pro','sta')
  dsclrprod = 0.d0

       ! This code doesn't work when frag1 /= frag2
       vstart = v1s  !    basestart(frag1)
       vstop  = v1e  !    basestop (frag1)

  select case (vchar)
          case ('n')
            dvec1 => vec1
            dvec2 => vec2

          case ('r')
            dvec1 => vec2
            dvec2 => vec1

        case default
           stop 1   ! prevent uninitialized compiler msg
  end select

!$omp parallel do private(il, d1, d2), shared(vstart,vstop, dvec1, dvec2), reduction(+ : dsclrprod)
  do il =vstart,vstop
        d1 = real(dvec1(il),kind=8)
        d2 = real(dvec2(il),kind=8)
        dsclrprod = dsclrprod+d1*d2
  end do
!$omp end parallel do

!$omp parallel do private(il, d1, d2), shared(vstart,vstop, dvec1, dvec2,dsclrprod)
  do il = vstart,vstop
        d1 = real(dvec1(il),kind=8)
        d2 = real(dvec2(il),kind=8)
        dvec1(il) = real(d1 - d2*dsclrprod,kind=lanc_prec)
  end do
!$omp end parallel do
  call clocker('pro','end')
  return
end subroutine dvecproj_p

!===================================================
      subroutine dvecdot_p(vchar,dsclrprod)
!
!  double-precision projection
!
!  vchar: if = 'n' vec1 = initial, vec2 = final
!         if = 'r' vec2 = initial, vec1 = final
!
! dsclrprod  = dvec1*dvec2
!
!
      use precisions
	use nodeinfo
      use localvectors
      use fragments
      implicit none
      character(1) :: vchar
      real(kind=lanc_prec), pointer :: dvec1(:),dvec2(:)
      real(kind=8) :: dsclrprod

!------------------ INTERMEDIATE -------------
      real(kind=8) :: d1,d2
      integer      :: i

      integer(kind=basis_prec) :: il,vstart,vstop

   if(nfragments > 1) then  ! this nfragments reference is ok
      if(iproc == 0) print *, "dvecdot_p: not supported with nfragments>1"
      stop 1
   end if

	! This code didn't work before when frag1 /= frag2
	vstart = v1s  !    basestart(frag1)
	vstop  = v1e  !    basestop (frag1)

   select case (vchar)
    case ('n')
      dvec1 => vec1
      dvec2 => vec2

    case ('r')
      dvec1 => vec2
      dvec2 => vec1

    case default
      stop 1   ! prevent uninitialized compiler msg

   end select
   dsclrprod = 0.d0
!$OMP PARALLEL SHARED(dvec1,dvec2,vstart,vstop), PRIVATE(i,d1,d2), REDUCTION(+:dsclrprod)
!$OMP DO SCHEDULE(STATIC)
   do il = vstart,vstop
         d1 = dvec1(il)
         d2 = dvec2(il)
         dsclrprod = dsclrprod+d1*d2
   enddo
!$OMP END DO
!$OMP END PARALLEL

   return
end subroutine dvecdot_p

!================================================================
      subroutine doverlapvec(dsclrprod,dotflag)
!
!  double-precision overlap
!  added in 7.7.0
!
! dsclrprod  = dvec1*dvec2
!
! CALLED BY overlap
!
      use precisions
	  use nodeinfo
      use localvectors
      use fragments
	  use bmpi_mod
      implicit none
      real(kind=8) :: dsclrprod,tmpprod
	  logical :: dotflag

!------------------ INTERMEDIATE -------------
      real(kind=8) :: d1,d2
      integer      :: i

      integer(kind=basis_prec) :: il,vstart,vstop
	  integer :: ierr


	vstart = v1s  !    basestart(frag1)
	vstop  = v1e  !    basestop (frag1)

   dsclrprod = 0.d0
   if(isfragroot .or. nproc==1)then
	   
	   if(dotflag)then
!$OMP PARALLEL SHARED(vec1,vec2,vstart,vstop), PRIVATE(i,d1,d2), REDUCTION(+:dsclrprod)
!$OMP DO SCHEDULE(STATIC)
   do il = vstart,vstop
         d1 = vec1(il)
         d2 = vec2(il)
         dsclrprod = dsclrprod+d1*d2
   enddo
!$OMP END DO
!$OMP END PARALLEL
     else ! compute the relative entropy also called the Kullbeck-Leibler divergence
!$OMP PARALLEL SHARED(vec1,vec2,vstart,vstop), PRIVATE(i,d1,d2), REDUCTION(+:dsclrprod)
!$OMP DO SCHEDULE(STATIC)
		    do il = vstart,vstop
		          d1 = vec1(il)
				  d1 = d1*d1
		          d2 = vec2(il)
				  d2 = d2*d2
				  if(d2 < 1.e-9 .or. d1 < 1.e-9)cycle
		          dsclrprod = dsclrprod + d1*log(d1/d2)
		    enddo
!$OMP END DO
!$OMP END PARALLEL		 
		 
	 end if

end if

#ifdef _MPI			  
call BMPI_REDUCE(dsclrprod,1,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif
!call BMPI_BARRIER(MPI_COMM_WORLD,ierr)
!print*,iproc,vstart,vstop,dsclrprod

!dsclprod=tmpprod
   return
end subroutine doverlapvec
!================================================================

! WRITE LANCZOS VECTOR_A
!
!  SHOULD ONLY BE USED TO WRITE TO DISK  
!
! CALLED BY:
!  exactdiag_p      in blamczosmain.f90
!  exactdiag_MPI              ''
!  thick_restart_sub_p   in blanczoslib.f90
!
!================================================================
subroutine write_lanczos_vector_a(vchar,ipoint,i,iunit)
  use lanczos_info
  use precisions
  use basis
  use nodeinfo
  use localvectors
  use fragments
  use flagger
  use bmpi_mod
  implicit none

  character(1) :: vchar,ipoint
  integer(4) :: i,iunit

  real(kind=lanc_prec),pointer   :: v(:)
  integer (kind=basis_prec) :: jl,vstart,vstop

!............................................................ 
#ifdef _MPI	  
type(MPI_status) :: status
!  integer(4)                    :: status(MPI_STATUS_SIZE)
#endif

!.................................................................

!..... choose which vector to normalize
!      complicated because during lanczos 
!      we switch between vec1 and vec2 (stored in module localvectors
!      also can normalize initial or final vector

  if(i==0)then
    print*,' oops bad index (W)  ',i,iproc
    stop
  end if
  if( (vchar == 'n' .and. ipoint == 'i' ) .or. &
        (vchar == 'r' .and. ipoint == 'f') )then
          v => vec1
          vstart = v1s
          vstop  = v1e
  else
          v => vec2
          vstart = v2s
          vstop  = v2e
  end if

!................... INTERNAL STORAGE.............


  if( storelanczosincore1 .or. storelanczosincoreMPI)then
        do jl = vstart,vstop
             lvec(jl,i) = v(jl+Lstart-1)      ! Lstart is set in routine distribute_lanczos_pieces
        end do
        if(.not.forcewritelanczos .or. storelanczosincoreMPI)return
  end if

!.................... WRITE TO DISK....................
  if ( iproc == 0 ) then
        if ( i == 1 ) rewind(iunit)
        write( iunit ) ( v(jl), jl = vstart,vstop )
  endif

  return
end subroutine write_lanczos_vector_a

!================================================================
! READ LANCZOS VECTOR_A  (A = all parallelisms)
! 
! for new parallel scheme
!
!================================================================
subroutine read_lanczos_vector_a(vchar,ipoint,i,iunit)
  use lanczos_info
  use precisions
  use basis
  use nodeinfo
  use localvectors
  use fragments
  use flagger
  use bmpi_mod
  implicit none

  character(1) :: vchar,ipoint
  integer(4) :: i,iunit

  real(kind=lanc_prec),pointer   :: v(:)
  integer (kind=basis_prec) :: jl,vstart,vstop

!........................................................................... 
#ifdef _MPI	  
type(MPI_status) :: status
!  integer(4)                    :: status(MPI_STATUS_SIZE)
#endif  
  integer ierr

!....................................................................

	if(useNewReorthog) then
		print *, "read_lanczos_vector_a not supported yet for new reorthog"
		stop 1
	end if

!..... choose which vector to normalize
!      complicated because during lanczos 
!      we switch between vec1 and vec2 (stored in module localvectors
!      also can normalize initial or final vector

  if(i==0)then
    print*,' oops bad index (R)  ',i,iproc
    stop
  end if

  if( (vchar == 'n' .and. ipoint == 'i' ) .or. & 
        (vchar == 'r' .and. ipoint == 'f') )then
          v => vec1
          vstart = v1s
          vstop  = v1e
  else
          v => vec2
          vstart = v2s
          vstop  = v2e
  end if


!................... INTERNAL STORAGE.............
  if(storelanczosincoreMPI)then
          vstart = 1
          vstop = Ldim      ! Ldim is set in routine distribute_lanczos_pieces
  end if

  if( storelanczosincore1 .or. storelanczosincoreMPI )then
           do jl =  vstart,vstop
             v(jl +Lstart- 1) = Lvec(jl,i)  ! Lstart is set in routine distribute_lanczos_pieces
           end do
           return
  end if

  if ( iproc == 0 ) then
        if ( i == 1 ) rewind(iunit)
        read( iunit,end=1 ) (v(jl),jl = vstart,vstop )  

  endif

  return
1 continue
  if(iproc==0)then
     write(6,*)' Error in reading lanczos vector ',i
     write(6,*)vchar,ipoint
  end if
  stop

end subroutine read_lanczos_vector_a

!================================================================
!
! read in lanczos vectors when restarting; place into core if needed
!
subroutine read_lanczos_vector_restart_a(vchar,ipoint,i,iunit)
  use lanczos_info
  use precisions
  use basis
  use nodeinfo
  use localvectors
  use fragments
  use flagger
  use bmpi_mod
  implicit none

  character(1) :: vchar,ipoint
  integer(4) :: i,iunit

  real(kind=lanc_prec),pointer   :: v(:)
  integer (kind=basis_prec) :: jl,vstart,vstop

!........................................................................... 
#ifdef _MPI	  
type(MPI_status) :: status
!  integer(4)                    :: status(MPI_STATUS_SIZE)
#endif
  integer ierr

!....................................................................
	if(useNewReorthog) then
		print *, "read_lanczos_vector_restart_a not supported yet with new reorthog"
		stop 1
	end if

!  if(storelanczosincore)then
!     call lanczos_piece_extract(vchar,ipoint,i,.false.,1.0d0)
!     return
!  end if

!..... choose which vector to normalize
!      complicated because during lanczos 
!      we switch between vec1 and vec2 (stored in module localvectors
!      also can normalize initial or final vector

  if( (vchar == 'n' .and. ipoint == 'i' ) .or. & 
        (vchar == 'r' .and. ipoint == 'f') )then
          v => vec1
          vstart = v1s
          vstop  = v1e
  else
          v => vec2
          vstart = v2s
          vstop  = v2e
  end if
  if ( iproc == 0 ) then
        if ( i == 1 ) rewind(iunit)
        read( iunit,end=1 ) (v(jl),jl = vstart,vstop )  
  endif
  if( storelanczosincore1 )then
           do jl =  vstart,vstop
             lvec(jl,i) = v(jl)
           end do
           return
  end if


  return
1 continue
  if(iproc==0)then
     write(6,*)' Error in reading lanczos vector ',i
     write(6,*)vchar,ipoint
  end if
  stop

end subroutine read_lanczos_vector_restart_a

!=================================================================
!
!  CALLS
!     lanczos_seq_orthog_a
!     lanczos_simul_orthog_a
!     read_lanczos_vector_a
!     clocker
!     dvecproj_p
!
subroutine reorthogonalize_a(iter,vchar,dot)

  use nodeinfo
  use localvectors
  use precisions
  use verbosity
  use lanczos_info
  use basis
  use flagger
  use fragments
  implicit none

  integer(4):: iter
  character(1) :: vchar
  integer(4) :: j
  real(8) :: dot,ddot
  integer ierr



  if(.not.reorthog)return
  if(iproc==0)call clocker('ort','sta')

  if(useNewReorthog) then
  	print *, "reorthogonalize_a:  use breorthog.f90"
	stop 1
  end if

  if( (storelanczosincoreMPI .or. storelanczosincore1))then  
     ! REORTHGONALIZE SIMULTANEOUSLY USING LANCZOS STORED IN CORE
     if(iproc==0) call clocker('ror','sta')

     if( orthog_sequential)then
     call lanczos_seq_orthog_a(vchar,iter,dot)

     else
     call lanczos_simul_orthog_a(vchar,iter,dot)
     endif
     if(iproc==0) call clocker('ror','end')
     if(iproc==0)  call clocker('ort','end')
     return
  end if

  if ( verbose_orthog .and. iproc == 0 ) write(6,*) ' Starting reorthog ',iter

  if( iproc ==0)then 
     do j = 1, iter
           if(iproc ==0)call clocker('ror','sta')
           call read_lanczos_vector_a(vchar,'f',j,lvec_file)   ! WILL READ EITHER FROM FILE OR RECREATE FROM MEMORY
           if(iproc==0)call clocker('ror','end')
           call dvecproj_p(vchar,dot)
     end do
  end if
  if(iproc==0)call clocker('ort','end')
  if ( verbose_orthog .and. iproc == 0) write(6,*) ' Finished reorthog '

  return

end subroutine reorthogonalize_a

!=======================================================
!
! orthogonalizes lanczos vector iter+1 
! (stored either in vec1 or vec2, labeled by vchar)
! against the previous iter vectors
!
!  HOWEVER: Need to orthogonalize sequentially against vectors with
!  substantial overlaps (e.g., iter and iter-1)
!
  subroutine lanczos_simul_orthog_a(vchar,iter,dalpha)
  use localvectors
  use nodeinfo
  use lanczos_info
  use basis
  use precisions
  use flagger
  use fragments
  use bmpi_mod
  implicit none

  character(1) :: vchar
  integer(4) :: ierr
  integer iter,ivec
  integer i
  integer (kind=8):: vshift,vstart,vend,vlength
  real(kind=8) :: dv1,dv2,dalpha, dovlp

  integer(kind=8) :: istate, istate1,istate2
 
  real(kind=8) :: pover(niter),sumover(niter),otmp

  real(kind=lanc_prec), pointer :: v(:)
  real(kind=8):: vtmp
  integer myproc  ! for testing only

  ! KSM - supposed to use breorthog.f90
  if(useNewReorthog) then
	print *, "lanczos_simul_orthog_a: should not call with new reorthog"
	stop 1
  end if

  if(iter==0)then
    print*,' oops bad index (O)  ',iter,iproc
    stop
  end if

  if(vchar=='n')then
     v => vec1
  else
     v => vec2
  end if


if(.not. alpha_before_orthog)then
!--------- ORTHOGONALIZE AGAINST LAST EIGENVECTOR

       otmp = 0.0d0
!$omp parallel do private(istate, dv1, dv2), shared(ivec,lstart,ldim)reduction(+ : otmp)
     do istate = 1, Ldim
        dv1 = real(Lvec(istate,iter), kind=8)
        dv2 = real(v(Lstart + istate - 1), kind=8)
        otmp = otmp + dv1*dv2
     end do
!$omp end parallel do
if(nproc > 1) then
#ifdef _MPI			  	
    call BMPI_ALLREDUCE(otmp, dalpha, 1,  MPI_SUM, MPI_COMM_WORLD, ierr) 
#endif	
else
     dalpha = otmp
end if
istate1 = Lstart - 1
!$omp parallel do private(istate,ivec,vtmp), shared(istate1, sumover, Lvec)
do istate = 1, Ldim
  vtmp = real(v(istate1+istate),kind=8)
      vtmp = vtmp - dalpha*real(Lvec(istate,iter),kind=8)
  v(istate1+istate)=real(vtmp,kind=lanc_prec)
end do
!$omp end parallel do

end if
!--------- ORTHOGONALIZE AGAINST NEXT-TO-LAST EIGENVECTOR

if(iter==1)return
       otmp = 0.0d0
!$omp parallel do private(istate, dv1, dv2), shared(ivec,lstart,ldim)reduction(+ : otmp)
     do istate = 1, Ldim
        dv1 = real(Lvec(istate,iter-1), kind=8)
        dv2 = real(v(Lstart + istate - 1), kind=8)
        otmp = otmp + dv1*dv2
     end do
!$omp end parallel do
if(nproc > 1) then
#ifdef _MPI			  	
    call BMPI_ALLREDUCE(otmp, dovlp, 1, MPI_SUM, MPI_COMM_WORLD, ierr) 
#endif	
else
     dovlp = otmp
end if
istate1 = Lstart - 1
!$omp parallel do private(istate,ivec,vtmp), shared(istate1, sumover, Lvec)
do istate = 1, Ldim
  vtmp = real(v(istate1+istate),kind=8)
      vtmp = vtmp - dovlp*real(Lvec(istate,iter-1),kind=8)
  v(istate1+istate)=real(vtmp,kind=lanc_prec)
end do
!$omp end parallel do

!.......... NOW ORTHOGONALIZE AGAINST THE REST

  do ivec = 1,iter
  pover(ivec) = 0.d0
  sumover(ivec) = 0.d0
  end do
  if(storelanczosincore1)then
     Ldim = dimbasis
  end if

  do ivec = 1, iter
     otmp = 0.0d0
!$omp parallel do private(istate, dv1, dv2), shared(ivec,lstart,ldim)reduction(+ : otmp)
     do istate = 1, Ldim
        dv1 = real(Lvec(istate,ivec), kind=8)
        dv2 = real(v(Lstart + istate - 1), kind=8)
        otmp = otmp + dv1*dv2
     end do
!$omp end parallel do
     pover(ivec) = otmp
  end do
 ! print*,pover(1:iter)

!----- REDUCE------ !
!  NOTE: right now I do this clumsily--I have all nodes reduce even though some may not contribute
!  This can be finessed later
!
  if(nproc > 1)then
#ifdef _MPI			  	  
    call BMPI_ALLREDUCE(pover, sumover, niter,  MPI_SUM, MPI_COMM_WORLD, ierr)   
#endif	
  else
    do ivec=1,iter
       sumover(ivec)=pover(ivec)
    end do
  end if
!----  Zero the vector where it doesn't exist on this processor
  do istate = 1, dimbasis
     if(istate < Lstart .or. istate >= Lstart + Ldim) v(istate) = 0.0
  end do


!----  Subtract the overlap of the previous Lanczos vectors
! NOTE: switched order of loops so can do subtraction in double precision -- CWJ @LLNL, 6/2013

istate1 = Lstart - 1
!$omp parallel do private(istate,ivec,vtmp), shared(istate1, sumover, Lvec)
do istate = 1, Ldim
  vtmp = real(v(istate1+istate),kind=8)
  do ivec = 1, iter
      vtmp = vtmp - sumover(ivec)*real(Lvec(istate,ivec),kind=8)
  end do
  v(istate1+istate)=real(vtmp,kind=lanc_prec)
end do
!$omp end parallel do

!  dalpha = sumover(iter)

  if(nproc < 2)return

#ifdef _MPI			  	 
  if( vchar == 'n' )then
     ! call block_reduce(dimbasis,vec1)
     call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, MPI_COMM_WORLD, ierr)
  else
     ! call block_reduce(dimbasis,vec2)
     call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, MPI_COMM_WORLD, ierr)
  end if
#endif	 

  return
  end subroutine lanczos_simul_orthog_a

!=======================================================
!
! orthogonalizes lanczos vector iter+1 
! (stored either in vec1 or vec2, labeled by vchar)
! against the previous iter vectors
!
  subroutine lanczos_seq_orthog_a(vchar,iter,dalpha)
  use localvectors
  use nodeinfo
  use lanczos_info
  use basis
  use precisions
  use flagger
  use fragments
  use bmpi_mod
  implicit none

  character(1) :: vchar
  integer(4) :: ierr
  integer iter,ivec
  integer i
  integer (kind=8):: vshift,vstart,vend,vlength
  real(kind=8) :: dv1,dv2,dalpha

  integer(kind=8) :: istate, istate1,istate2
 
  real(kind=8) :: pover,sumover,otmp

  real(kind=lanc_prec), pointer :: v(:)
  real(kind=8):: vtmp
  integer myproc  ! for testing only

  ! KSM - supposed to use breorthog.f90
  if(useNewReorthog) then
      print *, "lanczos_seq_orthog_a: should not call with new reorthog"
      stop 1
  end if

  if(vchar=='n')then
     v => vec1
  else
     v => vec2
  end if

  pover = 0.d0
  sumover = 0.d0

  if(storelanczosincore1)then
     Ldim = dimbasis
  end if

istate1 = Lstart - 1

  do ivec = 1, iter

     otmp = 0.0d0
!$omp parallel do private(istate, istate1, dv1, dv2), reduction(+ : otmp)
     do istate = 1, Ldim
        dv1 = real(Lvec(istate,ivec), kind=8)
        dv2 = real(v(Lstart + istate - 1), kind=8)
        otmp = otmp + dv1*dv2
     end do
!$omp end parallel do
          
!----- REDUCE------ !
!  NOTE: right now I do this clumsily--I have all nodes reduce even though some may not contribute
!  This can be finessed later
!
  if(nproc > 1)then 
#ifdef _MPI	  
    call BMPI_ALLREDUCE(otmp, sumover, 1,  MPI_SUM, MPI_COMM_WORLD, ierr) 
#endif
  else
       sumover=otmp
  end if

! KSM:  made ivec shared.  It was private, which leads to it not being initialized
!$omp parallel do private(istate,vtmp), shared(istate1, sumover, Lvec, ivec)
do istate = 1, Ldim
  vtmp =   real(v(istate1+istate),kind=8)
  vtmp = vtmp -sumover*real(Lvec(istate,ivec),kind=8)
  v(istate1+istate)= real(vtmp,kind=lanc_prec)
end do
!$omp end parallel do

  end do

  dalpha = sumover  ! after last iteration

  if(nproc < 2)return

!----  Zero the vector where it doesn't exist on this processor
  do istate = 1, dimbasis
     if(istate < Lstart .or. istate >= Lstart + Ldim) v(istate) = 0.0
  end do

#ifdef _MPI	  

  if( vchar == 'n' )then
     ! call block_reduce(dimbasis,vec1)
     call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, MPI_COMM_WORLD, ierr)
  else
     ! call block_reduce(dimbasis,vec2)
     call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, MPI_COMM_WORLD, ierr)
  end if
#endif
  return
  end subroutine lanczos_seq_orthog_a
  
!

!================================================
subroutine initialize_final(vchar)
   use nodeinfo
   use fragments
   use localvectors
   use precisions
   implicit none
   character(1) :: vchar

   integer(kind=basis_prec) :: i
   integer :: tid
   real(kind=lanc_prec), pointer :: v2p(:)
   ! OpenMP functions
   integer(kind=4) :: omp_get_thread_num, omp_get_num_threads

   if(useVec2Thread) then
!$omp parallel private(i, tid, v2p)  &
!$omp       shared(v2s, v2e)         &
!$omp       shared(vec2threadchunk, vec2threadchunkm1, vec2threadflat)
      tid = omp_get_thread_num();
      i = tid * vec2threadchunk
      v2p(v2s:v2e) => vec2threadflat(i: i + vec2threadchunkm1) 
      v2p(:) = 0.0
!$omp end parallel
      ! fall through and initialize vec2 as well
      ! needed anyway for final reduce
   end if

   select case(vchar)
      case('n')

      do i = basestart(frag2), basestop(frag2)
         vec2(i) = 0.0e0_lanc_prec
      end do

      case('r')

      do i = basestart(frag1), basestop(frag1)
         vec1(i) = 0.0e0_lanc_prec
      end do
   end select

   return
end subroutine initialize_final
!  

end module bvectorlib_mod
