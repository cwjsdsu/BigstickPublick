module switch_mod
contains

!
!==============================================================

!  subroutine sort_bundled1bops
!  by CWJ @ SDSU 10/2011
!
!  based on subroutine sort_blocksPN_states
!  by WEO @ LLNL - 9/2011
!
!  sorts jumps for OpenMP application of PN
!  renamed from sort_bundle1bops in 7.7.9
!
!  CALLED BY:
!       main routine
!
!  CALLS: subroutine sort1bjumps
!
subroutine sort_bundled_jumps
   use sectors
   use jumpNbody
   use jump3body
   use precisions
   use basis
   use nodeinfo
!   use distrinfo
   use verbosity
   use opbundles
   use nodeinfo
   use flagger
   use bmpi_mod
   implicit none

   integer(8) :: j

   integer(8) :: npjmps, nnjmps, startp, startn
   integer    :: ibundle 
   integer    :: iprocs, procstart,procstop

!........ OPTION TO SIMULATE MPI ON A SINGLE "CORE".....

  if((simulateMPI .and. nproc ==1) .or. .not.restrictjumps)then
      procstart = 0
      procstop  = nprocs -1
  else
      procstart = iproc
      procstop  = iproc
  end if

  do iprocs = procstart,procstop
     do ibundle = opbundlestart(iprocs), opbundleend(iprocs)

		 select case (opbundle(ibundle)%optype )
 !     if(opbundle(ibundle)%optype == 'PN')then
 !............................................................................. PN
            case ('PN')
               startp = opbundle(ibundle)%pxstart-1
               startn = opbundle(ibundle)%nsortstart-1
               npjmps = opbundle(ibundle)%pxend -startp 
               nnjmps = opbundle(ibundle)%nsortend -startn 
!--------------------   Sort proton jump states on final SD, p1b_fsd
!--------------------   Order them to enable OpenMP for the forward application of H_pn
               if(dosortjumps)call sort1bjumps(npjmps,p1b_isd(startp+1),p1b_fsd(startp+1),              &
                       p1b_cop(startp+1),p1b_dop(startp+1),                     &
                       p1b_phase(startp+1))

!--------------------   Sort neutron jump states on INITIAL SD, n1b_isd
!--------------------   Order them to enable OpenMP for the backward application of H_pn;
!
! NOTE BY CWJ: However, do we need to sort these every time?
! 
               if(dosortjumps)call sort1bjumps(nnjmps,n1b_fsd(startn+1),n1b_isd(startn+1),              &
                       n1b_cop(startn+1),n1b_dop(startn+1),                     &
                       n1b_phase(startn+1))   
!              end if
!............................................................................. PPN
            case('PPN')
!               if(opbundle(ibundle)%optype == 'PPN')then
                   startp = opbundle(ibundle)%pxstart-1
                startn = opbundle(ibundle)%nsortstart-1
                npjmps = opbundle(ibundle)%pxend -startp 
                nnjmps = opbundle(ibundle)%nsortend -startn 
!--------------------   Sort proton jump states on final SD, p1b_fsd
!--------------------   Order them to enable OpenMP for the forward application of H_pn
                if(dosortjumps)call sort2bjumps(npjmps,p2b_isd(startp+1),p2b_fsd(startp+1),              &
                       p2b_cop(startp+1),p2b_dop(startp+1),                     &
                       p2b_phase(startp+1))

!--------------------   Sort neutron jump states on final SD, n1b_isd
!--------------------   Order them to enable OpenMP for the backward application of H_pn;
!
! NOTE BY CWJ: However, do I need to sort these every time?
! 
                if(dosortjumps)call sort1bjumps(nnjmps,n1b_fsd(startn+1),n1b_isd(startn+1),              &
                       n1b_cop(startn+1),n1b_dop(startn+1),                     &
                       n1b_phase(startn+1))   
!      end if
!............................................................................. PNN
            case('PNN')
!      if(opbundle(ibundle)%optype == 'PNN')then
                startp = opbundle(ibundle)%pxstart-1
                startn = opbundle(ibundle)%nsortstart-1
                npjmps = opbundle(ibundle)%pxend -startp 
                nnjmps = opbundle(ibundle)%nsortend -startn 
!--------------------   Sort proton jump states on final SD, p1b_fsd
!--------------------   Order them to enable OpenMP for the forward application of H_pn
                if(dosortjumps)call sort1bjumps(npjmps,p1b_isd(startp+1),p1b_fsd(startp+1),              &
                       p1b_cop(startp+1),p1b_dop(startp+1),                     &
                       p1b_phase(startp+1))

!--------------------   Sort neutron jump states on final SD, n1b_isd
!--------------------   Order them to enable OpenMP for the backward application of H_pn;
!
! NOTE BY CWJ: However, do I need to sort these every time?
! 

                if(dosortjumps)call sort2bjumps(nnjmps,n2b_fsd(startn+1),n2b_isd(startn+1),              &
                       n2b_cop(startn+1),n2b_dop(startn+1),                     &
                       n2b_phase(startn+1))   

!      end if
!............................................................................. NN; added in 7.6.1
            case('NX')            ! I'm turning this off with a dummy optype, doesn't speed up
                startn = opbundle(ibundle)%nsortstart-1
                nnjmps = opbundle(ibundle)%nsortend -startn 

!--------------------   Sort neutron jump states on final SD, n1b_isd
!--------------------   Order them to enable OpenMP for the backward application of H_pn;
!
! NOTE BY CWJ: However, do I need to sort these every time?
! 

                if(dosortjumps)call sort2bjumps1op(nnjmps,n2b_fsd(startn+1),n2b_isd(startn+1),              &
                       n2b_op(startn+1), n2b_phase(startn+1))   

!............................................................................. PP; added in 7.6.1
			case('PX')  ! I'm turning this off with a dummy optype, doesn't speed up
				startp = opbundle(ibundle)%pxstart-1
			    npjmps = opbundle(ibundle)%pxend -startp 
!--------------------   Sort proton jump states on final SD, p1b_fsd
!--------------------   Order them to enable OpenMP for the forward application of H_pn
			    if(dosortjumps)call sort2bjumps1op(npjmps,p2b_isd(startp+1),p2b_fsd(startp+1),              &
					                          p2b_op(startp+1),  p2b_phase(startp+1))
					 
			case default    ! do nothing

        end select
   end do  ! ibundle
 end do  ! iprocs

end subroutine sort_bundled_jumps
   	
!
!
! Use shell sort to sort jumps.  The shell sort has
! the nice property that it is easy to thread
! not quite as fast for large arrays as qsort, but bulletproof.
! Previous qsort attempt was going out of bounds in rare circumstances.
!
! sorted on "x1b_fsd"
!
! CALLED BY: sort_bundled1bops
!
! CALLS: sort1bjumps_ins
!
subroutine sort1bjumps(ndim,x1b_isd,x1b_fsd,x1b_cop,x1b_dop,x1b_phase)
   use nodeinfo
   use precisions
   implicit none
   ! arguments
   integer(kind=8), intent(in) :: ndim
   integer(kind=basis_prec), intent(inout) :: x1b_isd(ndim), x1b_fsd(ndim)
   integer(kind=2), intent(inout) :: x1b_cop(ndim), x1b_dop(ndim)
   integer(kind=1), intent(inout) :: x1b_phase(ndim)
   ! locals
   integer(kind=basis_prec) :: off, h, i, k
   logical :: sorted
   ! For reasonable sized arrays a shell sort is  quite fast and is very simple
   ! We can also easily implement OMP on the inner loop
   ! shell sort starts with big strides and finishes with swaps
   ! should generate more strides or do dynamically.  See wikipedia article on shell sort
   ! I extended Sedgwick's table with Floor[2.25 a_{n-1}] + 1
   integer(kind=8) :: strides(28) = (/ &
      23422578909_8, 10410035071_8, 4626682254_8, 2056303224_8, &
      913912544_8, 406183353_8, 180525935_8, 80233749_8, &
      35659463_8, 15848650_8, 7043844_8, 3130597_8, &
      1391376_8, 463792_8, 198768_8, 86961_8, &
      33936_8, 13776_8, 4592_8, 1968_8, &
      861_8, 336_8, 112_8, 48_8, &
      21_8, 7_8, 3_8, 1_8/)
   integer(kind=8) :: stride

   do k=1, SIZE(strides,1)
      stride = strides(k)
      if(stride > ndim) cycle
      ! only off should be private
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(off)
!$OMP DO SCHEDULE(STATIC)
      do off=0, stride-1
         call sort1bjumps_ins(off,stride,ndim,x1b_isd,x1b_fsd,x1b_cop,x1b_dop,x1b_phase)
      end do
!$OMP END DO
!$OMP END PARALLEL
   end do

   ! check that sort worked
   sorted = .true.
   do i = 2, ndim
      if(x1b_fsd(i-1) > x1b_fsd(i))then
         sorted=.false.
         exit
      end if
   end do
   if(.not. sorted) then
      print *, "iproc=", iproc, ", sort1bjumps failed"
      stop 1 ! has msg
   end if
end subroutine sort1bjumps
!
!
! Perform an insertion sort with a stride and an offset
! Initial strides will be large
!
! CALLED BY: sort1bjumps
!
subroutine sort1bjumps_ins(off,stride, ndim,x1b_isd,x1b_fsd,x1b_cop,x1b_dop,x1b_phase)
   use nodeinfo
   use precisions
   implicit none
   ! arguments
   integer(kind=8), intent(in) :: off, stride, ndim
   integer(kind=basis_prec), intent(inout) :: x1b_isd(ndim), x1b_fsd(ndim)
   integer(kind=2), intent(inout) :: x1b_cop(ndim), x1b_dop(ndim)
   integer(kind=1), intent(inout) :: x1b_phase(ndim)
   ! locals
   integer(kind=basis_prec) :: s_isd, s_fsd
   integer(kind=2) :: s_cop, s_dop
   integer(kind=1) :: s_phase
   integer(kind=8) :: i, j, k

   ! scan array with lookback of stride for comparison
   do j = off+stride+1, ndim, stride
      s_isd = x1b_isd(j)
      s_fsd = x1b_fsd(j) ! sorting key
      s_cop = x1b_cop(j)
      s_dop = x1b_dop(j)
      s_phase = x1b_phase(j)
      k = j
      i = k - stride
      do while(i > 0) 
         if(x1b_fsd(i) <= s_fsd) exit
         x1b_isd(k) = x1b_isd(i)
         x1b_fsd(k) = x1b_fsd(i)
         x1b_cop(k) = x1b_cop(i)
         x1b_dop(k) = x1b_dop(i)
         x1b_phase(k) = x1b_phase(i)
         k = i
         i = i - stride
      end do
      x1b_isd(k) = s_isd
      x1b_fsd(k) = s_fsd
      x1b_cop(k) = s_cop
      x1b_dop(k) = s_dop
      x1b_phase(k) = s_phase
   end do
end subroutine sort1bjumps_ins

!==============================================================
!  subroutine switch
!  by WEO @ LLNL - 9/2011
!
!  SUBROUTINES CALLED:
!   switchi4
!   switchi2
!   switchi1
!
subroutine switch(idim,i1,i2,iarr1,iarr2,iarr3,iarr4,iarr5)
   use precisions
   implicit none

   integer(kind=8) :: idim, i1, i2
   integer(kind=basis_prec) :: iarr1(idim),iarr2(idim)
   integer(kind=2) :: iarr3(idim),iarr4(idim)
   integer(kind=1) :: iarr5(idim)

   call switchibp(idim,i1,i2,iarr1)
   call switchibp(idim,i1,i2,iarr2)
   call switchi2(idim,i1,i2,iarr3)
   call switchi2(idim,i1,i2,iarr4)
   call switchi1(idim,i1,i2,iarr5)

   return

end subroutine switch   

!==============================================================
!  subroutine switch1
!  by WEO @ LLNL - 9/2011
!
!  SUBROUTINES CALLED:
!  switchi4
!  switchi1

subroutine switch2(idim,i1,i2,iarr1,iarr2,iarr3,iarr5)
   implicit none
   integer(kind=8) :: idim, i1, i2
   integer(kind=4) :: iarr1(idim),iarr2(idim)
   integer(kind=4) :: iarr3(idim)
!   real(kind=4)    :: xarr4(idim)
   integer(kind=1) :: iarr5(idim)

   call switchi4(idim,i1,i2,iarr1)
   call switchi4(idim,i1,i2,iarr2)
   call switchi4(idim,i1,i2,iarr3)
!   call switchr4(idim,i1,i2,xarr4)
   call switchi1(idim,i1,i2,iarr5)

   return

end subroutine switch2  

subroutine switchi4(idim,i1,i2,iarray)
   use precisions
   implicit none
   integer(kind=8) :: idim, i1, i2
   integer(kind=4) :: iarray(idim)
   integer(kind=4) :: itemp

   itemp = iarray(i2)
   iarray(i2) = iarray(i1)
   iarray(i1) = itemp
   return
end subroutine switchi4

subroutine switchibp(idim,i1,i2,iarray)
   use precisions
   implicit none
   integer(kind=8) :: idim, i1, i2
   integer(kind=basis_prec) :: iarray(idim)
   integer(kind=8) :: itemp

   itemp = iarray(i2)
   iarray(i2) = iarray(i1)
   iarray(i1) = itemp
   return
end subroutine switchibp

subroutine switchr4(idim,i1,i2,xarray)
   implicit none
   integer(kind=8) :: idim, i1, i2
   real(kind=4)    :: xarray(idim)
   real(kind=4)    :: xtemp

   xtemp = xarray(i2)
   xarray(i2) = xarray(i1)
   xarray(i1) = xtemp
   return
end subroutine switchr4

subroutine switchi2(idim,i1,i2,iarray)
   implicit none
   integer(kind=8) :: idim, i1, i2
   integer(kind=2) :: iarray(idim)
   integer(kind=2) :: itemp

   itemp = iarray(i2)
   iarray(i2) = iarray(i1)
   iarray(i1) = itemp
   return
end subroutine switchi2

subroutine switchi1(idim,i1,i2,iarray)
   implicit none
   integer(kind=8) :: idim, i1, i2
   integer(kind=1) :: iarray(idim)
   integer(kind=1) :: itemp

   itemp = iarray(i2)
   iarray(i2) = iarray(i1)
   iarray(i1) = itemp
   return
end subroutine switchi1


!==============================================================
!   set_bundledXY_threadstart
!  modified by CWJ @ SDSU - 11/2011
!
!  based on subroutine set_blocksPN_start
!  by WEO @ LLNL - 9/2011
!
!  sets start (and stop) for OpenMP threads for PN application
!
! CALLED BY: main routine
!
subroutine set_bundledXY_threadstart

   use sectors
   use jumpNbody
   use precisions
   use basis
   use nodeinfo
 ! use distrinfo
   use verbosity
   use opbundles
   use nodeinfo
   use flagger
   use jumplimits
   use bmpi_mod
   implicit none

   integer(kind=8) :: i, j
   integer :: k  ! for threads
   integer(kind=8) npjmps, nnjmps, startp, startn
   integer(kind=8) :: npjmpsf, npjmpsb
   integer(kind=8) :: startp1
   
   integer(8) :: num_sd
   integer(8) ::  min_num, max_num
   integer(8) :: iavg
   integer    :: ithread
   integer(8) :: sumj
   real(kind=8) :: avg, sdavg, sdvar
   integer(8) :: tot_jmps
   integer ibundle
   real eff, eff_avg, tot
   
   integer :: num_threads

   integer :: omp_get_num_threads
   integer(basis_prec),pointer :: pXb_fsd(:), nYb_isd(:)

   integer :: procstart,procstop,iprocs
   integer :: aerr
   if(iproc==0)print*,' SORTING THREADS '
   
!$omp parallel shared(num_threads)  
   num_threads = omp_get_num_threads()
!$omp end parallel

   eff_avg = 0.0
   tot_jmps = 0

!........ OPTION TO SIMULATE MPI ON A SINGLE "CORE".....
  if((simulateMPI .and. nproc ==1) .or. .not. restrictjumps)then
      procstart = 0
      procstop  = nprocs -1
  else
      procstart = iproc
      procstop  = iproc
  end if

  do iprocs = procstart,procstop
     do ibundle = opbundlestart(iprocs), opbundleend(iprocs)

      if(opbundle(ibundle)%optype == 'PN' .or. & 
          opbundle(ibundle)%optype == 'PPN' & 
         .or. opbundle(ibundle)%optype == 'PNN')then

          select case( opbundle(ibundle)%optype )

             case ('PN')
!			 print*,' sorting ',iproc, ibundle,opbundle(ibundle)%pxstart-1,p1bjumplength(iproc)

                pXb_fsd => p1b_fsd
                nYb_isd => n1b_isd

             case ('PPN')
                pXb_fsd => p2b_fsd
                nYb_isd => n1b_isd

             case ('PNN')
                pXb_fsd => p1b_fsd
                nYb_isd => n2b_isd
             case default
                print *, "bad bun type ", opbundle(ibundle)%optype
                call printstacktrace
                stop 1 ! has msg
          end select

          startp = opbundle(ibundle)%pxstart-1
          startn = opbundle(ibundle)%nxstart-1
          npjmps = opbundle(ibundle)%pxend -startp 
          nnjmps = opbundle(ibundle)%nxend -startn 
          if(.not. allocated(opbundle(ibundle)%startp_thread)) then 
             allocate(opbundle(ibundle)%startp_thread(0:num_threads), stat=aerr)
             if(aerr /= 0) call memerror("set_bundledXY_threadstart 1")
          end if
          if(.not. allocated(opbundle(ibundle)%startn_thread)) then
             allocate(opbundle(ibundle)%startn_thread(0:num_threads), stat=aerr)
             if(aerr /= 0) call memerror("set_bundledXY_threadstart 2")
          end if

!--------   Set up blocks for protons - based on initial p1b_fsd
!--------   Used for backward application of H_pn

!---  Get data for averages, average num/thread, etc
          num_sd = 1
          sumj = 0
          sdavg = 0.0
          sdvar = 0.0
          do j = startp + 2, startp + npjmps
               if(pXb_fsd(j) == pXb_fsd(j-1))then
                   sumj = sumj + 1
               else
                   sdavg = sdavg + real(sumj)
                   sdvar = sdvar + real(sumj)**2
                   sumj = 0
                   num_sd = num_sd + 1
               end if      
           end do

           sdavg = sdavg/real(num_sd)
           sdvar = sdvar/real(num_sd)
           sdvar = sqrt(sdvar - sdavg**2)
           avg = real(npjmps)/real(num_threads)
           iavg = int(avg, 8)
           opbundle(ibundle)%startp_thread(0:num_threads-1) = 0
           ithread = 0
           sumj = 0
           opbundle(ibundle)%startp_thread(ithread) = startp

!  ----   Do appropriate break up for forward direction

           do j = startp + 2, startp + npjmps
               if(pXb_fsd(j) == pXb_fsd(j-1))then
                   sumj = sumj + 1
               else
                   if(ithread == num_threads - 1) then
                      sumj = sumj +1
                   else
                       if(sumj + sdavg/2 < iavg)then   ! See if current sum is too close to the average
               	            sumj = sumj + 1
                       else                    ! if exceeds, so let's stop
                           ithread = ithread + 1
                           sumj = 0
                           opbundle(ibundle)%startp_thread(ithread) = j - 1
                       end if
                   end if
               end if      
            end do  ! j

            do k = ithread+1, num_threads
                opbundle(ibundle)%startp_thread(k) = startp + npjmps
            end do   ! k

!------ Set up blocks for neutrons - based on initial n1b_isd
!------  Used for backward application of H_pn

             num_sd = 1
             sumj = 0
             sdavg = 0.0
             sdvar = 0.0
             do j = startn + 2, startn + nnjmps
                  if(nYb_isd(j) == nYb_isd(j-1))then
                    sumj = sumj + 1
                  else
                    sdavg = sdavg + real(sumj)
                    sdvar = sdvar + real(sumj)**2
                    sumj = 0
                    num_sd = num_sd + 1
                  end if      
              end do  ! j

              sdavg = sdavg/real(num_sd)
              sdvar = sdvar/real(num_sd)
              sdvar = sqrt(sdvar - sdavg**2)
              avg = real(nnjmps)/real(num_threads)
              iavg = int(avg, 8)
              opbundle(ibundle)%startn_thread(0:num_threads-1) = 0
              ithread = 0
              sumj = 0
              opbundle(ibundle)%startn_thread(ithread) = startn
               do j = startn + 2, startn + nnjmps

                   if(nYb_isd(j) == nYb_isd(j-1))then
                       sumj = sumj + 1
                   else
                       if(ithread == num_threads - 1) then
                           sumj = sumj + 1
                       else
                           if(sumj + sdavg/2 < iavg)then   ! See if current sum is too close to the average
               	              sumj = sumj + 1
                            else             ! if exceeds, so let's stop
                              ithread = ithread + 1
                              sumj = 0
                              opbundle(ibundle)%startn_thread(ithread) = j - 1
                            end if
                       end if
                    end if      
              end do ! j

              do k = ithread+1, num_threads
                   opbundle(ibundle)%startn_thread(k) = startn + nnjmps
              end do  !k
              k = 0
              
              do j = 1,num_threads
                  if( opbundle(ibundle)%startp_thread(j) > opbundle(ibundle)%startp_thread(j-1) )k=k+1
              end do
              opbundle(ibundle)%numpthreads = k
              k = 0
              
              do j = 1,num_threads
                  if( opbundle(ibundle)%startn_thread(j) > opbundle(ibundle)%startn_thread(j-1) )k=k+1
              end do
              opbundle(ibundle)%numnthreads = k

      tot_jmps = tot_jmps + npjmps

      end if

      
   end do  ! ibundle
  end do   ! iprocs

!   eff_avg = eff_avg/real(tot_jmps)
!   write(6,*)'Average efficiency =',eff_avg
   
end subroutine set_bundledXY_threadstart

!==============================================================

!  subroutine sort2bjumps
!  by WEO @ LLNL - 9/2011
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Try quick-sort algorithm isted of straight insertion        +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    ALSO needed to introduce "switch2b" to get this right
!
!  SUBROUTINES CALLED
!     switch2b

subroutine sort2bjumps(ndim, x1b_isd, x1b_fsd, x1b_cop,x1b_dop,x1b_phase)
   use bmpi_mod
   use precisions
   implicit none
!---   Data passed in   
   integer(kind=8) :: ndim
   integer(kind=basis_prec) :: x1b_isd(ndim), x1b_fsd(ndim)
   integer(kind=4) :: x1b_cop(ndim)
   integer(kind=4) :: x1b_dop(ndim)
   integer(kind=1) :: x1b_phase(ndim)

!----  Internal data
   integer(kind=8) :: m, nstack
   parameter (m=7,nstack=50)
   integer(kind=8) :: i, ir, j, k, l, istack(nstack)
   integer(kind=4) :: jstack
   integer(kind=basis_prec) :: isd, fsd
   integer(kind=4) :: cop, dop
   integer(kind=1) :: phase

   logical sorted
!.......... TEST IF ALREADY SORTED......
   sorted = .true.
   do i = 2,ndim
      if(x1b_fsd(i-1) > x1b_fsd(i))then
         sorted=.false.
         exit
      end if
   end do
   if(sorted)return
!...................................

   jstack=0
   l=1
   ir=ndim
1  if(ir-l.lt.M)then
      do 12 j=l+1,ir
!        a=x1b_fsd(j)
         isd = x1b_isd(j)
         fsd = x1b_fsd(j)
         cop = x1b_cop(j)
         dop = x1b_dop(j)
         phase = x1b_phase(j)
         do 11 i=j-1,1,-1
            if(x1b_fsd(i).le.fsd)goto 2
            x1b_fsd(i+1)=x1b_fsd(i)
            x1b_isd(i+1)=x1b_isd(i)
            x1b_cop(i+1)=x1b_cop(i)
            x1b_dop(i+1)=x1b_dop(i)
            x1b_phase(i+1)=x1b_phase(i)
11       continue
         i=0
2        x1b_fsd(i+1)=fsd
         x1b_isd(i+1)=isd
         x1b_cop(i+1)=cop
         x1b_dop(i+1)=dop
         x1b_phase(i+1)=phase
12    continue
      if(jstack.eq.0)return
      ir=istack(jstack)
      l=istack(jstack-1)
      jstack=jstack-2
   else
      k=(l+ir)/2

      call switch2b(ndim,k,l+1,x1b_isd,x1b_fsd,x1b_cop,x1b_dop,x1b_phase)

      if(x1b_fsd(l+1).gt.x1b_fsd(ir))then
         call switch2b(ndim,ir,l+1,x1b_isd,x1b_fsd,x1b_cop,x1b_dop,x1b_phase)
      end if
      if(x1b_fsd(l).gt.x1b_fsd(ir))then
         call switch2b(ndim,ir,l,x1b_isd,x1b_fsd,x1b_cop,x1b_dop,x1b_phase)
      end if
      if(x1b_fsd(l+1).gt.x1b_fsd(l))then
         call switch2b(ndim,l,l+1,x1b_isd,x1b_fsd,x1b_cop,x1b_dop,x1b_phase)
      end if
      i=l+1
      j=ir
      isd = x1b_isd(l)
      fsd = x1b_fsd(l)
      cop = x1b_cop(l)
      dop = x1b_dop(l)
      phase = x1b_phase(l)
!      a=x1b_fsd(l)
3     continue
      i=i+1
      if(x1b_fsd(i).lt.fsd)goto 3
!       if(x1b_fsd(i).lt.a)goto 3
4     continue
      j=j-1
      if(x1b_fsd(j).gt.fsd)goto 4
!        if(x1b_fsd(j).gt.a)goto 4
      if(j.lt.i)goto 5
      call switch2b(ndim,j,i,x1b_isd,x1b_fsd,x1b_cop,x1b_dop,x1b_phase)
      goto 3
5     x1b_fsd(l)=x1b_fsd(j)
      x1b_isd(l)=x1b_isd(j)
      x1b_cop(l)=x1b_cop(j)
      x1b_dop(l)=x1b_dop(j)
      x1b_phase(l)=x1b_phase(j)
      x1b_fsd(j)=fsd
      x1b_isd(j)=isd
      x1b_cop(j)=cop
      x1b_dop(j)=dop
      x1b_phase(j)=phase
!         x1b_fsd(j)=a
      jstack=jstack+2
      if(jstack.gt.NSTACK) then
         print *, 'NSTACK too small in sort'
         stop 1
      end if
      if(ir-i+1.ge.j-l)then
         istack(jstack)=ir
         istack(jstack-1)=i
         ir=j-1
      else
         istack(jstack)=j-1
         istack(jstack-1)=l
         l=i
      endif
   end if
   goto 1
   return
end subroutine sort2bjumps
!==============================================================

!  subroutine sort2bjumps
!  by WEO @ LLNL - 9/2011
!  modified in 7.6.1 by CWJ@SDSU 
!
!  only 1 kind of operator
!
!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Try quick-sort algorithm isted of straight insertion        +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!    ALSO needed to introduce "switch2b" to get this right
!
! CALLED BY:
!
!  SUBROUTINES CALLED
!     switch2b

subroutine sort2bjumps1op(ndim, x2b_isd, x2b_fsd, x2b_op,x2b_phase)
   use bmpi_mod
   use precisions
   implicit none
!---   Data passed in   
   integer(kind=8) :: ndim
   integer(kind=basis_prec) :: x2b_isd(ndim), x2b_fsd(ndim)
   integer(kind=4) :: x2b_op(ndim)
   integer(kind=1) :: x2b_phase(ndim)

!----  Internal data
   integer(kind=8) :: m, nstack
   parameter (m=7,nstack=50)
   integer(kind=8) :: i, ir, j, k, l, istack(nstack)
   integer(kind=4) :: jstack
   integer(kind=basis_prec) :: isd, fsd
   integer(kind=4) :: op
   integer(kind=1) :: phase

   logical sorted
!.......... TEST IF ALREADY SORTED......
   sorted = .true.
   do i = 2,ndim
      if(x2b_fsd(i-1) > x2b_fsd(i))then
         sorted=.false.
         exit
      end if
   end do
   if(sorted)return
!...................................

   jstack=0
   l=1
   ir=ndim
1  if(ir-l.lt.M)then
      do 12 j=l+1,ir
!        a=x1b_fsd(j)
         isd = x2b_isd(j)
         fsd = x2b_fsd(j)
         op = x2b_op(j)
         phase = x2b_phase(j)
         do 11 i=j-1,1,-1
            if(x2b_fsd(i).le.fsd)goto 2
            x2b_fsd(i+1)=x2b_fsd(i)
            x2b_isd(i+1)=x2b_isd(i)
            x2b_op(i+1)=x2b_op(i)
            x2b_phase(i+1)=x2b_phase(i)
11       continue
         i=0
2        x2b_fsd(i+1)=fsd
         x2b_isd(i+1)=isd
         x2b_op(i+1)=op
         x2b_phase(i+1)=phase
12    continue
      if(jstack.eq.0)return
      ir=istack(jstack)
      l=istack(jstack-1)
      jstack=jstack-2
   else
      k=(l+ir)/2

      call switch2b1op(ndim,k,l+1,x2b_isd,x2b_fsd,x2b_op,x2b_phase)

      if(x2b_fsd(l+1).gt.x2b_fsd(ir))then
         call switch2b1op(ndim,ir,l+1,x2b_isd,x2b_fsd,x2b_op,x2b_phase)
      end if
      if(x2b_fsd(l).gt.x2b_fsd(ir))then
         call switch2b1op(ndim,ir,l,x2b_isd,x2b_fsd,x2b_op,x2b_phase)
      end if
      if(x2b_fsd(l+1).gt.x2b_fsd(l))then
         call switch2b1op(ndim,l,l+1,x2b_isd,x2b_fsd,x2b_op,x2b_phase)
      end if
      i=l+1
      j=ir
      isd = x2b_isd(l)
      fsd = x2b_fsd(l)
      op =  x2b_op(l)
      phase = x2b_phase(l)
!      a=x1b_fsd(l)
3     continue
      i=i+1
      if(x2b_fsd(i).lt.fsd)goto 3
!       if(x1b_fsd(i).lt.a)goto 3
4     continue
      j=j-1
      if(x2b_fsd(j).gt.fsd)goto 4
!        if(x1b_fsd(j).gt.a)goto 4
      if(j.lt.i)goto 5
      call switch2b1op(ndim,j,i,x2b_isd,x2b_fsd,x2b_op,x2b_phase)
      goto 3
5     x2b_fsd(l)=x2b_fsd(j)
      x2b_isd(l)=x2b_isd(j)
      x2b_op(l)=x2b_op(j)
      x2b_phase(l)=x2b_phase(j)
      x2b_fsd(j)=fsd
      x2b_isd(j)=isd
      x2b_op(j)=op
      x2b_phase(j)=phase
!         x1b_fsd(j)=a
      jstack=jstack+2
      if(jstack.gt.NSTACK) then
         print *, 'NSTACK too small in sort'
         stop 1
      end if
      if(ir-i+1.ge.j-l)then
         istack(jstack)=ir
         istack(jstack-1)=i
         ir=j-1
      else
         istack(jstack)=j-1
         istack(jstack-1)=l
         l=i
      endif
   end if
   goto 1
   return
end subroutine sort2bjumps1op

!==============================================================
!  subroutine switch
!  by WEO @ LLNL - 9/2011
!
subroutine switch2b(idim,i1,i2,iarr1,iarr2,iarr3,iarr4,iarr5)
   use precisions
   implicit none
   integer(kind=8) :: idim, i1, i2
   integer(kind=basis_prec) :: iarr1(idim),iarr2(idim)
   integer(kind=4) :: iarr3(idim),iarr4(idim)
   integer(kind=1) :: iarr5(idim)

   call switchibp(idim,i1,i2,iarr1)
   call switchibp(idim,i1,i2,iarr2)
   call switchi4(idim,i1,i2,iarr3)
   call switchi4(idim,i1,i2,iarr4)
   call switchi1(idim,i1,i2,iarr5)

   return

end subroutine switch2b

subroutine switch2b1op(idim,i1,i2,iarr1,iarr2,iarr3,iarr5)
   use precisions
   implicit none
   integer(kind=8) :: idim, i1, i2
   integer(kind=basis_prec) :: iarr1(idim),iarr2(idim)
   integer(kind=4) :: iarr3(idim)
   integer(kind=1) :: iarr5(idim)

   call switchibp(idim,i1,i2,iarr1)
   call switchibp(idim,i1,i2,iarr2)
   call switchi4(idim,i1,i2,iarr3)
   call switchi1(idim,i1,i2,iarr5)

   return

end subroutine switch2b1op

end module switch_mod
