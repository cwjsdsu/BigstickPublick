!===================================================================
!
!  file BAPPLYHLIBNEW.f90
!
!  shan: using destination distribution as hint for work load distribution
!  07/05/2016
!
!  routines to control application of the Hamiltonian
!  for parallel calculations with (optional) fragmentation of the basis and vectors
!
!  started 9/2011 by CWJ
!
!  NOTE: As of version 7.0.0, "obsolete" libraries of apply H have been discarded
!
!===================================================================

module apply_ham_omp
	
	implicit none
	
	logical :: verbose_ompham = .false.
    
	integer(kind=8), allocatable :: bucket(:)
    real, allocatable :: bucketw(:)
    integer(kind=8), allocatable :: mybasisstart(:), mybasisstop(:)
    real, allocatable :: mywork(:)

    ! moved here from opbundle by shan 9/3/2016
    integer(kind=8), allocatable :: threadStart(:,:), threadStop(:,:)

    integer :: balanceAlg = 0
contains

subroutine getLoadDist()
  use nodeinfo
  use localvectors
  use opbundles
  use jumpNbody
  use operation_stats
  use basis
  use flagger
  use contigpointervectors
  use bmpi_mod
  use sectors
  use coupledmatrixelements, only:call_spe
  
  
  implicit none
  integer :: startbundle, endbundle
  integer :: i, j
  integer :: ibundle
  integer(kind=8) :: xjmpstart, xjmpend, csdstart, csdend, cstride
  integer(kind=8) :: ncstates, csd_index, xmp, csd, xjmp
  integer(kind=8) :: psdf, nsdf, nsd, psdi, nsdi, oldpsdf
  integer(kind=8) :: statef
  integer(kind=8) :: pstride
  integer(kind=8) :: statef_start, statef_end, psdstart, psdend, nsdstart, nsdend
  integer(kind=8) :: pjmpstart, pjmpend, njmpstart, njmpend, pjmp, njmp
  integer(kind=8) :: ip, in
  integer :: which, which0
  integer :: nn
  integer(kind=8) :: bksize, total, length, average
  real :: localwt
  real :: totalw, avg, sum, remain
  integer :: p, ierr, lastbucket, adjustlen
  integer :: nprint
  integer :: curt
  integer :: is, fs, ncxsd

  integer :: count(0:15), flag, tmpi
  integer(kind=8) :: acc0, acc1

  integer(kind=8) :: myoffmin, myoffmax
  logical :: needInit
  
!  if (iproc == 0) verbose_ompham = .true.

  if (verbose_ompham) print *, "MPI Time (getLoadDist start) : ", BMPI_Wtime() - mpiStartTime

  nn = 240
  count = 0
  allocate(mybasisstart(0:ompNumThreads))
  allocate(mybasisstop(0:ompNumThreads))
  allocate(mywork(0:ompNumThreads))
  allocate(bucket(0:ompNumThreads*nn))
  allocate(bucketw(0:ompNumThreads*nn))

  startbundle = opbundlestart(iproc)
  endbundle = opbundleend(iproc)

  allocate(threadStart(startbundle:endbundle, 0:ompNumThreads))
  allocate(threadStop(startbundle:endbundle, 0:ompNumThreads))
    
  if (ompNumThreads <= 1) then
     mybasisstart(0) = v2s
     mybasisstop(0)  = v2e

     do ibundle = startbundle, endbundle
!		 if(opbundle(ibundle)%annexed)cycle
        xjmpstart = opbundle(ibundle)%pxstart
        xjmpend   = opbundle(ibundle)%pxend
        threadStart(ibundle, 0) = xjmpstart
        threadStop(ibundle, 0)  = xjmpend
     end do
     return
  end if

  ! shan estimate the range of vec2 for process iproc
  myoffmin = v2e
  myoffmax = v2s
  needInit = .true.
  
!$omp parallel do schedule(dynamic, 1) &
!$omp private(xjmpstart, xjmpend, psdstart, psdend, pjmpstart, pjmpend) &
!$omp private(ibundle, xjmp, pjmp, psdf, psdi) &
!$omp private(p2b_1sd, p2b_2sd, n2b_1sd, n2b_2sd) &
!$omp firstprivate(needInit) &
!$omp private(is, fs, ncxsd) &  
!$omp reduction(max : myoffmax) &
!$omp reduction(min : myoffmin)  
  do ibundle = startbundle, endbundle
!	 if(opbundle(ibundle)%annexed)cycle
	  

     if (needInit) then
        myoffmin = v2e
        myoffmax = v2s
        needInit = .false.
     end if

     is = opbundle(ibundle)%isector
     fs = opbundle(ibundle)%fsector
     ncxsd  = xsd(1)%sector(fs)%ncxsd
     
     if( opbundle(ibundle)%hchar == 'b')then
        p2b_1sd => p2b_fsd
        p2b_2sd => p2b_isd
        n2b_1sd => n2b_fsd
        n2b_2sd => n2b_isd
     else
        p2b_1sd => p2b_isd
        p2b_2sd => p2b_fsd
        n2b_1sd => n2b_isd
        n2b_2sd => n2b_fsd
     endif

     select case(opbundle(ibundle)%optype)
     case ('PP')
        xjmpstart = opbundle(ibundle)%pxstart
        xjmpend   = opbundle(ibundle)%pxend

        do xjmp = xjmpstart,xjmpend
           psdf = p2b_2sd(xjmp)        
           if (psdf + ncxsd > myoffmax) myoffmax = psdf + ncxsd
           if (psdf + 1 < myoffmin) myoffmin = psdf + 1
        end do

     case ('NN')
        psdstart = opbundle(ibundle)%pxstart
        psdend   = opbundle(ibundle)%pxend
        psdf = pstart(psdstart)
        if (psdf + 1 < myoffmin) myoffmin = psdf + 1
        psdf = pstart(psdend)
        if (psdf + ncxsd > myoffmax) myoffmax = psdf + ncxsd
              
     case ('SPE')
	 if(call_spe)then
        psdstart = opbundle(ibundle)%pxstart
        psdend   = opbundle(ibundle)%pxend
        psdf = pstart(psdstart)
        if (psdf + 1 < myoffmin) myoffmin = psdf + 1
        psdf = pstart(psdend)
        if (psdf + ncxsd > myoffmax) myoffmax = psdf + ncxsd
	end if
     case ('PN')
                
        pjmpstart = opbundle(ibundle)%pxstart
        pjmpend   = opbundle(ibundle)%pxend
        
        if (opbundle(ibundle)%hchar /= 'b') then
           do pjmp = pjmpstart, pjmpend
              psdf = p1b_fsd(pjmp)       ! final proton SD
              if (psdf + ncxsd > myoffmax) myoffmax = psdf + ncxsd
              if (psdf + 1 < myoffmin) myoffmin = psdf + 1
           end do  ! pjmp        
        else
            do pjmp = pjmpstart,pjmpend
               psdi = p1b_isd(pjmp)       ! initial proton slater determinant
              if (psdi + ncxsd > myoffmax) myoffmax = psdi + ncxsd
              if (psdi + 1 < myoffmin) myoffmin = psdi + 1
            end do  ! pjmp
        end if
        
     case default
        print *, "wrong opbundle type ", opbundle(ibundle)%optype
        stop 30
  
     end select
     
  end do   ! ibundle
!$omp end parallel do

  if (verbose_ompham) print *, "MPI Time (range estimation) : ", BMPI_Wtime() - mpiStartTime
  if (verbose_ompham) print *, "Range for ", iproc, " is ", myoffmin, myoffmax, v2s, v2e, (myoffmax - myoffmin) * 1.0 / (v2e - v2s)
  
  bksize = (myoffmax - myoffmin) / (nn * ompNumThreads)
  bksize = (bksize + 16) / 16 * 16
  if (verbose_ompham)  print *, "Bksize is ", iproc, bksize
  
  bucket = 0
  bucketw = 0.0
  mybasisstart = 0
  mybasisstop = 0

if (verbose_ompham) write(*, 120) "p1b_isd ", iproc, size(p1b_isd), size(p1b_fsd), size(n1b_isd), size(n1b_fsd)
120 format(a8, i4, 4i12)

!$omp parallel do schedule(dynamic, 1) &
!$omp private(localwt, xjmpstart, xjmpend, csdstart, csdend, cstride, ncstates, csd_index) &
!$omp private(ibundle, xjmp, psdf, csd, nsd, statef, which) &
!$omp private(psdstart, psdend, pstride, ip)  &
!$omp private(nsdstart, nsdend, in) &
!$omp private(pjmpstart, pjmpend, njmpstart, njmpend, pjmp, njmp, nsdf, psdi, nsdi) &
!$omp private(which0) &
!$omp private(p2b_1sd, p2b_2sd, n2b_1sd, n2b_2sd) &
!$omp reduction(+ : bucket) &
!$omp reduction(+ : bucketw)
  do ibundle = startbundle, endbundle
!	 if(opbundle(ibundle)%annexed)cycle
	  
     
     if( opbundle(ibundle)%hchar == 'b')then
        p2b_1sd => p2b_fsd
        p2b_2sd => p2b_isd
        n2b_1sd => n2b_fsd
        n2b_2sd => n2b_isd
     else
        p2b_1sd => p2b_isd
        p2b_2sd => p2b_fsd
        n2b_1sd => n2b_isd
        n2b_2sd => n2b_fsd
     endif

     select case(opbundle(ibundle)%optype)
     case ('PP')
        localwt = opwtPP
        if (opbundle(ibundle)%hchar == 'b') localwt = opwtPPb
                
        xjmpstart = opbundle(ibundle)%pxstart
        xjmpend   = opbundle(ibundle)%pxend
        csdstart = opbundle(ibundle)%nxstart
        csdend   = opbundle(ibundle)%nxend
        cstride   = opbundle(ibundle)%cstride

        ncstates = (csdend +cstride -csdstart)/cstride
        csd_index = csdstart - cstride

        do xjmp = xjmpstart,xjmpend
           psdf = p2b_2sd(xjmp)        
           which0 = (psdf - myoffmin)/bksize
           bucket(which0) = bucket(which0) + ncstates
           bucketw(which0) = bucketw(which0) + ncstates * localwt
        end do

     case ('NN')
        localwt = opwtNN
        
        psdstart = opbundle(ibundle)%pxstart
        psdend   = opbundle(ibundle)%pxend
        pstride  = opbundle(ibundle)%cstride   !
        xjmpstart = opbundle(ibundle)%nxstart
        xjmpend   = opbundle(ibundle)%nxend

        do ip = psdstart, psdend, pstride
           psdf = pstart(ip)    
           which0 = (psdf- myoffmin)/bksize
           bucket(which0) = bucket(which0) + xjmpend - xjmpstart + 1
           bucketw(which0) = bucketw(which0) + (xjmpend - xjmpstart + 1) * localwt
        end do  ! xjmp
              
     case ('SPE')
	 
	 if(call_spe)then
        localwt = opwtSPE
       
        psdstart = opbundle(ibundle)%pxstart
        psdend   = opbundle(ibundle)%pxend
        nsdstart = opbundle(ibundle)%nxstart
        nsdend   = opbundle(ibundle)%nxend
        do ip = psdstart,psdend
           statef = pstart(ip) + nstart(nsdstart)
           which0 = (statef - myoffmin)/bksize
           bucket(which0) = bucket(which0) + nsdend - nsdstart + 1
           bucketw(which0) = bucketw(which0) + (nsdend - nsdstart + 1) * localwt
        end do  !ip
	endif
        
     case ('PN')
        localwt = opwtPN
        if (opbundle(ibundle)%hchar == 'b') localwt = opwtPNb
                
        pjmpstart = opbundle(ibundle)%pxstart
        pjmpend   = opbundle(ibundle)%pxend
        njmpstart = opbundle(ibundle)%nxstart
        njmpend   = opbundle(ibundle)%nxend
        
        if (opbundle(ibundle)%hchar /= 'b') then
           do pjmp = pjmpstart, pjmpend
              psdf = p1b_fsd(pjmp)       ! final proton SD
           which0 = (psdf - myoffmin)/bksize
           bucket(which0) = bucket(which0) + njmpend - njmpstart + 1
           bucketw(which0) = bucketw(which0) + (njmpend - njmpstart + 1)*localwt
           end do  ! pjmp        
        else
            do pjmp = pjmpstart,pjmpend
               psdi = p1b_isd(pjmp)       ! initial proton slater determinant
               which0 = (psdi - myoffmin)/bksize
               bucket(which0) = bucket(which0) + njmpend - njmpstart + 1
               bucketw(which0) = bucketw(which0) + (njmpend - njmpstart + 1)*localwt
            end do  ! pjmp
        end if
        
     case default
        print *, "wrong opbundle type ", opbundle(ibundle)%optype
        stop 30
        
     end select
     
  end do   ! ibundle
!$omp end parallel do
  
  if (verbose_ompham) print *, "MPI Time (bucket) : ", BMPI_Wtime() - mpiStartTime
  if (.false.) then
      open(unit=60,file='bucket',status='unknown')
      do i = 0, nn * ompNumThreads
        write(60, 10) i, bucket(i), bucketw(i)
      end do  
10  format(i8, i12, f14.2)

    close(60)
  endif

  total = 0
  totalw = 0.0
  do i = 0, ompNumThreads*nn -1
     total = total + bucket(i)
     totalw = totalw + bucketw(i)
  end do
  avg = totalw / ompNumThreads
!  if (iproc < 1 ) write(*, 100) iproc, startbundle, endbundle, total, totalw, avg
100 format('Load Distribution', 3i8, i12, f24.2, f24.2)
!  if (iproc < 1) print *, "Bucket ", bucket, " wb ", bucketw
!  if (total == 0) print *, "ZERO ", iproc, total, startbundle, endbundle, v2s, v2e

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  i = 0
  p = 0
  sum = 0.0
  lastbucket = 0
  adjustlen = bksize
  mybasisstart(p) = myoffmin

if (balanceAlg == 0) then
  do while (i < nn * ompNumThreads)
     do while (sum < avg .and. i < nn * ompNumThreads)
        sum = sum + bucketw(i)
        i = i + 1
     end do
     if (i == nn * ompNumThreads) then
        ! distribute evenly among remaining threads
        length = myoffmax - mybasisstart(p)
        average = length / (ompNumThreads - p)
        mybasisstop(p) = mybasisstart(p) + average
        do j = p+1, ompNumThreads-1
           mybasisstart(j) = mybasisstart(p) + average * (j - p) + 1
           mybasisstop(j) = mybasisstart(j) + average
        end do
        mybasisstop(ompNumThreads-1) = v2e   ! probably not needed
        exit
     end if

     if (sum < avg * 1.03) then
        length = bksize * i
     else
        remain = sum - avg
        i = i -1
        if (i /= lastbucket) adjustlen = bksize
        length = bksize * i + (1.0 - remain / bucketw(i)) * adjustlen
        bucketw(i) = remain
        if (i == lastbucket) adjustlen = remain / bucketw(i) * bksize
     end if
     lastbucket = i
     mybasisstop(p) = myoffmin + length -1
     if (verbose_ompham) print *, "Assign ", iproc, p, i, mybasisstart(p), mybasisstop(p), & 
	                  myoffmin, myoffmax, bucketw(i), avg, sum, remain, adjustlen
     p = p + 1
     if (p < ompNumThreads) then
        mybasisstart(p) = myoffmin + length
     else
        mybasisstop(p-1) = v2e
        do j = i, nn * ompNumThreads
           if (bucketw(j) > 0) print *, "Warning not count all buckets ", iproc, i, j, avg, sum, bucket(j)
        end do
        exit
     end if
     sum = 0.0
  end do

else   ! old algorithm

  do while (i < nn * ompNumThreads)
     do while (sum < avg .and. i < nn * ompNumThreads)
        sum = sum + bucketw(i)
        i = i + 1
     end do
    
     if (sum >= avg .and. i /= nn * ompNumThreads) then
        if ((i > 0) .and. (sum/avg-1.0 > 1.0 - (sum-bucketw(i-1))/avg)) then
           i = i - 1
           sum = sum - bucketw(i)
        end if
        length = bksize * i
        mybasisstop(p) = myoffmin + length - 1
        mywork(p) = sum - totalw / ompNumThreads
        if (verbose_ompham) print *, "Assign ", iproc, p, i, sum, avg, sum/avg
        p = p + 1
        if (p < ompNumThreads) then
           mybasisstart(p) = myoffmin + length
        else
           mybasisstop(p-1) = v2e
           if (i /= nn * ompNumThreads) then
              do j = i, nn * ompNumThreads
                if (bucketw(j) > 0) print *, "Warning not count all buckets ", iproc, i, j, avg, sum, bucket(j)
              end do
           end if
           exit
        end if
        totalw = totalw - sum
        avg = totalw / (ompNumThreads - p)
        sum = 0.0
     else
        if (i /= nn * ompNumThreads) write(*, 103) "Should not be here ", iproc, p, ompNumThreads, i, sum, avg
103 format(a19, i4, i4, i4, i6, f18.6, f18.6)
        mybasisstop(p) = v2e
        do j = p+1, ompNumThreads-1
           mybasisstart(j) = v2e
           mybasisstop(j) = v2e
        end do
     end if
     
  end do

end if

  if (verbose_ompham) write(*, 104) iproc, " is ", v2s, v2e, ompNumThreads, ' start ', mybasisstart(0:ompNumThreads)
  if (verbose_ompham) write(*, 104) iproc, " is ", v2s, v2e, ompNumThreads, ' stop  ', mybasisstop(0:ompNumThreads)
104 format('DIST ', i6, a4, 2i12, i6, a7, 257i12)

  ! adjust to pstart position
  j = 1
  do i = 0, ompNumThreads-2
     do while (pstart(j) < mybasisstop(i) .and. j < size(pstart))
        j = j + 1
     end do
     if (j == size(pstart)) then
        mybasisstop(i) = pstart(j)
        print *, "Warning: rare case A happened for proc ", iproc, ompNumThreads
        print *, "Info A pos: ", i, j, size(pstart)
        print *, "Info A offset: ", pstart(j), mybasisstop(i)
        print *, "Info A basis start : ", mybasisstart
        print *, "Info A basis stop  : ", mybasisstop
        print *, "increase the nn size ", nn, " may help to avoid such case"
        exit
     end if
     if (pstart(j) == mybasisstop(i)) then
        j = j + 1
        cycle
     end if
     if (mywork(i) < 0) then
        mybasisstop(i) = pstart(j)
        mybasisstart(i+1) = pstart(j) + 1
        j = j + 1
     else
        mybasisstop(i) = pstart(j-1)
        mybasisstart(i+1) = pstart(j-1) +1
     end if
  end do

  if (verbose_ompham) write(*, 105) iproc, " is ", v2s, v2e, ompNumThreads, ' start ', mybasisstart(0:ompNumThreads)
  if (verbose_ompham) write(*, 105) iproc, " is ", v2s, v2e, ompNumThreads, ' stop  ', mybasisstop(0:ompNumThreads)
105 format('DIST1', i6, a4, 2i12, i6, a7, 257i12)

  if (verbose_ompham) print *, "MPI Time (basis) : ", BMPI_Wtime() - mpiStartTime

!    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! now we can assign opbundle to threads

    count = 0

    startbundle = opbundlestart(iproc)
    endbundle = opbundleend(iproc)

!$omp parallel do schedule(dynamic, 1) &
!$omp private(ibundle, is, fs, ncxsd) &
!$omp private(p2b_1sd, p2b_2sd, n2b_1sd, n2b_2sd) &
!$omp private(psdstart, psdend, pstride, nsdstart, nsdend) &
!$omp private(curt, ip, statef)  &
!$omp private(xjmpstart, xjmpend, csdstart, csdend, cstride) &
!$omp private(flag, i, psdf, xjmp, j, tmpi) &
!$omp private(pjmpstart, pjmpend, njmpstart, njmpend, pjmp, psdi) &    
!$omp reduction(+ : count)    
    do ibundle = startbundle, endbundle
!	   if(opbundle(ibundle)%annexed)cycle	
       
       if( opbundle(ibundle)%hchar == 'b')then
          p2b_1sd => p2b_fsd
          p2b_2sd => p2b_isd
          n2b_1sd => n2b_fsd
          n2b_2sd => n2b_isd         
       else
          p2b_1sd => p2b_isd
          p2b_2sd => p2b_fsd
          n2b_1sd => n2b_isd
          n2b_2sd => n2b_fsd
       endif

       threadStart(ibundle, :) = 0
       threadStop(ibundle, :)  = -1

       is = opbundle(ibundle)%isector
       fs = opbundle(ibundle)%fsector
       ncxsd  = xsd(1)%sector(fs)%ncxsd
     
       select case(opbundle(ibundle)%optype)
       case ('SPE')
          
	   if(call_spe)then
          psdstart = opbundle(ibundle)%pxstart
          psdend   = opbundle(ibundle)%pxend
          nsdstart = opbundle(ibundle)%nxstart
          nsdend   = opbundle(ibundle)%nxend

          curt = 0
          ip = psdstart
          statef = pstart(ip) + 1
          
          do while (ip <= psdend)
             do while (statef > mybasisstop(curt))
                curt = curt + 1
             end do
             if (curt >= ompNumThreads) then
                print *, "Wrong out of range SPE ", iproc, curt, ompNumThreads, ibundle, ip, statef
                stop 3
             end if
          
             threadStart(ibundle, curt) = ip
             
             do while (ip < psdend .and. statef < mybasisstop(curt))  
                ip = ip + 1
                statef = pstart(ip) + 1
             end do

             if (statef == mybasisstop(curt) .and. ncxsd /= 1) then
                print *, "AAA SPE split ", iproc, curt, statef, ip, psdstart, psdend, ncxsd, ompNumThreads
             end if
                   
             if (ip == psdend) then
                statef = pstart(ip) + 1
                if (statef <= mybasisstop(curt)) then
                   threadStop(ibundle, curt) = ip
                else
                   threadStop(ibundle, curt) = psdend - 1
                   threadStart(ibundle, curt+1) = psdend
                   threadStop(ibundle, curt+1)  = psdend
                end if
                exit
             else
                if (statef == mybasisstop(curt)) then
                   threadStop(ibundle, curt) = ip
                   ip = ip + 1
                else
                   threadStop(ibundle, curt) = ip -1
             end if
             end if
             
             curt = curt + 1
          end do
          
!          write(*,110) ibundle, iproc, psdstart, psdend, threadStart(ibundle, 0:3), threadStop(ibundle, 0:3)
110 format("Bundle ", i6, i6, ' psd ', 2i12, ' offset ', 4i12, 4i12)

          end if

          
       case ('PP')   
                
        xjmpstart = opbundle(ibundle)%pxstart
        xjmpend   = opbundle(ibundle)%pxend
        csdstart = opbundle(ibundle)%nxstart
        csdend   = opbundle(ibundle)%nxend
        cstride   = opbundle(ibundle)%cstride

        flag = 0
        
        i = 0
        psdf = p2b_2sd(xjmpstart)
        do while (psdf >= mybasisstop(i))
           i = i + 1
        end do
        threadStart(ibundle, i) = xjmpstart
        threadStop(ibundle, i) = xjmpstart
        
        do xjmp = xjmpstart+1,xjmpend
           psdf = p2b_2sd(xjmp)

           j = 0
           do while (psdf >= mybasisstop(j))
              j = j + 1
           end do
           
           if (j == i) then
              threadStop(ibundle, i) = xjmp
              cycle
           end if

           if (j > i) then
              threadStart(ibundle, j) = xjmp
              threadStop(ibundle, j) = xjmp
              i = j
              cycle
           end if

           if (j == i - 1) then
!              print *, "Warning PP ", iproc, ibundle, xjmp, xjmpstart, xjmpend, psdf, ncxsd, i, j, mybasisstart(i), mybasisstop(i), &
!                   " thread ", threadStart(ibundle, 0:i), threadStop(ibundle, 0:i)
              ! reset j to prevent working
              threadStart(ibundle, i) = xjmp + 1
              threadStop(ibundle, i) = xjmp 
              if (flag == 0) flag = 1
           else
              flag = 2
              do tmpi = j + 1, i-1
                 threadStart(ibundle, tmpi) = 0
                 threadStop(ibundle, tmpi) = -1
              end do
              threadStart(ibundle, i) = xjmp + 1
              threadStop(ibundle, i) = xjmp
           end if

        end do

        count(flag+3) = count(flag+3) + 1

!            print *, "Warning PP ", iproc, ibundle, xjmp, xjmpstart, xjmpend, psdf, ncxsd,  &
!                   " thread ", threadStart(ibundle, 0:6), threadStop(ibundle, 0:6)

       case ('NN')

          psdstart = opbundle(ibundle)%pxstart
          psdend   = opbundle(ibundle)%pxend
          pstride  = opbundle(ibundle)%cstride   
          xjmpstart = opbundle(ibundle)%nxstart
          xjmpend   = opbundle(ibundle)%nxend

          curt = 0
          ip = psdstart
          statef = pstart(ip) + 1

          do while (ip <= psdend)
             do while (statef > mybasisstop(curt))
                curt = curt + 1
             end do
             if (curt >= ompNumThreads) then
                print *, "Wrong out of range NN ", iproc, curt, ompNumThreads,  ip, statef
                stop 3
             end if

             threadStart(ibundle, curt) = ip

             do while (ip < psdend .and. statef < mybasisstop(curt))
                ip = ip + pstride
                statef = pstart(ip) + 1
             end do

             if (statef == mybasisstop(curt) .and. ncxsd /= 1) then
                print *, "AAA NN split ", iproc, curt, statef, ip, psdstart, psdend, ncxsd, ompNumThreads
             end if             

             if (ip == psdend) then
                statef = pstart(ip) + 1
                if (statef <= mybasisstop(curt)) then
                   threadStop(ibundle, curt) = ip
                else
                   threadStop(ibundle, curt) = psdend - pstride
                   threadStart(ibundle, curt+1) = psdend
                   threadStop(ibundle, curt+1)  = psdend
                end if
                exit
             else
                if (statef == mybasisstop(curt)) then
                   threadStop(ibundle, curt) = ip
                   ip = ip + pstride
                else
                   threadStop(ibundle, curt) = ip - pstride
             end if
             end if

             curt = curt + 1
          end do

!          write(*,112) ibundle, iproc, psdstart, psdend, opbundle(ibundle)%threadStart(0:3), opbundle(ibundle)%threadStop(0:3)
!                       pstart(opbundle(ibundle)%threadStart(0)), pstart(opbundle(ibundle)%threadStop(0))
112 format("Bundle ", i6, i6, ' nn  ', 2i12, ' offset ', 4i12, 4i12)
        
       case ('PN')

          pjmpstart = opbundle(ibundle)%pxstart
          pjmpend   = opbundle(ibundle)%pxend
          njmpstart = opbundle(ibundle)%nxstart
          njmpend   = opbundle(ibundle)%nxend
          
          if (opbundle(ibundle)%hchar /= 'b') then
             curt = 0
             ip = pjmpstart
             statef = p1b_fsd(ip) + 1

             do while (ip <= pjmpend)
                do while (statef > mybasisstop(curt))
                   curt = curt + 1
                end do
                if (curt >= ompNumThreads) then
                   print *, "Wrong out of range PNf ", iproc, curt, ompNumThreads,ip,statef,pjmpstart,pjmpend,njmpstart,njmpend
                   print *, "Wrong mybasisStart ", iproc, mybasisstart(0:ompNumThreads)
                   print *, "Wrong mybasisStop  ", iproc, mybasisstop(0:ompNumThreads)
                   print *, "Wrong pstart ", iproc, size(pstart), pstart(0:3), pstart(size(pstart)-10 : size(pstart))
                   print *, "Wrong p1b_fsd ", iproc, size(p1b_fsd), p1b_fsd(pjmpstart:pjmpend)
                   print *, "Wrong p1b_fsd ", iproc, size(n1b_fsd), n1b_fsd(njmpstart:njmpend)
                   print *, "Wrong thread start ", iproc, ibundle, threadStart(ibundle, 0:ompNumThreads)
                   print *, "Wrong thread stop ", iproc, ibundle, threadStop(ibundle, 0:ompNumThreads)
                   print *, "Wrong bucketw ", iproc, bucketw(3800:3840)
                   stop 3
                end if

                threadStart(ibundle, curt) = ip

                do while (ip < pjmpend .and. statef < mybasisstop(curt))
                   ip = ip + 1
                   statef = p1b_fsd(ip) + 1
                end do

                if (statef == mybasisstop(curt) .and. ncxsd /= 1) then
                   print *, "AAA PNf split ", iproc, curt, statef, ip, pjmpstart, pjmpend, ncxsd, ompNumThreads
                end if

                if (ip == pjmpend) then
                   statef = p1b_fsd(ip) + 1
                   if (statef <= mybasisstop(curt)) then
                      threadStop(ibundle, curt) = ip
                   else
                      threadStop(ibundle, curt) = pjmpend - 1
                      threadStart(ibundle, curt+1) = pjmpend
                      threadStop(ibundle, curt+1)  = pjmpend
                   end if
                   exit
                else
                   if (statef == mybasisstop(curt)) then
                      threadStop(ibundle, curt) = ip
                      ip = ip + 1
                      if (curt == ompNumThreads-1) then
                         threadStop(ibundle, curt) = pjmpend 
                         do pjmp = ip, pjmpend
                           if (p1b_fsd(pjmp) /= p1b_fsd(ip-1)) then
                              print *, "Wrong out of range PNf ", iproc, curt, ompNumThreads, &
                                         pjmp, statef, pjmpstart, pjmpend, njmpstart, njmpend
                              stop 3
                           end if
                        end do
                        ip = pjmpend + 1
                      end if
                   else
                      threadStop(ibundle, curt) = ip -1
                   end if
                end if

                curt = curt + 1
             end do

!          write(*,113) ibundle, iproc, pjmpstart, pjmpend, opbundle(ibundle)%threadStart(0:3), opbundle(ibundle)%threadStop(0:3)
113 format("Bundle ", i6, i6, ' pnf ', 2i12, ' offset ', 4i12, 4i12)
             
          else

             flag = 0

             i = 0
             psdi = p1b_isd(pjmpstart)
             do while (psdi >= mybasisstop(i))
                i = i + 1
             end do

             threadStart(ibundle, i) = pjmpstart
             threadStop(ibundle, i) = pjmpstart

             do pjmp = pjmpstart+1,pjmpend
                psdi = p1b_isd(pjmp)       ! initial proton slater determinant

                j = 0
                do while (psdi >= mybasisstop(j))
                   j = j + 1
                end do

                if (j == i) then
                   threadStop(ibundle, i) = pjmp
                   cycle
                end if

                if (j > i) then
                   threadStart(ibundle, j) = pjmp
                   threadStop(ibundle, j) = pjmp
                   i = j
                   cycle
                end if

                if (j == i - 1) then
                   ! reset j to prevent working
                   threadStart(ibundle, i) = pjmp + 1
                   threadStop(ibundle, i) = pjmp 
                   if (flag == 0) flag = 1
                else
                   flag = 2
                   do tmpi = j + 1, i-1
                      threadStart(ibundle, tmpi) = 0
                      threadStop(ibundle, tmpi) = -1
                   end do
                   threadStart(ibundle, i) = pjmp + 1
                   threadStop(ibundle, i) = pjmp

                end if
             end do

             count(flag) = count(flag) + 1

          end if

!          write(*,114) ibundle, iproc, pjmpstart, pjmpend, opbundle(ibundle)%threadStart(0:3), opbundle(ibundle)%threadStop(0:3)
114 format("Bundle ", i6, i6, ' pnb ', 2i12, ' offset ', 4i12, 4i12)


       case default
          print *, "wrong opbundle type ", opbundle(ibundle)%optype
          stop 30
          
       end select
       
    end do   ! ibundle
!$omp end parallel do
    
    if(verbose_ompham)then      
       print *, "PNB ratio ", count(0:2), " for proc ", iproc
       print *, "PP  ratio ", count(3:5), " for proc ", iproc
       if (verbose_ompham) print *, "MPI Time (getLoadDist finish) : ", BMPI_Wtime() - mpiStartTime
   end if
   return
end subroutine getLoadDist

!
! subroutine applyHbundled
!
! note: default is going from vecin to vecout but this can be reversed depending on hchar
!
! INPUT:
!   ibundle : which "bundle" of operations (e.g., PP between two (sub) sectors, etc)
!   vchar = 'n' (normal), 'r' (reverse)
!      (fragments) of lanczos vectors stored in module localvectors
!        in vec1 and vec2; if normal  H vec1 = vec2
!                          if reverse H vec2 = vec1

subroutine applyHbundled_omp(vchar)
  use nodeinfo
  use flagger
  use flags3body
  use precisions
  use opbundles
  use fragments
!  use shampoo
  use interaction
  use basis
  use localvectors
  use bmpi_mod
  use lanczos_info
  use contigpointervectors
  use sectors
  use diagh
  use jumpNbody
  use timing_parallel
  use timing
  use coupledmatrixelements, only:call_spe
 
  implicit none

  character(1) :: vchar, hchar 
  real(kind=8) :: sum
  integer :: tid
  integer(kind=basis_prec) :: i, j
  integer :: ierr

  integer :: ibundle
  integer :: mythread
  integer :: fs

  integer(kind=basis_prec) :: bleft, bright
  integer(kind=basis_prec) :: psdstart, psdend, nsdstart, nsdend, csdstart, csdend
  integer(kind=basis_prec) :: ip, in, statef, statei

  integer(kind=basis_prec) :: cstride, xjmpstart, xjmpend, xjmp
  integer(kind=basis_prec) :: istart, iend, csd, psd
  integer(kind=basis_prec) :: Xoplabel
  
  integer(kind=basis_prec) :: ncstates, csd_index, nsd, psdi, psdf, nsdf, nsdi

  integer(kind=basis_prec) :: pjmpstart, pjmpend, njmpstart, njmpend, pjmp, njmp
  integer :: phasep, phasen, a, b, c, d
  integer(kind=basis_prec) :: coplabel, doplabel
        
  real(kind=4) :: xme
  real(4) :: pspe, nspe
  
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads

  integer, save :: numiter = 0
  double precision :: timer1, t12

  double precision :: timers(8)

  double precision :: thrtimers(0:256), elapsed
  
  if (vchar /= 'n') then
     print *, "Vchar NN can only be n not ", vchar
     stop 4
  end if
  
  call proc_clock(iproc,'sta')

  vecin  => vec1
  vecout => vec2

  timers = 0.0
  thrtimers = 0.0
  
  numiter = numiter + 1
  timer1 = BMPI_Wtime()
  if (verbose_ompham) print *, "Entering OMP      ", timer1, " iter ", numiter

  mythread = 0

!$omp parallel &
!$omp private(ibundle, mythread, fs, bleft, bright)  &
!$omp private(psdstart, psdend, nsdstart, nsdend, csdstart, csdend)  &
!$omp private(ip, pspe, statef, in, nspe, statei) &
!$omp private(cstride, xjmpstart, xjmpend, xjmp) &
!$omp private(istart, iend, csd, psd) &
!$omp private(Xoplabel, xme) &
!$omp private(hchar) & 
!$omp private(ncstates, csd_index, nsd, psdi, psdf, nsdf, nsdi) &
!$omp private(pjmpstart, pjmpend, njmpstart, njmpend, pjmp, njmp) &
!$omp private(phasep, phasen, a, b, c, d) &
!$omp private(coplabel, doplabel) &
!$omp private(p2b_1sd, p2b_2sd, n2b_1sd, n2b_2sd) &
!$omp private(i, t12, elapsed) &
!$omp reduction(+ : timers) 
        
  mythread = omp_get_thread_num()
  do ibundle = opbundlestart(iproc), opbundleend(iproc)
!	 if(opbundle(ibundle)%annexed)cycle
	  

     t12 = BMPI_Wtime()

     fs = opbundle(ibundle)%fsector
     bleft = xsd(1)%sector(fs)%basisstart
     bright = xsd(1)%sector(fs)%basisend

     select case(opbundle(ibundle)%optype)
     case ('SPE')
        
        psdstart = opbundle(ibundle)%pxstart
        psdend   = opbundle(ibundle)%pxend
        nsdstart = opbundle(ibundle)%nxstart
        nsdend   = opbundle(ibundle)%nxend

        do ip = threadStart(ibundle, mythread), threadStop(ibundle, mythread)
           pspe = pspe_h(ip)   ! the proton contribution to the single-particle energies
           statef = pstart(ip) + nstart(nsdstart)
           do in = nsdstart,nsdend
              nspe = nspe_h(in)  ! neutron contributions to s.p.e.
              ! if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) then
              !   print *, "Wrong exe SPE ", iproc, ibundle, statef, in, mythread, mybasisstart(mythread), mybasisstop(mythread)
              !   stop 11
              ! end if
              vecout(statef) = vecout(statef) + vecin(statef)*( pspe + nspe )  ! add spes
              statef = statef+1   ! neutron SDs are contiguous
           end do  !in
        end do  !ip

        elapsed = BMPI_Wtime() - t12
        timers(1) = timers(1) + elapsed
        thrtimers(mythread) = thrtimers(mythread) + elapsed
        
     case ('NN')

        hchar = opbundle(ibundle)%hchar
        if( hchar == 'b' )then
           n2b_1sd => n2b_fsd
           n2b_2sd => n2b_isd
        else
           n2b_1sd => n2b_isd
           n2b_2sd => n2b_fsd
        endif

        psdstart = opbundle(ibundle)%pxstart
        psdend   = opbundle(ibundle)%pxend
        cstride  = opbundle(ibundle)%cstride   !
        xjmpstart = opbundle(ibundle)%nxstart
        xjmpend   = opbundle(ibundle)%nxend
        
        istart = threadStart(ibundle, mythread)
        iend   = threadStop(ibundle, mythread)
        if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO NN MATRIX ELEMENTS
           do csd = istart, iend, cstride
              psd = pstart(csd)
              do xjmp = xjmpstart,xjmpend
                 Xoplabel = n2b_op(xjmp)
                 xme = hmatnn(Xoplabel)
                 xme = xme*n2b_phase(xjmp)
                 statei = n2b_1sd(xjmp) + psd
                 statef = n2b_2sd(xjmp) + psd
                 vecout(statef) = vecout(statef) +  xme*vecin(statei)
              end do  !
           end do  ! xjmp
           
        else
           !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
           do csd = istart, iend
              psd = pstart(csd)
              do xjmp = xjmpstart,xjmpend
                 !--------- FETCH MATRIX ELEMENT...............................
                 xme = n2b_me(xjmp)
                 !--------- GET PHASE.........................................
                 xme = xme*n2b_phase(xjmp)
                 statei = n2b_1sd(xjmp) + psd
                 statef = n2b_2sd(xjmp) + psd
!                 if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) then
!                    print *, "Wrong exe NN1 ", iproc, ibundle, statef, n2b_2sd(xjmp), mythread, mybasisstart(mythread), mybasisstop(mythread)
!                    stop 11
!                 end if
                 vecout(statef) = vecout(statef) +  xme*vecin(statei)
              end do  !
           end do  ! xjmp
        end if

        elapsed = BMPI_Wtime() - t12
        timers(2) = timers(2) + elapsed
        thrtimers(mythread) = thrtimers(mythread) + elapsed

     case ('PP')

        hchar = opbundle(ibundle)%hchar
        if( hchar == 'b')then
           p2b_1sd => p2b_fsd
           p2b_2sd => p2b_isd
        else
           p2b_1sd => p2b_isd
           p2b_2sd => p2b_fsd
        endif

        csdstart = opbundle(ibundle)%nxstart
        csdend   = opbundle(ibundle)%nxend
        cstride   = opbundle(ibundle)%cstride
        xjmpstart = opbundle(ibundle)%pxstart
        xjmpend   = opbundle(ibundle)%pxend

        istart = threadStart(ibundle, mythread)
        iend   = threadStop(ibundle, mythread)
        
        if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO PP MATRIX ELEMENTS

           do xjmp = istart, iend
              !------------------ FETCH MATRIX ELEMENT...............................
              Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
              xme = hmatpp(Xoplabel)
              !--------- GET PHASE.........................................
              xme = xme*p2b_phase(xjmp)
              psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
              psdf = p2b_2sd(xjmp)   ! KSM: final P SD

              do csd = csdstart, csdend, cstride
                 nsd = nstart(csd)
                 statef = psdf + nsd
                 statei = psdi + nsd
                 vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do ! csd
           end do  ! xjmp

           if (ompNumThreads <= 1) then
              elapsed = BMPI_Wtime() - t12
              thrtimers(mythread) = thrtimers(mythread) + elapsed

              if(hchar == 'b')then
                  timers(4) = timers(4) + elapsed
              else
                  timers(3) = timers(3) + elapsed
              end if
              cycle     ! may not necessary
           end if

           ! do irreglar: left
           istart = xjmpstart
           do i = mythread -1, 0, -1
              if (threadStart(ibundle, i) > 0) then
                 istart = threadStop(ibundle, i) + 1
                 exit
              end if
           end do
           if (threadStart(ibundle, mythread) > 0) then
              iend = threadStart(ibundle, mythread) -1
           else
              iend = xjmpend
              do i = mythread + 1, ompNumThreads-1
                 if (threadStart(ibundle, i) > 0) then
                    iend = threadStart(ibundle, i) -1
                    exit
                 end if
              end do
           end if
               
           do xjmp = istart, iend
              
              psdf = p2b_2sd(xjmp)   ! KSM: final P SD

              if (psdf < mybasisstart(mythread) -1) cycle
              if (psdf >= mybasisstop(mythread)) cycle 

!              i = 0
!              do while (psdf >= mybasisstop(i))
!                 i = i + 1
!              end do
!              if (i /= mythread) cycle

              Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
              xme = hmatpp(Xoplabel)
              !--------- GET PHASE.........................................
              xme = xme*p2b_phase(xjmp)
              psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
              
              do csd = csdstart, csdend, cstride
                 nsd = nstart(csd)
                 statef = psdf + nsd
                 statei = psdi + nsd
                 vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do ! csd
           end do  ! xjmp

           ! do irregular: right

           if (threadStart(ibundle, mythread) > 0 ) then
              istart = threadStop(ibundle, mythread) +1
              iend = xjmpend
              do i = mythread + 1, ompNumThreads -1
                 if (threadStart(ibundle, i) > 0) then
                    iend = threadStart(ibundle, i) -1 
                 end if
              end do
           else
              istart = 0
              iend   = -1
           end if
           
           do xjmp = istart, iend
              
              psdf = p2b_2sd(xjmp)   ! KSM: final P SD

              if (psdf < mybasisstart(mythread) -1) cycle
              if (psdf >= mybasisstop(mythread) ) cycle
              
!              i = 0
!              do while (psdf >= mybasisstop(i))
!                 i = i + 1
!              end do
!              if (i /= mythread) cycle

              Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
              xme = hmatpp(Xoplabel)
              !--------- GET PHASE.........................................
              xme = xme*p2b_phase(xjmp)
              psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
              
              do csd = csdstart, csdend, cstride
                 nsd = nstart(csd)
                 statef = psdf + nsd
                 statei = psdi + nsd
                 vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do ! csd
           end do  ! xjmp

           
        else    ! USE STORED PP MATRIX ELEMENTS DIRECTLY........

           do xjmp = xjmpstart,xjmpend
              psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
              psdf = p2b_2sd(xjmp)   ! KSM: final P SD
              xme = p2b_me(xjmp)
              xme = xme*p2b_phase(xjmp)
                               
i = 0
do while (psdf >= mybasisstop(i))
  i = i + 1
end do
if (i /= mythread) cycle

              do csd = csdstart, csdend, cstride
                 nsd = nstart(csd)
                 statef = psdf+nsd !csd_index
!!!                 if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
                 
                 statei = psdi+ nsd !csd_index
                 vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do  ! xjmp
           end do  ! csd
        end if

        elapsed = BMPI_Wtime() - t12
        thrtimers(mythread) = thrtimers(mythread) + elapsed

        if( hchar == 'b')then
           timers(4) = timers(4) + elapsed
        else
           timers(3) = timers(3) + elapsed
        end if
     case ('PN')

        pjmpstart = opbundle(ibundle)%pxstart
        pjmpend   = opbundle(ibundle)%pxend
        njmpstart = opbundle(ibundle)%nxstart
        njmpend   = opbundle(ibundle)%nxend

        hchar = opbundle(ibundle)%hchar
        istart = threadStart(ibundle, mythread)
        iend   = threadStop(ibundle, mythread)
        if (hchar /= 'b') then
           do pjmp = istart, iend
              psdi = p1b_isd(pjmp)       ! initial proton slater determinant
              psdf = p1b_fsd(pjmp)       ! final proton SD
              phasep = p1b_phase(pjmp)   ! phase of proton jumps
              a = p1b_cop(pjmp)     ! KSM: Proton 1-body creation label
              c = p1b_dop(pjmp)     ! KSM: Proton 1-body destruction label
              !--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
              do njmp = njmpstart,njmpend
                 nsdf = n1b_fsd(njmp)
                 statef = nsdf + psdf                ! final state in combined basis
!                 if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle

                 !----------- FIND MATRIX ELEMTN --------------------------------------------
                 b = n1b_cop(njmp)  ! KSM: Neutron 1-body creation label
                 d = n1b_dop(njmp)  ! KSM: Neutron 1-body destruction label
                 phasen = n1b_phase(njmp)
                 coplabel = cpnpair(b,a)
                 doplabel = dpnpair(d,c)
                 xme = hmatpn(coplabel + doplabel)   ! get matrix element
                 xme = xme*phasep*phasen             ! multiply matrix element by jump phases
                 nsdi = n1b_isd(njmp)
                 statei = nsdi + psdi                ! initial state in combined basis
                 vecout(statef) = vecout(statef) + xme*vecin(statei)
              end do  ! njmp
           end do  ! pjmp

          elapsed = BMPI_Wtime() - t12
          thrtimers(mythread) = thrtimers(mythread) + elapsed
          timers(5) = timers(5) + elapsed
           
      else
         !---- Backward direction using hermiticity -------------------
         ! do regular first

         istart = threadStart(ibundle, mythread)
         iend   = threadStop(ibundle, mythread)
         
         do pjmp = istart, iend
            psdf = p1b_fsd(pjmp)       ! final proton SD
            psdi = p1b_isd(pjmp)       ! initial proton slater determinant
            phasep = p1b_phase(pjmp)   ! phase of proton jumps
            a  = p1b_cop(pjmp)
            c  = p1b_dop(pjmp)

            do njmp = njmpstart,njmpend
               nsdi = n1b_isd(njmp)
               nsdf = n1b_fsd(njmp)
               b  = n1b_cop(njmp)
               d  = n1b_dop(njmp)
               phasen = n1b_phase(njmp)

               statef = nsdf + psdf                  ! final state in combined basis
               statei = nsdi + psdi                  ! initial state in combined basis
!!!               if (statei < mybasisstart(mythread) .or. statei > mybasisstop(mythread)) cycle
               
               coplabel = cpnpair(b,a)
               doplabel = dpnpair(d,c)
               xme = hmatpn(coplabel + doplabel)     ! get matrix element
               xme = xme*phasep*phasen               ! multiply matrix element by jump phases
               vecout(statei) = vecout(statei) + xme*vecin(statef)
            end do  ! pjmp
         end do  ! njmp

         if (ompNumThreads <= 1) then 
           elapsed = BMPI_Wtime() - t12
           thrtimers(mythread) = thrtimers(mythread) + elapsed
           timers(6) = timers(6) + elapsed
           cycle ! may not necessary
         end if

         ! do irreglar: left
         istart = pjmpstart
         do i = mythread -1, 0, -1
            if (threadStart(ibundle, i) > 0) then
               istart = threadStop(ibundle, i) + 1
               exit
            end if
         end do
         if (threadStart(ibundle, mythread) > 0) then
            iend = threadStart(ibundle, mythread) -1
         else
            iend = pjmpend
            do i = mythread + 1, ompNumThreads-1
               if (threadStart(ibundle, i) > 0) then
                  iend = threadStart(ibundle, i) -1
                  exit
               end if
            end do
         end if
         
!          write(*,117) ibundle, iproc, ' thread ', mythread, '  left ', istart, iend
117 format("Bundle ", i6, i6, a8, i4, a7, 2i12)

         do pjmp = istart, iend

            psdi = p1b_isd(pjmp)       ! initial proton slater determinant
            if (psdi < mybasisstart(mythread) -1) cycle
            if (psdi >= mybasisstop(mythread)) cycle

! i = 0
! do while (psdi >= mybasisstop(i))
!   i = i + 1
! end do
! if (i /= mythread) cycle

            psdf = p1b_fsd(pjmp)       ! final proton SD
            phasep = p1b_phase(pjmp)   ! phase of proton jumps
            a  = p1b_cop(pjmp)
            c  = p1b_dop(pjmp)

            do njmp = njmpstart,njmpend
               nsdi = n1b_isd(njmp)
               nsdf = n1b_fsd(njmp)
               b  = n1b_cop(njmp)
               d  = n1b_dop(njmp)
               phasen = n1b_phase(njmp)

               statef = nsdf + psdf                  ! final state in combined basis
               statei = nsdi + psdi                  ! initial state in combined basis
!!!               if (statei < mybasisstart(mythread) .or. statei > mybasisstop(mythread)) cycle
               
               coplabel = cpnpair(b,a)
               doplabel = dpnpair(d,c)
               xme = hmatpn(coplabel + doplabel)     ! get matrix element
               xme = xme*phasep*phasen               ! multiply matrix element by jump phases
               vecout(statei) = vecout(statei) + xme*vecin(statef)
            end do  ! pjmp
         end do  ! njmp

         ! do irregular: right

         if (threadStart(ibundle, mythread) > 0 ) then
            istart = threadStop(ibundle, mythread) +1
            iend = pjmpend
            do i = mythread + 1, ompNumThreads -1
               if (threadStart(ibundle, i) > 0) then
                  iend = threadStart(ibundle, i) -1 
               end if
            end do
         else
            istart = 0
            iend   = -1
         end if

!          write(*,118) ibundle, iproc, ' thread ', mythread, ' right ', istart, iend
118 format("Bundle ", i6, i6, a8, i4, a7, 2i12)

         do pjmp = istart, iend

            psdi = p1b_isd(pjmp)       ! initial proton slater determinant
            if (psdi < mybasisstart(mythread) -1) cycle
            if (psdi >= mybasisstop(mythread)) cycle

! i = 0
! do while (psdi >= mybasisstop(i))
!   i = i + 1
! end do
! if (i /= mythread) cycle

            
            psdf = p1b_fsd(pjmp)       ! final proton SD
            phasep = p1b_phase(pjmp)   ! phase of proton jumps
            a  = p1b_cop(pjmp)
            c  = p1b_dop(pjmp)

            do njmp = njmpstart,njmpend
               nsdi = n1b_isd(njmp)
               nsdf = n1b_fsd(njmp)
               b  = n1b_cop(njmp)
               d  = n1b_dop(njmp)
               phasen = n1b_phase(njmp)

               statef = nsdf + psdf                  ! final state in combined basis
               statei = nsdi + psdi                  ! initial state in combined basis
!!!               if (statei < mybasisstart(mythread) .or. statei > mybasisstop(mythread)) cycle
               
               coplabel = cpnpair(b,a)
               doplabel = dpnpair(d,c)
               xme = hmatpn(coplabel + doplabel)     ! get matrix element
               xme = xme*phasep*phasen               ! multiply matrix element by jump phases
               vecout(statei) = vecout(statei) + xme*vecin(statef)
            end do  ! pjmp
         end do  ! njmp

         
        elapsed = BMPI_Wtime() - t12
        thrtimers(mythread) = thrtimers(mythread) + elapsed
        timers(6) = timers(6) + elapsed
         
      end if
      
     case default
        print *, "wrong opbundle type exe ", opbundle(ibundle)%optype
        stop 30
        
     end select

  end do

!$omp end parallel
  
  time_procSPE(iproc) = time_procSPE(iproc) + timers(1)
  time_procNN(iproc)  = time_procNN(iproc)  + timers(2)
  time_procPP(iproc)  = time_procPP(iproc)  + timers(3)
  time_procPPb(iproc) = time_procPPb(iproc) + timers(4)
  time_procPN(iproc)  = time_procPN(iproc)  + timers(5)
  time_procPNb(iproc)  = time_procPNb(iproc) + timers(6)
  
  time_spe = time_spe + timers(1)
  time_nn  = time_nn  + timers(2)
  time_pp  = time_pp  + timers(3) + timers(4)
  time1body = time1body + timers(5) + timers(6)
  
  if (verbose_ompham) write(*, 230) "Finishing OMP     ", iproc, BMPI_Wtime() - timer1, " iter ", numiter
 230 format(a19, i6, f12.4, a6, i6)
  if (verbose_ompham) print *, iproc, "Thread timing dist ", thrtimers(0:ompNumThreads)

if (.false.) then
  !............SPE......................................
  if(call_spe)then
     call clocker('spe','sta')
     call procOP_clock(iproc,'sta','SPE')
     if(noisy0) print *, "Starting applyHbundled_g SPE"
     call applySPEbundled_omp(vchar,opbundlestart(iproc), opbundleend(iproc))
     if(applypntrace)call applyPNavgdiag(vchar,opbundlestart(iproc), opbundleend(iproc))
     call clocker('spe','end')
     call procOP_clock(iproc,'end','SPE')
  end if

  if(.not.threebody)then
     !........... PP .................
     call clocker('ppo','sta')
     call procOP_clock(iproc,'sta','PPO')
     if(noisy0) print *, "Starting applyHbundled_g PP"
     call applyhPPbundled_omp(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
     call procOP_clock(iproc,'end','PPO')
     
     call procOP_clock(iproc,'sta','PPB')
     call applyhPPbundled_omp(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     call procOP_clock(iproc,'end','PPB')

     call clocker('ppo','end')
     !............NN ..................
     call clocker('nno','sta')
     call procOP_clock(iproc,'sta','NNO')

     if(noisy0) print *, "Starting applyHbundled_g NN"
     call applyhNNbundled_omp(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
     call applyhNNbundled_omp(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     call clocker('nno','end')
     call procOP_clock(iproc,'end','NNO')

     !........... PN....................................
     call clocker('one','sta')
     
     if(noisy0) print *, "Starting applyHbundled_g PN"
     if(useTR)then
        call procOP_clock(iproc,'sta','PNO')
        call applyhPNbundledTR_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPNbundledTR_g(vchar,'h',opbundlestart(iproc), opbundleend(iproc))
        call procOP_clock(iproc,'end','PNO')
        
        call procOP_clock(iproc,'sta','PNB')
        call applyhPNbundledTR_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        call procOP_clock(iproc,'end','PNB')
     else
        call procOP_clock(iproc,'sta','PNO')
        call applyhPNbundled_omp(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPNbundled_omp(vchar,'h',opbundlestart(iproc), opbundleend(iproc))
        call procOP_clock(iproc,'end','PNO')

        call procOP_clock(iproc,'sta','PNB') 
        call applyhPNbundled_omp(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        call procOP_clock(iproc,'end','PNB')
     end if
     call clocker('one','end')

  else
     !............PPP.......................................
     call clocker('ppp','sta')
     call procOP_clock(iproc,'sta','PPP')
     
     call applyhPPPbundled_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
     call applyhPPPbundled_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     call clocker('ppp','end')
     call procOP_clock(iproc,'end','PPP')

     !............PPN.......................................
     call clocker('ppn','sta')
     call procOP_clock(iproc,'sta','PPN')
     
     if(useTR)then
        call applyhPPNbundledTR_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPPNbundledTR_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     else
        call applyhPPNbundled_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPPNbundled_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        
     endif
     call procOP_clock(iproc,'end','PPN')
     
     call clocker('ppn','end')
     !............PNN.......................................
     call clocker('pnn','sta')
     call procOP_clock(iproc,'sta','PNN')
     
     if(useTR)then
        call applyhPNNbundledTR_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPNNbundledTR_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        
     else
        call applyhPNNbundled_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
        call applyhPNNbundled_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
        
     endif
     call clocker('pnn','end')
     call procOP_clock(iproc,'end','PNN')
     
     !...........NNN......................................
     call clocker('nnn','sta')
     call procOP_clock(iproc,'sta','NNN')
     
     call applyhNNNbundled_g(vchar,'f',opbundlestart(iproc), opbundleend(iproc))
     call applyhNNNbundled_g(vchar,'b',opbundlestart(iproc), opbundleend(iproc))
     call clocker('nnn','end')
     call procOP_clock(iproc,'end','NNN')
     
     
  end if  ! if threebody

end if

  ! If using thread local output storage we have to reduce to vec2
  ! unconverted bundle types are still dumping in vec2 so we just add to it
  ! from vec2thread
  if(useVec2Thread) then
     ! the schedule here is static with a chunksize.  This make sense
     ! because each chunk will have predictable runtime
!$omp parallel do                      &
!$omp    private(i, j, tid, sum)          &
!$omp    shared(vec2threadflat, v2s, v2e)  &
!$omp    schedule(static, 1024)
     do i = v2s, v2e
        sum = vec2(i)
        !! do tid = 0, ompNumThreads-1
        !!    sum = sum + vec2thread(i, tid)
        !!    vec2thread(i, tid) = 0.0
        !! end do
        
        ! sum the contributions to vec2(i)
        do j = (i - v2s), vec2threadend, vec2threadchunk
           sum = sum + vec2threadflat(j)
           vec2threadflat(j) = 0.0;
        end do
        vec2(i) = real(sum,kind=lanc_prec)
     end do
     !$omp end parallel do
  end if
  
  call proc_clock(iproc,'end')

  !........... NOW REDUCE....................
  !            WHEN LANCZOS VECTOR BROKEN, 
  !            THIS WILL GET MORE COMPLICATED
  if(noisy0) print *, "applyhbundled_g: Doing Reduce"
  if(useNewReorthog) then
     if(vchar == 'n') then
        ! 'n' corresponds with vec2 here because we are looking at the output vector
        ! call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
        ! Do reduce only onto root node.  We will be sending this data to 
        ! the slices from the isfragroot nodes (rank=0 in fcomm1, fcomm2, hcomm),
        ! so other nodes don't need it.
        call clocker('blr', 'sta')
        call BMPI_REDUCE(vec2, size(vec2), MPI_SUM, 0, fcomm2, ierr) ! in place
        call clocker('blr', 'end')
     else
        ! call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, fcomm1, ierr) ! in place
        call BMPI_REDUCE(vec1, size(vec1), MPI_SUM, 0, fcomm1, ierr) ! in place
     end if
  else
     if(nproc > 1 .and. vchar == 'n')then
        call BMPI_ALLREDUCE(vec2, size(vec2), MPI_SUM, fcomm2, ierr) ! in place
     endif
     if(nproc > 1 .and. vchar == 'r')then
        call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, fcomm1, ierr) ! in place
     end if
  end if
  
  if (numiter == 12.and.verbose_ompham) write(*, 220)  "Finishing Reduce   ", iproc, BMPI_Wtime() - timer1, " iter ", numiter
220 format(a19, i6, f12.4, a6, i6)

  if(noisy0) print *, "applyhbundled_g: Returning"
  return

end subroutine applyHbundled_omp

subroutine applyhPPbundled_omp (vchar,hchar,startbundle,endbundle )
  use nodeinfo
  use localvectors
  use system_parameters
  use jumpNbody
  use precisions
  use interaction
  use opbundles
  use fragments
  use basis
  use lanczos_info
  use flagger
  use bmpi_mod
  use contigpointervectors

  use sectors

  implicit none

  integer :: ibundle
  character(1) :: hchar,vchar
  integer :: startbundle,endbundle

  integer :: fs
  integer(kind=8) :: bleft, bright

!------------------------------------------------------------

  integer(kind=8) csdstart, csdend, csd,cstride,ncstates, csd_index
  integer(kind=8) xjmp,xjmpstart,xjmpend
  integer(kind=8):: Xoplabel
  real(kind=4)   xme, prod
  integer(kind=8) :: statei, statef,nsd,psdi,psdf
  integer(kind=basis_prec) :: statefoff, stateistart, stateistop
!-------- OpenMP functions ---------------------------------
  integer(kind=4) :: omp_get_thread_num, omp_get_num_threads
  integer(kind=4) :: num_threads
  integer(kind=8) :: istart, iend, chunk
  integer(kind=8) :: vs
  integer(4) :: mythread
  logical, save :: yes = .true.

  if (vchar /= 'n') then
     print *, "Vchar can only be n not ", vchar
     stop 4
  end if
  
  vecin  => vec1
  vecout => vec2

  if( hchar == 'f')then
     p2b_1sd => p2b_isd
     p2b_2sd => p2b_fsd
  else
     p2b_1sd => p2b_fsd
     p2b_2sd => p2b_isd
  endif

!$omp parallel private(ibundle, fs, bleft, bright)  &
!$omp          private(xjmp, Xoplabel, xme, nsd)    &
!$omp          private(mythread) &
!$omp          private(istart, iend, csd, csd_index, statef, statei)  &
!$omp          private(psdi, psdf)  &
!$omp          private(csdstart, csdend, xjmpstart, xjmpend, cstride) &
!$omp          private(ncstates) 
  do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'PP')cycle
      if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
      
      mythread = omp_get_thread_num()

      fs = opbundle(ibundle)%fsector
      bleft = xsd(1)%sector(fs)%basisstart
      bright = xsd(1)%sector(fs)%basisend
!!      if (bleft > mybasisstop(mythread) .or. bright < mybasisstart(mythread)) cycle
      
      !...... EXTRACT INFORMATION FROM OPBUNDLE ........
      csdstart = opbundle(ibundle)%nxstart
      csdend   = opbundle(ibundle)%nxend
      xjmpstart = opbundle(ibundle)%pxstart
      xjmpend   = opbundle(ibundle)%pxend
      cstride   = opbundle(ibundle)%cstride
      
      ncstates = (csdend +cstride -csdstart)/cstride
      istart = 1
      iend = ncstates
      csd_index = csdstart + (istart - 1)*cstride - cstride

      if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO PP MATRIX ELEMENTS
         ! KSM - Test both ways!!
         if(.false.) then
            do csd = istart, iend
               csd_index = csd_index + cstride
               nsd = nstart(csd_index)

               !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
               do xjmp = xjmpstart,xjmpend
                  statef = p2b_2sd(xjmp)+ nsd !csd_index
                  if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
                  !--------- FETCH MATRIX ELEMENT...............................
                  Xoplabel = p2b_op(xjmp)
                  xme = hmatpp(Xoplabel)
                  !--------- GET PHASE.........................................
                  xme = xme*p2b_phase(xjmp)
                  !---------- GET INITIAL, FINAL SDs and place in basis..............
                  statei = p2b_1sd(xjmp)+ nsd !csd_index
                  vecout(statef) = vecout(statef) + xme*vecin(statei)
               end do  ! xjmp
            end do  ! csd
         else
            ! Trial speedup
            ! result so for is mysteriously slower
            do xjmp = xjmpstart,xjmpend
               !------------------ FETCH MATRIX ELEMENT...............................
               Xoplabel = p2b_op(xjmp)    ! KSM:  index to matrix element
               xme = hmatpp(Xoplabel)
               !--------- GET PHASE.........................................
               xme = xme*p2b_phase(xjmp)
               psdi = p2b_1sd(xjmp)   ! KSM: initial P SD
               psdf = p2b_2sd(xjmp)   ! KSM: final P SD
               
               !============================================
               csd_index = csdstart + (istart - 1)*cstride - cstride
               do csd = istart, iend
                  csd_index = csd_index + cstride
                  nsd = nstart(csd_index)
                  statef = psdf + nsd
                  if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
                  
                  statei = psdi + nsd 
                  vecout(statef) = vecout(statef) + xme*vecin(statei)
               end do ! csd

               !============================================
               ! KSM speedup.  Everything out of the loop that can go
               !!               statefoff = psdf - psdi;
               !!               stateistart = csd_index + psdi + cstride;
               !!               stateistop = stateistart + (iend - istart)*cstride
               !!               ! really iterating over chunk of neutrons.   Input
               !!               ! and output have the same state offset due to neutrons
               !!               ! because they are unaffected by PP
               !!               ! There is no way to generate overlap
               !!               do statei=stateistart, stateistop, cstride
               !!                  statef = statei + statefoff;
               !!                  prod = xme * vecin(statei)
               !!                  voutp(statef) = voutp(statef) + prod
               !!               end do
            end do  ! xjmp
         end if
      else    ! USE STORED PP MATRIX ELEMENTS DIRECTLY........
         do csd = istart, iend
            csd_index = csd_index + cstride
            nsd = nstart(csd_index)
            
            !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
            do xjmp = xjmpstart,xjmpend
               statef = p2b_2sd(xjmp)+nsd !csd_index
               if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle

               !--------- FETCH MATRIX ELEMENT...............................
               xme = p2b_me(xjmp)
               !                 Xoplabel = p2b_op(xjmp)
               !                 xme = hmatpp(Xoplabel)
               !--------- GET PHASE.........................................
               xme = xme*p2b_phase(xjmp)
               !---------- GET INITIAL, FINAL SDs and place in basis..............
               statei = p2b_1sd(xjmp)+ nsd !csd_index
               vecout(statef) = vecout(statef) + xme*vecin(statei)
            end do  ! xjmp
         end do  ! csd
      end if
      !--------------OR DO HERMITIAN/BACKWARDS APPLICATION----------
   end do ! ibundle
   !$omp end parallel
   return
 end subroutine applyhPPbundled_omp


! This version has each thread write into a different output
! buffer that is reduced with other threads at the end.
subroutine applyhNNbundled_omp (vchar,hchar,startbundle,endbundle )
   use localvectors
   use nodeinfo
   use system_parameters
   use jumpNbody
   use precisions
   use interaction
   use opbundles
   use fragments
   use basis
   use lanczos_info
   use flagger
   use bmpi_mod
   use contigpointervectors

   use sectors

   implicit none

   ! arguments
   character(1),intent(in) :: hchar,vchar
   integer,intent(in) :: startbundle,endbundle
   
   ! --------- NOTE: basestart, basestop stored in module fragments
   
   !------------------------------------------------------------
   
   integer :: ibundle
   integer(kind=8) csdstart, csdend,csd, csd_index,cstride,ncstates, pstride
   integer(kind=8) xjmp,xjmpstart,xjmpend
   integer(kind=8):: Xoplabel
   real(kind=4)   xme
   integer(kind=8) :: statei, statef,psd,nsdi,nsdf,statef_start,statef_end
   integer(kind=8) :: istart, iend
   integer(kind=8) :: vs

!-------- OpenMP functions ---------------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: mythread
   real(kind=lanc_prec), pointer :: voutp(:)

   integer :: fs
   integer(kind=8) :: bleft, bright

   if (vchar /= 'n') then
      print *, "Vchar NN can only be n not ", vchar
      stop 4
   end if

   vecin  => vec1
   vecout => vec2

   if( hchar == 'f' )then
      n2b_1sd => n2b_isd
      n2b_2sd => n2b_fsd
   else
      n2b_1sd => n2b_fsd
      n2b_2sd => n2b_isd
   endif

! firstprivate gives each thread its own copy, but initializes it
! this is   better than shared for read-only vars

!$omp parallel private(ibundle, fs, bleft, bright)        &
!$omp         private(mythread)       &
!$omp         private(csdstart, csdend, cstride)          &
!$omp         private(xjmpstart, xjmpend, ncstates)       &
!$omp         private(istart, iend, csd, csd_index, psd)            &
!$omp         private(pstride)                            &
!$omp         private(xjmp, Xoplabel, xme)                &
!$omp         private(statei, statef)                     &
!$omp         private(statef_start, statef_end)           &
!$omp         firstprivate(startbundle, endbundle)        

   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'NN')cycle
      if(opbundle(ibundle)%hchar /= hchar )cycle
!	  if(opbundle(ibundle)%annexed)cycle
	  
      !....... Thread specific setup
      mythread = omp_get_thread_num()

      fs = opbundle(ibundle)%fsector
      bleft = xsd(1)%sector(fs)%basisstart
      bright = xsd(1)%sector(fs)%basisend
!!!      if (bleft > mybasisstop(mythread) .or. bright < mybasisstart(mythread)) cycle

      !...... EXTRACT INFORMATION FROM OPBUNDLE ........
      csdstart = opbundle(ibundle)%pxstart
      csdend   = opbundle(ibundle)%pxend
      cstride  = opbundle(ibundle)%cstride   !
      xjmpstart = opbundle(ibundle)%nxstart
      xjmpend   = opbundle(ibundle)%nxend

      istart = threadStart(ibundle, mythread)
      iend   = threadStop(ibundle, mythread)
      if ( .not. storeXXmesjumps ) then   ! USED INDEX TO GET TO NN MATRIX ELEMENTS
         !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
         do csd = istart, iend, cstride 
            psd = pstart(csd)
            do xjmp = xjmpstart,xjmpend
              Xoplabel = n2b_op(xjmp)
              xme = hmatnn(Xoplabel)
              xme = xme*n2b_phase(xjmp)
              statei = n2b_1sd(xjmp) + psd
              statef = n2b_2sd(xjmp) + psd
              !              if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
!              if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) then
!                 print *, "Wrong exe NN ", iproc, ibundle, statef, istart, iend, cstride, csd, n2b_2sd(xjmp), mythread, mybasisstart(mythread), mybasisstop(mythread)
!                 print *, "Wrong exe NN2 ", pstart(116627), pstart(116840)
!                 stop 11
!              end if
              vecout(statef) = vecout(statef) +  xme*vecin(statei)
            end do  !
         end do  ! xjmp
      else
         !--------- LOOP OVER 2-BODY JUMPS IN THIS SECTOR JUMPS.............
         do csd = istart, iend
            psd = pstart(csd)
            do xjmp = xjmpstart,xjmpend
              !--------- FETCH MATRIX ELEMENT...............................
              xme = n2b_me(xjmp)
              !--------- GET PHASE.........................................
              xme = xme*n2b_phase(xjmp)
              statei = n2b_1sd(xjmp) + psd
              statef = n2b_2sd(xjmp) + psd
              !---------- GET INITIAL, FINAL SDs and place in basis..............
              ! if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
              if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) then
                 print *, "Wrong exe NN1 ", iproc, ibundle, statef, n2b_2sd(xjmp), mythread, &
				   mybasisstart(mythread), mybasisstop(mythread)
                 stop 11
              end if
              vecout(statef) = vecout(statef) +  xme*vecin(statei)
            end do  !
         end do  ! xjmp
      end if
   end do ! ibundle
!$omp end parallel
   return
 end subroutine applyhNNbundled_omp

 
subroutine applyhPNbundled_omp (vchar,hchar,startbundle,endbundle )

  use localvectors
   use nodeinfo
   use system_parameters
   use jumpNbody
   use precisions
   use interaction
   use opbundles
   use fragments
   use lanczos_info
   use flagger
   use bmpi_mod
   use contigpointervectors

   use sectors

   implicit none

!-- Arguments -----------------------------------------------
   integer, intent(in) :: startbundle,endbundle
   character(1),intent(in) :: hchar,vchar

!------------------------------------------------------------
   integer :: ibundle
   integer(kind=basis_prec) :: psdi,psdf,nsdi,nsdf

   integer(kind=8) pjmp,pjmpstart,pjmpend
   integer(kind=8) njmp,njmpstart,njmpend
   integer :: a,b,c,d
   integer(kind=8) :: coplabel,doplabel
   integer :: phasep,phasen
   real(kind=4) ::   xme
   integer(kind=basis_prec) :: statei, statef
   integer(kind=4) num_threads
   integer(kind=basis_prec) :: vs

!-------- OpenMP functions/vars -----------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: mythread,numpthreads,numnthreads
   integer(8) :: startp_thread, npjmps_thread
   integer(8) :: startn_thread, nnjmps_thread

   integer(kind=8) :: bleft, bright
   integer :: fs

   if (vchar /= 'n') then
      print *, "Vchar PN can only be n not ", vchar
      stop 4
   end if
   
   vecin  => vec1
   vecout => vec2

! schedule is dynamic because runtime per step is variable
!$omp parallel private(ibundle, fs, bleft, bright)           &
!$omp     private(mythread)                          &
!$omp     private(pjmpstart, pjmpend, njmpstart, njmpend)       &
!$omp     private(numpthreads, numnthreads)                     &
!$omp     private(startp_thread, npjmps_thread)                 &
!$omp     private(pjmp, njmp, psdi, psdf)                       &
!$omp     private(a, b, c, d)                                   &
!$omp     private(phasep, phasen, coplabel, doplabel)           &
!$omp     private(xme, nsdi, nsdf, statei, statef)              

   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'PN')cycle
      if(opbundle(ibundle)%hchar /= hchar )cycle
!	 if(opbundle(ibundle)%annexed)cycle
	  
      mythread = omp_get_thread_num()
      
      fs = opbundle(ibundle)%fsector
      bleft = xsd(1)%sector(fs)%basisstart
      bright = xsd(1)%sector(fs)%basisend
!!!      if (bleft > mybasisstop(mythread) .or. bright < mybasisstart(mythread)) cycle
      
      !...... EXTRACT INFORMATION FROM OPBUNDLE ........
      !
      pjmpstart = opbundle(ibundle)%pxstart
      pjmpend   = opbundle(ibundle)%pxend
      njmpstart = opbundle(ibundle)%nxstart
      njmpend   = opbundle(ibundle)%nxend
      numpthreads = opbundle(ibundle)%numpthreads
      numnthreads = opbundle(ibundle)%numnthreads
      !
      if (hchar /= 'b') then

         !  starting position for proton 1-body jumps for this thread
         ! bit of a hack here to keep this routine similar to _orig version
         ! we take the complete thread range and ignore the divisions
         startp_thread = opbundle(ibundle)%startp_thread(0)     
         npjmps_thread = opbundle(ibundle)%startp_thread(numpthreads) - startp_thread

         ! KSM:  start/stop set up so that each proton final state appears on only one thread
         ! KSM:  prevents collison over update of voutp(statef) below
         !---------   Forward direction ------------------------------
         do pjmp = pjmpstart, pjmpend
            psdi = p1b_isd(pjmp)       ! initial proton slater determinant
            psdf = p1b_fsd(pjmp)       ! final proton SD
            phasep = p1b_phase(pjmp)   ! phase of proton jumps
            a = p1b_cop(pjmp)     ! KSM: Proton 1-body creation label
            c = p1b_dop(pjmp)     ! KSM: Proton 1-body destruction label
            !--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
            do njmp = njmpstart,njmpend
               nsdf = n1b_fsd(njmp)
               statef = nsdf + psdf                ! final state in combined basis
               if (statef < mybasisstart(mythread) .or. statef > mybasisstop(mythread)) cycle
               
               !----------- FIND MATRIX ELEMTN --------------------------------------------
               b = n1b_cop(njmp)  ! KSM: Neutron 1-body creation label
               d = n1b_dop(njmp)  ! KSM: Neutron 1-body destruction label
               phasen = n1b_phase(njmp)
               coplabel = cpnpair(b,a)
               doplabel = dpnpair(d,c)
               xme = hmatpn(coplabel + doplabel)   ! get matrix element
               xme = xme*phasep*phasen             ! multiply matrix element by jump phases
               nsdi = n1b_isd(njmp)
               statei = nsdi + psdi                ! initial state in combined basis
               vecout(statef) = vecout(statef) + xme*vecin(statei)
            end do  ! njmp
         end do  ! pjmp        
      else
         !---- Backward direction using hermiticity ------------------- 
         
         !  starting position for proton 1-body jumps for this thread
         ! bit of a hack here to keep this routine similar to _orig version
         ! we take the complete thread range and ignore the divisions
         startn_thread = opbundle(ibundle)%startn_thread(0)     
         nnjmps_thread = opbundle(ibundle)%startn_thread(numnthreads) - startn_thread
         
         !...... OPTION TO SWITCH ORDER OF LOOPS WHEN NO OpenMP.......
         
         do njmp = njmpstart,njmpend     
               nsdi = n1b_isd(njmp)
               nsdf = n1b_fsd(njmp)
               b  = n1b_cop(njmp)
               d  = n1b_dop(njmp)
               phasen = n1b_phase(njmp)
               do pjmp = pjmpstart,pjmpend
                  psdf = p1b_fsd(pjmp)       ! final proton SD
                  psdi = p1b_isd(pjmp)       ! initial proton slater determinant
                  statef = nsdf + psdf                  ! final state in combined basis
                  statei = nsdi + psdi                  ! initial state in combined basis
                  if (statei < mybasisstart(mythread) .or. statei > mybasisstop(mythread)) cycle

                  phasep = p1b_phase(pjmp)   ! phase of proton jumps
                  a  = p1b_cop(pjmp) 
                  c  = p1b_dop(pjmp)
                  !--------- LOOP OVER NEUTRON JUMPS -----------------------------------------
                  !----------- FIND MATRIX ELEMENT -------------------------------------------
                  coplabel = cpnpair(b,a)
                  doplabel = dpnpair(d,c)	   
                  xme = hmatpn(coplabel + doplabel)     ! get matrix element
                  xme = xme*phasep*phasen               ! multiply matrix element by jump phases
                  vecout(statei) = vecout(statei) + xme*vecin(statef)
               end do  ! pjmp
         end do  ! njmp
            
      end if
   end do  ! ibundle
   !$omp end parallel
   
   return
 end subroutine applyhPNbundled_omp

! Alternate version that threads across bundles, not slices
subroutine applySPEbundled_omp(vchar,startbundle,endbundle )
   use basis
   use sectors
   use diagh
   use precisions
   use lanczos_info
   use localvectors
   use nodeinfo
   use system_parameters
   use opbundles
   use fragments
   use contigpointervectors
   implicit none

   integer :: ibundle,startbundle,endbundle
   character(1) :: vchar
!------------------------------------------------------------

   integer(kind=8) nsdstart,nsdend,psdstart,psdend
   integer(kind=8) xjmp,xjmpstart,xjmpend

   integer(kind=8) :: statef
!-------- OpenMP functions ---------------------------------
   integer :: omp_get_thread_num, omp_get_num_threads
   integer :: num_threads
   integer :: mythread
   integer(kind=basis_prec) :: vs

   real(4)    ::  pspe,nspe
   integer(kind=8) :: ip,in

   integer :: fs
   integer(kind=8) :: bleft, bright

   if (vchar /= 'n') then
      print *, "Vchar NN can only be n not ", vchar
      stop 4
   end if
   
   if(vchar == 'n')then
      vecin  => vec1
      vecout => vec2
   else
      vecin  => vec2
      vecout => vec1
   end if

! note: first private says that each such private copy is initialized with 
! the value before the parallel pragma

!$omp parallel private(ibundle,fs, bleft, bright, ip,in,pspe,nspe,statef) & 
!$omp  private(psdstart, psdend, nsdstart, nsdend)        &
!$omp  private(mythread)                              
   do ibundle = startbundle,endbundle
      if(opbundle(ibundle)%optype /= 'SPE')cycle
!	  if(opbundle(ibundle)%annexed)cycle
	  
      mythread = omp_get_thread_num()
      if (threadStart(ibundle, mythread) > threadStop(ibundle, mythread)) cycle

      fs = opbundle(ibundle)%fsector
      bleft = xsd(1)%sector(fs)%basisstart
      bright = xsd(1)%sector(fs)%basisend
!!!      if (bleft > mybasisstop(mythread) .or. bright < mybasisstart(mythread)) cycle

      psdstart = opbundle(ibundle)%pxstart
      psdend   = opbundle(ibundle)%pxend
      nsdstart = opbundle(ibundle)%nxstart
      nsdend   = opbundle(ibundle)%nxend

!.... LOOP OVER PROTON SDs in that sector..........
!      do ip = psdstart,psdend
      do ip = threadStart(ibundle, mythread), threadStop(ibundle, mythread)
           pspe = pspe_h(ip)   ! the proton contribution to the single-particle energies
           statef = pstart(ip) + nstart(nsdstart)
!......... LOOP OVER NEUTRON SDS in sector jsc........................
           do in = nsdstart,nsdend
              nspe = nspe_h(in)  ! neutron contributions to s.p.e.
              vecout(statef) = vecout(statef) + vecin(statef)*( pspe + nspe )  ! add spes
              statef = statef+1   ! neutron SDs are contiguous
           end do  !in
      end do  !ip
   end do ! ibundle
!$omp end parallel 
   return
end subroutine applySPEbundled_omp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module apply_ham_omp


