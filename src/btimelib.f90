!============================================================2
!
! collected routines for timing 
!
!  SUBROUTINES:
!    clocker
!    clockout
!    clockops
!    proc_clock
!    proc_clock_out
!    procOP_clock
!    procOP_clock_out
!    timeperopmaster
!    bundle_clock
!    bundle_clock_out
!
!===================================================================
! clocker: clock to measure elapsed time
!
! SUBROUTINES CALLED:
!    cpu_time         ! native fortran time
!    data-and_time    ! native fortran time
!    MPI_Wtime        ! MPI time
! ====================================================================
subroutine clocker(timer,starter)
  use timing
  use nodeinfo
  use bmpi_mod
  implicit none

  integer(4)        :: ierr
  character(len=3)  :: timer,starter
  real(8)           :: timenow
  character(len=12) :: real_clock(3)
  integer(kind=4)   :: time_values(8)
  real(kind=8)      :: wall_time                  !  Wall time in seconds
  real(kind=8)      :: smin, shour, sday

  smin = 60.0d0
  shour = smin*60.0d0
  sday = shour*24.0d0

!--------- THE FOLLOWING MIGHT BE PLATFORM DEPENDENT -----------------------
  if ( nproc == 1 ) then
!     call BMPI_BARRIER(icomm,ierr)
     call cpu_time(timenow)     ! sequential clock
     call date_and_time(real_clock(1), real_clock(2), real_clock(3),           &
                        time_values)
     wall_time = dfloat(time_values(8))/1000.0d0 + dfloat(time_values(7)) +    &
                 dfloat(time_values(6))*smin + dfloat(time_values(5))*shour +  &
                 dfloat((time_values(3)-1))*sday
     timenow = wall_time
  else
!     call BMPI_BARRIER(icomm,ierr)
     timenow = BMPI_Wtime()  ! MPI clock
     wall_time = BMPI_Wtime()
  end if
!----------- THE REST SHOULD BE PLATFORM INDEPENDENT -----------------------
  select case (timer)

  case ('all')
     if ( starter == 'sta' ) then
        startall = timenow

        time1body        = 0.0d0
		time_pn          = 0.0d0
        time2body        = 0.0d0
        timereorthog     = 0.0d0
        time_writeeigvec = 0.0d0
        time_reduce      = 0.0d0
        time_distr       = 0.0d0
        time_applobs     = 0.0d0
		time_writeobs    = 0.0d0
        time_pp          = 0.0d0
		time_ppb         = 0.0d0
		time_pnb         = 0.0d0
        time_nn          = 0.0d0
        time_pro         = 0.0d0
        time_dot         = 0.0d0
        time_swp         = 0.0d0
        time_spe         = 0.0d0
        time_rre         = 0.0d0
        time_ror         = 0.0d0
        time_sort1b      = 0.0d0
        time_sort2b      = 0.0d0
        time_hops        = 0.0d0
        time_ppp         = 0.0d0
        time_ppn         = 0.0d0
        time_pnn         = 0.0d0
        time_nnn         = 0.0d0
        time_jmpcnt      = 0.0d0
        time_meset       = 0.0d0  
        time_hmult       = 0.0d0
        time_trdens      = 0.0d0
        time_eigvec      = 0.0d0
        time_eig         = 0.0d0
        time_intron      = 0.0d0
		  time_pivot       = 0.0d0
        
        t2               = 0.0d0
        t3               = 0.0d0
        t4               = 0.0d0
        t5               = 0.0d0
        t6               = 0.0d0
        t7               = 0.0d0
     else
        endall = timenow
     end if

  case ('bas')
     if ( starter == 'sta' ) then
        startbasis = timenow
     else
        endbasis   = timenow
     end if

  case ('ham')
     if ( starter == 'sta' ) then
        startham = timenow
     else
        endham   = timenow
     end if

  case ('lan')
     if ( starter == 'sta' ) then
        startlanczos = timenow
     else
        endlanczos   = timenow
     end if

  case ('s1b')
     if ( starter == 'sta' ) then
        starts1b = timenow
     else
        ends1b   = timenow
     end if

  case ('p2b')
     if ( starter == 'sta' ) then
        startp2b = timenow
     else
        endp2b   = timenow
     end if

  case ('n2b')
     if ( starter == 'sta' ) then
        startn2b = timenow
     else
        endn2b   = timenow
     end if


  case ('obs')
     if ( starter == 'sta' ) then
        startobs = timenow
     else
        endobs   = timenow
     end if
     case ('den')
        if ( starter == 'sta' ) then
           startdens = timenow
        else
           enddens   = timenow
        end if	 
     case ('obw')
        if ( starter == 'sta' ) then
           startwriteobs = timenow
        else
           endwriteobs   = timenow
        end if

  case ('one')
     if ( starter == 'sta' ) then
        timelast1 = timenow
     else
        time1body = time1body + timenow - timelast1
     end if
  case ('pno')
     if ( starter == 'sta' ) then
           time_pn_last = timenow
     else
           time_pn = time_pn + timenow - time_pn_last
     end if
  case ('pnb')
     if ( starter == 'sta' ) then
              time_pnb_last = timenow
     else
              time_pnb = time_pnb + timenow - time_pnb_last
     end if
  case ('two')
     if ( starter == 'sta' ) then
        timelast2 = timenow
     else
        time2body = time2body + timenow - timelast2
     end if

  case ('ppp')
     if ( starter == 'sta' ) then
        time_ppp_last = timenow
     else
        time_ppp = time_ppp + timenow - time_ppp_last
     end if

  case ('ppn')
     if ( starter == 'sta' ) then
        time_ppn_last = timenow
     else
        time_ppn = time_ppn + timenow - time_ppn_last
     end if

  case ('pnn')
     if ( starter == 'sta' ) then
        time_pnn_last = timenow
     else
        time_pnn = time_pnn + timenow - time_pnn_last
     end if

  case ('nnn')
     if ( starter == 'sta' ) then
        time_nnn_last = timenow
     else
        time_nnn = time_nnn + timenow - time_nnn_last
     end if

  case ('hmu')    ! hamiltonian multiplication
     if ( starter == 'sta' ) then
        time_hmult_last = timenow
     else
        time_hmult = time_hmult + timenow - time_hmult_last
     end if

  case ('ort')
     if ( starter == 'sta' ) then
        timelast3 = timenow
     else
        timereorthog = timereorthog + timenow - timelast3
     endif

  case ('eig')
     if ( starter == 'sta' ) then
        time_eigvec_last = timenow
     else
        time_eigvec = time_eigvec + timenow - time_eigvec_last
     end if
  case ('egv')
     if ( starter == 'sta' ) then
        time_eig_last = timenow
     else
        time_eig = time_eig + timenow - time_eig_last
     end if

  case ('trd')
     if ( starter == 'sta' ) then
        time_trd_last = timenow
     else
        time_trdens = time_trdens + timenow - time_trd_last
     end if


  case ('hop')
     if ( starter == 'sta' ) then
        time_hops = timenow
     else
        time_hops = timenow - time_hops
     endif
! Extra timing ( added by P.G.K. )..........................................
! Timing for applyobs.......................................................
  case ('aob')
     if ( starter == 'sta' ) then
        timelast4 = timenow
     else
        time_applobs = time_applobs + timenow - timelast4
     end if
! Timing for distribution...................................................
  case ('dis')
     if ( starter == 'sta' ) then
        timelast6 = timenow
     else
        time_distr = time_distr + timenow - timelast6
     end if
! Timing for reduce (MPI_ALLREDUCE) in applyh only..........................
  case ('blr')
     if ( starter == 'sta' ) then
        timelast5 = timenow
     else
        time_reduce = time_reduce + timenow - timelast5
     end if
! Timing for writeeigenvec..................................................
  case ('wev')
     if ( starter == 'sta' ) then
        timelast7 = timenow
     else
        time_writeeigvec = time_writeeigvec + timenow - timelast7
     end if
! Timing for pp-operations..................................................
  case('ppo')
     if ( starter == 'sta') then
        time_pp_last = timenow
!        time_pp_last = wall_time
     else
        time_pp = time_pp + timenow - time_pp_last
     end if

! Timing for pp-operations 'backwards'..................................................
	case('ppb')
	      if ( starter == 'sta') then
	         time_ppb_last = timenow
	      else
	         time_ppb = time_ppb + timenow - time_ppb_last
	      end if
! Timing for nn-operations..................................................
  case('nno')
     if ( starter == 'sta' ) then
        time_nn_last = timenow
     else
        time_nn = time_nn + timenow - time_nn_last
     end if
! Timing for dot product....................................................
  case('dot')
     if ( starter == 'sta' ) then
        time_dot_last = timenow
     else
        time_dot = time_dot + timenow - time_dot_last
     end if
! Timing for projecting Lnaczos vector......................................
  case('pro')
     if ( starter == 'sta' ) then
        time_pro_last = timenow
     else
        time_pro = time_pro + timenow - time_pro_last
     end if
! Timing for vector swap....................................................
  case('swp')
     if ( starter == 'sta' ) then
        time_swp_last = timenow
     else
        time_swp = time_swp + timenow - time_swp_last
     end if
! Timing for applying sps energies..........................................
  case('spe')
     if ( starter == 'sta' ) then
        time_spe_last = timenow
     else
        time_spe = time_spe + timenow - time_spe_last
     end if
! Timing for reduce in reorthogonalization..................................
  case('rre')
     if ( starter == 'sta' ) then
        time_rre_last = timenow
     else
        time_rre = time_rre + timenow - time_rre_last
     end if
! Timing for read in reorthogonalization....................................
  case('ror')
     if ( starter == 'sta' ) then
        time_ror_last = timenow
     else
        time_ror = time_ror + timenow - time_ror_last
     end if

  case('jmc')  !count jumps
     if ( starter == 'sta' ) then
        time_jmpcnt_last = timenow
     else
        time_jmpcnt = time_jmpcnt + timenow - time_jmpcnt_last
     end if
  case('des')  ! create descendents
     if ( starter == 'sta' ) then
        time_desc_last = timenow
     else
        time_desc = time_desc + timenow - time_desc_last
     end if


  case('mes')  ! set up matrix elements
     if ( starter == 'sta' ) then
        time_meset_last = timenow
     else
        time_meset = time_meset + timenow - time_meset_last
     end if

  case('mun')  ! uncouple matrix elements
     if ( starter == 'sta' ) then
        time_munc_last = timenow
     else
        time_munc = time_munc + timenow - time_munc_last
     end if

     case('int')  ! intron mapping and deletion
        if ( starter == 'sta' ) then
           time_intron_last = timenow
        else
           time_intron = time_intron + timenow - time_intron_last
        end if	 

     case('piv')  ! pivot set up
        if ( starter == 'sta' ) then
           time_pivot_last = timenow
        else
           time_pivot = time_pivot + timenow - time_pivot_last
        end if	
! Temporary timing..........................................................
  case ('tt1')
     if ( starter == 'sta' ) then
        t1start = timenow
     else
        t1end = timenow
     end if

  case ('tt8')
     if ( starter == 'sta' ) then
        t8start = timenow
     else
        t8end = timenow
     end if

  case('tt2')
     if ( starter == 'sta' ) then
        t2last = timenow
     else
        t2 = t2 + timenow - t2last
     end if

  case('tt3')
     if ( starter == 'sta' ) then
        t3last = timenow
     else
        t3 = t3 + timenow - t3last
     end if

  case('tt4')
     if ( starter == 'sta' ) then
        t4last = timenow
     else
        t4 = t4 + timenow - t4last
     end if

  case('tt5')
     if ( starter == 'sta' ) then
        t5last = timenow
     else
        t5 = t5 + timenow - t5last
     end if

  case('tt6')
     if ( starter == 'sta' ) then
        t6last = timenow
     else
        t6 = t6 + timenow - t6last
     end if

  case('tt7')
     if ( starter == 'sta' ) then
        t7last = timenow
     else
        t7 = t7 + timenow - t7last
     end if
!...........................................................................

  case default
     if ( iproc == 0 ) then
        print*,' Sorry, wrong flag for timing '
        print*,timer
     end if
     call BMPI_ABORT(icomm,101,ierr)
     stop
  end select
  return
end subroutine clocker

!=====================================================================
! CLOCKOUT
!=====================================================================
subroutine clockout(timer)
  use timing
  use io
  use nodeinfo
  use bmpi_mod
  use flagger
  implicit none

  integer(4)       :: ierr
  character(len=3) :: timer

  select case (timer)
  
  case ('tmp')   ! time so far
  if ( iproc == 0 ) then
     print*,' Time to run so far : ',endall-startall
!     if(writeout)write(resultfile,*)' Time to run after Lanczos : ',endall-startall
     if(writeout)write(logfile,*)' Time to run after lanczos: ',endall-startall		
  end if

  case ('all')
     if ( iproc == 0 ) then
		print*,'  '
        print*,' Total time to run : ',endall-startall
        if(writeout)write(resultfile,*)' Total time to run : ',endall-startall
        if(writeout)write(logfile,*)' Total time to run : ',endall-startall		
     end if

  case ('bas')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to compute basis : ',endbasis-startbasis

        if(writeout)write(logfile,*)' Time to compute basis : ',endbasis-startbasis
     end if

  case ('jmc')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to count up jumps : ',time_jmpcnt
        if(writeout)write(logfile,*)' Time to count up jumps : ',time_jmpcnt
     end if
  case ('des')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to create descendents : ',time_desc
        if(writeout)write(logfile,*)' Time to create descendents : ',time_desc
     end if

  case ('mes')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to set up for matrix elements : ',time_meset
        if(writeout)write(logfile,*)' Time to setup for matrix elements : ',time_meset
     end if
  case ('mun')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to decouple matrix elements : ',time_munc
        if(writeout)write(logfile,*)' Time to decouple matrix elements : ',time_munc
     end if
  case ('ham')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to compute jumps : ',endham-startham
        if(writeout)write(logfile,*)' Time to compute jumps : ',endham-startham
     end if

  case ('s1b')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to sort 1b-jumps : ',ends1b-starts1b
        if(writeout)write(logfile,*)' Time to sort 1b-jumps : ',ends1b-starts1b
     end if

  case ('p2b')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to sort proton 2b-jumps : ',endp2b-startp2b
        if(writeout)write(logfile,*)' Time to proton sort 2b-jumps : ',endp2b-startp2b
     end if

  case ('n2b')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to sort neutron 2b-jumps : ',endn2b-startn2b
        if(writeout)write(logfile,*)' Time to neutron sort 2b-jumps : ',endn2b-startn2b
     end if

  case ('lan')
     if ( iproc == 0 ) then
        print*,' Time to compute lanczos : ',endlanczos-startlanczos
        if(writeout)write(logfile,*)' Time to compute lanczos : ',endlanczos-startlanczos
     
     end if

  case ('ort')
     if(iproc==0)then
!        print*,' Time in 2-body : ',time_pp+time_nn
!        if(writeout)write(resultfile,*)' Time in 2-body : ',time_pp+time_nn !time2body

        if(timingdetails)print*,' Time in reorthogonalization : ',timereorthog
        if(writeout)write(logfile,*)' Time in reorthogonalization : ',timereorthog
     end if

  case ('obs')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to compute J^2, T^2 : ',endobs -startobs
        if(writeout)write(logfile,*)' Time to compute J^2, T^2 : ',endobs -startobs
     end if
     case ('den')
        if ( iproc == 0 ) then
           if(timingdetails)print*,' Time to compute 1b densities : ',enddens -startdens
           if(writeout)write(logfile,*)' Time to compute 1b densities: ',enddens -startdens
        end if	 
     case ('obw')
        if ( iproc == 0 ) then
           if(timingdetails)print*,' Time to write 1b densities : ',endwriteobs -startwriteobs
           if(writeout)write(logfile,*)' Time to write 1b densities : ',endwriteobs -startwriteobs
        end if
  case ('hop')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to compute hops ',time_hops
        if(writeout)write(logfile,*)' Time to compute hops : ',time_hops
     end if

! Added by PGK..............................................................
  case ('aob')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in applyobs : ',time_applobs
        if(writeout)write(logfile,*)' Time in applyobs : ',time_applobs
     end if
  case ('dis')
     if ( iproc == 0 ) then
        print*,' Time to calculate distribution : ',time_distr
        if(writeout)write(logfile,*)' Time to calculate distribution : ',time_distr
     end if
  case ('blr')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in reduce (in applyh only) : ',time_reduce
        if(writeout)write(logfile,*)' Time in reduce (in applyh only) : ',time_reduce
     end if
   case ('wev')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in writeeigenvec : ',time_writeeigvec
        if(writeout)write(logfile,*)' Time in writeeigenvec : ',time_writeeigvec
     end if

  case ('one')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in pn : ',time1body
        if(writeout)write(logfile,*)' Time in pn : ',time1body
     end if
  case ('pno')
     if ( iproc == 0 ) then
          if(timingdetails) print*,' Time in pn : ',time_pn
           if(writeout)write(logfile,*)' Time in pn : ',time_pn
     end if
  case ('pnb')
     if ( iproc == 0 ) then
         if(timingdetails)  print*,' Time in pn(back) : ',time_pnb
            if(writeout)write(logfile,*)' Time in pn(back) : ',time_pnb
     end if	 
  case ('hmu')
     if ( iproc == 0 ) then
        print*,' Time total in H mat-vec multiply  : ',time_hmult
        if(writeout)write(logfile,*)' Time total in H mat-vec multiply  : ',time_hmult
     end if

  case ('ppo')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in 2-body (pp) : ',time_pp
        if(writeout)write(logfile,*)' Time in 2-body (pp) : ',time_pp
     end if
  case ('ppb')
        if ( iproc == 0 ) then
         if(timingdetails)  print*,' Time in 2-body (pp)(back) : ',time_ppb
           if(writeout)write(resultfile,*)' Time in 2-body (pp)(backwards) : ',time_ppb
        end if	 
  case ('nno')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in 2-body (nn) : ',time_nn
        if(writeout)write(logfile,*)' Time in 2-body (nn) : ',time_nn
     end if

  case ('ppp')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in 3-body (ppp) : ',time_ppp
        if(writeout)write(logfile,*)' Time in 3-body (ppp) : ',time_ppp
     end if

  case ('ppn')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in 3-body (ppn) : ',time_ppn
        if(writeout)write(logfile,*)' Time in 3-body (ppn) : ',time_ppn
     end if

  case ('pnn')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in 3-body (pnn) : ',time_pnn
        if(writeout)write(logfile,*)' Time in 3-body (pnn) : ',time_pnn
     end if

  case ('nnn')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in 3-body (nnn) : ',time_nnn
        if(writeout)write(logfile,*)' Time in 3-body (nnn) : ',time_nnn
     end if

  case ('dot')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in normalization : ',time_dot
        if(writeout)write(logfile,*)' Time in normalization : ',time_dot
     end if
  case ('pro')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in projection : ',time_pro
        if(writeout)write(logfile,*)' Time in projection : ',time_pro
     end if
  case ('swp')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time for vector swap : ',time_swp
        if(writeout)write(logfile,*)' Time for vector swap : ',time_swp
     end if
  case ('spe')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to apply sp energies : ',time_spe
        if(writeout)write(logfile,*)' Time to apply sp energies : ',time_spe
     end if
  case ('rre')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time in reduce (reorthog. only) : ',time_rre
        if(writeout)write(logfile,*)' Time in reduce (reorthog. only) : ',time_rre
     end if
  case ('ror')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time for read in reorthog. : ',time_ror
        if(writeout)write(logfile,*)' Time for read in reorthog. : ',time_ror
     end if
  case ('egv')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time spent diagonalizing. : ',time_eig
        if(writeout)write(logfile,*)' Time spent diagonalizing. : ',time_eig
     end if
  case ('eig')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to build eigenvectors. : ',time_eigvec
        if(writeout)write(logfile,*)' Time to build eigenvectors. : ',time_eigvec
     end if
  case ('trd')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to build trdens wave functions : ',time_trdens
        if(writeout)write(logfile,*)' Time to build trdens wave functions : ',time_trdens
     end if

  case ('int')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to compute and delete unused jump storage : ',time_intron
        if(writeout)write(logfile,*)' Time to compute and delete unused jump storage : ',time_intron
     end if

  case ('piv')
     if ( iproc == 0 ) then
        if(timingdetails)print*,' Time to set up pivot : ',time_pivot
        if(writeout)write(logfile,*)' Time to set up pivot : ',time_pivot
     end if

! Temporary.................................................................
  case ('tt1')
     if ( iproc == 0 ) then
        print*,' T1 = ',t1end - t1start
     end if
  case ('tt2')
     if ( iproc == 0 ) then
        print*,' T2 = ',t2
     end if
  case ('tt3')
     if ( iproc == 0 ) then
        print*,' T3 = ',t3
     end if
  case ('tt4')
     if ( iproc == 0 ) then
        print*,' T4 = ',t4
     end if
  case ('tt5')
     if ( iproc == 0 ) then
        print*,' T5 = ',t5
     end if
  case ('tt6')
     if ( iproc == 0 ) then
        print*,' T6 = ',t6
     end if
  case ('tt7')
     if ( iproc == 0 ) then
        print*,' T7 = ',t7
     end if
  case ('tt8')
     if ( iproc == 0 ) then
        print*,' T8 = ',t8end - t8start
     end if
!...........................................................................

  case default
      if ( iproc == 0 ) then
         print*,' Sorry, wrong flag for timing '
         print*,timer
      end if
      call BMPI_ABORT(icomm,101,ierr)
      stop
  end select

  return
end subroutine clockout


!===================================================================
! Timer: clock to measure elapsed time
! ====================================================================
subroutine proc_clock(iprocs,starter)
  use timing_parallel
  use nodeinfo
  use bmpi_mod
  implicit none

  integer(4)        :: ierr
  character(len=3)  :: starter
  real(8)           :: timenow
  character(len=12) :: real_clock(3)
  integer(kind=4)   :: time_values(8)
  real(kind=8)      :: wall_time                  !  Wall time in seconds
  real(kind=8)      :: smin, shour, sday
  integer           :: iprocs    ! dummy for which processor
  integer           :: aerr
  smin = 60.0d0
  shour = smin*60.0d0
  sday = shour*24.0d0

!--------- THE FOLLOWING MIGHT BE PLATFORM DEPENDENT -----------------------
  if ( nproc == 1 ) then
!     call BMPI_BARRIER(icomm,ierr)
     call cpu_time(timenow)     ! sequential clock
     call date_and_time(real_clock(1), real_clock(2), real_clock(3),           &
                        time_values)
     wall_time = dfloat(time_values(8))/1000.0d0 + dfloat(time_values(7)) +    &
                 dfloat(time_values(6))*smin + dfloat(time_values(5))*shour +  &
                 dfloat((time_values(3)-1))*sday
     timenow = wall_time
  else
!     call BMPI_BARRIER(icomm,ierr)
     timenow = BMPI_Wtime()  ! MPI clock
     wall_time = BMPI_Wtime()
  end if

  select case (starter)

  case ('set')
     allocate( time_Ham_MPI(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror("proc_clock")
     time_Ham_MPI(:) = 0.00
  case ('sta')
     timelast = timenow

  case ('end')
     time_Ham_MPI(iprocs) = time_Ham_MPI(iprocs) + timenow-timelast

  end select

  return
  end subroutine proc_clock
!===============================================
!
! called by: lanczos_p after all lanczos iterations
!
! AFTER THIS procOP_clock_out is called (see below)
!

  subroutine proc_clock_out

  use timing_parallel
  use nodeinfo
  use bmpi_mod
  use butil_mod
  use io
  use flags3body
  use operation_stats
  implicit none

  integer iprocs
  real(8)  :: avg_time
  real(8)  :: max_time
  integer tag
  integer(4) :: ierr
  integer stat(MPI_STATUS_SIZE)

!------ EVERY PROC SEND THEIR DATA TO PROC ZERO----

  if(iproc >0)then
    tag = iproc
    call BMPI_Send(time_Ham_MPI(iproc),1,0,tag,icomm,ierr) 

  end if
  if(iproc == 0)then
    open(unit=39,file='proctimes.bigstick',status='unknown')
!....... WRITE A HEADER TO proctimes.bigstick........
    write(39,*)'# timing weights assumed '
	if(.not.threebody)then
		write(39,1002)'SPE',opwtSPE
		write(39,1002)'PP ',opwtPP
		write(39,1002)'PPb',opwtPPb
		write(39,1002)'NN ',opwtNN
		write(39,1002)'PN ',opwtPN
		write(39,1002)'PNb',opwtPNb		
	else
		write(39,1002)'SPE',opwtSPE
		write(39,1002)'PPP',opwtPPP
		write(39,1002)'NNN',opwtNNN
		write(39,1002)'PPN',opwtPPN
		write(39,1002)'PNN',opwtPNN			
	endif 
    write(39,*)'#  proc    actual time     '
	write(39,*)'#  optype  total time     time per op (ns)   expected total time'	
!....... END OF HEADER................	
    max_time = 0.0
    avg_time = 0.0
    do iprocs = 0,nprocs-1
      if(iprocs > 0)then
        tag = iprocs
        call BMPI_Recv(time_Ham_MPI(iprocs),1,iprocs,tag,icomm,stat,  ierr) 
      end if
      avg_time = avg_time + time_Ham_MPI(iprocs)
      max_time = bmax(max_time, time_Ham_MPI(iprocs))
      
    end do

    avg_time =avg_time / float(nprocs)
    print*,' Average time per processor is ',avg_time, " max ", max_time
    print*,' Overall efficiency is ',avg_time/max_time
    print*,' (Individual processor data written to file proctimes.bigstick )'
    write(logfile,*)' Average time per processor is ',avg_time, " max ", max_time
    write(logfile,*)' Overall efficiency is ',avg_time/max_time
  end if
  return
  1002 format(1x,a3,2x,f7.2)
end subroutine proc_clock_out

!===================================================================
! Timer: clock to measure elapsed time
! ====================================================================
subroutine procOP_clock(iprocs,starter,optype)
  use timing_parallel
  use nodeinfo
  use bmpi_mod
  implicit none

  integer(4)        :: ierr
  character(len=3)  :: starter
  real(8)           :: timenow
  character(3)      ::optype
  character(len=12) :: real_clock(3)
  integer(kind=4)   :: time_values(8)
  real(kind=8)      :: wall_time                  !  Wall time in seconds
  real(kind=8)      :: smin, shour, sday
  integer iprocs    ! dummy for which processor
  integer :: aerr
  character (len=*), parameter :: memmsg = "procOP_clock"
  smin = 60.0d0
  shour = smin*60.0d0
  sday = shour*24.0d0

!--------- THE FOLLOWING MIGHT BE PLATFORM DEPENDENT -----------------------
  if ( nproc == 1 ) then
     call cpu_time(timenow)     ! sequential clock
     call date_and_time(real_clock(1), real_clock(2), real_clock(3),           &
                        time_values)
     wall_time = dfloat(time_values(8))/1000.0d0 + dfloat(time_values(7)) +    &
                 dfloat(time_values(6))*smin + dfloat(time_values(5))*shour +  &
                 dfloat((time_values(3)-1))*sday
     timenow = wall_time
  else
     timenow = BMPI_Wtime()  ! MPI clock
     wall_time = BMPI_Wtime()
  end if
  aerr = 0  ! make sure initialized for error detection.
  select case (starter)

  case ('set')
     if(.not.allocated(time_procSPE))allocate( time_procSPE(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procSPE(:) = 0.00
	  if(allocated(time_procPP)) deallocate(time_procPP)
	  allocate(time_procPP(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procPP(:) = 0.00	 
     if(.not.allocated(time_procNN)) allocate( time_procNN(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procNN(:) = 0.00	 
     if(.not.allocated(time_procPN)) allocate( time_procPN(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procPN(:) = 0.00

     !added by shan:
     allocate( time_procPPb(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procPPb(:) = 0.00
     allocate( time_procPNb(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procPNb(:) = 0.00

     if(.not.allocated(time_procPPP)) allocate( time_procPPP(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procPPP(:) = 0.00
     if(.not.allocated(time_procPPN)) allocate( time_procPPN(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procPPN(:) = 0.00
     if(.not.allocated(time_procPNN)) allocate( time_procPNN(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procPNN(:) = 0.00
     if(.not.allocated(time_procNNN)) allocate( time_procNNN(0:nprocs-1), stat=aerr)
     if(aerr /= 0) call memerror(memmsg)
     time_procNNN(:) = 0.00

  case ('sta')
     timelastprocop = timenow

  case ('end')
      select case (optype)
      case ('SPE')
          time_procSPE(iprocs) = time_procSPE(iprocs) + timenow-timelastprocop
      case ('PPO')
          time_procPP(iprocs) = time_procPP(iprocs) + timenow-timelastprocop
      case ('PNO')
          time_procPN(iprocs) = time_procPN(iprocs) + timenow-timelastprocop
      case ('NNO')
          time_procNN(iprocs) = time_procNN(iprocs) + timenow-timelastprocop

      !added by shan:
       case ('PPB')
          time_procPPb(iprocs) = time_procPPb(iprocs) + timenow-timelastprocop
      case ('PNB')
          time_procPNb(iprocs) = time_procPNb(iprocs) + timenow-timelastprocop
          
      case ('PPP')
          time_procPPP(iprocs) = time_procPPP(iprocs) + timenow-timelastprocop
      case ('PPN')
          time_procPPN(iprocs) = time_procPPN(iprocs) + timenow-timelastprocop
      case ('PNN')
          time_procPNN(iprocs) = time_procPNN(iprocs) + timenow-timelastprocop
      case ('NNN')
          time_procNNN(iprocs) = time_procNNN(iprocs) + timenow-timelastprocop


      end select
  end select

  return
  end subroutine procOP_clock
!===============================================
!
  subroutine procOP_clock_out

  use timing_parallel
  use nodeinfo
  use opbundles
  use fragments
  use lanczos_info
  use flags3body
  use bmpi_mod
  use coupledmatrixelements,only:call_spe

  implicit none

  integer iprocs
  integer(8) :: localops
  integer  :: ibundle
  integer tag
  integer(4) :: ierr
  integer stat(MPI_STATUS_SIZE)
  real(8) :: tmp(10)
  real(8) ::expectedtime
  
  call BMPI_BARRIER(icomm,ierr) 
  if(iproc >0)then
    tag = iproc
     tmp(1) = time_procSPE(iproc)
     tmp(2) = time_procPP(iproc)
     tmp(3) = time_procNN(iproc)
     tmp(4) = time_procPN(iproc)
     
     tmp(5) = time_procPPb(iproc)
     tmp(6) = time_procPNb(iproc)

     tmp(7) = time_procPPP(iproc)
     tmp(8) = time_procPPN(iproc)
     tmp(9) = time_procPNN(iproc)
     tmp(10) = time_procNNN(iproc)
     
     call BMPI_SEND(tmp,10,0,tag,icomm,ierr) 
  end if

  if(iproc == 0)then
    do iprocs = 1,nprocs-1
        tag = iprocs
        call BMPI_Recv(tmp,10,iprocs,tag,icomm,stat,ierr)
        time_procSPE(iprocs) = tmp(1)

        if(.not.threebody)then
           time_procPP(iprocs)  = tmp(2)
           time_procNN(iprocs)  = tmp(3)
           time_procPN(iprocs)  = tmp(4)

           time_procPPb(iprocs) = tmp(5)
           time_procPNb(iprocs) = tmp(6)
        else
           time_procPPP(iprocs) = tmp(7)
           time_procPPN(iprocs) = tmp(8)
           time_procPNN(iprocs) = tmp(9)
           time_procNNN(iprocs) = tmp(10)
        endif
       
    end do
    do iprocs = 0,nprocs-1
		expectedtime = 0.d0
       write(39,*)' processor, actual time ',iprocs, time_Ham_MPI(iprocs)
       localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'SPE')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        if(call_spe)write(39,*)'SPE',0.0
        else
      write(39,*)'SPE',time_procSPE(iprocs),time_procSPE(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtSPE*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtSPE
        end if

        if(.not.threebody)then
       localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'PP')cycle
            if(opbundle(ibundle)%hchar == 'b')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        write(39,*)'PP',0.0
        else
      write(39,*)'PP',time_procPP(iprocs),time_procPP(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPP  *1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPP
				 
        end if

!added by shan        
       localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
            if(opbundle(ibundle)%optype /= 'PP')cycle
            if(opbundle(ibundle)%hchar /= 'b')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        write(39,*)'PPb',0.0
        else
      write(39,*)'PPb',time_procPPb(iprocs),time_procPPb(iprocs)/dfloat(actual_iterations*localops)*1.e9 , & 
	             dfloat(localops*actual_iterations)*opwtPPb *1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPPb
        end if
  
       localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'NN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        write(39,*)'NN',0.0
        else
      write(39,*)'NN',time_procNN(iprocs),time_procNN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtNN*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtNN
        end if

       localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'PN')cycle
            if(opbundle(ibundle)%hchar == 'b')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        write(39,*)'PN',0.0
        else
      write(39,*)'PN',time_procPN(iprocs),time_procPN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPN*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPN
        endif

        !added by shan:
       localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
            if(opbundle(ibundle)%optype /= 'PN')cycle
            if(opbundle(ibundle)%hchar /= 'b')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        write(39,*)'PNb',0.0
        else
      write(39,*)'PNb',time_procPNb(iprocs),time_procPNb(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPNb*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPNb
        endif
        
        else
     localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'PPP')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        write(39,*)'PPP',0.0
        else
      write(39,*)'PPP',time_procPPP(iprocs),time_procPPP(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPPP *1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPPP
        end if
       localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'NNN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        write(39,*)'NNN',0.0
        else
      write(39,*)'NNN',time_procNNN(iprocs),time_procNNN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtNNN*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtNNN
        end if

       localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'PPN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        write(39,*)'PPN',0.0
        else
      write(39,*)'PPN',time_procPPN(iprocs),time_procPPN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPPN*1e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPPN
        endif

       localops = 0
         do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'PNN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
        write(39,*)'PNN',0.0
        else
      write(39,*)'PNN',time_procPNN(iprocs),time_procPNN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPNN*1e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPNN
        endif

        end if
		write(39,*)' Expected total time = ',expectedtime*1e-9

    end do

    close(39)
  end if
      call BMPI_BARRIER(icomm,ierr) 

  return

  end subroutine procOP_clock_out


!!=============================================================
!
! Subroutine to read in and write out time per operation
! Can use weightings as needed to help improve distributions
!
! replaces old subroutine "clockops" from above
!
! The distribution of operations depends on how much an operation "costs"
! BIGSTICK assumes a default "weight" for each operation but this 
! can be modified based upon a prior run, in order to improve efficiency
! This data is written to the file "timinginfo.bigstick"
!
! At startup, BIGSTICK looks for this file; if it does not exist, 
! or if flag modifytimeweights = F, then uses defaults
! Otherwise, the weights are modified. After the run, both prior and current
! weights are printed out. In particular the file looks like
!  OP         DEFAULT      PRIOR     CURRENT
!  SPE        .....
!  PP         ....
!  etc
! 
!  All times are in ns (nanoseconds) per operation
!
! IMPORTANT:  Must be called AFTER other clocks have been written out
!
subroutine timeperopmaster(starter)
  use flags3body
  use timing
  use timeinfo
  use distrinfo2
  use lanczos_info
  use nodeinfo
  use basis
  use timing_parallel
  use opbundles
  use fragments
  use program_info
  use operation_stats
  use bmpi_mod
  use coupledmatrixelements,only:call_spe
  implicit none

  character(3):: starter
  character(20) :: line
  real(8)  :: tdefault
  real(8)  :: tottimeX,tpoX
  integer(8) :: localops
  integer ibundle
  integer ip,ierr
  call BMPI_BARRIER(icomm,ierr) 

  if(starter=='set')then
     if(iproc == 0)then
!........ ATTEMPT TO OPEN FILE 'timinginfo.bigstick '.......
     open(unit=timefile,file='timinginfo.bigstick',status='old',err=1111)
     priortimeinfo = .true.
!....... can read in here
     read(timefile,'(a)',err=1112,end=1112)line
     read(timefile,1001,err=1112,end=1112)tdefault,tpoSPEprior, tpoSPE
     read(timefile,1001,err=1112,end=1112)tdefault,tpoPPprior, tpoPP
     read(timefile,1001,err=1112,end=1112)tdefault,tpoPPbprior, tpoPPb
     read(timefile,1001,err=1112,end=1112)tdefault,tpoNNprior, tpoNN
     read(timefile,1001,err=1112,end=1112)tdefault,tpoPNprior, tpoPN
     read(timefile,1001,err=1112,end=1112)tdefault,tpoPNbprior, tpoPNb
     read(timefile,1001,err=1112,end=1112)tdefault,tpoPPPprior, tpoPPP
     read(timefile,1001,err=1112,end=1112)tdefault,tpoNNNprior, tpoNNN
     read(timefile,1001,err=1112,end=1112)tdefault,tpoPPNprior, tpoPPN
     read(timefile,1001,err=1112,end=1112)tdefault,tpoPNNprior, tpoPNN

     if (modifytimewts)then
!................. CHECK FOR ERRORS......................
     if(tpoSPE < 1.0 .or. tpoSPE > 50.)then
         print*,' You might want to check for errors in or delete file timinginfo.bigstick '
     end if
     if(tpoPP < 1.0 .or. tpoPP > 50.)then
         print*,' You might want to check for errors in or delete file timinginfo.bigstick '
     end if
     if(tpoNN < 1.0 .or. tpoNN > 50.)then
         print*,' You might want to check for errors in or delete file timinginfo.bigstick '
     end if
     if(tpoPN < 1.0 .or. tpoPN > 50.)then
         print*,' You might want to check for errors in or delete file timinginfo.bigstick '
     end if
     if(tpoPPP < 1.0 .or. tpoPPP > 50.)then
         print*,' You might want to check for errors in or delete file timinginfo.bigstick '
     end if
     if(tpoPPN < 1.0 .or. tpoPPN > 50.)then
         print*,' You might want to check for errors in or delete file timinginfo.bigstick '
     end if
     if(tpoPNN < 1.0 .or. tpoPNN > 50.)then
         print*,' You might want to check for errors in or delete file timinginfo.bigstick '
     end if
     if(tpoNNN < 1.0 .or. tpoNNN > 50.)then
         print*,' You might want to check for errors in or delete file timinginfo.bigstick '
     end if
!.......................................................................

      opwtSPE= real(tpoSPE, kind(opwtSPE))
      opwtPP = real(tpoPP, kind(opwtPP))
      opwtPPb = real(tpoPPb, kind(opwtPPb))
      opwtNN = real(tpoNN, kind(opwtNN))
      opwtPN = real(tpoPN, kind(opwtPN))
      opwtPNb = real(tpoPNb, kind(opwtPNb))
      opwtPPP= real(tpoPPP, kind(opwtPPP))
      opwtNNN= real(tpoNNN, kind(opwtNNN))
      opwtPPN= real(tpoPPN, kind(opwtPPN))
      opwtPNN= real(tpoPNN, kind(opwtPNN))

      if (opwtSPE == 0.0) opwtSPE = tpoSPEdefault
      if (opwtPP == 0.0)  opwtPP  = tpoPPdefault
      if (opwtNN == 0.0)  opwtNN  = tpoNNdefault
      if (opwtPN == 0.0)  opwtPN  = tpoPNdefault
      if (opwtPPb == 0.0)  opwtPPb  = tpoPPbdefault
      if (opwtPNb == 0.0)  opwtPNb  = tpoPNbdefault
     else
      opwtSPE= tpoSPEdefault
      opwtPP = tpoPPdefault
      opwtPPb = tpoPPbdefault
      opwtNN = tpoNNdefault
      opwtPN = tpoPNdefault
      opwtPNb = tpoPNbdefault
      opwtPPP= tpoPPPdefault
      opwtNNN= tpoNNNdefault
      opwtPPN= tpoPPNdefault
      opwtPNN= tpoPNNdefault

     end if


1001 format(6x,3f7.2)
1002 format(1x,a3,2x,3f7.2)
1003 format(1x,a3,2x,3f7.2,' (NOT MODIFIED) ')
     goto 2020

1111 continue
     open(unit=timefile,file='timinginfo.bigstick',status='new')
!....... USE DEFAULTS.........
1112 continue
     priortimeinfo = .false.

     tpoSPEprior = tpoSPEdefault
     tpoSPE      = tpoSPEdefault
     tpoPPprior  = tpoPPdefault
     tpoPP       = tpoPPdefault

     tpoPPbprior  = tpoPPbdefault
     tpoPPb       = tpoPPbdefault
          
     tpoNNprior  = tpoNNdefault
     tpoNN       = tpoNNdefault
     tpoPNprior  = tpoPNdefault
     tpoPN       = tpoPNdefault

     tpoPNbprior  = tpoPNbdefault
     tpoPNb       = tpoPNbdefault
     
     tpoPPPprior = tpoPPPdefault
     tpoPPP      = tpoPPPdefault
     tpoNNNprior = tpoNNNdefault
     tpoNNN      = tpoNNNdefault
     tpoPPNprior = tpoPPNdefault
     tpoPPN      = tpoPPNdefault
     tpoPNNprior = tpoPNNdefault
     tpoPNN      = tpoPNNdefault

     opwtSPE= tpoSPEdefault
     opwtPP = tpoPPdefault

     opwtPPb = tpoPPbdefault
     
     opwtNN = tpoNNdefault
     opwtPN = tpoPNdefault

     opwtPNb = tpoPNbdefault
     
     opwtPPP= tpoPPPdefault
     opwtNNN= tpoNNNdefault
     opwtPPN= tpoPPNdefault
     opwtPNN= tpoPNNdefault

     end if 
2020 continue

      call BMPI_BCAST(opwtSPE,1,0,icomm,ierr)
      call BMPI_BCAST(opwtPP,1,0,icomm,ierr)
      call BMPI_BCAST(opwtPPb,1,0,icomm,ierr)
      call BMPI_BCAST(opwtNN,1,0,icomm,ierr)
      call BMPI_BCAST(opwtPN,1,0,icomm,ierr)
      call BMPI_BCAST(opwtPNb,1,0,icomm,ierr)
      
      call BMPI_BCAST(opwtPPP,1,0,icomm,ierr)
      call BMPI_BCAST(opwtPPN,1,0,icomm,ierr)
      call BMPI_BCAST(opwtPNN,1,0,icomm,ierr)
      call BMPI_BCAST(opwtNNN,1,0,icomm,ierr)
      call BMPI_BARRIER(icomm,ierr) 

      if (iproc == 0) print *, "WEIGHTS for SPE, PP, NN, PN ", opwtSPe, opwtPP, opwtNN, opwtPN 
     return
  end if

  if(starter=='out' .and. iproc==0)then
     rewind(timefile)
     write(timefile,*)'OP  DEFAULT   PRIOR  CURRENT  '
!............ WRITE OUT SPE............
     localops = 0
     do ibundle = 1,nopbundles
             if(opbundle(ibundle)%optype /= 'SPE')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
     end do  ! ibundle
     tottimeX = 0.d0
     do ip = 0,nprocs-1
       tottimeX = tottimeX + time_procSPE(ip)
     end do
     tpoX = tottimeX/(actual_iterations*localops)*1.e9
	 if(localops > 0)then
         write(timefile,1002)'SPE',tpoSPEdefault,tpoSPE,tpoX
	 else
         write(timefile,1002)'SPE',tpoSPEdefault,tpoSPE,tpoSPEdefault
	 end if
		 

     print*,' - - - - - - - - - - - - - - - - - - - '
       if(call_spe)write(6,'(" Time per SPE operation = ",f7.2," ns")')tpoX


     if(.not.threebody)then
!............ WRITE OUT PP............
       localops = 0
       do ibundle = 1,nopbundles
             if(opbundle(ibundle)%optype /= 'PP')cycle
          if(opbundle(ibundle)%hchar == 'b')cycle          
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
       end do  ! ibundle
       tottimeX = 0.d0
       do ip = 0,nprocs-1
          tottimeX = tottimeX + time_procPP(ip)
        end do
       tpoX = tottimeX/(actual_iterations*localops)*1.e9
       if(localops==0)then
         write(timefile,1003)'PP ',tpoPPdefault, tpoPPprior, tpoPPdefault
       else
         write(timefile,1002)'PP ',tpoPPdefault,tpoPP,tpoX
         write(6,'(" Time per PP operation = ",f7.2," ns")')tpoX
       endif

!............ WRITE OUT PPb............
       localops = 0
       do ibundle = 1,nopbundles
          if(opbundle(ibundle)%optype /= 'PP')cycle
          if(opbundle(ibundle)%hchar /= 'b')cycle          
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
       end do  ! ibundle
       tottimeX = 0.d0
       do ip = 0,nprocs-1
          tottimeX = tottimeX + time_procPPb(ip)
        end do
       tpoX = tottimeX/(actual_iterations*localops)*1.e9
       if(localops==0)then
         write(timefile,1003)'PPb',tpoPPbdefault, tpoPPbprior, tpoPPbdefault
       else
         write(timefile,1002)'PPb',tpoPPbdefault,tpoPPb,tpoX
         write(6,'(" Time per PPb operation = ",f7.2," ns")')tpoX
       endif       
      
!............ WRITE OUT NN............
       localops = 0
       do ibundle = 1,nopbundles
             if(opbundle(ibundle)%optype /= 'NN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
       end do  ! ibundle
       tottimeX = 0.d0
       do ip = 0,nprocs-1
          tottimeX = tottimeX + time_procNN(ip)
        end do
       tpoX = tottimeX/(actual_iterations*localops)*1.e9
       if(localops==0)then
         write(timefile,1003)'NN ',tpoNNdefault, tpoNNprior, tpoNNdefault
       else
         write(timefile,1002)'NN ',tpoNNdefault,tpoNN,tpoX
         write(6,'(" Time per NN operation = ",f7.2," ns")')tpoX
       endif
!............ WRITE OUT PN............
       localops = 0
       do ibundle = 1,nopbundles
             if(opbundle(ibundle)%optype /= 'PN')cycle
             if(opbundle(ibundle)%hchar == 'b')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
       end do  ! ibundle
       tottimeX = 0.d0
       do ip = 0,nprocs-1
          tottimeX = tottimeX + time_procPN(ip)
        end do
       tpoX = tottimeX/(actual_iterations*localops)*1.e9
       if(localops==0)then
         write(timefile,1003)'PN ',tpoPNdefault, tpoPNprior, tpoPNdefault
       else
         write(timefile,1002)'PN ',tpoPNdefault,tpoPN,tpoX
         write(6,'(" Time per PN operation = ",f7.2," ns")')tpoX
       endif
!     print*,' time per PN operation = ',tpoX,' ns '

       !............ WRITE OUT PNb............
       localops = 0
       do ibundle = 1,nopbundles
          if(opbundle(ibundle)%optype /= 'PN')cycle
          if(opbundle(ibundle)%hchar /= 'b')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
       end do  ! ibundle
       tottimeX = 0.d0
       do ip = 0,nprocs-1
          tottimeX = tottimeX + time_procPNb(ip)
        end do
       tpoX = tottimeX/(actual_iterations*localops)*1.e9
       if(localops==0)then
         write(timefile,1003)'PNb',tpoPNbdefault, tpoPNbprior, tpoPNb
       else
         write(timefile,1002)'PNb',tpoPNbdefault,tpoPNb,tpoX
         write(6,'(" Time per PNb operation = ",f7.2," ns")')tpoX
       endif
       

!-------- 3-BODY TIMES UNCHANGED-------------

       write(timefile,1003)'PPP ',tpoPPPdefault, tpoPPPprior, tpoPPP
       write(timefile,1003)'NNN ',tpoNNNdefault, tpoNNNprior, tpoNNN
       write(timefile,1003)'PPN ',tpoPPNdefault, tpoPPNprior, tpoPPN
       write(timefile,1003)'PPN ',tpoPNNdefault, tpoPNNprior, tpoPNN

     else
!-------- 2-BODY TIMES UNCHANGED-------------

       write(timefile,1003)'PP ',tpoPPdefault, tpoPPprior, tpoPP
        write(timefile,1003)'PPb',tpoPPbdefault, tpoPPbprior, tpoPPb
       write(timefile,1003)'NN ',tpoNNdefault, tpoNNprior, tpoNN
       write(timefile,1003)'PN ',tpoPNdefault, tpoPNprior, tpoPN
       write(timefile,1003)'PNb',tpoPNbdefault, tpoPNbprior, tpoPNb
!............ WRITE OUT PPP............
       localops = 0
       do ibundle = 1,nopbundles
             if(opbundle(ibundle)%optype /= 'PPP')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
       end do  ! ibundle
       tottimeX = 0.d0
       do ip = 0,nprocs-1
          tottimeX = tottimeX + time_procPPP(ip)
        end do
       tpoX = tottimeX/(actual_iterations*localops)*1.e9
       if(localops==0)then
         write(timefile,1003)'PPP',tpoPPPdefault, tpoPPPprior, tpoPPP
       else
         write(timefile,1002)'PPP',tpoPPPdefault,tpoPPP,tpoX
         write(6,'(" Time per PPP operation = ",f7.2," ns")')tpoX
       endif
!............ WRITE OUT NNN............
       localops = 0
       do ibundle = 1,nopbundles
             if(opbundle(ibundle)%optype /= 'NNN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
       end do  ! ibundle
       tottimeX = 0.d0
       do ip = 0,nprocs-1
          tottimeX = tottimeX + time_procNNN(ip)
        end do
       tpoX = tottimeX/(actual_iterations*localops)*1.e9
       if(localops==0)then
         write(timefile,1003)'NNN',tpoNNNdefault, tpoNNNprior, tpoNNN
       else
       write(timefile,1002)'NNN',tpoNNNdefault,tpoNNN,tpoX
       write(6,'(" Time per NNN operation = ",f7.2," ns")')tpoX
       endif
!............ WRITE OUT PPN............
       localops = 0
       do ibundle = 1,nopbundles
             if(opbundle(ibundle)%optype /= 'PPN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
       end do  ! ibundle
       tottimeX = 0.d0
       do ip = 0,nprocs-1
          tottimeX = tottimeX + time_procPPN(ip)
        end do
       tpoX = tottimeX/(actual_iterations*localops)*1.e9
       if(localops==0)then
         write(timefile,1003)'PPN',tpoPPNdefault, tpoPPNprior, tpoPPN
       else
       write(timefile,1002)'PPN',tpoPPNdefault,tpoPPN,tpoX
       write(6,'(" Time per PPN operation = ",f7.2," ns")')tpoX
       endif
!............ WRITE OUT PNN............
       localops = 0
       do ibundle = 1,nopbundles
             if(opbundle(ibundle)%optype /= 'PNN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
       end do  ! ibundle
       tottimeX = 0.d0
       do ip = 0,nprocs-1
          tottimeX = tottimeX + time_procPNN(ip)
        end do
       tpoX = tottimeX/(actual_iterations*localops)*1.e9
       if(localops==0)then
         write(timefile,1003)'PNN',tpoPNNdefault, tpoPNNprior, tpoPNN
       else
       write(timefile,1002)'PNN',tpoPNNdefault,tpoPNN,tpoX
       write(6,'(" Time per PNN operation = ",f7.2," ns")')tpoX
       endif
  endif

      write(timefile,*)' '
      write(timefile,*)' All times in ns (nanoseconds) '
      write(timefile,*)' ' 
      write(timefile,*)' BIGSTICK Version ',version,lastmodified
      write(timefile,*)' '
      close(timefile)
     print*,' - - - - - - - - - - - - - - - - - - - '
  end if
  return

end subroutine timeperopmaster
!===================================================================
!===================================================================
! Timer: clock to measure elapsed time
! ====================================================================
subroutine bundle_clock(ibundle,starter)
  use timing_parallel
  use opbundles
  use nodeinfo
  use bmpi_mod
  implicit none

  integer(4)        :: ierr
  character(len=3)  :: starter
  real(8)           :: timenow
  character(len=12) :: real_clock(3)
  integer(kind=4)   :: time_values(8)
  real(kind=8)      :: wall_time                  !  Wall time in seconds
  real(kind=8)      :: smin, shour, sday
  integer           :: ibundle    ! dummy for which bundle
  integer           :: aerr

  smin = 60.0d0
  shour = smin*60.0d0
  sday = shour*24.0d0

!--------- THE FOLLOWING MIGHT BE PLATFORM DEPENDENT -----------------------
  if ( nproc == 1 ) then
!     call BMPI_BARRIER(icomm,ierr)
     call cpu_time(timenow)     ! sequential clock
     call date_and_time(real_clock(1), real_clock(2), real_clock(3),           &
                        time_values)
     wall_time = dfloat(time_values(8))/1000.0d0 + dfloat(time_values(7)) +    &
                 dfloat(time_values(6))*smin + dfloat(time_values(5))*shour +  &
                 dfloat((time_values(3)-1))*sday
     timenow = wall_time
  else
!     call BMPI_BARRIER(icomm,ierr)
     timenow = BMPI_Wtime()  ! MPI clock
     wall_time = BMPI_Wtime()
  end if

  select case (starter)

  case ('set')
     allocate( time_bundle(nopbundles), stat=aerr)
     if(aerr /= 0) call memerror("bundle_clock")
     time_bundle(:) = 0.d00
  case ('sta')
     timelastbundle = timenow

  case ('end')
     time_bundle(ibundle) = time_bundle(ibundle) + timenow-timelastbundle

  end select

  return
  end subroutine bundle_clock
!===============================================

  subroutine bundle_clock_out

  use timing_parallel
  use opbundles
  use nodeinfo
  use lanczos_info

  use bmpi_mod
  use para_bundles_mod
  
  implicit none

  integer iob
  real(8)  :: avg_time
  real(8)  :: max_time,jmpstor
  integer(8)  :: localops
  integer ierr

  if (iproc == 0) then
      call MPI_Reduce(MPI_IN_PLACE, time_bundle, nopbundles, MPI_REAL8, MPI_SUM, 0, icomm, ierr)
  else
      call MPI_Reduce(time_bundle,  time_bundle, nopbundles, MPI_REAL8, MPI_SUM, 0, icomm, ierr)
  endif

  if(iproc == 0)then
    open(unit=49,file='bundletimes.bigstick',status='unknown')
	write(49,*)' bundle/type/dir/       # ops / tot time(s)    time/operation(ns)  '
    max_time = 0.0
    avg_time = 0.0
    do iob=1,nopbundles
		if(nproc == 1)call analyze_opbundle(iob,.false.,jmpstor)
      localops =  int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)* & 
           int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8)
      write(49,'(i8,2a4,1i12,g20.7,2f20.4)')iob,  & 
            opbundle(iob)%optype,opbundle(iob)%hchar, & 
!           int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8), & 
!           int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8), &
           localops,  time_bundle(iob),  &  !opbundle(iob)%njumps, &
           time_bundle(iob)*1.0e9/dfloat(localops*actual_iterations) !, &
!		   time_bundle(iob)
       !    , time_bundle(iob)*1.0e9/dfloat(actual_iterations)
    end do

    print*,' (Individual bundle data written to file bundletimes.bigstick )'
    close(49)
  end if
  return

  end subroutine bundle_clock_out
!=====================================================
!
! added in 7.7.9, routine to summarize all output timing
!
! replaces all or in part the following routines:
! proc_clock_out
! procOp_clock_out
! bundle_clock_out
!

subroutine all_clocks_out

  use timing_parallel
  use nodeinfo
  use bmpi_mod
  use butil_mod
  use io
  use flags3body
  use operation_stats
  use program_info
  
  use opbundles
  use lanczos_info

  use para_bundles_mod
  use localvectors,only:useHZSomp
  use sectors
  implicit none

  integer iprocs
  real(8)  :: avg_time
  real(8)  :: max_time
  integer tag
  integer(4) :: ierr
  integer stat(MPI_STATUS_SIZE)
  integer iob
  real(8)  :: jmpstor
  integer(8)  :: localops
  integer(8) :: length_ini, length_fin  ! length of initial, final sectors (p+n)
  integer  :: ips,fps,ins,fns
  integer  :: ibundle

  real(8) :: tmp(10)
  real(8) ::expectedtime
  character(10) :: date,time,zone
  integer :: datetimeval(8)
  integer omp_get_num_threads
  
  integer :: num_threads
  
  !$omp parallel shared(num_threads,iproc)  
     num_threads = omp_get_num_threads()
  !$omp end parallel
  
!------ EVERY PROC SEND THEIR DATA TO PROC ZERO----
  if(iproc >0)then
    tag = iproc
    call BMPI_Send(time_Ham_MPI(iproc),1,0,tag,icomm,ierr) 
  else
    do iprocs = 0,nprocs-1
        if(iprocs > 0)then
          tag = iprocs
          call BMPI_Recv(time_Ham_MPI(iprocs),1,iprocs,tag,icomm,stat,  ierr) 
        end if
	end do
  end if
  
call BMPI_BARRIER(icomm,ierr) 

if(iproc >0)then
  tag = iproc
   tmp(1) = time_procSPE(iproc)
   tmp(2) = time_procPP(iproc)
   tmp(3) = time_procNN(iproc)
   tmp(4) = time_procPN(iproc)
   
   tmp(5) = time_procPPb(iproc)
   tmp(6) = time_procPNb(iproc)

   tmp(7) = time_procPPP(iproc)
   tmp(8) = time_procPPN(iproc)
   tmp(9) = time_procPNN(iproc)
   tmp(10) = time_procNNN(iproc)
   
   call BMPI_SEND(tmp,10,0,tag,icomm,ierr) 
end if

if(iproc == 0)then
  do iprocs = 1,nprocs-1
      tag = iprocs
      call BMPI_Recv(tmp,10,iprocs,tag,icomm,stat,ierr)
      time_procSPE(iprocs) = tmp(1)

      if(.not.threebody)then
         time_procPP(iprocs)  = tmp(2)
         time_procNN(iprocs)  = tmp(3)
         time_procPN(iprocs)  = tmp(4)

         time_procPPb(iprocs) = tmp(5)
         time_procPNb(iprocs) = tmp(6)
      else
         time_procPPP(iprocs) = tmp(7)
         time_procPPN(iprocs) = tmp(8)
         time_procPNN(iprocs) = tmp(9)
         time_procNNN(iprocs) = tmp(10)
      endif
     
  end do
end if
!----------- BUNDLE DATA----------------


if (iproc == 0) then
    call MPI_Reduce(MPI_IN_PLACE, time_bundle, nopbundles, MPI_REAL8, MPI_SUM, 0, icomm, ierr)
else
    call MPI_Reduce(time_bundle,  time_bundle, nopbundles, MPI_REAL8, MPI_SUM, 0, icomm, ierr)
endif

!--------- NOW WRITE TO A FILE -----------------------------
    
if(iproc == 0)then
    open(unit=91,file='timingdata.bigstick',status='unknown')
    write(91,*)' BIGSTICK Version ',version,lastmodified
    call date_and_time(date,time,zone,datetimeval)
    write(91,*)' Run date: ',date(1:4),'-',date(5:6),'-',date(7:8)
	write(91,'(i5,"  iterations total")')actual_iterations
!....... WRITE A HEADER TO proctimes.bigstick........
    write(91,*)'# timing weights assumed '
	if(.not.threebody)then
		write(91,1002)'SPE',opwtSPE
		write(91,1002)'PP ',opwtPP
		write(91,1002)'PPb',opwtPPb
		write(91,1002)'NN ',opwtNN
		write(91,1002)'PN ',opwtPN
		write(91,1002)'PNb',opwtPNb		
	else
		write(91,1002)'SPE',opwtSPE
		write(91,1002)'PPP',opwtPPP
		write(91,1002)'NNN',opwtNNN
		write(91,1002)'PPN',opwtPPN
		write(91,1002)'PNN',opwtPNN			
	endif 
!    write(91,*)'#  proc    actual time     '
!....... END OF HEADER................	
    max_time = 0.0
    avg_time = 0.0
    do iprocs = 0,nprocs-1
      avg_time = avg_time + time_Ham_MPI(iprocs)
      max_time = bmax(max_time, time_Ham_MPI(iprocs))
      
    end do

    avg_time =avg_time / float(nprocs)
    print*,' Average time per processor is ',avg_time, " max ", max_time
    print*,' Overall efficiency is ',avg_time/max_time
    print*,' (Individual processor data written to file proctimes.bigstick )'
    write(logfile,*)' Average time per processor is ',avg_time, " max ", max_time
    write(logfile,*)' Overall efficiency is ',avg_time/max_time
  1002 format(1x,a3,2x,f7.2)
  write(91,*)' useHZSomp = ',useHZSomp
  write(91,*)num_threads,' OpenMP threads'

!----------- NOW LOOP OVER MPI PROCESSES ----------------------
  do iprocs = 0,nprocs-1
	 expectedtime = 0.d0
!     write(91,*)' processor, actual time ',iprocs, time_Ham_MPI(iprocs)
     write(91,'("MPI process   :",i7)')iprocs
!	 write(91,'("Total time (s):",f12.3)')time_Ham_MPI(iprocs)
 	 write(91,*)'#  optype  total time (s)  time per op (ns) expected total time (s)'	

     localops = 0
     do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'SPE')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
     end do
     if(localops ==0)then
          write(91,999)'SPE',0.0
     else
          write(91,999)'SPE',time_procSPE(iprocs),time_procSPE(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtSPE*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtSPE
     end if
!998  format(a3,2x,1f12.3)
999  format(a3,2x,f12.3,5x,f12.3,5x,f12.3 )

     if(.not.threebody)then
        localops = 0
        do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
           if(opbundle(ibundle)%optype /= 'PP')cycle
           if(opbundle(ibundle)%hchar == 'b')cycle
           localops = localops + &
             int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
             int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
           write(91,999)'PP ',0.0
        else
          write(91,999)'PP ',time_procPP(iprocs),time_procPP(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPP  *1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPP			 
        end if

!added by shan        
        localops = 0
        do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
            if(opbundle(ibundle)%optype /= 'PP')cycle
            if(opbundle(ibundle)%hchar /= 'b')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
           write(91,999)'PPb',0.0
        else
           write(91,999)'PPb',time_procPPb(iprocs),time_procPPb(iprocs)/dfloat(actual_iterations*localops)*1.e9 , & 
	             dfloat(localops*actual_iterations)*opwtPPb *1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPPb
        end if
  
        localops = 0
        do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'NN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
           write(91,999)'NN ',0.0
        else
           write(91,999)'NN ',time_procNN(iprocs),time_procNN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtNN*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtNN
        end if

        localops = 0
        do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'PN')cycle
            if(opbundle(ibundle)%hchar == 'b')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
           write(91,999)'PN ',0.0
        else
          write(91,999)'PN ',time_procPN(iprocs),time_procPN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPN*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPN
        endif

 !added by shan:
        localops = 0
        do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
            if(opbundle(ibundle)%optype /= 'PN')cycle
            if(opbundle(ibundle)%hchar /= 'b')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
           write(91,999)'PNb',0.0
        else
           write(91,999)'PNb',time_procPNb(iprocs),time_procPNb(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPNb*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPNb
        endif
        
     else  ! 3-body
        localops = 0
        do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'PPP')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
           write(91,999)'PPP',0.0
        else
           write(91,*)'PPP',time_procPPP(iprocs),time_procPPP(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPPP *1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPPP
        end if
        localops = 0
        do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'NNN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
           write(91,999)'NNN',0.0
        else
           write(91,999)'NNN',time_procNNN(iprocs),time_procNNN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtNNN*1.e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtNNN
        end if

        localops = 0
        do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'PPN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
           write(91,999)'PPN',0.0
        else
          write(91,999)'PPN',time_procPPN(iprocs),time_procPPN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPPN*1e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPPN
        endif

        localops = 0
        do ibundle = opbundlestart(iprocs), opbundleend(iprocs)
             if(opbundle(ibundle)%optype /= 'PNN')cycle
               localops = localops + &
           int( opbundle(ibundle)%pxend -opbundle(ibundle)%pxstart+1,8)* & 
           int( opbundle(ibundle)%nxend -opbundle(ibundle)%nxstart+1,8)
        end do
        if(localops ==0)then
           write(91,999)'PNN',0.0
        else
           write(91,999)'PNN',time_procPNN(iprocs),time_procPNN(iprocs)/dfloat(actual_iterations*localops)*1.e9, & 
	             dfloat(localops*actual_iterations)*opwtPNN*1e-9
				 expectedtime =expectedtime+  dfloat(localops*actual_iterations)*opwtPNN
        endif

     end if
	 write(91,'(" Total time (s): Actual ",f12.3,"   ; expected ",f12.3)')time_Ham_MPI(iprocs),expectedtime*1e-9

!     write(91,'(" Expected total time (s): ",f12.3)')expectedtime*1e-9
!	 write(91,'(" Total time (s):          ",f12.3)')time_Ham_MPI(iprocs)
	 
!.......... PRINT OUT BUNDLE INFO..............
     write(91,'("# bundle /optype  # ops  # p-jumps  # n-jumps  initial basis,  final basis  total time(s) time/op (ns)	")')
!	 '# bundle /optype     # ops    # p-jumps   # n-jumps  initial basis,  final basis     total time(s) time/op (ns)'	 
     do iob=opbundlestart(iprocs), opbundleend(iprocs)
 		if(nproc == 1)call analyze_opbundle(iob,.false.,jmpstor)
		
		call bundle_basis(iob,length_ini,length_fin)

!		length_ini = xsd(1)%sector(ips)
!		length_ini = int(xsd(1)%sector( ips)%nxsd,8)* int(xsd(2)%sector( ins)%nxsd,8)
!		length_fin = int(xsd(1)%sector( fps)%nxsd,8)* int(xsd(2)%sector( fns)%nxsd,8)

       localops =  int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)* & 
            int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8)
       write(91,'(i8,x,a3,a1,i14,4i12,f20.4,f11.3)')iob,  & 
             opbundle(iob)%optype,opbundle(iob)%hchar, & 
            localops, int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8),  & 
			int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8), & 
			length_ini,length_fin,time_bundle(iob),&  
            time_bundle(iob)*1.0e9/dfloat(localops*actual_iterations) 
!			write(92,*)iob,ips,ins,fps,fns,xsd(1)%sector( ips)%nxsd,xsd(2)%sector( ins)%nxsd
     end do

    end do  !iprocs

    close(91)
	print*,' '
	print*,' Detailed data on timing found in file timingdata.bigstick '
  end if
!  call BMPI_BARRIER(icomm,ierr) 

  return
	
	return
end subroutine all_clocks_out
!=====================================================
!
! added in 7.3.8
!
! creates a report on memory usage by MPI process
!
!
  subroutine memreport

  use nodeinfo
  use precisions
  use fragments
!  use tribution
  use mod_reorthog

  implicit none
  integer :: jproc

  integer fragi,fragf
  real  :: fragmem  ! memory for fragments
  real  :: lancmem  ! memory for lanczos storage
  real  :: jumpmem  ! memory for jumps

  if(iproc/= 0)return
  open(unit=64,file='procmemload.bigstick',status='unknown')
  write(64,*)' All memory in Gb '
  write(64,*)' Proc     Lanczos    Vec storage  ' 
  do jproc= 0, nproc-1
! . . . . . . . . . COMPUTE FRAGMENT MEMORY .......     
   fragi = nodal(jproc)%ifragment
   fragf = nodal(jproc)%ffragment
   fragmem = real(basestop(fragi)+basestop(fragf) - basestart(fragi)-basestart(fragf)+2)
   fragmem = fragmem*lanc_prec*1.0e-9    ! in Gb

!. . . . . . . . . .COMPUTE LANCZOS STORAGE MEMORY. . . . . . .
   lancmem = real(br_ostoplist(jproc) -br_ostartlist(jproc)+1)
   lancmem = lancmem*lanc_prec*1.0e-9*br_histmax 

! ; ; ; ; ; ; ; ; ; WRITE OUT TO FILE ; ; ; ; ; ; ; ; ; ; ; ; ; ; ;

    write(64,111)jproc,fragmem,lancmem

  end do   ! jproc
111 format(i4,5x,2f10.5)
  print*,' A report on memory usage can be found in file procmemload.bigstick '
  close(unit=64)
  return
  end subroutine memreport
