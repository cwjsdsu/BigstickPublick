!=================================================================
!
! routines to report on distributon
!
!
!===========================================================

module para_report_mod
	use para_bundles_mod
	use para_util_mod
contains

!===============================================================
!
! COMPUTES THE NUMBER OF OPERATIONS AS WELL THE JUMPS OVER WHICH THEY ARE DISTRIBUTED
!
!  ADDED 7.7.9
!   replaces print_ops_distro, print_jump_distro, print_opbundle_distro
!   NEW OUTPUT:  distrodata.bigstick 
!       REPLACES jumpdistro.bigstick,procjumpdistro.bigstick, opsdistro.bigstick,  
!
!  CALLED BY: 
!     master_para_distribute
!
!  SUBROUTINES CALLED:
!   count_jumps_on_node
!
subroutine print_all_distro
   use nodeinfo
   use opbundles
   use operation_stats
   use butil_mod
   use program_info
   use fragments
   use jumplimits,only:maxjumpmemory
   use io,only:logfile
   use menu_choices,only:menu_char
   
   implicit none

   integer :: iob,jnode
   integer(8) :: localops,maxlocalops,localjumps
   integer(8) :: startz,stopz,laststart
   integer :: alldistrofile=97
   real    :: localwt
   real(8) :: maxwtops, totwtops,wtlocalops
   character(10) :: date,time,zone
   integer :: datetimeval(8)
   integer(8) :: sumjumps,maxlocaljumps,avgjumps
   integer :: maxproc
   integer :: bytesPerJump = 21
   integer :: nlocalgreedy
   real(8) :: localgreedymemory,localjumpstorage
   integer :: ibundle

   if(iproc /= 0)return

   open(unit=alldistrofile,file='distrodata.bigstick',status='unknown')
   write(alldistrofile,*)' BIGSTICK Version ',version,lastmodified
   call date_and_time(date,time,zone,datetimeval)
   write(alldistrofile,*)' Run date: ',date(1:4),'-',date(5:6),'-',date(7:8)
   write(alldistrofile,*)'# bundles:    proton    x    neutron       total jumps'
   maxlocalops = 0
   maxwtops = 0.0
   totwtops = 0.0 
   avgjumps = 0
   maxlocaljumps = 0
   maxjumpmemory=0.d0
   maxproc = -1
   do jnode = 0,nprocs-1	
!	   print*,' JNODE ',jnode 
      localops = 0
      if(opbundleend(jnode) < 1)cycle
      wtlocalops = 0.0
      do iob = opbundlestart(jnode),opbundleend(jnode)
!		  if(iproc==0 .and. opbundle(iob)%greedy)print*,' greedy! ',opbundle(iob)%jumpstore
!		  print*,' IOB ',iob
         if(usewts)then
			 localwt = opbundle(iob)%opwt
!             select case (opbundle(iob)%optype)
!                case('SPE')
!                   localwt = opwtSPE
!                case('PP ')
!                    localwt = opwtPP
!                   if (opbundle(iob)%hchar == 'b') localwt = opwtPPb
!                case('NN ')
!                    localwt = opwtNN
!                case('PN ')
!                    localwt = opwtPN
!                   if (opbundle(iob)%hchar == 'b') localwt = opwtPNb
!                case('PPP')
!                    localwt = opwtPPP
!                case('NNN')
!                    localwt = opwtNNN
!                case('PPN')
!                    localwt = opwtPPN
!                    if (opbundle(iob)%hchar == 'b') localwt = opwtPPNb					
!                case('PNN')
!                    localwt = opwtPNN
!                    if (opbundle(iob)%hchar == 'b') localwt = opwtPNNb					
!                case default
!                    if(iproc == 0)print*,' OOPS wrong optype ',draft_opbundle(iob)%optype
!                    stop
!             end select
         else
             localwt = 1.d0
         end if
         select case ( opbundle(iob)%optype )
            case('PP ', 'PPP','PPN', 'SPE')
               startz = opbundle(iob)%pxstart
               stopz =  opbundle(iob)%pxend

            case('NN ','PN ','NNN','PNN')
               startz = opbundle(iob)%nxstart
               stopz =  opbundle(iob)%nxend

            case default
               startz = 0 ! prevent compiler warnings
               stopz = 0
         end select
         localops = localops + int(stopz-startz+1,8)*int(opbundle(iob)%min_nop,8)
         wtlocalops = wtlocalops+localwt*real(stopz-startz+1,8)*real(opbundle(iob)%min_nop,8)
         totwtops = totwtops + localwt*real(stopz-startz+1,8)*real(opbundle(iob)%min_nop,8)
      end do
      maxlocalops = bmax(maxlocalops, localops)
      maxwtops = bmax(maxwtops,wtlocalops)
      write(alldistrofile,*)' processor ',jnode !,' total operations ',localops,', time (ns) for ops ',wtlocalops
  
      write(alldistrofile,*)' initial, final fragments : ',nodal(jnode)%ifragment, nodal(jnode)%ffragment

      do iob = opbundlestart(jnode),opbundleend(jnode)

         select case ( opbundle(iob)%optype )
             case('SPE')
                localjumps = 0
             case('PP','PPP')
                localjumps = int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)
             case('NN','NNN')
                localjumps = int( opbundle(iob)%nsortend -opbundle(iob)%nsortstart+1,8)
             case('PN','PPN')
                localjumps = int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8)    & 
               + int( opbundle(iob)%nsortend -opbundle(iob)%nsortstart+1,8)

         end select


		 write(alldistrofile,999)iob,opbundle(iob)%optype,opbundle(iob)%hchar, & 
         int( opbundle(iob)%pxend -opbundle(iob)%pxstart+1,8), & 
         int( opbundle(iob)%nxend -opbundle(iob)%nxstart+1,8),  & 
		  localjumps
      end do
	  
!--------------- JUMPS STORAGE --------------------------------	  
      sumjumps = 0
      call count_jumps_on_node(jnode,0,sumjumps,localjumpstorage)
      write(alldistrofile,*)' Node ',jnode,', jump storage = ', &
         1.0e-6*localjumpstorage,' Mb '!
!         if(sumjumps > maxlocaljumps)maxproc = jnode
      if(localjumpstorage > maxjumpmemory)maxproc = jnode
      maxlocaljumps = bmax(maxlocaljumps,sumjumps)
	  maxjumpmemory=bmax(maxjumpmemory,localjumpstorage)
      avgjumps = avgjumps + sumjumps	
999 format(i8,2x,a3,a1,3i15)	  


   end do   ! jnode

   totwtops = totwtops/float(nprocs)
!   maxjumpmemory = 1.0e-6*bytesPerJump*maxlocaljumps	
   avgjumps = int(avgjumps/real(nprocs,kind=8),8)
   
   if(nproc > 1 .or. menu_char=='m')then
!      print*,' '
!      print*,' Countup of total weighted ops = ',totwtops*float(nprocs)
!	  print*,' '
      write(logfile,*)' Estimated efficiency of MPI distribution is approx ', & 
              totwtops/maxwtops
      write(logfile,*)' (Detailed information on distribution across MPI in distrodata.bigstick)'
!      print*,' '

      write(6,*)' max local storage of jumps ~', &
       1.0e-6*maxjumpmemory,' Mb, on process ',maxproc
       write(logfile,*)' max local storage of jumps ~', &
        1.0e-6*maxjumpmemory,' Mb, on process ',maxproc 
!       1.0e-6*bytesPerJump*maxlocaljumps,' Mb, on process ',maxproc

!............. INVESTIGATE GREEDY OPBUNDLES.....................

       if(ngreedy > 0)then
		   nlocalgreedy=0
		   localgreedymemory=0.
		   do ibundle=opbundlestart(maxproc),opbundleend(maxproc)
			   if( opbundle(ibundle)%greedy)then
				   nlocalgreedy=nlocalgreedy+1
				   call analyze_opbundle(ibundle,.false.,localjumpstorage)
				   localgreedymemory=localgreedymemory+ opbundle(ibundle)%jumpstore
				   
			
			   end if
			   
			   
		   end do
		   write(logfile,*)' This MPI rank has ',nlocalgreedy,' greedy opbundles requiring ',localgreedymemory*1.0e-6,' Mb memory' 
		   write(logfile,*)' (Note, however, this memory may overcount )'
		   
		   
	   end if  
	   
	   
      write(6,*)' avg storage of jumps is ',1.0e-6*bytesPerJump*avgjumps,' Mb '

      write(logfile,*)' max local storage of jumps = ', &
       1.0e-6*maxjumpmemory,' Mb, on process ',maxproc
!       1.0e-6*bytesPerJump*maxlocaljumps,' Mb, on process ',maxproc

      write(logfile,*)' avg storage of jumps is ',1.0e-6*bytesPerJump*avgjumps,' Mb '

   endif
	
   
   close(alldistrofile)
   return
 end subroutine print_all_distro 

!========================================================

end module para_report_mod
