!=================================================================
!
! file BPARALLEL_OPBUNDLES.F90
!
! routines to compute MPI distribution for multiple nodes; 
! originally only for one "fragment" (that is, lanczos vectors not broken up)
! all-in-one routine to be broken up in 7.2.13/14...
!
! started 4/2012 by CWJ
! original routines used when Lanczos vectors are NOT broken into fragments
! (supplanting routines handling fragments are found in bparallel_lib2.f90)
!
!===========================================================

module para_bundles_mod
	use opbundles
	
	implicit none
	
	real(kind=8):: excess_memory_needed
	
contains

!=================================================================
!
! semi-master routine for creating draft opbundles
!
!
! CALLED BY:
!    distro_opbundles_over_fragments (above)
!
! SUBROUTINES CALLED:
!    count_create_draft_opbundles_g
!    analyze_opbundle

subroutine draft_opbundles(nob_draft)

   use fragments
   use verbosity
   use nodeinfo
   use operation_stats
   use coupledmatrixelements,only:call_spe,dens2bflag
   use flagger
   use jumpstart
   use flags3body,only:applycentroids
   use io
   implicit none
 
   integer ifrag,ffrag
   integer nob,nob_draft
   real(4) :: localwt
   integer(8) :: localjumps,maxlocaljumps
   integer iprocs
   integer :: nnodesoverbooked, nnodesexcess
   logical :: print_bundle_stats=.true.
   integer :: aerr
   integer(8) :: allmyops
   real(8) :: opjumpstorage,maxopjumpstorage

   if(allocated(frag_draftopbundlestart))deallocate(frag_draftopbundlestart)
   allocate(frag_draftopbundlestart(nfragments,nfragments), stat=aerr )
   if(aerr /= 0) call memerror("draft_opbundles 1")
   
   if(allocated(frag_draftopbundleend))deallocate(frag_draftopbundleend)
   allocate(frag_draftopbundleend(nfragments,nfragments), stat=aerr )
   if(aerr /= 0) call memerror("draft_opbundles 2")


   nob_draft = 0
   do ifrag = 1,nfragments
       do ffrag = 1,nfragments
           if(call_spe)call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'SPE','h',nob_draft)
           if(applycentroids .and. ifrag==ffrag)call & 
		      count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'CEN','h',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PP0','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'NN0','f',nob_draft)
            call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PN0','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PN0','h',nob_draft)

           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PPP','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'NNN','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PPN','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PNN','f',nob_draft)
!           if(ifrag <= ffrag .and. .not.dens2bflag)then
		   
           if(ifrag <= ffrag)then
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.false.,'PP0','b',nob_draft)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.false.,'NN0','b',nob_draft)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.false.,'PN0','b',nob_draft)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.false.,'PPP','b',nob_draft)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.false.,'NNN','b',nob_draft)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.false.,'PPN','b',nob_draft)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.false.,'PNN','b',nob_draft)
           end if

       end do

   end do

   if(verbose_bundles .and. iproc == 0)then
     print*,nob_draft, ' opbundles  in first draft'
   end if
   if(allocated(draft_opbundle))deallocate(draft_opbundle)
   allocate(draft_opbundle(nob_draft), stat=aerr )
   if(aerr /= 0) call memerror("draft_opbundles 10")
   nob = 0
   do ifrag = 1,nfragments
       do ffrag = 1,nfragments
           frag_draftopbundlestart(ifrag,ffrag) = nob+1
           if(call_spe)call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'SPE','h',nob)
           if(applycentroids .and. ifrag==ffrag)call & 
		      count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'CEN','h',nob)
		   
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PP0','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'NN0','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PN0','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PN0','h',nob)

           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PPP','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'NNN','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PPN','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PNN','f',nob)
!           if(ifrag <= ffrag .and. .not. dens2bflag)then
	       if(ifrag <= ffrag )then
                  call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.true.,'PP0','b',nob)
                  call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.true.,'NN0','b',nob)
				  call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.true.,'PN0','b',nob)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.true.,'PPP','b',nob)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.true.,'NNN','b',nob)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.true.,'PPN','b',nob)
              call count_create_draft_opbundles_g(.true.,ffrag,ifrag,.true.,'PNN','b',nob)
           end if
           frag_draftopbundleend(ifrag,ffrag) = nob
       end do

   end do
   if(allocated(frag2fragops))deallocate(frag2fragops)
   allocate( frag2fragops(nfragments,nfragments), stat=aerr )
   if(aerr /= 0) call memerror("draft_opbundles 20")
   totalops = 0.d0

   frag2fragops(:,:) = 0.d0
   draft_opbundle(:)%min_nop=0
   draft_opbundle(:)%nsetops=0
   if(print_bundle_stats .and. iproc==0)then
	   open(unit=87,file='nnopsvsjumps.bigstick',status='unknown')
	   
	   open(unit=88,file='ppopsvsjumps.bigstick',status='unknown')
	   open(unit=89,file='pnopsvsjumps.bigstick',status='unknown')
   end if
   maxopjumpstorage=0.d0
   
   !........ ADDED in 7.9.6: can read in 'jumpstart' fine-grained timing data...
	 
   	 if(usewts .and. readjumpstart)then
        call readtimeinfo4jumpstart(nob_draft) ! read_opwt4jumpstart(nob_draft)
   	 end if
   
   excess_memory_needed = 0.0
   do nob = 1,nob_draft
     call analyze_opbundle(nob,.true.,opjumpstorage)
	 maxopjumpstorage=max(maxopjumpstorage,opjumpstorage)



     if(usewts .and. .not. readjumpstart)then
       select case (draft_opbundle(nob)%optype)
       case('SPE')
          localwt = opwtSPE
       case('PP ')
          localwt = opwtPP
          if (draft_opbundle(nob)%hchar == 'b') localwt = opwtPPb
	 	 if(print_bundle_stats .and. iproc==0)write(88,*)draft_opbundle(nob)%njumps,draft_opbundle(nob)%nops
		  
       case('NN ')
          localwt = opwtNN
 	 	 if(print_bundle_stats .and. iproc==0)write(87,*)draft_opbundle(nob)%njumps,draft_opbundle(nob)%nops
		  
       case('PN ')
          localwt = opwtPN
          if (draft_opbundle(nob)%hchar == 'b') localwt = opwtPNb
 	 	 if(print_bundle_stats .and. iproc==0)write(89,*) draft_opbundle(nob)%njumps,draft_opbundle(nob)%nops
		  
       case('PPP')
          localwt = opwtPPP
       case('NNN')
          localwt = opwtNNN
       case('PPN')
          localwt = opwtPPN
          if (draft_opbundle(nob)%hchar == 'b') localwt = opwtPPNb
		  
       case('PNN')
          localwt = opwtPNN
          if (draft_opbundle(nob)%hchar == 'b') localwt = opwtPNNb
		  
		  case('cen','CEN')
		  localwt = 1.	  
		  
       case default
          if(iproc == 0)print*,' OOPS ',nob,draft_opbundle(nob)%optype
          stop
       end select
  	   draft_opbundle(nob)%opwt=localwt
	   
     else
       localwt = 1.d0
     end if
	 if(readjumpstart)localwt=draft_opbundle(nob)%opwt
	 	 
     ifrag = draft_opbundle(nob)%ifragment
     ffrag = draft_opbundle(nob)%ffragment
!	 write(78,*)nob, draft_opbundle(nob)%nops,localwt
     frag2fragops(ifrag,ffrag) = frag2fragops(ifrag,ffrag) + dfloat(draft_opbundle(nob)%nops)*localwt
!	 if(iproc==0)print*,draft_opbundle(nob)%nops, localwt, ' weighting '
     totalops = totalops + dfloat(draft_opbundle(nob)%nops)*localwt
	 	 
   end do   
   call count_all_operations(.true.,nob_draft,allmyops)
   
   if(iproc==0 .and. excess_memory_needed > 0.0)then
	   print*,'  !!!!!! '
	   print*,' May need to increase memory for jumpbs by ',excess_memory_needed,' Gb or approx ', &
	   excess_memory_needed/maxjumpmemory_default,' MPI nodes '
	   print*,'  !!!!!! '
	   
	   write(logfile,*)'  !!!!!! '
	   write(logfile,*)' May need to increase memory for jumpbs by ',excess_memory_needed,' Gb or approx ', &
	   excess_memory_needed/maxjumpmemory_default,' MPI nodes '
	   write(logfile,*)'  !!!!!! '
	   
   end if

   nopbundles = nob_draft
   avgopstore = real(nprocs * maxjumpmemory_default*1e9,kind=8)/ real(allmyops,kind=8)
   avgwtopstore = real(nprocs * maxjumpmemory_default*1e9,kind=8)/ real(totalops,kind=8)

!........... LOOK FOR GREEDY BUNDLES AND RECOMPUTE WEIGHTS............
!             added 7.9.12
!      
   call greedy_bundles(nob_draft)
   call count_all_operations_new(.true.,nob_draft,allmyops)
!   call greedy_bundles(nob_draft)     ! re-iterate to balance better
!   call count_all_operations_new(.true.,nob_draft,allmyops)   

   if(iproc==0)then
	   if(maxopjumpstorage > 0.0)print*,' Max storage for jumps = ',maxopjumpstorage*1.0e-9,' Gb '
         if(usewts)then
              print*,' Estimated time per mat-vec multiply is ',totalops,' ns'
         else 
              print*,totalops,' operations  '
         end if
		 nopbundles = nob_draft
!		 call count_all_operations(.true.,allmyops)
		 print*,' Total operations = ',allmyops
   end if    

   call print_fragops  ! only for data gathering
   return

end subroutine draft_opbundles
!=====================================================================


! creates a draft set of bundles
!
! this is really the same routine as count_create_opbundles
! found in bparallel_lib2.f90 and the two can/should be merged
!
! Modified so as to group opbundles together by optype and by hchar
!
! REVISED 7.9.0 to loop over sectorjumps, not just search for one
!  this allows for multiple opbundles to be created, from split sectorjumps
!
! INPUT:
!   draft:  logical flag; if T then fill draft_opbundle, else fill opbundle
!   ifrag, ffrag: initial, final fragments of lanczos vectors
!   create: logical flag, if T then fill, else just count
!   optype: controls operation: PP, NN, PN, etc.
!   hchar : Hermitian flow control: F, B, or H
!
! OUTPUT
!   nob   : # of opbundles
!
! CALLED BY:
!   distro_opbundles_over_fragments
!   distro_opbundles_1frag_g
!   draft_opbundles
!
subroutine count_create_draft_opbundles_g(draft,ifrag,ffrag,create,optype,hchar,nob)
  use fragments
  use sectors
  use operation_stats
  use jumpNbody
  use jump3body
  use basis
  use nodeinfo
  use butil_mod
  use flagger
  use coupledmatrixelements,only:dens2bflag
  implicit none

  logical, intent(in) :: draft
  logical, intent(in) :: create
  integer, intent(inout) :: nob  !   temporary for # of opbundles
  integer, intent(in) :: ifrag,ffrag
  character(3), intent(in) :: optype
  character(1), intent(in) :: hchar

  integer :: is,fs  ! indices for sectors
  integer(8):: xstartx,xendx,ystarty,yendy
  integer :: xsj, ysj, isj
  integer :: cs,ncs,ics,csi,csf
  integer :: cstride
  integer :: insector,fnsector
  logical okay
  integer(8) :: localops
  integer(kind=basis_prec) :: ibasis

  type (bund), pointer :: bundle(:)
  type (bund), target :: bundummy(1)
	  
  integer(8) :: maxjumpsallowed,localjumps
  logical :: foundjumps
	  

  bundle => bundummy  ! prevent compiler warnings
  if(draft .and. create)bundle => draft_opbundle
  if(.not.draft .and. create)bundle => opbundle
  
  maxjumpsallowed = (maxjumpmemory_default*maxfrac_sectjump_store * 1e9)/bytesper2Bjump
  

  do fs = fragmentlist(ffrag)%ssectorstart, fragmentlist(ffrag)%ssectorend   ! order switched in 7.6.1
     do is = fragmentlist(ifrag)%ssectorstart, fragmentlist(ifrag)%ssectorend

!         if(is < fs .and. .not. dens2bflag)cycle   ! to prevent double counting; and this is the way
          if(is < fs)cycle   ! to prevent double counting; and this is the way
                                                   ! sector jumps are organized
												   ! modified 7.9.2 for two-body densities

         select case(optype)

!........... ASSIGN PP .......................

         case('PP0')
         if(opstat(is,fs)%nopPP > 0)then

!..........FIND START, STOP FOR PP JUMPS..........
!          ASSUME SUBSECTORS SAME AS SECTORS
               xsj = -1
               do isj = 1,x2bjump(1)%nsectjumps
                   if( is == x2bjump(1)%isector(isj) .and.  & 
                        fs == x2bjump(1)%fsector(isj)  ) then
                      xsj = isj

                      xstartx = x2bjump(1)%sjmp(xsj)%nstart+1
                      xendx   = xstartx + x2bjump(1)%sjmp(xsj)%njumps-1

!..........FIND START, STOP FOR NEUTRON SDs.........
                      csi = x2bjump(1)%csjmp(xsj)%cjump(1)  ! this is the first conjugate
                                                     ! neutron sector
                      ystarty = xsd(2)%sector(csi)%xsdstart
                      ncs = x2bjump(1)%csjmp(xsj)%ncjmps
                      csf  = x2bjump(1)%csjmp(xsj)%cjump(ncs)  
                      yendy = xsd(2)%sector(csf)%xsdend
                      localops = int(yendy-ystarty+1,8)*int(xendx-xstartx+1,8)
!..........CONVERT TO STATE INDICES

                      cstride = 1
                      nob = nob+1
                      if(create) then  ! set up for PP operations
                          bundle(nob)%optype     = 'PP '
                          bundle(nob)%hchar      =  hchar
                          bundle(nob)%pxstart    = xstartx
                          bundle(nob)%pxend      = xendx
                          bundle(nob)%nxstart    = ystarty
                          bundle(nob)%nxend      = yendy
                          bundle(nob)%cstride    = cstride
                          bundle(nob)%insector= csi
                          bundle(nob)%fnsector=csf
						  
!............. 'forward'.............
                          if(hchar=='f')then
                              bundle(nob)%ifragment  = ifrag
                              bundle(nob)%ffragment  = ffrag
                              bundle(nob)%isector = is
                              bundle(nob)%fsector = fs
                          endif
!............. 'backward'.................
                         if(hchar=='b')then
                              bundle(nob )%ifragment   = ffrag
                              bundle(nob  )%ffragment  = ifrag
                              bundle(nob  )%isector = fs
                             bundle(nob  )%fsector = is
                         endif
                      end if   ! create
                   end if ! found xsj

               end do  ! isj
!............. ERROR TRAP..................
               if(xsj == -1)then
                  print*,' Problem cannot find sector jumps PP '
                  print*,is,fs
                  print*,x2bjump(1)%nsectjumps
                  print*,xsd(1)%sector(is)%jzX, xsd(1)%sector(fs)%jzX
                  do isj = 1,x2bjump(1)%nsectjumps
                     print*,isj,x2bjump(1)%isector(isj) ,x2bjump(1)%fsector(isj) , & 
      xsd(1)%sector(x2bjump(1)%isector(isj) )%jzX, xsd(1)%sector(x2bjump(1)%fsector(isj))%jzX
                 end do
                 stop
               endif
!............. END ERROR TRAP..................
			   
         end if
!........... ASSIGN NN .......................
         case('NN0')

         if(opstat(is,fs)%nopNN > 0)then
            if(is /= fs)then
               print*,' Huh? NN op bundle ',is,fs
               stop
            end if
!..........FIND START, STOP FOR NN JUMPS..........
!          ASSUME SUBSECTORS SAME AS SECTORS
!  THIS CAN COMBINE TOGETHER SEVERAL NEUTRON 2-body SECTOR JUMPS
!  IN LARGE CASES, HOWEVER, THIS CAN BECOME TOO LARGE, 
!  SO WE NEED TO KEEP TRACK OF THE # OF JUMPS
            insector = -1000002 
            fnsector = -1000001
            xstartx = -99
            xendx   = 0
			localjumps = 0
			foundjumps = .false.
						
            do xsj = 1,x2bjump(2)%nsectjumps   ! loop over neutron sector jumps
               
                ncs =  x2bjump(2)%csjmp(xsj)%ncjmps 
                okay = .false.
                do cs = 1,ncs     ! loop over conjugate (proton) sectors
                  if( is == x2bjump(2)%csjmp(xsj)%cjump(cs)) then   ! if initial proton sector matches
                      okay =.true.
                      ysj = cs     ! this is conjugate proton sector
					  foundjumps = .true.
                      exit
                   end if 
                end do
                if(okay)then
!....................... CHECK WHETHER EXCEEDING JUMPS STORAGE LIMIT...............	
!  NB this was added in 7.9.0  early 2019 (cf 7.8.9) 
		
                    if(localjumps + x2bjump(2)%sjmp(xsj)%njumps > maxjumpsallowed)then  ! FILL CURRENT OPBUNDLE AND CONTINUE
						
						if(localjumps ==0 .and. iproc==0)print*,' Some error in loading jumps '
						
		                xendx = xendx+xstartx-1
!..........FIND START, STOP FOR PROTON SDs.........
		                ystarty = xsd(1)%sector(is)%xsdstart
		                yendy =   xsd(1)%sector(is)%xsdend
		                localops = int(yendy-ystarty+1,8)*int(xendx-xstartx+1,8)
!..........CONVERT TO STATE INDICES
		                cstride = 1
		                nob = nob + 1
		                if(create) then  ! set up for NN operations
		                      bundle(nob)%optype     = 'NN '
		                      bundle(nob)%hchar      = hchar
		                      bundle(nob)%pxstart    = ystarty
		                      bundle(nob)%pxend      = yendy
		                      bundle(nob)%nxstart    = xstartx
		                      bundle(nob)%nxend      = xendx
		                      bundle(nob)%nsortstart = xstartx
		                      bundle(nob)%nsortend   = xendx
		                      bundle(nob)%cstride    = cstride
		                      bundle(nob)%insector   = insector   !check if we actually use this
		                      bundle(nob)%fnsector   = fnsector   ! this is a bit of a kludge, but less wrong than before
!		   				      print*,'ntest a',nob,xstartx,xendx,ystarty,yendy
							  
!............. 'forward'.............
		                      if(hchar=='f')then
		                        bundle(nob)%ifragment  = ifrag
		                        bundle(nob)%ffragment  = ffrag
		                        bundle(nob)%isector = is
		                        bundle(nob)%fsector = fs
		                     endif
                
!............. 'backward'.................
		                     if(hchar=='b')then
		                        bundle(nob )%ifragment   = ffrag
		                        bundle(nob  )%ffragment  = ifrag
		                        bundle(nob  )%isector = fs
		                        bundle(nob  )%fsector = is
		                     end if
		                end if		! create				
!.......... NOW 'RESET'......................................	
                        xstartx = x2bjump(2)%sjmp(xsj)%nstart+1
                        insector =x2bjump(2)%isector(xsj) !xsj						
	                    xendx   = x2bjump(2)%sjmp(xsj)%njumps
	                    fnsector = x2bjump(2)%fsector(xsj) !xsj
					    localjumps = x2bjump(2)%sjmp(xsj)%njumps		
					else     ! localjumps are not too large
		
	                    if(xstartx==-99)then
	                       xstartx = x2bjump(2)%sjmp(xsj)%nstart+1
	                       insector =x2bjump(2)%isector(xsj) !xsj
	                    end if

						
	                    xendx   = xendx+x2bjump(2)%sjmp(xsj)%njumps
!						if(create .and. (is==2 .or. is==3))print*,xsj,x2bjump(2)%sjmp(xsj)%nstart,x2bjump(2)%sjmp(xsj)%njumps, &
!						xstartx,xendx
	                    fnsector = x2bjump(2)%fsector(xsj) !xsj
					    localjumps = localjumps + x2bjump(2)%sjmp(xsj)%njumps	
!						print*,' following ',xstartx,xendx,insector,fnsector									
					end if

			    end if  ! okay
				
				
             end do  ! xsj
             if(.not.foundjumps) then
               print *, "did not find NN jumps in count_create_draft_opbundles_g"
               stop 1
             end if
			 
!................... STILL SHOULD HAVE ONE MORE OPBUNDLE TO FILL	
             xendx = xendx+xstartx-1
!..........FIND START, STOP FOR PROTON SDs.........
!          assume subsectors = sectors!
             ystarty = xsd(1)%sector(is)%xsdstart
             yendy =   xsd(1)%sector(is)%xsdend
             localops = int(yendy-ystarty+1,8)*int(xendx-xstartx+1,8)
!..........CONVERT TO STATE INDICES
             cstride = 1
             nob = nob + 1
             if(create) then  ! set up for NN operations
                   bundle(nob)%optype     = 'NN '
                   bundle(nob)%hchar      = hchar
                   bundle(nob)%pxstart    = ystarty
                   bundle(nob)%pxend      = yendy
                   bundle(nob)%nxstart    = xstartx
                   bundle(nob)%nxend      = xendx
                   bundle(nob)%nsortstart = xstartx
                   bundle(nob)%nsortend   = xendx
                   bundle(nob)%cstride    = cstride
                   bundle(nob)%insector   = insector   !check if we actually use this
                   bundle(nob)%fnsector   = fnsector   ! this is a bit of a kludge, but less wrong than before
				   
!............. 'forward'.............
                   if(hchar=='f')then
                     bundle(nob)%ifragment  = ifrag
                     bundle(nob)%ffragment  = ffrag
                     bundle(nob)%isector = is
                     bundle(nob)%fsector = fs
                  endif
                
!............. 'backward'.................
                  if(hchar=='b')then
                     bundle(nob )%ifragment   = ffrag
                     bundle(nob  )%ffragment  = ifrag
                     bundle(nob  )%isector = fs
                     bundle(nob  )%fsector = is
                  end if
             end if
         end if

!........... ASSIGN PN.......................
         case('PN0')
		 
!		 print*,' measuring ',opstat(is,fs)%nopPN

         if(opstat(is,fs)%nopPN > 0)then
!..........FIND START, STOP FOR P JUMPS..........
!          ASSUME SUBSECTORS SAME AS SECTORS
             xsj = -1
             do isj = 1,x1bjump(1)%nsectjumps
                   if( is == x1bjump(1)%isector(isj) .and.  & 
                        fs == x1bjump(1)%fsector(isj)  ) then
                      xsj = isj
                      exit
                   end if 
             end do
             if(xsj == -1)then
                  print*,' Problem cannot find sector jumps PN '
                  print*,is,fs,hchar
                  stop
             endif

             xstartx = x1bjump(1)%sjmp(xsj)%nstart+1
             xendx = xstartx+ x1bjump(1)%sjmp(xsj)%njumps-1

             if(xendx < xstartx)cycle
!........LOOP OVER CONJUGATE N JUMPS

             ncs =  x1bjump(1)%csjmp(xsj)%ncjmps 
!
!.... NOTE THAT EACH "PN" OPBUNDLE IS DEFINED 
! FOR A SINGLE PROTON SECTOR JUMP and A SINGLE NEUTRON SECTOR JUMP
! WE DO THIS SO WE CAN SORT THE NEUTRON (AND PROTON) JUMPS WITHIN
! A SECTOR JUMP. 
!
!if(create)print*,xsj,' proton sector jump'
             do cs = 1,ncs

                 ysj = x1bjump(1)%csjmp(xsj)%cjump(cs)
				 
!				 if(create)print*,xsj,ysj, ' jump creation '
                 ystarty = x1bjump(2)%sjmp(ysj)%nstart+1
                 yendy   = ystarty+ x1bjump(2)%sjmp(ysj)%njumps-1
                 if(yendy < ystarty)cycle
                 localops = int(yendy-ystarty+1,8)*int(xendx-xstartx+1,8)
				 
!				 if(create)print*,ysj,' neutron sector jump '				 
!				 if(create)print*,' check ',hchar,is,fs,x1bjump(1)%sjmp(xsj)%diag,x1bjump(2)%sjmp(ysj)%diag

                 if( (hchar == 'h' .and.( is == fs .and. x1bjump(1)%sjmp(xsj)%diag & 
                                                  .and. x1bjump(2)%sjmp(ysj)%diag) )  & 
                    .or. ( (hchar == 'f' .or. hchar=='b') .and. (is /= fs  & 
           .or..not.x1bjump(1)%sjmp(xsj)%diag .or. .not. x1bjump(2)%sjmp(ysj)%diag)))then
                 nob = nob + 1
                 if(create)then
!					 print*,' going forward'

                   bundle(nob  )%optype     = 'PN '
                   bundle(nob  )%hchar      = hchar

                   bundle(nob  )%pxstart    = xstartx
                   bundle(nob  )%pxend      = xendx
                   bundle(nob  )%nxstart    = ystarty
                   bundle(nob  )%nxend      = yendy
                   bundle(nob  )%nsortstart = ystarty
                   bundle(nob  )%nsortend   = yendy
                   if(hchar == 'f')then
                         bundle(nob )%ifragment   = ifrag
                         bundle(nob  )%ffragment  = ffrag
                         bundle(nob  )%isector = is
                         bundle(nob  )%fsector = fs
                         bundle(nob  )%insector   = x1bjump(2)%isector(ysj) 
                         bundle(nob  )%fnsector   = x1bjump(2)%fsector(ysj) 
                   endif
                   if(hchar == 'b' .or. hchar == 'h')then
                         bundle(nob )%ifragment   = ffrag
                         bundle(nob  )%ffragment  = ifrag
                         bundle(nob  )%isector = fs
                         bundle(nob  )%fsector = is
                         bundle(nob  )%insector   = x1bjump(2)%fsector(ysj) 
                         bundle(nob  )%fnsector   = x1bjump(2)%isector(ysj) 
                   endif

                 endif
                 endif
             end do  ! ics

         end if
!........... ASSIGN SPE .......................
!     NOTE the SPE opbundle has been de-emphasized, with single particle potentials
!     absorbed into two-body operators
!     It is still needed, however, when we have a single particle or hole of either species
!
         case('SPE')

         if( is == fs)then
            nob = nob + 1

!..........FIND START, STOP FOR PSDs..........

            xstartx = xsd(1)%sector(is)%xsdstart !subsectorlist(is)%pSDstart
            xendx   = xsd(1)%sector(is)%xsdend  ! subsectorlist(is)%pSDend

!..........FIND START, STOP FOR NEUTRON SDs.........

            xsj = is ! subsectorlist(is)%parentsector
            cs  = xsd(1)%sector(xsj)%csector(1)
            ystarty = xsd(2)%sector(cs)%xsdstart
            ncs = xsd(1)%sector(xsj)%ncsectors
            cs  = xsd(1)%sector(xsj)%csector(ncs) 
            yendy = xsd(2)%sector(cs)%xsdend
            ! simple check that bundle will stay in bounds
            ibasis = pstart(xstartx) + nstart(ystarty)
            if(ibasis > fragmentlist(ifrag)%basisend) then
               write(6,"(A,I0)") "start ibasis=", ibasis
               write(6,"(A,I0,A,I0)") "pxstart=", xstartx, ", pxend=", xendx
               write(6,"(A,I0,A,I0)") "nxstart=", ystarty, ", nxend=", yendy
               call errstop("SPE out of bounds in count_create_draft_opbundles_g")
            end if
            ibasis = pstart(xendx) + nstart(yendy)
            if(ibasis > fragmentlist(ifrag)%basisend) then
               write(6,"(A,I0)") "end ibasis=", ibasis
               write(6,"(A,I0,A,I0)") "pxstart=", xstartx, ", pxend=", xendx
               write(6,"(A,I0,A,I0)") "nxstart=", ystarty, ", nxend=", yendy
               call errstop("SPE out of bounds in count_create_draft_opbundles_g")
            end if
            localops = int(yendy-ystarty+1,8)*int(xendx-xstartx+1,8)
!............. DIAGONAL so no directions .............
            if(create) then  ! set up for SPE operations

               bundle(nob  )%optype     = 'SPE'
               bundle(nob  )%hchar      = 'h'   ! not needed
               bundle(nob )%ifragment   = ffrag
               bundle(nob  )%ffragment  = ifrag
               bundle(nob  )%isector = fs
               bundle(nob  )%fsector = is
               bundle(nob)%insector= cs
               bundle(nob)%fnsector=cs
               bundle(nob  )%pxstart    = xstartx
               bundle(nob  )%pxend      = xendx
               bundle(nob  )%nxstart    = ystarty
               bundle(nob  )%nxend      = yendy
             end if
          end if

!........... ASSIGN PPP .......................
         case('PPP')

         if(opstat(is,fs)%nopPPP > 0)then
!..........FIND START, STOP FOR PPP JUMPS..........
!          ASSUME SUBSECTORS SAME AS SECTORS
            xsj = -1
			if(ifrag==0 .or. ffrag==0)then
				print*,' error fragment 0',ifrag,ffrag
			end if
            do isj = 1,x3bjump(1)%nsectjumps
                   if( is == x3bjump(1)%isector(isj) .and.  & 
                        fs == x3bjump(1)%fsector(isj)  ) then
                      xsj = isj
                      exit
                   end if 

            end do

            if(xsj == -1)then
                  print*,' Problem cannot find sector jumps  PPP'
                  print*,is,fs
                  stop
            endif

            xstartx = x3bjump(1)%sjmp(xsj)%nstart+1
            xendx   = xstartx + x3bjump(1)%sjmp(xsj)%njumps-1

!..........FIND START, STOP FOR NEUTRON SDs.........
            csi = x3bjump(1)%csjmp(xsj)%cjump(1)  ! this is the first conjugate
                                                     ! neutron sector
            ystarty = xsd(2)%sector(csi)%xsdstart
            ncs = x3bjump(1)%csjmp(xsj)%ncjmps
            csf  = x3bjump(1)%csjmp(xsj)%cjump(ncs)  
            yendy = xsd(2)%sector(csf)%xsdend
            cstride = 1

!..........CONVERT TO STATE INDICES

            localops = int(yendy-ystarty+1,8)*int(xendx-xstartx+1,8)
            nob = nob + 1
            if(create) then  ! set up for PPP operations
               bundle(nob)%optype     = 'PPP'
               bundle(nob)%hchar      = hchar
               bundle(nob)%pxstart    = xstartx
               bundle(nob)%pxend      = xendx
               bundle(nob)%nxstart    = ystarty
               bundle(nob)%nxend      = yendy
               bundle(nob)%cstride    = cstride
               bundle(nob)%insector= csi
               bundle(nob)%fnsector=csf
!............. 'forward'.............
               if(hchar == 'f')then
                  bundle(nob)%ifragment  = ifrag
                  bundle(nob)%ffragment  = ffrag
                  bundle(nob)%isector = is
                  bundle(nob)%fsector = fs
               endif
               if(hchar == 'b')then
!............. 'backward'.................
                  bundle(nob )%ifragment   = ffrag
                  bundle(nob  )%ffragment  = ifrag
                  bundle(nob  )%isector = fs
                  bundle(nob  )%fsector = is
               endif
            end if
         end if

!........... ASSIGN NNN .......................
         case('NNN')

         if(opstat(is,fs)%nopNNN > 0)then
            if(is /= fs)then
               print*,' Huh? NNN op bundle ',is,fs
               stop
            end if
!..........FIND START, STOP FOR NN JUMPS..........
            xstartx = -99
            xendx = 0
            do xsj = 1,x3bjump(2)%nsectjumps
               
                ncs =  x3bjump(2)%csjmp(xsj)%ncjmps 
                okay = .false.
                do cs = 1,ncs
                  if( is == x3bjump(2)%csjmp(xsj)%cjump(cs)) then
                      okay =.true.
                      ysj = cs
                      exit
                   end if 
                end do
                if(.not.okay)cycle
                if(xstartx==-99)xstartx = x3bjump(2)%sjmp(xsj)%nstart+1
                xendx   = xendx+x3bjump(2)%sjmp(xsj)%njumps 
             end do
             xendx = xendx +xstartx-1

!..........FIND START, STOP FOR PROTON SDs.........
                ystarty = xsd(1)%sector(is)%xsdstart
                yendy =   xsd(1)%sector(is)%xsdend
                localops = int(yendy-ystarty+1,8)*int(xendx-xstartx+1,8)
                if(localops == 0)then
                    print*,' NNN should have no operations '
                    stop
                end if
!..........CONVERT TO STATE INDICES

                   cstride = 1
                nob = nob + 1
                if(create) then  ! set up for NN operations

                  bundle(nob)%optype     = 'NNN'
                  bundle(nob)%hchar      = hchar
                  bundle(nob)%pxstart    = ystarty
                  bundle(nob)%pxend      = yendy
                  bundle(nob)%nxstart    = xstartx
                  bundle(nob)%nxend      = xendx
                  bundle(nob)%cstride    = cstride
                  bundle(nob)%insector= cs
                  bundle(nob)%fnsector=cs
!............. 'forward'.............
                  if(hchar=='f')then
                     bundle(nob)%ifragment  = ifrag
                     bundle(nob)%ffragment  = ffrag
                     bundle(nob)%isector = is
                     bundle(nob)%fsector = fs
                  endif
!............. 'backward'.................
                  if(hchar=='b')then
                     bundle(nob )%ifragment   = ffrag
                     bundle(nob  )%ffragment  = ifrag
                     bundle(nob  )%isector = fs
                     bundle(nob  )%fsector = is
                  endif
                end if
         end if

!........... ASSIGN PPN.......................
         case('PPN')

         if(opstat(is,fs)%nopPPN > 0)then

!..........FIND START, STOP FOR P JUMPS..........
!          ASSUME SUBSECTORS SAME AS SECTORS
             do ysj = 1,x1bjump(2)%nsectjumps
                ncs =  x1bjump(2)%csjmp(ysj)%ncjmps 

                do cs = 1,ncs
                   xsj = x1bjump(2)%csjmp(ysj)%cjump(cs)
                   if( (is == x2bjump(1)%isector(xsj) .and.  & 
                        fs == x2bjump(1)%fsector(xsj)) .or.   & 
                       (fs == x2bjump(1)%isector(xsj) .and.  & 
                        is == x2bjump(1)%fsector(xsj) )) then
     

                 xstartx = x2bjump(1)%sjmp(xsj)%nstart+1
                 xendx = xstartx+ x2bjump(1)%sjmp(xsj)%njumps-1
                 if(xendx < xstartx)cycle

!........LOOP OVER CONJUGATE N JUMPS

                 ystarty = x1bjump(2)%sjmp(ysj)%nstart+1
                 yendy   = ystarty+ x1bjump(2)%sjmp(ysj)%njumps-1
                 if(yendy < ystarty)cycle

                 nob = nob + 1
                 if(create)then

                   bundle(nob  )%optype     = 'PPN'
                   bundle(nob  )%hchar      = hchar 
                   bundle(nob )%ifragment   = ifrag
                   bundle(nob  )%ffragment  = ffrag
                   bundle(nob  )%isector = is
                   bundle(nob  )%fsector = fs
                   bundle(nob  )%pxstart    = xstartx
                   bundle(nob  )%pxend      = xendx
                   bundle(nob  )%nxstart    = ystarty
                   bundle(nob  )%nxend      = yendy
                   bundle(nob  )%nsortstart = ystarty
                   bundle(nob  )%nsortend   = yendy
                   if(hchar == 'f')then
                      bundle(nob )%ifragment   = ifrag
                      bundle(nob  )%ffragment  = ffrag
                      bundle(nob  )%isector = is
                      bundle(nob  )%fsector = fs
                      bundle(nob  )%insector   = x1bjump(2)%isector(ysj) 
                      bundle(nob  )%fnsector   = x1bjump(2)%fsector(ysj) 
                   end if
                   if(hchar == 'b')then
                      bundle(nob )%ifragment   = ffrag
                      bundle(nob  )%ffragment  = ifrag
                      bundle(nob  )%isector = fs
                      bundle(nob  )%fsector = is
                      bundle(nob  )%insector   = x1bjump(2)%fsector(ysj) 
                      bundle(nob  )%fnsector   = x1bjump(2)%isector(ysj) 
                   end if
                 end if
                 end if

               end do  ! cs
             end do ! ysj

         end if
!........... ASSIGN PNN.......................
         case('PNN')

         if(opstat(is,fs)%nopPNN > 0)then
                   
!..........FIND START, STOP FOR P JUMPS..........
!          ASSUME SUBSECTORS SAME AS SECTORS
             do xsj = 1,x1bjump(1)%nsectjumps

                   if( (is == x1bjump(1)%isector(xsj) .and.  & 
                        fs == x1bjump(1)%fsector(xsj))  .or.  & 
                          (fs == x1bjump(1)%isector(xsj) .and.  & 
                        is == x1bjump(1)%fsector(xsj))) then

             xstartx = x1bjump(1)%sjmp(xsj)%nstart+1
             xendx = xstartx+ x1bjump(1)%sjmp(xsj)%njumps-1
             if(xendx < xstartx)cycle

!........LOOP OVER CONJUGATE N JUMPS

             ncs =  x1bjump(1)%csjmp(xsj)%ncjmps 

             do cs = 1,ncs
                 ysj = x1bjump(1)%csjmp(xsj)%cjump(cs)
                 ystarty = x2bjump(2)%sjmp(ysj)%nstart+1
                 yendy   = ystarty+ x2bjump(2)%sjmp(ysj)%njumps-1
                 if(yendy < ystarty)cycle

                 nob = nob + 1

                 if(create)then

                   bundle(nob  )%optype     = 'PNN'
                   bundle(nob  )%hchar      = hchar
                   bundle(nob  )%pxstart    = xstartx
                   bundle(nob  )%pxend      = xendx
                   bundle(nob  )%nxstart    = ystarty
                   bundle(nob  )%nxend      = yendy
                   bundle(nob  )%nsortstart = ystarty
                   bundle(nob  )%nsortend   = yendy
                   if(hchar=='f')then
                      bundle(nob )%ifragment   = ifrag
                      bundle(nob  )%ffragment  = ffrag
                      bundle(nob  )%isector = is
                      bundle(nob  )%fsector = fs
                      bundle(nob  )%insector   = x2bjump(2)%isector(ysj) 
                      bundle(nob  )%fnsector   = x2bjump(2)%fsector(ysj) 
                   endif
                   if(hchar=='b')then
                      bundle(nob )%ifragment   = ffrag
                      bundle(nob  )%ffragment  = ifrag
                      bundle(nob  )%isector = fs
                      bundle(nob  )%fsector = is
                      bundle(nob  )%insector   = x2bjump(2)%fsector(ysj) 
                      bundle(nob  )%fnsector   = x2bjump(2)%isector(ysj) 
                   endif
                 end if

               end do  ! cs
             end if
           end do   ! xsj
         end if
!............... ADDED 7.9.7..................		 
		 case ('CEN')
		 if(ifrag/=ffrag)cycle
		 if(is /= fs)cycle
		 
		 if(hchar/='h')then
			 print*,' mis-setup of centroid opbundle '
			 stop
		 end if
		 
		 xstartx=xsd(1)%sector(is)%xsdstart
		 xendx  =xsd(1)%sector(is)%xsdend
		 
		 csi    =xsd(1)%sector(is)%csector(1)
		 ncs    =xsd(1)%sector(is)%ncsectors
		 csf    =xsd(1)%sector(is)%csector(ncs)
         ystarty = xsd(2)%sector(csi)%xsdstart
         yendy = xsd(2)%sector(csf)%xsdend
         localops = int(yendy-ystarty+1,8)*int(xendx-xstartx+1,8)		 

         cstride = 1
         nob = nob+1
         if(create) then  ! set up for PP operations
                          bundle(nob)%optype     = 'CEN'
                          bundle(nob)%hchar      =  hchar
                          bundle(nob)%pxstart    = xstartx
                          bundle(nob)%pxend      = xendx
                          bundle(nob)%nxstart    = ystarty
                          bundle(nob)%nxend      = yendy
                          bundle(nob)%nsortstart    = ystarty
                          bundle(nob)%nsortend      = yendy
                          bundle(nob)%cstride    = cstride
                          bundle(nob)%insector= csi
                          bundle(nob)%fnsector=csf
                              bundle(nob)%ifragment  = ifrag
                              bundle(nob)%ffragment  = ffrag
                              bundle(nob)%isector = is
                              bundle(nob)%fsector = fs
							  
	!						  print*,' CEN opbundle ',nob, xstartx,xendx,' sector ',is,ystarty,yendy

         end if   ! create
		 
		 
		 case default
		 
		 print*,' bad optype ', optype
		 stop
         end select
      end do ! fs
  end do  ! is
  return
end subroutine count_create_draft_opbundles_g

!===============================================================
!
! key routine to distribute operations across multiple compute nodes
! 
! Target # of operations is NOPSNODE.  Loop over "draft" opbundles
! starting where we left off 
!    -- NODESTARTOB is starting draft opbundle
!    -- and NODESTARTZ is where in sets of operations started
! then looping over remaining opbundles, count up operations until
! close to NOPSNODE.  This leaves NODESTOPOB (end of op bundles)
! and NODESTOPSZ (end of sets of operations)
! 
! INPUT:
!   nob_draft: # of "draft" opbundles
!   nodestartob:  which opbundle to start on
!   nodestartz:   which set of operations to start on
!
! OUTPUT:
!   nodestopob: which obbundle to stop on
!   nodestopz : which set of operations to stop on
!
! CALLED BY:
!   split_draft_opbundles_frag in bparallel_main.f90
!
subroutine search_draft_opbundles4split_g ( nob_draft, nodestartob, nodestopob, nodestartz, nodestopz)

   use nodeinfo
   use operation_stats
   implicit none
   integer :: nob_draft
   integer :: nodestartob, nodestopob
   integer(8) :: nodestartz,nodestopz
   integer(8) :: opsofar,diff_nop,opsleft

   integer(8) :: nset
   real (8) :: localwt

   integer iob

   opsofar = 0
   if(nodestartob > nob_draft)then
       print*,' error in starting search_draft_opbundles4split_g ',nodestartob,nob_draft
       stop
   end if
   do iob = nodestartob,nob_draft    ! LOOP OVER DRAFT OPBUNDLES
      if(draft_opbundle(iob)%pxstart < 0 .and. iproc ==0)then
           print*,iob,draft_opbundle(iob)%optype,' pxstart ',draft_opbundle(iob)%pxstart
      end if
      if(draft_opbundle(iob)%pxend < 0 .and. iproc ==0)then
           print*,iob,draft_opbundle(iob)%optype,' pxend ',draft_opbundle(iob)%pxend
      end if
      if(draft_opbundle(iob)%nxstart < 0.and. iproc ==0)then
           print*,iob,draft_opbundle(iob)%optype,' nxstart ',draft_opbundle(iob)%nxstart
      end if
      if(draft_opbundle(iob)%nxend < 0.and. iproc ==0)then
           print*,iob,draft_opbundle(iob)%optype,' nxend ',draft_opbundle(iob)%nxend
      end if

      if(draft_opbundle(iob)%nops == 0)cycle   ! HAS NO OPERATIONS, SKIP
      if(usewts)then
		  localwt = draft_opbundle(iob)%opwt
!      select case (draft_opbundle(iob)%optype)
!         case('PP ')
!             localwt = opwtPP
!             if (draft_opbundle(iob)%hchar == 'b') localwt = opwtPPb
!         case('NN ')
!             localwt = opwtNN
!         case('PN ')
!             localwt = opwtPN
!             if (draft_opbundle(iob)%hchar == 'b') localwt = opwtPNb
!         case('SPE')
!             localwt = opwtSPE
!         case('PPP')
!             localwt = opwtPPP
!         case('NNN')
!             localwt = opwtNNN
!         case('PPN')
!             localwt = opwtPPN
!             if (draft_opbundle(iob)%hchar == 'b') localwt = opwtPPNb
!         case('PNN')
!             localwt = opwtPNN
!             if (draft_opbundle(iob)%hchar == 'b') localwt = opwtPNNb
!         case default
!          if(iproc == 0)print*,' OOPS bad optype ',draft_opbundle(iob)%optype
!          stop
!      end select

      else
        localwt = 1.d0
      endif
!	  if(iproc==0)print*,' local weight ',localwt,nopsnode,nopsnode*32
!.......... COMPUTE # OF OPS AVAILABLE FROM THIS OPBUNDLE.......
      opsleft = int(draft_opbundle(iob)%nops*localwt,8)
      if(iob == nodestartob)then      ! CORRECT FOR OPS ALREADY USED 
         opsleft = opsleft-int(int(nodestartz-1,8)*int(draft_opbundle(iob)%min_nop,8)*localwt,8)
      end if
!.............. HOW MANY OPS NEEDED .......
      diff_nop = nopsnode-opsofar
      if(opsleft <= diff_nop)then  !  DON'T SPLIT;  3 cases


!............ LOOK AHEAD...........

          if(iob == nob_draft)then   ! CASE 1: must be at the end anyway
             nodestopob = iob
             nodestopz  = draft_opbundle(iob)%nsetops    
             opsofar = opsofar + opsleft   ! not needed but for debugging
             return
          else  ! SEE WHICH IS CLOSER TO IDEAL; CLOSING OUT HERE OR GOING TO NEXT OPBUNDLE

             if(abs( opsleft - diff_nop) < abs(opsleft-diff_nop+ & 
                   int( draft_opbundle(iob+1)%min_nop*localwt,8) ) )then 
                 nodestopob = iob
                 nodestopz  = draft_opbundle(iob)%nsetops    
                 opsofar = opsofar + opsleft   ! not needed but for debugging
                 return
             else
                 opsofar = opsofar + opsleft
                 cycle    ! keep looping
             end if
          end if

      else    ! must split
          nset = nint( real(diff_nop,8) / (real( draft_opbundle(iob)%min_nop, 8)*localwt), 8)
          if(iob == nodestartob)then
              nodestopz = nodestartz +nset-1
          else
              nodestopz = nset
          endif
          if(nodestopz > draft_opbundle(iob)%nsetops)then   ! ERROR TRAP
                 print *,' WENT TOO FAR! YIKES! '
                 print *, 'iob=', iob
                 print *, "nodestartob=", nodestartob
                 print *, "nodestopz=", nodestopz
                 print *, "draft_opbundle(iob)%nsetops=", draft_opbundle(iob)%nsetops
                 print *, "diff_nop=", diff_nop
                 print *, "draft_opbundle(iob)%nops=", draft_opbundle(iob)%nops
                 print *, "draft_opbundle(iob)%min_nop=", draft_opbundle(iob)%min_nop
                 print *, "draft_opbundle(iob)%optype=", draft_opbundle(iob)%optype
                 print *, "opsleft=", opsleft
                 print *,  "draft_opbundle(iob)%nops-int(nodestartz-1,8)*int(draft_opbundle(iob)%min_nop,8)", &
                     draft_opbundle(iob)%nops-int(nodestartz-1,8)*int(draft_opbundle(iob)%min_nop,8)
                 print *, "nset=", nset
                 print *, "localwt=", localwt
                 stop
          end if

          nodestopob = iob

 !............ THE FOLLOWING FOR DEBUGGING ONLY...........
      opsleft = int(nodestopz* draft_opbundle(iob)%min_nop*localwt,8)
      if(iob == nodestartob)then      ! CORRECT FOR OPS ALREADY USED 
         opsleft = opsleft-int(int(nodestartz-1,8)*int(draft_opbundle(iob)%min_nop,8)*localwt,8)
      end if         
          opsofar = opsofar + opsleft
          return

      end if

   end do ! iob
   return
end subroutine search_draft_opbundles4split_g

!==============================================================
! Debugging routine to dump the specification of a bundle
!
!  CALLED BY: 
!     check_opbundle_bounds
!
subroutine print_opbundle(fn, iob)
   use precisions
!   use opbundles
   use butil_mod
   implicit none
   integer, intent(in) :: fn ! file number
   integer, intent(in) :: iob ! index of bundle
   type(bund) :: b
   b = opbundle(iob)
   write(fn, "(A,I0,A)") "opbundle(", iob, ")"
   write(fn, "(A,I0)")   "   node=", b%node
   write(fn, "(A,A)")    "   optype=", b%optype
   write(fn, "(A,A)")    "   hchar=", b%hchar
   write(fn, "(A,I0,A,I0)") "   isector=", b%isector, ", fsector=", b%fsector
   write(fn, "(A,I0,A,I0)") "   insector=", b%insector, ", fnsector=", b%fnsector
   write(fn, "(A,I0,A,I0)") "   ifrag=", b%ifragment, "ofrag=", b%ffragment
   write(fn, "(A,I0,A,I0)") "   pxstart=", b%pxstart, ", pxend=", b%pxend
   write(fn, "(A,I0,A,I0)") "   nxstart=", b%nxstart, ", nxend=", b%nxend
   write(fn, "(A,I0)") "   cstride=",b%cstride
   write(fn, "(A,I0)") "   njumps=",b%njumps
   flush(fn)
end subroutine print_opbundle

!===============================================================
!
! Verifies that when the opbundles are evaluated, that the vec1/vec2 indices
! will stay in bounds for the fragments on the MPI processes they are
! assigned to.
!
! ADDED in 7.9.0: checks for "monster" opbundles which have large number of jumps
! and storage needs to be pre-reserved 
!
!  CALLED BY: 
!     MAIN
!
subroutine check_opbundle_bounds
   use precisions
   use nodeinfo
   use basis
   use fragments
!   use opbundles
   use butil_mod
   implicit none
   integer :: jnode  ! current simulated mpi rank
   integer :: iob    ! opbundle index on jnode
   integer :: frag1, frag2 ! fragments for vec1 and vec2
   integer(kind=basis_prec) :: f1s, f1e, f2s, f2e  ! fragment bounds
   integer(kind=8) :: psdstart, psdend   ! proton slater determinant bounds
   integer(kind=8) :: nsdstart, nsdend   ! neutron slater determinant bounds
   integer(kind=8) :: cstride
   integer(kind=8) :: imax, imin
   type(bund),save :: b

   if(iproc /= 0) return; ! only test in root
   do jnode = 0,nprocs-1
      ! get fragment bounds for rank jnode
      frag1 = nodal(jnode)%ifragment
      frag2 = nodal(jnode)%ffragment
      ! get bounds for vec1 and vec2
      f1s = basestart(frag1)
      f1e = basestop(frag1)
      f2s = basestart(frag2)
      f2e = basestop(frag2)
      do iob = opbundlestart(jnode),opbundleend(jnode)
         psdstart = opbundle(iob)%pxstart
         psdend   = opbundle(iob)%pxend
         nsdstart = opbundle(iob)%nxstart
         nsdend   = opbundle(iob)%nxend
         select case (opbundle(iob)%optype)
         case('SPE')
            imin = pstart(psdstart) + nstart(nsdstart)
            imax = pstart(psdend) + nstart(nsdend)
            if(imin < f1s .or. imin < f2s .or. imax > f1e .or. imax > f2e) then
               b = opbundle(iob)
               write(6,"(A,I0,A,I0)")  "opbundle generates out of bounds references ", &
                  imin, " or ", imax
               call print_node_bounds(6, jnode)
               call print_opbundle(6, iob)
               call errstop("check_opbundles_bounds SPE")
            end if
         case('PP ')
         case('NN ')
         case('PN ')
         case default
            ! also want to check 3 body some day
         end select
      end do ! iob
   end do ! jnode
   flush(6)
end subroutine check_opbundle_bounds

!===============================================================
!
! for a given opbundle, computes the total # operations in it, 
! and finds the minimal split; 
! this will be important for refining the distribution
!
! IN ADDITION (added 7.7.5) LOOK FOR OPBUNDLES WITH LARGE STORAGE REQUIREMENTS
!!
! INPUT:
!   iob : which opbundle
!   draft: logical flag; if T, then draft_opbundle, else opbundle
!
! OUTPUT: localjumpstorage,  how much memory is used
!
! CALLED BY:
!  draft_opbundles
!
subroutine analyze_opbundle(iob,draft,localjumpstorage)

    use jumpNbody
	use jump3body
	use operation_stats
	use fragments
	use nodeinfo
	use flagger
	use io
   implicit none
   integer,intent(in) :: iob   ! which opbundle
   logical,intent(in) :: draft  ! if draft or final opbundle
   real(8)    :: localjumpstorage        ! how much storage is required
   
   type (bund),pointer :: mybundle
   integer(8) :: xstartx,xendx,ystarty,yendy,nlocaljumps

!.................... TEMPORARY: FOR CONFIRMING COUNTS OF OPERATIONS............   
   integer :: is,fs   ! initial, final sectors; for testing
   integer(8) :: noptest,mynops,secondtest
   
   if(draft)then
	   mybundle=>draft_opbundle(iob)
   else
	   mybundle=>opbundle(iob)
   end if

  
!...... COMPUTE TOTAL # OF OPERATIONS IN THIS BUNDLE .....
   xstartx = mybundle%pxstart
   xendx   = mybundle%pxend
   ystarty = mybundle%nxstart
   yendy =   mybundle%nxend

   mybundle%nops = int(xendx-xstartx+1,8)*int(yendy-ystarty+1,8)
!....COMPUTE MIN SET OF OPS WHEN I SPLIT UP AN OPBUNDLE....
!
! ... represent # of operations as min_nop x nsetops
!     where min_nop is a block or set of operations
!     and nsetops is the # of such sets
!     E.G. # of PP operations is # neutrons SD (here = min_nop)
!            x # of PP jumps
!     IMPORTANT: In downstream routines, will split up
!     nsetops 
!
! *  note: for PP, NN, PPP, NNN this is easy: it's the # of conjugate SDs
! *  trickier for PN, PNN, PPN, because no obvious split
! *  choose one-body jumps because fewer than two-body
! *  for PN, chose proton jumps because (?) fewer for most cases
!
!  NOTE THIS IS ONLY NEEDED FOR DRAFT OPBUNDLES
!
   select case ( mybundle%optype)

       case('PP ','PPP','SPE')
       mybundle%min_nop = yendy -ystarty+1  ! governed by dimension of conjugate neutron SDs
        mybundle%nsetops = xendx -xstartx + 1

       case('NN ','NNN')
        mybundle%min_nop = xendx -xstartx+1  ! dimension of conjugate proton SDs
        mybundle%nsetops = yendy -ystarty + 1

       case('PN ', 'PNN')
       mybundle%min_nop = xendx-xstartx+1  ! dimension of proton 1-body jumps
        mybundle%nsetops = yendy -ystarty + 1

       case('PPN')
        mybundle%min_nop = yendy -ystarty+1  ! dimension of neutron 1-body jumps
        mybundle%nsetops = xendx -xstartx + 1
	   
		case('cen','CEN')
        mybundle%min_nop = xendx-xstartx+1  ! dimension of proton 1-body jumps
         mybundle%nsetops = yendy -ystarty + 1
		 
	   case default
	   print*,' OH NOES I could not find the optype (analyze_opbundles )', mybundle%optype
	   stop

   end select
!............ COUNT # OF JUMPS STORED......
   nlocaljumps = 0   
   select case(mybundle%optype)
	   
      case('PP0','PPP','PP')
         nlocaljumps = nlocaljumps + int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1

      case('NN0','NNN','NN')
	  !  NOTE nsortstart,nsortend not defined for draft opbundle
	  if(draft)then
         nlocaljumps = nlocaljumps + int(mybundle%nxend+1,8)-int(mybundle%nxstart,8)+1
	  else
         nlocaljumps = nlocaljumps + int(mybundle%nsortend+1,8)-int(mybundle%nsortstart,8)+1	 
	  end if

      case('PN0','PPN','PNN','PN','SPE','CEN')
         nlocaljumps = nlocaljumps + int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1
		 if(draft)then
             nlocaljumps = nlocaljumps + int(mybundle%nxend+1,8)-int(mybundle%nxstart,8)+1
	     else
             nlocaljumps = nlocaljumps + int(mybundle%nsortend+1,8)-int(mybundle%nsortstart,8)+1
	     end if
   end select
   mybundle%njumps = nlocaljumps   
   
   if(nlocaljumps < 1)then
	   print*,' hey some issue with # of jumps ',iob,draft, mybundle%optype
   end if
   
!........... DETERMINE STORAGE REQUIREMENTS AND IF IT GOES OVER....................   
    localjumpstorage = 0.d0
    if ( .not.restrictjumps .or. nprocs ==1)return
    select case(mybundle%optype)
	   case('PP')
  	      localjumpstorage = real(int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1,8)*bytesper2Bjump
   	   case('NN')
	      if(draft)then
     	    localjumpstorage = real(int(mybundle%nxend+1,8)-int(mybundle%nxstart,8)+1,8)*bytesper2Bjump 
		  else
       	    localjumpstorage = real(int(mybundle%nsortend+1,8)-int(mybundle%nsortstart,8)+1,8)*bytesper2Bjump 
		  endif 
   	   case('PN')
	     if(draft)then
     	   localjumpstorage = real(int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1,8)*bytesper1Bjump + & 
		    real(int(mybundle%nxend+1,8)-int(mybundle%nxstart,8)+1,8)*bytesper1Bjump 
		else
      	   localjumpstorage = real(int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1,8)*bytesper1Bjump + & 
 		    real(int(mybundle%nsortend+1,8)-int(mybundle%nsortstart,8)+1,8)*bytesper1Bjump 			
		end if
	   case('PPP')
	   	   localjumpstorage = real(int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1,8)*bytesper3Bjump
	   case('NNN')
	   if(draft)then
	      localjumpstorage = real(int(mybundle%nxend+1,8)-int(mybundle%nxstart,8)+1,8)*bytesper3Bjump 
	  else
	      localjumpstorage = real(int(mybundle%nsortend+1,8)-int(mybundle%nsortstart,8)+1,8)*bytesper3Bjump 
		  
	  endif
	   case('PPN')
	   if(draft)then
 	       localjumpstorage = real(int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1,8)*bytesper2Bjump + & 
	    real(int(mybundle%nxend+1,8)-int(mybundle%nxstart,8)+1,8)*bytesper1Bjump
	  else
          localjumpstorage = real(int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1,8)*bytesper2Bjump + & 
        real(int(mybundle%nsortend+1,8)-int(mybundle%nsortstart,8)+1,8)*bytesper1Bjump
	   endif 	   
 	   case('PNN')
	   if(draft)then
 	       localjumpstorage = real(int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1,8)*bytesper1Bjump + & 
	    real(int(mybundle%nxend+1,8)-int(mybundle%nxstart,8)+1,8)*bytesper2Bjump
	  else
          localjumpstorage = real(int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1,8)*bytesper1Bjump + & 
        real(int(mybundle%nsortend+1,8)-int(mybundle%nsortstart,8)+1,8)*bytesper2Bjump
	   endif 	  	
    end select
!    print*,iob, mybundle%optype, localjumpstorage , maxjumpmemory_default*1.0e9 , 'jump storage'
	if(localjumpstorage > maxjumpmemory_default*1.0e9 )then  ! flag for possible splitting
		if(iproc == 0)then
			excess_memory_needed = excess_memory_needed + localjumpstorage*10e-9 - maxjumpmemory_default
			write(logfile,*)' draft opbundle ',iob,' of type ',mybundle%optype,' requires memory ',localjumpstorage*1.0e-9, & 
			', greater than limit of ',maxjumpmemory_default
			if(nfragments > 1)write(logfile,*)' Located on fragments ',mybundle%ifragment,' -> ',mybundle%ffragment
		end if
	end if	
	mybundle%jumpstore = localjumpstorage  !*1.0e-9
	
!   print*,' storage per op  ',localjumpstorage/real(mybundle%nops,kind=8 )
   
   return

end subroutine analyze_opbundle

!===============================================================

subroutine count_all_operations(draft,nbundles,sumofallops)

    use nodeinfo
	implicit none
	logical, intent(IN) :: draft  
	integer,intent(IN) :: nbundles
	integer(8),intent(OUT) :: sumofallops
	
	integer :: ibundle
	
	type (bund),pointer :: mybundle(:)


    if(nbundles < 1)then
		print*,' oops have not set # of opbundles '
	end if		
	if(draft)then
		mybundle=>draft_opbundle
	else
		mybundle=>opbundle
	end if
	
	sumofallops = 0
	do ibundle = 1,nbundles		
		sumofallops = sumofallops + mybundle(ibundle)%nops
		
!		print*,ibundle,mybundle(ibundle)%pxend,mybundle(ibundle)%pxstart,mybundle(ibundle)%nxend,mybundle(ibundle)%nxstart
		
	end do
	
	return	

end subroutine count_all_operations

!===============================================================

subroutine count_all_operations_new(draft,nbundles,sumofallops)

    use nodeinfo
	use fragments
	use flagger
	use jumpstart
	use io
	implicit none
	logical, intent(IN) :: draft  
	integer,intent(IN) :: nbundles
	integer(8),intent(OUT) :: sumofallops
	
	integer :: ibundle
	integer :: ifrag,ffrag
	
    real(4) :: localwt
    integer(8) :: localjumps,maxlocaljumps
    integer iprocs
    integer :: nnodesoverbooked, nnodesexcess

    real(8) :: opjumpstorage,maxopjumpstorage
	
	
	type (bund),pointer :: mybundle(:)


    if(nbundles < 1)then
		print*,' oops have not set # of opbundles '
	end if		
	if(draft)then
		mybundle=>draft_opbundle
	else
		mybundle=>opbundle
	end if
	
	sumofallops = 0
	do ibundle = 1,nbundles		
		sumofallops = sumofallops + mybundle(ibundle)%nops
		
!		print*,ibundle,mybundle(ibundle)%pxend,mybundle(ibundle)%pxstart,mybundle(ibundle)%nxend,mybundle(ibundle)%nxstart
		
	end do
	
	totalops = 0.0
	frag2fragops = 0.0
	excess_memory_needed = 0.0
   do ibundle = 1,nbundles
     call analyze_opbundle(ibundle,.true.,opjumpstorage)
	 maxopjumpstorage=max(maxopjumpstorage,opjumpstorage)



     if(usewts .and. .not. readjumpstart)then

  	   localwt= mybundle(ibundle)%opwt
	   
     else
       localwt = 1.d0
     end if
!	 if(readjumpstart)localwt=mybundle(ibundle)%opwt
	 	 
     ifrag = mybundle(ibundle)%ifragment
     ffrag = mybundle(ibundle)%ffragment
!	 write(78,*)nob, draft_opbundle(nob)%nops,localwt
     frag2fragops(ifrag,ffrag) = frag2fragops(ifrag,ffrag) + dfloat(mybundle(ibundle)%nops)*localwt
!	 if(iproc==0)print*,draft_opbundle(nob)%nops, localwt, ' weighting '
     totalops = totalops + dfloat(mybundle(ibundle)%nops)*localwt
	 	 
   end do   
   if(iproc==0 .and. excess_memory_needed > 0.0)then
	   print*,'  !!!!!! 2 '
	   print*,' May need to increase memory for jumpbs by ',excess_memory_needed,' Gb or approx ', &
	   excess_memory_needed/maxjumpmemory_default,' MPI nodes '
	   print*,'  !!!!!! 2'
	   
	   write(logfile,*)'  !!!!!! 2'
	   write(logfile,*)' May need to increase memory for jumpbs by ',excess_memory_needed,' Gb or approx ', &
	   excess_memory_needed/maxjumpmemory_default,' MPI nodes '
	   write(logfile,*)'  !!!!!! 2'
	   
   end if
   
!   nopbundles = nob_draft
   avgopstore = real(nprocs * maxjumpmemory_default*1e9,kind=8)/ real(sumofallops,kind=8)
   avgwtopstore = real(nprocs * maxjumpmemory_default*1e9,kind=8)/ real(totalops,kind=8)
   return
   
end subroutine count_all_operations_new

!===============================================================
!
! diagnostic subroutine added in 7.4.9 to check for missing jumps
!
subroutine print_opbundle_distro
	use nodeinfo
	
	implicit none
	integer jproc  ! dummy
	integer ibundle
	
	if(iproc/=0)return
	
	open(unit=87,file='procjumpdistro.bigstick',status='unknown')
	
	do jproc = 0,nprocs-1  ! loop over dummy processes
		
		write(87,*)jproc,' process; bundles:  ',opbundlestart(jproc),opbundleend(jproc)
		do ibundle = opbundlestart(jproc),opbundleend(jproc)
			select case(opbundle(ibundle)%optype)
			case('PP')
  			    write(87,*)ibundle,opbundle(ibundle)%optype,opbundle(ibundle)%hchar,opbundle(ibundle)%pxstart, & 
				            opbundle(ibundle)%pxend,opbundle(ibundle)%nxstart,opbundle(ibundle)%nxend
			case('NN')
				write(87,*)ibundle,opbundle(ibundle)%optype,opbundle(ibundle)%hchar,opbundle(ibundle)%Nxstart, & 
								            opbundle(ibundle)%nxend		
			case('PN')
				write(87,*)ibundle,opbundle(ibundle)%optype,opbundle(ibundle)%hchar,opbundle(ibundle)%pxstart, & 
					opbundle(ibundle)%pxend,opbundle(ibundle)%nxstart,opbundle(ibundle)%nxend											
			end select
		
		end do
		write(87,*)
	end do ! jproc
	
	
	close(87)
	return
	
end subroutine print_opbundle_distro

!===============================================================
!
!  diagnostic subroutine added in 7.7.3
!
!  prints out # of operations between fragments
!
subroutine print_fragops
	use operation_stats
	use fragments
	use nodeinfo
	use menu_choices,only:menu_char
	use io
	implicit none
	integer :: ifrag,ffrag,fragcount
	
	real(8) :: avgopsperfrag,maxopsperfrag,minopsperfrag
	
	if(iproc /= 0)return
	
	avgopsperfrag = 0.0
	maxopsperfrag = 0.0
	minopsperfrag = frag2fragops(1,1)
	fragcount = 0
	open(unit=91,file='frag2fragops.bigstick',status='unknown')
	do ifrag = 1,nfragments
		do ffrag = 1,nfragments
			if(frag2fragops(ifrag,ffrag)==0.0)cycle
			fragcount = fragcount +1
			avgopsperfrag = avgopsperfrag + frag2fragops(ifrag,ffrag)
			maxopsperfrag = max(maxopsperfrag,frag2fragops(ifrag,ffrag))
			minopsperfrag = min(minopsperfrag,frag2fragops(ifrag,ffrag))
			write(91,*)frag2fragops(ifrag,ffrag)
			
		end do
		
	end do
	if(nproc > 1 .or. menu_char=='m')write(logfile,*)fragcount,' out of ',nfragments**2,' frag 2 frag have ops '
!	print*,' weighted average # of operations = ',avgopsperfrag/float(nfragments)**2
	if(nproc > 1 .or. menu_char=='m')write(logfile,*)' min/max weighted ops on fragments = ',minopsperfrag,maxopsperfrag
!	print*,' NOTE: weighted average depends on file timinginfo.bigstick !'
	close(91)
	
	return
	
end subroutine print_fragops

!=================================
!
! a subroutine to compute the size of the initial and final basis covered by an opbundle
!
!
! INPUT:  iob = index of opbundle
!
! OUTPUT:
!   ini_basis = length of the basis addressed by this opbundle for initial states
!   fin_basis = length of the basis addressed by this opbundle for final states
!
subroutine bundle_basis(iob,ini_basis,fin_basis)
	use sectors
	use nodeinfo
	implicit none
	integer,intent(in) :: iob
	integer(8),intent(out) :: ini_basis,fin_basis
	integer(8) :: npsdi,npsdf,nnsdi,nnsdf   !  # of initial, final proton/neutron SDs
	integer :: ps,ns
	integer :: sstart,sstop
	
	
	ini_basis=0
	fin_basis=0
	npsdi = 0
	nnsdi = 0
	npsdf = 0
	nnsdf = 0
	
	select case (opbundle(iob)%optype)
	
	
	case('PP','PP0','PPP')
	
	ps = opbundle(iob)%isector
	npsdi = xsd(1)%sector(ps)%nxsd
	ps = opbundle(iob)%fsector
	npsdf = xsd(1)%sector(ps)%nxsd
	
	if( opbundle(iob)%insector<= opbundle(iob)%fnsector)then
		sstart = opbundle(iob)%insector
		sstop  = opbundle(iob)%fnsector
	else
		sstart = opbundle(iob)%fnsector
		sstop  = opbundle(iob)%insector
	end if
	
	do ns = sstart, sstop
		nnsdi = nnsdi + xsd(2)%sector(ns)%nxsd
	end do
	nnsdf = nnsdi
	
	case('NN','NN0','NNN')
	
	ns = opbundle(iob)%insector
	nnsdi = xsd(2)%sector(ns)%nxsd
	ns = opbundle(iob)%fnsector
	nnsdf = xsd(2)%sector(ns)%nxsd
	
	if( opbundle(iob)%isector<= opbundle(iob)%fsector)then
		sstart = opbundle(iob)%isector
		sstop  = opbundle(iob)%fsector
	else
		sstart = opbundle(iob)%fsector
		sstop  = opbundle(iob)%isector
	end if
	
	do ps = sstart, sstop
		npsdi = npsdi + xsd(1)%sector(ps)%nxsd
	end do
	npsdf = npsdi
	
	case('PN','PN0','PPN','PNN')
	ps = opbundle(iob)%isector
	npsdi = xsd(1)%sector(ps)%nxsd
	ps = opbundle(iob)%fsector
	npsdf = xsd(1)%sector(ps)%nxsd	
	
	ns = opbundle(iob)%insector
	nnsdi = xsd(2)%sector(ns)%nxsd
	ns = opbundle(iob)%fnsector
	nnsdf = xsd(2)%sector(ns)%nxsd
	
	case('SPE','CEN')
	ps = opbundle(iob)%isector
	npsdi = xsd(1)%sector(ps)%nxsd
	ps = opbundle(iob)%fsector
	npsdf = xsd(1)%sector(ps)%nxsd	
	
	ns = opbundle(iob)%insector
	nnsdi = xsd(2)%sector(ns)%nxsd
	ns = opbundle(iob)%fnsector
	nnsdf = xsd(2)%sector(ns)%nxsd
	case default
	
	print*,' THIS SELECTION NOT AVAILABLE '
	print*,' ERROR IN bundle_basis '
	print*,iob,opbundle(iob)%optype
	stop
	
    end select
	
	ini_basis = npsdi*nnsdi
	fin_basis = npsdf*nnsdf
	
!	if(iproc==0)write(93,'(9i8)')iob,opbundle(iob)%isector,opbundle(iob)%fsector,opbundle(iob)%insector,opbundle(iob)%fnsector,& 
!	npsdi,nnsdi,npsdf,nnsdf
	
	
	return
	
end subroutine bundle_basis
	
!===============================================================
!
! ADDED in 7.9.3 for advanced distro 
!
! A GREEDY bundle is one whose operations take more memory that is
! available on average
!
!  NOTE: avgopstore =  nprocs * max jump memory storage / total # of ops
!      this is the avg memory AVAILALBE for storing jumps

subroutine greedy_bundles_old(nob_draft)
	use flagger
	use fragments
	use nodeinfo
	use io
	use opbundles
	use verbosity
	implicit none
	
	integer, intent(in) :: nob_draft
	
	integer :: iob

	integer :: ifrag,ffrag
	
	real(8) :: memrequired ,memrequiredperop  ! how much memory is used to store

	ngreedy = 0
	greedymemorytot = 0.0
	if(.not.allocated(greedymemoryf2f))allocate(greedymemoryf2f(nfragments,nfragments))
	greedymemoryf2f=0.0
    if(iproc==0)then
		print*,' '
		print*,' EFFICIENCY OF MEMORY PER OPERATION '
		print*,' avg storage per operation available '   ,avgopstore   ,' bytes'
		write(logfile,*)' avg storage per operation available '   ,avgopstore   ,' bytes'

	end if
	
    if ( .not.restrictjumps .or. nprocs ==1)return
	
	do iob = 1,nob_draft
		draft_opbundle(iob)%greedy = .false.  ! default
		
		memrequired = draft_opbundle(iob)%jumpstore
		memrequiredperop = memrequired/ real(draft_opbundle(iob)%nops,8)
		
!		print*,iob,memrequiredperop,avgopstore
		if(memrequiredperop > avgopstore)then
			
!			print*,' testing ',draft_opbundle(iob)%jumpstore, draft_opbundle(iob)%nops
			ngreedy  = ngreedy+1
			draft_opbundle(iob)%greedy = .true.
					
			ifrag = draft_opbundle(iob)%ifragment
			ffrag = draft_opbundle(iob)%ffragment
			greedymemorytot = greedymemorytot + memrequired
			greedymemoryf2f(ifrag,ffrag)= greedymemoryf2f(ifrag,ffrag)+memrequired
			
		end if

	end do

    if(iproc==0)then
		if(verbose_distro_report)then
			if(ngreedy> 0)then 
				print*,' Total amount of greedy opbundle memory = ',greedymemorytot*1e-9,' Gb from '
				print*,ngreedy,' greedy opbundles '
				print*,' This requires roughly ', int(greedymemorytot*1e-9/maxjumpmemory_default),' dedicated MPI processes '
				print*,' Note there should be about ',nprocs*maxjumpmemory_default,' Gb total across all processes '
			else
				print*,' No greedy opbundles '
			end if
		end if
		write(logfile,*)' Total amount of greedy opbundle memory = ',greedymemorytot*1e-9,' Gb'
		if(ngreedy> 0)then
			write(logfile,*)ngreedy,' greedy opbundles '
			write(logfile,*)' This requires roughly ', int(greedymemorytot*1e-9/maxjumpmemory_default),' dedicated MPI processes '
			
			write(logfile,*)' Note there should be about ',nprocs*maxjumpmemory_default,' Gb total across all processes '
		else
			write(logfile,*)' No greedy opbundles '
		end if
		print*,' '
	end if

	return
	
end subroutine greedy_bundles_old
	
!===============================================================
!===============================================================
!
! ADDED in 7.9.3 for advanced distro 
!
! A GREEDY bundle is one whose operations take more memory that is
! available on average
!
!  NOTE: avgopstore =  nprocs * max jump memory storage / total # of ops
!      this is the avg memory AVAILALBE for storing jumps
!
!  CALLED BY:
!   draft_opbundles

subroutine greedy_bundles(nob_draft)
	use flagger
	use fragments
	use nodeinfo
	use io
	use opbundles
	use BMPI_mod
	use verbosity
	implicit none
	
	integer, intent(in) :: nob_draft
	
	integer :: iob

	integer :: ifrag,ffrag
	
	real(8) :: memrequired ,memrequiredperop  ! how much memory is used to store
	real(8) :: wtmemrequired ,wtmemrequiredperop  ! how much memory is used to store
	
	real(4) :: greedy_frac_thresh = 0.1  ! what fraction of memory triggers a greedy opbundle;
	           ! small greedy opbundles can be (?) ignored
    real(4) :: greedy_NN_frac_thresh= 0.5   ! special case for NN opbundles
	! this is because NN ops can require more memory because of sorting requirement, so harder to split
	logical :: greedy_verbose=.false.		   
	integer :: ierr
	
	real :: inflate  ! multiplicative factor for increasing weights, added 7.9.12

	ngreedy = 0
	greedymemorytot = 0.0
	if(.not.allocated(greedymemoryf2f))allocate(greedymemoryf2f(nfragments,nfragments))
	greedymemoryf2f=0.0
    if(iproc==0)then
		print*,' '
		print*,' EFFICIENCY OF MEMORY PER OPERATION '
		print*,' avg storage per operation available '   ,avgopstore   ,' bytes'
		print*,' avg storage per weighted operation available '   ,avgwtopstore   ,' bytes/ns'
		
		write(logfile,*)' avg storage per operation available '   ,avgopstore   ,' bytes'
		write(logfile,*)' avg storage per weighted operation avail: '   ,avgwtopstore   ,' bytes/ns'

	end if
	do iob = 1,nob_draft
		draft_opbundle(iob)%greedy=.false.
	end do
		
    if ( .not.restrictjumps .or. nprocs ==1)return
	
	do iob = 1,nob_draft
		draft_opbundle(iob)%greedy = .false.  ! default
		
		memrequired = draft_opbundle(iob)%jumpstore
		wtmemrequired = draft_opbundle(iob)%jumpstore*draft_opbundle(iob)%opwt
		
		memrequiredperop = memrequired/ real(draft_opbundle(iob)%nops,8)
		wtmemrequiredperop = wtmemrequired/ real(draft_opbundle(iob)%nops*draft_opbundle(iob)%opwt,8)

		if(( memrequiredperop > avgopstore .and. memrequired*1e-9 > greedy_frac_thresh*maxjumpmemory_default ) .or. & 
		( wtmemrequiredperop > avgwtopstore .and. memrequired*1e-9 > greedy_frac_thresh*maxjumpmemory_default ) )then

			ngreedy  = ngreedy+1
			draft_opbundle(iob)%greedy = .true.
					
			ifrag = draft_opbundle(iob)%ifragment
			ffrag = draft_opbundle(iob)%ffragment
			greedymemorytot = greedymemorytot + memrequired
			greedymemoryf2f(ifrag,ffrag)= greedymemoryf2f(ifrag,ffrag)+memrequired
			if(greedy_verbose .and. iproc==0)print*,' Greedy draft opbundle : ',iob
			
			if(use_greedy)then  ! modify the weights; added in 7.9.12
				inflate = memrequiredperop/ avgopstore
				inflate = max (inflate, wtmemrequiredperop / avgwtopstore)
				draft_opbundle(iob)%opwt = draft_opbundle(iob)%opwt*inflate
!				if( iproc==0)print*,' Inflating weight of greedy draft opbundle : ',inflate
				
			end if
			
		end if

	end do

    if(iproc==0)then
		if(verbose_distro_report)then
			if(ngreedy> 0 .and. greedy_verbose)then 
				print*,' Total amount of greedy opbundle memory = ',greedymemorytot*1e-9,' Gb from '
				print*,ngreedy,' greedy opbundles '
			
				print*,' This requires roughly ', int(greedymemorytot*1e-9/maxjumpmemory_default),' dedicated MPI processes '
			
				print*,' Note there should be about ',nprocs*maxjumpmemory_default,' Gb total across all processes '
			else
				print*,' No greedy opbundles '
			end if
		end if
		write(logfile,*)' Total amount of greedy opbundle memory = ',greedymemorytot*1e-9,' Gb'
		if(ngreedy> 0)then
			write(logfile,*)ngreedy,' greedy opbundles '
			write(logfile,*)' Note there should be about ',nprocs*maxjumpmemory_default,' Gb total across all processes '
		else
			write(logfile,*)' No greedy opbundles '
		end if
		print*,' '
	end if

#ifdef _MPI	
    call BMPI_Barrier(MPI_COMM_WORLD,ierr)  ! for debugging purposes only
#endif
	return
	
end subroutine greedy_bundles
	
!===============================================================
! added in 7.10.2 to try to diagnose two-body densities

subroutine count_opbundle_types(optype,hchar,numob,numops)
	
	character*3 :: optype
	character :: hchar
	
	integer ::numob,iob
	integer(8) :: numops
	real(8) :: opjumpstorage
	
	
	numob = 0
	numops = 0
	do iob = 1,nopbundles
		
		if(opbundle(iob)%optype ==optype .and. opbundle(iob)%hchar==hchar)then
			numob=numob+1
			call analyze_opbundle(iob,.false.,opjumpstorage)

			numops = numops + opbundle(iob)%nops
		end if
		
	end do
	return
	
end subroutine count_opbundle_types

!===============================================================
	

end module para_bundles_mod
