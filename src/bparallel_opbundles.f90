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
   use coupledmatrixelements,only:call_spe
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
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PP0','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'NN0','f',nob_draft)
            call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PN0','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PN0','h',nob_draft)

           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PPP','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'NNN','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PPN','f',nob_draft)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.false.,'PNN','f',nob_draft)
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
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PP0','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'NN0','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PN0','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PN0','h',nob)

           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PPP','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'NNN','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PPN','f',nob)
           call count_create_draft_opbundles_g(.true.,ifrag,ffrag,.true.,'PNN','f',nob)
           if(ifrag <= ffrag)then
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
   do nob = 1,nob_draft
     call analyze_opbundle(nob,.true.,opjumpstorage)
	 maxopjumpstorage=max(maxopjumpstorage,opjumpstorage)
     if(usewts)then
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
       case('PNN')
          localwt = opwtPNN
       case default
          if(iproc == 0)print*,' OOPS ',nob,draft_opbundle(nob)%optype
          stop
       end select
     else
       localwt = 1.d0
     end if
     ifrag = draft_opbundle(nob)%ifragment
     ffrag = draft_opbundle(nob)%ffragment
!	 write(78,*)nob, draft_opbundle(nob)%nops,localwt
     frag2fragops(ifrag,ffrag) = frag2fragops(ifrag,ffrag) + dfloat(draft_opbundle(nob)%nops)*localwt
     totalops = totalops + dfloat(draft_opbundle(nob)%nops)*localwt
   end do

   if(iproc==0)then
	   if(maxopjumpstorage > 0.0)print*,' Max storage for jumps = ',maxopjumpstorage*1.0e-9,' Gb '
         if(usewts)then
              print*,' Estimated time per mat-vec multiply is ',totalops,' ns'
         else 
              print*,totalops,' operations  '
         end if
		 nopbundles = nob_draft
		 call count_all_operations(.true.,allmyops)
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
!  use opbundles
  use operation_stats
  use jumpNbody
  use jump3body
  use basis
  use nodeinfo
  use butil_mod
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

  bundle => bundummy  ! prevent compiler warnings
  if(draft .and. create)bundle => draft_opbundle
  if(.not.draft .and. create)bundle => opbundle

  do fs = fragmentlist(ffrag)%ssectorstart, fragmentlist(ffrag)%ssectorend   ! order switched in 7.6.1
     do is = fragmentlist(ifrag)%ssectorstart, fragmentlist(ifrag)%ssectorend

         if(is < fs )cycle   ! to prevent double counting; and this is the way
                               ! sector jumps are organized

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
                      exit
                   end if 

               end do

               if(xsj == -1)then
                  print*,' Problem cannot find sector jumps PP '
                  print*,is,fs
!                  print*,opsubstat(is,fs)%nopPP, opstat(iss,fss)%nopPP
                  print*,x2bjump(1)%nsectjumps
                  print*,xsd(1)%sector(is)%jzX, xsd(1)%sector(fs)%jzX
               do isj = 1,x2bjump(1)%nsectjumps
                     print*,isj,x2bjump(1)%isector(isj) ,x2bjump(1)%fsector(isj) , & 
      xsd(1)%sector(x2bjump(1)%isector(isj) )%jzX, xsd(1)%sector(x2bjump(1)%fsector(isj))%jzX
               end do
                  stop
               endif

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
!				if(iproc==0)print*,' test PP ',nob,is,fs,csi,csf

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
               end if
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
            insector = -1000002 
            fnsector = -1000001
            xstartx = -99
            xendx   = 0
            do xsj = 1,x2bjump(2)%nsectjumps
               
                ncs =  x2bjump(2)%csjmp(xsj)%ncjmps 
                okay = .false.
                do cs = 1,ncs
                  if( is == x2bjump(2)%csjmp(xsj)%cjump(cs)) then
                      okay =.true.
                      ysj = cs
                      exit
                   end if 
                end do
                if(.not.okay)cycle
                if(xstartx==-99)then
                   xstartx = x2bjump(2)%sjmp(xsj)%nstart+1
                   insector =x2bjump(2)%isector(xsj) !xsj
                end if
                xendx   = xendx+x2bjump(2)%sjmp(xsj)%njumps
                fnsector = x2bjump(2)%fsector(xsj) !xsj
             end do  ! xsj
             if(insector < 0) then
               print *, "insector < 0 in count_create_draft_opbundles_g"
               stop 1
             end if
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
!					if(iproc==0)write(95,*)nob,is,fs,insector,fnsector

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
             do cs = 1,ncs

                 ysj = x1bjump(1)%csjmp(xsj)%cjump(cs)
                 ystarty = x1bjump(2)%sjmp(ysj)%nstart+1
                 yendy   = ystarty+ x1bjump(2)%sjmp(ysj)%njumps-1
                 if(yendy < ystarty)cycle
                 localops = int(yendy-ystarty+1,8)*int(xendx-xstartx+1,8)

                 if( (hchar == 'h' .and.( is == fs .and. x1bjump(1)%sjmp(xsj)%diag & 
                                                  .and. x1bjump(2)%sjmp(ysj)%diag) )  & 
                    .or. ( (hchar == 'f' .or. hchar=='b') .and. (is /= fs  & 
           .or..not.x1bjump(1)%sjmp(xsj)%diag .or. .not. x1bjump(2)%sjmp(ysj)%diag)))then
                 nob = nob + 1
                 if(create)then

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
!   nopsnode: how many operations to put on (this) node
!   nob_draft: # of "draft" opbundles
!   nodestartob:  which opbundle to start on
!   nodestartz:   which set of operations to start on
!
! OUTPUT:
!   nodestopob: which obbundle to stop on
!   nodestopz : which set of operations to stop on
!
! CALLED BY:
!   distro_opbundles_over_fragments
!
subroutine search_draft_opbundles4split_g ( nob_draft, nodestartob, nodestopob, nodestartz, nodestopz)

!   use opbundles
   use nodeinfo
   use operation_stats
   implicit none
   integer :: nob_draft
!   integer(8) :: nopsnode
   integer :: nodestartob, nodestopob
   integer(8) :: nodestartz,nodestopz
   integer(8) :: opsofar,diff_nop,opsleft
!   integer(8) :: jumpsofar,diff_njump,jumpsleft

   integer(8) :: nset
   real (8) :: localwt

   integer iob

   opsofar = 0
!   jumpsofar = 0 
   if(nodestartob > nob_draft)then
       print*,' error in starting ',nodestartob,nob_draft
       stop
   end if
   do iob = nodestartob,nob_draft    ! LOOP OVER DRAFT OPBUNDLES
      if(draft_opbundle(iob)%pxstart < 0)then
           print*,iproc,iob,' pxstart ',draft_opbundle(iob)%pxstart
      end if
      if(draft_opbundle(iob)%pxend < 0)then
           print*,iproc,iob,' pxend ',draft_opbundle(iob)%pxend
      end if
      if(draft_opbundle(iob)%nxstart < 0)then
           print*,iproc,iob,' nxstart ',draft_opbundle(iob)%nxstart
      end if
      if(draft_opbundle(iob)%nxend < 0)then
           print*,iproc,iob,' nxend ',draft_opbundle(iob)%nxend
      end if

      if(draft_opbundle(iob)%nops == 0)cycle   ! HAS NO OPERATIONS, SKIP
      if(usewts)then
      select case (draft_opbundle(iob)%optype)
         case('PP ')
             localwt = opwtPP
             if (draft_opbundle(iob)%hchar == 'b') localwt = opwtPPb
         case('NN ')
             localwt = opwtNN
         case('PN ')
             localwt = opwtPN
             if (draft_opbundle(iob)%hchar == 'b') localwt = opwtPNb
         case('SPE')
             localwt = opwtSPE

         case('PPP')
             localwt = opwtPPP
         case('NNN')
             localwt = opwtNNN
         case('PPN')
             localwt = opwtPPN
         case('PNN')
             localwt = opwtPNN
         case default
          if(iproc == 0)print*,' OOPS bad optype ',draft_opbundle(iob)%optype
          stop

      end select

      else
        localwt = 1.d0
      endif
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
!          if(iob == nob_draft)then   ! CASE 1: must be at the end anyway
!             nodestopz  = draft_opbundle(iob)%nsetops    
!          end if
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
   print *, "Checking opbundles, only SPE for now"
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
   print *, "Checking opbundles, complete"
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
	   
	   case default
	   print*,' OH NOES I could not find the optype (analyze_opbundles_ops )', mybundle%optype
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

      case('PN0','PPN','PNN','PN','SPE')
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
	   case('PP0')
  	      localjumpstorage = real(int(mybundle%pxend+1,8)-int(mybundle%pxstart,8)+1,8)*bytesper2Bjump
   	   case('NN0')
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
    
	if(localjumpstorage > maxjumpmemory_default*1.0e9 )then  ! flag for possible splitting
		if(iproc == 0)then
			print*,' draft opbundle ',iob,' of type ',mybundle%optype,' requires memory ',localjumpstorage*1.0e-9, & 
			', greater than limit of ',maxjumpmemory_default
			if(nfragments > 1)print*,' Located on fragments ',mybundle%ifragment,' -> ',mybundle%ffragment
		end if
	end if
   
   return

end subroutine analyze_opbundle

!===============================================================

subroutine count_all_operations(draft,sumofallops)

	implicit none
	logical, intent(IN) :: draft  
	integer(8),intent(OUT) :: sumofallops
	
	integer :: ibundle
	
	type (bund),pointer :: mybundle(:)


    if(nopbundles < 1)then
		print*,' oops have not set # of opbundles '
	end if		
	if(draft)then
		mybundle=>draft_opbundle
	else
		mybundle=>opbundle
	end if
	
	sumofallops = 0
	do ibundle = 1,nopbundles
		sumofallops = sumofallops + mybundle(ibundle)%nops
!		print*,ibundle,mybundle(ibundle)%pxend,mybundle(ibundle)%pxstart,mybundle(ibundle)%nxend,mybundle(ibundle)%nxstart
		
	end do
	
	return	

end subroutine count_all_operations
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
	print*,fragcount,' out of ',nfragments**2,' frag 2 frag have ops '
	print*,' weighted average # of operations = ',avgopsperfrag/float(nfragments)**2
	print*,' min/max = ',minopsperfrag,maxopsperfrag
	print*,' NOTE: weighted average depends on file timinginfo.bigstick !'
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
	
	case default
	
	print*,' THIS SELECTION NOT AVAILABLE '
	print*,' ERROR IN bundle_basis '
	print*,iob,opbundle(iob)%optype
	stop
	
    end select
	
	ini_basis = npsdi*nnsdi
	fin_basis = npsdf*nnsdf
	
	if(iproc==0)write(93,'(5i8)')iob,opbundle(iob)%isector,opbundle(iob)%fsector,opbundle(iob)%insector,opbundle(iob)%fnsector,& 
	npsdi,nnsdi,npsdf,nnsdf
	
	
	return
	
end subroutine bundle_basis
	
	
	
	

end module para_bundles_mod
