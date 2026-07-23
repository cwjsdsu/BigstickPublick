!=========================================================================== 
!  Program BIGSTICK 
! 
!  Based upon code REDSTICK by W. E. Ormand 
!  BDEFORMEDLIB1.f90
!
!  optional subroutines for dealing with "deformed" systems, usually with
!  non-rotationally invariant single particle energies.
 
!=============================================== 
!
!  master routine for deformed (non-spherically-symmetric) single particle energies
!
! CALLED BY:
!    master_readin_tbmes
!
!  SUBROUTINES CALLED:
!    readdeformedSPEs
!
subroutine deformedSPEs 
 
use interaction_deformed
use spstate 
use haiku_info 
use io
 
implicit none 
integer it 
logical foundit 
integer :: aerr 
 
   allocate(speXm(-2:2, nhspsmax), stat=aerr) 
   if(aerr /= 0) call memerror("deformedSPEs") 
   speXm(:,:) = 0.0 
 
!   open(unit=22,file='deformed.spe',status='old',err=101) 
  write(logfile,*)' Reading in deformed single-particle energies '
  do it = 1,2 
    call readdeformedSPEs(it,foundit) 
    if(.not.foundit)then 
      print*,' Could not find file deformed.int ' 
      write(logfile,*)' Could not find file deformed.int ' 
    end if 
 
  end do 
 
return 
 
end subroutine deformedSPEs 
!=======================================================
!
!  INPUT:
!   it = species
!  OUTPUT
!   foundit = logical flag to signify success
!            in opening file deformed.int
!
!  CALLED BY:
!    deformedSPEs
!
subroutine readdeformedSPEs(it,foundit) 
   use verbosity 
   use spstate 
   use haiku_info 
   use system_parameters 
   use interaction_deformed
 
   implicit none 
   integer it      ! which species    
   logical foundit 
   integer ttmp, mtmp,jtmp,ntmp,ltmp 
   integer mtot,partot 
   integer, allocatable :: mapsps(:),mapm(:),mappar(:) 
   integer nspslocal 
   integer i 
   integer ntbme  ! # of uncoupled TBMEs in file 
   integer a,b,c,d 
   integer isps,jsps,ksps,lsps 
   integer cpair,dpair 
   real hfactor 
   real vxx 
   logical foundsps 
   integer iref, istart 
   integer cphase,dphase 
   integer itbme 
   real,pointer :: v(:) 
   integer,pointer :: map(:) 
   integer jzmin,jzmax 
   integer :: aerr 
 
   print*,' XX deformed ',it,nhspsmax 
   foundit = .false. 
   open(unit=23,file='deformed.int',status='old',err=101) 
   foundit = .true. 
!...... FIRST READ PAST AND GET MAPPING OF SINGLE PARTICLE STATES..... 
 
   read(23,*)nspslocal 
   allocate(mapsps(nspslocal),mapm(nspslocal),mappar(nspslocal), stat=aerr ) 
   if(aerr /= 0) then 
      call memerror("readdeformedSPEs") 
      stop 5 
   end if 
   mapsps(:) = -1 
   do i = 1,nspslocal 
      read(23,*)ttmp, ntmp, ltmp,jtmp,mtmp 
      if(ttmp /= it)then 
         mapsps(i) = -1 
         cycle 
      end if 
      mapm(i) = mtmp 
      mappar(i) = (-1)**ltmp 
      foundsps = .false. 
!............ SEARCH FOR SINGLE-PARTICLE STATES..... 
      if(mtmp > 0)then 
        do isps = 1,nhsps(it) 
          if( mtmp == hspsqn(it,isps)%m .and. jtmp == hspsqn(it,isps)%j  &  
          .and. ntmp == hspsqn(it,isps)%nr .and. ltmp == hspsqn(it,isps)%l )then 
              foundsps = .true. 
              mapsps(i) = isps+nhsps(-it) 
              exit 
          end if 
 
        end do 
 
      else 
        do isps = 1,nhsps(-it) 
          if( mtmp == hspsqn(-it,isps)%m .and. jtmp == hspsqn(-it,isps)%j  &  
          .and. ntmp == hspsqn(-it,isps)%nr .and. ltmp == hspsqn(-it,isps)%l )then 
              foundsps = .true. 
              mapsps(i) = isps 
              exit 
          end if 
 
        end do 
 
      end if 
      if(.not.foundsps)then   ! ERROR TRAP   
         print*,' DID NOT FIND CORRESPONDING STATE ' 
         print*,i,ntmp,ltmp,jtmp,mtmp 
         stop 
 
 
      end if 
 
   end do 
!................. READ SINGLE-PARTICLE ENERGIES  (added in V7.2.8).... 
 
   do a = 1,nspslocal 
      read(23,*)vxx 
      if(mapsps(a) < 0)cycle  ! wrong species of nucleon 
      isps = mapsps(a) 
      if(isps > nhsps(-it))then 
        speXm(-it,isps-nhsps(-it) ) = vxx 
      else 
        speXm(it,isps) = vxx 
      end if 
 
   end do 
 
   close(23) 
   deallocate(mapsps,mapm,mappar ) 
   return 
101 continue 
   print*,' file deformed.int does not exist ' 
   foundit = .false. 
   return 
end subroutine readdeformedSPEs 
!
!
!=======================================================
!
!  intermediate routine to create single-particle energies for haikus
!
!  alternate to routine makehaikuspes 
!
!  CALLED BY: 
!    makespeme
!
   subroutine makehaikuspes_deformed(ith,bspe)

   use bitstuff
   use diagh
   use blocks
   use haiku_info
   use spstate
   use interaction_deformed

   implicit none
   integer ith
   type (blockhspe), pointer :: bspe(:)
!   real, pointer :: xspe(:)

!------------------------------------
   integer iblock
   integer i
   integer isps
   integer j,iword,imask
   integer iorb
   integer, pointer :: hsd(:,:)
   integer :: aerr

   allocate(bspe( hblock(ith)%nhblocks ), stat=aerr)
   if(aerr /= 0) call memerror("makehaikuspes_deformed 1")
   do iblock = 1,hblock(ith)%nhblocks
     allocate( bspe( iblock)%hspe ( hblock(ith)%list(iblock)%nhsd ), stat=aerr)
     if(aerr /= 0) call memerror("makehaikuspes_deformed 10")
     bspe(iblock)%hspe = 0.0
     do i = 1,hblock(ith)%list(iblock)%nhsd
        hsd => hblock(ith)%list(iblock)%hsd
        do isps = 1,nhsps(ith)
          iword = (isps-1)/max_bit_word+1
          j = mod(isps-1,max_bit_word)+1 
          imask = ishft(1,j-1)
          if ( iand(imask,hsd(iword,i)) /=0 ) then
!              iorb = hspsqn(ith,isps)%orb
              bspe(iblock)%hspe(i) = bspe(iblock)%hspe(i) + speXm(ith,isps)
          end if

        enddo ! iop
     enddo ! i
   enddo  ! iblock
   return
   end subroutine makehaikuspes_deformed
!
!=======================================================
 
 
! 
! this version is mostly obsolete;  
! 
 
subroutine deformedSPEsOLD 
 
use interaction_deformed
use spstate 
use haiku_info 
 
implicit none 
integer ith 
integer isps 
integer m,j,l,nr 
integer i,jtmp,ltmp,mtmp,ttmp,nrtmp 
real spetmp 
integer :: aerr 
 
   allocate(speXm(-2:2, nhspsmax), stat=aerr) 
   if(aerr /= 0) call memerror("deformedSPEsOLD") 
   speXm(:,:) = 0.0 
 
   open(unit=22,file='deformed.spe',status='old',err=101) 
 
   do ith = -2,2 
       if(ith == 0)cycle 
       do isps = 1,nhsps(ith) 
          j = hspsqn(ith,isps)%j 
          m = hspsqn(ith,isps)%m 
          l = hspsqn(ith,isps)%l 
          nr= hspsqn(ith,isps)%nr 
          rewind(22) 
          do i=1,1000 
             read(22,*,end=102)ttmp,nrtmp, ltmp,jtmp,mtmp,spetmp 
             if(ttmp /= abs(ith))cycle 
             if(ltmp == l .and. jtmp == j .and. mtmp == m .and. nrtmp == nr)then 
                 speXm(ith,isps) = spetmp 
                 exit 
             end if 
 
          end do 
 
       end do 
 
   end do 
   close(22) 
 
return 
 
101  continue 
   print*,' Looking for file deformed.sps, does not exist ' 
   stop 
 
102 continue 
   print*,' ran out of states ' 
   print*,ith,l,j,m 
end subroutine deformedSPEsOLD 
