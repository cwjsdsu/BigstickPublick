!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  BIGSTICK configuration-interaction shell-model code
!
!  file BDIAGH.F
      module diagh

      implicit none
! diagonal interaction from single-particle energies

      type blockhspe
         real, pointer :: hspe(:)
      end type blockhspe

      real, pointer :: pspeME(:),nspeME(:)

      real, allocatable, target :: pspe_obs(:),nspe_obs(:)
      real, allocatable, target :: pspe_h(:),nspe_h(:)

contains
	
!
!  Diagonal part of Hamiltonian:
!  Hamiltonian routines for BIGSTICK
!
!  CALLED BY:
!      MAIN
!      lanczos_output
!      setup4obsmaster
!      pocc_compute_spocc
!      particle_occupation_p_orig
!
   subroutine makespeme(it,obschar)

!   use diagh
   use basis
   use coupledmatrixelements
   use interaction_deformed
   implicit none
!   include 'binterfaces.inc'

   integer it
   character(1) obschar
   real, pointer :: speME(:)
   real, pointer :: xspe(:)
   type (blockhspe), pointer :: lbspe(:), rbspe(:)
   integer :: aerr
   
   if(.not.call_spe)return

   select case (obschar)

   case ('H','J')

   if(it == 1)then
     if(.not.allocated(pspe_h)) then
        allocate(pspe_h( nXsd(it) ), stat=aerr)
        if(aerr /= 0) call memerror("makespeme 1")
     end if
     speME => pspe_h
     xspe  => pspe
   else
     if(.not.allocated(nspe_h)) then
        allocate(nspe_h( nXsd(it) ), stat=aerr)
        if(aerr /= 0) call memerror("makespeme 2")
     end if
     speME => nspe_h
     xspe  => nspe
   endif

   case ('T')
   if(it == 1)then
     if(.not.allocated(pspe_obs)) then
        allocate(pspe_obs( nXsd(it) ), stat=aerr)
        if(aerr /= 0) call memerror("makespeme 3")
     end if
     speME => pspe_obs
     xspe  => pspe
   else
     if(.not.allocated(nspe_obs)) then
        allocate(nspe_obs( nXsd(it) ), stat=aerr)
        if(aerr /= 0) call memerror("makespeme 4")
     end if
     speME => nspe_obs
     xspe  => nspe
   endif

   end select
   speME = 0.0
   if(deformed .and. obschar=='H')then
      call makehaikuspes_deformed(-it,lbspe)
      call makehaikuspes_deformed(it,rbspe)
   else
      call makehaikuspes(-it,lbspe,xspe)
      call makehaikuspes(it,rbspe,xspe)
   end if
   call combinehaikuspes(it,lbspe,rbspe,speME)
   deallocate(lbspe,rbspe)
   return
   end subroutine makespeme
!========================================================
!
!
   subroutine makehaikuspes(ith,bspe,xspe)

   use bitstuff
!   use diagh
   use blocks
   use haiku_info
   use spstate

   implicit none
   integer ith
   type (blockhspe), pointer :: bspe(:)
   real, pointer :: xspe(:)

!------------------------------------
   integer iblock
   integer i
   integer isps
   integer j,iword,imask
   integer iorb
   integer, pointer :: hsd(:,:)
   integer :: aerr

   allocate(bspe( hblock(ith)%nhblocks ), stat=aerr)
   if(aerr /= 0) call memerror("makehaikuspes 1")
   do iblock = 1,hblock(ith)%nhblocks
     allocate( bspe( iblock)%hspe ( hblock(ith)%list(iblock)%nhsd ), stat=aerr)
     if(aerr /= 0) call memerror("makehaikuspes 10")
     bspe(iblock)%hspe = 0.0
     do i = 1,hblock(ith)%list(iblock)%nhsd
        hsd => hblock(ith)%list(iblock)%hsd
        do isps = 1,nhsps(ith)
          iword = (isps-1)/max_bit_word+1
          j = mod(isps-1,max_bit_word)+1 
          imask = ishft(1,j-1)
          if ( iand(imask,hsd(iword,i)) /=0 ) then
              iorb = hspsqn(ith,isps)%orb
              bspe(iblock)%hspe(i) = bspe(iblock)%hspe(i) + xspe(iorb)
          end if

        enddo ! iop
     enddo ! i
   enddo  ! iblock
   return
   end subroutine makehaikuspes
!
!========================================================

   subroutine combinehaikuspes(it,lbspe,rbspe,speME)

   use sectors
!   use diagh
   use blocks

   implicit none
   integer it
   type (blockhspe), pointer :: lbspe(:),rbspe(:)
   real, pointer :: speME(:)

!------------------------------------
   integer is
   integer iblock
   integer rblock,lblock
   integer(8) :: xsdstart
   integer nradd,nladd
   integer ladd,radd
   integer(8) :: ix 
   do is = 1,nsectors(it) 
     do iblock = 1,xsd(it)%sector(is)%nhblocks
        rblock = xsd(it)%sector(is)%rhblock(iblock)
        lblock = xsd(it)%sector(is)%lhblock(iblock)
        xsdstart = xsd(it)%sector(is)%blockstart(iblock)-1

!---------------- always make left address the outermost loop (a convention)
        nradd = hblock(it)%list(rblock)%nhsd
        nladd = hblock(-it)%list(lblock)%nhsd

        do ladd = 1,nladd
            do radd = 1,nradd
                ix = xsdstart + radd+(ladd-1)*nradd
                speME(ix) = lbspe(lblock)%hspe(ladd)+rbspe(rblock)%hspe(radd)
            enddo  ! radd
        enddo  ! ladd

     enddo ! iblock

  enddo ! is
   
  return
  end subroutine combinehaikuspes
!=======================================================
end module diagh
