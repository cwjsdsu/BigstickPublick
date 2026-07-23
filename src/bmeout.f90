!
! ADDED April 2021 by CWJ, 7.9.12
!
! ROUTINES TO WRITE OUT uncoupled 1+2 body matrix elements to a file
! used rarely, but as a kludge for input into other codes, e.g., quantum computing
!


module mesout
	
contains
	
subroutine writeout_1and2bmes
    use system_parameters
    use io
    use menu_choices
    use verbosity
    use flags3body
    use nodeinfo
    use timing
    use btbme_mod
    use ntuple_info
	use diagh
	use flagger,only:subsume_SPE
	use sp_potter
	use hoppy
	implicit none
	integer :: itmax
	integer :: nme
	
	if(np(2)>1)then
		itmax=2   ! protons only
	else
		itmax=1   ! proton and neutrons
	end if
		
    call hopmaker  ! I appear to need this to generate the wfn
	
!
!----------------- RESTRICTIONS
! In order to save space, one needs to restrict the uncoupled matrix elements 
! created.  In order to do this, one needs the hops beforehand
  call master_q_ntuple(2,1)
  call master_q_ntuple(2,2)
  call master_cross_ntuple(1,1)   ! DON'T GET RID OF THIS... YET

!--------------- JUMPS------------------------------
!  NOTE: It is important to create 2-body jumps before 1-body
!  This is because certain linked lists are reused to plug memory leaks;
!  and an array contained in 2-body lists must be larger than for 1-body 
!
!  if(iproc==0)print*,' '

!----------------- COUNT UP JUMPS/OPERATIONS -------------
!........ SET UP NEEDED ARRAYS.....
!  if(iproc == 0)then
!     print*,' .... Setting up arrays for matrix elements ... '
!  end if
  call clocker('mes','sta')

!........ THIS WILL HAVE TO BE PUT ELSEWHERE....

  call prepare_to_uncoupleXXtbme(1)
  call prepare_to_uncoupleXXtbme(2)  
print*,' (BB)'
  
!  if(menu_char=='m')    call prepare_to_uncouplePNtbme(.false.)
   call prepare_to_uncouplePNtbme(.true.)

!  call prepare_to_uncouplePNtbme
  call clocker('mes','end')	
  
!  call threebodymaster
  call master_readin_tbmes(.false.)

 if(ham_readin)then
      call pandamaster
!     if(subsume_SPE)call subsume_sp_pot_into_2body
 end if

  
  call makespeme(1,'H')
  call makespeme(2,'H')
  call uncoupleXXtbme(1,'H')
  call uncoupleXXtbme(2,'H') 
  call delayedPNmatrixelements
  call count_uncouplepn('H',.true.)
  if(iproc /= 0)return
!.......... NOW TO WRITE OUT.....  
   open(unit=81,file='h12me.dat',status='unknown')	
   
   call write_sps2mesfile(81,itmax)
   call write_spe2mesfile(81,itmax,.true.,nme)
   call write_spe2mesfile(81,itmax,.false.,nme)
!.............. COUNT  UP TBMES..............   
   nme = 0
   call write_XXtbmes2file(81,itmax,.true.,nme)
   
   
!............... WRITE OUT TBMEs.............   
   write(81,'(i10)')nme
   nme = 0
   call write_XXtbmes2file(81,itmax,.false.,nme)

   close(81)
  return
	
end subroutine writeout_1and2bmes

!====================================================================
!
! subroutine to write out s.p.s information to the file h12me.dat
!
subroutine write_sps2mesfile(ifile,itmax)
	use spstate
	use haiku_info
	implicit none
	integer :: ifile 
	integer :: itmax
	integer :: it
	integer :: i,isps,nshift
	
	
	write(ifile,'(i3)')itmax  ! if =1, protons only, if = 2, protons +neutrons
!............... WRITE OUT # OF S.P. STATES............	
	do it = 1,itmax
		write(ifile,'(i5)')nsps(it)
	end do
	nshift = 0
	do it = 1,itmax
		do i = 1,nhsps(it)+nhsps(-it)
			if(i > nhsps(-it))then
				isps = i-nhsps(-it)
				
			else
				isps=i
				
			end if

			write(ifile,'(7i4)')isps+nshift,spsqn(it,isps)%orb,spsqn(it,isps)%nr, & 
			        spsqn(it,isps)%l,spsqn(it,isps)%j,spsqn(it,isps)%m,3-2*it
		end do ! isps
		nshift = nsps(it)
	end do ! it
	return
	
end subroutine write_sps2mesfile
!====================================================================
!
! subroutine to write out single particle potentials to the file h12me.dat
!
!  the potentials are stored by orbitals (psppot and nsppot) 
!  this routine writes them out explicitly in single-particle states (m-scheme)
!
subroutine write_spe2mesfile(ifile,itmax,count,nme)
	use spstate
	use coupledmatrixelements
	use sporbit
	implicit none
	integer :: ifile 
	integer :: itmax
	integer :: it
	integer :: isps,jsps
	integer :: iorb,jorb
	integer :: nshift
	logical count
	integer :: nme
	real,pointer :: potorb(:,:)
	nshift = 0
	if(.not.count)then
		write(ifile,'(i6)')nme  ! write how many matrix elements to expect
	end if
	
	nme = 0
	do it=1,itmax
		if(it==1)then
			potorb => psppot

		else
			potorb => nsppot
			
		end if
		do isps = 1,nsps(it)
			iorb = spsqn(it,isps)%orb
			do jsps = isps,nsps(it)
				jorb=spsqn(it,jsps)%orb

				if(potorb(iorb,jorb)==0.0)cycle
				nme = nme +1
				if(.not.count)then
					write(ifile,'(2i5,f10.5)')isps+nshift,jsps+nshift,potorb(iorb,jorb)
				end if
			end do
		end do
		nshift = nsps(1)
	end do
	return
		
end subroutine write_spe2mesfile
!====================================================================
!
! subroutine to write out uncoupled two-body matrix elements  to the file h12me.dat
! based upon uncoupleXXtbme in file btbme.f90
!
subroutine write_XXtbmes2file(ifile,itmax,count,nme)
	use spstate
	use interaction
	use coupledmatrixelements
	use haiku_info

	use btbme_mod
	implicit none
	integer :: ifile 
	integer :: itmax
	integer :: nme
	logical :: count
	integer :: it

    type (vjs), pointer :: tbme(:)
    integer pcref, pcstart, ppar
    integer,pointer :: xxmap(:)
    real,pointer :: v(:)
	integer :: nshift
	integer:: cpair,dpair
	integer :: m, par
	integer :: iref,istart,itbme
	integer :: i,j,k,l
	integer :: isps,jsps,ksps,lsps
	integer :: kth,lth
	real :: hfactor
	nshift = 0
	
	it =1

	do it=1, itmax
        if(it == 1)then
          tbme => ppme
          xxmap => ppcouplemap
          v => hmatpp
		  

        else
          tbme =>   nnme
          xxmap => nncouplemap
          v => hmatnn

        endif
        do dpair = 1,npairXX(it)        ! loop over destruction pair
			

          m = XX2(it)%pair(dpair)%m     ! find Jz, parity of pair
		
          par = XX2(it)%pair(dpair)%par

          iref = XX2(it)%meref(m,par)      ! iref, istart used to encode index for 
          istart = XX2(it)%mestart(m,par)  ! the uncoupled matrix element


  !---------------------- GET Q#s FOR PARTICLES IN DESTRUCTION PAIR --------
          k = XX2(it)%pair(dpair)%ia
          l = XX2(it)%pair(dpair)%ib
          if( k < 0)then
            ksps = -k
			kth = -it
          else
            ksps = k+nhsps(-it)
			kth = it
          endif

          if( l < 0)then
            lsps = -l
			lth = -it
          else
            lsps = l+nhsps(-it)
			lth = it
          endif
		  
!		  print*,m,hspsqn(kth,ksps)%m,hspsqn(lth,lsps)%m

          do cpair = 1,dpair        ! loop over creation pair  ! note triangular storage
            if(cpair == dpair)then
                hfactor = 0.5
            else
                hfactor = 1.0
            endif
            if(XX2(it)%pair(cpair)%m /=m)cycle       ! enforce quantum numbers
            if(XX2(it)%pair(cpair)%par /= par)cycle

  !--------------- HOW INDEX FOR UNCOUPLED (m-scheme) tbme IS COMPUTED---------

            itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref
			if(v(itbme)==0.)cycle


  !--------------- GET Q#s FOR CREATION PAIRS ---------------

            i = XX2(it)%pair(cpair)%ia
            j = XX2(it)%pair(cpair)%ib

            if( i < 0)then
              isps = -i
            else
              isps = i+nhsps(-it)
            endif

            if( j < 0)then
              jsps = -j
            else
              jsps = j+nhsps(-it)
            endif
			
			nme =nme+1
            if(.not.count)write(ifile,'(4i5,f12.5)')isps+nshift,jsps+nshift,ksps+nshift,lsps+nshift,v(itbme)

          enddo  ! loop over cpair

        enddo  ! loop over dpair
		nshift = nsps(1)
		
	end do ! loop over it
	
	return
end subroutine write_XXtbmes2file

!	
end module mesout