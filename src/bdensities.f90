!
! coupled one-body density matrices calculated FROM wavefunctions
!  
module densities
  use precisions
  implicit none

  real(kind=obs_prec), allocatable,target   :: p1bopme(:,:), n1bopme(:,:)
  logical :: pndensities    ! flag to denote printing out densities in proton-neutron formalism
                            ! added March 2017 to 7.7.5 
							
!.... ADDED in 7.8.2 to make using time-reversal easier, also writing out

 type densitypack
	 integer :: Jmin,Jmax
	 logical,allocatable :: zeroflag(:)
	 logical :: filled

	 real(kind=obs_prec),allocatable :: denmat(:,:,:,:) ! denmat(J,a,b,T)
	 
 end type densitypack
 
! logical :: usedensitybag = .true.     ! IF TRUE then use new density master routines
 										! automatically uses symmetry and writes out
										! DEPRECATED in 7.10.9
 logical, parameter :: usesymmetry = .false.   ! signals to only do i->f for i >= f FOR OLD ROUTINE
 logical :: writeunformatted1bden = .false.   ! if TRUE then write out as unformatted TO BE ADDED
  
 type (densitypack), allocatable,target :: densitybag(:,:)   ! densitybag(istate,fstate) ! note this is slightly backwards
 
 							
   contains
!
! Master routine for coupling density matrices
!
! from Edmonds "Angular momentum in quantum mechanics," Eq. (5.4.1)
!
!(Jf || O_K || Ji) = [Jf] ( Ji Mi, KM | Jf Mf)^-1 ( Jf Mf | O_KM | Ji Mi)
!
!For the density matrices I used the definition e.g. of Haxton and Donnelly:

!rho_K(a^+b) =  [K]^-1  (Jf || (a^+ b)_K || Ji)
!
!the advantage of  this definition is that the expectation value is just 
!the density matrix x the reduced matrix elements.
!
!Special case: number operator:  One can easily see that 

!rho(a^+a) =  sqrt( 2 Jf +1 )/ sqrt(2 j_a + 1) n(a)

!If one includes isospin then

!rho(a^+a) =  [ sqrt(2 Jf +1) sqrt(2 Tf +1)] / sqrt(2*(2 j_a +1 )) * n(a)
!
!
! NOTE: To save time/space we usually (starting in 7.6.0) compute i->f with i>=f
! and can use the symmetry
!
!  rho_K(a^+b,Ji->Jf) = rho_K(b^+a, Jf->Ji) x (-1)^(ja-jb + Ji-Jf)
!                       if isospin included, then also factor of (Ti-Tf)

!======================================================
! rewritten in 7.3.8 to allow for different particle spaces for protons, neutrons
!
!  CALLED BY
!  density1b_from_oldwfn
!  density1b_output 

  subroutine coupled_densities(Jt,Ji,Ti,Jf,Tf,Jzz,Tz,zeroflag, ndim,denmat)
  use system_parameters
  use sporbit
  use spstate
!  use densities
  use haiku_info
  use precisions
  use nodeinfo
  use bmpi_mod
  
  implicit none
  integer Jt, Tt      ! J, T of transition operator
  integer Ji,Ti       ! initial wfn J, T 
  integer Jf,Tf       ! final wfn J, T
  integer Jzz          ! of wfns
  integer Tz          ! = (Z -N )/2
  integer it
  logical zeroflag     ! flag to indicate no matrix elements possible
  integer:: ndim        ! dimensions of coupled density matrices
  real(8) denmat(ndim,ndim,0:1)
  real cleb
  real(8) cg,cgt
  integer a,b,ja,ma,jb,mb
  integer :: mya,myb
  integer ia,ib
  integer asps,ath,bsps,bth,asgn,bsgn
  integer iphase,tsign
  real, parameter :: cgtol = 5.0e-5 ! 1.0e-2 ! changed in 7.10.8 due to tiny clebsch (7 1, 10 0 | 7 1) previously 5.0e-6
  logical altflag
 real(kind=obs_prec), pointer :: x1bopme(:,:)
 integer :: ierr
 !....................................................................


  denmat(:,:,:) = 0.0

  cg = cleb(Ji,Jzz,2*Jt,0,Jf,Jzz)
  cg = cg*sqrt(float(2*Jt+1)) /sqrt(float( (jf+1)))
  if(abs(cg) < cgtol)then
     zeroflag = .true.
     return
  else
    zeroflag = .false.
  endif
do it = 1,2
  if(npeff(it) ==0)cycle
  if(it==1)then
      x1bopme => p1bopme
      tsign = 1
  else
      x1bopme => n1bopme
      tsign =-1
  end if
  do asps = 1, nhsps(it)+nhsps(-it)
     if(asps > nhsps(-it))then
         ia = asps - nhsps(-it)
         ath = it
         asgn = 1
     else
         ia = asps
         ath = -it
         asgn = -1
     endif

     a = hspsqn(ath,ia)%orb
     ja= hspsqn(ath,ia)%j
     ma= hspsqn(ath,ia)%m
     do bsps = 1, nhsps(it)+nhsps(-it)
        if(bsps > nhsps(-it))then
           ib = bsps - nhsps(-it)
            bth = it
           bsgn = 1
        else
           ib = bsps
           bth = -it
           bsgn = -1
        endif
        b = hspsqn(bth,ib)%orb
		
!		if(phconj(it))then
!			mya=b
!			myb=a
!		else
			mya=a
			myb=b
!		end if

        jb= hspsqn(bth,ib)%j
        mb= hspsqn(bth,ib)%m
        if( ma /= mb) cycle
        if(Jt > (ja + jb)/2 )cycle
        if(Jt < abs(ja - jb)/2 ) cycle
        iphase = (-1)**( (jb -mb)/2)
        if(isoflag .and. .not. pndensities)then
           do tt = 0,1
               cgt = cleb(Ti,Tz,2*Tt,0,Tf,Tz)
			   cgt = cgt* sqrt(float(2*Tt+1)/float(Tf+1))
               if(abs(cgt) < cgtol)cycle
        denmat(mya,myb,tt) = denmat(mya,myb,tt) + real(cleb(ja,ma,jb,-mb,2*Jt,0)*iphase,8) & 
            * tsign* x1bopme(ia*asgn,ib*bsgn)*real(cleb(1,tsign,1,-tsign,2*Tt,0),8) /(cgt)
           end do
         else
            denmat(mya,myb,it-1) = denmat(mya,myb,it-1)+ real(cleb(ja,ma,jb,-mb,2*Jt,0)*iphase,8) & 
            *  x1bopme(ia*asgn,ib*bsgn)
!			print*,it,x1bopme(ia*asgn,ib*bsgn),ja,ma,jb,mb
         endif

      end do  ! bsps
   enddo  ! asps
end do
denmat =denmat/cg
   return
  end subroutine coupled_densities
  
!=====================================================      
end module densities
