! routines for applying 2-body non-scalar operator

module trans2b
	implicit none
	
	integer :: J2btrans   ! J of 2-body transition operator
	integer :: T2btrans   ! optional T of 2-body transition operator
	

!........ FOR ANY GIVEN J, we can couple up Jab, Jcd to J	
    type gen2bop
		integer :: nvalues   ! combinations of Jab,Jcd
		integer, allocatable,target :: Jab(:),Jcd(:)
		real, allocatable, target :: op2bme(:)
	end type gen2bop	
	type (gen2bop), allocatable :: pp2bopme(:), nn2bopme(:),pn2bopme(:)


contains
!============================================================ 
! 
! set up storage of coupled PP,NN TBMEs for nonscalar operators
!  
!  started 5/2017 by CWJ @SDSU 
!  based upon routine set_tbme_XXarray in binput_tbme.f90 
! 
 
	subroutine set_trans_tbme_XXarray(it) 
 
	  use sporbit 
	  use coupledmatrixelements 
	  use butil_mod
	  implicit none 
 
	  integer it  ! species, 1 = P, 2 =N 
 
	  integer cpair,dpair,ddpair,ccpair 
	  integer pair1,pair2,tmp 
	  integer dparity 
	  integer iref,istart 
	  integer itbme 
	  integer a,b,c,d 
	  integer ja,jb,jc,jd 
	  integer jabmin,jabmax,jcdmin,jcdmax 
	  integer nv 
	  integer ncpair 
	  type (vjs), pointer :: xxme(:) 
	  integer, pointer :: cmap(:) 
	  integer :: aerr 
 
	!.................... COMPUTE # OF MATRIX ELEMENTS 
	  nv = 0 
	  ncpair = XXcouples(it)%meref(2) 
	  nv = ncpair*(ncpair+1)/2 
	  ncpair = ncouplesXX(it) - XXcouples(it)%meref(2) 
	  nv = nv + ncpair*(ncpair+1)/2 
	  nv2bmedim(it) = nv 
	  if(nv == 0)return 
	  if(it == 1)then 
	     allocate( ppme( nv2bmedim(it) ), stat=aerr) 
	     if(aerr /= 0) call memerror("set_tbme_XXarray 1") 
	     xxme => ppme 
	     cmap => PPcouplemap 
	  else 
	     allocate( nnme( nv2bmedim(it) ), stat=aerr) 
	     if(aerr /= 0) call memerror("set_tbme_XXarray 2") 
	     xxme => nnme 
	     cmap => NNcouplemap 
	  endif 
 
	  do dpair = 1,ncouplesXX(it) 
 
	     dparity = XXcouples(it)%pairc(dpair)%par 
	     iref = XXcouples(it)%meref(dparity) 
	     istart = XXcouples(it)%mestart(dparity) 
	     c = XXcouples(it)%pairc(dpair)%ia 
	     d = XXcouples(it)%pairc(dpair)%ib 
	     jc = orbqn(it,c)%j 
	     jd = orbqn(it,d)%j 
		 jcdmin = abs(jc-jd)/2
		 jcdmax = (jc+jd)/2
 
	     do cpair = 1, dpair 
 
	          if( XXcouples(it)%pairc(cpair)%par /= dparity )cycle 
	          a = XXcouples(it)%pairc(cpair)%ia 
	          b = XXcouples(it)%pairc(cpair)%ib 
 
	          ja = orbqn(it,a)%j 
	          jb = orbqn(it,b)%j 
	          jabmin =  abs(ja-jb)/2 
	          jabmax = (ja+jb)/2 
           
	          itbme = istart + (dpair-iref)*(dpair-iref-1)/2+cpair-iref 
	          xxme(itbme)%jmin = jmin 
	          xxme(itbme)%jmax = jmax 
	          if(jmax >= jmin) then 
 
	              allocate( xxme(itbme)%v(jmin:jmax), stat=aerr )  
	              if(aerr /= 0) call memerror("set_tbme_XXarray 3") 
	              xxme(itbme)%v(:) = 0.0 
 
	          endif       
	     end do  ! cpair 
 
	  end do ! dpair 
 
 
	  return 
	end subroutine set_trans_tbme_XXarray 
 
	!============================================================ 
end module dens2b
