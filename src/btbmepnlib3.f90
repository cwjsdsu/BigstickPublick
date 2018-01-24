!
!  file btbmepnlib3.f90     started 7/2017 by CWJ @ SDSU in V7.7.8
!
!  new storage of uncoupled pn matrix elements
!
!  the central problem is as follows. In uncoupling any matrix elements, 
!  we have creation and annihilation operators, either in pairs or (for 3-body) triplets.
!  We can exactly separate into blocks pairs by Jz and parity.  
!  We have, however, an additional constraint on W, the weighting. This is not conserved 
!  but is an upper limit. So if we have a creation pair, a^+ b^+, there is a maximal Wa+Wb.
!
!  In proton-neutron interactions, if we have mixed proton-neutron creation/destruction pairs, 
!  e.g., ap^+ bn^+, this still works. This makes for a complex fetching of the matrix element, 
!  however, as the one-body jumps are of the form ap^+ cp, a one-body density operator.
!
!  To appropriately reduce the number of stored matrix elements, we still need the constraint 
!  of a maximal W on the creation/annihilation pairs. To do this, we try the following:
!  given a neutron density operator bn^+ dn, compute max (Wb,Wd). This gives a constraint on the 
!  max (Wa,Wc) on any proton density operator ap^+ cp we can couple to. This takes the place of 
!  double constraint on (Wa+Wb), and on (Wc+Wd). 
!
! 
!  
!  (Old storage: create pn pairs, sorted by Jz and parity, with max W)
!  (             to retrieve uncoupled matrix elements, from each 1-body jump, extract 
!                the creation and destruction labels, e.g., pc and pd from protons,
!                nc, nd from neutrons, and then combine pc and nc to get a pn pair,
!                as well as pd,nd, and from those pairs look up the matrix element.)
!
!  New storage: create proton and neutron density operators (denop), that is, a^+_i a_j. 
!               for each denop, compute the max and min W, that is, max (Wi,Wj) and min(Wi,Wj)
!               These are used for partial selection rules, as max Wp + max Wn must be capped.
!               Denops are, like pn pairs, sorted by Jz and parity. Within those we sort, 
!               first on descending maxW, and then within each of those, descending minW. 
!
!  Another issue of storage is one of time reversal. (Actually there are multiple time-reversal issues)
!  Namely, there is no need to store both V(ab,cd) and V(cd,ab). For pn densities, we rewrite 
!  V_J(ab,cd) implicitly as W_K(ac,bd) where a,c are proton labels and b,d neutron labels, so that
!  \sum_abcd V_j(ab,cd) a^+ b+ d c  becomes W_K(ac,bd) (a^+ c) (b^+d). Time reversal symmetry means 
!  we shouldn't need to store both W_K(ac,bd) and W_K(ca,db). Can we do this without resorting 
!  to IF statements?
! 
!  If we order proton densities, a >= c  and then allow all possible b,d for neutrons (that is, 
!  both b,d and d,b) then we halve our storage requirements. 
! 
!  Related to this are issues of parallelization. Our aim is to store on a node only the decoupled 
!  matrix elements we need, or at least not much more than that. In that case, if we have a change 
!  in Jz this can already account for time reversal. The main issue is when we have delta Jz = 0, 
!  in which case we tend to ned a lot of uncouple matrix elements.  
!
!  Some of this information we can extract from "generations" used to control the creation of jumps.
!  This has information on the "groups", which are sets of creation/annihilation operators
!  combined by Jz, parity, and W.  Furthermore, we can find out what groups are combined together, 
!  either in pairs or in denops.  This will help limit the pairs/denops we need and thus the matrix 
!  element storage. Using max W as a further cap may help.


module newXYtbme_mod
	use pairdef
	implicit none
	integer :: nX1bdenops(2)
	type (pairinfo), allocatable, target :: P1bdenop(:),N1bdenop(:)
	
	
	contains    ! subroutines
	
	subroutine count_create_pndenops
		
		
	end subroutine count_create_pndenops
	
end module newXYtbme_mod
	
	
	