!
! library of routines for reading in single-particle potential
! added in 7.8.0 August 2017  by CWJ @ SDSU
!

module sp_potter
	use coupledmatrixelements
	use sporbit
	use system_parameters
	implicit none
	
	
contains
	
!......... CHECK IF SUBSUME...NEEDED TO SET UP BUNDLES....
!
!  CALLED BY: main routine
!

   subroutine do_I_subsume_spes
	use flagger,only: subsume_spe,spoccflag
	implicit none
	   
	call_spe = .true.
	
	if(.not.subsume_spe .or. spoccflag)return
	
	if(np(1)+np(2)<=1)then
		return
	end if
	call_spe = .false.	   
	return
	   
   end subroutine do_I_subsume_spes	
	
!...... SUBSUME S.P. ENERGIES + POTENTIALS INTO 2-BODY 
!       if there are more than 1 particle
	
	subroutine subsume_sp_pot_into_2body
!		use btbme_mod
        use butil_mod
		use flagger,only: subsume_spe,spoccflag
		implicit none
		integer :: a,b,c,d
		integer :: dfinish
		real    :: Afactor
		integer :: pair1,pair2,indx,tmp
		integer :: iref,istart, ipar
		integer :: jj,jjmin,jjmax
		integer :: ja,jb,jc,jd
		real :: vvv
		integer :: kdelta
		integer :: phasecd,phaseab
		
		if(call_spe)return
		
		Afactor = 1.0/float(np(1)+np(2)-1)
		
		if(npeff(2) > 0)then
			do a = 1,numorb(2)
				nsppot(a,a)=nsppot(a,a)+nspe(a)
			end do
		end if
		if(npeff(2) > 1)then
			do a = 1,numorb(2)
				ja = (orbqn(2,a)%j+1)/2
				do c = 1,a   ! ADD TO NN  
					jc = (orbqn(2,c)%j+1)/2
				    do b = 1,a
					    jb = (orbqn(2,b)%j+1)/2

						
						if(a==c)then
							dfinish = b
						else
							dfinish = c
						end if
						do d = 1,dfinish
							jd = (orbqn(2,d)%j+1)/2
							
							if(b/=d .and. a/=c .and. b/=c .and. a/=d)cycle

						
						! NOTE I do not think I have this looping correctly; I think I am leaving out 
						! stuff
						
				            pair1 = a*(a-1)/2 + b 
				            pair2 = c*(c-1)/2 + d 
						    pair1 = NNcouplemap(pair1)
						    pair2 = NNcouplemap(pair2)
						    if(pair1==-1 .or. pair2==-1)cycle
				            ipar = XXcouples(2)%pairc(pair1)%par 
				            if(ipar /= XXcouples(2)%pairc(pair2)%par)cycle
!				               print*,' problem with parity, boss (NN sp spot) ',a,b,c,b 
!				               stop 
!				            endif 
				            if(pair1 < pair2)then 
				              tmp = pair1 
				              pair1 = pair2 
				              pair2 = tmp 
				            endif 
				            iref = XXcouples(2)%meref(ipar) 
				            istart = XXcouples(2)%mestart(ipar)  
      
				            indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart 
						    jjmin = nnme(indx)%jmin
						    jjmax = nnme(indx)%jmax

						    do jj = jjmin,jjmax
							   phaseab = (-1)**( (ja-jb)+jj)
							   phasecd = (-1)**(( jc-jd)+jj)
						       vvv=nsppot(a,c)*kron_delta(b,d)+ nsppot(a,d)*kron_delta(b,c)*phasecd & 
							+   nsppot(b,c)*kron_delta(a,d)*phaseab + nsppot(b,d)*kron_delta(a,c)
						       vvv = Afactor/zeta(a,b)/zeta(d,c)*vvv
							   nnme(indx)%v(jj)= nnme(indx)%v(jj)+ vvv
						    end do   		! JJ
						end do ! d			
					
					end do ! b
										
				end do ! c
				
			end do ! a
		
		end if			
				
		if(npeff(1) > 0)then
			
			do a = 1,numorb(1)
				psppot(a,a)=psppot(a,a)+pspe(a)
			end do
			
			do a = 1,numorb(1)
				ja = (orbqn(1,a)%j+1)/2
				do c = 1,a   ! ADD TO PP  
					jc = (orbqn(1,c)%j+1)/2
					if(npeff(1)>1)then
				    do b = 1,a
					    jb = (orbqn(1,b)%j+1)/2

						
						if(a==c)then
							dfinish = b
						else
							dfinish = c
						end if
						do d = 1,dfinish
							jd = (orbqn(1,d)%j+1)/2
							
							if(b/=d .and. a/=c .and. b/=c .and. a/=d)cycle

						
						! NOTE I do not think I have this looping correctly; I think I am leaving out 
						! stuff
						
				            pair1 = a*(a-1)/2 + b 
				            pair2 = c*(c-1)/2 + d 
						    pair1 = PPcouplemap(pair1)
						    pair2 = PPcouplemap(pair2)
						    if(pair1==-1 .or. pair2==-1)cycle
				            ipar = XXcouples(1)%pairc(pair1)%par 
				            if(ipar /= XXcouples(1)%pairc(pair2)%par)cycle
!				               print*,' problem with parity, boss (PP) ',a,b,c,b 
!				               stop 
!				            endif 
				            if(pair1 < pair2)then 
				              tmp = pair1 
				              pair1 = pair2 
				              pair2 = tmp 
				            endif 
				            iref = XXcouples(1)%meref(ipar) 
				            istart = XXcouples(1)%mestart(ipar)  
      
				            indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart 
						    jjmin = ppme(indx)%jmin
						    jjmax = ppme(indx)%jmax

						    do jj = jjmin,jjmax
							   phaseab = (-1)**( (ja-jb)+jj)
							   phasecd = (-1)**(( jc-jd)+jj)
						       vvv=psppot(a,c)*kron_delta(b,d)+ psppot(a,d)*kron_delta(b,c)*phasecd & 
							+   psppot(b,c)*kron_delta(a,d)*phaseab + psppot(b,d)*kron_delta(a,c)
							   if(vvv==0.00)cycle
						       vvv = Afactor/zeta(a,b)/zeta(d,c)*vvv
						!	   print*,psppot(a,c)*kron_delta(b,d), psppot(a,d)*kron_delta(b,c)*phasecd & 
						!	,  psppot(b,c)*kron_delta(a,d)*phaseab , psppot(b,d)*kron_delta(a,c)
						!	   print*,a,b,c,d,vvv
							   ppme(indx)%v(jj)= ppme(indx)%v(jj)+ vvv
						    end do   		! JJ
						end do ! d			
					
					end do ! b
  				    end if
					
				    if(npeff(2) < 1)cycle
					
				
				    do b = 1,numorb(2)  ! ADD TO PN						
						if(a==c)then
							dfinish = b
						else
							dfinish = numorb(2)
						end if
						do d = 1,dfinish
					
						   pair1 = numorb(2)*(a-1) + b  						   
						   pair2 = numorb(2)*(c-1) + d
						   pair1 = PNcouplemap(pair1)
						   pair2 = PNcouplemap(pair2)
						   if(pair1==-1 .or. pair2==-1)cycle
						   ipar = PNcouples%pairc(pair1)%par
			               if(ipar /= PNcouples%pairc(pair2)%par)cycle
						   
				           if(pair1 < pair2)then 
				              tmp = pair1 
				              pair1 = pair2 
				              pair2 = tmp 
				           endif 
						   iref = PNcouples%meref(ipar) 
						   istart = PNcouples%mestart(ipar) 
						   indx = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart  
						   vvv = psppot(a,c)*kron_delta(b,d)+nsppot(b,d)*kron_delta(a,c)
						   vvv = Afactor*vvv!/zeta(a,b)/zeta(c,d)
						
						   jjmin = pnme(indx)%jmin
						   jjmax = pnme(indx)%jmax						
						   do jj = jjmin,jjmax
							  pnme(indx)%v(jj)= pnme(indx)%v(jj)+ vvv
						   end do   	! JJ
					   end do  ! d
					end do  ! b
				end do ! c
			end do ! a
					
		end if
		
		return
	end subroutine subsume_sp_pot_into_2body

!.......... Kronecker Delta function..........	
	integer function kron_delta(i,j)
		integer, intent(in) :: i,j
		kron_delta=0
		if(i==j)kron_delta=1
		return
	end function kron_delta

	
end module sp_potter