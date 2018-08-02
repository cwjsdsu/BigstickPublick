!=============================================================

subroutine r2scalarop(t)

   use phonon1b
   use sporbit
   implicit none
   integer t
   integer i1,i2
   integer n1,n2,l1,l2,j
   real r2me
   real tj2i   ! 3-j symbol from libra.f
   real facs, jx
   real tau ! reduced matrix element of tau operator
   real bosc

   allocate ( t1bme( numorb(1), numorb(1) ) )

   t1bme(:,:) = 0.0

!  reduced matrix elements of tau operator
   if (t==0) then
      tau = sqrt(2.)
   else
      tau = sqrt(6.)
   endif
   
   print*,' Enter value of oscillator length b in fm (enter 1 if not desired )'
   read*,bosc

   do i1 = 1,numorb(1)
      j  = orbqn(1,i1)%j
      n1 = orbqn(1,i1)%nr
      l1 = orbqn(1,i1)%l
      do i2 = 1,numorb(1)
          if( j /= orbqn(1,i2)%j)cycle
          l2 = orbqn(1,i2)%l
          if(l1 /= l2)cycle
          n2 = orbqn(1,i2)%nr
          t1bme(i1,i2) = r2me(n1,l1,n2,l2)

          if(spinless)then
          t1bme(i1,i2) = t1bme(i1,i2)*(2*l1+1)/sqrt(4.*3.1415926)  & 
          *tj2i(l1,0,l1,0,0,0)*(-1)**(l1)

          else
!       old code        
!          t1bme(i1,i2) = t1bme(i1,i2)*(2*j+1)/sqrt(4.*3.1415926)  & 
!          *tj2i(j,0,j,-1,0,1)*(-1)**((j+3)/2)

!       MKGK code
!          j is really 2j in this code
          jx = j/2.0
          facs=sqrt((2*jx+1)/(4.0*3.1415926))*tau
          t1bme(i1,i2)=t1bme(i1,i2)*facs*(-1)**(2*jx+1)*bosc**2
!          write(6,*) 't1bme = ',i1,i2,t1bme(i1,i2), r2me(n1,l1,n2,l2),&
!          foobar, j,l1
          endif

      end do
   end do ! il

   return

end subroutine r2scalarop
!=============================================================

subroutine quadop(t)

   use phonon1b
   use sporbit
   implicit none
   integer i1,i2
   integer n1,n2,l1,l2,j,j2
   real r2me
   integer t
   real tfact
   real tj2i   ! 3-j symbol from libra.f; with integer input 2xj
   real sj2i   ! 6-j symbol from libra.f 
   logical crossshell,allowelliott
   real peffcharge, neffcharge
   real tmp1b
   real bosc
   
   character :: crosschar
   allocate ( t1bme( numorb(1), numorb(1) ) )
   allowelliott = .true.    ! flag for allowing elimination of cross-shell matrix elements
                          ! to produce Elliott's SU(3) operator
   t1bme(:,:) = 0.0
   tfact = sqrt(2.*(2*T+1))
   
   if(allowelliott)then
   print*,' Include cross-shell matrix elements (n; default is yes)'
   read*,crosschar
   if(crosschar=='N' .or. crosschar=='n')then
	   crossshell=.false.
	   print*,' No cross-shell matrix elements to be included '
   else
	   crossshell=.true.
	   print*,' Including all allowed matrix elements '
   end if
   else
	crossshell=.true.
   end if
   
   
   print*,' Enter value of oscillator length b in fm (enter 1 if not desired )'
   read*,bosc
   
   if(t==0 .or. t==1)then

   do i1 = 1,numorb(1)
      j  = orbqn(1,i1)%j
      n1 = orbqn(1,i1)%nr
      l1 = orbqn(1,i1)%l
      do i2 = 1,numorb(1)
          l2 = orbqn(1,i2)%l
          j2  = orbqn(1,i2)%j

          if(abs(l1- l2)> 2)cycle
          if(abs(j-j2) > 4)cycle
          if(abs(l1+l2) < 2)cycle
          if(abs(j+j2) < 4)cycle
          n2 = orbqn(1,i2)%nr
		  if(.not.crossshell .and. 2*n2+l2 /= 2*n1+l1)cycle
		  
          t1bme(i1,i2) = r2me(n1,l1,n2,l2)*sqrt(5.*(2*l1+1)*(2*l2+1)) & 
          /sqrt(4.*3.1415926) *tj2i(l1*2,4,l2*2,0,0,0)*(-1)**(l1)
          if(.not.spinless)then
          t1bme(i1,i2) = t1bme(i1,i2)*sqrt((j+1.)*(j2+1.))*(-1)**(l1+(j2+1)/2) & 
                         *sj2i(2*l1,j,1,j2,2*l2,4)*tfact*bosc**2

 !   print*,i1,i2,tj2i(l1*2,4,l2*2,0,0,0),sj2i(2*l1,j,1,j2,2*l2,4),r2me(n1,l1,n2,l2)

          endif
      end do
   end do ! il
   
else
	pnformal = .true.
	print*,' Enter effective charges for proton, neutron '
	read*,peffcharge,neffcharge
    allocate ( p1bme( numorb(1), numorb(1) ),  n1bme( numorb(2), numorb(2) ) )
    p1bme(:,:) = 0.0	
	n1bme(:,:) =0.0
    do i1 = 1,numorb(1)
       j  = orbqn(1,i1)%j
       n1 = orbqn(1,i1)%nr
       l1 = orbqn(1,i1)%l
       do i2 = 1,numorb(1)
           l2 = orbqn(1,i2)%l
           j2  = orbqn(1,i2)%j

           if(abs(l1- l2)> 2)cycle
           if(abs(j-j2) > 4)cycle
           if(abs(l1+l2) < 2)cycle
           if(abs(j+j2) < 4)cycle
           n2 = orbqn(1,i2)%nr
 		  if(.not.crossshell .and. 2*n2+l2 /= 2*n1+l1)cycle
		  
           tmp1b = r2me(n1,l1,n2,l2)*sqrt(5.*(2*l1+1)*(2*l2+1)) & 
           /sqrt(4.*3.1415926) *tj2i(l1*2,4,l2*2,0,0,0)*(-1)**(l1)
           tmp1b = tmp1b*sqrt((j+1.)*(j2+1.))*(-1)**(l1+(j2+1)/2) & 
                          *sj2i(2*l1,j,1,j2,2*l2,4)*bosc**2
		
 		  p1bme(i1,i2)= peffcharge*tmp1b  
 		  n1bme(i1,i2)= neffcharge*tmp1b  				  
						  
       end do
    end do ! il	
	
end if

   return

end subroutine quadop

!---------------------------------------------
!
! note phase changed in Version 2, Nov 2013
!
  subroutine spinop(t)

   use phonon1b
   use sporbit
   implicit none
   integer i1,i2
   integer n1,n2,l1,l2,j1,j2
   real sj2i   ! 6-j symbol from libra.f
   integer t
   real tfact
   
   pnformal = .false.
   if(t== 0 .or. t == 1)then
      allocate ( t1bme( numorb(1), numorb(1) ) )

      t1bme(:,:) = 0.0

      if(t ==0)then
          tfact =sqrt(2.)      ! because reduced < 1/2 || 1 || 1/2 > = sqrt(2)
      else
         tfact = sqrt(6.)     ! because reduced < 1/2 || tau || 1/2 > = sqrt(6)
      endif
   

   do i1 = 1,numorb(1)
      j1  = orbqn(1,i1)%j
      n1 = orbqn(1,i1)%nr
      l1 = orbqn(1,i1)%l
      do i2 = 1,numorb(1)
          l2 = orbqn(1,i2)%l
          if(l1 /= l2)cycle
           
          n2 = orbqn(1,i2)%nr
          if(n1/= n2)cycle
          j2  = orbqn(1,i2)%j
          t1bme(i1,i2) = sj2i(1,j1,2*l1,j2,1,2)*sqrt(6.*(j1+1.)*(j2+1.)) & 
              *(-1)**(( 2*l1 + 3 + j1)/2)*tfact   !phase here changed by -1
      end do
   end do ! il
   return
else
	pnformal = .true.
	allocate(p1bme(numorb(1),numorb(1)))
	allocate(n1bme(numorb(2),numorb(2)))
	p1bme = 0.0
	n1bme = 0.0
	
	if(t==2)then
		print*,' Proton spin only '
	    do i1 = 1,numorb(1)
	       j1  = orbqn(1,i1)%j
	       n1 = orbqn(1,i1)%nr
	       l1 = orbqn(1,i1)%l
	       do i2 = 1,numorb(1)
	           l2 = orbqn(1,i2)%l
	           if(l1 /= l2)cycle
           
	           n2 = orbqn(1,i2)%nr
	           if(n1/= n2)cycle
	           j2  = orbqn(1,i2)%j
                
	           p1bme(i1,i2) = sj2i(1,j1,2*l1,j2,1,2)*sqrt(1.5*(j1+1.)*(j2+1.)) & 
	               *(-1)**(( 2*l1 + 3 + j1)/2) !phase here changed by -1
	       end do
	    end do ! il		
		
	else
		print*,' Neutron spin only '
	    do i1 = 1,numorb(2)
	       j1  = orbqn(2,i1)%j
	       n1 = orbqn(2,i1)%nr
	       l1 = orbqn(2,i1)%l
	       do i2 = 1,numorb(2)
	           l2 = orbqn(2,i2)%l
	           if(l1 /= l2)cycle
           
	           n2 = orbqn(2,i2)%nr
	           if(n1/= n2)cycle
	           j2  = orbqn(2,i2)%j
          
	           n1bme(i1,i2) = sj2i(1,j1,2*l1,j2,1,2)*sqrt(1.5*(j1+1.)*(j2+1.)) & 
	               *(-1)**(( 2*l1 + 3 + j1)/2) !phase here changed by -1
	       end do
	    end do ! il			
		
	end if
	
	
end if
return
  end subroutine spinop

!---------------------------------------------
!
! note phase changed in Version 2, Nov 2013
!
  subroutine gamowteller

   use phonon1b
   use sporbit
   implicit none
   integer i1,i2
   integer n1,n2,l1,l2,j1,j2
   real sj2i   ! 6-j symbol from libra.f
   integer t
   real tfact
   real gA
   allocate ( t1bme( numorb(1), numorb(1) ) )

   t1bme(:,:) = 0.0

   tfact = sqrt(3.)     ! because reduced < 1/2 || tau || 1/2 > = sqrt(6)
                            ! but divide by sqrt(2) because of raising/lowering
							
   print*,' Enter axial coupling gA if desired (enter 1 if no scaling)'
   read*,gA
   do i1 = 1,numorb(1)
      j1  = orbqn(1,i1)%j
      n1 = orbqn(1,i1)%nr
      l1 = orbqn(1,i1)%l
      do i2 = 1,numorb(1)
          l2 = orbqn(1,i2)%l
          if(l1 /= l2)cycle
           
          n2 = orbqn(1,i2)%nr
          if(n1/= n2)cycle
          j2  = orbqn(1,i2)%j
          
          t1bme(i1,i2) = sj2i(1,j1,2*l1,j2,1,2)*sqrt(6.*(j1+1.)*(j2+1.)) & 
              *(-1)**(( 2*l1 + 3 + j1)/2)*tfact *gA  !phase here changed by -1
      end do
   end do ! il

  end subroutine gamowteller

!----------------------------------
!
! compute reduce matrix elements of ang mom \vec{J}
!
  subroutine jvec

   use phonon1b
   use sporbit
   implicit none
   integer i1,i2
   integer n1,n2,l1,l2,j1,j2
   real sj2i   ! 6-j symbol from libra.f
   integer t
   real tfact
   allocate ( t1bme( numorb(1), numorb(1) ) )

   t1bme(:,:) = 0.0

   tfact = sqrt(2.)     
   do i1 = 1,numorb(1)
      j1  = orbqn(1,i1)%j
      n1 = orbqn(1,i1)%nr
      l1 = orbqn(1,i1)%l
      do i2 = 1,numorb(1)
          l2 = orbqn(1,i2)%l
          if(l1 /= l2)cycle
           
          n2 = orbqn(1,i2)%nr
          if(n1/= n2)cycle
          j2  = orbqn(1,i2)%j
          if(j1 /= j2)cycle
          t1bme(i1,i2) = 0.5*sqrt(float(j1*(j1+1)*(j1+2)))*tfact
      end do
   end do ! il

  end subroutine jvec
!----------------------------------
!
! compute reduce matrix elements of Fermi operator raising/lower tau_+/-
!
  subroutine fermi

   use phonon1b
   use sporbit
   implicit none
   integer i1,i2
   integer n1,n2,l1,l2,j1,j2
   real sj2i   ! 6-j symbol from libra.f
   integer t
   real tfact
   allocate ( t1bme( numorb(1), numorb(1) ) )

   t1bme(:,:) = 0.0

   tfact = sqrt(3.)     
   do i1 = 1,numorb(1)
      j1  = orbqn(1,i1)%j
      n1 = orbqn(1,i1)%nr
      l1 = orbqn(1,i1)%l
      do i2 = 1,numorb(1)
          l2 = orbqn(1,i2)%l
          if(l1 /= l2)cycle
           
          n2 = orbqn(1,i2)%nr
          if(n1/= n2)cycle
          j2  = orbqn(1,i2)%j
          if(j1 /= j2)cycle
          t1bme(i1,i2) = sqrt(float(j1+1))*tfact
      end do
   end do ! il

  end subroutine fermi

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function r2me(n1,l1,n2,l2)
!
!  matrix elements of r2 between h.o. wfns
!   < n1, l1 | r^2 | n2, l2 >
!  INPUT:
!    n1,l1,n2,l2  : n = nodal quantum number (0, 1, 2,...)
!
      implicit none
      integer n1,l1,n2,l2

      r2me = 0.0

      if(n1 == n2 .and. l1 == l2)then
        r2me =(2*n1 +l1 + 1.5)
      endif

      if(n1 == n2 .and. l2 == l1 +2)then  ! ?? not confirmed
        r2me = sqrt((n1+l1 + 1.5)*(n1+l1+2.5))
      endif
      if(n1 == n2 .and. l2 == l1 -2)then
        r2me = sqrt((n1+l2 + 1.5)*(n1+l2+2.5))
      endif

      if(l1 == l2 .and. n2 == n1+1)then
        r2me = -sqrt(n2*(n2+l1+0.5))
      endif
      if(l1 == l2 .and. n2 == n1-1)then
        r2me = -sqrt(n1*(n1+l1+0.5))
      endif

      if(l2 == l1+2 .and. n2 == n1 -1)then  ! ?? not confirmed
        r2me = -sqrt(2*n1*(2*n1+2*l1+3.0))
      endif
      if(l2 == l1-2 .and. n2 == n1 +1)then
        r2me = -sqrt(2*n2*(2*n2+2*l2+3.0))
      endif
!      write(17,*)n1,n2,r2me
      return
      end function r2me

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      real function p2me(n1,l1,n2,l2)
!
!  matrix elements of r2 between h.o. wfns
!   < n1, l1 | p^2 | n2, l2 >
!  INPUT:
!    n1,l1,n2,l2  : n = nodal quantum number (0, 1, 2,...)
!
      implicit none
      integer n1,l1,n2,l2

      p2me = 0.0
      if(n1 == n2 .and. l1 == l2)then
        p2me =(2*n1 +l1 + 1.5)
      endif
      if(n1 == n2 .and. l2 == l1 +2)then ! ?? not confirmed
        p2me = -sqrt((n1+l2 + 1.5)*(n1+l2+2.5))
      endif
      if(n1 == n2 .and. l2 == l1 -2)then
        p2me = -sqrt((n1+l1 + 1.5)*(n1+l1+2.5))
      endif

      if(l1 == l2 .and. n2 == n1+1)then
        p2me = sqrt(n2*(n2+l1+0.5))
      endif
      if(l1 == l2 .and. n2 == n1-1)then
        p2me = sqrt(n1*(n1+l1+0.5))
      endif

      if(l2 == l1+2 .and. n2 == n1 -1)then ! ?? not confirmed
        p2me = sqrt(2*n1+(2*n1+2*l1+3.0))
      endif
      if(l2 == l1-2 .and. n2 == n1 +1)then
        p2me = sqrt(2*n2+(2*n2+2*l2+3.0))
      endif

      return
      end function p2me


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine dipole(t)

   use phonon1b
   use sporbit
   implicit none
   integer i1,i2
   integer n1,n2,l1,l2,j1,j2
   real sj2i   ! 6-j symbol from libra.f
   real tj2i  ,tji ! 3-j symbol
   integer t
   real tfact
   real r1me
   real peffcharge, neffcharge
   real tmp1b
   real bosc
   
   print*,' Enter value of oscillator length b in fm (enter 1 if not desired )'
   read*,bosc
   
   if(t==0 .or. t==1)then
	   
	   pnformal = .false.
   allocate ( t1bme( numorb(1), numorb(1) ) )
   t1bme(:,:) = 0.0

   select case (t)
   
   case (0)
       tfact = sqrt(2.)
	   case(1)
       tfact = sqrt(6.)
	   
   end select
   
   do i1 = 1,numorb(1)
      j1  = orbqn(1,i1)%j
      n1 = orbqn(1,i1)%nr
      l1 = orbqn(1,i1)%l
      do i2 = 1,numorb(1)
          l2 = orbqn(1,i2)%l
           
          n2 = orbqn(1,i2)%nr
          j2  = orbqn(1,i2)%j

          if( abs(l1-l2) /= 1)cycle
          if( abs(n1-n2) > 1)cycle
         r1me = 0.0
         if( n1 == n2 )then
             if(l1 == l2+1)then
                r1me = sqrt( (n2+l2+1.5) )
             else
                r1me = sqrt((n2+l2+0.5))
             endif
         else
              if(n1 == n2 -1 .and.l1 == l2+1)then
                r1me = -sqrt(float(n2))
              endif
              if(n1 == n2 +1 .and.l1 == l2-1)then
                r1me = -sqrt(float(n2+1))
              endif

         endif

          t1bme(i1,i2) = sj2i(2*l1,j1,1,j2,2*l2,2)*sqrt((j1+1.)*(j2+1.)) & 
              *(-1)**(1+( j2 +1)/2)*r1me*tfact*sqrt( float ( (2*l1+1)*(2*l2+1)*3)/(4*3.141592)) & 
       * tj2i(2*l1,2,2*l2,0,0,0) *bosc 
      end do
   end do ! il
   return
   
else
	pnformal = .true.
	print*,' Enter effective charges for proton, neutron '
	print*,' (Standard values are +N/Z for protons, -Z/A for neutron)'
	read*,peffcharge,neffcharge
    allocate ( p1bme( numorb(1), numorb(1) ),  n1bme( numorb(2), numorb(2) ) )
    p1bme(:,:) = 0.0	
	n1bme(:,:) =0.0

    do i1 = 1,numorb(1)
       j1  = orbqn(1,i1)%j
       n1 = orbqn(1,i1)%nr
       l1 = orbqn(1,i1)%l
       do i2 = 1,numorb(1)
           l2 = orbqn(1,i2)%l
           
           n2 = orbqn(1,i2)%nr
           j2  = orbqn(1,i2)%j

           if( abs(l1-l2) /= 1)cycle
           if( abs(n1-n2) > 1)cycle
          r1me = 0.0
          if( n1 == n2 )then
              if(l1 == l2+1)then
                 r1me = sqrt( (n2+l2+1.5) )
              else
                 r1me = sqrt((n2+l2+0.5))
              endif
          else
               if(n1 == n2 -1 .and.l1 == l2+1)then
                 r1me = -sqrt(float(n2))
               endif
               if(n1 == n2 +1 .and.l1 == l2-1)then
                 r1me = -sqrt(float(n2+1))
               endif

          endif

          tmp1b = sj2i(2*l1,j1,1,j2,2*l2,2)*sqrt((j1+1.)*(j2+1.)) & 
               *(-1)**(1+( j2 +1)/2)*r1me*sqrt( float ( (2*l1+1)*(2*l2+1))) & 
			   *tj2i(2*l1,2,2*l2,0,0,0)*bosc *sqrt(0.75/3.1415926)
		  p1bme(i1,i2)= peffcharge*tmp1b  
		  n1bme(i1,i2)= neffcharge*tmp1b  
			   
			   
       end do
    end do ! il
	
	
end if

  end subroutine dipole


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!==========================================================
! M1 operator is based on Brussaard and Glaudemans
! pg 218 Equation 10.64
! output was checked by comparing to Table 10.5 pg 219
! in Brussaard.
! the output is in nuclear magnetons \mu_N
! (by MKGK @ LLNL April 2014)

  subroutine m1(t)
    use phonon1b
    use sporbit
    implicit none
    integer i1,i2
    integer n1,n2,l1,l2,j1,j2
    real sj2i   ! 6-j symbol from libra.f
    integer t
    real tfact
    real, parameter :: M_PI = 3.141592654
    real gpl, gnl, gps, gns ! orbital and spin g-factors
    real lgs, sgs
    real prefac, orb_part, spin_part
	integer :: it

    gpl = 1.0
    gnl = 0.0
    gps = 5.5857
    gns = -3.8263

    print*,' g-factors used: '
	print*,' p/n g_l : ',gpl,gnl
	print*,' p/n g_s : ',gps,gns
	
    if(t/=2)then
	  allocate ( t1bme( numorb(1), numorb(1) ) )

!   orbital g combination (read as orbital g's)
    lgs = 0.5*(gpl + (-1)**t * gnl)
!   spin g combination (read as spin g's)
    sgs = 0.5*(gps + (-1)**t * gns)
    t1bme= 0.0
    ! this factor is the sqrt(2I+1) in Eq. 10.64   
    tfact = sqrt(2*t+1.)

    ! note j1 and j2 are really 2*j1 and 2*j2

    do i1 = 1,numorb(1)
       j1  = orbqn(1,i1)%j
       n1 = orbqn(1,i1)%nr
       l1 = orbqn(1,i1)%l
       do i2 = 1,numorb(1)
          l2 = orbqn(1,i2)%l
          if(l1 /= l2)cycle

          n2 = orbqn(1,i2)%nr
          if(n1/= n2)cycle
          j2  = orbqn(1,i2)%j

          prefac = (-1)**l1 * sqrt(3.0*(j1+1.)*(j2+1.)/(4.0*M_PI)) * tfact

          orb_part = lgs*(-1)**((j2+3)/2) * sqrt(2.*l1*(l1+1.)*(2*l1+1.)) &
               *sj2i(2*l1,2*l2,2,j2,j1,1)

          spin_part = sgs*(-1)**((j1+3)/2) * sqrt(3.) * sj2i(1,1,2,j2,j1,2*l1)

          t1bme(i1,i2) = prefac*(orb_part + spin_part)

       end do
    end do ! i1


    else
	   allocate(p1bme(numorb(1),numorb(1)))
	   allocate(n1bme(numorb(2),numorb(2)))
	   p1bme = 0.0
	   n1bme = 0.0
	   pnformal = .true.
	   
	   do it = 1,2
		   if(it==1)then
			   lgs = gpl
			   sgs = gps
			   
		   else
			   lgs = gnl
			   sgs = gns
			   
		   endif
		   
	       do i1 = 1,numorb(it)
	          j1  = orbqn(it,i1)%j
	          n1 = orbqn(it,i1)%nr
	          l1 = orbqn(it,i1)%l
	          do i2 = 1,numorb(it)
	             l2 = orbqn(it,i2)%l
	             if(l1 /= l2)cycle

	             n2 = orbqn(it,i2)%nr
	             if(n1/= n2)cycle
	             j2  = orbqn(it,i2)%j

	             prefac = (-1)**l1 * sqrt(3.0*(j1+1.)*(j2+1.)/(4.0*M_PI))   

	             orb_part = lgs*(-1)**((j2+3)/2) * sqrt(l1*(l1+1.)*(2*l1+1.)) &
	                  *sj2i(2*l1,2*l2,2,j2,j1,1)

	             spin_part = sgs*(-1)**((j1+3)/2) * sqrt(1.5) * sj2i(1,1,2,j2,j1,2*l1)

				 if(it==1)then

	                p1bme(i1,i2) = prefac*(orb_part + spin_part)
				else
	                n1bme(i1,i2) = prefac*(orb_part + spin_part)
					
				endif

	          end do ! i2
		   
		  end do  ! i1
	   end do  ! it
 
    end if

    return
  end subroutine m1

 
    
