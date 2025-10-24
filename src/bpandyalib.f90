!=================================================================
!
! BPANDYALIB.f90
!
! routines for particle-hole conjugation
!
! initiated 11/2013 by CWJ @ SDSU
!
!=================================================================

subroutine pandamaster
   use system_parameters
   use coupledmatrixelements

   implicit none
   integer it
   ephshift = 0.0
   do it = 1,2
         call phshift(it)
   end do
   do it = 1,2
         call pandaspes(it)
   end do
   if((phconj(1) .and. .not.phconj(2)) .or. (phconj(2) .and. .not.phconj(1)))call pandacross

   return
end subroutine pandamaster

!============================================================

subroutine phshift(it)
   use system_parameters
   use coupledmatrixelements
   use sporbit

   implicit none
!.....INPUT........
   integer it   ! species
!.....INTERNAL......
   integer a,b    ! orbit labels
   integer :: J, Jmin, Jmax   ! angular momentum of two-body matrix elements
   type (vjs), pointer :: tbme(:)   ! points to XX two-body matrix elements
   integer,pointer :: xxmap(:)   
   real, pointer :: xspe(:)   ! points to single-particle energies
   
   integer pair1,iref,ipar,istart   ! elements for constructing the index

   integer indx

   if(.not.phconj(it))return

   if(it == 1)then
        tbme => ppme
        xxmap => ppcouplemap
        xspe  => pspe
   else
        tbme =>   nnme
        xxmap => nncouplemap
        xspe  => nspe

   endif

!........ SINGLE-PARTICLE-CONTRIBUTION........
   do a = 1,numorb(it)
	   if(orbqn(it,a)%w==99)cycle
         ephshift = ephshift+ xspe(a)*(orbqn(it,a)%j +1)
   end do
!.........TWO-BODY CONTRIBUTION FROM XX...............
   do a = 1,numorb(it)
	   if(orbqn(it,a)%w==99)cycle
	   
      do b = 1,a
   	   if(orbqn(it,b)%w==99)cycle
		  
          pair1 = a*(a-1)/2 + b
          if(xxmap(pair1)==-1)cycle
          pair1 = xxmap(pair1)
          ipar = XXcouples(it)%pairc(pair1)%par
          iref = XXcouples(it)%meref(ipar)
          istart = XXcouples(it)%mestart(ipar)
          indx = (pair1-iref)*(pair1-1-iref)/2+pair1-iref+istart   
          jmax = tbme(indx)%jmax
          jmin = tbme(indx)%jmin
          do j = jmin,jmax
             ephshift = ephshift + (2*j+1)*tbme(indx)%v(j)
          end do
      end do ! b
   end do  ! a
!..................... ALSO NEED CONTRIBUTION FROM PN......
   if(it ==2 .and. phconj(1) .and. phconj(2) )then
   do a = 1,numorb(1)
	   if(orbqn(1,a)%w==99)cycle
	   
      do b = 1,numorb(2)
   	   if(orbqn(2,b)%w==99)cycle
		  
          pair1 = numorb(2)*(a-1) + b
          pair1 = PNcouplemap(pair1)
          ipar = PNcouples%pairc(pair1)%par
          iref = PNcouples%meref(ipar)
          istart = PNcouples%mestart(ipar)
          indx = (pair1-iref)*(pair1-1-iref)/2+pair1-iref+istart   
          jmax = pnme(indx)%jmax
          jmin = pnme(indx)%jmin
          do j = jmin,jmax
             ephshift = ephshift + (2*j+1)*pnme(indx)%v(j)
          end do
      end do ! b
   end do  ! a

   end if
   return
end subroutine phshift
!============================================================
!
!  
!
subroutine pandaspes(it) 
   use system_parameters
   use coupledmatrixelements
   use sporbit
   use butil_mod

   implicit none
!.....INPUT........
   integer it   ! species
!.....INTERNAL......
   integer a,b    ! orbit labels
   integer :: J, Jmin, Jmax   ! angular momentum of two-body matrix elements
   integer :: nb   
   type (vjs), pointer :: tbme(:)   ! points to XX two-body matrix elements
   integer,pointer :: xxmap(:)   
   real, pointer :: xspe(:)   ! points to single-particle energies
   ! function in module:  real zeta   ! zeta(i,j) = sqrt(1+delta(i,j))
   real :: speshift
   integer pair1,iref,ipar,istart   ! elements for constructing the index

   integer indx

   if(it == 1)then
        tbme => ppme
        xxmap => PPcouplemap
        xspe  => pspe
   else
        tbme =>   nnme
        xxmap => NNcouplemap
        xspe  => nspe
   endif
   if(phconj(it))then
   do b = 1,numorb(it)
	   if(orbqn(it,b)%w==99)cycle
	   
     speshift = 0.0
     nb = orbqn(it,b)%j+1    ! max occupancy 
     do a = 1,numorb(it)
  	   if(orbqn(it,a)%w==99)cycle
		 
        if(a >= b)then
           pair1 =  a*(a-1)/2 + b
        else
           pair1 =  b*(b-1)/2 + a
        end if
        if(xxmap(pair1)==-1)cycle
        pair1 = xxmap(pair1)
        ipar = XXcouples(it)%pairc(pair1)%par
        iref = XXcouples(it)%meref(ipar)
        istart = XXcouples(it)%mestart(ipar)
        indx = (pair1-iref)*(pair1-1-iref)/2+pair1-iref+istart   
        jmax = tbme(indx)%jmax
        jmin = tbme(indx)%jmin
        do j = jmin,jmax
             speshift = speshift - zeta(a,b)**2* (2*j+1)*tbme(indx)%v(j)/float(nb)
        end do
     end do  ! a
!...........
     xspe(b) = -xspe(b) + speshift 
   end do   ! b
   end if
!..... SHIFT IN S.P.E.s DUE TO P-N INTERACTION.....
!   if(phconj(3-it) .and. np(1)*np(2) > 0)then
if(phconj(3-it) )then
   do b = 1,numorb(it)
	   if(orbqn(it,b)%w==99)cycle
	   
     speshift = 0.0
     nb = orbqn(it,b)%j+1    ! max occupancy 

        do a = 1,numorb(3-it)
	 	   if(orbqn(3-it,a)%w==99)cycle
			
            if(it==1)then
               pair1 = numorb(2)*(b-1)+a
            else
               pair1 = numorb(2)*(a-1)+b
            end if
            pair1 = PNcouplemap(pair1)
            ipar = PNcouples%pairc(pair1)%par
            iref = PNcouples%meref(ipar)
            istart = PNcouples%mestart(ipar)
            indx = (pair1-iref)*(pair1-1-iref)/2+pair1-iref+istart   
            jmax = pnme(indx)%jmax
            jmin = pnme(indx)%jmin
            do j = jmin,jmax
               speshift = speshift - (2*j+1)*pnme(indx)%v(j)/float(nb)
            end do

        end do

!...........
     if(phconj(it))then
         xspe(b) = xspe(b) + speshift 
     else
         xspe(b) = xspe(b) - speshift
     end if
   end do   ! b
   
   end if
   
!   print*,' spe shift ',it, xspe
   do a = 1,numorb(it)
!	   print*,it,orbqn(it,a)%w,ephshift,xspe(a)
	   if(orbqn(it,a)%w==99)cycle
	   
       xspe(a) = xspe(a) + ephshift/float(np(1)+np(2)) 
   end do 
   return
end subroutine pandaspes

!-----------------------------------------------------------
!
! if one species is p-h conjugated and the other is not
!
subroutine pandacross

   use system_parameters
   use coupledmatrixelements
   use sporbit
   use butil_mod

   implicit none
!.....INPUT........
   integer it   ! species
!....INTERNAL......
   type (vjs), pointer :: tbme(:)

   integer itc  ! conjugate species
   integer ax,by,cx,dy    ! orbit labels
   integer :: dystart
   integer :: nb   
!   real zeta   ! zeta(i,j) = sqrt(1+delta(i,j))
   real :: speshift
   integer pair1,pair2,iref,ipar,istart   ! elements for constructing the index

   integer indx1,indx2,tmp
   integer j,jmin1,jmax1,jmin2,jmax2,maxj
   integer :: k,kmin,kmax
   real, allocatable :: vtmp(:,:)
   real xja, xjb, xjc,xjd
   real :: sj  ! 6-j symbol, called from LIBRA.f
   integer :: aerr

   it = 1
   itc = 3 - it
   if(np(1)*np(2) ==0 )return   

!.......... FIND MAX J............
!           in order to reserve arrays
   maxj = 0
   do ax = 1,numorb(1)
	   if(orbqn(1,ax)%w==99)cycle
	   
      maxj = bmax(maxj,orbqn(1,ax)%j )
   end do
   do ax = 1,numorb(2)
	   if(orbqn(2,ax)%w==99)cycle
	   
      maxj = bmax(maxj,orbqn(2,ax)%j )
   end do

   allocate( vtmp(4, 0:maxj), stat=aerr )
   if(aerr /= 0) then
      call memerror("pandacross")
      stop 5
   end if
!........
   do ax = 1,numorb(it)
	   if(orbqn(it,ax)%w==99)cycle
	   
      xja = 0.5*orbqn(it,ax)%j
      do by = 1,numorb(itc)
   	   if(orbqn(itc,by)%w==99)cycle
		  
         xjb = 0.5*orbqn(itc,by)%j

         do cx = 1,ax
	  	   if(orbqn(it,cx)%w==99)cycle
			 
            xjc = 0.5*orbqn(it,cx)%j

            dystart = by
            do dy = dystart,numorb(itc)
		    	   if(orbqn(itc,dy)%w==99)cycle
				
                xjd = 0.5*orbqn(itc,dy)%j

!................. FIND INDICES AND CONJUGATES...................
                pair1 = numorb(2)*(ax-1) + by
                pair2 = numorb(2)*(cx-1) + dy
!     if( PNcouplemap(pair1) == -1) goto 1003
!     if( PNcouplemap(pair2) == -1) goto 1003
                pair1 = PNcouplemap(pair1)
                pair2 = PNcouplemap(pair2)

                ipar = PNcouples%pairc(pair1)%par
                if(ipar /= PNcouples%pairc(pair2)%par)cycle
                if(pair1 < pair2)then
                   tmp = pair1
                   pair1 = pair2
                   pair2 = tmp
                endif
                iref = PNcouples%meref(ipar)
                istart = PNcouples%mestart(ipar)
     
                indx1 = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart  
!.............. SWAP
                pair1 = numorb(2)*(ax-1) + dy
                pair2 = numorb(2)*(cx-1) + by
!     if( PNcouplemap(pair1) == -1) goto 1003
!     if( PNcouplemap(pair2) == -1) goto 1003
                pair1 = PNcouplemap(pair1)
                pair2 = PNcouplemap(pair2)

                ipar = PNcouples%pairc(pair1)%par

                if(pair1 < pair2)then
                   tmp = pair1
                   pair1 = pair2
                   pair2 = tmp
                endif
                iref = PNcouples%meref(ipar)
                istart = PNcouples%mestart(ipar)
                indx2 = (pair1-iref)*(pair1-1-iref)/2+pair2-iref+istart  
!                print*,ax,dy,cx,by !,pnme(indx2)%v(:)

!.................. FILL TEMPORARY ARRAYS
                vtmp = 0.0
                jmin1 = pnme(indx1)%jmin
                jmax1 = pnme(indx1)%jmax
                if(jmax1 > maxj)then   ! ERROR TRAP
                    print*,'(1) error in max js ',jmax1,maxj
                    stop
                end if
                do j = jmin1,jmax1
                   vtmp(1,j) = pnme(indx1)%v(j)
                end do ! j
                jmin2 = pnme(indx2)%jmin
                jmax2 = pnme(indx2)%jmax
                if(jmax2 > maxj)then   ! ERROR TRAP
                    print*,'(2) error in max js ',jmax2,maxj
                    stop
                end if
                do j = jmin2,jmax2
                   vtmp(2,j) = pnme(indx2)%v(j)
                end do ! j
!.................... TRANSFORM....................................       
                do j = jmin1,jmax1
                   vtmp(3,j) = 0.0
                   do k = jmin2,jmax2
                      vtmp(3,j) = vtmp(3,j) - (2*k+1) & 
             * sj(xja,xjd,float(k),xjc,xjb,float(j))*vtmp(2,k)
                   end do
                end do        
                do j = jmin2,jmax2
                   vtmp(4,j) = 0.0
                   do k = jmin1,jmax1
                      vtmp(4,j) = vtmp(4,j) - (2*k+1) & 
              *sj(xja,xjb,float(k),xjc,xjd,float(J) )*vtmp(1,k)
                   end do
                end do 
!................... NOW BACK FILL.................
                do j = jmin1,jmax1
                    pnme(indx1)%v(j) = vtmp(3,j)
                end do
                do j = jmin2,jmax2
                    pnme(indx2)%v(j) = vtmp(4,j)
                end do
            end do ! dy
         end do ! cx
      end do ! by
   end do  ! ax

   deallocate( vtmp)
   return
end subroutine pandacross
!=============================================================
!
! diagnostic routine to check on sum rule for PN part
! FOR TESTING ONLY
!
subroutine pnsumrule
   use system_parameters
   use coupledmatrixelements
   use sporbit

   implicit none
!.....INPUT........
   integer it   ! species
!....INTERNAL......
   integer itc  ! conjugate species
   integer ax,by   ! orbit labels
   
   integer pair1,iref,ipar,istart   ! elements for constructing the index

   integer indx1
   integer j,jmin1,jmax1
   real pnsum


   pnsum = 0.0
   do ax = 1,numorb(1)
      do by = 1,numorb(2)
!................. FIND INDICES AND CONJUGATES...................
                pair1 = numorb(2)*(ax-1) + by
!     if( PNcouplemap(pair1) == -1) goto 1003
!     if( PNcouplemap(pair2) == -1) goto 1003
                pair1 = PNcouplemap(pair1)

                ipar = PNcouples%pairc(pair1)%par

                iref = PNcouples%meref(ipar)
                istart = PNcouples%mestart(ipar)
     
                indx1 = (pair1-iref)*(pair1-1-iref)/2+pair1-iref+istart 
                jmin1 = pnme(indx1)%jmin
                jmax1 = pnme(indx1)%jmax
                
                do j = jmin1,jmax1
                   pnsum = pnsum + pnme(indx1)%v(j)*(2*j+1)
                end do ! j
      end do ! by
   end do  ! ax
   print*,' PN sum rule = ',pnsum

end subroutine pnsumrule

!=============================================================!
! a subroutine to write the transformed hamiltonian to file.
!
! normally not used
!
subroutine pandaprint
   use system_parameters
   use coupledmatrixelements
   use sporbit
   use nodeinfo
   use butil_mod
   implicit none
   integer i
   integer a,b,c,d,J,T
   integer pair1pn,pair2pn,pair1pp,pair2pp
   integer ipar,iref,istart,tmp,indxpp,indxpn
   integer jmin,jmax
   real xme
   real factor
   logical ident

   if(iproc/=0)return

   open(unit=41,file="pandya.int",status='unknown')
   write(41,*)'# energy shift = ',ephshift
   write(6,*)'# energy shift = ',ephshift

!......... WRITE OUT SINGLE PARTICLE ENERGIES with shift removed
    write(41,*)(pspe(i)-ephshift/float(np(1)+np(2)) ,i=1,numorb(1)),(nspe(i)-ephshift/float(np(1)+np(2)) ,i=1,numorb(2))
    write(6,*)(pspe(i)-ephshift/float(np(1)+np(2)) ,i=1,numorb(1))

!......... WRITE OUT T=0 MATRIX ELEMENTS
    do a = 1,numorb(1)
       do b = 1, a
          pair1pp = a*(a-1)/2 + b
          pair1pp = PPcouplemap(pair1pp)
          pair1pn = numorb(2)*(a-1) + b
          pair1pn = PNcouplemap(pair1pn)
          if(pair1pp == -1 .or. pair1pn == -1)then
            print*,' ACK! A PROBLEM '
          end if
          do c = 1,a
             do d = 1,c
                 if(c == a .and. d > b)cycle
                 factor = 1.0
                 if(a==b)factor = factor*sqrt(2.)
                 if(c==d)factor = factor*sqrt(2.)
                 pair2pp = c*(c-1)/2 + d
                 pair2pp = PPcouplemap(pair2pp)
                 pair2pn = numorb(2)*(c-1) + d
                 pair2pn = PNcouplemap(pair2pn)
                 if(pair2pp == -1 .or. pair2pn == -1)then
                    print*,' ACK! A PROBLEM '
                 end if
!................. COMPUTE PN INDEX
                
                 ipar = PNcouples%pairc(pair1pn)%par
                 if(ipar /= PNcouples%pairc(pair2pn)%par)cycle
!                      print*,' problem with parity, boss (pandya pn) ',a,b,c,d
 !                    stop
!                 endif
                 iref = PNcouples%meref(ipar)
                 istart = PNcouples%mestart(ipar)
                 if(pair1pn < pair2pn)then
                    indxpn = (pair2pn-iref)*(pair2pn-1-iref)/2+pair1pn-iref+istart   
                 else
                    indxpn = (pair1pn-iref)*(pair1pn-1-iref)/2+pair2pn-iref+istart   
                 end if
!................... COMPUTE PP INDEX...........

                 ipar = XXcouples(1)%pairc(pair1pp)%par
                 if(ipar /= XXcouples(1)%pairc(pair2pp)%par)then
                     cycle
!                     print*,' problem with parity, boss (pandya pp) ',a,b,c,d
!                     stop
                 endif
                 iref = XXcouples(1)%meref(ipar)
                 istart = XXcouples(1)%mestart(ipar)
                 if(pair1pp < pair2pp)then
                    indxpp = (pair2pp-iref)*(pair2pp-1-iref)/2+pair1pp-iref+istart  
                 else 
                    indxpp = (pair1pp-iref)*(pair1pp-1-iref)/2+pair2pp-iref+istart 
                 end if
!.................. WRITE OUT T=1 MATRIX ELEMENTS
                 ident = .false.
                 if(a==b .or. c==d)ident=.true.
                 do j = ppme(indxpp)%jmin,ppme(indxpp)%jmax
!... check if one can have T = 1; if identical particles need to have J even
                    if(ident .and. (-1)**j == -1)  cycle                  
                    write(41,'(6i4,f10.5)')a,b,c,d,j,1,ppme(indxpp)%v(j) !*4./factor
                 end do
				 
				 
!..........WRITE OUT T=0 MATRIX ELEMENTS
                 jmin = bmin(pnme(indxpn)%jmin,ppme(indxpp)%jmin)
                 jmax = bmax(pnme(indxpn)%jmax,ppme(indxpp)%jmax)

!                 print*,a,b,c,d,jmin,jmax
                 if(jmin>jmax)cycle
                 do j = jmin,jmax
                     xme = 0.0
                     if(ident .and. (-1)**j /=-1)cycle
                     if(j >= pnme(indxpn)%jmin .and. j <= pnme(indxpn)%jmax)then
                        xme = xme + pnme(indxpn)%v(j)*2/factor
                     end if
                     if(j >= ppme(indxpp)%jmin .and. j <= ppme(indxpp)%jmax)then
                        xme = xme - ppme(indxpp)%v(j)
                     end if
                     if(xme /= 0.0)write(41,'(6i4,f10.5)')a,b,c,d,j,0,xme!*2./factor
                 end do
             end do
          end do
       end do ! b
    end do ! a



   close(41)
   return
end subroutine pandaprint


