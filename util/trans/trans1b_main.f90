!==============================================================
!
! TRANS1B
!
! makes (doubly) reduced matrix elements for common 1-body transition operators
!
! started 9/2010 by CWJ @ SDSU
!
!===================================================================
!
!  MAIN CALLING PROGRAM
!

  print*,' WELCOME TO THE ONE-BODY TROPICS '
  print*,' (reduced matrix elements for 1-body transition operators )'
  print*, ' Version 6 June 2018'


  call get_orbit_info
  call menu
  call phoneout

  end
!=================================================================

  subroutine phoneout
  use phonon1b
  use sporbit
  implicit none
  integer i1,i2,nme
  integer ilast
  character :: ychar
  

  character(45) filename

  print*,' ' 
  print*,' Enter name of output file (.opme) '
  read(5,'(a)')filename
  ilast= index(filename,' ')-1
  open(unit = 3,file=filename(1:ilast)//'.opme',status='unknown')
  
  if(pnformal .and. .not.xpn)then
	  print*,' Do you wish to be in explicit, one-column proton-neutron format? (y/n)'
	  read(5,'(a)')ychar
	  if(ychar=='y'.or.ychar=='Y')then
		  xpn= .true.
	  else
		  print*,' Writing in two-column proton-neutron format (pns)'
	  end if
  end if
  
  if(pnformal)then
	  if(xpn)then
		  write(3,'(a3)')'xpn'
	  else
	     write(3,'(a3)')'pns'
	 end if
  else
     write(3,'(a3)')'iso'
 end if
 
 if(xpn)then
     write(3,*)numorb(1)
     do i1 = 1,numorb(1)
         write(3,29)i1,i1+numorb(1),orbqn(1,i1)%nr, orbqn(1,i1)%l, float(orbqn(1,i1)%j)/2.
 29 format(2i4,2i5,f5.1)
     end do	 
	 
     nme = 0
     do i1 =1, numorb(1)*2
       do i2 = 1,numorb(1)*2
          if(p1bme(i1,i2)/= 0.0 .or. n1bme(i1,i2)/=0.0) nme = nme + 1
       end do
     end do
	 
 else
  write(3,*)numorb(1)
  do i1 = 1,numorb(1)
      write(3,99)i1,orbqn(1,i1)%nr, orbqn(1,i1)%l, float(orbqn(1,i1)%j)/2.
99 format(i4,2i5,f5.1)
  end do

  nme = 0
  do i1 =1, numorb(1)
    do i2 = 1,numorb(1)
       if(t1bme(i1,i2)/= 0.0) nme = nme + 1
    end do
  end do
  
end if
  write(3,*)jtrans,ttrans
!  write(3,*)nme
  
if(.not.pnformal)then
  do i1 =1, numorb(1)
    do i2 = 1,numorb(1)
       if(t1bme(i1,i2)/= 0.0)then
         write(3,100)i1,i2,t1bme(i1,i2)
       endif
100  format(2i5,2f14.7)
    end do
  end do
else
	
	if(xpn)then
		
	    do i1 =1, numorb(1)
	      do i2 = 1,numorb(1)
	         if(p1bme(i1,i2)/= 0.0 )then
	           write(3,100)i1,i2,p1bme(i1,i2)
	         endif
	      end do
	    end do
		
	    do i1 =1, numorb(2)
	      do i2 = 1,numorb(2)
	         if(n1bme(i1,i2)/= 0.0 )then
	           write(3,100)i1+numorb(1),i2+numorb(1),n1bme(i1,i2)
	         endif
	      end do
	    end do
		
		
	else
		
    do i1 =1, numorb(1)
      do i2 = 1,numorb(1)
         if(p1bme(i1,i2)/= 0.0 .or. n1bme(i1,i2)/=0.0)then
           write(3,100)i1,i2,p1bme(i1,i2),n1bme(i1,i2)
         endif
      end do
    end do
	
   end if
end if
  close(3)
  return
  end subroutine phoneout



