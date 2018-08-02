!
! code to convert isospin-format phonon to pn format
!

implicit none

character*80 :: filename
integer      :: ilast
character*80 :: dummline
character    :: choicechar

logical      :: success,itsthere
integer      :: numorb
integer, allocatable :: nrad(:),l(:)
real, allocatable    :: j(:),op(:,:)
integer       :: k,i,a,b
real          :: x
integer       :: jtrans,ttrans
real          :: cgiso  ! isospin CG
integer       :: icc   ! denotes whether charge-changing or charge conserving 
real   :: tj  ! compile with libra.f


print*,' '
print*,' This code converts a transition operator from isospin format to full pn format '
print*,' '
success = .false.

do while(.not.success)
	print*,' Enter name of input .opme file (must be in iso format)'
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	inquire(file=filename(1:ilast)//".opme",exist=success)
	if(.not.success)then
		print*,filename(1:ilast),'.opme does not exist '
		cycle
	end if
	
	open(unit=1,file=filename(1:ilast)//".opme",status='old')
	exit
end do ! while not success

success=.false.
do while(.not.success)
     read(1,'(a)')dummline
	 if(dummline(1:1)/='#' .and. dummline(1:1)/='!')then
		 success = .true.
	 else
		 cycle
	 end if
	 if(dummline(1:3)/='iso' .and. dummline(2:4)/='iso')then
		 print*,' Does not appear to be in isospin format'
		 print*,dummline
		 stop
	 end if
end do

read(1,*)numorb
print*,numorb,' orbits '
allocate(nrad(numorb),l(numorb),j(numorb))

do i = 1,numorb
	read(1,*)k,nrad(k),l(k),x
	j(k)=nint(2*x)
	
end do

read(1,*)jtrans,ttrans

print*,' J, T of transition = ',jtrans,ttrans

if(ttrans /=0 .and. ttrans /=1)then
	print*,' T must be 0 or 1 '
	stop
end if

allocate(op(numorb,numorb))

op = 0.0

do i = 1,numorb**2
	read(1,*,end=303)a,b,op(a,b)
	
	
end do

303 continue

close(1)

print*,' '


success = .true.

do while(success)
print*,' Enter name of output .opme file (will be  in xpn format)'
read(5,'(a)')filename
ilast = index(filename,' ')-1
inquire(file=filename(1:ilast)//".opme",exist=success)
if(success)then
	print*,' that file already exists, choose a different name '	
	cycle
end if

end do

open(unit=2,file=filename(1:ilast)//".opme",status='new')

write(2,'(a3)')'xpn'
write(2,*)numorb
do i = 1,numorb
	write(2,'(4i4,f6.1)')i,i+numorb,nrad(i),l(i),0.5*j(i)
end do
write(2,*)jtrans,ttrans


if(ttrans == 1)then
	print*,' Do you want to include (1) Only charge-conserving transitions '
	print*,' (2) Only charge-changing transitions  or '
	print*,' (3) Both charge-changing and charge-conserving? '
	read*,icc
end if


! FACTOR cgiso COMES FROM WIGNER-ECKART THEOREM
!  = (-1)^{1/2-mtf} ( 1/2  mtf,  1/2 -mti | 1 mtf-mti)/sqrt(3)


!print out pp

if(icc /=2)then
!cgiso = 1.0/sqrt(6.)
cgiso = tj(0.5,float(ttrans),0.5, -0.5,0.0, 0.5)
do a = 1,numorb
	do b = 1,numorb
		if(op(a,b)/=0.0)write(2,111)a,b,cgiso*op(a,b)
	end do
end do


!cgiso = 1.0/sqrt(6.)
cgiso = -tj(0.5,float(ttrans),0.5, 0.5,0.0, -0.5)

! print out nn

do a = 1,numorb
	do b = 1,numorb
		if(op(a,b)/=0.0)write(2,111)a+numorb,b+numorb,cgiso*op(a,b)
	end do
end do

end if


if(icc /= 1)then
!  print out pn
!cgiso = 1.0/sqrt(3.)
cgiso = tj(0.5,float(ttrans),0.5, -0.5,1.0, -0.5)

do a = 1,numorb
	do b = 1,numorb
		if(op(a,b)/=0.0)write(2,111)a,b+numorb,cgiso*op(a,b)
	end do
end do


! print out np
!cgiso = -1.0/sqrt(3.)
cgiso = -tj(0.5,float(ttrans),0.5, 0.5,-1.0, 0.5)
do a = 1,numorb
	do b = 1,numorb
		if(op(a,b)/=0.0)write(2,111)a+numorb,b,cgiso*op(a,b)
	end do
end do

end if
111 format(2i4,f10.5) 
	 
close(2)
end
