!
! code to generate B(GT) from output of strength function results using Bigstick
! started 11/2020 by CWJ @ SDSU
!
!  FILE FORMAT
!  HEADER COMMETNT
!

implicit none

integer :: Zi,Zf,Ni,Nf
integer :: Zx,Nx
integer :: M2,Tz2, Ji2,Jf2,Ti2,Tf2
integer :: Tzi2,Tzf2,Tzx2

character(80) :: filenamein,filenameout
character(80) :: title
integer :: ilast
logical finished,foundit,first,firstinfile
integer :: iline,i,nlevels,ilevel,whichstrength,istrength
real :: Egsi,Egsf,Ei, ex,jx,tx
character :: ychar
real :: factor,ee,bgt
real :: cleb
real :: gtscale

print*,' '
print*,' Conversion of Gamow-Teller results from strength function'
print*,' '
print*,' Here is a list of what you will need: '
print*,' A .res file containing the ground state energy for the initial state space '
print*,' A .res file containing the ground state energy for the final state space '
print*,' A .res file containing the energy for the initial state  '
print*,' You may need to enter in by hand the initial J, T '
print*,' You can then enter in multiple final .res strength files '
print*,' You will need to enter in by hand the final J, T'
print*,' '
print*,' * * * * * * * * * '
print*,' '
print*,' Enter any scaling such as g_A, g_A x quench (this will get squared)'
read*,gtscale

print*,' Enter name of output file for B(GT) (leave off .bgt extension) '
read(5,'(a)')filenameout
ilast = index(filenameout,' ')-1
open(unit=1,file=filenameout(1:ilast)//".bgt",status='unknown')

finished = .false.
do while( .not. finished)
   print*,' Enter any header information you want to add  (multiple lines okay; enter "end" to stop)'
   read(5,'(a)')title
   if(title(1:3)=='end' .or. title(1:3)=='END')exit
	   
   write(1,"('# ',a80)")title
end do

write(1,"('# scaling (before squaring) = ',f10.5)")gtscale
gtscale = gtscale**2

print*,' '
print*,' Enter name of .res including the ground state of the initial state space '
print*,' (this is needed to establish relative energies )'
read(5,'(a)')filenamein
ilast = index(filenamein,' ')-1
open(unit=2,file=filenamein(1:ilast)//".res",status='old')
read(2,'(a)')title
print*,title
read(2,'(a)')title
print*,title
read(2,*)Zi,Ni
print*,'(valence) Z, N = ',Zi,Ni
write(1,*)Zi,Ni,' ! Initial space Z, N '
finished = .false.
do while(.not.finished)
	read(2,'(a)')title
	if(title(3:5)=="Sta")exit
	print*,title
end do
read(2,*)i,Egsi

if(i /= 1)then
	print*,' Did not find ground state '
	print*, i,Egsi
	stop
	
end if
print*,' Ground state energy initial = ',Egsi
write(1,*)Egsi, ' ! Initial space g.s. energy '
close(2)

!print*,' Enter 2 x initial J, T '
!print*,' (2 x T is probably ',abs(Tzi2),' )'
!read*,Ji2,Ti2
!print*,' '
!write(1,*)Ji2,Ti2

!print*,' Now you can enter in multiple files '

print*,' '
print*,' Enter name of .res file containing the ground state of the final state space '
print*,' (this is needed to establish relative energies )'
read(5,'(a)')filenamein
ilast = index(filenamein,' ')-1
open(unit=2,file=filenamein(1:ilast)//".res",status='old')
read(2,'(a)')title
print*,title
read(2,'(a)')title
print*,title
read(2,*)Zf,Nf
print*,'(valence) Z, N = ',Zf,Nf
write(1,*)Zf,Nf," ! Final space  Z, N"
finished = .false.
do while(.not.finished)
	read(2,'(a)')title
	if(title(3:5)=="Sta")exit
	print*,title
end do
read(2,*)i,Egsf

if(i /= 1)then
	print*,' Did not find ground state '
	print*, i,Egsf
	stop
	
end if
print*,' Ground state energy final = ',Egsf
write(1,*)Egsf,' ! Final space g.s. energy '
print*,' '
close(2)

!.............. NOW READ IN A STATE..........

2020 continue

print*,' Enter name of .res file with initial state energy '
print*,' (This may or may not be a file already read in )'
read(5,'(a)')filenamein
ilast = index(filenamein,' ')-1
open(unit=2,file=filenamein(1:ilast)//".res",status='old')

print*,filenamein(1:ilast),'.res. successfully opened! '

read(2,'(a)')title
print*,title
read(2,'(a)')title
print*,title
read(2,*)Zx,Nx
print*,'(valence) Z, N = ',Zx,Nx
finished = .false.
do while(.not.finished)
	read(2,'(a)')title
	if(title(3:5)=="Sta")exit
	print*,title
end do

do iline =1, 10000
	read(2,*,err=901,end=901)i,Ei,Ex,Jx,Tx
	nlevels = i
	if(i/=iline)then
		print*,' Some error in matching ',i,iline
		stop
	end if
	write(6,'(i4,2f10.5,2x,2f6.3)')i,Ei,Ex,Jx,Tx
end do
print*,' wait, I ran out '
stop
901 continue


print*,' Which state do you want? '
read*,ilevel
if(ilevel > nlevels) goto 901
rewind(2)
finished = .false.
do while(.not.finished)
	read(2,'(a)')title
	if(title(3:5)=="Sta")exit
	print*,title
end do
do iline = 1,ilevel
	read(2,*)i,Ei,Ex,Jx,Tx
end do
print*,' Initial energy is ', Ei,' which is an excitation energy of ',Ei-Egsi
write(1,*)Ei,' ! Energy of initial state '

print*,' '
print*,' Initial 2 x J, 2 x T appear to be : ',nint(2*Jx),nint(2*Tx)
print*,' Is this correct (y/n)?'
read(5,'(a)')ychar
if(ychar/='y' .and. ychar /='Y')then
	print*,' Enter correct initial 2 x J, 2 x T'
	read*,Ji2,Ti2
	print*,' You chose J of ',0.5*Ji2,' and T of ',0.5*Ti2
else
	Ji2 = Jx*2
	Ti2 = Tx*2
end if

write(1,*)Ji2,Ti2,' ! Initial state J x 2 , T x 2 '

!............... READ IN STRENGTH FUNCTION......

finished = .false.
firstinfile = .true.
do while(.not.finished)
print*,' '
print*,' Enter name of strength function file (leave off .res extension) '
if(.not.firstinfile)then
	print*,' (Enter "same" if you want to read more from this same .res file )'
end if

read(5,'(a)')filenamein
if(.not.firstinfile .and. filenamein(1:4)=='same')then
	rewind(2)
	write(1,'("# Sourced from same file")') 
	goto 333
	elseif(.not.firstinfile)then
		close(2)
end if
ilast = index(filenamein,' ')-1
open(unit=2,file=filenamein(1:ilast)//".res",status='old',err=331)

go to 332
331 continue
print*,' That file does not exist '
cycle


332 continue
print*,filenamein(1:ilast),'.res. successfully opened! '
print*,' '
write(1,'("# Sourced from ",a40,".res")')filenamein(1:ilast)

!write(1,*)"# Sourced from ",filenamein(1:ilast),'.res '
333 continue
read(2,'(a)')title
print*,title
read(2,*)Zx,Nx
print*,'(valence) Z, N = ',Zx,Nx
read(2,*)M2
print*,' 2 x Jz = ',M2
print*,' '
44 continue
print*,' Enter (by hand ) J x 2 , T x 2  of final state '
read*,Jf2,Tf2

!......... CHECK THIS IS CONSISTENT.....

if(mod(Jf2-M2,2)/=0 .or. mod(Tf2-Zx +Nx, 2)/=0)then
	print*,' Be sure to enter TWICE J, T '
	go to 44
end if
write(1,*)Jf2,Tf2,' ! Final state J x 2, T x 2 '
print*,' Enter which strength function you want '
print*,' (There may be multiple strength functions in a file;' 
print*,'  if only, enter "1" )'
read*,whichstrength


do istrength = 1,whichstrength
	foundit = .false.
	do while(.not.foundit)
		read(2,'(a)',end=111)title
		if(title(3:5)=='Ene')then
			foundit=.true.
			exit
		else
!			write(6,'(a)')title
		end if
	end do
111 continue	
	if(.not.foundit)then
		print*,' never found it'
		stop 123
	end if
	

end do
read(2,'(a)')title
print*,title


!............... CONMPUTE THE FACTOR HERE......
! The factor is given by eqn 5.34 in section 5.3.4 of the Bigstick manual

factor = 1.0

Tzx2 = Nx-Zx
Tzi2 = Ni-Zi
Tzf2 = Nf-Zf

!print*,Ti2,Tzi2,Tzf2-Tzi2,Tf2,Tzf2,cleb(Ti2,Tzi2,2,Tzf2-Tzi2 ,Tf2,Tzf2)


factor = cleb(Ti2,Tzi2,2,Tzf2-Tzi2 ,Tf2,Tzf2)
factor = factor / cleb(Ti2,Tzx2,2,0,Tf2,Tzx2)
factor = factor / cleb(Ji2,M2,2,0,Jf2,M2)
factor = factor**2
factor = factor*(Jf2+1.)/(Ji2+1.)

if(abs(factor) < 0.00001)then
	print*,' Some quantum numbers must have been entered incorrectly '
	stop
end if

!.................
first = .true.
do i = 1,10000
	read(2,*,err=321)ee,bgt
	if(bgt > 0.00001)then
		
		if(first)then
		    write(1,'(2f12.6,"  ! final energy, B(GT) ")')ee,bgt*factor
			first = .false.
		else

	       write(1,'(2f12.6)')ee,bgt*factor*gtscale
		   
	   end if
	end if
end do

321 continue


print*,' do you want to read in another strength for same initial state (y/n)? '
print*,' must be for the same initial energy '
read(5,'(a)')ychar

if(ychar=='y' .or. ychar=='Y')then
	finished = .false.
	firstinfile = .false.
else
	finished = .true.
end if


end do
close(2)

print*,' do you want to read in a different initial \ state (y/n)? '

read(5,'(a)')ychar
if(ychar=='y' .or. ychar=='Y')then
	write(1,'("# New initial state ")')
	goto 2020
end if

print*, ' All done !'


close(1)


end
