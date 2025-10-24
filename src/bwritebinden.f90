! routines to write out 1-body densities as binary files
! using routines written by O. Gorton for PANASh
! added 7.10.9, December 2023, by CWJ @ SDSU
!

module writebinden1b
	use sporbit
	use densities
	implicit none
	
    type momlist
      real (kind=8), allocatable :: rho(:,:)
    end type momlist
    type onebdenmat
      integer (kind=8) :: Kmin,Kmax ! min/max in change of angular momentum
      integer (kind=8) :: dpar     ! change in parity, for convenience
      type (momlist), allocatable :: klist(:)
    end type onebdenmat

    type level
      integer (kind=8) :: state
      real (kind=8) :: E
      real (kind=8) :: Ex
      integer (kind=8) :: jx2
      integer (kind=8) :: pi
      integer (kind=8) :: Tx2
    end type
	
contains
	
!==========================================================		
	subroutine binary_density_boss(ilevel,flevel)
		use io
		implicit none
		integer,intent(in) :: ilevel,flevel  ! 
        type (level),allocatable :: states(:)
		integer (kind=8) :: binfile
	    character(len=100) :: resfilename	
		
		resfilename = trim(outfile)//".dres.bin"
		allocate(states(flevel-ilevel+1))
		
		call fill_levels(ilevel,flevel,states)
		
	    open(newunit = binfile, file=resfilename, form='unformatted', access='stream')
		
		call write_binary_density(binfile,states)
		
		print*,' All done '
		
		return
		
	end subroutine binary_density_boss
!==========================================================	
	subroutine fill_levels(ilevel,flevel,states)
		use obs
		use system_parameters
	    use lanczos_util
		use sporbit
		
		implicit none
		integer,intent(in) :: ilevel,flevel  ! 
        type (level) :: states(:)
		integer :: nstates
		integer :: i,indx
	    logical :: evenAJ,evenAT    ! needed to force "correct" J, T
		integer :: myparity
		
	    nstates = size(states)
	    if( mod(np(1)+np(2),2) == 1 )then
	         evenAJ = .false.
	       evenAT = .false.
	    else
	       evenAJ = .true.
	       evenAT = .true.
	    end if
		
!		if(allsameparity)then
!			myparity = 1
			
!			if(np(1) > 0)myparity=myparity*(-1)**(np(1)* orbqn(1,1)%l)
!			if(np(2) > 0)myparity=myparity*(-1)**(np(2)* orbqn(2,1)%l)
!		end if
		indx = 0
		do i = ilevel,flevel
			indx = i -ilevel+1
			states(indx)%state = i
			states(indx)%E = energy(i)
			states(indx)%Ex = energy(i)-energy(1)
			states(indx)%jx2 = closest2J(evenAJ,xjlist(i))
			states(indx)%tx2 = closest2J(evenAT,xtlist(i))
!            if(allsameparity)then
!			    states(indx)%pi = myparity
!			else
			    states(indx)%pi = parity(i)
!			end if
		end do
		
		return
		
	end subroutine fill_levels
!==========================================================	
	
! original routine by OG (lightly edited by CWJ to conform to BIGSTICK internal arrays)================	
    subroutine write_binary_density(binfile, states)
      ! writes the proton and neutron one-body density matrices, prho and nrho,
      ! to a file in binary format:
      !   nstates, norbits, ncombo
      !   states(counter)%E, states(counter)%Ex, states(counter)%J, states(counter)%T <nstates lines>
      !   orbits(counter)%i, orbits(counter)%n, orbits(counter)%l, orbits(counter)%jx2 <norbits lines>
      !   rec_p
      !   rec_n
      ! where rec_p and rec_n are linear arrays containing the nominally-
      ! nonzero elements of the density matrices. The order of the elements 
      ! follows from the indices being ordered as:
      !   istate > fstate > a > b > k
      ! where k follows the triangle rule t(ji, jf, k) and t(ja, jb, k).
      ! binfile is the (integer (kind=8)) file identifier for the binary (unformatted,
      ! stream-access) file. 
      implicit none
      integer (kind=8), intent(in) :: binfile 
!      type (orbit), intent(in) :: orbits(:)
      type (level), intent(in) :: states(:)
!      type (onebdenmat), dimension(:,:), intent(in) :: prho, nrho

      real (kind=8), dimension(:), allocatable :: rec_p, rec_n
      integer (kind=8) :: istate, fstate, ji, jf
      integer (kind=8) :: aa, ja, bb, jb
      integer (kind=8) :: kmin, kmax, kk
      integer (kind=8) :: ncombo, counter
      integer (kind=8) :: norbits, nstates
      integer (kind=8) :: dpar_state, dpar_orbit

      norbits = max(numorb(1),numorb(2)) !size(orbits)
      nstates = size(states)

      print '("")'
      print '("Writing binary densities...")'
!      call system_clock(ta)
      print '("n orbits: ",i9)',norbits
      do aa = 1, norbits
 !         print '(5(i4))',orbits(aa)%i,orbits(aa)%n,orbits(aa)%l,orbits(aa)%jx2,&
  !            orbits(aa)%pi
        print '(5(i4))',aa,orbqn(1,aa)%nr,orbqn(1,aa)%l,orbqn(1,aa)%j,&
            orbqn(1,aa)%par
      end do
      print '("n states: ",i9)',nstates
      ncombo = 0
      do istate = 1, nstates
        do fstate = 1, nstates
           ji = states(istate)%jx2
           jf = states(fstate)%jx2
           dpar_state = states(istate)%pi * states(fstate)%pi
           do aa = 1, norbits
             do bb = 1, norbits
               dpar_orbit = (-1)**(orbqn(1,aa)%l  + orbqn(1,bb)%l)
               if (dpar_state /= dpar_orbit) cycle
               ja = orbqn(1,aa)%j
               jb = orbqn(1,bb)%j
               kmin = max(abs(ji-jf)/2, abs(ja-jb)/2)
               kmax = min(   (ji+jf)/2,    (ja+jb)/2)
               if (kmax < kmin) cycle
               ncombo = ncombo + (kmax-kmin+1)
             end do
           end do
         end do
      end do
      print '("n allowed: ",i9)',ncombo
      print '("Min. mem.: ",f10.3," MB")',real(2*ncombo)*sizeof(1.0)*1e-6 
      allocate(rec_p(ncombo), rec_n(ncombo))
      rec_p = 0.d0
      rec_n = 0.d0
      ncombo = 0
      do istate = 1, nstates
        do fstate = 1, nstates
           ji = states(istate)%jx2
           jf = states(fstate)%jx2
           dpar_state = states(istate)%pi * states(fstate)%pi
           do aa = 1, norbits
             do bb = 1, norbits
                 dpar_orbit = (-1)**(orbqn(1,aa)%l  + orbqn(1,bb)%l)
               !dpar_orbit = orbqn(1,aa)%par * orbqn(1,bb)%par
               if (dpar_state /= dpar_orbit) cycle
               ja = orbqn(1,aa)%j
               jb = orbqn(1,bb)%j
               kmin = max(abs(ji-jf)/2, abs(ja-jb)/2)
               kmax = min(   (ji+jf)/2,    (ja+jb)/2)
               if (kmax < kmin) cycle
               do kk = kmin, kmax
                 ncombo = ncombo + 1
                 rec_p(ncombo) = densitybag(istate, fstate)%denmat(kk,aa,bb,0) 
                 rec_n(ncombo) = densitybag(istate, fstate)%denmat(kk,aa,bb,1) 
!                 rec_p(ncombo) = prho(fstate, istate)%klist(kk)%rho(aa, bb)
!                 rec_n(ncombo) = nrho(fstate, istate)%klist(kk)%rho(aa, bb)
                 !print '(5i4,4f12.7)',istate,fstate,aa,bb,kk,rec_p(ncombo),rec_n(ncombo),&
                 !     prho(fstate, istate)%klist(kk)%rho(aa, bb), nrho(fstate, istate)%klist(kk)%rho(aa, bb)
               end do
             end do
           end do
         end do
      end do
!	  print*,' TEST ',nstates,norbits,ncombo
      write(binfile) nstates, norbits, ncombo
      do counter = 1, nstates
        write(binfile) states(counter)%E, states(counter)%Ex, states(counter)%Jx2, &
             states(counter)%Tx2, states(counter)%pi
      end do
      do counter = 1, norbits
        write(binfile) counter, int(orbqn(1,counter)%nr,kind=8), int(orbqn(1,counter)%l,kind=8), int(orbqn(1,counter)%j,kind=8)
      end do
      write(binfile)rec_p
      write(binfile)rec_n
      deallocate(rec_p)
      deallocate(rec_n)
  
      print '("Wrote ",I12," m.e..")',ncombo
  end subroutine write_binary_density
  
end module writebinden1b