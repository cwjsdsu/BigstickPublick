!=============================================================================
! bpartitionslib1.f90
!
! introduced into BIGSTICK July 2019 in version 7.9.1
!
! data arrays and routines taken from TRACER version 6 (from April 2018)
!
!---------------------------------------------------------------------------------
!
! this module has information on partitioning particles
!
! each partition/configuration, for N particles, with Norb possible orbits
! is a binary word that is at least N x Nb long
! where 2^Nb-1 >= Norb
! for example, if Nb = 3, then can have up to 7 orbits
!
! Could probably be more compact, but this ought to be fine; either have a few particles and a lot of orbits, 
! or many particles but only a few orbits
!
! -- CWJ  June 2010
!


module configtools
	
contains
!	
!  ADDITIONAL FUNCTIONS
!=====================================================================
!
!  function to compute W of a configuration array
!
! CALLED BY: count_create_configs, sortconfigX
!
function wconfig(it,nsize,nparticles,array)
use sporbit
implicit none
integer wconfig
integer it
integer nsize
integer nparticles
integer array(nsize)
integer ip,iorb

wconfig = 0
if(nparticles == 0)return
do ip = 1,nparticles
   iorb = array(ip)
   wconfig = wconfig + orbqn(it,iorb)%w
end do

return
end function wconfig
!=====================================================================
!
!  function to compute parity of a configuration array
!
!  called by: sortconfigX
!
function parconfig(it,nsize,nparticles,array)
use sporbit
implicit none
integer parconfig
integer it
integer nsize
integer nparticles
integer array(nsize)
integer ip,iorb

parconfig = 1
if(nparticles == 0)return
do ip = 1,nparticles
   iorb = array(ip)
   parconfig = parconfig *orbqn(it,iorb)%par
end do

return
end function parconfig
!=====================================================================
!
!  function to compute occupation of a configuration array
!
function occconfig(it,nsize,nparticles,array,iorb)
use sporbit
implicit none
integer occconfig
integer it
integer nsize
integer nparticles
integer array(nsize)
integer ip,iorb

occconfig = 0
if(nparticles == 0)return
do ip = 1,nparticles
   if(iorb ==array(ip))occconfig = occconfig +1
end do

return
end function occconfig
!=====================================================================
!
! convert a bit rep of a config to an array
!
! INPUT:
!  nword = # of integers to represent the configuration
!  config(i = 1,nword) = integer array, with configuration encoded
!  nparticles = # of particles in this configuration
!  mask0 = basic mask needed to extract the orbit information
!  bitorb = # of bits needed to label orbits
!  bitword = # of bits in a word
!  nsize = size of array
!
! OUTPUT:
!  array(i =1,nsize) = array of orbital occupation; array(i) = label of orbit of ith particle
!
subroutine convertconfig2array(nword,config,nparticles,mask0,bitorb,bitword,nsize,array)

	implicit none
integer nword
integer nparticles
integer bitorb,bitword
integer config(nword)
integer nsize
integer array(nsize)
integer :: masktmp,tmp,mask0
integer iword, icount, maxcount
integer :: i

masktmp = mask0
iword = 1

maxcount = bitword/bitorb
icount = 0
array(:) = 0
do i = 1,nparticles
    icount = icount + 1
    if(icount > maxcount)then  ! reset
        iword = iword + 1
        masktmp = mask0
        icount = 1
    end if
    tmp = iand(config(iword),masktmp)
    array(i) = ishft(tmp,-(icount-1)*bitorb)   ! shift back

    masktmp = ishft(masktmp,bitorb)
end do 

return
end subroutine convertconfig2array
!=====================================================================
!
! convert an array to a bit rep of a config
!
! INPUT:
!  nword = # of integers to represent the configuration
!  nparticles = # of particles in this configuration
!  mask0 = basic mask needed to extract the orbit information
!  bitorb = # of bits needed to label orbits
!  bitword = # of bits in a word
!  nsize = size of array
!  array(i =1,nsize) = array of orbital occupation; array(i) = label of orbit of ith particle
!
! OUTPUT:
!  config(i = 1,nword) = integer array, with configuration encoded

subroutine convertarray2config(nword,config,nparticles,bitorb,bitword,nsize,array)

	implicit none
integer nword
integer nparticles
integer bitorb,bitword
integer config(nword)
integer nsize
integer array(nsize)
integer iword, icount, maxcount
integer :: i
integer :: tmp

iword = 1

maxcount = bitword/bitorb
icount = 0
config(:) = 0
do i = 1,nparticles
    icount = icount + 1
    if(icount > maxcount)then  ! reset
        iword = iword + 1
        icount = 1
    end if
    
    tmp  = ishft( array(i) ,(icount-1)*bitorb)   ! shift back
    config(iword) = config(iword) + tmp
end do 

return
end subroutine convertarray2config
!=====================================================================	
!
! convert an array to a bit rep of a config
!
! INPUT:
!  nword = # of integers to represent the configuration
!  nparticles = # of particles in this configuration
!  mask0 = basic mask needed to extract the orbit information
!  bitorb = # of bits needed to label orbits
!  bitword = # of bits in a word
!  nsize = size of array
!  array(i =1,nsize) = array of orbital occupation; array(i) = label of orbit of ith particle
!
! OUTPUT:
!  config(i = 1,nword) = integer array, with configuration encoded

subroutine convertorbarray2config(nword,config,nparticles,bitorb,bitword,nsize,orbarray)

	implicit none
integer nword
integer nparticles
integer bitorb,bitword
integer config(nword)
integer nsize
integer orbarray(nsize)
integer array(nparticles)
integer iword, icount, maxcount
integer :: i,j,ip
integer :: tmp

iword = 1
maxcount = bitword/bitorb

!...... First convert orbarray to array

ip = 0
do i = 1,nsize
	if(orbarray(i) > 0)then
		do j = 1,orbarray(i)
			ip = ip + 1
			if(ip > nparticles)then
				print*,' problem with converting orbarray to array ',ip,nparticles
				stop
			end if
			array(ip)=i
			
		end do
		
	end if
end do

icount = 0
config(:) = 0

do i = 1,nparticles
    icount = icount + 1
    if(icount > maxcount)then  ! reset
        iword = iword + 1
        icount = 1
    end if
    
    tmp  = ishft( array(i) ,(icount-1)*bitorb)   ! shift back
    config(iword) = config(iword) + tmp
end do 

return
end subroutine convertorbarray2config
!=====================================================================	

subroutine convertarray2orbarray(it, nsize, nparticles, array, orbarray)

use sporbit
implicit none
integer it
integer nsize
integer nparticles
integer array(nsize)
integer orbarray( numorb(it) )

integer ip,iorb

orbarray(:) = 0

do ip = 1,nparticles
    iorb = array(ip)
!    if(iorb ==0)then
!      print*,' huh ',array(:)
!    end if
    orbarray(iorb) = orbarray(iorb) + 1
enddo

return
end subroutine convertarray2orbarray
	
end module configtools

!=====================================================================	
!=====================================================================	
!=====================================================================	


module configurations
	use configtools
  implicit none
  integer :: bitword = 31  ! probably larger but need to check
  integer :: bitorb(2)     ! #  of bits needed to capture the orbits
  integer :: mask0(2)
  integer :: nwordpart(2)  ! # of 8-byte words needed to capture
  integer numconfig(2), numtotconfig
!............. PROTON, NEUTRON CONFIGURATIONS STORED AS INTEGER WORDS
! combined configurations are factorized into proton, neutron parts
  integer, allocatable, target :: configp(:,:), confign(:,:)

  integer, allocatable :: pcfindx(:)  ! proton neutron configuration indices
                                                  ! as with basis states
                                                  ! pcfindx(ip)+in=config index

!..........DERIVED TYPES FOR STORING INFORMATION ABOUT CONFIGURATIONS

  type conjugateinfo     ! information on conjugate configurations
     integer :: cwmin,cwmax,cnw,cwliststart
     integer :: nconjcf   ! # of conjcf that can be joined up
  end type conjugateinfo
  type pinfo
     integer:: wmin,wmax, nw
     integer, allocatable :: wlist(:)
     integer, allocatable :: wstart(:), wstop(:)
     type (conjugateinfo), allocatable :: wconj(:)
  end type pinfo

  type configinfo
     integer, allocatable :: nstates(:)  ! dimension
     real, allocatable :: e0(:)  ! centroid for a particular configuration
     integer, allocatable :: W(:)
     integer, allocatable :: parity(:)
     integer :: parstart(2),parstop(2)
     type (pinfo) :: parinfo(2)

  end type configinfo

  type (configinfo),target  :: Xconfiginfo(2)
! map proton, neutron SDs to proton, neutron configurations
  integer(kind=4), allocatable, target :: pmapSD2config(:),nmapSD2config(:)
  
  logical  :: configout

!.............. DATA FOR CENTROIDS............. 
!  added in 7.9.7 


  contains
!
!
! master subroutine to generate configurations
!
! generate configurations recursively
!
! configurations are factorized into proton and neutron configurations
!
!  IMPORTANT: Proton/Neutron configurations stored as integer words
!
! CALLED BY:
!    main routine
!
! SUBROUTINES CALLED:
!   configsetup
!   recurse_configs
!   sortconfigX
!   findconjugateconfigs
!   combinedconfigs
!   masxSD2xconfig
!
subroutine bossconfigurator(setupcentroids)

	use nodeinfo
   implicit none
   logical :: setupcentroids ! only to set up centroids
   integer it   ! species

   call configsetup(1)
   call configsetup(2)   
   call recurse_configs(1)
   call recurse_configs(2)   
   call sortconfigX(1)   
   call sortconfigX(2)
   
   call findconjugateconfigs
   call combinedconfigs(.false.)    ! combine all into total configurations
   if(iproc==0)print*,' There are ',numtotconfig,' configurations ',numconfig(:)
   
   if(setupcentroids)call combinedconfigs(.true.)   
   call mapxSD2xconfig(1)
   call mapxSD2xconfig(2)
   
   if(iproc==0)print*,' All mapped to configurations'

   return 

end subroutine bossconfigurator
!=====================================================================
!
!  sets up information needed to store configurations
!
!  IMPORTANT: Proton/Neutron configurations stored as integer words
!  for each particle of species it (=1 or 2), 
!  we need bitorb(it) bits to store what orbit it is in
!  So if we have 3 orbits, we need 2 bits
!  If we have between 4 and 7 orbits, we need 3 bits
!
!  Then each particle's orbit is encoded using that number of bits
!
!  For the configuration, which is an integer word, we split the bits
!  up. Suppose we have 7 orbits, so for each particle we have 3 bits.
!  If we have 4 particles, we need 4 x 3 bits = 12. 
!  The first 3 bits are for the 1st particle, the next 3 for the 
!  second particle, and so on. 
!
!  To decode this, we use a mask (mask0(it)) which is 2^bitorb-1
!  So for 7 orbits, bitorb = 3 and mask0 = 111. 
!
!  The reason for this encoding is that it only takes 1 or 2 integers
!  to encode a configuration; if one uses an array, it will often 
!  require many more integers and one can run out of space.
!
!  The subroutine CONVERTCONFIG2ARRAY takes an integer representation
!  and converts it to an array of dimension np(it) =# of particles
!  with each element telling  what orbit is occupied.
!
!  The subroutine CONVERTARRAY2CONFIG reverses this process.
!
!  These configurations are created recursively: 
!  starting from an empty occupation, take an array
!  of occupations and then add one more particle. The rule is, add
!  only to the highest occupied orbit, or higher; this creates 
!  unique configurations. 
!
!  The subroutine CONVERTARRAY2ORBARRAY takes an array of particles
!  and their occupation (dimension np(it) ) and the converts it to 
!  an array of orbits (dimension numorb(it) ) 
!  and the # of occupying particles. Actually, I could probably 
!  skip the middle man. 
!
!
!  If we have 12 particles, we need 36 bits. The parameter
!  bitword is set to the # of available bits per word, e.g. 31 or 63. 
!  So one may have multiple words to encode a configuration. 
!
! CALLED BY: bossconfigurator
!
subroutine configsetup(it)
!use configurations
use sporbit
use system_parameters

implicit none
integer it

!......find min # of bits needed to label orbits.....

bitorb(it) = 0
do while(  2**bitorb(it) <= numorb(it) )
   bitorb(it) = bitorb(it) +1
end do

mask0(it) = 2**bitorb(it)-1

! now figure # of words needed

nwordpart(it) = bitorb(it)*np(it)/bitword+1


return
end subroutine configsetup
!=====================================================================
!
!  recursively creates configurations of species it
!
!  CALLED BY: bossconfigurator
!
!  SUBROUTINES CALLED: count_create_configs
!
subroutine recurse_configs(it)

use system_parameters
!use configurations
use W_info
use nodeinfo
implicit none

integer it  ! species

integer, allocatable :: array0(:), array1(:)
integer, pointer :: config0(:,:), config1(:,:)
! note: I am actually trying to minimize the number of arrays allocated and deallocated
integer n
integer nold,nnew
integer wlim

allocate(array0(np(it)), array1(np(it)) )  ! form temporary arrays for creating configurations

nold = 1
allocate( config0(nwordpart(it),1) )
config0(nwordpart(it),1) = 0

if(np(it) == 0)then
   numconfig(it) = 1
   if(it == 1)then
      allocate(configp(nwordpart(it),1) )
      configp(:,1) = 0
   else
      allocate(confign(nwordpart(it),1) )
      confign(:,1) = 0
   endif
   return
endif
wlim = maxW(it)

do n = 1, np(it)  
	
   call count_create_configs(it,np(it),array0,array1,n,nold,config0,nnew,config1,wlim,.false.)
   
   if(n == np(it) ) then
      if(it == 1)then
          allocate( configp(nwordpart(it), nnew) )
          config1 => configp
      else
          allocate( confign(nwordpart(it), nnew) )
          config1 => confign
      endif
   else
      allocate( config1(nwordpart(it), nnew) )
   endif
   
   call count_create_configs(it,np(it),array0,array1,n,nold,config0,nnew,config1,wlim,.true.)
   
   nold = nnew
   if( n < np(it))then
      deallocate(config0)
      allocate(config0(nwordpart(it),nnew) )
      config0(:,:) = config1(:,:)
      deallocate(config1)
   endif
end do  
numconfig(it) = nnew
if(iproc==0)then
    print*,' '
    print*,it,' there are ',numconfig(it),' configurations '
    print*,' '
endif
deallocate(config0)
deallocate(array0,array1)
return
end subroutine recurse_configs

!=====================================================================
!
! creates a new round of configurations with n particles from n-1 particles
!
! loop over nold configurations with nparticles-1 particles
!     convert the bit rep of the config to an array
!     (calculate the current W of the config)
!     find last orbit occupied and # particles in it
!     loop over available orbits ( >= last orbit occupied)
!        check if can add and not exceed W limit
!        if okay, then add a particle to this orbit and count
!            (if createflag, then actually create the configuration)
!
!  CALLED BY: recurse_configs
! 
!  SUBROUTINES CALLED:
!    convertconfig2array
!    convertarray2config
!
!  FUNCTIONS CALLED:
!    occconfig
!
subroutine count_create_configs(it,nsize,array0,array1,nparticles,nold,oldconfig,nnew,newconfig,wlim,createflag)

!use configurations
use sporbit
implicit none

integer it  ! species
integer nsize  ! dimensioning size of array0 and array1
integer :: array0(nsize),  array1(nsize)
integer nparticles ! # of particles
integer nold     ! # of old partitions
integer nnew     ! # of new partitions
integer, pointer :: oldconfig(:,:), newconfig(:,:)
integer wlim  ! limits on w for this species
logical createflag

integer iconfig
integer lastorb
integer iorb
integer orbstart
integer lastocc
!integer wconfig  ! function to compute W
!integer occconfig  ! function to compute occupation
integer w0
integer k

nnew = 0
do iconfig = 1,nold  ! loop over old configurations and find last occupied configuration
   if(nparticles == 1) then
      array0(:) = 0
      orbstart = 1
      w0 = 0
   else
      call convertconfig2array( nwordpart(it), oldconfig(:,iconfig), nparticles-1, mask0(it), bitorb(it), & 
                                bitword, nsize, array0)
      lastorb = array0(nparticles -1)
      if(lastorb == 0)then
		  print*,' some problem with orbits '
         print*,iconfig, nparticles
         print*, oldconfig(:,iconfig)
         print*, array0
         stop
      endif
!......... get # of particles occupying last orbit...........
      lastocc = occconfig(it,nsize,nparticles-1,array0,lastorb)

!......... determine W for this configuration................
      w0 = wconfig(it,nsize,nparticles-1,array0)

      if( lastocc < orbqn(it,lastorb)%j+1)then
         orbstart = lastorb      ! still room
      else
         orbstart = lastorb + 1  ! no room, go to next orbit
      endif
   endif
   do iorb = orbstart,numorb(it)
            if( orbqn(it,iorb)%w + w0 > wlim)cycle   ! w too large, go on
            nnew = nnew + 1
            if(createflag)then
                array1(:) = array0(:)
                array1(nparticles) = iorb
                call convertarray2config( nwordpart(it), newconfig(:,nnew), nparticles,  & 
                      bitorb(it), bitword, nsize, array1)
            end if
   end do ! iorb

end do  

return
end subroutine count_create_configs
!=====================================================================
! 
! sorts configurations of species it
! first by parity
! then by W
! sets up where all configurations of same parity and W 
! start and stop
! 
!  CALLED BY: bossconfigurator
!  
!  SUBROUTINES CALLED:
!     convertconfig2array
!
!  FUNCTIONS CALLED:
!   wconfig, parconfig
!
subroutine sortconfigX(it)

use system_parameters
!use configurations
use sporbit
implicit none

integer it
integer, pointer :: configX(:,:)
integer, allocatable :: array(:)
integer i,j,k
integer numodd
integer :: tmp
integer, pointer :: qinfo(:)
integer ipar
integer wstart,w0,nw,nwsteps

if(it ==1)then
    configX => configp
else
    configX => confign
endif

allocate( Xconfiginfo(it)%w(numconfig(it) ) )
allocate( Xconfiginfo(it)%parity(numconfig(it) ) )

!...... simple case
if(np(it) == 0)then
    xconfiginfo(it)%parstart(1) = 1
    xconfiginfo(it)%parstop(1) = 1
    xconfiginfo(it)%parinfo(1)%nw = 1
    allocate( Xconfiginfo(it)%parinfo(1)%wlist(1) )
    xconfiginfo(it)%parinfo(1)%wlist(1) = 0
    allocate( Xconfiginfo(it)%parinfo(1)%wstart(1) )
    allocate( Xconfiginfo(it)%parinfo(1)%wstop(1) )
    xconfiginfo(it)%parinfo(1)%wstart(1) = 1
    xconfiginfo(it)%parinfo(1)%wstop(1) = 1

    xconfiginfo(it)%parstart(2) = 0
    xconfiginfo(it)%parstop(2) = -1
    xconfiginfo(it)%parinfo(2)%nw = 0
    allocate( Xconfiginfo(it)%parinfo(2)%wlist(1) )
    xconfiginfo(it)%parinfo(2)%wlist(1) = 0
    allocate( Xconfiginfo(it)%parinfo(2)%wstart(1) )
    allocate( Xconfiginfo(it)%parinfo(2)%wstop(1) )
    xconfiginfo(it)%parinfo(2)%wstart(1) = 0
    xconfiginfo(it)%parinfo(2)%wstop(1) = -1

    return
endif

!............. ASSIGN W, PARITY TO EACH CONFIGURATION................
allocate( array( np(it) ) )
do i = 1,numconfig(it)
    call convertconfig2array(nwordpart(it),configX(:,i),np(it),mask0(it),bitorb(it),bitword,np(it),array)

    Xconfiginfo(it)%w(i) = wconfig(it,np(it),np(it),array)
    Xconfiginfo(it)%parity(i) = parconfig(it,np(it),np(it),array)
!    print*,array, Xconfiginfo(it)%parity(i)
end do

deallocate(array)

!............... NOW DO MODIFIED BUBBLE SORT.................
!                BEGIN BY SORTING ON PARITY

if ( .not. allsameparity) then  
   numodd = 0
   do i = 1,numconfig(it)  ! count up odd parity states
      if( Xconfiginfo(it)%parity(i) == -1)numodd = numodd + 1
   end do
   if( numodd == 0 .or. numodd == numconfig(it))goto 12

!..... PUT ODD PARITIES UP TOP.....
   i = 1
   j = numconfig(it)
   qinfo => Xconfiginfo(it)%parity
   do while (i < j)
      if( qinfo(i) == 1)then
          i = i+ 1
          if(i == j)exit
      else
          if( qinfo(j) == -1)then
             j = j - 1
             if(j == i)exit
          else  ! swap
             do k = 1,nwordpart(it)
                 tmp = configX(k,i)
                 configX(k,i) = configX(k,j)
                 configX(k,j) = tmp
             end do
             tmp = Xconfiginfo(it)%w(i)
             Xconfiginfo(it)%w(i) = Xconfiginfo(it)%w(j)
             Xconfiginfo(it)%w(j) = tmp
             tmp = Xconfiginfo(it)%parity(i)
             Xconfiginfo(it)%parity(i) = Xconfiginfo(it)%parity(j)
             Xconfiginfo(it)%parity(j) = tmp
             i = i + 1
             j = j -1
          endif
      endif
   end do

end if

!........
12 continue

if( Xconfiginfo(it)%parity(1) == 1)then
    Xconfiginfo(it)%parstart(1) = 1
    Xconfiginfo(it)%parstop(1) = numconfig(it)  ! default
    do i = 1,numconfig(it)
        if(Xconfiginfo(it)%parity(i) == -1)then
             Xconfiginfo(it)%parstop(1) = i-1
             exit
        endif
    enddo
    if( Xconfiginfo(it)%parstop(1) == numconfig(it) )then
        Xconfiginfo(it)%parstart(2) =0 
        Xconfiginfo(it)%parstop(2) = -1
    else
        Xconfiginfo(it)%parstart(2) = Xconfiginfo(it)%parstop(1)+1
        Xconfiginfo(it)%parstop(2) = numconfig(it)
    endif
else
    Xconfiginfo(it)%parstart(1) = 0
    Xconfiginfo(it)%parstop(1) = -1
    Xconfiginfo(it)%parstart(2) = 1
    Xconfiginfo(it)%parstop(2) = numconfig(it)
endif

!................ SORT ON W ..............................
!................ SORT NEUTRONS IN REVERSE ORDER......

do ipar = 1,2
    if ( .not. allsameW )then
       if(Xconfiginfo(it)%parstart(ipar) == 0)cycle

         nwsteps = 0
         wstart = Xconfiginfo(it)%parstart(ipar)
         do while( (it==1.and. wstart <= Xconfiginfo(it)%parstop(ipar)) .or. & 
             (it==2.and. wstart <= Xconfiginfo(it)%parstop(ipar))           )
            w0 = 1000
            do i = wstart,  Xconfiginfo(it)%parstop(ipar)
                w0 = min(w0, Xconfiginfo(it)%w(i) )   ! find the next lowest value of W
            end do
            
            nw = 0
            do i = wstart,  Xconfiginfo(it)%parstop(ipar)
                if ( w0 == Xconfiginfo(it)%w(i) )nw = nw+1  ! count how many configs have w = w0
            end do           
            j = nw+wstart
            nwsteps = nwsteps +1
!            print*,w0
            if( j > Xconfiginfo(it)%parstop(ipar) )exit
            do i = wstart, wstart+nw-1
                if( Xconfiginfo(it)%w(i) /= w0)then
                    j = nw+wstart   ! for some reason I need to start over from the bottom
                    do while( Xconfiginfo(it)%w(j) /= w0)
                        j = j + 1
                 
                        if(j > Xconfiginfo(it)%parstop(ipar) ) then
							print*,' Error in setting up configurations '
                            print*,wstart,nw,i
                            print*,j
                            do k = Xconfiginfo(it)%parstart(ipar), Xconfiginfo(it)%parstop(ipar)
                               print*,k,Xconfiginfo(it)%w(k)
                            enddo
                            stop
                        end if
                    end do
!................... NOW SWAP.................
                    do k = 1,nwordpart(it)
                       tmp = configX(k,i)
                       configX(k,i) = configX(k,j)
                       configX(k,j) = tmp
                    end do
                    tmp = Xconfiginfo(it)%w(i)
                    Xconfiginfo(it)%w(i) = Xconfiginfo(it)%w(j)
                    Xconfiginfo(it)%w(j) = tmp                    
                endif
                j = j+1
            end do
            wstart = wstart+nw
         end do
         Xconfiginfo(it)%parinfo(ipar)%nw = nwsteps
         if(nwsteps == 0) cycle
         allocate ( Xconfiginfo(it)%parinfo(ipar)%wlist( nwsteps) )    
         allocate ( Xconfiginfo(it)%parinfo(ipar)%wstart( nwsteps) )    
         allocate ( Xconfiginfo(it)%parinfo(ipar)%wstop( nwsteps) )    

         Xconfiginfo(it)%parinfo(ipar)%wlist(1) = Xconfiginfo(it)%w( Xconfiginfo(it)%parstart(ipar) )
         Xconfiginfo(it)%parinfo(ipar)%wstart(1) =  Xconfiginfo(it)%parstart(ipar) 

         k = 1
         do i =  Xconfiginfo(it)%parstart(ipar),  Xconfiginfo(it)%parstop(ipar)
            if( Xconfiginfo(it)%w(i) > Xconfiginfo(it)%parinfo(ipar)%wlist(k) )then
                Xconfiginfo(it)%parinfo(ipar)%wstop(k) =  i-1                
                k = k + 1
                Xconfiginfo(it)%parinfo(ipar)%wstart(k) =  i

                if( k > nwsteps)then
                     print*,' w problem ',k,nwsteps
                     print*, Xconfiginfo(it)%w(i) , Xconfiginfo(it)%parinfo(ipar)%wlist(k-1) 
                     stop
                endif
                Xconfiginfo(it)%parinfo(ipar)%wlist(k) = Xconfiginfo(it)%w(i)
            endif
         enddo 
         Xconfiginfo(it)%parinfo(ipar)%wstop(nwsteps) =  Xconfiginfo(it)%parstop(ipar)               

!         print*, Xconfiginfo(it)%parinfo(ipar)%wlist(1:nwsteps)
   else
      nwsteps = 1 
      allocate ( Xconfiginfo(it)%parinfo(ipar)%wlist( nwsteps) )    
      allocate ( Xconfiginfo(it)%parinfo(ipar)%wstart( nwsteps) )    
      allocate ( Xconfiginfo(it)%parinfo(ipar)%wstop( nwsteps) )   
      Xconfiginfo(it)%parinfo(ipar)%nw = 1
      Xconfiginfo(it)%parinfo(ipar)%wlist( nwsteps) = 0
      Xconfiginfo(it)%parinfo(ipar)%wstart( nwsteps) = Xconfiginfo(it)%parstart(ipar)
      Xconfiginfo(it)%parinfo(ipar)%wstop( nwsteps) = Xconfiginfo(it)%parstop(ipar)

   end if
end do

return
end subroutine sortconfigX

!---------------------------------------------------------------------!
! combine proton and neutron configurations; 
! compute total centroid and total dimension
!
! If there is no cut on parity and/or W, 
! then all proton configurations combine with all neutron configs.
! If there are cuts, then one must match, for example,
! If parity is positive, then must match positive parity proton 
! configs with +parity neutron configs, and - parity with -parity.
! If a W-cut, check that combined W of proton and neutron
! do not exceed maximum. 
!
! also: set up array pcfindx  allows one to find the final configuration index
!       from the proton and neutron configurations
!  
! CALLED BY:
!   bossconfigurator
!
!
subroutine combinedconfigs(calc_centroid)

use io
use system_parameters
use sporbit
use W_info
use threebodycentroids
use nodeinfo
implicit none

logical calc_centroid

integer(kind=8) :: totaldim,tmpdim,alldim,localdim

integer :: icminp,icminn
integer :: parp, parn
integer :: iwp, iwn
integer :: wp, wn
integer ip,in
integer :: ii,jj,localni,localnj
integer ntot
integer, allocatable :: parray(:), narray(:), porbarray(:), norbarray(:)

icminp = 0
totaldim = 0
alldim = 0
ntot= 0

if(np(1) > 0) allocate( parray(np(1)))
if(np(2) > 0) allocate( narray(np(2)))

allocate(porbarray(numorb(1)) )
allocate(norbarray(numorb(2)) )
porbarray(:) = 0
norbarray(:) = 0

if(calc_centroid)then
	if(numtotconfig < 1)then
		print*,' Need to set up configs first '
		stop
	end if
	allocate(centroid(numtotconfig)) 
	centroid = 0.0
end if

do parp = 1,2
    if( iparity == 1)then
        parn = parp
    else
        parn = 3-parp
    endif
    do iwp = 1,Xconfiginfo(1)%parinfo(parp)%nw
       wp = Xconfiginfo(1)%parinfo(parp)%wlist(iwp)
!       localdim=0
!.......... LOOP OVER PROTON CONFIGURATIONS
!!$OMP parallel do private(ip,iwn,wn,in,tmpdim,e0tmp)  &
!!$OMP firstprivate(porbarray,norbarray,narray,parray) reduction(+:e0tot,localdim)
       do ip = Xconfiginfo(1)%parinfo(parp)%wstart(iwp), Xconfiginfo(1)%parinfo(parp)%wstop(iwp)
          if(calc_centroid)then  ! unravel occupations
              if(np(1) >0)then
              call convertconfig2array(nwordpart(1),configp(:,ip),np(1),mask0(1),bitorb(1),bitword,np(1),parray)
              call convertarray2orbarray(1, np(1), np(1), parray, porbarray)
              endif
          endif
          do iwn = 1,Xconfiginfo(2)%parinfo(parn)%nw
             wn = Xconfiginfo(2)%parinfo(parn)%wlist(iwn)
             if( wp+wn > maxWtot .or. wp+wn < minWtot) cycle
!............... LOOP OVER NEUTRON CONFIGURATIONS .................

             do in = Xconfiginfo(2)%parinfo(parn)%wstart(iwn), Xconfiginfo(2)%parinfo(parn)%wstop(iwn)
                ntot = ntot + 1
	            if(calc_centroid)then  ! unravel occupations
	                if(np(2) >0)then
	                call convertconfig2array(nwordpart(2),confign(:,in),np(2),mask0(2),bitorb(2),bitword,np(2),narray)
	                call convertarray2orbarray(2, np(2), np(2), narray, norbarray)					
	                endif
					
!....................... CHECK CONFIGURATION MAPPING..........
					if(ntot /= pcfindx(ip)+in)then
						print*,' mismatch in configuration mapping '
						print*,ntot,ip,in,pcfindx(ip)+in
						stop
					end if

!............ COMPUTE 3-BODY CONFIGURATION ......................
					
					if(np(1)>2)call e0config3b(1,porbarray,centroid(ntot))
					if(np(2)>2)call e0config3b(2,norbarray,centroid(ntot))
					if(np(1)*np(2) > 0)call e0pn3b(porbarray,norbarray,centroid(ntot))
	            endif
				
				

222 format(2i7,2x,i12,f12.6)
333 format(f12.6,i12)
334 format(25i3)

             end do
   
          end do  ! iwn
       end do ! ip
!!$OMP end parallel do
!       alldim = alldim + localdim
       end do ! iwp

end do  ! parp

    numtotconfig = ntot
	if(iproc==0)then

        print*,ntot,' configurations'

        if(writeout)then
            write(resultfile,*)ntot,' configurations'

        endif
	end if

deallocate(porbarray,norbarray)
return
end subroutine combinedconfigs
!----------------------------------------------------------------------
!
! find "conjugate" configs whose quantum numbers parity, W add up acceptably
!
! CALLED BY:
!   bossconfigurator
!

subroutine findconjugateconfigs
use io
use system_parameters
use sporbit
use W_info
use nodeinfo

implicit none
integer :: icminp,icminn
integer :: parp, parn
integer :: iwp, iwn
integer :: wp, wn
integer ip,in
integer nconjug,ntot
integer :: wnmax,wnmin,wcount
integer :: nref

allocate ( pcfindx( numconfig(1)) )

pcfindx = 0
ntot = 0
do parp = 1,2
    if( iparity == 1)then
        parn = parp
    else
        parn = 3-parp
    endif
    allocate ( Xconfiginfo(1)%parinfo(parp)%wconj( Xconfiginfo(1)%parinfo(parp)%nw ) )
    do iwp = 1,Xconfiginfo(1)%parinfo(parp)%nw
       wp = Xconfiginfo(1)%parinfo(parp)%wlist(iwp)
       nconjug = 0

       wnmin = 10000
       wnmax = 0  
       wcount = 0
       nref = -1
       do iwn = 1,Xconfiginfo(2)%parinfo(parn)%nw
             wn = Xconfiginfo(2)%parinfo(parn)%wlist(iwn)
             if( wp+wn > maxWtot .or. wp+wn < minWtot) cycle
             wnmin = min(wnmin,wn)
             wnmax = max(wnmax,wn)
             wcount = wcount+1
             if(nref < 0)nref = Xconfiginfo(2)%parinfo(parn)%wstart(iwn)  ! starting index for neutron states
!............... LOOP OVER NEUTRON CONFIGURATIONS .................
             nconjug = nconjug + Xconfiginfo(2)%parinfo(parn)%wstop(iwn) -Xconfiginfo(2)%parinfo(parn)%wstart(iwn)+1   
!................ SET UP NEUTRON INDICES.................
             
       end do  ! iwn
       Xconfiginfo(1)%parinfo(parp)%wconj(iwp)%nconjcf = nconjug 
!................ SET UP PROTON INDICES..................
       do ip = Xconfiginfo(1)%parinfo(parp)%wstart(iwp), Xconfiginfo(1)%parinfo(parp)%wstop(iwp)
           pcfindx(ip) = ntot-nref+1
           ntot = ntot + nconjug
       end do
!       ntot = ntot + nconjug*(Xconfiginfo(1)%parinfo(parp)%wstop(iwp) -Xconfiginfo(1)%parinfo(parp)%wstart(iwp)+1) 
       Xconfiginfo(1)%parinfo(parp)%wconj(iwp)%cwmin    = wnmin
       Xconfiginfo(1)%parinfo(parp)%wconj(iwp)%cwmax    = wnmax  
       Xconfiginfo(1)%parinfo(parp)%wconj(iwp)%cnw      = wcount  
 
    end do ! iwp

end do  ! parp

if(iproc==0)print*,' expect ',ntot,' configurations '
return
end subroutine findconjugateconfigs

!----------------------------- MAP SDs TO CONFIGURATIONS --------------
!
!  CALLS: convertorbarray2config
! 
!  CALLED BY: bossconfigurator
!
subroutine mapxSD2xconfig(it)
	
	use sectors
	use blocks
	use haiku_info
	use haikus
	use spstate
	use sporbit
	use system_parameters
	use basis
	implicit none
	
	integer, intent(in) :: it ! species: 1,2 = proton, neutron
	
	integer :: xs   ! sector index
	integer :: ix
	integer :: xblock,xrblock,xlblock, xsdstart
	integer :: nxradd,nxladd,xradd,xladd
	integer :: n,i
    integer,pointer  :: hsd(:)
	
	integer :: parx,Wx
	integer :: nw, iw
	integer :: configstart,configend
	integer :: iorb
    integer, allocatable :: xocctmp(:)
	integer, pointer :: configx(:,:)
	logical :: success
	integer, allocatable :: xorbarray(:)
	integer, allocatable :: tmpconfig(:)
	integer, pointer :: xmapSD2config(:)
	integer :: iconfig,iword
	
!	print*,' setting up partition maps ',nXsd(it)
	if(it==1)then
		configx => configp
		allocate(pmapSD2config(nXsd(it)))
		xmapSD2config=>pmapSD2config
		
	else
		configx => confign
		allocate(nmapSD2config(nXsd(it)))
		xmapSD2config=>nmapSD2config
	end if
	xmapSD2config = 0
	allocate(xorbarray(numorb(it)))
	allocate(xocctmp(nhsps(it)))
	allocate(tmpconfig(nwordpart(it)))
	
    do xs = 1,nsectors(it)  ! loop over  sectors
		parx = xsd(it)%sector(xs)%parx
		Wx   = xsd(it)%sector(xs)%Wx
		
!............ FIND CORRESPONDING START, END FOR CONFIGURATIONS.....
		
        nw = Xconfiginfo(it)%parinfo(parx)%nw
		success = .false.
		do iw = 1,nw
			if(Wx == Xconfiginfo(it)%parinfo(parX)%wlist(iw))then
				configstart = Xconfiginfo(it)%parinfo(parX)%wstart(iw)
				configend   = Xconfiginfo(it)%parinfo(parX)%wstop(iw)
				success = .true.
			end if
			
		end do
		if(.not.success)then
			print*,' whoops did not find W value when searching through configs '
			stop
		end if
!------------------- NOW LOOP OVER Slater dets ------------------------		
		
        do xblock = 1,xsd(it)%sector(xs)%nhblocks
           xrblock = xsd(it)%sector(xs)%rhblock(xblock)
           xlblock = xsd(it)%sector(xs)%lhblock(xblock)
           xsdstart= xsd(it)%sector(xs)%blockstart(xblock)-1


 !---------------- always make left address the outermost loop (a convention)
           nxradd = hblock(it)%list(xrblock)%nhsd
           nxladd = hblock(-it)%list(xlblock)%nhsd
           do xladd = 1,nxladd
              do xradd = 1,nxradd
                 ix = xsdstart + xradd+(xladd-1)*nxradd   ! this is the index of the Slater det
		
 !--------------------- now construct occupation -- left haiku
                 hsd => hblock(-it)%list(xlblock)%hsd(:,xladd)
                 call convert_haiku_array(-1,hsd,xocctmp)
				 xorbarray = 0
                 n = 0
                 do i = 1,nhsps(-it)
                   if(xocctmp(i) == 1)then
					   iorb = hspsqn(-it,i)%orb
					   xorbarray(iorb)=xorbarray(iorb)+1
!                      n = n+1
 !                     xocc(n) = i
                   endif
                 end do  ! i
 !--------------------- now construct occupation -- right haiku

                 hsd => hblock(it)%list(xrblock)%hsd(:,xradd)
                 call convert_haiku_array(1,hsd,xocctmp)
                
                 do i = 1,nhsps(it)
                   if(xocctmp(i) == 1)then
					   iorb = hspsqn(it,i)%orb
					   xorbarray(iorb)=xorbarray(iorb)+1
					   
!                      n = n+1
!                      xocc(n) = i+nhsps(-1)
                   endif
                 end do  ! i
				 tmpconfig=0
                 call convertorbarray2config( nwordpart(it), tmpconfig, np(it),  & 
                       bitorb(it), bitword, numorb(it), xorbarray)
					   
					   
				 
!-------------------- FINALLY, SEARCH AMONG CONFIGURATIONS, LIMITED BY QUANTUM NUMBER ------	
				 do iconfig = configstart,configend
					 success = .true.
					 do iword = 1,nwordpart(it)
						 if(tmpconfig(iword)/= configx(iword,iconfig))then
							 success =.false.
							 exit
						 end if
						 
					 end do
					 if(success)then
	!					 print*,it,' map IT ',ix,iconfig
						 xmapSD2config(ix)=iconfig						 
						 exit
					 end if
					 
					 
				 end do
				 if(.not.success)then
					 print*,' Did not find the configuration in the list, some problem ',ix
					 print*,xorbarray,' is array '
					 print*,tmpconfig,' is encoded config '
					 print*,configstart,configend
					 print*,configx
					 stop
				 end if			 
				 
			 end do  ! xradd
		 end do ! xladd
	 end do ! xblock
 end do ! xs
 
 deallocate(tmpconfig,xorbarray,xocctmp)
	
return	
	
end subroutine mapxSD2xconfig

!----------------------------------------------------------------------
!
!  routine to compute distribution of configurations in final wavefunctions
!  based upon routine Wcounter (in blanczoslib.f90)
!  added in 7.9.1 by CWJ SDSU  7/2019
!  
!  CALLS SUBROUTINES: writeconfigsall
! 
!  CALLED BY:  writeoutconfigocc
!  
subroutine configcounter
   use W_info
   use sporbit
   use basis
   use sectors
   use io
   use nodeinfo
   use lanczos_info
   use localvectors
   use flagger
   use fragments
   use mod_reorthog
   use wfn_mod
   use bmpi_mod
   use butil_mod
   use bvectorlib_mod
   use menu_choices
   implicit none

   integer(4) :: ierr
   integer is,isc,jsc
   integer(8) :: ip,in
   integer(kind=basis_prec) :: ibasis
   integer Wp, Wn,W
   integer :: nWvals,n
   real(8), allocatable :: configfrac(:)
   real(8) ::  ftmp,dv
   integer istate
   integer :: idummy
   real(4) :: xj,xt,ei,ef,xtt,xjj
   integer :: aerr
   
   integer(8) :: ipartp,ipartn, ipartall

   allocate(configfrac(numtotconfig), stat=aerr )
   if(aerr /= 0) call memerror("config counter 1")

!...................LOOP OVER WFNS..............................

  if(.not.storelanczosincore1 .and. .not.storelanczosincoreMPI)then
	  if(menu_char=='cx')then
	      call overlaptribution
		  
  		  call setup_localvectors
	      call wfn_rewind(oldwfnfile)   
	      call read_wfn_header(oldwfnfile,.false.)
	      call wfn_read_nkeep(oldwfnfile, nkeep)  ! dummy reading in nkeep	  
	  else	  
        call wfn_rewind(wfnfile)   
        call read_wfn_header(wfnfile,.false.)
        call wfn_read_nkeep(wfnfile, n)  ! dummy reading in nkeep
	end if
     ! read(wfnfile)n   ! dummy reading in nkeep
  endif
  do istate = 1,nkeep
     if((storelanczosincore1 .or. storelanczosincoreMPI).and. menu_char/='cx')then
        if(useNewReorthog) then
           call br_retrieve_hist(istate)
           ! It turns out that we need the vector loaded into both vec1 and vec2
           ! of course, they are different slices on each node
           call br_restore_vec1()
        else
           vec1 = 0.0 ! all ranks
           call read_lanczos_vector_a('n','i',istate,lvec_file) ! read in rank=0
           if(storelanczosincoreMPI) then
              ! call block_reduce(dimbasis,vec1)
              ! each mpi process reads one slice.  allreduce is overkill but works
#ifdef _MPI			  
              call BMPI_ALLREDUCE(vec1, size(vec1), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place
#endif
           end if
        end if
     else
        ! new interface, we say which vector to read and it checks
		
		if(menu_char=='cx')then
            call wfn_readeigenvec(oldwfnfile,frag1, fcomm1_index, vec1,istate,ei,xj,xtt)
			
		else
			
            call wfn_readeigenvec(wfnfile,frag1, fcomm1_index, vec1,istate,ei,xj,xtt)
		
	    end if
     end if
     configfrac=0.d0

     ! figure out if this node is the first node in ifragment.
     if(nodal(iproc)%ifirst) then
        do is = 1,nsectors(1)
           wp = xsd(1)%sector(is)%Wx
   
           do isc = 1,xsd(1)%sector(is)%ncsectors
              jsc = xsd(1)%sector(is)%csector(isc)
              do ip = xsd(1)%sector(is)%xsdstart,xsd(1)%sector(is)%xsdend
                 do in = xsd(2)%sector(jsc)%xsdstart,xsd(2)%sector(jsc)%xsdend
                    ibasis = pstart(ip)+nstart(in)
					
					ipartp = pmapSD2config(ip)
					ipartn = nmapSD2config(in)
!					print*,ipartall,ibasis,vec1(ibasis)
					ipartall = pcfindx(ipartp)+ipartn
					
!					if(ip ==4942 .and. in==4942)then
!						print*,' testing partitions ',ip,in,ibasis,ipartp,ipartn,pcfindx(ipartp)
!					endif
					
                    if(ibasis .ge. v1s .and. ibasis .le. v1e) then
                       dv = vec1(ibasis)
 !                      ftmp = ftmp + dv*dv
		               configfrac(ipartall) = configfrac(ipartall)+dv*dv
					   
                    end if
                 end do
              end do
           end do ! isc
        end do   ! is
     end if
	 
     if(useNewReorthog) then
        ! have to reduce.  Note that condition %ifirst above suppresses nodes that 
        ! are not "first" in their fragment.   This is not very efficient, but it works.
#ifdef _MPI
        call BMPI_ALLREDUCE(configfrac(:), SIZE(configfrac), MPI_SUM, MPI_COMM_WORLD, ierr) ! in place reduce
#endif
     end if
     if(iproc==0)then
		 write(configfile,"('# State ')")
		 write(configfile,'(i5)')istate  !,' energy ',ei
		 call writeconfigsall(configfrac)
     end if
  end do  ! istate

  return
end subroutine configcounter

!===================================================================
!
! called by: configcounter
!
subroutine writeconfigsall(configfrac)
	
	use system_parameters
	use io,only:configfile
	use W_info
	implicit none
	real(8) :: configfrac(numtotconfig)
	integer :: parp,parn
	integer :: iwp,iwn,wp,wn
	integer :: ip,in
	integer :: ipart
	
	write(configfile,'("#config   p config   n config      fraction")')
	do parp = 1,2
	    if( iparity == 1)then
	        parn = parp
	    else
	        parn = 3-parp
	    endif
	    do iwp = 1,Xconfiginfo(1)%parinfo(parp)%nw
	       wp = Xconfiginfo(1)%parinfo(parp)%wlist(iwp)
!.......... LOOP OVER PROTON CONFIGURATIONS

	       do ip = Xconfiginfo(1)%parinfo(parp)%wstart(iwp), Xconfiginfo(1)%parinfo(parp)%wstop(iwp)
	          do iwn = 1,Xconfiginfo(2)%parinfo(parn)%nw
	             wn = Xconfiginfo(2)%parinfo(parn)%wlist(iwn)
	             if( wp+wn > maxWtot .or. wp+wn < minWtot) cycle
!............... LOOP OVER NEUTRON CONFIGURATIONS .................

	             do in = Xconfiginfo(2)%parinfo(parn)%wstart(iwn), Xconfiginfo(2)%parinfo(parn)%wstop(iwn)
					 ipart = in+pcfindx(ip)
					 write(configfile,'(3i8,5x,f15.10)')ipart,ip,in,configfrac(ipart)
					 
				 end do ! in
			 end do ! iwn
		   end do ! ip
	   end do ! iwp
   end do ! parp
	
	
	
	return
end subroutine writeconfigsall

!======================================================================
!
!  CALLED BY: writeoutconfigocc
!  revised in 7.10.5 to make easier to read by machine
subroutine startoutputconfig(ifile)

use nodeinfo
use io
use system_parameters
use sporbit
use W_info

implicit none
integer ifile
integer iorb
integer it

if(iproc > 0)return
if(.not.writeout)return
write(ifile,'("#   Z    N ")')
write(ifile,'(2i5)')np(1),np(2)
!101 format(' Z = ',i4,' N = ',i4)

it = 1
write(ifile,"('# Proton orbits ')")
do iorb = 1,numorb(it)
   write(ifile,102)orbqn(it,iorb)%nr, orbqn(it,iorb)%l, orbqn(it,iorb)%j, orbqn(it,iorb)%w
102 format(4i5)
enddo
it = 2
write(ifile,"('# Neutron orbits ')")
do iorb = 1,numorb(it)
   write(ifile,102)orbqn(it,iorb)%nr, orbqn(it,iorb)%l, orbqn(it,iorb)%j, orbqn(it,iorb)%w
enddo
write(ifile,"('# Parity  ')")
write(ifile,'(i2)')iparity
if(allsameW)then
    write(ifile,'("# No W truncation ")')
else
    write(ifile,"('# Max W-excitation ')")
	write(ifile,*)maxWtot -minWtot
endif

return
end subroutine startoutputconfig

!======================================================
!
! CALLED BY:  writeoutconfigocc
!
! SUBROUTINES CALLED:
!     convertconfig2array
!     convertarray2orbarray

subroutine writeXconfigs(it)
use io
use system_parameters
use sporbit
use nodeinfo
implicit none


integer it

integer iconfig
integer, allocatable :: array(:),orbarray(:)
integer, pointer :: config0(:,:)
integer(kind=8) :: totaldim
integer i
integer :: maxorbused


allocate( array(np(it) ) )
allocate( orbarray( numorb(it) ) )

if(it ==1)then
   config0 => configp
   if(configout)write(configfile,"('# Proton configurations ')")
else
   config0 => confign
   if(configout)write(configfile,"('# Neutron configurations ')")
endif


do iconfig = 1,numconfig(it)
    call convertconfig2array(nwordpart(it),config0(:,iconfig),np(it),mask0(it),  & 
                             bitorb(it),bitword,np(it),array)
    call convertarray2orbarray(it, np(it), np(it), array, orbarray)
!    write(17,*)orbarray,Xconfiginfo(it)%nstates(iconfig), Xconfiginfo(it)%e0(iconfig)
    if(configout)then
		! this was changed in version 7.11.4, August 25, 2025, to address larger spaces
        write(configfile,*)iconfig
		write(configfile,111)(orbarray(i),i= 1,numorb(it))
111     format(20i4)
    endif
end do
!print*,' Total dimension = ',totaldim
return
end subroutine writeXconfigs


!===============================================
!
!  CALLS
!     startoutputconfig
!     writeXconfigs
!     configcounter
subroutine writeoutconfigocc
	use nodeinfo
	use io
	implicit none
	integer :: ilast
	if(iproc/=0)return

	ilast = index(outfile,' ')-1
	
	
		open(unit=configfile,file=outfile(1:ilast)//".cfo",status='unknown')
		write(6,*)" Configuration occupations written to :",outfile(1:ilast),".cfo" ,configfile
		write(logfile,*)" Configuration occupations written to : ",outfile(1:ilast),".cfo "
	
	call startoutputconfig(configfile)
	call writeXconfigs(1)
	call writeXconfigs(2)
	call configcounter
	
	if(iproc==0)print*,' Finished writing configurations'
	close(configfile)
	
	return
	
end subroutine writeoutconfigocc

!.......... ...........  ............ ..............  .............. .............  .............

end module configurations


