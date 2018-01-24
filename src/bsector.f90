! bsector.f90
! routines moved here to remove loops in
! module dependencies.
!
module bsector_mod
contains

!==================================================================================
!
!  write out information on sectors needed to model distribution of workload on parallel nodes
!
! CALLED BY main routine
!
subroutine writesectors4modeling

 use sectors
 use verbosity
 use io
 use nodeinfo
 use flagger
 use butil_mod

 implicit none
 integer it,itc

 integer is
 integer ic,ics
 integer ilast
 integer(8) nxsd, ncxsd
 integer(8) ntot 
 integer(8) maxsector, maxsubsector
 integer :: nbigsectors

 if(iproc /= 0)return
 it =1
 itc = 3 -it
 if(writeout)then
    ilast = index(outfile,' ')-1
    open(unit=modelinfo_file,file=outfile(1:ilast)//'.modelinfo',status='unknown')
 else
    open(unit=modelinfo_file,file='modelinfo.bigstick',status='unknown')
 endif
 ntot = 0
!...... LOOP OVER PROTON SECTORS
 write(modelinfo_file,101)it,nsectors(it)
 maxsector = 0
! maxsubsector = 0
 nbigsectors = 0
 do is = 1,nsectors(it)
    ncxsd = 0
    nxsd = xsd(it)%sector(is)%nxsd
    do ic = 1,xsd(it)%sector(is)%ncsectors
      ics = xsd(it)%sector(is)%csector(ic)
      ncxsd = ncxsd + xsd(itc)%sector(ics)%nxsd
!      maxsubsector = bmax(maxsubsector, nxsd*xsd(itc)%sector(ics)%nxsd)
    end do
    ntot = ntot + ncxsd*nxsd
    if( ncxsd*nxsd > maxfragmentsize)then
        nbigsectors = nbigsectors +1
!        print*,is, ncxsd*nxsd
    end if
    maxsector = bmax(maxsector,ncxsd*nxsd)
    write(modelinfo_file,101)is,nxsd*ncxsd,nxsd,ncxsd

101 format(i7,i12,6x,2i8)
 enddo  !is
 write(modelinfo_file,*)ntot,' total states '
 print*,' '
 print*,' sector information written to file '

 print*,' '
 print*,' Total of ',nsectors(1),' sectors '
 print*, ' with max of ',maxsector,' states in a sector '
 print*,' With ',nbigsectors,' greater than limit of ',maxfragmentsize,' (var: maxfragmentsize) '
! print*, ' (Max size of subsector = ',maxsubsector,')'

!   call count_create_fragments_noopt(.false.)
!   call count_create_fragments_noopt(.true.) 
return
end subroutine writesectors4modeling

!===============================================================
!
! routine to extract change in quantum numbers from the sectors
!
! INPUT:
!    it: species
!    is,fs  : initial, final sectors
! OUTPUT:
!    dM, dpar,dW
!
!  CALLED BY:
!    survey_geneologies
!
subroutine get_delta_q(it,is,fs,dm,dpar,dW)
   use sectors
   implicit none
!....... INPUT VARIABLES.....
   integer :: it ! species
   integer :: is,fs  ! initial, final sector
!....... OUTPUT VARIABLES...
   integer :: dM, dpar, dW
   dM = (xsd(it)%sector(fs)%jzX-xsd(it)%sector(is)%jzX)/2
   dpar = abs(xsd(it)%sector(fs)%parX-xsd(it)%sector(is)%parX)
   dW = xsd(it)%sector(fs)%wX-xsd(it)%sector(is)%wX

   return
end subroutine get_delta_q

end module bsector_mod
