!
! BAPPLYCENTROID.f90
!
! added in 7.9.7
!
! routine to apply centroids to a vector
!



subroutine applycentroid(startbundle,endbundle)
	
    use nodeinfo
    use localvectors
	use threebodycentroids
!    use system_parameters
    use precisions
!    use interaction
    use opbundles
    use fragments
    use basis
!    use lanczos_info
    use flagger
!    use bmpi_mod
!    use butil_mod
	use configurations
!    use contigpointervectors, only : vecin,vecout
    implicit none

    integer :: ibundle
    integer :: startbundle,endbundle
	integer(kind=8) :: state,nsd,psd,pstate
	integer(kind=4) :: ip,in
	real(kind=4) :: vmx
	do ibundle= startbundle,endbundle
		if(opbundle(ibundle)%optype/='CEN')cycle
		
!		print*,ibundle, opbundle(ibundle)%pxstart,opbundle(ibundle)%pxend
		do psd = opbundle(ibundle)%pxstart,opbundle(ibundle)%pxend
			ip = pcfindx(pmapSD2config(psd))
			pstate =pstart(psd)
			do nsd = opbundle(ibundle)%nxstart,opbundle(ibundle)%nxend
				in = nmapSD2config(nsd)
				
!				if(in+ip > numtotconfig)then
!					print*,psd,nsd,in,pmapSD2config(psd),in+ip
!					print*,ibundle,opbundle(ibundle)%pxstart,opbundle(ibundle)%pxend, opbundle(ibundle)%nxstart,opbundle(ibundle)%nxend
!					stop
!				end if
!				if(nstart(nsd)/=nsd)then
!					print*,' huh ',nsd,nstart(nsd)
!					stop
!				end if
				state = nstart(nsd)+pstate
				vmx = centroid(ip+in)
				vec2(state)=vec2(state)+ vec1(state)*vmx
			end do ! nsd
		end do ! psd
	end do
	
	return
end subroutine applycentroid