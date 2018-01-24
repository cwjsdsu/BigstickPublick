!===========================================================================
!  This code uses LAPACK ROUTINES
!
!  LAPACK copyright statements and license
!
!Copyright (c) 1992-2013 The University of Tennessee and The University
!                        of Tennessee Research Foundation.  All rights
!                        reserved.
!Copyright (c) 2000-2013 The University of California Berkeley. All
!                        rights reserved.
!Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
!                        reserved.

!Additional copyrights may follow

!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are
!met:
!
!- Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!
!- Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer listed
!  in this license in the documentation and/or other materials
!  provided with the distribution.
!
!- Neither the name of the copyright holders nor the names of its
!  contributors may be used to endorse or promote products derived from
!  this software without specific prior written permission.
!
!The copyright holders provide no reassurances that the source code
!provided does not infringe any patent, copyright, or any other
!intellectual property rights of third parties.  The copyright holders
!disclaim any liability to any recipient for claims brought against
!recipient by any third party for infringement of that parties
!intellectual property rights.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===== END OF LACK COPYRIGHT STATEMENT ========================


! for preconditioning of Lanczos pivot via self-consistent mean-field
! added Dec 2011 CWJ @ 2011
! NOTE: preconditioning removed in 7.5.5
!       restored in 7.7.2, as potentially part of LOBPCG
!
module shampoo
  use precisions
  implicit none
  logical, parameter :: ask_precondition = .true.
  logical :: precondition
  logical :: fasth     ! does only part of the Hamiltonian 
  integer :: nprecond    ! # of preconditioning iterations
  real, parameter :: betamin = 1.0e-4    ! minimum size for beta; stop if this is reached
  
  
!..... SOME TEMPORARY ARRAYS FOR LANCZOS..... in 7.7.6...
!  eventually these should be replaced by the structures and routines in breorthog.f90
!  However I cannot seem to make those work

  real(kind=lanc_prec), allocatable :: templancvecs(:,:)  


contains
!
      subroutine makepnHFpot

      use verbosity
      use spstate
      use sporbit
      use haiku_info
      use system_parameters
      use interaction
      use coupledmatrixelements
      use nodeinfo
      use ntuple_info
      use onebodypot
      use densities
      use butil_mod

      implicit none


      integer dpair,cpair
      integer Jtot,jmin,jmax
      integer T
      integer M
      integer par
      integer ik,il
      integer ii,ij

      integer i,j,k,l
      integer ith,jth,kth,lth
      integer isps,jsps,ksps,lsps
      integer itr,jtr,ktr,ltr

      integer nk,nl,ni,nj
      integer jk,jl,ji,jj
      integer mk,ml,mi,mj

      integer iref,istart
      integer itbme,itbmetr
      integer a,b,c,d
      integer ia,ib,ic,id
      integer pair1,pair2
      integer pair11,pair22
      integer ppar,pcref,pcstart
      integer phasekl,phaseij

      integer indx
      integer itmp
      real vtmp

      !! function in btbme_mod:  real zeta   ! zeta(i,j) = sqrt(1+delta(i,j))
      real cleb
      logical first
      integer nsection
      integer ipair,pairstart
      integer cindx,dindx  ! creation/destruction indices
      integer gi,gj,gk,gl ! group indices
      integer n0,n0j,n0k,nsumk
      integer kmin,kmax

      n0 = 0
      n0j = 0
      n0k = 0
      nsumk = 0
      if(Np(1)*Np(2) < 1)return

      vmatpn => hmatpn

	  ! initialize to prevent compiler warnings
      m = -9999
      par = 0
      nsection = 0
      pairstart = -1

      do dpair = 1,npairpn       ! loop over destruction pairs
!---------- FIND OUT IF THERE IS A CHANGE IN M OR PARITY
        if(par /= PN2%pair(dpair)%par .or. m /= PN2%pair(dpair)%m)then

!------------- ADD TO TOTAL THOSE FROM LAST TIME

!---------IF SO, THEN COUNT UP HOW MANY HAVE THE SAME M, PARITY
          m = PN2%pair(dpair)%m
          par = PN2%pair(dpair)%par
          pairstart = dpair


          do ipair = dpair,npairpn  ! loop upwards until M, PARITY change
                                      ! pairs have previously been sorted
             if(par == PN2%pair(ipair)%par .and. m == PN2%pair(ipair)%m)then
!------------ NSECTION IS HOW MANY PAIRS HAVE SAME M, PARITY
                nsection = ipair -dpair+1
             else
                exit
             endif
          enddo
          pairstart = dpair
          do ipair = dpair,1,-1  ! loop downwards until M, PARITY change
                                      ! pairs have previously been sorted
             if(par == PN2%pair(ipair)%par .and. m == PN2%pair(ipair)%m)then
!------------ NSECTION IS HOW MANY PAIRS HAVE SAME M, PARITY
                pairstart = ipair
             else
                exit
             endif
          enddo

        endif

!---------------- FROM NOW ON, UNCOUPLING -------------------------

!-------------- EXTRACT QUANTUM NUMBERS OF DESTRUCTION PAIR--------
        k = PN2%pair(dpair)%ia
        l = PN2%pair(dpair)%ib
        if( k < 0)then
          ksps = -k
          kth = -1
        else
          ksps = k
          kth  = 1
        endif
        ia = hspsqn(kth,ksps)%orb
        jk = hspsqn(kth,ksps)%j
        mk = hspsqn(kth,ksps)%m
        gk = hspsqn(kth,ksps)%group

        if( l < 0)then
          lsps = -l
          lth = -2
        else
          lsps = l
          lth  = 2
        endif
        ib = hspsqn(lth,lsps)%orb
        jl = hspsqn(lth,lsps)%j
        ml = hspsqn(lth,lsps)%m
        gl = hspsqn(lth,lsps)%group
        phasekl = 1

!-------------- FIND INDEX "PAIR1" used to find coupled tbme

        pair1 = numorb(2)*(ia-1)+ib
        pair1 = PNcouplemap(pair1)
        if(pair1 == -1)then
            print*,' uh problem pn map boss 1 '
            stop
        endif

!----------------- LOOP OVER CREATION PAIRS--------------
        cindx = 0
        do cpair = pairstart,pairstart+nsection-1
     
!-------------- EXTRACT QUANTUM NUMBERS OF CREATION PAIR--------

           i = PN2%pair(cpair)%ia
           j = PN2%pair(cpair)%ib

           if( i < 0)then
             isps = -i
             ith = -1
           else
             isps = i
             ith  = 1
           endif
           ic = hspsqn(ith,isps)%orb
           ji = hspsqn(ith,isps)%j
           mi = hspsqn(ith,isps)%m
           gi = hspsqn(ith,isps)%group
          if( j < 0)then
             jsps = -j
            jth = -2
          else
            jsps = j
            jth  = 2
          endif
          id = hspsqn(jth,jsps)%orb
          jj = hspsqn(jth,jsps)%j
          mj = hspsqn(jth,jsps)%m
          gj = hspsqn(jth,jsps)%group


!------------- PUT INTO "STANDARD" ORDER; MAY GET PHASE -------

          phaseij = 1

          pair2 = numorb(2)*(ic-1)+id
          pair2 = PNcouplemap(pair2)
          if(pair2==-1)then
           print*,' oh boy should not have gotten that '
           stop
          endif

!-------------- ALSO MUST HAVE STANDARD ORDERING OF PAIRS -----------

          if(pair1 < pair2)then
            pair11 = pair2
            pair22 = pair1
          else
           pair11 = pair1
           pair22 = pair2
          endif

          ppar = PNcouples%pairc(pair1)%par
          pcref = PNcouples%meref(ppar)
          pcstart = PNcouples%mestart(ppar)

          indx = (pair11-pcref)*(pair11-pcref-1)/2+pair22-pcref +pcstart
          if(indx <= 0)then
             print*,' problem indx ',indx,pair11,pair22
             print*,ia,ib,ic,id
             print*,cpair,dpair
             stop
          endif

!-------------FIND JMIN,JMAX ALLOWED ----------------

          jmax = MIN( jk+jl,ji+jj) /2
          jmin = MAX( abs(jk-jl),abs(ji-jj))/2
          jmin = MAX(jmin,abs(m))

!-------------NOW EXTRACT FROM TBMEs ---------------
!-------      NOTE FACTOR ZETA = SQRT(1+DELTA(A,B)) INCLUDED
          vtmp = 0.

          do jtot = jmin,jmax
                vtmp = vtmp+ pnme(indx)%v(jtot) & 
              * cleb(jk,mk,jl,ml,2*jtot,mk+ml)*cleb(ji,mi,jj,mj,2*jtot,mi+mj)
          enddo    ! jtot
!          print*, i,k,' ... ',j,l
          ppot_h(i,k) = ppot_h(i,k) + vtmp*n1bopme(j,l)
          npot_h(j,l) = npot_h(j,l) + vtmp*p1bopme(i,k)

        enddo  ! loop over cpair

      enddo  ! loop over dpair

      return
      end subroutine makepnHFpot

!==========================================================
!
!      subroutine uncoupleXXtbme
!
!  uncouples the T =1 matrix elements for pp, nn
!  includes the factor zeta = sqrt(1+delta(a,b)) etc.
!
!  IMPORTANT: How pp/nn matrix elements are encoded; 
!  needed for setting up the operator arrays for 
!  two-body "jumps"
!
!  INPUT:
!   it = species, 1 = proton, 2 = neutrons

!
      subroutine makeXXHFpot(it)

      use verbosity
      use spstate
      use haiku_info
      use system_parameters
      use interaction
      use coupledmatrixelements
      use onebodypot
      use densities
      use butil_mod
      use btbme_mod

      implicit none
      integer it      ! which species

      integer(kind=8) :: dpair,cpair
      integer Jtot,jmin,jmax
      integer T
      integer M
      integer par
      integer ik,il
      integer ii,ij

      integer i,j,k,l
      integer ith,jth,kth,lth
      integer isps,jsps,ksps,lsps

      integer nk,nl,ni,nj
      integer jk,jl,ji,jj
      integer mk,ml,mi,mj

      integer(kind=8) :: iref,istart
      integer(kind=8) :: itbme

      integer ia,ib,ic,id
      integer(kind=8) :: pair1,pair2
      integer(kind=8) :: pair11,pair22
      integer phasekl,phaseij

      integer(kind=8) :: indx
      integer itmp
      real vtmp

      ! function in module btbme_mod:   real zeta   ! zeta(i,j) = sqrt(1+delta(i,j))
      real cleb
      real hfactor  ! FOR HERMICIITY
      real,pointer :: v(:)

      type (vjs), pointer :: tbme(:)
      integer pcref, pcstart, ppar
      integer,pointer :: xxmap(:)

      if(Np(it) < 2)return

      if(it == 1)then
        tbme => ppme
        xxmap => ppcouplemap

      else
        tbme =>   nnme
        xxmap => nncouplemap

      endif

      T = 1        ! for pp, nn matrix elements

      do dpair = 1,npairXX(it)        ! loop over destruction pair

        m = XX2(it)%pair(dpair)%m     ! find Jz, parity of pair
        par = XX2(it)%pair(dpair)%par

        iref = XX2(it)%meref(m,par)      ! iref, istart used to encode index for 
        istart = XX2(it)%mestart(m,par)  ! the uncoupled matrix element

        if(iref > dpair)then
          print*,' problem with XX iref (d) ',it
          print*,dpair,iref
          print*,m,par
          print*,xx2(it)%mestart(:,1)
          stop
        endif

!---------------------- GET Q#s FOR PARTICLES IN DESTRUCTION PAIR --------
        k = XX2(it)%pair(dpair)%ia
        l = XX2(it)%pair(dpair)%ib
        if( k < 0)then
          ksps = -k
          kth = -it
        else
          ksps = k
          kth  = it
        endif
        ia = hspsqn(kth,ksps)%orb
        jk = hspsqn(kth,ksps)%j
        mk = hspsqn(kth,ksps)%m

        if( l < 0)then
          lsps = -l
          lth = -it
        else
          lsps = l
          lth  = it
        endif
        ib = hspsqn(lth,lsps)%orb
        jl = hspsqn(lth,lsps)%j
        ml = hspsqn(lth,lsps)%m

        phasekl = 1
!-------------- MUST HAVE ia > ib, else swap and pick up a phase ---

        if(ia < ib)then
             itmp = ia
             ia = ib
             ib = itmp

             itmp = jk
             jk = jl
             jl = itmp

             itmp = mk
             mk   = ml
             ml   = itmp

             phasekl = -1  ! the rest of the phase comes from the Clebsch Gordan

        endif
 

        do cpair = 1,dpair !,npairXX(it)        ! loop over creation pair
          if(cpair == dpair)then
              hfactor = 1.0  ! 0.5  I don't think I am using Hermiticity here
          else
              hfactor = 1.0
          endif
          if(XX2(it)%pair(cpair)%m /=m)cycle       ! enforce quantum numbers
          if(XX2(it)%pair(cpair)%par /= par)cycle

          if(iref > cpair)then
            print*,' problem with XX iref (c) '
            print*,cpair,iref
            stop
          endif

!--------------- GET Q#s FOR CREATION PAIRS ---------------

          i = XX2(it)%pair(cpair)%ia
          j = XX2(it)%pair(cpair)%ib

          if( i < 0)then
            isps = -i
            ith = -it
          else
            isps = i
            ith  = it
          endif
          ic = hspsqn(ith,isps)%orb
          ji = hspsqn(ith,isps)%j
          mi = hspsqn(ith,isps)%m

          if( j < 0)then
            jsps = -j
            jth = -it
          else
            jsps = j
            jth  = it
          endif
          id = hspsqn(jth,jsps)%orb
          jj = hspsqn(jth,jsps)%j
          mj = hspsqn(jth,jsps)%m
!------------- PUT INTO "STANDARD" ORDER; MAY GET PHASE -------

          phaseij = 1

          if(ic < id) then

             itmp = ic
             ic = id
             id = itmp

             itmp = ji
             ji = jj
             jj = itmp
             itmp = mi
             mi   = mj
             mj   = itmp

             phaseij = -1  ! the rest of the phase comes from the Clebsch Gordan

          endif
          pair1 = ia*(ia-1)/2 + ib  ! used to find index of coupled TBME

          pair2 = ic*(ic-1)/2 + id



!-------------- ALSO MUST HAVE STANDARD ORDERING OF PAIRS -----------

          pair1 = xxmap(pair1)
          pair2 = xxmap(pair2)
          if(pair1==-1 .or. pair2 == -1)then
             print*,' hmm strange ...',it
             print*,pair1,pair2 
             print*,ia,ib,ic,id
             stop
          endif
          ppar = XXcouples(it)%pairc(pair1)%par
          if(ppar /= XXcouples(it)%pairc(pair2)%par)then
              print*,' sigh parity problem '
              print*,ia,ib,ic,id,ppar,par
              print*,pair1,pair2
              stop
          endif
          if(pair1 < pair2)then
            pair11 = pair2
            pair22 = pair1
          else
           pair11 = pair1
           pair22 = pair2
          endif
          pcref = XXcouples(it)%meref(ppar)
          pcstart = XXcouples(it)%mestart(ppar)

          indx = (pair11-pcref)*(pair11-1-pcref)/2+pair22-pcref+pcstart
          if(indx <= 0)then
             print*,' problem indx ',indx
             stop
          endif

!-------------FIND JMIN,JMAX ALLOWED ----------------

          jmax = MIN( jk+jl,ji+jj) /2
          jmin = MAX( abs(jk-jl),abs(ji-jj))/2
          jmin = MAX(jmin,abs(m))

!-------------NOW EXTRACT FROM TBMEs ---------------
!-------      NOTE FACTOR ZETA = SQRT(1+DELTA(A,B)) INCLUDED
          vtmp = 0.

          do jtot = jmin,jmax

             vtmp = vtmp+ tbme(indx)%v(jtot) * zeta(ia,ib)*zeta(ic,id)   * &
                    cleb(jk,mk,jl,ml,2*jtot,2*m)*cleb(ji,mi,jj,mj,2*jtot,2*m)
          enddo           

          if(it == 1)then
             ppot_h(i,k) = ppot_h(i,k) + vtmp*phasekl*phaseij*p1bopme(j,l)
             ppot_h(j,l) = ppot_h(j,l) + vtmp*phasekl*phaseij*p1bopme(i,k)

             if(cpair < dpair)then
                ppot_h(k,i) = ppot_h(k,i) + vtmp*phasekl*phaseij*p1bopme(l,j)
                ppot_h(l,j) = ppot_h(l,j) + vtmp*phasekl*phaseij*p1bopme(k,i)

             end if
          else
             npot_h(i,k) = npot_h(i,k) + vtmp*phasekl*phaseij*n1bopme(j,l)
             npot_h(j,l) = npot_h(j,l) + vtmp*phasekl*phaseij*n1bopme(i,k)

             if(cpair < dpair)then
                npot_h(k,i) = npot_h(k,i) + vtmp*phasekl*phaseij*n1bopme(l,j)
                npot_h(l,j) = npot_h(l,j) + vtmp*phasekl*phaseij*n1bopme(k,i)

             end if
          end if
        enddo  ! loop over cpair

      enddo  ! loop over dpair
!      if(it == 1)then
!        do i = -3,3
!          write(6,'(7f8.3)')(ppot_h(i,j),j=-3,3)
!        enddo
!      end if
      return
      end subroutine makeXXHFpot
!===================================================
subroutine makespHFpot(it)

!    use interaction
    use coupledmatrixelements
	
    use onebodypot
    use spstate
    use haiku_info
    use system_parameters

    implicit none
    integer it
    integer isps
    integer i
    integer ith
    integer iorb

    if(np(it) < 0)return

    do isps = -nhsps(-it), nhsps(it)
       if(isps == 0)cycle
       if( isps < 0)then
          i = -isps
          ith = -it
       else
          i = isps
          ith = it
       endif
       iorb = hspsqn(ith,i)%orb
       if(it == 1 )then
          pPOT_h(isps,isps) = pPOT_h(isps,isps) + pspe(iorb)
       else
          nPOT_h(isps,isps) = nPOT_h(isps,isps) + nspe(iorb)

       end if

    end do

    return
end subroutine makespHFpot
!===================================================
!
! sets up arrays for self-consistent mean-field
!


subroutine setup4MF
	use densities
	use onebodypot
	use haiku_info
	use spstate
	implicit none
	integer :: aerr

    allocate( p1bopme(-nhspsmax:nhspsmax, -nhspsmax:nhspsmax), stat=aerr )
    if(aerr /= 0) then
       call memerror("setupMF 1")
       stop 5
    end if
    allocate( n1bopme(-nhspsmax:nhspsmax, -nhspsmax:nhspsmax), stat=aerr )
    if(aerr /= 0) then
       call memerror("setupMF 2")
       stop 5
    end if
    allocate(pPOT_h(-nhsps(-1):nhsps(1), -nhsps(-1):nhsps(1)), stat=aerr)
    if(aerr /= 0) call memerror("setupMF 3")
    allocate(nPOT_h(-nhsps(-2):nhsps(2), -nhsps(-2):nhsps(2)), stat=aerr)
    if(aerr /= 0) call memerror("setupMF 4")

	return
end subroutine setup4MF
!===================================================
!  added 7.7.4
!
! routine to create and apply mean-field potential
!

subroutine make_MF
	use densities
	use onebodypot
	use mod_reorthog
	implicit none
	
!.......... COPY vec1 TO vec2 FOR COMPUTING DENSITY
     call br_load_vec2_from_vec1()
!..........	COMPUTE DENSITY...............
     call master_density1b_bundled

!...........COMPUTE ONE-BODY POTs...........
     pPOT_h(:,:) = 0.0
     nPOT_h(:,:) = 0.0

     call  makepnHFpot
     call  makexxHFpot(1)
     call  makexxHFpot(2)
     call  makespHFpot(1)
     call  makespHFpot(2)
	 
!	 print*,' potential ',nPOT_h(:,:)

    return	
end subroutine make_MF

!===================================================
!
!  MASTER subroutine for invoking self-consistent mean-field
!
subroutine scmf_master
	
	use nodeinfo
	use lanczos_info,only: nkeep
	use lanczos_util
	use bmpi_mod
	use basis, only : dimbasis
	use io
	use fragments
	use localvectors
	use wfn_mod
	use mod_reorthog
	
	implicit none
	integer :: nloops
	integer :: ierr
	real,allocatable  :: emf(:)
	integer :: ikeep
	
	logical :: write_out_vectors_special  ! 
	
	
	write_out_vectors_special = .false.
	
!.............. MENU....................................
    if(iproc==0)then
        print*,' '
        print*,' / ------------------------------------------------------------------------\ ' 
        print*,' |                                                                         | '
        print*,' |    SELFCONSISTENT MEAN-FIELD OPTIONS (choose one)                       | '
        print*,' |                                                                         | '		
		print*,' |    (lf)  Default Lanczos iteration with mean-field potential (standard) | '
		print*,' |    (ex)  Full diagonalization of mean-field (testing only)              | '
        print*,' |                                                                         | '
        print*,' \ ------------------------------------------------------------------------/ ' 
        print*,' '
        read(5,'(a)')lanczchar
  
        select case (lanczchar)
		    case('ex','Ex','EX')	
			lanczchar = 'ex'
			
			niter = dimbasis   ! not needed, but this initializes the (unused) variable to prevent warning flags
		
		    case default
			
			lanczchar = 'lf'
			
			print*,' Enter number of Lanczos iterations per selfconsistent loop, number of vectors to print out '
			read*,niter,nkeep
		
   	    end select	
		
		print*,' Enter # of selfconsistent loops to use '
		read*,nloops
		
	end if	
    call BMPI_BCAST(lanczchar,2,0,icomm,ierr)
    call BMPI_BCAST(niter,1,0,icomm,ierr)
    call BMPI_BCAST(nloops,1,0,icomm,ierr)
	
	call setup4MF
	call setup_for_lanczos
	allocate(emf(nkeep))
	if(lanczchar == 'ex')then
		call fulldiag_scmf(nloops,emf)
	else
		call lanczos_scmf(nloops,emf)
	end if
	
!........... WRITE OUT FINAL VECTOR TO A .WFN FILE.................	
    call wfn_write_nkeep(nkeep) ! write number of vectors to wfn file
	
	if(write_out_vectors_special)then
		open(unit=37,file='SCMF_vectors.dat',status='unknown')
		write(37,*)dimbasis,nkeep

	end if

    if(writeout)then
		do ikeep = 1,nkeep
			
            call br_retrieve_hist(ikeep)
            call br_restore_vec1()
			call wfn_writeeigenvec(wfnfile,frag1, vec1,ikeep,emf(ikeep),0.0,0.0)
			write(6,*)vec1(1),vec1(dimbasis)
			if(write_out_vectors_special)then
				write(37,*)ikeep
				write(37,*)vec1
			end if
							 
		 end do
	 end if
	if(write_out_vectors_special)close(37)
	
	return

end subroutine scmf_master
!===================================================
!
! modified Lanczos routine for Self-consistent Mean-Field
!
! In this routine:
!   niter = # of Lanczos applications on each loop
!   nloops= # of loops of diagonalization to apply
!

subroutine lanczos_scmf(nloops,emf)
	use nodeinfo
	use mod_reorthog
	use flagger
	use bvectorlib_mod
	use lanczos_util
	use localvectors
	use basis
	use flagger
	use lanczos_info,only:nkeep
	implicit none
	
	integer :: nloops,iloop
    integer(4) :: ierr

    integer(4) :: iter

    real(kind=8) :: da,db
    integer(4) :: i,j,k
    real(kind=8) :: dnorm, dnorm0
	integer :: aerr
    real(kind=egv_prec), allocatable :: h(:,:)
    real(kind=lanc_prec), allocatable :: vi(:), vf(:)
    real(kind=egv_prec), allocatable :: efull(:), u(:,:), work(:),fullvec(:,:)
	logical :: testham,testlanc,altlanc
	logical :: printham_at_end   ! flag to write SCMF to a file at the end
	integer :: info
	real :: emf(nkeep)
	testham = .false.
    testlanc = .false.
	printham_at_end = .false.
	if(testham)then
	
	    allocate(h(dimbasis,dimbasis), stat=aerr)
	    if(aerr /= 0) call memerror("fulldiag_scmf 1")	

	    allocate( efull(dimbasis), work(dimbasis*3),fullvec(dimbasis,dimbasis), stat=aerr)
	    if(aerr /= 0) call memerror("fulldiag_scmf 2")	
	end if

!-----------------------------------------------------------------
	if(.not.allocated(alpha)) allocate( alpha(niter),beta(niter), stat=aerr )

	if(aerr /= 0) then
	     call memerror("lanczos_scmf 1")
	     stop 5
    end if
	if(.not.allocated(e)) allocate( eiglvec(niter,niter),e(niter), stat=aerr )
	if(aerr /= 0) then
	     call memerror("lanczos_scmf 1")
	     stop 5
	end if
!.............. FOR TESTING........	
	if(altlanc)then
	   allocate (templancvecs(dimbasis,niter+1))	
	   call pushvec1(1)
	   call normalizetemp(1,da)
    end if
!.........................

	if(.not.UseNewReorthog)then
		if(iproc==0)print*,' Not using new reorthog? Something wrong '
		stop
		
	end if
    call br_grab_vec1()          ! pushes vec1 to br_reg
    call br_normalize(dnorm0)    ! normalizes br_reg
    call br_restore_vec1()       ! pulls (normalized) vec1 from br_reg
    call br_add2hist(1)         ! adds br_reg to history at position 1
	call make_MF

	do iloop = 1,nloops
		
		call make_MF
!	    call br_normalize(dnorm0)
!		print*,' initial norm ',dnorm0,br_histpos

  	    alpha = 0.0
  	    beta  = 0.0
!............ MORE TESTING....................		
		if (testham)then
		    do i = 1, dimbasis
		       vec1(:) = 0.d0
		       vec1(i) = 1.d0
		       vec2(:) = 0.d0
		       if(useVec2Thread) vec2threadflat(:) = 0.0
			   call master_apply1bop
		       do j = 1,dimbasis
		          h(i,j) = real(vec2(j),kind(egv_prec))   ! store as real(4); this can be changed

		       end do  !J 
		   end do
		end if	
!.............................
		do iter = 1,niter			

!.............. TESTING.................
            if(altlanc)then
			   call normalizetemp(iter,db)   ! this seems to be necessary... but why?	
			   call pullvec1(iter)
		    end if
!................................

			call initialize_final('n')  ! prepares vec2
			
			call master_apply1bop  ! H vec1 = vec2
			 
!.............. TESTING.................			
			if(altlanc)then
  			   call pushvec2(iter+1)			
			   call orthofoot(iter+1,iter,da)
!			print*,iloop,iter,' alpha = ',da ,', norm = ',db
			   if(iter > 1)then
			      do j = iter-1,1,-1
			 	     call orthofoot(iter+1,j,db)
!				write(99,*)iloop,iter+1,j,db
			      end do
   		       end if
			   db = 0.0
			   call normalizetemp(iter+1,db)
		    end if
!.................. 
		    call br_grab_vec2()   ! pushes vec2 to br_reg
!	        if(alpha_before_orthog) ASSUME THIS
	        call br_remove_prev_overlap(da)   ! projects br_reg against previous in history 
	        call br_orthogonalize(1) ! 1 says ignore last history entry (prev)
	        call br_normalize(db)   ! normalize br_reg
			

	        alpha(iter) = real(da,kind(0.0e0))
	        beta(iter) = real(db, kind(0.0e0))
			
!..................... TESTING.........			
			if(altlanc)call normalizetemp(iter+1,db)
			
!			write(98,*)iloop,iter,alpha(iter),beta(iter)
			
!			print*,' alpha beta ',iter, alpha(iter),beta(iter)
			
!............. IMPORTANT .... TEST FOR RANDOM RESTART................

            if(beta(iter) < restart_tol*5 .and. iter < niter)then
					if(iproc==0)print*,iter, ' restarting, beta = ',beta(iter)
					beta(iter)=0.0
					call random_restart_with_history
				
				
			end if
	        ! Finally, add br_reg to history
	        if(iter < niter)then
				 call br_add2hist(iter+1)	! stores br_reg in history as iter+1	
	 	        ! Restore back to vec1 for next iteration
	 	        ! This is instead of swap_vchar in older code
			
	 	        call br_restore_vec1()    ! pulls vec1 from br_reg
			end if
				 	
		end do  ! iter
		
		
!........ TEST BY CHECKING ORTHONORMALITY....
!print*,' testing orthonormality ',niter
!do iter = 1,niter

!	call br_retrieve_hist(iter)
!    call br_normalize(db)
!	print*,iter,db
!end do
!	call br_restore_vec1()
 !   print*,vec1(1)
!	do j = 1,iter
!		call br_retrieve_hist(j)
		
!		call br_restore_vec2()
!		dnorm = 0.0
!		do i = 1,dimbasis
!			dnorm = dnorm+vec1(i)*vec2(i)
!		end do
!		print*,iter,j,dnorm
!	end do
	
!end do		
!......... TEST BY CREATING ENTIRE HAMILTONIAN........		


if (testham)then
    do i = 1, dimbasis
       vec1(:) = 0.d0
       vec1(i) = 1.d0
       vec2(:) = 0.d0
       if(useVec2Thread) vec2threadflat(:) = 0.0
	   call master_apply1bop
	 !         if(nproc > 1) then
           ! leaves data only on isfragroot nodes
           ! In normal use, would move to vector space (breorthog.f90)
           ! need to broadcast
              ! output in vec2 with input 'n'
!              call BMPI_BCAST(vec2, size(vec2), 0, fcomm2, ierr)
!       end if
       do j = 1,dimbasis
          h(i,j) = real(vec2(j),kind(egv_prec))   ! store as real(4); this can be changed

       end do  !J 
   end do
	
	
end if


if(testlanc)then
	print*,'testing'
	do i =1,iter
		call pullvec1(i)
		call pullvec2(i+1)
		da = 0.d0
		db = 0.d0
		do j = 1,dimbasis
			
			do k = 1,dimbasis
				da = da + vec1(j)*h(j,k)*vec1(k)
				db = db + vec1(j)*h(j,k)*vec2(k)
				
			end do
			
		end do
		write(99,*)i,alpha(i),da,beta(i),db
		
	end do
	
end if
!------------------- RESTART WITH INITIAL VECTOR-----------		
!  invoke thick-restart routine to accomplish this
!
!print*,'alphas ',alpha(1),alpha(2),alpha(3),alpha(4)

        call find_lanczos_eigenvalues(niter,.true.,0)
		do i = 1,nkeep
		   emf(i) = e(i)
	    end do
        if(iproc==0)print*,iloop,' lanczos ', e(1)

		if(testham)then
           call DSYEV( 'V','U', dimbasis, h, dimbasis, efull, WORK, 3*dimbasis, INFO )
           if(iproc==0)print*,'exact 2', efull(1),efull(2),efull(3),efull(4),efull(5),efull(6),efull(7),efull(8)
	    end if
        dnorm = 0.d0
		
		if(altlanc)then
		   call constructvec1
		   call pushvec1(1) 
	    endif
	    if(iloop==nloops)then
	        call br_transform_basis(eiglvec,nkeep, niter)    ! only transform initial vector
			
		else
           call br_transform_basis(eiglvec, 1, niter)    ! only transform initial vector
           call br_retrieve_hist(1)   
		
		   call br_normalize(dnorm0)
 
        ! init for next iteration. 
        ! load vec1 from the top of the history.
!        call br_retrieve_hist(1)   

	       call br_restore_vec1()
		   call br_set_histpos(0)
           call br_add2hist(1) ! br_retrieve_hist left in br_reg, put at new head
	   endif
 		
	end do ! iloop
!............... ONLY FOR WRITING OUT SCMF.............	
	if(printham_at_end)then
		open(unit=48,file='scmf.dat',status='unknown')
		print*,' opened file scmf.dat '
		write(48,*)dimbasis
	    do i = 1, dimbasis
	       vec1(:) = 0.d0
	       vec1(i) = 1.d0
	       vec2(:) = 0.d0
	       if(useVec2Thread) vec2threadflat(:) = 0.0
		   call master_apply1bop
		 !         if(nproc > 1) then
	           ! leaves data only on isfragroot nodes
	           ! In normal use, would move to vector space (breorthog.f90)
	           ! need to broadcast
	              ! output in vec2 with input 'n'
	!              call BMPI_BCAST(vec2, size(vec2), 0, fcomm2, ierr)
	!       end if
	
	 
	       do j = 1,i
			   if(vec2(j) /= 0.d0)write(48,*)i,j,vec2(j)
	       end do  !J 
	   end do		
	   
		
		close(48)
		print*,' finished writing scmf.dat '
		
	end if
	
	return
end subroutine lanczos_scmf

!------------------ SUBROUTINES FOR TESTING LANCZOS-----------
!                  TO BE ELIMINATED IN LATER VERSIONS


subroutine pushvec1(pos)
	use precisions
	use basis,only: dimbasis
	use localvectors
	
	implicit none
	
	integer(kind=basis_prec) :: i
	integer pos
	
	do i = 1,dimbasis
		templancvecs(i,pos)=vec1(i)
	end do
	return
	
end subroutine pushvec1
!===================================================
subroutine pushvec2(pos)
	use precisions
	use basis,only: dimbasis
	use localvectors
	
	implicit none
	
	integer(kind=basis_prec) :: i
	integer pos
	
	do i = 1,dimbasis
		templancvecs(i,pos)=vec2(i)
	end do
	return
	
end subroutine pushvec2
!===================================================
subroutine pullvec1(pos)
	use precisions
	use basis,only: dimbasis
	use localvectors
	
	implicit none
	
	integer(kind=basis_prec) :: i
	integer pos
	
	do i = 1,dimbasis
		vec1(i)=templancvecs(i,pos)
	end do
	return
	
end subroutine pullvec1
!===================================================
subroutine pullvec2(pos)
	use precisions
	use basis,only: dimbasis
	use localvectors
	
	implicit none
	
	integer(kind=basis_prec) :: i
	integer pos
	
	do i = 1,dimbasis
		vec2(i) = templancvecs(i,pos)
	end do
	return
	
end subroutine pullvec2

!===================================================
subroutine normalizetemp(pos,dnorm)
	use precisions
	use basis,only: dimbasis
	implicit none
	
	integer(kind=basis_prec) :: i
	integer pos
	real(8) :: dnorm,scale
	
	dnorm = 0.d0
	
	do i = 1,dimbasis
		dnorm = dnorm + templancvecs(i,pos)**2
		
	end do
	dnorm = sqrt(dnorm)
	scale = 1.d0/dnorm
	
	do i = 1,dimbasis
		templancvecs(i,pos)=templancvecs(i,pos)*scale
	end do
	
	return
end subroutine normalizetemp

!===================================================
subroutine orthofoot(pos1,pos2,dnorm)
	use precisions
	use basis,only: dimbasis
	implicit none
	
	integer(kind=basis_prec) :: i
	integer pos1,pos2
	real(8) :: dnorm,scale
	
	dnorm = 0.d0
!	print*,' dot of ',pos1,pos2
	do i = 1,dimbasis
		dnorm = dnorm + templancvecs(i,pos1)*templancvecs(i,pos2)
!		print*,i,templancvecs(i,pos1),templancvecs(i,pos2)
		
	end do

	
	do i = 1,dimbasis
		templancvecs(i,pos1)=templancvecs(i,pos1)-dnorm*templancvecs(i,pos2)
	end do
	
	return
end subroutine orthofoot

!===================================================
subroutine constructvec1
	use precisions
	use basis,only: dimbasis
	use lanczos_util
	use localvectors
	implicit none
	integer(kind=basis_prec) :: i
	real(8) :: dtmp
	integer :: pos
	do i = 1,dimbasis
		dtmp = 0.d0
		do pos = 1,niter
			dtmp = dtmp + templancvecs(i,pos)*eiglvec(pos,1)
			
		end do
		vec1(i)=real(dtmp,kind=lanc_prec)
	
	
	end do
	return
	
end subroutine constructvec1
!===================================================
subroutine ersatz_random_restart(pos)
	
	use precisions
	use basis,only: dimbasis
	use lanczos_util
	use localvectors
	implicit none
	integer(kind=basis_prec) :: i
	real :: rtmp,scale
	
	integer :: pos

	call random_seed

    scale = 1.0/sqrt(real(dimbasis))
	do i = 1,dimbasis	
		call random_number(rtmp)
		vec1(i)=(rtmp-0.5)*scale
	end do
	
	
	return
end subroutine ersatz_random_restart


!===================================================

!
!  FOR TESTING ONLY
!
!  version of SCMF with full diagonalization
!
subroutine fulldiag_scmf(nloops,emf)
	use nodeinfo
	use precisions
	use basis, only : dimbasis
	use localvectors
	use bmpi_mod
	use apply_ham
	use lanczos_info,only:nkeep
	
	implicit none
	
	integer(4) :: nloops 
	
	integer(4) :: iloops
    integer(4) :: ierr
    real(kind=egv_prec), allocatable :: h(:,:)
    real(kind=lanc_prec), allocatable :: vi(:), vf(:)
    real(kind=egv_prec), allocatable :: e(:), u(:,:), work(:),eiglvec(:,:)
    integer(kind=basis_prec) :: i, j, k
    integer :: vidx
    real :: xj,xt	
	integer(4) :: aerr,info
	real :: emf(nkeep)
	
    allocate(h(dimbasis,dimbasis), stat=aerr)
    if(aerr /= 0) call memerror("fulldiag_scmf 1")	

    allocate( e(dimbasis), work(dimbasis*3), stat=aerr)
    if(aerr /= 0) call memerror("fulldiag_scmf 2")	
!    call   procOP_clock(0,'set','all')
	
!---------- INITIAL VECTOR FOR CREATING DENSITIES----------------------	
		call make_MF

	do iloops = 1,nloops
		call make_MF
	    do i = 1, dimbasis
	       vec1(:) = 0.d0
	       vec1(i) = 1.d0
	       vec2(:) = 0.d0
	       if(useVec2Thread) vec2threadflat(:) = 0.0
   		   call master_apply1bop
	       if(nproc > 1) then
	           ! leaves data only on isfragroot nodes
	           ! In normal use, would move to vector space (breorthog.f90)
	           ! need to broadcast
	              ! output in vec2 with input 'n'
	              call BMPI_BCAST(vec2, size(vec2), 0, fcomm2, ierr)
	       end if
	       do j = 1,dimbasis
	          h(i,j) = real(vec2(j),kind(egv_prec))   ! store as real(4); this can be changed

	       end do  !J 
		   
!		   print*,h(i,:)

	    end do  !i
	    print*,' Finished with creating the Hamiltonian '
		
!	    close(60)
!------------------ DIAGONALIZE VIA HOUSEHOLDER -------------------

	    if(egv_prec==4)then
	        call SSYEV( 'V','U', dimbasis, h, dimbasis, e, WORK, 3*dimbasis, INFO )
	    else
	        call DSYEV( 'V','U', dimbasis, h, dimbasis, e, WORK, 3*dimbasis, INFO )
	    end if		
		
!------------------- RESTART WITH INITIAL VECTOR-----------		
!       
        if(iproc == 0 ) print*,e(1),e(2),e(3),e(4)
!		print*,h(:,1)
        do i = 1, dimbasis
			vec1(i)= h(i,1)   ! 
		end do
!		print*,'vec1 :',vec1(:)
		
	end do  ! iloops
	emf(1) = e(1)
	return
	
end subroutine fulldiag_scmf

!===================================================

end module shampoo







