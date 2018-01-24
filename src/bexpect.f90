
!=============================================================
!   subroutine expectator_p
!       master subroutine to control calculation of expectation values
!
!  NOTE: There will have to be modifications when sections of the lanczos vector are broken into pieces
!
!  CALLS
!   setup_localvectors
!   read_wfn_nkeep
!   dnormvec_p
!   wfn_readeigenvec
!    dnormvec_p
!  br_load_vec2_from_vec1
!   applyobsbundled
subroutine expectator_p

  use precisions
  use basis
  use io
  use localvectors
  use obs
  use fragments
  use nodeinfo
  use flags3body
  use mod_reorthog
  use wfn_mod
  use bmpi_mod
  use bvectorlib_mod
  use lanczos_util
  use apply_obs_mod
  implicit none

  integer :: ikeep
  integer :: i,j,n
  real :: xe,xj,xt
  real (kind = lanc_prec) :: exph
  real (kind = 8) :: dnorm
  logical :: zeroflag
  integer :: ierr

  !.....NOTE: LIMITED TO 2-BODY OBSERVABLES....CAN BE FIXED BUT MUST ADD ROUTINES TO applybsbundled 

  if(iproc==0 .and. threebody)print*,' Expectation values not yet set up for 3-body forces'
!........... TEMP...... WILL NEED TO FIX LATER....
 
  call setup_localvectors
  call wfn_read_nkeep(oldwfnfile, nkeep) ! does BCAST
  if(iproc==0)then
       print*,nkeep,' states '
       if(writeout)write(resultfile,'(''nkeep ='',1x,i10)')nkeep
       print*,' '
       print*,' STATE       E          J           <H >       (norm)'
       if(writeout)write(resultfile,*)'  STATE      E            J           T^2           <H >      (norm)'
  end if
  twoobsflag = .false.  ! may change
  do ikeep = 1,nkeep
     i = ikeep
     ! new interface, we say which vector to read - it checks
     call wfn_readeigenvec(oldwfnfile, frag1, fcomm1, vec1,i,xe,xj,xt) ! KSM: updated
     call dnormvec_p('n','i',dnorm,zeroflag)
     call br_load_vec2_from_vec1()
     xj2 = 0.0e0_lanc_prec
     call applyobsbundled(1)
     if(iproc==0)then
         write(6,101)i,xe,xj,xj2, dnorm
!         if(writeout)write(resultfile,101)i,e,xj,xj2,dnorm
         if(writeout)write(resultfile,101)i,xe,xj,xt,xj2,dnorm
101  format(i5,2x,3(1x,f11.4),1x,f13.6,1x,f10.5)
     end if
  end do ! ikeep
  
  return
end subroutine expectator_p

!==================================================================
