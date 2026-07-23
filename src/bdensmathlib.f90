
!   ADDED BY KSM 7.4.3

module dm_mod
   type dmtab
      integer(kind=4) :: statei, statej
      integer(kind=4) :: jt
      real(kind=4), pointer :: tab(:,:,:)  ! tab(orba, orbb, species)
      type(dmtab), pointer :: next
   end type dmtab

   type(dmtab), pointer :: dmhead => NULL()
contains

! clean up linked list
subroutine dm_reset()
   implicit none
   type(dmtab), pointer :: dp

   do while(associated(dmhead))
      dp => dmhead
      dmhead => dp%next
      deallocate(dp%tab)
      deallocate(dp)
   end do
end subroutine dm_reset

! push new record onto linked list
subroutine dm_add_tab(statei, statej, jt, numorb)
   implicit none
   ! args
   integer :: species, statei, statej, jt,numorb
   ! locals
   real(kind=4), pointer :: tab(:,:,:)
   type(dmtab),  pointer :: dnp

   allocate(dnp)
   allocate(dnp%tab(1:numorb, 1:numorb, 0:1))
   dnp%tab = 0.0
   dnp%statei = statei
   dnp%statej = statej
   dnp%jt = jt
   dnp%next => dmhead
   dmhead => dnp
end subroutine dm_add_tab

!======================================================
! Write out density matrix elements to Mathematica format
! KSM
subroutine density1b_write_math_me(mfn, nkeep, numorbmax, dp)
   implicit none
   ! arguments
   integer :: mfn  ! file
   integer :: species, nkeep, numorbmax
   real(kind=4) :: tab(1:numorbmax,1:numorbmax,0:1)
   type(dmtab), pointer :: dp

   ! locals
   integer :: m, n
   real(kind=4) :: me0, me1
   logical :: wcomma

   write(mfn, '(A)')  "(* densMatElem[statei, statej, jt, tt, {table}]        *)"
   write(mfn, '(A)')  "(* Table entries   {orba, orbb, me0, me1}              *)"
   write(mfn, '(A,I4,A,I4,A,I4,A)') "densMatElem[", &
      dp%statei, ",", dp%statej, ", ", dp%jt,  ", {"
   wcomma = .false.
   do m=1,numorbmax
      do n=1,numorbmax
         me0 = dp%tab(m,n,0)
         me1 = dp%tab(m,n,1)
         if(me0 /= 0.0 .or. me1 /= 0.0) then
            if(wcomma) write(mfn, '(A)') ","
            wcomma = .true.
            write(mfn, '(A,I3,A,I3,A,F12.8,A,F12.8,A)', advance='NO') "   {", m, ", ", n, ", ", me0,", ", me1, "}"
         end if
      end do
   end do
   write(mfn, '(A)') " " ! newline
   write(mfn, '(A)') "}];"
end subroutine density1b_write_math_me

subroutine density1b_write_math(nkeep, stateE, stateJ, stateT)
   use nodeinfo
   use sporbit
   use io
   use menu_choices
   use pocc_mod
   use system_parameters
   use interaction
   implicit none
   ! arguments
   integer :: nkeep
   real(kind=4) :: dm(0:1, 1:nkeep, 1:nkeep, 1:numorbmax, 1:numorbmax)
   real(kind=4) :: stateE(1:nkeep), stateJ(1:nkeep), stateT(1:nkeep)
   ! locals
   integer :: i, j, m, n
   logical :: wcomma
   integer :: mfn ! file
   character(len=1024) :: mdenspath
   type(dmtab), pointer :: dp

  ! open mathematica format output file if requested
  mfn = -1  ! test  mfn > 0 for writing below
  if(iproc == 0 .and. menu_dx_omathematica) then
     i = index(outfile, ' ') - 1
     mdenspath = outfile(1:i) // '_dres.m'
     open(unit=mdensfile, file=mdenspath,  action='WRITE', status='REPLACE', err=101)
     mfn = mdensfile
101  continue
  end if
  if(mfn > 0) then
     write(mfn, "(A)")  "(* Bigstick density matrix file *)"
     if(isoflag) then
        write(mfn, '(A)') "densSetIso[True];"
     else
        write(mfn, '(A)') "densSetIso[False];"
     endif
     call pocc_write_orbits_m(mfn)
     write(mfn, '(A,I4,A,I4,A)') "densSetNumParticles[", np(1), ",", np(2), "];"
     write(mfn, '(A)')  "densSetStates[{ (* {#, E, J*2, T*2} *)"
     wcomma = .false.
     do i=1,nkeep
        if(wcomma) write(mfn, '(A)') ","
        wcomma = .true.
        write(mfn, '(A,I3,A,F14.8,A,F12.5,A,F12.5,A)', advance='NO') &
          "   {", i, ",", stateE(i), ",", stateJ(i), ",", stateT(i), "}"
     end do
     write(mfn, '(A)')  " "  ! final newline
     write(mfn, '(A)')  "}];"
     dp => dmhead
     do while(associated(dp))
        call density1b_write_math_me(mfn, nkeep, numorbmax, dp)
        dp => dp%next
     end do
  end if
  if(mfn > 0) close(mfn)

end subroutine density1b_write_math

end module dm_mod
