!
! routines for reading in J-coupled 3-body forces
! modified from routines supplied by P. Navratil @ TRIUMF 9/2012
!
!
      module v3b
      real(kind(0.0)),allocatable :: v3b_cJ(:)   ! must be defined in v3b module
      integer,allocatable :: index_abc(:,:,:),index_abcdef(:,:), &
          start_abcdef(:),nlj_orb(:) ! must be defined in v3b module
      real(kind=kind(0.0)),allocatable :: cgj12(:,:,:,:,:),  &
          cgj123(:,:,:,:,:)
      real(kind=kind(0.0)) :: cgt12(0:1,0:1,0:1), & 
          cgt123(0:1,0:1,-1:1,0:1)
      end module v3b

      module spbasis_3BJC
      integer :: N1_max,isp1ntot,isp3ntot,dim_abc,dim_abcdef
      integer,allocatable :: isp1n_cJ(:,:),isp3n_cJ(:,:),sp3ntot(:,:,:)
      integer,allocatable :: map2orb(:), mapfromorb(:)

      integer, allocatable,target :: map2psps(:), map2nsps(:)   ! maps to "standard" bigstick s.p. states
                                                           ! protons and neutrons
      integer, allocatable,target :: cmap2psps(:), cmap2nsps(:)  ! maps to s.p. states with opposite m

      integer, allocatable :: mapXXX(:), phaseXXX(:)
      integer, allocatable :: cmapXXX(:), cphaseXXX(:)   ! maps to 3-body states with all M's reversed
      integer, allocatable :: tphaseXXX(:)  ! phase that arises from time-reversal of all M's
                                              ! (-1)**(j1+j2+j3 + 1)

      end module spbasis_3BJC

!===================================================================================
      subroutine fetch3bodyinputcoupled

      use nodeinfo
      use flags3body
      use v3b
      use io
      use interactions3body
	  use bmpi_mod
      implicit none

      character(len=80) :: input_file_name
      integer ierr

      if(.not.threebodycheck .or. .not.threebody) return
      scale3body = 0.0
      if(iproc == 0)then
         write(6,*)' Enter filename for 3-body force '
         write(6,*)' (if none then enter "none"; will convert 2-body to 3-body)'
         write(6,*)' (If filename is long, store in file FILENAME3BODY.TXT '
         write(6,*)'  and it will be auto read from there; enter "auto") '
         if(auto_input)then
            read(autoinputfile,'(a)')input_file_name
            print*,input_file_name
         else
            read(5,'(a)')input_file_name
            print*,input_file_name

         endif
      end if
  
      call BMPI_BARRIER(icomm,ierr)

      call BMPI_BCAST(input_file_name,12,0,icomm,ierr)
      select case (input_file_name(1:4))

      case ('none','NONE')   ! no 3-body
       if(iproc==0)then 
         print*,' no 3-body input '
         if(.not.auto_input)write(autoinputfile,'(a)')input_file_name
       endif
       return

      case ('auto','AUTO')
       threebody = .true.
       open(unit=4,file='FILENAME3BODY.TXT',status='old',err=111)
       read(4,'(a)')input_file_name

       goto 112

111   continue
       if(iproc==0)print*,' You have to create FILENAME3BODY.TXT '
       stop

      case default
       threebody = .true.
   
      end select
112   continue
      if(.not.auto_input .and. iproc==0)write(autoinputfile,'(a)')input_file_name

      if(iproc == 0)then
       if(auto_input)then
           read(autoinputfile,*)scale3body
       else
          write(6,*)' Enter scaling strength for 3-body '
          read(5,*)scale3body
          write(autoinputfile,*)scale3body
       endif
      end if 

      call threebodysetup_cJ(input_file_name)

      return
      end subroutine fetch3bodyinputcoupled
!===================================================================================
      subroutine master3bodyuncoupler

      use interactions3body
      implicit none
      if(scale3body == 0.0)return
      call uncoupleXXX(1)
      call uncoupleXXX(2)
      call uncoupleXXY(1)
      call uncoupleXXY(2)

      return
      end subroutine master3bodyuncoupler

!===================================================================================
!
! CALLs function ham3b_cJ
!
  subroutine uncoupleXXX(it)

  use system_parameters
  use haiku_info
  use interactions3body
  use ntuple_info
  use spstate
  implicit none
  integer it
  integer m, mstart,mend
  integer par
  integer tripstart,tripend,itrip,ftrip
  integer iref3,istart3
  integer a,b,c,d,e,f
  integer ath,bth,cth,dth,eth,fth
  integer i3bme
  real, pointer :: hmatXXX(:)
  real v3me
  real diagfactor
  real ham3b_cJ    ! a function

  if(np(it) < 3)return

  if(it ==1)then
      hmatXXX => hmatPPP  
  else
      hmatXXX => hmatNNN
  end if

  mstart = XXX3(it)%jzstart !trip(1)%m
  mend   = XXX3(it)%jzend   ! trip(nXXX(it) )%m

  do m = mstart,mend,2
     do par = 1,2
        tripstart = nXXX(it)+1
        tripend   = 0
        do itrip = 1,nXXX(it)
           if( XXX3(it)%trip(itrip)%m == m .and. XXX3(it)%trip(itrip)%par==par)then
              tripstart = MIN(tripstart,itrip)
              tripend   = MAX(tripend,itrip)
           end if
        end do   ! itrip
        if(tripstart > tripend)cycle
        iref3 = XXX3(it)%meref(m,par)
        istart3 = XXX3(it)%mestart(m,par)
        do itrip =tripstart,tripend
           a = XXX3(it)%trip(itrip)%ia
           b = XXX3(it)%trip(itrip)%ib
           c = XXX3(it)%trip(itrip)%ic
           if(a < 0) then
                a = -a
                ath = -it
           else
                ath = it
           end if
           if(b < 0) then

                b = -b
                bth = -it
           else
                bth = it
           end if
           if(c < 0) then
                c = -c
                cth = -it
           else
                cth = it
           end if
           do ftrip = itrip,tripend
!................... BECAUSE HERMITICITY EITHER WAY, DIAGONALS APPLIED TWICE SO MULTIPLY BY 1/2 .....
                if( ftrip == itrip )then
                    diagfactor = 0.5
                else
                    diagfactor = 1.0
                end if
                d = XXX3(it)%trip(ftrip)%ia
                if(d < 0) then
                     d = -d
                    dth = -it
                else
                    dth = it
                end if
                e = XXX3(it)%trip(ftrip)%ib
                if(e < 0) then
                    e = -e
                    eth = -it
                else
                    eth = it
                end if
                f = XXX3(it)%trip(ftrip)%ic
                if(f < 0) then
                     f = -f
                     fth = -it
                else
                    fth = it
                end if

                i3bme = istart3 + (ftrip-iref3)*(ftrip-iref3-1)/2+itrip-iref3
                v3me = ham3b_cJ(a,b,c,d,e,f,ath,bth,cth,dth,eth,fth)
                hmatXXX(i3bme) = hmatXXX(i3bme) + v3me*diagfactor*scale3body

           end do
        end do  ! itrip


     end do
  end do ! m



  return
  end subroutine uncoupleXXX

!===================================================================================

  subroutine uncoupleXXY(itx)

  use system_parameters
  use haiku_info
  use interaction

  use interactions3body
  use ntuple_info
  implicit none
  integer itx,ity
  integer m, mstart,mend
  integer istart,ifinish
  integer i,iend
  integer parstart,parend
  integer par
  integer iref3,istart3
  integer a,b,c,d,e,f
  integer ath,bth,cth,dth,eth,fth 
  real diagfactor
  logical semidiagflag  ! if a, b = d,e but c /= f
  integer itrip,ftrip
  integer iXXY, fXXY
  integer iXXpair,fXXpair
   integer,pointer :: cXXYtriplet(:,:), dXXYtriplet(:,:)

  integer i3bme
  real, pointer :: hmatXXY(:)
  real v3me
  real ham3b_cJ    ! a function

  integer(kind=8):: nxyzmat

  ifinish = -1
  iend = -1
  ity = 3-itx

  if(np(itx) < 2 .or. np(ity) < 1 )return

  if(itx ==1)then
      hmatXXY => hmatPPN  
      cXXYtriplet => cppntriplet
      dXXYtriplet => dppntriplet
  else
      hmatXXY => hmatPNN
      cXXYtriplet => cpnntriplet
      dXXYtriplet => dpnntriplet
  end if


  mstart = XXY(itx)%trip(1)%m
  mend   = XXY(itx)%trip(nXXY(itx))%m

   istart  = 1
   do m = mstart,mend, 2   
      if(XXY(itx)%trip(istart)%m /= m)then
           print*,' mismatch ',m,XXY(itx)%trip(istart)%m
           stop
      endif
!......given m, search for start,finish of parity
      if(istart == nXXY(itx))then
         ifinish = nXXY(itx)
      else
         do i = istart+1,nXXY(itx)
            ifinish = i
            if(XXY(itx)%trip(i)%m /= m)then
	       ifinish = i-1
               exit
            endif
         enddo
      endif
      parstart = XXY(itx)%trip(istart)%par
      parend   = XXY(itx)%trip(ifinish)%par
      
      do par = parstart,parend
!............ NOW FIND ACTUAL START, FINISH FOR THIS M, PARITY
         if(istart > ifinish)cycle
         if(istart == ifinish)then
	    iend = ifinish
          else
            do i = istart+1,ifinish
               iend = i
               if(XXY(itx)%trip(i)%par /= par)then
                  iend = i-1
                  exit
               endif
            enddo
         endif
!............ NOW LOOP OVER XXY TRIPLETS
         do itrip = istart,iend
             iXXpair = XXY(itx)%trip(itrip)%ia
             c = XXY(itx)%trip(itrip)%ib
             iXXY =  cXXYtriplet(iXXpair,c) 
             if(c < 0) then
                c = -c
                cth = -ity
             else
                cth = ity
             end if
             a = XX2(itx)%pair(iXXpair)%ia
             if(a < 0)then
                 a = -a
                 ath = -itx
             else
                 ath = itx
             end if
             b = XX2(itx)%pair(iXXpair)%ib
             if(b < 0)then
                 b = -b
                 bth = -itx
             else
                 bth = itx
             end if
             do ftrip = istart,iend


                fXXpair = XXY(itx)%trip(ftrip)%ia
                f = XXY(itx)%trip(ftrip)%ib
                fXXY =  dXXYtriplet(fXXpair,f) 
                if(f < 0) then
                   f = -f
                   fth = -ity
                else
                   fth = ity
                end if

                d = XX2(itx)%pair(fXXpair)%ia
                if(d < 0)then
                   d = -d
                   dth = -itx
                else
                   dth = itx
                end if
                e = XX2(itx)%pair(fXXpair)%ib
                if(e < 0)then
                    e = -e
                    eth = -itx
                else
                    eth = itx
                end if
                i3bme = iXXY+fXXY

                semidiagflag = .false.
                if( iXXpair == fXXpair .and.  XXY(itx)%trip(ftrip)%ib/=  XXY(itx)%trip(itrip)%ib)then
                    semidiagflag = .true.
                end if
                if(itrip==ftrip .or. semidiagflag)then
                     diagfactor = 0.5
                else
                     diagfactor = 1.0
                end if
                v3me = ham3b_cJ(a,b,c,d,e,f,ath,bth,cth,dth,eth,fth)
                hmatXXY(i3bme) = hmatXXY(i3bme) + v3me*diagfactor*scale3body
             end do ! ftrip
        end do !  itrip

!.......................
         istart = iend + 1
      end do ! par
  end do ! m

  return
  end subroutine uncoupleXXY

!===================================================================================
      real function ham3b_cJ(a,b,c,d,e,f,ath,bth,cth,dth,eth,fth)


!      real function ham3b_cJ(ispcr_1_in,ispcr_2_in,ispcr_3_in, &
!          ispan_1_in,ispan_2_in,ispan_3_in,ith1c,ith2c,ith3c,ith1a,ith2a,ith3a)
!      use parameters
!      use spb
!      use spsdata
!      use spbasis_3BJC, only: isp1ntot
      use spstate
      use haiku_info
!      use hamc
      use v3b
      use spbasis_3BJC
      implicit none
!      integer,intent(IN) :: ispcr_1_in,ispcr_2_in,  &
!                ispcr_3_in,ispan_1_in,ispan_2_in,ispan_3_in
      integer :: ith1a,ith1c,ith2a,ith2c,ith3a,ith3c ! species/handedness for particle
      integer :: it1,it2,it3  ! species for particles 1,2,3
      integer :: mt1,mt2,mt3  ! isospin for particles 1,2,3
      integer :: ispcr_1,ispcr_2, &
                ispcr_3,ispan_1,ispan_2,ispan_3  
      integer :: iphase,itemp,mjtot2,mttot2,itemp_i !,ipar
      integer:: i_a,i_b,i_c,i_d,i_e,i_f,i_abc,i_def,i_abcdef,iii
!     $     ,order_abc,order_def
      real(kind(0.d0)) :: sum,cg1,cg2,cg3,cg4,cgt1(0:1),cgt2(0:1), &
          cgt3(0:1,0:1),cgt4(0:1,0:1),clebd,sum3
      integer :: j2_a,j2_b,j2_c,j2_d,j2_e,j2_f,j_ab,j_de,t_ab,t_de, &
          J_3,T_3,m2_a,m2_b,m2_c,m2_d,m2_e,m2_f,mt2_a,mt2_b,mt2_c,  &
          mt2_d,mt2_e,mt2_f
      integer :: isp1n_i,isp1n_j,isp1n_l,isp1n_m,isp1n_n,isp1n_k
      integer :: a,b,c,d,e,f,ath,bth,cth,dth,eth,fth
      real(kind=kind(0.0)) :: ham

      ham3b_cJ = 0.0
!.......... COPY OVER TO KEEP ORIGINALS FROM GETTING SWAPPED
      ispcr_1=a
      ispcr_2=b
      ispcr_3=c
      ispan_1=d
      ispan_2=e
      ispan_3=f
      ith1c = ath
      ith2c = bth
      ith3c = cth
      ith1a = dth
      ith2a = eth
      ith3a = fth
      it1 = abs(ith1a)
      it2 = abs(ith2a)
      it3 = abs(ith3a)
      mt1 = (3-it1*2)
      mt2 = (3-it2*2)
      mt3 = (3-it3*2)
      mttot2 = mt1 + mt2 + mt3 

      mjtot2 = hspsqn(ith1a,ispan_1)%m + hspsqn(ith2a,ispan_2)%m + hspsqn(ith3a,ispan_3)%m

      i_a = mapfromorb( hspsqn(ith1c,ispcr_1)%orb )
      i_b = mapfromorb( hspsqn(ith2c,ispcr_2)%orb )
      i_c = mapfromorb( hspsqn(ith3c,ispcr_3)%orb )
      i_d = mapfromorb( hspsqn(ith1a,ispan_1)%orb )
      i_e = mapfromorb( hspsqn(ith2a,ispan_2)%orb )
      i_f = mapfromorb( hspsqn(ith3a,ispan_3)%orb )

      iphase=1
      if (i_b>i_a) then
         call swap3(ispcr_1,ith1c,i_a,ispcr_2,ith2c,i_b)
         iphase=-iphase
      endif
      if (i_c>i_b) then

         call swap3(ispcr_2,ith2c,i_b,ispcr_3,ith3c,i_c)

         iphase=-iphase
         if (i_b>i_a) then
            call swap3(ispcr_1,ith1c,i_a,ispcr_2,ith2c,i_b)
            iphase=-iphase
         endif
      endif
      if (i_e>i_d) then
         call swap3(ispan_1,ith1a,i_d,ispan_2,ith2a,i_e)
         iphase=-iphase
      endif
      if (i_f>i_e) then
         call swap3(ispan_2,ith2a,i_e,ispan_3,ith3a,i_f)
         iphase=-iphase
         if (i_e>i_d) then
            call swap3(ispan_1,ith1a,i_d,ispan_2,ith2a,i_e)
            iphase=-iphase
         endif
      endif

      i_abc=index_abc(i_a,i_b,i_c)
      i_def=index_abc(i_d,i_e,i_f)
      if (i_def>i_abc) then
         j2_a =hspsqn(ith1a,ispan_1)%j
         j2_b =hspsqn(ith2a,ispan_2)%j
         j2_c =hspsqn(ith3a,ispan_3)%j
         j2_d =hspsqn(ith1c,ispcr_1)%j
         j2_e =hspsqn(ith2c,ispcr_2)%j
         j2_f =hspsqn(ith3c,ispcr_3)%j

         m2_a =hspsqn(ith1a,ispan_1)%m
         m2_b =hspsqn(ith2a,ispan_2)%m
         m2_c =hspsqn(ith3a,ispan_3)%m
         m2_d =hspsqn(ith1c,ispcr_1)%m
         m2_e =hspsqn(ith2c,ispcr_2)%m
         m2_f =hspsqn(ith3c,ispcr_3)%m

         mt2_a = (3- abs(ith1a)*2)  !  mt1
         mt2_b = (3- abs(ith2a)*2)  !mt2
         mt2_c =(3- abs(ith3a)*2)  !mt3
         mt2_d = (3- abs(ith1c)*2)  !mt1
         mt2_e = (3- abs(ith2c)*2)  !mt2
         mt2_f =(3- abs(ith3c)*2)  !mt3
         i_abcdef=index_abcdef(i_def,i_abc)
      else
         j2_d =hspsqn(ith1a,ispan_1)%j
         j2_e =hspsqn(ith2a,ispan_2)%j
         j2_f =hspsqn(ith3a,ispan_3)%j
         j2_a =hspsqn(ith1c,ispcr_1)%j
         j2_b =hspsqn(ith2c,ispcr_2)%j
         j2_c =hspsqn(ith3c,ispcr_3)%j

         m2_d =hspsqn(ith1a,ispan_1)%m
         m2_e =hspsqn(ith2a,ispan_2)%m
         m2_f =hspsqn(ith3a,ispan_3)%m
         m2_a =hspsqn(ith1c,ispcr_1)%m
         m2_b =hspsqn(ith2c,ispcr_2)%m
         m2_c =hspsqn(ith3c,ispcr_3)%m

         mt2_a = (3- abs(ith1c)*2)  !  mt1
         mt2_b = (3- abs(ith2c)*2)  !mt2
         mt2_c =(3- abs(ith3c)*2)  !mt3
         mt2_d = (3- abs(ith1a)*2)  !mt1
         mt2_e = (3- abs(ith2a)*2)  !mt2
         mt2_f =(3- abs(ith3a)*2)  !mt3

         i_abcdef=index_abcdef(i_abc,i_def)
      endif

      cgt1(0:1)=cgt12(0:1,(mt2_a+1)/2,(mt2_b+1)/2)
      cgt2(0:1)=cgt12(0:1,(mt2_d+1)/2,(mt2_e+1)/2)

      cgt3(0:1,0:1)=cgt123(0:1,0:1,(mt2_a+mt2_b)/2,(mt2_c+1)/2)
      cgt4(0:1,0:1)=cgt123(0:1,0:1,(mt2_d+mt2_e)/2,(mt2_f+1)/2)

      sum3=0.d0
      iii=start_abcdef(i_abcdef)
      if(m2_a + m2_b + m2_c /= m2_d + m2_e + m2_f) then
        print*,' mismatch in ms '
        print*, m2_a, m2_b,m2_c, m2_d, m2_e,m2_f 
        print*,hspsqn(ath,a)%m, hspsqn(bth,b)%m,hspsqn(cth,c)%m, hspsqn(dth,d)%m,hspsqn(eth,e)%m,hspsqn(fth,f)%m
        write(6,'(12i3)')ath,a,bth,b,cth,c,dth,d,eth,e,fth,f
        write(6,'(12i3)')ith1c,ispcr_1,ith2c,ispcr_2,ith3c,ispcr_3,ith1a,ispan_1,ith2a,ispan_2,ith3a,ispan_3
        stop
      end if
      do j_ab=abs(j2_a-j2_b),(j2_a+j2_b),2
            cg1=cgj12(j_ab/2,j2_a/2,(j2_a+m2_a)/2,j2_b/2,(j2_b+m2_b)/2)
         do j_de=abs(j2_d-j2_e),(j2_d+j2_e),2
               cg2=cgj12(j_de/2,j2_d/2,(j2_d+m2_d)/2, &
                   j2_e/2,(j2_e+m2_e)/2)
            do J_3=MAX(abs(j_ab-j2_c),  &
                abs(j_de-j2_f)),        &
                MIN(j_ab+j2_c,j_de+j2_f),2  

                  cg3=cgj123(J_3/2,j_ab/2,(m2_a+m2_b)/2, &
                      j2_c/2,(j2_c+m2_c)/2)
                  cg4=cgj123(J_3/2,j_de/2,(m2_d+m2_e)/2, &
                      j2_f/2,(j2_f+m2_f)/2)
  
               sum3=sum3+cg1*cg2*cg3*cg4*(                             &
                   v3b_cJ(iii)*cgt1(0)*cgt2(0)*cgt3(0,0)*cgt4(0,0)    &
                   +v3b_cJ(iii+1)*cgt1(0)*cgt2(1)*cgt3(0,0)*cgt4(1,0) &
                   +v3b_cJ(iii+2)*cgt1(1)*cgt2(0)*cgt3(1,0)*cgt4(0,0) &
                   +v3b_cJ(iii+3)*cgt1(1)*cgt2(1)*cgt3(1,0)*cgt4(1,0) &
                   +v3b_cJ(iii+4)*cgt1(1)*cgt2(1)*cgt3(1,1)*cgt4(1,1))
               iii=iii+5

            end do
         end do
      end do

      ham3b_cJ = sum3*real(iphase,kind(0.d0))

      end

!-------------------------------------------------------------
     subroutine swap3(n1,n2,n3,m1,m2,m3)

     implicit none
     integer n1,n2,n3,m1,m2,m3
     integer temp

     temp = m1
     m1 = n1
     n1=temp

     temp = m2
     m2 = n2
     n2=temp

     temp = m3
     m3 = n3
     n3=temp

     return
     end subroutine swap3

!-------------------------------------------------------------
!     any N_a <=  nhom1sp
!     any N_a + N_b <= nhom2sp
!     any N_a + N_b + N_c <= nhom

      subroutine threebodysetup_cJ(v3intfile)
!** Setup for a three-body interaction ***
!** Petr Navratil, 30 June 2011, TRIUMF *****
!
!  SUBROUTINES CALLED: 
!      sp1nbas_cJ
!      cg_init
!
!      use parameters
      use nodeinfo
!      use spb
!      use spsdata
      use spbasis_3BJC
      use v3b
      use system_parameters
      use W_info
      use sporbit
      use spstate
!      use hamc
	  use bmpi_mod
      implicit none
      integer :: nhom2sp
      character(len=80),intent(IN) :: v3intfile
      integer :: J_3,T_3,tot_num_of_3bmatel,ini,fin
      integer :: minz,nspsmin,nhom,ierr

      integer :: ntot_a,n_a,l_a,j2_a,isp1n_a
      integer :: ntot_b,n_b,l_b,j2_b,isp1n_b
      integer :: ntot_c,n_c,l_c,j2_c,isp1n_c
      integer :: ntot_d,n_d,l_d,j2_d,isp1n_d
      integer :: ntot_e,n_e,l_e,j2_e,isp1n_e,isp1n_e_max
      integer :: ntot_f,n_f,l_f,j2_f,isp1n_f,isp1n_f_max
      integer :: ii,j_ab,t_ab,j_de,t_de,i_abc,i_def,i_abcdef
      integer(4) :: ih3,ih3oibuf,iix
      integer :: i,ibuf,pi_i,pi_f
      integer iii
      real vx
      integer :: Nmin_HO    ! must be a function
      integer :: aerr

      minz=Nmin_HO(np(1))+Nmin_HO(np(2))
      nspsmin=MAX(minz-3,0)
      if (iproc==0) print *,' nspsmin=',nspsmin
      nspsmin=MIN(Nmin_HO(np(1)-3)+Nmin_HO(np(2) ), &
          Nmin_HO(np(1)-2)+Nmin_HO(np(2)-1),       &
          Nmin_HO(np(1)-1)+Nmin_HO(np(2)-2),       &
          Nmin_HO(np(1))+Nmin_HO(np(2)-3))
      if (iproc==0) print *,' nspsmin=',nspsmin
!      nhom=nhw-nspsmin                             !nhw = excitation???
      
      nhom = maxWtot-nspsmin
      nspsmin=MIN(Nmin_HO(np(1)-2)+Nmin_HO(np(2) ), &
          Nmin_HO(np(1)-1)+Nmin_HO(np(2)-1),       &
          Nmin_HO(np(1))+Nmin_HO(np(2)-2) )
      nhom2sp = maxWtot - nspsmin
      if(iproc==0)print*,' N12_max = ',nhom2sp

      if (iproc==0) print *,' N123_max=',nhom
!*** test
!      nhom=MIN(nhom,3*(2*n_sp(nasps)+l_sp(nasps)))
!      nhom=MAX(nhom,nhom2sp) ! good only for truncating 6Li,4He inter. use rather nhom=MAX(nhom,N123max)
!*** test

!      N1_max=nshll-1
!      N1_max = numorb(1) -1 
      nspsmin = MIN( Nmin_HO(np(1)-1) + Nmin_HO(np(2)), Nmin_HO(np(1))+Nmin_HO(np(2) -1) )
      N1_max = maxWtot-nspsmin
      call sp1nbas_cJ(N1_max,9)

      allocate(index_abc(isp1ntot,isp1ntot,isp1ntot), stat=aerr)
      if(aerr /= 0) call memerror("threebodysetup_cJ 1")
      ii=0
      do isp1n_a=1,isp1ntot
         n_a=isp1n_cJ(1,isp1n_a)
         l_a=isp1n_cJ(2,isp1n_a)
         j2_a=isp1n_cJ(3,isp1n_a)
         ntot_a=2*n_a+l_a
         if (ntot_a>nhom) cycle
         do isp1n_b=1,isp1n_a
            n_b=isp1n_cJ(1,isp1n_b)
            l_b=isp1n_cJ(2,isp1n_b)
            j2_b=isp1n_cJ(3,isp1n_b)
            ntot_b=2*n_b+l_b
            if (ntot_a+ntot_b>nhom2sp) cycle
            do isp1n_c=1,isp1n_b
               n_c=isp1n_cJ(1,isp1n_c)
               l_c=isp1n_cJ(2,isp1n_c)
               j2_c=isp1n_cJ(3,isp1n_c)
               ntot_c=2*n_c+l_c
               if (ntot_a+ntot_c>nhom2sp) cycle
               if (ntot_b+ntot_c>nhom2sp) cycle
               if (ntot_a+ntot_b+ntot_c>nhom) cycle
               ii=ii+1
               index_abc(isp1n_a,isp1n_b,isp1n_c)=ii
            end do
         end do
      end do
      dim_abc=ii
!      print *,' dim_abc=',dim_abc
      allocate(index_abcdef(dim_abc,dim_abc), stat=aerr)
      if(aerr /= 0) call memerror("threebodysetup_cJ 2")
      index_abcdef=0
      ii=0
      do isp1n_a=1,isp1ntot
         n_a=isp1n_cJ(1,isp1n_a)
         l_a=isp1n_cJ(2,isp1n_a)
         j2_a=isp1n_cJ(3,isp1n_a)
         ntot_a=2*n_a+l_a
         if (ntot_a>nhom) cycle
         do isp1n_b=1,isp1n_a
            n_b=isp1n_cJ(1,isp1n_b)
            l_b=isp1n_cJ(2,isp1n_b)
            j2_b=isp1n_cJ(3,isp1n_b)
            ntot_b=2*n_b+l_b
            if (ntot_a+ntot_b>nhom2sp) cycle
            do isp1n_c=1,isp1n_b
               n_c=isp1n_cJ(1,isp1n_c)
               l_c=isp1n_cJ(2,isp1n_c)
               j2_c=isp1n_cJ(3,isp1n_c)
               ntot_c=2*n_c+l_c
               if (ntot_a+ntot_c>nhom2sp) cycle
               if (ntot_b+ntot_c>nhom2sp) cycle
               if (ntot_a+ntot_b+ntot_c>nhom) cycle
               i_abc=index_abc(isp1n_a,isp1n_b,isp1n_c)
!               order_abc=(isp1n_a-1)*((isp1ntot-1)*(isp1ntot+1)+1)
!     $              +(isp1n_b-1)*isp1ntot+isp1n_c-1

!               do isp1n_d=1,isp1ntot
               do isp1n_d=1,isp1n_a
                  n_d=isp1n_cJ(1,isp1n_d)
                  l_d=isp1n_cJ(2,isp1n_d)
                  j2_d=isp1n_cJ(3,isp1n_d)
                  ntot_d=2*n_d+l_d
                  if (ntot_d>nhom) cycle
                  if (isp1n_d==isp1n_a) then
                     isp1n_e_max=isp1n_b
                  else
                     isp1n_e_max=isp1n_d
                  endif
                  do isp1n_e=1,isp1n_e_max
                     n_e=isp1n_cJ(1,isp1n_e)
                     l_e=isp1n_cJ(2,isp1n_e)
                     j2_e=isp1n_cJ(3,isp1n_e)
                     ntot_e=2*n_e+l_e
                     if (ntot_d+ntot_e>nhom2sp) cycle
                     if (isp1n_d==isp1n_a.and.isp1n_e==isp1n_b) then
                        isp1n_f_max=isp1n_c
                     else
                        isp1n_f_max=isp1n_e
                     endif
                     do isp1n_f=1,isp1n_f_max
                        n_f=isp1n_cJ(1,isp1n_f)
                        l_f=isp1n_cJ(2,isp1n_f)
                        j2_f=isp1n_cJ(3,isp1n_f)
                        ntot_f=2*n_f+l_f
                        if (ntot_d+ntot_f>nhom2sp) cycle
                        if (ntot_e+ntot_f>nhom2sp) cycle
                        if (ntot_d+ntot_e+ntot_f>nhom) cycle
                        if (mod(l_a+l_b+l_c+l_d+l_e+l_f,2)==1) cycle
                        i_def=index_abc(isp1n_d,isp1n_e,isp1n_f)

                        if (i_def>i_abc) then
                           print *,'a,b,c=',isp1n_a,isp1n_b,isp1n_c
                           print *,'e,d,f=',isp1n_e,isp1n_d,isp1n_f

                           print *,' i_abc=',i_abc
                           print *,' i_def=',i_def
                           stop
                        endif
                        ii=ii+1
                        index_abcdef(i_abc,i_def)=ii
                     end do
                  end do
               end do
            end do
         end do
      end do
      dim_abcdef=ii
!      print *,' dim_abcdef=',dim_abcdef
      allocate(start_abcdef(dim_abcdef), stat=aerr)
      if(aerr /= 0) call memerror("threebodysetup_cJ 3")
      start_abcdef=0
      ii=0
      do isp1n_a=1,isp1ntot
         n_a=isp1n_cJ(1,isp1n_a)
         l_a=isp1n_cJ(2,isp1n_a)
         j2_a=isp1n_cJ(3,isp1n_a)
         ntot_a=2*n_a+l_a
         if (ntot_a>nhom) cycle
         do isp1n_b=1,isp1n_a
            n_b=isp1n_cJ(1,isp1n_b)
            l_b=isp1n_cJ(2,isp1n_b)
            j2_b=isp1n_cJ(3,isp1n_b)
            ntot_b=2*n_b+l_b
            if (ntot_a+ntot_b>nhom2sp) cycle
            do isp1n_c=1,isp1n_b
               n_c=isp1n_cJ(1,isp1n_c)
               l_c=isp1n_cJ(2,isp1n_c)
               j2_c=isp1n_cJ(3,isp1n_c)
               ntot_c=2*n_c+l_c
               if (ntot_a+ntot_c>nhom2sp) cycle
               if (ntot_b+ntot_c>nhom2sp) cycle
               if (ntot_a+ntot_b+ntot_c>nhom) cycle
               i_abc=index_abc(isp1n_a,isp1n_b,isp1n_c)
!               do isp1n_d=1,isp1ntot
               do isp1n_d=1,isp1n_a
                  n_d=isp1n_cJ(1,isp1n_d)
                  l_d=isp1n_cJ(2,isp1n_d)
                  j2_d=isp1n_cJ(3,isp1n_d)
                  ntot_d=2*n_d+l_d
                  if (ntot_d>nhom) cycle
                  if (isp1n_d==isp1n_a) then
                     isp1n_e_max=isp1n_b
                  else
                     isp1n_e_max=isp1n_d
                  endif
!                  do isp1n_e=1,isp1n_d
                  do isp1n_e=1,isp1n_e_max
                     n_e=isp1n_cJ(1,isp1n_e)
                     l_e=isp1n_cJ(2,isp1n_e)
                     j2_e=isp1n_cJ(3,isp1n_e)
                     ntot_e=2*n_e+l_e
                     if (ntot_d+ntot_e>nhom2sp) cycle
                     if (isp1n_d==isp1n_a.and.isp1n_e==isp1n_b) then
                        isp1n_f_max=isp1n_c
                     else
                        isp1n_f_max=isp1n_e
                     endif
                     do isp1n_f=1,isp1n_f_max
                        n_f=isp1n_cJ(1,isp1n_f)
                        l_f=isp1n_cJ(2,isp1n_f)
                        j2_f=isp1n_cJ(3,isp1n_f)
                        ntot_f=2*n_f+l_f
                        if (ntot_d+ntot_f>nhom2sp) cycle
                        if (ntot_e+ntot_f>nhom2sp) cycle
                        if (ntot_d+ntot_e+ntot_f>nhom) cycle
                        if (mod(l_a+l_b+l_c+l_d+l_e+l_f,2)==1) cycle
                        i_def=index_abc(isp1n_d,isp1n_e,isp1n_f)
                        i_abcdef=index_abcdef(i_abc,i_def)
                        start_abcdef(i_abcdef)=ii+1
                        do j_ab=abs(j2_a-j2_b)/2,(j2_a+j2_b)/2
                           do j_de=abs(j2_d-j2_e)/2,(j2_d+j2_e)/2
                              do J_3=MAX(abs(2*j_ab-j2_c),  &
                                  abs(2*j_de-j2_f)), &
                                  MIN(2*j_ab+j2_c,2*j_de+j2_f),2
                                 do t_ab=0,1

                                    do t_de=0,1


                                       do T_3=MAX(abs(2*t_ab-1),   &
                                           abs(2*t_de-1)),         &
                                           MIN(2*t_ab+1,2*t_de+1),2
                                          ii=ii+1
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      tot_num_of_3bmatel=ii
      if (iproc==0) then
         print *,' tot_num_of_3bmatel=',tot_num_of_3bmatel
!         write(2,1000) tot_num_of_3bmatel
! 1000    format(/,' Number of 3N matrix elements in threebodysetup_cJ:',
!     $        i8)
      endif
      allocate(v3b_cJ(tot_num_of_3bmatel), stat=aerr)
      if(aerr /= 0) call memerror("threebodysetup_cJ 2")
      if (iproc==0) then
         open(27,file=trim(v3intfile),status='old',form='unformatted',  &
             action='read')
         ii=0
         do isp1n_a=1,isp1ntot
            n_a=isp1n_cJ(1,isp1n_a)
            l_a=isp1n_cJ(2,isp1n_a)
            j2_a=isp1n_cJ(3,isp1n_a)
            ntot_a=2*n_a+l_a
            if (ntot_a>nhom) cycle
            do isp1n_b=1,isp1n_a
               n_b=isp1n_cJ(1,isp1n_b)
               l_b=isp1n_cJ(2,isp1n_b)
               j2_b=isp1n_cJ(3,isp1n_b)
               ntot_b=2*n_b+l_b
               if (ntot_a+ntot_b>nhom2sp) cycle
               do isp1n_c=1,isp1n_b
                  n_c=isp1n_cJ(1,isp1n_c)
                  l_c=isp1n_cJ(2,isp1n_c)
                  j2_c=isp1n_cJ(3,isp1n_c)
                  ntot_c=2*n_c+l_c
                  if (ntot_a+ntot_c>nhom2sp) cycle
                  if (ntot_b+ntot_c>nhom2sp) cycle
                  if (ntot_a+ntot_b+ntot_c>nhom) cycle
                  i_abc=index_abc(isp1n_a,isp1n_b,isp1n_c)
                  pi_i=(-1)**(l_a+l_b+l_c)
                  do isp1n_d=1,isp1n_a
                     n_d=isp1n_cJ(1,isp1n_d)
                     l_d=isp1n_cJ(2,isp1n_d)
                     j2_d=isp1n_cJ(3,isp1n_d)
                     ntot_d=2*n_d+l_d
                     if (ntot_d>nhom) cycle
                     if (isp1n_d==isp1n_a) then
                        isp1n_e_max=isp1n_b
                     else
                        isp1n_e_max=isp1n_d
                     endif
                     do isp1n_e=1,isp1n_e_max
                        n_e=isp1n_cJ(1,isp1n_e)
                        l_e=isp1n_cJ(2,isp1n_e)
                        j2_e=isp1n_cJ(3,isp1n_e)
                        ntot_e=2*n_e+l_e
                        if (ntot_d+ntot_e>nhom2sp) cycle
                        if (isp1n_d==isp1n_a.and.isp1n_e==isp1n_b) then
                           isp1n_f_max=isp1n_c
                        else
                           isp1n_f_max=isp1n_e
                        endif
                        do isp1n_f=1,isp1n_f_max
                           n_f=isp1n_cJ(1,isp1n_f)
                           l_f=isp1n_cJ(2,isp1n_f)
                           j2_f=isp1n_cJ(3,isp1n_f)
                           ntot_f=2*n_f+l_f
                           if (ntot_d+ntot_f>nhom2sp) cycle
                           if (ntot_e+ntot_f>nhom2sp) cycle
                           if (ntot_d+ntot_e+ntot_f>nhom) cycle
                           pi_f=(-1)**(l_d+l_e+l_f)
                           if (pi_i/=pi_f) cycle
                           i_def=index_abc(isp1n_d,isp1n_e,isp1n_f)
!                           if (i_def<i_abc) cycle
!     i_abcdef=index_abcdef(i_abc,i_def)
!                        start_abcdef(i_abcdef)=ii+1
                           do j_ab=abs(j2_a-j2_b)/2,(j2_a+j2_b)/2
                              do j_de=abs(j2_d-j2_e)/2,(j2_d+j2_e)/2
                                 do J_3=MAX(abs(2*j_ab-j2_c),  &
                                     abs(2*j_de-j2_f)),       &
                                     MIN(2*j_ab+j2_c,2*j_de+j2_f),2
                                    do t_ab=0,1

                                       do t_de=0,1
                                          
                                          do T_3=MAX(abs(2*t_ab-1), & 
                                              abs(2*t_de-1)),      &
                                              MIN(2*t_ab+1,2*t_de+1),2
                                             ii=ii+1
                                             read(27) v3b_cJ(ii)
!     print *, isp1n_a,isp1n_b,
!     $                                         isp1n_c,isp1n_d,isp1n_e,
!     $                                         isp1n_f,j_ab,t_ab,
!     $                                         j_de,t_de,J_3,T_3,
!     $                                         v3b_cJ(ii)
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      endif

      if (iproc==0) print *,' v3b_cJ read, size ',tot_num_of_3bmatel
!------------------ check
      rewind(27)
      iii = 0
      do ii = 1,tot_num_of_3bmatel*4
         read(27,end=321)vx
         iii = iii+1
      end do
      if(iproc==0)print*,' hey ran out of room '
321   continue
      if(iii /= tot_num_of_3bmatel)then
        if(iproc==0)print*,' actual size = ',iii
        stop
      end if
      call BMPI_Barrier(icomm,ierr)
      ibuf=3000000
      if (tot_num_of_3bmatel<=ibuf) then
         i=tot_num_of_3bmatel
         call BMPI_BCAST(v3b_cJ(1),i,0,icomm,ierr) 
      else
         ih3oibuf=tot_num_of_3bmatel/ibuf
         iix=1
         do i=1,ih3oibuf
            call BMPI_Bcast(v3b_cJ(iix),ibuf,0,icomm,ierr) 
            if (iproc==0) print *,' Broadcasted iix=',iix
            iix=iix+ibuf
         end do
         i=mod(tot_num_of_3bmatel,int(ibuf,kind(4)))
         if (i/=0) then
            call BMPI_BCASt(v3b_cJ(iix),i,0,icomm,ierr)
            if (iproc==0) print *,' Broadcasted iix=',iix
         endif
      endif
      if (iproc==0) print *,' H3 broadcasted'
      call cg_init(N1_max,nhom2sp,nhom)
      if (iproc==0) print *,' cg_init called'
      return
      end
!=========================================================================================
      subroutine sp1nbas_cJ(nhom,iunitout)
!      use parameters
!      use spsdata
      use spbasis_3BJC
      use v3b
      use sporbit
      use spstate
      use nodeinfo
	  use bmpi_mod
      implicit none

      integer,intent(IN) :: nhom,iunitout
      integer :: n,l,j2
      integer :: ii,ntot
      integer :: mxsps2  ! = 2 x # of single-particles, rounded up to multiple of 64 -- PN
      integer :: iorb
      integer :: aerr

      mxsps2 = 2*nsps(2)
      mxsps2 = (mxsps2/64)*64
      if(mxsps2 < 2*nsps(2))mxsps2 = mxsps2+64

      ii=0
      do ntot=0,nhom
         do l=mod(ntot,2),ntot,2
            n=(ntot-l)/2
            do j2=iabs(2*l-1),2*l+1,2
               ii=ii+1
            end do
         end do
      end do
      isp1ntot=ii
      if (allocated(isp1n_cJ)) deallocate(isp1n_cJ)
      allocate(isp1n_cJ(3,isp1ntot), stat=aerr)
      if(aerr /= 0) call memerror("sp1nbas_cJ")
      if (iproc==0) then
         write(iunitout,1000) isp1ntot
!      write(6,1000) isp1ntot
 1000    format(' Number of nlj states:',i5)
      endif
      ii=0

!....... ALLOCATE MAPS TO STANDARD BIGSTICK S.P. STATES
! !  if(allocated(map2psps))deallocate(map2psps,cmap2psps) 
!   if(allocated(map2nsps))deallocate(map2nsps,cmap2nsps) 
!   allocate( map2psps( isp1ntot ),map2nsps( isp1ntot+1:2*isp1ntot ) , stat=aerr)
!     if(aerr /= 0) call memerror("sp1nbas_cJ")
!   allocate( cmap2psps( isp1ntot ),cmap2nsps( isp1ntot+1:2*isp1ntot ) , stat=aerr)
!     if(aerr /= 0) call memerror("sp1nbas_cJ")
!   map2psps(:) = -1
!   map2nsps(:) = -1
!   cmap2psps(:)= -1
!   cmap2nsps(:)= -1
    if(allocated(map2orb))deallocate(map2orb, stat=aerr)
    if(aerr /= 0) call memerror("sp1nbas_cJ 10")
    if(allocated(mapfromorb))deallocate(mapfromorb)
    allocate( map2orb(isp1ntot), stat=aerr)
    if(aerr /= 0) call memerror("sp1nbas_cJ 11")
    allocate( mapfromorb(numorb(2) ), stat=aerr)
    if(aerr /= 0) call memerror("sp1nbas_cJ 12")
    map2orb(:) = -1
    mapfromorb(:) = -1
    do ntot=0,nhom
         do l=mod(ntot,2),ntot,2
            n=(ntot-l)/2
            do j2=iabs(2*l-1),2*l+1,2
               ii=ii+1
               isp1n_cJ(1,ii)=n
               isp1n_cJ(2,ii)=l
               isp1n_cJ(3,ii)=j2
!......... MAP THIS STATE TO NATIVE BIGSTICK ORBIT...................
               do iorb = 1,numorb(2)
                   if( n == orbqn(2,iorb)%nr .and. l == orbqn(2,iorb)%l .and. j2 == orbqn(2,iorb)%j)then
                       map2orb(ii) = iorb
                       if(iorb> numorb(2) )then
                             print*,' problem with iorb ',iorb
                             stop
                       endif
                       mapfromorb(iorb) = ii
                       exit
                   endif
               end do ! iorb

!               write(2,1100) ii,n,l,j2
!               write(6,1100) ii,n,l,j2
! 1100          format(' #',i5,'  n=',i3,'  l=',i3,
!     +              '  j=',i3,'/2')
            end do
         end do
      end do
!      allocate(nlj_orb(mxsps2), stat=aerr)
!      if(aerr /= 0) call memerror(sp1nbas_cJ")
!      nlj_orb=0
!      do ii=1,mxsps2
!         do l=1,isp1ntot
!            if (n_sp(ii)==isp1n_cJ(1,l).and.l_sp(ii)==isp1n_cJ(2,l) &
!                .and.j2_sp(ii)==isp1n_cJ(3,l)) then
!               nlj_orb(ii)=l
!               exit
!            endif
!         end do
!      end do
      return
      end


!=========================================================================================
      subroutine cg_init(N1max,N12max,N123max)
      use v3b
      use spbasis_3BJC
      implicit none
      integer,intent(IN) :: N1max,N12max,N123max
      integer :: mt1,mt2,t12,mt12,mt3,T3
      integer :: j1,m1,j2,m2,j12,m12,j3,m3,j123,j1max,j12max,j123max
!      real(kind(0.d0)) :: clebd
      real             :: cleb
      integer :: aerr

      cgt12=0.d0
      cgt123=0.d0
      do mt1=-1,1,2
         do mt2=-1,1,2
            do t12=abs(mt1+mt2),2,2
               cgt12(t12/2,(mt1+1)/2,(mt2+1)/2)= &
                   cleb(1,mt1,1,mt2,t12,mt1+mt2)
            end do
         end do
      end do
      do t12=0,2,2
         do mt12=-t12,t12,2
            do mt3=-1,1,2
               do T3=MAX(abs(t12-1),abs(mt12+mt3)),t12+1,2
                  cgt123(t12/2,T3/2,mt12/2,(mt3+1)/2)= &
                      cleb(t12,mt12,1,mt3,T3,mt12+mt3)
               end do
            end do
         end do
      end do

      j1max=2*N1max+1
      j12max=2*(N12max+1)
      j123max=2*N123max+3

      allocate(cgj12(0:j12max/2,0:j1max/2,0:j1max,0:j1max/2,0:j1max), stat=aerr)
      if(aerr /= 0) call memerror("cg_init 1")
      allocate(cgj123(0:j123max/2,0:j12max/2,-j12max/2:j12max/2, 0:j1max/2,0:j1max), stat=aerr)
      if(aerr /= 0) call memerror("cg_init 2")
      cgj12=0.d0
      cgj123= 0.d0
      do j1=1,j1max,2
         do m1=-j1,j1,2
            do j2=1,j1max,2
               do m2=-j2,j2,2
                  do j12=MAX(abs(j1-j2),abs(m1+m2)),MIN(j1+j2,j12max),2
                     cgj12(j12/2,j1/2,(j1+m1)/2,j2/2,(j2+m2)/2)= &
                         cleb(j1,m1,j2,m2,j12,m1+m2)
                  end do
               end do
            end do
         end do
      end do

      do j12=0,j12max,2
         do m12=-j12,j12,2
            do j3=1,j1max,2
               do m3=-j3,j3,2
                  do j123=MAX(abs(j12-j3),abs(m12+m3)), &
                      MIN(j12+j3,j123max),2
                     cgj123(j123/2,j12/2,m12/2,j3/2,(j3+m3)/2)= &
                         cleb(j12,m12,j3,m3,j123,m12+m3)
                  end do
               end do
            end do
         end do
      end do
      return
      end

!======================================================================
!
!compute min HO for n particle
! a bit klugdy but should be sufficient
!

      function nmin_ho(n)

      implicit none
      integer n
      integer nmin_ho
      integer i,j,k
      integer w(40)

      nmin_ho = 0
      if(n < 1)return
      if(n > 40)then
         print*,' Problem with Nmin_HO -- must raise size ',n
         stop
      end if
! .... set up array
      j = 1
      do i = 0,3
        do k = 1,(i+1)*(i+2)/2
           w(j) = i
           j = j+1
           w(j) = i
           j = j + 1
        end do
      end do
      do i = 1,n
        Nmin_HO = Nmin_HO + w(i)
      end do
     
      return
      end function nmin_ho

!==========================================================================

