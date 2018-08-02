C**********************************************************************
C   FINDEX      function to find index of s.p. state if QN are known  *
C   SP_INDEX    subroutine to setup array of s.p. states with QN      *
C   RLI         matrix elements of r**L with H.O. wave functions      *
C               input is orbit label setup in SP_INDEX                *
C   RL          matrix elements of r**L with H.O. wave functions      *
C   TR          checks triangle relations for adding angular momentum *
C   ITR         same as TR with integer input                         *
C   XFAC        factorial                                             *
C   DFAC        double factorial                                      *
C   BINOM       binomial coefficients                                 *
C   PAR         (-1.)**X (where x is real)                            *
C   TH_FACINIT  log of factorial- stored for use in CLEB routines     *
C   XGAMMA      gamma function of X                                   *
C   HATI        sqrt(2*I+1)                                           *
C   HAT         sqrt(2*X+1)   entry into HATI                         *
C   CLEBR       Clebsch-Gordon (CG) coefficients with real input      *
C   CLEBI       CG coefficients (integer input - entry inot CLEBR)    *
C   CLEB        CG coefficients (2*integer input - entry inot CLEBR)  *
C   TJ2I        three-J with integer twice normal value (uses CLEBR)  *
C   TJI         three-J with integer input (entry into TJ2I           *
C   TJ          three-J with real input (entry into TJ2I)             *
C   RACAH       Racah coefficient (2*integer input)                   *
C   RACAHI      Racah coefficient (integer input - entry into RACAH)  *
C   RACAHR      Racah coefficient (real input - entry into RACAH)     *
C   SJ          six-J coefficient (real input) uses subroutine RACAH  *
C   SJ2I        six-J with integer input of twice normal value        *
C   COEF9J      nine-J coefficient                                    *
C   AFAC        AFAC(N,L)=(-1.)**N/SQRT((2*N)!!*(2*N+2*L+1)!!)        * 
C   FFAC        factor used in BMOSH                                  *
C   BMOSH       Moshinsky brackets                                    *
C   DELTAI      Kronecker Delta with integer arguments                *
C   DELTAR      Kronecker Delta with real arguments                   *
C**********************************************************************
C      QN   stands for quantum numbers. s.p. for single-particle      *
C**********************************************************************
C     Function to find index for single-particle state                *
C**********************************************************************
      FUNCTION FINDEX(XN,XL,XJ)
      FINDEX = 0
      IF(XN.EQ.0.) RETURN
      IF(TR(XL,XJ,0.5).EQ.-1.) RETURN
      N = XN
      L = XL
      J2 = 2.*XJ
      FINDEX = ((2*N+L)*(2*N+L+3)-(J2-1)+2)/2
      RETURN
      END
C************************************************************************
C                      Subroutine to set up quantum numbers for all     *
C                      shell-model single particle states               *
C                      Note that n starts with 0 and not 1              *
C                      Also, note change in name to avoid conflict      *
C                      with intrinsic function index (same for common   *
C                      block index1                                     *
C************************************************************************
      SUBROUTINE SP_INDEX 
      CHARACTER*1 LABEL
      COMMON/CLABEL/LABEL(20)
      COMMON/INDEX1/XL(3,200)
      COMMON/INDEX2/IORB(50,50),INUM(50)
      I = 0
      DO 100 I1 = 1,1000
      INUM(I1)=0
      J=1
      N = I1-1
      LN = 0
      L = N-2*LN
200   IF(L.LT.0) GO TO 100
      I = I + 1
      IF(I.GT.200) GO TO 101
      XJ2 = 2*L+1
      IORB(I1,J)=I
      J=J+1
      INUM(I1)=INUM(I1)+1
C      XL(1,I) = LN+1
      XL(1,I) = LN
      XL(2,I) = L
      XL(3,I) = XJ2/2.
      IF(L.EQ.0) GO TO 201
      I = I + 1
      IF(I.GT.200) GO TO 101
      IORB(I1,J)=I
      J=J+1
      INUM(I1)=INUM(I1)+1
      XJ2 = 2*L-1
C      XL(1,I) = LN+1
      XL(1,I) = LN
      XL(2,I) = L
      XL(3,I) = XJ2/2.
201   L = L-2
      LN = LN + 1
      GO TO 200
100   CONTINUE
101   CONTINUE
      LABEL(1) = 's'
      LABEL(2) = 'p'
      LABEL(3) = 'd'
      LABEL(4) = 'f'
      LABEL(5) = 'g'
      LABEL(6) = 'h'
      LABEL(7) = 'i'
      LABEL(8) = 'j'
      LABEL(9) = 'k'
      LABEL(10) = 'l'
      LABEL(11) = 'm'
      LABEL(12) = 'n'
      LABEL(13) = 'o'
      LABEL(14) = 'p'
      LABEL(15) = 'q'
      LABEL(16) = 'r'
      LABEL(17) = 's'
      LABEL(18) = 't'
      LABEL(19) = 'u'
      LABEL(20) = 'v'
      END
C**********************************************************************
C                        Matrix elements <n1,l1|r**L| n2,l2>          *
C                        with H.O. wave functions                     *
C**********************************************************************
      FUNCTION RLI(i1,i2,L)
      COMMON /INDEX1/ XL(3,200)
      PI = 2.0*asin(1.0)
      XN1 = XL(1,I1)
      XL1 = XL(2,I1)
      XN2 = XL(1,I2)
      XL2 = XL(2,I2)
      XX = L
      X2 = 2.0*XL1+2.0*XN1+1.0
      X3 = 2.0*XL2+2.0*XN2+1.0
      IM = XN1+1.0
      JM = XN2 + 1.0
      S = 0.0
      DO 13 I = 1,IM
      DO 13 J = 1,JM
      XI = I-1
      XJ = J-1
      IF(PAR(XL1+XL2+XX)) 20,20,21
20    X6P = XI+XJ+(XL1+XL2+XX+1.)/2.
21    CONTINUE
      X6 = XL1+XL2+2.0*XI+2.0*XJ+XX+1.0
      X7 = 2.0*XL1+2.0*XI+1.0
      X8 = 2.0*XL2+2.0*XJ+1.0
      T = PAR(XI+XJ)/(XFAC(XI)*XFAC(XJ))
      IF(PAR(XL1+XL2+XX)) 30,30,31
30    T = T*XFAC(X6P)
      XMM = (X6+1.)/2.
      GO TO 32
31    T = T*DFAC(X6)
32    CONTINUE
      T = T/(DFAC(X7)*DFAC(X8)*XFAC(XN1-XI)*XFAC(XN2-XJ))
      IF(PAR(XL1+XL2+XX)) 40,40,13
40    T = T*(2.**XMM)/SQRT(PI)
13    S = S + T
      S = S*SQRT(DFAC(X2)*DFAC(X3))
      S = S*SQRT(XFAC(XN1)*XFAC(XN2))
      S = S/SQRT(2.0**(XN1+XN2+XX))
      RLI = S
      RETURN
      END
      function RL(xn1,xl1,xn2,xl2,L)
      PI = 2.0*asin(1.0)
      XX = L
      X2 = 2.0*XL1+2.0*XN1+1.0
      X3 = 2.0*XL2+2.0*XN2+1.0
      IM = XN1+1.0
      JM = XN2 + 1.0
      S = 0.0
      DO 13 I = 1,IM
      DO 13 J = 1,JM
      XI = I-1
      XJ = J-1
      IF(PAR(XL1+XL2+XX)) 20,20,21
20    X6P = XI+XJ+(XL1+XL2+XX+1.)/2.
21    CONTINUE
      X6 = XL1+XL2+2.0*XI+2.0*XJ+XX+1.0
      X7 = 2.0*XL1+2.0*XI+1.0
      X8 = 2.0*XL2+2.0*XJ+1.0
      T = PAR(XI+XJ)/(XFAC(XI)*XFAC(XJ))
      IF(PAR(XL1+XL2+XX)) 30,30,31
30    T = T*XFAC(X6P)
      XMM = (X6+1.)/2.
      GO TO 32
31    T = T*DFAC(X6)
32    CONTINUE
      T = T/(DFAC(X7)*DFAC(X8)*XFAC(XN1-XI)*XFAC(XN2-XJ))
      IF(PAR(XL1+XL2+XX)) 40,40,13
40    T = T*(2.**XMM)/SQRT(PI)
13    S = S + T
      S = S*SQRT(DFAC(X2)*DFAC(X3))
      S = S*SQRT(XFAC(XN1)*XFAC(XN2))
      S = S/SQRT(2.0**(XN1+XN2+XX))
      RL = S
      RETURN
      END
C**********************************************************************
C                                  Kroneker delta                     *
C**********************************************************************
      FUNCTION DELTAK(ia,ib)
      deltak=0.0
      if(ia.eq.ib)deltak=1.0
      RETURN
      END
C**********************************************************************
C                          subroutine to check triangle relation for  *
C                          angular momentum coupling returns -1 if    *
C                          not possible  (used in CLEB routines)      *
C**********************************************************************
      FUNCTION TR(A1,A2,A3)
      IF(A1.LT.0.) GO TO 2
      IF(A2.LT.0.) GO TO 2
      IF(A3.LT.0.) GO TO 2
      AMIN = ABS(A1-A2)
      AMAX = A1+A2
      XA = AMIN
10    CONTINUE
      IF(A3.EQ.XA) GO TO 20
      XA = XA+1.
      IF(XA.GT.AMAX) GO TO 2
      GO TO 10
20    TR=1.
      RETURN
2     TR=-1.
      RETURN
      END

C**********************************************************************
C                        subroutine to check triangle relation for    *
C                        angular momentum coupling --- integer input  *
C                        returns -1 if not possible                   *
C**********************************************************************
      FUNCTION ITR(IA1,IA2,IA3)
      IF(IA1.LT.0.) GO TO 2
      IF(IA2.LT.0.) GO TO 2
      IF(IA3.LT.0.) GO TO 2
      IAMIN = IABS(IA1-IA2)
      IAMAX = IA1+IA2
      IXA = IAMIN
10    CONTINUE
      IF(IA3.EQ.IXA) GO TO 20
      IXA = IXA+1
      IF(IXA.GT.IAMAX) GO TO 2
      GO TO 10
20    ITR=1
      RETURN
2     ITR=-1
      RETURN
      END

C**********************************************************************
C                                             X  factorial            *
C**********************************************************************
      FUNCTION XFAC(X)
      IF(X) 3,4,4
3     XFAC = 0.0
      RETURN
4     S = 1.0
      DO 1 I = 1,1000
      Y = I
      IF(X-Y) 2,1,1
1     S = S*Y
2     XFAC = S
      RETURN
      END

C**********************************************************************
C                                      Double factorial               *
C**********************************************************************
        FUNCTION DFAC(A)
        IF(A.LT.0.)THEN
          DFAC=0.
          RETURN
        END IF
        IF(ABS(A).LT.0.000001)then   !  IT IS ZERO
           DFAC=1.0
           RETURN
        END IF
8       N=A
        if(par(a).eq.1.)write(6,7777)
7777    FORMAT(1X,'ERROR IN DFAC: A MUST BE ODD')
        IF(PAR(A).EQ.1.) STOP
        DFAC=1.
        DO 9 I=1,N,2
        AI=I
9       DFAC=AI*DFAC
        RETURN
        END
C**********************************************************************
C                                         binomial coefficient        *
C**********************************************************************
        FUNCTION BINOM(A,B)
        IF(A.EQ.B)GO TO 12
        IF(B.EQ.0.)GO TO 12
        IF(A-B.GT.B)GO TO 9
        M=A-B
        GO TO 11
9       M=B
11      BINOM=1.
        X=A-FLOAT(M)
        DO 10 I=1,M
        AI=I
10      BINOM=BINOM*(X+AI)/AI
        RETURN
12      BINOM=1.
        RETURN
        END
C**********************************************************************
C                                      parity, i.e. (-1.)**a          *
C**********************************************************************
      FUNCTION PAR(A)
      I = A
      I = I/2
      P = A/2.0
      Q = I
      IF(P-Q) 1,2,1
1     PAR = -1.0
      RETURN
2     PAR = 1.0
      RETURN
      END
	BLOCK DATA INITIAL
	real*8 faclog
	LOGICAL FIRST
c	PARAMETER (LFACTC=200)
	PARAMETER (LFACTC=500)
	DATA FIRST/.TRUE./,LFACT/LFACTC/
	COMMON / LOGFAC / FIRST,LFACT,FACLOG(LFACTC)
	END 
C**********************************************************************
C                               log of factorials - used in cleb      *
C**********************************************************************
      SUBROUTINE TH_FACINIT
C ...  SET UP LOG OF FACTORIALS
      PARAMETER (LFACTC=500)
c      PARAMETER (LFACTC=200)
      LOGICAL FIRST
      COMMON / LOGFAC / FIRST,LFACT,FACLOG(LFACTC)
      REAL*8 FACLOG,FN
      FIRST=.FALSE.
      FACLOG(1)=0.0
      FACLOG(2)=0.0
      FN=1.0
      DO 10 I=3,LFACTC
      FN=FN+1.0
      FACLOG(I)=FACLOG(I-1)+LOG(FN)
   10 CONTINUE
      RETURN
      END
C**********************************************************************
C             GAMMA OF X, WHERE X IS EITHER INTEGER OR HALF INTEGER   *
C             FUNCTION RETURNS ZERO IF THE ARGUEMENT IS NOT INTEGER   *
C             OR HALF INTEGER                                         *
C**********************************************************************
      FUNCTION XGAMMA(X)
      IX=INT(2.0*X)
      XGAMMA=0.0
      IF(FLOAT(IX).NE.(2.*X))RETURN
      IF(X.EQ.0.5)GOTO 11
      IF(IAND(IX,1).EQ.1)GOTO 10
      XGAMMA=XFAC(X-1.)
      RETURN
   10 XGAMMA=DFAC(2.0*(X-0.5)-1.0)/(2.0**(X-0.5))*SQRT(3.14159265)
      RETURN
   11 XGAMMA=SQRT(3.14159265)
      RETURN
      END
C**********************************************************************
C        RETURNS THE VALUE HI=SQRT(2*I+1), WHERE I IS AN INTEGER      *
C**********************************************************************
      FUNCTION HAT2I(I)
      x=i
      hat2i=sqrt(x+1.)
      return
      end
      FUNCTION HATI(I)
      X=FLOAT(I)
      HATI=SQRT(2.0*X+1.)
      return
      end      
      function HAT(x)
      HAT=SQRT(2.0*X+1.)
      RETURN
      END
C                               Subroutines CLEBR, CLEBI, CLEB        *
C**********************************************************************
      REAL FUNCTION CLEBR(A,B,C,D,E,F)
C
C      ARGUMENTS ARE REAL AND OF TRUE VALUE; J1,M1,J2,M2,J3,M3
C
CWEO      COMMON / LOGFAC / FIRST,LFACT,FACLOG(1)
      REAL*8 FACLOG
      PARAMETER (LFACTC=500)
      LOGICAL FIRST
      COMMON / LOGFAC / FIRST,LFACT,FACLOG(LFACTC)
CCCCCC      IA=2*J1,ID=2*M1 ETC.(J1 IS OF TRUE VALUE)
      IA=NINT(2.*A)
      IB=NINT(2.*C)
      IC=NINT(2.*E)
      ID=NINT(2.*B)
      IE=NINT(2.*D)
      IF=NINT(2.*F)
      GOTO 7000
C ...............  CLEBI  ...........................
C
C      ARGUMENTS ARE INTEGER AND OF TRUE VALUE
C
      ENTRY CLEBI(LL1,LM1,LL2,LM2,LL3,LM3)
      IA=2*LL1
      IB=2*LL2
      IC=2*LL3
      ID=2*LM1
      IE=2*LM2
      IF=2*LM3
      GOTO 7000
C ..............  CLEB  ..............................
C
C      ARGUMENTS ARE INTEGER AND REPRESENT TWICE THE REAL VALUE
C
      ENTRY CLEB(I2J1,I2M1,I2J2,I2M2,I2J3,I2M3)
      IA=I2J1
      IB=I2J2
      IC=I2J3
      ID=I2M1
      IE=I2M2
      IF=I2M3
 7000 IF (FIRST) CALL TH_FACINIT
      RAC=0.0
      IF(ID+IE-IF) 1000,105,1000
  105 K1=IA+IB+IC
      IF((-1)**K1) 1000,110,110
  110 K1=IA+IB-IC
      K2=IC+IA-IB
      K3=IB+IC-IA
      K4=IA-IABS (IB-IC)
      K5=IB-IABS (IC-IA)
      K6=IC-IABS (IA-IB)
      K7= MIN0 (K1,K2,K3,K4,K5,K6)
      IF(K7) 1000,120,120
  120 IF((-1)**(IA+ID)) 1000,1000,130
  130 IF((-1)**(IB+IE)) 1000,1000,140
  140 IF((-1)**(IC+IF)) 1000,1000,150
  150 IF(IA-IABS (ID)) 1000,152,152
  152 IF(IB-IABS (IE)) 1000,154,154
  154 IF(IC-IABS (IF)) 1000,160,160
  160 SIGNFC=1.0
      IAM=IA
      IBM=IB
      ICM=IC
      IDM=ID
      IEM=IE
      IFM=IF
      IF(IA-IB) 210,220,220
  210 IF(IA-IC) 215,225,225
  215 IT=IA
      IA=IB
      IB=IT
      IT=ID
      ID=IE
      IE=IT
      SIGNFC=(-1.0)**((IA+IB-IC)/2)
      GO TO 235
  220 IF(IC-IB) 225,235,235
  225 IT=IC
      IC=IB
      IB=IT
      IT=IF
      IF=-IE
      IE=-IT
      FIBM=IBM+1
      FICM=ICM+1
      SIGNFC=(-1.)**((IAM-IDM)/2)*SQRT (FICM/FIBM)
  235 IF(IB) 237,236,237
  236 RAC=SIGNFC
      GO TO 900
  237 IF(IE) 250,250,240
  240 SIGNFC=SIGNFC*((-1.0)**((IA+IB-IC)/2))
      ID=-ID
      IE=-IE
      IF=-IF
  250 FC2=IC+1
      IABCP=(IA+IB+IC)/2+1
      IABC=IABCP-IC
      ICAB=IABCP-IB
      IBCA=IABCP-IA
      IAPD=(IA+ID)/2+1
      IAMD=IAPD-ID
      IBPE=(IB+IE)/2+1
      IBME=IBPE-IE
      ICPF=(IC+IF)/2+1
      ICMF=ICPF-IF
      SQFCLG=0.5*(ALOG(FC2)-FACLOG(IABCP+1)
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))
      NZMIC2=(IB-IC-ID)/2
      NZMIC3=(IA-IC+IE)/2
      NZMI= MAX0 (0,NZMIC2,NZMIC3)+1
      NZMX= MIN0 (IABC,IAMD,IBPE)
      IF(NZMI-NZMX) 310,310,900
  310 SS=0.0
      S1=(-1.0)**(NZMI-1)
      DO 400 NZ=NZMI,NZMX
      NZM1=NZ-1
      NZT1=IABC-NZM1
      NZT2=IAMD-NZM1
      NZT3=IBPE-NZM1
      NZT4=NZ-NZMIC2
      NZT5=NZ-NZMIC3
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)
      SSTERM=S1*EXP (TERMLG)
      SS=SS+SSTERM
  400 S1=-S1
      RAC=SIGNFC*SS
  900 IA=IAM
      IB=IBM
      IC=ICM
      ID=IDM
      IE=IEM
      IF=IFM
 1000 CLEB=RAC
      RETURN
      END
C**********************************************************************
C                             Three-J symbol                          *
C                             given in terms of a clebsh              *
C**********************************************************************
      FUNCTION TJ2I(J21,J22,J23,M21,M22,M23)
      A1=FLOAT(J21)/2.
      B1=FLOAT(M21)/2.
      A2=FLOAT(J22)/2.
      B2=FLOAT(M22)/2.
      A3=FLOAT(J23)/2.
      B3=FLOAT(M23)/2.
      goto 9
      entry TJI(J1,J2,J3,M1,M2,M3)
      A1=FLOAT(J1)
      B1=FLOAT(M1)
      A2=FLOAT(J2)
      B2=FLOAT(M2)
      A3=FLOAT(J3)
      B3=FLOAT(M3)
      entry TJ(AA1,AA2,AA3,BB1,BB2,BB3)
      A1=AA1
      B1=BB1
      A2=AA2
      B2=BB2
      A3=AA3
      B3=BB3
 9    TJ = PAR(A1-A2-B3)*CLEBR(A1,B1,A2,B2,A3,-B3)/SQRT(2.*A3+1.)
      RETURN
      END
C**********************************************************************
C  CALCULATES THE PRODUCT OF 4 CLEBSCH GORDON COEFFICIENTS
C    WHERE EACH IS OF THE FORM
C
C        (K,J,0,0|L,0)
C**********************************************************************
      FUNCTION CLEBSCH(K1,J1,L1,K2,J2,L2,KK1,KK2,L,JJ1,JJ2,LL)
      CLEBSCH=CLEBI(K1,0,J1,0,L1,0)*CLEBI(K2,0,J2,0,L2,0)*
     $CLEBI(KK1,0,KK2,0,L,0)*CLEBI(JJ1,0,JJ2,0,LL,0)
      RETURN
      END
C**********************************************************************
C                           Computes Racah coefficient                *
C**********************************************************************
      FUNCTION RACAH(JAD,JBD,JCD,JDD,JED,JFD)
C
C        CALCULATES RACAH COEFFICIENTS
C
C        SOURCE : UNKNOWN
C        MODIFIED : MARCH 1982 , OLAF SCHOLTEN
C              RUN TIME OPTIMIZED FOR VAX780 MACHINE
C
C        ENTRIES : RACAH , RACAHI , RACAHR
C            RACAH  : INTEGER ARGUMENTS = 2*J
C            RACAHI : ARGUMENTS = TRUE INTEGER VALUE
C            RACAHR : ARGUMENTS = TRUE REAL VALUE
C        EXTERNAL : TH_FACINIT , GENERATES FACTORIAL TABLE
C
      DIMENSION I(16)
      REAL*8 G,S
      PARAMETER (LFACTC=500)
      LOGICAL FIRST
      COMMON / LOGFAC / FIRST,LFACT,G(LFACTC)
      EQUIVALENCE(I(1),I1),(I(2),I2),(I(3),I3),(I(4),I4),(I(5),I5),
     1 (I(6),I6),(I(7),I7),(I(8),I8),(I(9),I9),(I(10),I10),(I(11),I11),
     2 (I(12),I12),(I(13),I13),(I(14),I14),(I(15),I15),(I(16),I16)
C        MAKE USEFULL COMBINATIONS
      K=JAD+JBD-JED+2
      I1=K/2
      IF((2*I1).NE.K) GOTO 300
      K=JCD+JDD-JED+2
      I4=K/2
      IF((2*I4).NE.K) GOTO 300
      K=JAD+JCD-JFD+2
      I7=K/2
      IF((2*I7).NE.K) GOTO 300
      K=JBD+JDD-JFD+2
      I10=K/2
      IF((2*I10).NE.K) GOTO 300
      I13=I1+JED
      I14=I4+JED
      I15=I7+JFD
      I16=I10+JFD
      I2=I13-JAD
      I3=I13-JBD
      I5=I14-JCD
      I6=I14-JDD
      I8=I15-JAD
      I9=I15-JCD
      I11=I16-JBD
      I12=I16-JDD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,2,2
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    2 IL=MAX(I13,I14,I15,I16)
      IF(MIN(JAD,JBD,JCD,JDD,JED,JFD)) 300,20,1
C    ..............
      ENTRY RACAHI(JA1,JB1,JC1,JD1,JE1,JF1)
C        MAKE USEFULL COMBINATIONS
      I13=JA1+JB1+JE1+1
      I14=JC1+JD1+JE1+1
      I15=JA1+JC1+JF1+1
      I16=JB1+JD1+JF1+1
      I1=I13-JE1*2
      I2=I13-JA1*2
      I3=I13-JB1*2
      I4=I14-JE1*2
      I5=I14-JC1*2
      I6=I14-JD1*2
      I7=I15-JF1*2
      I8=I15-JA1*2
      I9=I15-JC1*2
      I10=I16-JF1*2
      I11=I16-JB1*2
      I12=I16-JD1*2
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,4,4
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    4 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA1,JB1,JC1,JD1,JE1,JF1)
      IF(LMIN)300,20,1
C     ............
      ENTRY RACAHR(A,B,C,D,E,F)
C     CONVERT ARGUMENTS TO INTEGER 
      JA=NINT(2.*A)
      JB=NINT(2.*B)
      JC=NINT(2.*C)
      JD=NINT(2.*D)
      JE=NINT(2.*E)
      JF=NINT(2.*F)
C        MAKE USEFULL COMBINATIONS
      K=JA+JB-JE+2
      I1=K/2
      IF((2*I1-K).NE.0) GOTO 300
      K=JC+JD-JE+2
      I4=K/2
      IF((2*I4-K).NE.0) GOTO 300
      K=JA+JC-JF+2
      I7=K/2
      IF((2*I7-K).NE.0) GOTO 300
      K=JB+JD-JF+2
      I10=K/2
      IF((2*I10-K).NE.0) GOTO 300
      I13=I1+JE
      I14=I4+JE
      I15=I7+JF
      I16=I10+JF
      I2=I13-JA
      I3=I13-JB
      I5=I14-JC
      I6=I14-JD
      I8=I15-JA
      I9=I15-JC
      I11=I16-JB
      I12=I16-JD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,3,3
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    3 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA,JB,JC,JD,JE,JF)
      IF(LMIN)300,20,1
C      ------------
    1 IF(FIRST) CALL TH_FACINIT
      IF(IL.GE.LFACT) STOP 'RACAH: LENGTH FACTORIAL TABLE INSUFFICIENT'
      J1=IL-I13+1 
      J2=IL-I14+1 
      J3=IL-I15+1 
      J4=IL-I16+1
      J5=I13+I4-IL 
      J6=I15+I5-IL 
      J7=I16+I6-IL
      PH=1.
      IF(2*(J5/2).EQ.J5) PH=-1.
      H=PH*EXP ((G(I1)+G(I2)+G(I3)-G(I13+1)+G(I4)+G(I5)+G(I6)-
     1G(I14+1)+G(I7)+G(I8)+G(I9)-G(I15+1)+G(I10)+G(I11)+G(I12)-G(I16+1))
     2*.5+G(IL+1)-G(J1)-G(J2)-G(J3)-G(J4)-G(J5)-G(J6)-G(J7))
      IF(N)300,110,120
C
  110 RACAH=H 
      RETURN
C
  120 S=1.
      K=N-1
      KL=IL+1
      J5=J5-1
      J6=J6-1
      J7=J7-1
      DO 130 J=1,N   ! K=N-J
      S=1.-((KL+K)*(J5-K)*(J6-K)*(J7-K))*S/((J1+K)*(J2+K)*(J3+K)*(J4+K))
      K=K-1
  130 CONTINUE  
      RACAH=H*S
      RETURN
C
C      ONE OF THE ARGUMENTS =0
   20 IAD=IL
      IBD=IL
      DO 21 J=13,16
      IF(IAD.LT.I(J)) GOTO 22
      IF(IAD.LT.IBD) IBD=IAD
      IAD=I(J)
      GOTO 21
   22 IF(IBD.GT.I(J)) IBD=I(J)
   21 CONTINUE
      J5=I13+I4-IL 
      PH=1.
      IF(2*(J5/2).EQ.J5) PH=-1.
      RACAH=PH/SQRT(FLOAT(IAD*IBD))
      RETURN
C
C      IMPOSSIBLE COMBINATION OF ARGUMENTS
  300 RACAH=0.
      RETURN
      END
C**********************************************************************
C                               sixj routine                          *
C                               given in terms of a racah coeff.      *
C**********************************************************************
      FUNCTION SJ(A1,A2,B2,B1,A3,B3)
      SJ = PAR(A1+A2+A3+B1)*RACAHR(A1,A2,A3,B1,B2,B3)
      RETURN
      END
C**********************************************************************
C                               sixj with integer input (2*I)         *
C                               given in terms of a racah coeff.      *
C**********************************************************************
      FUNCTION SJ2I(J1,J2,K2,K1,J3,K3)
      SJ2I = PAR(FLOAT(J1+J2+J3+K1)/2.)*RACAH(J1,J2,J3,K1,K2,K3)
      RETURN
      END
C**********************************************************************
C                              Subroutine  COEF9J                     *
C                              Computes nine-J coefficient            *
C**********************************************************************
      Function COEF9J(J1,J2,J3,J4,J5,J6,J7,J8,J9)
CCCC  TAMURAS NUMBERING CONVENTION IS USED HERE.
CCCCC THE ARGUMENTS OF COEF9J, WHEN NUMBERED SEQUENTIALLY 1 THROUGH 9,
CCCC    CORRESPOND TO THE ARRAY
CCCC                               1  2  5
CCCC                               3  4  6
CCCC                               7  8  9
C
      DIMENSION LT(9)
CCCCCC
CCCCCC      ALL L9 MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
CCCCCC
C
C               CHANGED FOR THEORY LIBRARY 4/8/82   HK
C
      U9=0.0
      LT(1)=J1
      LT(2)=J2
      LT(3)=J3
      LT(4)=J4
      LT(5)=J5
      LT(6)=J6
      LT(7)=J7
      LT(8)=J8
      LT(9)=J9
      LMIN=LT(1)
      IMIN=1
      DO 20 I=2,9
      IF(LT(I)-LMIN) 15,20,20
   15 LMIN=LT(I)
      IMIN=I
   20 CONTINUE
      KEX=0
      GO TO (110,110,110,110,150,150,170,170,190),IMIN
  110 MM=(IMIN-1)/2+1
      M1=MM+MM-1
      M2=M1+1
      M3=MM+4
      L1=LT(7)
      LT(7)=LT(M1)
      LT(M1)=L1
      L1=LT(8)
      LT(8)=LT(M2)
      LT(M2)=L1
      L1=LT(9)
      LT(9)=LT(M3)
      LT(M3)=L1
      IMIN=IMIN+(7-M1)
      GO TO 175
  150 KEX=1
      M1=7
      M2=8
      M3=IMIN+IMIN-9
      M4=M3+1
      GO TO 180
  170 KEX=1
  175 M1=5
      M2=6
      M3=IMIN-6
      M4=M3+2
  180 L1=LT(M1)
      L1=LT(M1)
      LT(M1)=LT(M3)
      LT(M3)=L1
      L1=LT(M2)
      LT(M2)=LT(M4)
      LT(M4)=L1
      L1=LT(9)
      LT(9)=LT(IMIN)
      LT(IMIN)=L1
  190 IF(LT(9)) 200,200,300
  200 IF(LT(5)-LT(6)) 1000,210,1000
  210 IF(LT(7)-LT(8)) 1000,220,1000
  220 RT=(LT(5)+1)*(LT(7)+1)
      K=(LT(5)+LT(7)-LT(1)-LT(4))/2
      RAC= RACAH(LT(1),LT(2),LT(3),LT(4),LT(5),LT(7))
      PH=1.
      IF (2*(K/2) .NE. K) PH=-1.
      U9=(RAC/SQRT(RT))*PH
      GO TO 370
  300 K1=IABS(LT(2)-LT(7))
      K2=IABS(LT(3)-LT(5))
      K3=IABS(LT(4)-LT(9))
      NMIN=MAX0(K1,K2,K3)
      K1=LT(2)+LT(7)
      K2=LT(3)+LT(5)
      K3=LT(4)+LT(9)
      NMAX=MIN0(K1,K2,K3)
      IF (NMIN-NMAX) 320, 320, 1000
  320 DO 350 N=NMIN,NMAX,2
      W1=N+1
      RAC= RACAH(LT(2),LT(5),LT(7),LT(3),LT(1),N)
      IF (RAC) 321, 350, 321
  321 W1=W1*RAC
      RAC= RACAH(LT(2),LT(4),LT(7),LT(9),LT(8),N)
      IF (RAC) 322, 350, 322
  322 W1=W1*RAC
      RAC= RACAH(LT(3),LT(4),LT(5),LT(9),LT(6),N)
      IF (RAC) 323, 350, 323
  323 U9=U9+W1*RAC
  350 CONTINUE
  370 IF(KEX) 400,1000,400
  400 KP=0
      DO 410 I=1,9
  410 KP=KP+LT(I)
      K=KP/2
      PH=1.
      IF (2*(K/2) .NE. K) PH=-1.
      U9=U9*PH
 1000 COEF9J=U9
      RETURN
      END
C**********************************************************************
C                   COEFR9J  -  nine-J coefficient with real input    *
C                               calls COEF9J
C**********************************************************************
      Function COEFR9J(XJ1,XJ2,XJ3,XJ4,XJ5,XJ6,XJ7,XJ8,XJ9)
CCCCC THE ARGUMENTS OF COEF9J, WHEN NUMBERED SEQUENTIALLY 1 THROUGH 9,
CCCC    CORRESPOND TO THE ARRAY
CCCC                               1  2  5
CCCC                               3  4  6
CCCC                               7  8  9
      j1=int(2.*xj1)
      j2=int(2.*xj2)
      j3=int(2.*xj3)
      j4=int(2.*xj4)
      j5=int(2.*xj5)
      j6=int(2.*xj6)
      j7=int(2.*xj7)
      j8=int(2.*xj8)
      j9=int(2.*xj9)
      coefr9j=COEF9J(J1,J2,J3,J4,J5,J6,J7,J8,J9)
      return
      end
C**********************************************************************
C                   AFAC(N,L)=(-1.)**N/SQRT((2*N)!!*(2*N+2*L+1)!!)    *
C                   USED IN BMOSH                                     *
C**********************************************************************
      FUNCTION AFAC(N,L)
      AFAC=PAR(FLOAT(N))/(SQRT(DFAC(FLOAT(2*N))*DFAC(FLOAT(2*N+2*L+1))))
      RETURN
      END
C**********************************************************************
C                                FACTOR USED BMOSH                    *
C**********************************************************************
      FUNCTION FFAC(NU,N,L,K,J)
      UP=SQRT(3.14159265)*XFAC(FLOAT(N))*XGAMMA(FLOAT(N+L)+1.5)
      DOWN=XFAC(FLOAT(NU))*(SQRT(2.)**(FLOAT(2*N+L)))*
     $XGAMMA(FLOAT(2*N+L-K-J)/2.-FLOAT(NU)+1.0)*
     $XGAMMA(FLOAT(2*N+L+3+J-K)/2.-FLOAT(NU))*
     $(XGAMMA(FLOAT(K+NU)+1.5))
      FFAC=UP/DOWN
      RETURN
      END
C*********************************************************************
C   CALCULTES MOSHINSKY BRACKETS USING THE FORMULA (A6.19)           *
C    IN THEORY OF THE NUCLEAR SHELL MODEL, R.D.LAWSON P.497          *
C*********************************************************************
      FUNCTION BMOSH(NN1,L1,NN2,L2,LB,N,L,NN,LL)
      N1=NN1-1
      N2=NN2-1
      BMOSH=0.0
      IF(LB.GT.(L+LL).OR.(LB.GT.(L1+L2)))RETURN
      IF(LB.LT.IABS(L-LL).OR.(LB.LT.IABS(L1-L2)))RETURN
      IF((2*N1+L1+2*N2+L2).NE.(2*N+L+2*NN+LL))RETURN
      IF(IAND((L+LL+L1+L2),1).EQ.1)RETURN
      IF((N1.EQ.N2).AND.(L1.EQ.L2).AND.(IAND((LB+LL),1).EQ.1))RETURN
      IF((N.EQ.NN).AND.(L.EQ.LL).AND.(IAND((LB+L1),1).EQ.1))RETURN
      ANORM=AFAC(N1,L1)*AFAC(N2,L2)/(AFAC(N,L)*AFAC(NN,LL)*4.0)
      MAX=2*N+L
      NU1MAX=MAX/2
      DO 10 I=1,NU1MAX+1
      NU1=I-1
      NU2MAX=(MAX-2*NU1)/2
      DO 12 J=1,NU2MAX+1
      NU2=J-1
      K1MAX=MAX-2*NU1-2*NU2
      DO 13 K=1,K1MAX+1
      K1=K-1
      LEFT=K1MAX-K1
      K2=LEFT
      J1MIN=IABS(K1-L1)
      J1MAX=K1+L1
      J2MIN=IABS(K2-L2)
      J2MAX=K2+L2
      IJ1=J1MAX-J1MIN+1
      IJ2=J2MAX-J2MIN+1
      DO 14 M=1,IJ1
      J1=J1MIN+M-1
      DO 15 NP=1,IJ2
      J2=J2MIN+NP-1
      ISUM=2*NU1+2*NU2+K1+K2
      X1=FLOAT(2*N1+L1+3+J1-K1)/2.-FLOAT(NU1)
      Y1=FLOAT(2*N1+L1-K1-J1)/2.-FLOAT(NU1)+1
      X2=FLOAT(2*N2+L2+3+J2-K2)/2.-FLOAT(NU2)
      Y2=FLOAT(2*N2+L2-K2-J2)/2.-FLOAT(NU2)+1
      IF((X1.LT.-1.).OR.(X2.LT.-1.)) GOTO 15
      IF((Y1.LT.1.).OR.(Y2.LT.1.)) GOTO 15
      JTL=2*L
      JTLL=2*LL
      JTK1=2*K1
      JTK2=2*K2
      JTJ1=2*J1
      JTJ2=2*J2
      JTLB=2*LB
      JTL1=2*L1
      JTL2=2*L2
      BMOSH=BMOSH+PAR(FLOAT(K2))*
     &((HATI(K1)*HATI(K2)*HATI(J1)*HATI(J2))**2.)*
     $COEF9J(JTK1,JTK2,JTJ1,JTJ2,JTL,JTLL,JTL1,JTL2,JTLB)*
     $CLEBSCH(K1,J1,L1,K2,J2,L2,K1,K2,L,J1,J2,LL)*
     $FFAC(NU1,N1,L1,K1,J1)*FFAC(NU2,N2,L2,K2,J2)
  15  CONTINUE
  14  CONTINUE
  13  CONTINUE
  12  CONTINUE
  11  CONTINUE
  10  CONTINUE
      BMOSH=ANORM*BMOSH
      RETURN
      END
      real function DELTAI(I1,I2)
      INTEGER I1,I2
      DELTAI=0.0E0
      IF(I1.EQ.I2)DELTAI=1.0E0
      RETURN
      END
      real function DELTAR(X1,X2)
      REAL X1,X2
      DELTAR=0.0E0
      IF(ABS(X1-X2).LT.1.0E-5)DELTAR=1.0E0
      RETURN
      END
      FUNCTION CHARM(N,L,K)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      XL=2*L+1
      XN=N
      XK=K
      ACHARM=PAR(XK)*XFAC(XN)*DFAC(XL)*2.**(XK)
      BCHARM=XFAC(XK)*XFAC(XN-XK)*DFAC(XL+2.*XK)
      CHARM=ACHARM/BCHARM
      RETURN
      END
      FUNCTION rme(n1,l1,n2,l2,lam)
      pi=2.*asin(1.0)
      xn1=float(n1)
      xl1=float(l1)
      xn2=float(n2)
      xl2=float(l2)
      ANORM=2.**(XL1-XN1+2)*DFAC(2.*XL1+2.*XN1+1)
      BNORM=XFAC(XN1)*(DFAC(2.*XL1+1)**2.)*SQRT(pi)
      XNOR1=SQRT(ANORM/BNORM)
      ANORM=2.**(XL2-XN2+2)*DFAC(2.*XL2+2.*XN2+1)
      BNORM=XFAC(XN2)*(DFAC(2.*XL2+1)**2.)*SQRT(pi)
      XNOR2=SQRT(ANORM/BNORM)
      rme=0.0
      I1=INT(XN1)+1
      I2=INT(XN2)+1
      rm=0.0
      do 10 k1=1,I1
      do 10 k2=1,I2
      xk1=FLOAT(k1-1)
      xk2=FLOAT(k2-1)
      rm=rm+CHARM(N1,L1,K1-1)*CHARM(N2,L2,K2-1)*
     1XGAMMA((2.*xk1+2.*xk2+xl1+xl2+lam+3.)/2.)
10    CONTINUE
      rm=0.5*xnor1*xnor2*rm
      rme=rm
      return
      end
      
