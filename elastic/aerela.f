      PROGRAM AERELA
!     0PROGRAM AERELA(INPUT=201,OUTPUT=201,TAPE5=INPUT,TAPE6=OUTPUT,     ARLA0010
!     . TAPE20=201)                                                      ARLA0020
C***********************************************************************ARLA0030
C                                                                       ARLA0040
C  AERELA USES THE GEOMETRY AND STRUCTURAL DATA FILES CREATED IN AEREAD ARLA0050
C  AND AERLAS TO COMPUTE THE ELASTIC AND RIGID STABILITY DERIVATIVES    ARLA0060
C  AND INDUCED DRAG FOR A THIN ELASTIC AIRPLANE.                        ARLA0070
C                                                                       ARLA0080
C                                                                       ARLA0100
      COMMON IND(1)                                                     ARLA0110
C***********************************************************************ARLA0090
      COMMON/CNTROL/N,NEAMAX,NCNTL,HALFSW,CREF,AM,CLDES,                ARLA0120
     . NCP2,NCP3,NCP4,WNGPNL,HT1PNL,HT2PNL,NCPHT1,NCPHT2,NWNGP,         ARLA0130
     . NHT1P,NHT2P,HALFB,ALPH,HALFB1,HALFB2,NCPSWL,NCPWNG,MAX,IALPHA    ARLA0140
      DIMENSION D(1)                                                    ARLA0150
      EQUIVALENCE (IND(1),D(1))                                         ARLA0160
      COMMON/INDMS/ IA,IB,IC                                            ARLA0170
      INTEGER WNGPNL,HT1PNL,HT2PNL                                      ARLA0180
      CALL OPENMS(20,IND,14,0)                                          ARLA0190
      CALL READMS(20,N,26,1)                                            ARLA0200
C                                                                       ARLA0210
C  PARTITION COMMON FOR PROBLEM ARRAYS                                  ARLA0220
C                                                                       ARLA0230
C  DISC ADDRESSES OF MATRICES                                           ARLA0240
      IA=1                                                              ARLA0250
      IB=IA+N                                                           ARLA0260
      IC=IB+N                                                           ARLA0270
C                                                                       ARLA0280
C  CORE STORAGE FOR N X N MATRIX                                        ARLA0290
      IBASE=15                                                          ARLA0300
      ICBASE=IC+IBASE                                                   ARLA0310
      IMM=IBASE+3*N+1                                                   ARLA0320
      NN=N*N                                                            ARLA0330
C                                                                       ARLA0340
C  CONTROL POINTS                                                       ARLA0350
      IXP=IMM+NN                                                        ARLA0360
      IYP=IXP+N                                                         ARLA0370
      IZP=IYP+N                                                         ARLA0380
C                                                                       ARLA0390
C  PANEL CG, AREA, MASS, CP-FLAT                                        ARLA0400
      IXG=IZP+N                                                         ARLA0410
      IAR=IXG+N                                                         ARLA0420
      IAM=IAR+N                                                         ARLA0430
      ICP=IAM+N                                                         ARLA0440
C                                                                       ARLA0450
C  PANEL CORNER COORDINATES, SLOPES                                     ARLA0460
      N4=4*N                                                            ARLA0470
      IXN=ICP+N                                                         ARLA0480
      IYN=IXN+N4                                                        ARLA0490
      IZN=IYN+N4                                                        ARLA0500
      IBP=IZN+N4                                                        ARLA0510
C                                                                       ARLA0520
C  MULTIPLE ALLOCATIONS (EQUIVALENCING).FOR WORKING                     ARLA0530
C  STORAGE, LOCAL ALPHAS, CP-CAMBERED                                   ARLA0540
      IPT=IXG                                                           ARLA0550
      IWA=IZN                                                           ARLA0560
      IWB=IWA+N                                                         ARLA0570
      IWC=IWB+N                                                         ARLA0580
      IWI=IWC+N                                                         ARLA0590
      IAL=IWI+N                                                         ARLA0600
      ICL=IAL+N                                                         ARLA0610
      IWS=ICL+N                                                         ARLA0620
C                                                                       ARLA0630
C  READ DATA NEEDED FOR DOWNWASH MATRIX FROM DISC.                      ARLA0640
      CALL DISCRD(IND(1),D(IXP),D(IYP),D(IZP),D(IXG),D(IAR),D(IAM),     ARLA0650
     .            D(IXN),D(IYN),D(IZN),D(IBP),D(IMM),IND(IPT),N)        ARLA0660
C                                                                       ARLA0670
C  READ C-MATRIX DISC INDEX FOR ELASTIC PROBLEM                         ARLA0680
      IF(NCNTL.NE.1) CALL READMS(20,IND(ICBASE),N,12)                   ARLA0690
C                                                                       ARLA0700
C                                                                       ARLA0710
C  COMPUTE AERODYNAMIC DOWNWASH MATRIX AND RIGID DERIVATIVES.           ARLA0720
C  COMPUTE RIGID FORCE COEFFICIENTS FOR FLAT AND CAMBERED SURFACES.     ARLA0730
C                                                                       ARLA0740
      CALL LINK1L(IND(IBASE),IND(IBASE),IND(IPT),D(IMM),D(IXP),D(IYP),  ARLA0750
     .            D(IZP),D(IXG),D(IAR),D(IAM),D(ICP),D(IXN),D(IYN),     ARLA0760
     .            D(IZN),D(IBP),D(IWA),D(IWB),D(IWC),D(IWI),D(IAL),     ARLA0770
     .            D(ICL),N)                                             ARLA0780
      IF(NCNTL.EQ.1) GO TO 3                                            ARLA0790
C                                                                       ARLA0800
C                                                                       ARLA0810
C  COMPUTE ELASTIC DERIVATIVES.  COMPUTE ELASTIC FORCE COEFFICIENTS     ARLA0820
C  FOR FLAT AND CAMBERED SURFACES.                                      ARLA0830
C                                                                       ARLA0840
      CALL LINK3L(IND(IBASE),D(IMM),D(IXP),D(IYP),D(IZP),D(IXG),D(IAR), ARLA0850
     .            D(IAM),D(ICP),D(IXN),D(IYN),D(IWA),D(IWB),D(IWC),     ARLA0860
     .            D(IWI),D(IAL),D(ICL),D(IWS),N)                        ARLA0870
    3 WRITE (6,7)                                                       ARLA0880
    7 FORMAT (1H1,55X,21HTHIS CASE IS COMPLETE,/,56X,21H****************ARLA0890
     .*****)                                                            ARLA0900
C  CHANGE INDEX KEY TO SET UP DISC FOR NEXT CASE                        ARLA0910
      CALL STINDX(20,IND(1),12)                                         ARLA0920
      STOP                                                              ARLA0930
      END                                                               ARLA0940
      SUBROUTINE DISCRD(IND,XCP,YCP,ZCP,XCG,AREA,AMASS,XN,YN,ZN,BPRIM,  DISC0010
     . DUM,IPNT,N)                                                      DISC0020
C.......................................................................DISC0030
C     SUBROUTINE DISCRD READS THE DATA STORED ON TAPE20 INTO CORE       DISC0040
C.......................................................................DISC0050
      DIMENSION XCP(N),YCP(N),ZCP(N),XCG(N),AREA(N),AMASS(N),           DISC0060
     . XN(N,4),YN(N,4),ZN(N,4),BPRIM(N,4),IPNT(N,4)                     DISC0070
      DIMENSION IND(12),DUM(4,N)                                        DISC0080
      N4 = 4 *N                                                         DISC0090
      CALL READMS(20,XCP(1),N,2)                                        DISC0100
      CALL READMS(20,YCP(1),N,3)                                        DISC0110
      CALL READMS(20,ZCP(1),N,4)                                        DISC0120
      CALL READMS(20,IPNT,N4,13)                                        DISC0130
      CALL RESTORE(DUM,IND,N,N4,8,XN)                                   DISC0210
      CALL RESTORE(DUM,IND,N,N4,9,YN)                                   DISC0220
      CALL RESTORE(DUM,IND,N,N4,10,ZN)                                  DISC0230
      CALL RESTORE(DUM,IND,N,N4,11,BPRIM)                               DISC0240
      RETURN                                                            DISC0250
      END                                                               DISC0260
C
      SUBROUTINE RESTORE(DUM,IND,N,N4,NREC,STORE)                       REST0010
      DIMENSION DUM(4,N),IND(1),STORE(N,4)                              REST0020
C                                                                       REST0030
      CALL READMS(20,DUM,N4,NREC)                                       REST0040
      DO 1 I=1,N                                                        REST0050
      DO 1 J=1,4                                                        REST0060
    1 STORE(I,J)=DUM(J,I)                                               REST0070
      RETURN                                                            REST0080
      END                                                               REST0090
      SUBROUTINE LINK1L(IND,KIND,IPNT,AW,XCP,YCP,ZCP,XCG,AREA,AMASS,CP, LN1L0010
     .                  XN,YN,ZN,BPRIM,WA,WB,WC,WI,ALPHAL,CPL,N)        LN1L0020
C                                                                       LN1L0030
C***********************************************************************LN1L0040
C                                                                       LN1L0050
C  LINK1L DIRECTS COMPUTATION OF THE AERODYNAMIC INFLUENCE MATRIX,      LN1L0060
C  RIGID DERIVATIVES AND INDUCED DRAG, AND PRESSURE COEFFICIENTS.       LN1L0070
C                                                                       LN1L0080
C***********************************************************************LN1L0090
C                                                                       LN1L0100
      COMMON/INDMS/IA,IB,IC                                             LN1L0110
      COMMON/CNTROL/LPANEL,NEAMAX,NCNTL,HALFSW,CREF,AM,CLDES,           LN1L0120
     . NCP2,NCP3,NCP4,WNGPNL,HT1PNL,HT2PNL,NCPHT1,NCPHT2,NWNGP,         LN1L0130
     . NHT1P,NHT2P,HALFB,ALPHA,HALF1,HALF2,NCPSWL,NCPWNG,MAX,IALPHA     LN1L0140
C                                                                       LN1L0150
      DIMENSION IND(1),KIND(1),IPNT(N,4),AW(N,N),XCP(N),YCP(N),ZCP(N),  LN1L0160
     .          XCG(N),AREA(N),AMASS(N),CP(N),XN(N,4),YN(N,4),ZN(N,4),  LN1L0170
     .          BPRIM(N,4),WA(N),WB(N),WC(N),WI(N),ALPHAL(N),CPL(N)     LN1L0180
C                                                                       LN1L0190
      INTEGER WNGPNL,HT1PNL,HT2PNL                                      LN1L0200
C                                                                       LN1L0210
C                                                                       LN1L0220
C  COMPUTE DOWNWASH MATRIX.                                             LN1L0230
C                                                                       LN1L0240
      CALL ADICMX(AM,XN,YN,ZN,XCP,YCP,ZCP,N,BPRIM,AW, N ,IPNT,IPNT)     LN1L0250
C                                                                       LN1L0260
C                                                                       LN1L0270
C  RETRIEVE C.G., AREA AND MASS ARRAYS FROM DISC, THEN CHANGE           LN1L0280
C  DISC FILE KEY TO MAKE RECORD NO. ONE THE FIRST ROW OF A.             LN1L0290
C                                                                       LN1L0300
      IBASE=15                                                          LN1L0310
      N3=3*N+1                                                          LN1L0320
      CALL READMS(20,XCG(1),N,5)
      CALL READMS(20,AREA(1),N,6)
      CALL READMS(20,AMASS(1),N,7)
      CALL STINDX(20,KIND,N3)                                           LN1L0340
C                                                                       LN1L0350
C                                                                       LN1L0360
C                                                                       LN1L0370
C  COMPUTE RIGID PRESSURE, LIFT AND MOMENT COEFFICIENTS FOR             LN1L0380
C    (1) ZERO INCIDENCE ANGLE                                           LN1L0390
C    (2) INCIDENCE ANGLE YIELDING DESIRED LIFT COEFFICIENT (CLDES)      LN1L0400
C  IF CLDES NOT SPECIFIED, INCIDENCE ANGLE ALPHA IS SET TO 1.0 RADIANS. LN1L0410
C                                                                       LN1L0420
C  INITIALIZE                                                           LN1L0430
C                                                                       LN1L0440
      ALPHA=1.0                                                         LN1L0450
      DO 10  JJ = 1,N                                                   LN1L0460
      CPL(JJ) = 0.                                                      LN1L0470
      WI(JJ) =  1.0                                                     LN1L0480
      ALPHAL(JJ) = 0.0                                                  LN1L0490
   10 CONTINUE                                                          LN1L0500
C                                                                       LN1L0510
C                                                                       LN1L0520
C  INVERT DOWNWASH MATRIX AND MULTIPLY ONTO AREA MATRIX, * (-2.0)       LN1L0530
C  TO OBTAIN THE AERODYNAMIC INFLUENCE COEFFICIENT MATRIX A.            LN1L0540
C                                                                       LN1L0550
      CALL MINV(AW,N)                                                   LN1L0560
      DO 15  J = 1,N                                                    LN1L0570
      DO 15  I = 1,N                                                    LN1L0580
   15 AW(I,J) = -2. * AREA(I) * AW(I,J)                                 LN1L0590
C                                                                       LN1L0600
C  STORE MATRIX A ON DISC                                               LN1L0610
      DO 30 I=1,N                                                       LN1L0620
      IAA=IA+I-1                                                        LN1L0630
      DO 25 J=1,N                                                       LN1L0640
   25 WA(J)=AW(I,J)                                                     LN1L0650
   30 CALL WRITMS(20,WA,N,IAA)                                          LN1L0660
C                                                                       LN1L0670
C  COMPUTE RIGID STABILITY DERIVATIVES.                                 LN1L0680
C                                                                       LN1L0690
      DO 35 I=1,N                                                       LN1L0700
   35 XCG(I)=-XCG(I)                                                    LN1L0710
      CALL DERIVR(HALFSW,CREF,N,AW,XCG,WA,WB,CLALPA)                    LN1L0720
C                                                                       LN1L0730
C  CAMBER EFFECTS. IF IALPHA = 0, NO CAMBER.                            LN1L0740
C                     IALPHA = 1, LOCAL ALPHAS INPUT IN RADIANS.        LN1L0750
C                     IALPHA =-1, SLOPES INPUT, CONVERT TO RADIANS.     LN1L0760
C                                                                       LN1L0770
      IF(IALPHA.EQ.0) GO TO 45                                          LN1L0780
C***********************************************************************LN1L0790
      READ(5,910)(ALPHAL(I),I=1,N)                                      LN1L0800
C***********************************************************************LN1L0810
      IF(IALPHA.NE.-1)GO TO 45                                          LN1L0820
      DO 40 J=1,N                                                       LN1L0830
   40 ALPHAL(J)=ATAN(-ALPHAL(J))                                        LN1L0840
C                                                                       LN1L0850
C                                                                       LN1L0860
C  COMPUTE FORCE COEFFICIENTS FOR ZERO INCIDENCE ANGLE.                 LN1L0870
C                                                                       LN1L0880
C  PRESSURE COEFFICIENTS                                                LN1L0890
   45 CALL MATMPY(AW,ALPHAL,WA,N,N,1)                                   LN1L0900
      DO 50 J = 1,N                                                     LN1L0910
   50 CPL(J)=WA(J)/AREA(J)                                              LN1L0920
C                                                                       LN1L0930
C  LIFT COEFFICIENT                                                     LN1L0940
      CALL TRNPRD(WI,WA,CLO,N)                                          LN1L0950
      CLO = CLO/HALFSW                                                  LN1L0960
C                                                                       LN1L0970
C  MOMENT COEFFICIENT                                                   LN1L0980
      CALL TRNPRD(XCG, WA, CMO, N)                                      LN1L0990
      CMO = CMO/(CREF * HALFSW )                                        LN1L1000
C                                                                       LN1L1010
C  COMPUTE INCIDENCE ANGLE CORRESPONDING TO DESIRED LIFT                LN1L1020
C  COEFFICIENT (CLDES).  IF CLDES = 0, SET INCIDENCE ANGLE              LN1L1030
C  TO 1.0.  ADD INCIDENCE ANGLE TO LOCAL ALPHAS.                        LN1L1040
C                                                                       LN1L1050
   60 IF(CLDES .EQ.0.)GO TO 65                                          LN1L1060
      ALPHA = (CLDES - CLO) / CLALPA                                    LN1L1070
   65 DO 70 I =1,N                                                      LN1L1080
      WA(I) = ALPHAL(I) +ALPHA                                          LN1L1090
   70 CONTINUE                                                          LN1L1100
C                                                                       LN1L1110
C  COMPUTE FORCE COEFFICIENTS FOR INCIDENCE ANGLE ALPHA.                LN1L1120
C  PLUS LOCAL ALPHAS.                                                   LN1L1130
C                                                                       LN1L1140
C  PRESSURE COEFFICIENTS                                                LN1L1150
      CALL MATMPY(AW,WA,CP,N,N,1)                                       LN1L1160
C                                                                       LN1L1170
C  MOMENT COEFFICIENT                                                   LN1L1180
      CALL TRNPRD(XCG,CP,CMD,N)                                         LN1L1190
      CMD=CMD/(CREF*HALFSW)                                             LN1L1200
C                                                                       LN1L1210
      DO 75 J =1,N                                                      LN1L1220
      CP(J) = CP(J)/AREA(J)                                             LN1L1230
   75 CONTINUE                                                          LN1L1240
C                                                                       LN1L1250
C  CL DESIRED, IF NOT INPUT                                             LN1L1260
      IF(CLDES .NE.0.)GO TO 80                                          LN1L1270
      CLDES = CLO + CLALPA                                              LN1L1280
   80 CONTINUE                                                          LN1L1290
C                                                                       LN1L1300
      DO 90 I=1,N                                                       LN1L1310
   90 XCG(I)=-XCG(I)                                                    LN1L1320
C                                                                       LN1L1330
C  OUTPUT                                                               LN1L1340
      WRITE(6,900)                                                      LN1L1350
      IF(IALPHA) 110,100,110                                            LN1L1360
  100 WRITE(6,901)                                                      LN1L1370
      GO TO 120                                                         LN1L1380
  110 WRITE(6,902)                                                      LN1L1390
  120 CONTINUE                                                          LN1L1400
      WRITE(6,903)AM                                                    LN1L1410
      WRITE(6,904)CLDES                                                 LN1L1420
      WRITE(6,908)CLO,CMO                                               LN1L1430
      WRITE(6,905) ALPHA                                                LN1L1440
      WRITE(6,909) CMD                                                  LN1L1450
      WRITE(6,907)                                                      LN1L1460
      WRITE(6,906)(I,XCP(I),YCP(I),ZCP(I),XCG(I),AREA(I),AMASS(I),      LN1L1470
     . ALPHAL(I),CP(I),I = 1,N )                                        LN1L1480
C                                                                       LN1L1490
C  COMPUTE SECTIONAL COEFFICIENTS, INDUCED DRAG                         LN1L1500
C                                                                       LN1L1510
      CL = CLDES                                                        LN1L1520
C                                                                       LN1L1530
C                                                                       LN1L1540
      CALL SETDRAG(XN,YN,YCP,AREA,CP,XCG,WA,ALPHAL,CL,AW,N)             LN1L1550
 900  FORMAT(1H1)                                                       LN1L1560
 901  FORMAT(42X,40HRIGID AIRCRAFT WITH FLAT LIFTING SURFACE///)        LN1L1570
 902  FORMAT(40X,44HRIGID AIRCRAFT WITH CAMBERED LIFTING SURFACE,///)   LN1L1580
 903  FORMAT(55X,13HMACH NUMBER =,F5.2,/55X,18H******************,//)   LN1L1590
 904  FORMAT(48X,25HDESIRED LIFT COEFFICIENT=,F10.5,/,48X,35H***********LN1L1600
     .************************,/)                                       LN1L1610
 905  FORMAT(35X,40HCP VALUES ARE FOR AN INCIDENCE ANGLE OF ,F10.6, 8H RLN1L1620
     .ADIANS,/,35X,58H**************************************************LN1L1630
     .********/)                                                        LN1L1640
 906  FORMAT((7X,I3,5(5X,F10.5),5X,F10.7,5X,F10.5,5X,E12.5))            LN1L1650
 907  FORMAT(7X,5HPANEL,8X,3HXCP,12X,3HYCP,10X,3HZCP,15X,3HXCG,11X,4HARELN1L1660
     .A,10X,4HMASS,10X,6HALPHAL,10X,2HCP,/,7X,6HNUMBER )                LN1L1670
 908  FORMAT(34X,25HFOR INCIDENCE ANGLE = 0  ,5HCL = ,F10.5,5X,5HCM = , LN1L1680
     .F10.5,/34X,60(1H*),/)                                             LN1L1690
 909  FORMAT(45X,"AT THIS INCIDENCE ANGLE  CM =",F10.5/45X,39(1H*)///)  LN1L1700
 910  FORMAT(8F10.5)                                                    LN1L1710
      RETURN                                                            LN1L1720
      END                                                               LN1L1730
      SUBROUTINE ADICMX(AM,XN,YN,ZN,XCP,YCP,ZCP,N,BPRIM,AW,             ADIC0010
     .                  MAX,IPNT,WSTUFF)                                ADIC0020
C                                                                       ADIC0030
C***********************************************************************ADIC0040
C                                                                       ADIC0050
C  ADICMX COMPUTES AERODYNAMIC DOWNWASH INFLUENCE COEFFICIENT MATRIX    ADIC0060
C  USING WOODWARD'S METHOD.  AN OPTIMIZATION ALGORITHM ELIMINATES       ADIC0070
C  REDUNDANT EVALUATIONS OF THE LIFTING SURFACE INTEGRAL AT CORNERS     ADIC0080
C  COMMON TO TWO OR MORE PANELS                                         ADIC0090
C                                                                       ADIC0100
C***********************************************************************ADIC0110
C                                                                       ADIC0120
      DIMENSION XN(N,4),YN(N,4),ZN(N,4),BPRIM(N,4),IPNT(N,4),           ADIC0130
     .          WSTUFF(N,4),XCP(N),YCP(N),ZCP(N),AW(N,N),               ADIC0140
     .          W(4),TWO(2),CON(2),A(2)                                 ADIC0150
C                                                                       ADIC0160
      DATA TWO/1.,-1./,CON/0.25,0.50/,PAI/3.14159265/,EPSLON/1.E-8/     ADIC0170
C                                                                       ADIC0180
C  MACH NUMBER REGIME                                                   ADIC0190
      IJ=1                                                              ADIC0200
      IF(AM.GE.1.0) IJ=2                                                ADIC0210
      BETA1=1.0/SQRT(ABS(1.0-AM*AM))                                    ADIC0220
C                                                                       ADIC0230
C  MULTIPLY PANEL SLOPES BY TANGENT OF MACH ANGLE                       ADIC0240
      PAIBET=1.0/(BETA1*PAI)                                            ADIC0250
      DO 10 I=1,N                                                       ADIC0260
      DO 10 J=1,4                                                       ADIC0270
 10   BPRIM(I,J)=BETA1*BPRIM(I,J)                                       ADIC0280
C                                                                       ADIC0290
C  DOWNWASH MATRIX COMPUTATION                                          ADIC0300
C       I IS INFLUENCED PANEL                                           ADIC0310
C       J IS INFLUENCING PANEL                                          ADIC0320
C       K IS CORNER OF PANEL J                                          ADIC0330
C       II IS THE SIDE OF SYMMETRIC AIRPLANE BEING CONSIDERED           ADIC0340
C                                                                       ADIC0350
      F5=0.0                                                            ADIC0360
      DO 600 I=1,N                                                      ADIC0370
      DO 600 J=1,N                                                      ADIC0380
      DO 500 K=1,4                                                      ADIC0390
C                                                                       ADIC0400
C  PRESENCE OF POINTER IN IPNT MEANS W(K) HAS BEEN COMPUTED             ADIC0410
C  FOR A PRIOR J.  UNPACK POINTER AND ACCESS W(K) IN WSTUFF.            ADIC0420
C                                                                       ADIC0430
      ITAB=IPNT(J,K)                                                    ADIC0440
      IF(ITAB.GT.3004.OR.ITAB.LT.1) GO TO 100                           ADIC0450
      JJ=ITAB/10                                                        ADIC0460
      KK=MOD(ITAB,10)                                                   ADIC0470
      W(K)=WSTUFF(JJ,KK)                                                ADIC0480
      GO TO 500                                                         ADIC0490
C                                                                       ADIC0500
C  COMPUTATIONS COMMON TO SUBSONIC AND SUPERSONIC REGIMES.              ADIC0510
 100  CONTINUE                                                          ADIC0520
      BPM=BPRIM(J,K)                                                    ADIC0530
      XIPRIM=(XCP(I)-XN(J,K))*BETA1                                     ADIC0540
      ZPRIM=ZCP(I)-ZN(J,K)                                              ADIC0550
C                                                                       ADIC0560
C  KK IS PANEL SLOPE SIGN FLAG                                          ADIC0570
      KK=1                                                              ADIC0580
      IF(BPM.GE.0.0) GO TO 110                                          ADIC0590
      KK=2                                                              ADIC0600
      BPM=-BPM                                                          ADIC0610
 110  IF(BPM.LE.EPSLON) BPM=0.0                                         ADIC0620
      IF(ABS(BPM-1.).LE.EPSLON) BPM=1.0                                 ADIC0630
      BPM2=BPM*BPM                                                      ADIC0640
      BB=ABS(1.-BPM2)                                                   ADIC0650
      IF(IJ.EQ.2) BB=SQRT(BB)                                           ADIC0660
      IF(BB.GT.EPSLON) GO TO 20                                         ADIC0670
      BPM=1.0                                                           ADIC0680
      BB=0.0                                                            ADIC0690
 20   B1=XIPRIM*XIPRIM                                                  ADIC0700
C                                                                       ADIC0710
C  Y-PRIME DEPENDS ON II                                                ADIC0720
      DO 400 II=1,2                                                     ADIC0730
      YPRIM=TWO(II)*YCP(I)-YN(J,K)                                      ADIC0740
      IF(KK.EQ.2) YPRIM=-YPRIM                                          ADIC0750
      B2=YPRIM*YPRIM+ZPRIM*ZPRIM                                        ADIC0760
      GO TO (200,300),IJ                                                ADIC0770
C                                                                       ADIC0780
C  SUBSONIC REGIME                                                      ADIC0790
C                                                                       ADIC0800
C  FOR COMPUTERS OTHER THAN THE CDC 6000 SERIES, THE FOLLOWING          ADIC0810
C  COMPUTATIONS MAY REQUIRE DOUBLE PRECISION EVALUATION IF              ADIC0820
C  EITHER A11 OR A4 IS LESS THAN (-1000.)                               ADIC0830
C                                                                       ADIC0840
 200  A1=XIPRIM+SQRT(B1+B2)                                             ADIC0850
      A11=XIPRIM/SQRT(B2)                                               ADIC0860
      A2=BPM*XIPRIM+YPRIM                                               ADIC0870
      A3=XIPRIM-BPM*YPRIM                                               ADIC0880
      A33=(BPM2+1.0)*ZPRIM*ZPRIM                                        ADIC0890
      SQA3=1./SQRT(A3*A3+A33)                                           ADIC0900
      A4=A2*SQA3                                                        ADIC0910
      F1=ALOG(A11+SQRT(A11*A11+1.0))                                    ADIC0920
      A7=SQRT(A4*A4+1)                                                  ADIC0930
      A8=A4+A7                                                          ADIC0940
      A6=ALOG(A8)                                                       ADIC0950
      F2=A6/SQRT(1.+BPM2)                                               ADIC0960
C  END DOUBLE PRECISION SECTION                                         ADIC0970
C                                                                       ADIC0980
      A5=BPM*SQRT(B2)*SQA3                                              ADIC0990
      F5=0.0                                                            ADIC1000
      IF(A5.GT.EPSLON) F5=ALOG(A5)                                      ADIC1010
      F6=A1/B2                                                          ADIC1020
      GO TO 400                                                         ADIC1030
C                                                                       ADIC1040
C  SUPERSONIC REGIME                                                    ADIC1050
 300  A1=B2                                                             ADIC1060
      A2=SQRT(A1)                                                       ADIC1070
      IF(XIPRIM.GT.A2) GO TO 340                                        ADIC1080
C                                                                       ADIC1090
      F1=0.0                                                            ADIC1100
      F2=0.0                                                            ADIC1110
      F6=0.0                                                            ADIC1120
      IF(BPM2.GT.1.0) GO TO 400                                         ADIC1130
C                                                                       ADIC1140
      TEST=YPRIM-BPM*XIPRIM                                             ADIC1150
      IF(XIPRIM.EQ.A2) GO TO 310                                        ADIC1160
      IF(YPRIM.LE.0.0) GO TO 400                                        ADIC1170
C                                                                       ADIC1180
      CONTL=BPM*YPRIM+BB*ABS(ZPRIM)                                     ADIC1190
      IF(XIPRIM-CONTL) 400,320,310                                      ADIC1200
 310  IF(TEST) 400,320,330                                              ADIC1210
 320  F2=1.57079633/BB                                                  ADIC1220
      GO TO 400                                                         ADIC1230
 330  F2=PAI/BB                                                         ADIC1240
      GO TO 400                                                         ADIC1250
C                                                                       ADIC1260
 340  A3=XIPRIM/A1                                                      ADIC1270
      SQXI=SQRT(XIPRIM*XIPRIM-A1)                                       ADIC1280
      F6=SQXI/A1                                                        ADIC1290
      A11=XIPRIM/A2                                                     ADIC1300
      F1=ALOG(A11+SQRT(A11*A11-1.0))                                    ADIC1310
      IF(BPM2.EQ.1.0) GO TO 350                                         ADIC1320
C                                                                       ADIC1330
      A4=XIPRIM-BPM*YPRIM                                               ADIC1340
      A5=(BPM2-1.0)*ZPRIM*ZPRIM                                         ADIC1350
      A6=SQRT(A4*A4+A5)                                                 ADIC1360
      A7=BPM*XIPRIM-YPRIM                                               ADIC1370
      IF(BPM2.GT.1.0) GO TO 360                                         ADIC1380
      F2=(1./BB)*ACOS(A7/A6)                                            ADIC1390
      GO TO 400                                                         ADIC1400
C                                                                       ADIC1410
 350  F2=SQXI/(XIPRIM-YPRIM)                                            ADIC1420
      GO TO 400                                                         ADIC1430
C                                                                       ADIC1440
 360  A8=A7/A6                                                          ADIC1450
      F2=(1./BB)*ALOG(A8+SQRT(A8*A8-1.0))                               ADIC1460
C                                                                       ADIC1470
 400  A(II)=CON(IJ)*PAIBET*((BPM2+TWO(IJ))*F2-BPM*(F1-F5)               ADIC1480
     .     -YPRIM*F6)*TWO(KK)                                           ADIC1490
      W(K)=A(1)+A(2)                                                    ADIC1500
C                                                                       ADIC1510
C  STORE W(K) IN WSTUFF TO AVOID RECOMPUTATION LATER.                   ADIC1520
      WSTUFF(J,K)=W(K)                                                  ADIC1530
 500  CONTINUE                                                          ADIC1540
 600  AW(I,J)=W(1)-W(2)-W(3)+W(4)                                       ADIC1550
      RETURN                                                            ADIC1560
      END                                                               ADIC1570
      SUBROUTINE DERIVR(AREA,CREF,LPANEL,A,X,A1,UNIT1,CLALPA)           DRVR0010
C                                                                       DRVR0020
C***********************************************************************DRVR0030
C  DERIVR COMPUTES THE RIGID STABILITY DERIVATIVES                      DRVR0040
C***********************************************************************DRVR0050
      DIMENSION A(LPANEL,LPANEL),X(LPANEL,1),A1(LPANEL,1),UNIT1(LPANEL,1DRVR0060
     .)                                                                 DRVR0070
C                                                                       DRVR0080
      WRITE(6,6)                                                        DRVR0090
    6 FORMAT(1H1)                                                       DRVR0100
      DO 1 I=1,LPANEL                                                   DRVR0110
    1 UNIT1(I,1)=1.                                                     DRVR0120
C                                                                       DRVR0130
C  CLALPA COMPUTATION                                                   DRVR0140
C  (A1)=(A)*(UNIT1)                                                     DRVR0150
C                                                                       DRVR0160
C  MATMPY COMPUTES THE PRODUCT OF TWO MATRICES.                         DRVR0170
C                                                                       DRVR0180
      CALL MATMPY(A,UNIT1,A1,LPANEL,LPANEL,1)                           DRVR0190
C  (A3)=(UNIT1)T*(A1)                                                   DRVR0200
C                                                                       DRVR0210
C  TRNPRD COMPUTES THE DOT PRODUCT OF TWO VECTORS.                      DRVR0220
C                                                                       DRVR0230
      CALL TRNPRD(UNIT1,A1,A3,LPANEL)                                   DRVR0240
      CLALPA=A3/AREA                                                    DRVR0250
      WRITE(6,2) CLALPA                                                 DRVR0260
C                                                                       DRVR0270
C  CMALPA COMPUTATION                                                   DRVR0280
C  (A4)=(X)T*A1                                                         DRVR0290
C                                                                       DRVR0300
      CALL TRNPRD(X,A1,A4,LPANEL)                                       DRVR0310
      CMALPA=A4/(AREA*CREF)                                             DRVR0320
      WRITE(6,3) CMALPA                                                 DRVR0330
C                                                                       DRVR0340
C  CLQ COMPUTATION                                                      DRVR0350
C  A1=A*X(I)                                                            DRVR0360
C                                                                       DRVR0370
      CALL MATMPY(A,X,A1,LPANEL,LPANEL,1)                               DRVR0380
C  (A3)=(UNIT1)T*(A1)                                                   DRVR0390
      CALL TRNPRD(UNIT1,A1,A3,LPANEL)                                   DRVR0400
      CLQ=-2.*A3/(AREA*CREF)                                            DRVR0410
      WRITE(6,4) CLQ                                                    DRVR0420
C                                                                       DRVR0430
C  CMQ COMPUTATION                                                      DRVR0440
C  (A4)=(X)T*(A1)                                                       DRVR0450
C                                                                       DRVR0460
      CALL TRNPRD(X,A1,A4,LPANEL)                                       DRVR0470
      CMQ=-2.*A4/(AREA*CREF*CREF)                                       DRVR0480
      WRITE(6,5) CMQ                                                    DRVR0490
    2 FORMAT (10X,9HCLALPR  =,E12.5)                                    DRVR0500
    3 FORMAT (10X,9HCMALPR  =,E12.5)                                    DRVR0510
    4 FORMAT (10X,9HCLQ     =,E12.5)                                    DRVR0520
    5 FORMAT (10X,9HCMQ     =,E12.5)                                    DRVR0530
      RETURN                                                            DRVR0540
      END                                                               DRVR0550
      SUBROUTINE LINK3L(IND,P,XCP,YCP,THETAE,XCG,AREA,AMASS,CP,XN,YN,   LN3L0010
     .                  WA,WB,WC,WI,ALPHAL,CPL,WS,N)                    LN3L0020
C                                                                       LN3L0030
C.......................................................................LN3L0040
C                                                                       LN3L0050
C  LINK3L DIRECTS COMPUTATION OF THE ELASTIC STABILITY DERIVATIVES,     LN3L0060
C  INDUCED DRAG AND PRESSURE COEFFICIENTS.                              LN3L0070
C                                                                       LN3L0080
C  THE P-ARRAY IS USED TO HOLD THE A AND B-MATRICES AS NEEDED.          LN3L0090
C                                                                       LN3L0100
C.......................................................................LN3L0110
      COMMON/INDMS/IA,IB,IC                                             LN3L0120
      DIMENSION IND(1),P(N,N),XCP(N),YCP(N), THETAE(N),XCG(N),AREA(N),  LN3L0130
     .          AMASS(N),CP(N),XN(N,4),YN(N,4),WA(N),WB(N),WC(N),       LN3L0140
     .          WI(N),ALPHAL(N),CPL(N),WS(N)                            LN3L0150
      COMMON/CNTROL/LPANEL,NEAMAX,NCNTL,HALFSW,CREF,AMACH,CLDES,        LN3L0160
     . NCP2,NCP3,NCP4,WNGPNL,HT1PNL,HT2PNL,NCPHT1,NCPHT2,NWNGP,         LN3L0170
     . NHT1P,NHT2P,HALFB,ALPHA,HALFB1,HALFB2,NCPSWL,NCPWNG,MAX,IALPHA   LN3L0180
      INTEGER WNGPNL,HT1PNL,HT2PNL                                      LN3L0190
      EQUIVALENCE (AM,AMACH)                                            LN3L0200
      EQUIVALENCE(CL,CLDES)                                             LN3L0210
C                                                                       LN3L0220
C***********************************************************************LN3L0230
      READ (5,19) W,ALTNUM                                              LN3L0240
C***********************************************************************LN3L0250
      NALT=ALTNUM                                                       LN3L0260
C                                                                       LN3L0270
C  ENTER DYNAMIC PRESSURE (Q) LOOP.                                     LN3L0280
C                                                                       LN3L0290
      DO 200 K=1,NALT                                                   LN3L0300
      WRITE(6,25) AMACH                                                 LN3L0310
      CLO=0.0                                                           LN3L0320
      WRITE (6,20) W                                                    LN3L0330
      DO 1 I = 1,N                                                      LN3L0340
    1 XCG(I)=-XCG(I)                                                    LN3L0350
C***********************************************************************LN3L0360
      READ (5,19) RHO,SOS,ALT                                           LN3L0370
      READ (5,19) Q                                                     LN3L0380
C***********************************************************************LN3L0390
      V=AMACH*SOS                                                       LN3L0400
      WRITE (6,21) RHO,SOS,ALT                                          LN3L0410
C                                                                       LN3L0420
      CLTRIM=W/(Q*2.*HALFSW)                                            LN3L0430
      G=32.2                                                            LN3L0440
      WRITE (6,23) Q                                                    LN3L0450
      WRITE (6,24) CLTRIM                                               LN3L0460
      NUM = N                                                           LN3L0470
C                                                                       LN3L0480
C  COMPUTE ELASTIC STABILITY DERIVATIVES.                               LN3L0490
C                                                                       LN3L0500
      CALL DERIV(AMASS,HALFSW,CREF,XCG,Q,G,V,CLTRIM,CLALPE,P,           LN3L0510
     . WA,WB,WC,WI,N,IND)                                               LN3L0520
      DO 10 L = 1,N                                                     LN3L0530
   10 WA(L) = 0.0                                                       LN3L0540
C                                                                       LN3L0550
C  COMPUTE ELASTIC FORCE COEFFICIENTS FOR LOCAL ALPHAS AT               LN3L0560
C  ZERO INCIDENCE ANGLE.                                                LN3L0570
C                                                                       LN3L0580
C  PRESSURE COEFFICIENTS.  PLACE A MATRIX IN P                          LN3L0590
      DO 15 I = 1,N                                                     LN3L0600
      IAA = IA + I -1                                                   LN3L0610
      CALL READMS(20,WB,N     ,IAA)                                     LN3L0620
      DO 15 J = 1,N                                                     LN3L0630
   15 P(I,J) = WB(J)                                                    LN3L0640
C                                                                       LN3L0650
C  INERTIAL = (C)(MASS)*G                                               LN3L0660
      CALL DISMPY(WS,AMASS,WC,N,N,1,20,IC,0,IND)                        LN3L0670
      DO 30 J=1,N                                                       LN3L0680
      WC(J)=G*WC(J)                                                     LN3L0690
C                                                                       LN3L0700
C  AERODYNAMIC = (ALPHAL)                                               LN3L0710
      WA(J)=ALPHAL(J)-WC(J)                                             LN3L0720
   30 CONTINUE                                                          LN3L0730
C                                                                       LN3L0740
C  (CP) = (B)(A)(AERODYNAMIC - INERTIAL)/(AREA)                         LN3L0750
      CALL  MATMPY(P,WA,WB,N,N,1)                                       LN3L0760
      CALL DISMPY(WS,WB,WC,N,N,1,20,IB,0,IND)                           LN3L0770
C                                                                       LN3L0780
C  CL, CM AT ZERO INCIDENCE ANGLE                                       LN3L0790
      CALL TRNPRD(WI,WC,CLO,N)                                          LN3L0800
      CLO =CLO /HALFSW                                                  LN3L0810
      CALL TRNPRD(XCG,WC,CMO,N)                                         LN3L0820
      CMO = CMO/(CREF* HALFSW)                                          LN3L0830
C                                                                       LN3L0840
      IF(IALPHA) 45,40,45                                               LN3L0850
   40 WRITE(6,28)                                                       LN3L0860
      GO TO 50                                                          LN3L0870
   45 WRITE(6,29)                                                       LN3L0880
   50 CONTINUE                                                          LN3L0890
      WRITE(6,31)AM                                                     LN3L0900
      WRITE(6,33)CLO,CMO                                                LN3L0910
C                                                                       LN3L0920
C  COMPUTE INCIDENCE ANGLE YIELDING DESIRED CL FOR ELASTIC              LN3L0930
C  AIRPLANE.                                                            LN3L0940
C                                                                       LN3L0950
      ALPHA=(CLDES-CLO)/CLALPE                                          LN3L0960
   60 DO 70 KL = 1,N                                                    LN3L0970
      WA(KL) = WA(KL) + ALPHA                                           LN3L0980
   70 CONTINUE                                                          LN3L0990
C                                                                       LN3L1000
C  PRESSURE COEFFICIENTS                                                LN3L1010
      CALL MATMPY(P,WA,WB,N,N,1)                                        LN3L1020
      CALL DISMPY(WS,WB,WC,N,N,1,20,IB,0,IND)                           LN3L1030
      DO 80 L =1,N                                                      LN3L1040
      CP(L) = WC(L)/AREA(L)                                             LN3L1050
   80 CONTINUE                                                          LN3L1060
C                                                                       LN3L1070
C  ELASTIC CM                                                           LN3L1080
C                                                                       LN3L1090
      CALL TRNPRD(XCG,WC,CMR,N)                                         LN3L1100
      CMR = CMR/(CREF * HALFSW)                                         LN3L1110
      WRITE(6,27) ALPHA,CLDES,CMR                                       LN3L1120
      DO 90 I =1,N                                                      LN3L1130
      WC(I) = Q * (ALPHAL(I)+ALPHA)                                     LN3L1140
   90 CONTINUE                                                          LN3L1150
      CALL MATMPY(P,WC,WA,N,N,1)                                        LN3L1160
      DO 100 I =1,N                                                     LN3L1170
      WA(I)=WA(I) - G * AMASS(I)                                        LN3L1180
  100 CONTINUE                                                          LN3L1190
C                                                                       LN3L1200
C  COMPUTATION OF DEFORMED ANGLE OF PANELS (THETAE)                     LN3L1210
C                                                                       LN3L1220
      CALL DISMPY(WC,WA,THETAE,N,N,1,20,IC,0,IND)                       LN3L1230
      CALL DISMPY(WC,P,WB,N,N,N,20,IC,IB,IND)                           LN3L1240
C                                                                       LN3L1250
C  READ MATRIX B FROM DISC AND STORE IT IN MATRIX P                     LN3L1260
C                                                                       LN3L1270
      DO 110 I =1,N                                                     LN3L1280
      IBB=IB+I-1                                                        LN3L1290
      CALL READMS(20,WB,N,IBB)                                          LN3L1300
      DO 110 J =1,N                                                     LN3L1310
  110 P(I,J) = WB(J)                                                    LN3L1320
C                                                                       LN3L1330
      DO 120 I =1,N                                                     LN3L1340
      DO 120 J =1,N                                                     LN3L1350
      P(I,J)=-Q*P(I,J)                                                  LN3L1360
      IF(I.EQ.J) P(I,J) = 1.0 + P(I,J)                                  LN3L1370
  120 CONTINUE                                                          LN3L1380
      CALL MINV(P,N)                                                    LN3L1390
      CALL MATMPY(P,THETAE,WA,N,N,1)                                    LN3L1400
C                                                                       LN3L1410
C  (NOTE THAT DEFORMED ANGLES ARE NOW STORED IN WA)                     LN3L1420
C                                                                       LN3L1430
      WRITE(6,26)                                                       LN3L1440
      WRITE(6,32)(I,ALPHAL(I),WA(I),CP(I),I =1,N)                       LN3L1450
C                                                                       LN3L1460
C                                                                       LN3L1470
C  COMPUTE SECTIONAL COEFFICIENTS, INDUCED DRAG.                        LN3L1480
C                                                                       LN3L1490
      DO 130 I = 1,N                                                    LN3L1500
      WA(I)=WA(I)+ALPHAL(I)+ALPHA                                       LN3L1510
  130 XCG(I)=-XCG(I)                                                    LN3L1520
      CALL SETDRAG(XN,YN,YCP,AREA,CP,XCG,WA,ALPHAL,CL,P,N)              LN3L1530
C                                                                       LN3L1540
  200 CONTINUE                                                          LN3L1550
      RETURN                                                            LN3L1560
   19 FORMAT (3F10.5)                                                   LN3L1570
   20 FORMAT (10X,7HWEIGHT=,F10.0)                                      LN3L1580
   21 FORMAT (10X,4HRHO=,F12.8,5X,7HS.O.S.=,F10.3,5X,9HALTITUDE=,F6.0)  LN3L1590
   22 FORMAT (10X,9HVELOCITY=,E12.5)                                    LN3L1600
   23 FORMAT (10X,17HDYNAMIC PRESSURE=,E12.5)                           LN3L1610
   24 FORMAT (10X,7HCLTRIM=,F10.5)                                      LN3L1620
   25 FORMAT (1H1,9X,36HELASTIC DERIVATIVES FOR MACH NUMBER=,F5.2,///)  LN3L1630
   26 FORMAT(25X,5HPANEL,15X,5HALPHA,15X,8HTHETA(E),15X,5HCP(E),/, 25X, LN3L1640
     . 6HNUMBER,14X,5HLOCAL,15X,6HALPHAE,17X,6HALPHAE,/)                LN3L1650
   27 FORMAT(24X,22HFOR INCIDENCE ANGLE = ,F10.6,3X,5HCL = ,F10.5,5X,   LN3L1660
     .5HCM = , F10.5,/24X,70(1H*),//)                                   LN3L1670
   28 FORMAT(1H1,38X,42HELASTIC AIRCRAFT WITH FLAT LIFTING SURFACE/)    LN3L1680
   29 FORMAT(1H1,36X,46HELASTIC AIRCRAFT WITH CAMBERED LIFTING SURFACE/)LN3L1690
   31 FORMAT(51X,13HMACH NUMBER =,F5.2,/51X,18(1H*),/)                  LN3L1700
   32 FORMAT(25X,I3,14X,E12.5, 8X,E12.5, 8X, E12.5)                     LN3L1710
   33 FORMAT(30X,25HFOR INCIDENCE ANGLE = 0  ,5HCL = ,F10.5,5X,5HCM = , LN3L1720
     . F10.5,/ 30X,60(1H*),/)                                           LN3L1730
      END                                                               LN3L1740
      SUBROUTINE DERIV(AMASS,AREA,CREF,X,Q,G,V,CLTRIM,CLALPE,P,         DERI0010
     . WA,WB,WC,WI,N,IND)                                               DERI0020
C.......................................................................DERI0030
C                             D*E*R*I*V                                 DERI0040
C      SUBROUTINE  DERIV  CALCULATES ALL THE ELASTIC DERIVATIVES        DERI0050
C.......................................................................DERI0060
      COMMON/INDMS/IA,IB,IC                                             DERI0070
      DIMENSION IND(1),THETAE(1),AMASS(1),CP(1),X(1),WA(1),WB(1),WC(1), DERI0080
     .          WI(1), P(N,N)                                           DERI0090
      DO 1 I=1,N                                                        DERI0100
    1 WI(I)=1.                                                          DERI0110
C                                                                       DERI0120
C  GENERATION OF (B) MATRIX                                             DERI0130
C                                                                       DERI0140
      DO 19 I=1,N                                                       DERI0150
      ICC=IC+I-1                                                        DERI0160
      CALL READMS(20,WB,N,ICC)                                          DERI0170
      DO 19 J=1,N                                                       DERI0180
   19 P(I,J) = WB(J)                                                    DERI0190
      CALL DISMPY(WA,P,WB,N,N,N, 20,IA,IB,IND)                          DERI0200
      DO 20 I=1,N                                                       DERI0210
      IBB=IB+I-1                                                        DERI0220
      CALL READMS(20,WB,N,IBB)                                          DERI0230
      DO 20 J=1,N                                                       DERI0240
   20 P(I,J) = WB(J)                                                    DERI0250
      DO 3 I=1,N                                                        DERI0260
      DO 3 J=1,N                                                        DERI0270
      P(I,J) = -Q * P(I,J)                                              DERI0280
      IF(I.EQ.J) P(I,J) = 1.0 + P(I,J)                                  DERI0290
    3 CONTINUE                                                          DERI0300
C                                                                       DERI0310
C  MINV INVERTS A WELL-CONDITIONED MATRIX                               DERI0320
C                                                                       DERI0330
      CALL MINV(P,N)                                                    DERI0340
      DO 21 I=1,N                                                       DERI0350
      IBB=IB+I-1                                                        DERI0360
      DO 22 J=1,N                                                       DERI0370
   22 WB(J) = P(I,J)                                                    DERI0380
   21 CALL WRITIN(20,WB,N,IBB)                                          DERI0390
C  THE B MATRIX HAS BEEN GENERATED AND THE INVERSE OF B                 DERI0400
C  HAS BEEN STORED ON DISC.                                             DERI0410
C                                                                       DERI0420
C  CLALPA COMPUTATION                                                   DERI0430
C                                                                       DERI0440
      CALL DISMPY(WB,WI,WA,N,N,1,20,IA,0,IND)                           DERI0450
      CALL MATMPY(P,WA,WB,N,N,1)                                        DERI0460
C                                                                       DERI0470
      CALL TRNPRD(WI,WB,A3,N)                                           DERI0480
      CLALPA=A3/AREA                                                    DERI0490
      WRITE(6,5) CLALPA                                                 DERI0500
C                                                                       DERI0510
C  CMALPA COMPUTATION                                                   DERI0520
C                                                                       DERI0530
      CALL TRNPRD(X,WB,A4,N)                                            DERI0540
      CMALPA=A4/(AREA*CREF)                                             DERI0550
      WRITE(6,6) CMALPA                                                 DERI0560
C                                                                       DERI0570
C  CLQEBR COMPUTATION                                                   DERI0580
C                                                                       DERI0590
      CALL DISMPY(WB,X,WA,N,N,1,20,IA,0,IND)                            DERI0600
      CALL MATMPY(P,WA,WB,N,N,1)                                        DERI0610
      CALL TRNPRD(WI,WB,A3,N)                                           DERI0620
      CLQEBR=-2.*A3/(AREA*CREF)                                         DERI0630
      WRITE(6,7) CLQEBR                                                 DERI0640
C                                                                       DERI0650
C  CMQEBR COMPUTATION                                                   DERI0660
C                                                                       DERI0670
      CALL TRNPRD(X,WB,A4,N)                                            DERI0680
      CMQEBR=-2.*A4/(AREA*CREF*CREF)                                    DERI0690
      WRITE(6,8) CMQEBR                                                 DERI0700
C                                                                       DERI0710
C  STABILITY DERIVATIVES AT VARYING LOAD FACTOR                         DERI0720
C                                                                       DERI0730
      AC=(AREA*CREF)/(2.*V*V)                                           DERI0740
C                                                                       DERI0750
C  CLQI COMPUTATION                                                     DERI0760
C                                                                       DERI0770
      CALL DISMPY(WB,AMASS,WA,N,N,1,20,IC,0,IND)                        DERI0780
      CALL DISMPY(WC,WA,WB,N,N,1,20,IA,0,IND)                           DERI0790
      CALL MATMPY(P,WB,WC,N,N,1)                                        DERI0800
      CALL TRNPRD(WI,WC,ACLQI,N)                                        DERI0810
      CLQI=ACLQI/AC                                                     DERI0820
      WRITE(6,9) CLQI                                                   DERI0830
C                                                                       DERI0840
C  CMQI COMPUTATION                                                     DERI0850
C                                                                       DERI0860
      CALL TRNPRD(X,WC,ACMQI,N)                                         DERI0870
      CMQI=ACMQI/(AC*CREF)                                              DERI0880
      WRITE(6,10) CMQI                                                  DERI0890
C                                                                       DERI0900
C  CLWDI COMPUTATION                                                    DERI0910
C                                                                       DERI0920
      CLWDI=-ACLQI/AREA                                                 DERI0930
      WRITE(6,11) CLWDI                                                 DERI0940
C                                                                       DERI0950
C  CMWDI COMPUTATION                                                    DERI0960
C                                                                       DERI0970
      CMWDI=-ACMQI/(AREA*CREF)                                          DERI0980
      WRITE(6,12) CMWDI                                                 DERI0990
      DO 4 I=1,N                                                        DERI1000
    4 WB(I)=AMASS(I)*X(I)                                               DERI1010
C                                                                       DERI1020
C  CLTDDI COMPUTATION                                                   DERI1030
C                                                                       DERI1040
      CALL DISMPY(WC,WB,WA,N,N,1,20,IC,0,IND)                           DERI1050
      CALL DISMPY(WC,WA,WB,N,N,1,20,IA,0,IND)                           DERI1060
      CALL MATMPY(P,WB,WC,N,N,1)                                        DERI1070
      CALL TRNPRD(WI,WC,ALTDDI,N)                                       DERI1080
      CLTDDI=ALTDDI/AREA                                                DERI1090
      WRITE(6,13) CLTDDI                                                DERI1100
C                                                                       DERI1110
C  CMTDDI COMPUTATION                                                   DERI1120
C                                                                       DERI1130
      CALL TRNPRD( X,WC,AMTDDI,N)                                       DERI1140
      CMTDDI=AMTDDI/(AREA*CREF)                                         DERI1150
      WRITE(6,14) CMTDDI                                                DERI1160
C                                                                       DERI1170
C  DCLDN COMPUTATION                                                    DERI1180
C                                                                       DERI1190
      DCLDN=-32.2*CLWDI                                                 DERI1200
      WRITE(6,15) DCLDN                                                 DERI1210
C                                                                       DERI1220
C  DCMDN COMPUTATION                                                    DERI1230
C                                                                       DERI1240
      DCMDN=-32.2*CMWDI                                                 DERI1250
      WRITE(6,16) DCMDN                                                 DERI1260
C                                                                       DERI1270
C  CLALPAE COMPUTATION                                                  DERI1280
C                                                                       DERI1290
      CLALPE=CLALPA/(1.-DCLDN/CLTRIM)                                   DERI1300
      WRITE(6,17) CLALPE                                                DERI1310
C                                                                       DERI1320
C  CMALPAE COMPUTATION                                                  DERI1330
C                                                                       DERI1340
      CMALPE=CMALPA+DCMDN*CLALPA/(CLTRIM-DCLDN)                         DERI1350
      WRITE(6,18) CMALPE                                                DERI1360
    5 FORMAT (10X,7HCLALPA=,E12.5)                                      DERI1370
    6 FORMAT (10X,7HCMALPA=,E12.5)                                      DERI1380
    7 FORMAT (10X,7HCLQEBR=,E12.5)                                      DERI1390
    8 FORMAT (10X,7HCMQEBR=,E12.5)                                      DERI1400
    9 FORMAT (10X,7HCLQI  =,E12.5)                                      DERI1410
   10 FORMAT (10X,7HCMQI  =,E12.5)                                      DERI1420
   11 FORMAT (10X,7HCLWDI =,E12.5)                                      DERI1430
   12 FORMAT (10X,7HCMWDI =,E12.5)                                      DERI1440
   13 FORMAT (10X,7HCLTDDI=,E12.5)                                      DERI1450
   14 FORMAT (10X,7HCMTDDI=,E12.5)                                      DERI1460
   15 FORMAT (10X,7HDCLDN =,E12.5)                                      DERI1470
   16 FORMAT (10X,7HDCMDN =,E12.5)                                      DERI1480
   17 FORMAT (10X,7HCLALPE=,E12.5)                                      DERI1490
   18 FORMAT (10X,7HCMALPE=,E12.5)                                      DERI1500
      RETURN                                                            DERI1510
      END                                                               DERI1520
      SUBROUTINE SETDRAG(XN,YN,YCP,AREA,CP,XCG,WA,ALPHAL,CL,P,N)        SETD0010
C...................................................................... SETD0020
C   SETDRAG COMPUTES THE LEADING EDGE THRUST COEFFICIENT                SETD0030
C   AND THE INDUCED DRAG COEFFICIENT                                    SETD0040
C...................................................................... SETD0050
      COMMON/CNTROL/ LPANEL,NEAMAX,NCNTL,HALFSW,CREF,AM,CLDES,          SETD0060
     . NCP2,NCP3,NCP4,WNGPNL,HT1PNL,HT2PNL,NCPHT1,NCPHT2,NWNGP,         SETD0070
     . NHT1P,NHT2P,HALFB,ALPH,HALFB1,HALFB2,NCPSWL,NCPWNG,MAX           SETD0080
      INTEGER WNGPNL,HT1PNL,HT2PNL                                      SETD0090
      DIMENSION XN(N,4),YN(N,4),YCP(1),AREA(1),CP(1),XCG(1),WA(1)       SETD0100
      DIMENSION HEAD(6),P(1),ALPHAL(1)                                  SETD0110
      DATA HEAD/5HWING ,5HHORIZ,5HONTAL,5H TAIL,5H(ONE),5H(TWO)/        SETD0120
C                                                                       SETD0130
      LFUS=WNGPNL-1                                                     SETD0140
      CLF=0.                                                            SETD0150
      CDF=0.0                                                           SETD0160
      IF(LFUS.LE.0) GO TO 20                                            SETD0170
C  WA IS THE SUM OF THE (LOCAL +DEFORMED+INCIDENCE)ANGLE                SETD0180
      DO 10 I=1,LFUS                                                    SETD0190
      DCL=CP(I)*AREA(I)                                                 SETD0200
      CLF=CLF+DCL                                                       SETD0210
 10   CDF=CDF+DCL*WA(I)                                                 SETD0220
 20   CONTINUE                                                          SETD0230
C                                                                       SETD0240
      WRITE(6,15) HEAD(1)                                               SETD0250
      CTT=0.                                                            SETD0260
      CD=CDF                                                            SETD0270
C                                                                       SETD0280
C  SET UP FOR FIRST CALL TO DRAG.                                       SETD0290
      NA=NCP2                                                           SETD0300
      NAM=NCP2-1                                                        SETD0310
      NP=NWNGP                                                          SETD0320
      NW=WNGPNL                                                         SETD0330
      HB=HALFB                                                          SETD0340
      ASSIGN 30 TO JBACK                                                SETD0350
      GO TO 100                                                         SETD0360
 30   IF(HT1PNL.LE.WNGPNL) GO TO 200                                    SETD0370
C  SET UP FOR SECOND CALL AND WRITE HEADING.                            SETD0380
      NA=NCP3                                                           SETD0390
      NAM=NCP3-1                                                        SETD0400
      NP=NHT1P                                                          SETD0410
      NW=HT1PNL                                                         SETD0420
      HB=HALFB1                                                         SETD0430
C                                                                       SETD0440
      IH=4                                                              SETD0450
      IF(HT2PNL.GT.HT1PNL) IH=5                                         SETD0460
      WRITE(6,15) (HEAD(I),I=2,IH)                                      SETD0470
      ASSIGN 40 TO JBACK                                                SETD0480
      GO TO 100                                                         SETD0490
C                                                                       SETD0500
C  SET UP FOR THIRD CALL                                                SETD0510
 40   IF(HT2PNL.LE.HT1PNL) GO TO 200                                    SETD0520
      NA=NCP4                                                           SETD0530
      NAM=NCP4-1                                                        SETD0540
      NW=HT2PNL                                                         SETD0550
      HB=HALFB2                                                         SETD0560
      WRITE(6,15) (HEAD(I),I=2,4),HEAD(6)                               SETD0570
      ASSIGN 200 TO JBACK                                               SETD0580
C  CALL DRAG AND INCREMENT SUMS                                         SETD0590
C                                                                       SETD0600
 100 0CALL DRAG(XN,YN,XCG,YCP,AREA,CP,WA,P(711),N,NA,NAM,NP,NW,HB,AM,CL,SETD0610
     . CTA,CDA,P(1),P(36),P(71),P(106),P(141),P(176),P(211),P(246),     SETD0620
     . P(281),P(316),P(351),P(386))                                     SETD0630
      CTT=CTT+CTA                                                       SETD0640
      CD=CD+CDA                                                         SETD0650
      GO TO JBACK,(30,40,200)                                           SETD0660
C                                                                       SETD0670
C  COMPLETE DRAG COMPUTATIONS AND RETURN                                SETD0680
 200  CD=CD/HALFSW                                                      SETD0690
      CTT=CTT/HALFSW                                                    SETD0700
      CDBCL=CD/(CL*CL)                                                  SETD0710
      WRITE(6,17) CTT,CD,CDBCL                                          SETD0720
      RETURN                                                            SETD0730
 15   FORMAT(//,5X,23HSECTIONAL CD/CL**2 FOR ,6A5)                      SETD0740
 17  0FORMAT(//,5X,38HTOTAL LEADING EDGE THRUST COEFFICIENT=,E12.5/5X,  SETD0750
     . 31HTOTAL INDUCED DRAG COEFFICIENT=,E12.5/5X, 9HCD/CL**2=,E12.5)  SETD0760
      END                                                               SETD0770
     0SUBROUTINE DRAG(XN,YN,XCG,YCP,AREA,CP,ALPHA,A,K,NCPCWL,NN,NCPSWL, DRAG0010
     . IPANEL,HALFB,AM,CL,CTT,CD,DELTL,DELTY,XCGL,DELTC,ZBAR,CT,CLA,CDL,DRAG0020
     . GAMMA,CAK,Y,AK)                                                  DRAG0030
C                                                                       DRAG0040
C.......................................................................DRAG0050
C  DRAG COMPUTES THE INDUCED DRAG DISTRIBUTION                          DRAG0060
C.......................................................................DRAG0070
      DIMENSION XN(K,4),YN(K,4),XCG(1),YCP(1),AREA(1),CP(1),ALPHA(1)    DRAG0080
      DIMENSION A(NN,1),DELTL(1),DELTY(1),XCGL(1),DELTC(1),ZBAR(1)      DRAG0090
      DIMENSION CT(1),CLA(1),CDL(1),GAMMA(1),CAK(1),Y(1),AK(1)          DRAG0100
      NCC=NCPCWL-1                                                      DRAG0110
      NCS=NCPSWL                                                        DRAG0120
      K2=IPANEL-1                                                       DRAG0130
      DO 1 I=1,NCS                                                      DRAG0140
      K4=K2+NCC                                                         DRAG0150
      K3=K2+1                                                           DRAG0160
      DELTL(I)=SQRT((XN(K3,2)-XN(K3,1))**2+(YN(K3,2)-YN(K3,1))**2)      DRAG0170
      DELTY(I)=YN(K3,2)-YN(K3,1)                                        DRAG0180
      XCGL(I)=XN(K3,1)+(XN(K3,2)-XN(K3,1))*(YCP(K3)-YN(K3,1))/DELTY(I)  DRAG0190
      XCGT=XN(K4,3)+(XN(K4,4)-XN(K4,3))*(YCP(K3)-YN(K4,3))/DELTY(I)     DRAG0200
      DELTC(I)=XCGT-XCGL(I)                                             DRAG0210
      K2=K4                                                             DRAG0220
    1 CONTINUE                                                          DRAG0230
      N=IPANEL                                                          DRAG0240
      DO 2 I=1,NCC                                                      DRAG0250
      II=N+I-1                                                          DRAG0260
    2 ZBAR(I)=(XCG(II)-XCGL(1))/DELTC(1)                                DRAG0270
      NCC1=NCC+1                                                        DRAG0280
      KK=1                                                              DRAG0290
C                                                                       DRAG0300
C  PRESSURE DISTRIBUTION INTERPOLATION IN STREAMWISE DIRECTION          DRAG0310
C                                                                       DRAG0320
    3 CONTINUE                                                          DRAG0330
      AB=AM*AM*(DELTY(KK)*DELTY(KK)/(DELTL(KK)*DELTL(KK)))              DRAG0340
      IF (AB .GE. 1.0) GO TO 6                                          DRAG0350
      DO 5 J=1,NCC                                                      DRAG0360
      JN=J+N-1                                                          DRAG0370
      CAK(J)=-CP(JN)*SQRT(ZBAR(J)/(1.-ZBAR(J)))                         DRAG0380
      DO 4 I=1,NCC                                                      DRAG0390
      A(J,I)=ZBAR(J)**(I-1)                                             DRAG0400
    4 CONTINUE                                                          DRAG0410
    5 CONTINUE                                                          DRAG0420
C                                                                       DRAG0430
      CALL MINV(A,NCC)                                                  DRAG0440
      CALL MATMPY(A,CAK,AK,NCC,NCC,1)                                   DRAG0450
      AC = SQRT(1.0 - AB)                                               DRAG0460
      CTL=(3.14159265*DELTL(KK)*AC)/(8.0*DELTY(KK))                     DRAG0470
      CT(KK)=-CTL*AK(1)*AK(1)                                           DRAG0480
      GO TO 7                                                           DRAG0490
    6 CT(KK)=0.0                                                        DRAG0500
    7 N=N+NCC                                                           DRAG0510
      KK = KK+1                                                         DRAG0520
      IF (KK.GT.NCS) GO TO 8                                            DRAG0530
      GO TO 3                                                           DRAG0540
    8 CONTINUE                                                          DRAG0550
      K2=IPANEL-1                                                       DRAG0560
      CD = 0.0                                                          DRAG0570
      CTT = 0.0                                                         DRAG0580
      III=K2                                                            DRAG0590
      DO 10 I=1,NCS                                                     DRAG0600
      CLA(I) =0.0                                                       DRAG0610
      CDL(I)=0.                                                         DRAG0620
      III=(I-1)*NCC+K2                                                  DRAG0630
      DO 9 J=1,NCC                                                      DRAG0640
      II = III + J                                                      DRAG0650
      DCY=DELTC(I)*DELTY(I)                                             DRAG0660
      TT=(CP(II)*AREA(II))/DCY                                          DRAG0670
      CDL(I)=CDL(I)+TT*ALPHA(II)                                        DRAG0680
    9 CLA(I) = CLA(I) + TT                                              DRAG0690
      CDL(I)=CDL(I)+CT(I)                                               DRAG0700
      CTT=CTT+CT(I)*DCY                                                 DRAG0710
      GAMMA(I)=CLA(I)*DELTC(I)/(4.0*HALFB)                              DRAG0720
      CD=CD+CDL(I)*DCY                                                  DRAG0730
   10 III=III+NCC                                                       DRAG0740
      N=IPANEL                                                          DRAG0750
      WRITE (6,13)                                                      DRAG0760
      J=N                                                               DRAG0770
      DO 11 I=1,NCS                                                     DRAG0780
      CDL2=CDL(I)/(CL*CL)                                               DRAG0790
      Y(I)=YCP(J)/HALFB                                                 DRAG0800
      J=J+NCC                                                           DRAG0810
   11 WRITE (6,12) Y(I),CDL(I),CLA(I),GAMMA(I),CT(I),CDL2               DRAG0820
   12 FORMAT (6(5X,F10.5))                                              DRAG0830
   13 FORMAT (11X,1HY,13X,3HCDI,12X,3HCLI,11X,5HGAMMA,11X,2HCT,10X,8HCD/DRAG0840
     .CL**2)                                                            DRAG0850
      RETURN                                                            DRAG0860
      END                                                               DRAG0870
      SUBROUTINE MINV(A,NN)                                             MINV0010
C...................................................................... MINV0020
C**GAUSSIAN ELIMINATION ROUTINE TO INVERT A WELL-CONDITIONED,           MINV0030
C  DIAGONALLY DOMINANT MATRIX USING THE DIAGONAL ELEMENTS               MINV0040
C  FOR PIVOTING.                                                        MINV0050
C...................................................................... MINV0060
      DIMENSION A(NN,NN)                                                MINV0070
      DO 30 I=1,NN                                                      MINV0080
      PIV=1./A(I,I)                                                     MINV0090
      A(I,I)=1.0                                                        MINV0100
      DO 31 L=1,NN                                                      MINV0110
 31   A(I,L)=A(I,L)*PIV                                                 MINV0120
C                                                                       MINV0130
      DO 30 M=1,NN                                                      MINV0140
      IF(M-I) 33,30,33                                                  MINV0150
 33   TT=A(M,I)                                                         MINV0160
      A(M,I)=0.0                                                        MINV0170
C                                                                       MINV0180
      DO 32 L=1,NN                                                      MINV0190
 32   A(M,L)=A(M,L)-A(I,L)*TT                                           MINV0200
 30   CONTINUE                                                          MINV0210
      RETURN                                                            MINV0220
      END                                                               MINV0230
      SUBROUTINE MATMPY(A,B,C,L,M,N)                                    MATM0010
C                                                                       MATM0020
C  PRODUCT OF TWO MATRICES                                              MATM0030
C                                                                       MATM0040
      DIMENSION A(L,M),B(M,N),C(L,N)                                    MATM0050
      DO 1 I=1,L                                                        MATM0060
      DO 1 J=1,N                                                        MATM0070
      C(I,J)=0.                                                         MATM0080
      DO 1 K=1,M                                                        MATM0090
    1 C(I,J)=C(I,J)+A(I,K)*B(K,J)                                       MATM0100
      RETURN                                                            MATM0110
      END                                                               MATM0120
      SUBROUTINE DISMPY(A,B,C,L,M,N,LU,IA,IC,IND)                       DISM0010
      DIMENSION A(M),B(M,N),C(N),IND(1)                                 DISM0020
C                                                                       DISM0030
C  DISMPY MULTIPLIES THE L*M MATRIX A STORED BY ROWS STARTING           DISM0040
C  AT DISC INDEX IA ONTO THE M*N MATRIX B RESIDING IN CORE.  THE        DISM0050
C  N*L PRODUCT MATRIX C IS STORED BY ROWS STARTING AT DISC INDEX        DISM0060
C  IC                                                                   DISM0070
C                                                                       DISM0080
C  A AND C ARE SINGLY-SUBSCRIPTED ARRAYS USED FOR DISC I/O.             DISM0090
C  IF N IS 1, THE PRODUCT IS RETURNED IN C WITHOUT WRITING TO DISC.     DISM0100
C                                                                       DISM0110
C  IF IC IS NEGATIVE WRITIN IS USED RATHER THAN WRITMS.                 DISM0120
C  IND IS THE DISC INDEX ARRAY.  LU IS THE DISC LOGICAL UNIT NUMBER.    DISM0130
C                                                                       DISM0140
      JWRIT=1                                                           DISM0150
      IF(IC.GT.0) GO TO 10                                              DISM0160
      JWRIT=2                                                           DISM0170
      IC=-IC                                                            DISM0180
10    IF(N.EQ.1) JWRIT=3                                                DISM0190
      IAA=IA                                                            DISM0200
      ICC=IC                                                            DISM0210
C                                                                       DISM0220
      DO 100 I=1,L                                                      DISM0230
      CALL READMS(LU,A,M,IAA)                                           DISM0240
      DO 50 J=1,N                                                       DISM0250
      C(J)=0.0                                                          DISM0260
      DO 50 K=1,M                                                       DISM0270
 50   C(J)=C(J)+A(K)*B(K,J)                                             DISM0280
C                                                                       DISM0290
C  GO TO WRITMS, WRITIN OR NO DISC WRITE                                DISM0300
      GO TO (60,70,80), JWRIT                                           DISM0310
 60   CALL WRITMS(LU,C,N,ICC)                                           DISM0320
      GO TO 90                                                          DISM0330
 70   CALL WRITIN(LU,C,N,ICC)                                           DISM0340
      GO TO 90                                                          DISM0350
C                                                                       DISM0360
C  FOR N=1, STORE RESULT TEMPORARILY IN C(I+1)                          DISM0370
 80   IF(I.LT.L) C(I+1)=C(1)                                            DISM0380
C  INCREMENT DISC INDICES                                               DISM0390
 90   IAA=IAA+1                                                         DISM0400
      ICC=ICC+1                                                         DISM0410
 100  CONTINUE                                                          DISM0420
      IF(N.GT.1) RETURN                                                 DISM0430
C  RESTORE C-VECTOR FOR N=1                                             DISM0440
      A(1)=C(1)                                                         DISM0450
      DO 120 I=2,L                                                      DISM0460
 120  C(I-1)=C(I)                                                       DISM0470
      C(L)=A(1)                                                         DISM0480
      RETURN                                                            DISM0490
      END                                                               DISM0500
      SUBROUTINE TRNPRD(X,Y,Z,LPANEL)                                   TRNP0010
C                                                                       TRNP0020
C  DOT PRODUCT OF TWO VECTORS                                           TRNP0030
C                                                                       TRNP0040
      DIMENSION X(LPANEL,1),Y(LPANEL,1)                                 TRNP0050
      Z=0.                                                              TRNP0060
      DO 1 I=1,LPANEL                                                   TRNP0070
    1 Z=Z+X(I,1)*Y(I,1)                                                 TRNP0080
      RETURN                                                            TRNP0090
      END                                                               TRNP0100
      SUBROUTINE WRITIN(LUN,A,N,IREC)                                   WRTN0010
      DIMENSION A(N)                                                    WRTN0020
      CALL WRITMS(LUN,A(1),N,IREC)                                      WRTN0030
      RETURN                                                            WRTN0040
      END                                                               WRTN0050
