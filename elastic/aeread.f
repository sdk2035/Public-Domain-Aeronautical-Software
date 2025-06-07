      PROGRAM AEREAD
!!!      PROGRAM AEREAD(INPUT=201,OUTPUT=201,TAPE5=INPUT,TAPE6=OUTPUT,     AERE0010
!!!     . TAPE20=401)                                                      AERE0020
C***********************************************************************AERE0030
C                                                                       AERE0040
C  THE UNIVERSITY OF KANSAS AEROELASTIC PROGRAM COMPUTES ALPHA AND      AERE0050
C  Q STABILITY DERIVATIVES FOR A THIN ELASTIC AIRPLANE AT SUBSONIC      AERE0060
C  AND SUPERSONIC SPEEDS.  THE AIRPLANE MAY HAVE CAMBER.  MAXIMUM       AERE0070
C  PROBLEM SIZE IS 300 PANELS.                                          AERE0080
C                                                                       AERE0090
C  THE PROGRAM CONSISTS OF THREE JOB STEPS                              AERE0100
C                                                                       AERE0110
C    JOB STEP AEREAD PROCESSES INPUT DATA AND PANELS THE AIRPLANE.      AERE0120
C    JOB STEP AERLAS COMPUTES THE STRUCTURAL (C) MATRIX.                AERE0130
C    JOB STEP AERELA COMPUTES THE RIGID AND ELASTIC AERODYNAMIC         AERE0140
C       MATRICES AND STABILITY DERIVATIVES                              AERE0150
C                                                                       AERE0160
C  COMMUNICATION BETWEEN JOB STEPS IS VIA RANDOM ACCESS FILE TAPE20.    AERE0170
C  THE RANDOM ACCESS FILE ROUTINES ARE CDC-6000 DEPENDENT AND MUST      AERE0180
C  BE MODIFIED FOR USE ON OTHER COMPUTERS.                              AERE0190
C                                                                       AERE0200
C  PROGRAM MODIFICATIONS BY COMPUTER SCIENCES CORPORATION FOR           AERE0210
C  LANGLEY RESEARCH CENTER INCLUDE                                      AERE0220
C     FIELD LENGTH REDUCTION                 CSC 45351  JULY 1972       AERE0230
C     C-MATRIX TAPE OPTION                   CSC 45258  SEPT 1972       AERE0240
C     DOWNWASH MATRIX OPTIMIZATION           CSC 45378  JAN. 1973       AERE0250
C     INCLUSION OF CAMBER EFFECTS            CSC 45383  JAN. 1973       AERE0260
C                                                                       AERE0270
C  RESPONSIBLE INDIVIDUALS.  DR. JOHN E. LAMAR,  LRC  827-3711          AERE0280
C                            CHARLES W. BOLZ,JR, CSC  826-1725          AERE0290
C                                                                       AERE0300
C***********************************************************************AERE0310
C                                                                       AERE0320
C  INPUT FORMAT CONTROL.                                                AERE0330
C                                                                       AERE0340
C                                                                       AERE0350
C                                                                       AERE0360
C     IF  NCONTL=0  USE INPUT FORMATS OF NASA CR-112229                 AERE0370
C     IF  NCONTL=1  USE INPUT FORMATS OF NASA TM-X 2074                 AERE0380
C                                                                       AERE0390
C.......................................................................AERE0400
      COMMON/CNTROL/ LPANEL,NEAMAX,NCNTL,HALFSW,CREF,AM,CLDES,          AERE0410
     . NCP2,NCP3,NCP4,WNGPNL,HT1PNL,HT2PNL,NCPHT1,NCPHT2,NWNGP,         AERE0420
     . NHT1P,NHT2P,HALFB,ALPH,HALFB1,HALFB2,NCPSWL,NCPWNG,MAX,IALPHA    AERE0430
C                                                                       AERE0440
C  AEREAD BLANK COMMON IS SET UP FOR A MAXIMUM                          AERE0450
C  PROBLEM SIZE OF 300 PANELS.                                          AERE0460
C                                                                       AERE0470
      COMMON XCG(300),XCP(300),YCP(300),ZCP(300),AREA(300),AMASS(300),  AERE0480
     . XN(300,4),YN(300,4),ZN(300,4),BPRIM(300,4)                       AERE0490
      INTEGER WNGPNL,HT1PNL,HT2PNL                                      AERE0500
C                                                                       AERE0510
      DIMENSION DUM(25),IPNT(300,4),IRRAY(1200),RANDOM(4,300),TITLE(16) AERE0520
     .         ,IND(14)                                                 AERE0530
      EQUIVALENCE (LPANEL,N,NUM,DUM),(NEAMAX,M),(XCG,RANDOM,IPNT,IRRAY) AERE0540
C                                                                       AERE0550
      DO 10 I=1,25                                                      AERE0560
   10 DUM(I)=0.                                                         AERE0570
      MAX=300                                                           AERE0580
C***********************************************************************AERE0590
      READ (5,4) TITLE                                                  AERE0600
      READ (5,5) N,M,NCNTL,IALPHA                                       AERE0610
      READ (5,7) NCONTL                                                 AERE0620
C***********************************************************************AERE0630
      WRITE (6,6) TITLE                                                 AERE0640
C                                                                       AERE0650
C  SUBROUTINE DRGDTA READS THE GEOMETRY DATA ACCORDING TO THE           AERE0660
C  FORMAT GIVEN IN NASA-TMX-2074 AND MAKES AERODYNAMIC PANELS           AERE0670
C  CONFORMAL TO THIS PROGRAM                                            AERE0680
C                                                                       AERE0690
      IF(NCONTL.GT.0)CALL DRGDTA                                        AERE0700
C                                                                       AERE0710
C  SUBROUTINE NEWPNL USES THE INPUT DATA FORMAT OF NASA CR-112229       AERE0720
C  TO DETERMINE AERODYNAMIC PANELS.                                     AERE0730
C                                                                       AERE0740
      IF(NCONTL.EQ.0)CALL NEWPNL                                        AERE0750
C                                                                       AERE0760
C  FLCHEK COMPUTES THE FIELD LENGTH REQUIRED TO EXECUTE THIS PROBLEM    AERE0770
C  AND TERMINATES EXECUTION IF DECLARED FIELD LENGTH IS INSUFFICIENT.   AERE0780
C                                                                       AERE0790
      CALL FLCHEK(N)                                                    AERE0800
C***********************************************************************AERE0810
      READ (5,8) AM,CLDES                                               AERE0820
      IF(NCNTL.NE.1)READ(5,8)(AMASS(I),I=1,LPANEL)                      AERE0830
C***********************************************************************AERE0840
      N4=4*N                                                            AERE0850
C                                                                       AERE0860
      CALL OPENMS(20,IND,14,0)                                          AERE0870
      CALL WRITMS(20,LPANEL,26,1)                                       AERE0880
      CALL WRITMS(20,  XCP(1),N, 2)                                     AERE0890
      CALL WRITMS(20,  YCP(1),N, 3)                                     AERE0900
      CALL WRITMS(20,  ZCP(1),N, 4)                                     AERE0910
      CALL WRITMS(20,  XCG(1),N, 5)                                     AERE0920
      CALL WRITMS(20, AREA(1),N, 6)                                     AERE0930
      CALL WRITMS(20,AMASS(1),N, 7)                                     AERE0940
C                                                                       AERE0950
C  ARRAYS XN,YN,ZN AND BPRIM ARE STORED WITH ROWS AND COLUMNS           AERE0960
C  REVERSED TO FIT DATA STRUCTURE OF JOBSTEP AERELA.                    AERE0970
C                                                                       AERE0980
      DO 20 I=1,300                                                     AERE0990
      DO 20 J=1,4                                                       AERE1000
   20 RANDOM(J,I)=   XN(I,J)                                            AERE1010
      CALL WRITMS(20,RANDOM(1,1),N4, 8)                                 AERE1020
      DO 21 I=1,300                                                     AERE1030
      DO 21 J=1,4                                                       AERE1040
   21 RANDOM(J,I)=   YN(I,J)                                            AERE1050
      CALL WRITMS(20,RANDOM(1,1),N4, 9)                                 AERE1060
      DO 22 I=1,300                                                     AERE1070
      DO 22 J=1,4                                                       AERE1080
   22 RANDOM(J,I)=   ZN(I,J)                                            AERE1090
      CALL WRITMS(20,RANDOM(1,1),N4,10)                                 AERE1100
      DO 23 I=1,300                                                     AERE1110
      DO 23 J=1,4                                                       AERE1120
   23 RANDOM(J,I)=BPRIM(I,J)                                            AERE1130
      CALL WRITMS(20,RANDOM(1,1),N4,11)                                 AERE1140
C                                                                       AERE1150
C  CREATE POINTER TABLE FOR DOWNWASH MATRIX OPTIMIZATION                AERE1160
C                                                                       AERE1170
      CALL POINTER(IPNT,N)                                              AERE1180
      K=0                                                               AERE1190
      DO 24  J =1,4                                                     AERE1200
      DO 24 I = 1,LPANEL                                                AERE1210
      K = K+1                                                           AERE1220
      IRRAY(K)=IPNT(I,J)                                                AERE1230
   24 CONTINUE                                                          AERE1240
      CALL WRITMS(20,IRRAY,N4,13)                                       AERE1250
   25 CONTINUE                                                          AERE1260
C                                                                       AERE1270
C                                                                       AERE1280
    4 FORMAT (16A5)                                                     AERE1290
    5 FORMAT(4I3)                                                       AERE1300
    6 FORMAT (1H1,/////,25X,16A5,//////)                                AERE1310
    7 FORMAT (I2)                                                       AERE1320
    8 FORMAT(8F10.5)                                                    AERE1330
      STOP                                                              AERE1340
      END                                                               AERE1350
      SUBROUTINE POINTER(IPNT,N)                                        POIN0010
C                                                                       POIN0020
C  POINTER CREATES THE REFERENCE TABLE USED TO ELIMINATE                POIN0030
C  REDUNDANT PANEL CORNER COMPUTATIONS IN ADICMX.                       POIN0040
C                                                                       POIN0050
      COMMON DUMMY(1800),XN(300,4),YN(300,4),ZN(300,4),BP(300,4)        POIN0060
      DIMENSION IPNT(300,4)                                             POIN0070
C                                                                       POIN0080
      DO 50 K=1,4                                                       POIN0090
      DO 50 J=1,300                                                     POIN0100
 50   IPNT(J,K)=0                                                       POIN0110
C                                                                       POIN0120
      DO 100 JP=1,N                                                     POIN0130
      DO 100 KP=1,4                                                     POIN0140
C                                                                       POIN0150
C  HAS A POINTER BEEN PLACED IN THIS CORNER                             POIN0160
      IF(IPNT(JP,KP).NE.0) GO TO 100                                    POIN0170
C                                                                       POIN0180
C  CORNER 4 CANNOT HAVE A POINTER.                                      POIN0190
      DO 100 K=1,3                                                      POIN0200
C  TWO CORNERS WITH THE SAME NUMBER CANNOT COINCIDE.                    POIN0210
      IF(KP.EQ.K) GO TO 100                                             POIN0220
      DO 100 J=JP,N                                                     POIN0230
      IF(IPNT(J,K).NE.0) GO TO 100                                      POIN0240
C                                                                       POIN0250
C  DOES CORNER (J,K) MATCH CORNER (JP,KP)                               POIN0260
      IF(XN(J,K).NE.XN(JP,KP)) GO TO 100                                POIN0270
      IF(YN(J,K).NE.YN(JP,KP)) GO TO 100                                POIN0280
      IF(ZN(J,K).NE.ZN(JP,KP)) GO TO 100                                POIN0290
      IF(BP(J,K).NE.BP(JP,KP)) GO TO 100                                POIN0300
C                                                                       POIN0310
C  CORNERS MATCH.  SET IPNT(J,K) TO POINT TO (JP,KP)                    POIN0320
      IPNT(J,K)=10*JP+KP                                                POIN0330
 100  CONTINUE                                                          POIN0340
      RETURN                                                            POIN0350
      END                                                               POIN0360
      SUBROUTINE FLCHEK(N)                                              FLCH0010
C...................................................................... FLCH0020
C  FLCHEK DETERMINES PROBLEM FIELD LENGTH REQUIREMENT AND STOPS         FLCH0030
C  EXECUTION IF FIELD LENGTH IS INSUFFICIENT.                           FLCH0040
C...................................................................... FLCH0050
      DIMENSION NFL(1),O(4),L(6)                                        FLCH0060
      INTEGER CFL,RFL,FL,O                                              FLCH0070
      DATA CFL/10835/,O/32768,4096,512,64/,L(5),L(6)/2*0/               FLCH0080
      CFL=8116                                                          FLCH0090
C  CFL IS OBJECT CODE LENGTH.  RFL IS FL REQUIRED BY N-PANEL PROBLEM.   FLCH0100
C  FL IS USER-ASSIGNED FIELD LENGTH.                                    FLCH0110
C                                                                       FLCH0120
      RFL=CFL+26*N+N*N                                                  FLCH0130
      RFL=MAX0(RFL,19264)                                               FLCH0140
      K=LOCF(NFL)                                                       FLCH0150
      FL=NFL(63-K)                                                      FLCH0160
      RFL=((RFL-1)/64 + 1)*64                                           FLCH0170
      IF(FL.EQ.RFL) RETURN                                              FLCH0180
C                                                                       FLCH0190
C  CONVERT FIELD LENGTH TO OCTAL                                        FLCH0200
      L(1)=RFL/O(1)                                                     FLCH0210
      DO 10 I=2,4                                                       FLCH0220
      MFL=MOD(RFL,O(I-1))                                               FLCH0230
 10   L(I)=MFL/O(I)                                                     FLCH0240
 900 0FORMAT(1H0,"THE MINIMUM FIELD LENGTH REQUIREMENT FOR THIS PROBLEM FLCH0250
     .IS ",6I1," OCTAL WORDS.")                                         FLCH0260
      WRITE(6,900) L                                                    FLCH0270
C                                                                       FLCH0280
C  STOP JOB IF FIELD LENGTH IS TOO SMALL FOR PROBLEM.                   FLCH0290
      IF(FL.GE.RFL) RETURN                                              FLCH0300
      WRITE(6,901)                                                      FLCH0310
 901 0FORMAT(1H0,"YOUR FIELD LENGTH IS TOO SMALL FOR THIS PROBLEM.  EXECFLCH0320
     .UTION TERMINATED")                                                FLCH0330
      CALL EXIT                                                         FLCH0340
      END                                                               FLCH0350
      SUBROUTINE BREAK(X,Y,Z,N,YY,NN,NNN)                               BRAK0010
C                                                                       BRAK0020
C.......................................................................BRAK0030
C                                                                       BRAK0040
C  BREAK GENERATES BREAK LINES ON WING DUE TO HORIZONTAL TAIL           BRAK0050
C  AND VICE VERSA                                                       BRAK0060
C                                                                       BRAK0070
C.......................................................................BRAK0080
      DIMENSION X(1),Y(1),Z(1),YY(1)                                    BRAK0090
      NNN=N                                                             BRAK0100
      DO 2 I=1,NN                                                       BRAK0110
      YYY=YY(I)                                                         BRAK0120
      DO 1 J=1,N                                                        BRAK0130
      IF (Y(J).GE.YYY) GO TO 2                                          BRAK0140
      IF (Y(J).LT.YYY.AND.Y(J+1).GT.YYY) GO TO 3                        BRAK0150
      GO TO 1                                                           BRAK0160
    3 NNN=NNN+1                                                         BRAK0170
C                                                                       BRAK0180
C  SUBROUTINE NEWBRK DETERMINES THE BREAK LINE COORDINATES              BRAK0190
C                                                                       BRAK0200
      CALL NEWBRK(X(J),X(J+1),Y(J),Y(J+1),Z(J),Z(J+1),YYY,X(NNN),Y(NNN),BRAK0210
     .Z(NNN))                                                           BRAK0220
    1 CONTINUE                                                          BRAK0230
    2 CONTINUE                                                          BRAK0240
      RETURN                                                            BRAK0250
      END                                                               BRAK0260
      SUBROUTINE DRGDTA                                                 DRGD0010
C                                                                       DRGD0020
C.......................................................................DRGD0030
C  DRGDTA ACCEPTS GEOMETRY DATA IN THE FORMAT OF NASA TM-X 2074 AND     DRGD0040
C  PANELS THE AIRPLANE IN CONFORMANCE WITH THIS PROGRAM.                DRGD0050
C                                                                       DRGD0060
C  THIS COMPUTER PROGRAM CAN REDUCE DATA ONLY FOR FOUR CASES            DRGD0070
C      (1) WING, FUSELAGE AND HORIZONTAL TAIL                           DRGD0080
C      (2) WING AND FUSELAGE                                            DRGD0090
C      (3) WING AND HORIZONTAL TAIL                                     DRGD0100
C      (4) WING ONLY                                                    DRGD0110
C                                                                       DRGD0120
C.......................................................................DRGD0130
      COMMON/CNTROL/ LPANEL,NEAMAX,NCNTL,HALFSW,CREF,AM,CLDES,          DRGD0140
     . NCP2,NCP3,NCP4,WNGPNL,HT1PNL,HT2PNL,NCPHT1,NCPHT2,NWNGP,         DRGD0150
     . NHT1P,NHT2P,HALFB,ALPH,HALFB1,HALFB2,NCPSWL,NCPWNG,MAX           DRGD0160
      COMMON XCG(300),XCP(300),YCP(300),ZCP(300),AREA(300),AMASS(300),  DRGD0170
     . XN(300,4),YN(300,4),ZN(300,4),BPRIM(300,4)                       DRGD0180
      COMMON/PNL/CPCWL(35),CPSWL(35),IPANEL                             DRGD0190
      DIMENSION XHT1LE(35),XHT1TE(35),YHT1(35),XHLE(2),YHLE(2),XHTE(2)  DRGD0200
     . ,XWLE1(2),XWTE1(2),YL1(2),XHT2LE(35),XHT2TE(35),YHT2(35)         DRGD0210
      DIMENSION X(20),Y(20),Z(20),CHORD(20),XFUS(120),YFUS(120),ZFUS(120DRGD0220
     .),NRADX(4),NFORX(4),DUMMY(30),ZHT(2,2),SLOPLE(19),SLOPTE(19),ZAVHTDRGD0230
     .(2),XHT(2,20),YHT(2,20),CHRDHT(2,20),XFL(2),XFT(2),YF(2),YL(30),XWDRGD0240
     .LE(30),XWTE(30),XFUS1(25),YFUS1(25),ZFUS1(25)                     DRGD0250
      INTEGER WNGPNL,HT1PNL,HT2PNL                                      DRGD0260
      PAI=3.141592654                                                   DRGD0270
C***********************************************************************DRGD0280
      READ (5,57)J0,J1,J2,J3,J4,J5,J6,NWAF,NWAFOR,NFUS,(NRADX(I),NFORX(IDRGD0290
     .),I=1,4),NP,NPODOR,NF,NFINOR,NCAN,NCANOR                          DRGD0300
C***********************************************************************DRGD0310
      IF (J1.NE.0.AND.J2.NE.0.AND.J5.NE.0) NOPTON=1                     DRGD0320
      IF (J1.NE.0.AND.J2.NE.0.AND.J5.EQ.0) NOPTON=2                     DRGD0330
      IF (J1.NE.0.AND.J2.EQ.0.AND.J5.NE.0) NOPTON=3                     DRGD0340
      IF (J1.NE.0.AND.J2.EQ.0.AND.J5.EQ.0) NOPTON=4                     DRGD0350
      IF (J1.EQ.0.AND.J2.NE.0.AND.J5.NE.0) GO TO 53                     DRGD0360
      IF (J1.EQ.0.AND.J2.NE.0.AND.J5.EQ.0) GO TO 54                     DRGD0370
      IF (J1.EQ.0.AND.J2.EQ.0.AND.J5.NE.0) GO TO 55                     DRGD0380
      IF (J1.EQ.0.AND.J2.EQ.0.AND.J5.EQ.0) GO TO 56                     DRGD0390
C                                                                       DRGD0400
C  J0  (REFERENCE AREA)                                                 DRGD0410
C                                                                       DRGD0420
C***********************************************************************DRGD0430
      IF (J0.EQ.1) READ (5,58) REFA                                     DRGD0440
C***********************************************************************DRGD0450
C                                                                       DRGD0460
C  J1  (WING DATA)                                                      DRGD0470
C                                                                       DRGD0480
      IF (J1.EQ.0) GO TO 6                                              DRGD0490
C***********************************************************************DRGD0500
      READ (5,58)(GARB,I=1,NWAFOR)                                      DRGD0510
      DO 1 I=1,NWAF                                                     DRGD0520
      READ (5,58)X(I),Y(I),Z(I),CHORD(I)                                DRGD0530
    1 CHORD(I)=CHORD(I)+X(I)                                            DRGD0540
      IF (J1.EQ.(-1)) GO TO 3                                           DRGD0550
      DO 2 I=1,NWAF                                                     DRGD0560
    2 READ (5,58)(GARB,K=1,NWAFOR)                                      DRGD0570
    3 DO 4 I=1,NWAF                                                     DRGD0580
    4 READ (5,58)(GARB,K=1,NWAFOR)                                      DRGD0590
C***********************************************************************DRGD0600
C                                                                       DRGD0610
C  ZAVWNG   (AVERAGE HEIGHT OF WING)                                    DRGD0620
C                                                                       DRGD0630
      SUM1=0.                                                           DRGD0640
      DO 5 I=1,NWAF                                                     DRGD0650
    5 SUM1=SUM1+Z(I)                                                    DRGD0660
      ZAVWNG=SUM1/(FLOAT(NWAF))                                         DRGD0670
      NEWAF=NWAF                                                        DRGD0680
C                                                                       DRGD0690
C  J2  (FUSELAGE DATA)                                                  DRGD0700
C  J6  (FUSELAGE DATA)                                                  DRGD0710
C                                                                       DRGD0720
      FUSWTH=0.                                                         DRGD0730
    6 IF (J2.EQ.0) GO TO 20                                             DRGD0740
C***********************************************************************DRGD0750
      NFUSXI=1                                                          DRGD0760
      DO 8 I=1,NFUS                                                     DRGD0770
      NRAD=NRADX(I)                                                     DRGD0780
      NFUSXF=NFUSXI+NFORX(I)-1                                          DRGD0790
      READ (5,58)(XFUS(K),K=NFUSXI,NFUSXF)                              DRGD0800
      IF (J2.EQ.1) GO TO 9                                              DRGD0810
      GO TO 100                                                         DRGD0820
    9 DO 10 K=NFUSXI,NFUSXF                                             DRGD0830
      READ (5,58)(DUMMY(L),L=1,NRAD)                                    DRGD0840
C                                                                       DRGD0850
C  SUBROUTINE  MAXVLU  FINDS THE VALUE OF THE LARGEST ELEMENT           DRGD0860
C  OF ARRAY DUMMY                                                       DRGD0870
C                                                                       DRGD0880
      CALL MAXVLU(DUMMY,NRAD,YMAX,KK)                                   DRGD0890
      YFUS(K)=YMAX                                                      DRGD0900
      READ (5,58)(DUMMY(L),L=1,NRAD)                                    DRGD0910
      ZFUS(K)=DUMMY(KK)                                                 DRGD0920
   10 CONTINUE                                                          DRGD0930
      GO TO 8                                                           DRGD0940
  100 IF (J6.EQ.0) READ (5,58)(ZFUS(K),K=NFUSXI,NFUSXF)                 DRGD0950
      READ (5,58)(YFUS(K),K=NFUSXI,NFUSXF)                              DRGD0960
      DO 7 K=NFUSXI,NFUSXF                                              DRGD0970
    7 YFUS(K)=SQRT(YFUS(K)/PAI)                                         DRGD0980
    8 NFUSXI=NFUSXF                                                     DRGD0990
C***********************************************************************DRGD1000
C                                                                       DRGD1010
C  ZAVFUS   (AVERAGE HEIGHT OF FUSELAGE)                                DRGD1020
C                                                                       DRGD1030
      SUM2=0.                                                           DRGD1040
      DO 12 K=1,NFUSXF                                                  DRGD1050
   12 SUM2=SUM2+ZFUS(K)                                                 DRGD1060
      ZAVFUS=SUM2/FLOAT(NFUSXF)                                         DRGD1070
      XLE=X(1)                                                          DRGD1080
      XTE=CHORD(1)                                                      DRGD1090
C                                                                       DRGD1100
C  FIND WIDTH OF FUSELAGE                                               DRGD1110
C                                                                       DRGD1120
C                                                                       DRGD1130
C  SUBROUTINE WIDTH FINDS OUT THE AVERAGE FUSELAGE WIDTH BETWEEN        DRGD1140
C  THE GIVEN TWO STATIONS                                               DRGD1150
C                                                                       DRGD1160
      CALL WIDTH(XFUS,YFUS,NFUSXF,XLE,XTE,FUSWTH)                       DRGD1170
      RAD2=FUSWTH/2.                                                    DRGD1180
      DIFF=ABS(ZAVWNG-ZAVFUS)                                           DRGD1190
      IF (DIFF.GT.RAD2) GO TO 13                                        DRGD1200
C                                                                       DRGD1210
C  WING AND FUSELAGE ARE IN THE SAME PLANE                              DRGD1220
C                                                                       DRGD1230
      ZAVWNG=ZAVFUS                                                     DRGD1240
      FUSWTH=Y(1)                                                       DRGD1250
      GO TO 14                                                          DRGD1260
C                                                                       DRGD1270
C  CHECK FOR WING INTERSECTION WITH CENTER LINE                         DRGD1280
C                                                                       DRGD1290
   13 NEWAF=NEWAF+1                                                     DRGD1300
      Y(NEWAF)=0.                                                       DRGD1310
      X(NEWAF)=X(1)-Y(1)*SIN(SLOPLE(1))/COS(SLOPLE(1))                  DRGD1320
      CHORD(NEWAF)=CHORD(1)-Y(1)*SIN(SLOPTE(1))/COS(SLOPTE(1))          DRGD1330
      X(1)=X(1)-(Y(1)-FUSWTH)*SIN(SLOPLE(1))/COS(SLOPLE(1))             DRGD1340
      CHORD(1)=CHORD(1)-(Y(1)-FUSWTH)*SIN(SLOPTE(1))/COS(SLOPTE(1))     DRGD1350
      Y(1)=FUSWTH                                                       DRGD1360
   14 CONTINUE                                                          DRGD1370
      NCOUNT=1                                                          DRGD1380
      NFUSX1=NFUSXF-1                                                   DRGD1390
C                                                                       DRGD1400
C  FRONT AND REAR FUSELAGE MODIFICATION (STARTS)                        DRGD1410
C                                                                       DRGD1420
      AREAR=0.                                                          DRGD1430
      AREAF=0.                                                          DRGD1440
      DO 19 I=1,NFUSX1                                                  DRGD1450
      IF (XFUS(I).LT.XLE) GO TO 15                                      DRGD1460
      IF (XFUS(I).GE.XTE) GO TO 16                                      DRGD1470
      GO TO 19                                                          DRGD1480
   15 AREAF=AREAF+(YFUS(I)+YFUS(I+1))*(XFUS(I+1)-XFUS(I))/2.            DRGD1490
      XF=XFUS(I+1)-XFUS(1)                                              DRGD1500
      GO TO 19                                                          DRGD1510
   16 AREAR=AREAR+(YFUS(I)+YFUS(I+1))*(XFUS(I+1)-XFUS(I))/2.            DRGD1520
      IF (NCOUNT.EQ.1) GO TO 17                                         DRGD1530
      GO TO 19                                                          DRGD1540
   17 IF (XTE.GT.XFUS(I-1).AND.XTE.LT.XFUS(I)) GO TO 18                 DRGD1550
      XR=XFUS(I)                                                        DRGD1560
      NCOUNT=2                                                          DRGD1570
      GO TO 19                                                          DRGD1580
   18 AREAR=AREAR+(YFUS(I-1)+YFUS(I))*(XFUS(I)-XFUS(I-1))/2.            DRGD1590
      XR=XFUS(I-1)                                                      DRGD1600
      NCOUNT=2                                                          DRGD1610
   19 CONTINUE                                                          DRGD1620
      A=2.*XF-2.*AREAF/FUSWTH                                           DRGD1630
      XFL(2)=XFUS(1)+A                                                  DRGD1640
      XX=XFUS(NFUSXF)-XR                                                DRGD1650
      A=2.*XX-2.*AREAR/FUSWTH                                           DRGD1660
      XFT(2)=XFUS(NFUSXF)-A                                             DRGD1670
      XFL(1)=XFUS(1)                                                    DRGD1680
      XFT(1)=XFUS(NFUSXF)                                               DRGD1690
      YF(1)=0.                                                          DRGD1700
      YF(2)=FUSWTH                                                      DRGD1710
C                                                                       DRGD1720
C  FRONT AND REAR FUSELAGE MODIFICATION (ENDS)                          DRGD1730
C                                                                       DRGD1740
C                                                                       DRGD1750
C  J3  (POD DATA)                                                       DRGD1760
C                                                                       DRGD1770
   20 IF (J3.EQ.0) GO TO 22                                             DRGD1780
C***********************************************************************DRGD1790
      DO 21 K=1,NP                                                      DRGD1800
      READ (5,58)(GARB,I=1,3)                                           DRGD1810
      READ (5,58)(GARB,I=1,NPODOR)                                      DRGD1820
   21 READ (5,58)(GARB,I=1,NPODOR)                                      DRGD1830
C***********************************************************************DRGD1840
C                                                                       DRGD1850
C  J4  (FIN DATA)                                                       DRGD1860
C                                                                       DRGD1870
   22 IF (J4.EQ.0) GO TO 24                                             DRGD1880
C***********************************************************************DRGD1890
      DO 23 K=1,NF                                                      DRGD1900
      READ (5,58)(GARB,I=1,8)                                           DRGD1910
      READ (5,58)(GARB,I=1,NFINOR)                                      DRGD1920
   23 READ (5,58)(GARB,I=1,NFINOR)                                      DRGD1930
C***********************************************************************DRGD1940
C                                                                       DRGD1950
C  J5  (CANARD DATA)                                                    DRGD1960
C                                                                       DRGD1970
   24 IF (J5.EQ.0) GO TO 37                                             DRGD1980
C***********************************************************************DRGD1990
      DO 26 K=1,NCAN                                                    DRGD2000
      READ (5,58)(XHT(K,I),YHT(K,I),ZHT(K,I),CHRDHT(K,I),I=1,2)         DRGD2010
      CHRDHT(K,1)=CHRDHT(K,1)+XHT(K,1)                                  DRGD2020
      CHRDHT(K,2)=CHRDHT(K,2)+XHT(K,2)                                  DRGD2030
      IF (NCANOR.LT.0) GO TO 25                                         DRGD2040
      READ (5,58)(GARB,I=1,NCANOR)                                      DRGD2050
      READ (5,58)(GARB,I=1,NCANOR)                                      DRGD2060
      GO TO 26                                                          DRGD2070
   25 NCANO1=-NCANOR                                                    DRGD2080
      READ (5,58)(GARB,I=1,NCANO1)                                      DRGD2090
      READ (5,58)(GARB,I=1,NCANO1)                                      DRGD2100
      READ (5,58)(GARB,I=1,NCANO1)                                      DRGD2110
   26 CONTINUE                                                          DRGD2120
C***********************************************************************DRGD2130
      WIDTH1=0.                                                         DRGD2140
      WIDTH2=0.                                                         DRGD2150
      IF (J2.EQ.0) GO TO 27                                             DRGD2160
      XLE=XHT(1,1)                                                      DRGD2170
      XTE=CHRDHT(1,1)                                                   DRGD2180
      CALL WIDTH(XFUS,YFUS,NFUSXF,XLE,XTE,WDTH1)                        DRGD2190
      IF (NCAN.EQ.1) GO TO 27                                           DRGD2200
      XLE=XHT(2,1)                                                      DRGD2210
      XTE=CHRDHT(2,1)                                                   DRGD2220
      CALL WIDTH(XFUS,YFUS,NFUSXF,XLE,XTE,WDTH2)                        DRGD2230
   27 CONTINUE                                                          DRGD2240
      MAA=2                                                             DRGD2250
      DO 28 I=1,2                                                       DRGD2260
      XFUS(I)=XHT(1,I)                                                  DRGD2270
      YFUS(I)=YHT(1,I)                                                  DRGD2280
   28 ZFUS(I)=CHRDHT(1,I)                                               DRGD2290
      ZAVHT(1)=(ZHT(1,1)+ZHT(1,2))/2.                                   DRGD2300
      IF (NCAN.EQ.1) GO TO 29                                           DRGD2310
      NCONTL=0                                                          DRGD2320
      IF (XHT(1,2).EQ.XHT(2,1).AND.YHT(1,2).EQ.YHT(2,1).AND.ZHT(1,2).EQ.DRGD2330
     .ZHT(2,1).AND.CHRDHT(1,2).EQ.CHRDHT(2,1)) NCONTL=1                 DRGD2340
      IF (NCONTL.EQ.0) GO TO 30                                         DRGD2350
   29 XFUS(1)=XFUS(1)-(XFUS(2)-XFUS(1))*(YFUS(1)-FUSWTH)/(YFUS(2)-YFUS(1DRGD2360
     .))                                                                DRGD2370
      ZFUS(1)=ZFUS(1)-(ZFUS(2)-ZFUS(1))*(YFUS(1)-FUSWTH)/(YFUS(2)-YFUS(1DRGD2380
     .))                                                                DRGD2390
      YFUS(1)=FUSWTH                                                    DRGD2400
      IF (NCAN.EQ.1) GO TO 33                                           DRGD2410
      XFUS(3)=XHT(2,2)                                                  DRGD2420
      YFUS(3)=YHT(2,2)                                                  DRGD2430
      ZFUS(3)=CHRDHT(2,2)                                               DRGD2440
      ZAVHT(1)=(ZHT(1,1)+ZHT(1,2)+ZHT(2,2))/3.                          DRGD2450
      NCAN=1                                                            DRGD2460
      MAA=3                                                             DRGD2470
      GO TO 33                                                          DRGD2480
   30 DO 31 I=1,2                                                       DRGD2490
      XFUS1(I)=XHT(2,I)                                                 DRGD2500
      YFUS1(I)=YHT(2,I)                                                 DRGD2510
   31 ZFUS1(I)=CHRDHT(2,I)                                              DRGD2520
      XFUS1(1)=XFUS1(1)-(XFUS1(2)-XFUS1(1))*(YFUS1(1)-FUSWTH)/(YFUS1(2)-DRGD2530
     .YFUS1(1))                                                         DRGD2540
      ZFUS1(1)=ZFUS1(1)-(ZFUS1(2)-ZFUS1(1))*(YFUS1(1)-FUSWTH)/(YFUS1(2)-DRGD2550
     .YFUS1(1))                                                         DRGD2560
      YFUS1(1)=FUSWTH                                                   DRGD2570
      ZAVHT(2)=(ZHT(2,1)+ZHT(2,2))/2.                                   DRGD2580
C                                                                       DRGD2590
C  SUBROUTINE BREAK LOCATES BREAK LINES ON ONE LIFTING                  DRGD2600
C  SURFACE DUE TO THE BREAK LINES ON THE OTHER LIFTING SURFACE          DRGD2610
C                                                                       DRGD2620
C                                                                       DRGD2630
C  BREAK LINES ON TAIL 1 DUE TO TAIL 2                                  DRGD2640
C                                                                       DRGD2650
      CALL BREAK(XFUS,YFUS,ZFUS,2,YFUS1,2,NN1)                          DRGD2660
C                                                                       DRGD2670
C  BREAK LINES ON TAIL 2 DUE TO TAIL 1                                  DRGD2680
C                                                                       DRGD2690
      CALL BREAK(XFUS1,YFUS1,ZFUS1,2,YFUS,NN1,NN2)                      DRGD2700
C                                                                       DRGD2710
C  BREAK LINES ON WING DUE TO TAIL 1                                    DRGD2720
C                                                                       DRGD2730
      CALL BREAK(X,Y,CHORD,NEWAF,YFUS,NN1,NEWAF1)                       DRGD2740
C                                                                       DRGD2750
C  BREAK LINES ON WING DUE TO TAIL 2                                    DRGD2760
C                                                                       DRGD2770
      CALL BREAK(X,Y,CHORD,NEWAF1,YFUS1,NN2,NEWAF2)                     DRGD2780
C                                                                       DRGD2790
C  BREAK LINES ON TAIL 1 DUE TO WING                                    DRGD2800
C                                                                       DRGD2810
      CALL BREAK(XFUS,YFUS,ZFUS,NN1,Y,NEWAF2,NN3)                       DRGD2820
C                                                                       DRGD2830
C  BREAK LINES ON TAIL 2 DUE TO WING                                    DRGD2840
C                                                                       DRGD2850
      CALL BREAK(XFUS1,YFUS1,ZFUS1,NN2,Y,NEWAF2,NN4)                    DRGD2860
      NN2=NN3                                                           DRGD2870
      NEWAF=NEWAF2                                                      DRGD2880
      WDTH2=WDTH2/2.                                                    DRGD2890
      DIFF=ABS(ZAVFUS-ZAVHT(2))                                         DRGD2900
      IF (DIFF.LT.WDTH2) GO TO 32                                       DRGD2910
      NN4=NN4+1                                                         DRGD2920
      YFUS1(NN4)=0.                                                     DRGD2930
      XFUS1(NN4)=XFUS1(1)-(XFUS1(2)-XFUS1(1))*YFUS1(1)/(YFUS1(2)-YFUS1(1DRGD2940
     .))                                                                DRGD2950
      ZFUS1(NN4)=ZFUS1(1)-(ZFUS1(2)-ZFUS1(1))*YFUS1(1)/(YFUS1(2)-YFUS1(1DRGD2960
     .))                                                                DRGD2970
      NN4=NN4+1                                                         DRGD2980
      YFUS1(NN4)=FUSWTH                                                 DRGD2990
      YFUS11=YFUS1(1)-FUSWTH                                            DRGD3000
      XFUS1(NN4)=XFUS1(1)-(XFUS1(2)-XFUS1(1))*YFUS11/(YFUS1(2)-YFUS1(1))DRGD3010
      ZFUS1(NN4)=ZFUS1(1)-(ZFUS1(2)-ZFUS1(1))*YFUS11/(YFUS1(2)-YFUS1(1))DRGD3020
      NHT2=NN4                                                          DRGD3030
      GO TO 34                                                          DRGD3040
   32 IF (J2.NE.0) ZAVHT(2)=ZAVFUS                                      DRGD3050
      GO TO 34                                                          DRGD3060
   33 CONTINUE                                                          DRGD3070
C                                                                       DRGD3080
C  BREAK LINES ON WING DUE TO TAIL                                      DRGD3090
C                                                                       DRGD3100
      CALL BREAK(X,Y,CHORD,NEWAF,YFUS,MAA,NEWAF2)                       DRGD3110
C                                                                       DRGD3120
C  BREAK LINES ON TAIL DUE TO WING                                      DRGD3130
C                                                                       DRGD3140
      CALL BREAK(XFUS,YFUS,ZFUS,MAA,Y,NEWAF2,NN2)                       DRGD3150
      NEWAF=NEWAF2                                                      DRGD3160
   34 CONTINUE                                                          DRGD3170
      WDTH1=WDTH1/2.                                                    DRGD3180
      DIFF=ABS(ZAVFUS-ZAVHT(1))                                         DRGD3190
      IF (DIFF.LT.WDTH1) GO TO 35                                       DRGD3200
      NN2=NN2+1                                                         DRGD3210
      YFUS(NN2)=0.                                                      DRGD3220
      XFUS(NN2)=XFUS(1)-(XFUS(2)-XFUS(1))*YFUS(1)/(YFUS(2)-YFUS(1))     DRGD3230
      ZFUS(NN2)=ZFUS(1)-(ZFUS(2)-ZFUS(1))*YFUS(1)/(YFUS(2)-YFUS(1))     DRGD3240
      NN2=NN2+1                                                         DRGD3250
      YFUS(NN2)=FUSWTH                                                  DRGD3260
      YFUS11=YFUS(1)-FUSWTH                                             DRGD3270
      XFUS(NN2)=XFUS(1)-(XFUS(2)-XFUS(1))*YFUS11/(YFUS(2)-YFUS(1))      DRGD3280
      ZFUS(NN2)=ZFUS(1)-(ZFUS(2)-ZFUS(1))*YFUS11/(YFUS(2)-YFUS(1))      DRGD3290
      NHT1=NN2                                                          DRGD3300
      GO TO 36                                                          DRGD3310
   35 IF (J2.NE.0) ZAVHT(1)=ZAVFUS                                      DRGD3320
   36 CONTINUE                                                          DRGD3330
C                                                                       DRGD3340
C  SUBROUTINE  ORDER  ARRANGES ALL THE BREAK LINES IN ORDER             DRGD3350
C  FROM INBOARD TO OUTBOARD OF WING                                     DRGD3360
C                                                                       DRGD3370
C  COMPARE FOR ALL THE OPTIONS                                          DRGD3380
   37 CONTINUE                                                          DRGD3390
      CALL ORDER(X,Y,CHORD,NEWAF,XWLE,XWTE,YL,IWNG)                     DRGD3400
      IF (NOPTON.EQ.2.OR.NOPTON.EQ.4) GO TO 38                          DRGD3410
      IF (NCAN.NE.0) CALL ORDER(XFUS,YFUS,ZFUS,NN2,XHT1LE,XHT1TE,YHT1,IHDRGD3420
     .T1)                                                               DRGD3430
      IF (NCAN.EQ.2) CALL ORDER(XFUS1,YFUS1,ZFUS1,NN4,XHT2LE,XHT2TE,YHT2DRGD3440
     .,IHT2)                                                            DRGD3450
   38 CONTINUE                                                          DRGD3460
      NEWAF=IWNG                                                        DRGD3470
      NC=NEWAF-1                                                        DRGD3480
      IF (J2.NE.0) NC=NC+1                                              DRGD3490
      IF (IHT1.NE.0) NC=NC+IHT1-1                                       DRGD3500
      IF (IHT2.NE.0) NC=NC+IHT2-1                                       DRGD3510
      WRITE (6,65) NC                                                   DRGD3520
      KKK=0                                                             DRGD3530
      IPANEL=1                                                          DRGD3540
C***********************************************************************DRGD3550
      READ (5,59) NCP1,NCP2,NCP3,NCP4                                   DRGD3560
C***********************************************************************DRGD3570
      NPNLF=NCP1                                                        DRGD3580
      NPNLW=NCP2                                                        DRGD3590
      NPNLT1=NCP3                                                       DRGD3600
      NPNLT2=NCP4                                                       DRGD3610
      NWNGP=0                                                           DRGD3620
      NHT1P=0                                                           DRGD3630
      NHT2P=0                                                           DRGD3640
      NPNLF1=NPNLF                                                      DRGD3650
      NPNLW1=NPNLW                                                      DRGD3660
      NPNL11=NPNLT1                                                     DRGD3670
      NPNL22=NPNLT2                                                     DRGD3680
      IF (NOPTON.GT.2) GO TO 39                                         DRGD3690
C                                                                       DRGD3700
C  FUSELAGE ANALYSIS                                                    DRGD3710
C                                                                       DRGD3720
C***********************************************************************DRGD3730
      READ (5,60) (CPCWL(I),I=1,NPNLF1)                                 DRGD3740
C***********************************************************************DRGD3750
      CPSWL(1)=0.                                                       DRGD3760
      CPSWL(2)=100.                                                     DRGD3770
C                                                                       DRGD3780
C  SUBROUTINE  RITXYZ  WRITES THE MODIFIED GEOMETRY DATA                DRGD3790
C                                                                       DRGD3800
      CALL RITXYZ(KKK,XFL,XFT,YF,CPCWL,NPNLF1,CPSWL,2,ZAVFUS)           DRGD3810
C                                                                       DRGD3820
C  SUBROUTINE PANEL DIVIDES A LIFTING SURFACE INTO AERODYNAMIC          DRGD3830
C  PANELS                                                               DRGD3840
C                                                                       DRGD3850
      CALL PANEL(XFL,YF,XFT,NPNLF1,2,ZAVFUS)                            DRGD3860
      IPANEL=LPANEL+1                                                   DRGD3870
C                                                                       DRGD3880
C  WING ANALYSIS                                                        DRGD3890
C                                                                       DRGD3900
      CREF=XWTE(1)-XWLE(1)                                              DRGD3910
   39 IST=1                                                             DRGD3920
      WNGPNL=IPANEL                                                     DRGD3930
      NEWAF1=NEWAF-1                                                    DRGD3940
C***********************************************************************DRGD3950
      READ (5,60) (CPCWL(I),I=1,NPNLW1)                                 DRGD3960
C***********************************************************************DRGD3970
      HALFB=YL(NEWAF)                                                   DRGD3980
      W1=(YL(NEWAF)-YL(1))/10.                                          DRGD3990
      IF (ZAVFUS.NE.ZAVWNG) GO TO 42                                    DRGD4000
   40 DO 41 I=IST,NEWAF1                                                DRGD4010
C                                                                       DRGD4020
C  SUBROUTINE  NEWCOR  OUTPUTS THE INBOARD AND OUTBOARD CHORD           DRGD4030
C  OF LIFTING SURFACE UNDER CONSIDERATION                               DRGD4040
C                                                                       DRGD4050
      CALL NEWCOR(XWLE,XWTE,YL,I,XWLE1,XWTE1,YL1)                       DRGD4060
      INT=(YL(I+1)-YL(I))/W1+0.5                                        DRGD4070
      IF (INT.EQ.0) INT=1                                               DRGD4080
C                                                                       DRGD4090
C  SUBROUTINE NEWCOR LOCATES THE CONSTANT PERCENT STREAMWISE            DRGD4100
C  LINES (CPSWL) FOR THE LIFTING SURFACE                                DRGD4110
C                                                                       DRGD4120
      CALL GNCCWL(INT,CPSWL,1.)                                         DRGD4130
      INT1=INT+1                                                        DRGD4140
      CALL RITXYZ(KKK,XWLE1,XWTE1,YL1,CPCWL,NPNLW1,CPSWL,INT1,ZAVWNG)   DRGD4150
      CALL PANEL(XWLE1,YL1,XWTE1,NPNLW1,INT1,ZAVWNG)                    DRGD4160
      IPANEL=LPANEL+1                                                   DRGD4170
      NWNGP=NWNGP+INT1-1                                                DRGD4180
   41 CONTINUE                                                          DRGD4190
      GO TO 43                                                          DRGD4200
   42 CPSWL(1)=0.                                                       DRGD4210
      CPSWL(2)=100.                                                     DRGD4220
      CALL NEWCOR(XWLE,XWTE,YL,1,XWLE1,XWTE1,YL1)                       DRGD4230
      CALL RITXYZ(KKK,XWLE1,XWTE1,YL1,CPCWL,NPNLW1,CPSWL,2,ZAVWNG)      DRGD4240
      CALL PANEL(XWLE1,YL1,XWTE1,NPNLW1,2,ZAVWNG)                       DRGD4250
      IPANEL=LPANEL+1                                                   DRGD4260
      IST=2                                                             DRGD4270
      NWNGP=1                                                           DRGD4280
      GO TO 40                                                          DRGD4290
   43 CONTINUE                                                          DRGD4300
      IF (NOPTON.EQ.2.OR.NOPTON.EQ.4) GO TO 48                          DRGD4310
C                                                                       DRGD4320
C  HORIZONTAL TAIL ANALYSIS                                             DRGD4330
C                                                                       DRGD4340
      IST=1                                                             DRGD4350
      HALFB1=YHT1(IHT1)                                                 DRGD4360
      IHT11=IHT1-1                                                      DRGD4370
      HT1PNL=IPANEL                                                     DRGD4380
C***********************************************************************DRGD4390
      READ (5,60) (CPCWL(I),I=1,NPNL11)                                 DRGD4400
C***********************************************************************DRGD4410
      IF (NOPTON.LE.2.AND.ZAVFUS.NE.ZAVHT(1)) GO TO 51                  DRGD4420
   44 DO 45 I=IST,IHT11                                                 DRGD4430
      INT=(YHT1(I+1)-YHT1(I))/W1+0.5                                    DRGD4440
      CALL NEWCOR(XHT1LE,XHT1TE,YHT1,I,XHLE,XHTE,YHLE)                  DRGD4450
      IF (INT.EQ.0) INT=1                                               DRGD4460
      INT1=INT+1                                                        DRGD4470
      CALL GNCCWL(INT,CPSWL,1.)                                         DRGD4480
      CALL RITXYZ(KKK,XHLE,XHTE,YHLE,CPCWL,NPNL11,CPSWL,INT1,ZAVHT(1))  DRGD4490
      CALL PANEL(XHLE,YHLE,XHTE,NPNL11,INT1,ZAVHT(1))                   DRGD4500
      NHT1P=NHT1P+INT1-1                                                DRGD4510
   45 IPANEL=LPANEL+1                                                   DRGD4520
      IF (NCAN.NE.2) GO TO 48                                           DRGD4530
      IST=1                                                             DRGD4540
      HALFB2=YHT2(IHT2)                                                 DRGD4550
      HT2PNL=IPANEL                                                     DRGD4560
C***********************************************************************DRGD4570
      READ (5,60) (CPCWL(I),I=1,NPNL22)                                 DRGD4580
C***********************************************************************DRGD4590
      IHT21=IHT2-1                                                      DRGD4600
      IF (NOPTON.LE.2.AND.ZAVFUS.NE.ZAVHT(2)) GO TO 52                  DRGD4610
   46 DO 47 I=IST,IHT21                                                 DRGD4620
      INT=(YHT2(I+1)-YHT2(I))/W1+0.5                                    DRGD4630
      CALL NEWCOR(XHT2LE,XHT2TE,YHT2,I,XHLE,XHTE,YHLE)                  DRGD4640
      IF (INT.EQ.0) INT=1                                               DRGD4650
      INT1=INT+1                                                        DRGD4660
      CALL GNCCWL(INT,CPSWL,1.)                                         DRGD4670
      CALL RITXYZ(KKK,XHLE,XHTE,YHLE,CPCWL,NPNL22,CPSWL,INT1,ZAVHT(2))  DRGD4680
      CALL PANEL(XHLE,YHLE,XHTE,NPNL22,INT1,ZAVHT(2))                   DRGD4690
      NHT2P=NHT2P+INT1-1                                                DRGD4700
   47 IPANEL=LPANEL+1                                                   DRGD4710
   48 CONTINUE                                                          DRGD4720
      LSTPNL=LPANEL                                                     DRGD4730
C***********************************************************************DRGD4740
      READ (5,60) CCREF                                                 DRGD4750
C***********************************************************************DRGD4760
      IF (CCREF.NE.0.) CREF=CCREF                                       DRGD4770
      IF (REFA.NE.0.) GO TO 50                                          DRGD4780
      IM=HT1PNL-1                                                       DRGD4790
      IF (HT1PNL.LE.0) IM=LPANEL                                        DRGD4800
      HALFSW=0.                                                         DRGD4810
      DO 49 I=WNGPNL,IM                                                 DRGD4820
   49 HALFSW=HALFSW+AREA(I)                                             DRGD4830
      REFA=2.*HALFSW                                                    DRGD4840
      WRITE (6,66) REFA,CREF                                            DRGD4850
      RETURN                                                            DRGD4860
   50 HALFSW=REFA/2.                                                    DRGD4870
      WRITE (6,66) REFA,CREF                                            DRGD4880
      RETURN                                                            DRGD4890
   51 CPSWL(1)=0.                                                       DRGD4900
      CPSWL(2)=100.                                                     DRGD4910
      CALL NEWCOR(XHT1LE,XHT1TE,YHT1,1,XHLE,XHTE,YHLE)                  DRGD4920
      CALL RITXYZ(KKK,XHLE,XHTE,YHLE,CPCWL,NPNL11,CPSWL,2,ZAVHT(1))     DRGD4930
      CALL PANEL(XHLE,YHLE,XHTE,NPNL11,2,ZAVHT(1))                      DRGD4940
      IST=2                                                             DRGD4950
      NHT1P=1                                                           DRGD4960
      GO TO 44                                                          DRGD4970
   52 CPSWL(1)=0.                                                       DRGD4980
      CPSWL(2)=100.                                                     DRGD4990
      CALL NEWCOR(XHT2LE,XHT2TE,YHT2,1,XHLE,XHTE,YHLE)                  DRGD5000
      CALL RITXYZ(KKK,XHLE,XHTE,YHLE,CPCWL,NPNL22,CPSWL,2,ZAVHT(2))     DRGD5010
      CALL PANEL(XHLE,YHLE,XHTE,NPNL22,2,ZAVHT(2))                      DRGD5020
      IST=2                                                             DRGD5030
      NHT2P=1                                                           DRGD5040
      GO TO 46                                                          DRGD5050
   53 WRITE (6,61)                                                      DRGD5060
      RETURN                                                            DRGD5070
   54 WRITE (6,62)                                                      DRGD5080
      RETURN                                                            DRGD5090
   55 WRITE (6,63)                                                      DRGD5100
      RETURN                                                            DRGD5110
   56 WRITE (6,64)                                                      DRGD5120
      RETURN                                                            DRGD5130
   57 FORMAT (24I3)                                                     DRGD5140
   58 FORMAT (10F7.3)                                                   DRGD5150
   59 FORMAT (4I2)                                                      DRGD5160
   60 FORMAT (8F10.5)                                                   DRGD5170
   61 FORMAT (90HTHIS TEST CASE HAS FUSELAGE AND HORIZONTAL TAIL. CHANGEDRGD5180
     . HORIZONTAL TAIL DATA TO WING DATA.   )                           DRGD5190
   62 FORMAT (89HTHIS TEST CASE HAS ONLY FUSELAGE WHICH DOES NOT MAKE ANDRGD5200
     .Y SENSE. FIRST DATA CARD IS WRONG.    )                           DRGD5210
   63 FORMAT (82HTHIS TEST CASE HAS ONLY HORIZONTAL TAIL. CHANGE HORIZONDRGD5220
     .TAL TAIL DATA TO WING DATA.    )                                  DRGD5230
   64 FORMAT (79HTHIS TEST CASE DOES NOT HAVE ANY AERODYNAMIC SURFACE. FDRGD5240
     .IRST DATA CARD IS WORNG.    )                                     DRGD5250
   65 FORMAT (1H1,//,11X,45H********************************************DRGD5260
     .*,/,11X,32H* GEOMETRY DATA IS DIVIDED IN TO,I2,11H SECTIONS *,/,11DRGD5270
     .X,45H*********************************************,///)           DRGD5280
   66 FORMAT (10X,10HREF. AREA=,E12.5,/,10X,11HREF. CHORD=,E12.5)       DRGD5290
      END                                                               DRGD5300
      SUBROUTINE GNCCWL(NPNL,CPCWL,A)                                   GNCC0010
C                                                                       GNCC0020
C.......................................................................GNCC0030
C                                                                       GNCC0040
C  GNCCWL LOCATES THE CONSTANT PERCENT STREAMWISE LINES                 GNCC0050
C  FOR THE LIFTING SURFACE                                              GNCC0060
C                                                                       GNCC0070
C.......................................................................GNCC0080
      DIMENSION CPCWL(1)                                                GNCC0090
      CPCWLL=100./FLOAT(NPNL)                                           GNCC0100
      CPCWL(1)=0.                                                       GNCC0110
      CPCWL(2)=CPCWLL/A                                                 GNCC0120
      CPCWL(NPNL+1)=100.                                                GNCC0130
      DO 1 I=3,NPNL                                                     GNCC0140
    1 CPCWL(I)=CPCWL(I-1)+CPCWLL                                        GNCC0150
      RETURN                                                            GNCC0160
      END                                                               GNCC0170
      SUBROUTINE MAXVLU(Y,N,YMAX,KK)                                    MAXV0010
C                                                                       MAXV0020
C.......................................................................MAXV0030
C                                                                       MAXV0040
C  MAXVLU  RETURNS THE NUMBER AND VALUE OF THE LARGEST                  MAXV0050
C  ELEMENT OF ARRAY Y                                                   MAXV0060
C                                                                       MAXV0070
C.......................................................................MAXV0080
      DIMENSION Y(1)                                                    MAXV0090
      I=1                                                               MAXV0100
      NN=N+1                                                            MAXV0110
    1 IF (Y(I).GT.Y(I+1)) GO TO 2                                       MAXV0120
      I=I+1                                                             MAXV0130
      IF (I.EQ.NN) GO TO 2                                              MAXV0140
      GO TO 1                                                           MAXV0150
    2 YMAX=Y(I)                                                         MAXV0160
      KK=I                                                              MAXV0170
      RETURN                                                            MAXV0180
      END                                                               MAXV0190
      SUBROUTINE NEWBRK(X1,X2,Y1,Y2,C1,C2,YY,XNEW,YNEW,CNEW)            NBRK0010
C                                                                       NBRK0020
C.....                                                                  NBRK0030
C                                                                       NBRK0040
C  NEWBRK COMPUTES THE BREAK LINE COORDINATES                           NBRK0050
C                                                                       NBRK0060
C.......................................................................NBRK0070
      YY1=(YY-Y1)/(Y2-Y1)                                               NBRK0080
      X21=X2-X1                                                         NBRK0090
      XNEW=X1+X21*YY1                                                   NBRK0100
      C21=C2-C1                                                         NBRK0110
      CNEW=C1+C21*YY1                                                   NBRK0120
      YNEW=YY                                                           NBRK0130
      RETURN                                                            NBRK0140
      END                                                               NBRK0150
      SUBROUTINE NEWCOR(XWLE,XWTE,YL,I,XWLE1,XWTE1,YL1)                 NCOR0010
C.......................................................................NCOR0020
C                                                                       NCOR0030
C  NEWCOR RETURNS THE INBOARD AND OUTBOARD COORDINATE OF THE            NCOR0040
C  LIFTING SURFACE UNDER CONSIDERATION                                  NCOR0050
C                                                                       NCOR0060
C.......................................................................NCOR0070
      DIMENSION XWLE(1),XWTE(1),YL(1),XWLE1(1),XWTE1(1),YL1(1)          NCOR0080
      YL1(1)=YL(I)                                                      NCOR0090
      YL1(2)=YL(I+1)                                                    NCOR0100
      XWLE1(1)=XWLE(I)                                                  NCOR0110
      XWLE1(2)=XWLE(I+1)                                                NCOR0120
      XWTE1(1)=XWTE(I)                                                  NCOR0130
      XWTE1(2)=XWTE(I+1)                                                NCOR0140
      RETURN                                                            NCOR0150
      END                                                               NCOR0160
      SUBROUTINE NEWPNL                                                 NPNL0010
C                                                                       NPNL0020
C.......................................................................NPNL0030
C                                                                       NPNL0040
C  NEWPNL PANELS THE AIRPLANE USING THE INPUT FORMATS OF NASA-CR-112229.NPNL0050
C                                                                       NPNL0060
C  KONT=1  FUSELAGE DATA                                                NPNL0070
C  KONT=2  WING DATA                                                    NPNL0080
C  KONT=3  HORIZONTAL TAIL 1 DATA                                       NPNL0090
C  KONT=4  HORIZONTAL TAIL 2 DATA                                       NPNL0100
C                                                                       NPNL0110
C.......................................................................NPNL0120
      COMMON/CNTROL/ LPANEL,NEAMAX,NCNTL,HALFSW,CREF,AM,CLDES,          NPNL0130
     . NCP2,NCP3,NCP4,WNGPNL,HT1PNL,HT2PNL,NCPHT1,NCPHT2,NWNGP,         NPNL0140
     . NHT1P,NHT2P,HALFB,ALPH,HALFB1,HALFB2,NCPSWL,NCPWNG,MAX           NPNL0150
      COMMON XCG(300),XCP(300),YCP(300),ZCP(300),AREA(300),AMASS(300),  NPNL0160
     . XN(300,4),YN(300,4),ZN(300,4),BPRIM(300,4)                       NPNL0170
      COMMON/PNL/CPCWL(35),CPSWL(35),IPANEL                             NPNL0180
      DIMENSION XXL(2),YL(2),XXT(2)                                     NPNL0190
      INTEGER WNGPNL,HT1PNL,HT2PNL                                      NPNL0200
      KKK=0                                                             NPNL0210
      IPANEL=1                                                          NPNL0220
      NWNGP=0                                                           NPNL0230
      NHT1P=0                                                           NPNL0240
      NHT2P=0                                                           NPNL0250
C***********************************************************************NPNL0260
      READ (5,8) NCP1,NCP2,NCP3,NCP4                                    NPNL0270
      READ (5,8) NC                                                     NPNL0280
C***********************************************************************NPNL0290
      NCPFUS=NCP1-1                                                     NPNL0300
      NCPWNG=NCP2-1                                                     NPNL0310
      NCPHT1=NCP3-1                                                     NPNL0320
      NCPHT2=NCP4-1                                                     NPNL0330
      WRITE (6,9) NC                                                    NPNL0340
      KONT=0                                                            NPNL0350
      DO 6 KK=1,NC                                                      NPNL0360
C***********************************************************************NPNL0370
      READ (5,8) NCPSWL,KONTL                                           NPNL0380
      IF (KONTL.EQ.1) NCPCWL=NCP1                                       NPNL0390
      IF (KONTL.EQ.2) NCPCWL=NCP2                                       NPNL0400
      IF (KONTL.EQ.3) NCPCWL=NCP3                                       NPNL0410
      IF (KONTL.EQ.4) NCPCWL=NCP4                                       NPNL0420
      READ (5,10) (XXL(I),XXT(I),YL(I),I=1,2),Z                         NPNL0430
      IF (KONT.NE.KONTL) READ (5,10) (CPCWL(I),I=1,NCPCWL)              NPNL0440
      READ (5,10) (CPSWL(K),K=1,NCPSWL)                                 NPNL0450
C***********************************************************************NPNL0460
      IF (KONTL.EQ.2) GO TO 1                                           NPNL0470
      GO TO 2                                                           NPNL0480
    1 NWNGP=NWNGP+NCPSWL-1                                              NPNL0490
      HALFB=YL(2)                                                       NPNL0500
    2 IF (KONTL.EQ.3) NHT1P=NHT1P+NCPSWL-1                              NPNL0510
      IF (KONTL.EQ.4) NHT2P=NHT2P+NCPSWL-1                              NPNL0520
      IF (KONTL.EQ.3) HALFB1=YL(2)                                      NPNL0530
      IF (KONTL.EQ.4) HALFB2=YL(2)                                      NPNL0540
C                                                                       NPNL0550
C  SUBROUTINE RITXYZ WRITES THE GEOMETRY DATA                           NPNL0560
C                                                                       NPNL0570
      CALL RITXYZ(KKK,XXL,XXT,YL,CPCWL,NCPCWL,CPSWL,NCPSWL,Z)           NPNL0580
      IF (KONT.LT.2.AND.KONTL.EQ.2) GO TO 3                             NPNL0590
      GO TO 4                                                           NPNL0600
    3 CREF=XXT(1)-XXL(1)                                                NPNL0610
      WNGPNL=IPANEL                                                     NPNL0620
      GO TO 5                                                           NPNL0630
    4 IF (KONT.LT.3.AND.KONTL.EQ.3) HT1PNL=IPANEL                       NPNL0640
      IF (KONT.LT.4.AND.KONTL.EQ.4) HT2PNL=IPANEL                       NPNL0650
C                                                                       NPNL0660
C  SUBROUTINE PANEL DIVIDES A LIFTING SURFACE INTO AERODYNAMIC          NPNL0670
C  PANELS                                                               NPNL0680
C                                                                       NPNL0690
    5 CALL PANEL(XXL,YL,XXT,NCPCWL,NCPSWL,Z)                            NPNL0700
      IPANEL=LPANEL+1                                                   NPNL0710
      KONT=KONTL                                                        NPNL0720
    6 CONTINUE                                                          NPNL0730
      IM=HT1PNL-1                                                       NPNL0740
      IF (HT1PNL.LE.0) IM=LPANEL                                        NPNL0750
      HALFSW=0.                                                         NPNL0760
      DO 7 I=WNGPNL,IM                                                  NPNL0770
    7 HALFSW=HALFSW+AREA(I)                                             NPNL0780
C***********************************************************************NPNL0790
      READ (5,10) SREF,CCREF                                            NPNL0800
C***********************************************************************NPNL0810
      IF (SREF.NE.0.) HALFSW=SREF/2.                                    NPNL0820
      IF (CCREF.NE.0.) CREF=CCREF                                       NPNL0830
      LSTPNL=LPANEL                                                     NPNL0840
      SW=2.*HALFSW                                                      NPNL0850
      WRITE (6,11) SW                                                   NPNL0860
      WRITE (6,12) CREF                                                 NPNL0870
      RETURN                                                            NPNL0880
    8 FORMAT (4I2)                                                      NPNL0890
    9 FORMAT (1H1,//,11X,45H********************************************NPNL0900
     .*,/,11X,32H* GEOMETRY DATA IS DIVIDED IN TO,I2,11H SECTIONS *,/,11NPNL0910
     .X,45H*********************************************,///)           NPNL0920
   10 FORMAT (8F10.5)                                                   NPNL0930
   11 FORMAT (10X,10HREF. AREA=,E12.5)                                  NPNL0940
   12 FORMAT (10X,11HREF. CHORD=,E12.5)                                 NPNL0950
      END                                                               NPNL0960
      SUBROUTINE ORDER(X,Y,CHORD,INTL,XWLE,XWTE,YL,II)                  ORDR0010
C.......................................................................ORDR0020
C                                                                       ORDR0030
C  ORDER ORDERS THE BREAK LINES FROM INBOARD TO OUTBOARD.               ORDR0040
C                                                                       ORDR0050
C.......................................................................ORDR0060
      DIMENSION X(1),Y(1),CHORD(1),XWLE(1),XWTE(1),YL(1)                ORDR0070
      INT1=INTL+1                                                       ORDR0080
      DO 2 I=1,INT1                                                     ORDR0090
      JJ=1                                                              ORDR0100
      YL(I)=Y(1)                                                        ORDR0110
      XWLE(I)=X(1)                                                      ORDR0120
      XWTE(I)=CHORD(1)                                                  ORDR0130
      DO 1 J=2,INTL                                                     ORDR0140
      IF (YL(I).LE.Y(J)) GO TO 1                                        ORDR0150
      YL(I)=Y(J)                                                        ORDR0160
      XWLE(I)=X(J)                                                      ORDR0170
      XWTE(I)=CHORD(J)                                                  ORDR0180
      JJ=J                                                              ORDR0190
    1 CONTINUE                                                          ORDR0200
      Y(JJ)=1.E06                                                       ORDR0210
    2 CONTINUE                                                          ORDR0220
      I=0                                                               ORDR0230
      DO 3 J=1,INTL                                                     ORDR0240
      IF (XWLE(J).EQ.XWLE(J+1).AND.YL(J).EQ.YL(J+1).AND.XWTE(J).EQ.XWTE(ORDR0250
     .J+1)) GO TO 3                                                     ORDR0260
      I=I+1                                                             ORDR0270
      X(I)=XWLE(J)                                                      ORDR0280
      Y(I)=YL(J)                                                        ORDR0290
      CHORD(I)=XWTE(J)                                                  ORDR0300
    3 CONTINUE                                                          ORDR0310
      II=I                                                              ORDR0320
      DO 4 J=1,II                                                       ORDR0330
      XWLE(J)=X(J)                                                      ORDR0340
      YL(J)=Y(J)                                                        ORDR0350
    4 XWTE(J)=CHORD(J)                                                  ORDR0360
      RETURN                                                            ORDR0370
      END                                                               ORDR0380
      SUBROUTINE PANEL(XXL,YL,XXT,NCW,NSW,Z)                            PANL0010
C.......................................................................PANL0020
C                                                                       PANL0030
C  PANEL DIVIDES A LIFTING SURFACE INTO AERODYNAMIC PANELS.             PANL0040
C                                                                       PANL0050
C.......................................................................PANL0060
      COMMON/CNTROL/ LPANEL,NEAMAX,NCNTL,HALFSW,CREF,AM,CLDES,          PANL0070
     . NCP2,NCP3,NCP4,WNGPNL,HT1PNL,HT2PNL,NCPHT1,NCPHT2,NWNGP,         PANL0080
     . NHT1P,NHT2P,HALFB,ALPH,HALFB1,HALFB2,NCPSWL,NCPWNG,MAX           PANL0090
      COMMON XCG(300),XCP(300),YCP(300),ZCP(300),AREA(300),AMASS(300),  PANL0100
     . XN(300,4),YN(300,4),ZN(300,4),BPRIM(300,4)                       PANL0110
      COMMON/PNL/CPCWL(35),CPSWL(35),IPANEL                             PANL0120
      DIMENSION XXL(2),YL(2),XXT(2)                                     PANL0130
      INTEGER WNGPNL,HT1PNL,HT2PNL                                      PANL0140
      DIMENSION C(2),X(35,35),Y(35,35),SLOPE(35),XL(2,35)               PANL0150
      NSW1=NSW-1                                                        PANL0160
      NCW1=NCW-1                                                        PANL0170
      DO 1 I=1,2                                                        PANL0180
      C(I)=XXT(I)-XXL(I)                                                PANL0190
      DO 1 J=1,NCW                                                      PANL0200
    1 XL(I,J)=XXL(I)+CPCWL(J)*C(I)/100.                                 PANL0210
      SPAN=YL(2)-YL(1)                                                  PANL0220
      DO 2 J=1,NCW                                                      PANL0230
    2 SLOPE(J)=(XL(2,J)-XL(1,J))/SPAN                                   PANL0240
      DO 3 K=1,NSW                                                      PANL0250
      YK=CPSWL(K)*SPAN/100.                                             PANL0260
      DO 3 J=1,NCW                                                      PANL0270
      Y(J,K)=YK+YL(1)                                                   PANL0280
      X(J,K)=XL(1,J)+SLOPE(J)*(Y(J,K)-YL(1))                            PANL0290
    3 CONTINUE                                                          PANL0300
      DO 6 K=1,NSW1                                                     PANL0310
      DO 6 J=1,NCW1                                                     PANL0320
      NPANEL=(K-1)*NCW1+J-1+IPANEL                                      PANL0330
      DO 5 I=1,4                                                        PANL0340
      ZN(NPANEL,I)=Z                                                    PANL0350
      KI1=K+I-1                                                         PANL0360
      KI3=K+I-3                                                         PANL0370
      IF (I.LE.2) GO TO 4                                               PANL0380
      XN(NPANEL,I)=X(J+1,KI3)                                           PANL0390
      YN(NPANEL,I)=Y(J+1,KI3)                                           PANL0400
      BPRIM(NPANEL,I)=SLOPE(J+1)                                        PANL0410
      GO TO 5                                                           PANL0420
      YN(NPANEL,I)=Y(J,KI1)                                             PANL0440
      BPRIM(NPANEL,I)=SLOPE(J)                                          PANL0450
    4 XN(NPANEL,I)=X(J,KI1)                                             PANL0430
    5 CONTINUE                                                          PANL0460
      A1=X(J+1,K)-X(J,K)                                                PANL0470
      A2=X(J+1,K+1)-X(J,K+1)                                            PANL0480
      B=Y(J,K+1)-Y(J,K)                                                 PANL0490
      YBAR=B*(A1+2.*A2)/(3.*(A1+A2))                                    PANL0500
      CBAR=A1+(A2-A1)*YBAR/B                                            PANL0510
      XBAR=X(J,K)+YBAR*(X(J,K+1)-X(J,K))/B                              PANL0520
      FE=X(J+1,K+1)-X(J+1,K)                                            PANL0530
      AE=A1+FE                                                          PANL0540
      BC=X(J,K+1)-X(J,K)                                                PANL0550
      XK=(AE*AE+2.*AE*A1-A1*FE-BC*BC)/(3.*(A1+A2))                      PANL0560
      XCG(NPANEL)=X(J,K)+XK                                             PANL0570
      XCP(NPANEL)=XBAR+0.95*CBAR                                        PANL0580
      YCP(NPANEL)=Y(J,K)+YBAR                                           PANL0590
      ZCP(NPANEL)=Z                                                     PANL0600
      AREA(NPANEL)=(A1+A2)/2.*B                                         PANL0610
    6 CONTINUE                                                          PANL0620
      LPANEL=NPANEL                                                     PANL0630
      RETURN                                                            PANL0640
      END                                                               PANL0650
      SUBROUTINE RITXYZ(KK,XXL,XXT,YL,CPCWL,NCPCWL,CPSWL,NCPSWL,Z)      RITX0010
C                                                                       RITX0020
C.......................................................................RITX0030
C                                                                       RITX0040
C  RITXYZ PRINTS THE GEOMETRY DATA                                      RITX0050
C                                                                       RITX0060
C.......................................................................RITX0070
      DIMENSION XXL(1),XXT(1),YL(1),CPCWL(1),CPSWL(1)                   RITX0080
      KK=KK+1                                                           RITX0090
      WRITE (6,1) KK                                                    RITX0100
      WRITE (6,2)                                                       RITX0110
      WRITE (6,3) (XXL(I),XXT(I),YL(I),Z,I=1,2)                         RITX0120
      WRITE (6,4) NCPCWL,(CPCWL(I),I=1,NCPCWL)                          RITX0130
      WRITE (6,5) NCPSWL,(CPSWL(I),I=1,NCPSWL)                          RITX0140
    1 FORMAT (1H0,10X,21HDEFINITION OF SECTION,I2,/,11X,23H*************RITX0150
     .**********)                                                       RITX0160
    2 FORMAT (1H0,10X,14HX-LEADING EDGE,5X,15HX-TRAILING EDGE,5X,14HY-LERITX0170
     .ADING EDGE,12X,1HZ)                                               RITX0180
    3 FORMAT (4(10X,F11.5))                                             RITX0190
    4 FORMAT (10X,10HTHERE ARE ,I2,36H CONSTANT PERCENT CHORDWISE LINES RITX0200
     .AT,/,10X,(12F10.3))                                               RITX0210
    5 FORMAT (10X,10HTHERE ARE ,I2,37H CONSTANT PERCENT STREAMWISE LINESRITX0220
     . AT,/,10X,(12F10.3))                                              RITX0230
      RETURN                                                            RITX0240
      END                                                               RITX0250
      SUBROUTINE WIDTH(XFUS,YFUS,NFUSXF,XLE,XTE,FUSWTH)                 WIDT0010
C                                                                       WIDT0020
C.......................................................................WIDT0030
C                                                                       WIDT0040
C  WIDTH COMPUTES THE AVERAGE FUSELAGE WIDTH BETWEEN TWO GIVEN STATIONS.WIDT0050
C                                                                       WIDT0060
C.......................................................................WIDT0070
      DIMENSION XFUS(1),YFUS(1)                                         WIDT0080
      K=0                                                               WIDT0090
      SUM=0.                                                            WIDT0100
      DO 1 I=1,NFUSXF                                                   WIDT0110
      IF (XFUS(I).LT.XLE.OR.XFUS(I).GT.XTE) GO TO 1                     WIDT0120
      K=K+1                                                             WIDT0130
      SUM=SUM+YFUS(I)                                                   WIDT0140
    1 CONTINUE                                                          WIDT0150
      IF (SUM.EQ.0.) GO TO 2                                            WIDT0160
      FUSWTH=SUM/FLOAT(K)                                               WIDT0170
      RETURN                                                            WIDT0180
    2 FUSWTH=0.                                                         WIDT0190
      RETURN                                                            WIDT0200
      END                                                               WIDT0210
