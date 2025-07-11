      PROGRAM AERELAS
!!!     0PROGRAM AERELAS(INPUT=201,OUTPUT=201,TAPE5=INPUT,TAPE6=OUTPUT,    ARLS0010
!!!     . TAPE7=1001,TAPE20=401)                                           ARLS0020
C.......................................................................ARLS0030
C                                                                       ARLS0040
C  AERELAS GENERATES OR READS FROM TAPE7 THE STRUCTURAL INFLUENCE       ARLS0050
C  COEFFICIENT MATRIX (C) USED IN AERELA.                               ARLS0060
C                                                                       ARLS0070
C  IF NCNTL = 0, (C) IS COMPUTED.                                       ARLS0080
C  IF NCNTL =+3, (C) IS READ FROM TAPE7 BY ROWS.                        ARLS0090
C  IF NCNTL =-3, (C) IS READ FROM TAPE7 BY COLUMNS.                     ARLS0100
C                                                                       ARLS0110
C  THIS PROGRAM CANNOT HANDLE A CONFIGURATION WITH CANARD SURFAC        ARLS0120
C                                                                       ARLS0130
C.......................................................................ARLS0140
      COMMON IND(1)                                                     ARLS0150
      COMMON/CNTROL/ N,NEAMAX,NCNTL,HALFSW,CREF,AM                      ARLS0160
      COMMON/INDMS/ IA,IB,IC                                            ARLS0170
      DIMENSION D(1)                                                    ARLS0180
      EQUIVALENCE(IND,D),(NEAMAX,M)                                     ARLS0190
C                                                                       ARLS0200
C  OPEN RANDOM-ACCESS FILE, READ FIRST RECORD                           ARLS0210
      CALL OPENMS(20,IND,14,0)                                          ARLS0220
      CALL READMS(20,N,6,1)                                             ARLS0230
      NN=N*N                                                            ARLS0240
      IC=2*N+1                                                          ARLS0250
C                                                                       ARLS0260
C  PARTITION COMMON FOR PROBLEM ARRAYS                                  ARLS0270
      IBASE = 15                                                        ARLS0280
      IXP=IBASE+3*N+1                                                   ARLS0290
      IYP=IXP+N                                                         ARLS0300
      IZP=IYP+N                                                         ARLS0310
      IXG=IZP+N                                                         ARLS0320
      IAR=IXG+N                                                         ARLS0330
      IAM=IAR+N                                                         ARLS0340
      IXN=IAM+N                                                         ARLS0350
C                                                                       ARLS0360
      N4=4*N                                                            ARLS0370
      IYN=IXN+N4                                                        ARLS0380
      IZN=IYN+N4                                                        ARLS0390
      IBP=IZN+N4                                                        ARLS0400
      ICP=IBP+N4                                                        ARLS0410
      IMM=ICP+N                                                         ARLS0420
      IWA=IMM+NN                                                        ARLS0430
      IWB=IWA+N                                                         ARLS0440
C                                                                       ARLS0450
C  READ PROBLEM REDUCED INPUT DATA FROM DISK AND REKEY.                 ARLS0460
     0CALL DISCRD(IND(1),D(IXP),D(IYP),D(IZP),D(IXG),D(IAR),D(IAM),     ARLS0470
     . D(IXN),D(IYN),D(IZN),D(IBP),D(ICP),N)                            ARLS0480
C                                                                       ARLS0490
      N3=3*N+1                                                          ARLS0500
      CALL STINDX(20,IND(IBASE),N3)                                     ARLS0510
    2 M1=M-1                                                            ARLS0520
      MM=M*M                                                            ARLS0530
      INC=1                                                             ARLS0540
      INPHUP=IMM                                                        ARLS0550
      INPHUM=INPHUP+MM                                                  ARLS0560
      INPHUT=INPHUM+MM                                                  ARLS0570
      INXEA=INPHUT+MM                                                   ARLS0580
      INYEA=INXEA+M                                                     ARLS0590
      INPHIP=INYEA+M                                                    ARLS0600
      INPHIM=INPHIP+M                                                   ARLS0610
      INPHIT=INPHIM+M                                                   ARLS0620
      INEI=INPHIT+M                                                     ARLS0630
      INGJ=INEI+M1                                                      ARLS0640
      INCA=INGJ+M1                                                      ARLS0650
      INGAMA=INCA+M1                                                    ARLS0660
      ING=INGAMA+M1                                                     ARLS0670
      INSTRP=ING+M1*M1                                                  ARLS0680
      INSTRM=INSTRP+M                                                   ARLS0690
C                                                                       ARLS0700
C  LINK  LINK22  CALCULATES THE STRUCTURAL INFLUENCE COEFFICIENT        ARLS0710
C  MATRIX                                                               ARLS0720
C                                                                       ARLS0730
C                                                                       ARLS0740
C  STRUCTURAL MATRIX                                                    ARLS0750
     0CALL LINK2L(D(IXG),D(IYP),D(INPHUP),D(INPHUM),D(INPHUT),D(INXEA), ARLS0760
     . D(INYEA),D(INPHIP),D(INPHIM),D(INPHIT),D(INEI),D(INGJ),D(INCA),  ARLS0770
     . D(INGAMA),D(ING),D(INSTRP),D(INSTRM),M,M1,N  ,D(IWA),D(IWB),     ARLS0780
     . IND(IBASE),D(IMM))                                               ARLS0790
      IF(NCNTL.EQ.2) STOP                                               ARLS0800
C                                                                       ARLS0810
C  RESET DISC INDEX KEY FOR AERELA                                      ARLS0820
      CALL STINDX(20,IND(1),14)                                         ARLS0830
      ICBASE=IC+IBASE                                                   ARLS0840
      CALL WRITMS(20,IND(ICBASE),N,12)                                  ARLS0850
      STOP                                                              ARLS0860
      END                                                               ARLS0870
      SUBROUTINE DISCRD(IND,XCP,YCP,ZCP,XCG,AREA,AMASS,XN,YN,ZN,BPRIM,  DISC0010
     . DUM,N)                                                           DISC0020
C.......................................................................DISC0030
C     SUBROUTINE DISCRD READS THE DATA STORED ON TAPE20 INTO CORE       DISC0040
C.......................................................................DISC0050
      DIMENSION XCP(N),YCP(N),ZCP(N),XCG(N),AREA(N),AMASS(N),           DISC0060
     . XN(N,4),YN(N,4),ZN(N,4),BPRIM(N,4)                               DISC0070
      DIMENSION IND(12),DUM(4,N)                                        DISC0080
      CALL READMS(20,XCP(1),N,2)                                        DISC0090
      CALL READMS(20,YCP(1),N,3)                                        DISC0100
      CALL READMS(20,ZCP(1),N,4)                                        DISC0110
      CALL READMS(20,XCG(1),N,5)                                        DISC0120
      CALL READMS(20,AREA(1),N,6)                                       DISC0130
      CALL READMS(20,AMASS(1),N,7)                                      DISC0140
      N4=4*N                                                            DISC0150
      CALL RESTORE(DUM,IND,N,N4,8,XN)                                   DISC0160
      CALL RESTORE(DUM,IND,N,N4,9,YN)                                   DISC0170
      CALL RESTORE(DUM,IND,N,N4,10,ZN)                                  DISC0180
      CALL RESTORE(DUM,IND,N,N4,11,BPRIM)                               DISC0190
      RETURN                                                            DISC0200
      END                                                               DISC0210
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
      SUBROUTINE LINK2L(XCG,YCG,PHUP,PHUM,PHUT,XEA,YEA,PHIUP,PHIUM,     LN2L0010
     . PHIUT,EI,GJ,A,GAMMA,G,STORP,STORM,M1,M2,NUMPNL,IASIGN,C,IND,CC)  LN2L0020
C                                                                       LN2L0030
C.......................................................................LN2L0040
C                                                                       LN2L0050
C  LINKZL CONTROLS COMPUTATION OF THE STRUCTURAL INFLUENCE MATRIX.      LN2L0060
C                                                                       LN2L0070
C.......................................................................LN2L0080
C                                                                       LN2L0090
      DIMENSION C(1),PHUP(M1,M1),PHUM(M1,M1),PHUT(M1,M1),G(M2,M2)       LN2L0100
      DIMENSION XCG(1),YCG(1),IASIGN(1),XEA(1),YEA(1),PHIUP(1),PHIUM(1),LN2L0110
     .PHIUT(1),EI(1),GJ(1),A(1),GAMMA(1),STORP(1),STORM(1)              LN2L0120
      DIMENSION HEAD(12),LPNL(4)                                        LN2L0130
      DIMENSION CC(NUMPNL,NUMPNL)                                       LN2L0140
      INTEGER FPNLFF,FPNLRF,FPNLMW,FPNLHT                               LN2L0150
      DIMENSION IND(1)                                                  LN2L0160
      COMMON/INDMS/IA,IB,IC                                             LN2L0170
      INTEGER WNGPNL,HT1PNL,HT2PNL                                      LN2L0180
      COMMON /PNLNUM/ FPNLFF,LPNLFF,FPNLRF,LPNLRF,FPNLMW,LPNLMW,FPNLHT,LLN2L0190
     .PNLHT,MOVTEL                                                      LN2L0200
      COMMON/CNTROL/LPANEL,NEAMAX,NCNTL,HALFSW,CREF,AM                  LN2L0210
      DATA (HEAD(I),I=1,12)/5HFRONT,5H FUSE,5HLAGE ,5H REAR,5H FUSE,5HLALN2L0220
     .GE ,5H   MA,5HIN WI,5HNG   ,5HHORIZ,5HONTAL,5H TAIL/              LN2L0230
C.......................................................................LN2L0240
C                                                                       LN2L0250
C  NCONTL=1     FRONT FUSELAGE ONLY                                     LN2L0260
C  NCONTL=2     REAR FUSELAGE ONLY                                      LN2L0270
C  NCONTL=3     MAIN WING ONLY                                          LN2L0280
C  NCONTL=4     HORIZONTAL TAIL ONLY                                    LN2L0290
C  N1= ELASTIC AXIS POINT ATTACHED TO HORIZONTAL TAIL                   LN2L0300
C  MOVTEL=0 NON-MOVABLE TAIL                                            LN2L0310
C  MOVTEL=1 ALL MOVABLE TAIL                                            LN2L0320
C                                                                       LN2L0330
C.......................................................................LN2L0340
C                                                                       LN2L0350
C  SPECIAL CODE FOR READING THE STRUCTURAL MATRIX FROM TAPE7.           LN2L0360
C                                                                       LN2L0370
C  IF NCNTL = +3, READ MATRIX (C) FROM TAPE7 BY ROWS AND RETURN.        LN2L0380
C  IF NCNTL = -3, READ MATRIX (C) FROM TAPE7 BY COLUMNS, STORE BY       LN2L0390
C  ROWS, AND RETURN                                                     LN2L0400
C  NOTE.  AEREL EXPECTS EACH ROW (OR COLUMN) TO BE A LOGICAL RECORD.    LN2L0410
C  HOWEVER, IF AN END OF FILE IS ENCOUNTERED AFTER THE FIRST READ,      LN2L0420
C  AN ATTEMPT WILL BE MADE TO READ THE ENTIRE MATRIX FROM TAPE7.        LN2L0430
C                                                                       LN2L0440
      NCP=NUMPNL                                                        LN2L0450
      IF(IABS(NCNTL).NE.3) GO TO 99                                     LN2L0460
      REWIND 7                                                          LN2L0470
      DO 40 I=1,NCP                                                     LN2L0480
      ICC=IC+I-1                                                        LN2L0490
      READ(7) (C(J),J=1,NCP)                                            LN2L0500
      IF(EOF,7) 50,40                                                   LN2L0510
   40 CALL WRITMS(20,C,NCP,ICC)                                         LN2L0520
      IF(NCNTL) 80,80,90                                                LN2L0530
C                                                                       LN2L0540
C  EOF ENCOUNTERED                                                      LN2L0550
   50 NREC=I-1                                                          LN2L0560
      IF(NREC) 55,55,60                                                 LN2L0570
   55 WRITE(6,56)                                                       LN2L0580
   56 FORMAT(1H1,"NO DATA FOUND ON TAPE7.  EXECUTION TERMINATED")       LN2L0590
      CALL EXIT                                                         LN2L0600
C                                                                       LN2L0610
   60 WRITE(6,61) NREC                                                  LN2L0620
   610FORMAT(1H1,I5," RECORDS FOUND ON TAPE7"/"     MATRIX READ OF FIRSTLN2L0630
     . RECORD BEING ATTEMPTED")                                         LN2L0640
      REWIND 7                                                          LN2L0650
      READ(7) ((CC(I,J),J=1,NCP),I=1,NCP)                               LN2L0660
C                                                                       LN2L0670
      DO 70 I=1,NCP                                                     LN2L0680
      ICC=IC+I-1                                                        LN2L0690
      DO 65 J=1,NCP                                                     LN2L0700
   65 C(J)=CC(I,J)                                                      LN2L0710
   70 CALL WRITMS(20,C,NCP,ICC)                                         LN2L0720
      IF(NCNTL) 80,80,90                                                LN2L0730
C                                                                       LN2L0740
C  MATRIX WAS READ BY COLUMNS.  STORE BY ROWS.                          LN2L0750
   80 CALL CWRITE(PHUP(1),C,NCP,IND)                                    LN2L0760
   90 REWIND 7                                                          LN2L0770
      RETURN                                                            LN2L0780
C.......................................................................LN2L0790
C                                                                       LN2L0800
   99 CONTINUE                                                          LN2L0810
      DO 29 I=1,NCP                                                     LN2L0820
   29 C(I)=0.0                                                          LN2L0830
      DO 31 I=1,NCP                                                     LN2L0840
      ICC=IC+I-1                                                        LN2L0850
   31 CALL WRITMS(20,C,NCP,ICC)                                         LN2L0860
C***********************************************************************LN2L0870
      READ(5,12)(LPNL(I),I=1,4),N1,MOVTEL,NEIGJ                         LN2L0880
C***********************************************************************LN2L0890
      FPNLFF=1                                                          LN2L0900
      LPNLFF=LPNL(1)                                                    LN2L0910
      FPNLRF=LPNLFF+1                                                   LN2L0920
      LPNLRF=LPNLFF+LPNL(2)                                             LN2L0930
      FPNLMW=LPNLRF+1                                                   LN2L0940
      LPNLMW=LPNLRF+LPNL(3)                                             LN2L0950
      FPNLHT=LPNLMW+1                                                   LN2L0960
      LPNLHT=LPNLMW+LPNL(4)                                             LN2L0970
      IF (LPNLFF.GE.FPNLFF) WRITE (6,16) FPNLFF,LPNLFF                  LN2L0980
      IF (LPNLRF.GE.FPNLRF) WRITE (6,17) FPNLRF,LPNLRF                  LN2L0990
      IF (LPNLMW.GE.FPNLMW) WRITE (6,18) FPNLMW,LPNLMW                  LN2L1000
      IF (LPNLHT.GE.FPNLHT) WRITE (6,19) FPNLHT,LPNLHT                  LN2L1010
C***********************************************************************LN2L1020
      READ (5,9) (IASIGN(I),I=1,NCP)                                    LN2L1030
      READ (5,10) CSMW,CSHT,WDTHMW,WDTHHT,XREF                          LN2L1040
C***********************************************************************LN2L1050
      IF (CSMW.NE.0.) WRITE (6,20) CSMW                                 LN2L1060
      IF (CSHT.NE.0.) WRITE (6,21) CSHT                                 LN2L1070
      IF (WDTHMW.NE.0.) WRITE (6,22) WDTHMW                             LN2L1080
      IF (WDTHHT.NE.0.) WRITE (6,23) WDTHHT                             LN2L1090
      WRITE (6,24) XREF                                                 LN2L1100
C***********************************************************************LN2L1110
      IF (NCNTL.EQ.2) READ (5,10) (XCG(I),YCG(I),I=1,NUMPNL)            LN2L1120
C***********************************************************************LN2L1130
      WRITE (6,13) (I,XCG(I),YCG(I),IASIGN(I),I=1,NCP)                  LN2L1140
    1 CONTINUE                                                          LN2L1150
C***********************************************************************LN2L1160
      READ (5,9) NCONTL,NEA                                             LN2L1170
C***********************************************************************LN2L1180
      IF (NCONTL.GT.4) GO TO 7                                          LN2L1190
      IF (NCONTL.EQ.0) GO TO 6                                          LN2L1200
      IF (NEA.GT.M1) GO TO 8                                            LN2L1210
      NN=NEA-1                                                          LN2L1220
C***********************************************************************LN2L1230
      READ (5,10) (XEA(I),YEA(I),I=1,NEA)                               LN2L1240
      READ (5,11) (EI(I),GJ(I),I=1,NN)                                  LN2L1250
C                                                                       LN2L1260
C  SUBROUTINE  CNVERT  IS USED TO CONVERT EI AND GJ VALUES FROM         LN2L1270
C  (LB-INCH**2) UNITS TO (LB-FEET**2) UNITS                             LN2L1280
C                                                                       LN2L1290
      IF (NEIGJ.EQ.0) CALL CNVERT(EI,GJ,NN)                             LN2L1300
C***********************************************************************LN2L1310
      NCCC2=3*NCONTL                                                    LN2L1320
      NCCC1=NCCC2-2                                                     LN2L1330
      WRITE (6,27) (HEAD(I),I=NCCC1,NCCC2),NEA                          LN2L1340
      WRITE (6,14) (I,XEA(I),YEA(I),I=1,NEA)                            LN2L1350
      WRITE (6,15) (I,EI(I),GJ(I),I=1,NN)                               LN2L1360
      GO TO (2,3,4,5),NCONTL                                            LN2L1370
C                                                                       LN2L1380
C  USE OF  CTHETA  FOR FRONT FUSELAGE                                   LN2L1390
C                                                                       LN2L1400
    2 NCP1=FPNLFF                                                       LN2L1410
      NCP2=LPNLFF                                                       LN2L1420
      NW1=FPNLMW                                                        LN2L1430
      NW2=LPNLMW                                                        LN2L1440
      GO TO 30                                                          LN2L1450
C                                                                       LN2L1460
C  USE OF  CTHETA  FOR REAR FUSELAGE                                    LN2L1470
C                                                                       LN2L1480
    3 NCP1=FPNLRF                                                       LN2L1490
      NCP2=LPNLRF                                                       LN2L1500
      NW1=FPNLMW                                                        LN2L1510
      NW2=LPNLMW                                                        LN2L1520
      NH1=FPNLHT                                                        LN2L1530
      NH2=LPNLHT                                                        LN2L1540
      GO TO 30                                                          LN2L1550
C                                                                       LN2L1560
C  USE OF  CTHETA  FOR WING                                             LN2L1570
C                                                                       LN2L1580
    4 NCP1=FPNLMW                                                       LN2L1590
      NCP2=LPNLMW                                                       LN2L1600
      GO TO 30                                                          LN2L1610
C                                                                       LN2L1620
C  USE OF  CTHETA  FOR HORIZONTAL TAIL                                  LN2L1630
C                                                                       LN2L1640
    5 NCP1=FPNLHT                                                       LN2L1650
      NCP2=LPNLHT                                                       LN2L1660
      NH1=FPNLHT                                                        LN2L1670
      NH2=LPNLHT                                                        LN2L1680
      NC1=FPNLRF                                                        LN2L1690
      NC2=LPNLRF                                                        LN2L1700
C                                                                       LN2L1710
C  SUBROUTINE  CTHETA  IS THE ACTUAL SUBROUTINE WHICH CALCULATES        LN2L1720
C  STRUCTURAL MATRIX                                                    LN2L1730
C                                                                       LN2L1740
   30 CALL CTHETA(NCONTL,XEA,YEA,XCG,YCG,EI,GJ,IASIGN,CSMW,WDTHMW,CSHT,WLN2L1750
     .DTHHT,NCP1,NCP2,NW1,NW2,N1,NH1,NH2,NC1,NC2,NEA,NN,PHIUP,PHIUM,PHIULN2L1760
     .T,PHUP,PHUM,PHUT,A,STORP,STORM,GAMMA,G,C,NCP,XREF,M1,M2,IND)      LN2L1770
      GO TO 1                                                           LN2L1780
C  CWRITE STORES THE C MATRIX ON DISC BY ROWS                           LN2L1790
    6 CALL CWRITE(PHUP(1),C,NUMPNL,IND)                                 LN2L1800
C                                                                       LN2L1810
C  WRITEC PRINTS OUT THE INFLUENCE COEFFICIENT MATRIX                   LN2L1820
      IF(NCNTL.EQ.2.OR.NCNTL.EQ.4) CALL WRITEC(PHUP(1),NCP)             LN2L1830
      RETURN                                                            LN2L1840
    7 IERROR=2                                                          LN2L1850
      WRITE (6,28) IERROR                                               LN2L1860
      WRITE (6,26)                                                      LN2L1870
      RETURN                                                            LN2L1880
    8 IERROR=3                                                          LN2L1890
      WRITE (6,28) IERROR                                               LN2L1900
      WRITE (6,25) NEA,M1                                               LN2L1910
      RETURN                                                            LN2L1920
    9 FORMAT (40I2)                                                     LN2L1930
   10 FORMAT (8F10.5)                                                   LN2L1940
   11 FORMAT (4E15.8)                                                   LN2L1950
   12 FORMAT (10I3)                                                     LN2L1960
   13 FORMAT(1H1,65H COORDINATES OF UNIT LOADING POINTS AND LOADING POINLN2L1970
     .T ASSIGNMENTS,//33H    I       XCG       YCG  IASIGN//(I5,2F10.5,ILN2L1980
     .8))                                                               LN2L1990
   14 FORMAT(///,37H COORDINATES OF ELASTIC AXIS SEGMENTS//25H    I     LN2L2000
     .  XEA       YEA//(I5,2F10.5))                                     LN2L2010
   15 FORMAT(///,46H ELASTIC AXIS TORSIONAL AND BENDING STIFFNESS /     LN2L2020
     .4X,1HI,10X,2HEI,9X,3HKGJ//(I5,2E12.4))                            LN2L2030
   16 FORMAT (10X,38HFIRST PANEL NUMBER ON FRONT FUSELAGE =,I3,/,10X,38HLN2L2040
     .LAST PANEL NUMBER ON FRONT FUSELAGE  =,I3)                        LN2L2050
   17 FORMAT (10X,38HFIRST PANEL NUMBER ON REAR FUSELAGE  =,I3,/,10X,38HLN2L2060
     .LAST PANEL NUMBER ON REAR FUSELAGE   =,I3)                        LN2L2070
   18 FORMAT (10X,38HFIRST PANEL NUMBER ON MAIN WING      =,I3,/,10X,38HLN2L2080
     .LAST PANEL NUMBER ON MAIN WING       =,I3)                        LN2L2090
   19 FORMAT (10X,38HFIRST PANEL NUMBER ON HORIZONTAL TAIL=,I3,/,10X,38HLN2L2100
     .LAST PANEL NUMBER ON HORIZONTAL TAIL =,I3)                        LN2L2110
   20 FORMAT (10X,42H     STRUCTURAL CHORD FOR MAIN WING      =,F10.5)  LN2L2120
   21 FORMAT (10X,42H     STRUCTURAL CHORD FOR HORIZONTAL TAIL=,F10.5)  LN2L2130
   22 FORMAT (10X,64HHALF THE WIDTH OF FUSELAGE AT ROOT CHORD OF THE MAILN2L2140
     .NWING       =,F10.5)                                              LN2L2150
   23 FORMAT (10X,64HHALF THE WIDTH OF FUSELAGE AT ROOT CHORD OF THE HORLN2L2160
     .IZONTAL TAIL=,F10.5)                                              LN2L2170
   24 FORMAT (10X,28HX-COORDINATE OF FIXED POINT=,F10.5)                LN2L2180
   25 FORMAT (10X,45HNEA IS GREATER THAN MAXIMUM ASSUMED NEA VALUE,/,10XLN2L2190
     .,30HFIRST INPUT DATA CARD IS WRONG,/,10X,26HSECOND VARIABLE SHOULDLN2L2200
     . BE ,I2,9H AND NOT ,I2)                                           LN2L2210
   26 FORMAT (10X,31HNCONTL CANNOT BE GREATER THAN 4)                   LN2L2220
   27 FORMAT (1H1,33HNUMBER OF ELASTIC AXIS POINTS ON ,3A5,1H=,I2,///)  LN2L2230
   28 FORMAT (1H1,10X,7HIERROR=,I2)                                     LN2L2240
      END                                                               LN2L2250
      SUBROUTINE CTHETA(NCONTL,XEA,YEA,XCG,YCG,EI,GJ,IASIGN,CSMW,WDTHMW,CTHE0010
     .CSHT,WDTHHT,NCP1,NCP2,NW1,NW2,N1,NH1,NH2,NC1,NC2,NEA,NN,PHIUP,PHIUCTHE0020
     .M,PHIUT,PHUP,PHUM,PHUT,A,STORP,STORM,GAMMA,G,C,NCP,XREF,M1,M2,IND)CTHE0030
C                                                                       CTHE0040
C.......................................................................CTHE0050
C                                                                       CTHE0060
C  CTHETA COMPUTES THE STRUCTURAL MATRIX C.  C RESIDES ON DISC AND      CTHE0070
C  INDIVIDUAL ROWS ARE ACCESSED AS REQUIRED.                            CTHE0080
C                                                                       CTHE0090
C                                                                       CTHE0100
C  NCONTL=1     FRONT FUSELAGE ONLY                                     CTHE0110
C  NCONTL=2     REAR FUSELAGE ONLY                                      CTHE0120
C  NCONTL=3     MAIN WING ONLY                                          CTHE0130
C  NCONTL=4     HORIZONTAL TAIL ONLY                                    CTHE0140
C                                                                       CTHE0150
C.......................................................................CTHE0160
      DIMENSION A(1),XEA(1),YEA(1),XCG(1),YCG(1),EI(1),GJ(1),GAMMA(1)   CTHE0170
      DIMENSION PHIUP(1),PHIUM(1),PHIUT(1)                              CTHE0180
      DIMENSION STORP(1),STORM(1),IASIGN(1)                             CTHE0190
      DIMENSION PHUP(M1,M1),PHUM(M1,M1),PHUT(M1,M1),G(M2,M2)            CTHE0200
      DIMENSION IND(1),C(1)                                             CTHE0210
      COMMON/INDMS/IA,IB,IC                                             CTHE0220
      INTEGER FPNLFF,FPNLRF,FPNLMW,FPNLHT                               CTHE0230
      COMMON /PNLNUM/ FPNLFF,LPNLFF,FPNLRF,LPNLRF,FPNLMW,LPNLMW,FPNLHT,LCTHE0240
     .PNLHT,MOVTEL                                                      CTHE0250
      N=NCP                                                             CTHE0260
      IF (NCONTL.GT.2) GO TO 24                                         CTHE0270
      DO 1 I=2,NEA                                                      CTHE0280
    1 A(I-1)=ABS(XEA(I)-XEA(I-1))                                       CTHE0290
      IF (NCONTL.EQ.1) CK=1.                                            CTHE0300
      IF (NCONTL.EQ.2) CK=-1.                                           CTHE0310
      DO 5 I=1,NN                                                       CTHE0320
      J=NEA-I                                                           CTHE0330
      JP=J+1                                                            CTHE0340
    2 JJ=J+1                                                            CTHE0350
      XM=(XEA(JJ)-XEA(JP))*CK                                           CTHE0360
      PHIUP(JJ)=(XM*A(J)/EI(J)+A(J)*A(J)/(2.*EI(J)))*CK                 CTHE0370
      PHIUM(JJ)=(A(J)/EI(J))*CK                                         CTHE0380
      J=J-1                                                             CTHE0390
      IF (JJ.GT.2) GO TO 2                                              CTHE0400
      PHUP(1,JP)=0.                                                     CTHE0410
      PHUM(1,JP)=0.                                                     CTHE0420
      DO 4 II=2,NEA                                                     CTHE0430
      IF (II.GT.JP) GO TO 3                                             CTHE0440
      PHUP(II,JP)=PHUP(II-1,JP)+PHIUP(II)                               CTHE0450
      PHUM(II,JP)=PHUM(II-1,JP)+PHIUM(II)                               CTHE0460
      GO TO 4                                                           CTHE0470
    3 PHUP(II,JP)=PHUP(JP,JP)                                           CTHE0480
      PHUM(II,JP)=PHUM(JP,JP)                                           CTHE0490
    4 CONTINUE                                                          CTHE0500
    5 CONTINUE                                                          CTHE0510
C                                                                       CTHE0520
C  EFFECT OF FRONT FUSELAGE ON FRONT FUSELAGE                           CTHE0530
C                                                                       CTHE0540
      IF (NCONTL.EQ.1) GO TO 7                                          CTHE0550
      N2=IASIGN(N1)                                                     CTHE0560
      DO 6 II=1,NEA                                                     CTHE0570
      STORP(II)=PHUP(II,N2)                                             CTHE0580
    6 STORM(II)=PHUM(II,N2)                                             CTHE0590
C                                                                       CTHE0600
C  EFFECT OF REAR FUSELAGE ON REAR FUSELAGE                             CTHE0610
C                                                                       CTHE0620
    7 DO 12 K=NCP1,NCP2                                                 CTHE0630
      ICC=IC+K-1                                                        CTHE0640
C                                                                       CTHE0650
C  READ C MATRIX FROM DISC                                              CTHE0660
      CALL READMS(20,C,N,ICC)                                           CTHE0670
      I=IASIGN(K)                                                       CTHE0680
      IF (I.EQ.1) GO TO 8                                               CTHE0690
      XM=(XEA(I)-XCG(K))*CK                                             CTHE0700
      GO TO 10                                                          CTHE0710
    8 DO 9 J=NCP1,NCP2                                                  CTHE0720
      L=IASIGN(J)                                                       CTHE0730
    9 C(L)=0.                                                           CTHE0740
      GO TO 12                                                          CTHE0750
   10 DO 11 J=NCP1,NCP2                                                 CTHE0760
      L=IASIGN(J)                                                       CTHE0770
   11 C(J)=PHUP(L,I)+PHUM(L,I)*XM                                       CTHE0780
C                                                                       CTHE0790
C  REPLACE C MATRIX WITH REAR FUSELAGE ON REAR FUSELAGE EFFECTS         CTHE0800
   12 CALL WRITIN(20,C,N,ICC)                                           CTHE0810
C                                                                       CTHE0820
C  EFFECT OF FRONT FUSELAGE ON MAIN WING                                CTHE0830
C  EFFECT OF REAR FUSELAGE ON MAIN WING                                 CTHE0840
C                                                                       CTHE0850
      NHH1=NCP1                                                         CTHE0860
      NHH2=NCP2                                                         CTHE0870
      NCP11=NCP1                                                        CTHE0880
      IF (NCONTL.EQ.2) NCP11=NCP1-1                                     CTHE0890
      NCP22=NCP2                                                        CTHE0900
   13 CSMW2=CSMW/2.                                                     CTHE0910
      YLIMIT=WDTHMW+CSMW2                                               CTHE0920
      DO 17 K=NHH1,NHH2                                                 CTHE0930
      ICC=IC+K-1                                                        CTHE0940
C  READ C MATRIX FROM DISC                                              CTHE0950
      CALL READMS(20,C,N,ICC)                                           CTHE0960
      DO 16 I=NW1,NW2                                                   CTHE0970
      IF(YCG(I).GT.YLIMIT) GO TO 165                                    CTHE0980
      XCONT=(XCG(I)-XREF)*CK                                            CTHE0990
      IF(XCONT.GE.0.) GO TO 16                                          CTHE1000
C                                                                       CTHE1010
      DO 14 J=NCP11,NCP22                                               CTHE1020
 14   IF(XCG(I).GT.XCG(J).AND.XCG(I).LE.XCG(J+1)) GO TO 15              CTHE1030
C                                                                       CTHE1040
 15   Y=1.-(YCG(I)-WDTHMW)/CSMW2                                        CTHE1050
      X=(XCG(J+1)-XCG(I))/(XCG(J+1)-XCG(J))                             CTHE1060
      CMAX=C(J+1)  +(C(J)  -C(J+1)  )*X                                 CTHE1070
 16   C(I)  =CMAX*Y                                                     CTHE1080
 165  CONTINUE                                                          CTHE1090
C                                                                       CTHE1100
C  REPLACE C MATRIX WITH FRONT AND REAR FUSELAGE ON MAIN WING EFFEC     CTHE1110
 17   CALL WRITIN(20,C,N,ICC)                                           CTHE1120
   18 CONTINUE                                                          CTHE1130
      IF (NCONTL.EQ.4) GO TO 36                                         CTHE1140
      IF (NCONTL.EQ.1) RETURN                                           CTHE1150
C                                                                       CTHE1160
C  EFFECT OF REAR FUSELAGE ON HORIZONTAL TAIL                           CTHE1170
C                                                                       CTHE1180
      CSHT2=CSHT/2.                                                     CTHE1190
      YLIMIT=WDTHHT+CSHT2                                               CTHE1200
      DO 23 I=NH1,NH2                                                   CTHE1210
      IF (YCG(I).GT.YLIMIT) GO TO 21                                    CTHE1220
      DO 19 J=NCP1,NCP2                                                 CTHE1230
      IF (XCG(I).LE.XCG(J+1).AND.XCG(I).GT.XCG(J)) GO TO 20             CTHE1240
   19 CONTINUE                                                          CTHE1250
   20 Y=1.-(YCG(I)-WDTHHT)/CSHT2                                        CTHE1260
      X=(XCG(J)-XCG(I))/(XCG(J+1)-XCG(J))                               CTHE1270
      GO TO 22                                                          CTHE1280
   21 Y=0.                                                              CTHE1290
   22 DO 23 K=NCP1,NCP2                                                 CTHE1300
      ICC=IC+K-1                                                        CTHE1310
C  READ C MATRIX FROM DISC                                              CTHE1320
      CALL READMS(20,C,N,ICC)                                           CTHE1330
      CMAX=C(J)+(C(J)-C(J+1))*X                                         CTHE1340
      C(I)=C(N1)*(1.-Y)+CMAX*Y                                          CTHE1350
C  REPLACE C MATRIX WITH REAR FUSELAGE ON HORIZONTAL TAIL EFFECTS       CTHE1360
   23 CALL WRITIN(20,C,N,ICC)                                           CTHE1370
      RETURN                                                            CTHE1380
   24 CONTINUE                                                          CTHE1390
C                                                                       CTHE1400
C  EFFECT OF MAIN WING ON MAIN WING                                     CTHE1410
C  EFFECT OF HORIZONTAL TAIL ON HORIZONTAL TAIL                         CTHE1420
C  (WITHOUT REAR FUSELAGE BENDING)                                      CTHE1430
C                                                                       CTHE1440
      DO 26 I=2,NEA                                                     CTHE1450
      A(I-1)=SQRT((XEA(I)-XEA(I-1))**2+(YEA(I)-YEA(I-1))**2)            CTHE1460
      IF (YEA(I).EQ.YEA(I-1)) GO TO 25                                  CTHE1470
      GAMMA(I-1)=ATAN((XEA(I)-XEA(I-1))/(YEA(I)-YEA(I-1)))              CTHE1480
      GO TO 26                                                          CTHE1490
   25 GAMMA(I-1)=1.570796325                                            CTHE1500
   26 CONTINUE                                                          CTHE1510
      DO 27 I=1,NN                                                      CTHE1520
      DO 27 J=1,NN                                                      CTHE1530
   27 G(I,J)=GAMMA(J)-GAMMA(I)                                          CTHE1540
      DO 31 I=1,NN                                                      CTHE1550
      J=NEA-I                                                           CTHE1560
      JP=J+1                                                            CTHE1570
   28 JJ=J+1                                                            CTHE1580
      SINE=SIN(GAMMA(J))                                                CTHE1590
      COSINE=COS(GAMMA(J))                                              CTHE1600
      T=(YEA(JP)-YEA(JJ))*SINE-(XEA(JP)-XEA(JJ))*COSINE                 CTHE1610
      XM=(YEA(JP)-YEA(JJ))*COSINE+(XEA(JP)-XEA(JJ))*SINE                CTHE1620
      SINJJP=SIN(G(J,JP-1))                                             CTHE1630
      COSJJP=COS(G(J,JP-1))                                             CTHE1640
      THETP=XM*A(J)/EI(J)+A(J)*A(J)/(2.*EI(J))                          CTHE1650
      PHIP=T*A(J)/GJ(J)                                                 CTHE1660
      THETM=COSJJP*A(J)/EI(J)                                           CTHE1670
      PHIM=-SINJJP*A(J)/GJ(J)                                           CTHE1680
      THETT=SINJJP*A(J)/EI(J)                                           CTHE1690
      PHIT=COSJJP*A(J)/GJ(J)                                            CTHE1700
      PHIUP(JJ)=PHIP*COSINE-THETP*SINE                                  CTHE1710
      PHIUM(JJ)=PHIM*COSINE-THETM*SINE                                  CTHE1720
      PHIUT(JJ)=PHIT*COSINE-THETT*SINE                                  CTHE1730
      J=J-1                                                             CTHE1740
      IF (JJ.GT.2) GO TO 28                                             CTHE1750
      PHUP(1,JP)=0.                                                     CTHE1760
      PHUT(1,JP)=0.                                                     CTHE1770
      PHUM(1,JP)=0.                                                     CTHE1780
      DO 30 II=2,NEA                                                    CTHE1790
      IF (II.GT.JP) GO TO 29                                            CTHE1800
      PHUP(II,JP)=PHUP(II-1,JP)+PHIUP(II)                               CTHE1810
      PHUM(II,JP)=PHUM(II-1,JP)+PHIUM(II)                               CTHE1820
      PHUT(II,JP)=PHUT(II-1,JP)+PHIUT(II)                               CTHE1830
      GO TO 30                                                          CTHE1840
   29 PHUP(II,JP)=PHUP(JP,JP)                                           CTHE1850
      PHUT(II,JP)=PHUT(JP,JP)                                           CTHE1860
      PHUM(II,JP)=PHUM(JP,JP)                                           CTHE1870
   30 CONTINUE                                                          CTHE1880
   31 CONTINUE                                                          CTHE1890
      DO 33 K=NCP1,NCP2                                                 CTHE1900
      ICC=IC+K-1                                                        CTHE1910
C  READ C MATRIX FROM DISC                                              CTHE1920
      CALL READMS(20,C,N,ICC)                                           CTHE1930
      I=IASIGN(K)                                                       CTHE1940
      SING=SIN(GAMMA(I-1))                                              CTHE1950
      COSG=COS(GAMMA(I-1))                                              CTHE1960
      T=(YCG(K)-YEA(I))*SING-(XCG(K)-XEA(I))*COSG                       CTHE1970
      XM=(YCG(K)-YEA(I))*COSG+(XCG(K)-XEA(I))*SING                      CTHE1980
      DO 32 J=NCP1,NCP2                                                 CTHE1990
      L=IASIGN(J)                                                       CTHE2000
   32 C(J)=PHUP(L,I)+PHUM(L,I)*XM+PHUT(L,I)*T                           CTHE2010
C                                                                       CTHE2020
C  REPLACE C MATRIX WITH MAIN WING ON MAIN WING EFFECTS                 CTHE2030
   33 CALL WRITIN(20,C,N,ICC)                                           CTHE2040
      IF (NCONTL.EQ.3) RETURN                                           CTHE2050
C                                                                       CTHE2060
C  EFFECT OF HORIZONTAL TAIL ON REAR FUSELAGE                           CTHE2070
C                                                                       CTHE2080
      DO 35 K=NCP1,NCP2                                                 CTHE2090
      ICC=IC+K-1                                                        CTHE2100
C  READ C MATRIX FROM DISC                                              CTHE2110
      CALL READMS(20,C,N,ICC)                                           CTHE2120
      XM=XCG(K)-XCG(N1)                                                 CTHE2130
      DO 34 J=NC1,NC2                                                   CTHE2140
      L=IASIGN(J)                                                       CTHE2150
   34 C(J)=STORP(L)+STORM(L)*XM                                         CTHE2160
C                                                                       CTHE2170
C  REPLACE C MATRIX WITH HORIZONTAL TAIL ON REAR FUSELAGE EFFECTS       CTHE2180
   35 CALL WRITIN(20,C,N,ICC)                                           CTHE2190
C                                                                       CTHE2200
C  EFFECT OF HORIZONTAL TAIL ON MAIN WING                               CTHE2210
C                                                                       CTHE2220
      NH1=FPNLHT                                                        CTHE2230
      NH2=LPNLHT                                                        CTHE2240
      NHH1=FPNLHT                                                       CTHE2250
      NHH2=LPNLHT                                                       CTHE2260
      NW1=FPNLMW                                                        CTHE2270
      NW2=LPNLMW                                                        CTHE2280
      NCP1=FPNLRF                                                       CTHE2290
      NCP2=LPNLRF                                                       CTHE2300
      CK=-1.                                                            CTHE2310
      NCP11=LPNLFF                                                      CTHE2320
      NCP22=LPNLRF                                                      CTHE2330
      GO TO 13                                                          CTHE2340
   36 CONTINUE                                                          CTHE2350
C                                                                       CTHE2360
C  EFFECT OF HORIZONTAL TAIL ON HORIZONTAL TAIL                         CTHE2370
C  (WITH REAR FUSELAGE BENDING)                                         CTHE2380
C                                                                       CTHE2390
      IF (MOVTEL.EQ.1) GO TO 43                                         CTHE2400
      CSHT2=CSHT/2.                                                     CTHE2410
      YLIMIT=WDTHHT+CSHT2                                               CTHE2420
      DO 42 K=NH1,NH2                                                   CTHE2430
      ICC=IC+K-1                                                        CTHE2440
C  READ C MATRIX FROM DISC                                              CTHE2450
      CALL READMS(20,C,N,ICC)                                           CTHE2460
C                                                                       CTHE2470
      DO 41 I=NH1,NH2                                                   CTHE2480
      IF(YCG(I).GT.YLIMIT) GO TO 39                                     CTHE2490
      DO 37 J=NC1,NC2                                                   CTHE2500
 37   IF(XCG(I).GT.XCG(J).AND.XCG(I).LE.XCG(J+1)) GO TO 38              CTHE2510
 38   Y=1.-(YCG(I)-WDTHHT)/CSHT2                                        CTHE2520
      X=(XCG(J)-XCG(I))/(XCG(J+1)-XCG(J))                               CTHE2530
      GO TO 40                                                          CTHE2540
 39   Y=0.                                                              CTHE2550
      X=0.                                                              CTHE2560
 40   CMAX=C(J)+(C(J)-C(J+1))*X                                         CTHE2570
      C(I)=C(I)+C(N1)*(1.-Y)+CMAX*Y                                     CTHE2580
 41   CONTINUE                                                          CTHE2590
C                                                                       CTHE2600
C  REPLACE C MATRIX WITH HORIZONTAL TAIL ON HORIZONTAL TAIL EFFECTS     CTHE2610
 42   CALL WRITIN(20,C,N,ICC)                                           CTHE2620
      RETURN                                                            CTHE2630
   43 DO 44 J=NH1,NH2                                                   CTHE2640
      ICC=IC+K-1                                                        CTHE2650
C  READ C MATRIX FROM DISC                                              CTHE2660
      CALL READMS(20,C,N,ICC)                                           CTHE2670
      DO 44 I=NH1,NH2                                                   CTHE2680
      C(I)=C(I)+C(N1)                                                   CTHE2690
C                                                                       CTHE2700
C  REPLACE C MATRIX WITH ALL-MOVABLE HORIZONTAL TAIL EFFECTS            CTHE2710
   44 CALL WRITIN(20,C,N,ICC)                                           CTHE2720
      WRITE (6,45)                                                      CTHE2730
   45 FORMAT (///1X,"THIS AIRCRAFT HAS ALL-MOVABLE HORIZONTAL TAIL"///)  CTHE2740
      RETURN                                                            CTHE2750
      END                                                               CTHE2760
      SUBROUTINE CWRITE(C,DUM,N,IND)                                    CWRT0010
C...................................................................... CWRT0020
C  CTHETA CREATED AND STORED MATRIX C BY COLUMNS.  CWRITE READS         CWRT0030
C  C BACK INTO CORE AND STORES ON DISC BY ROWS FOR LATER USAGE.         CWRT0040
C...................................................................... CWRT0050
C                                                                       CWRT0060
      DIMENSION C(N,N),DUM(N)                                           CWRT0070
      DIMENSION IND(1)                                                  CWRT0080
      COMMON/INDMS/IA,IB,IC                                             CWRT0090
      DO 1 I=1,N                                                        CWRT0100
      ICC=IC+I-1                                                        CWRT0110
      CALL READMS(20,DUM(1),N,ICC)                                      CWRT0120
      DO 1 J=1,N                                                        CWRT0130
    1 C(J,I)=DUM(J)                                                     CWRT0140
      DO 2 J=1,N                                                        CWRT0150
      ICC=IC+J-1                                                        CWRT0160
      DO 3 I=1,N                                                        CWRT0170
    3 DUM(I)=C(J,I)                                                     CWRT0180
    2 CALL WRITIN(20,DUM(1),N,ICC)                                      CWRT0190
      RETURN                                                            CWRT0200
      END                                                               CWRT0210
      SUBROUTINE WRITEC(C,NCP)                                          WRTC0010
C                                                                       WRTC0020
C.......................................................................WRTC0030
C  WRITEC WRITES THE C-MATRIX TO TAPE7 AND TO THE OUTPUT FILE.          WRTC0040
C.......................................................................WRTC0050
      DIMENSION C(NCP,NCP)                                              WRTC0060
    1 FORMAT(1H1,47X,36HINFLUENCE COEFFICIENT MATRIX  C(J,K)//          WRTC0070
     .           6X,88HNOTATION FOR MATRIX C(J,K). J IS THE DEFLECTED PAWRTC0080
     .NEL. K IS THE PANEL WITH THE UNIT LOAD.//4X,7HJ    K=,6X,I3,7(12X,WRTC0090
     .I3))                                                              WRTC0100
    2 FORMAT(I5,8E15.5)                                                 WRTC0110
      REWIND 7                                                          WRTC0120
      K2=0                                                              WRTC0130
      NPAGE=NCP/8                                                       WRTC0140
      IF(MOD(NCP,8).NE.0) NPAGE=NPAGE+1                                 WRTC0150
      DO 3 I=1,NPAGE                                                    WRTC0160
      K1=K2+1                                                           WRTC0170
      K2=K2+8                                                           WRTC0180
      IF(K2.GT.NCP) K2=NCP                                              WRTC0190
      WRITE (6,1) (K,K=K1,K2)                                           WRTC0200
      DO 3 J=1,NCP                                                      WRTC0210
      WRITE (6,2) J,(C(J,K),K=K1,K2)                                    WRTC0220
    3 CONTINUE                                                          WRTC0230
C                                                                       WRTC0240
      DO 4 I=1,NCP                                                      WRTC0250
    4 WRITE(7) (C(I,J),J=1,NCP)                                         WRTC0260
      RETURN                                                            WRTC0270
      END                                                               WRTC0280
      SUBROUTINE CNVERT(EI,GJ,NN)                                       CNVT0010
C.......................................................................CNVT0020
C  CNVERT CONVERTS EI AND GJ VALUES FROM LB/IN**2 TO LB/FT**2.          CNVT0030
C.......................................................................CNVT0040
      DIMENSION EI(NN),GJ(NN)                                           CNVT0050
      C=1./144.                                                         CNVT0060
      DO 1 I=1,NN                                                       CNVT0070
      EI(I)=EI(I)*C                                                     CNVT0080
    1 GJ(I)=GJ(I)*C                                                     CNVT0090
      RETURN                                                            CNVT0100
      END                                                               CNVT0110
      SUBROUTINE WRITIN(LUN,A,N,IREC)                                   WRTN0010
      DIMENSION A(N)                                                    WRTN0020
      CALL WRITMS(LUN,A(1),N,IREC)                                      WRTN0030
      RETURN                                                            WRTN0040
      END                                                               WRTN0050
