
! POTENTIAL PROBLEMS WITH THIS PROGRAM
! Module subprogram name(ABLATE)
!  2005-W: "ablate.f90", line 1378: EMV is used but never set.
! Module subprogram name(DEE)
!  2005-W: "ablate.f90", line 2127: TEMATA is used but never set.
! Module subprogram name(EXEC2)
!  2005-W: "ablate.f90", line 2952: ALLNEW is used but never set.
! Module subprogram name(PBETA)
!  2005-W: "ablate.f90", line 4784: CG is used but never set.
!  2005-W: "ablate.f90", line 4784: GVAR is used but never set.




MODULE HeadG
IMPLICIT NONE  
  CHARACTER(LEN=80):: head
END Module HeadG   ! ========================================================

MODULE TimeDt
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
  REAL:: ctime
  INTEGER:: ktime   ! added by RLC  12 Nov 2005
  REAL:: time
END Module TimeDt   ! =======================================================

MODULE RevisedProcedures

CONTAINS

!
SUBROUTINE CRATIO(a,b) 
!     THIS IS A DUMMY SUBROUTINE AND SHOULD NOT BE ENTERED              
  REAL,INTENT(IN):: a,b
  WRITE(*,*) 'Stopping in Cratio, a,b=', a,b
  STOP 
END Subroutine Cratio

!
SUBROUTINE DEPOLY(a) 
!     THIS IS A DUMMY SUBROUTINE AND SHOULD NOT BE ENTERED
  REAL,INTENT(IN):: a
!----------------------------------------------------------------------------
  WRITE(*,*) 'Stopping in Depoly, a=', a  
  STOP 
END Subroutine Depoly   ! ---------------------------------------------------

!
SUBROUTINE MDUMP(k)
  INTEGER,INTENT(IN):: k
!     THIS IS A DUMMY SUBROUTINE AND SHOULD NOT BE ENTERED              
  WRITE(*,*) 'Stopping in Mdump, k=', k
  STOP 
END Subroutine Mdump

!
SUBROUTINE ErrorMessage(a)
! ---------------------------------------------------------------------------
  CHARACTER(LEN=*),INTENT(IN):: a
  REAL:: X = -1.0
!----------------------------------------------------------------------------
  WRITE(6,'(3X,A)') a
!     TO GET A SUBROUTINE BACK TRACE WE DO THIS                         
  X = SQRT(X) 
  STOP 
END Subroutine ErrorMessage   ! ---------------------------------------------


!+
SUBROUTINE InitializeCommon()
! ---------------------------------------------------------------------------
! PURPOSE - Initialize the values of certain COMMON variables
USE TimeDt

!$IBFTC INITBD
!!!      BLOCK DATA
!****                         INPUT     COMMON

      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)

      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)

      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25)

      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)

      COMMON  /INPTC /   PVAR ( 50)

      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)

      COMMON  /INPTC /   HORA ( 50),    TTABL(510)

      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)


!****                              COMPUTATION COMMON

      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4)

      COMMON  /COUNTR/   COUNT ( 25)

      COMMON  /HEATDP/   Q    (4,13)

      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11)

      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50)

      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6)

      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)

!!!      COMMON  /TIMEDT/   CTIME     ,    TIME

      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25)
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10)

      COMMON  /OUTDPT/  AUS(25)

!      COMMON  /HEADG /  HEAD(14)
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4)
!----------------------------------------------------------------------------
  eptab(1:7) = (/.2,.01,0.0001,0.0001,0.0001,.01,.0001/)
  agptb(1:9) = 0.0
  agptb(10:15) = (/.11,2.E5,.2,.071,.25,.33/)
  htconv(1:3) = (/1.9858,.476E-12,1.74532925E-2/)
END Subroutine InitializeCommon   ! -----------------------------------------

!
SUBROUTINE QBLK(QC,DHRHW,QBLOCK,EMV) 
! ---------------------------------------------------------------------------
USE TimeDt
  REAL,INTENT(IN):: qc,dhrhw
  REAL,INTENT(OUT):: qblock
  REAL,INTENT(IN):: emv

!****                         INPUT     COMMON

      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     

      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      

      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 

      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)

      COMMON  /INPTC /   PVAR ( 50) 

      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)

      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 

      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)


!****                              COMPUTATION COMMON
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 

      COMMON  /OUTDPT/  AUS(25) 

!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
!----------------------------------------------------------------------------
!****         THE EQUATIONS ARE DEFINED AND DISCUSSED IN SECTION II.C   
      IF ( CHART(1) ) 1, 100, 1 
    1 QDC = QC 
      HRHWD = DHRHW 
      EMDOT = EMV 
      II = KONVAR(6) * KONVAR(1) + 1 
      IF ( HRHWD ) 2, 3, 2 
    2 IF(QDC) 4, 3, 4 
    3 QBL = 0.0 
      GO TO 98 
    4 QCDH = QDC/HRHWD 
      IF(ABS(CHART(4))-CTIME) 7,7,5 
!****              IF CHART(4) IS MINUS, THEN WE HAVE POWERED FLIGHT,
!****                   SO TURBULENT COMES BEFORE LAMINAR.
    5 IF(CHART(4)) 50, 50, 10 
!****              IF CHART(4) IS PLUS, THEN WE HAVE RE-ENTRY,
!****                   SO LAMINAR COMES BEFORE TURBULENT.
    7 IF(CHART(4)) 10, 50, 50 
!****              LAMINAR BLOCKING ACTION STARTS AT EFN 10.
   10 EXPONL = .33333333 
      TERM = 0.69*(CHART(9)**EXPONL)*EMDOT/QCDH/(CHART(10)**EXPONL) 
      IF(TERM - .8) 30, 30, 20 
   20 TERM = .8 
   30 QBL = QDC*TERM 
!****              TURBULENT BLOCKING ACTION STARTS AT EFN 50.
      GO TO 98 
   50 EXPONT = -.38*CHART(11)*EMDOT/QCDH 
      QBL = QDC*(1. - EXP(EXPONT)) 
   98 QBLOCK = QBL 
  100 RETURN 
END Subroutine Qblk   ! -----------------------------------------------------


!+
SUBROUTINE QNET()
! ---------------------------------------------------------------------------
! PURPOSE -
!****         THE EQUATIONS USED IN THE QNET ROUTINE ARE DISCUSSED AND
!****         DEFINED IN SECTION II.C

! NOTES (RLC) - The original source code used qnet as both the name of the
!  subroutine and a local variable. This gets flagged as an error today.
!  I renamed the local variable qnetLocal.
USE TimeDt
!****                         INPUT     COMMON

      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)

      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)

      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25)

      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)

      COMMON  /INPTC /   PVAR ( 50)

      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)

      COMMON  /INPTC /   HORA ( 50),    TTABL(510)

      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)


!****                              COMPUTATION COMMON

      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4)

      COMMON  /COUNTR/   COUNT ( 25)

      COMMON  /HEATDP/   Q    (4,13)

      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11)

      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50)

      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6)

      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)

!!!      COMMON  /TIMEDT/   CTIME     ,    TIME

      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25)
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10)

      COMMON  /OUTDPT/  AUS(25)

!      COMMON  /HEADG /  HEAD(14)
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4)

!      DIMENSION THR(15),S(16)
!      DIMENSION HWTEMP(12),HWPRES(7),HWRTAB(12,7)


  REAL,PARAMETER,DIMENSION(15):: THR = (/ 0.,2000.,3000.,4000.,5000., &
    6000.,7000.,8000.,9000.,10000.,11000.,12000.,13000.,14000.,160000./)
  REAL,DIMENSION(16):: S = (/ 14.8,10.8,8.75,6.90,5.60,4.60, &
    3.95,3.50,3.10,2.85,2.65,2.52,2.42,2.38,2.38,1.0/)
  REAL,PARAMETER,DIMENSION(12):: HWTEMP = (/0.,3500.,3600.,4000.,4500., &
    5000.,5500.,6000.,6500.,7000.,7500.,8000./)
  REAL,PARAMETER,DIMENSION(7):: HWPRES = (/0.,.001,.01,.1,1.,10.,100./)
  REAL,PARAMETER,DIMENSION(12,7):: HWRTAB = RESHAPE( (/ &
    910.,910.,1020.,1290.,1900.,2670.,3100.,3330.,3550.,3900.,4500.,5500., &
    910.,910.,1020.,1290.,1900.,2670.,3100.,3330.,3550.,3900.,4500.,5500., &
    910.,910.,1000.,1170.,1500.,2100.,2770.,3210.,3450.,3700.,4000.,4480., &
    910.,910., 990.,1130.,1350.,1700.,2190.,2770.,3230.,3550.,3820.,4100., &
    910.,910., 990.,1120.,1310.,1540.,1840.,2240.,2680.,3140.,3540.,3840., &
    910.,910., 990.,1110.,1280.,1470.,1690.,1960.,2250.,2630.,3020.,3440., &
    910.,910., 990.,1110.,1280.,1450.,1650.,1860.,2070.,2310.,2590.,2900./), &
    (/12,7/) )
!****         ENTRY 16 IN THE S TABLE IS JUST FOR THE USE OF TABUP.
!****         THE 1.0 INDICATES THAT THE TABLE IS VARIABLE.
!----------------------------------------------------------------------------
      TEMPO = CTIME
      II = KONVAR(6)*KONVAR(1)+1
      TW = T(II)
      CALL TABUP(HORA,TTABL,TEMPO,TIVAR,INDTL(5),-8)
      KMELT = KONVAR(8)
      IF (GENRL(5)) 5,30,30
!****              COMPUTATION FOR QNET WITH FRONT FACE
!****                   T READ IN TTABL(1).
    5 CALL TabUpSingle(WARM,PROPT,TW,CAY,INDTL(4),2)
      TX = -3.*T(II)+4.*T(II+1)-T(II+2)
      qnetLocal = -CAY*SPACVR(3)/THIKN(1)*TX/2./DELETA(1)
      QRR = 0.0
      QHGR = 0.0
      QC = 0.0
      QBLOCK = 0.0
      GO TO 900
!****              QC OR QC/DH IS IN TTABL(1)
   30 CALL TABUP(WARM,HTABL,TW,HEVAR,INDTL(4),-3)
      QRR = HEVAR(2)*HTCONV(2)*TW**4
      IF (KMELT-4) 35,400,35
   35 QHGR = TIVAR(3)
      QC = TIVAR(1)
      HW = HEVAR(1)*TW
      HRMHW = TIVAR(2)-HW
      IF (GENRL(5)) 40,45,40
   40 QC = QC*HRMHW
   45 IF (CHART(1)) 65,60,65
   60 QBLOCK = 0.
      GO TO 90
   65 CALL QBLK(QC,HRMHW,QBLOCK,EMDG(1))
   90 qnetLocal = QC+QHGR-QRR-QBLOCK

      IF(KMELT) 95, 100, 95
   95 GO TO (100,200,300,400,500,600,700),KMELT
!****              EFN 100 STARTS THE FIXED CHAR LENGTH
!****                   OPTION.
!****         QNET IN THIS OPTION IS DEFINED IN II.C.1
  100 CONTINUE
      GO TO 900


!****              EFN 200 STARTS THE GRAPHITE SUBLIMATION
!****                   OPTION.
!****         QNET IN THIS OPTION IS DEFINED IN II.C.2
  200 HTCONV(6) = QC
      EXP1 = -11.05E4/TW
      QSTR = 3.96E8*TIVAR(4)**(-.67)*EXP(EXP1)
      HTCONV(10) = QSTR
!****              IF KONVAR(9) = 1, THE ROUTINE WAS CALLED
!****                   FROM THE ABLATE ROUTINE FOR QC.
      IF (KONVAR(9)) 998,220,998
  220 IF(ABS(CHART(4)) - TEMPO) 230, 230, 225
  225 IF(CHART(4)) 250, 250, 240
  230 IF(CHART(4)) 240, 250, 250
  240 CALL TabUpSingle(THR,S,TIVAR(2),HRS,15,1)
      GO TO 260
  250 HR = TIVAR(2)
      HRS = 26.781 + HR*(-7.9902E-3 + HR*(1.1199E-6 + HR*(-7.4989E-11 + &
     &HR*1.9445E-15)))
  260 QCP = QC*(1.-HRS*QSTR)
      Q(2,13) = QCP
      CALL QBLK(QCP,HRMHW,QBLOCK,EMDG(1))
      qnetLocal = QCP+QHGR-QRR-QBLOCK
      GO TO 900


!****              EFN 300 STARTS THE REFRASIL OPTION
!****         QNET IN THIS OPTION IS DEFINED IN II.C.3
  300 qnetLocal = QC+QHGR-QRR-QBLOCK-HTCONV(5)*TIVAR(6)
      GO TO 900


!****         EFN 400 STARTS THE GRAPHITE AT STAGNATION PT. OPTION.
!****         ASSUME CORRECT NAMES FOR CLARITY.
  400 HE = TIVAR(1)
      PE = TIVAR(2)

!*** changed S to s(1) in line below. RLC 7Sep05
  410 S(1) = 18.68-4.418E-3*HE + 3.945E-7*HE**2 + 1.146E-12*HE**3          &
     &-2.057E-15*HE**4 + 8.333E-20*HE**5
!****         TEST TO DETERMINE HOW  HW(AIR)  WILLBE COMPUTED.
  415 IF(TW.LE.3500.) GO TO 435
!****         HERE GET  HW(AIR)  FROM TABLE LOOK UP.
      PEINAT = PE/2116.
      CALL TabUpSingle(HWTEMP,HWRTAB(:,1),T(II),AIR1,12,1)   ! modified by RLC
      CALL TabUpSingle(HWPRES,HWRTAB(:,1),PEINAT,HWAIR,7,1)  ! modified by RLC
      GO TO 436
!****         HERE GET  HW(AIR)  DIRECTLY.
  435 HWAIR = .26*TW
  436 IF(TMELT(8))  460,437,460
  437 RNOSE = TMELT(6)
      IF(TMELT(7))  440,450,440
!****         TEST TO SEE HOW  QBL  WILL BE COMPUTED.
  440 CHLR = TIVAR(4)
      GO TO 455
!****         HERE COMPUTE  QBL  DIRECTLY BUT FIRST GET  CHLR.
  450 CALL CRATIO(TRAJV(18),CHLR)
  455 QBL = .037*(HE-HWAIR)*CHLR*TMELT(5)*SQRT(PE/2116.*1./RNOSE)
      GO TO 465
!****         GET  QBL  FROM TIVAR(3) WHICH CAME FROM T-TABLE  3.
  460 QBL = TIVAR(3)
  465 QSUB = QBL*S(1)*2.4E6*(1./(PE/2116.)**.67)*EXP(-1.11E5/TW) ! s(1) was s
  475 qnetLocal = QBL - QSUB - QRR
      QHGR = QBL
      QC = QSUB
!****         THAT ENDS THIS OPTION.
      GO TO 900


!****              EFN 500 STARTS THE FIXED MELTING T OPTION.
!****              CHECK IN PTX1 FOR RHO*SPM TERM
!****         QNET IN THIS OPTION IS DEFINED IN II.C.1.A
  500 GO TO 100


  600 CALL TabUpSingle(WARM,GASPR,TW,HCG,INDTL(4),3)
      CALL TabUpSingle(WARM,HTABL,TW,HEVAR(1),INDTL(4),1)   ! modified by RLC
!****          HTCONV(5) = SPM
  625 RTVAP = HTCONV(5)*HCG
  650 DMDG =  DENSE(1)*HTCONV(5)
      QHGR = TIVAR(3)
      QC = TIVAR(1)
  655 HW = HEVAR(1)*TW
      HRMHW = TIVAR(2) - HW
  660 IF(GENRL(5)) 665,670,665
  665 QC = QC*HRMHW
!****          PLACE DUMMY NUM. IN CHART(1) FOR SUB. QBLK.
  670 CHART(1) = 1.0
  675 CALL QBLK(QC,HRMHW,QBLOCK,DMDG)
!****          REMOVE DUMMY NUM. FROM CHART(1) MIGHT GET SCREWED OTHERWISE
      CHART(1) = 0.0
      qnetLocal = qnetLocal-RTVAP-QBLOCK
      GO TO 900


!****    EFN 700 STARTS A COMBINATION GRAPHITE SUBLIMATION AND
!****    FIXED MELTING OPTION
  700 IF(TW-TMELT(3)) 200,200,500


  900 KT = KONVAR(2)+1
      Q(KT,1) = qnetLocal
      Q(KT,2) = QRR
      Q(KT,3) = QHGR
      Q(KT,4) = QC
      Q(KT,5) = QBLOCK
  998 RETURN
END Subroutine Qnet   ! -----------------------------------

!+
      SUBROUTINE TableOutput()
! ---------------------------------------------------------------------------
! PURPOSE -!     THIS IS THE HEADING SECTION

USE HeadG,ONLY: head
USE TimeDt
!****                         INPUT     COMMON

      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)

      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)

      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25)

      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)

      COMMON  /INPTC /   PVAR ( 50)

      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)

      COMMON  /INPTC /   HORA ( 50),    TTABL(510)

      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)


!****                              COMPUTATION COMMON

      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4)

      COMMON  /COUNTR/   COUNT ( 25)

      COMMON  /HEATDP/   Q    (4,13)

      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11)

      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50)

      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6)

      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)

!!!      COMMON  /TIMEDT/   CTIME     ,    TIME

      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25)
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10)

      COMMON  /OUTDPT/  AUS(25)

!      COMMON  /HEADG /  HEAD(14)
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4)

!************COMMON FOR VARIABLE T-P-D INPUT*******FOR BETA-4 TABLE****
      DIMENSION TMDENS(1),TMT(1),TMP(1),XFORT(1),XFORP(1),XFORDN(1)
!     ABOVE IS DUMMY DIMENSION -JUST TO ALLOW SUBSCRIBTING
      EQUIVALENCE(TMDENS(1),DENS(301)),(TMT(1),T(301)),(TMP(1),P(301))
      EQUIVALENCE(XFORT(1),T(151)),(XFORP(1),P(151)),(XFORDN(1),DENS(151&
     &))
      EQUIVALENCE(NT,T(449)),(NXT,T(450)),(NP,P(449)),(NXP,P(450))
      EQUIVALENCE(ND,DENS(449)),(NXD,DENS(450))
!     THE T-P-AND DENS INPUTS ARE STORED IN THE
!     LAST THIRD OF THEIR RESPECTIVE ARRAYS.
!     NT-NP-NDENS ETC ARE STORED IN LOC
!     NOT USED -DUE TO THE IMPOSSIBILITY OF HAVING 150 POINTS
!     THE X TABLE FOR T(XFORT)ETC ARE STORED IN THE
!     MIDDLE THIRD OF THEIR RESPECTIVE ARRAYS
!     FOR THE BETA 4 TABLE HAVING PSEUDO - VARIABLE DIMENSION

      COMMON TXQTAB(200),NOFTS,NOFXS,NTXQ(50),ICON,JCON,END
!!!      DIMENSION FMT(13)    ! 78 characters
      CHARACTER(LEN=6),DIMENSION(13):: fmt

  CHARACTER(LEN=6),PARAMETER:: ANCT = 'CONST.', ANVAR = 'VARIB.'
  REAL,PARAMETER:: ZERO = 0.0
!----------------------------------------------------------------------------
  WRITE(6,*) ACHAR(12)
  WRITE(6,'(1X,A)') head
!      WRITE (6,1500) (HEAD(I), I = 2,13)
      LAYER = KONVAR(4)
      JGAP = KONVAR(5)
      NTIME = INDTL(5)
      NTEMP = INDTL(4)
      NDENS = INDTL(1)
      WRITE (6,1501)
!****         SECTION FOR THE PTSIN, DENSE, AND THIKN TABLES.
      WRITE (6,1542)
      DO K = 1,LAYER
        WRITE(6,1543)  K,PTSIN(K),THIKN(K),DENSE(K)
      END DO

      WRITE (6,1611)
!****         SECTION FOR TEMPERATURE DEPENDENT TABLES.
!****         WRITE OUT THE APPROPRIATE HEADINGS.
      WRITE (6,1608)
      IF(STVAR(1) <= ZERO) THEN   ! 12, 12, 15
     12 DO J = 1,LAYER
        WRITE (6,1502) J
        WRITE (6,1505)
        DO I = 1,NTEMP
          MN = (NTEMP+1)*(J-1)+I
          MN1 = (NTEMP+1)*(3*(J-1))+I
          MN2 = (NTEMP+1)*(3*(J-1)+1)+I
          MN3 = (NTEMP+1)*(3*(J-1)+2)+I
   14     WRITE(6,1506)WARM(I),PROPT(MN1),PROPT(MN2),PROPT(MN3),CTABL(MN)
        END DO
        END DO
      ELSE

   15   DO J = 1,LAYER ! TEMP. DEPENDENT MATERIAL AND STRESS PROPERTY TABLES.
          WRITE (6,1502) J
          WRITE(6,1563)
          NDGR1 = NDGREE + 1
          DO I = 1,NTEMP
            MN = (NTEMP+1)*(J-1)+I
            MN1 = (NTEMP+1)*(3*(J-1))+I
            MN2 = (NTEMP+1)*(3*(J-1)+1)+I
            MN3 = (NTEMP+1)*(3*(J-1)+2)+I
            MD = NDGR1 * (4*(J-1)) + I
            MD3 = NDGR1 * (4*(J-1)+2) + I
            WRITE(6,1564) WARM(I), PROPT(MN1),PROPT(MN2),PROPT(MN3), &
             CTABL(MN),ATETB(MD),ATETB(MD3)
          END DO
        END DO
        WRITE(6,1567)STVAR(1),STVAR(2),STVAR(3)
      END IF


!****         TEMP. DEPENDENT GAS PROPERTY TABLES.
   28 WRITE(6,1508)
      DO I = 1,NTEMP
      MN1 = I
      MN2 = NTEMP+1+I
      MN3 = MN2+NTEMP+1
      MN4 = MN3+NTEMP+1
      MN5=MN4+NTEMP+1
   39 WRITE(6,1509) WARM(I),GASPR(MN1),GASPR(MN2),GASPR(MN3),GASPR(MN4), &
       GASPR(MN5),HTABL(MN1),HTABL(MN2)
      END DO

!****         SECTION TO TEST AND WRITE OUT IF TABLES ARE
!****              CONSTANT OR VARIABLE.
      N = NTEMP+1
      DO 360 J = 1,LAYER
!!!      J = J
      MN = J*N
      MD = J * NDGR1
      NDG3 = 3 * MD
      DO I = 1,3
      MN1 = I*J*N
      IF (PROPT(MN1)==ZERO) THEN    ! 308,304,308
  304 FMT(I) = ANCT
      ELSE   ! GO TO 310
  308 FMT(I) = ANVAR
      END IF
      END DO

      DO I = 1,5
        MN1 = I*J*N
        IF (GASPR(MN1)==ZERO) THEN   ! 318,314,318
          FMT(I+3) = ANCT
        ELSE
          FMT(I+3) = ANVAR
        END IF
      END DO

      DO I = 1,2
        MN1 = I*J*N
        IF (HTABL(MN1)==ZERO) THEN   ! 326,322,326
          FMT(I+8) = ANCT
          JK2 = JK1+N
          JK3 = JK2+N
        ELSE
          FMT(I+8) = ANVAR
        END IF
      END DO

      IF (CTABL(MN)==ZERO) THEN   ! 336,332,336
        FMT(11) = ANCT
      ELSE
        FMT(11) = ANVAR
      END IF

  340 IF (STVAR(1) /= ZERO) THEN   ! 341,358,341
  341   IF (ATETB(MD)==ZERO) THEN   ! 346, 342, 346
      342 FMT(12) = ANCT
        ELSE
  346 FMT(12) = ANVAR
        END IF
  350 IF (ATETB(NDG3)==ZERO) THEN   ! 356, 352, 356
  352 FMT(13) = ANCT
      ELSE
  356 FMT(13) = ANVAR
      END IF
      END IF

  358 WRITE(6,1600)FMT(1),FMT(2),FMT(3),FMT(11)
      WRITE(6,1601)FMT(4),FMT(5),FMT(6),FMT(7)
      WRITE(6,1602)FMT(8),FMT(9),FMT(10)
      IF (STVAR(1)/=ZERO) WRITE(6,1603)FMT(12),FMT(13)
  360 END DO

      WRITE(6,1611)
      WRITE(6,1609)
!****         SECTION FOR TIME DEPENDENT TABLES.
!****         WRITE OUT THE APPROPRIATE HEADINGS.
      N = NTIME+1
  DO  K = 1,5,4   ! once with k=1; again with k=5
    IF (K < 5) THEN   ! 430, 431, 431
      WRITE(6,1531)
      KZ = 0
    ELSE
      WRITE(6,1532)
      KZ = 5*N
    END IF

    DO KL = 1,NTIME
      JK1 = KZ + KL
      JK2 = JK1 + N
      JK3 = JK2 + N
      JK4 = JK3+N
      JK5 = JK4+N
      WRITE(6,1533) HORA(KL), &
        TTABL(JK1),TTABL(JK2),TTABL(JK3),TTABL(JK4),TTABL(JK5)
    END DO
  END DO

!*****        SECTION TO TEST AND WRITE OUT IF THE TABLES ARE
!****              CONSTANT OR VARIABLE.
  N = NTIME+1
  DO I = 1,10
    MN1 = I*N
    IF (TTABL(MN1)==ZERO) THEN   ! 466,462,466
      FMT(I) = ANCT
    ELSE
      FMT(I) = ANVAR
    END IF
  END DO

  WRITE(6,1604) FMT(1),FMT(2),FMT(3),FMT(4)
  WRITE(6,1605) FMT(5),FMT(6),FMT(7),FMT(8)
  WRITE(6,1606) FMT(9),FMT(10)
  WRITE(6,1611)

!****         SECTION FOR THE DENSITY DEPENDENT TABLES.
!****         WRITE OUT THE APPROPRIATE HEADINGS.
  IF (CHART(1) /= ZERO) THEN   ! 512, 616, 512
    WRITE(6,1610)
    WRITE(6,1512)
    DO  L = 1,NDENS
      IZE2 = (NDENS+1)+L
      IZE3 = IZE2+NDENS+1
      IZE4 = IZE3+NDENS+1
      WRITE(6,1513)TDENS(L),ZORET(L),ZORET(IZE2),ZORET(IZE3),ZORET(IZE4)
    END DO
!*****        SECTION TO TEST AND WRITE OUT IF THE TABLES ARE
!****              CONSTANT OR VARIABLE.
    N = NDENS+1
    DO I = 1,2
      MN1 = I*N
      IF (ZORET(MN1)==ZERO) THEN   ! 576,572,576
        FMT(I) = ANCT
      ELSE
        FMT(I) = ANVAR
      END IF
    END DO
    WRITE(6,1607)FMT(1),FMT(2)
    WRITE(6,1611)
  END IF


!****         SECTION TO WRITE OUT THE GENRL TABLE.
  616 WRITE(6,1546)( GENRL(I), I = 1,15)
      WRITE(6,1611)
!****         SECTION TO WRITE OUT THE DTIME TABLE.
  639 WRITE(6,1539)(DTIME(K),K = 1,20)
      WRITE(6,1611)
!****         SECTION TO WRITE OUT THE EPTAB TABLE.
      WRITE(6,1550)(EPTAB(I),I = 1,10)
      WRITE(6,1611)
!****         SECTION TO WRITE OUT THE CHART TABLE.
  649 WRITE(6,1549)(CHART(J),J = 1,20)
      WRITE(6,1611)
!****         SECTION TO WRITE OUT THE TMELT TABLE.
      WRITE(6,1569)(TMELT(J),J = 1,18)
      WRITE(6,1611)
!****         SECTION FOR THE OUTTB TABLE.
  653 WRITE(6,1553)(OUTTB(J),J = 1,25)
      WRITE(6,1611)
!****         SECTION TO PRINT OUT THE SPACE OPTION.
  652 WRITE(6,1552)(SPACE(J),J=1,10)
      WRITE(6,1611)

!****         SECTION TO WRITE OUT THE AGPTB TABLE.
  IF (JGAP /= ZERO) THEN   ! 690,700,690
    WRITE(6,1559)
    WRITE(6,1560)JGAP,(AGPTB(J),J=2,6)
    WRITE(6,1561)(AGPTB(J),J=10,15)
    WRITE(6,1562)(AGPTB(J),J = 7,9)
    WRITE(6,1611)
  END IF

!****         SECTION TO WRITE OUT BETA 4 CORRECTION TABLE.
  IF (NOFTS /= 0) THEN   ! 710,725,710
    WRITE(6,1620)
    WRITE(6,1625) (TXQTAB(JD), JD = 1,NOFTS)
    WRITE(6,1611)
    LIMOUT = NOFTS + 1
    LIMOT1 = NOFTS + NOFXS
    WRITE(6,1630) (TXQTAB(JD), JD = LIMOUT,LIMOT1)
    WRITE(6,1611)
    DO JD = 2,NOFXS
      LIMOUT = (JD-1)*NOFXS +NOFTS + 1
      LIMOT1 = JD*NOFXS + NOFTS
      WRITE(6,1635) (TXQTAB(JDH) , JDH = LIMOUT,LIMOT1)
    END DO
    WRITE(6,1611)
  END IF

!****         SECTION TO WRITE OUT VARIABLE T'P'DENSITY INPUT)
  IF (GENRL(4)==ZERO) THEN   ! 733,730,733
    WRITE(6,1640)
    WRITE(6,1645) (XFORT(JD) , JD = 1,NXT)
    WRITE(6,1611)
    WRITE(6,1650) (TMT(JD) , JD = 1,NT)
    WRITE(6,1611)
    WRITE (6,1655) (XFORP(JD) ,JD = 1,NXP)
    WRITE(6,1611)
    WRITE(6,1660) (TMP(JD) , JD = 1,NP)
    WRITE(6,1611)
    WRITE (6,1665) (XFORDN(JD) , JD = 1,ND)
    WRITE(6,1611)
    WRITE (6,1670) (TMDENS(JD) ,JD = 1,NXD)
    WRITE(6,1611)
  END IF


  IF(AMPLCT(2) /= ZERO) THEN   ! 735,740,735
    WRITE (6,1611)
    WRITE (6,1675) (AMPLCT(I),I = 1,10)
  END IF

  IF(BACKFC(1) /= ZERO) THEN   ! 745,1000,745
    WRITE(6,1611)
    WRITE (6,1680) (BACKFC(I) ,I= 1,4)
  END IF

  RETURN

!****         SECTION FOR FORMAT STATEMENTS.
 1500 FORMAT(1H1,12A6)
 1501 FORMAT (/,55X,'INPUT TABLES')
 1502 FORMAT (/,33X,'LAYER ',I2)
 1505 FORMAT (/,3X, &
   'TEMP',9X,'CP TAB',7X,'K TAB',8X,'DKDT TAB',5X,'C TAB')
 1506 FORMAT(5E13.4)
 1508 FORMAT (/,3X,'TEMP',9X,'HGF',10X,'CPG',10X,'HCG',9X, &
   'MBAR',10X, 'MBARP',8X,'CPBL',9X,'EFF')
 1509 FORMAT(8E13.4)
 1512 FORMAT (/,3X,'DENSITY TABLE',4X,'Z',13X,'E',13X,'Z1',12X,'E1')
 1513 FORMAT(5E13.4)
 1531 FORMAT (/,3X,'TIME',9X,'TTABL(1)',6X,'TTABL(2)',6X, &
    'TTABL(3)',6X, 'TTABL(4)',6X,'TTABL(5)')
 1532 FORMAT (/,3X,'TIME',9X,'TTABL(6)',6X,'TTABL(7)',6X, &
   'TTABL(8)',6X, 'TTABL(9)',6X,'TTABL(10)')
 1533 FORMAT(7E13.4)
 1539 FORMAT(/,3X,'DTIME TABLE = ',5E18.5)
 1542 FORMAT (/,3X, &
   'LAYER',5X,'INTERIOR PTS',4X,'THICKNESS',7X,'DENSITYS')
 1543 FORMAT (I6,F17.0,E19.3,E15.3)
 1546 FORMAT(/,3X,'GENRL TABLE = ',5E18.5)
 1549 FORMAT (/,3X,'CHART TABLE = ',5E18.5)
 1550 FORMAT (/,3X,'EPTAB TABLE = ',5E18.5)
 1552 FORMAT(/,3X,'SPACE TABLE = ',5E18.5)
 1553 FORMAT (/,3X,'OUTTB TABLE = ',5E18.5)
 1559 FORMAT (/,45X,'AIRGAP PROPERTIES')
 1560 FORMAT(1X,'LAYER=',I2,3X,'FE=',E10.3,3X,'FA=',E10.3,3X, &
       'EFF COND=',E10.3,3X,'VISC=',E10.3,3X,'HEIGHT=',E10.3)
 1561 FORMAT(' EXP. OF X/L=',E10.3,2X,'TRIP NO.=',E10.3,2X,'COEF(-,+)=(', &
       E10.3,',',E10.3,')',2X,'EXP(-,+)=(',E10.3,',',E10.3,')')
 1562 FORMAT(2X,'G1 = ',E15.8,3X,'G2 = ',E15.8,3X,'G3 = ',E15.8)
 1563 FORMAT (/,3X,'TEMP',9X,'CP TAB',7X,'K TAB',8X,'DKDT TAB',5X, &
      'C TAB',8X,'ALPHA T',7X,'E-MOD(IN.)')
 1564 FORMAT(7E13.4)
 1567 FORMAT(/,3X,'STVAR(1) = ',E15.8,4X,'T(REF) = ',E15.8,4X, &
       'T(AB) = ',E15.8)
 1568 FORMAT(/,3X,'TRAJT = ',5E18.5)
 1569 FORMAT (/,3X,'TMELT TABLE = ',5E18.5)
 1600 FORMAT (/,3X,'PROPT(1) = ',1A6,3X,'PROPT(2) = ',1A6,3X, &
      'PROPT(3) = ',1A6,3X,'CTABL(1) = ',1A6)
 1601 FORMAT (/,3X,'GASPR(1) = ',1A6,3X,'GASPR(2) = ',1A6,3X, &
       'GASPR(3) = ',1A6,3X,'GASPR(4) = ',1A6)
 1602 FORMAT (/,3X,'GASPR(5) = ',1A6,3X,'HTABL(1) = ',1A6,3X, &
       'HTABL(2) = ',1A6)
 1603 FORMAT (/,3X,'ALPHT(1) = ',1A6,3X,'EMOD(1) = ',1A6)
 1604 FORMAT (/,3X,'TTABL(1) = ',1A6,3X,'TTABL(2) = ',1A6,3X, &
      'TTABL(3) = ',1A6,3X,'TTABL(4) = ',1A6)
 1605 FORMAT (/,3X,'TTABL(5) = ',1A6,3X,'TTABL(6) = ',1A6,3X, &
      'TTABL(7) = ',1A6,3X,'TTABL(8) = ',1A6)
 1606 FORMAT (/,3X,'TTABL(9) = ',1A6,3X,'TTABL(10) = ',1A6)
 1607 FORMAT (/,3X,'ZORET(1) = ',1A6,3X,'ZORET(2) = ',1A6)
 1608 FORMAT(/,'***************',3X,'TEMP. DEP. TABLES',3X,'***************')
 1609 FORMAT(/,'***************',3X,'TIME  DEP. TABLES',3X,'***************')
 1610 FORMAT(/,'***************',3X,'DENS. DEP. TABLES',3X,'***************')
 1611 FORMAT(/,'*************************')
 1620 FORMAT(/,'***************',3X,'BETA 4 .Q. TABLES',3X,'***************')
 1625 FORMAT(/,3X,'IND T TABLE =',5E18.5)
 1630 FORMAT(/,3X,'IND X TABLE =',5E18.5)
 1635 FORMAT(/,3X,'DEP Q TABLE =',5E18.5)
 1640 FORMAT(/,'***************',3X,'T  P  DENS TABLES',3X,'***************')
 1645 FORMAT(/,3X,'XFORT TABLE =',5E18.5)
 1650 FORMAT(/,3X,'TEMPS TABLE =',5E18.5)
 1655 FORMAT(/,3X,'XFORP TABLE =',5E18.5)
 1660 FORMAT(/,3X,'PRESS TABLE =',5E18.5)
 1665 FORMAT(/,3X,'XFORD TABLE =',5E18.5)
 1670 FORMAT(/,3X,'DENS. TABLE =',5E18.5)
 1675 FORMAT(/,3X,'AMPLCT TABLE=',5E18.5)
 1680 FORMAT(/,'BACKFC TABLE=',4E18.5)
END Subroutine TableOutput   ! ----------------------------------------------

!+
SUBROUTINE TABUP(TABIND,TABDEP,VARIND,VARDEP,NENTRY,KIN)
! ---------------------------------------------------------------------------
! PURPOSE -
USE TimeDt
  REAL,INTENT(IN),DIMENSION(:):: tabind,tabdep
  REAL,INTENT(IN):: varind
  REAL,INTENT(OUT),DIMENSION(:):: vardep
  INTEGER,INTENT(IN):: nentry
  INTEGER,INTENT(IN):: kin



!****                         INPUT     COMMON

      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)

      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)

      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25)

      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)

      COMMON  /INPTC /   PVAR ( 50)

      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)

      COMMON  /INPTC /   HORA ( 50),    TTABL(510)

      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)


!****                              COMPUTATION COMMON

      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4)

      COMMON  /COUNTR/   COUNT ( 25)

      COMMON  /HEATDP/   Q    (4,13)

      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11)

      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50)

      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6)

      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)

!!!      COMMON  /TIMEDT/   CTIME     ,    TIME

      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25)
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10)

      COMMON  /OUTDPT/  AUS(25)

!      COMMON  /HEADG /  HEAD(14)
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4)

!                    THIS SUBROUTINE PERFORMS A TABLE LOOKUP
!                    ON ANY ONE OR ANY NUMBER OF TABLES IN A GIVEN ARRAY
!                    EXPLANATION OF CALLING SEQUENCE
!                    TABIND IS THE INDEPENDENT TABLE ARRAY
!                    TABDEP IS THE DEPENDENT TABLE ARRAY
!                    VARIND IS THE INDEPENDENT VARIABLE
!                    VARDEP WILL CONTAIN THE ANSWER OR ANSWERS
!                    RESULTING FROM THE TABLE LOOK'UP
!                    NENTRY CONTAINS THE NUMBER OF VALUES IN THE TABLE
!                    KIN IS THE CONTROL FOR THE NUMBER
!                    OF TABLE LOOKUPS
!                    A NEGATIVE SIGN MAY MEAN MORE THAN ONE TABLE
!                    LOOK'UP.  NUMERICAL VALUE OF KIN WILL GIVE THE
!                    NUMBER OF TABLES GOING FROM THE FIRST
!                    ON WHICH TO PERFORM A TABLE LOOK'UP
!                    A POSITIVE SIGN MEANS ONE TABLE LOOK'UP AND THE
!                    NUMERICAL VALUE OF KIN WILL GIVE WHICH
!                    TABLE IN THE ARRAY THAT IS TO BE USED
!      DIMENSION VARDEP(1),TABIND(1),TABDEP(1)
!----------------------------------------------------------------------------
      KAY=KIN
!                    TEST FOR NUMBER OF TABLE LOOK'UPS
      IF(KAY) 1,6,2
!                    MORE THAN ONE TABLE LOOK'UP DESIRED
    1 LL=-KAY
      KAY=1
      GO TO 3
!                    ONLY ONE TABLE LOOK'UP
    2 LL=1
!                    TEST IF BELOW TABLE RANGE
    3 IF(VARIND-TABIND(1)) 5,4,7
!                    EXACT VALUE
    4 JC=1
      GO TO 9
    5 CALL ErrorMessage('TABLE LOOKUP OUT OF RANGE')
    6 CALL ErrorMessage('KIN HAS A ZERO ARGUMENT')
    7 DO 8 JC=2,NENTRY
!!!      JC=JC
!                    FIND LOCATION OF VARIABLE IN THE
!                    INDEPENDENT TABLE
      IF(VARIND-TABIND(JC)) 10,9,8
    8 END DO
      GO TO 5
    9 MM=1
      GO TO 11
!                    INTERPOLATION NECESSARY
   10 MM=2
   11 DO 16 J=1,LL
      LLL=(NENTRY+1)*J*KAY
!                    TEST IF TABLE IS CONSTANT
      IF(TABDEP(LLL)) 13,12,13
!                    CONSTANT TABLE
   12 VARDEP(J)=TABDEP(LLL-1)
      GO TO 16
!                    N IS SUBSCRIPT TO LOCATE TABLE
   13 N=(1+NENTRY)*(KAY*J-1)+JC
      GO TO (14,15),MM
!                    EXACT VALUE
   14 VARDEP(J)=TABDEP(N)
      GO TO 16
!                    INTERPOLATION ROUTINE
   15 M=N-1
      DIV=TABIND(JC)-TABIND(JC-1)
      TERM1=VARIND*(TABDEP(N)-TABDEP(M))
      TERM2=TABIND(JC)*TABDEP(M)-TABIND(JC-1)*TABDEP(N)
      VARDEP(J)=(TERM1+TERM2)/DIV
   16 END DO
   20 RETURN
END Subroutine Tabup

!+
SUBROUTINE TabUpSingle(TABIND,TABDEP,VARIND,VARDEP,NENTRY,KIN)
! ---------------------------------------------------------------------------
! PURPOSE -
USE TimeDt
  REAL,INTENT(IN),DIMENSION(:):: tabind,tabdep
  REAL,INTENT(IN):: varind
  REAL,INTENT(OUT):: vardep
  INTEGER,INTENT(IN):: nentry
  INTEGER,INTENT(IN):: kin



!****                         INPUT     COMMON

      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)

      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)

      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25)

      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)

      COMMON  /INPTC /   PVAR ( 50)

      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)

      COMMON  /INPTC /   HORA ( 50),    TTABL(510)

      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)


!****                              COMPUTATION COMMON

      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4)

      COMMON  /COUNTR/   COUNT ( 25)

      COMMON  /HEATDP/   Q    (4,13)

      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11)

      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50)

      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6)

      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)

!!!      COMMON  /TIMEDT/   CTIME     ,    TIME

      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25)
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10)

      COMMON  /OUTDPT/  AUS(25)

!      COMMON  /HEADG /  HEAD(14)
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4)

!                    THIS SUBROUTINE PERFORMS A TABLE LOOKUP
!                    ON ANY ONE OR ANY NUMBER OF TABLES IN A GIVEN ARRAY
!                    EXPLANATION OF CALLING SEQUENCE
!                    TABIND IS THE INDEPENDENT TABLE ARRAY
!                    TABDEP IS THE DEPENDENT TABLE ARRAY
!                    VARIND IS THE INDEPENDENT VARIABLE
!                    VARDEP WILL CONTAIN THE ANSWER OR ANSWERS
!                    RESULTING FROM THE TABLE LOOK'UP
!                    NENTRY CONTAINS THE NUMBER OF VALUES IN THE TABLE
!                    KIN IS THE CONTROL FOR THE NUMBER
!                    OF TABLE LOOKUPS
!                    A NEGATIVE SIGN MAY MEAN MORE THAN ONE TABLE
!                    LOOK'UP.  NUMERICAL VALUE OF KIN WILL GIVE THE
!                    NUMBER OF TABLES GOING FROM THE FIRST
!                    ON WHICH TO PERFORM A TABLE LOOK'UP
!                    A POSITIVE SIGN MEANS ONE TABLE LOOK'UP AND THE
!                    NUMERICAL VALUE OF KIN WILL GIVE WHICH
!                    TABLE IN THE ARRAY THAT IS TO BE USED
!      DIMENSION VARDEP(1),TABIND(1),TABDEP(1)
!----------------------------------------------------------------------------
      KAY=KIN
!                    TEST FOR NUMBER OF TABLE LOOK'UPS
      IF(KAY) 1,6,2
!                    MORE THAN ONE TABLE LOOK'UP DESIRED
    1 LL=-KAY
      KAY=1
      GO TO 3
!                    ONLY ONE TABLE LOOK'UP
    2 LL=1
!                    TEST IF BELOW TABLE RANGE
    3 IF(VARIND-TABIND(1)) 5,4,7
!                    EXACT VALUE
    4 JC=1
      GO TO 9
    5 CALL ErrorMessage('TABLE LOOKUP OUT OF RANGE')
    6 CALL ErrorMessage('KIN HAS A ZERO ARGUMENT')
    7 DO 8 JC=2,NENTRY
!!!      JC=JC
!                    FIND LOCATION OF VARIABLE IN THE
!                    INDEPENDENT TABLE
      IF(VARIND-TABIND(JC)) 10,9,8
    8 END DO
      GO TO 5
    9 MM=1
      GO TO 11
!                    INTERPOLATION NECESSARY
   10 MM=2
   11 DO 16 J=1,LL
      LLL=(NENTRY+1)*J*KAY
!                    TEST IF TABLE IS CONSTANT
      IF(TABDEP(LLL)) 13,12,13
!                    CONSTANT TABLE
   12 VARDEP=TABDEP(LLL-1)
      GO TO 16
!                    N IS SUBSCRIPT TO LOCATE TABLE
   13 N=(1+NENTRY)*(KAY*J-1)+JC
      GO TO (14,15),MM
!                    EXACT VALUE
   14 VARDEP=TABDEP(N)
      GO TO 16
!                    INTERPOLATION ROUTINE
   15 M=N-1
      DIV=TABIND(JC)-TABIND(JC-1)
      TERM1=VARIND*(TABDEP(N)-TABDEP(M))
      TERM2=TABIND(JC)*TABDEP(M)-TABIND(JC-1)*TABDEP(N)
      VARDEP=(TERM1+TERM2)/DIV
   16 END DO
   20 RETURN
END Subroutine TabupSingle




END Module RevisedProcedures   ! ============================================


!+
MODULE AblateProcedures

USE HeadG,ONLY: head
USE RevisedProcedures

CONTAINS


      SUBROUTINE ABLATE (H1) 
USE TimeDT
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATIONS USED FOR THE RATES OF MELT ARE DEFINED IN   
!****         SECTION II, PART C.                                       
      H = H1 
!****              IF H = 0.,CALCULATE SPC AND SPM.                     
!****              IF H = NUMBER, COMPUTE SC AND SM.                    
      KT = KONVAR(2) 
      IF (TMELT(1)) 5,4,5 
    4 SPM = 0.0 
      SPC = 0.0 
      GO TO 900 
    5 KMELT = KONVAR(8) 
      CTIME = CTIME 
      II = KONVAR(6)*KONVAR(1)+1 
      TW = T(II) 
      IF (H) 800,10,800 
!****         SPC IS DEFINED IN SECTION II, PARTC, SUB-PART 1,B.        
   10 SPC = EMDG(1)/(DENSE(1)-CHART(2)) 
      IF(GENRL(9)) 20,30,20 
   20 SPC = TRAJV(50)*SPC/(TRAJV(50) -HTCONV(8)) 
   30 GO TO (100,200,300,400,500,600,700),KMELT 
!****              EFN 100 STARTS THE FIXED CHAR                        
!****                   LENGTH OPTION.                                  
  100 CALL TabUpSingle(HORA,TTABL,CTIME,SCM,INDTL(5),7) 
      SX = (SCHAR(KT+1)-SMELT(KT+1)) - SCM 
      IF (SX) 115,115,118 
  115 SPM = 0. 
      GO TO 900 
  118 CALL TabUpSingle(HORA,TTABL,CTIME,SPCM,INDTL(5),8) 
      SPM = SPC - SPCM 
      GO TO 900 
!****              EFN 200 STARTS THE GRAPHITE                          
!****                   SUBLIMATION OPTION.                             
  200 IF (TW-1800.) 210,215,215 
  210 SPM = 0. 
      GO TO 900 
  215 KONVAR(9) = 1 
      CALL QNET 
      KONVAR(9) = 0 
!****              TIVAR VALUES ARE COMPUTED IN QNET.                   
      CALL TabUpSingle(WARM,HTABL,TW,CPBL,INDTL(4),1) 
      DMD = HTCONV(6)/(TIVAR(7)+TIVAR( 8)*(TIVAR(2)-CPBL*TW)) 
      EXP1 = -11.05E4/TW 
  220 IF(ABS(CHART(4)) - CTIME) 230, 230, 225 
  225 IF(CHART(4)) 250, 250, 240 
  230 IF(CHART(4)) 240, 250, 250 
  240 DMT = (1.+2.64E9*TIVAR(4)**(-.67)*EXP(EXP1))*DMD 
      GO TO 260 
  250 DMT = (1. + HTCONV(10)*( 6.7 + HTCONV(10)*11.))*DMD 
  260 SPM = DMT/DENS(II) 
      GO TO 900 
!****              EFN 300 STARTS THE REFRASIL OPTION.                  
  300 EXP1 = -TMELT(4)/TW 
      SPM = TMELT(2)*TW**TMELT(3)*EXP (EXP1) 
      GO TO 900 
!****         EFN 400 STARTS THE GRAPHITE AT STAGNATION PT. OPTION.     
!****         GET  PE  FROM THE T-TABLE                                 
  400 CALL TabUpSingle(HORA,TTABL ,CTIME,PE,INDTL(5),2) 
!****         GET  CHLR.                                                
      RNOSE = TMELT(6) 
      IF (TMELT(7))  402,403,402 
  402 CALL TabUpSingle(HORA,TTABL ,CTIME,CHLR,INDTL(5),4) 
      GO TO 405 
  403 CALL CRATIO (TRAJV(18),CHLR) 
  405 RMDOT = TMELT(3)*SQRT(.21* PE/2116.) *EXP(-TMELT(2)/(1.104*TW)) 
  410 DMDOT = TMELT(4)*CHLR*SQRT(PE/2116.*1./RNOSE) 
  415 CSUBCW = .15*   (RMDOT**2/(RMDOT**2+DMDOT**2))**.5 + 2.4E6*1./    &
     &(PE/2116.)**.67 *EXP(-11.1E4 /TW)                                 
  420 WMDOT = .0425*CSUBCW*CHLR*TMELT(5)*SQRT(PE/2116.*1./RNOSE) 
  425 SPM = WMDOT/DENS(II) 
      AUS(1) = CHLR 
      AUS(2) = CSUBCW 
      AUS(3) = 1./(TMELT(3)*SQRT(RNOSE)/100.)**4 
  430 IF(AUS(3).GE.PE/2116.) GO TO 435 
      AUS(3) = -AUS(3) 
  435 GO TO 900 
!****         THAT ENDS THIS OPTION.                                    
!****              EFN 500 STARTS THE FIXED MELTING T                   
!****                   OPTION.                                         
  500 IF(TW-TMELT(2)+ABS(EPTAB(3))) 510,515,514 
  510 SPM = 0. 
      GO TO 900 
  514 T(II) = TMELT(2) 
  515 TD = -3.*T(II)+4.*T(II+1)-T(II+2) 
      TX = TD/2./DELETA(1) 
      CALL CHARP(II,1) 
      CALL TabUpSingle(WARM,PROPT,TMELT(2),CAY,INDTL(4),2) 
      KT = KONVAR(2)+1 
      CALL QNET 
      TQ = -CHARV(1)*CAY*TX*SPACVR(3)/THIKN(1) 
      EPS = -EPTAB(6) 
      TSPM = SPM 
!****         GAMMA AND RHOL ARE IN TIVAR(6) AND (7)                    
!***               AND WERE LOOKED UP IN QNET.                          
      GRHOL = TIVAR(6) * TIVAR(7) 
      HRMHW = TIVAR(2) - HEVAR(1) * T(II) 
      QC = Q(KT,4) 
  520 SMVMG = TIVAR(7) * DENS(II) * SPM + EMDG(1) 
      CALL QBLK(QC,HRMHW,QBKSTR,SMVMG) 
      QVAP = GRHOL * SPM 
      IF(TIVAR(7)) 523, 522, 523 
  522 GRHOL = TIVAR(6) 
  523 SPM = (Q(KT,1) + Q(KT,5) - TQ - QBKSTR)/GRHOL 
      IF(TIVAR(7)) 524,540,524 
  524 CALL TNEST(SPM,TSPM,EPS) 
      KSPM = TSPM 
      EPS = ABS(EPS) 
      GO TO (540,527,530),KSPM 
  527 SPM = AMAX1 (0.0,SPM) 
      GO TO 520 
  530 WRITE(6,1200) EMV,SPM,TQ,QBKSTR,QC,HRMHW    !!! emv never set ???
      GO TO 520 
  540 Q(KT,1) = TQ 
      Q(1,13) = QVAP 
      Q(KT,5) = QBKSTR 
      GO TO 900 
!****    EFN 600 STARTS THE DEPOLYMERIZATION OPTION                     
  600 CALL DEPOLY(SPM) 
      GO TO 900 
!****    EFN 700 STARTS A COMBINATION GRAPHITE SUBLIMATION AND          
!****    FIXED MELTING OPTION                                           
  700 IF(TW-TMELT(3)) 200,200,500 
!****         THE EQUATIONS FOR THE AMOUNTS OF MELT AND CHAR            
!****         ARE DISCUSSED IN SECTION II, PART C, SUB-PART 1,B.        
!****              EFN 800 STARTS THE COMPUTATIONS FOR THE              
!****                   AMOUNT OF CHAR AND MELT.                        
  800 KT = KONVAR(2) 
      SCHAR(KT+1) = SCHAR(KT) + H*HTCONV(4) 
      RHOB = .97*DENSE(1) 
      INDTL(6) = PTSIN(1) + 1. 
      IF(DENS(II)-RHOB) 809, 809, 804 
  804 SSTAR = SCHAR(KT) 
      GO TO 815 
  809 IF(CHART(12)) 810,815,810 
  810 CALL TabUpSingle(DENS(II:),XLOC,RHOB,SSTAR,INDTL(6),1)   ! modified by RLC
  815 HTCONV(8) = SCHAR(KT+1)+CHART(12)*(SSTAR-SCHAR(KT+1)) 
      SMELT(KT+1) = SMELT(KT)+H*HTCONV(5) 
      HTCONV(7) = SMELT(KT+1) 
      THIKN(1) = HTCONV(9)-HTCONV(7) 
      GO TO 998 
!**** EFN 900 SETS UP STORAGE FOR SPM AND SPC                           
  900 HTCONV(5) = SPM 
      HTCONV(4) = SPC 
  998 RETURN 
 1200 FORMAT(/,4X,'MV = ',E15.8,3X,'SPM = ',E15.8,3X,'TQ = ',E15.8,3X, &
        'QBKSTR = ',E15.8//5X,'QC = ',E15.8,3X,'HRMHW = ',E15.8//)               
      STOP 
END Subroutine Ablate   ! ---------------------------------------------
                                          

!C$IBFTC BEE                                                            
      SUBROUTINE BEE (II,NM,N,ITTEMP) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!************COMMON FOR VARIABLE T-P-D INPUT*******FOR BETA-4 TABLE**** 
      DIMENSION TMDENS(1),TMT(1),TMP(1),XFORT(1),XFORP(1),XFORDN(1) 
!     ABOVE IS DUMMY DIMENSION -JUST TO ALLOW SUBSCRIBTING              
      EQUIVALENCE(TMDENS(1),DENS(301)),(TMT(1),T(301)),(TMP(1),P(301)) 
      EQUIVALENCE(XFORT(1),T(151)),(XFORP(1),P(151)),(XFORDN(1),DENS(151&
     &))                                                                
      EQUIVALENCE(NT,T(449)),(NXT,T(450)),(NP,P(449)),(NXP,P(450)) 
      EQUIVALENCE(ND,DENS(449)),(NXD,DENS(450)) 
!     THE T-P-AND DENS INPUTS ARE STORED IN THE                         
!     LAST THIRD OF THEIR RESPECTIVE ARRAYS.                            
!     NT-NP-NDENS ETC ARE STORED IN LOC                                 
!     NOT USED -DUE TO THE IMPOSSIBILITY OF HAVING 150 POINTS           
!     THE X TABLE FOR T(XFORT)ETC ARE STORED IN THE                     
!     MIDDLE THIRD OF THEIR RESPECTIVE ARRAYS                           
!     FOR THE BETA 4 TABLE HAVING PSEUDO - VARIABLE DIMENSION           
      COMMON TXQTAB(200),NOFTS,NOFXS,NTXQ(50),ICON,JCON,END 
      DIMENSION     GVAR(6),XTEMP(150) 
!****         THE BEE ROUTINE ASSUMES THAT THE CHARP SUBROUTINE HAS     
!****         BEEN PREVIOSLY CALLED ,FOR THIS POINT TO BE COMPUTED.     
!                              THIS SUBROUTINE CALCULATES               
!                              BETA(1-4) FOR THE DIFFERENCE             
!                              EQUATIONS.                               
!                              SET UP COUNTERS                          
!****         COMPUTE THE BETA COEFFICIENTS IN CARTESIAN COODINATES.    
      I = II 
      MN = NM 
      J = N 
      TDN = 2. * DELETA(J) 
      NPTS = PTSIN(J) + 1. 
!                       TEST IF THIS IS A PRESSURE OPTION.              
      IF(KONVAR(7)-2) 2,1,2 
    1 CALL PBETA (I,MN,J,ITTEMP) 
      GO TO 400 
!                              GET CORRECT TEMPERATURE,                 
!                              DENSITY, AND PRESSURE.                   
    2 TE = T(I) 
!             SET UP THE TRANSFORMATION TERMS.                          
    8 IF(J-1) 10, 9, 10 
    9 DEN1 = DENS(I) 
      EX = ETASBX(MN) 
      EXX = ETAXX(MN) 
      ET = ETASBT(MN) 
      GO TO 201 
   10 DEN1 = DENSE(J) 
      EX = 1./THIKN(J) 
      EXX = 0.0 
      ET = 0.0 
!                              CALCULATE BETA(1-4)                      
!                       TABLE LOOK-UP FOR CK,DK/DT,CP,                  
  201 MN1 = (INDTL(4)+1)*(J-1)*3+1 
      CALL TABUP(WARM,PROPT(MN1:),TE,PROP,INDTL(4),-3)   ! modified by RLC
      SS = PROP(1) 
      CK = PROP(2) 
      DKDTUP = PROP(3) 
      BETA1 = CHARV(1)*CK/(DEN1*SS*CHARV(3)) 
      BETA2 = BETA1*DKDTUP/CK 
      IF(J-1) 210, 206, 210 
  206 IF(CHART(1)) 214, 210, 214 
  210 BETA3 = 0.0 
      BETA4 = 0.0 
      GO TO 249 
!****              THERE IS CHARRING SO BETA-S 3 AND 4                  
!****                   HAVE VALUES.                                    
  214 IF(MN-(JLAYER(J)+1)) 225, 230, 225 
  225 IF(MN-(JLAYER(J) + NPTS + 1)) 239, 230, 239 
!****         COMPUTE THE DIFFERENCES FOR FRONT OR BACKFACES.           
  230 CALL DFINTC(DENS,I,MN,J,RHON,SD) 
      RHON = RHON/TDN 
      GO TO 247 
  239 RHON = (DENS(I+1) - DENS(I-1)) / TDN 
  247 CALL TABUP(WARM,GASPR,TE,GVAR,INDTL(4),-5) 
      CG = GVAR(2)+GVAR(3) 
      BETA3 = BETA1*(CHARV(2)/CHARV(1)*EX*RHON+EMDG(MN)*CG/CHARV(1)/CK) 
      GEE = GVAR(1)*WDOTP(MN) 
      BETA4 = BETA1*GEE/CHARV(1)/CK 
!****         SEE IF ANOTHER COORDINATE SYSTEM IS DESIRED.              
  249 IF(GENRL(9)) 250, 301, 260 
!****         EFN 250 MODIFIES BETA(3,MN) FOR SPHERICAL COORDINATES.    
  250 IF(GENRL(10)) 254, 252, 254 
!****         Q TRAVELS IN NORMAL DIRECTION.                            
  252 CTRM = 2./ XLOC(MN) 
      GO TO 270 
!****         Q TRAVELS FROM INSIDE TO OUTSIDE.                         
  254 CTRM = -2./ XLOC(MN) 
      GO TO 270 
!****         EFN 260 MODIFIES BETA(3,MN) FOR CYLINDRICAL COORDINATES.  
  260 IF(GENRL(10)) 264, 262, 264 
!****         Q TRAVELS IN NORMAL DIRECTION.                            
  262 CTRM = 1. / XLOC(MN) 
      GO TO 270 
!****         Q TRAVELS FROM INSIDE TO OUTSIDE.                         
  264 CTRM = -1. / XLOC(MN) 
!****         NOW MODIFY BETA(3,MN).                                    
  270 BETA3 = BETA3 - BETA1 * CTRM 
!                              NOW COMPLETE BETA(1-4).                  
  301 ETAXSQ = EX**2 
      BETA(1,MN) = BETA1*ETAXSQ 
      BETA(2,MN) = BETA2*ETAXSQ 
      BETA(3,MN) = BETA3*EX-ET+BETA1*EXX 
      BETA(4,MN) = BETA4 
!                              TEST QUANTITIES TO CHECK                 
!                              IF ITERATION IS NECESSARY                
!                              IN INTPTS                                
      IF (CHART(1)) 375,304,375 
  304 IF (BETA2) 375,305,375 
  305 IF (ABS (WDOTP(MN))) 310,315,310 
  310 IF (ABS(CG))  375,315,375 
  315 NT1 = INDTL(4)+1 
      J1 = (J-1)*NT1*3+NT1 
      IF(ABS(PROPT(J1))) 375,380,375 
  375 ITTEMP = 1 
!                       ITTEMP = 1 MEANS ITERATION IN MIDPTS.           
      GO TO 400 
  380 ITTEMP = 0 
!                       ITTEMP = 0. MEANS NO ITERATION IN MIDPTS.       
!     SEE IF B4 CORRECTION IS NEEDED FOR THIS RUN.                      
  400 IF(NOFXS) 425,600,425 
!     B4 CORRECTION IS NEEDED SO BEGIN COMPUTATION AT 425               
!     IS THIS THE FIRST TIME STEP.                                      
  425 IF (KONVAR(2).NE.1) GO TO 550 
!     DOES THE FIRST TIME FALL BELOW THE TIME TABLE                     
  430 IF (CTIME.GE.TXQTAB(1)) GO TO 440 
  435 CALL ErrorMessage('ARG OUT OF RANGE TOO LOW ON TIME TABLE FOR B4') 
!     DOES THE FIRST XLOC FALL BELOW THE TABLE                          
  440 NUMOFC = 1 
      C = XLOC(1) 
      NTNX = NOFTS + NOFXS 
      NPTOTL = KONVAR(1) 
  443 IF(C - TXQTAB(NTNX))  445,459,450 
  445 IF (TXQTAB(NOFTS+1) -C)  459,459,446 
  446 CALL ErrorMessage('ARG. OUT OF RANGE ON X TABLE')
  450 IF(TXQTAB(NOFTS+1) -C)  446,459,459 
  459 IF(NUMOFC .NE. 1)  GO TO 460 
      C = XLOC(NPTOTL) 
      NUMOFC = 2 
      GO TO 443 
!     INITIALISE TIME TABLE REMEMBERING DEVICE                          
  460 LTIME = 1 
  470 DO 490 JD = LTIME,NOFTS 
!     USE LESS THAN BECAUSE TIME TABLE IS  ASCENDING                    
  480 IF (CTIME.LT.TXQTAB(JD)) GO TO 510 
  490 END DO 
  500 CALL ErrorMessage('ARG OUT OF RANGE  TOO HIGH IN TIME FOR  B4')
!     AT 510 BEGIN TO WEIGHT TVALUES - GENERATE 1WAY TABLE IN XTCOMP.   
  510 JTU = NOFTS+(JD  )*NOFXS 
!     SAVE JD IN LTIME FOR NEXT TIME STEP.                              
      JTL = NOFTS+(JD-1)*NOFXS 
      LTIME = JD 
!     GENERATE XTEMP FROM THE X-ARRAYS STARTING AT JXU AND JXL          
!     FIRST FIND WEIGHTING FACTOR FROM TIME TABLE                       
  525 TFACT = (CTIME - TXQTAB(JD-1))/(TXQTAB(JD) - TXQTAB(JD-1)) 
!     NOW WEIGHT ALL T VALUES.                                          
  530 DO 535 JD = 1,NOFXS 
      JTUALT = JTU+JD 
      JTLALT = JTL+JD 
      XTEMP(JD) = TXQTAB(JTLALT)+(TXQTAB(JTUALT)-TXQTAB(JTLALT))*TFACT 
  535 END DO 
!     SET UP TIME SAVE TO DECIDE WHEN TO COMP. NEW T-RATIO.             
      TIMLST = CTIME 
      GO TO 560 
  550 IF (TIMLST.NE.CTIME) GO TO 470 
      IF(GENRL(9))  555,552,555 
  552 NTLST = NOFTS + JLAST 
      IF(XLOC(MN) .LT. TXQTAB(NTLST)) GO TO 565 
      GO TO 560 
!****         TABLE CAN BE ASCENDING OR DESCENDING IF CART. OR CYL.     
  555 NTLST = NOFTS + JLAST 
      IF (XLOC(MN).GT. TXQTAB(NTLST))  GO TO 565 
  560 JLAST = 1 
!     START X SEARCH FROM WHEREVER YOU LEFT OFF LAST TIME.              
  565 DO 580 JD = JLAST,NOFXS 
      JDNT = NOFTS + JD 
      IF (GENRL(9)) 576,577,576 
  576 IF (GENRL(10)) 577,575,577 
  575 IF (XLOC(MN).GT.TXQTAB(JDNT)) GO TO 585 
      GO TO 580 
  577 IF (XLOC(MN) .LT. TXQTAB(JDNT)) GO TO 585 
  580 END DO 
!     AT 585 GENERATE B4 FROM JD AND JD-1                               
  585 B4ALT = XTEMP(JD-1) + ((XLOC(MN) - TXQTAB(JDNT  -1))/(TXQTAB(JDNT &
     &  ) - TXQTAB(JDNT  -1)))*(XTEMP(JD)-XTEMP(JD-1))                  
      B4ALT = B4ALT/(DEN1*SS*CHARV(3)) 
      JLAST = JD 
!     ADD CORRECTION TERM IN                                            
  590 BETA(4,MN)=BETA(4,MN)+B4ALT 
  600 RETURN 
END Subroutine Bee ! --------------------------------------------------

!$IBFTC CHARP                                                           
      SUBROUTINE CHARP(II,J1) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATIONS FOR THE CHARP ROUTINE ARE DISCUSSED AND     
!****         DEFINED IN SECTION II.B.2( AND IN APPENDIX A              
!                       THIS SUBROUTINE CALCULATES THE                  
!                              CHAR PROPERTIES K(A), K(A) PRIME,        
!                              K BAR,AND CPA.                           
!                       CHAR TABLE ENTRIES.                             
!                       NON-PRESSURE DEPENDENT                          
!****              CHART(1)  = CONTROL                                  
!****                   1 MEANS NORMAL CHARRING OPTION.                 
!****                   2 MEANS TRANSIENT PRESSURE OPTION.              
!****              CHART(2) = DENSITY OF CHAR                           
!****              CHART(3) = ORDER OF REACTION                         
!****              CHART(4) = TIME TO CHANGE FROM LAMINAR TO            
!****                   TURBULENT OR VICE-VERSA.                        
!****                   MINUS MEANS POWERED FLIGHT TRAJECTORY.          
!****                   PLUS MEANS RE-ENTRY TRAJECTORY.                 
!****              CHART(5) = PARAMETER FOR AK                          
!****              CHART(6) = PARAMETER FOR AK                          
!****              CHART(7) = PARAMETER FOR CPA                         
!****              CHART(8) = PARAMETER FOR CPA                         
!****              CHART(9) = CONSTANT FOR LAMINAR FLOW EQUATION        
!****                   EQUATION.                                       
!****              CHART(10) CONSTANT FOR LAMINAR FLOW IN QBLOCK        
!****                   EQUATION.                                       
!****              CHART(11) = CONSTANT FOR TURBULENT FLOW IN QBLOCK    
!****                   EQUATION.                                       
!****              CHART(12) CONSTANT FOR AMOUNT OF CHAR EQUATION.      
      I = II 
      J = J1 
      HARD = DENS(I) 
      DENSC = CHART(2) 
      DDENS = HARD-DENSC 
      IF(CHART(1)) 1,11,1 
    1 IF (J-1) 21,21,11 
!                              PROPERTIES FOR ZERO SMAX AND             
!                              ALL LAYERS GREATER THAN 1.               
   11 AK = 1. 
      CPA = 1. 
      AKP = 0. 
      GO TO 71 
!                              TEST FOR PRESSURE OPTION.                
   21 IF(KONVAR(7)-2) 31, 11, 31 
   31 IF (DDENS) 51,41,51 
   41 AK = 1.-CHART(6) 
      CPA = 1.-CHART(8) 
      AKP = 0. 
      GO TO 71 
!                       COMPUTE (RHO-RHO(C))/(RHO(P)-RHO(C))            
!                       = DELRHO                                        
   51 DELRHO = DDENS/(DENSE(J)-DENSC) 
   54 AK = (1.-CHART(6))+CHART(6)*DELRHO**CHART(5) 
      CPA = (1.-CHART(8))+CHART(8)*DELRHO**CHART(7) 
      AKP = CHART(5)*(AK-1.+CHART(6))/DDENS 
   71 CHARV(1) = AK 
      CHARV(2) = AKP 
      CHARV(3) = CPA 
   91 RETURN 
END Subroutine Charp

!$IBFTC CLCINT                                                          
      FUNCTION CLCINT(IND,DARG,FARG,STOARG) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!                             INTEGRATION SUBROUTINE                    
      DIMENSION STOARG(5) 
    1 DX=DARG 
      FX=FARG 
      IF(DX.NE.0.0) GO TO 41 
   21 STOARG(1)=0.0 
      GO TO 401 
   41 IF(IND.EQ.0) GO TO 101 
      IF(DX.EQ.STOARG(5)) GO TO 201 
!     TRAPEZOIDAL RULE                                                  
  101 STOARG(2)=DX/2.0*(FX+STOARG(3)) 
      GO TO 301 
!     SIMPSON,S RULE                                                    
  201 STOARG(2)=DX/3.0*(FX+4.0*STOARG(3)+STOARG(4))-STOARG(2) 
  301 STOARG(1)=STOARG(1)+STOARG(2) 
  401 CLCINT=STOARG(1) 
      STOARG(4)=STOARG(3) 
      STOARG(3)=FX 
      STOARG(5)=DX 
 5001 RETURN 
END Function Clcint

!$IBFTC COORD                                                           
      SUBROUTINE COORD (NM,JAY) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      RZERO = 0.0 
!****         CHECK TO SEE IF COORDINATE SYSTEM OTHER THAN CARTESIAN    
!****         IS TO BE USED.                                            
      IF(GENRL(9)) 4,270,4 
!****         COMPUTE THE SUM OF THE LAYER WIDTHES.                     
    4 L = KONVAR(4) 
      IF(KONVAR(2) - 1) 5,5,250 
    5 THIKLR = 0.0 
      DO 6 J = 2,L 
    6 THIKLR = THIKLR + THIKN(J) 
!****         EFN 250 STARTS THE MODIFICATION OF THE XLOC(MN) FOR       
!****              SPHERICAL OR CYLINDRICAL COORDINATES.                
  250 MN = NM 
      IF( GENRL(10)) 254, 252, 254 
!****         Q TRAVELS IN THE NORMAL DIRECTION                         
  252 RZERO = ABS(GENRL(9)) + THIKN(1) + THIKLR + HTCONV(7) 
      XLOC(MN) = RZERO - XLOC(MN) 
      GO TO 270 
!****         Q TRAVELS FROM THE INSIDE TO THE OUTSIDE WALL.            
  254 RZERO = ABS(GENRL(9)) 
      XLOC(MN) = RZERO + XLOC(MN) 
  270 TRAJV(50) = RZERO 
      RETURN 
END Subroutine Coord

                                           
!$IBFTC DEE                                                             
      SUBROUTINE DEE(III,MNI) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
!****             CHART TABLE FOR TRANSIENT PRESSURES.                  
!****              CHART(1) = 2                                         
!****              CHART(2) = DENSITY OF CHAR,RHO(C)                    
!****              CHART(3) = ORDER OF REACTION.                        
!****              CHART(4) =RE  TRIP TIME                              
!****              CHART(5) = EMISSIVITY-E(F)                           
!****              CHART(6) = GAMMA                                     
!****              CHART(7) = N FOR K/E COMPUTATION                     
!****              CHART(8) = RHO(C)-TILDA                              
!****              CHART(9) = MA/MC                                     
!****              CHART(10) = PR                                       
!****              CHART(11) = CT                                       
!****              CHART(12) = DELTA                                    
!****              CHART(13) = RHO - TILDA(P)                           
!****              COMPUTE THE EMISSIVITY AND STORE IN                  
!****                TABLD(9).                                          
      II = III 
      MN = MNI 
      JJ = II - KONVAR(1) 
      KT = KONVAR(2) 
      RCFP = CHART(8)*(1.-CHART(5)) 
      RCP = (RCFP+(CHART(6)-1.)*DENS(II))/CHART(6) 
      RPP = DENS(II)-RCP 
      TABLD(9) = 1.-RCP/CHART(8)-RPP/CHART(13) 
!****              COMPUTE THE DERIVATIVE OF EMMISSIVITY                
!****                WRT. ETA AND STORE IN TABLD(10).                   
!****              LOGIC TO COMPUTE 3 POINT DIFFERENCES FOR FRONT       
!****                   AND BACK FACES.                                 
      NPTS = PTSIN(1) + 2. 
      IF(MN-1) 6,8,6 
    6 IF(MN-NPTS) 12, 10, 12 
    8 A = -3. 
      B = 4.0 
      C = -1.0 
      IP1 = II +2 
      IP2 = II 
      IP3 = II +1 
      GO TO 11 
   10 A = 3. 
      B = -4.0 
      C = 1.0 
      EMM = 1.-RCP/CHART(8)-RPP/CHART(13) 
      CALL TabUpSingle(HORA,TTABL,CTIME,TIVAR(4),INDTL(5),4) 
      P(II) = EMM*TIVAR(4)*GRAPH(9)/GRAPH(8)/T(II) 
      IP1 = II-2 
      IP2 = II 
      IP3 = II-1 
   11 JP1 = IP1 - KONVAR(1) 
      JP2 = IP2 - KONVAR(1) 
      JP3 = IP3 - KONVAR(1) 
      TRHON = A * DENS(IP2) + B * DENS(IP3) + C * DENS(IP1) 
      TPORJ = A * DENS(JP2) + B * DENS(JP3) + C * DENS(JP1) 
      TPORI = TRHON 
      TRHONN = DENS(IP1) + DENS(IP2) - 2. * DENS(IP3) 
      RHON = TRHON/2./DELETA(1) 
      RHONN = TRHONN/DELETA(1)**2 
      TPORS = TPORI - TPORJ 
      GO TO 14 
   12 RHON = (DENS(II+1)-DENS(II-1))/2./DELETA(1) 
      RHONN = (DENS(II+1)+DENS(II-1)-2.*DENS(II))/DELETA(1)**2 
      TPORS = DENS(II+1)-DENS(II-1) -DENS(JJ+1) + DENS(JJ-1) 
   14 TEMETA = -((CHART(6)-1.)/CHART(8)+1./CHART(13))/CHART(6) 
      TABLD(10) = TEMETA*RHON 
!****              COMPUTE THE 2ND DERIV. OF EMISSIVITY                 
!****                WRT. ETA AND STORE IN TABLD(11).                   
      TABLD(11) = TEMETA*RHONN 
!****              COMPUTE D(1-4) AND STORE IN TABLD(1-4)               
!****              CALL KFS FOR K(F),K(F)-ETA/K(F),                     
!****                K(F)-ETA-ETA/K(F),K(F)-T/K(F),WHICH ARE            
!****                STORED IN SCDTB(12-15) INC.                        
      CALL KFS(II,MN) 
      EN1 = CHART(7)-1. 
      CAYDE = TABLD(9)**CHART(7)*SCDTB(12) 
      TRM2 = (SCDTB(14)-SCDTB(13)**2)*CAYDE 
   15 IF(CHART(7)) 30, 20, 30 
   20 CAYDEN = CAYDE * SCDTB(13) 
      TRM1 = CAYDEN*SCDTB(13) 
      CAYENN = TRM1 + TRM2 
      GO TO 60 
   30 IF(CHART(7)-1.) 15, 35, 55 
   35 CAYDEN = CAYDE*SCDTB(13)+CHART(7)*TABLD(9)**EN1*SCDTB(12)*        &
     &TABLD(10)                                                         
      TRM1 = CAYDEN*SCDTB(13) 
      TRM4 = CHART(7)*TABLD(9)**EN1*(TABLD(10)*SCDTB(13)*SCDTB(12)+     &
     &SCDTB(12)*TABLD(11))                                              
      CAYENN = TRM1 + TRM2 + TRM4 
      GO TO 60 
   55 CAYDEN = CAYDE*SCDTB(13)+CHART(7)*TABLD(9)**EN1*SCDTB(12)*        &
     &TABLD(10)                                                         
      TRM1 = CAYDEN*SCDTB(13) 
      TRM3 = CHART(7)*EN1*TABLD(9)**(EN1-1.)*SCDTB(12)*TABLD(10)**2 
      TRM4 = CHART(7)*TABLD(9)**EN1*(TABLD(10)*SCDTB(13)*SCDTB(12)+     &
     &SCDTB(12)*TABLD(11))                                              
      CAYENN = TRM1+TRM2+TRM3+TRM4 
   60 TABLD(4) = -WDOTP(MN) 
      TABLD(12) = TEMATA * TPORS /2./ DELETA(1) / DELTO(KT)  !!! temata never set ???
      TABLD(13) = TEMETA * (DENS(II)-DENS(JJ))/DELTO(KT) 
!****              D(1)-T IS STORED IN TABLD(8).                        
      TABLD(18) = CHART(7)*TABLD(9)**EN1*TABLD(13)*SCDTB(12)+           &
     &TABLD(9)**CHART(7)*SCDTB(15)*SCDTB(12)                            
      TABLD(14) = RCFP 
      TABLD(15) = CAYDE 
      TABLD(16) = CAYDEN 
      TRM1NT = TABLD(18)*SCDTB(13) 
      TRM2NT = TABLD(15)*SCDTB(16) 
      TRM3NT =CHART(7)*EN1*TABLD(9)**(EN1-1.)*SCDTB(12)*TABLD(12) 
      TRM4NT = CHART(7)*TABLD(9)**EN1*TABLD(10)*SCDTB(15)*SCDTB(12) 
      TABLD(17) = TRM1NT + TRM2NT + TRM3NT + TRM4NT 
      TABLD(19) = CAYENN 
!****              FIRST CALL PHI FOR PHI(1),PHI(1)-ETA,                
!****                PHI(1)-ETA-ETA,PHI(1)-T,WHICH ARE                  
!****                STORED IN SCDTB(6-10) INC.                         
      CALL PHI(II,MN) 
      TABLD(1) = SCDTB(6)*CAYDE 
      D1N = SCDTB(7)*CAYDE + SCDTB(6)*CAYDEN 
      D1NN = SCDTB(8)*CAYDE +2.*SCDTB(7)*CAYDEN+SCDTB(6)*CAYENN 
      TABLD(2) = 3.*ETASBX(MN)*D1N 
      TABLD(3) = D1NN*ETASBX(MN)**2+D1N*ETAXX(MN) 
      TABLD(5) = D1N 
      TABLD(6) = D1NN 
      TABLD(8) = SCDTB(10)*CAYDE + SCDTB(6)*TABLD(18) 
      TD1NT1 = SCDTB(9)*TABLD(15) 
      TD1NT2 = SCDTB(7)*TABLD(18) 
      TD1NT3 = SCDTB(10)*TABLD(16) 
      TD1NT4 = SCDTB(6)*TABLD(17) 
      D1NT = TD1NT1 + TD1NT2 + TD1NT3 + TD1NT4 
      TABLD(7) = D1NT 
  201 RETURN 
END Subroutine Dee
                                           
!$IBFTC DENSIT                                                          
      SUBROUTINE DENSIT(ISTR,M,K) 
      USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION DBAR(1),CBAR(1) 
      DIMENSION AI(150),BI(150),CI(150),DI(150) 
      EQUIVALENCE  (DBAR,DI) , (CBAR,CI) 
!****       THE EQUATIONS IN THIS SUBROUTINE ARE DISCUSSED IN SECTION V 
!                       I OR II, MEANS POINT WRT. TIME-STEP.            
!                       M OR MN, MEANS POINT WRT. POSITION.             
      II = ISTR 
      KT = KONVAR(2) 
      DT = DELTO(KT) 
      MN = M 
!                       K = ZERO, MEANS USE IMPLICIT SCHEME.            
!                       K = 1, MEANS USE EXPLICIT SCHEME.               
      KAY = K 
      NPTS = PTSIN(1)+2. 
      IT = KONVAR(1) 
      JJ = II - IT 
      KNG = KONVAR(6) 
      IF(CHART(1)) 5,60,5 
    5 IF (IMEX) 100, 7,14 
    7 IF (KAY) 30,30,10 
!                       EXPLICIT SCHEME STARTS AT STATEMENT 10.         
   10 IF (CHART(1)) 11,60,11 
   11 IF(KT-1) 12, 14, 12 
   12 DTR = DT/DELTO(KT-1) 
      DO 13 I = 1,NPTS 
      II = 2*IT + I 
      JJ = II-IT 
!****         DENS(II) IS DEFINED IN EQUATION V.14                      
      DENS(II) = (DTR+1.)*DENS(JJ)-DTR*DENS(I) 
      IF(DENS(II)-CHART(2)*1.000001) 113, 113, 114 
  113 DENS(II) = CHART(2) 
      GO TO 13 
  114 IF(ABS(DENS(II)-DENS(JJ))-ABS(DENS(JJ)*1.E-6)) 115,13,13 
  115 DENS(II) = DENS(JJ) 
   13 END DO 
      GO TO 500 
   14 DO 27 N = 1,NPTS 
      II = KONVAR(1)*KONVAR(6) +N 
      JJ = II-KONVAR(1) 
      MN = N 
      WDPMN = WDOTP(MN) 
      ETATMN = ETASBT(MN) 
      IF (ETATMN) 15,20,20 
   15 DELRHO = DENS(JJ+1)-DENS(JJ) 
      GO TO 21 
   20 DELRHO = DENS(JJ)-DENS(JJ-1) 
!****         DENS(II) IS DEFINED IN EQUATION V.10                      
   21 DENS(II) = DENS(JJ)+DT*(WDPMN-ETATMN*DELRHO/DELETA(1)) 
      IF(DENS(II)-CHART(2)*1.000001) 23, 23, 24 
   23 DENS(II) = CHART(2) 
      GO TO 27 
   24 IF(ABS(DENS(II)-DENS(JJ))-ABS(DENS(JJ)*1.E-6)) 25, 27, 27 
   25 DENS(II) = DENS(JJ) 
   27 END DO 
      GO TO 500 
!                       IMPLICIT SCHEME STARTS AT STATEMENT 30.         
   30 D = DENS(II) 
      EPS = -EPTAB(7) 
      ETATMN = ETASBT(MN) 
   33 CALL WDP(II,MN) 
      WDPMN = WDOTP(MN) 
      IF (ETATMN) 36,35,35 
   35 RHOB = DENS(II-1) 
      GO TO 37 
   36 RHOB = DENS(II+1) 
   37 ALPHA = ABS(DT*ETATMN/DELETA(1)) 
!****         DENS(II) IS DEFINED IN EQUATION V.11                      
      DENS(II) = (DENS(JJ)+DT*WDPMN+ALPHA*RHOB)/(1.+ALPHA) 
      IF(DENS(II)-CHART(2)*1.000001) 40, 40, 45 
   40 DENS(II) = CHART(2) 
      GO TO 50 
   45 IF(ABS(DENS(II)-DENS(JJ))-ABS(DENS(JJ)*1.E-6)) 46, 50, 50 
   46 DENS(II) = DENS(JJ) 
   50 CALL TNEST(DENS(II),D,EPS) 
      EPS = ABS(EPS) 
      COUNT(11) = COUNT(11) + 1. 
      KD = D 
!                       D IS CHANGED IN TNEST( 0.= NO CONV.,            
!                       NON-ZERO = CONVERGENCE).                        
      GO TO (500,33,55),KD 
   55 WRITE(6,550)   ALPHA,RHOB,WDPMN,DENS(II) 
      GO TO 33 
!                       COMPUTATIONS FOR CONSTANT DENSITY.              
   60 NPSI = KONVAR(1) 
      DO 65 L = 1,NPSI 
      II = KONVAR(6)*KONVAR(1)+L 
      JJ = II - KONVAR(1) 
   65 DENS(II) = DENS(JJ) 
      GO TO 500 
!****          AT 100 THE STRICTLY IMPLICIT SCHEME HAS BEEN REQUESTED.  
  100 IF(ETASBT(1) )  105,105,103 
  103 CALL ErrorMessage('A(1) NOT ZERO. CANNOT START SCHEME.') 
  105 RZERO = DELTO(KTIME)/DELETA(1)*ETASBT(MN) 
      DO 150 MN = 1,NPTS 
      JJ = (KNG-1)*IT +MN 
      II = JJ + IT 
      IF (ETASBT(MN))  110,120,120 
  110 ALPHA = 1. 
      GO TO 125 
  120 ALPHA = 0.0 
  125 IF(KT-1) 130,135,130 
  130 DTR = DT/DELTO(KT-1) 
      GO TO 140 
  135 DTR = 0.0 
  140 DENS(II) = (DTR+1.)*DENS(JJ) - DTR*DENS(MN) 
      CALL WDP(II,MN) 
      AI(MN) = RZERO*(ALPHA-1.) 
      BI(MN) = 1.+ ABS(RZERO) 
      CI(MN) = RZERO*ALPHA 
      DI(MN) = DELTO(KTIME)*WDOTP(MN) + DENS(JJ) 
!****          COMPUTE INTERMEDIATE VALUES.                             
      CBAR(MN) = CI(MN)/(BI(MN) -CBAR(MN-1)*AI(MN)) 
      DBAR(MN) = (DI(MN)-DBAR(MN-1)*AI(MN))/(BI(MN)-CBAR(MN-1)*AI(MN)) 
  150 END DO 
!****          NOW COMPUTE THE DENSITIES COMPLETELY IMPLICITLY'''       
!****          DO THE BOUNDARY LAYER FIRST.                             
      NBNLYR = (KNG+1)*IT 
  160 DENS(NBNLYR) = DBAR(NPTS) 
      NPTS1 = NPTS -1 
      DO 200 I = 1,NPTS1 
      MN = NPTS - I 
      II = KNG*IT + MN 
      DENS(II) = DBAR(MN) -CBAR(MN)*DENS(II+1) 
  200 END DO 
  500 RETURN 
  550 FORMAT(/,3X,'ALPHA = ',E18.8,3X,'RHOB = ',E18.8,3X,'WDP = ',E18.8, &
        3X,'RHO-I = ',E18.8)                                               
      STOP 
END Subroutine Densit
                                           
                                           
!$IBFTC DFINTC                                                          
      SUBROUTINE DFINTC (DFFC, III, MNI, JL, FD, SD) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION DFFC(1) 
      MN = MNI 
      J = JL 
      II = III 
      NPT = PTSIN(J) + 2. 
      NPTS=JLAYER(J) +NPT 
      IF(MN-1) 5,10,5 
    5 IF(MN-NPTS) 20, 15, 20 
   10 FDI = -3.*DFFC(II) + 4.*DFFC(II+1) - DFFC(II+2) 
      SDI = DFFC(II+2) + DFFC(II) -2.*DFFC(II+1) 
      GO TO 30 
   15 FDI = 3.*DFFC(II) - 4.* DFFC(II-1) + DFFC(II-2) 
      SDI = DFFC(II-2) + DFFC(II) -2.*DFFC(II-1) 
      GO TO 30 
   20 FDI = DFFC(II+1) - DFFC(II-1) 
      SDI = DFFC(II+1) + DFFC(II-1) - 2.*DFFC(II) 
   30 FD = FDI 
      SD = SDI 
  100 RETURN 
END Subroutine Dfintc
                                          
!$IBFTC DPRESS                                                          
      SUBROUTINE DPRESS                                                 &
     &             (III,                                                &
     &             MN1,                                                 &
     &             IMP)                                                 
!****              WHERE III MEANS THE SUBSCRIPT OF THE POINT           
!****                   TO BE COMPUTED.                                 
!****              SUBSCRIPT RELATED TO III, BUT FOR ONLY               
!****                   ONE TIME STEP.                                  
!****              ONE MEANS EXPLICIT POINTS ARE TO BE                  
!****                   CALCULATED.                                     
!****              ZERO MEANS IMPLICIT POINTS ARE TO BE                 
!****                   CALCULATED.                                     
!****              - 1 MEANS, THE ROUTINE WAS CALLED FROM STARTT.       
!****                         INPUT     COMMON                          
   
USE TimeDt                                                                     
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
      IM = IMP 
      II = III 
      MN = MN1 
      JJ = II - KONVAR(1) 
      CTIME = CTIME 
      NPTS = PTSIN(1)+2. 
      KT = KONVAR(2) 
      IF(IM) 200, 300, 10 
   10 IF (KT-1) 60,15,60 
!****              EFN 10 STARTS THE EXPLICIT SCHEME                    
   15 NPTS = NPTS - 2 
      TIVAR(9) = 0.0 
      TIVAR(10) = 0.0 
      NPTS1 = NPTS +1 
      KN = 1 
      KM = 0 
      ITIME = KONVAR(1) 
      P(ITIME+1) = P(1) 
      DO 50 IR = 2,NPTS1 
      II = KONVAR(6) *KONVAR(1) + IR 
      JJ = II-KONVAR(1) 
      MN = IR 
!****              COMPUTE THE D VARIABLES FOR EACH POINT.              
      CALL DEE(JJ,MN) 
      DTN = DELTO(KT)/DELETA(1) 
      RZ = DTN/2. 
      B1 = TABLD(1)*ETASBX(MN)**2 
      R1 = DTN/DELETA(1) *B1 
      R2 = P(JJ+1)+P(JJ-1) 
      R3 = P(JJ+1)-P(JJ-1) 
      TRM = R1*R2 
      TRM1 = RZ*R3 
      RN = R3/2./DELETA(1) 
      B2 = TABLD(2)*ETASBX(MN)+TABLD(1)*ETAXX(MN) 
      B3 = B1*RN 
      B4 = B3-ETASBT(MN) 
      TRM3 = (TABLD(3)*DELTO(KT) - 2.*R1)*P(JJ) 
      P(II) = P(JJ)*(1.+TRM+B2*TRM1+TRM3)+B4*TRM1+DELTO(KT)*TABLD(4) 
   50 END DO 
      KEND = NPTS +2 +KONVAR(1) 
   51 P(KEND) = P(NPTS+2) 
      GO TO 400 
!****              EFN 60 STARTS EXPLICIT SCHEME AFTER 1ST TIME STEP.   
   60 DTR = DELTO(KT)/DELTO(KT-1) 
      DO 100 IR = 1,NPTS 
      II = 2*KONVAR(1) + IR 
      JJ = II - KONVAR(1) 
      P(II) = (DTR+1.)*P(JJ)-DTR*P(IR) 
  100 END DO 
      GO TO 400 
!****              COMPUTE THE FRONT FACE DENSITY, AND THEN             
!****                   SET ALL THE DENSITIES EQUAL TO IT.              
  200 CALL RHOFFC(1,1,-1) 
      DO  250 IR = 1,NPTS 
      P(IR) = P(1) 
  250 END DO 
      GO TO 400 
!****              EFN 300 STARTS THE IMPLICIT SCHEME.                  
  300 CALL DEE(II,MN) 
      IF(MN-1) 310, 312, 310 
  310 IF(MN-NPTS) 315, 313, 315 
  312 NI = 0 
      GO TO 314 
  313 NI = 1 
  314 CALL RHOFFC(II,MN,NI) 
      GO TO 400 
  315 B1 = TABLD(1)*ETASBX(MN)**2 
      B2 = TABLD(2)*ETASBX(MN)+TABLD(1)*ETAXX(MN) 
      B3 = B1*(P(II+1)-P(II-1))/2./DELETA(1) 
      B4 = B3 - ETASBT(MN) 
      Q(KT,8) = AMAX1(Q(KT,8),B2) 
      Q(KT,9) = AMAX1(Q(KT,9),B4) 
      DTN = DELTO(KT)/DELETA(1) 
      RZ = DTN/DELETA(1) 
      R1 = DTN 
      DNB2 = DELETA(1) * B2 
      IF(B2) 318, 318, 320 
  318 A1 = 0.0 
      GO TO 322 
  320 A1 = 1.0 
  322 IF(B4) 324, 324, 326 
  324 A2 = 0.0 
      GO TO 328 
  326 A2 = 1.0 
  328 G1 = RZ*(B1 + A1*DNB2) 
      G2 = RZ*(B1 + (A1 - 1.0)*DNB2) 
      G3 = R1 * B4 * A2 
      G4 = R1 * B4*  ( A2 - 1.0) 
      A = G1 + G2 -TABLD(3) * DELTO(KT) 
      B = 1.0 + G3 + G4 - G1 * P(II+1) - G2 * P(II-1) 
      C = G3*P(II+1) + G4*P(II-1) + P(JJ) + TABLD(4)*DELTO(KT) 
      D = 4.*A*C/(B**2) 
      IF(A) 390,340,340 
  340 IF(C) 390, 342, 342 
  342 IF(B) 350, 360, 352 
  350 A3 = 1.0 
      GO TO 360 
  352 A3 = -1.0 
  356 IF(.02-D) 360, 358, 358 
  358 IF(1.E-8-ABS(D)) 1358, 359, 359 
 1358 TRM = -D/2.0*(1.0 -D/4.0 +D**2/8.-.078125*D**3) 
      GO TO 362 
  359 TRM = -D/2. 
      GO TO 362 
  360 TRM = 1.0 + A3*SQRT(1.+D) 
  362 P(II) = -B/2./A*TRM 
      GO TO 400 
  390 CALL ErrorMessage('A OR C IS MINUS')
  400 RETURN 
END Subroutine Dpress
                                          
!$IBFTC EMG                                                             
      SUBROUTINE EMG(II) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATION FOR THE MASS FLUX OF THE GAS, M-DOT-G,       
!****         IS DEFINED IN SECTION II.B.2                              
      K2 = II 
      IF(KONVAR(7)-2) 1,121,1 
    1 K1 = KONVAR(6)*KONVAR(1) 
      NPTS = PTSIN(1)+1. 
      EMDG(NPTS+1) = 0. 
      IF (K2) 3,11,3 
    3 NPTS = 1 
      K = K2-K1 
      GO TO 21 
   11 K = 1 
   21 DO 107 I1 = 1,NPTS 
      N = (NPTS-I1+1)*K 
      I = N+K1 
      IF(N) 45, 45, 55 
   45 N = I 
   55 CALL WDP(I,N) 
!****              USING TRAPEZOIDAL RULE, SOLVE FOR                    
!****                   M-DOT-G OF EACH POINT.                          
      IF(GENRL(9)) 65, 60, 75 
!****         EFN 60 COMPUTES EMDG(N) FOR CARTESIAN COORDINATES.        
   60 EMDG(N) = EMDG(N+1)-.5*DELETA(1)*(WDOTP(N)/ETASBX(N)+WDOTP(N+1)/ET&
     &ASBX(N+1))                                                        
      GO TO 107 
!****         EFN 65 COMPUTES EMDG(N) FOR SPHERICAL COORDINATES.        
   65 CALL ErrorMessage('OPTION FOR EMDG SPERICAL COORDINATES NOT DETERMINED')
   75 DRS1 = XLOC(N) 
      DRS2 = XLOC(N+1) 
   85 EMDG(N) = EMDG(N+1) - .5*DELETA(1)/TRAJV(50)    *(WDOTP(N)*DRS1/  &
     & ETASBX(N) + WDOTP(N+1) * DRS2/ETASBX(N+1))                       
  107 END DO 
      GO TO 200 
  121 CALL PEMDG(K2) 
  200 RETURN 
END Subroutine Emg
                                           

!$IBFTC EXEC2                                                           
      SUBROUTINE EXEC2 
      USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!************COMMON FOR VARIABLE T-P-D INPUT*******FOR BETA-4 TABLE**** 
      DIMENSION TMDENS(1),TMT(1),TMP(1),XFORT(1),XFORP(1),XFORDN(1) 
!     ABOVE IS DUMMY DIMENSION -JUST TO ALLOW SUBSCRIBTING              
      EQUIVALENCE(TMDENS(1),DENS(301)),(TMT(1),T(301)),(TMP(1),P(301)) 
      EQUIVALENCE(XFORT(1),T(151)),(XFORP(1),P(151)),(XFORDN(1),DENS(151&
     &))                                                                
      EQUIVALENCE(NT,T(449)),(NXT,T(450)),(NP,P(449)),(NXP,P(450)) 
      EQUIVALENCE(ND,DENS(449)),(NXD,DENS(450)) 
!     THE T-P-AND DENS INPUTS ARE STORED IN THE                         
!     LAST THIRD OF THEIR RESPECTIVE ARRAYS.                            
                                                                        
!     NT-NP-NDENS ETC ARE STORED IN LOC                                 
!     NOT USED -DUE TO THE IMPOSSIBILITY OF HAVING 150 POINTS           
!     THE X TABLE FOR T(XFORT)ETC ARE STORED IN THE                     
!     MIDDLE THIRD OF THEIR RESPECTIVE ARRAYS                           
!     FOR THE BETA 4 TABLE HAVING PSEUDO - VARIABLE DIMENSION           
      COMMON TXQTAB(200),NOFTS,NOFXS,NTXQ(50),ICON,JCON,END 
      COMMON MK 
      DIMENSION TABLN(120) 
   18 GO TO (19,69),JCON 
   19 CALL DELTAT(0) 
      ZERO = 0.0 
!****              MAKE SURE TNEST IS CLEARED.                          
      CALL TNEST(ZERO,ZERO,ZERO) 
      KT = KONVAR(2) 
      CALL GAMA(-DELTO(KT)) 
!                       CHECK NZERO TO SEE                              
!                       WHICH POINT TO START THE                        
!                       EXPLICIT EQUATION.                              
      NZERO1 = 5 - KONVAR(3) 
!                       CHECK IF PAST THE STOP TIME,                    
!                       IF NOT, CALCULATE THE TIME FOR                  
!                       THE TIME STEP.                                  
      IF(IMEX.EQ.-1.AND.CHART(1).NE.0.0)  CALL TMPIMP(-1) 
      IF(CTIME+DELTO(KT) - GENRL(2)) 22, 75, 75 
   22 CTIME = CTIME + DELTO(KT) 
!                       DELTA T IS STORED IN DELTO(KTIME)               
      NPTS = PTSIN(1) 
      IF(NPTS-1) 24, 23, 24 
!                  THIN SKIN COMPUTATIONS.                              
   23 CALL PTX12 
      GO TO 69 
!                  CALCULATE THE DENSITIES EXPLICITLY.                  
   24 IF(IMEX.EQ.-1)  GO TO 30 
      CALL DENSIT(1,1,1) 
      IF(KONVAR(7)-2) 26, 25, 26 
!                  PRESSURE COMPUTATIONS.                               
   25 CALL PRESSM(0,0,0,1) 
!                       CALCULATE THE TEMPERATURES                      
!                       OF THE INTERIOR POINTS BY                       
!                       USING THE EXPLICIT SCHEME.                      
   26 CALL INTPTS(0,NZERO1) 
!                       SET UP SPACING FOR TIME STEP                    
!                       BEING COMPUTED.                                 
!****              COMPUTE ONLY R VALUES                                
   30 CALL RITER(2) 
!                       COMPUTE THE AMOUNT MELTED                       
!                       AND CHARRED.                                    
   35 CALL ABLATE (DELTO(KT)) 
!****              COMPUTE ETA DERIVATIVES AND X LOCATIONS.             
      CALL RITER(3) 
!****         M-DOT-G NEEDS TO BE COMPUTED IF THERE IS CHARRING         
      IF(CHART(1)) 40, 50, 40 
   40 CONTINUE 
      CALL EMG(0) 
!                       COMPUTE THE RATE OF MELTING                     
!                       AND CHARRING.                                   
   50 CALL ABLATE(0.0) 
!                       SET UP CONSTANTS FOR ITERATION                  
      EP = -ABS(EPTAB(6)) 
      RMELT = HTCONV(5) 
      RMELT1 = RMELT 
!                       COMPUTE DELTA T FOR ITERATION                   
!                       DEPENDING ON KD.                                
      KD = -1 
      CALL DELTAT(KD) 
      IF(KD) 54, 53, 53 
   53 CALL GAMA(-DELTO(KT)) 
      GO TO 24 
!****                   COMPUTE THE ETA(T)-S.                           
   54 CALL NSUBT 
      IF (IMEX) 249,249,276 
  249 IF (ABS(RMELT1*ETASBX(1))-.001) 275,275,250 
  250 IF(EPTAB(6))  254,275,275 
!                  SECTION FOR ITERATION ON EXPLICIT POINTS ONLY,       
!                       DUE TO RATE OF MELT.                            
  254 IF(KONVAR(3) -3) 260,256,260 
  256 N4 = PTSIN(1) 
      CALL INTPTS(1,N4) 
  260 CALL PTX1 
      CALL ABLATE(0.0) 
      CALL TNEST(HTCONV(5),RMELT,EP(1))    ! modified by RLC
      EP = ABS(EP) 
      KRM = RMELT 
      COUNT(3) = COUNT(3) + 1. 
      GO TO (275,54,265),KRM 
  265 II = KONVAR(6)*KONVAR(1) + 1 
      WRITE(6,500) HTCONV(5),T(II),T(II+1),T(II+2) 
      GO TO 54 
!                       NOW CALL INTPTS TO COMPUTE THE                  
!                       VALUES IMPLICITLY.                              
  275 CALL INTPTS(-1,KONVAR(3)) 
      IF(IMEX) 156,276,276 
!                       COMPUTE THE VALUES FOR THE INTER-FACE.          
  276 CALL INTFC 
!                       CALCULATE THE VALUES FOR                        
!                       THE AIR-GAP.                                    
      IF(KONVAR(5)) 55, 56, 55 
   55 CALL AIRGAP 
!                       CALCULATE THE VALUES AT THE                     
!                       FRONT-FACE.                                     
   56 CALL PTX1 
!                       COMPUTE THE RATE OF MELTING                     
!                       AND CHARRING.                                   
  156 CALL ABLATE(0.0) 
      IF(IMEX) 57,57,60 
!                       TEST RATE OF MELT.                              
   57 IF(ABS(RMELT1*ETASBX(1))-.001) 60, 60, 58 
   58 IF(EPTAB(6)) 60, 60, 59 
!                       USING TNEST, ITERATE ON THE                     
!                       RATE OF MELTING.                                
   59 CALL TNEST(HTCONV(5),RMELT,EP(1))    ! modified by RLC
      EP = ABS(EP) 
      COUNT(3) = COUNT(3) +1. 
      KRM = RMELT 
      GO TO (60, 54, 159),KRM 
  159 II = KONVAR(6)*KONVAR(1) +1 
      WRITE (6,500) HTCONV(5),T(II),T(II+1),T(II+2) 
      GO TO 54 
!                       COMPUTE THE HEAT BALANCE.                       
   60 CALL HTBLNS(1) 
!                       PRINT OUTPUT                                    
      CALL OUTPUT(0) 
!****         CHECK TO SEE IF IT IS TIME FOR THE TAPE INTERRUPT.        
      IF( GENRL(8) .EQ. 0.0 .OR. CTIME .LE. GENRL(8)) GO TO 69 
      CALL MDUMP (MK) 
      GENRL(8) = 0.0 
!****         WHEN MK =                                                 
!****              0, MEANS CHECKPOINT WAS SUCCESSFUL.                  
!****              +, MEANS THE RESTART WAS SUCCESSFUL.                 
!****              -, MEANS THE RESTART TAPE WAS NOT GENERATED.         
      IF(MK) 69,69,68 
   68 ICON = 2 
      GO TO 110 
!                       CHECK TO MAKE SURE THAT                         
!                       ENTIRE LENGTH IS NOT CHARRED                    
!                       OR MELTED AWAY.                                 
   69 IF(SMELT(KT+1)-HTCONV(9)) 70,71,71 
!                       INITIALIZE FOR NEXT TIME STEP.                  
   70 CALL INIT(2) 
      GO TO 18 
   71 WRITE (6,501) 
!                       RUN IS COMPLETED, CALL INIT                     
!                       TO START NEXT CASE                              
   75 CALL INIT(3) 
      GO TO 90 
      CALL INIT(3) 
!****          IF GENRL 4  = 0  ALLNEW MUST = 1 TO RUN MULT. CASES.     
   90 IF(GENRL(4)) 97,95,97 
   95 IF(ALLNEW) 97,100,97    !!! allnew is never set ???
!****          IF END = 1 NEXT CASE USES ALL NEW INPUT.                 
   97 IF(END) 102,103,102 
  102 ICON = 1 
      GO TO 110 
  103 ICON = 2 
      GO TO 110 
  100 WRITE (6,105) 
  105 FORMAT(/' CANNOT USE VARIABLE TEMP WITHOUT ALLNEW')
  110 RETURN 
  500 FORMAT(/,3X,'SPM = ',E18.8,3X,'T(1-3) = ',3E20.8) 
  501 FORMAT(/,4X,'ENTIRE LAYER MELTED AWAY') 
END Subroutine Exec2
                                          
!$IBFTC HTBLNS                                                          
      SUBROUTINE HTBLNS(IKBLN) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATIONS FOR THIS ROUTINE ARE DEFINED AND DISCUSSED  
!****         IN SECTION IX, ACCURACY OF THE PROGRAM.                   
!                       THIS SUBROUTINE CALCULATES THE HEAT BALANCE     
!                       FOR THE PROGRAM.  IT IS STORED IN QI.           
      DIMENSION QI(2,10),QIDOT(2,10),TMP(5) 
      LAYER = KONVAR(4) 
!****         SEE IF AIRGAP IS LAST LAYER.                              
      IF (IKBLN) 1,5,1 
    1 IF(KONVAR(5) -LAYER) 2,4,2 
    2 MN = KONVAR(1) 
      II = KONVAR(1)*(KONVAR(6)+1) 
      KTIME = KONVAR(2) 
      TERM = 3.*T(II)-4.*T(II-1)+T(II-2) 
      MN1 = (INDTL(4)+1)*(3*(KONVAR(4)-1)+1)+1 
      CALL TabUpSingle(WARM,PROPT(MN1:),T(II),CAYUP,INDTL(4),1)   ! modified by RLC
      Q(KTIME+1,6) = -CAYUP*ETASBX(MN)*TERM/2./DELETA(LAYER) 
      GO TO 6 
!****         IF LAST IS AIRGAP     DO TABLE LOOK UP ON BACKFC          
    4 CALL TabUpSingle (HORA,TTABL,CTIME,TIVAR(5),INDTL(5),5) 
      Q(KTIME+1,6) = TIVAR(5) 
    5 IF(OUTTB(19)) 600,6,600 
!                       IF KBLNS = 0.  THE ROUTINE                      
!                       WAS CALLED FROM STARTT.                         
!                       IF KBLNS = 1.  THE ROUTINE WAS                  
!                       CALLED FROM EXEC.                               
    6 QITOT = 0.0 
      N0 = 5 - KONVAR(3) 
      KTM = KONVAR(2) 
      KT1 = KTM + 1 
      IF(IKBLN) 50, 8, 50 
    8 DT = ABS(GENRL(3)) 
      LAYER = KONVAR(4) 
      IF(GENRL(3)) 10, 10, 15 
!****              CONSTANT HEAT FLUX                                   
   10 DO 11 K = 1,6 
   11 Q(KT1,K+6) = DT*Q(KT1,K) 
      QITOT = Q(KT1,7) 
      GO TO 20 
!****              LINEAR HEAT FLUX                                     
   15 DO 16 K = 1,6 
   16 Q(KT1,K+6) = DT*Q(KT1,K)/2. 
      QITOT = Q(KT1,7) 
   20 Q(4,13) = QITOT 
      QIDOT(2,1) = Q(KT1,1) 
      QI(2,1) = QITOT 
      DO 25 I = 2,10 
      QIDOT(2,I) = 0.0 
   25 QI(2,I) = 0.0 
      GO TO 550 
!****              TRAPEZOIDAL RULE FOR Q-S.                            
   50 DO 55 J = 7,12 
   55 Q(KT1,J) = Q(KTM,J) + .5*DELTO(KTM)*(Q(KT1,J-6)+Q(KTM,J-6)) 
      DO 500 J = 1,LAYER 
!****              SKIP FOR THE AIRGAP                                  
      IF(J-KONVAR(5)) 56,70,56 
   56 CONTINUE 
      MNP = 1 + (J-1)*3*(INDTL(4)+1) 
      NPTS = PTSIN(J)+1. 
      TDN = 2. * DELETA(J) 
      IF (IMEX .EQ. -1) GO TO 61 
!****              CALL BEE FOR EXPLICIT POINTS ONLY.                   
      DO 60 N = N0,NPTS,2 
      II = KONVAR(6) * KONVAR(1) + JLAYER(J) + N 
      MN = JLAYER(J) + N 
      IF(KONVAR(7) - 2) 58, 57, 58 
   57 CALL DEE(II, MN) 
   58 CALL CHARP(II,J) 
      CALL BEE(II,MN,J,ITTEMP) 
   60 END DO 
   61 NFFC = KONVAR(6) * KONVAR(1) + JLAYER(J) + 1 
      NBFC = KONVAR(6) * KONVAR(1) + JLAYER(J) + NPTS +1 
      MNFFC = JLAYER(J) + 1 
      MNBFC = JLAYER(J)  + NPTS +1 
      IF(KONVAR(7) - 2) 63, 62, 63 
   62 CALL DEE(NFFC, MNFFC) 
   63 CALL CHARP(NFFC,J) 
      CALL BEE(NFFC,MNFFC,J,ITTEMP) 
      IF(KONVAR(7) - 2) 65, 64, 65 
   64 CALL DEE ( NBFC, MNBFC) 
   65 CALL CHARP(NBFC,J) 
!****         IF OTHER THAN CARTESIAN COORDINATES ARE USED WE  ASSUME   
!****              THAT THE LAST POINT HAS THE SAME BETA-S AS THE       
!****              PREVIOUS POINT.                                      
   66 IF(GENRL(9)) 67, 69, 67 
   67 DO 68 I = 1,4 
      BETA(I,MNBFC) = BETA(I,MNBFC-1) 
   68 END DO 
      GO TO 70 
   69 CALL BEE( NBFC,MNBFC,J,ITTEMP) 
   70 NPTSS = NPTS+1 
      DO 400 I = 1,NPTSS 
      II = KONVAR(6)*KONVAR(1)+I + JLAYER(J) 
      JJ = II-KONVAR(1) 
      MN = JLAYER(J)+I 
      CALL CHARP(II,J) 
      IF(KONVAR(7) - 2) 75, 80, 75 
   75 CALL TABUP(WARM,PROPT(MNP:),T(II),PROP,INDTL(4),-2)    ! modified by RLC
      CAY = PROP(2)*CHARV(1) 
      GO TO 90 
   80 CALL KEFFS(II,MN) 
      CAY = HTCONV(15) 
   90 IF(J-1) 101, 103, 101 
  101 EXX = 0.0 
      EX = 1./THIKN(J) 
      GO TO 109 
  103 EXX = ETAXX(MN) 
      EX = ETASBX(MN) 
  109 CALL DFINTC(DENS,II,MN,J,RHON,SD) 
      CALL DFINTC(T,II,MN,J,TN,SD) 
      RHON = RHON/TDN 
      TN = TN/TDN 
  130 TT = (T(II)-T(JJ))/DELTO(KTM) 
      S4BAR = BETA(1,MN)*RHON*CHARV(2)/CHARV(1) 
      TRM1 = CAY*EX/BETA(1,MN) 
      TRM2 = BETA(1,MN)*EXX/EX**2-BETA(3,MN) + S4BAR 
      IF (J-KONVAR(5)) 150,140,150 
  140 TTI = (T(II+1)-T(JJ+1))/DELTO(KTM) 
!****         QINT IS DEFINED IN EQUATION IX.6.                         
      QINT = THIKN(J)*DENSE(J)*PROP(1)*(TT+TTI)/2. 
      TT = TTI 
      JI = II + 1 
      GO TO 410 
  150 IF(I-1) 160,155,160 
  155 DX = 0.0 
      GO TO 170 
  160 IF(I-NPTS) 165,168,165 
  165 DX = DELETA(J) 
      GO TO 170 
  168 DX = ABS(XLOC(MN) - XLOC(MN-1)) 
      DX = DX * ETASBX(MN) 
!****         QINT IS DEFINED IN EQUATION IX.4                          
  170 QINT = CLCINT(0,DX,TRM1*(TT+TRM2*TN-BETA(4,MN)),TMP) 
  400 END DO 
      JI = II 
  410 IF(J-LAYER) 411,412,411 
  411 CALL TabUpSingle(WARM,CTABL,T(JI),CI,INDTL(4),J) 
      GO TO 416 
  412 IF(J-KONVAR(5)) 415,411,415 
  415 CI = 0.0 
  416 QIDOT(2,J) = QINT + CI*TT 
!****         TRAPEZOIDAL RULE FOR QI TOTAL.                            
      QI(2,J) = QI(1,J)+DELTO(KTM)*(QIDOT(2,J)+QIDOT(1,J))/2. 
      QITOT = QITOT+QI(2,J) 
  500 END DO 
      QITOT = QITOT + Q(KTM,12) 
      Q(4,13) = QITOT 
!****              SET UP VALUES FOR NEXT TIME STEP.                    
  550 DO 560 J = 1,10 
      QIDOT(1,J) = QIDOT(2,J) 
      QI(1,J) = QI(2,J) 
  560 END DO 
  600 RETURN 
END Subroutine Htblns
                                          
!$IBFTC INIT                                                            
      SUBROUTINE INIT(LL) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!************COMMON FOR VARIABLE T-P-D INPUT*******FOR BETA-4 TABLE**** 
      DIMENSION TMDENS(1),TMT(1),TMP(1),XFORT(1),XFORP(1),XFORDN(1) 
!     ABOVE IS DUMMY DIMENSION -JUST TO ALLOW SUBSCRIBTING              
      EQUIVALENCE(TMDENS(1),DENS(301)),(TMT(1),T(301)),(TMP(1),P(301)) 
      EQUIVALENCE(XFORT(1),T(151)),(XFORP(1),P(151)),(XFORDN(1),DENS(151&
     &))                                                                
      EQUIVALENCE(NT,T(449)),(NXT,T(450)),(NP,P(449)),(NXP,P(450)) 
      EQUIVALENCE(ND,DENS(449)),(NXD,DENS(450)) 
!     THE T-P-AND DENS INPUTS ARE STORED IN THE                         
!     LAST THIRD OF THEIR RESPECTIVE ARRAYS.                            
!     NT-NP-NDENS ETC ARE STORED IN LOC                                 
!     NOT USED -DUE TO THE IMPOSSIBILITY OF HAVING 150 POINTS           
!     THE X TABLE FOR T(XFORT)ETC ARE STORED IN THE                     
!     MIDDLE THIRD OF THEIR RESPECTIVE ARRAYS                           
!     FOR THE BETA 4 TABLE HAVING PSEUDO - VARIABLE DIMENSION           
      COMMON TXQTAB(200),NOFTS,NOFXS,NTXQ(50),ICON,JCON,END 
!****              THIS SUBROUTINE INITIALIZES THE INPUT.               
!****              IT SETS UP CONSTANTS, COMPUTES THE X-LOCATIONS,      
!****              ETA-S,DELETA-S,ETC. LL = 1                           
!****              IT ALSO INITIALIZES FOR THE NEXT TIME STEP           
!***               LL = 2, AND FOR THE NEXT CASE LL = 3                 
!****              LL = 4 REINITIALIZES AFTER AN MDUMP RESTART.         
      L = LL 
      GO TO ( 10, 150, 190, 191, 15 ),L 
!****              PRIMARY INITIALIZATION                               
   10 KONVAR(5) = AGPTB(1)*1.01 
      NOFTS = NTXQ(1) 
      NOFXS = NTXQ(2) 
      TIME = GENRL(1) 
      CTIME = TIME + ABS(GENRL(3)) 
      TMELT(8) = TMELT(8)*HTCONV(3) 
      TMELT(9) = TMELT(9)*HTCONV(3) 
      KONVAR(2) = 1 
      KONVAR(3) = 3 
      LAYER = ABS(GENRL(6)) * 1.01 
      KONVAR(4) = LAYER 
      KONVAR(6) = 1 
      KONVAR(7) = CHART(1)*1.01 
      KONVAR(8) = TMELT(1)*1.01 
      COUNT(1) = 1. 
!****              POINT TO RESTART FROM MDUMP.                         
   15 HTCONV(9) = THIKN(1) 
      ZERO = 0.0 
!****              CALL ENLARG TO FILL UP THE TABLES                    
!****                        IF NECESSARY.                              
      CALL ENLARG 
!****              COMPUTE THE DK/DT TABLE BY CALLING DIFTAB.           
!****              COMPUTE SUBSCRIPTS TO USE FOR FIRST LOCATION         
!****              OF THE CK TABLE AND DK/DT TABLE.                     
   30 NT1 = INDTL(4)+1 
      DO 35 MN1 = 1,LAYER 
      MN11 = MN1-1 
      MN2 = NT1*(3*MN11+1)+1 
      MN3 = NT1*(3*MN11+2)+1 
      IF(PROPT(MN3))  35, 34, 35 
   34 CALL DIFTAB(PROPT(MN2),WARM,PROPT(MN3),INDTL(4),1) 
   35 END DO 
!****              TEST TO SEE IF CHARRING OPTION IS USED.              
      IF (KONVAR(7)) 50,50,40 
   40 I = 3*NT1+1 
      J = 4*NT1+1 
!****              COMPUTE THE DIFFERENTIAL OF M-BAR                    
      CALL DIFTAB(GASPR(I),WARM,GASPR(J),INDTL(4),1) 
   50 IF (KONVAR(8)-1) 60,55,60 
   55 L = 6*(INDTL(5)+1)+1 
      K = L+INDTL(5)+1 
!****              COMPUTE THE DERIVATIVE OF S-CHAR-MAX.                
      CALL DIFTAB(TTABL(L),HORA,TTABL(K),INDTL(5),1) 
   60 IF(KONVAR(7)) 61, 64, 61 
!****              THERE IS CHAR, SO TEST FOR A DENSITY                 
!****                        OF CHAR                                    
   61 IF (CHART(2)) 62,62,64 
   62 CALL ErrorMessage('THE DENSITY IS ZERO OR NEGATIVE')
   64 IF (KONVAR(5)-1) 68,65,68 
   65 CALL ErrorMessage ('THE AIR GAP CANT BE FIRST')
   68 DO 80 J = 1,LAYER 
      IF (THIKN(J)) 70,78,70 
   70 IF (DENSE(J)) 72,78,72 
   72 IF (PTSIN(J)) 80,74,80 
   74 IF (J-KONVAR(5)) 78,80,78 
   78 CALL ErrorMessage('EITHER DENSE, PTSIN, OR THIKN IS ZERO')
   80 END DO 
!****              PTSIN MUST BE ODD                                    
      JLAYER(1) = 0 
      DO 87 J = 1,LAYER 
      EX = PTSIN(J)/2. 
      KX = IFIX(EX) 
      FX = FLOAT(KX) 
      IF( J .EQ. KONVAR(5)) GO TO 84 
      IF(EX .NE. FX) GO TO 85 
   83 PTSIN(J) = PTSIN(J)+1. 
      GO TO 85 
   84 PTSIN(J) = 0. 
   85 DELETA(J) = 1./(PTSIN(J)+1.) 
      NPTS = PTSIN(J-1) 
      IF (J-1) 86,87,86 
   86 JLAYER(J) = JLAYER(J-1)+NPTS+1 
   87 END DO 
      LAYER1 = LAYER + 1 
      NPTS = PTSIN(LAYER) 
      JLAYER(LAYER1) = JLAYER(LAYER) + NPTS + 1 
!****              BLOCK TO COMPUTE ETASBX(MN)                          
      DO 95 J = 1,LAYER 
      J1 = LAYER+1-J 
      NPTS2 = PTSIN(J1)+2. 
      DO 90 K = 1,NPTS2 
      MN = JLAYER(J1)+K 
   90 ETASBX(MN) = 1./THIKN(J1) 
   95 END DO 
!****              TEST FOR ROUND OFF IN THE MACHINE.                   
      IF (GENRL(4)-WARM(1)) 97,97,98 
   97 WARM(1) = WARM(1)*.99999 
   98 NPTS3 = PTSIN(LAYER)+2. 
      KONVAR(1) = JLAYER(LAYER)+NPTS3 
!****              MAKE SURE TNEST IS CLEARED                           
      CALL TNEST(ZERO,ZERO,ZERO) 
      GO TO 300 
!****              BLOCK TO PREPARE FOR NEXT TIME STEP.                 
!****              ALTERNATE KONVAR(3) BETWEEN 2 AND 3.                 
  150 KONVAR(3) = 5-KONVAR(3) 
!****              UPDATE OR RE-INITIALIZE COUNTERS.                    
      COUNT(1) = COUNT(1)+1. 
      COUNT(2) = COUNT(2)+COUNT(3) 
      COUNT(3) = 0. 
      COUNT(4) = COUNT(4)+COUNT(5) 
      COUNT(5) = 0 
      COUNT(6) = COUNT(6)+COUNT(7) 
      COUNT(7) = 0. 
      COUNT(8) = COUNT(8)+COUNT(9) 
      COUNT(9) = 0. 
      COUNT(10) = COUNT(10)+COUNT(11) 
      COUNT(11) = 0. 
      COUNT(12) = COUNT(12)+COUNT(13) 
      COUNT(13) = 0. 
  155 IF( KONVAR(2) .NE. 1) GO TO 159 
  156 KONVAR(2) = KONVAR(2)+1 
      GO TO 300 
  159 KPTS = 2*KONVAR(1) 
      DO 160 JTRANS = 1, KPTS 
      MTIME = KONVAR(1)+JTRANS 
      DENS(JTRANS) = DENS(MTIME) 
      P(JTRANS) = P(MTIME) 
      T(JTRANS) = T(MTIME) 
  160 END DO 
      IF(KONVAR(2) .NE. 3) GO TO 175 
  161 DO 164 J = 2,3 
      K = J-1 
      DELTO(K) = DELTO(J) 
      SMELT(K) = SMELT(J) 
  164 SCHAR(K) = SCHAR(J) 
      SCHAR(3) = SCHAR(4) 
      SMELT(3) = SMELT(4) 
  167 DO 168 I = 1,12 
      DO 168 J = 2,4 
      K = J-1 
      Q(K,I) = Q(J,I) 
  168 CONTINUE 
      GO TO 300 
  175 KONVAR(2) = 3 
      GO TO 300 
!****              BLOCK TO PREPARE NEXT CASE.                          
  190 CALL OUTPUT(1) 
!****              FIRST PRINTOUT OUTPUT.                               
      THIKN(1) = HTCONV(9) 
  191 DO 192 I = 1,5 
      CHARV(I) = 0.0 
  192 SPACVR(I) = 0.0 
      DO 195 I = 1,4 
      SMELT(I) = 0.0 
      DELTO(I) = 0.0 
      SCHAR(I) = 0.0 
      DO 193 J = 1,13 
  193 Q(I,J) = 0.0 
      DO 194 K = 1,150 
  194 BETA(I,K) = 0.0 
  195 END DO 
      DO 196 I = 1,25 
      AUS(I) = 0.0 
      COUNT(I) = 0.0 
  196 KONVAR(I) = 0 
      DO 197 I = 1,150 
      EMDG(I) = 0.0 
      ETASBT(I) = 0.0 
      ETASBX(I) = 0.0 
      ETAXX(I) = 0.0 
      WDOTP(I) = 0.0 
  197 XLOC(I) = 0.0 
      DO 198 I = 1,450 
      DENS(I) = 0.0 
      P(I) = 0.0 
  198 T(I) = 0.0 
      DO 199 I = 1,10 
      DELETA(I) = 0.0 
      TIVAR(I) = 0.0 
      HEVAR(I) = 0.0 
  199 JLAYER(I) = 0 
      DO 203 I = 4,25 
  203 HTCONV(I) = 0.0 
      IMPTME = 0 
      IF ( L .NE. 4) GO TO 302 
      ALLNEW = 0.0 
      GENRL(4) = 0.0 
      NTXQ(1) = 0 
      NTXQ(2) = 0 
      NOFXS = 0.0 
      NOFTS = 0.0 
      END = 0. 
      DO 204 I = 1,404 
  204 CHART(I) = 0.0 
      DO 205 I = 16,45 
  205 AGPTB(I) = 0.0 
      DO 206 I = 8,3378 
  206 EPTAB(I) = 0.0 
      DO 207 I = 1,9 
  207 CHARV(I) = 0.0 
      DO 208 I = 1,25 
  208 COUNT(I) = 0.0 
      DO 209 I = 1,4 
      DO 209 JD = 1,13 
  209 Q(I,JD) = 0.0 
      DO 210 I = 1,21 
  210 DELETA(I) = 0.0 
      DO 211 I = 1,54 
  211 SMELT(I) = 0.0 
      DO 212 I = 1,15 
  212 DELTO(I) = 0.0 
      DO 213 I = 1,4 
      DO 213 JD = 1,150 
  213 BETA(I,JD) = 0.0 
      DO 1214 I = 1,2250 
 1214 DENS(I) = 0.0 
      CTIME = 0.0 
      TIME = 0.0 
      DO 214 I = 4,50 
  214 HTCONV(I) = 0.0 
      DO 215 I = 1,23 
  215 TIVAR(I) = 0.0 
      DO 216 I = 1,25 
  216 AUS(I) = 0.0 
      DO 217 I = 1,16 
  217 AMPLCT(I) = 0.0 
      WRITE ( 6,1001) 
  300 KONVAR(6) = MIN0(2,KONVAR(2)) 
      RETURN 
 1001 FORMAT('1') 
  302 RETURN 
END Subroutine Init          
                                 
                                           
!$IBFTC INTPTS                                                          
      SUBROUTINE INTPTS(KT,NT) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE INTERIOR POINTS OF EACH LAYER ARE DISCUSSED IN        
!****         SECTION IV OF THE ANALYSIS BY P. GORDON                   
      J2 = KT 
      KEX = KT 
      N1 = NT 
      IT = KONVAR(1) 
      KTIM = KONVAR(2) 
!****          WHEN                                                     
!****              KEX = 0.,USE EXPLICIT EQUATIONS.                     
!****              KEX = -1,USE IMPLICIT EQUATIONS.                     
!****              KEX IS GREATER THAN ZERO, EVALUATE ONE               
!****                   OR MORE POINTS BY THE IMPLICIT                  
!****                   EQUATIONS.                                      
      IF (IMEX) 300,4,10 
    4 IF(KEX) 100,5,191 
    5 IF (KONVAR(2)-1) 50,10,50 
!****              NOW USE THE EXPLICIT SCHEME FOR                      
!****                   THE FIRST TIME STEP.                            
   10 K = 1 
      L = KONVAR(4) 
      DO 40 J = 1,L 
!****         RZ IS DEFINED BY EQUATION IV.4                            
   12 RZ = DELTO(KTIM)/DELETA(J)**2 
!****              STORES THE FIRST GUESS OF THE INTER-FACE TEMPERATURES
      NPTS = PTSIN(J)+1. 
      JJ1 = JLAYER(J) + 1 
      IJ1 = JJ1 + KONVAR(1)*KONVAR(6) 
      T(IJ1) = T(JJ1) 
      IJ2 = IJ1 + NPTS 
      JJ2 = JJ1 + NPTS 
      IF(J-KONVAR(5)) 15, 39, 15 
   15 DO 30 I = 2,NPTS 
      MN = JLAYER(J) + I 
      II = KONVAR(1)*KONVAR(6) + I + JLAYER(J) 
      JJ = II-KONVAR(1) 
!***          R1 IS DEFINED BY EQUATION IV.5                            
      R1 = (T(JJ+1)-T(JJ-1))/2./DELETA(J) 
!***          R2 IS DEFINED BY EQUATION IV.56                           
      R2 = T(JJ+1)+T(JJ-1) 
!****              CALL CHARP TO CALCULATE CPA,AK,AKP.                  
      CALL CHARP(JJ,J) 
!****         IF CHARRING IN THE FIRST LAYER COMPUTE M-DOT-G            
      IF (CHART(1)) 17,21,17 
   17 IF(J-1) 21, 20, 21 
   20 CALL EMG(JJ) 
!****         COMPUTE THE BETA COEFFICIENTS OF THE T EQUATION           
   21 CALL BEE(JJ,MN,J,ITTEMP) 
      BB3 = BETA(3,MN) + BETA(2,MN)*R1 
      IF (BB3) 25,23,25 
   23 ALPHA = .5 
      GO TO 26 
   25 RBDN = BETA(1,MN)/(ABS(BB3))/DELETA(J) 
      ALPHA = AMIN1(.5,RBDN) 
   26 IF (BB3) 27,28,28 
   27 ALF1 = ALPHA 
      ALF2 = 1.-ALPHA 
      GO TO 29 
   28 ALF1 = 1.-ALPHA 
      ALF2 = ALPHA 
!***          R3 IS DEFINED BY EQUATION IV.8                            
   29 R3 = T(JJ)+RZ*((BETA(1,MN)+ALF1*DELETA(J)*BB3)*T(JJ+1)+(BETA(1,MN)&
     &-ALF2*DELETA(J)*BB3)*T(JJ-1))+DELTO(KTIM)*BETA(4,MN)              
!****              NOW, COMPUTE THE TEMPERATURE.                        
!***          T(II) IS DEFINED IN EQUATION IV.9                         
      T(II) = RZ*((ALF2-ALF1)*DELETA(J)*BB3-2.*BETA(1,MN))*T(JJ)+R3 
   30 END DO 
!****              STORES THE 1ST GUESS FOR THE BACK-FACE TEMPERATURE   
   39 T(IJ2) = T(JJ2) 
   40 END DO 
      GO TO 200 
   50 DELTTR = DELTO(KTIM)/DELTO(KTIM-1) 
      K = 1 
      DO 60 I = 1,IT 
      II = 2*KONVAR(1)+I 
      JJ = II-KONVAR(1) 
!****              EXPLICIT EQUATION FOR AFTER ONE                      
!****                   TIME STEP.                                      
!***          T(II) IS DEFINED IN EQUATION IV.10                        
      T(II) = (DELTTR+1.)*T(JJ)-DELTTR*T(I) 
   60 END DO 
      GO TO 200 
!****              EFN 100 STARTS THE IMPLICIT SCHEME.                  
  100 J1 = 1 
      K = 0 
      J2 = KONVAR(4) 
      N2 = 2 
      IF (N1-2) 102,101,102 
  101 N1 = N1+2 
  102 DO 190 J = J1,J2 
      IF (J-KONVAR(5)) 103,190,103 
!***          RZ IS DEFINED IN EQUATION IV.4                            
  103 RZ = DELTO(KTIM)/DELETA(J)**2 
      NPTS1 = PTSIN(J) 
      NPTS = NPTS1+4-N1 
  105 DO189 I = N1,NPTS,N2 
      IF (KEX) 106,106,110 
  106 MN = JLAYER(J)+NPTS1+3-I 
      II = KONVAR(6)*KONVAR(1)+MN 
      GO TO 115 
  110 IF (J-J1) 115,106,114 
  114 II = II+N6 
      MN = MN+N6 
  115 JJ = II-KONVAR(1) 
!***          R1 IS DEFINED IN EQUATION IV.5                            
      R1 = (T(II+1)-T(II-1))/2./DELETA(J) 
!***          R2 IS DEFINED IN EQUATION IV.6                            
      R2 = T(II+1)+T(II-1) 
      T1 = T(II) 
      EPS = -EPTAB(5) 
!****              SAVE TEMPERATURE FOR TNEST IN T1.                    
!****              CALL CHARP FOR CPA,AK,AKP                            
      CALL CHARP(II,J) 
!****              TEST TO SEE IF CHARRING.                             
  116 IF (CHART(1)) 117,121,117 
  117 IF (J-1) 121,118,121 
  118 CALL DENSIT(II,MN,0) 
      IF(KONVAR(7)-2) 120, 119, 120 
  119 CALL PRESSM(II,MN,J,0) 
!****         IF CHARRING COMPUTE M-DOT-G                               
  120 CALL EMG(II) 
!****              CALL BEE TO COMPUTE BETA-S.                          
  121 CALL BEE(II,MN,J,ITTEMP) 
      BB3 = BETA(3,MN)+BETA(2,MN)*R1 
      IF (BB3) 169,168,169 
  168 ALPHA = .5 
      GO TO 170 
  169 RBDN = BETA(1,MN)/(ABS(BB3))/DELETA(J) 
      ALPHA = AMIN1(.5,RBDN) 
  170 IF (BB3) 171,172,172 
  171 ALF1 = ALPHA 
      ALF2 = 1.-ALPHA 
      GO TO 174 
  172 ALF1 = 1.-ALPHA 
      ALF2 = ALPHA 
!***          R3 IS DEFINED IN EQUATION IV.8                            
  174 R3 = T(JJ)+RZ*((BETA(1,MN)+ALF1*DELETA(J)*BB3)*T(II+1)+(BETA(1,MN)&
     &-ALF2*DELETA(J)*BB3)*T(II-1))+DELTO(KTIM)*BETA(4,MN)              
!****              NOW COMPUTE THE TEMPERATURE.                         
!***          T(II) IS DEFINED BY EQUATION IV.9                         
      T(II) = R3/(1.+RZ*(2.*BETA(1,MN)+(ALF1-ALF2)*DELETA(J)*BB3)) 
!****              USE TNEST TO ITERATE FOR THE                         
!****              FINAL TEMP. VALUE. ITTEMP IS                         
!****              A CONTROL SET IN BEE TO DETERMINE                    
!****              IF THE ITERATION IS NECESSARY.                       
!****         STORE ITTEMP IN ANOTHER CELL, IN ORDER TO KEEP FROM       
!****         OSCILLATING ON ITERATIONS, DUE TO BETA2 IN BEE.           
  175 IF(EPS) 176,2176, 2176 
  176 ITTP1 = ITTEMP 
 2176 IF(ITTP1) 177, 183, 177 
  177 CALL TNEST(T(II),T1,EPS) 
      KT1 = T1 
      GO TO (183,178,179),KT1 
  178 EPS = ABS(EPS) 
      T(II) = AMAX1(WARM(1),T(II)) 
      COUNT (7) = COUNT(7)+1. 
      GO TO 116 
  179 WRITE(6,500)   T(II),ALF1,ALF2,BETA(1,MN),BETA(2,MN),BETA(3,MN),  &
     &BETA(4,MN),MN                                                     
      GO TO 116 
!****              T1 IS CHANGED IN TNEST                               
!****                   1.CONVERGENCE                                   
!****                   2.NO-CONVERGENCE                                
!****                   3.ERROR RETRACE                                 
!****         COMPUTE M-DOT-G IF CHARRING                               
  183 IF (CHART(1)) 184,189,184 
  184 IF (J-1) 189,185,189 
!****              GO TO EMG TO GET VALUES FOR POINTS                   
!****                   SKIPPED IN IMPLICIT SCHEME. IF NO               
!****                   ITERATION ON T.                                 
!****              COMPUTE CORRECT M-DOT-G.                             
  185 IF(ITTP1) 187, 186, 187 
  186 K7 = II 
      GO TO 188 
  187 K7 = II-1 
  188 CALL EMG(K7) 
      IF (K7-II+1) 187,189,187 
  189 END DO 
  190 END DO 
      GO TO 200 
!****              COMPUTE INDICES FOR COMPUTATION                      
!****              OF ONE OR TWO INTERFACE POINTS.                      
!****                   BLOCK TO COMPUTE LAYER INDICES.                 
  191 K = 0 
      IF (KEX-1) 193,192,193 
  192 J1 = J2 
!****              FOR FRONT FACE, POINT N1 IS INPUT                    
!****                   TO THIS ROUTINE AND EQUALS                      
!****                   PTSIN(J)+1.                                     
      N6 = 2 
      N2 = PTSIN(J1)+4. 
      GO TO 102 
  193 IF(KEX-1-KONVAR(4)) 195, 194, 195 
  194 J1 = J2-1 
      J2 = J1 
!****              FOR THE BACKFACE, POINT N1 IS INPUT                  
!****                   TO THIS ROUTINE AND = 2.                        
      N2 = PTSIN(J1)+4. 
      N6 = 2 
      GO TO 102 
  195 IF (KEX-KONVAR(5)) 196,197,196 
  196 J1 = J2-1 
!****              FOR INTERIOR POINTS N1 = 2.                          
      N2 = PTSIN(J1)+4.+PTSIN(J2) 
      N6 = 2 
      GO TO 102 
  197 N6 = 3 
!****              FOR AIR-GAP,KEX = KONVAR(5)                          
!****                   AND N1 = 2.                                     
      J1 = KONVAR(5)-1 
      N2 = PTSIN(J1)+4. 
      IF (KONVAR(5)-KONVAR(4)) 198,102,198 
  198 J2 = J2+1 
      N2 = PTSIN(J1)+PTSIN(J2)+4. 
      GO TO 102 
  200 RETURN 
  300 CALL TMPIMP(1) 
      GO TO 200 
  500 FORMAT (' BEE EFN 177 T = ',E20.8,3X,'A1 = ',E20.8,3X,'A2 = ',E20.8// &
        3X,'B1 = ',E20.8,3X,'B2 = ',E20.8,3X,'B3 = ',E20.8,3X,'B4 = ',E20.8,'I = ',I3)
      STOP 
END Subroutine Intpts
                                          
!$IBFTC KEFFS                                                           
      SUBROUTINE KEFFS(III,MNI) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
      II = III 
!****              THIS ROUTINE COMPUTES K(EFF)                         
!****                K(EFF)-T,K(EFF)-RHO(S)-PRIME AND                   
!****                STORES THEM IN SCDTB(1-3)                          
      MN = MNI 
      RR = TABLD(14)/DENS(II) 
      RR2 = RR/DENS(II) 
      RTRM = RR-1. 
      SCDTB(1) = GRAPH(1) + GRAPH(2) + GRAPH(3)*RTRM 
      SCDTB(2) = 0.0 
      SCDTB(3) = -GRAPH(3)*RR2 
      HTCONV(15) = SCDTB(1) 
  201 RETURN 
END Subroutine Keffs
                                          
!$IBFTC KFS                                                             
      SUBROUTINE KFS(III,MNI) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
      SCDTB(12) = GRAPH(7) 
      SCDTB(13) = 0.0 
      SCDTB(14) = 0.0 
      SCDTB(15) = 0.0 
      SCDTB(16) = 0.0 
!****              K(F) IS STORED IN SCDTB(12)                          
!****              K(F)-ETA IS STORED IN SCDTB(13)                      
!****              K(F)-ETA-ETA IS STORED IN SCDTB(14).                 
      RETURN 
END Subroutine Kfs     


                                          
!$IBFTC MOMENT                                                          
      SUBROUTINE MOMENT(XARAY) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION XARAY(1) 
!                       THIS SUBROUTINE CALCULATES THE                  
!                       MOMENTS OF THE BODY.                            
!                       INPUT ARRAY COMMON.                             
!                       XSTAR IS THE DISTANCE, FROM THE                 
!                       FIRST POINT OF THE BODY, TO THE                 
!                       TEMPERATURE WHERE T= OR IS GREATER              
!                       THAN T(AB).                                     
!                       X IS THE DISTANCE, FROM THE FIRST               
!                       POINT OF THE BODY, TO THE POINT IN              
!                       QUESTION.                                       
!                       THE LOCATIONS STRESS (21-26) ARE                
!                       WHERE THE INTEGRATED VALUES ARE                 
!                       STORED.                                         
!                       THE LOCATIONS STRESS (11-16) ARE                
!                       WHERE THE SUM OF THE INTEGRANDS                 
!                       ARE STORED.                                     
!                       THE LOCATIONS STRESS (1-7) ARE                  
!                       WHERE THE INDIVIDUAL VALUES OF THE              
!                       INTEGRAND ARE STORED.                           
!                       EMOD(I) UNITS ARE IN INCHES.                    
!                       ALPHA T IS DIMENSIONLESS.                       
!                      T(REF) IS STORED IN STVAR(2).                    
!                       T(AB) IS STORED IN STVAR(3).                    
      NTEMP = INDTL(4) 
      ITIME = KONVAR(1) 
      LAYER = KONVAR(4) 
      IF (GENRL(11)) 1,3,1 
    1 DO 2 I = 1,ITIME 
    2 XARAY(I) = XARAY(I)/12. 

    3 N3 = KONVAR(6) 
      XLEN = 0. 
!****              TOTAL LENGTH OF LAYERS IS IN TLEN.                   
      TLEN = 0. 
      DO 4 J = 1,LAYER 
      TLEN = TLEN + THIKN(J) 
    4 END DO 

      DO 5 I = 1,30 
    5 STRESS(I) = 0. 

      XSTAR = 0. 
!                       BLOCK TO COMPUTE XSTAR.                         
!                       (WORK FROM BACK-FACE FORWARD.)                  
      DO 50 I = 1,ITIME 
      IT = (N3+1)*ITIME +1 -I 
      IF(T(IT) - STVAR(3)) 50, 60, 60 
   50 END DO 

      XSTAR = 0. 
      GO TO 115 
   60 IT = IT - N3*ITIME 
!                       THE POINT IS STORED IN IT.                      
!                       NOW FIND THE LAYER.                             
      DO 70 J = 1,LAYER 
!!!      J = J 
      IF (JLAYER(J)-IT) 70,80,80 
   70 END DO 

      JT = LAYER 
      GO TO 81 
!                       THE LAYER IS STORED IN JT.                      
!                       NOW, COMPUTE ETA.                               
   80 JT = J-1 
   81 IT2 = IT-JLAYER(JT) 
      IF(JT-1) 90, 85, 90 
   85 XSTAR = XARAY(IT) 
      GO TO 115 
   90 JT1 = JT-1 
      DO 92 J = 1,JT1 
      XLEN = XLEN + THIKN(J) 
   92 END DO 
      XSTAR = XLEN + XARAY(IT) 
  115 IF (XSTAR) 119,118,119 
  118 IT2 = 1 
      JT = 1 
  119 XLEN = 0.0 
      DO 210 J4 = JT,LAYER 
      NPTS = PTSIN(J4) +2. 
      IF(J4-KONVAR(5)) 120, 210, 120 
  120 IF(J4-JT-1) 122,121,121 
  121 IT2 = 1 
  122 DO 199 I4 = IT2,NPTS 
      XB = 0.0 
      MN = JLAYER(J4)+I4 
      II = N3*ITIME+MN 
      XB = (XLEN + XARAY(MN)-XSTAR)*12. 
      IF(J4-LAYER) 124,124,129 
  124 IF(J4-1) 129,129,125 
  125 IF(I4-1) 129,126,129 
  126 XLEN = XLEN + THIKN(J4-1) 
  129 NLA = (J4-1) * 4 * (NDGREE + 1) + 1 
!                       NOW COMPUTE THE MOMENT INTEGRALS.               
      NLA = (J4-1) * 4 * (NDGREE + 1) + 1 
      NDGRE2 = NLA + 2 * (NDGREE + 1) 
      CALL TabUpSingle(DGREE,ATETB(NLA:), STVAR(2), ATR, NDGREE,1)  ! modified by RLC
      CALL TabUpSingle(DGREE,ATETB(NLA:),T(II),ATTAB,NDGREE,1)   ! modified by RLC
      CALL TabUpSingle(DGREE,ATETB(NDGRE2:),T(II),TABE,NDGREE,1)   ! modified by RLC
      IF(J4-1) 135, 130, 135 
  130 STRESS(1) = 1./ETASBX(MN) 
      GO TO 137 
  135 STRESS(1) = THIKN(J4) 
  137 STRESS(2) = TABE*STRESS(1) 
      STRESS(3) = STRESS(2)*XB 
      STRESS(4) = STRESS(3)*XB 
      STRESS(7) = ATTAB-ATR 
      STRESS(5) = STRESS(2)*STRESS(7) 
      STRESS(6) = STRESS(5)*XB 
      IF (I4-IT2) 150,151,150 
  150 IF(I4-NPTS) 155, 151, 155 
  151 PT4 = 1. 
      GO TO 160 
  155 PT4 = 2. 
!                       MAKE THE TRAPEZOIDAL RULE SUM.                  
  160 STRESS(12) = STRESS(2)*PT4+STRESS(12) 
      STRESS(13) = STRESS(3)*PT4+STRESS(13) 
      STRESS(14) = STRESS(4)*PT4+STRESS(14) 
      STRESS(15) = STRESS(5)*PT4+STRESS(15) 
      STRESS(16) = STRESS(6)*PT4+STRESS(16) 
  199 END DO 
!****              NOW DIVIDE BY 2.0 AND MULT. BY DELETA(J4)            
!****              COMPUTE TOTAL SUM AND CONVERT TO INCHES.             
      DEN = DELETA(J4) * 12./2. 
  200 STRESS(22) = STRESS(12)*DEN+STRESS(22) 
      STRESS(23) = STRESS(13)*DEN+STRESS(23) 
      STRESS(24) = STRESS(14)*DEN+STRESS(24) 
      STRESS(25) = STRESS(15)*DEN+STRESS(25) 
      STRESS(26) = STRESS(16)*DEN+STRESS(26) 
!                       CLEAN OUT STRESS(11-16) FOR                     
!                       NEXT LAYER.                                     
      DO 205 J8 =  1,16 
      STRESS(J8) = 0. 
  205 END DO 
  210 END DO 
!                       NOW SET UP OUTPUT                               
      STRESS(31) = STRESS(25) 
      STRESS(32) = STRESS(22) 
      STRESS(33) = STRESS(26)-(STRESS(23)*STRESS(25))/STRESS(22) 
      STRESS(34) = STRESS(24)-(STRESS(23)**2)/STRESS(22) 
      STRESS(35) = (TLEN - XSTAR) * 12. 
      STRESS(36) = STRESS(23)/STRESS(22) 
      STRESS(37) = STRESS(33)/STRESS(34) 
      STRESS(38) = STRESS(25)/STRESS(22) 
      STRESS(39) = SQRT(12.*STRESS(34)/STRESS(22)) 
      STRESS(40) = STRESS(37)/2. 
      STRESS(41) = STRESS(38)+STRESS(40)*STRESS(39) 
      STRESS(42) = STRESS(38)-STRESS(40)*STRESS(39) 
      STRESS(43) = STRESS(22)/STRESS(39) 
!                       NOW WRITE OUT OUTPUT                            
      WRITE(6,500)STRESS(31),STRESS(32),STRESS(33),STRESS(34),STRESS(35)&
     &,STRESS(36),STRESS(37),STRESS(38),STRESS(41),STRESS(42),STRESS(39)&
     &,STRESS(43)                                                       
!                       FORMAT STATEMENTS
  500 FORMAT(/,  &
   4X,'I(EAT) = ',E15.8,3X,'I(E) =   ',E15.8,3X,'I(EATZ) =  ',E15.8// &
   5X,'I(EZZ) = ',E15.8,4X,'T(ACT) = ',E15.8,3X,'CENTROID = ',E15.8// &
   5X,'MBAR =   ',E15.8,3X,'PBAR =   ',E15.8,4X,'ALPHA T(OUTER) = ',E15.8// &
   5X,'ALPHA T(INNER) = ',E15.8,3X,'T(EFF) = ',E15.8,3X,'E(EFF) = ',E15.8)
  216 RETURN 
      STOP 
      END Subroutine Moment
                                          
!$IBFTC NSUBT                                                           
      SUBROUTINE NSUBT 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATION FOR THE PARTIAL DERIVATIVE OF ETA WITH       
!****         RESPECT TO TIME IS DEFINED IN APPENDIX B                  
      LT = PTSIN(1)+2. 
      C = SPACE(6) 
      LAYER = KONVAR(4) 
      ETASBT(1) = -ETASBX(1)* HTCONV(5) 
      DO 100 J = 1,LAYER 
      NPTS = PTSIN(J)+1. 
      DO 80 I = 1,NPTS 
      MN = JLAYER(J)+I+1 
      IF (J-2) 10,75,75 
   10 TEN = I-1 
      EN = TEN*DELETA(J) 
!****         XBAR IS DEFINED IN EQUATION B.1                           
      XBAR = ABS(XLOC(MN) - XLOC(1))/THIKN(1) 
      IF (SPACVR(1)) 20,15,20 
   15 ENR = EN*(1.-EN)/2. 
      GO TO 50 
   20 AOR = THIKN(1)/SPACVR(1) 
!****         ENR, ETA-SUB-R, IS DEFINED IN EQUATION B.6.2              
      ENR = AOR*(XBAR*ETASBX(MN)/(1.+C*XBAR*(2.*XBAR-1.))-EN*ETASBX(LT)/&
     &(1.+C))                                                           
!****         ETASBT, ETA-SUB-T, IS DEFINED IN EQUATION B.6.1           
   50 ETASBT(MN) = ENR*SPACVR(2)-ETASBX(MN)*(1.-XBAR)*HTCONV(5) 
      GO TO 80 
   75 ETASBT(MN) = 0. 
   80 END DO 
  100 END DO 
      ETASBT(LT) = 0.0 
  200 RETURN 
END Subroutine Nsubt                   
                       
!$IBFTC OUTPUT                                                          
      SUBROUTINE OUTPUT (KERR) 

USE HeadG,ONLY: head
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TRAY(150) 
      DIMENSION XRAY(175) 
      IF (CTIME .EQ. 0.0  .OR.  CTIME .NE. PRVTME) GO TO 8 
      RETURN 

    8 PRVTME = CTIME 
      N3 = KONVAR(6) 
      KTIME = KONVAR(2) 
      LAYER = KONVAR(4) 
      ITIME = KONVAR(1) 
      IF(GENRL(11))  13,  14,  13 
   13 FTORIN = 12.0 
      GO TO 15 
   14 FTORIN = 1. 
   15 IF (KERR) 39, 20, 39 
   20 IF(KTIME -1) 25, 24, 25 
   24 FRSTM = GENRL(1) 
      KPAGE=0 
   25 DO 30 I = 2,18,2 
      IJ = I 
      IF(CTIME-OUTTB(I)) 33, 33, 30 
   30 END DO 
   33 IJ = IJ 
      L = 0 
      IF((CTIME-FRSTM) - OUTTB(IJ-1)) 350, 37, 37 
   37 FRSTM = CTIME 
      L=1 
   39 WRITE (6,'(1X,A)') head
      KPAGE=KPAGE+1 
      WRITE(6,501)KPAGE 
      WRITE(6,502) CTIME,DELTO(KTIME) 
      WRITE(6,513) 
      IF(CHART(1)) 42, 40, 42 
   40 SCHAR(KTIME + 1) = 0.0 
   42 SC  = FTORIN * HTCONV(8) 
      SPC = FTORIN * HTCONV(4) 
      SM  = FTORIN * HTCONV(7) 
      SPM = FTORIN * HTCONV(5) 
      WRITE(6,503)  SC,SPC,SM,SPM,EMDG(1),                              &
     &SPACVR(1),SPACVR(2),SPACVR(3)                                     
      WRITE(6,514) 
      WRITE(6,504)(TIVAR(I),I=1,10) 
   53 IF(OUTTB(19)) 57, 55, 57 
   55 WRITE(6,515) 
      WRITE(6,505)(Q(KTIME+1,I),I=7,12),Q(4,13) 
   57 WRITE(6,516) 
      WRITE(6,505) (Q(KTIME+1,I),I = 1,6),Q(1,13),Q(2,13) 
      NXS = 1 
      XRAY(ITIME + 1) = 1. 
      NPSL1 = PTSIN(1) + 2. 
   62 DO 149 J = 1,LAYER 
      NPSI = PTSIN(J) + 2. 
      NXS = NXS + NPSI-1 
!     WRITE OUT THE HEADING FOR THE COORDINATES                         
      JL1 = JLAYER(J) + 1 
      IF(JL1 - 1) 113, 112, 113 
  112 JPLA = 0 
      GO TO 114 
  113 XRAY(JL1) = 0.0 
      JPLA = NXS-NPSI 
  114 JPLA1 = JPLA + 1 
      DO 115 I = JPLA1,NXS 
      K = I + 1 
      XRAY(K) = 0.0 
      IF(J-1) 2116,2115,2116 
 2115 XRAY(I) = FTORIN * ABS(XLOC(I) - TRAJV(50)) 
      GO TO 115 
 2116 XRAY(I) = FTORIN * ABS (XLOC(I) - XLOC(JPLA1)) 
  115 END DO 
      IF(J-1) 116, 116, 121 
  116 XLAST = XRAY(NPSL1) 
  121 WRITE(6,563) J 
      NJL1 = JL1 + NPSI -1 
      WRITE(6,566)(XRAY(JJ),JJ = JL1,NJL1) 
!     WRITE OUT HEADING FOR TEMPERATURES                                
      WRITE(6,565) J 
      IK = ITIME*N3+JLAYER(J) + 1 
      KI=IK+NPSI-1 
      IF(J.EQ.1) JDS=IK 
      WRITE(6,566)(T(K),K = IK,KI) 
      IF (J-1) 149, 126, 149 
  126 IF(CHART(1)) 127, 146, 127 
!     WRITE OUT THE DENSITY AND HEADING IF DESIRED                      
  127 WRITE(6,567) J 
      WRITE(6,566)(DENS(K),K=IK,KI) 
  146 IF(KONVAR(7)-2) 149, 145, 149 
  145 WRITE(6,568) J 
      WRITE(6,566)(P(K),K=IK,KI) 
  149 END DO 
      NPTS=NXS 
      JE=JDS+NPTS-1 
      IF(OUTTB(20)) 310, 200, 310 
  200 NX = PTSIN(1) + 2.0 
      IKI = ITIME*N3 + 1 
      IT = IKI + ITIME 
      T(IT) = 1. 
      XXX=0.0 
      XRAY(NPSL1) = XLAST 
      NXM1 = NX - 1 
      DO 220 K = 1,NXM1 
      IF(XXX-XRAY(1)) 210,215,215 
  210 TRAY(K)=0.0 
      GO TO 220 
  215 IF(K-NX) 216,222,222 
  216 IF(XXX-XRAY(NX)) 218, 218, 217 
  217 STOP   ! changed from CALL EXIT by RLC 7Sep05
  218 CALL TabUpSingle(XRAY,T(IKI:),XXX,TRAY(K),ITIME,1) 
  220 XXX = XXX + DELETA(1) * HTCONV(9) * FTORIN 
  222 LASTT = ITIME*N3 + NX 
      TRAY(NX)=T(LASTT) 
      WRITE(6,570) 
      WRITE(6,566)(TRAY(JKL),JKL=1,NX) 
      WRITE(6,552)(COUNT(N),N = 1,18) 
  310 IF(STVAR(1)) 315, 330, 315 
  315 CALL MOMENT(XRAY) 
!****         TEST IF INEQUALITY HAS BEEN VIOLATED.                     
  330 IF(AUS(3)) 335,340,340 
  335 WRITE(6,585) 
!****         RESET AUS(3) FOR OUTPUT.                                  
      AUS(3) = -AUS(3) 
  340 WRITE(6,590) 
      WRITE(6,566) (AUS(JD),JD = 1,8) 
      WRITE(6,600) IMEX 
!****          SEE IF A PLOTTING TAPE HAS BEEN REQUESTED.               
  350 IF(GENRL(13)) 400,400,360 
!****          AT EFN 360 A TAPE FOR THE S.C.--4020 PLOTTING            
!****          ROUTINES HAS BEEN REQUESTED.                             
!****          THE TAPE FORMAT IS FORTRAN 4  BINARY.  THE FIRST RECORD  
!****          IS WRITTEN ONCE AT THE START OF THE RUN-----FOLLOWING    
!****          THIS AT EVERY NORMAL OUTPUT INTERVAL THERE WILL          
!****          BE WRITTEN 1 RECORD OF INFORMATION.                      
!****          THIS RECORD WILL CONTAIN ALL THE DATA NECESSARY FOR      
!****          CONSTRUCTING ANY OF THE PLOTS NORMALLY REQUESTED.        
  360 IF(L.EQ.0) GO TO 400 
      IF(KPAGE.GT.1) GO TO 380 
!****          THIS IS THE SINGLE RECORD FOR THE FIRST TIME.            
  370 WRITE( 9) LAYER,GENRL(11),GENRL(1),GENRL(2),(THIKN(JD),JD=1,LAYER)&
     &,(DENSE(JD),JD=1,LAYER)                                           
!****          THIS IS FOR THE SUBSEQUENT DATA RECORDS.                 
  380 WRITE( 9)  CTIME,Q(4,13),SM,SC,(Q(KTIME+1,JD),JD=1,6),Q(1,13),    &
     &Q(2,13),EMDG(1),NPTS,Q(KTIME+1,7),(PTSIN(JD),JD=1,LAYER),(XRAY(JD)&
     &,JD=1,NPTS),(T(JD),JD=JDS,JE),(DENS(JD),JD=JDS,JE)                
      L = 0 
  400 RETURN 
  501 FORMAT(112X,'PAGE',I3) 
  502 FORMAT('     TIME=',F14.8,'  DELTA TIME=',F13.8) 
  503 FORMAT(1X,8E12.4) 
  504 FORMAT(1X,5E12.4//) 
  505 FORMAT(1X,8E12.4) 
  513 FORMAT(/,2X,'AMT CHAR',4X,'RATE CHAR',3X,'AMT MELT',4X,'RATE MELT', &
   4X,'MDG(1)',4X,'AR1(1)',6X,'AR1P(1)',5X,'ETA0(1)'  )           
  514 FORMAT(/,13X,'TIME TABLE ENTRIES 1 THROUGH 10')
  515 FORMAT(/, &
   'INT QNET     INT QRR    INT QGHR     INT QC      INT QBLK   INT QBKFAC    QI TOTAL')
  516 FORMAT(/,4X,'QDNET',7X,'QDRR',8X,'QGHR',9X,'QC',9X,'QDBLK',9X, &
     'QBKFAC',4X,'QVAP',7X,'QC *')                                     
  552 FORMAT(/,2X,'COUNT = ',1P6E13.5/(11X,ES13.5) ) 
  554 FORMAT(4E15.8) 
  563 FORMAT(/,4X,'COORDINATE',1X,'VALUES',3X,'LAYER=',I2) 
  565 FORMAT(/,4X,'TEMPERATURES',4X,'LAYER=',I3) 
  566 FORMAT(/(1PE16.5,7E14.5)) 
  567 FORMAT(/,4X,'DENSITYS',8X,'LAYER=',I2) 
  568 FORMAT(/,4X,'GAS DENSITIES',6X,'LAYER = ',I2) 
  570 FORMAT(/,4X,'EQUALLY SPACED TEMPERATURES   FIRST LAYER') 
  580 FORMAT(1H1,12A6) 
  585 FORMAT(/,4X, &
   '**********CONDITION VIOLATED RUN CONTINUED**********'///)
  590 FORMAT(/,4X,'AUS  TABLE') 
  600 FORMAT(/,10X,'*****',2X,'METHOD OF COMPUTATION WAS ',I2)
      STOP 
END Subroutine Output
                                          
!$IBFTC PBETA                                                           
      SUBROUTINE PBETA (III,MNI,JJJ,ITM) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      DIMENSION GVAR(6) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
!                              THIS SUBROUTINE CALCULATES               
!                              BETA(1-4) FOR THE DIFFERENCE             
!                              EQUATIONS.                               
!                              SET UP COUNTERS                          
      RBAR = GRAPH(8)/GRAPH(9) 
      II = III 
      J = JJJ 
      MN = MNI 
      GJ = 25051.6 
      KT = KONVAR(2) 
      TDN = 2.*DELETA(1) 
      IF(KT-1) 6, 5, 6 
    5 GRMU = GRAPH(8) *  32.2 / GRAPH(10) / GRAPH(9) 
      ALPH = .33333333 
    6 CTIME = CTIME 
      JJ = II-KONVAR(1) 
      TE = T(II) 
      RHOGP = P(II) 
      WDS = WDOTP(MN) 
!****         ETVGN HAS NOT BEEN DEFINED AS YET, SO THE VALUE IS ASSUMED
!****         TO BE ZERO                                                
      ETVGN = SCDTB(18) 
      CALL TABUP(WARM,PROPT,T(II),PROP,INDTL(4),-3) 
!****         THE ARRANGEMENT OF STORING HG, EC, EP IN PROPT(1-3)       
!****         ALLOWS ONLY ONE LAYER TO BE USED)  LATER THIS WILL HAVE   
!****         TO BE CHANGED.                                            
      HG = PROP(1) 
      EC = PROP(2) 
      EP(1) = PROP(3)    ! modified by RLC
      EBAR = (EP(1)-EC)/CHART(6)+EC    ! modified by RLC
      DEL1 = (T(II) - T(II-1))/DELETA(1) 
      DEL2 = (T(II+1) - T(II))/DELETA(1) 
      IF (J-1) 200,10,200 
!****              CALL FOR K(EFF) AND RHO(CP)-EFF                      
   10 CALL KEFFS(II,MN) 
      CALL RCPEFS(II,MN) 
      NPTS = PTSIN(1) + 2. 
   25 CALL DFINTC(DENS,II,MN,J,TPTFMR,SD) 
      CALL DFINTC(P,II,MN,J,TPTFMI,RHOGNN) 
      CALL DFINTC(P,JJ,MN,J,TPTFMJ,SD) 
      CALL DFINTC(T,II,MN,J,TNFI,SDT) 
      CALL DFINTC(T,JJ,MN,J,TNFJ,SD) 
      RHON = TPTFMR/TDN 
      RHOGN = TPTFMI/TDN 
      TN = TNFI/TDN 
      TNN = SDT/DELETA(1)**2 
      TNT = (TNFI-TNFJ)/TDN/DELTO(KT) 
      RHOGNT = (TPTFMI-TPTFMJ)/TDN/DELTO(KT) 
      RHOGNN = RHOGNN/DELETA(1)**2 
   55 RHOGT = (P(II)-P(JJ))/DELTO(KT) 
      TT = (T(II)-T(JJ))/DELTO(KT) 
      RMBART = GRAPH(11)/GRAPH(9) 
      RMBRTT = GRAPH(12)/GRAPH(9) 
      RMBRSQ = RMBART ** 2 
!****         COMPUTE F TERMS                                           
      F1 = -RMBART + ALPH/TE 
      F1T = -RMBRTT + RMBRSQ - ALPH/TE**2 
      F2 = GRMU * TE ** ALPH 
      F1F2 = F1 * F2 
      RKN1 = RHOGN * TABLD(15) 
      RKN2 = TABLD(16) * RHOGP 
      RHOKN = RKN1 + RKN2 
      RKNN = RHOGNN * TABLD(15) + 2. * RHOGN * TABLD(16) + TABLD(19) *  &
     & RHOGP                                                            
      ETATX =( ETASBT(MN+1) - ETASBT(MN-1))/2./DELETA(1) 
      RHOKT = T ABLD(15) * RHOGT + TABLD(18) * RHOGP 
      RKNT = TABLD(17)*RHOGP + TABLD(16)*RHOGT + TABLD(18)*RHOGN +      &
     & TABLD(15)*RHOGNT                                                 
!****         COMPUTE THE SET OF A TERMS.                               
      A1 = F1F2 * TABLD(15) * RHOGP 
      A2 = F2 * ETASBX(MN) * RHOKN 
      A1N = A1 * ETASBX(MN) 
      A1NTN = A1N*TN 
      A1TRM = A1NTN + 2.*A2 
      A3 = F1F2 * RHOKT + F1*A2*ETASBT(MN)/ETASBX(MN) 
      A4 = F2*(RKNN*ETASBX(MN)*ETASBT(MN) + RHOKN*ETATX +               &
     & ETASBX(MN)*RKNT)                                                 
      A5 = F2*(RKNN*ETASBX(MN)**2 + RHOKN*ETAXX(MN)) 
!****         COMPUTE BETA-S.                                           
      BETA1 = (SCDTB(1) + RHOGP/GJ*A1*A2**2)/SCDTB(4) 
      BETA2 = SCDTB(2)/SCDTB(4) 
      BETA3 = SCDTB(3)*RHON*ETASBX(MN)/SCDTB(4) 
!****         COMPUTE THE TERMS FOR BETA 4.                             
      B41A = -A2**2/2./GJ - HG-EBAR 
      B41B = RHOGP * A2 * (A2*A5-A4)/GJ 
      B41C = RBAR/778. * TE * (RHOGT + ETASBT(MN) * RHOGN ) 
      B41 = WDS * B41A + B41B + B41C 
      B42 = GRAPH(4)*RHOGP*ETASBX(MN)*(A1NTN + A2) 
      IF(B42) 70, 70, 75 
   70 TETA = DEL1 
      GO TO 80 
   75 TETA = DEL2 
   80 B42 = B42*TETA 
      B43 = -WDS*  A1N* A1TRM/2. 
      IF(B43) 82, 82, 85 
   82 TETA = DEL1 
      GO TO 87 
   85 TETA = DEL2 
   87 B43 = B43 * TETA / GJ 
      B44 = RHOGP * A1N *  A5 * A1TRM 
      IF(B44) 90, 90, 92 
   90 TETA = DEL1 
      GO TO 95 
   92 TETA = DEL2 
   95 B44 = B44 * TETA / GJ 
      B451 = RHOGP *  A2**2 * ETASBX(MN) 
      B452 = 2. * A2 * F1 
      B453 = 4. * F1 + F1T + F1 ** 2 
      B454 = A1 * TN * ETASBX(MN) / F1 
      B45 = B451 * (B452 + B453 * B454) 
      IF(B45) 97, 97, 99 
   97 TETA = DEL1 
      GO TO 100 
   99 TETA = DEL2 
  100 B45 = B45 * TETA / GJ 
      RA1N = RHOGP * ETASBX(MN) * A1 ** 2 
      B461 = RA1N * ETASBX(MN) ** 2 
      B462 = A1TRM-ETASBT(MN)/ ETASBX(MN) 
      B46 = B461 * B462 * TNN 
      IF(B46) 102, 102, 105 
  102 TETA = DEL1 
      GO TO 107 
  105 TETA = DEL2 
  107 B46 = B46 * TETA / GJ 
      B471 = RA1N * TN 
      B472 = A1TRM * ETAXX(MN) - ETATX * ETASBX(MN) 
      B47 = B471 * B472 
      IF(B47) 110, 110, 112 
  110 TETA = DEL1 
      GO TO 115 
  112 TETA = DEL2 
  115 B47 = B47 * TETA / GJ 
      B481 = -RHOGP * A3 * ETASBX(MN) 
      B48 = B481 * (A1NTN + A2) 
      IF(B48) 117, 117, 119 
  117 TETA = DEL1 
      GO TO 120 
  119 TETA = DEL2 
  120 B48 = B48 * TETA / GJ 
      B49 = -RHOGP * A1 * A4 * ETASBX(MN) 
      IF(B49) 122, 122, 124 
  122 TETA = DEL1 
      GO TO 126 
  124 TETA = DEL2 
  126 B49 = B49 * TETA / GJ 
      B4101 = B461 * TN ** 2 / F1 
      F1F1T = 2. * F1 ** 2 + F1T 
      B4102 = 2. * A2 * F1F1T 
      B4103 = A1 * (F1**2 + F1T) * ETASBX(MN) * TN 
      B410 = B4101 * (B4102 + B4103) 
      IF(B410) 130, 130, 132 
  130 TETA = DEL1 
      GO TO 135 
  132 TETA = DEL2 
  135 B410 = B410 * TETA  / GJ 
      B4111 = -RHOGP * A1 * A2 / F1 
      B4112 = F1F1T * ETASBX(MN) * ETASBT(MN) * TN 
      B4113 = F1 * ETATX * ETASBX(MN) 
      B4114 = - B471*ETASBX(MN) * ETASBT(MN) 
      B411 = B4111 * (B4112 + B4113) + B4114 
      IF(B411) 137, 137, 139 
  137 TETA = DEL1 
      GO TO 140 
  139 TETA = DEL2 
  140 B411 = B411 * TETA  / GJ 
      B412 = -RHOGP * A1 * A2 * ETASBX(MN) * ET ASBT(MN) * TNN 
      B412 = B412/GJ 
      B4131 = B4111 * ETASBX(MN) 
      B4132 = B4131 * TN * TT 
      B4133 = B4131 * F1 
      B4134 = B4133 * TNT 
      B4135 = - RHOGP * A1N**2 * TN * TNT 
      B4136 = -RA1N * (2.*F1**2 + F1T) * TN * ETASBT(MN)/F1 
      IF(B4136) 150,150,152 
  150 TETA = DEL1 
      GO TO 154 
  152 TETA = DEL2 
  154 B4136 = B4136 * TETA 
      B413 = B4132 + B4134 + B4135  + B4136 
      B413 = B413/GJ 
      B414 = ETVGN * ET ASBX(MN) 
      BETA4 = (B41 + B42 + B43+ B44 + B45 + B46 + B47 + B48 + B49 +     &
     & B410 + B411 + B412 + B413 + B414) / SCDTB(4)                     
      EX = ETASBX(MN) 
      EXX = ETAXX(MN) 
      ET = ETASBT(MN) 
      GO TO 301 
  200 DEN1 = DENSE(J) 
      NPTS = PTSIN(J) + 1. 
!                              CALCULATE BETA(1-4)                      
!                       TABLE LOOK-UP FOR CK,DK/DT,CP,                  
  201 MN1 = (INDTL(4)+1)*(J-1)*3+1 
      CALL TABUP(WARM,PROPT(MN1:),TE,PROP,INDTL(4),-3) 
      SS = PROP(1) 
      CK = PROP(2) 
      DKDTUP = PROP(3) 
      BETA1 = CHARV(1)*CK/(DEN1*SS*CHARV(3)) 
      BETA2 = BETA1*DKDTUP/CK 
  206 EX = 1./THIKN(J) 
      EXX = 0.0 
      ET = 0.0 
  301 ETAXSQ = EX**2 
      BETA(1,MN) = BETA1*ETAXSQ 
      BETA(2,MN) = BETA2*ETAXSQ 
      BETA(3,MN) = BETA3*EX-ET+BETA1*EXX 
      BETA(4,MN) = BETA4 
      TIVAR(10) = AMAX1(ABS(BETA(2,MN)*DENS(II)),ABS(BETA(4,MN)),       &
     &  TIVAR(10))                                                      
!                              TEST QUANTITIES TO CHECK                 
!                              IF ITERATION IS NECESSARY                
!                              IN INTPTS                                
      IF (CHART(1)) 375,304,375 
  304 IF (BETA2) 375,305,375 
  305 IF (ABS(WDOTP(MN))) 310,315,310 
  310 IF (ABS(CG)+ABS(GVAR(1))) 375,315,375 
  315 NT1 = INDTL(4)+1 
      J1 = (J-1)*NT1*3+NT1 
      IF(ABS(PROPT(J1))) 375,380,375 
  375 ITTEMP = 1 
!                       ITTEMP = 1 MEANS ITERATION IN MIDPTS.           
      GO TO 400 
  380 ITTEMP = 0 
!                       ITTEMP = 0. MEANS NO ITERATION IN MIDPTS.       
  400 ITM = ITTEMP 
  501 RETURN 
END Subroutine Pbeta
                                          
!$IBFTC PEMDG                                                           
      SUBROUTINE PEMDG(KI) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
      K2 = KI 
    1 K1 = KONVAR(6)*KONVAR(1) 
      NPTS = PTSIN(1)+1. 
      EMDG(NPTS+1) = 0. 
      IF (K2) 3,11,3 
    3 NPTS = 1 
      K = K2-K1 
      GO TO 21 
   11 K = 1 
   21 DO 107 I1 = 1,NPTS 
      N = (NPTS-I1+1)*K 
      I = N+K1 
      IF (N) 45,45,50 
   45 N = I 
   50 IF(K2) 51, 53, 51 
   51 IF(K-1) 54, 53, 54 
   53 CALL DEE(I,N) 
   54 IF( N-1) 65,55,65 
   55 RETA = (-3.*P(I) + 4.*P(I+1) - P(I+2))/2./DELETA(1) 
      GO TO 70 
   65 RETA = (P(I+1)-P(I-1))/2./DELETA(1) 
   70 VG = -(TABLD(2)*P(I) + 3.*TABLD(1)*ETASBX(N)*RETA)/3. 
      EMDG(N) =-P(I)*VG 
  107 END DO 
  202 RETURN 
END Subroutine Pemdg
                                          
!$IBFTC PHI                                                             
      SUBROUTINE PHI(III,MNI) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
      II = III 
      MN = MNI 
      CTIME = CTIME 
      JJ = II - KONVAR(1) 
      KT = KONVAR(2) 
      J = 1 
      CALL DFINTC(T,II,MN,J,TNP,TNNP) 
      CALL DFINTC(T,JJ,MN,J,TNPJ,TNNPJ) 
      TN = TNP/2./DELETA(1) 
      TNT = (TNP-TNPJ)/2./DELETA(1)/DELTO(KT) 
      TNN = TNNP/DELETA(1)**2 
      TT = (T(II)-T(JJ))/DELTO(KT) 
      C = GRAPH(8)*32.2/GRAPH(10)/GRAPH(9) 
      PHI1 =  C  * T(II) ** .33333333 
      PHI1N = PHI1* TN/3./T(II) 
      PHI1T = PHI1* TT/3./T(II) 
      PHI1NN = PHI1 * (3. * T(II)* TNT - 2. * TN * TT)/9./T(II)**2 
      PHI1NT = PHI1 * (3. *T(II) * TNN - 2. * TN ** 2)/9./T(II) ** 2 
      SCDTB(6) = PHI1 
      SCDTB(7) = PHI1N 
      SCDTB(8) = PHI1NN 
      SCDTB(9) = PHI1NT 
      SCDTB(10) = PHI1T 
  300 RETURN 
END Subroutine Phi
                                          
!$IBFTC PRESSM                                                          
      SUBROUTINE PRESSM(III,JJJ,N,KKK) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
                                                                        
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
      II = III 
      JJ = JJJ 
      MN = II - KONVAR(6)*KONVAR(1) 
      K = KKK 
      IF(MN) 10, 15, 15 
   10 MN = 0 
   15 CALL DPRESS(II,MN,K) 
      RETURN 
END Subroutine Pressm
                                          
                                          
!$IBFTC RCPEFS                                                          
      SUBROUTINE RCPEFS(III,MNI) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
      II = III 
!****              THIS ROUTINE COMPUTES RHO-CP-EFF                     
!****                AND STORES IT IN SCDTB(4).                         
      MN = MNI 
      CV = (GRAPH(5)-GRAPH(6))/CHART(6) 
      ALPH = .33333333 
      RBAR = GRAPH(8)/GRAPH(9) 
      F1 = -GRAPH(11)/GRAPH(9) + ALPH/T(II) 
      RHOKN = TABLD(16) * P(II) + TABLD(15) *(P(II+1) - P(II-1))/2. /   &
     &DELETA(1)                                                         
      A2 = 32.2 * RBAR / GRAPH(10) * T(II) ** ALPH * ETASBX(MN) * RHOKN 
      TRM1 = P(II)*GRAPH(4) 
      TRM2 = DENS(II)*(CV + GRAPH(6)) 
      TRM3 = -CV*TABLD(14) 
      TRM4 = P(II)/778.*((A2**2*F1)/32.2-RBAR*(.66666667+T(II) * F1)) 
      SCDTB(4) = TRM1 + TRM2 + TRM3 + TRM4 
      IF(SCDTB(4)) 198, 199, 199 
  198 CALL ErrorMessage('R-CP-EFF IS -')
  199 RETURN 
END Subroutine RCPEFS
                                          
!$IBFTC RHOFFC                                                          
      SUBROUTINE RHOFFC(N,M,IM) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
                                                                        
      IMP = IM 
      KN = N 
      KM = M 
!****              KN IS THE SAME AS II.                                
!****              KM IS THE SAME AS MN.                                
!****              COMPUTE CONTROL FOR COMPUTED GO TO BY ADDING 2       
!****                   TO INPUT NUMBER IMP.                            
      KIMP = IMP + 2 
      GO TO (100, 200, 300), KIMP 
!****              FRONT FACE GAS DENSITY CALLED FROM STARTT.           
  100 RCFP = CHART(8)*(1.-CHART(5)) 
      RCP = (RCFP+(CHART(6)-1.)*DENS(KN))/CHART(6) 
      RPP = DENS(KN) - RCP 
      EMM = 1.-RCP/CHART(8)-RPP/CHART(13) 
      CALL TabUpSingle(HORA,TTABL,CTIME,TIVAR(4),INDTL(5),4) 
      P(KN) = EMM*TIVAR(4)*GRAPH(9)/GRAPH(8)/T(KN) 
      GO TO 998 
!****              FRONT FACE DENSITY OF GAS CALLED FROM PTX1.          
  200 GO TO 100 
!****              BACK FACE RHO-G CALLED FROM INTFC.                   
  300 P(KN) = (4.*P(KN-1)-P(KN-2))/3. 
  998 RETURN 
END Subroutine Rhoffc
                                          
!$IBFTC RITER                                                           
      SUBROUTINE RITER(KON) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!                       THIS ROUTINE SETS UP THE                        
!                       SPACING OF THE POINTS ALONG                     
!                       THE X AXIS.                                     
!****              KLOG = 1, COMPUTE THE R VALUES,THE X VALUES, AND     
!****                   THE ETA DERIVATIVES.                            
!****              KLOG = 2, COMPUTE THE R VALUES.                      
!****             KLOG = 3, COMPUTE THE X VALUES AND THE ETA DERIVATIVES
      KLOG = KON 
      CTIME = CTIME 
      LAYER = KONVAR(4) 
      ITIME = KONVAR(1) 
      GO TO (800, 800, 1000),KLOG 
  800 STOPT = GENRL(2) 
      IF (SPACE(1)) 120,110,120 
!                       INITIALIZE SPACING.                             
  110 ETAO = 1.0 
      AR = 0.0 
      ARP = 0.0 
      GO TO 246 
  120 IF (KONVAR(2)-1) 129,121,129 
!                       FIRST TIME THROUGH THE ROUTINE.                 
!                       BLOCK FOR INITIALIZATION FOR                    
!                       FIRST TIME STEP.                                
  121 IF(SPACE(3)) 1121, 124, 1121 
 1121 ETAMAX = AMAX1(SPACE(2),SPACE(3),SPACE(4)) 
      DO 122 I=2,4 
      IF(ETAMAX-SPACE(I)) 122,123,122 
  122 END DO 
  123 XI=I-1 
      XJ1=7.-8.*XI+2.*XI**2 
      XJ2=5.-5.*XI+XI**2 
      XJ3=22.-27.*XI+7.*XI**2 
      J4=XJ2+2. 
      J5=4.-XJ2-XI 
      XJ5=J5 
      XK1=ALOG(ETAMAX/SPACE(J4+1)) 
      XK2=ALOG(ETAMAX/SPACE(J5+1)) 
      AA=ETAMAX 
      IF(XK1*XK2) 125,126,127 
  125 CALL ErrorMessage('PRODUCT K1,K2 IS NEG.')
  126 FPRIME = 0. 
      XK6 = 1. 
      XK4=1. 
      GO TO 129 
  124 ETAO = SPACE(2) 
      BB = 0. 
      EFE = 1.0 
      GO TO 131 
  127 XK3=SQRT(XK2/XK1) 
      BB=4.*XK1 
      TAUBAR=1.-XJ2*SQRT(1.-XJ1*XK3) 
      XK4=(-XJ1*TAUBAR+XJ3/2.)/((XJ5-2.)*TAUBAR+XJ1)*STOPT/(STOPT-SPACE(&
     &5))                                                               
      FPRIME=(1.+XJ1)*TAUBAR/2.+(1.-XJ1)+(XJ2-1.)*XK4/2. 
      XK6 = (1.-XJ2)/2.*XK4*SPACE(5) 
!                       END OF INITIALIZATION BLOCK.                    
!                       COMPUTE F,TAU,AND ETAO(1).                      
  129 IF(SPACE(3)) 130, 246, 130 
  130 EFE = FPRIME*CTIME+XK6 
      TAU = EFE/((1.-XK4)*CTIME+XK4*SPACE(5)) 
      ETAO = AA*EXP(-BB*(TAU*(1.-TAU/2.))**2) 
  131 IF (ETAO) 138,144,144 
  138 CALL ErrorMessage('ETAO IS NEGATIVE')
!****              END OF ETAO BLOCK.                                   
!****              BLOCK FOR ITERATION ON AR.                           
  144 IF (ETAO-1.) 155,180,146 
  146 AR = ETAO 
  148 RI = ETAO*(1.-EXP(-AR)) 
      IF (ABS(AR-RI)-.0005) 225,225,150 
  150 AR = RI 
      GO TO 148 
  155 C1 = -ETAO*(1.+(2.*(2.-EXP(ETAO))) /((EXP(ETAO) -1.)**2)) 
      IF ((-2.*ETAO)-C1) 157,159,159 
  157 AR = -2.*ETAO 
      GO TO 160 
  159 AR = C1 
  160 RI = -ALOG(1.-AR/ETAO) 
      IF (ABS(AR-RI)-.0005) 225,225,161 
  161 AR = RI 
      GO TO 160 
  180 AR = 0. 
!****              END OF AR BLOCK.                                     
!****              BLOCK TO COMPUTE ARP.                                
  225 IF (AR) 232,227,232 
  227 FSUBR = 2. 
      GO TO 240 
  232 FSUBR = AR/(1.0-ETAO+AR) 
      IF (EFE) 240,235,240 
  235 ARP = 0.0 
      GO TO 246 
  240 TRMATA = -BB*TAU*(TAU-1.)*(TAU-2.)*TAU*(FPRIME-TAU*(1.-XK4))/EFE 
  245 ARP = FSUBR*TRMATA 
!****              NEW TRANSFORMATION                                   
  246 SPACVR(1) = AR 
      SPACVR(2) = ARP 
      SPACVR(3) = ETAO 
!****              END OF ARP BLOCK.                                    
      GO TO (1000, 600, 1000),KLOG 
 1000 DO 599 J = 1,LAYER 
      NPTS = PTSIN(J)+2. 
      IF(J-1) 247, 247, 248 
  247 N1 = 1 
      GO TO 249 
  248 N1 = 2 
  249 DO 598 I = N1,NPTS 
      K = 0 
      MN = JLAYER(J)+I 
      EI = I-1 
      ETAI = EI*DELETA(J) 
      IF( I - NPTS) 1248, 1247, 1248 
 1247 SQUIGL = 1.0 
      IF(J-1) 265, 265, 560 
 1248 IF(J-1) 600, 1249, 560 
 1249 IF(AR) 260, 250, 260 
  250 SQUIGL = ETAI 
      GO TO 265 
  260 SQUIGL = ALOG((1.-ETAI)+ETAI*EXP(-AR)) 
      SQUIGL = SQUIGL/(-AR) 
!****         IF ( C .GT. 8 ) CALL ErrorMessage MESSAGE AT 390                 
  265 IF(SPACE(6) - 8.) 267, 390, 390 
  267 IF(SPACE(6)) 280, 270, 280 
  270 XBAR = SQUIGL 
      GO TO 400 
  280 SQSTR = .5*EXP(-SPACE(6)/4.) 
      IF (SQUIGL-SQSTR) 290,285,300 
  285 XBAR = .5 
      GO TO 400 
  290 X = SQUIGL 
  293 X1 = SQUIGL*EXP(SPACE(6)*X*(1.-X)) 
      IF(ABS(X1-X)-.0005*X) 297, 297, 295 
  295 X = X1 
      K = K + 1 
      IF(K-100) 293,293,296 
  296 CALL ErrorMessage('NO CONVERGENCE IN RITER EFN 293')
  297 XBAR = X1 
      GO TO 400 
  300 IF(SPACE(6) -2.) 320, 320, 301 
  301 ARG1 = SQUIGL * EXP( SPACE(6) / 4. ) 
      X = AMIN1(ARG1,1.) 
  303 FX = SQUIGL*EXP(SPACE(6)*X*(1.-X)) 
      FP = SPACE(6)*(1.-2.*X)*FX 
      X1 = (FX-X*FP)/(1.-FP) 
      IF(ABS(X1-X)-.0005*X) 310, 310, 305 
  305 X = X1 
      K = K + 1 
      IF( K-100) 303, 303, 306 
  306 CALL ErrorMessage('NO CONVERGENCE IN RITER EFN 303')
  310 XBAR = X1 
      GO TO 400 
  320 X  = (1. + SQRT(2./SPACE(6)))/2. 
  322 FX = SQUIGL * EXP(SPACE(6) * X * (1.-X)) 
      X1 = (FX + 4. * X) / 5. 
      IF(ABS(X1-X) - .0005 * X) 310, 310, 325 
  325 X = X1 
      K = K + 1 
      IF(K - 100) 322, 322, 336 
  336 CALL ErrorMessage('NO CONVERGENCE IN RITER EFN 336')
  390 CALL ErrorMessage('SPACE(6) IS GREATER THAN 8')
  400 XLOC(MN) = THIKN(J)*XBAR+HTCONV(7) 
      IF (MN-1) 410,405,410 
  405 ETASBX(MN) = ETAO/THIKN(J) 
      ETAXO = ETAO/THIKN(J) 
  410 TSQXB = 1. + SPACE(6)*XBAR*(2.*XBAR-1.) 
      IF(XBAR) 420, 415, 420 
  415 SQXB = TSQXB 
      GO TO 425 
  420 SQXB = SQUIGL*TSQXB/XBAR 
  425 ETASBX(MN) = ETAXO*((1.-ETAI)+ETAI*EXP(-AR))*SQXB 
!****              SECTION TO COMPUTE 2ND PARTIAL OF ETA                
!****                   WITH RESPECT TO X.                              
      SXX = SPACE(6)*((6.*XBAR-2.)+SPACE(6)*XBAR*(2.*XBAR-1.)**2) 
      SXXDSX = SXX/TSQXB 
      ETAXX(MN) = ETASBX(MN)*(-AR*SQXB+SXXDSX)/THIKN(J) 
      GO TO 598 
  560 NPTS1 = JLAYER(J) + 1 
      XLOC(MN) = THIKN(J)*ETAI + XLOC(NPTS1) 
      ETASBX(MN) = 1./THIKN(J) 
      ETAXX(MN) = 0.0 
  598 END DO 
  599 END DO 
  600 DO 640 J = 1,LAYER 
      J1 = J+1 
      NPTS = PTSIN(J) + 2. 
      IF(J-1)  612,612,615 
  612 N1 = 1 
      GO TO 617 
  615 N1 = 2 
  617 DO 630 I= N1,NPTS 
      MN = JLAYER(J) +I 
      CALL COORD(MN,J1) 
  630 END DO 
  640 END DO 
  652 RETURN 
END Subroutine Riter
                                          
                                          
!$IBFTC TNEST                                                           
      SUBROUTINE TNEST(TNEW,OUT,EPSS1) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATIONS IN THIS ROUTINE ARE DEFINED AND DISCUSSED   
!****         IN APPENDIX C                                             
      DIMENSION ITCNT(10), J1(10),L(10),T1(10), T2(10), T3(10), TAVG(10) 
      DIMENSION OVALIN(10), OVALOT(10) 
!                  TNEW IS THE NEW GUESS AND OUTPUTOF THE ROUTINE.      
!                  TNINT IS THE NEW GUESS (INTERNAL TO THE ROUTINE).    
!                  TGUESS IS THE LATEST GUESS.                          
      EPSS = EPSS1 
!                  FIRST TIME THROUGH,  OUT1 IS INITIAL GUESS.          
      OUTI = OUT 
      TGUESS = TNEW 
!                  J1(N) = 1       NORMAL OPERATION.                    
!                  J1(N) = 2       ITERATIVE OPERATION TO GET           
!                       INFORMATION,  BECAUSE OF NON-CONVERGENCE.       
      J1N = J1(N) 
      IF(EPSS) 1,108,2 
!                       IF EPSS =                                       
!                       -  IT IS THE FIRST TIME THROUGH.                
!                       0. IT IS INITIALIZATION.                        
!                       + 2 OR 3 POINTS ARE IN CORE                     
!                       SO INTERPOLATION AND NESTING CAN BE DONE.       
    1 N = N+1 
      IF(N-10) 600, 600, 106 
  600 OVALIN(N) = OUTI 
      OVALOT(N) = TNEW 
      L(N) = 1 
      OUTI = 2.01 
      ITCNT(N) = 0 
      TNINT = TGUESS 
      T1(N) = OVALIN(N) 
      T2(N) = TGUESS 
      IF(ABS(T1(N)-T2(N))-ABS(EPSS*T2(N))) 105, 105, 110 
!     SECOND TIME THROUGH                                               
    2 N = N 
      LL = L(N) 
!                  L(N) = 1        FUNCTION IS MONOTONIC.               
!                  L(N) = 2        FUNCTION IS NON-MONOTONIC.           
      GO TO (3,42),LL 
    3 T3(N) = TGUESS 
!                       TEST FOR CONVERGENCE                            
      IF(ABS(T2(N)-T3(N))-ABS(EPSS*T2(N))) 105, 105, 6 
    6 TMAX =AMAX1(T1(N),T2(N),T3(N)) 
      TMIN =AMIN1(T1(N),T2(N),T3(N)) 
!                  TMAX IS THE MAXIMUM OF THE 3 GUESSES.                
!                  TMIN IS THE MINIMUM OF THE 3 GUESSES.                
!                       TEST IF MONOTONE                                
      IF(T1(N)- TMAX)  30, 27, 30 
   27 IF(T3(N)- TMIN)  37,35,37 
   30 IF(T1(N)- TMIN ) 37,32,37 
   32 IF(T3(N)- TMAX ) 37,35,37 
!                       NUMBERS ARE MONOTONIC                           
!                  T1(N), T2(N), T3(N)  ARE CONSECUTIVE GUESSES         
!                       IN MONOTONIC ITERATION.                         
   35 T1(N) = T2(N) 
      T2(N) = T3(N) 
      TNINT = T2(N) 
!                       NO CONVERGENCE, SO SET ITERATION INDICATOR      
      OUTI = 2.01 
      GO TO 88 
!                       NON-MONOTONIC OR NESTING                        
   37 T1(N) = TMIN 
      T2(N) = TMAX 
!                  TAVG(N) IS THE AVERAGE OF T1(N) AND T2(N).           
!                       IN NESTING ITERATION,  NEXT GUESS.              
   38 TAVG(N)  = (T1(N) + T2(N))/2.0 
      TNINT = TAVG(N) 
      L(N) = 2 
!                       NO CONVERGENCE, SO SET ITERATION INDICATOR      
      OUTI = 2.01 
      GO TO 88 
!                       NON-MONTONIC OR NESTING                         
   42 IF(ABS(TAVG(N)-TGUESS)-ABS(EPSS*TAVG(N))) 105,105,43 
   43 IF(TAVG(N)-TGUESS) 44,44,46 
!                  T1(N), T2(N)  ARE THE BOUNDS IN THE                  
!                       NON-MONOTONIC ITERATION.                        
   44 T1(N) = TAVG(N) 
      GO TO 38 
   46 T2(N) = TAVG(N) 
      GO TO 38 
!                       ITERATION TEST AND FINIS                        
   88 ITCNT(N) = ITCNT(N) + 1 
   95 IF(ITCNT(N)-50) 110,110,99 
   99 J1N = J1N 
      GO TO (100,200),J1N 
  100 J1(N) = 2 
      WRITE(6,500) 
  112 OUTI = 3.01 
      TNEW = OVALOT(N) 
      T1(N) = OVALIN(N) 
      T2(N) = OVALOT(N) 
      WRITE(6,506) T1(N),T2(N) 
      L(N) = 1 
      GO TO 115 
  200 WRITE(6,501) N, L(N), TGUESS, TNINT, T1(N), T2(N), T3(N), TAVG(N),&
     &TMAX, TMIN                                                        
      IF(ITCNT(N) - 66) 201, 205, 205 
  201 OUTI = 3.01 
      GO TO 110 
  205 CALL ErrorMessage('15 RETRIES COMPLETED. PDUMP FOLLOWS')
!                       CONVERGED, SO SET INDICATOR                     
  105 OUTI = 1.01 
      TNINT = TGUESS 
      TAVG(N) = 0. 
      J1(N) = 1 
      N = N-1 
      IF(N) 106, 110, 110 
  106 WRITE(6,505) N 
      CALL ErrorMessage('INIT. OF EPSILON IS INCORRECT')
  108 DO 109 I=1,90 
  109 ITCNT(I) = 0 
      EPSS = 0.0 
      J1N = 0 
      LL = 0 
      OUTI = 1.0 
      N = 0 
      TMAX = 0.0 
      TMIN = 0.0 
      DO 1110 J = 1,10 
 1110 J1(J) = 1 
  110 TNEW = TNINT 
      OUT = OUTI 
!                  OUT = 1    MEANS CONVERGENCE.                        
!                  OUT = 2    MEANS NO CONVERGENCE  PASS ONE            
!                  OUT = 3    MEANS NO CONVERGENCE  PASS TWO.           
  115 RETURN 
  500 FORMAT(////' NO CONVERGENCE IN TNEST, INFORMATION FOLLOWS.')
  501 FORMAT(/,3X,'N = ',I3,3X,'L = ',I3,3X,'GUESS = ',E18.8,3X, &
      'NEW GUESS = ',E18.8// &
      3X,'T1(N) = ',E18.8,3X,'T2(N) = ',E18.8,3X,'T3(N) = ',E18.8// &
      3X,'TAVG(N) = ',E18.8,3X,'TMAX = ',E18.8,3X,'TMIN = ',E18.8// &
      1X,'*******' )                                                     
  505 FORMAT(/,3X,'N = ',I3) 
  506 FORMAT(/,'ORIG. GUESS = ',E20.8,4X,'ORIG. ANS. =', E20.8) 
END Subroutine Tnest
                                          
!$IBFTC WDP                                                             
      SUBROUTINE WDP(II,L) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION ZORE(4) 
!****         THE EQUATION FOR THE RATE OF CHANGE OF DENSITY WITH       
!****         RESPECT TO TIME, W-DOT-P, IS DEFINED IN SECTION II.B.2    
!****       THIS SUBROUTINE COMPUTES THE W-DOT-P FOR EACH POINT.        
      I = II 
      TE = T(I) 
!****         SECTION FOR SINGLE REACTION.                              
      D = DENS(I) 
      TRM = DENSE(1)*((D-CHART(2))/DENSE(1))**CHART(3) 
      IF (GENRL(12)) 9, 1, 9 
    1 NDENS = INDTL(1) 
!                              TABLE LOOK-UP FROM Z OR E TABLE.         
      CALL TABUP(TDENS,ZORET,D,ZORE,NDENS,-4) 
      EX1 = -ZORE(2)/HTCONV(1)/TE 
      WDOTP(L) = -ZORE(1)*EXP(EX1)*TRM 
!                              TABLE LOOK-UP FROM Z OR E TABLE.         
    5 RETURN 
!****         SECTION FOR MULTIPLE REACTION.                            
    9 NZ = GENRL(12) 
      DRDT = 0.0 
      DO 10 N = 1,NZ 
      NPNZ = 2 * N 
      EX1 = -ZORET(NPNZ)/HTCONV(1) / TE 
      N1 = 2 * (N-1) + 1 
      DRDT = DRDT + ZORET(N1) * EXP(EX1) 
   10 END DO 
      WDOTP(L) = -DRDT*TRM 
      GO TO 5 
END Subroutine Wdp
                                          
!$ORIGIN        ALPHA                                                   
!$IBFTC CLPOLY                                                          
      FUNCTION CLPOLY(ARG,COFARG,NARG) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!     POLYNOMIAL EVALUATION                                             
      DIMENSION COFARG(1) 
    1 CLPOLY=COFARG(NARG+1) 
      IF(NARG.EQ.0) GO TO 5001 
      DO 91 J1=1,NARG 
      N1=NARG-J1+1 
      CLPOLY=CLPOLY*ARG+COFARG(N1) 
   91 END DO 
 5001 RETURN 
END Function ClPoly
                                          
!$IBFTC DIFTAB                                                          
      SUBROUTINE DIFTAB(Y,X,DYDX,NX,NL) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
      DIMENSION X(1),Y(1),DYDX(1) 
!     NX = INDEPENDENT TABLE LENGTH,NL = NUMBER OF LAYERS.              
!     DETERMINE STARTING POINT OF Y TABLE INTERVAL JK1.                 
      NX = NX 
      NL = NL 
      NL1 = NL-1 
      NX1 = NX+1 
      NX2 = NX-2 
      DO 100 JL = 1,NL 
    2 JK = NX1*JL 
      JK1 = JK-NX 
    4 IF (Y(JK)) 5,89,5 
    5 IF (NX2) 50,50,10 
!     CALCULATE THE DIFFERENCES FOR THE FIRST POINT.                    
   10 IF(ABS(Y(JK1+1)-Y(JK1))-ABS(Y(JK1+1)*1.E-7)) 15, 15, 11 
   11 D1 = X(2) - X(1) 
      D2 = X(3)-X(2) 
      A = -(2.*D1+D2)/D1/(D1+D2) 
      B = (D1+D2)/D1/D2 
      C = -D1/D2/(D1+D2) 
!     CALCULATE THE DERIVATIVE OF THE FIRST POINT.                      
      DYDX(JK1)  = A*Y(JK1)+B*Y(JK1+1)+C*Y(JK1+2) 
      GO TO 16 
   15 DYDX(JK1) = 0. 
!     CALCULATE THE DIFFERENCES FOR THE MIDDLE POINTS.                  
   16 DO 20 N = 1,NX2 
      N1 = N+1 
      D1 = X(N1)-X(N) 
      D2 = X(N1+1)-X(N1) 
      A = -D2/D1/(D1+D2) 
      B = (D2-D1)/D1/D2 
      C = D1/D2/(D1+D2) 
!     CALCULATE THE DERIVATIVES FOR THE MIDDLE POINTS.                  
      JK2 = JK1+N 
   19 DYDX(JK2) = A*Y(JK2-1)+B*Y(JK2)+C*Y(JK2+1) 
   20 END DO 
!     CALCULATE THE DIFFERENCES FOR THE END POINT.                      
      JK3 = JK1+NX-1 
      IF(ABS(Y(JK3)-Y(JK3-1))-ABS( Y(JK3)*1.E-7)) 29,29,25 
   25 D1 = X(NX-1)-X(NX2) 
      D2 = X(NX)-X(NX-1) 
      A = D2/D1/(D1+D2) 
      B = -(D1+D2)/D1/D2 
      C = (D1+2.*D2)/D2/(D1+D2) 
!     CALCULATE THE DERIVATIVE OF THE LAST POINT.                       
   30 DYDX(JK3) = A*Y(JK3-2)+B*Y(JK3-1)+C*Y(JK3) 
      GO TO 60 
   29 DYDX(JK3) = 0. 
      GO TO 60 
   50 D1 = X(NX)-X(NX-1) 
      JK3 = JK1+1 
      DYDX(JK1) = (Y(JK3)-Y(JK1))/D1 
      DYDX(JK3) = DYDX(JK1) 
   60 NX3 = NX-1 
      DO 87 J = 1,NX3 
      JK5 = JK1 +J - 1 
      IF (ABS(DYDX(JK5+1)-DYDX(JK5))-ABS(DYDX(JK5+1))*1.E-7) 75,80,80 
   75 DYDX(JK) = 0. 
      GO TO 87 
   80 DYDX(JK)  = 1. 
      GO TO 100 
   87 END DO 
      GO TO 100 
   89 DO 95 K = JK1,JK 
   95 DYDX(K) = 0. 
  100 END DO 
  200 RETURN 
END Subroutine Diftab
                                          
!$IBFTC ENLARG                                                          
      SUBROUTINE ENLARG 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!************COMMON FOR VARIABLE T-P-D INPUT*******FOR BETA-4 TABLE**** 
      DIMENSION TMDENS(1),TMT(1),TMP(1),XFORT(1),XFORP(1),XFORDN(1) 
!     ABOVE IS DUMMY DIMENSION -JUST TO ALLOW SUBSCRIBTING              
      EQUIVALENCE(TMDENS(1),DENS(301)),(TMT(1),T(301)),(TMP(1),P(301)) 
      EQUIVALENCE(XFORT(1),T(151)),(XFORP(1),P(151)),(XFORDN(1),DENS(151))
      EQUIVALENCE(NT,T(449)),(NXT,T(450)),(NP,P(449)),(NXP,P(450)) 
      EQUIVALENCE(ND,DENS(449)),(NXD,DENS(450)) 
!     THE T-P-AND DENS INPUTS ARE STORED IN THE                         
!     LAST THIRD OF THEIR RESPECTIVE ARRAYS.                            
!     NT-NP-NDENS ETC ARE STORED IN LOC                                 
!     NOT USED -DUE TO THE IMPOSSIBILITY OF HAVING 150 POINTS           
!     THE X TABLE FOR T(XFORT)ETC ARE STORED IN THE                     
!     MIDDLE THIRD OF THEIR RESPECTIVE ARRAYS                           
!     FOR THE BETA 4 TABLE HAVING PSEUDO - VARIABLE DIMENSION           
      COMMON TXQTAB(200),NOFTS,NOFXS,NTXQ(50),ICON,JCON,END 
      CALL EXPAND(NWARM,PROPT,630,WARM,NPROPT,30) 
      CALL EXPAND(NWARM,GASPR,105,WARM,NGASPR,5) 
      CALL EXPAND(NWARM,HTABL,42,WARM,NHTABL,2) 
      CALL EXPAND(NWARM,CTABL,210,WARM,NCTABL,5) 
      CALL EXPAND(NHORA,TTABL,510,HORA,NTTABL,10) 
      CALL EXPAND( NTDENS,ZORET,84,TDENS,NZORET,2) 
      CALL EXPAND( NDGREE,ATETB,420,DGREE,NATETB,20) 
      INDTL(1) = NTDENS 
      INDTL(2) = NDGREE 
      INDTL(4) = NWARM 
      INDTL(5) = NHORA 
!     TEST DATA FOR VARIABLE T-P-DENS DISTRIBUTIONS.                    
!     NO. OF ENTRIES FOR ASSOCIATED TABLES MUST BE EQUAL.               
!     CHECK DENSITY INPUT                                               
      IF (GENRL(4)) 275,200,275 
  200 IF (ND.EQ.NXD) GO TO 215 
      WRITE (6,205)ND,NXD 
  205 FORMAT (/,'BAD COUNT ON VARIABLE DENSITY TABLE   ND = ',I4,5X,'ND = ',I4)
      CALL ErrorMessage('1') 
!     CHECK PRESSURE INPUT                                              
  215 IF (NP.EQ.NXP) GO TO 230 
      WRITE (6,220) NP,NXP 
  220 FORMAT (/,'BAD COUNT ON VARIABLE PRESSURE TABLE   NP = ',I4,5X,'NXP = ',I4)
      CALL ErrorMessage('2') 
  230 IF (NT.EQ.NXT) GO TO 275 
!     CHECK TEMP INPUT                                                  
      WRITE (6,240)NT,NXT 
  240 FORMAT (/,'BAD COUNT ON VARIABLE TEMP. TABLE   NT = ',I4,5X,'NXT = ',I4)
      CALL ErrorMessage('3') 
!     THIS ENDS THE TESTING FOR VARIABLE DISTRIBUTION INPUT.            
!     BEGIN HERE TO TEST INPUT FOR INTERNAL HEAT GENERATION             
!     BETA-4CORRECTION.-IF NEEDED                                       
  275 IF (NOFTS) 300,350,300 
!     CHECK NUMBER OF TEMP TABLES.                                      
  300 IF (NOFTS.LT.148) GO TO 315 
      WRITE (6,305)NOFTS 
  305 FORMAT (/,'MAXIMUM NO OF TEMP TABLES = 148 YOU HAVE ',I4)
      CALL ErrorMessage('4')
!     CHECK NUMBER OF ENTRIES PER TABLE.                                
  315 DO 310 JD = 3,NOFTS 
      IF (NOFXS.EQ.NTXQ(JD)) GO TO 310 
      WRITE (6,307)JD,NOFXS,NTXQ(JD) 

  307 FORMAT (/,'BAD COUNT ON B4 TEMP TABLE',I4,3X,'SHOULD BE= ',I4,4X,'YOU HAVE ',I4)
      CALL ErrorMessage('5') 
  310 END DO 
!     THIS FINISHES CHECKING THE B4 INPUT. PROCEED WITH THE PROGRAM.    
  350 RETURN 
END Subroutine Enlarg

!$IBFTC EXEC1
      SUBROUTINE EXEC1()

USE HeadG,ONLY: head
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!************COMMON FOR VARIABLE T-P-D INPUT*******FOR BETA-4 TABLE**** 
      DIMENSION TMDENS(1),TMT(1),TMP(1),XFORT(1),XFORP(1),XFORDN(1) 
!     ABOVE IS DUMMY DIMENSION -JUST TO ALLOW SUBSCRIBTING              
      EQUIVALENCE(TMDENS(1),DENS(301)),(TMT(1),T(301)),(TMP(1),P(301)) 
      EQUIVALENCE(XFORT(1),T(151)),(XFORP(1),P(151)),(XFORDN(1),DENS(151&
     &))                                                                
      EQUIVALENCE(NT,T(449)),(NXT,T(450)),(NP,P(449)),(NXP,P(450)) 
      EQUIVALENCE(ND,DENS(449)),(NXD,DENS(450)) 
!     THE T-P-AND DENS INPUTS ARE STORED IN THE                         
!     LAST THIRD OF THEIR RESPECTIVE ARRAYS.                            
                                                                        
!     NT-NP-NDENS ETC ARE STORED IN LOC                                 
!     NOT USED -DUE TO THE IMPOSSIBILITY OF HAVING 150 POINTS           
!     THE X TABLE FOR T(XFORT)ETC ARE STORED IN THE                     
!     MIDDLE THIRD OF THEIR RESPECTIVE ARRAYS                           
!     FOR THE BETA 4 TABLE HAVING PSEUDO - VARIABLE DIMENSION           
      COMMON TXQTAB(200),NOFTS,NOFXS,NTXQ(50),ICON,JCON,END 
      COMMON MK 
      DIMENSION TABLN(120) 
      INTEGER:: errCode

!----------------------------------------------------------------------------
!                       SET UP THE SPACING IN THE FIRST LAYER.          
!****              COMPUTE R VALUES,ETA DERIVATIVES AND X LOCATIONS.    
!                       THIS SUBROUTINE CONTROLS                        
!                       THE MAIN FLOW OF THE                            
!                       PROGRAM. IT CALLS THE                           
!                       SUBROUTINES IN THEIR CORRECT                    
!                       AND LOGICAL ORDER, AND CAN                      
!                       BE CONSIDERED AS A GENERAL                      
!                       TYPE OF FLOW CHART.                             
    2 GO TO (3,4),ICON 
    3 CALL INIT(4) 
       MK = 0 
!                       READ IN HEADING                                 
    4 READ(5,'(A)',IOSTAT=errCode) head   ! (HEAD(I), I = 2,13) 
      IF (errCode < 0) THEN
        WRITE(*,*) 'End of input data'
        STOP
      END IF
!    5 FORMAT(12A6) 
!                       READ IN INPUT VALUES BY CALLING INPUT           
      CALL INPUT(1,END,ALLNEW) 
!                  INITIALIZE TABLES, ETC.                              
!****         CHECK MK TO SEE IF A RESTART TAPE WAS USED.               
      IF (MK) 8, 10, 8 
    8 CALL INIT (5) 
      CALL TableOutput()
      JCON = 2 
      GO TO 25 
   10 CALL INIT(1) 
!                       PRINT OUT THE INPUT TABLES.                     
   11 CALL TableOutput()
   12 CALL RITER(1) 
!****                   COMPUTE THE ETA(T)-S INITIALLY.                 
   13 CALL NSUBT 
!                       COMPUTE THE STARTING TEMPERATURES,              
!                       DENSITIES, AND PRESSURES.                       
   14 CALL STARTT 
!                       COMPUTE THE DELTA T.                            
      JCON = 1 
   25 RETURN 
   97 IF(END)3,98,3 
   98 CALL INIT(3) 
      GO TO 4 
END Subroutine Exec1
                                     
!$IBFTC EXPAND                                                          
      SUBROUTINE EXPAND(MIND,XTABL,MARAY,ENDTBL,NVLUS,NMTB) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!                    THIS SUBROUTINE IS USED TO EXPAND A                
!                    DEPENDENT TABLE ACCORDING TO THE NUMBER OF VALUES  
!                    IN THE CORRESPONDING INDEPENDENT TABLE             
!                    THE SUBROUTINE FOLLOWS THE PROCEDURE OF PLACING    
!                    THE TABLES  (GOING FROM THE LAST TO THE FIRST)     
!                    AT THE END OF THE ARRAY READING FROM THE LAST      
!                    ENTRY TO THE FIRST IN THE TABLE.                   
!                    ONCE ALL THE TABLES HAVE BEEN EXPANDED             
!                    AND MOVED TO THE BACK OF THE ARRAY,THE BLOCK       
!                    IS SHIFTED FORWARD SO THAT THE TABLES ARE          
!                    IN THEIR PROPER LOCATIONS                          
!                    THE NPLCE ARRAY STORES THE LOCATION                
!                    FROM THE START OF THE DEPENDENT TABLE ARRAY        
!                    TO THE LAST ENTRY OF A PARTICULAR TABLE            
!                    IN THE ARRAY                                       
!                    NSLSH IS THE LOCATION WHICH STORES                 
!                    THE NUMBER OF TABLES READ INTO THE ARRAY           
      DIMENSION XTABL(1),NVLUS(1),NPLCE(31),ENDTBL(1) 
      IF(NVLUS(1).EQ.0) GO TO 60 
      NIND=MIND 
      IF(NIND .EQ. 0 ) GO TO 60 
      NARAY=MARAY 
      MXTB=NMTB 
      NSLSH=0 
      NPLCE(1)=0 
      NADD=0 
      DO 7 I=1,MXTB 
      NADD=NADD+NVLUS(I) 
      IF(NVLUS(I)) 8,9,1 
    1 NSLSH=NSLSH +1 
      NPLCE(NSLSH+1)=NADD 
    7 END DO 
      GO TO 9 
    8 CALL ErrorMessage('THERE IS A TABLE INPUT ERROR')
!                            TEST IF A TABLE WAS READ IN                
    9 IF(NSLSH) 8,60,10 
   10 DO 49 I = 1,NSLSH 
      KSLSH=NSLSH-I+1 
      NUMB=NVLUS(KSLSH) 
!                            TEST TO FIND THE NUMBER OF                 
!                            VALUES IN EACH TABLE                       
      IF(NUMB-2) 20,30,18 
!                    TEST IF ALL VALUES ARE IN THE TABLE                
   18 IF(NUMB-NIND) 8,40,8 
!                            THERE IS ONE VALUE SO THE                  
!                            TABLE IS CONSTANT                          
   20 JK=NARAY-(I-1)*(NIND+1) 
      JK1=NPLCE(KSLSH) 
      TEMVL=XTABL(JK1+1) 
      XTABL(JK)=0. 
      DO 25 J=1,NIND 
      KJ = JK-J 
   25 XTABL(KJ) = TEMVL 
      GO TO 49 
!                            THERE ARE TWO VALUES SO THE                
!                            TABLE MUST BE EXPANDED BY                  
!                            INTERPOLATION                              
   30 KK = NPLCE(KSLSH) 
      DEP1 = XTABL(KK+1) 
      DEP2 = XTABL(KK+2) 
      JK = NARAY-(I-1)*(NIND+1) 
      XTABL(JK) = 1. 
      DO 35 J = 1,NIND 
      KJ = JK+J-NIND-1 
   35 XTABL(KJ) = (ENDTBL(J)*(DEP2-DEP1)/(ENDTBL(NIND)-ENDTBL(1)))+((   &
     &ENDTBL(NIND)*DEP1)-(ENDTBL(1)*DEP2))/(ENDTBL(NIND)-ENDTBL(1))     
      GO TO 49 
!                            THE WHOLE TABLE IS IN AND                  
!                            MUST BE TRANSFERRED AND                    
!                            EXPANDED                                   
   40 LL = NPLCE(KSLSH+1)+1 
      JK = NARAY-(I-1)*(NIND+1) 
      XTABL(JK) = 1. 
      DO 45 J = 1,NIND 
      KJ = JK-J 
      LL1=LL-J 
   45 XTABL(KJ) = XTABL(LL1) 
   49 END DO 
!                            NOW THE TABLES MUST BE                     
!                            MOVED TO THEIR PROPER LOCATIONS            
!                            IN THE ARRAY                               
   50 KSTART = NSLSH*(NIND+1) 
      JSTART = NARAY-KSTART+1 
      KCONT = 0 
      DO 55 J = JSTART,NARAY 
      KCONT = KCONT+1 
   55 XTABL(KCONT) = XTABL(J) 
      KCONT = KCONT+1 
      IF(KCONT.GT.NARAY) GO TO 57 
      DO 56 L = KCONT,NARAY 
   56 XTABL(L) = 0. 
   57 DO 58 IK=1,NSLSH 
   58 NVLUS(IK)=0 
   60 RETURN 
END Subroutine Expand          
                                
!$IBFTC INPUT                                                           
      SUBROUTINE INPUT(K,X,Y) 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!************COMMON FOR VARIABLE T-P-D INPUT*******FOR BETA-4 TABLE**** 
      DIMENSION TMDENS(150),TMT(150),TMP(150),XFORT(300),XFORP(300),    &
     &XFORDN(300)                                                       
      EQUIVALENCE(TMDENS(1),DENS(301)),(TMT(1),T(301)),(TMP(1),P(301)) 
      EQUIVALENCE(XFORT(1),T(151)),(XFORP(1),P(151)),(XFORDN(1),DENS(151&
     &))                                                                
      EQUIVALENCE(NT,T(449)),(NXT,T(450)),(NP,P(449)),(NXP,P(450)) 
      EQUIVALENCE(ND,DENS(449)),(NXD,DENS(450)) 
!     THE T-P-AND DENS INPUTS ARE STORED IN THE                         
!     LAST THIRD OF THEIR RESPECTIVE ARRAYS.                            
!     NT-NP-NDENS ETC ARE STORED IN LOC                                 
!     NOT USED -DUE TO THE IMPOSSIBILITY OF HAVING 150 POINTS           
!     THE X TABLE FOR T(XFORT)ETC ARE STORED IN THE                     
!     MIDDLE THIRD OF THEIR RESPECTIVE ARRAYS                           
!     FOR THE BETA 4 TABLE HAVING PSEUDO - VARIABLE DIMENSION           
      COMMON TXQTAB(200),NOFTS,NOFXS,NTXQ(50),ICON,JCON,END 
!!!      DIMENSION X(1),Y(1) 
  REAL:: x,y      
  INTEGER:: errCode
      DIMENSION TABLD(20),SCDTB(20),SCITB(20),GRAPH(20) 
      EQUIVALENCE (PVAR,GRAPH) , (PVAR(21),SCITB) , (TRAJV,SCDTB) ,     &
     &             (TRAJV(21),TABLD)                                    
!****         NAMELIST STATEMENTS FOR GE-STC.       

      NAMELIST /REKAP/ AGPTB,ALLNEW,AMPLCT,ATETB,BACKFC,CHART,CTABL,     &
     & DENSE,DGREE,DTIME,END,EPTAB,GASPR,GENRL,HORA,HTABL,NATETB,        &
     & NCTABL,ND,NDGREE,NGASPR,NHORA,NHTABL,NP,NPHI,NPROPT,NT,NTDENS,    &
     & NTTABL,NTXQ,NU,NWARM,NXD,NXP,NXT,NZORET,OUTTB,PHII,PROPT,PTSIN,   &
     & PVAR,SPACE,STVAR,TDENS,THIKN,TMDENS,TMELT,TMP,TMT,TRAJT,TTABL,    &
     & TXQTAB,WARM,XFORDN,XFORP,XFORT,ZORET

                    
!      NAMELIST /REKAP /  AGPTB,    DENSE,    DTIME,    EPTAB,    GENRL, &
!     &                   OUTTB,    PTSIN,    SPACE,    THIKN,           &
!     &                   CTABL,    HTABL,    NCTABL,   NHTABL,   NPROPT,&
!     &                   NWARM,    PROPT,    WARM,                      &
!     &                   CHART,    GASPR,    NGASPR,   NTDENS,   NZORET,&
!     &                   HORA ,    NHORA,    NTTABL,   TTABL,           &
!     &                   TDENS,    ZORET,    TMELT ,   TRAJT ,          &
!     &                   ATETB,    DGREE ,   NDGREE,   NATETB,  NPHI   ,&
!     &                   NU   ,    PHII  ,   STVAR ,   PVAR,            &
!     &                   TXQTAB,  NTXQ,  TMT,  NT,  XFORT,  NXT,        &
!     &                   TMP,  NP,  XFORP,  NXP,  TMDENS,  XFORDN,  NXD,&
!     &                   ALLNEW,  END,  AMPLCT,  BACKFC,  ND            
!****         ADD TO NAMELIST.                                          
!                                                                       
      READ (5,REKAP,IOSTAT=errCode) 
      IF (errCode < 0) THEN
        WRITE(*,*) 'End of file reading /REKAP/ '
        STOP
      END IF
      INDTL(1) = NTDENS 
      INDTL(2) = NDGREE 
      INDTL(4) = NHORA 
      INDTL(5) = NWARM 
      X = END 
      Y = ALLNEW 
      RETURN 
END Subroutine Input
                                          
!$IBFTC STARTT                                                          
      SUBROUTINE STARTT 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!************COMMON FOR VARIABLE T-P-D INPUT*******FOR BETA-4 TABLE**** 
      DIMENSION TMDENS(1),TMT(1),TMP(1),XFORT(1),XFORP(1),XFORDN(1) 
!     ABOVE IS DUMMY DIMENSION -JUST TO ALLOW SUBSCRIBTING              
      EQUIVALENCE(TMDENS(1),DENS(301)),(TMT(1),T(301)),(TMP(1),P(301)) 
      EQUIVALENCE(XFORT(1),T(151)),(XFORP(1),P(151)),(XFORDN(1),DENS(151&
     &))                                                                
      EQUIVALENCE(NT,T(449)),(NXT,T(450)),(NP,P(449)),(NXP,P(450)) 
      EQUIVALENCE(ND,DENS(449)),(NXD,DENS(450)) 
!     THE T-P-AND DENS INPUTS ARE STORED IN THE                         
!     LAST THIRD OF THEIR RESPECTIVE ARRAYS.                            
!     NT-NP-NDENS ETC ARE STORED IN LOC                                 
!     NOT USED -DUE TO THE IMPOSSIBILITY OF HAVING 150 POINTS           
!     THE X TABLE FOR T(XFORT)ETC ARE STORED IN THE                     
!     MIDDLE THIRD OF THEIR RESPECTIVE ARRAYS                           
!     FOR THE BETA 4 TABLE HAVING PSEUDO - VARIABLE DIMENSION           
      COMMON TXQTAB(200),NOFTS,NOFXS,NTXQ(50),ICON,JCON,END 
!****         THIS ROUTINE IS A CONDUCTION SOLUTION DISCUSSED IN        
!****         SECTION VI OF THE ANALYSIS.                               
!****         IT DETERMINES THE STARTING TEMPERATURES FOR THE RUN.      
      DIMENSION COF(7) 
!****         COF IS USED IN THE ERROR FUNCTION COMPUTATION.            
      DATA COF/1.0,7.05230784E-2,4.22820123E-2,9.2705272E-3,1.520143E-4,&
     &2.765672E-4,4.30638E-5/                                           
      IF (GENRL(4))  1,80,1 
    1 LAYER = KONVAR(4) 
      KTIME = KONVAR(2) 
      NTEMP = INDTL(4) 
      NTIME = INDTL(5) 
      ITIME = KONVAR(1) 
      KT1 = KONVAR(2)+1 
      IT1 = KONVAR(1)+1 
      T(1) = GENRL(4) 
      T1 = T(1) 
      TIMEO = ABS(GENRL(3)) + GENRL(1) 
      DENS(1) = DENSE(1) 
      DO 9 J = 1,LAYER 
      NPTS = PTSIN(J) + 1. 
!****         STORE INITIAL T-S, AND DENSITIES IN THEIR PROPER ARRAYS   
      DO 8 I = 1,NPTS 
      MN = JLAYER(J) + I + 1 
    5 T(MN) = GENRL(4) 
    7 DENS(MN)=DENSE(J) 
    8 END DO 
    9 END DO 
!****              CLEAR Q-S FOR HTBLNS ROUTINE.                        
      DO 11 N = 1,12 
      Q(1,N) = 0.0 
   11 END DO 
   15 IF (GENRL(3)) 17,64,16 
!****         COMPUTE R1 FOR LINEAR HEAT FLUX                           
   16 R1 = 1. 
      GO TO 18 
!****         COMPUTE R1 FOR JUMP HEAT FLUX                             
   17 R1 = 1.5 
   18 CALL CHARP(1,1) 
      CALL TABUP(WARM,PROPT,T(1),PROP,NTEMP,-2) 
      CCPP = PROP(1) 
      CKAY = PROP(2) 
      CPA = CHARV(3) 
      AK = CHARV(1) 
      BARC1 = 4./3.*SQRT(TIMEO/3.1415927)*(1./SQRT(DENSE(1)*CCPP*CPA*AK*&
     &CKAY))*R1                                                         
      BARC2 = GENRL(4) 
!                       Q ITERATION BLOCK.                              
!                       IF GENRL(5) = -  LOOK UP FRONT FACE TEMP.       
      IF(GENRL(5)) 20,30,30 
   20 CALL TabUpSingle (HORA,TTABL,CTIME,TIVAR(1),NTIME,1) 
      T(1) = TIVAR(1) 
      Q(KT1,1) = (T(1)-GENRL(4))/BARC1 
      GO TO 45 
!                       IF GENRL(5) = +  CALCULATE THE FRONT FACE TEMP. 
   30 TOH = T(1) 
      T(IT1) = T(1) 
      EP = -EPTAB(3) 
!****         COMPUTE THE Q OF THE FRONT FACE                           
   31 CALL QNET 
!****         T(II) IS DEFINED IN EQUATION VI.6                         
      T(IT1) = GENRL(4)+BARC1*Q(KT1,1) 
      T2 = T(IT1) 
      IF(T2-WARM(NTEMP)) 32,32,33 
   32 IF(T2-WARM(1)) 34,38,38 
   33 T2 = WARM(NTEMP) 
      GO TO 38 
   34 T2 = WARM(1) 
   38 CALL TNEST(T(ITIME+1),T1,EP(1))    ! modified by RLC
      EP = ABS(EP) 
      NT1 = T1 
      GO TO (40,31,39),NT1 
   39 WRITE(6,500)  T(IT1),Q(KT1,1) 
      GO TO 31 
   40 T(1) = T(IT1) 
!                       END OF Q ITERATION BLOCK.                       
   45 DO 47 N = 1,6 
   47 Q(1,N)=Q(2,N) 
      NPTS = PTSIN(1)+2. 
      QPDNET = Q(KT1,1)/TIMEO 
!                       SET UP LOOP TO ACTUALLY CALCULATE               
!                       THE STARTING TEMPERATURES.                      
      DO 63 I = 1,NPTS 
      XLOCI = ABS(XLOC(1) - XLOC(I)) 
      ZEE = XLOCI/(2.*SQRT((CKAY*AK)/(DENSE(1)*CCPP*CPA)*TIMEO)) 
!****         ERFCZ IS DEFINED IN EQUATION VI.1                         
      IF(ZEE .EQ. 0.0) GO TO 147 
      CLCLCL=CLPOLY(ABS(ZEE),COF,6) 
!****         CHECK IF CLCLCL WILL OVERFLOW WHEN                        
!****         RAISED TO THE SIXTEENTH POWER                             
      IF(CLCLCL.GE.128.) GO TO 146 
      ERFCZ=SIGN(1.0,ZEE)*(1.0-(CLCLCL)**(-16)) 
      GO TO 160 
  146 ERFCZ=1.0 
      GO TO 160 
  147 ERFCZ = 0.0 
  160 PRT1 = (4./3.*R1*QPDNET*TIMEO**1.5)/SQRT(CKAY*AK*DENSE(1)*CCPP*CPA&
     &)                                                                 
      IF((ZEE**2)-80.) 60,60,59 
   59 PRT2=0. 
      PRT3=0. 
      GO TO 63 
   60 PRT2=(1.+ZEE**2)**(3.-2.*R1)/(SQRT(3.1415927))*EXP(-(ZEE**2)) 
      PRT3=(ZEE*((3.+2.*ZEE**2)/2.)**(3.-2.*R1))*(1.-ERFCZ) 
!****         T(I) IS DEFINED IN EQUATION VI.7                          
   63 T(I) = PRT1*(PRT2-PRT3)+GENRL(4) 
   64 CALL HTBLNS(0) 
      IF(KONVAR(7) - 2) 70, 65, 70 
!                       COMPUTE INTIAL PRESSURE                         
!                       DISTRIBUTION.                                   
!****              INITIAL TEMPERATURE AND RHO(G) MUST BE IN FIRST.     
   65 CALL DPRESS(1,1,-1) 
   70 RETURN 
!     AT 80 A VARIABLE TEMP-PRESS-DENS DISTRIBUTION IS                  
!     BEING READ IN- SO PUT VALUES IN                                   
!     CORRECT LOCATIONS                                                 
!     PLACE A 1 IN LAST ELEMENT OF EACH DEP TABLE                       
   80 TMT(NT+1) = 1. 
      TMP(NP+1) = 1. 
      TMDENS(ND+1) = 1. 
      KON2L = KONVAR(1) 
      DO 100 JD = 1,KON2L 
      CALL TabUpSingle(XFORT,TMT,XLOC(JD),T(JD),NT,1) 
      CALL TabUpSingle(XFORP,TMP,XLOC(JD),P(JD),NP,1) 
      CALL TabUpSingle(XFORDN,TMDENS,XLOC(JD),DENS(JD),ND,1) 
  100 END DO 
!     ERASE LAST TWO THIRDS OF T-P AND DENS ARRAYS                      
      DO 110 JD = 1,150 
      TMDENS(JD) = 0. 
      XFORDN(JD) = 0. 
      TMP(JD) = 0. 
      XFORP(JD) = 0. 
      TMT(JD) = 0. 
      XFORT(JD) = 0. 
  110 END DO 
      GO TO 70 
  500 FORMAT (/,3X,'TNEST IN STARTT',4X,'T(F.F.) =', E20.8,4X,'QNET = ',E20.8)
      STOP 
END Subroutine Startt                                          


!$ORIGIN        ALPHA                                                   
!$IBFTC AIRGAP                                                          
      SUBROUTINE AIRGAP 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE AIR-GAP EQUATIONS ARE DEFINED AND DISCUSSED IN        
!****         SECTION IV OF THE ANALYSIS.                               
!                       TEST IF KTIME = 1.                              
!                       SET UP COUNTERS FOR                             
!                       FIRST TIME THROUGH                              
      KTIME = KONVAR(2) 
      IF (KTIME-1) 10,5,10 
    5 JGAP = KONVAR(5) 
      J1 = JGAP+1 
      NTEMP = INDTL(4) 
      LAYER = KONVAR(4) 
      ITIME = KONVAR(1) 
      J2 = JGAP-1 
      MN = JLAYER(JGAP)+1 
   10 II = KONVAR(6)*KONVAR(1)+MN 
      JJ = II-ITIME 
      KK = JJ-ITIME 
!                       COMPUTE TERMS IN THE C COEFFICIENTS             
!                       THAT REMAIN CONSTANT THROUGOUT                  
!                       THE ITERATIONS.                                 
   11 AGTAV = (T(II)+T(II+1))/2. 
!                       USE TABLE LOOK-UP TO GET THE FOLLOWING          
!                       1.  THE K FOR THE (J-1) LAYER = CAYL.           
!                       2.  THE K FOR THE (J+1) LAYER = CAYR.           
!                       3.  THE C(P) FOR THE AIR-GAP = AGCP.            
!                       4.  THE RHO-C(PL) OF THE (J-1) LAYER = RCPLL.   
!                       5.  THE RHO-C(PL) OF THE AIR-GAP = RCPLAG.      
!                       6.  THE K FOR THE AIR-GAP = CAYAG               
      MN1 = (NTEMP+1)*(3*(J2-1)+1)+1 
      CALL TabUpSingle(WARM,PROPT(MN1:),T(II),CAYL,NTEMP,1) 
!                       TEST TO SEE IF AIR-GAP IS LAST LAYER.           
      IF (JGAP-LAYER) 18,16,18 
   16 B2 = 0. 
      GO TO 19 
   18 MN4 = (NTEMP+1)*(3*(J1-1)+1)+1 
      CALL TabUpSingle(WARM,PROPT(MN4:),T(II+1),CAYR,NTEMP,1) 
      B2 = -CAYR/THIKN(J1)/2./DELETA(J1) 
   19 CALL TabUpSingle(WARM,CTABL,T(II),RCPLL,NTEMP,J2) 
      CALL TabUpSingle(WARM,CTABL,T(II+1),RCPLAG,NTEMP,JGAP) 
      B1 = -CAYL*ETASBX(MN)/2./DELETA(J2) 
      MN3 = (NTEMP+1)*(3*(JGAP-1))+1 
      CALL TabUpSingle(WARM,PROPT(MN3:),AGTAV,AGCP,NTEMP,1) 
      RCPA2 = DENSE(JGAP)*AGCP*THIKN(JGAP)/2. 
      TMP1 = RCPLAG+RCPA2 
      TMP2 = RCPLL+RCPA2 
!****          SEE IF 2 PT. OR 3 PT. FORMULA IS TO BE USED.             
!****          GENRL(6) = - , USE 2 PT.     IF =  +  USE 3 PT.          
      IF(GENRL(6)) 21,23,23 
   21 ALF3 = 0.0 
      GO TO 24 
   23 ALF3 = 1.0 
   24 BC1 = (2.+ALF3)*B2 - GAMMA(6)*TMP1 
      BC2 = (2.+ALF3)*B1 - GAMMA(6)*TMP2 
      BC5 = 1./((2.+ALF3)*B1 - RCPLL*GAMMA(6)) 
      B3 = TMP1*(GAMMA(4)*T(KK+1)+GAMMA(5)*T(JJ+1)) 
      B4 = TMP2*(GAMMA(4)*T(KK)+GAMMA(5)*T(JJ)) 
      B5 = RCPLL*(GAMMA(4)*T(KK)+GAMMA(5)*T(JJ)) 
!                       ITERATION SCHEME WILL START HERE.               
      EP1 = -EPTAB(4) 
      TX5 = T(II) 
      TX5I = T(II) 
   26 TX6 = T(II+1) 
      EPS = -EPTAB(4) 
  261 CONTINUE 
      IF (KONVAR(3)-2) 28,27,28 
   27 CALL INTPTS(JGAP,2) 
   28 TMP3 = ((2.+2.*ALF3)*T(II-1)-ALF3*T(II-2))*B1 
      BC3 = B2*((2.+2.*ALF3) * T(II+2) -ALF3*T(II+3)) +TMP3 + B3 + B4 
!****         T(II+1) IS DEFINED IN EQUATION IV.22                      
      T(II+1) = (BC3-BC2*T(II))/BC1 
      IF (KONVAR(3)-2) 30,29,30 
   29 CALL TNEST(T(II+1),TX6,EPS) 
      EPS = ABS(EPS) 
      KTX6 = TX6 
      GO TO (30,261,127),KTX6 
  127 WRITE(6,500)  BC1,BC2,BC3,BC4,BC5,B1,B2,B3,B4,B5,T(II) 
      GO TO 261 
   30 BC4 = BC5*(B1*((2.+2.*ALF3)*T(II-1)-ALF3*T(II-2))+B5) 
!                       BLOCK FOR BAR-C-6.                              
      CALL QGAP(BC6,QDCOND,QR,QCR) 
      BC6 = BC6*(T(II)-T(II+1)) 
!****         T(II) IS DEFINED IN EQUATION IV.21                        
      T(II) = BC4+BC5*BC6 
      CALL TNEST(T(II),TX5,EP1) 
      EP1 = ABS(EP1) 
      COUNT(13) = COUNT(13) + 1. 
      KTX5 = TX5 
      GO TO (60,26,55),KTX5 
   55 WRITE(6,501)  QDCOND,QR,QCR,BC4,BC5,BC6,T(II+1) 
      GO TO 26 
   60 RETURN 
  500 FORMAT(/,4X,'C BARS = ',5E20.8//7X,'B-S = ',5E20.8//5X,'T(II) = ',E18.8)
  501 FORMAT (/,4X,'QC = ',E18.8,3X,'QR = ',E18.8,3X,'QCR = ',E18.8// &
       3X,'BC4 = ',E18.8,2X,'BC5 = ',E18.8,3X,'BC6 = ',E18.8//3X,'T(II+1) = ',E18.8)
      STOP 
END Subroutine Airgap
                                          
!$IBFTC DELTAT                                                          
      SUBROUTINE DELTAT(K) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!                       THIS ROUTINE CALCULATES THE DELTA               
!                       T FOR THE TIME STEP                             
!                       K IS A COUNTER WHICH DETERMINES                 
!                       IF DELTA T IS CHOSEN FROM TABLE                 
!                       (K=0) OR MUST BE MODIFIED (K=1).                
      KAY = K 
      KTIME = KONVAR(2) 
      NLAYER = KONVAR(4) 
      IF(KAY) 180, 2, 200 
!****     THIS COUNT IS RESET EACH TIME KAY IS ZERO                     
    2 KOUNTR = 0 
!                       GET CORRECT DELTA BAR FROM DTIME(X).            
      DO 4 I = 1,10 
      N = 2*I 
      N1 = N-1 
      IF(DTIME(N) - CTIME - DELTO(KTIME)) 4,4,3 
    3 DELBAR = DTIME(N1) 
      GO TO 30 
    4 END DO 
!                       SET DELTO(KTIME) = DELBAR                       
   30 TEST1 = ABS(ETASBT(1)) 
      LMSRCH = KONVAR(1) 
      DO 35 IJSRCH = 1, LMSRCH 
      IF (TEST1-ABS(ETASBT(IJSRCH))) 31,33,33 
   31 TEST1 = ABS(ETASBT(IJSRCH)) 
   33 IF (TEST1-ABS(BETA(3,IJSRCH))) 34,35,35 
   34 TEST1 = ABS(BETA(3,IJSRCH)) 
   35 END DO 
      DELTO(KTIME) = DELBAR 
      IF(KTIME-1) 56, 56, 50 
   50 E1 = (1.+EPTAB(1))*DELTO(KTIME-1) 
!                       KTIME IS GREATER THAN ONE.                      
!                       SELECT CORRECT DELTO(KTIME).                    
      DELZMN = AMIN1(DELBAR,E1) 
      IF (DELBAR-DELTO(KTIME-1))  56,55,55 
   55 DELTO(KTIME) = DELZMN 
   56 IF(ABS(TEST1*DELTO(KTIME)/DELETA(1))-1.) 180,58,58 
   58 DELTO(KTIME) = DELETA(1)*.95/ABS(TEST1) 
      COUNT(15) = COUNT(15) + 1. 
      IF(DELTO(KTIME) -.01*DELTO(KTIME-1)) 190,190,180 
  180 BB = -HTCONV(5) * ETASBX(1) 
      KOUNTR = KOUNTR +1 
      IF (KOUNTR.GT.10) GO TO 191 
!****     THE CONSTANT IN THE FOLLOWING STATEMENT                       
!****     MUST BE LARGER THAN THE ONE IN STATEMENT 186                  
      IF(ABS(DELTO(KTIME)*BB/DELETA(1))-.5) 200,200,184 
  184 IF(KAY.EQ.0) GO TO 186 
  185 CTIME = CTIME - DELTO(KTIME) 
      K = 0 
!****              COMPUTE R VALUES,ETA DERIVATIVES AND X LOCATIONS.    
      CALL RITER(1) 
  186 DELTO(KTIME) = .45*DELETA(1)/ABS(BB) 
!****     SEE PREVIOUS COMMENT CONCERNING CONSTANT                      
      IF(KAY.EQ.0) GOTO 189 
      CTIME = CTIME + DELTO(KTIME) 
  189 CONTINUE 
      COUNT(14) = COUNT(14) + 1. 
      IF(DELTO(KTIME) - .01 * DELTO(KTIME-1)) 190, 190, 200 
  190 WRITE (6,500) 
      STOP   ! changed from CALL EXIT by RLC 7Sep05 
  191 WRITE (6,501) 
      STOP   ! changed from CALL EXIT by RLC 7Sep05 
  200 IF (TIVAR(10)) 205,240,205 
!****     THE CONSTANT IN STATEMENT NUMBER 205 MUST BE GREATER          
!****     THAN THE ONE IN STATEMENT NUMBER 225                          
  205 IF (DELTO(KTIME) - .50*DELETA(1)/TIVAR(10)) 230, 230, 225 
  225 DELTO(KTIME) = .45 * DELETA(1)/ TIVAR(10) 
      COUNT(14) = COUNT(14) +1. 
  230 TIVAR(10) = 0.0 
!****         SEE IF IMPLICIT SCHEME IS REQUESTED FOR THIS TIME STEP.   
  240 CONTINUE 
      IMPTME = 0 
  245 IF (CTIME - AMPLCT(IMPTME+1)) 260,250,250 
  250 IF (CTIME - AMPLCT(IMPTME+2)) 255,258,258 
!****         AT 255 THE IMPLICIT SCHEME HAS BEEN REQUESTED             
!****         FOR THIS TIME STEP.                                       
  255 IMEX = -1 
      GO TO 400 
  258 IMPTME = IMPTME + 2 
      IF (IMPTME.LE.8 ) GO TO 245 
!****         AT 260 THE IMPLICIT SCHEME HAS NOT BEEN REQUESTED         
!****         FOR THIS TIME STEP SO DECIDE BETWEEN EXP. AND EXP-IMP.    
!****         J = THE LAYER.                                            
!****         MN = THE X LOCATION   MN IS THE TRUE MN.                  
  260 DO 268 J = 1,NLAYER 
      MNBGIN = JLAYER(J) + 2 
      TEST = BETA(1,MNBGIN) 
      MNBGIN = MNBGIN + 1 
      MNSTOP = JLAYER(J+1) 
      DO 263 MN = MNBGIN , MNSTOP 
      IF(TEST - BETA(1,MN))  262,263,263 
  262 TEST = BETA(1,MN) 
  263 END DO 
!****WE NOW HAVE THE LARGEST BETA FOR LAYER J.                          
      IF(J-1)  266,265,266 
  265 DELGES = DELETA(J)**2/(2.*TEST) 
      GO TO 268 
  266 TSTDEL = DELETA(J)**2/(2.*TEST) 
      IF(TSTDEL - DELGES)  267,268,268 
  267 DELGES = TSTDEL 
  268 END DO 
      IF(DELGES - DELTO(KTIME)) 275,270,270 
  270 IMEX = 1 
      KONVAR(3) = 3 
      GO TO 400 
  275 IMEX = 0 
  400 RETURN 
  500 FORMAT(/,4X,'DT CUT 2 ORDERS OF MAGNITUDE')
  501 FORMAT(/,4X,'ITERATING TOO MUCH IN DELTAT')
END Subroutine Deltat
             
                             
!$IBFTC GAMA                                                            
      SUBROUTINE GAMA(HH) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATIONS FOR THE APPROXIMATION COEFFICIENTS ARE      
!****         DEFINED AND DISCUSSED IN SECTION IV,THE INTERFACE         
!****         EQUATIONS.                                                
!                       THIS SUBROUTINE CALCULATES THE                  
!                       APPROXIMATION COEFFICIENTS CALLED               
!                       GAMMA(4-6).                                     
!****         GAMMA(1-3) ARE AVAILABLE FOR FUTURE MODIFICATION.         
      H=HH 
      KTIME = KONVAR(2) 
    1 H = -H 
      VRDEL2 = DELTO(KTIME) 
      GO TO (2,3,3),KTIME 
    2 GAMMA(4) = 0.0 
      GAMMA(5) = -1./VRDEL2 
      GAMMA(6) = 1./VRDEL2 
      GO TO 4 
    3 VRDEL1 = DELTO(KTIME-1) 
      GAMMA(4) = VRDEL2/(VRDEL1*(VRDEL1+VRDEL2)) 
      GAMMA(5) = -(VRDEL1+VRDEL2)/(VRDEL1*VRDEL2) 
      GAMMA(6) = (VRDEL1+2.*VRDEL2)/(VRDEL2*(VRDEL1+VRDEL2)) 
    4 RETURN 
END Subroutine Gama
                                          
!$IBFTC INTFC                                                           
      SUBROUTINE INTFC 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATIONS OF THIS ROUTINE ARE DEFINED AND             
!****         DISCUSSED IN SECTION IV, NUMERICAL SOLUTIONS OF ENERGY    
!****         EQUATIONS AND BOUNDARY CONDITIONS.                        
!                       THIS SUBROUTINE CALCULATES THE                  
!                       INTER-FACE AND BACK-FACE                        
!                       TEMPERATURES AND PRESSURES.                     
!                       SUBROUTINE GAMA,CALLED WITH A                   
!                       (-) ARGUEMENT, CALCULATES GAMMA(1-6).           
      ITIME = KONVAR(1) 
      KTIME = KONVAR(2) 
      LAYER = KONVAR(4) 
      JGAP = KONVAR(5) 
      N3 = KONVAR(6) 
      NTEMP = INDTL(4) 
      NTIME = INDTL(5) 
!****          SEE IF 2 PT. OR 3 PT. FORMULA IS TO BE USED.             
!****          GENRL(6) = - , USE 2 PT.     IF =  +  USE 3 PT.          
      IF(GENRL(6)) 2,3,3 
    2 ALF3 = 0.0 
      GO TO 4 
    3 ALF3 = 1.0 
    4 DO 30 J = 1,LAYER 
!                       LOOP COMPUTES INTERFACE TEMPERATURES.           
      L = J+1 
      TIVAR(5) = 0. 
      IF (L-JGAP) 5,30,5 
    5 IF (JGAP+1-L) 7,30,7 
    7 NPTS = PTSIN(J)+1. 
      MN = JLAYER(J)+NPTS+1 
      II = N3*ITIME+MN 
      JJ = II-ITIME 
      KK = JJ-ITIME 
      TY4 = T(II) 
      EPS = -EPTAB(5) 
    8 MN1 = (NTEMP+1)*(3*(J-1)+1)+1 
      CALL TabUpSingle(WARM,PROPT(MN1:),T(II),CAY1,NTEMP,1)   ! modified by RLC
      BBAR1 = -CAY1*ETASBX(MN)/DELETA(J)/2. 
      IF (J-LAYER) 22,10,22 
   10 BBAR2 = 0. 
      CALL TabUpSingle(HORA,TTABL,CTIME,TIVAR(5),NTIME,5) 
      IF (GENRL(7)) 20,23,20 
   20 T(II) = TIVAR(5) 
      GO TO 23 
   22 MN1 = (NTEMP+1)*(3*(L-1)+1)+1 
      CALL TabUpSingle(WARM,PROPT(MN1:),T(II),CAY2,NTEMP,1)   ! modified by RLC
      BBAR2 = -CAY2*ETASBX(MN+1)/DELETA(L)/2. 
!                       TEST IF ROUTINE MUST                            
!****                   CALL INTPTS.                                    
   23 IF (KONVAR(3)-2) 26,24,26 
   24 CALL INTPTS(L,2) 
!                       COMPUTE TERMS FOR TEMP. EQUATION.               
   26 TRM1 = BBAR1*((2.+2.*ALF3)*T(II-1) - ALF3*T(II-2)) 
   27 IF(J-LAYER) 28,127, 127 
  127 IF (GENRL(7)) 30,117,30 
  117 QBFC = TIVAR(5) 
      GO TO 128 
   28 QBFC = 0. 
  128 TRM2 = BBAR2*((2.+2.*ALF3)*T(II+1) - T(II+2)*ALF3) + QBFC 
      CALL TabUpSingle(WARM,CTABL,T(II),SEEUP,NTEMP,J) 
      TRM3 = SEEUP*(GAMMA(4)*T(KK) + GAMMA(5)*T(JJ)) 
!****         T(II) IS DEFINED IN EQUATION IV.14                        
!****         WHEN BBAR2 = 0.0, T(II) IS DEFINED BY EQUATION IV.16.     
      T(II) = (TRM1 + TRM2 + TRM3)/(BBAR1*(2.+ALF3) + BBAR2*(2.+ALF3)   &
     &-SEEUP*GAMMA(6))                                                  
!                       CALL TNEST FOR ITERATION.                       
      IF (KONVAR(3)-2) 30,29,30 
   29 CALL TNEST(T(II),TY4,EPS) 
      EPS = ABS(EPS) 
      KTY4 = TY4 
      GO TO (30,123,150),KTY4 
  123 COUNT(9) = COUNT(9) + 1. 
      GO TO 23 
  150 WRITE(6,500)  BBAR1,BBAR2,SEEUP,TRM1,TRM2,TRM3 
      GO TO 23 
   30 END DO 
  131 NPTS2 = PTSIN(1) + 2. 
      MND = JLAYER(1) + NPTS2 
      IID = KONVAR(6)*KONVAR(1) + MND 
  233 IF(KONVAR(7)-2) 139, 133, 139 
  133 CALL DPRESS(IID,MND,0) 
  139 RETURN 
  500 FORMAT(/,4X,'BBAR1 = ',E18.8,3X,'BBAR2 = ',E18.8,3X,'SEEUP = ',E18.8// &
       4X,'TRM1 = ',E18.8,3X,'TRM2 = ',E18.8,3X,'TRM3 = ',E18.8)        
      STOP 
END Subroutine Intfc  
                                        
!$IBFTC PTX1                                                            
      SUBROUTINE PTX1 
USE TimeDt
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATIONS IN THIS SUBROUTINE ARE DISCUSSED IN SECTION 
!****         IV, PART B OF THE ANALYSIS.                               
!                       THIS SUBROUTINE CALCULATES                      
!                       THE TEMPERATURE AND PRESSURE                    
!                       OF THE FRONT FACE.                              
!****         K(RHO) IS NOT INCLUDED IN THE TEMPERATURE ITERATION.      
      NTEMP = INDTL(4) 
      NTIME = INDTL(5) 
      II = KONVAR(6)*KONVAR(1) + 1 
      CALL CHARP(II,1) 
    5 BARC1T = DELETA(1)/1.5/ETASBX(1)/CHARV(1) 
    7 TTIME = CTIME 
!****              TEST TO SEE IF FRONT FACE T IS READ IN.              
      IF(GENRL(5)) 107, 8, 8 
  107 CALL TabUpSingle(HORA,TTABL,TTIME,TIVAR(1),NTIME,1) 
      T(II) = TIVAR(1) 
!****          SEE IF 2 PT. OR 3 PT. FORMULA IS TO BE USED.             
!****          GENRL(6) = - , USE 2 PT.     IF =  +  USE 3 PT.          
    8 IF(GENRL(6)) 108,109,109 
  108 ALF3 = 0.0 
      GO TO 110 
  109 ALF3 = 1.0 
  110 EPS = -EPTAB(3) 
      T1 = T(II) 
!                       TEST NZERO TO SEE IF                            
!                       SUBROUTINE MIDPTS NEEDS                         
!                       TO BE CALLED.                                   
    9 IF (KONVAR(3)-2) 15,10,15 
   10 N4 = PTSIN(1)+1. 
      CALL INTPTS(1,N4) 
   15 C3 = 2. + ALF3 
      C2 =(2. + 2.*ALF3)/C3*T(II+1) - ALF3/C3*T(II+2) 
      IF (CHART(1)) 20,27,20 
!****         COMPUTE M-DOT-G IF CHARRING                               
   20 IF (IMEX) 21,21,22 
   21 CALL DENSIT(II,1,0) 
   22 CALL EMG(II) 
!****         COMPUTE THE Q OF THE FRONT FACE                           
   27 CALL QNET 
      IF(GENRL(5)) 38,28,28 
   28 N5 = KONVAR(2)+1 
!****         PRESSURE OPTION                                           
      IF(KONVAR(7)-2) 18, 16, 18 
   16 CALL DPRESS(II,1,0) 
      CALL KEFFS(II,1) 
      CAYUP = HTCONV(15) 
      GO TO 19 
   18 CALL TabUpSingle(WARM,PROPT,T(II),CAYUP,NTEMP,2) 
   19 BARC1 = BARC1T/CAYUP 
      C1 = 3. * BARC1/C3 
!****         T(II) IS DEFINED IN EQUATION IV.11                        
      T(II) = C1 * Q(N5,1) + C2 
!****              TEST TO SEE IF FIXED MELTING T OPTION.               
      IF (KONVAR(8)-5) 32,29,32 
   29 IF (HTCONV(5)) 32,32,30 
   30 IF (T(II)-TMELT(2)) 32,32,31 
   31 IF (EPS)35,32,32 
   32 T(II) = AMAX1(T(II),WARM(1)) 
      CALL TNEST(T(II),T1,EPS) 
      EPS = ABS(EPS) 
      COUNT(5) = COUNT(5) + 1. 
      JT1 = T1 
      GO TO (134,9,234),JT1 
  234 WRITE(6,500)ALF3,C1,C2,Q(N5,1),DENS(II),DENS(II+1),T(II+1),EMDG(1) 
      GO TO 9 
  134 IF(T(II)-WARM(1)) 41, 41, 34 
   34 IF (KONVAR(8)-5) 38,35,38 
   35 T(II) = AMIN1 (T(II),TMELT(2)) 
      IF (TMELT(2)-T(II)) 38,36,38 
!                       RECOMPUTE QUANTITIES INVOLVED IN THE            
!                       FRONT FACE TEMP. CALCULATIONS.                  
   36 IF (KONVAR(3)-2) 38,37,38 
   37 CALL INTPTS(1,N4) 
   38 IF ( CHART(1) .EQ. 0.0 ) GO TO 139 
      CALL DENSIT(II,1,0) 
      IF(KONVAR(7)-2) 139,138,139 
  138 CALL DPRESS(II,1,0) 
  139 CALL QNET 
   40 RETURN 
   41 CALL ErrorMessage('T(II) = WARM(1)')
  500 FORMAT(/,4X,'ALF3 = ',E18.8,3X,'C1 = ',E18.8,3X,'C2 = ',E18.8,3X, &
       'QDNET = ',E18.8//4X,'DENS(1) = ',E18.8,3X,'DENS(2) = ',E18.8,3X, &
       'T(2) = ',E18.8,3X,'MDG(1) = ',E18.8)                                  
      STOP 
END Subroutine Ptx1        
                                  
!$IBFTC PTX12                                                           
      SUBROUTINE PTX12 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****         THE EQUATIONS USED IN THIS SUBROUTINE ARE DISCUSSED IN    
!****         APPENDIX E, THIN SKIN OPTION, OF THE ANALYSIS.            
      KTIME = KONVAR(2) 
      II = KONVAR(6)*KONVAR(1) + 1 
      JJ = II - KONVAR(1) 
      AK = 1. 
      CPA = 1. 
      AKP = 0. 
      IF(KTIME-1) 5, 5, 10 
    5 T(JJ) = T(1) 
      T(II) = T(1) 
      DENS(JJ) = DENS(1) 
      DENS(II) = DENS(1) 
      GO TO 15 
   10 DELTR = DELTO(KTIME)/DELTO(KTIME-1) 
      T(II) =(DELTR + 1.)*T(JJ)-DELTR*T(1) 
      DENS(II) = (DELTR + 1.)*DENS(JJ)-DELTR*DENS(1) 
      CALL ABLATE(DELTO(KTIME)) 
   15 T1 = T(II) 
      EPS = -EPTAB(3) 
   20 N5 = KTIME + 1 
      CALL TabUpSingle(WARM,PROPT,T(II),CP,INDTL(4),1) 
      BARC1 = DELTO(KTIME)/THIKN(1)/DENSE(1)/CP 
!****         COMPUTE THE Q OF THE FRONT FACE                           
      CALL QNET 
      T(II) = T(JJ)+BARC1*Q(N5,1) 
      T(II) = AMAX1(T(II),WARM(1)) 
      CALL TNEST(T(II),T1,EPS) 
      EPS = ABS(EPS) 
      COUNT(5) = COUNT(5) + 1. 
      KT1 = T1 
      GO TO (50,20,30),KT1 
   30 WRITE(6,500) T(II),T(JJ),BARC1,Q(N5,1) 
      GO TO 20 
   50 IF(T(II)-WARM(1)) 75, 75, 60 
   60 IF(KONVAR(8)-5)68,65,68 
   65 T(II) = AMIN1(T(II),TMELT(2)) 
      IF(TMELT(2)-T(II)) 68, 66, 68 
   66 CALL QNET 
   68 CALL ABLATE(0.0) 
      CALL OUTPUT(0) 
   72 RETURN 
   75 CALL ErrorMessage('T(II) = WARM(1)')
      GO TO 68 
  500 FORMAT(/,3X,'T(II) = ',E18.8,3X,'T(JJ) = ',E18.8,3X,'BARC1 = ', &
       E18.8,3X,'QNET = ',E18.8)
      STOP 
END Subroutine Ptx12
                                          
!$IBFTC QGAP                                                            
      SUBROUTINE  QGAP(BC6,QDCOND,QR,QCR) 
USE TimeDt

!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 
                                                                        
!****          THIS SUBROUTINE COMPUTES THE Q OF THEAIRGAP.             
!****          THE EXPRESSION (T(II)-T(II+1)) HAS BEEN FACTORED         
!****          OUT OF  QCR---QR---QDCOND.                               
      II = KONVAR(6)*KONVAR(1) + JLAYER(JGAP)+1 
      NTEMP = INDTL(4) 
      PT = CTIME-GENRL(1) 
      IF(KONVAR(2)-1)  10,5,10 
    5 JGAP = KONVAR(5) 
      AGTAV = (T(II) + T(II+1))/2. 
!                       COMPUTE COEFFICIENTS FOR                        
!                       FOR PARABOLIC CURVE FIT.                        
      G0 = AGPTB(7) 
      G1 = AGPTB(8) 
      G2 = AGPTB(9) 
      T1 = GENRL(2)-GENRL(1) 
      A1 = (G0-2.*G1+G2)*2./T1**2 
      A2 = (-3.*G0+4.*G1-G2)/T1 
      A3 = G0 
      MN2 = (NTEMP+1)*(3*(JGAP-1)+1)+1 
      CALL TabUpSingle(WARM,PROPT(MN2:),AGTAV,CAYAG,NTEMP,1)   ! modified by RLC
      MN3 = (NTEMP+1)*(3*(JGAP-1))+1 
      CALL TabUpSingle(WARM,PROPT(MN3:),AGTAV,AGCP,NTEMP,1)   ! modified by RLC
   10 G = A1*PT**2+A2*PT+A3 
   15 QR = HTCONV(2)*AGPTB(2)*AGPTB(3)*(T(II)+T(II+1))*                 &
     &(T(II)**2+T(II+1)**2)                                             
   20 QDCOND = AGPTB(4)/THIKN(JGAP) 
      BETDTM = 2.*(T(II)-T(II+1))/(T(II)+T(II+1)) 
      BETDTM = ABS(BETDTM) 
      IF(G)  40,30,40 
   30 QCR = 0.0 
      GO TO 100 
   40 IF(AGPTB(5)) 41,30,41 
   41 GR = THIKN(JGAP)**3*DENSE(JGAP)**2*G*BETDTM/AGPTB(5)**2 
      IF(GR - AGPTB(11))  50,50,45 
   45 C = AGPTB(13) 
      EN = AGPTB(15) 
      GO TO 70 
   50 C = AGPTB(12) 
      EN = AGPTB(14) 
   70 QCR = C*CAYAG/THIKN(JGAP)*(THIKN(JGAP)/AGPTB(6))**AGPTB(10) 
      QCR = QCR*(GR*((AGCP*AGPTB(5))/CAYAG))**EN 
  100 BC6 = QCR + QR + QDCOND 
      RETURN 
END Subroutine Qgap  
                                        
!$IBFTC TMPIMP                                                          
      SUBROUTINE TMPIMP(LEM) 
!     THIS IS A DUMMY SUBROUTINE AND SHOULD NOT BE ENTERED              
      STOP 
END Subroutine TmpImp

END Module AblateProcedures   ! =============================================





PROGRAM AblativeNozzleMaterials
! ---------------------------------------------------------------------------

USE AblateProcedures
USE TimeDt

!$IBFTC MAIN                                                            
      COMMON TXQTAB(200),NOFTS,NOFXS,NTXQ(50),ICON,JCON,END 
!****                         INPUT     COMMON                          
                                                                        
      COMMON  /INPTC/  CHART( 25),   GASPR(105),  TDENS( 20),           &
     &                   ZORET( 84)                                     
                                                                        
      COMMON  /INPTC /   INDTL( 25),   NATETB( 20),   NCTABL( 5),       &
     &                   NDGREE    ,   NGASPR(  5),   NHORA     ,       &
     &                  NHTABL( 2),    NNU   (  5),   NPHI      ,       &
     &                  NPROPT(30),    NTDENS     ,   NTTABL(10),       &
     &                  NWARM      ,   NZORET(  4)                      
                                                                        
      COMMON  /INPTC /   TMELT( 25),    TRAJT( 25) 
                                                                        
      COMMON  /INPTC /  AGPTB( 25),    DENSE( 10),    DTIME (25),       &
     &                  EPTAB( 10),    GENRL( 25),    OUTTB( 25),       &
     &                  PTSIN( 10),    SPACE ( 10),   THIKN (10)        
                                                                        
      COMMON  /INPTC /   PVAR ( 50) 
                                                                        
      COMMON  /INPTC /  CTABL(210),    HTABL ( 42),   PROPT(630),       &
     &                  WARM ( 20)                                      
                                                                        
      COMMON  /INPTC /   HORA ( 50),    TTABL(510) 
                                                                        
      COMMON  /INPTC /  AT    (150),   ATETB (420),   ATP  (150),       &
     &                  DGREE( 20),    E     (150),   EP   (150),       &
     &                  NU    (  5),   PHII  (105),   PPHIP(150),       &
     &                  PPHI2P(150),   SIGMA (150),   STRESS(150),      &
     &                  STVAR ( 25)                                     
                                                                        
                                                                        
!****                              COMPUTATION COMMON                   
                                                                        
      COMMON  /CHARDT/   CHARV (  5),  SCHAR (  4) 
                                                                        
      COMMON  /COUNTR/   COUNT ( 25) 
                                                                        
      COMMON  /HEATDP/   Q    (4,13) 
                                                                        
      COMMON  /LAYERD/   DELETA( 10),  JLAYER( 11) 
                                                                        
      COMMON  /MELTDT/   SMELT (  4),  TRAJV ( 50) 
                                                                        
      COMMON  /PHYSDT/   DELTO (  4),  SPACVR(  5),    GAMMA (  6) 
                                                                        
      COMMON  /POINTD/  BETA(4,150),   DENS  (450),   EMDG (150),       &
     &                  ETASBT(150),   ETASBX(150),   ETAXX(150),       &
     &                  P      (450),  T     (450),   WDOTP(150),       &
     &                  XLOC  (150)                                     
                                                                        
!!!      COMMON  /TIMEDT/   CTIME     ,    TIME 
                                                                        
      COMMON  /VARIBL/   HTCONV( 25),   KONVAR( 25) 
      COMMON  /VARIBL/  TIVAR(10)   ,  PROP    (3),   HEVAR  (10) 
                                                                        
      COMMON  /OUTDPT/  AUS(25) 
                                                                        
!      COMMON  /HEADG /  HEAD(14) 
      COMMON/EXPIMP/ AMPLCT(10),IMEX,IMPTME,BACKFC(4) 

  INTEGER:: errCode
  CHARACTER(LEN=80):: fileName
!----------------------------------------------------------------------------
  WRITE(*,*) 'ablate - Ablative nozzle materials'
  DO
    WRITE(*,*) 'Enter the name of the input file: '
    READ(*,'(A)') fileName
    IF (LEN_TRIM(fileName)==0) STOP
    OPEN(UNIT=5,FILE=fileName,STATUS='OLD',IOSTAT=errCode,ACTION='READ')
    IF (errCode==0) EXIT
    WRITE(*,*) 'Unable to open this file. Try again.'
  END DO
  OPEN(UNIT=6,FILE='ablate.out',STATUS='REPLACE',ACTION='WRITE')
  CALL InitializeCommon()
!      X=X 
      ICON = 1 
    1 CALL EXEC1()
      CALL EXEC2()
      GO TO 1 
      END Program AblativeNozzleMaterials   ! ===============================
