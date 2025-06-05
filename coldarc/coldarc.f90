      MODULE ColdArcProcedures

!      IMPLICIT NONE

      CONTAINS


!
      SUBROUTINE PHGETT (P,H,T)
!     S.N.GREENSCHLAG - RESPONSIBLE PROGRAMMER
!
      COMMON CR(200)
      EQUIVALENCE                                                       &
     &   (CR(   3),  SORC   ),                                          &
     &   (CR(  92),  IERR   )
      DATA EPS/1.0E-05/
!
!!!      IF (P) 2,2,1
!!!    1 IF (T) 2,2,5
  IF (p<=0.0 .OR. t<=0.0) THEN
    WRITE(6,*) 'Subroutine Phgett abort exit!'
    WRITE(6,*)'p=',P,'  T=',T
    IERR = 1
    RETURN
  END IF 

     IC = 0
  490 PC=P                          ! label 490 not used
      KC = 0
      TI=T
!
  500 IC=IC+1
      TC=TI
      IP=-1
!
  510 CALL TNP (PC,TC,HC)
  IF(IERR .NE. 0) THEN
    WRITE(6,*) 'Subroutine Phgett abort exit!'
    WRITE(6,*)'p=',P,'  T=',T
    RETURN
  END IF

      F=HC-H
!
!*****ITERATION VALUES CLOSEST TO THE ROOT ARE SAVED
      IF(IC.GT.1) GO TO 512
      FSAV1=0.0
      TSAV1=0.0
      FSAV2=0.0
      TSAV2=0.0
  512 IF(F.LT.0.0) GO TO 515
      IF(ABS(F).GE.ABS(FSAV1).AND.FSAV1.NE.0.0) GO TO 520
      FSAV1 = F
      TSAV1 = TC
      GO TO 520
  515 IF(ABS(F).GE.ABS(FSAV2).AND.FSAV2.NE.0.0) GO TO 520
      FSAV2 = F
      TSAV2 = TC
  520 IF(ABS(F/H)-EPS) 530,530,540
!
  530 SOR=SORC
      T=TC
      GO TO 1000
  540 IF(IC.GE.25) GO TO 800
!
      IF(IP) 580,680,680
  580 IF(IC-1)620,620,590
!
  590 IF(DFDT)600,620,600
!
  600 DT=-F/DFDT
      GO TO 630
!
  620 DT=-F
!
  630 IF( ABS(DT)-.01*TI)640,670,650
!
  640 DT= SIGN(.01*TI,DT)
      GO TO 670
!
  650 IF( ABS(DT)-.9*TI)670,670,660
  660 DT= SIGN(.9*TI,DT)
!
  670 TII=TI+DT
      FI=F
      IP=1
      TC=TII
      GO TO 510
!
  680 DFDT=(F-FI)/DT
!
      IF(DFDT)700,690,700
!
  690 TI=TII
      GO TO 500
!
  700 DT=-FI/DFDT
  710 IF( ABS(DT)-.9*TI)730,720,720      ! label 710 not used
!
  720 DT= SIGN(.9*TI,DT)
!
  730 TI=TI+DT
      GO TO 500
!*****SECOND ATTEMPT AT CONVERGING ON A ROOT USING A SEARCHING METHOD
  800 TU = TSAV1
      IF(TSAV1.EQ.0.0) TU = 27000.0
      TL = TSAV2
      GM = (SQRT(5.0)-1.0)/2.0
  805 KC = KC +1
      DELRT = (TU-TL)*GM
      TC = TL + DELRT
      CALL TNP (PC,TC,HC)
      F = HC - H
      IF(ABS(F/H).GT.EPS) GO TO 810
      GO TO 530
  810 IF(F) 815,815,820
  815 TL = TC
      GO TO 825
  820 TU = TC
  825 IF(KC.LT.45) GO TO 805
      FOH = F/H
      IKCT = IC + KC
      WRITE(6,850) IKCT,TC,F,FOH,P,H,SOR,T,EPS
  850 FORMAT(/,' PASSES IN PHGETT = ',I3, ' TC, F, FOH = ',3ES20.8/      &
     &  ' INITIAL ARGUMENT P,H,SOR,T,EPS =',5ES14.6)
      SOR = SORC
      T = TC
 1000 ICNTR = IC + KC
 9999 RETURN
      END Subroutine Phgett

!
SUBROUTINE TEMPE (RM,ET,GAMF,ZF)
! ---------------------------------------------------------------------------
IMPLICIT NONE
  REAL,INTENT(IN):: rm
  REAL,INTENT(OUT):: et
  REAL,INTENT(IN):: gamf
  REAL,INTENT(IN):: zf

  REAL:: expo,part1
!----------------------------------------------------------------------------
  expo   = 0.5 * (gamf + 1.0) / (gamf - 1.0)
  part1  = (2.0 + (gamf - 1.0) * (rm ** 2.0)) / (gamf + 1.0)
  et     = (SQRT(zf)/ rm) * (part1 ** expo)
  RETURN
END Subroutine Tempe   ! ----------------------------------------------------

!
      SUBROUTINE INTERP (RM,SI,SISTR)
  REAL,INTENT(IN):: rm
  REAL,INTENT(OUT):: si,sistr

      DIMENSION FSM(51),S1(51),S2(51)
      DATA S1/    1.286,     1.304,     1.358,     1.447,               &
     &            1.573,     1.734,     1.930,     2.163,               &
     &            2.430,     2.733,     3.071,     3.444,               &
     &            3.852,     4.294,     4.770,     5.284,               &
     &            5.831,     6.413,     7.029,     7.680,               &
     &            8.365,     9.085,     9.839,    10.627,               &
     &           11.451,    12.308,    13.200,    14.126,               &
     &           15.091,    16.082,    17.112,    18.175,               &
     &           19.273,    20.400,    21.573,    22.775,               &
     &           24.010,    25.279,    26.584,    27.922,               &
     &           29.295,    30.702,    32.144,    33.620,               &
     &           35.130,    36.675,    38.254,    39.867,               &
     &           41.515,    43.197,    44.914/
      DATA S2/   10.286,    10.314,    10.399,    10.540,               &
     &           10.735,    10.984,    11.282,    11.634,               &
     &           12.031,    12.475,    12.963,    13.494,               &
     &           14.067,    14.680,    15.328,    16.023,               &
     &           16.752,    17.517,    18.319,    19.157,               &
     &           20.030,    20.938,    21.881,    22.858,               &
     &           23.869,    24.914,    25.993,    27.105,               &
     &           28.259,    29.431,    30.644,    31.889,               &
     &           33.168,    34.471,    35.825,    37.205,               &
     &           38.615,    40.059,    41.536,    43.046,               &
     &           44.589,    46.165,    47.774,    49.416,               &
     &           51.091,    52.799,    54.539,    56.313,               &
     &           58.119,    59.959,    61.832/
!----------------------------------------------------------------------------
      A = 0.0
      DO I=1,51
      FSM(I) = A
      A = A + 0.2
      END DO
!
      DO I=1,51
      IF (RM.LE.FSM(I)) GO TO 25
      END DO

      I=50
   25 F = (RM - FSM(I)) / (FSM(I+1) - FSM(I))
      SI = F * (S1(I+1) - S1(I)) + S1(I)
      SISTR = F * (S2(I+1) - S2(I)) + S2(I)
      RETURN
      END Subroutine Interp
!
      SUBROUTINE TNP (P,T,H)
  REAL,INTENT(IN):: p,t
  REAL,INTENT(OUT):: h

!     ENTER WITH P IN CR(1) AND T IN CR(2)
!     EXIT WITH SOR IN CR(3), H IN CR(4), Z IN CR(5), RHO IN CR(6)

!!!      REAL LPF8,LPF10,LPF2,LPF24

      COMMON CR
      DIMENSION CR(200),PF(40)
      EQUIVALENCE                                                       &
     &   (CR(3 ),SOR ) ,                  (CR(5 ),Z   ) , (CR(6 ),RHO ) &
     & , (CR(7 ),GAM ) , (CR(8 ),VIS )  , (CR(9 ),PRTH) , (CR(10),A   ) &
     & , (CR(11),CON ) , (CR(12),CONTH) , (CR(13),CP  ) , (CR(14),CV  ) &
     & , (CR(15),PR  ) , (CR(16),CI  )  , (CR(17),SQCI) , (CR(18),CII ) &
     & , (CR(19),CIII) , (CR(20),SIG )  , (CR(21),TLN ) , (CR(22),ZE  ) &
     & , (CR(23),ZO  ) , (CR(24),ZOI )  , (CR(25),ZO2 ) , (CR(26),ZOPL) &
     & , (CR(27),ZN  ) , (CR(28),ZN2 )  , (CR(29),ZNII) , (CR(30),ZNPL) &
     & , (CR(31),PF(1)), (CR(68),DP5 )  , (CR(69),DP25) , (CR(70),DP26) &
     & , (CR(92),IERR)
!
      IF (p>0.0 .AND. t>0.0) GO TO 5
!      IF (P) 2,2,1
!    1 IF (T) 2,2,5
      WRITE(6,3)P,T
    3 FORMAT(/,' SUBROUTINE TNP ABORT EXIT, P=',ES13.5,'  T=',ES13.5)
      IERR = 1
      GO TO 9999

    5 DEL=P/2116.2169
      DLN= LOG(DEL)
      TLN= LOG(T)
      CPF11=.0671+.02584*TLN
      ARG =     8.1-55500./T-DLN*.5
      IF ( ABS(ARG) - 88.028) 20,25,25
   20 SQCI =  EXP(ARG)
      IF(SQCI-1.0E-20)25,30,30
   25 SQCI = 0.
      Z = 1.
      ZO = 0.
      ZN = 0.
      ZOPL = 0.
      ZNPL = 0.
      ZE = 0.
      CI = 0.
      CII = 0.
      CIII = 0.
      ZNII = 0.
      ZO2 = .2095
      ZN2 = .7808
      IF (300. - T) 26,27,27
   26 PF(4) = PF4 (T)
      PF(6) = PF6 (T)
      PF(8) = LPF8 (TLN)
      PF(10) = LPF10 (TLN)
      SOR = .520655 + .2095 * PF(8) + .7808 * PF(10) + CPF11 - DLN
      H = .0686 * T * (.2095 * PF(4) + .7808 * PF(6) + .024 )
      GO TO 230
   27 H = .24 * T
      SOR = 3.49 * TLN + 2. - DLN
      GO TO 230
!
   30 IF(SQCI-.244949E-1)40,40,50
!
   40 ZOI=.4575*SQCI
      GO TO 60
   50 CI=SQCI**2
      ZOI=2.*( SQRT((.3657*CI+.8378)*CI)-.3953*CI)/(4.+CI)
      IF(ZOI-.22)60,60,70
!
   60 Z=1.+.5*ZOI
      ZO=ZOI
      ZO2=.2095-.5*ZO
      ZN2=.7808
      ZN=0.
      ZOPL=0.
      ZNPL=0.
      ZE=0.
      CII=0.
      CIII=0.
      ZNII = 0.
      PF(3)=PF3 (T)
      PF(4)=PF4 (T)
      PF(6)=PF6 (T)
      PF(7)=3.42+2.54*TLN
      PF(8)=LPF8 (TLN)
      PF(10)=LPF10 (TLN)
      SOR=ZO*(PF(7)- ALOG(ZO))+ZO2*(PF(8)- ALOG(ZO2))+.7808*PF(10)+CPF11&
     &+.193198-Z* ALOG(DEL/Z)
      H=3665.*ZO+.0686*T*(ZO*PF(3)+ZO2*PF(4)+.7808*PF(6)+.024)
      GO TO 230
!
   70 PF(2)=LPF2 (TLN)
      CII= EXP(PF(2)-204000./T)/DEL
      IF(CII-.0001)80,80,90
!
   80 ZNII= SQRT(.7808*(1.+.5*ZOI)*CII)
      GO TO 100
   90 ZNII=2.*( SQRT((.9903*CII+3.778)*CII)-.2143*CII)/(4.+CII)
      IF(ZNII-.7)100,100,105
!
  100 ZO2=.2095-.5*ZOI
      PF(4)=PF4 (T)
      PF(8)=LPF8 (TLN)
      CIII=0.
      ZEII=9.55E-6* SQRT(ZOI*ZNII)* EXP(.784*TLN-28500./T)
      SZO2=ZO2*(PF(8)- ALOG(ZO2))
      HZO2=ZO2*PF(4)
      GO TO 130
!
  105 IF(ZNII-1.5616)120,120,110
  110 ZNII=1.5616
  120 PF(24)=LPF24 (TLN)
      ZO2=0.
      CIII= EXP(PF(24)-300000./T)/DEL
      ZEIII=1.99* SQRT(CIII/(1.+CIII))
      ZEII=9.55E-6* SQRT(ZOI*ZNII)* EXP(.784*TLN-28500./T)
      IF(ZEII-ZEIII)140,125,125
!
  125 SZO2=0.
      HZO2=0.
  130 PF(3)=PF3 (T)
      ZN2=.7808-.5*ZNII
      ZOPL=0.
      ZNPL=0.
      DP5=DELPF5 (T)
      ZO=ZOI
      ZE=ZEII
      Z=1.+.5*(ZOI+ZNII)+ZEII
      PF(5)=2.5+DP5
      ZN=ZNII
      PF(7)=3.42+2.54*TLN
      PF(9)=3.*TLN-1.5
      PF(10)=LPF10 (TLN)
      PF(33)=2.5*TLN-13.21
      PF(6)=PF6 (T)
      SOR=ZOI*(PF(7)- ALOG(ZOI))+SZO2+ZNII*(PF(9)- ALOG(ZNII))+ZN2*(PF(1&
     &0)- ALOG(ZN2))+ZEII*(PF(33)- ALOG(ZEII))+CPF11     -Z* ALOG(DEL/Z)
      H=3665.*ZOI+6993.*ZNII+14600.*ZEII+.0686*T*(ZOI*PF(3)+HZO2+ZNII*PF&
     &(5)+ZN2*PF(6)+2.5*ZEII+.024)
      GO TO 230
!
  140 PF(37)=PF37 (T)
      CIIIPR= EXP(PF(37)-19500./T)
      F1=1.-CIIIPR
      F2=CIIIPR*ZNII+.419+ZEIII*F1
!
      IF(CIIIPR-.98)170,160,150
  150 IF(CIIIPR-1.02)160,160,170
!
  160 ZOPL=.419/F2*ZEIII
      GO TO 180
!
  170 ZOPL=.5/F1*(F2- SQRT(F2**2-1.676*ZEIII*F1))
  180 IF(ZOPL-ZOI)200,200,190
  190 ZOPL=ZOI
  200 ZNPL=ZEIII-ZOPL
      IF(ZNPL-ZNII)220,220,210
  210 ZNPL=ZNII
  220 ZN2=.7808-.5*ZNII
      ZO=ZOI-ZOPL
      ZN=ZNII-ZNPL
      Z=1.+.5*(ZNII+ZOI)+ZEIII
      ZE=ZEIII
      PF(3)=PF3 (T)
      DP5=DELPF5 (T)
      PF(5)=2.5+DP5
      PF(6)=PF6 (T)
      PF(7)=3.42+2.54*TLN
      PF(9)=3.*TLN-1.5
      PF(10)=LPF10 (TLN)
      DP25=DEPF25 (T)
      PF(25)=2.5+DP25
      DP26=DEPF26 (T)
      PF(26)=2.5+DP26
      PF(31)=2.87*TLN-.47
      PF(32)=1.5+2.73*TLN
      PF(33)=2.5*TLN-13.21
      SOR = CPF11 -Z*ALOG(DEL/Z)
      IF (ZO    .GT. 0.)  SOR = SOR + ZO*(PF(7) - ALOG(ZO))
      IF (ZN    .GT. 0.)  SOR = SOR + ZN*(PF(9) - ALOG(ZN))
      IF (ZN2   .GT. 0.)  SOR = SOR + ZN2*(PF(10) - ALOG(ZN2))
      IF (ZOPL  .GT. 0.)  SOR = SOR + ZOPL*(PF(31) - ALOG(ZOPL))
      IF (ZNPL  .GT. 0.)  SOR = SOR + ZNPL*(PF(32)-ALOG(ZNPL))
      IF (ZEIII .GT. 0.)  SOR = SOR + ZEIII*(PF(33) - ALOG(ZEIII))
      H=3665.*ZOI+6993.*ZNII+20600.*ZEIII+.0686*T*(ZO*PF(3)+ZN*PF(5)+ZN2&
     &*PF(6)+ZOPL*PF(25)+ZNPL*PF(26)+2.5*ZEIII+.024)
!
  230 SIG=491.68/Z*DEL/T
      RHO=P/Z/T/1716.48272
 9999 RETURN
      END Subroutine Tnp

!+
      SUBROUTINE TVSPR (P,H,T,VIS,PR)
! PURPOSE - Calculate miscellaneous (remaining) air properties, mu, k,
!     c(p)/r, c(v)/r, gamma, a, pr
! NOTE - In addition to computing vis and pr, numerous elements of the
!  vector CR in blank common are changed
  REAL,INTENT(IN):: p    ! apparently not used ??? 
  REAL,INTENT(IN):: h    ! apparently not used ???
  REAL,INTENT(IN):: t
  REAL,INTENT(OUT):: vis
  REAL,INTENT(OUT):: pr

      COMMON CR
      DIMENSION CR(200),PF(40)
      EQUIVALENCE                                                       &
     &   (CR(3 ),SOR ) ,                  (CR(5 ),Z   ) , (CR(6 ),RHO ) &
     & , (CR(7 ),GAM ) ,                  (CR(9 ),PRTH) , (CR(10),A   ) &
     & , (CR(11),CON ) , (CR(12),CONTH) , (CR(13),CP  ) , (CR(14),CV  ) &
     & ,                 (CR(16),CI  )  , (CR(17),SQCI) , (CR(18),CII ) &
     & , (CR(19),CIII) , (CR(20),SIG )  , (CR(21),TLN ) , (CR(22),ZE  ) &
     & , (CR(23),ZO  ) , (CR(24),ZOI )  , (CR(25),ZO2 ) , (CR(26),ZOPL) &
     & , (CR(27),ZN  ) , (CR(28),ZN2 )  , (CR(29),ZNII) , (CR(30),ZNPL) &
     & , (CR(31),PF(1)), (CR(68),DP5 )  , (CR(69),DP25) , (CR(70),DP26) &
     & , (CR(72),ZOM ),(CR(89),INST)

      REAL:: delfe = 0.0                           ! added by RLC

      CONR= SQRT(T)/(1.+202./T)
      VISR=7.33E-7*CONR
      CONR=2.39E-7*CONR
      PF(21) = 2.15 - .179 * TLN
      F111T=111000./T
!
      CPII=0.
      CPIII=0.
      CVII=0.
      CVIII=0.
      CONII=0.
      CONIII=0.
      IF (SQCI - .0244949) 1,1,12
!
    1 IF (T - 300.) 2,2,3
    2 CVOR = 2.5
      CPOR = 3.5
      VISRAT = 1.
      CONTH = 1.
      CONRAT = 1.
      GAM = CPOR / CVOR
      GO TO 290
!
    3 F1=5600./T
      F2= EXP(F1)
      F3=(F1/(F2-1.))**2*F2
      CVTH=2.5+F3
      CPTH=CVTH+1.
      CVI=SQCI/T
      COT=13.6E8/T
      CPI=CVI*(COT+12700.)
      CVI=CVI*(COT-12230.)
      CPOR=CPTH+CPI
      CVOR=CVTH+CVI
      GAM=CPOR/CVOR
!
      IF (T - 1800.) 5,5,6
!
    5 PF(36) = 1.
      GO TO 7
    6 PF(36) = 2.05 - .14 * TLN
    7 VISRAT = 1.0 / PF(36)
      CONTH = (1.0+(.209*F3))/PF(36)
      CONI=.0358*SQCI/PF(21)*(F111T+1.)**2
      CONRAT=CONTH+CONI
      GO TO 290
!
   12 ZA=ZO+ZN
      IF(ZA)14,14,16
   14 ZM=0.
      GO TO 20
   16 ZM = 1.0 - (ZOI + ZNII) / 2.0
   20 PF(12)=PF12 (T)
!
      IF (T-1800.) 21,21,23
!
   21 PF(36) = 1.0
      GO TO 24
   23 PF(36) = 2.05 - .14 * TLN
   24 PF(4)=PF4 (T)
      PF(13)=PF13 (T)
      PF(38)=DELPF5 (T)
      PF(6)=PF6 (T)
      PF(14)=PF14 (T,PF(38))
      PF(3)=PF3 (T)
      PF(5)=2.5+PF(38)
      PF(15)=PF15 (T)
      PF(18)=2.*PF(3)-PF(4)
      PF(19)=2.*PF(5)-PF(6)
      PF(20) = .8 - 1.6 * T / 100000.0
      PF(22) = 1.0/3.0 + 1600.0/(T+2400.0)
      PF(23)=1.852-.1516*TLN
      PF(39)=DEPF25 (T)
      PF(25)=2.5+PF(39)
      PF(40)=DEPF26 (T)
      PF(26)=2.5+PF(40)
      PF(27)=PF27 (T,PF(39))
      PF(28)=PF28 (T)
      PF(30)=.21*(PF(25)-PF(3))+.79*(PF(26)-PF(5))
      PF(34) = 40.0* T **(-.39)
      IF (ZE) 300,300,301
!     SIMULATED INFINITY VALUE IS NOT USED
  300 PF(35) =-.9E38
      GO TO 302
  301 PF(35)=3.6E9*( ALOG(T**3/SIG/ZE)*.5-14.17)/(T**2+(202.*T))
  302 CVTH=ZO*PF(12)+ZO2*PF(13)+ZN*PF(14)+ZN2*PF(15)+ZOPL*PF(27)+ZNPL*PF&
     &(28)+1.5*ZE+.014
      CPTH=CVTH+Z
      FI=CI/(4.+CI)*((.4189+.3657*CI)/ SQRT((.3657*CI+.8378)*CI)-.3953-Z&
     &OI/2.)
      CPI=F111T**2*FI
      CG1=FI/(FI+Z)*(F111T-1.)
      CVI=Z*CG1*(F111T-1.)
      CG4=FI*F111T/Z
      ZAPL=ZOPL+ZNPL
!
      IF(ZM)22,22,25
!
   22 CU1=0.
      CK1=0.
      FKU5=0.
      FMM= 40.
      FKU11=0.
      GO TO 26
!
   25 FMM =(28.967-(14.01*ZNII)-(16.0*ZOI))/ ZM
      SQMM= SQRT(FMM)
      FU1=.2625*ZM*SQMM
      FKU2=1.414*ZM
      FKU5=PF(20)*ZM
      FK1=3.6*ZM+1.6*(ZO2*PF(13)+ZN2*PF(15))
      FKU11=1.264*ZM
!
   26 IF(ZA)30,30,40
!
   30 CU2=0.
      CK2=0.
      FMA=0.
      FKU8=0.
!
      IF(ZM)50,50,35
!
   35 CU1=FU1/FKU2
      CK1=FK1/SQMM/FKU2
      GO TO 50
!
   40 FMA=(16.*ZO+14.*ZN)/ZA
      SQMA= SQRT(FMA)
      FKU8=PF(22)*ZA
      FKU11=ZA+FKU11
!
      IF(ZM)44,44,42
!
   42 FKU5=FKU5* SQRT(1.+FMA/FMM)
      CU1 = FKU2*PF(36)+(ZA*PF(20)+ZAPL*PF(34))* SQRT(1.0+FMM/FMA)
      CK1 = (FK1+.0233) / (CU1*SQMM)
      CU1 = FU1/CU1
   44 CU2 = FKU5+1.414*(ZA*PF(22)+ZAPL*PF(34))
      CK2 =(3.6*ZA+1.6*(ZO*PF(12)+ZN*PF(14)))/(CU2*SQMA)
      CU2 = .2625*ZA*SQMA/CU2
   50 IF (ZE) 55,55,70
   55 CK3 = 0.
   60 CU3 = 0.
      GO TO 130
   70 CK3 = 255.9 * ZE/(.31*(ZA+ZM) + ZE*PF(35))
   80 IF (ZA) 60,60,90                            ! label 80 not used
   90 CU3 = ZAPL/(1.414 * ZA * PF(34) + PF(35)*ZAPL)
!
  130 VISRAT=CU1+CU2+CU3
      CONTH =CK1+CK2+CK3
      ZOZO2=ZO*ZO2
      ZNZN2=ZN*ZN2
      ZAZAPL=ZA*ZAPL
      CONI=0.
      F204TP=204000./T+PF(19)
      F3TP=300000./T+PF(30)
!
      IF(ZOZO2)160,160,150
  150 CONI=.179*ZOZO2/(PF(21)*(.101+1.8*ZO2)+.552*ZO)*F111T**2

  160 IF(ZNZN2)190,190,180
  180 CONII=.179*ZNZN2/(PF(21)*(1.41+.242*ZN)+.838*PF(23)*ZN2)*F204TP**2

  190 IF(ZAZAPL)220,220,210
  210 CONIII = 10.3385 * ZAZAPL / (115.*PF(34)) * (F3TP + 2.5)**2

  220 IF(CII)230,230,240
  230 CG2=0.
      CG5=0.
      GO TO 250

  240 FII=CII/(4.+CII)*((1.887+.9903*CII)/ SQRT((.9903*CII+3.774)*CII)-.&
     &2143-ZNII/2.)
      CPII=FII*F204TP
      CVII=F204TP-1.
      CG2=FII/(FII+Z)*CVII
      CVII=Z*CG2*CVII
      CG5=CPII/Z
      CPII=CPII*F204TP
!
  250 IF(CIII)260,260,270
  260 CG3=0.
      CG6=0.
      GO TO 280
  270 SQC9=.99* SQRT(CIII)
      SQC3=1.+CIII
      SQC3= SQRT(SQC3)*SQC3
      CG6=SQC9/SQC3*(F3TP+2.5)
      CPIII=CG6*(F3TP+2.5)
      CG6=CG6/Z
      CG3=SQC9/(Z*SQC3+SQC9)*(F3TP+1.5)
      CVIII=Z*CG3*(F3TP+1.5)
!
  280 CPOR=CPTH+CPI+CPII+CPIII
      CVOR=CVTH+CVI+CVII+CVIII
      CONRAT=CONTH+CONI+CONII+CONIII
      GAM=CPOR/CVOR*(1.+CG1+CG2+CG3)/(1.+CG4+CG5+CG6)
  290 VIS=VISRAT*VISR
      CON=CONRAT*CONR
      AA=1716.*GAM*T*Z
      A= SQRT(AA)
      CP=.0686*CPOR
      CV=.0686*CVOR
      PR=.210526*CPOR/CONRAT*VISRAT
      PRTH=.210526*CPTH/CONTH*VISRAT
!               DIMENSIONED FORM
      CONTH = CONTH * CONR
      IF(T .LT. 180.) RETURN
      TK = T/1.8
      DELFA = -3.030993 * (1.-ALOG(TK))
      DELFB = -1.309476E-4*TK
      DELFC = +1.725756E-9*TK**2
      DELFD = + .4235962E-13*TK**3
      DELFS = DELFA+DELFB+DELFC+DELFD+DELFE-14.75022   ! delfe never computed
      DELTAF = -1.9865*TK*DELFS-3.22788E+4
      EXPK = EXP(DELTAF/(1.9865*TK))
      ZCPL = 0.0
      SUMPL = ZNPL + ZOPL + ZCPL
      CONRHO = 0.5154/28.966*RHO
      ZOM = 82.05*TK*CONRHO*EXPK*SUMPL*ZO
      RETURN
      END Subroutine Tvspr



      FUNCTION Lpf2(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f
      f  = 9.17 * EXP (.06 * t)
      RETURN
      END Function Lpf2

      FUNCTION Pf3(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f
      REAL:: a,b
      A = EXP(-41000./T)
      B = EXP(-410./T)
      f = (190000.* A + 1500.* B) / (T *(5.*(A +1.)+ 3.* B))+ 2.5
      RETURN
      END Function Pf3

      FUNCTION Pf4(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

      f=3.5+((4090./(EXP(4090./T)-1.)+20000.*EXP(-21500./T))/T)
      RETURN
      END Function Pf4

      FUNCTION DelPf5(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

      IF (T < 4500.) THEN
        f    = 0.
      ELSE
        f    = 120000. * (EXP (-49000./ T) / T)
      END IF

      RETURN
      END Function Delpf5

      FUNCTION Pf6(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f
      REAL:: a
      A = 6100. / T
      f = 3.5 + A /(EXP(A)-1.)
      RETURN
      END Function Pf6

      FUNCTION Lpf8(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

      IF ( T   >= 7.2) THEN
        f  = 4.5 *( T -1.)
      ELSE
        f = 3.5 * T + 2.7
      END IF
      RETURN
      END Function Lpf8

      FUNCTION Lpf10(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

      IF (T>=7.32) THEN
        f = 4.325 *T-5.
      ELSE
        f  = 3.5 * T + 1.04
      END IF
      RETURN
      END Function Lpf10

      FUNCTION Pf12(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

!!!      IF (T.GE. 4100.) GO TO 25
!!!      PFS  = 1.5
!!!      GO TO 50
!!!   30 PFS  = 1.4 * EXP (T/60000.)
!!!      GO TO 50
!!!   25 IF(T.LT. 14400.) GO TO 30
!!!      PFS  = 1.78
!!!      GO TO 50

      IF (t < 4100.0) THEN
        f=1.5
      ELSEIF (t < 14400.0) THEN
        f=1.4*EXP(t/60000.0)
      ELSE
        f=1.78
      END IF

      RETURN
      END Function Pf12

      FUNCTION Pf13(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f
      REAL:: a

      A = EXP(4090. /T)
      f=2.5+((4090.0/(1.-A))**2*A + 4.3*10.E7*EXP(-21500./T))/ T**2
      RETURN
      END FUNCTION Pf13

      FUNCTION Pf15(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f
      REAL:: a

      A = EXP ( 6100. /T)
      f  =2.5+(6100.   /(T*(A-1.))) **2 * A
      RETURN
      END Function Pf15

      FUNCTION Lpf24(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

      f = .62 *EXP(.275*T)
      RETURN
      END Function Lpf24

      FUNCTION Depf25(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

      IF (T.LE. 7000.) THEN
        f=0.0
      ELSE
        f   = (190000.*EXP (-70000./T))/ T
      END IF
      RETURN
      END Function Depf25

      FUNCTION Depf26(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

      IF (T < 8100.) THEN
        f    = .04
      ELSE
        f   = (T/90000.) - .05
      END IF

      RETURN
      END Function Depf26

      FUNCTION Pf28(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

      IF (T > 5400.) THEN
        f  = 1.5 + .35 * (1.- EXP(-.00015 * T))
      ELSE
        f = 1.5
      END IF

      RETURN
      END Function Pf28

      FUNCTION Pf37(t) RESULT(f)
      REAL,INTENT(IN):: t
      REAL:: f

      f  = 1.58 - 4.2 * EXP(-70000. /T)
      RETURN
      END Function Pf37

      FUNCTION Pf27(t,d) RESULT(f)
      REAL,INTENT(IN):: t,d
      REAL:: f
      REAL:: a,b

      A = EXP(-69500. /T)
      B = EXP(-104900. /T)
      f=1.5+(((241.*A)+330.*B)/(2.+3.*B+5.*A))*10.E7/T**2-D**2
      RETURN
      END Function Pf27

      FUNCTION Pf14(t,d) RESULT(f)
      REAL,INTENT(IN):: t,d
      REAL:: f
      REAL:: a,b

      a = EXP(-49900./t)
      b = EXP(-74800./t)

      f = (  1.E+8 / T**2) * (124.1 * A + 167.5 * B ) / ( 5.* A + 3.    &
     & * B + 2.) -(D)**2 + 1.5
      RETURN
      END Function Pf14


      END Module ColdArcProcedures   !=======================================






!
!     ***************************************************************
!     *                                                             *
!     *          B. A. MILLER'S  PLASMA COLD ARC PROGRAM            *
!     *                                                             *
!     *          MARCH 13, 1981                                     *
!     *          PROGRAMMER : JAB                                   *
!     *                                                             *
!     ***************************************************************
!
!          PTOT      TOTAL PRESSURE                  PSF
!          HTOT      TOTAL ENTHALPY                  BTU/LB
!          E         A/A*
!          EMIS      EMISSIVITY
!          WDOT      FLOW RATE                       LB/SEC
!          AEX       AREA AT NOZZLE EXIT             FT2
!          TW        SURFACE TEMP OF TEST ARTICLE    DEGREES R
!          PL2       SURFACE PRES OF TEST ARTICLE    PSF
!          ETA       CATALYTICITY



      PROGRAM ColdArc
      USE ColdArcProcedures

      REAL NOZMOM,KEBL,LEQUIV,MRAT53,MRAT54,MRAT55,MRAT56,MRAT57,       &
     &     MRAT58,MRAT59,MOMEN,KESL,INTESL,MOMSL1,MOMSL2,MOMSL

  CHARACTER(LEN=80):: fileName 
      INTEGER:: errCode
!
      DATA CON467,CON530/0.467E-12,530.0/
!----------------------------------------------------------------------------
  WRITE(*,*) 'ColdArc - dissociated air flow effects during plasma arc testing'
  DO
    WRITE(*,*) 'Enter the name of the input file:'
    READ(*,'(A)') fileName
    IF (LEN_TRIM(fileName)==0) STOP
    OPEN(UNIT=5, FILE=fileName, STATUS='OLD', IOSTAT=errCode, ACTION='READ')
    IF (errCode==0) EXIT
    WRITE(*,*) 'Unable to open this file. Try again.'
  END DO
  OPEN(UNIT=7,FILE='coldarc.out',STATUS='REPLACE',ACTION='WRITE')

      NCASE = 1
    1 READ(5,*,IOSTAT=errCode) PTOT,HTOT,E,EMIS,WDOT,AEX,TW,PL2,ETA,PSI
      IF (errCode < 0) GO TO 999
!!!      IF (EOF(5)) 999,2
!
!   INPUT INITIAL VALUE OF TTOT
!
    2 TTOT = 8300.0
!
      CALL PHGETT (PTOT,HTOT,TTOT)
      ZF = (HTOT + 7380.0 - .2744 * TTOT) / (7380.0 + 0.0343 * TTOT)
      IF (ZF.LE.1.21) GO TO 5
!
      ZF = (HTOT - 0.2744 * TTOT + 15365.8) / (13980.0 + 0.0343 * TTOT)
    5 GAMF = (ZF + 8.0) / (8.0 - ZF)
!
!   FROZEN CP WITH VIB. DEGREE OF FREEDOM FULLY ACTIVE
!
      CPF = 0.0686 * (0.5 * ZF + 4.0)
!
      IF (ZF.GT.1.21) GO TO 10
      HD = 7380.0 * (ZF - 1.0)
      GO TO 20
!
   10 IF (ZF.GT.2.0) GO TO 15
      HD = 1550.0 + 13980.0 * (ZF - 1.21)
      GO TO 20
!
!  IONIZATION
!
   15 HD = 12595.0
!
!  FIND THE LOCAL MACH NO. AT NOZZLE EXIT
!
   20 RM = 0.5
   30 CALL TEMPE (RM,ET,GAMF,ZF)
      IF (E-ET) 45,40,35
!
   35 RM = RM + 0.5
      GO TO 30
!
   40 RML1 = RM
      GO TO 60
!
   45 RMOLD = RM - 0.5
      RMNEW = RM
   47 RM = (RMOLD + RMNEW) / 2.0
      CALL TEMPE (RM,ET,GAMF,ZF)
      IF ((ABS(E-ET)).LT.0.01) GO TO 60
      IF (E-ET) 50,40,55
!
   50 RMNEW = RM
      GO TO 47
!
   55 RMOLD = RM
      GO TO 47
!
   60 RML1 = RM
!
      PART1 = (1.28586 * WDOT) / (RML1 * AEX)
      PART2 = ((HTOT - HD) * ZF/GAMF) ** 0.5
      PART3 = (CPF + ZF * GAMF * 0.0342 * (RML1**2.0)) ** 0.5
      PL1 = PART1 * (PART2 / PART3)
!
      TL1 = (0.6048 * GAMF / ZF) * ((AEX * PL1 * RML1 / WDOT) ** 2.0)
      HL1 = HTOT - 0.0342 * GAMF * ZF * TL1 * (RML1 ** 2.0)
      UL1 = 223.84 * ((HTOT - HL1) ** 0.5)
      RHOL1 = PL1 / (1716.0 * ZF * TL1)
!
      PART1 = (PL2 - PL1) ** 2.0
      PART2 = 2.0 * (PL2 - PL1) -PL2
      PART3 = (UL1 ** 2.0) * RHOL1
      SIN2TH = (PART1 / PART2) / PART3
      SINTH = SQRT(SIN2TH)
      SHOCK = 57.3 * ASIN(SINTH)
!
      PART1 = TL1 * (PL2/PL1)
      PART2 = 1.0 - (PL2/PL1 - 1.0) / (RML1**2.0 * SIN2TH * GAMF)
      TL2 = PART1 * PART2
!
      PART1 = 0.03425 * GAMF * TL1 * ZF * (RML1 ** 2.0) * SIN2TH
      PART2 = 1.0 - (((TL2/TL1) / (PL2/PL1)) ** 2.0)
      HL2 = PART1 * PART2 + HL1
      HL2 = CPF * TL2 + HD
!
      UL2 =223.84 * (HTOT - HL2) ** 0.5
      RHOL2 = PL2 / (1716.0 * ZF * TL2)
      RML2 = UL2 / ((GAMF * 1716.0 * TL2 * ZF) ** 0.5)
!
      CALL TVSPR (PTOT,HTOT,TTOT,VISC,PRF)
!
      UC = VISC / 32.174
      OF = ((TL2/TTOT) ** 1.5) * ((TTOT + 202.0)/(TL2 + 202.0))
      UF = UC * OF
!
      CALL TNP (PL2,TW,HWEQ)
!
      HWNE = CPF * TW + ETA * HD
!
      HAWEQ = (UL2 ** 2.0) * (PRF ** 0.33) / 50103.0 + HL2
!
!   ENTHALPY DIFFERENCES
!
      DELHEQ = HAWEQ - HWEQ
      DELHNE = HAWEQ - HWNE
!
!   REYNOLDS ANALOGY FACTOR
!
      FAC = 32.174 / (UL2 * (PRF ** 0.67))
!
      MOMEN = RHOL2 * (UL2 ** 2.0)
!
!   TURBULENT BOUNDARY LAYER THICKNESS, FLAT PLATE
!
      REX = RHOL2 * UL2 * PSI / UF
      DELTA = 0.37 * PSI / (REX ** 0.2)
!
!   ECKERT'S REFERENCE ENTHALPY
!
      HSTARF = HL2 + 0.5 * (HWNE - HL2) + 0.22 * (HAWEQ - HL2)
      TSTARF = (HSTARF - HD) / CPF
!
!   ECKERT'S TURBULENT, FLAT PLATE, CF/2
!
      PART1 = ((HL2 - HD)/(HSTARF - HD)) ** 0.8
      PART2 = (TSTARF / TL2) ** 0.3
      PART3 = ((TL2 + 202.0)/(TSTARF + 202.0)) ** 0.2
      CFD2EC = 0.0296 /  (REX ** 0.2) * PART1 * PART2 * PART3
!
!   BARTZ'S TURBULENT, WITH PRESSURE GRADIENT, PRIMARILY FOR NOZZLES
!
      PART1 = (UC * TSTARF / (UL2 * DELTA * RHOL2 * TL2)) ** 0.25
      PART2 = (2.0 * TL2 / (TW + TL2)) ** 0.6
      PART3 = (TL2 / TTOT) ** 0.15
      CFD2BT = 0.0228 * PART1 * PART2 * PART3
!
!   H. SCHLICHTING'S TURBULENT, FLAT PLATE, CF/2
!   REPLACED BY SPALDING ) CHI    4/30/81
!
      A = TW/TL2
      G = (GAMF - 1.0) * RML2**2.0
      C1 = -0.5 * PRF ** (0.333 * G)
      B = 1.0 - C1 - A
      FC = (1.0 /(-C1)** 0.5 * (ASIN(-(2.0 * C1 + B)/                   &
     &      SQRT(ABS(B) ** 2.0 - 4.0 * A * C1)) - ASIN(-B /             &
     &      SQRT(ABS(B) ** 2.0 - 4.0 * A * C1))))
      IF (FC.NE.0.0) FC = 1.0 / FC**2.0
      IF (FC.EQ.0.0) FC = 1.0E+35
      TAW = (HAWEQ - HD) / CPF
      FRO = A ** (-0.702) * (TAW / TW) ** 0.772
      FRX = FRO / FC
      C = ALOG (FRX * REX)
      CFEXP = 9.280864 + C * (-4.734 + C * (6.685E-01 + C *             &
     &       (-4.1876E-02 + C * (-5.505E-04 + C * (2.837E-04 +          &
     &        C * (-2.1249E-05 + C * (8.0167E-07 + C * (-1.59E-08 +     &
     &        C * (1.323E-10)))))))))
      CFD2SH = 0.5 * (EXP(CFEXP) / FC)
!
!   TABLES IN INTERP DONE BY WILSON, TURBULENT FLAT PLATE
!
      CALL INTERP (RML2,DELSDO,DELDO)
!
!   SQUIRE ) YOUNG, TURBULENT FLAT PLATE
!
      THETA = DELTA / DELDO
      DELSTR = DELSDO * THETA
!
      CFD2SY = 0.0288 / ((ALOG10(4.075 * REX * THETA /PSI)) ** 2.0)
!
!   PRANDTL'S TURBULENT FLAT PLATE
!
      CFD2PL = 0.0128 / ((REX * THETA / PSI) ** 0.25)
!
!   FALKNER'S TURBULENT FLAT PLATE
!
      PART1 = ((2.0 * TL2) / (TL2 + TW)) ** 0.715
      CFD2FK = 0.0131 / (REX ** 0.143) * PART1
!
!   FRANKL ) VOISHEL TURBULENT FLAT PLATE
!
      PART1 = 1.0 - 1.12/ALOG10(REX)
      PART2 = (ALOG10(REX)) ** (-2.58)
      PART3 = (TL2/TW) ** 0.467
      CFD2FV = 0.236 * PART1 * PART2 * PART3
!
!  CALCULATE SHEAR = CF/2 * MOMENTUM
!
      TAUECK = MOMEN * CFD2EC
      TAUBTZ = MOMEN * CFD2BT
      TAUSHG = MOMEN * CFD2SH
      TAUSY  = MOMEN * CFD2SY
      TAUPTL = MOMEN * CFD2PL
      TAUFLK = MOMEN * CFD2FK
      TAUFKV = MOMEN * CFD2FV
!
!   USING SHEAR AND REYNOLD'S ANALOGY, COMPUTE HEAT FLUXES
!
      QNEECK = FAC * DELHNE * TAUECK
      QNEBTZ = FAC * DELHNE * TAUBTZ
      QNESHG = FAC * DELHNE * TAUSHG
      QNESY  = FAC * DELHNE * TAUSY
      QNEPTL = FAC * DELHNE * TAUPTL
      QNEFLK = FAC * DELHNE * TAUFLK
      QNEFKV = FAC * DELHNE * TAUFKV
!
!   HEAT FLUX INFERRED BY WALL TEMPERATURE MEASUREMENT
!
      QNERAD = CON467 * (TW ** 4.0 - CON530 ** 4.0) * EMIS
!
!   AN EFFECTIVE LENGTH USED TO COMPUTE AREAS
!
      LEQUIV = SQRT (AEX / 0.785)
!
!   BOUNDARY LAYER FLOW RATE
!
      WDOTBL = LEQUIV * (DELTA - DELSTR) * UL2 * RHOL2 * 32.174
!
!   BOUNDARY LAYER KINETIC ENERGY USING AN INTEGRAL AVERAGE OVER DELTA
!
      KEBL   = WDOTBL * ((0.83 * UL2) ** 2.0) / 50103.0
!
!   POTENTIAL ENERGY
!
      ENGINT = WDOTBL * HL2
!
!   TOTAL FLOW ENERGY PRODUCED BY ARC AT ENTRANCE TO TEST SECTION
!
      ENGNOZ = WDOT * HTOT
!
      BLMOM  = WDOTBL * UL2 * 0.83 / 32.174 + PL2 * DELTA * LEQUIV
      AREAS  = LEQUIV * PSI
!
!        SHOCK LAYER DATA
!
      PART1 = (GAMF-1.0) * (RML1**2.0) * SIN2TH + 2.0
      PART2 = 2.0 * GAMF * (RML1**2.0) * SIN2TH - (GAMF-1.0)
      SIN2DT = (PART1/PART2) / (RML2**2.0)
!
      DELANG = ATAN(ASIN(SQRT(SIN2DT)))
!
!   FLOW RATE - SL
!
!     WDOTSL = LEQUIV * (PSI * DELANG - DELTA) * RHOL2 * UL2 * 32.174
      WDOTSL = WDOT-WDOTBL
!
!   KINETIC ENERGY - SL
!
      KESL = WDOTSL * (UL2**2.0) / 50103.0
!
!   INTERNAL ENERGY - SL
!
      INTESL = HL2 * WDOTSL
      ENGSL  = INTESL + KESL
!
!   SHOCK LAYER MOMENTUM
!
      MOMSL1 = WDOTSL * UL2 / 32.174
      MOMSL2 = PL2 * LEQUIV * (PSI * DELANG - DELTA)
      MOMSL  = MOMSL1 + MOMSL2
!
!   CHECK ON SOL. ACC. W.R.T. TOTAL AVERAGE ENERGY
!
      ERAT61 = (KEBL + ENGINT + AREAS * QNEECK) / ENGNOZ + ENGSL/ENGNOZ
      ERAT62 = (KEBL + ENGINT + AREAS * QNEBTZ) / ENGNOZ + ENGSL/ENGNOZ
      ERAT63 = (KEBL + ENGINT + AREAS * QNESHG) / ENGNOZ + ENGSL/ENGNOZ
      ERAT64 = (KEBL + ENGINT + AREAS * QNESY ) / ENGNOZ + ENGSL/ENGNOZ
      ERAT65 = (KEBL + ENGINT + AREAS * QNEPTL) / ENGNOZ + ENGSL/ENGNOZ
      ERAT66 = (KEBL + ENGINT + AREAS * QNEFLK) / ENGNOZ + ENGSL/ENGNOZ
      ERAT67 = (KEBL + ENGINT + AREAS * QNEFKV) / ENGNOZ + ENGSL/ENGNOZ
!
      NOZMOM = WDOT * UL1 / 32.174 + PL1 * AEX
!
!   MOMENTUM CHECK
!
      MRAT53 = (BLMOM + AREAS * TAUECK) / NOZMOM + MOMSL/NOZMOM
      MRAT54 = (BLMOM + AREAS * TAUBTZ) / NOZMOM + MOMSL/NOZMOM
      MRAT55 = (BLMOM + AREAS * TAUSHG) / NOZMOM + MOMSL/NOZMOM
      MRAT56 = (BLMOM + AREAS * TAUSY ) / NOZMOM + MOMSL/NOZMOM
      MRAT57 = (BLMOM + AREAS * TAUPTL) / NOZMOM + MOMSL/NOZMOM
      MRAT58 = (BLMOM + AREAS * TAUFLK) / NOZMOM + MOMSL/NOZMOM
      MRAT59 = (BLMOM + AREAS * TAUFKV) / NOZMOM + MOMSL/NOZMOM
!
!        WRITE STATEMENTS WILL GO HERE
!
      WRITE(7,*) ACHAR(12)
      WRITE(7,1000)
      WRITE(7,1010) NCASE
      WRITE(7,1020) PTOT,HTOT,E,EMIS,WDOT,AEX,TW,PL2,ETA,PSI
      WRITE(7,1030)
      WRITE(7,1040) RML1,RML2,UL1,UL2,PL1,PL2,TL1,TL2,                   &
     &               RHOL1,RHOL2,HL1,HL2,HWEQ
      WRITE(7,1060) HWNE,DELTA,DELSTR,THETA,SHOCK
      WRITE(7,1070) TAUECK,TAUBTZ,TAUSHG,TAUSY,TAUPTL,TAUFLK,TAUFKV,     &
     &               QNEECK,QNEBTZ,QNESHG,QNESY,QNEPTL,QNEFLK,QNEFKV,    &
     &               MRAT53,MRAT54,MRAT55,MRAT56,MRAT57,MRAT58,MRAT59,   &
     &               ERAT61,ERAT62,ERAT63,ERAT64,ERAT65,ERAT66,ERAT67
      WRITE(7,1110) ZF,CPF,GAMF,PRF,HAWEQ,HD,TTOT
      WRITE(7,1120) WDOTBL,KEBL,ENGINT,QNERAD,WDOTSL,KESL,INTESL
      NCASE = NCASE + 1
      GO TO 1   ! keep looking for more cases until EOF

  999 STOP
 1000 FORMAT (33X,"***  PLASMA COLD ARC RESULTS  ***",33X,"PAGE 1"//)

 1010 FORMAT (39X,"INPUT DATA     CASE ",I2/)
 1020 FORMAT(1X,109("_")/1X,"!  PTOT   !   HTOT   !   A/A*   !   EMIS",  &
     &       "   !   WDOT   !    AEX   !    TW    !    PL2   !    ETA",  &
     &       "   !   PSI"/1X,"!  PSF    !  BTU/LB  !",3(10X,"!"),4X,     &
     &       "FT2   !  DEG. R  !    PSF   !",10X,"!"/1X,"!_________",    &
     &       9("!__________")/1X,"!",9X,9("!          ")/                &
     &       1X,"!",F7.0,2X,"!",F8.0,2X,"!",F7.0,3X,"!",F7.2,3X,"!",     &
     &       F7.2,3X,"!",F8.3,2X,"!",F8.0,2X,"!",F7.0,3X,"!",F8.3,2X,    &
     &       "!",F6.1/1X,"!_________",9("!__________")///)               


 1030 FORMAT (45X,"OUTPUT DATA"/4X,103("_")/4X,"!",46X,"!  PRESHOCK",    &
     &       "/NOZZLE EXIT",5X,"!",5X,"POSTSHOCK/B.L. EDGE  !"/4X,       &
     &       "!",46("_"),"!",27("_"),"!",26("_"),"!"/                    &
     &       4X,"!",46X,"!",27X,"]",26X,"]")                             

 1040 FORMAT (4X,"!LOCAL MACH NO.",21X,"---",8X,"!",5X,1PE11.4,11X,"!",    & 
     &       8X,1PE11.4,7X,"!"/                                            & 
     &       4X,"!LOCAL VELOCITY",21X,"FPS",8X,"!",5X,1PE11.4,11X,"!",     & 
     &       8X,1PE11.4,7X,"!"/                                            & 
     &       4X,"!LOCAL PRESSURE",21X,"PSF",8X,"!",5X,1PE11.4,11X,"!",     & 
     &       8X,1PE11.4,7X,"!"/                                            & 
     &       4X,"!LOCAL TEMPERATURE",18X,"DEG. R",5X,"!",5X,1PE11.4,11X,   & 
     &   "!",8X,1PE11.4,7X,"!"/                                            & 
     &       4X,"!LOCAL DENSITY",22X,"S/FT2",6X,"!",5X,1PE11.4,11X,"!",    & 
     &       8X,1PE11.4,7X,"!"/                                            & 
     &       4X,"!LOCAL ENTHALPY",21X,"BTU/LB",5X,"!",5X,1PE11.4,11X,      & 
     &       "!",8X,1PE11.4,7X,"!"/                                        & 
     &       4X,"!EQUILIBRIUM WALL ENTHALPY",10X,"BTU/LB",5X,"!",27X,      & 
     &       "!",8X,1PE11.4,7X,"!")

 1060 FORMAT(4X,"!NONEQUILIBRIUM WALL ENTHALPY       BTU/LB",5X,"!",       & 
     &       27X,"!",8X,1PE11.4,7X,"!"/                                    & 
     &       4X,"!DELTA",30X,"FT",9X,"!",27X,"!",8X,1PE11.4,7X,"!"/        & 
     &       4X,"!DELTA *",28X,"FT",9X,"!",27X,"!",8X,1PE11.4,7X,"!"/      & 
     &       4X,"!MOMENTUM THICKNESS",17X,"FT",9X,"!",27X,"!",8X,          & 
     &       1PE11.4,7X,"!"/                                               & 
     &       4X,"!SHOCK ANGLE",24X,"DEGREES",4X,"!",27X,"!",8X,            & 
     &       1PE11.4,7X,"!"/4X,"!",46("_"),"!",27("_"),"!",26("_"),"!")

 1070 FORMAT ("1",99X,"PAGE 2"//1X,103("_")/1X,"!",17X,"!  ECKERT   !",    &
     &       "   BARTZ   ! SPALDING  !   SQUIRE  !  PRANDTL  !",           &
     &       "  FALKNER  !  FRANKL   !"/1X,"!",17X,"!",11X,"!",11X,        &
     &       "!   CHI     ]   YOUNG   !",11X,"!",11X,"!  VOISHEL  !"/      &
     &       1X,"!",17("_"),"!",7("___________!")/                         &
     &       1X,"!",17X,"!",7(11X,"!")/                                    &
     &       1X,"!SHEAR ",8(11X,"!")/                                      &
     &       1X,"!  LB/FT2",9X,"!",7(1PE10.3," !")/                        &
     &       1X,"!",17X,"!",7(11X,"!")/                                    &
     &       1X,"!NONEQ. QDOT/A",4X,"!",7(11X,"!")/                        &
     &       1X,"!  BTU/FT2-SEC    !",7(1PE10.3," !")/                     &
     &       1X,"!",17X,"!",7(11X,"!")/                                    &
     &       1X,"!CONSERVATION EST.!",7(11X,"!")/                          &
     &       1X,"!  BLMOM/NOZMOM   !",7(1PE10.3," !")/                     &
     &       1X,"!",17X,"!",7(11X,"!")/                                    &
     &       1X,"!  BLENER/NOZENER !",7(1PE10.3," !")/                     &
     &       1X,"!",17("_"),"!",7("___________!")///)
                            
 1110 FORMAT (41X,"FROZEN PROPERTIES/MEAS. VALUES *"//1X,99("_")/          &
     &       1X,"!",7(13X,"!")/                                            &
     &       1X,"!",6X,"ZF",5X,"!     CPF     !",                          &
     &       5X,"GAMF    !     PRF     !     HAW     !      HD",           &
     &       5X,"!     TTOT    !"/1X,"!",7("_____________!")/              &
     &       1X,"!",7(13X,"!")/                                            &
     &       1X,"!",F10.4,3X,"!",F10.5,3X,"!",F10.4,3X,"!",F10.5,3X,       &
     &       "!",F10.1,3X,"!",F10.1,3X,"!",F10.1,3X,"!"/                   &
     &       1X,"!",7("_____________!"))                                   
 1120 FORMAT (///1X,99("_")/1X,"!",7(13X,"!")/                             &
     &       1X,"!    WDOTBL   !    KEBL     !    ENGINT   !    QNERAD",   &
     &       "   !    WDOTSL   !    KESL     !    INTESL   ]"/             &
     &       1X,"!",7("_____________!")/1X,"!",7(13X,"!")/                 &
     &       1X,"!",7(1PE12.5," !")/1X,"!",7("_____________!"))


END Program ColdArc   ! =====================================================
