!+
! PROGRAM SonicBoom
! ------------------------------------------------------------------------------
! PURPOSE - Compute the mean aerodynamic chord of a wing.
! AUTHORS
!        - Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 28Feb92  0.5   RLC   Original coding (Fortran 77) (from old fragments)
! 22Oct00  0.6   RLC   Recoded
! 28Dec02  0.7   RLC   Write output to a file
! 16Oct09  0.8   RLC   Write xte in output
! 21Oct09  0.81  RLC   Minor adjustments to output
! 03Nov09  0.85  RLC   All reals to double precision
!+
MODULE SonicBoomProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Replace the original blank common and encapsulate all of the 
!    subroutines of the original program.

  INTEGER,PARAMETER:: WP = SELECTED_REAL_KIND(6)

  REAL(WP),DIMENSION(9):: ay
  REAL(WP),DIMENSION(21,5,2):: data
  INTEGER:: jnr
  INTEGER:: jobs
  INTEGER:: kase
  INTEGER:: kend
  INTEGER:: kq
  INTEGER:: mp
  INTEGER:: n
  INTEGER:: nend
  INTEGER:: nn
  INTEGER:: nv
  REAL(WP),DIMENSION(21):: phi
  REAL(WP),DIMENSION(7):: pj
  REAL(WP),DIMENSION(9):: q,s,w,x,y
  REAL(WP),DIMENSION(9,100):: z      
 ! DIMENSION Z(9,100),S(9),W(9),X(9),Y(9),PJ(7),Q(9),AY(9), 
 !& DATA(21,5,2),PHI(21)

  REAL(WP):: acc     
  REAL(WP):: af      
  REAL(WP):: alt     
  REAL(WP):: apr     
  REAL(WP):: avs     
    
  REAL(WP):: b       
  REAL(WP):: bong    
  REAL(WP):: bsa     
  REAL(WP):: bsc     
  REAL(WP):: bsv     
  REAL(WP):: bvp     

  REAL(WP):: c       
  REAL(WP):: c1      
  REAL(WP):: c2      
    REAL(WP):: crv     
  REAL(WP):: cth,sth 

  
  INTEGER,PARAMETER:: DBG = 3
  REAL(WP):: dl      

  REAL(WP):: el      
  REAL(WP):: elh     
  REAL(WP):: em,en   
  REAL(WP):: ez      

  REAL(WP):: f       
    REAL(WP):: fl    
  REAL(WP):: g       

  REAL(WP):: h       
  REAL(WP):: hh      
  
  REAL(WP):: pr      
    REAL(WP):: ps      
  REAL(WP):: psi     
    
  REAL(WP):: r1      
  REAL(WP):: rc     
  
  REAL(WP):: t    
    REAL(WP):: tau     
  REAL(WP):: test      
  REAL(WP):: ts      
  
        
  REAL(WP):: u       
  REAL(WP):: vf   
  REAL(WP):: vs,vp   
  REAL(WP):: vt      

  REAL(WP):: wl 
  REAL(WP):: wt      



  REAL(WP):: xd,yd,zd



CONTAINS

!
SUBROUTINE ALTA()
! ------------------------------------------------------------------------------
! PURPOSE - Resequence all of the atmospheric data, starting at the aircraft
!  altitude and going down to the ground. Aircraft altitude must be less than or
! equal to the highest atmospheric data point.

!!!INCLUDE 'common.inc'  
!-------------------------------------------------------------------------------
  DO  K = 1, KEND  
   IF (ALT > Z (7, K) ) EXIT
  END DO  ! 50

  ka = k   
  r = (alt-z(7,ka)) / (z(7,ka-1)-z(7,ka))  

  DO J = 1, 9  
    z(j,1) = z(j,ka) + r*(z(j,ka-1) - z(j,ka))  
    DO  k = ka,kend  
      kk = k - ka + 2  
      z(j,kk) = z(j,k)
    END DO ! 53
  END DO  ! 52
  kend = kend-ka+2  
  g = (z(3,1) / z(2,1) )**0.25 / z(3,1)
  bsv = vf/t*bong**0.75  
  bsa = fl*SQRT(wt/z(2,1) / SQRT(wl))  
  RETURN  
END Subroutine Alta   ! --------------------------------------------------------

!+
SUBROUTINE ONE()
! ------------------------------------------------------------------------------
! PURPOSE - Initial conditions for all integrations are determined.
!     wind components are computed.
!     The variable I(z) in Eq.(III.4) is determined.

!!!INCLUDE'common.inc'  
!-------------------------------------------------------------------------------
  EZ = 0.  
  B = EM * EM - 1.  
  BQ = SQRT (B)  
  D = SIN (PHI (N) / 57.2958)  
  DQ = SQRT (1.0 - D * D)  
  BSC = SQRT (BSV**2 + BSA**2 * DQ * BQ / (B + 1.0) )  
  CPS = COS (PSI / 57.2958)  
  SPS = SIN (PSI / 57.2958)  
   11 AA = SQRT (D * D+ (CPS / BQ + DQ * SPS) **2)  
  STH = D / AA  
  CTH = SQRT (1.0 - STH**2)  
  ELH = - BQ * AA / EM  
  H = SQRT (1.0 - elh * elh)  
  hh = 1.0 / h  
  c = z (3, 1) / elh  
  bj = acc / (z(3,1)*em)**2 / b / elh / cth  
  bk = 1.0 - sps * (sps - cps * dq * bq) / h / h  
  bl = bj * bk  
  bm = crv * dq / em / h / h  
  bn = sps * dq - bq * cps + sth * d / cth  
  bp = bm * bn  
   12 c1 = bp - bl  
  c2 = 1.414 * bsc * (em**3 / b)**0.25  
   14 R1 = 1.0 / (Z(3,1)**1.5 * SQRT(z(2,1) * h) )  
  DO 15 M = 1, 6  
   Y (M) = 0.  
15 END DO  
   
  w(1) = 0.0  
  w(2) = - elh / h  
  w(3) = 0.0 
  w(4) = -1.0/(z(3,1)*h)  
  w(5) = c1/h/h  
  z(4,1) = 0.  
  z(5,1) = 0.  
  z(1,1) = 0.  
  z(6,1) = 1.0  
  DO 8 k = 2,kend  
    z(4,k) = cth*(z(8,k)-z(8,1))+sth*(z(9,k)-z(9,1))
    z(5,k) = cth*(z(9,k)-z(9,1))-sth*(z(8,k)-z(8,1))
    b = z(4,k)*z(3,k)**2 / (c-z(4,k))
    v = z(4,k)**2 + z(5,k)**2  
    z(1,k) = (b+v) / SQRT(z(3,k)**2+v+2.0*b)  
    DO 5 i = 1, 3  
      s(i) = 0.5 * (z(i,k) + z(i,k-1))
    5 END DO  
    b = z(1,k)-z(1,k-1)  
    v = (z(2,k)-z(2,k-1)) / s(2)
    f = 2.0*(z(3,k)-z(3,k-1)) / s(3)  
    z(6,k) = z(6,k-1)*EXP((b-0.2*s(1)*(v-f)) / (s(1)+s(3)))
  8 END DO
  q(1) = (z(3,2)+z(1,2)-z(3,1)) / (z(7,2)-z(7,1))  
13 w(6) = q(1)*elh/c/h/h  

  do 16 kz = 1, 9  
16 q(kz) = 0.  

  RETURN  
END Subroutine One   ! ---------------------------------------------------------

!+
SUBROUTINE MID()
! ------------------------------------------------------------------------------
! PURPOSE -  all integrations are carried out, using the trapezoidal method.

!!!INCLUDE'common.inc'  
!-------------------------------------------------------------------------------
  S (7) = S (7) + DL  
  IF (S (7) - Z (7, KEND) - 5.0) 65, 65, 61  
   65 DL = Z (7, KEND) + DL - S (7)  
  S(7) = Z(7,KEND)  


   61 DO 62 K = KQ, KEND  
!!!   IF (S (7) - Z (7, K) ) 62, 63, 63  
    IF (s(7) > z(7,k)) GO TO 63
   62 END DO  
   63 KQ = K - 1  

  R = (S (7) - Z (7, KQ) ) / (Z (7, KQ + 1) - Z (7, KQ) )  
  DO 64 J = 1, 6  
   S (J) = Z (J, KQ) + R * (Z (J, KQ + 1) - Z (J, KQ) )  
   64 END DO  
  RETURN  
END Subroutine Mid   ! ---------------------------------------------------------

!+
SUBROUTINE LINT()
! ------------------------------------------------------------------------------

!!!INCLUDE'common.inc'  
!-------------------------------------------------------------------------------
  DO 24 k = 2, kq  
    el = z(3,k) / (c-z(4,k))  
    IF (ABS (el) - .999) 20, 21, 21  
21    q (1) = 1.  
    s(7) = z(7,k)  

    GO TO 23  
20  en = -SQRT(1.0- el*el)  
    x(4) = 1. / (z(3,k) * en)  
    x(3) = z(5,k)*x(4)  
    x(2) = (fl*z(3,k)+z(4,k))*x(4)  
    r2 = SQRT(1.0+x(3)**2+x(2)**2)  
    x(5) = c1*r2  
    zdl = z(7,k)-z(7,k-1)  
    b = (z(3,k)-z(3,k-1))/zdl  

    ba = el*(z(4,k)-z(4,k-1))/zdl  
26  x(6) = -el*(b+ba)/c/en/h  
    b = 0.5 * zdl  
    DO 22 m = 2, 6  
      y(m) = y(m)+b*(w(m)+x(m))  
      w(m) = x(m)  
22 END DO   
   
25  f = (1.0+y(5)+y(6))/h  
    IF (f) 27, 27, 29  
27  f = f - y (5) / h  
    y(5) = 0.  
    w(5) = 0.  
    c1 = 0.0  
    q(2) = 1.0  
    q(4) = z(7, k)  
29  r4 = r2/((z(3,k)+z(1,k))**2*SQRT(f*z(2,k)/z(3,k))*z(6,k))
    x(1) = (r1-r4) / SQRT(z(7,1)-z(7,k))  
    y(1) = b * (w(1)+x(1)) + y(1)  
    w(1) = x(1)  
24 END DO  
   
   23 RETURN  
END Subroutine Lint   ! --------------------------------------------------------

!+
SUBROUTINE FIN()
! ------------------------------------------------------------------------------
! PURPOSE - 
!!!INCLUDE'common.inc'  
!-------------------------------------------------------------------------------

  y(1) = y(1) + 2.0*r1*SQRT(z(7,1)-s(7))  
  el = s(3) / (c - s(4) )  
  IF (ABS (el) - .999) 26, 27, 27  
   27 Q (1) = 1.  
   GOTO 30  
   
26 en = - SQRT(1.0 - el*el)  
  x(4) = 1. / (en * s(3) )  
  x(3) = s(5) * x(4)  
  x(2) = (el * s(3) + s(4) ) * x(4)  
  r2 = SQRT (1.0 + x(2)**2 + x(3)**2)  
  zdl = s(7) - z(7,kq)  
  b = (s(3) - z(3,kq) ) / zdl  
  ba = el * (s(4) - z(4,kq) ) / zdl  
  x(5) = c1 * r2  
23 x(6) = - el * (b + ba) / c / en / h  
  b = 0.5 * zdl  

  DO 25 m = 2,6  
   y(m) = y(m) + b*(x(m)+w(m))  
25 w(m) = x(m)  

   24 f = (1. + y(5) + y(6) ) / h  
  IF (f) 34, 34, 37  
   34 f = f - y(5)/h  
  c1 = 0.  
  q(2) = 1.  
  y(5) = 0.0  
  w(5) = 0.0  
  GO TO 37
    
37 r4 = r2/((s(3)+s(1))**2*s(6)*SQRT(f*s(2)/s(3)))
  x(1) = (r1-r4)/SQRT(z(7,1)-s(7))  
  y(1) = b*(w(1)+x(1))+y(1)  
  w(1) = x(1)  
  w(1) = w(1) - r1 / SQRT (z (7, 1) - s (7) )  
28 dl = z (7, kq)  
  kl = kend-1  
  DO 35 k = kq,kl  
    b = z(7,k)-z(7,k+1)
35 DL = MIN(dl,b)

  dl = -0.25*dl  
29 b = g*s(6)*s(3)*s(3)
  ha = SQRT(y(1)*s(2)*f*(z(7,1)-s(7))/s(3))  
  pr = c2 / b / ha  
  vs = s(3) * (1. + .4286 * pr)  
  vp = (vs - s(3) ) / dl  
33 f = f * (z (7, 1) - s (7) )  
30 RETURN  

END Subroutine Fin   ! ---------------------------------------------------------

!+
SUBROUTINE INTEG()
! ------------------------------------------------------------------------------
! PURPOSE - 
!!!INCLUDE'common.inc'  
!-------------------------------------------------------------------------------

!!  write(dbg,*) " entering subroutine Integ"

41 el = - vs / (s(4)-c)  
  IF (ABS (el)-0.999) 415, 43, 43  
43 q(1) = 1.0  

  GO TO 49  ! bail out, setting q(1)=1.0 as an error code
  
  
415 en = - SQRT (1.0 - el * el)  
  x(4) = 1.0 / en / vs  
  x(2) = (el * vs + s(4) ) * x(4)  
  x(3) = s(5) * x(4)  
  r2 = SQRT (1.0 + x(2)**2 + x(3)**2)  
  b1 = (s(4) - u) * el / dl  
  x(5) = c1 * r2  
44 x(6) = - el * (vp + b1) / c / en / h  
  b = 0.5 * dl  
  DO 42 m = 2, 6  
42 ay(m) = y(m) + b * (w(m) + x(m))  

  af = 1.0 + ay(5) + ay(6)  
  if (af - 0.01) 48, 417, 417  
48 af = af-ay(5)  
  c1 = 0.0  
  ay(5) = 0.0  
  y(5) = 0.0  
  w(5) = 0.0  
  q(5) = 3.0  
  417 f=(z(7,1)-s(7))*af/h  
  416 x(1) = r2/s(6)/((s(3)+s(1))**2 *SQRT(s(2)*f/s(3)))
  x(1) = - x(1)  
  ay(1) = y(1) + b * (w(1) + x(1))  
  b = g*s(6)*s(3)**2  
  pr = 0.5 * (c2 / b / SQRT (ay(1) * s(2) * f / s(3) ) + apr)  
  vs = 0.5 * (avs + s (3) * (1.0 + 0.4286 * pr) )  
  vp = (vs - avs) / dl  

49 RETURN  


END Subroutine Integ   ! -------------------------------------------------------

!+
SUBROUTINE CORR(OUT)
! ------------------------------------------------------------------------------
! PURPOSE -
  INTEGER,INTENT(IN):: OUT   ! external file number for output messages
!!!INCLUDE 'common.inc'  
!-------------------------------------------------------------------------------
!!! VALF is an arithmetic statement function
  VALF (XX) = SQRT (1.0 + (ALF + 2.0 * BETA * XX) **2)  

  write(dbg,*) " Entering Corr"
  CRV = CRV * 1.0E-06  
  PSI = PSI / 57.2957795  
  IF (PSI) 3, 4, 3  
    4 PSI = ABS (PSI)  
    3 TS = SIGN (TAU, PSI)  

  vt=SQRT(((em*z(3,1))*COS(psi)-z(8,1))**2 + ((em*z(3,1))*SIN(ps))**2)
  ALF = - TAN (PSI)  
  BETA = 0.5 * CRV * (1.0 + ALF**2) **1.5  
  XD = 0.0  
  QUAD = 0.  


  DO 1 I = 1, 10  
   IF (XD) 11, 12, 11  
   11    QUAD = (VALF(0.) + 4.0 * VALF (XD / 2.0) + VALF (XD) )  * XD / 6.0
   
   12    XI = XD  
   XD = XI - (QUAD-TS * VT) / VALF (XI)  
!!!   IF (ABS ( (XD-XI) / XD) - 0.01) 2, 2, 1  
    IF (ABS((xd-xi)/xd) <= 0.01) GO TO 2
1 END DO  


WRITE(OUT, 1100) XI, XD  
 1100 FORMAT(/," ERROR MESSAGE",2X, &
& "CORRECTION FOR CURVED FLIGHT PATH DID NOT CONVERGE."/)

2 ZD = (XD * BETA + ALF) * XD  
  YD = - TS * Z (9, 1)  
  IF (ZD) 5, 5, 44  
   44 IF (PSI) 45, 46, 45  
   45 XD = - XD  
  YD = - YD  
  ZD = - ZD  
  PSI = PSI - VT * CRV * TS  
  GO TO 5  
46 PSI = TAU * VT * CRV  
  ZD = - ZD  
5 PSI = PSI * 57.2957795  

  write(dbg,*) " Leaving Corr"
  RETURN  

END Subroutine Corr   ! --------------------------------------------------------

!+
SUBROUTINE SORT(OUT)
! ------------------------------------------------------------------------------
! PURPOSE - Determine the ground-shock intersection. Only used for diving or
!     climbing flight paths.
!!!DIMENSION J(21), SG(21,5)  
  INTEGER,INTENT(IN):: OUT   ! external file number for output messages
  INTEGER,DIMENSION(21):: j
  REAL(WP),DIMENSION(21,5):: sg
!!!INCLUDE 'common.inc'
!-------------------------------------------------------------------------------
  NSUM = 2 * NEND  
  TAVG = 0.0  
  DO 100 N = 1, NEND  
    J (N) = 2  
    IF (DATA(N,4,2)) 11, 11, 12  
11  DO 110 L = 1, 3  
  110    DATA(N,L,2) = 0.0  
  
   J (N) = 1  
   GOTO 10  
   
12 DATA(N,4,2) = DATA(N,4,2)-TS  
   DATA(N,3,2) = DATA(N,3,2)+YD  
   DATA(N,2,2) = DATA(N,2,2)+XD  

10 IF (DATA(N,4,1)) 13,13,14  
13 DO 130 L = 1,3  
130 DATA(N,L,1) = 0.0  
   J(N) = 1  
14 IND = J(N)  
   GO TO (15,16), IND  
15 NSUM = NSUM - 2  
   GO TO 101
16 TAVG = TAVG + DATA(N,4,1) + DATA(N,4,2)  

101 continue
100 END DO  
IF (NSUM) 17, 17, 18  
   17 K = 1  
GOTO 300  
   18 K = 2  

  TAVG = TAVG / FLOAT (NSUM)  
  DO 200 N = 1, NEND  
    IND = J (N)  
    GO TO (21, 22), IND  
   21    DO 210 L = 1, 3  
  210    SG (N, L) = 0.0  
   GOTO 201
   22    B = (DATA (N, 4, 1) - TAVG) / (DATA (N, 4, 1) - DATA (N, 4, 2) )
    DO 220 L = 1, 3  
  220    SG (N, L) = DATA (N, L, 1) - B * (DATA(N,L,1) - DATA(N,L,2))
  201 continue
  200 END DO  

  300 WRITE(OUT, 1100) KASE, EM  
 1100 FORMAT("1SONIC BOOM, CASE",I5,",  M=",F6.3/ &
& 10X,"ANGLE",7X,"PRESSURE JUMP",6X,"X",15X,"Y",10X,"TIME"/ &
& 10X,"(DEG)",10X,"(PSF)",9X,"(FT)",12X,"(FT)",9X,"(SEC)"/ &
& "   RAY-GROUND DATA")

 1101 FORMAT(/" ALTITUDE =",F7.0/(F15.2,F15.3,2F15.0,F13.2))  
WRITE(OUT, 1102)  
 1102 FORMAT(//" SHOCK-GROUND DATA"/)  
GOTO (301, 302), K  
  301 WRITE(OUT, 1103)  
 1103 FORMAT(/," NO SHOCK-GROUND DATA.")  
RETURN  
  302 WRITE(OUT, 1104) (PHI (N), (SG (N, L), L = 1, 3), TAVG, N = 1, NEND)  
 1104 FORMAT(F15.2,F15.3,2F15.0,F13.2)  
RETURN  
END Subroutine Sort   ! --------------------------------------------------------


END Module SonicBoomProcedures   ! =============================================







!+
PROGRAM Boom  
! ------------------------------------------------------------------------------

USE SonicBoomProcedures
IMPLICIT NONE
!!!INCLUDE 'common.inc'  
  REAL(WP):: adl
  REAL(WP):: bb
  REAL(WP):: dlt
  REAL(WP):: dst
  REAL(WP):: dvc
  INTEGER:: errCode
  CHARACTER(LEN=132):: fileName
  INTEGER:: i,j,k,m
  INTEGER,PARAMETER:: IDAT=1, OUT=2
  REAL(WP):: um
  REAL(WP):: v
  REAL(WP):: vg

  CHARACTER(LEN=*),PARAMETER:: GREETING = " Sonic Boom Calculator. NASA CR-157"
  CHARACTER(LEN=*),PARAMETER:: VERSION = " Version 1.01 (11 June 2017)"
  CHARACTER(LEN=*),PARAMETER:: FMT001 = "(1X,A,5F12.4/(10X,5F12.4))"

!-------------------------------------------------------------------------------
 WRITE(*,'(A/A)') GREETING,VERSION
 DO
    WRITE(*,*) "Enter the name of the input file: "
    READ(*,'(A)') fileName
    IF (fileName == ' ') STOP
    OPEN(UNIT=IDAT, FILE=fileName, IOSTAT=errCode, &
      STATUS='OLD', ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) "Unable to open this file. Try again."
  END DO
  INQUIRE(UNIT=IDAT, NAME=fileName)
  WRITE(*,*) "Reading from "//Trim(fileName)
  OPEN(UNIT=OUT,FILE="boom.out",IOSTAT=errCode,STATUS='REPLACE',ACTION='WRITE')
  IF (errCode /= 0) THEN
    WRITE(*,*) " Unable to open boom.out for output. Job terminated."
    STOP
  END IF
  
  OPEN(UNIT=DBG,FILE="boom.dbg",STATUS='REPLACE',ACTION='WRITE')

310 CONTINUE   ! keep coming back here to do additional cases
  READ(IDAT,200,IOSTAT=errCode) KEND, NEND, KASE, NV, NN, JOBS  
  IF (errCode /= 0) THEN
    WRITE(*,*) "End of data records (1). Job terminated. Output is on boom.out"
    STOP
  END IF
  READ(IDAT,201,IOSTAT=errCode) (z(7,k),z(2,k),z(1,k),z(8, k),z(9,k),k=1,kend)
  IF (errCode /= 0) THEN
    WRITE(*,*) "End of data records (2). Job terminated. Output is on boom.out"
    STOP
  END IF
  READ(IDAT,202,IOSTAT=errCode) (phi(n),n=1,nend)
  IF (errCode /= 0) THEN
    WRITE(*,*) "End of data records (3). Job terminated. Output is on boom.out"
    STOP
  END IF
  READ(IDAT,202,IOSTAT=errCode) acc,bong,rc,em,alt,vf,fl,wt,t,wl,crv,psi,tau
  IF (errCode /= 0) THEN
    WRITE(*,*) "End of data records (4). Job terminated. Output is on boom.out"
    STOP
  END IF
! all input data for this case has now been read. Proceed with computations.

  UM = ASIN(1.0/EM)*57.2958  
  IF (90.0 - UM - PSI) 11, 11, 13  
11 WRITE(OUT,204) KASE, EM, ALT  
  WRITE(OUT,218)
  218 FORMAT(/,"NO SHOCKS AT GROUND DUE TO CLIMB ANGLE")  
  GOTO 38  

13 B = (ALT - 100.0 * BONG) / 1000.  
  IF (B - Z(7,KEND) ) 12, 12, 39  
12 WRITE(OUT,204) KASE, EM, ALT  
  WRITE(OUT, 213)  

  GOTO 38  
39 DO 112 K = 1, KEND  
    Z(7,K) = Z(7,K)*1000.0  
112 Z(3,K) = 49.0 * SQRT(Z(1,K) + 459.6)  
    2 WRITE(OUT, 205) KASE,EM,ALT,ACC,RC,VF,FL,WT,BONG,T,WL,PSI,TAU,CRV
  WRITE(OUT, '(/A,6I5)') " kend,nend,kase,nv,nn=", kend,nend,kase,nv,nn  
  WRITE(OUT, 206)  
  WRITE(OUT, 207) (z(7,k), z(8,k), z(9,k), z(2,k), z(3,k), z(1,k), k=1,kend)
  DVC = 0.0  
  MP = 1  

  IF (ABS(psi) + ABS(crv)) 40, 40, 41  
41 DVC = 1.0  
   CALL CORR(OUT)  
   
40 CALL ALTA  
   
   WRITE(DBG,*) " In main program, at start of loop 111, nend=",nend
   
1  DO 111 N = 1, NEND  
     write(dbg,*)"starting do 111 with n= ", n
!!!      N=N
     DO J = 1, 5  
       DATA (N, J, MP) = 0.0  
     END DO
     CALL ONE  
     S(7) = Z(7,1) - 100.0*BONG*COS(phi(n)/57.3)  
!     ACOUSTIC INITIAL INTEGRATION
14   KQ = 1  
     DL = 0.  
     CALL MID  
     IF (KQ-1) 18,18,19  
19   CALL LINT  
     IF (Q(2)-1.0) 25,23,25  
23     Q (2) = 0.  
       WRITE(OUT, 216) Q (4), PHI (N)  
       
25   IF (Q (1) - 1.0) 18, 21, 18  
18     CALL FIN  
     IF (Q(2) - 1.0) 22, 24, 22  
24     Q(2) = 0.  
    WRITE(OUT, 216) S (7), PHI (N)  
       
22  IF (Q (1) - 1.0) 20, 21, 20  
21  WRITE(OUT, 215) S (7), PHI (N)  
    GOTO 1111  

!     SHOCK INTEGRATION STARTS HERE
20    ADL = .75 * DL  
   write(dbg,*) " in main, after statement 20"  
   DLT = DL  
80 U = S (4)  
   AVS = VS  
   BVP = VP  
   write(dbg,*) "After #80, u,avs,bvp= ", u,avs,bvp
   
81 CALL MID  
99 TEST = -2.0  
85 CALL INTEG  

98 IF (Q (1) - 1.) 3, 4, 3  
 4 Q(1) = 0.  

    IF (DL + 10.0) 78, 31, 31  
3   IF (TEST + 1.0) 82, 83, 83  
82  TEST = TEST + 1.0  
    APR = PR  
    write(dbg,*) "After 82, test,apr=", test,apr
    GO TO 85  
83  V = 2.0*ABS(PR-APR) / (PR+APR)  
    write(dbg,*) "In main, after #83, testing v= ", v
    IF (V - .01) 86, 86, 87  
87  TEST = TEST + 1.0  
    APR = PR  
    write(dbg,*) "After 87, test,apr=", test,apr
    IF (TEST - 10.0) 85, 85, 88  
88  IF (DL + 5.) 78, 77, 77  
78  S(7) = S(7)-DL  

    DL = 0.5*DL  
96  VS = AVS  
    VP = BVP  
    GO TO 81  
77  WRITE(OUT, 217) PHI (N)  

    GO TO 91  
31  WRITE(OUT, 214) PHI (N)  
91  WRITE(OUT, 211)  
    DATA(N,4,MP) = -1.0  

    GO TO 90  
86  IF (NN - N) 89, 37, 89  

37  continue
    write(out,*) " at #37, so nn=n= ", n
    IF (DVC - 1.0) 45, 45, 89  
    
45  IF (EZ) 36, 36, 90  
36  EZ = 1.0  
    WRITE(OUT, 210) (PHI (NN) )  
    WRITE(OUT, 211)  
90  PJ (6) = S (2)  
    PJ (4) = PR  
    PJ (5) = S (2) * PR  
    PJ (1) = S (7)  
    PJ (3) = AY (3) * CTH + AY (2) * STH + AY (4) * Z (9, 1)  
    PJ (2) = AY (2) * CTH - AY (3) * STH + AY (4) * Z (8, 1)  
    write(dbg,FMT001) " before write statement with 212", s,pj
    WRITE(OUT, 212) (PJ (I), I = 1, 6)  
    IF (DATA(N,4,MP)) 1111, 89, 89  
89  DO 71 M = 1, 6  
      W (M) = X(M)  
71    Y (M) = AY(M)  
      write(dbg,FMT001) "after 71, s(7),z(7,kend),q,af,dl=", s(7),z(7,kend),q(5),af,dl
    IF (DL+10.0) 72, 76, 76  
76  DL = ADL  
72  IF (S(7)-Z(7, KEND) ) 70, 70, 95  
95  IF (Q(5)-1.0) 92, 93, 94  
92  IF (AF-0.05) 97, 97, 80  
97  DL = MAX(-50.0, DL)  

    GO TO 80  
94  Q(5) = Q(5)-1.0  

    GO TO 80  
93  DL = DLT  
    Q(5) = 0.0  

    GO TO 80  
70  DATA(N,1,MP) = RC * Z(2,KEND) * PR  
    DATA(N,4,MP) = Y(4)  
    DATA(N,3,MP) = Y(4) * Z(9,1) + Y(2) * STH + Y(3)* CTH
    DATA(N,2,MP) = Y(2) * CTH - Y(3) * STH + Y(4) * Z(8,1)  

1111 continue
111 END DO  
  
   WRITE(DBG,*) "At the end of loop 111,   n=", n
  
   IF (DVC - 1.0) 42, 48, 49  
48 DVC = 2.0  
   ALT = ALT + ZD  
   PSI = PSI - VT * CRV * SIGN (TAU, PSI)  
   MP = 2  
   GO TO 40  
49 CALL SORT(OUT) 

   GO TO 38  
42 VG = SQRT (Z (9, 1) **2 + (EM * Z (3, 1) - Z (8, 1) ) **2)  
   B = Z (9, 1) / VG  
   BB = (EM * Z (3, 1) - Z (8, 1) ) / VG  

   DO 102 N = 1, NEND  
      IF (DATA(N,4,MP)) 101,101,103  
101   DATA (N, 2, MP) = 0.  
      DATA (N, 3, MP) = 0.  
      DATA (N, 1, MP) = 0.  
      GO TO 104
103   DST = VG * (DATA (N, 4, MP) - DATA (NV, 4, MP) )  
      DATA (N, 2, MP) = DATA (N, 2, MP) + DST * BB  
      DATA (N, 3, MP) = DATA (N, 3, MP) - DST * B  
104   continue
102 END DO  


    WRITE(OUT, 204) KASE, EM, ALT  
    WRITE(OUT, 203)
    WRITE(OUT, 208)  
    WRITE(OUT, 209) (PHI (N), (DATA (N, J, MP), J = 1, 4), N = 1, NEND)  
38  IF (JOBS) 311, 310, 310  
             !!! CALL EXIT
311 CONTINUE
    WRITE(*,*) "No more input data. Job terminated. Output is on boom.out"
    STOP  


  200 FORMAT(6I10)  
  201 FORMAT(F10.0,4F10.3)  
  202 FORMAT(7F10.2)  
  203 FORMAT(/," SHOCK-GROUND DATA" )  
  204 FORMAT("1SONIC BOOM, CASE ",I5, ", M=",F6.3,", ALTITUDE=",F7.0)  
  205 FORMAT("1SONIC BOOM, CASE ",I5, ", M=",F6.3,", ALIITUDE=",F7.0, &
& ", ACC=", F7.3, ", RC=", F4.2, ", VF=",F6.3," LF=",F6.3/ &
& "    WT=",F9.1, ", LENGTH=", F6.1,", FR=",F5.2,", WL=",F6.3, &
& ", ANGLE=",F7.2,", TAU=",F6.2, " SEC.", ", CURVATURE=",F7.3//)
  206 FORMAT (/,8X,"ALTITUDE",7X,"HEADWIND",8X,"SIDEWIND",8X,"PRESSURE", &
& 4X,"SOUND SPEED",4X,"TEMPERATURE"/11X,"(F1)",11X,"(FP5)",11X, &
& "(FPS)",11X, "(PSF)", 8X,"(FPS)",9X,"(DEG F)"/)
  207 FORMAT (F16.0,3F16.3,F13.3,F15.3)  
  208 FORMAT(1H0,9X,"ANGLE",7X,"PRESSURE JUMP", 6X,"X",15X,"Y",10X, &
& "TIME"/10X,"(DEG)",10X,"(PSF)",9X,"(FT)",12X,"(FT)",9X,"(SEC)"/ &
& " ")
  209 FORMAT(F15.2,F15.3,2F15.0,F13.2)  
  210 FORMAT(/,"HISTORY OF SHOCK STRENGTH VARIATION, ANGLE=", F6.2)  
  211 FORMAT("0",11X,"Z",15X,"X",15X,"Y",8X, &
& "PRESSURE RATIO     PRESSURE JUMP  PRESSURE"/ &
& 10X,"(FT)",12X,"(FT)",12X,"(FT)",31X,"(PSF)",10X,"(PSF)"/" ")
  212 FORMAT(1X,3F15.0,F19.7,F17.3,F15.3)  
  213 FORMAT(/," GROUND IS CLOSER THAN 100 BODY LENGTHS")  
  214 FORMAT(/," CUTOTF ALTITUDE, ANGLE=", F6.2)  
  215 FORMAT(/," CUTOFF BEFORE 100 BODY LENGTHS, ALTITUDE=", F10.2, &
& ", ANGLE=",F6.2)
  216 FORMAT(/, &
& "ACCELERATION EFFECTS BEFORE 100 BODY LENGTHS, ALTITUDE=", &
& F10.2,", ANGLE=",F6.2)
  217 FORMAT (/," COMPUTATION DOES NOT CONVERGE, ANGLE=",F6.2)  


END Program Boom   ! ===========================================================