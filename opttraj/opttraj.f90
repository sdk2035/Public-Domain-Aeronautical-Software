MODULE Extra

CONTAINS

!+
SUBROUTINE AT62(ZFT,ANS)
! ---------------------------------------------------------------------------
! PURPOSE - 1962 standard atmosphere
! called by Crutbl, Engepr, Trim1, UpDown, Voptrj

  REAL,INTENT(IN):: zft
  REAL,INTENT(OUT),DIMENSION(4):: ans  ! dimension is 4

      REAL PH,HZ,A,B,WA,WB,D1,D2,D3,PZ
!      DIMENSION ANS(4)
      DIMENSION HT(8),TH(8),THD(8),PH(8)
      DATA HT/0.,11.,20.,32.,47.,52.,61.,79./
      DATA TH/288.15,216.65,216.65,228.65,270.65,270.65,252.65,180.65/
      DATA THD/-6.5,0.,1.,2.8,0.,-2.,-4.,0./
      DATA PH/101325.,22632.0638,5474.88855,868.018647,110.906298,      &
     &59.0009367,18.2100724,1.03771164/
      DIMENSION ZT(13),TZ(13),TZD(13),HZ(13),A(13),B(13)
      DATA ZT/90.,100.,110.,120.,150.,160.,170.,190.,230.,300.,400.,    &
     &500.,600./
      DATA TZ/180.65, 210.65, 260.65, 360.65, 960.65, 1110.65, 1210.65, &
     & 1350.65, 1550.65, 1830.65, 2160.65, 2420.65, 2590.65/
      DATA TZD/3.,5.,10.,20.,15.,10.,7.,5.,4.,3.3,2.6,1.7,1.1/
      DATA HZ/88.7433565,98.4509829,108.128578,117.776280,146.541401,   &
     &156.070901,165.571187,184.484657,221.966870,286.476269,376.312415,&
     &463.526097,548.230014/
      DATA A/.99999916,.99999897,.99999877,.99999832,.99999776,.99999746&
     &,.99999698,.99999592,.99999355,.99998878,.99998131,.99997196,     &
     &.99996075/
      DATA B/.00015734766,.00015734953,.00015735140,.00015735513,       &
     &.00015735887,.00015736074,.00015736355,.00015736915,.00015737943, &
     &.00015739532,.00015741401,.00015743271,.00015745140/
      DIMENSION WA(13),WB(13),WC(13)
      DATA WA/21.998808, 15.798995, 31.044527, 40.387675, 29.538575,    &
     &32.268971, 27.789444, 32.166670, 30.241635, 34.561172,            &
     &36.099504, 38.195672, 18.258073/
      DATA WB/.15479092, .27878720, .0015957013, -.15412343,-.0094687678&
     &, -.043598715, .0091016009, -.036974463, -.020235026,             &
     &-.049031942, -.056723605, -.065108273, .0013503901/
      DATA WC/-.85994958D-3, -.14799309D-2, -.21996960D-3, .42886012D-3,&
     &-.53322091D-4, .53333994D-4, -.10166693D-3, .19585867D-4,         &
     &-.16804213D-4, .31190648D-4, .40805227D-4, .49189895D-4,          &
     &-.61923241D-5/
      DIMENSION D1(13),D2(13),D3(13)
      DATA D1/.0017834765,.0010654122,.00053055610,.00026454351,        &
     &.00035360997,.00053348782,.00076836496,.0010889831,.0013783559,   &
     &.0016975137,.0022189663,.0037023997,.0067578185/
      DATA D2/-11.281753,-6.7098914,-3.3278396,-1.6546388,-2.2171667,   &
     &-3.3643151,-4.8850055,-7.0083025,-8.9810162,-11.235530,-15.122423,&
     &-27.520411,-59.311259/
      DATA D3/.016920782,.024329051,.039545102,.057409044,.016199137,   &
     &.0093014845,.0059339235,.0037645169,.0026065966,.0018120459,      &
     &.0011923023,.00064736059,.00033627561/
      DIMENSION PZ(13)
      DATA PZ/.16438012,.030075034,.0073545270,.0025216927,             &
     &.00050617890,.00036943532,.00027926462,.00016852498,.69605367D-4, &
     &.18838777D-4,.40304321D-5,.10956964D-5,.34502614D-6/
!----------------------------------------------------------------------------
      ALT=ZFT*0.3048
      Z=ALT/1000.    ! kilometers

      IF (Z.LT.-5.)Z=-5.
      IF (Z.GT.700.)Z=700.
      IF (Z.GT.90.)GO TO 90

      DEN=1.0+0.00015733831D0*Z
      H=Z/DEN
      GMW=28.9644
      IF(H.GE.47.)GO TO 47
      IF(H.GE.20.)GO TO 20
      J=1
      IF(H.GE.11.)J=2
      GO TO 21
   20 J=3
      IF(H.GE.32.)J=4
      GO TO 21
   47 IF(H.GE.61.)GO TO 61
      J=5
      IF(H.GE.52.)J=6
      GO TO 21
   61 J=7
      IF(H.GE.79.)J=8
   21 TM=TH(J)+THD(J)*(H-HT(J))
      IF(THD(J).EQ.0.)GO TO 5
      PLOG=-34.163195D0*ALOG(TM/TH(J))/THD(J)
      GO TO 2
    5 PLOG=-34.163195D0*(H-HT(J))/TM
    2 PB=PH(J)
      GO TO 100
   90 IF(Z.LT.170.)GO TO 11
      IF(Z.LT.300.)GO TO 12
      IF(Z.LT.500.)GO TO 13
      J=13
      IF(Z.LT.600.)J=12
      GO TO 10
   13 J=11
      IF(Z.LT.400.)J=10
      GO TO 10
   12 J=9
      IF(Z.LT.230.)J=8
      IF(Z.LT.190.)J=7
      GO TO 10
   11 IF(Z.LT.120.)GO TO 14
      J=6
      IF(Z.LT.160.)J=5
      IF(Z.LT.150.)J=4
      GO TO 10
   14 J=3
      IF(Z.LT.110.)J=2
      IF(Z.LT.100.)J=1
   10 GMW=WA(J)+Z*(WB(J)+Z*WC(J))
      TM=TZ(J)+TZD(J)*(Z-ZT(J))
      DEN=A(J)+Z*B(J)
      H=Z/DEN
      DELTAH=H-HZ(J)
      PLOG=D1(J)*DELTAH+D2(J)*ALOG(1.0+D3(J)*DELTAH)
      PB=PZ(J)
  100 P=PB*EXP(PLOG)
      ANS(1)=6.75944794D-6*P/TM
      ANS(2)=P*0.020885434D0
      ANS(3)=GMW*TM/28.9644
      ANS(4)=894.50046D0
      ARG1=4325.73899D0*TM
      IF(ZFT.LT.300000.)ANS(4)=SQRT(ARG1)
      RETURN
END Subroutine AT62   ! -----------------------------------------------------

!+
SUBROUTINE ICLOCK(TIME, IHR, IMIN, ISEC)
! ---------------------------------------------------------------------------
! PURPOSE - break out the time in seconds (?) (REAL) into hours, min, sec
  REAL,INTENT(IN):: time
  INTEGER,INTENT(OUT):: ihr, imin, isec

  INTEGER:: itime
  INTEGER,PARAMETER:: SIXTY = 60
!----------------------------------------------------------------------------
  itime = time    ! truncates to integer
  imin = itime/SIXTY
  isec = itime - imin*SIXTY
  ihr = imin/SIXTY
  imin = imin - ihr*SIXTY
  RETURN
END Subroutine Iclock   ! ---------------------------------------------------

!+
FUNCTION Jtrunc(x,n) RESULT(j)                              ! a PURE function
! ---------------------------------------------------------------------------
! PURPOSE - Find the first element of the array X that is NOT less than
!   the following element.
! called by Cruzop and Wlefhv

  REAL,INTENT(IN),DIMENSION(:):: x
  INTEGER,INTENT(IN):: n

  INTEGER:: j

  INTEGER k
!----------------------------------------------------------------------------
  DO k=1,SIZE(x)-1
   IF (x(k)<=x(k+1)) THEN
     j=k
     RETURN
   END IF
  END DO
  j=SIZE(x)
  RETURN
END Function Jtrunc   ! -----------------------------------------------------

!+
SUBROUTINE Lsqpol(x,y,w,coeff)
! ---------------------------------------------------------------------------
! PURPOSE - Least square polynomial fit
! NOTE - the original program allowed multiple y (2-d)
!  complete restructure of arguments - RLC - 20Oct99

  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(IN),DIMENSION(:):: y
  REAL,INTENT(IN),DIMENSION(:):: w  ! weights
  REAL,INTENT(OUT),DIMENSION(:):: coeff

  REAL,ALLOCATABLE,DIMENSION(:,:):: a
  REAL,ALLOCATABLE,DIMENSION(:,:):: b
  REAL,ALLOCATABLE,DIMENSION(:,:):: c
  INTEGER:: i,j,k

  INTEGER:: m
  INTEGER:: n
!----------------------------------------------------------------------------
  n=SIZE(x)
  m=SIZE(coeff)
  ALLOCATE(a(m,m), b(m,1), c(n,m) )

  c(:,1)=1.0
  DO k=2,m
    c(:,k)=c(:,k-1)*x(:)
  END DO

  a(:,:)=0.0
  DO I=1,M
    DO J=1,M
      DO K=1,N
        a(i,j)=a(i,j)+c(k,i)*c(k,j)*w(k)
      END DO
    END DO
  END DO

  b=0.0
  DO i=1,m
    DO k=1,n
      b(i,1)=b(i,1)+c(k,i)*y(k)*w(k)
    END DO
  END DO

  CALL MATINV (a(1:m,1:m),b(1:m,1:1),determ)
  coeff(:)=b(:,1)

  DEALLOCATE(a,b,c)
  RETURN
END Subroutine Lsqpol   ! ---------------------------------------------------


!+
SUBROUTINE MATINV(A,B,DETERM)
! ---------------------------------------------------------------------------
! PURPOSE - MATRIX INVERSION WITH ACCOMPANYING SOLUTION OF LINEAR EQUATIONS
  REAL,INTENT(IN OUT),DIMENSION(:,:):: a
  REAL,INTENT(IN OUT),DIMENSION(:,:):: b
  REAL,INTENT(OUT):: determ
!
!      DIMENSION A(20,20),B(20,1)
!      COMMON/ANE206/PIVOT(20),INDEX(20,2),IPIVOT(20),DUMI(320)
      EQUIVALENCE (IROW,JROW), (ICOLUM,JCOLUM), (AMAX, T, SWAP)

  INTEGER:: i,j,k
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: index
  INTEGER,ALLOCATABLE,DIMENSION(:):: ipivot
  INTEGER:: m,n
  REAL,ALLOCATABLE,DIMENSION(:):: pivot
!----------------------------------------------------------------------------
  n=SIZE(a,1)
  m=SIZE(b,2)

  ALLOCATE(index(n,2), ipivot(n), pivot(n) )
!
!     INITIALIZATION
!
   10 DETERM=1.0
      ipivot(1:n)=0
!   15 DO 20 J=1,N
!      IPIVOT(J)=0
!   20 END DO

   30 DO 551 I=1,N
!
!     SEARCH FOR PIVOT ELEMENT
!
   40 AMAX=0.0
   45 DO 105 J=1,N
   50 IF (IPIVOT(J)-1) 60, 105, 60
   60 DO 100 K=1,N
   70 IF (IPIVOT(K)-1) 80, 100, 740
   80 IF (ABS(AMAX)-ABS(A(J,K))) 85, 100, 100
   85 IROW=J
   90 ICOLUM=K
   95 AMAX=A(J,K)
  100 END DO
  105 END DO
  110 IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1
!
!     INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL
!
  130 IF (IROW-ICOLUM) 140, 260, 140
  140 DETERM=-DETERM
  150 DO 200 L=1,N
  160 SWAP=A(IROW,L)
  170 A(IROW,L)=A(ICOLUM,L)
  200 A(ICOLUM,L)=SWAP
!  205 IF(M) 260, 260, 210

  210 DO 250 L=1, M
  220 SWAP=B(IROW,L)
  230 B(IROW,L)=B(ICOLUM,L)
      B(ICOLUM,L)=SWAP
  250 END DO

  260 INDEX(I,1)=IROW
  270 INDEX(I,2)=ICOLUM
  310 PIVOT(I)=A(ICOLUM,ICOLUM)
  320 DETERM=DETERM*PIVOT(I)
!
!     DIVIDE PIVOT ROW BY PIVOT ELEMENT
!
  330 A(ICOLUM,ICOLUM)=1.0
  340 DO 350 L=1,N
      A(ICOLUM,L)=A(ICOLUM,L)/PIVOT(I)
  350 END DO

!  355 IF(M) 380, 380, 360
  360 DO 370 L=1,M
      B(ICOLUM,L)=B(ICOLUM,L)/PIVOT(I)
  370 END DO
!
!     REDUCE NON-PIVOT ROWS
!
  380 DO 550 L1=1,N
  390 IF(L1-ICOLUM) 400, 550, 400
  400 T=A(L1,ICOLUM)
  420 A(L1,ICOLUM)=0.0
  430 DO 450 L=1,N
      A(L1,L)=A(L1,L)-A(ICOLUM,L)*T
  450 END DO
! 455 IF(M) 550, 550, 460
  460 DO 500 L=1,M
      B(L1,L)=B(L1,L)-B(ICOLUM,L)*T
  500 END DO
  550 END DO
  551 END DO

!
!     INTERCHANGE COLUMNS
!
  600 DO 710 I=1,N
  610 L=N+1-I
  620 IF (INDEX(L,1)-INDEX(L,2)) 630, 710, 630
  630 JROW=INDEX(L,1)
  640 JCOLUM=INDEX(L,2)
  650 DO 705 K=1,N
  660 SWAP=A(K,JROW)
  670 A(K,JROW)=A(K,JCOLUM)
  700 A(K,JCOLUM)=SWAP
  705 END DO
  710 END DO

  DEALLOCATE(index, ipivot, pivot)
  740 RETURN
END Subroutine Matinv   ! ---------------------------------------------------


!+
FUNCTION MINF2(AX,BX, F, X, IPRINT) RESULT(fOutput)
! ---------------------------------------------------------------------------
! PURPOSE - Find the minimum of F(X) OVER [AX, BX].
!   Fibonnaci search

  REAL,INTENT(IN):: ax,bx
  REAL,INTENT(OUT):: x
  INTEGER,INTENT(IN):: iprint   !      IPRINT = 0 NO PRINT

  REAL:: fOutput

INTERFACE
  FUNCTION f(x) RESULT(ff)
    REAL,INTENT(IN):: x
    REAL:: ff
  END Function f
END INTERFACE

  
      
  REAL,DIMENSION(8):: dex
  REAL:: dxab

  REAL,PARAMETER,DIMENSION(8):: FIBONO = (/ 0.3818182, 0.2363636, &
    0.1454545, 0.09090909, 0.05454545, 0.03636361, 0.01818182, 0.01818182 /)
  INTEGER:: i
  REAL:: x1,x2
  REAL:: xa,xb

!----------------------------------------------------------------------------
   90 XA = AX
      XB = BX
      DXAB = XB - XA
      dex(:)=dxab*FIBONO(:)
!      DO 100 I= 1,8
!  100 DEX(I) = FIBONO(I) *DXAB

      I= 1
      X1 = XA + DEX(I)
      X2 = XB - DEX(I)
      FX1 = F(X1)
      FX2 = F(X2)
!
!
      IF( IPRINT .EQ. 0) GO TO 101
      WRITE(6, 50)
   50 FORMAT(/ 'XA, X1, X2, XB        DEX(I),       FX1, FX2')
  101 I= I+ 1
      IF( IPRINT .EQ. 0) GO TO 151
      WRITE(6, 51) XA, X1, X2, XB,DEX(I -1), FX1, FX2
   51 FORMAT(/4F10.2, 5X, F10.4, 2F15.4)
  151 IF( FX2 .LT. FX1) GO TO 201
      XB = X2
      X2= X1
      X1 = XA + DEX(I)
      FX2 = FX1
      FX1 = F(X1)
      GO TO 203
  201 XA = X1
      X1 = X2
      X2 = XB - DEX(I)
      FX1 = FX2
      FX2 = F(X2)
  203 IF(I .GE. 8) GO TO 204   ! always do 8 iterations
      GO TO 101

  204 IF( IPRINT .EQ. 0) GO TO 152
      WRITE(6, 51) XA, X1, X2, XB,DEX(I -1), FX1, FX2

  152 X = 0.5*(X1 + X2)

  fOutput=F(X)
  RETURN
END Function MinF2   ! ------------------------------------------------------

!+
FUNCTION Minf(ax, bx, F, x, iprint) RESULT(fOutput)
! ---------------------------------------------------------------------------
! PURPOSE - minimize F(X) over the interval [AX, BX]
IMPLICIT NONE

  REAL,INTENT(IN):: ax,bx       ! limits of interval
  REAL,INTENT(OUT):: x          ! x-location of minimum
  INTEGER,INTENT(IN):: iprint   ! 0=no print; otherwise print

  REAL:: fOutput                ! F(x)

INTERFACE
  FUNCTION f(x) RESULT(ff)
    REAL,INTENT(IN):: x
    REAL:: ff
  END Function f
END INTERFACE

  REAL,DIMENSION(10):: dex
  REAL:: dxab
  REAL,PARAMETER,DIMENSION(10):: FIBONO = (/ &
    0.38194,  0.23611,  0.145833, 0.09028,   0.055555, &
    0.034722, 0.020833, 0.013889, 0.0069444, 0.006944 /)


!      DIMENSION FIBONO(10), DEX(10)
!      DATA FIBONO/.38194,.23611, .145833, .09028, .055555, .034722,     &
!     & .020833, .013889,  2*.0069444/
  CHARACTER(LEN=*),PARAMETER:: FMT50 = &
    "XA, X1, X2, XB        DEX(I),    FX1, FX2"
  CHARACTER(LEN=*),PARAMETER:: FMT51 = '(4F10.2, 5X, F10.4, 2F15.4)'
  REAL:: fx1,fx2
  INTEGER:: i
  REAL:: x1,x2
  REAL:: xa,xb


!----------------------------------------------------------------------------
! use 4 points:  xa,x1,x2,xb to enclose the minimum and close in on it
!   by successive iterations. This implementation always uses 10 iterations.

  IF (iprint/=0) WRITE(6,*) FMT50
!   50 FORMAT(/ 'XA, X1, X2, XB        DEX(I),       FX1, FX2')

  xa = ax
  xb = bx
  dex(:)=(xb-xa)*FIBONO(:)

  i= 1
  x1 = xa + dex(i)
  x2 = xb - dex(i)
  fx1 = f(x1)
  fx2 = f(x2)

!   51 FORMAT(/4F10.2, 5X, F10.4, 2F15.4)

  DO i=2,10
    IF (iprint/=0) WRITE(6,FMT51) xa, x1, x2, xb, dex(i-1), fx1, fx2

    IF (fx2 .lt. fx1) THEN
      xa = x1                 ! x can't be less than x1
      x1 = x2
      x2 = xb - dex(i)        ! but don't change xb
      fx1 = fx2
      fx2 = F(x2)
    ELSE
      xb = x2                 ! x can't be greater than x2
      x2= x1
      x1 = xa + dex(i)        ! but don't change xa
      fx2 = fx1
      fx1 = f(x1)
    END IF
  END DO

  IF (iprint/=0) WRITE(6,FMT51) xa, x1, x2, xb, dex(i-1), fx1, fx2
  x = 0.5*(x1+x2)
  fOutput= F(x)
  RETURN
END Function MinF   ! -------------------------------------------------------

!+
SUBROUTINE Page()
! ---------------------------------------------------------------------------
! PURPOSE - 
!----------------------------------------------------------------------------
  RETURN
END Subroutine Page   ! -----------------------------------------------------

!+
FUNCTION Polye1(x,m,b) RESULT(f)  ! a PURE function
! ---------------------------------------------------------------------------
! PURPOSE - Evaluate polynomial
!   f=b(1) + b(2)*x + b(3)*x**2 + b(4)*x**3 + ... + b(m)*x**(m-1)
! soon, drop m from argument list

  REAL,INTENT(IN):: x
  INTEGER,INTENT(IN):: m
  REAL,INTENT(IN),DIMENSION(m):: b

  REAL:: f
  INTEGER:: k

!----------------------------------------------------------------------------
  f=b(m)
  DO k=m-1,1,-1
    f = f*x + b(k)
  END DO
  RETURN
END Function Polye1   ! -----------------------------------------------------


!+
SUBROUTINE SERCH1(tx,x,nx,pf,l,limit)
! ---------------------------------------------------------------------------
! PURPOSE -

  REAL,INTENT(IN),DIMENSION(:):: tx
  REAL,INTENT(IN):: x
  INTEGER,INTENT(IN):: nx
  REAL,INTENT(OUT):: pf
  INTEGER,INTENT(OUT):: l
  INTEGER,INTENT(OUT):: limit

!!!      DIMENSION TX(1)
!XXXXX TX(N),N=1....NX  IS ASSUMED IN ASCENDING ORDER
!----------------------------------------------------------------------------
!
    7 LIMIT = 0
      L = NX/2
      IF(X.LT.TX(1).OR.X.GT.TX(NX)) GO TO 20
    8 A = TX(L)
      B = TX(L+1)
      IF (X .LT. A) GO TO 11
      IF( X .GT. B) GO TO 12
      PF=(X- A)/(B-A)
      RETURN
   11 L=L-1
      IF (L .LE. 0) GO TO 20
      GO TO 8
   12 L=L+1
      IF (L .GT. NX) GO TO 20
      IF( L .EQ. NX) GO TO 19
      GO TO 8
   19 PF = 1.
      L = NX - 1
      RETURN
   20 LIMIT= 1
      WRITE(6,100) L,NX,X, TX(1), TX(NX)
  100 FORMAT(/'INPUT TO SERCH1 OUTSIDE TABLE (L,NX,X,TX(L),TX(NX)='  &
     & 2I4, 4X, 3F15.7)
      RETURN
END Subroutine Serch1   ! ---------------------------------------------------


!+
SUBROUTINE SERCHD(TX,X,NX,PF,L,LIMIT)
! ---------------------------------------------------------------------------
! PURPOSE - 
  REAL,INTENT(IN),DIMENSION(:):: tx
  REAL,INTENT(IN):: x
  INTEGER,INTENT(IN):: nx
  REAL,INTENT(OUT):: pf
  INTEGER,INTENT(OUT):: l
  INTEGER,INTENT(OUT):: limit


!!!      DIMENSION TX(1)
!XXXXX TX(N),N=1....NX  ARE MONOTONICALLY DECREASING
!----------------------------------------------------------------------------
!
    7 LIMIT = 0
      L = NX/2
      IF(X.GT.TX(1).OR.X.LT.TX(NX)) GO TO 20
    8 B = TX(L)
      A = TX(L+1)
      IF (X .LT. A) GO TO 12
      IF( X .GT. B) GO TO 11
      PF=(B- X)/(B-A)
      RETURN

   11 L=L-1
      IF (L .LE. 0) GO TO 20
      GO TO 8

   12 L=L+1
      IF (L .GT. NX) GO TO 20
      IF( L .EQ. NX) GO TO 19
      GO TO 8

   19 PF = 1.
      L = NX - 1
      RETURN

   20 LIMIT= 1
      WRITE(6,100) L,NX,X, TX(1), TX(NX)
  100 FORMAT(/'INPUT TO SERCH1 OUTSIDE TABLE (L,NX,X,TX(L),TX(NX)='  &
     & 2I4, 4X, 3F15.7)
      RETURN
END Subroutine Serchd   ! ---------------------------------------------------


!+
FUNCTION TwoDinterp(xtab,ytab,ftab, x,y) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Two dimensional linear interpolation in rectangular grid

IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: xtab,ytab ! increasing
  REAL,INTENT(IN),DIMENSION(:,:):: ftab
  REAL,INTENT(IN):: x,y

  REAL:: f

  INTEGER:: j,k
!----------------------------------------------------------------------------
  f=0.0
  RETURN
END Function TwoDinterp   ! -------------------------------------------------

!+
FUNCTION VALUE2(TEMP, IT, ITMAX, TEMPA, ALT24, IALT, IALMAX, H, EPRMAX, I)
! ---------------------------------------------------------------------------
! PURPOSE -

!! USE Extra,ONLY:Serch1,Xtrpl1,Xtrpl2

  REAL,INTENT(IN),DIMENSION(:):: temp
  INTEGER,INTENT(IN):: it
  INTEGER,INTENT(IN):: itmax
  REAL,INTENT(IN):: tempa
  REAL,INTENT(IN),DIMENSION(:):: alt24
  INTEGER,INTENT(IN):: ialt
  INTEGER,INTENT(IN):: ialmax
  REAL,INTENT(IN):: h
  REAL,INTENT(IN),DIMENSION(:,:):: eprmax
  INTEGER,INTENT(IN):: i

!      DIMENSION TEMP(IT), ALT24(IALT), EPRMAX(IT, IALT)
!
!
!
!     IT, IALT, FOR DIMENSIONING, ITMAX, IALMAX FOR SEARCH1
!     ITEMP ROW INDEX, IH COLUMN INDEX
!----------------------------------------------------------------------------
!
   70 IF ( I .EQ. 0) GO TO 80
      IF( H .LE. ALT24(1)) GO TO 101
      IF( H .GE. ALT24(IALMAX)) GO TO 102
      IF ( TEMPA .LT. TEMP(1)) GO TO 106
      IF( TEMPA .GT. TEMP(ITMAX)) GO TO 107
   80 CALL SERCH1(ALT24, H,  IALMAX, PFH, IH, LIMIT)
      CALL SERCH1(TEMP, TEMPA, ITMAX, PFT, ITEMP, LIMIT)
      EPRMX1 = EPRMAX(ITEMP, IH) + PFT*(EPRMAX(ITEMP + 1, IH)           &
     & - EPRMAX(ITEMP    , IH))
      EPRMX2 = EPRMAX(ITEMP, IH+ 1) + PFT*( EPRMAX(ITEMP + 1, IH+ 1) -  &
     & EPRMAX(ITEMP    , IH+ 1))
  100 VALUE2 = EPRMX1 + PFH*(EPRMX2 - EPRMX1)
      RETURN
!
!
  101 IH = 1
      GO TO 103
  102 IH = IALMAX
  103 IF( TEMPA .LT. TEMP(1)) GO TO 104
      IF( TEMPA .GT. TEMP(ITMAX)) GO TO 105
      CALL SERCH1(TEMP, TEMPA, ITMAX, PFT, ITEMP, LIMIT)
      VALUE2 = EPRMAX(ITEMP, IH) + PFT*(EPRMAX(ITEMP + 1, IH)           &
     & - EPRMAX(ITEMP    , IH))
      RETURN
  104 VALUE2 = XTRPL1(TEMP(1), TEMP(2), TEMPA, EPRMAX(1,1), EPRMAX(2,1))
      RETURN
  105 N1 = ITMAX - 1
      VALUE2 = XTRPL2 (TEMP(N1), TEMP(ITMAX), TEMPA, EPRMAX(N1,IALMAX), &
     & EPRMAX(ITMAX, IALMAX))
      RETURN
  106 CALL SERCH1(ALT24, H,  IALMAX, PFH, IH, LIMIT)
      EPRMX2 = XTRPL1(TEMP(1), TEMP(2), TEMPA, EPRMAX(1,IH + 1),        &
     &EPRMAX(2, IH+1))
      EPRMX1 = XTRPL1(TEMP(1), TEMP(2), TEMPA, EPRMAX(1,IH),EPRMAX(2,IH)&
     &)
      GO TO 100
  107 CALL SERCH1(ALT24, H,  IALMAX, PFH, IH, LIMIT)
      N1 = ITMAX - 1
      EPRMX1= XTRPL2(TEMP(N1), TEMP(ITMAX), TEMPA, EPRMAX(N1,IH),       &
     & EPRMAX(ITMAX, IH))
      EPRMX2= XTRPL2(TEMP(N1), TEMP(ITMAX), TEMPA, EPRMAX(N1,IH + 1),   &
     & EPRMAX(ITMAX, IH+1))
      GO TO 100
END Function Value2   ! ---------------------------------------------------



!+
FUNCTION  XTRPL1(a1, a2, a, b1, b2) RESULT(f)  ! a PURE function
! ---------------------------------------------------------------------------
! PURPOSE - straight line thru (a1,b1) and (a2,b2). Evaluate at x=a
  REAL,INTENT(IN):: a1,a2, a, b1,b2
  REAL:: f
!----------------------------------------------------------------------------
  f = b2 - (a2-a)*(b2-b1)/(a2-a1)
  RETURN
END Function Xtrpl1   ! -----------------------------------------------------

!+
FUNCTION XTRPL2(an1, an, a, bn1, bn) RESULT(f)  ! a PURE function
! ---------------------------------------------------------------------------
! PURPOSE - straight line thru (an1,bn1) and (an,bn). Evaluate at x=a
  REAL,INTENT(IN):: an1,an,a,bn1,bn
  REAL:: f
!----------------------------------------------------------------------------
  f = bn1 + (a-an1)*(bn-bn1)/(an-an1)
  RETURN
END Function Xtrpl2   ! -----------------------------------------------------





END Module Extra   ! ========================================================

!+
MODULE Engine
! ---------------------------------------------------------------------------


!----------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE CPMEP1 (IPRINT)
! ---------------------------------------------------------------------------
! PURPOSE -
!   Reads from unit 9 (two places and in each call to EPRIO)
  INTEGER,INTENT(IN):: iprint


  INCLUDE 'boeing.inc'
!      REAL MACHNO
!      REAL MACH24,MACH27,MACHN1,MACHN2
!      COMMON/BOEING/ALT(10),MACHNO(10),FNIDL(10,10),WFIDL(10,10)        &
!     &,ALT24(7),TEMP(13),EPRMAX(13,7),ALTFF(10)                         &
!     &,EPR(14),FNMAX(14,10)                                             &
!     &,DECL(21),DECD(21,10),MACH24(10)                                  &
!     &,ALT27(9),THRUST(30,9),MACH27(8,9),TSFC(30,8,9)                   &
!     &,MACHN1(10),MACHN2(10)

!      DIMENSION HDING(20),LABLE(20),FMTINP(5),FMTOUT(10)

  CHARACTER(LEN=132):: dummyLabel
  CHARACTER(LEN=80):: fmtinp, fmtout
  INTEGER:: irow,icol
  INTEGER:: jprint,kprint
  INTEGER:: k
  REAL,DIMENSION(10,10):: x
!----------------------------------------------------------------------------
!
!     IDLE THRUST IN LBS. ENTERED
!
  write(3,*) "Entering cpmep1, iprint=", iprint

  icol=10
  irow=1
  fmtinp='(10F7.3)'
  fmtout='(" ALT/MACHNO    ", 10F10.3)'
  jprint=1
  kprint=1
  CALL EPRIO(MACHNO(1:10),X,icol,irow,fmtinp, fmtout, IPRINT,jprint,kprint)
!  READ(9,'(A)') dummyLabel
!  READ(9,'(A)') dummyLabel
!  READ(9,fmtinp) machno(1:10)

  fmtinp='(11F7.2)'
  fmtout='(F11.0,5X,10F10.0)'
  icol=10
  irow=10
  CALL EPRIO(ALT,FNIDL,icol,irow,fmtinp,fmtout,IPRINT,0,0)

  fmtinp='(10F7.3)'
  fmtout='(" ALT/MACHNO    ", 10F10.3)'
  CALL EPRIO(MACHN1,X,10,1,fmtinp,fmtout,IPRINT,1,0)

  fmtinp='(11F7.2)'
  fmtout='(F11.0,5X,10F10.0)'
  CALL EPRIO(ALTFF,WFIDL,10,10,fmtinp,fmtout,IPRINT,0,0)

  ! CONTINUE with these
  fmtinp='(8F7.0)'
  fmtout='(11X,5F10.0)'
  CALL EPRIO(ALT24,X,5,1,fmtinp,fmtout,IPRINT,1,1)

  fmtinp='(8F7.3)'
  fmtout='(8F10.3)'
  CALL EPRIO(TEMP,EPRMAX,5,13,fmtinp,fmtout,IPRINT,0,0)

  fmtinp='(10F7.2)'
  fmtout='(" EPR/ FN/DELTA(AM)", F7.3, 9F10.3)'
  CALL EPRIO(MACHN2,X,10,1,fmtinp,fmtout,IPRINT,1,0)

  fmtinp='(F7.2,10F7.0)'
  fmtout='(F11.2,5X,10F10.0)'
  CALL EPRIO(EPR,FNMAX,10,14,fmtinp,fmtout,IPRINT,0,0)

  fmtinp='(10F7.3)'
  fmtout='(" DELTA CL/MACHNO", F6.3, 9F10.3)'
  CALL EPRIO(MACH24,X,10,1,fmtinp,fmtout,IPRINT,1,1)

  fmtinp='(11F7.5)'
  fmtout='(11F10.5)'
  CALL EPRIO(DECL,DECD,10,21,fmtinp,fmtout,IPRINT,0,0)

      DO 107 I=1,9
      READ(9,72) ALT27(I),(MACH27(J,I),J=1,8)
   72 FORMAT(F7.0,8F7.4)
      IF (IPRINT .EQ. 0) GO TO 105
      WRITE(6,70)
   70 FORMAT('1CHART NO. 27273A')
      WRITE(6,71)
   71 FORMAT(/' WEIGHTED AVERAGE TSFC TABLE')
      WRITE(6,73)I,ALT27(I)
   73 FORMAT(' PAGE ',I1,5X,'ALT = ',F7.0)
      WRITE(6,74)(MACH27(J,I),J=1,8)
   74 FORMAT(/' THRUST/MACHNO',F9.5,7F10.4)

  105 DO 108 K = 1,30
      READ(9,75) THRUST(K,I),(TSFC(K,J,I),J=1,8)
   75 FORMAT(F7.0,8F7.4)
      IF (IPRINT .EQ. 0) GO TO 108
      WRITE(6,76) THRUST(K,I),(TSFC(K,J,I),J=1,8)
   76 FORMAT(F7.0,4X,8F10.4)
  108 END DO
  107 END DO

  RETURN
END Subroutine Cpmep1   ! ---------------------------------------------------

!+
SUBROUTINE ENGEPR(H, MAKNO, EPRX,INCRUZ, THRST, FDOT,MFGR)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: At62,PolyE1,Serch1,Value2
  REAL,INTENT(IN):: h
  REAL,INTENT(IN):: makno
  REAL,INTENT(IN OUT):: eprx
  INTEGER,INTENT(IN):: incruz
  REAL,INTENT(OUT):: thrst   ! really in out ??
  REAL,INTENT(OUT):: fdot
  INTEGER,INTENT(IN):: mfgr

!      REAL MACHNO,
  REAL:: KC
!      REAL MACH24,MACH27
!  REAL:: thrust
!  REAL:: tsfc

  INCLUDE 'boeing.inc'
!      COMMON/BOEING/ALT(10),MACHNO(10),FNIDL(10,10),WFIDL(10,10)        &
!     &,ALT24(7),TEMP(13),EPRMAX(13,7),ALTFF(10)                         &
!     &,EPR(14),FNMAX(14,10)                                             &
!     &,DECL(21),DECD(21,10),MACH24(10)                                  &
!     &,ALT27(9),THRUST(30,9),MACH27(8,9),TSFC(30,8,9)                   &
!     &,MACHN1(10),MACHN2(10)

      COMMON/COST/EPR1,ICOST,FC,TC,DTEMPK,W,FUELDT
 
!      DIMENSION ANS(4)   ! WF(6,5)
  INTEGER,PARAMETER,DIMENSION(9):: MAXMAK = (/6, 6, 7, 7, 7, 7, 6, 6, 6/)
  INTEGER,PARAMETER,DIMENSION(9):: MAXT = &
                                   (/30, 30, 30, 30, 30, 30, 28, 29, 29/)
!      DATA PZ , TZ/2116.2, 288./

!  REAL:: alt
  REAL,DIMENSION(4):: ans
  REAL:: delta, deltam
  REAL:: fnde, fnde3
  INTEGER:: i
  INTEGER:: ih,ih1
  INTEGER,PARAMETER:: itDummy=30, jtDummy=8

  REAL,PARAMETER:: PZ = 2116.2
  REAL,PARAMETER,DIMENSION(6,5):: WF = RESHAPE( (/ &
    -6053.414, 3158.016,5070.555,-1103.562,-889.9834, 332.9619,  &
    -4155.113, 2862.128,3424.648,-583.6367,-632.5012, 239.0175,  &
    -617.3311,-743.4917,3462.835, 190.6228,-750.6321, 214.6133,  &
    -230.7649,-428.9509,2437.024, 641.0098,-678.4924, 169.6867,  &
     994.3538,-1921.537,2557.642, 1251.383,-1011.685, 230.7005 /), (/6,5/) )

   REAL,PARAMETER:: TZ = 288.0

!----------------------------------------------------------------------------


  CALL AT62(h, ans)
  tempk = ans(3) + dtempk          ! static temperature
  t2=tempk *(1.0 + 0.2*makno**2)   ! total temperature
  tempa = t2 - 273.15
  deltam = ans(2)/PZ
  delta = deltam*(1.+.2*makno**2)**3.5
  theta = tempk /TZ

!     Table 24J011 given total air temp, alt, look up max epr
  eprmx = Value2(temp, 13, 13, tempa, alt24, 5, 5, h, eprmax,0)

  IF (incruz .EQ. 1) eprmx = eprmx - 0.1
  IF (eprx .GT. eprmx) eprx= eprmx

!     Table 18L001 given epr, mach no., look up fn/de
  fnde  = Value2(epr, 14, 14, eprx, machno, 10, 10, makno, fnmax,0)
  fnde3= 3.*fnde
  thrst = fnde3*deltam
!
! following statement replaced  -  RLC  -  27Oct99
!      GO TO (101, 102), MFGR

  SELECT CASE(mfgr)
    CASE(1)                        !     PRATT WHITNEY CURVES FOR FDOT
      IF (h.GT.35000.) GO TO 115
      IF (makno.LT..8) GO TO 111
      wfc = Polye1(eprx, 6, wf(1:6,3))
      GO TO 114
  111 IF(makno.LT. .4) GO TO 112
      l1 = 2
      l2 = 3
      pfm = (makno - .4) / .4
      GO TO 113
  112 l1 = 1
      l2 =2
      pfm = makno / .4
  113 wf1 = Polye1(eprx, 6, wf(1:6,l1))
      wf2 = Polye1(eprx, 6, wf(1:6,l2))
      wfc = wf1 + pfm * (wf2 - wf1)
  114 IF(h.LE.25000.) GO TO 121
      wf1 = wfc
      wf2 = Polye1(eprx, 6, wf(1:6,4))
      pfhh = (h - 25000.) / 10000.
      GO TO 120
  115 wf1 = Polye1(eprx, 6, wf(1:6,4))
      wf2 = Polye1(eprx, 6, wf(1:6,5))
      pfhh = (h - 35000.) / 10000.
  120 wfc = wf1 + pfhh * (wf2 - wf1)
  121 kc =.00223181*tempa + .9675897
   82 fdot= 3.*  wfc*delta*kc

    CASE(2)                            !     BOEING TABULATED DATA FOR FDOT
                     !     TABLE 27273A GIVEN THRUST, MACH NO, LOOK UP TSFC
      CALL Serch1(alt27, h, 9,    pfh, ih, limit)
      ih1 = ih + 1
!      J = MAXT(IH)
!      J1 = MAXT(IH1)
!      K= MAXMAK(IH)
!      K1 = MAXMAK(IH1)


    i=1
    tsfc1 =Value2(thrust(:,IH), itdummy, MAXT(ih), fnde3, &
        MACH27(:,ih), jtDummy, MAXMAK(IH), makno, tsfc(1:,1:,iH), i)
    tsfc2 =Value2(thrust(:,ih1), itDummy, MAXT(IH1), fnde3, &
        MACH27(:, ih1), jtDummy, MAXMAK(ih1), makno, tsfc(:, :, ih1), i)

    fdot = (tsfc1 + pfh*(tsfc2-tsfc1))*thrst*SQRT(theta)
  END SELECT
  RETURN
END Subroutine Engepr   ! ---------------------------------------------------





!+
SUBROUTINE Eprio(xvar,yvar,icol,irow,fmtinp,fmtout,iprint,kprint,lprint)
! ---------------------------------------------------------------------------
! PURPOSE -
!   Reads from 9 (two places)
! called by Cpmep1

IMPLICIT NONE
  REAL,INTENT(OUT),DIMENSION(:):: xvar
  REAL,INTENT(OUT),DIMENSION(:,:):: yvar
  INTEGER,INTENT(IN):: icol,irow
  CHARACTER(LEN=*),INTENT(IN):: fmtInp, fmtOut
  INTEGER,INTENT(IN):: iPrint,kPrint,lPrint

!!!      DIMENSION XVAR(IROW),YVAR(IROW,ICOL),HDING(20),LABLE(20)

  CHARACTER(LEN=80):: heading, label
  INTEGER:: i,j
!----------------------------------------------------------------------------
  write(*,*) "entering eprio"
  write(3,*) "entering eprio, icol, irow =", icol, irow
  write(3,*) "in eprio, print codes=", iprint, kprint, lprint

  IF (kprint==0) THEN
    DO i = 1,irow
      READ(9,FMTINP) xvar(i),(yvar(i,j),j=1,icol)
      IF (iprint/=0) WRITE(6,FMTOUT) xvar(i),(yvar(i,j),j=1,icol)
    END DO
    RETURN
  END IF

  IF (lprint/=0) READ(9,'(A)') heading
  READ(9,'(A)') label

  IF (IPRINT/=0) THEN
    IF (LPRINT/=0) WRITE(6,*) " "//heading
    WRITE(6,*) " "//label
  END IF

  READ(9,FMTINP) (xvar(j),j=1,icol)
  IF (IPRINT/=0) WRITE(6,FMTOUT)(xvar(j),j=1,icol)

  RETURN
END Subroutine Eprio   ! ----------------------------------------------------

!+
SUBROUTINE IdleCharacteristics(h,mach, t,ff)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the thrust and fuel flow at idle for a given
!   altitude-mach combination
USE Extra,ONLY: Value2
IMPLICIT NONE

  REAL,INTENT(IN):: h
  REAL,INTENT(IN):: mach
  REAL,INTENT(OUT):: t   ! idle thrust in pounds
  REAL,INTENT(OUT):: ff  ! idle fuel flow in pounds per hour

  INCLUDE 'boeing.inc'
!      COMMON/BOEING/ALT(10),MACHNO(10),FNIDL(10,10),WFIDL(10,10)        &
!     &,ALT24(7),TEMP(13),EPRMAX(13,7),ALTFF(10)                         &
!     &,EPR(14),FNMAX(14,10)                                             &
!     &,DECL(21),DECD(21,10),MACH24(10)                                  &
!     &,ALT27(9),THRUST(30,9),MACH27(8,9),TSFC(30,8,9)                   &
!     &,MACHN1(10),MACHN2(10)

!----------------------------------------------------------------------------

  t = Value2(alt, 10, 10, h, machno, 10, 10, mach, fnidl, 0)*3.
  ff = Value2(altff, 10, 10, h, machno, 10, 10, mach, wfidl,0)*3.

  RETURN
END Subroutine IdleCharacteristics   ! --------------------------------------

END Module Engine   ! =======================================================


!+
      PROGRAM OptimumTrajectory
! ---------------------------------------------------------------------------
! PURPOSE - 

! AUTHORS -  ?

!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!     ?    1.0    ?    Release of program ARC-11282 to COSMIC
! 12Oct99  1.1   RLC   Converted to F90 free format. Replaced alternate
!                      returns in routines PCCOMP and TRIM1 with returnCode
!                      variable. Format statements used as input to
!                      subroutine EPRIO changed to character variables.
! 13Oct99  1.11  RLC   Added variables to common /IO/ in subs Cruzop,Wlefhv
!                      Added variables to common /DRAGMN/ in sub UpDown
!                      FORMAT: no strings continues to new line; comma after
!                      strings, eliminated lots of Hollerith strings (not all)
!                      added log file on unit 3
! 18Oct99  1.12  RLC   New PolyE1, debugging statements
! 22Oct99  1.15  RLC   Continued debugging
! 27Oct99  1.16  RLC   Restructured a bit

USE Extra
USE Engine

      IMPLICIT REAL(M)


      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, THSMAX,MFGR

  REAL:: lambda,e,edot,ff
  INTEGER:: igrnd,iths
  REAL:: mnvtas,mxvtas, vtas2
      COMMON/ENERGY/LAMBDA,E,EDOT,FF,IGRND,ITHS,MNVTAS, MXVTAS,VTAS2
      COMMON/EPSILN/EPSIL1, EPSIL2, ISPLMT

  INTEGER:: iprint, idrag
      COMMON/III/IPRINT, IDRAG

      COMMON/IO/ ANS(4), WS(11), EOPTS(11), HSTARS(11), MSTARS(11),    &
     & PISTRS(11), LAMBS(11), VTASOP(11), FUELFL(11), CRDIST(11),       &
     &CRTIME(11), MNCOST(2), MOPIAS(2), MOPTAS(2), MACHOP(2), FDTOPT(2),&
     & OPTALT(2), HOPT(2), EPRS(2),IWMAX,IOPARM,HTO, VTO, HOLNDG,VOLNDG,&
     &  ETO
      COMMON/OPT/OPTMAK,IOPT,EDTMX
      COMMON/VTRCRU/WSS(10), JJCRUZ(10), LLCRUZ(50,10),EECRUZ(50,10),   &
     & DLLDEE(2, 10),IWMAXX, WTO, WCRUZ ,CRUZCT,HHCRUZ(50,10),          &
     & FFCRUZ(50,10),JLAST1, JLAST2
      COMMON/WINDY/IWIND, PSIA, VWA
      COMMON/CLIMB/MACH,D,ICLIMB,TDUMMY
      COMMON/COST/EPR,ICOST,FC,TC,DTEMPK,W,FUELDT
      COMMON/CRURNG/ECRUZ,FCRUZ,HCRUZ,IL1,IL2,IW,PFW

!      DATA FS2KNT, G, RD2DEG, RHOSL/1.68781,32.2,57.296,.0023769/
!
!
  REAL:: epsil1,epsil2
  INTEGER:: errCode
  CHARACTER(LEN=132):: fileName
!  REAL,PARAMETER:: FS2KNT = 1.68781   ! ft/sec to knots (never used)
!  REAL,PARAMETER:: G = 32.2   ! acceleration of gravity ( never used)
  INTEGER:: IBTABL
  INTEGER:: IOPARM
  INTEGER:: ISPLMT
  INTEGER:: ITAB
  INTEGER:: IWIND
  REAL:: LAMBS
  REAL:: LAMBDX
  REAL:: LLCRUZ
  INTEGER:: MFGR
!  REAL,PARAMETER:: RD2DEG = 57.296   ! radians to degrees (never used)
  INTEGER:: returnCode
!  REAL,PARAMETER:: RHOSL = 0.0023769  ! sea level density, slugs /cu.ft.
  CHARACTER(LEN=*),PARAMETER:: VERSION = "1.16 (27Oct99)"
!----------------------------------------------------------------------------

  WRITE(*,*) "Optimum Trajectory Program, version "//VERSION

! open data files
  DO
    WRITE(*,*) "Enter the name of the general input file: "
    READ(*,'(A)') fileName
!    fileName='case1.in5'
    IF (Len_Trim(fileName)==0) STOP
    OPEN(UNIT=1, FILE=fileName, STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) "Unable to open "//Trim(fileName)//" for input. Try again."
  END DO

  DO
    WRITE(*,*) "Enter the name of the input 7 file: "
    READ(*,'(A)') fileName
!    fileName='case1.in7'
    IF (Len_Trim(fileName)==0) STOP
    OPEN(UNIT=7, FILE=fileName, STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) "Unable to open "//Trim(fileName)//" for input. Try again."
  END DO

  DO
    WRITE(*,*) "Enter the name of the input 8 file: "
    READ(*,'(A)') fileName
!    fileName='case1.xxx'
    IF (Len_Trim(fileName)==0) STOP
    OPEN(UNIT=8, FILE=fileName, STATUS='OLD', &
      IOSTAT=errCode, ACTION='READWRITE', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) "Unable to open "//Trim(fileName)//" for input. Try again."
  END DO

  DO
    WRITE(*,*) "Enter the name of the input 9 file: "
    READ(*,'(A)') fileName
!    fileName='case1.dat'
    IF (Len_Trim(fileName)==0) STOP
    OPEN(UNIT=9, FILE=fileName, STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) "Unable to open "//Trim(fileName)//" for input. Try again."
  END DO

  DO
    OPEN(UNIT=6, FILE='opttraj.out', STATUS='REPLACE', &
      IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) "Unable to open opttraj.out for output. Aborting."
    STOP
  END DO

  DO
    OPEN(UNIT=3, FILE='opttraj.log', STATUS='REPLACE', &
      IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) "Unable to open opttraj.log for output. Aborting."
    STOP
  END DO


  READ(1,'(20I4)') itab,iprint,ioparm, iwind,mfgr,ibtabl,isplmt
  WRITE(3,*) itab,iprint,ioparm,iwind,mfgr,ibtabl,isplmt

  CALL Cpmep1(ibtabl)
  write(3,*) "back from cpmep1"

  IF(ioparm==0) THEN
    WRITE(6,*) ' OPTIMIZING OVER V ONLY'
  ELSE     
    WRITE(6,*) ' OPTIMIZING OVER BOTH V AND PI'
  END IF

  SELECT CASE(mfgr)
    CASE(1)
      WRITE(6,*) "FUEL FLOW RATE FROM PRATT WHITNEY CURVES"
    CASE(2)
      WRITE(6,*) "FUEL FLOW RATE FROM "// &
                  "BOEING TABULATED DATA TABLE CT NO 27273A"
  END SELECT

  IF (itab==0) CALL CRUZOP()
  write(3,*) "back from cruzop"

  CALL CRUTBL(itab,eoptmx,hstar, wcruz, mstar, thstar,lambda,range, &
    iprint, ecruz, ef)
  write(3,*) "back from crutbl"

  lambdx = lambda
  epsil2 = .162*lambdx
  epsil1 = epsil2*2.

  CALL PAGE()
  wto = w
  ipc = 1
  pcmax1 =(llcruz(10, iw)/lambs(iw) -1.)*100.
  pcmax2 =(llcruz(10, iw+ 1)/lambs(iw+1) - 1.)*100.
  pc = MIN(50.0, pcmax1, pcmax2)
!     PC = 50.
  x3 = 1. /range
  linear = 0
  ispliz = 1

  write(3,*) "starting the 190 loop in main"
! looks as IF you will never get out of this loop...
  190 wcruz = wto - fulest(range, ioparm, eto, fc, tc, pc, ecruz, ef)
      iclimb=1
      write(3,*) "calling updown 1"
      CALL UpDown(ef, 1, wto, wlndg, ecruz, wcruz)
      iclimb=2
      CALL Watest(ecruz, lambda, wlndg, 1, wcruzf, range,0,pc,ioparm)
      CALL UpDown(ef, 0, wto, wlndg, ecruz, wcruz)
  350 CALL Voptrj(wlndg, tdist,pc , lambda,ispliz, range,ioparm,ef,1)
      CALL UpDown(ef, 1, wto, wlndg, ecruz, wcruz)
      CALL Voptrj (wlndg, tdist,pc, lambda,ispliz, range,ioparm,ef,2)
      CALL PcComp(pc, ipc, x3,  range, tdist, linear, returnCode,       &
     & ioparm, ispliz)
      write(3,*) "return from PcComp, returnCode=", returnCode
      IF (returnCode==1) GO TO 350

      GO TO 190

      END Program OptimumTrajectory   ! -------------------------------------

!!CONTAINS

!+
SUBROUTINE Cdrag(mach, cl, gear, df,cd)
! ---------------------------------------------------------------------------
! PURPOSE - drag coeffucient

USE Extra, ONLY: PolyE1,Serch1
  REAL,INTENT(IN):: mach
  REAL,INTENT(IN):: cl
  INTEGER,INTENT(IN):: gear
  REAL,INTENT(IN):: df
  REAL,INTENT(OUT):: cd

  INTEGER:: iflap
  REAL:: fflaps
  COMMON/FLAP_SETTINGS/ iflap,FFLAPS

!      DIMENSION CD223(3,7) ! CD224(5, 9) !,CD10(3,6)
  REAL,PARAMETER,DIMENSION(3,7):: CD223 = RESHAPE( (/ &
    0.017681,  -.004785, .057098,   &
    0.05534,   -.06353,  .085691,   &
    0.084281,  -.083402, .083985,   &
    0.09862,   -.084673, .079563,   &
    0.133276,  -.080345, .071928,   &
    0.145686,  -.021176, .048164,   &
    0.17858,   .005597,  .041766/), (/3,7/) )

  REAL,PARAMETER,DIMENSION(5,9):: CD224 = RESHAPE( (/              &
    1.694372E-2, 1.251969E-2, 8.321553E-5, 6.333011E-2, 0.,        &
   -2.137362E-5, 2.122186E-3,-1.415352E-2, 3.065697E-2, 0.,        &
    1.434742E-3,-1.360926E-2, 3.045387E-2, 1.264659E-2, 0.,        &
   -2.396688E-3, 4.467819E-2,-2.118082E-1, 3.233827E-1, 0.,        &
   -3.356013E-3, 6.753743E-2,-3.372280E-1, 5.439879E-1, 0.,        &
   -1.555442E-3, 3.311260E-2,-4.547477E-2,-5.356485E-1, 1.475826,  &
   -2.261359E-3, 7.315081E-2,-4.401166E-1, 8.555356E-1, 0.,        &
   -3.643833E-3, 1.446852E-1,-9.749968E-1, 2.11356,     0.,        &
    9.698153E-3,-4.923731E-2,-9.285760E-2, 7.444338E-1, 0./), (/5,9/) )

  REAL,PARAMETER,DIMENSION(3,6):: CD10 = RESHAPE( (/ &
    0.027126,  -.008443, .003171,    &
    0.027126,  -.008443, .003171,    &
    0.023293,  -.006582, .002594,    &
    0.023345,  -.010857, .003687,    &
    0.020875,  -.01045,  .0035,      &
    0.010316,  -.002858, .00985 /), (/3,6/) )
  REAL,PARAMETER,DIMENSION(3):: CD11 = (/ 0.02781, -0.040903, 0.054867/)
  INTEGER:: i1,i2
  INTEGER,PARAMETER,DIMENSION(9):: IPOLY = (/ 4,4,4,4,4,5,4,4,4 /)
  REAL,PARAMETER,DIMENSION(9):: MVAL = &
    (/ 0.7, 0.76, 0.8, 0.82, 0.84, 0.85, 0.86, 0.88, 0.9 /)
!      DATA IPOLY/5*4, 5, 3*4/
!----------------------------------------------------------------------------

!     with flaps use 2.2-3; without flaps 2.2-4
   99 IF( df .EQ. 0.) GO TO 110
      i2 = iflap +1
      cdbsc1 = Polye1(cl, 3, CD223(1:3,iflap))
      cdbsc2 = Polye1(cl, 3, CD223(1:3, i2))
      cdbasc = cdbsc1 + fflaps*(cdbsc2 - cdbsc1)
      GO TO 113

!     no flaps: M .le. .7, de.cd=0, m gt .7 interpolate; 1st curve
!     CD(CL) basic; 2ND TO 9TH CURVE DE.CD(CL)
  110 IF (mach .LE. .7) GO TO 111
      IF (mach .GE. .9) GO TO 114
      CALL Serch1(mval, mach, 9, fmach, i, limit)
      i2 = i +1
      IF( mach .GT. .76) GO TO 115
      cdbsc1 = 0.
      GO TO 116
  115 cdbsc1 = Polye1(cl, ipoly(i),  CD224(1:5,i))
  116 cdbsc2 = Polye1(cl, ipoly(i2), CD224(1:5, i2))
      decd   = cdbsc1 +  fmach*(cdbsc2 - cdbsc1)
      GO TO 112

  111 decd = 0.
      GO TO 112

  114 decd = Polye1(cl, ipoly(9), CD224(1:5,9))

  112 cdbasc = Polye1(cl, ipoly(1), CD224(1:5, 1)) + decd


  113 IF( gear .EQ. 0) GO TO 119
      IF(df .EQ. 0.0) GO TO 117

      i1 = iflap
      IF( df .LT. 5.)i1 = 1
      i2 = i1 + 1
      cdgr1 = polye1(cl, 3, CD10(1:3, i1))
      cdgr2 = Polye1(cl, 3, CD11(:))
      cdgear = cdgr1 + fflaps*(cdgr2 -cdgr1)
      GO TO 118

  117 cdgear =Polye1(mach, 3, CD11(1:3))
      GO TO 118

  119 cdgear  = 0.0


  118 cd = cdbasc + cdgear

  RETURN
END Subroutine Cdrag   ! ----------------------------------------------------

!+
SUBROUTINE Clift(mach, h, alphap, df, gear, cl,itrim)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: PolyE1

  REAL,INTENT(IN):: mach
  REAL,INTENT(IN):: h
  REAL,INTENT(IN):: alphap
  REAL,INTENT(IN):: df      ! flap angle, degrees
  INTEGER,INTENT(IN):: gear
  REAL,INTENT(OUT):: cl
  INTEGER,INTENT(IN):: itrim
  
  INTEGER:: iflap
  REAL:: fflaps
  COMMON/FLAP_SETTINGS/ iflap, FFLAPS
      DIMENSION CL216(3,7), CL217(5), CL219(3,3)
      DIMENSION CL17(3,7),  CL218(7,3),L218(3)

      COMMON/TRIM2/CL0, CLDA,CLGEAR


! (params)
  REAL,PARAMETER,DIMENSION(7):: FLAPS = (/0., 2., 5., 15., 25., 30., 40./)
  REAL,PARAMETER,DIMENSION(3):: ALT = (/0.0, 20000.0, 40000.0/)
      DATA L218/7,5,4/
      DATA CL216/                       .012578 , .101778,  -.001426,   &
     &    .06282,   .105492,  -.000959, .133972,  .109409,  -.000737,   &
     &    .297618,  .117523,  -.000866, .604512,  .110909,  -.000483,   &
     &    .891152,  .108086,  -.000499, 1.207021, .106244,  -.000104/
      DATA CL217/   1.129635, 1.044612, -9.53392, 20.1947, -13.5972/
      DATA CL218/   .017578,  -.12377,  1.72797,  -7.56877, 6.78399,    &
     &    18.7373,  -28.9716,                                           &
     &    .017364,  .017688,  -.172465, .369741, -.292488,  2*0.,       &
     &    .017732,  .003925,  .02681,   -.065259,3*0./
      DATA CL219/                       .088103,  -.022472, -.002787,   &
     &    .088353,  -.018279, .019896,  .08963,   -.033409, .063451/
      DATA CL17/                        .036537,  -.000448, -.000122,   &
     &    .036537,  -.000448, -.000122, .018362,  .003372,  -.000329,   &
     &    .012889,  -.002635, .000025,  -.012501, -.00536,  .000615,    &
     &    -.013979, -.002833, .00007,   -.039666, -.002491, -.000002/

  REAL,PARAMETER,DIMENSION(3):: CL18 = (/ 0.127893,  -0.496562, 0.576844/)

  REAL:: clbasc
!----------------------------------------------------------------------------
  clbasc=0.0    ! added by RLC

  100 DO 200 i=1,7
      IF (dF .GE. FLAPS(I)) GO TO 200
      GO TO 201
  200 END DO

  201 i = i - 1
      i2 = i + 1
      fflaps= (df - FLAPS(I))/(FLAPS(I+1) - FLAPS(I))
      DO 202 j= 1,3
      IF(h .GE. ALT(J)) GO TO 202
      GO TO 203
  202 END DO

  203 j = j - 1
      j2 = j + 1
      falt = (h - ALT(j)) /(ALT(j2) - ALT(j))

      IF (itrim .EQ. 1) GO TO 204
      clbsc1 = Polye1(alphap, 3, CL216(1:3,I))
      clbsc2 = Polye1(alphap, 3, CL216(1:3,I2))
      clbasc = clbsc1 +fflaps *(clbsc2  - clbsc1)

      IF ( df .NE. 0. ) GO TO 204
      clbmax = Polye1(mach, 5, CL217(1:5))
      IF( clbasc .GT. clbmax) clbasc = clbmax

  204 cl01 = Polye1(mach, L218(J ), CL218(1:7,J))
      Cl02 = Polye1 (mach, L218(J2), CL218(1:7,J2))
      cl0 = cl01 + falt*(cl02-cl01) - .0175

      clda1= Polye1(mach, 3, CL219(1:3, j))
      clda2 = Polye1(mach, 3, CL219(1:3, j2))
      clda = clda1 + falt*(clda2 - clda1) -.088

      IF (gear .EQ. 0) GO TO 229
      IF( df .LT. 2.) GO TO 227
      clgr1 = Polye1(alphap,3,CL17(I:i+2,1))   ! can this be right??
      clgr2 = Polye1(alphap,3,CL17(1:3,I2))
      clgear = clgr1 +fflaps*(clgr2 - clgr1)
      GO TO 230

  227 IF (mach .GT. .4  ) GO TO 228
      clgear = .02
      GO TO 230
  228 clgear = POLYE1(mach, 3, CL18(1:3))
      GO TO 230
  229 clgear = 0.

  230 ialt = j                             ! never used ???
      iflap = i
      cl = clbasc + cl0 + clda*alphap + clgear

  RETURN
END Subroutine Clift   ! ----------------------------------------------------

!+
SUBROUTINE CRUTBL(ITAB, EOPT, HSTAR, WCRUZ, MSTAR, THSTAR,LAMBDA, &
     & RANGE,  IPRINT, ECRUZ, EF)
! ---------------------------------------------------------------------------
! PURPOSE -
!   reads from 8

USE Extra,ONLY: At62,Iclock,Page,Serch1,SerchD

  IMPLICIT REAL(m)
  INTEGER,INTENT(IN):: itab
  REAL,INTENT(IN OUT):: eopt
  REAL,INTENT(IN OUT):: hstar
  REAL,INTENT(IN OUT):: wcruz
  REAL,INTENT(IN OUT):: mstar
  REAL,INTENT(IN):: thstar
  REAL,INTENT(IN OUT):: lambda
  REAL,INTENT(IN OUT):: range
  INTEGER,INTENT(IN):: iprint
  REAL,INTENT(IN OUT):: ecruz
  REAL,INTENT(IN OUT):: ef

! NOTE - thstar not used ????
!      IMPLICIT REAL(M)
      REAL LAMBS,LLCRUZ
      COMMON/COST/EPR, ICOST, FC, TC, DTEMPK, W, FUELDT
      COMMON/IO/ ANS(4), WS(11), EOPTS(11), HSTARS(11), MSTARS(11),    &
     & PISTRS(11), LAMBS(11), VTASOP(11), FUELFL(11), CRDIST(11),       &
     &CRTIME(11), MNCOST(2), MOPIAS(2), MOPTAS(2), MACHOP(2), FDTOPT(2),&
     & OPTALT(2), HOPT(2), EPRS(2),IWMAX,IOPARM,HTO, VTO, HOLNDG, VOLNDG&
     & ,ETO
      COMMON/VTRCRU/WSS(10), JJCRUZ(10), LLCRUZ(50,10),EECRUZ(50,10),   &
     & DLLDEE(2, 10),IWMAXX, WTO, WCRUZ1,CRUZCT,HHCRUZ(50,10),          &
     & FFCRUZ(50,10),JLAST1, JLAST2
      COMMON/WINDY/IWIND, PSIA, VWA
!      DIMENSION T120(1), T136(1), T150(1)   ! not used ???

  REAL,PARAMETER:: FS2KNT=1.68781, G=32.2, RHOSL=0.0023769
!----------------------------------------------------------------------------
  write(3,*) "entering Crutbl"

      IF (itab .EQ. 0) GO TO 200
      write(3,*) "reading from 8 in crutbl 1"
      READ(8,400)  fc, tc, dtempk, psia
  400 FORMAT(8E15.7)
   20 FORMAT(20I4)
   26 IF (iwind==0) THEN
        WRITE(6,*) 'No wind run'
        vwa = 0.
      ELSE
        CALL WindIn()
        WRITE(6, 31) PSIA
   31 FORMAT(/'Aircraf heading = ', F10.0, 2X, 'deg')
      END IF

      WRITE(6, 25) fc, tc, dtempk
   25 FORMAT(/'FUEL COST($/#)', F9.4, 2X,'TIME COST($/HR)= ', F7.2, &
        2X, 'TEMP VAR (DEG K) =', F7.2)

      iw = 0

   85 write(3,*) "reading from 8 in crutbl 85"
      READ(8, 400) w
      IF( w .LT. 0.) GO TO 402          ! jumps out of the big loop
      iw = iw+ 1
      wss(iw) = w
      IF(iprint .EQ. 0) GO TO 92
      WRITE(6, 24) w
   24 FORMAT(/'Aircraft cruise wt = ' F10.0, 2X, '#S')
      WRITE(6, 91)
   91 FORMAT(/3X, 'ALT', 3X, 'MN DRAG  MAX', 9X, &
        '****  MINIMUM COST/DISTANCE  ****', 7X, &
        '**** MINIMUM FUEL/TIME   ****'/3X, 'FT' 6X, &
        'SPEED  SPEED', 10X, 'SPEED' 11X, 'PWR SETG FUEL', &
        14X, 'SPEED' 8X, 'PWR SETG  FUEL'/12X, 'KIAS  KIAS' 6X, &
        'KIAS  KTAS   MACH', 5X, ' EPR', 4X, '$/NM', 6X, &
        'KIAS   KTAS   MACH', 4X, ' EPR', 5X, '#/HR')
   92 fdtopz= 0.   ! never used ???
      ez = 0.      ! never used ???
      jjcruz(iw) = 0


  390 write(3,*) "reading from 8 in crutbl 390"
      READ(8, 400) h, makias, fbias, (mopias(i), moptas(i), machop(i), &
     & eprs(i), fdtopt(i), i=1,2),fueldt
      IF ( h .LT. -10.) GO TO 401    ! jumps out of this loop
      e = h + (moptas(1) *FS2KNT)**2/64.4
      IF(iprint .EQ. 0) GO TO 93
      WRITE(6, 100) h, makias, fbias, (mopias(i), moptas(i), machop(i), &
     & eprs(i), fdtopt(i), i=1,2),e
  100 FORMAT(1X, 3F7.0, 2X, 2F7.0, F7.3, 2X, F8.3, F7.2, 2X, 2F7.0, F7.3,&
     & 2X, F8.3, F7.0, F10.0)
   93 jjcruz(iw) = jjcruz(iw) + 1
      jj= jjcruz(iw)
      ffcruz(jj, iw) = fueldt
      hhcruz(jj, iw) =h
      llcruz(jj, iw) = fdtopt(1)
      eecruz(jj, iw) = e
      fdtopz= fdtopt(1)
      ez = e
      GO TO 390


  401 write(3,*) "reading from 8 in crutbl 401"
      READ(8, 400) hopt(1), optmak, optias, opttas, mncost(1), epr,eopt&
     & , fueld1
  102 IF(iprint .EQ. 0) GO TO 94
      WRITE(6, 105)
  105 FORMAT(/'MINIMIZING FUEL/DISTANCE:')
      WRITE(6, 101) hopt(1), optmak, optias, opttas, mncost(1), epr,eopt
  101 FORMAT (/' OPT ALT= ', F7.0, 1X,'FT, OPT SPEED = ', F7.4, 'MACH,',  &
       F7.0, ' KIAS, ', F7.0, 'KTAS, MIN (FDOT/V) =', F7.3, &
       ' $/NM, CRUISE POWER SETG = ', F7.4, ' EPR'/ &
       ' OPTIMUM CRUISE ENERGY = ', F8.0, 2X, 'FT')

   94 CALL Serch1 (hhcruz(:, iw),hopt(1),jj,pf,jjopt, limit)
      jjopt1 = jjopt + 1
      jsum = jj + jjopt1
      DO 391 i= jjopt1, jj
      j = jsum - i
      hhcruz(j+1, iw) = hhcruz(j, iw)
      ffcruz(j+1, iw) = ffcruz(j, iw)
      llcruz(j+1, iw) = llcruz(j, iw)
      eecruz(j+1, iw) = eecruz(j, iw)
  391 END DO
      llcruz(jjopt1, iw) = mncost(1)
      eecruz(jjopt1, iw) = eopt
      hhcruz(jjopt1, iw) = hopt(1)
      ffcruz(jjopt1, iw) = fueld1
      jjcruz(iw) = jjcruz(iw) + 1

      write(3,*) "reading from 8 in crutbl 400"
      READ(8, 400) hopt(2), optmak, optias, opttas, mncost(2), epr
  103 IF(iprint .EQ. 0) GO TO 85
      WRITE(6, 106)
  106 FORMAT(/'MINIMIZING FUEL/TIME:')
      WRITE(6, 104) hopt(2), optmak, optias, opttas, mncost(2), epr
  104 FORMAT(/'OPT ALT= ', F7.0, ' FT, OPT SPEED = ', F7.4, ' MACH,', &
        F7.0, ' KIAS, ', F7.0, 'KTAS, MIN (FDOT) = ', F7.0, &
        '#/HR, CRUISE POWER SETG =', F7.4, ' EPR')
      CALL PAGE
      GO TO 85     ! end of the big loop

  402 write(3,*) "reading from 8 in crutbl 402"
      READ(8, 20) iwmax
      iwmaxx= iwmax
      READ(8, 400) (dlldee(1,j), dlldee(2,j), j=1, iwmax)

      DO 197 i= 1, iwmax
      write(3,*) "reading from 8 in crutbl 197"
      READ(8, 400) ws(i), eopts(i), mstars(i), hstars(i), pistrs(i),   &
     & lambs(i), vtasop(i),fuelfl(i)
  197 END DO

      WRITE(6,191)
  191 FORMAT(/' D(LAMBDA)/DE = A E + B IN $/NM**2'/' CRUISE WT', 8X, &
     & 'A', 14X, 'B')

      WRITE(6, 192) (ws(i), dlldee(1,i), dlldee(2,i), i=1,iwmax)
  192 FORMAT(/F10.0, 2E15.7)

  200 READ(1, 21) w, range
      READ(1, 21) hto, vto, holndg, volndg
   21 FORMAT(8F10.2)
      wto = w
      CALL AT62(hto, ans)
      rho=ans(1)*(1.+ dtempk/ans(3))
      vo =vto*SQRT(RHOSL/rho)*FS2KNT
      eto=hto+ .5*  (vo)**2/G

      pcDummy=0.0
      wcruz = wto - Fulest(range, ioparm, eto, fc, tc, pcDummy, ecruz, ef)
      wcruz1 = wcruz
      CALL Serchd(ws, wcruz, iwmax, pf, iw, limit)
      iwnd1 = iw + 1
      iwcruz = iw
      eopt  = eopts (iw) + pf*(eopts (iw + 1) - eopts (iw))
      hstar = hstars(iw) + pf*(hstars(iw + 1) - hstars(iw))
      mstar = mstars(iw) + pf*(mstars(iw + 1) - mstars(iw))
      eprtar= pistrs(iw) + pf*(pistrs(iw + 1) - pistrs(iw))
      lambda= lambs (iw) + pf*(lambs (iw + 1) - lambs (iw))
      cruzct = lambda
      fueldt = fuelfl(iw) + pf*(fuelfl(iw+1) - fuelfl(iw))

      IF( pf .EQ. 1.) GO TO 195
      isum = iwmax + iwnd1
      DO 220 i= iwnd1,iwmax
      j= isum - i
      ws(j+ 1) = ws(j)
      eopts(j+ 1) = eopts(j)
      hstars(j+ 1) = hstars(j)
      mstars(j+1) = mstars(j)
      pistrs(j+1) = pistrs(j)
       lambs(j+1) = lambs(j)
      fuelfl(j+1) = fuelfl(j)
  220 END DO
      ws(iwnd1) = wcruz
      eopts(iwnd1) = eopt
      hstars(iwnd1) = hstar
      mstars(iwnd1) = mstar
      pistrs(iwnd1) = eprtar
      lambs(iwnd1) = lambda
      fuelfl (iwnd1) = fueldt
      iwmax = iwmax + 1
  195 CALL Page
      iwmaxz = iwmax - 1

      WRITE(6, 198)
  198 FORMAT(/'  CRUISE WT   OPT H      KTAS    OPT MACH    EPR     ', &
       'COST   FUEL FLOW   OPT E'/6X,'LBS',7X,'FT',26X,'SETTING',4X,    &
     &'$/NM',6X,'#/HR',6X,'FT')

  write(3,*) "297 loop with iwmax=", iwmax
      DO 297 i= 1, iwmax
      CALL AT62(hstars(i), ans)
      tempk = ans(3) + dtempk
      a = 65.76 * SQRT(tempk)
      vtasop(i) = mstars(i) * a / FS2KNT
      WRITE(6,196) ws(i),hstars(i),vtasop(i),mstars(i),pistrs(i),       &
     &lambs(i),fuelfl(i),eopts(i)
  297 END DO

  196 FORMAT (/F10.0,F10.0,F10.2,F10.4,F 9.4,F10.3,F10.2,F10.0)

      crdist(1) = 0.
      crtime(1) = 0.
      vwa1 = 0.
      vwa2 = 0.

  DO i= iwnd1, iwmaxz
    IF (iwind/=0) THEN
      CALL Wind(hstars(i),   psia, vwa1)
      CALL Wind(hstars(i+1), psia, vwa2)
    END IF
    vg1 = vtasop(i) + vwa1/FS2KNT
    vg2 = vtasop(i+1) + vwa2/FS2KNT
    avevel = vg1 + vg2
    cruzd =  (ws(i) -ws(i+1))*avevel/(fuelfl(i) + fuelfl(i+1))
    avevel = avevel/2.
    j = i+2 - iwnd1
    crdist(j) = crdist(j-1) + cruzd
    crtime(j) = crtime(j-1) + cruzd*3600./avevel
  END DO

  last = j
  ilast = 0
  IF (range  .GE. crdist(last)) GO TO 222
  CALL Serch1(crdist, range , last,pfc, ilast, limit)
  last = ilast + 1

  222 CALL Page
  WRITE(6,'(A)') ' CRUISE DIST  TIME      WEIGHT    ENERGY'// &
        '  ALTITUDE  MACH NO       KTAS   GRD SPEED   LAMBDA  PWR SETG'
  WRITE(6,'(A)') '      NM    HR:MN:SEC     #         FT         FT'// &
        '                         KNOT      $/NM       EPR'
  dedist = 100.
  cruzd = 0.
  cruzt = 0.
  CALL Iclock (cruzt, ihr, imin, isec)
  IF (iwind .EQ. 0) THEN
    gsprnt = vtasop(iwnd1)
  ELSE
    CALL Wind(hstars(iwnd1), psia, vwa)
    gsprnt = vtasop(iwnd1) + vwa/FS2KNT
  END IF
  
  WRITE(6, 225) cruzd, ihr, imin, isec, ws(iwnd1 ), eopts(iwnd1 ),  &
     & hstars(iwnd1 ), mstars(iwnd1 ), vtasop(iwnd1 ), gsprnt,          &
     & lambs(iwnd1 ), pistrs(iwnd1 )
  225 FORMAT(/F10.2, 2X, 2(I2, ':'), I2, 3F10.0, F10.4, 2F10.2,2F10.3)
  cruzd = cruzd + dedist


  224 CALL Serch1(crdist, cruzd, last, pf, idist, limit) ! returns pf,idist,limit
 2241 cruzt = crtime(idist) + pf *(crtime(idist +1) - crtime(idist))
      CALL Iclock (cruzt, ihr, imin, isec)
      iw = idist + iwcruz
      wt = ws(iw) + pf*(ws(iw+1) - ws(iw))
      eprnt = eopts(iw) + pf*(eopts(iw+1) - eopts(iw))
      hprnt = hstars(iw) + pf *(hstars(iw+1) - hstars(iw))
      mprnt = mstars(iw) + pf*(mstars(iw+1) - mstars(iw))
      vprnt = vtasop(iw) + pf*(vtasop(iw+1) - vtasop(iw))
      costnt = lambs(iw) + pf*(lambs(iw+1) - lambs(iw))
      piprnt = pistrs(iw) + pf*(pistrs(iw+1) - pistrs(iw))
      IF ( iwind .EQ. 0) GO TO 226
      CALL Wind(hprnt, psia, vwa)
      gsprnt = vprnt + vwa/FS2KNT
      GO TO 227
  226 gsprnt = vprnt
  227 WRITE(6, 225)  cruzd, ihr, imin, isec, wt, eprnt, hprnt, mprnt,   &
     & vprnt,gsprnt, costnt, piprnt
      IF( ilast .EQ. -1) RETURN
      cruzd = cruzd + dedist
      cruzd = cruzd + dedist

  IF (cruzd <= range) THEN
    IF (cruzd==range) ilast=-1
    GO TO 224
  END IF

!      IF( CRUZD - RANGE) 224, 2242, 2243           ! obscolescent feature
!     IF( CRUZD  .LT. CRDIST(LAST)) GO TO 224
!     IF ( ILAST .EQ. 0) GO TO 228
!     IDIST = ILAST
!     PF = PFC
! 2242 ILAST = -1
!      GO TO 224

 2243 IF ((cruzd - dedist) .GE. range) RETURN   ! the escape
      cruzd = range
  ilast=-1
  GO TO 224

! 228 IW = LAST + IWCRUZ
!     CALL ICLOCK (CRTIME(LAST), IHR, IMIN, ISEC)
!     IF ( IWIND .EQ. 0) GO TO 215
!     CALL WIND(HSTARS(IW), PSIA, VWA)
!     GSPRNT = VTASOP(IW) + VWA/FS2KNT
!     GO TO 216
! 215 GSPRNT = VTASOP(IW)
! 216 WRITE(6, 225) CRDIST(LAST), IHR, IMIN, ISEC, WS(IW), EOPTS(IW),
!    1 HSTARS(IW), MSTARS(IW), VTASOP(IW), GSPRNT,     LAMBS(IW),
!    2 PISTRS(IW)
!
!     RETURN
END Subroutine Crutbl   ! ---------------------------------------------------

!+
SUBROUTINE CRUZOP()
! ---------------------------------------------------------------------------
! PURPOSE -
!   writes to 8

USE Extra,ONLY: At62,Jtrunc,Lsqpol,MinF,MinF2,Page,Serch1
USE Engine,ONLY: EngEPR
      IMPLICIT REAL(M)
      REAL LAMBS,LLCRUZ
!!!      EXTERNAL Fbound, FCOST, FOPT,FCLIMB,FCLMB2


  REAL,EXTERNAL:: Fbound, Fcost, Fopt, Fclimb, Fclmb2


  REAL:: ans
      COMMON/IO/ ANS(4), WS(11), EOPTS(11), HSTARS(11), MSTARS(11),    &
     & EPRSTR(11), LAMBS(11), VTASOP(11), FUELFL(11), CRDIST(11),       &
     &CRTIME(11), MNCOST(2), MOPIAS(2), MOPTAS(2), MACHOP(2), FDTOPT(2),&
     & OPTALT(2), HOPT(2), EPRS(2),IWMAX,IOPARM,   &
        hto,vto,holndg,volndg,eto

  REAL:: rho,p,tempk,a,ratio,h,alpha,cl,eprmax
  INTEGER:: mfgr
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, EPRMAX,MFGR

  INTEGER:: iprint, idrag
      COMMON/III/IPRINT, IDRAG
      COMMON/COST/EPR, ICOST, FC, TC, DTEMPK, W, FUELDT
      COMMON/OPT/OPTMAK,IOPT,EDTMX

  REAL:: lambda,e,edot,ff
  INTEGER:: igrnd,iths
  REAL:: mnvtas,mxvtas, vtas2
      COMMON/ENERGY/LAMBDA,E,EDOT,FF,IGRND,ITHS,MNVTAS, MXVTAS,VTAS2

  REAL:: mach,d
  INTEGER:: iclimb
  REAL:: tdummy
      COMMON/CLIMB/MACH, D, ICLIMB,TDUMMY

  INTEGER:: iwind
  REAL:: psia,vwa
      COMMON/WINDY/IWIND, PSIA, VWA


      COMMON/VTRCRU/WSS(10), JJCRUZ(10), LLCRUZ(50,10),EECRUZ(50,10),   &
     & DLLDEE(2, 10),IWMAXX, WTO, WCRUZ1,CRUZCT,HHCRUZ(50,10),          &
     & FFCRUZ(50,10),JLAST1,JLAST2

!      DIMENSION WTGFC(50), AP(20,3), BP(20)
!      DATA WTGFC /50*1./

  REAL,DIMENSION(20):: bp
  REAL:: fa,fb
  REAL:: fac,fbc
  REAL,PARAMETER:: FS2KNT = 1.68781
  REAL,PARAMETER:: G = 32.2
!      DATA FS2KNT, G, RD2DEG/1.68781, 32.2, 57.296/  ! RD2DEG not used
      DATA ENDATA/-1.E6/

  REAL:: dragA
  REAL,DIMENSION(50):: wtgfc

!----------------------------------------------------------------------------
  write(3,*) "entering cruzop"

   88 IPRINT = 0    ! arbitrarily changes the global setting !!! ???
      READ(1, 21) FC, TC, DTEMPK, PSIA
      write(3,*) "reading #1 in cruzop", fc,tc,dtempk,psia
   20 FORMAT(20I4)
   21 FORMAT(8F10.2)
  write(3,*) "writing to 8 in cruzop 1"
      WRITE(8,400)  FC, TC, DTEMPK, PSIA
  400 FORMAT(8E15.7)


  IF (IWIND==0) THEN
    WRITE(6,*) "NO WIND RUN (in cruzop)"
    VWA = 0.0
  ELSE
    CALL WINDIN()
    WRITE(6, 31) PSIA
   31 FORMAT(/'AIRCRAFT HEADING = ' F10.0, '  DEG')
  END IF


  WRITE(6, 25) FC, TC, DTEMPK
   25 FORMAT(/'FUEL COST($/#)', F9.4, 2X,'TIME COST($/HR)= ', F7.2, 2X,&
     & 'TEMP VAR (DEG K) =', F7.2)
      READ(1, 21) W,WN, DEW


!     W LOOP
      IW = 1
      EPR = 1.5
   85 WRITE(6, 24) W    ! keep coming back here
      ICEILG = 0
      write(3,*) "writing to 8 in cruzop 2"
      WRITE(8, 400) W
   24 FORMAT(/' AIRCRAFT CRUISE WT = ', F10.0, 2X, '#S')
      FDTOPZ= 0.
      EZ = 0.
      JJCRUZ(IW) = 0
      H = 0.
      MNCOST(1) = 1.E6
      MNCOST(2) = 1.E6
      DEH = 1000.
      WRITE(6, 91)
   91 FORMAT(/3X, 'ALT', 3X, 'MN DRAG  MAX', 9X, &
       '****  MINIMUM COST/DISTANCE  ****', 7X, &
       '**** MINIMUM FUEL/TIME   ****'/3X, 'FT' 6X, &
       'SPEED  SPEED', 10X, 'SPEED', 11X, 'PWR SETG FUEL',14X, &
       'SPEED', 8X, 'PWR SETG  FUEL'/12X, 'KIAS  KIAS', 6X, &
       'KIAS  KTAS   MACH', 5X, 'EPR ', 4X, '$/NM' 6X, &
       'KIAS   KTAS   MACH', 4X, 'EPR ' 5X, '#/HR')
      ALPHA = 0.

   89 IF ( H .LE. 39999.) GO TO 92
      IF((H- 39999.) .GT. 50.) GO TO 135
      H = 39999.
   92 write(3,*) "at line 92 in cruzop, h=", h, ans
      CALL AT62(H, ANS)
      TEMPK = ANS(3) + DTEMPK
      P = ANS(2)
      RHO = P / (3092.40 * TEMPK)
      A = 65.76 * SQRT(TEMPK)
      RATIO = SQRT(RHO/.0023769)
   90  IDRAG = 1
      EPRMAX = 2.4
      MINDRG = MINF(0.1, 0.9, Fbound, MACH, IPRINT) 
      CALL ENGEPR(H, MACH, EPRMAX, 1, TMAX, FF, MFGR)
! following statement replaced   -  RLC 27Oct99
!      IF( TMAX - (MINDRG + 50.)) 129,128,120
  IF (tmax < mindrg+50.0) GO TO 129
  IF (tmax==mindrg+50.0) GO TO 128   ! bet this will never branch

  120 FA =  MACH - .1
      CALL ENGEPR(H, FA, EPRMAX, 1, TMAXA, FFA, MFGR)
      DRAGA = Fbound(FA)
      IF( TMAXA .GE. DRAGA)GO TO 127
      IDRAG = 2
      FAC = MINF(0.1, MACH, Fbound, FA, IPRINT)      ! not used ???
  127 FB = .9
      CALL ENGEPR(H, FB, EPRMAX, 1, TMAXB, FFB, MFGR )
      IDRAG = 1
      DRAGB= Fbound (FB)
      IF( TMAXB   .GE. DRAGB)GO TO 125
  124 IDRAG = 2
      FBC = MINF(MACH, .9, Fbound, FB, IPRINT)      ! not used ???
  125 MAKIAS = 29.*SQRT(P*((1.+.2*MACH*MACH)**3.5-1.))/FS2KNT
      FBIAS = 29.*SQRT(P*((1.+.2*FB*FB)**3.5-1.))/FS2KNT

      IF (IWIND  .EQ. 0) GO TO 123
      CALL WIND(H, PSIA, VWA)
  123 DO 122 I=1,2
      ICOST = I
      FDTOPT(I) = MINF(FA, FB, FCOST, MACHOP(I), IPRINT)
      IF( FDTOPT(I) .GE. MNCOST(I)) GO TO 121
      MNCOST(I) = FDTOPT(I)
      OPTALT(I) = H
  121 MOPTAS(I) =MACHOP(I)*A/FS2KNT
      MOPIAS(I)=29.*SQRT(P*((1.+.2*MACHOP(I)*MACHOP(I))**3.5-1.))/FS2KNT
      EPRS(I) = EPR
      IF( I .EQ. 1) FUELD1 = FUELDT
  122 END DO
      EE= H + (MOPTAS(1) *FS2KNT)**2/64.4
      JJCRUZ(IW) = JJCRUZ(IW) + 1
      JJ= JJCRUZ(IW)
      HHCRUZ(JJ, IW) =H
      FFCRUZ(JJ, IW) =  FUELD1
      LLCRUZ(JJ, IW) = FDTOPT(1)
      EECRUZ(JJ, IW) = EE
      FDTOPZ= FDTOPT(1)
      EZ = EE
      WRITE(6, 100) H, MAKIAS, FBIAS, (MOPIAS(I), MOPTAS(I), MACHOP(I), &
     & EPRS(I), FDTOPT(I), I=1,2),EE
  100 FORMAT(' ',3F7.0, 2X, 2F7.0, F7.3, 2X, F8.3,F7.2,   2X, 2F7.0, &
        F7.3, 2X, F8.3, F7.0,F10.0)
      write(3,*) "writing to 8 in cruzop 3"
      WRITE(8, 400) H, MAKIAS, FBIAS, (MOPIAS(I), MOPTAS(I), MACHOP(I), &
     & EPRS(I), FDTOPT(I), I=1,2),FUELD1
      IF (ICEILG.EQ.0) GO TO 126
      ICEILG = ICEILG + 1
      IF (ICEILG .GT. 5) GO TO 135
      DEH = DEH/2.
  126 H = H + DEH
      GO TO 89
  128 MACHOP(1) = MACH
      MACHOP(2) = MACH
      FDTOPT(2) = FF
      GO TO 135
  129 IF (ICEILG .NE. 0) GO TO 130
      ICEILG = 1
  130 DEH = DEH/2.
      H = H - DEH
      GO TO 89
  135 HCEILG = H
      IF ( HCEILG .GT. 39999.) HCEILG = 39999.
      write(3,*) "writing to 8 in cruzop 4"
      WRITE(8,400)ENDATA,MAKIAS,FBIAS,(MOPIAS(I), MOPTAS(I), MACHOP(I), &
     & EPRS(I), FDTOPT(I), I=1,2)
      DO 136 I=1,2
      ICOST = I
      H1 = OPTALT(I) - 1000.
      H2 = OPTALT(I)  + 1000.
      IF (H2 .GT. HCEILG) H2 = HCEILG
      MNCOST(I) = MINF2(H1, H2, FOPT, HOPT(I), IPRINT)
      H1 = HOPT(I) - 37.
      H2 = HOPT(I) + 37.
      IF (H2 .GT. HCEILG) H2 = HCEILG
      MNCOST(I) = MINF2( H1, H2, FOPT, HOPT(I), IPRINT)
      OPTTAS = OPTMAK*A/FS2KNT
      OPTIAS = 29.*SQRT(P*((1.+.2*OPTMAK*OPTMAK)**3.5-1.))/FS2KNT

! following statement replaced  -  RLC  27Oct99
!      GO TO (102, 103), ICOST
    SELECT CASE(icost)
      CASE(1)
        WRITE(6,*) " "
        WRITE(6,*) 'MINIMIZING FUEL/DISTANCE:'
        EOPT = HOPT(1) + .5*(OPTTAS*FS2KNT)**2/G
        WRITE(6, 101) HOPT(I), OPTMAK, OPTIAS, OPTTAS, MNCOST(I), EPR,EOPT
        write(3,*) "writing to 8 in cruzop 5"
        WRITE(8, 400) &
          HOPT(I), OPTMAK, OPTIAS, OPTTAS, MNCOST(I), EPR, EOPT, FUELDT
  101 FORMAT (/'OPT ALT= ' F7.0, 1X,'FT, OPT SPEED = ' F7.4, 'MACH,' &
     & F7.0, ' KIAS, ' F7.0, 'KTAS, MIN (FDOT/V) =' F7.3, ' $/NM, CRUISE&
     & POWER SETG = ' F7.4, 'EPR '/' OPTIMUM CRUISE ENERGY = ' F8.0,  &
     & 2X, 'FT')
        HSTAR  = HOPT(1)
        MSTAR = OPTMAK
        EPRTAR = EPR
        LAMBDA = MNCOST(1)
        FUELD1 = FUELDT
      CASE(2)
!      GO TO 136
        WRITE(6,*) " "
        WRITE(6,*) 'MINIMIZING FUEL/TIME:'
        WRITE(6, 104) HOPT(I), OPTMAK, OPTIAS, OPTTAS, MNCOST(I), EPR
        write(3,*) "writing to 8 in cruzop  6"
        WRITE(8, 400) HOPT(I), OPTMAK, OPTIAS, OPTTAS, MNCOST(I), EPR
  104 FORMAT (/'OPT ALT= ', F7.0, ' FT, OPT SPEED = ', F7.4, 'MACH,' &
     & F7.0, ' KIAS, ' F7.0, 'KTAS, MIN (FDOT) = ' F7.0, &
       '#/HR, CRUISE POWER SETG =', F7.4, 'EPR ')
    END SELECT
  136 END DO
      WS(IW) = W
      EOPTS(IW) = EOPT
      HSTARS (IW) = HSTAR
      MSTARS(IW) = MSTAR
      EPRSTR(IW) = EPRTAR
      LAMBS(IW) = LAMBDA
      FUELFL(IW) = FUELD1
      WSS(IW) = W
      JJOPT = JTRUNC(LLCRUZ(:, IW), JJ)
      JJOPT1 = JJOPT + 1
      JSUM = JJ + JJOPT1
      DO 391 I= JJOPT1, JJ
      J = JSUM - I
      HHCRUZ(J+1, IW) = HHCRUZ(J, IW)
      FFCRUZ(J+1, IW) = FFCRUZ(J, IW)
      LLCRUZ(J+1, IW) = LLCRUZ(J, IW)
      EECRUZ(J+1, IW) = EECRUZ(J, IW)
  391 END DO
      LLCRUZ(JJOPT1, IW) = MNCOST(1)
      EECRUZ(JJOPT1, IW) = EOPT
      FFCRUZ(JJOPT1, IW) = FUELD1
      HHCRUZ(JJOPT1, IW) = HOPT(1)
      JJCRUZ(IW) = JJCRUZ(IW) + 1
      JJ = JJCRUZ(IW)
      CALL SERCH1(HHCRUZ(:, IW), 20000., JJ, PF, ISTART, LIMIT)
      IF (PF .GT. 0.9) ISTART = ISTART + 1
  393 NH = JJOPT1 - ISTART + 1

! only call to LSQPOL:  (fitting a quadratic)

  wtgfc(:)=1.0
  WRITE(3,*) "Calling LsqPol, istart,jjopt,nh=", istart,jjopt,nh
  CALL LSQPOL(EECRUZ(ISTART:jjopt, IW), LLCRUZ(ISTART:jjopt, IW), wtgfc,bp(1:3))

      DLLDEE(1, IW) = 12160.*BP(3)
      BPRIME        = 6080.*BP(2)
      DLLDEE(2, IW) = - DLLDEE(1,IW)*(EOPT + 1000.)
      W = W - DEW
      IF ( W .LT. WN) GO TO 200
      CALL PAGE
      IW = IW + 1
      IF( IW .GT.10) GO TO 900
      GO TO 85
!
!
  900 WRITE(6, 901) IW
  901 FORMAT(/'DIMENSIONED FOR ONLY TEN DIFFERENT WEIGHTS'/ &
        ' COMPUTING THE ', I2, 'TH WEIGHT')
      GO TO 1000
  200 IWMAX = IW
      write(3,*) "writing to 8 in cruzop 7"
      WRITE(8,400) ENDATA
      WRITE(8, 20) IWMAX
      IWMAXX = IWMAX
      write(3,*) "writing to 8 in cruzop 8"
      WRITE(8, 400) (DLLDEE(1,J), DLLDEE(2,J), J=1, IWMAX)
      CALL PAGE
      WRITE(6, 198)
  198 FORMAT(/'  CRUISE WT   OPT H      KTAS    OPT'// &
       ' MACH    EPR      COST   FUEL FLOW   OPT E'/ &
       5X,'LBS',7X,'FT',26X,'SETTING',4X, '$/NM',6X,'#/HR',6X,'FT')
      DO 197 I= 1, IWMAX
      CALL AT62(HSTARS(I), ANS)
      A = 65.76 * SQRT(ANS(3) + DTEMPK)
      VTASOP(I) = MSTARS(I)* A / FS2KNT
      write(3,*) "writing to 8 in cruzop 9"
      WRITE(8, 400) WS(I), EOPTS(I), MSTARS(I), HSTARS(I), EPRSTR(I),   &
     & LAMBS(I), VTASOP(I),FUELFL(I)
      WRITE(6,196) WS(I),HSTARS(I),VTASOP(I),MSTARS(I),EPRSTR(I),       &
     &LAMBS(I),FUELFL(I),EOPTS(I)
  197 END DO
  196 FORMAT (/F10.0,F10.0,F10.2,F10.4,F 9.4,F10.3,F10.2,F10.0)
 1000 RETURN
END Subroutine CruzOp   ! ---------------------------------------------------


!+
FUNCTION Fbound(MACH) RESULT(fOutput)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Engine,ONLY:EngEPR

  REAL,INTENT(IN):: MACH

  REAL:: fOutput

  REAL:: eprmax
  REAL:: h
  INTEGER:: mfgr
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, EPRMAX,MFGR

  INTEGER:: iprint, idrag
      COMMON/III/IPRINT, IDRAG

      COMMON/COST/EPR, ICOST, FC, TC, DTEMPK,W,FUEL
  REAL:: ff
  INTEGER:: incruz

  REAL:: tmax
!----------------------------------------------------------------------------
  f=0.5*(mach*a)**2*1560.*rho     ! 1560 = wing area ???? always???
  cl = w/f
  CALL Cdrag(mach, cl, 0, 0., cd)
  d = f*cd

  IF (idrag==1) THEN
    fOutput = D
    RETURN
  END IF

  incruz=1
  CALL Engepr(h, mach, eprmax, incruz, tmax, ff, mfgr)
! 152 Fbound = ABS(TMAX - D)

  IF (tmax < d) THEN
    fOutput=1.E6*ABS(tmax - d)
  ELSE
    fOutput= tmax - D
  END IF
  IF (iprint==0) RETURN

  WRITE(6,100) eprmax, tmax, d, fOutput, mach
  100 FORMAT(/'EPRMAX, TMAX, D, ABS(T-D)', F10.2, 3F10.0, F10.5)
  RETURN
END Function Fbound   ! -----------------------------------------------------

!+
FUNCTION FCLIMB(VTAS) RESULT(fOutput)                           ! not a PURE function
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY:AT62,Value2
USE Engine,ONLY: EngEPR,IdleCharacteristics

  REAL,INTENT(IN):: vtas   ! true air speed (units ???)
  REAL:: fOutput

!!!      REAL MACH, MACHNO, MACH24, MACH27, MACHN1
  REAL:: mach

  REAL:: lambda,e,edot,ff
  INTEGER:: igrnd,iths
  REAL:: mnvtas,mxvtas, vdummy ! not vtas2 ???
     
COMMON/ENERGY/LAMBDA,E,EDOT,FF,IGRND,ITHS,MNVTAS,MXVTAS,VDUMMY
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, EPRSET,MFGR

  REAL:: epr
      COMMON/COST/EPR, ICOST, FC, TC, DTEMPK,W,FUEL

  REAL:: t
      COMMON/CLIMB/MACH, D, ICLIMB,T
      COMMON/OPT/OPTMAK,IOPT,EDTMX
      COMMON/WINDY/IWIND, PSIA, VWA
      COMMON/CRUZ/NEARCZ,IDLE

!  INCLUDE 'boeing.inc'

!      COMMON/BOEING/ALT(10),MACHNO(10),FNIDL(10,10),WFIDL(10,10)        &
!     &,ALT24(7),TEMP(13),EPRMX(13,7) ,ALTFF(10)                         &
!     &,EPRS(14),FNMAX(14,10)                                            &
!     &,DECL(21),DECD(21,10),MACH24(10)                                  &
!     &,ALT27(9),THRUST(30,9),MACH27(8,9),TSFC(30,8,9)                   &
!     &,MACHN1(10),MACHN2(10)
!      DIMENSION ANS(4)
!      DATA G, EPRMAX/32.2, 2.4/

  REAL,DIMENSION(4):: ans
  REAL,PARAMETER:: EPRMAX = 2.4
  REAL,PARAMETER:: G = 32.2
  REAL:: h
  INTEGER:: incruz
  INTEGER:: mfgr
  
!----------------------------------------------------------------------------


   80 H = E - .5*  VTAS*VTAS/G
   81 CALL AT62(H, ANS)
      P = ANS(2)
      TEMPK = ANS(3) + DTEMPK
      RHO = P / (3092.40 * TEMPK)
      A = 65.76 * SQRT(TEMPK)    ! speed of sound (units ?)
   82 MACH = VTAS/A
      IF (MACH .LT. .1) MACH = .1
      IF ( MACH .GT. .9) MACH = .9
      F = .5*VTAS**2*1560.*RHO                             ! 1560. is what?
      CL = W/F
      CALL CDRAG(MACH, CL, 0, 0., CD)   ! gear=0.0  df=0.0
      D = F*CD

      IF (IWIND .EQ. 0) GO TO 90
      CALL WIND (H, PSIA, VWA)
  90  IF (ICLIMB .EQ. 1 ) GO TO 91
      IF( IOPT .NE. 1) GO TO 93
!!!      T = VALUE2(ALT, 10, 10, H, MACHNO, 10, 10, MACH, FNIDL, 0)*3.
!!!      FF = VALUE2(ALTFF, 10, 10, H, MACHNO, 10, 10, MACH, WFIDL,0)*3.
      CALL IdleCharacteristics(h,mach, t,ff)   ! computes t and ff
      GO TO 92
   91 EPR = EPRSET


   93 incruz=0
      CALL ENGEPR( H, MACH, EPR,0, T, FF, MFGR)
   92 EDOT = (T-D)*(VTAS)/W
      fOutput= (FC*FF + TC - LAMBDA*(VTAS + VWA))/ ABS(EDOT)
      IF( ABS(EDOT) .LE. 5.) GO TO 95
      IF( NEARCZ .EQ. 1 .AND.  fOutput .GT. 0..AND.ICLIMB.EQ.2)GO TO 95
      IF(ICLIMB .NE. 1 .OR. IOPT .NE. 3) RETURN
      IF ( EDOT .GT. EDTMX) RETURN
   95 fOutput = ABS(fOutput)*1.E6

  RETURN
END FUNCTION FCLIMB   ! -----------------------------------------------------

!+
FUNCTION FCLMB2(EPR) RESULT(fOutput)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: MinF

  REAL,INTENT(IN OUT):: epr

  REAL:: fOutput

!!!      EXTERNAL FCLIMB, PILIMT, FTHRST, FDRAG

  REAL,EXTERNAL:: Fclimb, Fdrag, Fthrst, Pilimt

!      REAL MINF
      REAL:: MACH,MACHID
      COMMON/CLIMB/MACH, D, ICLIMB,T
      COMMON/COST/EPR1,ICOST, FC, TC, DTEMPK,W,FUEL

  INTEGER:: nearcz,idle
      COMMON/CRUZ/NEARCZ,IDLE
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, EPRSET,MFGR

  REAL:: lambda,e,edot,ff
  INTEGER:: igrnd,iths
  REAL:: mnvtas,mxvtas, vtas2
      COMMON/ENERGY/LAMBDA,E,EDOT,FF,IGRND,ITHS,MNVTAS,MXVTAS,VTAS

      COMMON/EPSILN/EPSIL1, EPSIL2, ISPLMT

  INTEGER:: iprint, idrag
      COMMON/III/IPRINT, IDRAG

      COMMON/OPT/OPTMAK,IOPT,EDTMX
      COMMON/WINDY/IWIND, PSIA, VWA
!!!      DATA G, EDTMAX/32.2, 5./

  REAL:: costpi
  REAL:: costv
  REAL,PARAMETER:: EDTMAX=5.0
  REAL,PARAMETER:: G=32.2
  INTEGER:: icount
  REAL:: vmin,vmax,vmid
!----------------------------------------------------------------------------


   68 COSTV  = 1.E6
      COSTPI = 1.E6
      IDLE = 0
      ICOUNT = 1

!     COMPUTE SPEED LIMITS
      IF ( ICLIMB .EQ. 1) GO TO 69
      IF( NEARCZ .EQ. 1) MNVTAS = VTAS   - 50.
   69 VMAX= SQRT(2.*G*(E-H))
      IF(ISPLMT .EQ. 0) GO TO 70
      IF( H .GT. 10000.) GO TO  70
      VLIMIT = 250.*1.68781*SQRT(.0023769/RHO)
      Vmax = min(vmax, vlimit)
   70 vmax1 = vmax
      Vmax2 = a*.89
      Vmax = min(vmax1, vmax2)
!     IF( nearcz .eq. 1) vmax = vmax2

  IF (e > 39999.0) THEN
    VMIN1 = SQRT(2.0*G*(E - 39999.0))
  ELSE
    Vmin1 = 0.
  END IF
  VMIN = MAX( VMIN1, MNVTAS)

! following statement replaced - RLC - 27Oct99
!      GO TO (121, 122), ICLIMB

  SELECT CASE(iclimb)
    CASE(1)                                    !     CLIMB MODE
      IF (IOPT .EQ. 1) GO TO 80
      VMID = VTAS
      GO TO 88
   80 idrag = 1
      Eprset = 2.4
      Drgmin = minf(vmin,  vmax,  fdrag,vmid,  iprint)
      Idrag = 0
   88 epr1 = 2.4
   81 eprset = epr1
   83 dummy = minf(vmin,  vmid,  fdrag, v1, iprint)
      Dummy = fdrag(vmax)
      IF( T  .LT.  D) GO TO 87
      V2 = vmax
      GO TO 85
   87 dummy = minf(vmid, vmax, fdrag,  v2,   iprint)
   85 costv = minf(  v1, v2,   fclimb, vtas, iprint)

! following statement replaced  -  RLC  -  27Oct99
!      GO TO(103, 84, 86), IOPT
!
    IF (iopt < 0) THEN
      fOutput=costv
      RETURN
    END IF

    IF (iopt > 0) THEN                        !     FOR  HIGHER ACCURACY RUN
      V1x = vtas - 10.
      V2x = vtas + 10.
      V11 = max(v1, v1x, vmin)
      V22 = min(v2, v2x)
      Costv = minf(v11, v22, fclimb, vtas, iprint)
    END IF

    IF (ABS(COSTV - COSTPI) .LE. EPSIL1) GO TO 101
    DUMMY = MINF(1.1, 2.4, PILIMT, EPRMIN, IPRINT)
    COSTPI = MINF(EPRMIN, 2.4, FTHRST, EPR1, IPRINT)

    IF (iopt==3) THEN                      !     FOR  HIGHER ACCURACY RUN
      Eprmx1 = epr1 + .01
      Eprmn1 = epr1 - .01
      Eprmx = amin1( eprmx1, 2.4 )
      Eprmn = amax1(eprmn1, eprmin)
      Costpi = minf(eprmn, eprmx, fthrst, epr1, iprint)
    END IF
!
    IF( ABS(COSTV - COSTPI) .LE. EPSIL2 .OR. ICOUNT .GT. 3) GO TO 101
    ICOUNT = ICOUNT + 1
    GO TO 81

 101 EPR = EPR1
    fOutput = MIN(COSTV, COSTPI)

  CASE(2)                                             !     DESCEND MODE
    IF (IOPT .EQ. 1) GO TO 95
      Ioptid = iopt
      Iopt = 1
      Costid = minf(vmin, vmax, fclimb, vtas, iprint)
      Tid = t
      Hid = h
      Machid = mach
      Vtasid = vtas
      Edotid = edot
      Ffid = ff
      Iopt = ioptid
   92 dummy = minf(1.1,  2.4, pilimt, eprmax, iprint)
      IF ( iwind .eq. 0) GO TO 94
      CALL WIND (h, psia, vwa)
   94 COSTPI = MINF(1.1,  EPRMAX, FTHRST, EPR1, IPRINT)
      IF ( ABS(COSTV - COSTPI) .GT. EPSIL2) GO TO 95

   96 IF( (COSTV .LE. COSTID) .OR. (COSTPI .LE. COSTID)) THEN
        epr=epr1
        fOutput=MIN(costv, costpi)
        RETURN
      END IF
      Idle = 1
      T = tid
      Mach = machid
      Edot = edotid
      Ff = ffid
      foutput = costid
      RETURN

   95 COSTV = MINF( VMIN,  VMAX, FCLIMB, VTAS, IPRINT)
      IF( IOPT .EQ. 1) GO TO 93
      IF( ABS(COSTV - COSTPI) .LE. EPSIL2 .OR. ICOUNT .GT. 3) GO TO 96
      Icount = icount + 1
      Idrag = 1
      D = fdrag(vtas)
      GO TO 92
   93 IDLE = 1
  103 fOutput = COSTV
  END SELECT
  RETURN
END Function Fclmb2   ! -----------------------------------------------------

!+
FUNCTION FCOST(MACH) RESULT(fOutput)
! ---------------------------------------------------------------------------
! PURPOSE -
IMPLICIT NONE

  REAL,INTENT(IN):: MACH

  REAL:: fOutput

  REAL:: rho,p,tempk,a,ratio,h,alpha,cl,eprmax
  INTEGER:: mfgr
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, EPRMAX,MFGR

  REAL:: eprf
  INTEGER:: icost
  REAL:: fc,tc,dtempk,w,fuel
      COMMON/COST/EPRF, ICOST,FC,TC, DTEMPK,W,FUEL

  INTEGER:: iwind
  REAL:: psia, vwa
      COMMON/WINDY/IWIND, PSIA, VWA

  REAL:: gamma
  INTEGER:: returnCode
  REAL:: vtas
!----------------------------------------------------------------------------
  VTAS = A*MACH

  gamma=0.0   ! added by RLC  30Oct99   undefined otherwise
  CALL TRIM1(EPRF,VTAS , H, W, 0., 0, 1, GAMMA, ALPHA, 1, 1, 2, 0., &
     & 1.,0., FUEL, returnCode)

  IF (returnCode==1) THEN
    WRITE(6,*) "0FUNCTION OUTSIDE RANGE"
    RETURN
  END IF

!     ICOST = 1 RETURN FUEL/VTAS; =2 RETURN FUEL
  IF (icost==1) THEN
    fOutput = ((FC*FUEL + TC)/(VTAS + VWA))/(3600./6080.)
  ELSE
    fOutput = FUEL
  END IF

  RETURN
END Function Fcost   ! ------------------------------------------------------

!+
FUNCTION FDRAG(VTAS) RESULT(fOutput)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: At62
USE Engine,ONLY: EngEPR
  REAL,INTENT(IN):: vtas

  REAL:: fOutput

  INTEGER:: iprint, idrag
  COMMON/III/IPRINT, IDRAG

  REAL:: lambda,e,edot,ff
  INTEGER:: igrnd,iths
  REAL:: mnvtas,mxvtas,vtas1
      COMMON/ENERGY/LAMBDA,E,EDOT,FF,IGRND,ITHS,MNVTAS,MXVTAS,VTAS1

  REAL:: h, eprset
  INTEGER:: mfgr
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, EPRSET,MFGR

  REAL:: epr
      COMMON/COST/EPR,    ICOST, FC, TC, DTEMPK, W,FF1

  REAL:: mach, d, t
      COMMON/CLIMB/MACH, D, ICLIMB,T

  REAL,DIMENSION(4):: ans
  REAL,PARAMETER:: G = 32.2
  INTEGER:: incruz
!----------------------------------------------------------------------------
  h = e - 0.5*vtas*vtas/g
  IF (h .lt. 0.0) h = 0.0
  CALL AT62(H,ANS)
  Tempk = ans(3) + dtempk
  rho = ans(2) / (3092.40 * tempk)
  F = 0.5*rho*vtas**2*1560.0
  cl = w/f
  call cdrag(mach, cl, 0, 0., cd)
  D = f*cd

  IF (IDRAG==1) THEN
    fOutput=D
    RETURN
  END IF

  Epr = eprset
  mach = vtas/ans(4)
  incruz=0
  call engepr(h, mach, epr, incruz, t, ff, mfgr)
  foutput = abs(t-d)
  IF (t < d) foutput = foutput*1.0e6
  IF (IPRINT==0) RETURN

  WRITE(6, 100) Epr,  t, d, foutput, vtas
  100 FORMAT(/'EPR, T, D, ABS(T-D)', F10.2, 4F10.0)

  RETURN
END FUNCTION FDRAG   ! ------------------------------------------------------

!+
FUNCTION FOPT(ALT) RESULT(fOutput)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: At62,MinF
USE Engine,ONLY: EngEPR

      IMPLICIT REAL(M)
  REAL,INTENT(IN):: alt

  REAL:: fOutput


!!!      EXTERNAL Fbound, FCOST

INTERFACE

!FUNCTION Fbound(MACH) RESULT(f)
!  REAL,INTENT(IN):: MACH
!  REAL:: f
!END Function Fbound

!FUNCTION FCOST(MACH) RESULT(f)
!  REAL,INTENT(IN):: MACH
!  REAL:: f
!END Function Fcost

END INTERFACE

  REAL,EXTERNAL:: Fbound
  REAL,EXTERNAL:: Fcost
  REAL:: mach
  REAL:: mindrg

  REAL:: rho,p,tempk,a,ratio,h,alpha,cl,eprmax
  INTEGER:: mfgr
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, EPRMAX,MFGR

  REAL:: machop
  INTEGER:: iopt
  REAL:: edtmx
      COMMON/OPT/MACHOP,IOPT,EDTMX

  INTEGER:: iprint, idrag
      COMMON/III/IPRINT, IDRAG

  REAL:: epr
  INTEGER:: icost
  REAL:: fc,tc,dtempk,w,fuel
      COMMON/COST/EPR, ICOST, FC, TC, DTEMPK,W,FUEL

  INTEGER:: iwind
  REAL:: psia, vwa
      COMMON/WINDY/IWIND, PSIA, VWA

  REAL,DIMENSION(4):: ans
  REAL:: fa,fb
!----------------------------------------------------------------------------
      H =ALT
   89 CALL AT62(ALT, ANS)
      Tempk = ans(3) + dtempk
      P = ans(2)
      Rho = p / (3092.40 * tempk)
      A = 65.76 * sqrt(tempk)
   90 eprmax = 2.4
      Idrag = 1
      Mindrg = minf(0., .9,fbound, mach, iprint)  ! look for min of fbound
      CALL ENGEPR(ALT, MACH, EPRMAX, 1, TMAX, FF, MFGR)

  IF (TMAX < MINDRG) THEN
    WRITE(6,'(A,2F10.0)') " DRAG EXCEEDS MAX THRUST", MINDRG, TMAX
    fOutput=0.0
    RETURN
  END IF

      IF (icost .eq. 1) GO TO 94
      Fa = mach - .1
      F =.5*(a*fa )**2*1560.*rho
      Cl = w/f
        Call cdrag (fa, cl, 0, 0., cd)
      D = f*cd
      Call engepr(alt, fa, eprmax, 1, t, ff, mfgr)
      IF( t .ge. D) GO TO 91
      Idrag = 2
      Fac = minf(.1, mach, fbound, fa, iprint)
      GO TO 91
   94 fa = mach
   91 fb = .9
      Call engepr(alt, fb, eprmax, 1, t, ff, mfgr)
      F =.5*(a*fb )**2*1560.*rho
      Cl = w/f
        Call cdrag (fb, cl, 0, 0., cd)
      D = f*cd
      IF( t .ge. D) GO TO 95
      Idrag = 2
      Fbc = minf(mach, .9, fbound, fb, iprint)
   95 IF ( iwind .ne. 0) call wind(h, psia, vwa)

  fOutput=MINF(FA, FB, Fcost, MACHOP, IPRINT)

  RETURN
END Function Fopt   ! -------------------------------------------------------

!+
FUNCTION FTHRST(EPR)   ! not PURE
! ---------------------------------------------------------------------------
! PURPOSE -

USE Engine,ONLY: EngEPR

  REAL,INTENT(IN OUT):: epr

  REAL:: lambda,e,edot,ff
  INTEGER:: igrnd,iths
  REAL:: mnvtas,mxvtas, vtas
      COMMON/ENERGY/LAMBDA,E,EDOT,FF,IGRND,ITHS,MNVTAS,MXVTAS,VTAS

  REAL:: rho,p,tempk,a,ratio,h,alpha,cl,eprset
  INTEGER:: mfgr
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, EPRSET,MFGR

      COMMON/COST/EPR1,ICOST, FC, TC, DTEMPK,W,FUEL

  REAL:: mach,d
  INTEGER:: iclimb
  REAL:: t
      COMMON/CLIMB/MACH, D,ICLIMB,T
      COMMON/OPT/OPTMAK,IOPT,EDTMX
      COMMON/WINDY/IWIND, PSIA, VWA
      COMMON/CRUZ/NEARCZ,IDLE

   INTEGER:: incruz
!----------------------------------------------------------------------------

  incruz=0
   90 call engepr(h, mach, epr, incruz, t, ff, mfgr)
      Edot = (t-d)*vtas/w
   93 fthrst= (fc*ff + tc - lambda*(vtas+vwa))/ abs(edot)
      IF( abs(edot) .le. 5.) GO TO 95
      IF( nearcz .eq. 1 .and.  Fthrst .gt. 0. .and. Iclimb .eq. 2) GO TO&
     &    95
      IF(iclimb .ne. 1 .or. Iopt .ne. 3) return
      IF ( edot .gt. Edtmx) return
   95 fthrst = abs(fthrst)*1.e6
      RETURN
END Function Fthrst   ! -----------------------------------------------------

!+
FUNCTION FULEST(RANGE, IOPARM, ETO, FC, TC, PC, ECRUZ, EF) ! not PURE
! ---------------------------------------------------------------------------
! PURPOSE - Compute climb fuel
  REAL,INTENT(IN):: range    ! range not used ???
  INTEGER,INTENT(IN):: ioParm
  REAL,INTENT(IN):: eto,fc,tc,pc
  REAL,INTENT(IN OUT):: ecruz
  REAL,INTENT(IN OUT):: ef

      REAL K1, K2
      COMMON/VTRCRU/WSS(10), JJCRUZ(10), LLCRUZ(50,10),EECRUZ(50,10),   &
     & DLLDEE(2, 10),IWMAXX, WTO, WCRUZ ,CRUZCT,HHCRUZ(50,10),          &
     & FFCRUZ(50,10),JLAST1, JLAST2
      DATA EVOPT1, VOPT1, EVOPT2, VOPT2/.1130, .1079, 6.1463E-6,        &
     & 4.6849E-6/
!----------------------------------------------------------------------------
  write(3,*) "entering Fulest"
   70 wcruz = wto
      Fulst1 = 0.
      IF (ioparm .ne. 0) GO TO 71
      K1 = vopt1
      K2 = 1. + (tc/fc)*vopt2
      GO TO 72
   71 k1 = evopt1
      K2 = 1. + (tc/fc)*evopt2

  write(3,*) "starting the 72 loop in fulest, k1,k2=", k1,k2
   72 call wlefhv(1, ecruz, ef, pc, vgknt, vcktas)
      Fulest = k1*(ecruz - eto)*k2*wto/136000.
      Wcruz = wto - fulest
      IF( abs(fulest - fulst1) .le. 100.) return
      Fulst1 = fulest
      GO TO 72
!
!
!     COMPUTE TOTAL FUEL USED
END Function Fulest   ! -----------------------------------------------------


!+
SUBROUTINE PCCOMP( PC, IPC, X3, CRUDST, TDIST,LINEAR,returnCode,  &
     & IOPARM, ISPLIZ)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: Serch1

  REAL,INTENT(IN OUT):: pc
  INTEGER,INTENT(IN OUT):: ipc
  REAL,INTENT(IN):: x3
  REAL,INTENT(IN):: crudst
  REAL,INTENT(IN):: tdist
  INTEGER,INTENT(IN):: linear
  INTEGER,INTENT(OUT):: returnCode
  INTEGER,INTENT(IN OUT):: ioparm
  INTEGER,INTENT(IN OUT):: ispliz

!!!      DIMENSION RANGE(10), COSTPC(10)
!!!      EQUIVALENCE(RMIN, RANGE(1)), (RMAX, RANGE(2))
!      DATA ONE, TWICE/1., 2./
  REAL,PARAMETER:: ONE=1.0, TWICE=2.0

  REAL,SAVE,DIMENSION(10):: range,costpc
!----------------------------------------------------------------------------
  write(3,*) "Entering PcComp, ipc=", ipc
      returnCode=0
!
      Iter = ipc - 2
      Write(6, 30) pc, iter
   30 format(/'cost (% over lambda) = ', 2x, f10.2, &
       'No of iterations = ', i5)
      IF ( ipc .ge. 10) stop
      IF( ispliz .eq. 0) stop
      I= ipc
      IF( I .gt. 3) i= 3

! following statement replaced  -  rlc  -  27oct99
!      GO TO (101, 102, 103), I
  select case(i)
    Case(1)
      Costpc(1) = pc
      Range(1) = tdist
      IF (crudst <= range(1)) stop   ! wow, that's a powerful error
      Ipc = 2
      Pc = 1.
      Costpc(2) = 1.
      Init = 0
      Ialter = 1
    Case(2)
      Range(2) = tdist
      IF (abs(crudst - range(2)) <=  5.0) stop
      IF (range(2)  <= (crudst+ 5.)) then
        Ispliz = 0
        IF( ioparm .eq.1) ioparm = 2
        returncode=1
        Return
      END IF
    Case(3)
      IF ( abs(tdist - crudst) .lt. 5.) stop
      Range(insrt1) = tdist
    Case default     ! don't think you should need this, but.....
      Costpc(1) = pc
      Range(1) = tdist
      IF ( crudst .le. range(1)) stop   ! wow, that's a powerful error
      Ipc = 2
      Pc = 1.
      Costpc(2) = 1.
      Init = 0
      Ialter = 1

  End select

!     Compute percentage change in lambda(*)

!     Spacefor next tdist such that r(1)...r(insert),tdist,r(insrt+1)
!     Will store to be computed tdist in r(insrt1)
  104 call serch1(range, crudst, ipc, pf, insert, limit)
      Insrt1  = insert + 1
  353 isum = ipc + insrt1
      write(3,*) "at 353 in pccomp, ipc,insert,isum=", ipc,insert,isum
      Do i= insrt1, ipc
        J = isum - I
        Range(j+ 1) = range(j)
        Costpc(j+1) = costpc(j)
      End do

!     Find best two points r(i1), r(i2)

      Dist2 = crudst - range(insert)
      Dist3 = range(insrt1) - crudst
      IF (dist3 .gt. 100. .and. Dist2 .gt. 100. .and. Init .eq. 0) GO TO 358
      Der= 20.
      IF( costpc(insert)  .le. 2.) der= 10.
      I1= insert
      I2 = insrt1
      IF( insert .lt. 2) GO TO 345
      Dist1 = crudst - range(insert - 1)
      IF (dist3 .le. (Dist1+ der)) GO TO 345
      I2 = insert - 1
      GO TO 354
  345 IF( ipc .lt. (Insrt1+1)) GO TO 354
      Dist4 = range(insrt1 + 2) - crudst
      IF( dist2 .le. (Dist4 + der)) GO TO 354
      I1 = insrt1 + 2

  354 pc1 = costpc(i1)
      Pc2 = costpc(i2)
        IF( abs(pc1 - pc2) .le. .2  .or. Ialter .eq. 1) GO TO  373
  355 x1 = one/range(i1)
      X2 = one/range(i2)
      IF( costpc(insert)  .gt. 2.) GO TO 370

!     Linear reciprocal fit
      Det = x1 - x2
      Y11 = one    /det
      Y12 = - one/det
      Y21 = -x2/det
      Y22 = x1/det
      GO TO 371

!     Quadratic reciprocal fit
  370 det = x1*x2*(x1 - x2)
      Y11 = x2/det
      Y12 = - x1/det
      Y21 = - x2*x2/det
      Y22 =  x1*x1/det
  371 a = y11*pc1 +  y12*pc2
      B = y21*pc1 + y22 *pc2
      IF( costpc(insert)  .le. 2.) GO TO 372
      Costpc (insrt1) = a*x3*x3 + b*x3
      GO TO 360
  372 costpc (insrt1) = a*x3 + b
      Ialter = 1
      GO TO 360

!     Linear interpolation or extrapolation

  373 costpc( insrt1) = costpc(i1) + (crudst - range(i1))*(costpc(i2) - &
     & costpc(i1))/(range(i2) - range(i1))
      IF( costpc(insrt1) .le. 0.) GO TO 355
      Ialter = 2
      GO TO 360

!     Linear interpolation

  358 costpc(insrt1) = costpc(insert) + pf*(costpc(insrt1) -            &
     & costpc(insert))
      Init = 1
  360 pc = costpc(insrt1)
      Ipc = ipc + 1

  RETURN
END Subroutine PcComp   ! ---------------------------------------------------

!+
FUNCTION PILIMT(EPR) RESULT(fOutput)   ! not a PURE function
! ---------------------------------------------------------------------------
! PURPOSE -
USE Engine,ONLY: EngEPR

  REAL,INTENT(IN OUT):: epr

  REAL:: fOutput

  REAL:: MACH,d
  INTEGER:: iclimb
  REAL:: t
      COMMON/CLIMB/MACH, D,ICLIMB, T

  REAL:: rho,p,tempk,a,ratio,h,alpha,cl,eprset
  INTEGER:: mfgr
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, EPRSET,MFGR

  INTEGER:: incruz
!----------------------------------------------------------------------------
  incruz=0
  CALL ENGEPR(H, MACH, EPR, incruz, T, FF, MFGR)
  fOutput=ABS(T-D)
  SELECT CASE(iclimb)    !!! GO TO (90, 91),  ICLIMB
    CASE(1)
      IF (t < d) fOutput = fOutput*1.E6
    CASE(2)
      IF (t > D) fOutput = fOutput*1.E6
  END SELECT

  RETURN
END Function PiLimt   ! -----------------------------------------------------


!+
SUBROUTINE TRIM1(EPR, VTAS, H, W, DF, GEAR,INCRUZ,GAMMA,ALPHA,    &
  ICNTRL, INIT, MODE, VDOT, COSPHI, PSI, FF, returnCode)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: At62
USE Engine,ONLY:EngEPR

  REAL,INTENT(OUT):: epr
  REAL,INTENT(IN):: vtas
  REAL,INTENT(IN):: h
  REAL,INTENT(IN):: w
  REAL,INTENT(IN):: df
  INTEGER,INTENT(IN):: gear
  INTEGER,INTENT(IN):: incruz
  REAL,INTENT(IN OUT):: gamma
  REAL,INTENT(OUT):: alpha
  INTEGER,INTENT(IN):: icntrl
  INTEGER,INTENT(IN):: init
  INTEGER,INTENT(IN):: mode   ! =1   =2   =3
  REAL,INTENT(IN):: vdot
  REAL,INTENT(IN):: cosphi
  REAL,INTENT(IN):: psi
  REAL,INTENT(OUT):: ff
  INTEGER,INTENT(OUT):: returnCode

!      INTEGER GEAR
      REAL MACH
!     EXTERNAL TRIMT
      COMMON/B1/WINDE(23,5),DWXDH, DWYDH, AKW, VW, PSIW
      COMMON/COST/EPRXX,ICOST,FC,TC,DTEMPK,WW,FUELDT

  REAL:: rho,p,tempk,asos,ratio,alt,alpha1,cl1,thsmax
  INTEGER:: mfgr
      COMMON/DRAGMN/RHO, P, TEMPK, ASOS,RATIO,ALT,ALPHA1,CL1,THSMAX,MFGR
      COMMON/ENGIN/RPMAX

  INTEGER:: iflap
  REAL:: fflaps
  COMMON/FLAP_SETTINGS/ iflap,FFLAPS

  REAL:: fn,fa
  INTEGER:: l
  REAL:: d,t
  COMMON/RITE/FN,FA,L,D,T,CL,CD,Q,ALPH1,TSINA,VDO1,GAMDOT,MACH,DUMMY

      COMMON/TTRIM1/H1, MACH1,T2, FF1
      COMMON/TRIM2/CL0, CLDA,CLGEAR
      DIMENSION ANS(4)
      DIMENSION CL216(3,7)
      DATA CL216/                       .012578 , .101778,  -.001426,   &
     &    .06282,   .105492,  -.000959, .133972,  .109409,  -.000737,   &
     &    .297618,  .117523,  -.000866, .604512,  .110909,  -.000483,   &
     &    .891152,  .108086,  -.000499, 1.207021, .106244,  -.000104/
      DATA ALFAMN,ALFAMX,RD2DEG/-5., 25., 57.296/
      DATA RHOZ/.0023769/
      DATA DTDEPZ,G/41000.,32.2/
!
!
!
!     ICNTRL = 1 CONSTANT VTAS, 2 CONSTANT MACH, 3 CONSTANT VIAS
!
!----------------------------------------------------------------------------
  returnCode=0

!      GO TO (97, 102,103), ICNTRL
!  102 AZERO = ASOS
!      GO TO 97
!  103 RHOLAS = RHO

  IF (icntrl==0) azero=asos   ! from /DRAGMN/
  IF (icntrl>0) rholas=rho    ! from /DRAGMN/
  CALL AT62(H, ANS)   ! was statement 97

  Tempk = ans(3) + dtempk
  rho = ans(2) / (3092.40 * tempk)
  Asos = 65.76 * sqrt(tempk)
  Dtdepr = dtdepz
  q = rho    *vtas*vtas/2.
  Qs = q*1560.
  Mach = vtas/asos
  alphap=0.0   ! added by rlc 30oct99. Otherwize undefined

  Call clift(mach, h, alphap, df, gear, cl,1)
  I2 = iflap + 1


  Call engepr(h, mach, epr, incruz, t, ff, mfgr )

   96 GAMMR = GAMMA/RD2DEG          ! gamma is dummy arg.

      GO TO (110, 111, 112), MODE          ! next 5 mod by RLC 2Nov99
  110 GAMMRZ = GAMMR
      GO TO 112

  111 FAZ = W*VDOT/G

!  SELECT CASE(mode)
!    CASE(1)
!      gammrz=gammr
!    CASE(2)
!      faz=w*vdot/G
!  END SELECT

! keep coming back here until ...

  112 B11= CL216(1, iflap) + FFLAPS*( CL216(1, I2) - CL216(1,iflap ))
      B2 = CL216(2, iflap) + FFLAPS*( CL216(2, I2) - CL216(2,iflap ))
      B3 = CL216(3, iflap) + FFLAPS*( CL216(3, I2) - CL216(3,iflap ))
      A = qs*b3*cosphi
      B =(qs*(b2 + clda) + t/rd2deg)*cosphi
      Psir = psi/rd2deg
!     Call winmod(h, psir)
      Akw = 0.
      C =qs*(b11+cl0 + clgear)*cosphi - w*(cos(gammr) - akw*vtas*       &
     & sin(gammr)**2/g)
      Disc = b*b - 4.*a*c
      IF (disc .lt. 0.) GO TO 990     ! error
      Rad = sqrt(disc)
      Alpha1 = (-b + rad)/(2.*a)
      IF (alpha1 .ge. Alfamn .and. Alpha1 .le. Alfamx) GO TO 210
      Alpha = (-b - rad)/(2.*a)
      IF( alpha .lt. Alfamn .or. Alpha .gt. Alfamx) GO TO 991   ! error
      GO TO 211
  210 alpha = alpha1
  211 alphar = alpha/rd2deg

      IF ( MODE .EQ. 3) GO TO 992
      Cl = b11+ b2*alpha + b3 * alpha*alpha + clda*alpha+cl0+clgear
      Call cdrag(mach, cl, gear, df,cd)
      D = qs*cd
      IF( mode .eq. 4) return
      Wsgam = w*sin(gammr)
      Cosa = cos(alphar)
      Ftd = t*cosa - d
      Ftdw = ftd/w
      GO TO (220, 221, 992), mode

  220 IF (icntrl .eq. 1 .or. Init .eq. 1   .or. Abs(h - hzero) .lt. 5.) &
     & GO TO 205
      Xx =   vtas /g
      GO TO (205, 201,202),icntrl
  201 f =    Xx*mach*( asos  - azero)/(h - hzero)
      GO TO 206

  202 f=-vtas*xx*(  rho  - rholas)/(2.*(h - hzero)* rho)
      GO TO 206

  205 f = 0.

  206 a1 = 1. + f
      A2 = akw*vtas/g
      Gammr = ftdw/(a1 + a2)
  207 twogam = 2.*gammr
      Fgamma = a1 *sin(gammr) + .5*a2*sin(twogam) - ftdw
      Den = a1*cos(gammr) + a2*cos(twogam)
      Def = .0018*abs(den)
      IF( abs(fgamma) .le. Def  ) GO TO 208
      Degam = fgamma/den
      Gammr = gammr - degam
      GO TO 207
  208 IF( abs(gammr - gammrz) .le. .0018) GO TO 992
      GO TO 110

  221 wsgam = w*sin(gammr)*(1.+ akw*vtas*cos(gammr)/g)
      Fa = ftd - wsgam
      IF( abs(fa - faz) .le. 4. .or. (Abs(fa-faz) .le. 25. .and. Mfgr   &
     &.eq. 2 .and. H .gt. 35000.))GO TO 992
      T2 = (d + faz + wsgam)/cosa      ! required thrust
      Tz = t
      Eprz=epr

  311 epr1 = eprz - (tz - t2)/dtdepr
      Call engepr(h, mach, epr1, incruz, t1, ff, mfgr)
      IF( abs(t2 - t1) .le. 4.) GO TO 312
      IF (epr1 .lt. 2.4) GO TO 313
      Call engepr(h, mach, epr1, incruz, t1, ff, mfgr)

      IF ( t1 <= t2) then
        Write(6, '(a,2f15.0)' ) ' required thrust, max thrust=', t2, t1
        returncode=1
        Return
      END IF

  313 IF (eprz .eq. Epr1) GO TO 993
      Dtdepr = (tz - t1)/(eprz - epr1)
      Eprz = epr1
      Tz = t1
      GO TO 311   ! loop until ...

  312 t = t1
      Epr = epr1
      GO TO 112     ! keep looping until something happens

  990 write(6, 50) disc
   50 format(/'discriminat = ', e15.6)
      returncode=1
      Return

  991 write(6, 51) alpha1, alpha, alfamn, alfamx
   51 format(/'alpha(1), alpha(2), alpha(min), alpha(max) = ', 4e15.6)
      returncode=1
      Return

  993 IF ( abs (fa - faz) .le. 20.)GO TO 992
      Ff = 1.e6

  992 gamma = gammr*rd2deg
      Hzero = h
      RETURN

END Subroutine Trim1   ! ----------------------------------------------------


!+
SUBROUTINE UPDOWN(EOPT, IPRINT, WTO, WLNDG, ECRUZ,WCRUZ)
! ---------------------------------------------------------------------------
! PURPOSE -
!   called by main program

USE Extra,ONLY: At62,Iclock,Page,Serch1

      IMPLICIT REAL(M)

  REAL,INTENT(IN):: eopt
  INTEGER,INTENT(IN):: iprint
  REAL,INTENT(IN):: wto
  REAL,INTENT(IN):: wlndg
  REAL,INTENT(IN):: ecruz
  REAL,INTENT(IN):: wcruz

      REAL:: LAMBS

!!!      EXTERNAL FCLIMB, FCLMB2
INTERFACE

FUNCTION FCLIMB(VTAS) RESULT(f)
  REAL,INTENT(IN):: vtas
  REAL:: f
END Function Fclimb

FUNCTION FCLMB2(EPR) RESULT(f)
  REAL,INTENT(IN):: epr
  REAL:: f
END Function Fclmb2

END INTERFACE



      COMMON/CLIMB/MACH, D, ICLIMB,TDUMMY
      COMMON/COST/EPR, ICOST, FC, TC, DTEMPK, W, FUELDT
      COMMON/CRUZ/NEARCZ,IDLE
      COMMON/DRAGMN/RHO, P, TEMPK, A, RATIO, H, ALPHA, CL, THSMAX, mfgr

  REAL:: lambda,e,edot,ff
  INTEGER:: igrnd,iths
  REAL:: mnvtas,mxvtas, vtas2
      COMMON/ENERGY/LAMBDA,E,EDOT,FF,IGRND,ITHS,MNVTAS, MXVTAS,VTAS2

      COMMON/IO/ ANS(4), WS(11), EOPTS(11), HSTARS(11), MSTARS(11),    &
     & PISTRS(11), LAMBS(11), VTASOP(11), FUELFL(11), CRDIST(11),       &
     &CRTIME(11), MNCOST(2), MOPIAS(2), MOPTAS(2), MACHOP(2), FDTOPT(2),&
     & OPTALT(2), HOPT(2), EPRS(2),IWMAX,IOPARM,HTO, VTO, HOLNDG, VOLNDG&
     & ,ETO

  REAL:: optmak
  INTEGER:: iopt
  REAL:: edtmx
      COMMON/OPT/OPTMAK,IOPT,EDTMX


      COMMON/VUPDWN/EECLMB(110), FFCLMB(110), DDCLMB(110), EEDOWN(110), &
     & FFDOWN(110), DDDOWN(110), SSCOST(110) ,JCLIMB,JDESCN,TTCLMB(110),&
     & TTDOWN(110)

  INTEGER:: iwind
  REAL:: psia,vwa
      COMMON/WINDY/IWIND, PSIA, VWA

!      DATA FS2KNT, G, RD2DEG, RHOSL/1.68781,32.2,57.296,.0023769/
  CHARACTER(LEN=*),PARAMETER:: FMT189= &
    "(/'AIRCRAFT TAKE OFF WT = ', F8.0,2X, '#S, INITIAL WT =',"//  &
    "F10.0, '#S, CRUISE ENERGY = ', F7.0,2X,'FT')"
  CHARACTER(LEN=*),PARAMETER:: FMT190= &
    "(/'INITIAL ALT (FT), SPEED (KIAS)', 2F10.0 )"

  CHARACTER(LEN=*),PARAMETER:: FMT191A= "CLIMB OPTIMIZATION:", &
    FMT191B=" ENERGY  ALTITUDE    MACH   VIAS  VTAS      EDOT    GAMMA"// &
     "     TIME    DIST     FUEL USED PWR SETG   COST/E", &
    FMT191C="    FT        FT       NO    KNOT  KNOT     FT/SEC"// &
     "     DEG     HR:MN:SEC  N MILE        #       EPR      $/ E FT"

  CHARACTER(LEN=*),PARAMETER:: FMT210= &
    "(2F10.0, F9.3, 2F6.0,2F10.2, 4X, 2(I2, ':'), I2, F10.3, F10.0, 3F10.3)"
  CHARACTER(LEN=*),PARAMETER:: FMT228= &
    "(2F10.0,  F9.3,    2F6.0,2F10.2, 4X, 2(I2, ':'),  I2,"//   &
    "F10.3, F10.0, 5X, 'IDLE ', 2F10.3)"

  CHARACTER(LEN=*),PARAMETER:: FMT289= &
    "(/'AIRCRAFT LANDING WT = ', F10.0, 2X, '#S, CRUISE ENERGY = ',"// &
    "F10.0, 2X, 'FT')"
  CHARACTER(LEN=*),PARAMETER:: FMT290= &
    "(/'FINAL ALT(FT), SPEED (KIAS) = ', 2F10.0)"
  CHARACTER(LEN=*),PARAMETER:: FMT291= &
    "/'DESCEND OPTIMIZATION:'/"//                               &
    "4X, 'ENERGY  ALTITUDE    MACH   VIAS  VTAS', 6X, 'EDOT', 5X,"//  &
    "'GAMMA', 6X,'TIME', 6X,'DIST', 5X,"// &
    "'FUEL USED PWR SETG   COST/E  SUM COST/E'/"// &
    "6X, 'FT', 9X, 'FT', 7X, 'NO', 4X, 'KNOT  KNOT', 5X, 'FT/SEC', 5X,"//  &
    "'DEG', 5X, 'HR:MN:SEC  N MILE', 8X, '#', 7X, ' EPR', 5X,"// &
    "'$/ E FT',  4X, '$/E FT')"


  REAL,PARAMETER:: FS2KNT=1.68781, G=32.2, RD2DEG=57.296, RHOSL=0.0023769
  REAL,DIMENSION(110):: clmcst  !    DIMENSION CLMCST(110)
  REAL:: gamma
  INTEGER:: icount
  INTEGER,SAVE:: iupmax
!----------------------------------------------------------------------------
   25 FORMAT(/'FUEL COST($/#)', F9.4, 2X,'TIME COST($/HR)= ', F7.2, 2X, &
       'TEMP VAR (DEG K) =', F7.2, 'LAMBDA =', F10.3, '$/NM')

  189 FORMAT(/'AIRCRAFT TAKE OFF WT = ', F8.0,2X, '#S, INITIAL WT =',  &
       F10.0, '#S, CRUISE ENERGY = ', F7.0,2X,'FT')
  190 FORMAT(/'INITIAL ALT (FT), SPEED (KIAS)', 2F10.0 )
  191 FORMAT(/'CLIMB OPTIMIZATION:'/                                 &
       4X, 'ENERGY  ALTITUDE    MACH   VIAS  VTAS', 6X, 'EDOT' 5X, &
      'GAMMA', 6X, 'TIME' 6X, 'DIST' 5X, 'FUEL USED PWR SETG   COST/E'/ &
       6X, 'FT', 9X, 'FT', 7X, 'NO', 4X, 'KNOT  KNOT', 5X, 'FT/SEC', 5X,  &
     &  'DEG', 5X, 'HR:MN:SEC  N MILE', 8X, '#', 7X, ' EPR', 5X, '$/ E FT')
  289 FORMAT(/'AIRCRAFT LANDING WT = ', F10.0, 2X, &
       '#S, CRUISE ENERGY = ', F10.0, 2X, 'FT')
  290 FORMAT(/'FINAL ALT(FT), SPEED (KIAS) = ', 2F10.0)
  291 FORMAT(/'DESCEND OPTIMIZATION:'/                               &
       4X, 'ENERGY  ALTITUDE    MACH   VIAS  VTAS', 6X, 'EDOT', 5X,  &
       'GAMMA', 6X,'TIME', 6X,'DIST', 5X, &
       'FUEL USED PWR SETG   COST/E  SUM COST/E'/ &
       6X, 'FT', 9X, 'FT', 7X, 'NO', 4X, 'KNOT  KNOT', 5X, 'FT/SEC', 5X,  &
       'DEG', 5X, 'HR:MN:SEC  N MILE', 8X, '#', 7X, ' EPR', 5X, &
       '$/ E FT',  4X, '$/E FT')



  clmcst(:)=0.0   ! added by RLC, but shouldn't need it


   80 nearcz = 0
!      GO TO (81, 300), iclimb

  IF (iclimb==1) then
    Ho = hto
    Vo = vto
    W = wto
    Epr = 2.4
    Jclimb= 0
  else
    Call page
    Ho = holndg
    Vo = volndg
    W = wlndg
    Jdescn = 0
  END IF

  H = ho
  call at62(ho, ans)
  Rho=ans(1)*(1.+ dtempk/ans(3))
  P = ans(2)
  Tempk = ans(3) + dtempk
  a = ans(4)
  Vo = vo*sqrt(rhosl/rho)*fs2knt
  e = ho + .5*  (Vo)**2/g


  !    GO TO (182, 183), iclimb
  select case(iclimb)
    Case(1)
      Write(6,189) w,wcruz,ecruz
      Write(6, 190) ho, vto
      Write(6, 25) fc, tc, dtempk, lambda
      Lambda = lambda*(3600./6080.)
      Write(6, 191)

    Case(2)
      Write(6, 289) w, ecruz
      IF( iprint .eq. 0) GO TO 192
      Write(6, 290) ho,volndg
      Write(6, 291)
  End select

192   denrgy = 500.
      Iopt = 1
      Mnvtas = 330.
      Mxvtas = 840.
      Eshift = eopt - (3000. - 1.)
      Ediff = eshift - e
      Ie500 = ediff/500.
      Deinit = ediff - ie500*500.
      Instep = 0

      IF (iclimb/=1) then
      Icsame = 1
      IF ( e .eq. Eto) icsame = 0
      END IF

  194 time = 0.
      Dist = 0.
      Gamma = 0.
      Init = 0
      Fueluz = 0.
      Ienry = 0
      Iemax = 0
      IF (iwind/=0) call wind(ho, psia, vwa)

  206 vgo = vo + vwa
      Icount = 1


! start of big loop.....
  201 ccost = fclmb2(epr)*6080./3600.
      IF (abs(edot) .ge. 5.0) GO TO 207
      IF (denrgy .lt. 10.0)GO TO 208
      Denrgy = denrgy / 2.
      E = e - denrgy
      GO TO 201

  208 icount = icount + 1
      GO TO 900

  207 vtas = vtas2
      IF ( ioparm .eq. 0) GO TO 202
      IF( iopt .eq. 3) GO TO 205
      IF ( iopt .eq. 2) GO TO 204
      IF   (E .ge. 10000.) iopt = 2
      GO TO 205

  204 IF( (eopt - e) .gt. 3000.) GO TO 205
      Denrgy = denrgy/2.
      Iopt = 3
      GO TO 205

  202 IF ( ienry .eq. 1) GO TO 205
      IF ( (eopt - e) .gt. 3000.) GO TO 205
      Denrgy = denrgy/2.
      Ienry = 1

  205 det = denrgy/edot
      W = w - (det*ff)/3600.
      Time = time +abs(det)
      Call iclock(TIME, IHR, IMIN, ISEC)
      IF( init .eq. 0) GO TO 203
      Gsprnt = vtas + vwa
      Singa =  2*(h - ho)/(det*(vgo + gsprnt))
      Gamma = asin(singa)    ! was arsin   (Rlc)
      Dex = .5*(vgo + gsprnt)*cos(gamma)*det
      Dist = dist + abs(dex/6080.)
      Gamma = gamma*rd2deg

  203 vgo = vtas + vwa
      Ho = h
      Fueluz = fueluz + ff*abs(det)/3600.
      Vtask = vtas/fs2knt
      Viask = 29.*sqrt(p*((1.+.2*mach*mach)**3.5-1.))/fs2knt
      IF ( (ecruz -(e + denrgy)) .lt. 5000.) nearcz = 1


      GO TO (220, 221), iclimb

  220 WRITE(6, 210) E, H, MACH, VIASK, VTASK, EDOT, GAMMA, IHR, IMIN,   &
     & ISEC, DIST, FUELUZ ,EPR, CCOST

      Clmcst(icount) = ccost   ! icount could get out of bounds
      Jclimb= jclimb + 1
      Eeclmb (jclimb) = e
      Ffclmb(jclimb) = fueluz
      Ddclmb(jclimb) = dist
      Ttclmb(jclimb) = time
  210  format(2f10.0,  f9.3,    2f6.0,2f10.2, 4x, 2(i2, ':'),  i2,   &
     & f10.3, f10.0, 3f10.3)
      GO TO 222


  221 IF ( icsame .ne. 0) GO TO 223
      Upcost = clmcst(icount)
      GO TO 225

  223 IF ( e .lt. Eeclmb(1)) GO TO 224
      IF( e .gt. Eeclmb(iupmax)) GO TO 224
      Call serch1(eeclmb, e, iupmax, pf, iup, limit)
      Upcost = clmcst(iup) + pf*(clmcst(iup+1) - clmcst(iup))
      GO TO 225

  224 upcost = 0.

  225 sumcst = ccost + upcost
      IF( iprint .eq. 0) GO TO 91
      IF (idle .eq. 1) GO TO 227
      Write(6, 210) e, h, mach, viask, vtask, edot, gamma, ihr, imin,   &
     & isec, dist, fueluz ,epr, ccost                    ,sumcst
      GO TO 91

  227 WRITE(6, 228) E, H, MACH, VIASK, VTASK, EDOT, GAMMA, IHR, IMIN,   &
     & ISEC, DIST, FUELUZ, CCOST,SUMCST
  228  FORMAT(2F10.0,  F9.3,    2F6.0,2F10.2, 4X, 2(I2, ':'),  I2,   &
     & F10.3, F10.0, 5X, 'IDLE ', 2F10.3)


   91 jdescn= jdescn +1
      Eedown(jdescn) = e
      Ffdown(jdescn) = fueluz
      Dddown(jdescn) = dist
      Ttdown(jdescn) = time
      Sscost(jdescn) = sumcst


  222 icount = icount + 1
      Init = 1
      IF( instep .eq. 0) GO TO 226
      E = e + denrgy
      IF (e .lt. Eopt - 100.) GO TO 201
      IF (iemax .eq. 1) GO TO 900   ! jump out of loop
      E = eopt - 50.
      Iemax = 1
      GO TO 201

  226 e = e + deinit
      Instep = 1
      GO TO 201

  900 IF( iclimb .ne. 1) return
      Iupmax = icount - 1
      RETURN
END Subroutine UpDown   ! ---------------------------------------------------






SUBROUTINE VOPTRJ(WLANDG,TDIST,PC,LAMBDX, ISPLIZ, RANGE,IOPARM,EF,ITER)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: At62,Page,Serch1

      REAL LLCRUZ, LAMBDI, LAMBDF,LAMBDX,MACHI, MACHF
      COMMON/CRURNG/ ECRUZ, FCRUZ, HCRUZ, IL1, IL2, IW, PFW
      COMMON/COST/THS, ICOST, FC, TC, DTEMPK,WX, FUELDT

  REAL:: lambda,e,edot,ff
  INTEGER:: igrnd,iths
  REAL:: mnvtas,mxvtas, vtas2
      COMMON/ENERGY/LAMBDA,E,EDOT,FF,IGRND,ITHS,MNVTAS, MXVTAS,VTAS2

      COMMON/VTRCRU/WSS(10), JJCRUZ(10), LLCRUZ(50,10),EECRUZ(50,10),   &
     & DLLDEE(2, 10),IWMAX , WTO, WCRUZ ,CRUZCT,HHCRUZ(50,10),          &
     & FFCRUZ(50,10),JLAST1, JLAST2
      COMMON/VUPDWN/EECLMB(110), FFCLMB(110), DDCLMB(110), EEDOWN(110), &
     & FFDOWN(110), DDDOWN(110), SSCOST(110) ,JCLIMB,JDESCN,TTCLMB(110),&
     & TTDOWN(110)
      COMMON/WINDY/IWIND, PSIA, VWA
      DIMENSION ANS(4)
      DATA T2G,FT2KNT,T3600,SQRRTG,SQR2RO/64.4,1.68781,3600.,65.76,29./
!----------------------------------------------------------------------------
!
!
!     Compute(climb, descend) (fuel, time, distance)
   79 IF( ispliz .eq. 0 )  GO TO 85
      IF ( ioparm .eq. 0) GO TO 840
   81 llast = min0(jclimb, jdescn) - 1
      Do 82 l = 1, llast
      Cost1 = sscost(l)/ abs(sscost(l))
      Cost2 = sscost(l+ 1)/abs(sscost(l+1))
      IF( abs(cost1 - cost2) .gt. .5) GO TO 83
   82 end do
      GO TO 85
   83 pf = -sscost(l)/(sscost(l+1) - sscost(l))
      Ecruz = eeclmb(l) + pf*(eeclmb(l+1) - eeclmb(l))
  840 IF( iter .gt. 1) GO TO 84
      Call serch1(eeclmb, ecruz, jclimb, pfc,lc,limit)
   84 call serch1(eedown,ecruz, jdescn, pfd,ld,limit)
      GO TO 86
   85 pfc = 0.
      Pfd = 0.
      Lc = jclimb
      Ld = jdescn
   86 IF( iter .gt. 1) GO TO 87
      Fclmb = ffclmb(lc) + pfc*(ffclmb(lc+1) - ffclmb(lc))
      Tclmb = ttclmb(lc) + pfc*(ttclmb(lc+1) - ttclmb(lc))
      Dclmb = ddclmb(lc) + pfc*(ddclmb(lc+1) - ddclmb(lc))
   87 fdown = ffdown(ld) + pfd*(ffdown(ld + 1) - ffdown(ld))
      Tdown = ttdown(ld) + pfd*(ttdown(ld + 1) - ttdown(ld))
      Ddown = dddown(ld) + pfd*(dddown(ld + 1) - dddown(ld))
!
!
!
!
!     Compute cruise distance, fuel, time
      IF( ispliz .eq. 0) GO TO 911
      IF( iter .gt. 1) GO TO 95
      IF( ioparm .ne. 0) GO TO 912
   80 dlde1 = dlldee(1, iw)*ecruz + dlldee(2, iw)
      Dlde2= dlldee(1, iw+1)*ecruz + dlldee(2, iw+1)
      Dlde=dlde1 + pfw*(dlde2 - dlde1)
      Scost = sscost(lc) + pfc*(sscost(lc+1) - sscost(lc))
   90 dcruz = abs(scost/dlde)
      GO TO 930
  911 dcruz = range - dclmb - ddown
      IF (iter .gt. 1) GO TO 931
      GO TO 930
  912 dcruz = 0.
      Fcrulb = 0.
      Tcruz = 0.
      Efcruz = 0.
!
!     Initial cruise weight or end condition of climb
  930 IF( iter .gt. 1) GO TO 931
      Wcruzi = wto -  fclmb
      Wcruz = wcruzi
      Ecruzi = ecruz
      Call wlefhv(3, ecruzi, ef, pc, vgknti, vktasi)
      Hcruzi = hcruz
      Fcruzi = fcruz
      Lambdi = lambda
      Call at62(hcruzi,ans)
      Machi = vktasi*ft2knt/(sqrrtg * sqrt(ans(3)+dtempk))
      Vkiasi=sqr2ro*sqrt(ans(2)*((1.+.2*machi*machi)**3.5-1.))/ft2knt
      IF( ioparm .eq. 1) GO TO 95
!
!     Average cruise weight for fuel comsumption  computation
  931 fcrulb= dcruz*fcruz/vgknti
      Wcruza= wto - fclmb - fcrulb/2.
      Wcruz = wcruza
      Call wlefhv(3, ecruza, ef, pc, vgknta, vktasa)
      Fcrulb = dcruz*fcruz  /vgknta
!
!     Final cruise weight or end condition for descend
   95 wcruzf = wcruzi - fcrulb
      Wlandg = wcruzf - fdown
      Wcruz = wcruzf
      Call wlefhv(3, ecruzf, ef, pc, vgkntf, vktasf)
      Hcruzf = hcruz
      Lambdf = lambda
      Lambda = lambda/ft2knt
      Call at62(hcruzf,ans)
      Machf = vktasf*ft2knt/(sqrrtg*sqrt(ans(3)+ dtempk))
      Vkiasf=sqr2ro*sqrt(ans(2)*((1.+.2*machf*machf)**3.5-1.))/ft2knt
!
!
!     PRINT OUT CRUISE DATA
      CALL PAGE
      WRITE(6, 106)
  106 FORMAT(/15X, 'INITIAL',6X,'FINAL', 20X,'INITIAL', 9X,'FINAL'/      &
     & 15X, 'CRUISE', 6X, 'CRUISE', 22X, 'CRUISE', 8X, 'CRUISE')
      WRITE(6, 107) WCRUZI, WCRUZF, VKTASI, VKTASF, LAMBDI, LAMBDF,     &
     & VKIASI, VKIASF, ECRUZI, ECRUZF, VGKNTI, VGKNTF,HCRUZI,HCRUZF,    &
     & MACHI, MACHF
  107 FORMAT(/'WEIGHT(LB)', F11.0, F12.0, 6X, 'TAS', F18.0, F14.0/     &
     & ' COST($/NM)', F11.3, F12.3, 6X, 'IAS', F18.2, F14.2/          &
     & ' ENERGY(FT)', F11.0, F12.0, 6X, 'GR SP KN', F13.2, F14.2/       &
     & ' ALTITUDE', F13.0,F12.0, 6X, 'MACH NO', F14.5, F14.5)
!
!
!
      IF( dcruz .eq. 0.) GO TO 950
      Efcruz =fcrulb/dcruz
      Tcruz = dcruz*t3600/vgknta
  950 tfuel = fclmb + fdown + fcrulb
      Tdist = dclmb + ddown + dcruz
      Ttime = tclmb + tdown + tcruz
      Effcnz = tfuel/tdist
      WRITE(6, 100)
  100 FORMAT(/10X, &
       'FUEL USED(#)    DISTANCE(N M),HR:MIN:SEC,  COST($), $/NM')
      CALL WRITE1(FCLMB, TCLMB, DCLMB, 'CLIMB   ')
      CALL WRITE1(FDOWN, TDOWN, DDOWN, 'DESCEND ')
      CALL WRITE1(FCRULB, TCRUZ, DCRUZ, 'CRUISE  ')
      CALL WRITE1(TFUEL, TTIME, TDIST, 'TOTAL   ')
      WRITE(6, 105) WLANDG
  105 FORMAT(/'LANDING WEIGHT = ', F10.0)
      WRITE(6, 108) EFCRUZ, EFFCNZ
  108 FORMAT(/'CRUISE & OVERALL EFFICIENCY', 2F10.3, '#/NM')
      RETURN
END SUBROUTINE VOPTRJ  ! ----------------------------------------------------



!+
SUBROUTINE WATEST(ECRUZ, LAMBDA, WLNDG, N, WCRUZF, RANGE, INIT, PC, IOPARM)
! ---------------------------------------------------------------------------
! PURPOSE -
!
!
! ESTIMATE CRUISE DISTANCE AND CRUISE WT, COMP FINAL WEIGHT
!     FOR CLIMB-CRUISE DESCEND TRAJ IF ioparm==0
! otherwise,     !     COMPUTE FINAL WEIGHT FOR CLIMB-DESCEND TRAJ

! called by main program

USE Extra,ONLY:PolyE1,Serch1

  REAL,INTENT(IN OUT):: ecruz
  REAL,INTENT(IN):: lambda    ! not used
  REAL,INTENT(OUT):: wlndg
  INTEGER,INTENT(IN):: n         ! not used
  REAL,INTENT(IN):: wcruzf       ! not used
  REAL,INTENT(IN):: range        ! not used
  INTEGER,INTENT(IN):: init      ! not used
  REAL,INTENT(IN):: pc
  INTEGER,INTENT(IN):: ioparm


!      REAL LAMBD1, LAMBD2    ! never used ???
  REAL:: LLCRUZ
!!!      DIMENSION DNFUEL(3), DCOEF(5)
      COMMON/CRURNG/ECRUZX, FCRUZ, HCRUZ, IL1, IL2, IW, PFW
      COMMON/VTRCRU/WSS(10), JJCRUZ(10), LLCRUZ(50,10),EECRUZ(50,10),   &
     & DLLDEE(2, 10),IWMAX , WTO, WCRUZ ,CRUZCT,HHCRUZ(50,10),          &
     & FFCRUZ(50,10),JLAST1, JLAST2
      COMMON/VUPDWN/EECLMB(110), FFCLMB(110), DDCLMB(110), EEDOWN(110), &
     & FFDOWN(110), DDDOWN(110), SSCOST(110) ,JCLIMB,JDESCN,TTCLMB(110),&
     & TTDOWN(110)

  REAL,PARAMETER,DIMENSION(5):: DCOEF = (/ &
    115.8777, -15.07233, 8.081622E-1, -1.865838E-2, 1.527675E-4/)
  REAL,PARAMETER,DIMENSION(3):: DNFUEL = (/689.89, -1.28974, -0.111423/)
!
!----------------------------------------------------------------------------

! Compute climb fuel, estimate descend fuel
  ecruzx = ecruz
  IF (ecruzx .gt. Eeclmb(jclimb)) ecruzx = eeclmb(jclimb)
  Call serch1(eeclmb, ecruzx,jclimb, pf, l, limit)
  Fclmb = ffclmb(l) + pf*(ffclmb(l+1) - ffclmb(l))
  Fdown = (dnfuel(3)*pc +dnfuel(2))*pc + dnfuel(1)

  IF (ioparm==0) then
    Dcruz = polye1(pc, 5, dcoef)
    Call wlefhv(2, ecruz, ef, pc, vgknt, vcktas) ! sets ecruz,ef,vgknt,vcktas
    Fcrulb = dcruz *fcruz/vgknt
    Wlndg = wto - fcrulb - fdown -fclmb
  else
    Wlndg = wto- fclmb - fdown
  END IF
  RETURN
END Subroutine Watest   ! ---------------------------------------------------

!+
SUBROUTINE WIND(H, PSIA, VWA)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: Serch1
  REAL,INTENT(IN):: h
  REAL,INTENT(IN):: psia
  REAL,INTENT(OUT):: vwa

      COMMON/WINDTD/HWIND(24), VW(24), PSIW(24)
!
!----------------------------------------------------------------------------
!
   80 call serch1(hwind, h, 24, pf, I, limit)
      Vwa = vw(i) + pf* (vw(i+1) - vw(i))
      Psiwa = psiw(i) + pf*( psiw(i+1) - psiw(i))
      A =(psia - psiwa)/57.29578
      Vwa =vwa*cos(a)
      RETURN
END Subroutine Wind   ! -----------------------------------------------------


!+
SUBROUTINE WINDIN()
! ---------------------------------------------------------------------------
! PURPOSE -

      COMMON/WINDTD/HWIND(24), VW(24), PSIW(24)

  INTEGER:: i
  REAL:: wind
!----------------------------------------------------------------------------
  hwind(1:24)=2000.0*real( (/ (i, i=0,23) /) )
!  Hwind(1) = 0.
!      Do 100 i= 1, 23
!  100 hwind(i+1) = hwind(i) + 2000.

  Read(7, 90) (psiw(i), vw(i), i= 1, 24)
   90 format(13(f2.0, 1x, f3.1))

  Write(6,*) "wind data"
  Write(6,*) "alt(ft),  vw(knots), vw(ft/sec), psiw(deg)"

  psiw(:)=180.0+10.0*psiw(:)
  psiw(:)=modulo(psiw(:),360.0)

  Do i= 1,24
    Wind = vw(i)
    Vw(i) = vw(i) *1.68781
  !    Psiw(i)=psiw(i)*10.+180.
  !    Psiw(i) = mod(psiw(i), 360.)
    Write(6,102) hwind(i), wind, vw(i), psiw(i)
  END DO
  102 FORMAT(F10.0, 2F10.2, F10.0)

  RETURN
END Subroutine WindIn   ! ---------------------------------------------------

!+
SUBROUTINE WLEFHV(I, ECRUZX, EF, PC, VGKNT, VCKTAS)
! ---------------------------------------------------------------------------
! PURPOSE -

USE Extra,ONLY: Jtrunc,SerchD
  INTEGER,INTENT(IN):: i
  REAL,INTENT(OUT):: ecruzx
  REAL,INTENT(OUT):: ef
  REAL,INTENT(IN):: pc
  REAL,INTENT(OUT):: vgknt
  REAL,INTENT(OUT):: vcktas

      REAL LLCRUZ, LAMBD1, LAMBD2,LAMBDX,LAMBS
      COMMON/CRURNG/ ECRUZ, FCRUZ, HCRUZ, IL1, IL2, IW, PFW

  REAL:: lambda,e,edot,ff
  INTEGER:: igrnd,iths
  REAL:: mnvtas,mxvtas, vtas2
      COMMON/ENERGY/LAMBDA,E,EDOT,FF,IGRND,ITHS,MNVTAS, MXVTAS,VTAS2

      COMMON/IO/ ANS(4), WS(11), EOPTS(11), HSTARS(11), MSTARS(11),    &
     & PISTRS(11), LAMBS(11), VTASOP(11), FUELFL(11), CRDIST(11),       &
     &CRTIME(11), MNCOST(2), MOPIAS(2), MOPTAS(2), MACHOP(2), FDTOPT(2),&
     & OPTALT(2), HOPT(2), EPRS(2),IWMAX,IOPARM, &
        hto,vto,holndg,volndg,eto
      COMMON/VTRCRU/WSS(10), JJCRUZ(10), LLCRUZ(50,10),EECRUZ(50,10),   &
     & DLLDEE(2, 10),IWMAXX, WTO, WCRUZ ,CRUZCT,HHCRUZ(50,10),          &
     & FFCRUZ(50,10),JLAST1, JLAST2
      COMMON/WINDY/IWIND, PSIA, VWA
      DATA T2G,FT2KNT/64.4, 1.68781/
!
!----------------------------------------------------------------------------


      Fpc = (100.0 + pc)/100.0
      IF (i==0) GO TO 102
!      GO TO (101, 102, 101), I


!     Compute new lambda(wcruz) and ecruz(lambda)
  101 call serchd(wss, wcruz, iwmax, pfw, iw, limit)
      Jlast1 = jtrunc(llcruz(1:,iw),   jjcruz(iw))
      Jlast2 = jtrunc(llcruz(1:,iw+1), jjcruz(iw+1))
   81 lambd1 = llcruz(jlast1, iw)*fpc
      Lambd2 = llcruz(jlast2, iw+1)*fpc
      Lambda = lambd1 + pfw*(lambd2 - lambd1)
   82 call serchd(llcruz(1:,iw),lambd1, jlast1, pf1,il1, limit)
      Ecruz1 = eecruz(il1, iw) +pf1*(eecruz(il1 +1, iw) -eecruz(il1,iw))
      Call serchd(llcruz(1:,iw+1), lambd2, jlast2,pf2,il2,limit)
      Ecruz2 = eecruz(il2,iw+1)+pf2*(eecruz(il2+1,iw+1) -eecruz(il2,    &
     & iw+1))
      Ecruz = ecruz1 + pfw*(ecruz2- ecruz1)
      Eoptmx = eopts(iw) + pfw*(eopts(iw+1) - eopts(iw))
      IF (ecruz.lt. (Eoptmx-100.)) GO TO 99
      Ef = ecruz
      GO TO 100
   99 ef = ecruz + 100.
  100 ecruzx = ecruz

      IF( I .EQ. 1) RETURN


!     Compute (fdot, h, vktas) as function of lambda
  102 fcruz1 = ffcruz(il1, iw) +pf1*(ffcruz(il1 +1, iw) -ffcruz(il1,iw))
      Hcruz1 = hhcruz(il1, iw) +pf1*(hhcruz(il1 +1, iw) -hhcruz(il1,iw))
      Fcruz2 = ffcruz(il2,iw+1)+pf2*(ffcruz(il2+1,iw+1) -ffcruz(il2,    &
     & iw+1))
      Hcruz2 = hhcruz(il2,iw+1)+pf2*(hhcruz(il2+1,iw+1) -hhcruz(il2,    &
     & iw+1))
      Fcruz = fcruz1 + pfw*(fcruz2- fcruz1)
      Hcruz = hcruz1 + pfw*(hcruz2- hcruz1)
!
   91 vcktas = sqrt(t2g *(ecruz - hcruz))/ft2knt
      Vgknt = vcktas
      IF (iwind .eq. 0) then
        vgknt=vcktas
      Else
        Call wind(hcruz, psia, vwa)
        Vgknt = vcktas + vwa/ft2knt
      END IF
     
  RETURN
END SUBROUTINE WLEFHV   ! ---------------------------------------------------

!+
SUBROUTINE WRITE1(FUEL, TIME, DIST, LABEL)
! ---------------------------------------------------------------------------
! PURPOSE - Print a message on unit 6

USE Extra,ONLY: Iclock
IMPLICIT NONE
  REAL,INTENT(IN):: fuel
  REAL,INTENT(IN):: time
  REAL,INTENT(IN):: dist
  CHARACTER(LEN=*),INTENT(IN):: label

  REAL:: ths
  INTEGER:: icost
  REAL:: fc, tc, dtempk, wx, fueldt
  COMMON/COST/THS, ICOST, FC, TC, DTEMPK,WX, FUELDT
  REAL:: dolrnm
  REAL:: fcost
  CHARACTER(LEN=*),PARAMETER:: FMT1 = &
    "(/1X, A, 2F10.2, 4X, 2(I2, ':'), I2, 2F10.2)"
  INTEGER:: ihr,imin,isec
  REAL,PARAMETER:: T3600 = 3600.0

!----------------------------------------------------------------------------

  fcost = fc*fuel + tc*time/t3600
  IF (dist==0.0) THEN
    dolrnm=0.0
  ELSE
    Dolrnm = fcost/dist
  END IF

  CALL ICLOCK(TIME, IHR, IMIN, ISEC)
  WRITE(6,FMT1) trim(label), fuel, dist, ihr, imin, isec, fcost, dolrnm

  RETURN
END Subroutine Write1   ! ---------------------------------------------------


!!END Program OptimumTrajectory   ! -------------------------------------

