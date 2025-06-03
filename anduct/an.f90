! PROGRAM AnnularDuct
! ------------------------------------------------------------------------------
! PURPOSE - Calculate the velocity field in an axisymmetric annular duct using
!  the velocity gradient method.
! AUTHORS - Theodore Katsanis, NASA Lewis (later Glenn) Research Center
!           William D. McNally, NASA Lewis (later Glenn) Research Center
!           Ralph L. Carmichael, Public Domain Aeronautical Software
!
! REVISION HISTORY
!   DATE  VERS RERSON  STATEMENT OF CHANGES
!   1977   1.0 TK&WDM  Original coding and release of NASA TN D-8430 and 8431
! Jul1982  1.1   TK    Publication of NASA TP 2029
!   1994   1.2   RLC   Acquired COSMIC program LEW-14000
!   2001   1.3   RLC   Converted to Fortran 90
! 30Nov09  1.4   RLC   Moved all procedures to a module
! 01Dec09  1.5   RLC   Allowed for single or double precision

! NOTE - Original program description from COSMIC:
! Turbomachinery components are often connected by ducts, which are usually
! annular. The configurations and aerodynamic characteristics of these ducts
! are crucial to the optimum performance of the turbomachinery blade rows.
! The ANDUCT computer program was developed to calculate the velocity
! distribution along an arbitrary line between the inner and outer walls of an
! annular duct with axisymmetric swirling flow. Although other programs are
! available for duct analysis, the use of the velocity gradient method makes
! the ANDUCT program fast and convenient while requiring only modest computer
! resources.

! A fast and easy method of analyzing the flow through a duct with axisymmetric
! flow is the velocity gradient method, also known as the stream filament or
! streamline curvature method. This method has been used extensively for blade
! passages but has not been widely used for ducts, except for the radial
! equilibrium equation. In ANDUCT, a velocity gradient equation derived from
! the momentum equation is used to determine the velocity variation along an
! arbitrary straight line between the inner and outer wall of an annular duct.
! The velocity gradient equation is used with an assumed variation of
! meridional streamline curvature. Upstream flow conditions may vary between
! the inner and outer walls, and an assumed total pressure distribution may be
! specified. ANDUCT works best for well-guided passages and where the curvature
! of the walls is small as compared to the width of the passage.





!+
MODULE AnnularDuctProcedures
! ---------------------------------------------------------------------------
! PURPOSE - Collect several widely used procedures for re-use in PDAS software.
! PURPOSE - Collect subroutines Contin, Pabc, and Splint used by AnnularDuct
!  and define some global constants.

IMPLICIT NONE

  INTEGER,PARAMETER:: NREAD=1, NWRIT=2, DBG=3


! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 21Oct01  0.1   RLC   Collected procedures from several sources



! G L O B A L   V A R I A B L E S :
  CHARACTER(LEN=*),PARAMETER:: VERSION = " Version 1.5 (1 November 2009)"

  INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6)    ! single precision
  INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(14)   ! double precision
  INTEGER,PARAMETER:: WP = DP
  
  REAL(WP),PARAMETER:: PI = 3.141592653589793238462643383279502884197_WP
  REAL(WP),PARAMETER:: HALFPI=PI/2, TWOPI=PI+PI
  REAL(WP),PARAMETER:: ZERO = 0.0_WP, HALF = 0.5_WP, ONE=1.0_WP, TWO=2.0_WP
  REAL(WP),PARAMETER:: SIX = 6.0_WP
  REAL(WP),PARAMETER:: DEG2RAD = PI/180,   RAD2DEG = 180/PI


  PUBLIC:: Contin
  PRIVATE:: Pabc
  PUBLIC:: Splint

!----------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE Assert(bool, message)
! ---------------------------------------------------------------------------
! PURPOSE - Display a message if the logical assertion is .FALSE.
!  Inspired by the function in the C programming language.

  LOGICAL,INTENT(IN):: bool   ! logical assertion
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: message
!----------------------------------------------------------------------------
  IF (.NOT. bool) THEN
    WRITE(*,*) "Assertion failed."
    IF (Present(message)) WRITE(*,*) message
  END IF
  RETURN
END Subroutine Assert   ! ---------------------------------------------------

!+
SUBROUTINE FillArray(startValue,endValue, array, spacingCode)
! ------------------------------------------------------------------------------
! PURPOSE - fill an array from start to end. The intermediate points are
!    computed according to various spacing rules.

! NOTE - This subroutine uses module private constants PI and HALFPI. 
!   If you extract this routine for other uses, add them back in.

  REAL(WP),INTENT(IN):: startValue,endValue
  REAL(WP),INTENT(OUT),DIMENSION(:):: array
  INTEGER,INTENT(IN),OPTIONAL:: spacingCode
                              ! =2 full cosine
                              ! =3 half cosine
                              ! =4 half sine
                              ! anything else (or nothing)= uniform spacing
  INTEGER:: k,n
  REAL(WP),ALLOCATABLE,DIMENSION(:):: temp
!-------------------------------------------------------------------------------
  n=SIZE(array)
  IF (n <= 0) RETURN

  array(n)=endValue
  array(1)=startValue
  IF (n <= 2) RETURN

  ALLOCATE(temp(n-2))
  temp= (/ (REAL(k), k=1,n-2) /) / REAL(n-1)

  IF (Present(spacingCode)) THEN
    SELECT CASE(spacingCode)
      CASE (2)
        temp=0.5*(1.0-COS(PI*temp))       ! full cosine, dense near both ends
      CASE (3)
        temp=1.0-COS(HALFPI*temp)             ! half cosine, dense near start
      CASE (4)
        temp=SIN(HALFPI*temp)                     ! half sine, dense near end
    END SELECT
  END IF

  array(2:n-1)=startValue + (endValue-startValue)*temp

  DEALLOCATE(temp)
  RETURN
END Subroutine FillArray   ! ---------------------------------------------------

!+
FUNCTION GetDateTimeStr() RESULT(s)
! ------------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time.
!   You can change the first I2.2 below to I2 if you don't like a leading
!   zero on early morning times, i.e., 6:45 instead of 06:45
  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: MONTH='JanFebMarAprMayJunJulAugSepOctNovDec'
  CHARACTER(LEN=*),PARAMETER:: FMT = '(I2.2,A1,I2.2,I3,A3,I4)'
  CHARACTER(LEN=15):: s
  INTEGER,DIMENSION(8):: v
  INTRINSIC:: DATE_AND_TIME
!-------------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)

  WRITE(s,FMT) v(5), ':', v(6), v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
  RETURN
END FUNCTION GetDateTimeStr   ! ------------------------------------------------

!+
SUBROUTINE OpenInputFile(efu,fname,ext1)
! ------------------------------------------------------------------------------
! PURPOSE - Open an existing file read-only. If the file does not exist,
!   make a new file name by adding the extention ext1
!   Inability to open this file is a fatal error.
  INTEGER,INTENT(IN):: efu              ! # of the external file unit
  CHARACTER(LEN=*),INTENT(IN):: fname   ! name of the file
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: ext1   ! possible extention
  INTEGER:: errCode
  INTRINSIC:: TRIM
!-------------------------------------------------------------------------------
  OPEN(UNIT=efu, FILE=fname, STATUS='OLD', &
    IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
  IF (errCode==0) RETURN

  IF (Present(ext1)) THEN
    OPEN(UNIT=efu, FILE=Trim(fname)//'.'//AdjustL(ext1), STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) RETURN
  END IF

  WRITE(*,*) "FATAL ERROR. Unable to open "//Trim(fname)
  STOP

END Subroutine OpenInputFile   ! -----------------------------------------------

!+
SUBROUTINE OpenOutputFile(efu,fileName)
! ------------------------------------------------------------------------------
! PURPOSE -
  INTEGER,INTENT(IN):: efu   ! external file unit to be opened
  CHARACTER(LEN=*),INTENT(IN):: fileName
  INTEGER:: errCode
!-------------------------------------------------------------------------------
  OPEN(UNIT=efu, FILE=fileName, STATUS='REPLACE', &
      IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode==0) RETURN
  WRITE(*,*) "FATAL ERROR! Unable to open "//Trim(fileName)
  STOP
END Subroutine OpenOutputFile   ! -----------------------------------------------

!+
SUBROUTINE PrintArrays(efu, x,y, a,b,c,d,e,f)
! ------------------------------------------------------------------------------
! PURPOSE - Print arrays, neatly formatted. X and Y are always printed.
!   The next 6 optional variables are printed when present.

  INTEGER,INTENT(IN):: efu   ! unit number of the external file
  REAL(WP),INTENT(IN),DIMENSION(:):: x,y
  REAL(WP),INTENT(IN),DIMENSION(:),OPTIONAL:: a,b,c,d,e,f

  CHARACTER(LEN=80):: buffer
  CHARACTER(LEN=*),PARAMETER:: FMT = '(2F10.4)'
  INTEGER:: k
  INTRINSIC:: SIZE,TRIM
!-------------------------------------------------------------------------------
  buffer=""
  DO k=1,SIZE(x)
    WRITE(buffer(1:20),FMT) x(k),y(k)
    IF (Present(a)) WRITE(buffer(21:30),FMT) a(k)
    IF (Present(b)) WRITE(buffer(31:40),FMT) b(k)
    IF (Present(c)) WRITE(buffer(41:50),FMT) c(k)
    IF (Present(d)) WRITE(buffer(51:60),FMT) d(k)
    IF (Present(e)) WRITE(buffer(61:70),FMT) e(k)
    IF (Present(f)) WRITE(buffer(71:80),FMT) f(k)
    WRITE(efu,*) Trim(buffer)
  END DO
  RETURN
END Subroutine PrintArrays   ! -------------------------------------------------

!+
SUBROUTINE PrintArraysNumbered(efu, x,y, a,b,c,d,e,f)
! ------------------------------------------------------------------------------
! PURPOSE - Print arrays, neatly formatted. X and Y are always printed.
!   The next 6 optional variables are printed when present.

  INTEGER,INTENT(IN):: efu   ! unit number of the external file
  REAL(WP),INTENT(IN),DIMENSION(:):: x,y
  REAL(WP),INTENT(IN),DIMENSION(:),OPTIONAL:: a,b,c,d,e,f

  CHARACTER(LEN=85):: buffer
  CHARACTER(LEN=*),PARAMETER:: FMT = '(2F10.4)'
  INTEGER:: k
  INTRINSIC:: SIZE,TRIM
!-------------------------------------------------------------------------------
  buffer=""
  DO k=1,SIZE(x)
    WRITE(buffer(1:5),'(I5)') k
    WRITE(buffer(6:25),FMT) x(k),y(k)
    IF (Present(a)) WRITE(buffer(26:35),FMT) a(k)
    IF (Present(b)) WRITE(buffer(36:45),FMT) b(k)
    IF (Present(c)) WRITE(buffer(46:55),FMT) c(k)
    IF (Present(d)) WRITE(buffer(56:65),FMT) d(k)
    IF (Present(e)) WRITE(buffer(66:75),FMT) e(k)
    IF (Present(f)) WRITE(buffer(76:85),FMT) f(k)
    WRITE(efu,*) Trim(buffer)
  END DO
  RETURN
END Subroutine PrintArraysNumbered   ! -----------------------------------------




!+
SUBROUTINE Contin(xest,ycalc,ind,jz,ygiv,xdel)
! ---------------------------------------------------------------------------
! PURPOSE - Calculate an estimate of the relative flow velocity for use in 
!  the velocity gradient equation
  REAL(WP),INTENT(IN OUT):: xest
  REAL(WP),INTENT(IN OUT):: ycalc
  INTEGER,INTENT(IN OUT):: ind
  INTEGER,INTENT(IN OUT):: jz
  REAL(WP),INTENT(IN OUT):: ygiv
  REAL(WP),INTENT(IN OUT):: xdel

  REAL(WP):: acb2
  REAL(WP):: apa,bpb,cpc
  REAL(WP):: discr
  REAL(WP),DIMENSION(3):: x,y
  REAL(WP):: xorig
  REAL(WP):: xoshft
  INTEGER,SAVE:: ncall = 0
!----------------------------------------------------------------------------

  ncall = ncall+1
  IF (ind.NE.1 .AND. ncall.GT.100) GO TO 160

      GO TO (10,30,40,50,60,110,150),ind

!--FIRST CALL
   10 ncall = 1
      xorig = xest
      IF (ycalc.GT.YGIV .AND. jz.EQ.1) GO TO 20
      ind = 2
      y(1) = ycalc
      x(1) = 0.
      xest = xest+xdel
      RETURN
   20 ind = 3
      y(3) = ycalc
      x(3) = 0.
      xest = xest-xdel
      RETURN

!--SECOND CALL
   30 ind = 4
      y(2) = ycalc
      x(2) = xest-xorig
      xest = xest+xdel
      RETURN

   40 ind = 5
      y(2) = ycalc
      x(2) = xest-xorig
      xest = xest-xdel
      RETURN

!--THIRD OR LATER CALL - FIND SUBSONIC OR SUPERSONIC SOLUTION
   50 y(3) = ycalc
      x(3) = xest-xorig
      GO TO 70

   60 y(1) = ycalc
      x(1) = xest-xorig
   70 IF (ygiv.LT.MIN(y(1),y(2),y(3))) GO TO (120,130),jz
      ind = 6   !!! was labelled statement 80
      CALL Pabc(x,y,apa,bpb,cpc)
      discr = bpb**2-4.*apa*(cpc-ygiv)
      IF (discr.LT.0.) GO TO 140
      IF (ABS(400.*apa*(cpc-ygiv)).LE.bpb**2) GO TO 90

      xest = -bpb-SIGN(SQRT(discr),apa)
      IF (jz.EQ.1 .AND. apa.GT.0. .AND. y(3).GT.y(1)) xest=-bpb+SQRT(discr)
      IF (jz.EQ.2 .AND. apa.LT.0.) xest = -bpb-sqrt(discr)
      xest = xest/2./apa
      GO TO 100
   90 IF (jz.EQ.2 .AND. bpb.GT.0.) GO TO 130
      acb2 = apa/bpb*(cpc-ygiv)/bpb
      IF (ABS(acb2).LE.1.E-8) acb2=0.
      xest = -(cpc-ygiv)/bpb*(1.+acb2+2.*acb2**2)
  100 IF (xest.GT.x(3)) GO TO 130
      IF (xest.LT.x(1)) GO TO 120
      xest = xest+xorig
      RETURN

!--FOURTH OR LATER CALL - NOT CHOKED
  110 IF(xest-xorig.GT.x(3)) GO TO 130
      IF(xest-xorig.LT.x(1)) GO TO 120
      y(2) = ycalc
      x(2) = xest-xorig
      GO TO 70

!--THIRD OR LATER CALL - SOLUTION EXISTS,
!--BUT RIGHT OR LEFT SHIFT REQUIRED
  120 ind = 5
!--LEFT SHIFT
      xest = x(1)-xdel+xorig
      xoshft = xest-xorig
      xorig = xest
      y(3) = y(2)
      x(3) = x(2)-xoshft
      y(2) = y(1)
      x(2) = x(1)-xoshft
      RETURN
  130 ind = 4
!--RIGHT SHIFT
      xest = x(3)+xdel+xorig
      xoshft = xest-xorig
      xorig = xest
      y(1) = y(2)
      x(1) = x(2)-xoshft
      y(2) = y(3)
      x(2) = x(3)-xoshft
      RETURN
!--THIRD OR LATER CALL - APPEARS TO BE CHOKED
  140 xest = -bpb/2./apa
      ind = 7
      IF (xest.LT.x(1)) GO TO 120
      IF (xest.GT.x(3)) GO TO 130
      xest = xest+xorig
      RETURN

!--FOURTH OR LATER CALL - PROBABLY CHOKED
  150 IF (ycalc.GE.ygiv) GO TO 110
      ind = 10
      RETURN

!--NO SOLUTION FOUND IN 100 ITERATIONS
  160 ind = 11

  RETURN
END Subroutine Contin   ! ------------------------------------------------------


!+
SUBROUTINE Pabc(x,y,a,b,c)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate the coefficients a,b,c of the parabola
!     y=a*x**2+b*x+c, passing through the given x,y points

IMPLICIT NONE
  REAL(WP),INTENT(IN),DIMENSION(3):: x,y   ! the three input points
  REAL(WP),INTENT(OUT):: a,b,c   ! the polynomial coefficients

  REAL(WP):: c1,c2
!----------------------------------------------------------------------------
  c1 = x(3)-x(1)
  c2 = (y(2)-y(1))/(x(2)-x(1))
  a = (c1*c2-y(3)+y(1))/c1/(x(2)-x(3))
  b = c2-(x(1)+x(2))*a
  c = y(1)-x(1)*b-x(1)**2*a
  RETURN
END Subroutine Pabc


!+
SUBROUTINE Splint (x,y,n,z,maxx,yint,dydx,d2ydx2)
! ---------------------------------------------------------------------------
! PURPOSE - Calculate the interpolated points and derivatives for a spline 
!  curve. The end condition is that the second derivatives at end points are
!  sdr1 and sdrn times the second derivatives at adjacent points.

  REAL(WP),INTENT(IN),DIMENSION(:):: x,y
  INTEGER,INTENT(IN):: n
  INTEGER,INTENT(IN):: maxx
  REAL(WP),INTENT(IN),DIMENSION(:):: z
  REAL(WP),INTENT(OUT),DIMENSION(:):: yint,dydx,d2ydx2

!!!      COMMON NREAD,NWRIT
!!!      DIMENSION X(N),Y(N),Z(MAX),YINT(MAX),DYDX(MAX),D2YDX2(MAX)
   
  REAL(WP):: a
  REAL(WP):: acb2
  REAL(WP):: apa
  REAL(WP):: bpb
  REAL(WP):: c
  REAL(WP):: cpc
  REAL(WP),DIMENSION(101):: em,g,sb
  REAL(WP):: f
  INTEGER:: i
  INTEGER:: ierr
  INTEGER:: k
  INTEGER:: no
  REAL(WP):: s2
  REAL(WP):: sdr1,sdrn
  REAL(WP):: sk,sn
  REAL(WP):: toler
  REAL(WP):: w
  REAL(WP):: y0,ynp1
   
!!!      DIMENSION G(101),SB(101),EM(101)
!----------------------------------------------------------------------------
  write(DBG,*) 'entering splint, n,maxx=', n,maxx

      ierr = 0
      sdr1 = .5
      sdrn = .5
!!!      write(*,*) 'n=', n
      toler= 1E-5*ABS(x(n)-x(1))/REAL(n)
      c = x(2)-x(1)
      IF (c.EQ.0.) GO TO 130
      
      
      sb(1) = -sdr1
      g(1) = 0.
      no = n-1
      IF (no.LE.0) GO TO 140     ! n=1 ERROR
      
!!!      IF (no.EQ.1) GO TO 20      ! n=2 TRIVIAL

  DO i=2,no
    a = c
    c = x(i+1)-x(i)
    IF (a*c.EQ.0.) GO TO 130
    IF (a*c.LT.0.) ierr = 1
    w = 2.*(a+c)-a*sb(i-1)
    sb(i) = c/w
    f = (y(i+1)-y(i))/c-(y(i)-y(i-1))/a
    g(i) = (6.0*f-a*g(i-1))/w
  END DO

!!!   20 em(n) = sdrn*g(n-1)/(1.+sdrn*sb(n-1))
  em(n) = sdrn*g(n-1)/(1.+sdrn*sb(n-1))

  DO i=2,n
    k = n+1-i
    em(k) = g(k)-sb(k)*em(k+1)
  END DO
  IF (maxx.LE.0) RETURN

!!!      ENTRY SPLENT (Z,MAX,YINT,DYDX,D2YDX2)

      DO 120 i=1,maxx
      k=2
      IF (ABS(z(i)-x(1)).LT.toler) GO TO 40
      IF (z(i).GT.2.0*x(1)-x(2)) GO TO 50
      GO TO 80

   40 yint(i) = y(1)
      sk = x(k)-x(k-1)
      GO TO 110

   50 IF (ABS(z(i)-x(k)).LT.toler) GO TO 60
      IF (z(i).GT.x(k)) GO TO 70
      GO TO 100

   60 yint(i) = y(k)
      sk = x(k)-x(k-1)
      GO TO 110

   70 IF (k.GE.n) GO TO 90
      k = k+1
      GO TO 50

   80 s2 = x(2)-x(1)
      y0 = em(1)*s2**2+2.*y(1)-y(2)
      dydx(i) = (y(2)-y(1))/s2-7.*em(1)/6.*s2
      yint(i) = y0+dydx(i)*(z(i)-x(1)+s2)
      d2ydx2(i) = 0.
      CYCLE   !!! GO TO 120

   90 IF (z(i).lt.2.*x(n)-x(n-1)) GO TO 100
      sn = x(n)-x(n-1)
      ynp1 = em(n)*sn**2+2.*y(n)-y(n-1)
      dydx(i) = (y(n)-y(n-1))/sn+7.*em(n)/6.*sn
      yint(i) = ynp1+dydx(i)*(z(i)-x(n)-sn)
      d2ydx2(i) = 0.
      CYCLE   !!! GO TO 120

  100 sk = x(k)-x(k-1)
      yint(i) = em(k-1)*(x(k)-z(i))**3/6./sk  + &
        em(k)*(z(i)-x(k-1))**3/6.0/sk +         &
        (y(k)/sk  -em(k)*sk/6.)*(z(i)-x(k-1)) + &
        (y(k-1)/sk-em(k-1)*sk/6.)*(x(k)-z(i))

  110 dydx(i)=-em(k-1)*(x(k)-z(i))**2/2.0/sk  + &
       em(k)*(x(k-1)-z(i))**2/2.0/sk + &
       (y(k)-y(k-1))/sk  - &
       (em(k)-em(k-1))*sk/6.
      d2ydx2(i) = em(k)-(x(k)-z(i))/sk*(em(k)-em(k-1))
  120 END DO


      IF (ierr.EQ.0) RETURN

  130 WRITE(NWRIT,*) "SPLINT ERROR -- ONE OF THREE POSSIBLE CAUSES"
      WRITE(NWRIT,*) "   ADJACENT X POINTS ARE DUPLICATES OF EACH OTHER"
      WRITE(NWRIT,*) "   SOME X POINTS ARE OUT OF SEQUENCE"
      WRITE(NWRIT,*) "   SOME X POINTS ARE UNDEFINED"
      WRITE(NWRIT,*) "NUMBER OF POINTS = ", n
      WRITE(NWRIT,1020) N,(X(I),Y(I),I=1,N)
      IF (IERR.EQ.0) STOP
      IERR = 0
      RETURN

  140 WRITE(NWRIT,1010)
      WRITE(NWRIT,*) "    NUMBER OF POINTS = ", n
      WRITE(NWRIT,1020) N,(X(I),Y(I),I=1,N)
      STOP

! 1000 FORMAT (1H1,10X,44HSPLINT ERROR -- ONE OF THREE POSSIBLE CAUSES/  &
!     &17X,51H1.  ADJACENT X POINTS ARE DUPLICATES OF EACH OTHER./       &
!     &17X,38H2.  SOME X POINTS ARE OUT OF SEQUENCE./                    &
!     &17X,32H3.  SOME X POINTS ARE UNDEFINED.)

 1010 FORMAT('SPLINT ERROR -- NUMBER OF SPLINE POINTS GIVEN  LESS THAN TWO',//)

 1020 FORMAT (//17X,"X  ARRAY",6X,"Y  ARRAY"/(17X,2ES13.5))

END Subroutine Splint   ! ------------------------------------------------------

END Module AnnularDuctProcedures   ! =========================================================



!+
PROGRAM AnnularDuct
! ---------------------------------------------------------------------------
! PURPOSE - Solve for the freestream velocity gradient equation for an
!   annular duct with variable whirl and stagnation pressure.


!
USE AnnularDuctProcedures

 IMPLICIT NONE

!!!      COMMON NREAD,NWRIT
  REAL(WP),DIMENSION(50):: strfn,qdist,zlamda,tip,rhoip
!!!      DIMENSION STRFN(50),QDIST(50),ZLAMDA(50),TIP(50),RHOIP(50)
      
  INTEGER,PARAMETER:: NDIM=101
  
  REAL(WP),DIMENSION(NDIM):: aaa,bbb

  REAL(WP),DIMENSION(NDIM):: r,zlmint,tint,rhoint,samp
  REAL(WP),DIMENSION(NDIM):: a1,a2,a3,e,f,g1
  REAL(WP),DIMENSION(NDIM):: delprs,dellam,deltip
  REAL(WP),DIMENSION(NDIM):: rcardq,wsubm,u,q
  REAL(WP),DIMENSION(NDIM):: al,curv,b1,w
  REAL(WP),DIMENSION(NDIM):: beta,pres
  REAL(WP),DIMENSION(NDIM):: unew,wwcr

  REAL(WP),DIMENSION(NDIM-1):: whub,wtip,wwcrh,wwcrt,wmh,wmt,beth,bett,prsh,prst,hdist,tdist

        
  
!      DIMENSION R(101),ZLMINT(101),TINT(101),RHOINT(101),SAMP(101),     &
!     &   A1(101),A2(101),A3(101),E(101),F(101),G1(101),DELPRS(101),     &
!     &   DELLAM(101),DELTIP(101),RCARDQ(101),WSUBM(101),U(101),Q(101),  &
!     &   AL(101),CURV(101),B1(101),W(101),BETA(101),PRES(101),UNEW(101),&
!     &   WWCR(101)
!      DIMENSION WHUB(100),WTIP(100),WWCRH(100),WWCRT(100),WMH(100),     &
!     &   WMT(100),BETH(100),BETT(100),PRSH(100),PRST(100),HDIST(100)
 
  REAL(WP):: acoef
  REAL(WP):: alh,alt
  REAL(WP):: ar
  REAL(WP):: camp
  REAL(WP):: cp
  REAL(WP):: curvh
  REAL(WP):: curvt
  CHARACTER(LEN=15):: dateTimeStr
  REAL(WP):: delal
  REAL(WP):: delman
  REAL(WP):: delmax
  REAL(WP):: delq
  REAL(WP):: denom
  REAL(WP):: errmax
  REAL(WP):: expon
  CHARACTER(LEN=80):: fileName
  REAL(WP):: gam
  REAL(WP):: grt

  INTEGER:: i
  INTEGER:: ind
  INTEGER:: iprint
  INTEGER:: istat
  INTEGER:: jz
  INTEGER:: lsfr
  INTEGER:: mson
  INTEGER:: msup
  INTEGER:: nadd
  INTEGER:: ncount
  INTEGER:: next
  INTEGER:: nht
  INTEGER:: nsub

  REAL(WP):: prsint
  REAL(WP):: prsnxt
  REAL(WP):: psi
  REAL(WP):: qht
  REAL(WP):: rhlast
  REAL(WP):: rhub
  REAL(WP):: rtip
  REAL(WP):: rtlast
  REAL(WP),PARAMETER:: RTOLER = 1E-3_WP
  REAL(WP):: rva
  REAL(WP):: rvas
  REAL(WP),PARAMETER:: SBAND = 0.01_WP
  CHARACTER(LEN=80):: stitle
  REAL(WP):: temp
  CHARACTER(LEN=80):: title
  REAL(WP):: tgrgp1
  REAL(WP):: ttip
  REAL(WP):: was
  REAL(WP):: wass
  REAL(WP):: wcr
  REAL(WP):: wm2
  REAL(WP):: wmhub
  REAL(WP):: wmmin = 0.0
  REAL(WP):: wmmax = 0.0
  REAL(WP):: wsq
  REAL(WP):: wth
  REAL(WP):: wth2
  REAL(WP):: zhlast
  REAL(WP):: zhub
  REAL(WP):: zmnflo
  REAL(WP):: zmsfl
  REAL(WP):: zmsfle
  REAL(WP):: zmxflo
  REAL(WP):: ztip
  REAL(WP):: ztlast

!----------------------------------------------------------------------------
  dateTimeStr=GetDateTimeStr()
  OPEN(UNIT=DBG,FILE='anduct.dbg',STATUS='REPLACE',ACTION='WRITE')
  WRITE(DBG,*) 'Annular Duct Program'
  WRITE(DBG,*) 'PDAS'//VERSION
  WRITE(DBG,*) 'Date/time of run: '//dateTimeStr
  OPEN(UNIT=NWRIT,FILE='anduct.out',STATUS='REPLACE',ACTION='WRITE')
!!  CALL OpenOutputFile(NWRIT,'anduct.out')
!!  CALL OpenOutputFile(DBG, 'anduct.dbg')

    5 ISTAT = 0
      HDIST(1) = 0.
      TDIST(1) = 0.

  WRITE(*,*) "Enter the name of the input file:"
  READ(*,'(A)') fileName
  CALL OpenInputFile(NREAD,fileName)
  WRITE(DBG,*) 'Reading from '//TRIM(fileName)


      READ(NREAD,'(A)') title

   10 WRITE(NWRIT,1100)
      istat = istat+1
      WRITE(NWRIT,*) title

      READ (NREAD,'(A)') stitle
      WRITE(NWRIT,*) stitle

      READ (NREAD,1010) nht,lsfr,iprint,next
      WRITE(NWRIT,1120) nht,lsfr,iprint,next

      READ (NREAD,1020) gam,ar,zmsfl
      WRITE(NWRIT,1130) gam,ar,zmsfl

      READ (NREAD,1020) rhub,rtip,zhub,ztip,curvh,curvt,alh,alt
      WRITE(NWRIT,1140) rhub,rtip,zhub,ztip,curvh,curvt,alh,alt

      IF(lsfr.NE.0) GO TO 20
      READ (NREAD,1020) (strfn(i),i=1,nht)
      WRITE(NWRIT,1150) (strfn(i),i=1,nht)
      GO TO 30

   20 READ (NREAD,1020) (qdist(i),i=1,nht)
      WRITE(NWRIT,1160) (qdist(i),i=1,nht)

   30 READ (NREAD,1020) (zlamda(i),i=1,nht)
      WRITE(NWRIT,1170) (zlamda(i),i=1,nht)
      READ (NREAD,1020) (tip(i),i=1,nht)
      WRITE(NWRIT,1180) (tip(i),i=1,nht)
      READ (NREAD,1020) (rhoip(i),i=1,nht)
      WRITE(NWRIT,1190) (rhoip(i),i=1,nht)
!
!--ALL INPUT HAS BEEN READ
!
!
!--CALCULATE CONSTANTS AND INTERPOLATE AS NECESSARY, CONVERT ALH + ALT
!     TO RADIANS
!
      cp = ar/(gam-1.)*gam
      expon = 1./(gam-1.)
      tgrgp1 = 2.*gam*ar/(gam+1.)
      psi = ATAN2(zhub-ztip,rtip-rhub)
      qht = SQRT((rtip-rhub)**2+(ztip-zhub)**2)
      delq = qht/100.
      alh = alh/180.*pi
      alt = alt/180.*pi
      delal = (alt-alh)/100.
      jz = 1
!
!--initialize and make initial guess for wmhub
!
      ncount = 0
      nadd = 0
      nsub = 0
      zmxflo = 0.
      zmnflo = 1.E10
      wmhub = zmsfl/rhoip(1)/pi/(rhub+rtip)/qht
      wcr = SQRT(tgrgp1*tip(1))
      wmhub = MIN(wmhub,wcr)
      delmax = wmhub/20.
      wmmin = -1.E30    ! was -1E50

      DO 40 i=1,101
      r(i) = rhub+REAL(i-1)/100.*(rtip-rhub)
      q(i) = REAL(i-1)*delq
      u(i) = REAL(i-1)/100.*(r(i)+rhub)/(rtip+rhub)
      al(i) = alh+REAL(i-1)*delal
      curv(i) = curvh+REAL(i-1)/100.*(curvt-curvh)
   40 END DO
     
      IF (LSFR.NE.0) GO TO 60
!
!--Stream function given as input(lsfr = 0)
!
   50 CALL Splint(strfn,zlamda,nht,u,101,zlmint,aaa,bbb)
      WRITE(DBG,*) '50: STRFN and ZLAMBDA'
      CALL PrintArraysNumbered(DBG,strfn(1:nht),zlamda(1:nht))
      WRITE(DBG,*) '50: U,ZLMINT,AAA,BBB'
      CALL PrintArraysNumbered(DBG,u,zlmint,aaa,bbb) 
   
      CALL Splint(strfn,tip   ,nht,u,101,tint  ,aaa,bbb)
      WRITE(DBG,*) '50: STRFN and TIP'
      CALL PrintArraysNumbered(DBG,strfn(1:nht),tip(1:nht))
      WRITE(DBG,*) '50: U,TINT,AAA,BBB'
      CALL PrintArraysNumbered(DBG,u,tint,aaa,bbb) 
   
      CALL Splint(strfn,rhoip ,nht,u,101,rhoint,aaa,bbb)
      WRITE(DBG,*) '50: STRFN and RHOIP'
      CALL PrintArraysNumbered(DBG,strfn(1:nht),rhoip(1:nht))
      WRITE(DBG,*) '50: U,RHOINT,AAA,BBB'
      CALL PrintArraysNumbered(DBG,u,rhoint,aaa,bbb) 
      GO TO 70
!
!--Radius given as input(lsfr = 1)
!
   60 CALL Splint(qdist,zlamda,nht,q,101,zlmint,aaa,bbb)
      WRITE(DBG,*) '60: QDIST and ZLAMBDA'
      CALL PrintArraysNumbered(DBG,qdist(1:nht),zlamda(1:nht))
      WRITE(DBG,*) '60: Q,ZLMINT,AAA,BBB'
      CALL PrintArraysNumbered(DBG,q,zlmint,aaa,bbb) 
  
      CALL Splint(qdist,tip,nht,q,101,tint,aaa,bbb)
      WRITE(DBG,*) '60: QDIST and TIP'
      CALL PrintArraysNumbered(DBG,qdist(1:nht),tip(1:nht))
      WRITE(DBG,*) '60: Q,TINT,AAA,BBB'
      CALL PrintArraysNumbered(DBG,q,tint,aaa,bbb) 
      
      CALL Splint(qdist,rhoip,nht,q,101,rhoint,aaa,bbb)
      WRITE(DBG,*) '60: QDIST and RHOIP'
      CALL PrintArraysNumbered(DBG,qdist(1:nht),rhoip(1:nht))
      WRITE(DBG,*) '60: Q,RHOINT,AAA,BBB'
      CALL PrintArraysNumbered(DBG,q,rhoint,aaa,bbb) 
!
!--Calculate quantities which are independent of w-sub-m
!
   70 prsint = rhoint(1)*ar*tint(1)

      DO 80 i=1,101
      samp(i) = SIN(al(i)-psi)
      camp = COS(al(i)-psi)
      a1(i) = camp*curv(i)
      a2(i) = samp(i)/camp*curv(i)
      a3(i) = -SIN(al(i))/r(i)
      b1(i) = -samp(i)/camp
      e(i) = -zlmint(i)/r(i)/r(i)
      f(i) = (zlmint(i)/r(i))**2/2./tint(i)
      g1(i) = ar/prsint
      rcardq(i) = rhoint(i)*camp*2.*pi*r(i)*delq
      IF (i.GE.101) CYCLE   !!! GO TO 80
      prsnxt = rhoint(i+1)*ar*tint(i+1)
      delprs(i) = prsnxt-prsint
      prsint = prsnxt
      dellam(i) = zlmint(i+1)-zlmint(i)
      deltip(i) = tint(i+1)-tint(i)
   80 END DO

!--Solve velocity gradient equation on given straight line from hub to tip

!--Start or restart iteration procedure with ind = 1

  100 ind = 1
!--statements 110 to 140 perform a single iteration
  110 wsubm(1) = wmhub
      ncount = ncount+1
      msup = 0
      mson = 0
      zmsfle = 0.
      wth2= (zlmint(1)/r(1))**2
      wm2 = wsubm(1)**2
      temp = tint(1)-(wth2+wm2)/2./cp
      IF(temp.LE.0. .AND. wsubm(1).LT.0.) GO TO 150
      IF(temp.LE.0.) GO TO 160
      ttip = temp/tint(1)
      rva = ttip**expon*wsubm(1)*rcardq(1)


      DO 140 I=1,100                  ! start of loop
      grt = gam*ar*temp
      if(wm2.gt.grt) msup = 1
      denom = 1.-wm2/grt
      IF(ABS(denom).GT.sband) GO TO 120
      mson = 1
      denom = sband
  120 acoef = a1(i)+samp(i)/denom*(a2(i)+a3(i)*(1.+wth2/grt))
      was = wsubm(i)*(1.+acoef*delq+b1(i)/denom*delal+deltip(i)/2./tint(i))
      was = was+(e(i)*dellam(i)+f(i)*deltip(i)+g1(i)*temp*delprs(i))/wsubm(i)
      IF (was.LT.ABS(wmhub/1000.)) GO TO 150
      wm2 = was**2
      wth2 = (zlmint(i+1)/r(i+1))**2
      temp = tint(i+1)-(wth2+wm2)/2./cp
      IF (temp.LE.0. .AND. was.LT.0.) GO TO 150
      IF (temp.LE.0.) GO TO 160
      grt = gam*ar*temp
      IF (wm2.GT.grt) msup = 1
      denom = 1.-wm2/grt
      IF (ABS(denom).GT.sband) GO TO 130
      mson = 1
      denom = sband
  130 acoef = a1(i+1)+samp(i+1)/denom*(a2(i+1)+a3(i+1)*(1.+wth2/grt))
      wass = wsubm(i)+was*(acoef*delq+b1(i+1)/denom*delal+deltip(i)/2./ &
     &   tint(i+1))
      wass = wass+(e(i+1)*dellam(i)+f(i+1)*deltip(i)+g1(i+1)*temp*      &
     &   delprs(i))/was
      IF (wass.LT.ABS(wmhub/1000.)) GO TO 150
      wsubm(i+1) = (was+wass)/2.
      wm2 = wsubm(i+1)**2
      temp = tint(i+1)-(wth2+wm2)/2./cp
      IF (temp.LE.0. .AND. wsubm(i+1).LT.0.) GO TO 150
      IF (temp.LE.0.) GO TO 160
      IF (wsubm(i+1).LT.0.) GO TO 150
      ttip = temp/tint(i+1)
      rvas = ttip**expon*wsubm(i+1)*rcardq(i+1)
      zmsfle = (rva+rvas)/2.+zmsfle
      unew(i+1) = zmsfle/zmsfl
  140 rva = rvas


      IF (ind.GE.6 .AND. ABS(zmsfl-zmsfle).LE.zmsfl*rtoler) GO TO 200
      zmxflo = MAX(zmsfle,zmxflo)
      zmnflo = MIN(zmsfle,zmnflo)

      IF (wmmin.LT.-1.E30) wmmax = wmhub    ! was E40
      IF (wmmin.LT.-1.E30) wmmin = wmhub    ! was E40
      wmmax = MAX(wmhub,wmmax)
      wmmin = MIN(wmhub,wmmin)
      CALL Contin(wmhub,zmsfle,ind,jz,zmsfl,delmax)
      IF (wmhub.EQ.0.) GO TO 300
      IF (ind.LT.10) GO TO 110
!--IND = 10 INDICATES CHOKED FLOW
      IF (ind.EQ.10) GO TO 200
!--IND = 11 INDICATES NO SOLUTION FOUND IN 100 ITERATIONS
      GO TO 300


  150 wmhub = wmhub+.501*delmax
      nadd = nadd+1
      IF (ncount.LT.1000) GO TO 100
      GO TO 300
  160 wmhub = wmhub-.499*delmax
      nsub = nsub+1
      IF (ncount.LT.1000) GO TO 100
      GO TO 300
!
!--CONVERGED INNER ITERATION
!
  200 CONTINUE
      IF(lsfr.NE.0) GO TO 220
      errmax = 0.
  DO i=2,101   ! was 201 loop
    IF (lsfr.EQ.0) errmax = MAX(errmax,ABS(unew(i)-u(i)))
    u(i) = unew(i)
  END DO
  IF (errmax.GT.rtoler .AND. ncount.LT.1000) GO TO 50
!
!--CONVERGED OUTER ITERATION - PRINT DETAILED OUTPUT IF REQUESTED
!
  220 IF (iprint.NE.0) WRITE(NWRIT,1100)
      IF (iprint.NE.0) WRITE(NWRIT,*) STITLE
      IF (errmax.GT.RTOLER) WRITE(NWRIT,2030) ERRMAX
      IF (ind.EQ.10) WRITE(NWRIT,2000) ZMSFLE
      IF (msup.EQ.1) WRITE(NWRIT,2010)
      IF (mson.EQ.1) WRITE(NWRIT,2020)

  DO I=1,101   ! was 230 loop
    wth = zlmint(i)/r(i)
    wsq = wsubm(i)**2+wth**2
    w(i) = SQRT(wsq)
    wwcr(i) = w(i)/SQRT(tgrgp1*tint(i))
    beta(i) = ATAN2(wth,wsubm(i))/PI*180.
    prsint = rhoint(i)*ar*tint(i)
    ttip = 1.-wsq/2./cp/tint(i)
    pres(i) = prsint*ttip**(gam*expon)
  END DO

  IF (istat.GT.1) hdist(istat)=hdist(istat-1)+SQRT((rhub-rhlast)**2+(zhub-zhlast)**2)
  IF (istat.GT.1) tdist(istat)=tdist(istat-1)+SQRT((rtip-rtlast)**2+(ztip-ztlast)**2)
  rhlast = rhub
  rtlast = rtip
  zhlast = zhub
  ztlast = ztip
  whub(istat) = w(1)
  wtip(istat) = w(101)
  wwcrh(istat) = wwcr(1)
  wwcrt(istat) = wwcr(101)
  wmh(istat) = wsubm(1)
  wmt(istat) = wsubm(101)
  beth(istat) = beta(1)
  bett(istat) = beta(101)
  prsh(istat) = pres(1)
  prst(istat) = pres(101)
  IF (iprint.EQ.0) GO TO 500
  WRITE(NWRIT,1200) (i,w(i),wwcr(i),wsubm(i),beta(i),pres(i),u(i),i=1,101,2)
  GO TO 500
!
!--Print error messages with information if a satisfactory solution
!     cannot be obtained
  300 CONTINUE
      WRITE(NWRIT,1100)
      WRITE(NWRIT,*) stitle
      IF (ind.EQ.11) WRITE(NWRIT,2040)
      IF (ncount.GE.1000) WRITE(NWRIT,2050)
      WRITE(NWRIT,2060) zmxflo,zmnflo,wmmax,wmmin,ncount
      IF (nsub.GT.0) WRITE(NWRIT,2070) nsub
      IF (nadd.GT.0) WRITE(NWRIT,2080) nadd
  500 IF (next.EQ.1 .AND. istat.LT.100) GO TO 10
!
!--Print summary for a given geometry
!
      IF(next.EQ.1 .AND. istat.GE.100) WRITE(NWRIT,2110)
      WRITE(NWRIT,1100)
      WRITE(NWRIT,*) title
      WRITE(NWRIT,2130) istat
      WRITE(NWRIT,2120) (hdist(i),whub(i),wwcrh(i),wmh(i),beth(i),      &
     &   prsh(i),i=1,istat)
      WRITE(NWRIT,2140) istat
      WRITE(NWRIT,2120) (tdist(i),wtip(i),wwcrt(i),wmt(i),bett(i),      &
     &   prst(i),i=1,istat)
!
!--Continue with next case if any
!
      IF (next.EQ.0) STOP
      IF (next.EQ.-1) GO TO 5   ! go back and read another complete case
      ISTAT = 0
      GO TO 10                  ! go back and read another station for this case

!--FORMAT STATEMENTS
!!! 1000 FORMAT(20A4)
 1010 FORMAT(16I5)
 1020 FORMAT(8F10.0)
 1100 FORMAT(1H1)
!!! 1110 FORMAT(20X,20A4)
 1120 FORMAT(6X,'NHT',3X,'LSFR',1X,'IPRINT',3X,'NEXT'/2X,16I7)
 1130 FORMAT(5X,'GAM',13X,'AR',12X,'ZMSFL'/1X,8G16.7)
 1140 FORMAT(5X,'RHUB',12X,'RTIP',12X,'ZHUB',12X,'ZTIP',11X,'CURVH',     &
     &   11X,'CURVT',12X,'ALH',13X,'ALT'/1X,8G16.7)
 1150 FORMAT(10X,'STRFN ARRAY'/(1X,8G16.7))
 1160 FORMAT(10X,'QDIST ARRAY'/(1X,8G16.7))
 1170 FORMAT(10X,'ZLAMDA ARRAY'/(1X,8G16.7))
 1180 FORMAT(10X,'TIP ARRAY'/(1X,8G16.7))
 1190 FORMAT(10X,'RHOIP ARRAY'/(1X,8G16.7))
 1200 FORMAT(3X,'I',7X,'V',13X,'V/VCR',11X,'VSUBM',12X,'BETA',6X,        &
     &   'STATIC PRESSURE',3X,'STREAM FUNCTION'/(2X,I3,6G16.7))
 2000 FORMAT(' ******************************************************'/  &
     & ' THE PASSAGE IS CHOKED AT THIS STATION'/                         &
     & ' THE CHOKING MASS FLOW IS ',G12.4/                               &
     &       ' ******************************************************')

 2010 FORMAT(10X,                                                        &
     & ' SUPERSONIC MERIDIONAL VELOCITY COMPONENT AT THIS STATION')
 2020 FORMAT(' ******************************************************'/  &
     & '  SONIC MERIDIONAL VELOCITY COMPONENT OCCURS AT THIS STATION'/   &
     & '  THIS MAY RESULT IN AN INACCURATE SOLUTION'/                    &
     &       ' ******************************************************')
 2030 FORMAT(' ******************************************************'/  &
     &       '  FULLY CONVERGED SOLUTION COULD NOT BE OBTAINED'/         &
     &       '  IN 1000 ITERATIONS AT THIS STATION'/                     &
     &       '  THE STREAM FUNCTION CHANGED BY',F6.3/                    &
     &       '  BETWEEN THE LAST TWO ITERATIONS'/                        &
     &       ' ******************************************************')

 2040 FORMAT('     NO SOLUTION COULD BE FOUND IN 100 ITERATIONS')
 2050 FORMAT(                                                                   &
    ' ITERATION PROCEDURE HAD TO BE RESTARTED TO AVOID NEGATIVE TEMPERATURE'/   &
    ' RESTART PROCEDURE WAS ABORTED AFTER 1000 TOTAL NUMBER OF ITERATIONS')
 2060 FORMAT(                                            &
    ' THE MAXIMUM MASS FLOW FOR WHICH A SOLUTION COULD BE OBTAINED WAS',G16.7/  &
    ' THE MINIMUM MASS FLOW FOR WHICH A SOLUTION COULD BE OBTAINED WAS',G16.7/     &
    ' THE MAXIMUM VALUE OF VSUBM AT THE HUB FOR WHICH A SOLUTION COULD BE OBTAINED WAS',G16.7/   &
    ' THE MINIMUM VALUE OF VSUBM AT THE HUB FOR WHICH A SOLUTION COULD BE OBTAINED WAS',G16.7/   &
    ' THE TOTAL NUMBER OF ITERATIONS WAS',I7)
 2070 FORMAT(' NSUB =',I7)
 2080 FORMAT(' NADD =',I7)
 2110 FORMAT(' ******************************************************'/  &
     & ' THE LIMIT OF 100 STATIONS PER CASE HAS BEEN EXCEEDED'/          &
     & ' OUTPUT IS GIVEN ONLY FOR THE FIRST 100 STATIONS'/               &
     &       ' ******************************************************')
 2120 FORMAT(6X,'DISTANCE',12X,'V',11X,'V/VCR',11X,'V-SUB-M',11X,   &
     &   'BETA',7X,'STATIC PRESSURE'/(2X,6G16.6))
 2130 FORMAT(' HUB OUTPUT DATA FOR ',I4,' STATIONS')
 2140 FORMAT(' SHROUD OUTPUT DATA FOR',I4,' STATIONS')



END Program AnnularDuct   ! ====================================================
