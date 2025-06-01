!+
! PROGRAM AIRSMO  (Airfoil Smoother)
! ------------------------------------------------------------------------------
! PURPOSE - Smooth a table af airfoil coordinates in preparation for
!   submission to a aerodynamic prediction code.

! AUTHORS - Harry L. Morgan, NASA Langley Research Center
!           Michael Selig, U. Illinois
!           Ralph L. Carmichael, Public Domain Aeronautical Software

!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1986   1.0   HLM   Original coding from NASA TM 84666
!   2000   1.1   MSS   Used by Michael Selig at U. Illinois with minor mods
!2021Jun10 1.2   RLC   PDAS version created from afsmo_2015.f from Selig
!2021Jun13 1.21  RLC   Subroutine array dimension of (1) -> (*)
!2021Jun13 1.22  RLC   Added padding to make common blocks uniform
!2021Jun14 1.23  RLC   Replaced constant -1 in call to CSDS with variable
!2021Jun24 1.24  RLC   Made new plot files JGNU1,JGNU2,JGNU3
!2021Jun27 1.25  RLC   write sqrt(curvature) to plot file
!2021Jul13 1.26  RLC   removed page eject from 1st and last pages
!2021Sep30 1.3   RLC   Added blank line in three places to help post processing




! REF: NASA TM 84666 by Harry Morgan


CC*DECK AIRSMO
      PROGRAM AIRSMO
C 
C     THIS PROGRAM PRESENTS A TECHNIQUE FOR SMOOTHING AIRFOIL 
C     COORDINATES USING LEAST SQUARES POLYNOMIAL AND CUBIC SPLINE 
C     METHODS. THIS IS FTN5 VERSION (JULY 1986).
C 
C        CODED BY -- HARRY MORGAN  NASA/LARC/AAD/SAB       1992 
C
C     Calls INPUT,SMOXY,PCARD,CAMTK,INTP
      IMPLICIT REAL*8(A-H,O-Z)
C 
      CHARACTER IFILE*30, TITLE(8)*10
C
      DIMENSION XINT(100), X(200), Y(200), W(200), YSMO(200), 
     1YPS(200), YPPS(202), THETA(202) 
C 
      COMMON /HLM/ DUMMX(2000)
C 
      COMMON /SMY/ DUMMY(2130)
C 
      COMMON /BLK1/ PI,PI2,RAD,CONS 
C 
      COMMON /INOUT/ JREAD,JWRITE,IPRINT,JGNU1,JGNU2,JGNU3     ! RLC  20210625
      
      JREAD=3 
      JWRITE=4
      IPRINT=1
      JGNU1=11     ! RLC  20210625
      JGNU2=12      ! RLC  20210625
      JGNU3=13      ! RLC 20210626
C 
C     OPEN FILES
C
      OPEN (5)
      OPEN (6)
      WRITE (6,500)
500   FORMAT (/'ENTER NAME OF INPUT FILE ===>'$)
      READ (5,'(A30)') IFILE
      OPEN (UNIT=JREAD,FILE=IFILE,STATUS='OLD',ACTION='READ')
      
      
      
      WRITE (6,501)
501   FORMAT (/'ENTER NAME OF OUTPUT FILE ===>'$)
      READ (5,'(A30)') IFILE
      OPEN (4,FILE=IFILE,STATUS='REPLACE',ACTION='WRITE')      
      WRITE (6,502)
502   FORMAT (/'ENTER NAME OF PUNCH FILE ===>'$)
      READ (5,'(A30)') IFILE
      OPEN (1,FILE=IFILE,STATUS='REPLACE',ACTION='WRITE')

!!!      OPEN(UNIT=4, FILE='afsmo.out', STATUS='REPLACE',ACTION='WRITE')
!!!      OPEN(UNIT=1, FILE='afsmo.pch', STATUS='REPLACE')      
      OPEN(UNIT=JGNU1, FILE='afsmo1.gnu', 
     X   STATUS='REPLACE',ACTION='WRITE')     ! RLC  20210625
      OPEN(UNIT=JGNU2, FILE='afsmo2.gnu', 
     X   STATUS='REPLACE',ACTION='WRITE')     ! RLC  20210625   
      OPEN(UNIT=JGNU3, FILE='afsmo3.gnu', 
     X   STATUS='REPLACE',ACTION='WRITE')     ! RLC  20210626   

C
C     INITIALIZE PROGRAM CONSTANTS
C 
      PI=DACOS(-1.D00)
      PI2=PI/2. 
      RAD=180./PI 
      CONS=1./(1.+DATAN(DSINH(PI2)))
 !!!     write(*,*) cons
      EPS=1.D-6 
      DF=1.D-4
      REWIND 1
C 
C     READ INPUT DATA 
C 
1     CALL INPUT (TITLE,ITER,IPLOT,IPUNCH,IOP,ICAMTK,INTR,YLTE,YNOSE,YUT
     1E,NINT,XINT,CNEW,NP,X,Y,W,THETA,YPS,YPPS,NOSE,CHORD,IERR) 
      WRITE(*,*) 'Input completed with ierr= ', ierr
      IF (IERR-1) 2,1,5 
C 
C     SMOOTH AIRFOIL COORDINATES
C 
2     CALL SMOXY (THETA,X,Y,W,YSMO,YPS,YPPS,NP,NOSE,YLTE,YNOSE,YUTE,EPS,
     1 DF,ITER,TITLE,IOP,IERR) 
      WRITE(*,*) 'Smoothing completed with ierr= ',ierr
      IF (IERR.NE.0) GO TO 1
C 
C     PUNCH OUTPUT DATA 
C 
      IF (IPUNCH.GE.1.AND.IPUNCH.LE.4) THEN
         CALL PCARD (IPUNCH,X,Y,W,THETA,YSMO,YPS,
     X    YPPS,NOSE,NP,CHORD,TITLE)
         WRITE(*,*) 'Punch out completed'
      END IF
      
      
C     COMPUTE THICKNESS AND CAMBER DISTRIBUTION 
C 
      IF (ICAMTK.EQ.1) THEN
         CALL CAMTK (THETA,YSMO,YPPS,NOSE,NP,EPS,KPLOT,IPUNCH,TITLE)
         WRITE(*,*) 'Camber and thickness completed'
      END IF   

C     INTERPOLATE NEW COORDINATES 
C 
      IF (INTR.GT.0) THEN
         CALL INTP (THETA,X,YSMO,YPPS,NP,NOSE,CHORD,TITLE,
     X      NINT,XINT,CNEW,INTR,IPUNCH) 
         WRITE(*,*) 'New coordinates interpolated'
      END IF   

C     RETURN AND READ NEXT CASE 
C 
      GO TO 1 
C 
5     WRITE (JWRITE,6)
      ENDFILE 1
      REWIND 1
      ENDFILE JWRITE
      REWIND JWRITE
      WRITE(*,*) 'All cases completed'
      STOP
C 
6     FORMAT (////48X,"-- THE LAST CASE HAS BEEN PROCESSED --") 
      END PROGRAM AIRSMO
      
C*DECK INTER 
      SUBROUTINE INTER (XINT,YINT,N,X,Y,JSTART,JEND,ICD)
C 
C     INTERPOLATION ROUTINE 
C 
C        ROUTINE SOURCE -- NORTH AMERICAN ROCKWELL L. A. DIVISION  1973 
C 
C        ICD=0  WEIGHTING METHOD USED 
C        ICD=1  LINEAR INTERPOLATION
C 
C CALLED BY BADPT AND SMOXY;   NO CALLS
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(N), Y(N)
C 
C     CHECK TO SEE IF XINT IS OUTSIDE BOUNDS OF X-ARRAY 
C 
      JEND=JSTART 
      IF (JSTART.EQ.N) GO TO 12 
C        CHECK TO SEE IF X ARRAY IS INCREASING OR DECREASING
      SGN=1.
      IF (X(N).LT.X(JSTART)) SGN=-1.
      D1=SGN*(XINT-X(N))
      IF (D1.GE.0.0) GO TO 12 
      
      D1=SGN*(XINT-X(JSTART)) 
      IF (D1.LE.0.0) GO TO 13 
      
      IF (ICD.EQ.1) GO TO 14
      
C     WEIGHTING METHOD REQUIRES AT LEAST 4 VALUES IN X AND Y ARRAYS 
      IF (N.LT.4) GO TO 14
C 
C     WEIGHTING METHOD
C 
C        DETERMINE X-ARRAY INDICES FOR TWO POINTS FORWARD (J,L) AND TWO 
C        POINTS AFT (K,M) OF XINT 
      DO 1 L=JSTART,N 
      J=L 
      D1=SGN*(X(J)-XINT)
      IF (D1) 1,2,3 
1     JEND=J


2     YINT=Y(J) 
      RETURN
3     IF (J.LE.2) GO TO 5 
      IF (J.EQ.N) GO TO 4 
      JJ=3
      GO TO 6 
      
4     JJ=2
      J=N-1 
      GO TO 6 
      
5     JJ=1
      J=3 
      
6     K=J-1 
      M=J-2 
      L=J+1 
C        INTERPOLATE A YINT VALUE (YSL) BY FITTING A STRAIGHT LINE
C        BETWEEN K AND J
      D1=XINT-X(M)
      D2=XINT-X(K)
      D3=XINT-X(J)
      D=(XINT-X(K))/(X(J)-X(K)) 
      YSL=D*Y(J)+(1.0-D)*Y(K) 
C        INTERPOLATE A YINT VALUE (YP1) BY FITTING A QUADRATIC BETWEEN
C        M, K, AND J
      C1=D3*D2/((X(M)-X(K))*(X(M)-X(J)))
      C2=D1*D3/((X(K)-X(M))*(X(K)-X(J)))
      C3=D2*D1/((X(J)-X(M))*(X(J)-X(K)))
      YP1=C1*Y(M)+C2*Y(K)+C3*Y(J) 
C        INTERPOLATE A YINT VALUE (YP2) BY FITTING A QUADRATIC BETWEEN
C        K, J, AND L
      D4=XINT-X(L)
      C1=D4*D3/((X(K)-X(J))*(X(K)-X(L)))
      C2=D2*D4/((X(J)-X(K))*(X(J)-X(L)))
      C3=D3*D2/((X(L)-X(K))*(X(L)-X(J)))
      YP2=C1*Y(K)+C2*Y(J)+C3*Y(L) 
C 
      IF (JJ-2) 7,8,9 
7     YP2=YP1 
      D=(XINT-X(1))/(X(2)-X(1)) 
      YSL=D*Y(2)+(1.0-D)*Y(1) 
      GO TO 9 
      
8     YP1=YP2 
      D=(XINT-X(N-1))/(X(N)-X(N-1)) 
      YSL=D*Y(N)+(1.0-D)*Y(N-1) 
      
C        COMPUTE DEVIATION BETWEEN LINEAR AND QUADRATIC YINT VALUES 
9     DEV1=DABS(YP1-YSL) 
      DEV2=DABS(YP2-YSL) 
      IF (DEV1+DEV2) 10,10,11 
      
10    YINT=YSL
      RETURN
C        COMPUTE WEIGHTING FACTORS

11    WT2=(DEV1*D)/(DEV1*D+(1.0-D)*DEV2)
      WT1=1.0-WT2 
C        COMPUTE FINAL YINT 
      YINT=WT2*YP2+WT1*YP1
      RETURN
      
12    YINT=Y(N) 
      JEND=N
      RETURN
      
13    YINT=Y(JSTART)
      RETURN
C 
C     LINEAR INTERPOLATION METHOD 
C 
14    DO 15 L=JSTART,N
      J=L 
      D1=SGN*(X(J)-XINT)
      IF (D1) 15,2,16 
15    JEND=J

16    YINT=Y(J-1)+(Y(J)-Y(J-1))*(XINT-X(J-1))/(X(J)-X(J-1)) 
      RETURN
      
      END SUBROUTINE INTER
      
      
C*DECK INPUT 
      SUBROUTINE INPUT (TITLE,ITER,IPLOT,IPUNCH,IOP,ICAMTK,INTR,YLTE,YNO
     1SE,YUTE,NINT,XINT,CNEW,NP,X,Y,W,THETA,YPS,YPPS,NOSE,CHORD,IERR) 
C 
C     ROUTINE TO READ INPUT DATA FOR AIRFOIL SMOOTHING PROGRAM
C 
C CALLED BY MAIN PROGRAM    CALLS TRNSRT, BADPT(twice) 

C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C 
C***********************************************************************
C*                                                                     *
C*          DESCRIPTION OF INPUT CARDS FOR SMOOTHING PROGRAM           *
C*                                                                     *
C* CARD NUMBER                  DESCRIPTION                            *
C*                                                                     *
C*.....................................................................*
C*   1   FORMAT(8A10)                                                  *
C*            TITLE CARD                                               *
C*.....................................................................*
C*   2   FORMAT(8F10.0)                                                *
C*            ITER - MAXIMUM NUMBER OF SMOOTHING ITERATIONS            *
C*            IPLOT - PLOTTING OPTION                                  *
C*                    0 - NO PLOTS                                     *
C*                    1 - PLOT SMOOTHED AND UNSMOOTHED Y/C, SMOOTHED   *
C*                        YPS, AND SMOOTHED YPPS VS THETA              *
C*                    2 - PLOT SMOOTHED AND UNSMOOTHED Y/C VS X/C      *
C*                    3 - PLOT SMOOTHED CURVATURE VS THETA             *
C*                    4 - PLOT CAMBER AND THICKNESS DISTRIBUTION       *
C*                    5 - PLOT OPTIONS 1 AND 2                         *
C*                    6 - PLOT OPTIONS 1 AND 3                         *
C*                    7 - PLOT OPTIONS 1, 2, AND 3                     *
C*                    8 - PLOT OPTIONS 1 AND 4                         *
C*                    9 - PLOT OPTIONS 1, 2, AND 4                     *
C*                   10 - PLOT OPTIONS 1, 2, 3, AND 4                  *
C*            IPUNCH - PUNCH OUTPUT OPTION                             *
C*                     0 - NO PUNCHED OUTPUT                           *
C*                     1 - SMOOTHED (X,Y,W) PUNCHED                    *
C*                     2 - SMOOTHED (THETA,Y/C,W) PUNCHED              *
C*                     3 - SMOOTHED (THETA,YPS,W) PUNCHED (YLTE,       *
C*                         YNOSE, AND YUTE ALSO PUNCHED)               *
C*                     4 - SMOOTHED (THETA,YPPS,W) PUNCHED (YLTE,      *
C*                         YNOSE, AND YUTE ALSO PUNCHED)               *
C*                     5 - THICKNESS AND CAMBER DISTRIBUTION (X/C,     *
C*                         Y/C, T/C/2, AND SLOPE) PUNCHED              *
C*                     6 - INTERPOLATED COORDINATES PUNCHED            *
C*            IOP - INPUT DATA OPTION                                  *
C*                  0 - (X,Y,W) INPUT                                  *
C*                  1 - (THETA,Y/C,W) INPUT                            *
C*                  2 - (THETA,YPS,W) INPUT                            *
C*                  3 - (THETA,YPPS,W) INPUT                           *
C*            ICAMTK - THICKNESS AND CAMBER DISTRIBUTION OPTION        *
C*                     0 - DO NOT COMPUTE THICKNESS AND CAMBER         *
C*                     1 - COMPUTE THICKNESS AND CAMBER                *
C*            IBAD - BAD COORDINATE CHECK OPTION                       *
C*                   0 - DO NOT CHECK FOR BAD COORDINATES              *
C*                   1 - CHECK FOR BAD COORDINATES                     *
C*            ITRN - INPUT COORDINATE TRANSLATION AND ROTATION OPTION  *
C*                   0 - DO NOT TRANSLATE AND ROTATE                   *
C*                   1 - TRANSLATE AND ROTATE SO THAT X-AXIS           *
C*                       CORRESPONDS TO THE LONGEST CHORDLINE          *
C*            INTR - COORDINATE INTERPOLATION OPTION                   *
C*                   0 - NO INTERPOLATION DESIRED                      *
C*                   1 - INTERPOLATE NEW COORDINATES USING STANDARD 57 *
C*                       X/C COORDINATES DEFINED IN SUBROUTINE INTP    *
C*                   2 - INTERPOLATE NEW COORDINATES AT INPUT X/C      *
C*                       VALUES  (0.0 .GE. X/C .LE. 1.0)               *
C*.....................................................................*
C* mss033195 took out formatted read
C*   3   FORMAT(10.0)                                                  *
C*            NU - NUMBER OF UPPER SURFACE INPUT COORDINATES           *
C*.....................................................................*
C* mss033195 took out formatted read, WL, WU disabled in process
C*   4   FORMAT(3F10.3)                                                *
C*            XU,YU,WU - UPPER SURFACE INPUT COORDINATES AND WEIGHTING *
C*                       (NU CARDS ARE INPUT)                          *
C*                    IF IOP=0, XU=X AND YU=Y COORDINATES              *
C*                    IF IOP=1, XU=THETA AND YU=Y/C                    *
C*                    IF IOP=2, XU=THETA AND YU=YPS                    *
C*                    IF IOP=3, XU=THETA AND YU=YPPS                   *
C*                    FOR ALL IOP, WU=WEIGHTING FACTOR                 *
C*.....................................................................*
C* mss033195 took out formatted read
C*   5   FORMAT(10.0)                                                  *
C*            NL - NUMBER OF LOWER SURFACE INPUT COORDINATES           *
C*.....................................................................*
C* mss033195 took out formatted read, WL, WU disabled in process
C*   6   FORMAT(3F10.3)                                                *
C*            XL,YL,WL - LOWER SURFACE INPUT COORDINATES AND WEIGHTING *
C*                       (NL CARDS ARE INPUT)                          *
C*                    IF IOP=0, XL=X AND YL=Y COORDINATES              *
C*                    IF IOP=1, XL=THETA AND YL=Y/C                    *
C*                    IF IOP=2, XL=THETA AND YL=YPS                    *
C*                    IF IOP=3, XL=THETA AND YL=YPPS                   *
C*                    FOR ALL IOP, WL=WEIGHTING FACTOR                 *
C*.....................................................................*
C*   7   FORMAT(3F10.0)      SKIP IF IOP=0 OR 1                        *
C*            YLTE,YNOSE,YUTE - LOWER SURFACE TRAILING-EDGE, NOSE,     *
C*                              AND UPPER SURFACE TRAILING-EDGE        *
C*                              Y/C COORDINATES                        *
C*.....................................................................*
C* mss 033195 took out formatted read
C*   8   FORMAT(F10.0)       SKIP IF INTR=0 OR 1                       *
C*            NINT - NUMBER OF INTERPOLATION X/C COORDINATES           *
C*.....................................................................*
C*  mss 033195
C*   9   FORMAT(*)      SKIP IF INTR=0 OR 1                            *
C*            XINT - INTERPOLATION X/C COORDINATES (NINT VALUES INPUT) *
C*.....................................................................*
C* mss 033195
C*  10   FORMAT(*)           SKIP IF INTR=0 OR 1                       *
C*            CNEW - DESIRED CHORD LENGTH OF INTERPOLATED COORDINATES  *
C*.....................................................................*
C*                                                                     *
C*       RESTRICTIONS:                                                 *
C*             ITER NOT GREATER THAN 300                               *
C*             NU OR NL NOT GREATER THAN 100                           *
C*             NINT NOT GREATER THAN 100                               *
C*                                                                     *
C***********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C 
      CHARACTER TITLE(8)*10
C
      DIMENSION VAR(8), XINT(*), X(*), Y(*), W(*), THETA(*), YPS(*), 
     1 YPPS(*)
C 
      COMMON /SMY/ XU(100),YU(100),WU(100),XL(100),YL(100),WL(100),
     X DUMMY2(1530) 
C 
      COMMON /BLK1/ PI,PI2,RAD,CONS 
C 
      COMMON /INOUT/ JREAD,JWRITE,IPRINT,JGNU1,JGNU2,JGNU3     ! RLC  20210625
C 
C     INITIALIZE ROUTINE CONSTANTS
C 
      ITRMAX=300
      NMAX=100
      TOLR=1.D-2
      IERR=0
C 
C     READ AND PRINT INPUT DATA 
C 
C        READ AND WRITE TITLE 
      READ (JREAD,27,END=25) TITLE
C 
      WRITE (JWRITE,28) TITLE 
C        READ AND WRITE OPTIONS 
      READ (JREAD,29) VAR 
      ITER   = IDINT(VAR(1))         ! IDINT converts REAL*8 to INTEGER
      IPLOT  = IDINT(VAR(2))
      IPUNCH = IDINT(VAR(3)) 
      IOP    = IDINT(VAR(4))
      ICAMTK = IDINT(VAR(5)) 
      IBAD   = IDINT(VAR(6)) 
      ITRN   = IDINT(VAR(7)) 
      INTR   = IDINT(VAR(8)) 
C             CHECK LIMITS OF OPTIONS 
      IF (ITER.GT.ITRMAX) ITER=ITRMAX 
C      IF (IPLOT.GT.10) IPLOT=0
      IPLOT = 0
C      IF (IPUNCH.GT.6) IPUNCH=0 
      IPUNCH = 1
      IF (IOP.GT.3) GO TO 23
      IF (ICAMTK.NE.0) ICAMTK=1 
      IF (IBAD.NE.0) IBAD=1 
      IF (ITRN.NE.0) ITRN=1 
      IF (INTR.GT.2) INTR=0 
      WRITE (JWRITE,30) ITER,IPLOT,IPUNCH,IOP,ICAMTK,IBAD,ITRN,INTR 
C        READ AND WRITE NUMBER OF UPPER SURFACE INPUT POINTS
ccc      READ (JREAD,29) VAR(1)  mss033195
      READ(JREAD,*) VAR(1)
      NU=IDINT(VAR(1)) 
      IF (NU.GT.NMAX) GO TO 22
      WRITE (JWRITE,31) NU
C        READ AND WRITE UPPER SURFACE INPUT POINTS AND WEIGHTING
ccc      READ (JREAD,32) (XU(I),YU(I),WU(I),I=1,NU) mss033195
      READ(JREAD,*)  (XU(I),YU(I),I=1,NU)
      DO 2 I=1,NU 
         WU(I) = 0.
         IF (WU(I).LT.1.0) WU(I)=1.0 
2     CONTINUE
      IF (IOP.EQ.0) WRITE (JWRITE,33) (XU(I),I=1,NU)
      IF (IOP.NE.0) WRITE (JWRITE,34) (XU(I),I=1,NU)
      IF (IOP.LT.2) WRITE (JWRITE,35) (YU(I),I=1,NU)
      IF (IOP.EQ.2) WRITE (JWRITE,36) (YU(I),I=1,NU)
      IF (IOP.EQ.3) WRITE (JWRITE,37) (YU(I),I=1,NU)
      WRITE (JWRITE,38) (WU(I),I=1,NU)
C        READ AND WRITE NUMBER OF LOWER SURFACE INPUT POINTS
ccc   READ (JREAD,29) VAR(1) mss033195
      READ(JREAD,*) VAR(1)
      NL=IDINT(VAR(1)) 
      IF (NL.GT.NMAX) GO TO 22
      WRITE (JWRITE,39) NL
C        READ AND WRITE LOWER SURFACE INPUT POINTS AND WEIGHTING
ccc   READ (JREAD,32) (XL(I),YL(I),WL(I),I=1,NL) mss033195
      READ (JREAD,*) (XL(I),YL(I),I=1,NL)
      DO 3 I=1,NL 
         WL(I) = 0.
         IF (WL(I).LT.1.0) WL(I)=1.0 
3     CONTINUE
      IF (IOP.EQ.0) WRITE (JWRITE,40) (XL(I),I=1,NL)
      IF (IOP.NE.0) WRITE (JWRITE,41) (XL(I),I=1,NL)
      IF (IOP.LT.2) WRITE (JWRITE,42) (YL(I),I=1,NL)
      IF (IOP.EQ.2) WRITE (JWRITE,43) (YL(I),I=1,NL)
      IF (IOP.EQ.3) WRITE (JWRITE,44) (YL(I),I=1,NL)
      WRITE (JWRITE,45) (WL(I),I=1,NL)
C        READ AND WRITE TRAILING-EDGE COORDINATES 
      IF (IOP.LE.1) GO TO 4 
      READ (JREAD,29) YLTE,YNOSE,YUTE 
      WRITE (JWRITE,46) YLTE,YNOSE,YUTE 
C        READ AND WRITE NUMBER OF INTERPOLATION COORDINATES 
4     IF (INTR.EQ.0) GO TO 6
      IF (INTR.NE.2) GO TO 5
ccc      READ (JREAD,29) VAR(1) mss033195
      READ(JREAD,*) VAR(1)
      NINT=IDINT(VAR(1)) 
      IF (NINT.GT.NMAX) GO TO 24
      WRITE (JWRITE,47) NINT
      WRITE(6,*) ' NINT =', NINT
C        READ AND WRITE INTERPOLATION COORDINATES 
C mss033195 made unformatted read
      READ (JREAD,*) (XINT(I),I=1,NINT)
      WRITE (JWRITE,48) (XINT(I),I=1,NINT)
C        READ AND WRITE NEW CHORD OF INTERPOLATED COORDINATES 
C mss 033195 made unformatted read
5     IF (INTR.EQ.2) READ (JREAD,*) CNEW 
      IF (INTR.EQ.2) WRITE (JWRITE,49) CNEW 
C 
C     CHECK UPPER SURFACE COORDINATES FOR BAD POINTS
C 
6     IF (IOP.NE.0) GO TO 7 
      IF (IBAD.EQ.1) CALL BADPT (XU,YU,NU,TOLR,1,IERR)
      IF (IERR.NE.0) GO TO 26 
C 
C     CHECK LOWER SURFACE COORDINATES FOR BAD POINTS
C 
      IF (IBAD.EQ.1) CALL BADPT (XL,YL,NL,TOLR,2,IERR)
      IF (IERR.NE.0) GO TO 26 
C 
C     TRANSLATE AND ROTATE THE INPUT COORDINATES SO THAT THE X-AXIS 
C     CORRESPONDS TO THE LONGEST CHORDLINE OF THE AIRFOIL 
C 
      IF (ITRN.EQ.1) CALL TRNSRT (XU,YU,WU,NU,XL,YL,WL,NL,TITLE)
C 
C     LOAD X, Y, THETA, YPS, AND YPPS ARRAYS
C 
7     IF (IOP) 8,8,15 
C        IF IOP=0, COMPUTE THETA FROM INPUT X 
C             COMPUTE THETA FOR LOWER SURFACE 
8     CHORD=XL(NL)-XL(1)
      DELTA=XU(NU)-XU(1)
      IF (DELTA.GT.CHORD) CHORD=DELTA 
      NP=0
      DO 11 I=1,NL
      NP=NP+1 
      J=NL+1-I
      W(NP)=WL(J) 
      DELTA=(XL(J)-XL(1))/CHORD 
      IF (DELTA.LE.CONS) GO TO 9
      DELTA=DTAN(DELTA/CONS-1.)
      THETA(NP)=-PI2-DLOG(DELTA+DSQRT(DELTA*DELTA+1.)) 
      GO TO 10
      
9     THETA(NP)=-DACOS(1.-DELTA/CONS)

10    X(NP)=XL(J)/CHORD 
11    Y(NP)=YL(J)/CHORD 
      NOSE=NP 
C             COMPUTE THETA FOR UPPER SURFACE 
      J=1 
      IF (XL(1).EQ.XU(1).AND.YL(1).EQ.YU(1)) J=2
      DO 14 I=J,NU
      NP=NP+1 
      W(NP)=WU(I) 
      DELTA=(XU(I)-XU(1))/CHORD 
      IF (DELTA.LE.CONS) GO TO 12 
      DELTA=DTAN(DELTA/CONS-1.)
      THETA(NP)=PI2+DLOG(DELTA+DSQRT(DELTA*DELTA+1.))
      GO TO 13
12    THETA(NP)=DACOS(1.-DELTA/CONS) 
13    X(NP)=XU(I)/CHORD 
14    Y(NP)=YU(I)/CHORD 
      GO TO 20
C        IF IOP=1, 2, OR 3, COMPUTE X/C FROM INPUT THETA
C             COMPUTE X/C FOR LOWER SURFACE 
15    CHORD=1.0 
      NP=0
      DO 17 I=1,NL
      NP=NP+1 
      J=NL+1-I
      W(NP)=WL(J) 
      IF (IOP.EQ.1) Y(NP)=YL(J) 
      IF (IOP.EQ.2) YPS(NP)=YL(J) 
      IF (IOP.EQ.3) YPPS(NP)=YL(J)
      THETA(NP)=XL(J)/RAD 
      DELTA=DABS(THETA(NP))
      IF (DELTA.GT.PI2) GO TO 16
      XL(J)=CONS*(1.-DCOS(DELTA))
      GO TO 17
16    XL(J)=CONS*(DATAN(DSINH(DELTA-PI2))+1.) 
17    X(NP)=XL(J) 


      NOSE=NP 
C             COMPUTE X/C FOR UPPER SURFACE 
      XU(1)=XL(1) 
      DO 19 I=2,NU
      NP=NP+1 
      W(NP)=WU(I) 
      IF (IOP.EQ.1) Y(NP)=YU(I) 
      IF (IOP.EQ.2) YPS(NP)=YU(I) 
      IF (IOP.EQ.3) YPPS(NP)=YU(I)
      THETA(NP)=XU(I)/RAD 
      DELTA=DABS(THETA(NP))
      IF (DELTA.GT.PI2) GO TO 18
      XU(I)=CONS*(1.-DCOS(DELTA))
      GO TO 19
18    XU(I)=CONS*(DATAN(DSINH(DELTA-PI2))+1.) 
19    X(NP)=XU(I) 
C 
C     PRINT SUMMARY OF INPUT DATA 
C 
      WRITE(*,*) 'Printing summary of input data, ip,ierr=', np, ierr
20    WRITE (JWRITE,50) TITLE 
      DO 21 I=1,NP
      DELTA=THETA(I)*RAD
      
      IF (IOP.LE.1) THEN
         WRITE (JWRITE,51) I,X(I),Y(I),DELTA,W(I)
         WRITE (JGNU1,51)  I,X(I),Y(I),DELTA,W(I)   ! RLC   2021-06
      END IF   
      
      IF (IOP.EQ.2) WRITE (JWRITE,52) I,X(I),DELTA,YPS(I),W(I)
      IF (IOP.EQ.3) WRITE (JWRITE,53) I,X(I),DELTA,YPPS(I),W(I) 
21    CONTINUE
      WRITE (JWRITE,54) CHORD 
      GO TO 26
C 
C     PRINT ERROR MESSAGES
C 
22    NN=IDINT(VAR(1)) 
      WRITE (JWRITE,55) NN
      GO TO 25
23    WRITE (JWRITE,56) IOP 
      GO TO 25
24    WRITE (JWRITE,57) NINT
C 
C     NO ADDITIONAL INPUT DATA
C 
25    IERR=2
C 
C     RETURN TO CALLING PROGRAM 
C 
26    RETURN
C 
27    FORMAT (8A10) 
28    FORMAT (57X,14H--INPUT DATA--//5X,7HTITLE--,2X,8A10)   ! RLC 20210713
29    FORMAT (8F10.5) 
30    FORMAT (/5X,6HITER =,I4,3X,7HIPLOT =,I3,3X,8HIPUNCH =,I3,3X,5HIOP 
     1=,I3,3X,8HICAMTK =,I3,3X,6HIBAD =,I3,3X,6HITRN =,I3,3X,6HINTR =,I3
     2) 
31    FORMAT (/5X,4HNU =,I4)
32    FORMAT (3F10.5) 
33    FORMAT (/5X,3HXU=,8E15.6/(8X,8E15.6)) 
34    FORMAT (/5X,3HTU=,8E15.6/(8X,8E15.6)) 
35    FORMAT (/5X,3HYU=,8E15.6/(8X,8E15.6)) 
36    FORMAT (/4X,4HYPU=,8E15.6/(8X,8E15.6))
37    FORMAT (/3X,5HYPPU=,8E15.6/(8X,8E15.6)) 
38    FORMAT (/5X,3HWU=,8E15.6/(8X,8E15.6)) 
39    FORMAT (/5X,4HNL =,I4)
40    FORMAT (/5X,3HXL=,8E15.6/(8X,8E15.6)) 
41    FORMAT (/5X,3HTL=,8E15.6/(8X,8E15.6)) 
42    FORMAT (/5X,3HYL=,8E15.6/(8X,8E15.6)) 
43    FORMAT (/4X,4HYPL=,8E15.6/(8X,8E15.6))
44    FORMAT (/3X,5HYPPL=,8E15.6/(8X,8E15.6)) 
45    FORMAT (/5X,3HWL=,8E15.6/(8X,8E15.6)) 
46    FORMAT (/3X,6HYLTE =,E15.6,5X,7HYNOSE =,E15.6,5X,6HYUTE =,E15.6)
47    FORMAT (/3X,6HNINT =,I4)
48    FORMAT (/3X,5HXINT=,8E15.6/(8X,8E15.6)) 
49    FORMAT (/3X,6HCNEW =,F10.3) 
50    FORMAT (1H1,29X,25H--SUMMARY OF INPUT DATA--//5X,9HTITLE--  ,8A10/
     1/9X,1HI,10X,3HX/C,12X,3HY/C,12X,5HTHETA,10X,3HYPS,12X,4HYPPS,14X,1
     2HW) 
51    FORMAT (I10,2F15.6,F15.2,30X,F15.2) 
52    FORMAT (I10,F15.6,15X,F15.2,F15.6,15X,F15.2)
53    FORMAT (I10,F15.6,15X,F15.2,15X,F15.6,F15.2)
54    FORMAT (/5X,7HCHORD =,F15.6)
55    FORMAT (//5X,28HINPUT CARD ERROR  NU OR NL =,I4)
56    FORMAT (//5X,23HINPUT CARD ERROR  IOP =,I4) 
57    FORMAT (//5X,24HINPUT CARD ERROR  NINT =,I5)
      END  SUBROUTINE INPUT
C*DECK TRNSRT
      SUBROUTINE TRNSRT (XU,YU,WU,NU,XL,YL,WL,NL,TITLE) 
C 
C     ROUTINE TO TRANSLATE AND ROTATE THE INPUT AIRFOIL COORDINATES SO
C     THAT THE X-AXIS CORRESPONDS TO THE LONGEST CHORDLINE
C 
C CALLED BY INPUT;   NO CALLS
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C
      IMPLICIT REAL*8(A-H,O-Z)
C 
      CHARACTER TITLE(8)*10
C
      DIMENSION XU(*), YU(*), WU(*), XL(*), YL(*), WL(*)   ! RLC  2021-06-13
C 
      COMMON /HLM/ X(200),Y(200),W(200), DUMMY6(1400)
C 
      COMMON /BLK1/ PI,PI2,RAD,CONS 
C 
      COMMON /INOUT/ JREAD,JWRITE,IPRINT,JGNU1,JGNU2,JGNU3     ! RLC  20210625
C 
C     PRINT INPUT COORDINATES 
C 
      WRITE (JWRITE,13) TITLE 
      J=NU
      IF (NL.GT.NU) J=NL
      DO 1 I=1,J
      IF (I.LE.NU.AND.I.LE.NL) WRITE (JWRITE,14) I,XU(I),YU(I),XL(I),YL(
     1I)
      IF (I.LE.NU.AND.I.GT.NL) WRITE (JWRITE,14) I,XU(I),YU(I)
      IF (I.GT.NU.AND.I.LE.NL) WRITE (JWRITE,15) I,XL(I),YL(I)
1     CONTINUE
C 
C     COMPUTE LONGEST CHORDLINE 
C 
C        LOAD LOWER SURFACE COORDINATES 
      N=0 
      DO 2 I=1,NL 
      J=NL+1-I
      N=N+1 
      W(N)=WL(J)
      X(N)=XL(J)
2     Y(N)=YL(J)
      J=1 
      IF (XL(1).EQ.XU(1).AND.YL(1).EQ.YU(1)) J=2
C        LOAD UPPER SURFACE COORDINATES 
      DO 3 I=J,NU 
      N=N+1 
      W(N)=WU(I)
      X(N)=XU(I)
3     Y(N)=YU(I)
C        COMPUTE MIDPOINT OF TRAILING-EDGE BASE 
      XTE=0.5*(X(1)+X(N)) 
      YTE=0.5*(Y(1)+Y(N)) 
C        FIND MOST FORWARD LEADING-EDGE POINT AND LONGEST CHORD 
      CHORD=0.0 
      DO 5 I=1,N
      DIST=DSQRT((X(I)-XTE)**2+(Y(I)-YTE)**2)
      IF (DIST-CHORD) 5,5,4 
4     CHORD=DIST
      NOSE=I
      XNOSE=X(I)
      YNOSE=Y(I)
5     CONTINUE
C 
C     TRANSLATE AND ROTATE AIRFOIL
C 
      IF (CHORD.LE.0.0) GO TO 6 
      COSA=(XTE-XNOSE)/CHORD
      SINA=(YTE-YNOSE)/CHORD
      ANGLE=DATAN(SINA/COSA)*RAD 
      GO TO 7 
6     COSA=0.0
      SINA=0.0
      ANGLE=0.0 
7     DO 8 I=1,N
      DIST=X(I) 
      X(I)=(DIST-XNOSE)*COSA+(Y(I)-YNOSE)*SINA
8     Y(I)=(Y(I)-YNOSE)*COSA-(DIST-XNOSE)*SINA
C 
C     REDEFINE LOWER AND UPPER SURFACE COORDINATES
C 
      DO 9 I=1,NOSE 
      J=NOSE+1-I
      WL(I)=W(J)
      XL(I)=X(J)
9     YL(I)=Y(J)
      NL=NOSE 
      DO 10 I=NOSE,N
      J=I+1-NOSE
      WU(J)=W(I)
      XU(J)=X(I)
10    YU(J)=Y(I)
      NU=J
C 
C     PRINT NEW AIRFOIL COORDINATES 
C 
      WRITE (JWRITE,16) TITLE 
      J=NU
      IF (NL.GT.NU) J=NL
      DO 11 I=1,J 
      IF (I.LE.NU.AND.I.LE.NL) WRITE (JWRITE,14) I,XU(I),YU(I),XL(I),YL(
     1I)
      IF (I.LE.NU.AND.I.GT.NL) WRITE (JWRITE,14) I,XU(I),YU(I)
      IF (I.GT.NU.AND.I.LE.NL) WRITE (JWRITE,15) I,XL(I),YL(I)
11    CONTINUE
      WRITE (JWRITE,12) XNOSE,YNOSE,ANGLE 
      RETURN
C 
C....leading slash on next line insures a blank line after writing the tables
12    FORMAT (/5X,7HXNOSE =,F15.6,5X,7HYNOSE =,F15.6,5X,7HANGLE =,F8.3) 
13    FORMAT (1H1,32X,21H--INPUT COORDINATES--//5X,7HTITLE--,2X,8A10//
     1 9X,1HI,11X,2HXU,13X,2HYU,13X,2HXL,13X,2HYL) 
14    FORMAT (5X,I5,4F15.6) 
15    FORMAT (5X,I5,30X,2F15.6) 
16    FORMAT (/1H1,21X,38H--TRANSLATED AND ROTATED COORDINATES--//5X,
     1 7HTITLE--,2X,8A10//9X,1HI,11X,2HXU,13X,2HYU,13X,2HXL,13X,2HYL)
      END  SUBROUTINE TRNSRT
C*DECK BADPT 
      SUBROUTINE BADPT (X,Y,NP,TOLR,ISURF,IERR) 
C 
C     ROUTINE TO EDIT BAD POINTS FROM X AND Y INPUT COORDINATES 

C CALLED BY INPUT;   CALLS INTER
C 
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C 
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(NP),Y(NP),SURF(2) 
C 
      COMMON /HLM/ TI(100),YI(100),YN(100),THETA(100),DUMMY7(1600)
C 
      COMMON /BLK1/ PI,PI2,RAD,CONS 
C 
      COMMON /INOUT/ JREAD,JWRITE,IPRINT,JGNU1,JGNU2,JGNU3     ! RLC  20210625
C 
      DATA SURF(1)/5HUPPER/,SURF(2)/5HLOWER/
C 
C     IF TOLERANCE IS ZERO OR NEGATIVE RETURN 
C 
      IERR=0
      IF (TOLR.LE.0.0) RETURN 
C 
C     COMPUTE LOCAL CHORD 
C 
      CHORD=X(NP)-X(1)
C 
C     INITIALIZE ITERATION PARAMETERS 
C 
      ICD=0 
      IPTP=0
      N1=NP-1 
      NMAX=0
      TOLC=TOLR*CHORD 
C 
C     COMPUTE THETA EQUIVALENT OF X 
C 
      DO 2 I=1,NP 
      DELTA=(X(I)-X(1))/CHORD 
      IF (DELTA.LE.CONS) GO TO 1
      DELTA=DTAN(DELTA/CONS-1.)
      THETA(I)=PI2+DLOG(DELTA+DSQRT(DELTA*DELTA+1.)) 
      GO TO 2 
1     THETA(I)=DACOS(1.-DELTA/CONS)
2     CONTINUE
C 
C     LOOP TO SEARCH FOR BAD POINTS 
C 
3     NMAX=NMAX+1 
      JSTART=1
C        COMPUTE NEW Y VALUE BY INTERPOLATION 
      DO 5 I=2,N1 
      K=0 
C             LOAD TI AND YI ARRAY - OMIT THE I(TH) INPUT DATA POINT
      DO 4 J=1,NP 
      IF (I.EQ.J) GO TO 4 
      K=K+1 
      TI(K)=THETA(J)
      YI(K)=Y(J)
4     CONTINUE
C             INTERPOLATE I(TH) DATA POINT
      CALL INTER (THETA(I),YN(I),K,TI,YI,JSTART,JEND,ICD) 
      JSTART=JEND 
5     CONTINUE
C        CHECK TOLERANCE OF INTERPOLATED POINTS 
      IPT=0 
      ERRMAX=0. 
      DO 7 I=2,N1 
      ERRMIN=0. 
      ERR=DABS(YN(I)-Y(I)) 
      IF (ERR.GE.TOLC) ERRMIN=ERR 
      IF (ERRMIN-ERRMAX) 7,7,6
6     IPT=I 
      ERRMAX=ERRMIN 
7     CONTINUE
      IF (IPT.EQ.0) RETURN
C        PRINT COORDINATES OF BAD POINTS
      IF (NMAX.EQ.1) WRITE (JWRITE,9) SURF(ISURF),TOLC
      WRITE (JWRITE,10) IPT,X(IPT),Y(IPT),YN(IPT) 
C        REPLACE BAD POINT WITH INTERPOLATED VALUE
      Y(IPT)=YN(IPT)
C        CHECK TO SEE IF THIS BAD POINT IS ADJACENT TO THE PREVIOUS BAD 
C        POINT -- IF IT IS, PRINT A WARNING MESSAGE AND TERMINATE 
C        PROGRAM EXECUTION
      IF ((IPTP.EQ.IPT-1).OR.(IPTP.EQ.IPT+1)) GO TO 8 
      IF (IPTP.EQ.IPT) GO TO 8
      IPTP=IPT
      IF (NMAX.GE.NP) RETURN
C 
C     RETURN TO START OF LOOP AND SEARCH FOR NEXT BAD POINT 
C 
      GO TO 3 
C 
C     WARNING MESSAGE PRINT STATEMENT 
C 
8     WRITE (JWRITE,11) 
      IERR=1
      RETURN
C 
9     FORMAT (1H1//1X,44HWARNING -- BAD POINTS HAVE BEEN FOUND ON THE,1X
     1,A5,1X,37HSURFACE BASED ON AN EDIT TOLERANCE OF,F10.6/) 
10    FORMAT (1X,15HBAD POINT AT I=,I4,5X,4HX = ,F10.6,5X,4HY = ,F10.6,5
     1X,18HREPLACED WITH Y = ,F10.6/) 
11    FORMAT (1X,93HADJACENT BAD POINTS HAVE BEEN FOUND -- PLEASE CORREC
     1T YOUR INPUT DATA AND RESUBMIT THIS CASE.)
      END  SUBROUTINE BADPT
      
      
C*DECK SMOXY 
      SUBROUTINE SMOXY (THETA,X,Y,W,YSMO,YPS,YPPS,NP,NOSE,YLTE,YNOSE,YUT
     1E,EPS,DF,ITER,TITLE,IOP,IERR) 
C 
C     THIS SUBROUTINE PRESENTS A TECHNIQUE FOR SMOOTHING Y INPUT
C     COORDINATES USING LEAST SQUARES POLYNOMIAL AND CUBIC SPLINE 
C     METHODS 
C 
C     IF IOP=0 OR 1, COMPUTE YPPU (UNSMOOTHED SECOND DERIVATIVES) FROM
C     LEAST SQUARES POLYNOMIAL FITTING OF Y VS THETA. THEN COMPUTE
C     YPPS (SMOOTHED SECOND DERIVATIVES) FROM LEAST SQUARES CUBIC 
C     SPLINE FITTING OF YPPU VS THETA. FINALLY COMPUTE YSMO (SMOOTHED Y 
C     COORDINATES) USING INVERSE CUBIC SPLINE METHOD. 
C 
C     IF IOP=2, COMPUTE SECOND DERIVATIVES FROM INPUT FIRST DERIVATIVES.
C     THEN COMPUTE UNSMOOTHED Y COORDINATES FROM SECOND DERIVATIVES AND 
C     FOLLOW SAME PROCEDURES AS OUTLINED ABOVE FOR IOP 0 OR 1.
C 
C     IF IOP=3, COMPUTE UNSMOOTHED Y COORDINATES FROM INPUT SECOND
C     DERIVATIVES. THEN FOLLOW SAME PROCEDURES AS OUTLINED ABOVE FOR
C     IOP 0 OR 1. 
C 

C CALLED BY MAIN PROGRAMCALLS INTER,CSDS,LSQSMO,YNEW
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C 
C        DIMENSION THETA, X, Y, W, YSMO, YPS, AND YPPS BY NP IN CALLING 
C        PROGRAM
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER TITLE(8)*10
C
      DIMENSION THETA(NP),X(NP),Y(NP),W(NP),YSMO(NP),YPS(NP),YPPS(NP) 
      
      DIMENSION WORKDUMMY(2000)   ! added by RLC  2021/06/14
C 
      COMMON /HLM/ WK(200,10)  
C 
      COMMON /SMY/ YPP(200),YUSMO(200),DUM(200),A(200,4),YN(200),
     X YPPU(200),SUMY(300),LTER(60) 
C 
      COMMON /BLK1/ PI,PI2,RAD,CONS 
C 
      COMMON /INOUT/ JREAD,JWRITE,IPRINT,JGNU1,JGNU2,JGNU3     ! RLC  20210625
C 
      DATA LMX/200/,WT/100./
C 

!!!      write(*,*) 'Entering SMOXY, np=', np

      IERR=0
      IF (IOP.EQ.0.OR.IOP.EQ.1) GO TO 13
      IF (IOP.EQ.2) GO TO 1 
      IF (IOP.EQ.3) GO TO 11
C 
C        IF IOP=2, COMPUTE SECOND DERIVATIVES FROM INPUT FIRST
C        DERIVATIVES. THEN COMPUTE INITIAL Y/C COORDINATES FROM SECOND
C        DERIVATIVES. 
C 
C             COMPUTE SECOND DERIVATIVES USING CSDS 
1     DO 2 I=1,NP 
2     DUM(I)=1.0
      T1=0.0
      CALL CSDS (LMX,NP,THETA,YPS,DUM,T1,-1,A,WK,IERR)
      IF (IERR.NE.0) GO TO 71 
      DO 4 I=1,NP 
      IF (I.EQ.NP) GO TO 3
      YPPS(I)=A(I,2)
      GO TO 4 
3     DELTA=THETA(I)-THETA(I-1) 
      YPPS(I)=(3.*A(I-1,4)*DELTA+2.*A(I-1,3))*DELTA+A(I-1,2)
4     CONTINUE
C             COMPUTE SECOND DERIVATIVES USING LSQSMO 
      DELTA=1.
      CALL LSQSMO (THETA,YPS,W,DUM,YPP,YUSMO,NP,1,NP,NOSE,DELTA,EPS,IERR
     1) 
      IF (IERR.NE.0) RETURN 
C             COMPUTE Y/C COORDINATES 
      CALL YNEW (THETA,YPPS,Y,NOSE,NP,YLTE,YNOSE,YUTE,EPS,DUM,WK,JWRITE,
     10)
      CALL YNEW (THETA,YPP,YUSMO,NOSE,NP,YLTE,YNOSE,YUTE,EPS,DUM,WK,JWRI
     1TE,0) 
C             COMPUTE NEW FIRST DERIVATIVES AND COMPARE WITH INPUT
C             FIRST DERIVATIVES 
      WRITE (JWRITE,73) TITLE 
      SUM1=0.0
      SUM2=0.0
      DO 7 I=1,NP 
      IF (I.EQ.1) GO TO 5 
      DELTA=THETA(I)-THETA(I-1) 
      YN(I)=YPPS(I-1)*DELTA/6.+YPPS(I)*DELTA/3.+(Y(I)-Y(I-1))/DELTA 
      DUM(I)=YPP(I-1)*DELTA/6.+YPP(I)*DELTA/3.+(YUSMO(I)-YUSMO(I-1))/DEL
     1TA
      GO TO 6 
5     DELTA=THETA(2)-THETA(1) 
      YN(1)=-YPPS(1)*DELTA/3.-YPPS(2)*DELTA/6.+(Y(2)-Y(1))/DELTA
      DUM(1)=-YPP(1)*DELTA/3.-YPP(2)*DELTA/6.+(YUSMO(2)-YUSMO(1))/DELTA 
6     T1=YPS(I)-YN(I) 
      T2=YPS(I)-DUM(I)
      SUM1=SUM1+T1*T1 
      SUM2=SUM2+T2*T2 
7     WRITE (JWRITE,74) I,YPS(I),YN(I),T1,DUM(I),T2 


      WRITE (JWRITE,75) SUM1,SUM2 
C             SELECT OUTPUT FROM EITHER CSDS OR LSQSMO
      DO 10 I=1,NP
      IF (SUM2.LT.SUM1) GO TO 8 
      YPP(I)=YPPS(I)
      GO TO 9 
8     Y(I)=YUSMO(I) 
      YN(I)=DUM(I)
9     YSMO(I)=Y(I)
10    YUSMO(I)=Y(I) 
      IF (SUM2.GE.SUM1) WRITE (JWRITE,76) 
      IF (SUM2.LT.SUM1) WRITE (JWRITE,77) 
      IF (ITER.EQ.0) GO TO 48 
      GO TO 13
C 
C        IF IOP=3, COMPUTE INITIAL Y/C FROM INPUT SECOND DERIVATIVES
C        AND Y/C COORDINATES AT THE UPPER AND LOWER SURFACE TRAILING
C        EDGE AND NOSE
C 
11    CALL YNEW (THETA,YPPS,Y,NOSE,NP,YLTE,YNOSE,YUTE,EPS,DUM,WK,JWRITE,
     10)
C             COMPUTE FIRST DERIVATIVES 
      DO 12 I=1,NP
      YSMO(I)=Y(I)
      YUSMO(I)=Y(I) 
      IF (I.EQ.1) GO TO 12
      DELTA=THETA(I)-THETA(I-1) 
      YN(I)=YPPS(I-1)*DELTA/6.+YPPS(I)*DELTA/3.+(Y(I)-Y(I-1))/DELTA 
12    YPP(I)=YPPS(I)
      DELTA=THETA(2)-THETA(1) 
      YN(1)=-YPPS(1)*DELTA/3.-YPPS(2)*DELTA/6.+(Y(2)-Y(1))/DELTA
      IF (ITER.EQ.0) GO TO 48 
C 
C        INITIALIZE ARRAYS
C 
13    DO 14 I=1,NP
      YUSMO(I)=Y(I) 
      IF (IOP.LT.2) YPP(I)=0.0
      YSMO(I)=THETA(I)*RAD
14    DUM(I)=1. 
      IF (ITER.GT.0) GO TO 17 
C 
C        IF IOP=0 OR 1 AND NO SMOOTHING DESIRED (I.E. ITER=0) , COMPUTE 
C        SECOND DERIVATIVE FROM CUBIC SPLINE SUBROUTINE 
C 
      IPTXXX=-1
      CALL CSDS (LMX,NP,THETA,Y,DUM,0.0D00,IPTXXX,A,WK,IERR)    ! RLC  2021-06-15
      IF (IERR.NE.0) GO TO 71 
C        COMPUTE Y AND SECOND DERIVATIVE
      DO 16 I=1,NP
      IF (I.EQ.NP) GO TO 15 
      YSMO(I)=A(I,1)
      YN(I)=A(I,2)
      YPP(I)=2.*A(I,3)
      GO TO 16
15    DELTA=THETA(I)-THETA(I-1) 
      YSMO(I)=((A(I-1,4)*DELTA+A(I-1,3))*DELTA+A(I-1,2))*DELTA+A(I-1,1) 
      YN(I)=(3.*A(I-1,4)*DELTA+2.*A(I-1,3))*DELTA+A(I-1,2)
      YPP(I)=6.*A(I-1,4)*DELTA+2.*A(I-1,3)
16    CONTINUE


      GO TO 48
C 
C        FIND MAXIMUM INPUT Y VALUE AND ITS LOCATION FOR UPPER AND
C        LOWER SURFACES 
C             LOWER SURFACE 
17    YMAX=0.0
      JMAXL=1 
      DO 19 I=1,NOSE
      J=NOSE+1-I
      IF (DABS(Y(J)).GT.YMAX) GO TO 18 
      GO TO 19
18    YMAX=DABS(Y(J))
      JMAXL=J 
19    CONTINUE
C             UPPER SURFACE 
      YMAX=0.0
      JMAXU=1 
      DO 21 I=NOSE,NP 
      IF (DABS(Y(I)).GT.YMAX) GO TO 20 
      GO TO 21
20    YMAX=DABS(Y(I))
      JMAXU=I 
21    CONTINUE
C 
C        COMPUTE UNSMOOTHED SECOND DERIVATIVE USING LEAST 
C        SQUARES POLYNOMIAL METHOD
C 
      J1=0
      ICON=0
      MTER=0
      J=ITER
      KTI=0 
      IF (IPRINT.NE.0) WRITE (JWRITE,78) TITLE
      
      
      DO 23 I=1,30
      KTI=KTI+1 
      LTER(I)=10
      J=J-10
      IF (J) 22,24,23 
22    LTER(I)=10+J
      GO TO 24
23    CONTINUE

24    CONTINUE

!!!      write(*,*) 'Iteration Loop,  KTI,ITER,J=', KTI,ITER,J
!!!      write(*,*) LTER(1:KTI)
      
      DO 39 LL=1,KTI
      N1=LTER(LL) 
      DO 34 I=1,N1
C             CALL LEAST SQUARES POLYNOMIAL SMOOTHING ROUTINE 
      CALL LSQSMO (THETA,YUSMO,W,YN,DUM,YPPU,NP,JMAXL,JMAXU,NOSE,WT,EPS,
     1IERR) 
      IF (IERR.NE.0) RETURN 
      J = 10*(LL-1) + I
      WRITE (6,500) J  
500   FORMAT ('  COMPLETED ITERATION ',I4)
C             COMPUTE ERROR TERM
      SUMY(I)=0.0 
      DO 25 J=1,NP
25    SUMY(I)=SUMY(I)+(YPPU(J)-YPP(J))**2 


      J1=J1+1 
      IF ((I.LE.3).AND.(LL.EQ.1)) GO TO 26
      IF (I.EQ.1) GO TO 26
C             CHECK FOR OSCILLATIONS IN CONVERGENCE OF ERROR TERM 
      IF (SUMY(I)-SUMY(I-1)) 26,26,32 
C             LOAD ARRAYS FOR NEXT ITERATION
26    DO 31 J=1,NP
      WK(J,I)=YPPU(J) 
      IF (LL.EQ.1.AND.I.EQ.1) YPPS(J)=YPPU(J) 
      YPP(J)=YPPU(J)
      CC=YUSMO(J) 
      IF (J1-2) 29,28,27
27    AA=YN(J)-YUSMO(J) 
      BB=A(J,1)-A(J,2)
      T1=DSIGN(1.D00,AA)
      T2=DSIGN(1.D00,BB)
      IF (T1.EQ.T2.OR.AA.EQ.BB) GO TO 28
      YUSMO(J)=A(J,2)-BB*(YUSMO(J)-A(J,2))/(AA-BB)
      GO TO 30
28    YUSMO(J)=0.5*(YUSMO(J)+YN(J)) 
      GO TO 30
29    YUSMO(J)=YN(J)
30    A(J,1)=YN(J)
31    A(J,2)=CC 
      GO TO 33
32    NTER=I-1
      ICON=2
      GO TO 36
33    NTER=I
C             CHECK FOR CONVERGENCE BASED ON INPUT EPS
      IF (SUMY(I).LE.EPS) GO TO 35
34    CONTINUE
      GO TO 36
35    ICON=1
C 
C        PRINT SECOND DERIVATIVES GENERATED DURING SMOOTHING PROCESS
C 
36    IF (IPRINT.NE.0) GO TO 38 
      WRITE (JWRITE,80) TITLE 
      DO 37 J=1,NP
37    WRITE (JWRITE,81) J,YSMO(J),(WK(J,I),I=1,NTER)
      WRITE (JWRITE,82) (SUMY(I),I=1,NTER)
38    IF (IPRINT.NE.0) WRITE (JWRITE,79) LL,(SUMY(I),I=1,NTER)
      MTER=MTER+NTER
      IF (ICON.NE.0) GO TO 40 
39    CONTINUE                                   ! end of iteration loop


40    IF (ICON.EQ.0) WRITE (JWRITE,83) MTER 
      IF (ICON.EQ.1) WRITE (JWRITE,84) MTER 
      IF (ICON.EQ.2) WRITE (JWRITE,85) MTER 
C 
C        COMPUTE SMOOTHED SECOND DERIVATIVE USING LEAST SQUARES 
C        CUBIC SPLINE 
C 
      DO 41 I=1,NP
41    DUM(I)=DF 
C             CALL LEAST SQUARES CUBIC SPLINE ROUTINE 
      IPTXXX=-1
      CALL CSDS (LMX,NP,THETA,YPPU,DUM,DFLOAT(NP),
     X  IPTXXX,A,WORKDUMMY,IERR)   ! modified RLC   2021/06/14
      IF (IERR.NE.0) GO TO 71 
C             COMPUTE SECOND DERIVATIVE 
      SUM=0.0 
      DO 44 I=1,NP
      IF (I.EQ.NP) GO TO 42 
      YPP(I)=A(I,1) 
      GO TO 43
42    DELTA=THETA(I)-THETA(I-1) 
      YPP(I)=((A(I-1,4)*DELTA+A(I-1,3))*DELTA+A(I-1,2))*DELTA+A(I-1,1)
43    SUM=SUM+(YPPU(I)-YPP(I))**2 
44    YPPU(I)=YPPS(I) 
      WRITE (JWRITE,88) SUM 
C 
C        COMPUTE NEW Y COORDINATES FROM SMOOTHED SECOND DERIVATIVES 
C 
      CALL YNEW (THETA,YPP,YSMO,NOSE,NP,YUSMO(1),YUSMO(NOSE),YUSMO(NP),E
     1PS,DUM,WK,JWRITE,1) 
C 
C        CHECK NEW Y COORDINATES FOR SMOOTHNESS 
C 
C             CALL LEAST SQUARES POLYNOMIAL ROUTINE 
      DO 45 I=1,NP
45    A(I,1)=1.0
      CALL LSQSMO (THETA,YSMO,A,YN,DUM,YPPS,NP,1,NP,NOSE,WT,EPS,IERR) 
      IF (IERR.NE.0) RETURN 
C             COMPUTE ERROR TERMS 
      SUM1=0.0
      SUM2=0.0
      DO 46 I=1,NP
      A(I,1)=YSMO(I)-YN(I)
      A(I,2)=YPP(I)-YPPS(I) 
      SUM1=SUM1+A(I,1)**2 
46    SUM2=SUM2+A(I,2)**2 
C 
C        COMPUTE FIRST DERIVATIVE FROM SMOOTHED SECOND DERIVATIVE 
C 
      N1=NP-1 
      DO 47 I=1,N1
      DELTA=THETA(I+1)-THETA(I) 
47    YN(I)=-YPP(I)*DELTA/3.-YPP(I+1)*DELTA/6.+(YSMO(I+1)-YSMO(I))/DELTA
      DELTA=THETA(NP)-THETA(N1) 
      YN(NP)=YPP(N1)*DELTA/6.+YPP(NP)*DELTA/3.+(YSMO(NP)-YSMO(N1))/DELTA
C 
C        PRINT SUMMARY OF SMOOTHED AND UNSMOOTHED DATA
C 
48    WRITE (JWRITE,86) TITLE 
      DO 53 I=1,NP
      YPS(I)=YN(I)
      IF (THETA(I).LE.0.) YN(I)=-YN(I)
      T1=DABS(THETA(I))
      IF (T1.GT.PI2) GO TO 49 
      GP=CONS*DSIN(T1) 
      GPP=CONS*DCOS(T1)
      GO TO 50
49    DIF=DCOSH(T1-PI2)
      DELTA=DSINH(T1-PI2)
      GP=CONS/DIF 
      GPP=-CONS*DELTA/(DIF*DIF) 
50    IF (I.EQ.NOSE) GO TO 51 
      DYDX=YN(I)/GP 
      DY2DX=(YPP(I)*GP-YN(I)*GPP)/(GP**3) 
      CURV=DABS(DY2DX)/(DSQRT(1.+DYDX**2)**3) 
      GO TO 52
51    DYDX=0.1D99 
      DY2DX=0.1D99
      CURV=CONS/(YN(I)**2)
      RLE=1./CURV 
52    DELTA=THETA(I)*RAD
      DIF=Y(I)-YSMO(I)
      YPPS(I)=YPP(I)
      WRITE (JWRITE,87) I,DELTA,X(I),Y(I),YUSMO(I),YSMO(I),DIF,YPS(I),
     1   YPP(I),DYDX,DY2DX,CURV
!!! plot data added RLC 2021-06-28     
      WRITE (JGNU2,87) I,DELTA,X(I),Y(I),YUSMO(I),YSMO(I),DIF,YPS(I),
     1   YPP(I),DYDX,DY2DX,SQRT(CURV)   ! note that sqrt(curv) to plot file, not curv
53    END DO     
      WRITE (JWRITE,89) RLE 
C 
C     CHECK FOR INTERSECTION OF UPPER AND LOWER SURFACES
C 
C        DEFINE ITERATION INTERVAL
      KRT=1001
      N1=2*KRT
      TE=THETA(NP)
      TN=-THETA(1)
      IF (TN.LT.TE) TE=TN 
      DIF=TE/DFLOAT(KRT-1) 
      BB=0.5*DIF
      AA=0.85*TE
      YL1=YSMO(NOSE)
      YU1=YSMO(NOSE)
      TP=0.0
      TN=0.0 
      J1=NOSE 
      J2=2
C        DO-LOOP TO SEARCH FOR INTERSECTION 
      DO 59 I=2,N1
      IF (TP.LE.AA) TN=TN+DIF 
      IF (TP.GT.AA) TN=TN+BB
      IF (TN.GT.TE) GO TO 61
      TI=TN 
C        FIND UPPER SURFACE Y-COORDINATE AT THETA = TN
      DO 54 K=J1,NP 
      J=K-1 
      IF (TI.GE.THETA(J).AND.TI.LE.THETA(J+1)) GO TO 55 
54    CONTINUE
55    DELTA=THETA(J+1)-THETA(J) 
      T2=THETA(J+1)-TI
      T1=TI-THETA(J)
      YU2=YPPS(J)*(T2**3/(6.*DELTA)-T2*DELTA/6.)+YPPS(J+1)*(T1**3/(6.*DE
     1LTA)-T1*DELTA/6.)+(YSMO(J)*T2+YSMO(J+1)*T1)/DELTA 
      J1=J
      IF (J1.LT.NOSE) J1=NOSE 
C        FIND LOWER SURFACE Y-COORDINATE AT THETA = TN
      TI=-TN
      DO 56 K=J2,NOSE 
      J=NOSE+1-K
      IF (TI.GE.THETA(J).AND.TI.LE.THETA(J+1)) GO TO 57 
56    CONTINUE
57    DELTA=THETA(J+1)-THETA(J) 
      T2=THETA(J+1)-TI
      T1=TI-THETA(J)
      YL2=YPPS(J)*(T2**3/(6.*DELTA)-T2*DELTA/6.)+YPPS(J+1)*(T1**3/(6.*DE
     1LTA)-T1*DELTA/6.)+(YSMO(J)*T2+YSMO(J+1)*T1)/DELTA 
      J2=NOSE+1-J 
      IF (J2.LT.2) J2=2 
C        COMPUTE THETA FOR INTERSECTION OF STRAIGHT LINE SEGMENTS THRU
C        LAST TWO POINTS ON EACH SURFACE
      CC=(YU2-YU1-YL2+YL1)/(TN-TP)
      IF (DABS(CC).LT.1.D-10) GO TO 58 
      T1=(YL1-YU1)/CC+TP
      IF (I.EQ.2) GO TO 58
C        CHECK TO SEE IF INTERSECTION THETA IS BETWEEN THIS TN-VALUE
C        AND THE PREVIOUS TN-VALUE
      IF (T1.GE.TP.AND.T1.LE.TN) GO TO 60 
C        CONTINUE TO NEXT TN-VALUE
58    YU1=YU2 
      YL1=YL2 
      TP=TN 
59    CONTINUE
      GO TO 61
60    IF (T1.GE.TE) GO TO 61
C        IF INTERSECTION OCCURS WRITE ERROR MESSAGE AND RETURN TO 
C        CALLING PROGRAM
      T1=T1*RAD 
      WRITE (JWRITE,72) T1
      IERR=1
      RETURN
C 
C        FIND LOCATIONS WHERE DY/DX=0.
61    KRT=0 
      N1=NP-1 
      DO 66 I=1,N1
      DELTA=THETA(I+1)-THETA(I) 
      AA=(YPP(I)-YPP(I+1))/(2.*DELTA) 
      IF (AA.EQ.0.0) GO TO 66 
      BB=(YPP(I+1)*THETA(I)-YPP(I)*THETA(I+1))/DELTA
      CC=(YPP(I)*THETA(I+1)**2-YPP(I+1)*THETA(I)**2)/(2.*DELTA)+(YPP(I+1
     1)-YPP(I))*DELTA/6.-(YSMO(I+1)-YSMO(I))/DELTA
      GP=BB*BB-4.*AA*CC 
      IF (GP) 66,62,62
62    GP=DSQRT(GP) 
      T1=(-BB+GP)/(2.*AA) 
      T2=(-BB-GP)/(2.*AA) 
      IF (T1.GE.THETA(I).AND.T1.LE.THETA(I+1)) GO TO 63 
      GO TO 64
63    KRT=KRT+1 
      WK(KRT,1)=T1
64    IF (T2.GE.THETA(I).AND.T2.LE.THETA(I+1)) GO TO 65 
      GO TO 66
65    KRT=KRT+1 
      WK(KRT,1)=T2
66    CONTINUE
      IF (KRT.EQ.0) GO TO 70
C             FIND X/C AND Y/C WHERE DY/DX=0. 
      DO 69 I=1,KRT 
      CALL INTER (WK(I,1),WK(I,2),NP,THETA,X,1,KTI,0) 
      DO 67 J=1,N1
      J1=J
      J2=J+1
      IF (WK(I,1).GE.THETA(J).AND.WK(I,1).LE.THETA(J+1)) GO TO 68 
67    CONTINUE
68    AA=THETA(J2)-WK(I,1)
      BB=WK(I,1)-THETA(J1)
      WK(I,1)=WK(I,1)*RAD 
      DELTA=THETA(J2)-THETA(J1) 
69    WK(I,3)=YPP(J1)*(AA**3/(6.*DELTA)-AA*DELTA/6.)+YPP(J2)*(BB**3/(6.*
     1DELTA)-BB*DELTA/6.)+(YSMO(J1)*AA+YSMO(J2)*BB)/DELTA 
70    CONTINUE
      IF (KRT.GT.0) WRITE (JWRITE,90) (WK(I,2),WK(I,3),WK(I,1),I=1,KRT) 
C 
C        PRINT RESULTS OF SMOOTHNESS CHECK
C 
      IF (ITER.EQ.0) RETURN 
      WRITE (JWRITE,91) TITLE,DF
      WRITE (JWRITE,92) (I,A(I,1),A(I,2),I=1,NP)
      WRITE (JWRITE,93) SUM1,SUM2 
      RETURN
C 
C        PRINT WARNING MESSAGE IF ERROR OCCURRED IN CALL TO CSDS
C 
71    WRITE (JWRITE,94) IERR
      RETURN
C 
72    FORMAT (/5X,108HERROR MESSAGE  ---  SMOOTHING PROCESS RESULTED IN 
     1AN INTERSECTION OF THE UPPER AND LOWER SURFACES AT THETA =,F10.3) 
73    FORMAT (1H1,1X,7HTITLE--,2X,8A10//12X,62H--CHECK OF FIRST DERIVATI
     1VES GENERATED FROM IOP=2 INPUT DATA--///9X,1HI,5X,12HDY/DT(INPUT),
     24X,11HDY/DT(CSDS),8X,3HDIF,6X,13HDY/DT(LSQSMO),8X,3HDIF/) 
74    FORMAT (5X,I5,5(5X,F10.6))
75    FORMAT (/25X,16HSUM OF SQUARES =,4X,F10.6,20X,F10.6)
76    FORMAT (/10X,25HOUTPUT FROM CSDS SELECTED)
77    FORMAT (/10X,27HOUTPUT FROM LSQSMO SELECTED)
78    FORMAT (1H1,1X,7HTITLE--,2X,8A10//30X,53H--SUM OF SQUARES GENERATE
     1D DURING SMOOTHING PROCESS--) 
79    FORMAT (/1X,I5,10F12.7) 
80    FORMAT (1H1,1X,7HTITLE--,2X,8A10//30X,67H--SECOND DERIVATIVES W/R 
     1THETA GENERATED DURING SMOOTHING PROCESS--/4X,1HI,5X,5HTHETA,10(5X
     2,6HDY2/DT)/)
81    FORMAT (I5,F10.2,10F11.6) 
82    FORMAT (/1X,14HSUM OF SQUARES,(10F11.6))
83    FORMAT (/3X,41HSMOOTHING PROCESS HAS NOT CONVERGED AFTER,I4,1X,10H
     1ITERATIONS) 
84    FORMAT (/3X,33HSMOOTHING PROCESS CONVERGED AFTER,I4,1X,10HITERATIO
     1NS) 
85    FORMAT (/3X,41HSMOOTHING PROCESS BEGAN OSCILLATING AFTER,I4,1X,10H
     1ITERATIONS) 
86    FORMAT (1H1,1X,7HTITLE--,2X,8A10//48X,28H--SMOOTHING OUTPUT SUMMAR
     1Y--//4X,1HI,5X,5HTHETA,5X,3HX/C,7X,3HY/C,7X,4HYT/C,5X,6HYSMO/C,4X,
     25HDELTA,7X,3HYPS,6X,4HYPPS,8X,5HDY/DX,7X,11HD(DY/DX)/DX,6X, 
     3  9HCURVATURE)   ! final slash omitted   RLC 20210928
87    FORMAT (I5,F10.2,7F10.6,3E15.6) 
88    FORMAT (/3X,58HSUM OF SQUARES FROM LEAST SQUARES CUBIC SPLINE SMOO
     1THING =,E12.4)
89    FORMAT (/3X,22HLEADING-EDGE RADIUS/C=,F10.6)
90    FORMAT (/3X,16HDY/DX=0. AT X/C=,F10.6,5X,4HY/C=,F10.6,5X,6HTHETA=,
     1F10.3)
91    FORMAT (1H1,1X,7HTITLE--,2X,8A10//12X,29HCHECK OF SMOOTHED COORDIN
     1ATES,3X,3HDF=,F10.6//9X,1HI,5X,20H(YSMO/C-CHECK VALUE),7X,
     218H(YPPS-CHECK VALUE)/) 
92    FORMAT (5X,I5,10X,F10.6,15X,F10.6)
93    FORMAT (/5X,15HSUM OF SQUARES=,F10.6,15X,F10.6) 
94    FORMAT (/3X,21HINPUT ERROR -- POINT ,I3,18H IS NOT INCREASING/) 
      END  SUBROUTINE SMOXY
C*DECK YNEW
      SUBROUTINE YNEW (THETA,YPP,Y,NOSE,NP,YLTE,YNOSE,YUTE,EPS,DUM,WK,JW
     1RITE,IPT) 
C 
C     ROUTINE TO COMPUTE NEW Y/C COORDINATES USING AN ITERATION 
C     PROCEDURE THAT INSURES A DESIRED Y/C COORDINATE AT THE NOSE 
C     (IPT=0) OR THAT INSURES CONTINUITY OF THE FIRST DERIVATIVE W/R TO 
C     THETA AT THE NOSE (IPT=1) 
C 
C CALLED BY SMOXY (4 PLACES); CALLS INVY (4 PLACES)
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C 
C        DIMENSION THETA, YPP, Y, AND DUM BY NP AND WK BY 2*NP IN 
C        CALLING PROGRAM
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION THETA(NP),YPP(NP),Y(NP),DUM(NP),WK(*) 
C 
C     INITIALIZE ITERATION PARAMETERS 
C 
      NMAX=20 
      N1=-1 
      DELTA=0.
      T1=THETA(NOSE)-THETA(NOSE-1)
      T2=THETA(NOSE+1)-THETA(NOSE)
      DO 1 I=1,NP 
1     DUM(I)=YPP(I) 
C 
C     ITERATION LOOP TO COMPUTE INCREMENTAL ADJUSTMENT TO SECOND
C     DERIVATIVE TO INSURE THAT THE DESIRED CONVERGENCE OPTION AT 
C     THE NOSE IS OBTAINED
C 
2     N1=N1+1 
      IF (N1.GT.NMAX) GO TO 11
      IF (IPT.EQ.1) GO TO 3 
C        IF IPT=0, COMPUTE UPPER AND LOWER SURFACE Y/C COORDINATES
C        CONCURRENTLY 
      CALL INVY (THETA,DUM,1,NP,Y,YLTE,YUTE,WK) 
C        COMPUTE DIFFERENCE BETWEEN OUTPUT AND DESIRED Y/C COORDINATE 
C        AT THE NOSE
      DIF=Y(NOSE)-YNOSE 
      GO TO 4 
C        IF IPT=1, COMPUTE UPPER AND LOWER SURFACE Y/C COORDINATES
C        CONSECUTIVELY
3     CALL INVY (THETA,DUM,NOSE,NP,Y,YNOSE,YUTE,WK) 
      CALL INVY (THETA,DUM,1,NOSE,Y,YLTE,YNOSE,WK)
C        COMPUTE DIFFERENCE BETWEEN UPPER AND LOWER SURFACE FIRST 
C        DERIVATIVES AT THE NOSE
      AA=-DUM(NOSE)*T2/3.-DUM(NOSE+1)*T2/6.+(Y(NOSE+1)-Y(NOSE))/T2
      BB=DUM(NOSE-1)*T1/6.+DUM(NOSE)*T1/3.+(Y(NOSE)-Y(NOSE-1))/T1 
      DIF=AA-BB 
C        CHECK FOR CONVERGENCE
4     IF (DABS(DIF).LE.EPS) GO TO 9
C        COMPUTE ADJUSTMENT VALUE TO SECOND DERIVATIVE
      IF (N1.EQ.0) GO TO 6
      IF (DIF.EQ.DIFP) GO TO 5
      SP=(DELTA-DELTAP)/(DIF-DIFP)
      DELTAP=DELTA
      DIFP=DIF
      DELTA=DELTA-DIF*SP
      GO TO 7 
5     DELTA=0.5*(DELTA+DELTAP)
      GO TO 7 
6     DELTAP=DELTA
      DIFP=DIF
      DELTA=DELTA+DIF 
C        ADD ADJUSTMENT VALUE TO SECOND DERIVATIVE
7     DO 8 I=1,NP 
8     DUM(I)=YPP(I)+DELTA 
C        CONTINUE TO ITERATE
      GO TO 2 
C 
C     PRINT CONVERGENCE MESSAGE 
C 
9     WRITE (JWRITE,14) N1,DELTA
C        REDEFINE THE SECOND DERIVATIVE 
      DO 10 I=1,NP
10    YPP(I)=DUM(I) 
      IF (IPT.EQ.1) GO TO 12
      GO TO 13
C 
C     PRINT NON-CONVERGENCE MESSAGE 
C 
11    N1=N1-1 
      WRITE (JWRITE,15) N1
C 
C     COMPUTE NEW UPPER AND LOWER SURFACE Y/C COORDINATES CONCURRENTLY
C 
12    CALL INVY (THETA,YPP,1,NP,Y,YLTE,YUTE,WK) 
C 
C     RETURN TO CALLING PROGRAM 
C 
13    RETURN
C 
14    FORMAT (/3X,88HITERATION PROCEDURE TO COMPUTE INCREMENTAL ADJUSTME
     1NT TO SECOND DERIVATIVE CONVERGED IN ,I3,23H ITERATIONS AND DELTA 
     2=,E12.4)
15    FORMAT (///10X,40HWARNING THE FOLLOWING ERROR HAS OCCURRED//3X,95H
     1ITERATION PROCEDURE TO COMPUTE INCREMENTAL ADJUSTMENT TO SECOND DE
     2RIVATIVE DID NOT CONVERGE IN ,I3,11H ITERATIONS)
      END  SUBROUTINE YNEW
C*DECK INVY
      SUBROUTINE INVY (X,YPP,NS,NE,Y,YSTART,YEND,A) 
C 
C     THIS ROUTINE COMPUTES Y VALUES FROM KNOWN SECOND DERIVATIVES AND
C     END CONDITIONS

C CALLED BY YNEW (4 PLACES);  NO CALLS
C 
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C 
C        IN CALLING PROGRAM DIMENSION X, YPP, AND Y BY NE AND A BY 2*NE 
C
      IMPLICIT REAL*8(A-H,O-Z)
C 
      DIMENSION X(NE),YPP(NE),Y(NE),A(NE,2) 
C 
C     SET END CONDITIONS
C 
      Y(NS)=YSTART
      Y(NE)=YEND
C 
C     PERFORM FORWARD ELIMINATION 
C 
      A(1,1)=YSTART 
      A(1,2)=0.0
      N=NE-NS+1 
      N1=N-1
      DO 1 I=2,N1 
      J=NS+I-1
      H1=X(J)-X(J-1)
      H2=X(J+1)-X(J)
      C=(H1*YPP(J-1)/6.+(H1+H2)*YPP(J)/3.+H2*YPP(J+1)/6.)*H1*H2 
      D=-H2*(A(I-1,2)+1.)-H1
      A(I,2)=H1/D 
1     A(I,1)=(C-H2*A(I-1,1))/D
C 
C     PERFORM BACK SUBSTITUTION 
C 
      J=NE
      DO 2 I=2,N1 
      J=J-1 
      N=N-1 
2     Y(J)=A(N,1)-A(N,2)*Y(J+1) 
C 
C     RETURN TO CALLING PROGRAM 
C 
      RETURN
      END  SUBROUTINE INVY
C*DECK LSQSMO
      SUBROUTINE LSQSMO (X,Y,W,YN,YP,YPP,N,IMAX,JMAX,NOSE,WT,EPS,IERR)
C 
C     THIS SUBROUTINE IS USED TO SMOOTH X AND Y BY CONSECUTIVELY FITTING
C     A LEAST SQUARES POLYNOMIAL OF DEGREE 4 THRU 7 POINTS AT A TIME
C 
C CALLED BY SMOXY (2 PLACES);  NO CALLS
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C 
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(N), Y(N), W(N), YN(N), YP(N), YPP(N)
C 
      DIMENSION XI(7), YI(7), WW(7), A(5,6), B(5) 
C 
      COMMON /INOUT/ JREAD,JWRITE,IPRINT,JGNU1,JGNU2,JGNU3     ! RLC  20210625
      
      
!!!      write(*,*) 'Entering LSQSMO with N=', n
      
C 
C     CHECK NOSE REGION FOR SYMMETRY
C 
      ISYM=1
      DO 1 I=1,3
      IF (DABS(X(NOSE-I)+X(NOSE+I)).GT.EPS) ISYM=0 
      IF (DABS(Y(NOSE-I)+Y(NOSE+I)).GT.EPS) ISYM=0 
1     CONTINUE
      IERR=0
C 
C     FIT A LEAST SQUARES POLYNOMIAL OF DEGREE 4 THRU 7 POINTS
C 
      DO 14 I=1,N 
C        LOAD 7 POINTS FOR LEAST SQUARES POLYNOMIAL FIT 
      IF (I.LT.4) GO TO 2 
      IF (I.GT.N-3) GO TO 3 
      J1=I-3
      J2=I+3
      GO TO 4 
2     J1=1
      J2=7
      GO TO 4 
3     J1=N-6
      J2=N
4     KK=0
      IF (ISYM.EQ.0) GO TO 7
      IF (I.GT.NOSE-3.AND.I.LE.NOSE) GO TO 5
      IF (I.LT.NOSE+3.AND.I.GT.NOSE) GO TO 6
      GO TO 7 
5     J1=NOSE-6 
      J2=NOSE 
      GO TO 7 
6     J1=NOSE 
      J2=NOSE+6 
7     DO 8 L=J1,J2
      J=L 
      IF (I.LE.NOSE) J=J1+J2-L
      KK=KK+1 
      WW(KK)=1.0
      IF (I.EQ.J) WW(KK)=W(I) 
      IF (J.EQ.IMAX.OR.J.EQ.JMAX) WW(KK)=WT*W(J)
      XI(KK)=X(J) 
8     YI(KK)=Y(J) 
      IF (I.LE.4) WW(7)=7.*W(1) 
      IF (I.GE.N-3) WW(7)=7.*W(N) 
C        COMPUTE LEAST SQUARES MATRIX 
      DO 9 L=1,5
      DO 9 J=1,6
9     A(L,J)=0. 
      DO 11 K=1,7 
      T1=1. 
      DO 11 J=1,5 
      T2=T1 
      DO 10 L=1,5 
      A(J,L)=A(J,L)+T2*WW(K)
10    T2=T2*XI(K) 
      A(J,6)=A(J,6)-YI(K)*T1*WW(K)
11    T1=T1*XI(K) 
C        SOLVE FOR COEFFICIENTS OF LEAST SQUARES POLYNOMIAL 
      DO 12 K=1,4 
      DO 12 J=K,4 
      T1=A(J+1,K)/A(K,K)
      DO 12 L=K,6 
12    A(J+1,L)=A(J+1,L)-A(K,L)*T1 
      B(5)=-A(5,6)/A(5,5) 
      
      DO 13 L=2,5 
      K=6-L 
      B(K)=-A(K,6)/A(K,K) 
      K1=K+1
      DO 13 J=K1,5
13    B(K)=B(K)-B(J)*A(K,J)/A(K,K)


C        COMPUTE NEW Y , FIRST , AND SECOND DERIVATIVE
      YN(I)=(((B(5)*X(I)+B(4))*X(I)+B(3))*X(I)+B(2))*X(I)+B(1)
      YP(I)=((4.*B(5)*X(I)+3.*B(4))*X(I)+2.*B(3))*X(I)+B(2) 
      YPP(I)=(12.*B(5)*X(I)+6.*B(4))*X(I)+2.*B(3) 
14    CONTINUE
      IF (ISYM.EQ.0) RETURN 
      YN(NOSE)=0.0
      YPP(NOSE)=0.0 
      YP(NOSE)=1.0
      RETURN
      END  SUBROUTINE LSQSMO
C*DECK CSDS
      SUBROUTINE CSDS (MAX,IX,X,F,DF,S,IPT,COEF,WK,IERR)
C***********************************************************************
C*                                                                     *
C*    PURPOSE:                                                         *
C*                 SUBROUTINE CSDS FITS A SMOOTH CUBIC SPLINE TO A     *
C*                 UNIVARIATE FUNCTION.  DATA MAY BE UNEQUALLY SPACED. *
C*                                                                     *
C*        USE:                                                         *
C*                 CALL CSDS(MAX,IX,X,F,DF,S,IPT,COEF,WK,IERR)         *
C*                                                                     *
C*         MAX     INPUT INTEGER SPECIFYING THE MAXIMUM NUMBER OF DATA *
C*                 POINTS FOR THE INDEPENDENT VARIABLE.                *
C*                                                                     *
C*          IX     INPUT INTEGER SPECIFYING THE ACTUAL NUMBER OF DATA  *
C*                 POINTS FOR THE INDEPENDENT VARIABLE.  IX@MAX.       *
C*                                                                     *
C*           X     ONE-DIMENSIONAL INPUT ARRAY DIMENSIONED AT LEAST    *
C*                 IX IN THE CALLING PROGRAM.  UPON ENTRY TO CSDS,     *
C*                 X(I) MUST CONTAIN THE VALUE OF THE INDEPENDENT      *
C*                 VARIABLE AT POINT I.                                *
C*                                                                     *
C*           F     ONE-DIMENSIONAL INPUT ARRAY DIMENSIONED AT LEAST    *
C*                 IX IN THE CALLING PROGRAM.  UPON ENTRY TO CSDS,     *
C*                 F(I) MUST CONTAIN THE VALUE OF THE FUNCTION AT      *
C*                 POINT X(I).                                         *
C*                                                                     *
C*          DF     ONE-DIMENSIONAL INPUT ARRAY DIMENSIONED AT LEAST    *
C*                 IX IN THE CALLING PROGRAM.  UPON ENTRY TO CSDS,     *
C*                 DF(I) MUST CONTAIN AN ESTIMATE OF THE STANDARD      *
C*                 DEVIATION OF F(I).                                  *
C*                                                                     *
C*           S     A NON-NEGATIVE INPUT PARAMETER WHICH CONTROLS THE   *
C*                 EXTENT OF SMOOTHING.  S SHOULD BE IN THE RANGE      *
C*                 (IX-(2*IX)**.5)<S<(IX+(2*IX)**.5).                  *
C*                                                                     *
C*         IPT     INPUT INITIALIZATION PARAMETER.  THE USER MUST      *
C*                 SPECIFY IPT=-1 WHENEVER A NEW X ARRAY IS            *
C*                 INPUT.  THE ROUTINE WILL ALSO CHECK TO INSURE THAT  *
C*                 THE X ARRAY IS IN STRICTLY INCREASING ORDER.        *
C*                                                                     *
C*        COEF     A TWO-DIMENSIONAL OUTPUT ARRAY DIMENSIONED (MAX,4)  *
C*                 IN THE CALLING PROGRAM.  UPON RETURN, COEF(I,J)     *
C*                 CONTAINS THE J-TH COEFFICIENT OF THE SPLINE FOR     *
C*                 THE INTERVAL BEGINNING AT POINT X(I).  THE          *
C*                 FUNCTIONAL VALUE OF THE SPLINE AT ABSCISSA X1,      *
C*                 WHERE X(I)<X1<X(I+1), IS GIVEN BY:                  *
C*                    F(X1)=((COEF(I,4)*H+COEF(I,3))*H+COEF(I,2))*H    *
C*                    +COEF(I,1)                                       *
C*                 WHERE H=X1-X(I)                                     *
C*                                                                     *
C*          WK     A ONE-DIMENSIONAL WORK AREA ARRAY DIMENSIONED AT    *
C*                 LEAST (7*IX+9) IN THE CALLING PROGRAM.              *
C*                                                                     *
C*        IERR     OUTPUT ERROR PARAMETER:                             *
C*                 =0   NORMAL RETURN.  NO ERROR DETECTED.             *
C*                 =J   THE J-TH ELEMENT OF THE X ARRAY IS NOT IN      *
C*                      STRICTLY INCREASING ORDER.                     *
C*                 =-1  THERE ARE LESS THAN FOUR VALUES IN THE X ARRAY.*
C*                                                                     *
C*                 UPON RETURN FROM CSDS, THIS PARAMETER SHOULD BE     *
C*                 TESTED IN THE CALLING PROGRAM.                      *
C*                                                                     *
C*                                                                     *
C*                                                                     *
C*    REQUIRED ROUTINES                -NONE                           *
C*                                                                     *
C*    SOURCE                           IMSL ROUTINE ICSSMU MODIFIED BY *
C*                                     COMPUTER SCIENCES CORPORATION   *
C*                                                                     *
C*    LANGUAGE                         -FORTRAN                        *
C*                                                                     *
C*    DATE RELEASED                    SEPTEMBER 5, 1973               *
C*                                                                     *
C*    LATEST REVISION                  MARCH 1975                      *
C***********************************************************************
C 
C CALLED BY SMOXY (2 PLACES)
C 
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(IX), F(IX), DF(IX), COEF(MAX,4), WK(2000)
C 
C                                  SET UP WORKING AREAS 
C 
!!!      write(*,*) 'Entering CSDS, ipt,ix,max=',ipt,ix,max 
      IERR=0
!!!      write(*,*) 'ierr set to zero', ierr
      IF (IPT.NE.-1) GO TO 4
      IPT=0 
!!!      write(*,*) 'ipt set to zero', ipt,ierr
      IF (IX.LT.4) GO TO 1
      GO TO 2 
 1    IERR=-1 
      write(*,*) 'Error return from CSDS=', ierr 
      RETURN
      
 2    IX1=IX-1
!!!      write(*,*) 'checking that x-array is increasing' 
      DO 3 I=1,IX1
      IF (X(I+1)-X(I).GT.0) GO TO 3 
      IERR=I+1
      write(*,*) 'Error return from CSDS=', ierr      
      RETURN
      
 3    CONTINUE
!!!      write(*,*) 'proceeding 3 continue in CSDS'
      NP1=IX+1
      IB1=NP1 
      IB2=IB1+NP1 
      IB3=IB2+NP1+1 
      IB4=IB3+NP1 
      IB5=IB4+NP1 
      IB6=IB5+NP1+1 
!!!      write(*,*) 'Set IB1...IB6', ib1,ib2,ib3,ib4,ib,ib6
      WK(1)=0.
      WK(2)=0.
      WK(IB2)=0.
      WK(IB3)=0.
      IJK2=IB2+NP1
      WK(IJK2)=0. 
      IJK5=IB5+1
      WK(IJK5)=0. 
      IJK5=IB5+2
      WK(IJK5)=0. 
      WK(IB6)=0.
      IJK5=IB5+NP1
      WK(IJK5)=0. 
!!!      write(*,*) 'Finished ipt=-1 setup in csds.', WK(1:ijk5)
      
4     CONTINUE
!!!      write(*,*) 'proceeding at 4 continue'
      P=0.
      H=X(2)-X(1) 
      F2=-S 
      FF=(F(2)-F(1))/H
      IF (IX.LT.3) GO TO 10 
      DO 5 I=3,IX 
      G=H 
      H=X(I)-X(I-1) 
      E=FF
      FF=(F(I)-F(I-1))/H
      COEF(I-1,1)=FF-E
      IJK3=IB3+I
      WK(IJK3)=(G+H)*.66666666666667
      IJK4=IB4+I
      WK(IJK4)=H/3. 
      IJK2=IB2+I
      WK(IJK2)=DF(I-2)/G
      WK(I)=DF(I)/H 
      IJK1=IB1+I
      WK(IJK1)=-DF(I-1)/G-DF(I-1)/H 
5     CONTINUE
      DO 6 I=3,IX 
      IJK1=IB1+I
      IJK2=IB2+I
      COEF(I-1,2)=WK(I)*WK(I)+WK(IJK1)*WK(IJK1)+WK(IJK2)*WK(IJK2) 
      COEF(I-1,3)=WK(I)*WK(IJK1+1)+WK(IJK1)*WK(IJK2+1)
      COEF(I-1,4)=WK(I)*WK(IJK2+2)
6     CONTINUE
C 
C                                  NEXT ITERATION 
C 
7     IF (IX.LT.3) GO TO 10 
      DO 8 I=3,IX 
      IJK1=IB1+I-1
      IJK0=I-1
      WK(IJK1)=FF*WK(IJK0)
      IJK2=IB2+I-2
      IJK0=I-2
      WK(IJK2)=G*WK(IJK0) 
      IJK0=I
      IJK3=IB3+I
      WK(IJK0)=1./(P*COEF(I-1,2)+WK(IJK3)-FF*WK(IJK1)-G*WK(IJK2)) 
      IJK5=IB5+I
      IJKN=IJK5-1 
      IJK0=IJKN-1 
      WK(IJK5)=COEF(I-1,1)-WK(IJK1)*WK(IJKN)-WK(IJK2)*WK(IJK0)
      IJK4=IB4+I
      FF=P*COEF(I-1,3)+WK(IJK4)-H*WK(IJK1)
      G=H 
      H=COEF(I-1,4)*P 
8     CONTINUE
      DO 9 I=3,IX 
      J=IX-I+3
      IJK5=IB5+J
      IJK6=IJK5+1 
      IJK7=IJK6+1 
      IJK1=IB1+J
      IJK2=IB2+J
      WK(IJK5)=WK(J)*WK(IJK5)-WK(IJK1)*WK(IJK6)-WK(IJK2)*WK(IJK7) 
9     CONTINUE
10    E=0 
      H=0 
C 
C                                  COMPUTE U AND ACCUMULATE E 
C 
      DO 11 I=2,IX
      G=H 
      IJK5=IB5+I
      H=(WK(IJK5+1)-WK(IJK5))/(X(I)-X(I-1)) 
      IJK6=IB6+I
      WK(IJK6)=(H-G)*DF(I-1)*DF(I-1)
      E=E+WK(IJK6)*(H-G)
11    CONTINUE
      G=-H*DF(IX)*DF(IX)
      IJK6=IB6+NP1
      WK(IJK6)=G
      E=E-G*H 
      G=F2
      F2=E*P*P
      IF (F2.GE.S.OR.F2.LE.G) GO TO 14
      FF=0. 
      IJK6=IB6+2
      H=(WK(IJK6+1)-WK(IJK6))/(X(2)-X(1)) 
      IF (IX.LT.3) GO TO 13 
      DO 12 I=3,IX
      G=H 
      IJK6=IB6+I
      H=(WK(IJK6+1)-WK(IJK6))/(X(I)-X(I-1)) 
      IJK1=IB1+I-1
      IJK2=IB2+I-2
      G=H-G-WK(IJK1)*WK(I-1)-WK(IJK2)*WK(I-2) 
      FF=FF+G*WK(I)*G 
      WK(I)=G 
12    CONTINUE
13    H=E-P*FF
      IF (H.LE.0) GO TO 14
C 
C                                  UPDATE THE LAGRANGE MULTIPLIER P 
C                                     FOR THE NEXT ITERATION
C 
      P=P+(S-F2)/((DSQRT(S/E)+P)*H)
      GO TO 7 
C 
C                                  IF E LESS THAN OR EQUAL TO S,
C                                  COMPUTE THE COEFFICIENTS AND RETURN. 
C 
14    DO 15 I=2,NP1 
      IJK6=IB6+I
      COEF(I-1,1)=F(I-1)-P*WK(IJK6) 
      IJK5=IB5+I
      COEF(I-1,3)=WK(IJK5)
15    CONTINUE
      DO 16 I=2,IX
      H=X(I)-X(I-1) 
      COEF(I-1,4)=(COEF(I,3)-COEF(I-1,3))/(3.*H)
      COEF(I-1,2)=(COEF(I,1)-COEF(I-1,1))/H-(H*COEF(I-1,4)+COEF(I-1,3))*
     1H 
16    CONTINUE

      write(*,*) 'CSDS completed, ierr=', ierr
      RETURN
      END  SUBROUTINE CSDS
C*DECK PCARD 
      SUBROUTINE PCARD (IPUNCH,X,Y,W,THETA,YSMO,YPS,YPPS,NOSE,NP,CHORD,T
     1ITLE) 
C 
C     ROUTINE TO PUNCH OUTPUT DATA   (TAPE 1 IS PUNCH FILE) 
C 
C CALLED BY MAIN PROGRAM; NO CALLS
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C 
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER TITLE(8)*10
C
      DIMENSION X(*), Y(*), W(*), THETA(*), YSMO(*), YPS(*), YPPS(*)
C 
      COMMON /HLM/ DX(200),DY(200),DW(200), DUMMY8(1400)
C 
      COMMON /BLK1/ PI,PI2,RAD,CONS 
C 
      COMMON /INOUT/ JREAD,JWRITE,IPRINT,JGNU1,JGNU2,JGNU3     ! RLC  20210625
C 
      IF (IPUNCH.LE.0.OR.IPUNCH.GE.5) RETURN
C 
C     PUNCH TITLE CARD
C 
      WRITE (JWRITE,10) IPUNCH,TITLE
      WRITE (1,11) TITLE
C 
C     DETERMINE OUTPUT PUNCH OPTION 
C 
      IOP=0 
      IF (IPUNCH.EQ.2) IOP=1
      IF (IPUNCH.EQ.3) IOP=2
      IF (IPUNCH.EQ.4) IOP=3
      WRITE (JWRITE,12) IOP 
      XI=DFLOAT(IOP) 
      WRITE (1,13) XI 
C 
C     PUNCH UPPER SURFACE QUANTITIES
C 
      J=0
      KP=0
      DO 1 I=NOSE,NP
      J=J+1 
      DW(J)=W(I)
      IF (W(I).GT.1.0) KP=1 
      IF (IOP.EQ.0) DX(J)=X(I)*CHORD
      IF (IOP.NE.0) DX(J)=THETA(I)*RAD
      IF (IOP.EQ.0) DY(J)=YSMO(I)*CHORD 
      IF (IOP.EQ.1) DY(J)=YSMO(I) 
      IF (IOP.EQ.2) DY(J)=YPS(I)
      IF (IOP.EQ.3) DY(J)=YPPS(I) 
1     CONTINUE
      WRITE (JWRITE,14) J 
      XI=DFLOAT(J) 
      WRITE (1,15) XI 
      IF (IOP.EQ.0) WRITE (JWRITE,16) (DX(I),I=1,J) 
      IF (IOP.NE.0) WRITE (JWRITE,7) (DX(I),I=1,J)
      WRITE (JWRITE,17) (DY(I),I=1,J) 
      IF (KP.EQ.1) WRITE (JWRITE,21) (DW(I),I=1,J)
      DO 3 I=1,J
      IF (IOP.NE.0) GO TO 2 
      IF (DW(I).GT.1.0) WRITE (1,22) DX(I),DY(I),DW(I)
      IF (DW(I).LE.1.0) WRITE (1,18) DX(I),DY(I)
      GO TO 3 
2     IF (DW(I).GT.1.0) WRITE (1,8) DX(I),DY(I),DW(I) 
      IF (DW(I).LE.1.0) WRITE (1,9) DX(I),DY(I) 
3     CONTINUE
C 
C     PUNCH LOWER SURFACE QUANTITIES
C 
      J=0
      KP=0
      DO 4 I=1,NOSE 
      J=J+1 
      K=NOSE+1-I
      DW(J)=W(K)
      IF (W(K).GT.1.0) KP=1 
      IF (IOP.EQ.0) DX(J)=X(K)*CHORD
      IF (IOP.NE.0) DX(J)=THETA(K)*RAD
      IF (IOP.EQ.0) DY(J)=YSMO(K)*CHORD 
      IF (IOP.EQ.1) DY(J)=YSMO(K) 
      IF (IOP.EQ.2) DY(J)=YPS(K)
      IF (IOP.EQ.3) DY(J)=YPPS(K) 
4     CONTINUE
      WRITE (JWRITE,19) J 
      XI=DFLOAT(J) 
      WRITE (1,15) XI 
      IF (IOP.EQ.0) WRITE (JWRITE,16) (DX(I),I=1,J) 
      IF (IOP.NE.0) WRITE (JWRITE,7) (DX(I),I=1,J)
      WRITE (JWRITE,17) (DY(I),I=1,J) 
      IF (KP.EQ.1) WRITE (JWRITE,21) (DW(I),I=1,J)
      DO 6 I=1,J
      IF (IOP.NE.0) GO TO 5 
      IF (DW(I).GT.1.0) WRITE (1,22) DX(I),DY(I),DW(I)
      IF (DW(I).LE.1.0) WRITE (1,18) DX(I),DY(I)
      GO TO 6 
5     IF (DW(I).GT.1.0) WRITE (1,8) DX(I),DY(I),DW(I) 
      IF (DW(I).LE.1.0) WRITE (1,9) DX(I),DY(I) 
6     CONTINUE
C 
C     PUNCH YLTE AND YUTE 
C 
      IF (IOP.LE.1) RETURN
      YLTE=YSMO(1)
      YNOSE=YSMO(NOSE)
      YUTE=YSMO(NP) 
      WRITE (JWRITE,20) YLTE,YNOSE,YUTE 
      WRITE (1,18) YLTE,YNOSE,YUTE
C 
C     RETURN TO CALLING PROGRAM 
C 
      RETURN
C 
7     FORMAT (/3X,4HTH =,8F10.5/(7X,8F10.5))
8     FORMAT (F10.5,F10.6,F10.2)
9     FORMAT (F10.5,F10.6)
10    FORMAT (1H1,10X,36HTHE FOLLOWING DATA HAVE BEEN PUNCHED,5X,7HIPUNC
     1H=,I4//3X,8A10) 
11    FORMAT (8A10) 
12    FORMAT (/3X,5HIOP =,I4) 
13    FORMAT (30X,F10.2)
14    FORMAT (/3X,4HNU =,I4)
15    FORMAT (F10.2)
16    FORMAT (/3X,4HDX =,8F10.6/(7X,8F10.6))
17    FORMAT (/3X,4HDY =,8F10.6/(7X,8F10.6))
18    FORMAT (3F10.6) 
19    FORMAT (/3X,4HNL =,I4)
20    FORMAT (/3X,6HYLTE =,F10.6,5X,7HYNOSE =,F10.6,5X,6HYUTE =,F10.6)
21    FORMAT (/3X,4HDW =,8F10.2/(7X,8F10.2))
22    FORMAT (2F10.6,F10.2) 

      END  SUBROUTINE PCARD
      
      
      FUNCTION F(X1,X2,X3,X4,X5,X6,X7,X8,X9)
      IMPLICIT REAL*8(A-H,O-Z)
      F=X1*(X5*X9-X6*X8)+X2*(X6*X7-X4*X9)+X3
     1*(X4*X8-X5*X7)
      RETURN
      END FUNCTION F
      
      
C*DECK CAMTK 
      SUBROUTINE CAMTK (THETA,YSMO,YPPS,NOSE,NP,EPS,KPLOT,IPUNCH,TITLE) 
C 
C     THIS SUBROUTINE COMPUTES THE THICKNESS AND CAMBER DISTRIBUTIONS 
C     OF THE SMOOTHED AIRFOIL 
C C CALLED BY MAIN PROGRAM ; NO CALLS
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C
      IMPLICIT REAL*8(A-H,O-Z)
C 
      CHARACTER TITLE(8)*10
C
      DIMENSION THETA(*), YSMO(*), YPPS(*)
C 
      COMMON /SMY/ TU(100),YPPU(100),TL(100),YPPL(100),DYXU(100),LX(102)
     1,XLS(101),YLS(101),TH(101),XU(102),YU(102),XL(102),YL(102),XC(103)
     2,YC(103),TK(103),DUMMY9(559)
C 
      COMMON /BLK1/ PI,PI2,RAD,CONS 
C 
      COMMON /INOUT/ JREAD,JWRITE,IPRINT,JGNU1,JGNU2,JGNU3     ! RLC  20210625
C 
      DATA NM/2001/,SIZ/.40/,ISIZ/3/
C 
C 
C        LOAD THETA, X/C, Y/C, AND SECOND DERIVATIVES INTO SEPARATE 
C        ARRAYS FOR UPPER AND LOWER SURFACES
C 
      J=0 
      NU=NP-NOSE+1
      DO 2 I=NOSE,NP
      J=J+1 
      TU(J)=THETA(I)
      YU(J)=YSMO(I) 
      TP=DABS(THETA(I))
      IF (TP.GT.PI2) GO TO 1
      XU(J)=CONS*(1.-DCOS(TP)) 
      GO TO 2 
1     XU(J)=CONS*(DATAN(DSINH(TP-PI2))+1.)
2     YPPU(J)=YPPS(I) 
      NL=NOSE 
      J=NOSE+1
      DO 4 I=1,NOSE 
      J=J-1 
      TL(J)=THETA(I)
      YL(J)=YSMO(I) 
      TP=DABS(THETA(I))
      IF (TP.GT.PI2) GO TO 3
      XL(J)=CONS*(1.-DCOS(TP)) 
      GO TO 4 
3     XL(J)=CONS*(DATAN(DSINH(TP-PI2))+1.)
4     YPPL(J)=YPPS(I) 
C             COMPUTE FIRST DERIVATIVES OF UPPER SURFACE
      DO 5 I=2,NU 
      DELTA=TU(I)-TU(I-1) 
      DYXU(I)=YPPU(I)*DELTA/3.+YPPU(I-1)*DELTA/6.+(YU(I)-YU(I-1))/DELTA 
      IF (TU(I).LE.PI2) DYXU(I)=DYXU(I)/(CONS*DSIN(TU(I))) 
      IF (TU(I).GT.PI2) DYXU(I)=DYXU(I)*DCOSH(TU(I)-PI2)/CONS
5     CONTINUE
      DYXU(1)=0.1D99
C 
C        COMPUTE THICKNESS AND CAMBER DISTRIBUTIONS BY FINDING LOWER
C        SURFACE COORDINATE (XLS,YLS) CORRESPONDING TO INPUT UPPER
C        SURFACE COORDINATE (XU,YU) 
C 
      NT=0
      KSAVE=1 
      NS=1
      NL1=NL-1
      NM1=NM-1
      A1=PI/DFLOAT(NM1)
      DEL=1./(DFLOAT(NM1)**2)
      DO 12 I=1,NU
C             LOAD XU AND YU
      IJ=NU+1-I 
      XXU=XU(IJ)
      YYU=YU(IJ)
      DYU=DYXU(IJ)
      NN=1
C             FIND XLS
      DO 9 K=NS,NM
      TP=A1*DFLOAT(NM-K) 
      IF (K.EQ.1) TP=DABS(TL(NL))
      IF (K.EQ.NM) TP=DABS(TL(1))
      IF (TP.LE.PI2) XXL=CONS*(1.-DCOS(TP))
      IF (TP.GT.PI2) XXL=CONS*(DATAN(DSINH(TP-PI2))+1.) 
      IF (NN.EQ.NL) NN=NL1
      DO 6 J=NN,NL1 
      J2=NL-J 
      J1=J2+1 
      IF (TP.GE.DABS(TL(J2)).AND.TP.LE.DABS(TL(J1))) GO TO 7
6     CONTINUE
7     DELTA=TL(J2)-TL(J1) 
      T1=-TP-TL(J1) 
      T2=TL(J2)+TP
      YYL=YPPL(J1)*(T2**3/(6.*DELTA)-T2*DELTA/6.)+YPPL(J2)*(T1**3/(6.*DE
     1LTA)-T1*DELTA/6.)+(YL(J1)*T2+YL(J2)*T1)/DELTA 
      DYL=YPPL(J1)*(DELTA/6.-T2*T2/(2.*DELTA))+YPPL(J2)*(T1*T1/(2.*DELTA
     1)-DELTA/6.)+(YL(J2)-YL(J1))/DELTA 
      IF (TP.LE.PI2) DELTA=CONS*DSIN(TP) 
      IF (TP.GT.PI2) DELTA=CONS/DCOSH(TP-PI2)
      IF (TP.LE.0.0) DYL=0.1D99 
      IF (TP.GT.0.0) DYL=-DYL/DELTA 
      NN=NL+1-J1
      D=DSQRT((XXL-XXU)**2+(YYL-YYU)**2) 
      IF (I.EQ.1.AND.D.LE.DEL) GO TO 10 
      IF (D.LE.DEL) GO TO 9 
      COST=(YYU-YYL)/D
      SINT=(XXL-XXU)/D
      IF (DYU.NE.0.1D99) DU=(COST*DYU-SINT)/(SINT*DYU+COST) 
      IF (DYU.EQ.0.1D99.AND.SINT.NE.0.0) DU=COST/SINT 
      IF (DYU.EQ.0.1D99.AND.SINT.EQ.0.0) DU=0.1D99
      IF (DYL.NE.0.1D99) DL=-(COST*DYL-SINT)/(SINT*DYL+COST)
      IF (DYL.EQ.0.1D99.AND.SINT.NE.0.0) DL=-COST/SINT
      IF (DYL.EQ.0.1D99.AND.SINT.EQ.0.0) DL=-0.1D99 
      IF (K.EQ.NS) GO TO 8
      DKL=(DL-DLP)/(XXL-XP) 
      DKU=(DU-DUP)/(XXL-XP) 
      IF (DKU.EQ.DKL) GO TO 8 
      XK=XP+(DLP-DUP)/(DKU-DKL) 
      IF (XK.LE.XP+DEL.AND.XK.GE.XXL-DEL) GO TO 11
8     KSAVE=K 
      XP=XXL
      DUP=DU
      DLP=DL
9     CONTINUE
      IF (I.GT.1) GO TO 12
10    XK=XL(NL) 
      KSAVE=NS
11    NT=NT+1 
      LX(NT)=IJ 
      XLS(NT)=XK
      NS=KSAVE
12    CONTINUE
C             COMPUTE YLS FOR EACH XLS AND PRINT RESULTS
      WRITE (JWRITE,44) TITLE 
      DO 19 I=1,NT
      IJ=LX(I)
      DELTA=XLS(I)
      IF (DELTA.GT.1.) DELTA=1. 
      IF (DELTA.LE.CONS) GO TO 13 
      DELTA=DTAN(DELTA/CONS-1.)
      TP=PI2+DLOG(DELTA+DSQRT(DELTA*DELTA+1.)) 
      GO TO 14
13    TP=DACOS(1.-DELTA/CONS)
14    DO 15 J=1,NL1 
      J2=NL-J 
      J1=J2+1 
      IF (TP.GE.DABS(TL(J2)).AND.TP.LE.DABS(TL(J1))) GO TO 16 
15    CONTINUE
16    DELTA=TL(J2)-TL(J1) 
      T1=-TP-TL(J1) 
      T2=TL(J2)+TP
      YYL=YPPL(J1)*(T2**3/(6.*DELTA)-T2*DELTA/6.)+YPPL(J2)*(T1**3/(6.*DE
     1LTA)-T1*DELTA/6.)+(YL(J1)*T2+YL(J2)*T1)/DELTA 
      YLS(I)=YYL
      XC(I)=(XU(IJ)+XLS(I))/2.
      YC(I)=(YU(IJ)+YYL)/2. 
      TK(I)=0.5*DSQRT((XU(IJ)-XLS(I))**2+(YU(IJ)-YYL)**2)
      IF (YU(IJ).EQ.YYL) TH(I)=0.0
      IF (YU(IJ).NE.YYL) TH(I)=DATAN((XLS(I)-XU(IJ))/(YU(IJ)-YYL)) 
      IF (TK(I).LE.0.0) GO TO 17
      DYL=YPPL(J1)*(DELTA/6.-T2*T2/(2.*DELTA))+YPPL(J2)*(T1*T1/(2.*DELTA
     1)-DELTA/6.)+(YL(J2)-YL(J1))/DELTA 
      IF (TP.LE.PI2) DELTA=CONS*DSIN(TP) 
      IF (TP.GT.PI2) DELTA=CONS/DCOSH(TP-PI2)
      IF (TP.LE.0.0) DYL=0.1D99 
      IF (TP.GT.0.0) DYL=-DYL/DELTA 
      COST=(YU(IJ)-YYL)/(2.*TK(I))
      SINT=(XLS(I)-XU(IJ))/(2.*TK(I)) 
      DU=(COST*DYXU(IJ)-SINT)/(SINT*DYXU(IJ)+COST)
      DL=(COST*DYL-SINT)/(SINT*DYL+COST)
      T2=DABS(DABS(DU)-DABS(DL)) 
      GO TO 18
17    T2=0.0
18    T1=TH(I)*RAD
      WRITE (JWRITE,45) I,XU(IJ),YU(IJ),XLS(I),YYL,XC(I),YC(I),TK(I),T1,
     1T2
!!! plot file added  RLC   2021-06-28     
      WRITE (JGNU3,45) I,XU(IJ),YU(IJ),XLS(I),YYL,XC(I),YC(I),TK(I),T1,
     1T2
     
19    CONTINUE
C 
C        COMPUTE STARTING LOCATION OF CAMBER DISTRIBUTION (I.E. 
C        THICKNESS = 0) BY FITTING SECOND ORDER CURVE TO LAST THREE 
C        COMPUTED CAMBER LINE COORDINATES AND THEN DETERMINING
C        INTERSECTION OF THAT CURVE WITH AIRFOIL SURFACE
C 
      ISYM=1
      DO 20 I=1,5 
      IF (DABS(XU(I)-XL(I)).GT.EPS) ISYM=0 
      IF (DABS(YU(I)+YL(I)).GT.EPS) ISYM=0 
20    CONTINUE

      IF (ISYM.EQ.1) GO TO 30 
      IF (XC(NT).LE.DEL) GO TO 31 
      X1=XC(NT)**2
      X2=XC(NT-1)**2
      X3=XC(NT-2)**2
      
CCC replaced 1. with ONE (set to 1.0) to avoid REAL*4-REAL*8 mismatch  RLC 2021Jun12
      ONE=1.0      
      D=F(X1,XC(NT),ONE,X2,XC(NT-1),ONE,X3,XC(NT-2),ONE) 
      A1=F(YC(NT),XC(NT),ONE, YC(NT-1),XC(NT-1),ONE, 
     X YC(NT-2),XC(NT-2),ONE)/D
      A2=F(X1,YC(NT),ONE, X2,YC(NT-1),ONE, X3,YC(NT-2),ONE)/D
      A3=YC(NT)-A1*X1-A2*XC(NT) 
      NM1=NM/4
      D=XC(NT)/DFLOAT(NM1) 
      X=0.0 
      XP=X
      YYUP=YU(1)
      YYLP=YL(1)
      YYCP=(A1*X+A2)*X+A3 
      NM1=NM1+1 
      DO 27 I=2,NM1 
      X=X+D 
      IF (X.GT.CONS) GO TO 27 
      TP=DACOS(1.-X/CONS)
      DO 21 K=2,NU
      K1=K-1
      K2=K
      IF (TP.GE.TU(K1).AND.TP.LE.TU(K2)) GO TO 22 
21    CONTINUE

22    DELTA=TU(K2)-TU(K1) 
      T1=TP-TU(K1)
      T2=TU(K2)-TP
      YYU=YPPU(K1)*(T2**3/(6.*DELTA)-T2*DELTA/6.)+YPPU(K2)*(T1**3/(6.*DE
     1LTA)-T1*DELTA/6.)+(YU(K2)*T1+YU(K1)*T2)/DELTA 
      DO 23 J=2,NL
      J2=J-1
      J1=J
      IF (TP.GE.DABS(TL(J2)).AND.TP.LE.DABS(TL(J1))) GO TO 24 
23    CONTINUE

24    DELTA=TL(J2)-TL(J1) 
      T1=-TP-TL(J1) 
      T2=TL(J2)+TP
      YYL=YPPL(J1)*(T2**3/(6.*DELTA)-T2*DELTA/6.)+YPPL(J2)*(T1**3/(6.*DE
     1LTA)-T1*DELTA/6.)+(YL(J1)*T2+YL(J2)*T1)/DELTA 
      YYC=(A1*X+A2)*X+A3
      DKC=(YYC-YYCP)/(X-XP) 
      DKU=(YYU-YYUP)/(X-XP) 
      IF (DKU.EQ.DKC) GO TO 25
      XKU=XP+(YYCP-YYUP)/(DKU-DKC)
      IF (XKU.GE.XP.AND.XKU.LE.X) GO TO 28
25    DKL=(YYL-YYLP)/(X-XP) 
      IF (DKL.EQ.DKC) GO TO 26
      XKL=XP+(YYCP-YYLP)/(DKL-DKC)
      IF (XKL.GE.XP.AND.XKL.LE.X) GO TO 29
26    XP=X
      YYLP=YYL
      YYUP=YYU
      YYCP=YYC
27    CONTINUE
      GO TO 31
      
28    NT=NT+1 
      LX(NT)=0
      XLS(NT)=XKU 
      XC(NT)=XKU
      DU=(A1*XKU+A2)*XKU+A3 
      TK(NT)=0. 
      TH(NT)=DATAN(2.*A1*XKU+A2) 
      TP=DACOS(1.-XKU/CONS)
      DELTA=TU(K2)-TU(K1) 
      T1=TP-TU(K1)
      T2=TU(K2)-TP
      YYU=YPPU(K1)*(T2**3/(6.*DELTA)-T2*DELTA/6.)+YPPU(K2)*(T1**3/(6.*DE
     1LTA)-T1*DELTA/6.)+(YU(K2)*T1+YU(K1)*T2)/DELTA 
      YLS(NT)=YYU 
      YC(NT)=YLS(NT)
      D=DABS(DABS(DU)-DABS(YC(NT)))
      T1=TH(NT)*RAD 
      WRITE (JWRITE,45) NT,XLS(NT),YLS(NT),XLS(NT),YLS(NT),XC(NT),YC(NT)
     1,TK(NT),T1,D
      GO TO 31
      
29    NT=NT+1 
      LX(NT)=0
      XLS(NT)=XKL 
      XC(NT)=XKL
      DL=(A1*XKL+A2)*XKL+A3 
      TK(NT)=0. 
      TH(NT)=DATAN(2.*A1*XKL+A2) 
      TP=DACOS(1.-XKL/CONS)
      DELTA=TL(J2)-TL(J1) 
      T1=-TP-TL(J1) 
      T2=TL(J2)+TP
      YYL=YPPL(J1)*(T2**3/(6.*DELTA)-T2*DELTA/6.)+YPPL(J2)*(T1**3/(6.*DE
     1LTA)-T1*DELTA/6.)+(YL(J1)*T2+YL(J2)*T1)/DELTA 
      YLS(NT)=YYL 
      YC(NT)=YLS(NT)
      D=DABS(DABS(DL)-DABS(YC(NT)))
      T1=TH(NT)*RAD 
      WRITE (JWRITE,45) NT,XLS(NT),YLS(NT),XLS(NT),YLS(NT),XC(NT),YC(NT)
     1,TK(NT),T1,D
      GO TO 31
      
30    IF (LX(NT).EQ.1) GO TO 31 
      NT=NT+1 
      LX(NT)=1
      XC(NT)=0.0
      YC(NT)=YU(1)
      XLS(NT)=0.0 
      YLS(NT)=YL(1) 
      TK(NT)=0.0
      TH(NT)=0.0
      D=0.0 
      WRITE (JWRITE,45) NT,XC(NT),YC(NT),XLS(NT),YLS(NT),XC(NT),YC(NT),T
     1K(NT),TH(NT),D
C 
C        PUNCH CAMBER AND THICKNESS DISTRIBUTIONS 
C 
31    IF (IPUNCH.NE.5) GO TO 33 
      WRITE (1,46) TITLE
      WRITE (JWRITE,41) IPUNCH,TITLE,NT 
C 
      D=DFLOAT(NT) 
      WRITE (1,42) D
C 
      DO 32 I=1,NT
      J=NT+1-I
      WRITE (JWRITE,43) XC(J),YC(J),TK(J),TH(J) 
      WRITE (1,47) XC(J),YC(J),TK(J),TH(J)
32    CONTINUE

   33 WRITE(JWRITE,*) " "     ! added 20210928 by RLC
      RETURN                  ! adds blank line after data
      
41    FORMAT (1H1,5X,47HTHE FOLLOWING CAMBERLINE DATA HAVE BEEN PUNCHED,
     15X,7HIPUNCH=,I4//5X,8A10//5X,4HNT =,I4//9X,3HX/C,7X,3HY/C,5X,5HT/C
     2/2,5X,5HSLOPE)
42    FORMAT (F10.2)
43    FORMAT (5X,4F10.6)
44    FORMAT (1H1,1X,7HTITLE--,2X,8A10//32X,37H--THICKNESS AND CAMBER DI
     1STRIBUTION--//4X,1HI,5X,4HXU/C,6X,4HYU/C,6X,4HXL/C,6X,4HYL/C,6X,3H
     2X/C,7X,3HY/C,6X,5HT/C/2,5X,5HSLOPE,10X,5HERROR)
!      final slash on line above omitted   RLC 20210928     
45    FORMAT (I5,7F10.6,F10.4,5X,F10.6) 
46    FORMAT (8A10) 
47    FORMAT (4F10.6) 
      END  SUBROUTINE CAMTK
C*DECK INTP
      SUBROUTINE INTP (THETA,X,YSMO,YPPS,NP,NOSE,CHORD,TITLE,NINT,XINT,C
     1NEW,INTR,IPUNCH)
C 
C     ROUTINE TO INTERPOLATE ADDITIONAL UPPER AND LOWER SURFACE 
C     COORDINATES 
C 
C CALLED BY MAIN PROGRAM;  CALLS COORD
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C 
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER TITLE(8)*10
C
      DIMENSION THETA(*), X(*), YSMO(*), YPPS(*), XINT(*) 
C 
      DIMENSION XSAV(57)
C 
      COMMON /INOUT/ JREAD,JWRITE,IPRINT,JGNU1,JGNU2,JGNU3     ! RLC  20210625
C 
      COMMON /BLK1/ PI,PI2,RAD,CONS 
C 
      COMMON /HLM/ XU(100),YU(100),XL(100),YL(100),TLS(100),DUMMY9(1500)
C 
C          STANDARD X/C COORDINATE INTERPOLATION VALUES 
      DATA (XSAV(I),I=1,57)/0.0,.00025,.0005,.00075,.001,.0015,.002,.002
     15,.005,.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.125,.15,.175,.2,.22
     25,.25,.275,.3,.325,.35,.375,.4,.425,.45,.475,.5,.525,.55,.575,.6,.
     3625,.65,.675,.7,.725,.75,.775,.8,.825,.85,.875,.9,.925,.95,.97,.98
     4,.99,1.0/ 
C 
C     IF INTR EQUAL 1, LOAD STANDARD X/C COORDINATE VALUES
C 
      IF (INTR.NE.2) CNEW=CHORD 
      IF (INTR.EQ.2) GO TO 2
      NINT=57 
      DO 1 I=1,NINT 
      XINT(I)=XSAV(I) 
1     END DO
C 
C     INTERPOLATE UPPER SURFACE COORDINATES 
C 
2     WRITE (JWRITE,7) TITLE
      XUP=X(NP)*CHORD 
      XNOSE=X(NOSE)*CHORD 
      XLO=X(1)*CHORD
      RATIO=CNEW/CHORD
      DO 5 I=1,NINT 
      XU(I)=XINT(I)*CHORD+XNOSE 
      XL(I)=XU(I) 
      IF (XU(I).GT.XUP) XU(I)=XUP 
      IF (XL(I).GT.XLO) XL(I)=XLO 
      XU(I)=(XU(I)-XNOSE)*RATIO 
      XL(I)=(XL(I)-XNOSE)*RATIO 
      DELTA=XINT(I) 
      IF (DELTA.LE.CONS) GO TO 3
      DELTA=DTAN(DELTA/CONS-1.)
      TU=PI2+DLOG(DELTA+DSQRT(DELTA*DELTA+1.)) 
      GO TO 4 
3     TU=DACOS(1.-DELTA/CONS)
4     TL=-TU
      IF (TL.LT.THETA(1)) TL=THETA(1) 
      IF (TU.GT.THETA(NP)) TU=THETA(NP) 
      TLS(I)=TL  
      CALL COORD (THETA,YPPS,YSMO,NP,TU,YU(I),DYDX,DY2DX,CURV)
      YU(I)=YU(I)*CNEW
      WRITE (JWRITE,8) I,XU(I),YU(I),DYDX,DY2DX,CURV
C     033195
C     mss added this line for interpolated coordinates
!!!      WRITE (30,1045) XU(I), YU(I)    !!! commented out RLC  20210928
 1045 FORMAT(2X,2F12.5)
5     CONTINUE
      WRITE (JWRITE,9) CNEW 
C 
C     INTERPOLATE LOWER SURFACE COORDINATES 
C 
      WRITE (JWRITE,10) TITLE 
      DO 6 I=1,NINT 
      TL=TLS(I) 
      CALL COORD (THETA,YPPS,YSMO,NP,TL,YL(I),DYDX,DY2DX,CURV)
      YL(I)=YL(I)*CNEW
      WRITE (JWRITE,8) I,XL(I),YL(I),DYDX,DY2DX,CURV
C     033195
C     mss added lines below
!!!      WRITE(30,1045) XL(I), YL(I)        !!! commented out  RLC 20210715
6     CONTINUE

C.... Write a blank line to indicate end of data
      WRITE(JWRITE,*) " "    ! added 20210928   RLC
C 
C     PUNCH COORDINATES 
C 
      IF (IPUNCH.NE.6) RETURN 
      WRITE (JWRITE,11) CNEW,TITLE
      WRITE (1,12) TITLE
      WRITE (JWRITE,13) NINT
      XNT=DFLOAT(NINT) 
      WRITE (1,14) XNT
      WRITE (JWRITE,15) (XU(I),I=1,NINT)
      WRITE (JWRITE,16) (YU(I),I=1,NINT)
      WRITE (1,17) (XU(I),YU(I),I=1,NINT) 
      WRITE (JWRITE,18) NINT
      WRITE (1,14) XNT
      WRITE (JWRITE,19) (XL(I),I=1,NINT)
      WRITE (JWRITE,20) (YL(I),I=1,NINT)
      WRITE (1,17) (XL(I),YL(I),I=1,NINT) 
C 
C     RETURN TO CALLING PROGRAM 
C 
      RETURN
C 
7     FORMAT (1H1,5X,9HTITLE--  ,8A10//26X,42H--UPPER SURFACE INTERPOLAT
     1ED COORDINATES--//9X,1HI,10X,2HXU,13X,2HYU,11X,5HDY/DX,6X,11HD(DY/
     2DX)/DX,6X,9HCURVATURE)
8     FORMAT (I10,2F15.6,3E15.6)
9     FORMAT (/10X,7HCHORD =,F10.6) 
10    FORMAT (1H1,5X,9HTITLE--  ,8A10//26X,42H--LOWER SURFACE INTERPOLAT
     1ED COORDINATES--//9X,1HI,10X,2HXL,13X,2HYL,11X,5HDY/DX,6X,11HD(DY/
     2DX)/DX,6X,9HCURVATURE)
11    FORMAT (1H1,10X,50HTHE FOLLOWING DATA HAVE BEEN PUNCHED FOR A CHOR
     1D =,F10.6//3X,9HTITLE--  ,8A10) 
12    FORMAT (8A10) 
13    FORMAT (5X,4HNU =,I4) 
14    FORMAT (F10.2)
15    FORMAT (5X,4HXU =,8F10.6/(9X,8F10.6)) 
16    FORMAT (5X,4HYU =,8F10.6/(9X,8F10.6)) 
17    FORMAT (2F10.6) 
18    FORMAT (5X,4HNL =,I4) 
19    FORMAT (5X,4HXL =,8F10.6/(9X,8F10.6)) 
20    FORMAT (5X,4HYL =,8F10.6/(9X,8F10.6)) 
      END  SUBROUTINE INTP
      
C*DECK COORD 
      SUBROUTINE COORD (THETA,YPPS,YSMO,NP,TI,YI,DYDX,DY2DX,CURV) 
C 
C     ROUTINE TO COMPUTE THE Y COORDINATE, DY/DX, D(DY/DX)/DX, AND
C     CURVATURE AT A GIVEN VALUE OF THETA 
C 
C CALLED BY INTP;  NO CALLS
C        CODED BY -- HARRY MORGAN  NASA/LARC/TAD/AAB       1982 
C 
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION THETA(*), YPPS(*), YSMO(*)
C 
      COMMON /BLK1/ PI,PI2,RAD,CONS 
C 
      DO 1 K=2,NP 
      J=K-1 
      IF (TI.GE.THETA(J).AND.TI.LE.THETA(K)) GO TO 2
1     CONTINUE

2     DELTA=THETA(J+1)-THETA(J) 
      T2=THETA(J+1)-TI
      T1=TI-THETA(J)
      YI=YPPS(J)*(T2**3/(6.*DELTA)-T2*DELTA/6.)+YPPS(J+1)*(T1**3/(6.*DEL
     1TA)-T1*DELTA/6.)+(YSMO(J)*T2+YSMO(J+1)*T1)/DELTA
      YPI=YPPS(J)*(DELTA/6.-T2*T2/(2.*DELTA))+YPPS(J+1)*(T1*T1/(2.*DELTA
     1)-DELTA/6.)+(YSMO(J+1)-YSMO(J))/DELTA 
      YPPI=(YPPS(J)*T2+YPPS(J+1)*T1)/DELTA
      DELTA=YPI 
      IF (TI.LE.0.0) DELTA=-DELTA 
      TP=DABS(TI)
      IF (TP.GT.PI2) GO TO 3
      GP=CONS*DSIN(TP) 
      GPP=CONS*DCOS(TP)
      GO TO 4 
3     T1=DCOSH(TP-PI2) 
      T2=DSINH(TP-PI2) 
      GP=CONS/T1
      GPP=-CONS*T2/(T1*T1)
4     IF (TP.LE.0.0.OR.GP.EQ.0.0) GO TO 5 
      DYDX=DELTA/GP 
      DY2DX=(YPPI*GP-DELTA*GPP)/(GP**3) 
      CURV=DABS(DY2DX)/(DSQRT(1.+DYDX**2)**3) 
      RETURN
5     DYDX=0.1D99 
      DY2DX=0.1D99
      CURV=CONS/(DELTA*DELTA) 
      RETURN
      END  SUBROUTINE COORD
      
