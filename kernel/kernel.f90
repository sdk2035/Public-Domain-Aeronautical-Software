!+
      PROGRAM Kernel

! NOTES (RLC) - Removed the overlay structure and replaced with subroutine
!  calls according to the following overlay scheme:
!     1,0   OMDLIB
!     2,0   OMANRD
!     2,1   OREDPR
!     2,2   OGEOTY
!     3,0   OMANAE
!     3,4   TRANUN
!     3,7   TRANST
!     4,0   AEROPR
!     4,1   STYCOF
!     4,2   UNSYCO
!     4,3   QGENF

! Replaced Break with Breakpoint to avoid conflict with Fortran intrinsic
! Replaced all EOF and EXIT calls with modern Fortran


!      OVERLAY(R2T,0,0)
!      PROGRAM MAIN(INPUT=66,OUTPUT=66,PUNCH=66,TAPE5=INPUT,TAPE6=
!     1OUTPUT,TAPE7=PUNCH,TAPE1,TAPE9,TAPE10,TAPE4,TAPE2,TAPE3)
!     MAIN PROGRAM TO COMPUTE STEADY AND UNSTEADY AERODYNAMICS FOR
!     FLUTTER AND DYNAMIC RESPONSE ANALYSES OF ARRAYS OF MULTIPLE
!     SURFACES IN SUBSONIC, TRANSONIC, AND SUPERSONIC FLOW.
!            BY ATLEE M. CUNNINGHAM, JR.
!            REFS. AIAA PAPERS 71-329, 73-670, 74-359  AND 75-99
!                  NASA CR-112264 (1973) AND NASA CR-144895 (1976)
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), HMACH
      COMMON /OPTION/ IOP1, IOP2, IOP3, IOP4, NUNIT, IQT, IOPLU, IND
!
      COMMON /GAMA/ GAMMA
!!!      COMMON /FLTLE/ TITLE(32)
      CHARACTER(LEN=132):: title
      COMMON /FLTLE/ title   !!! TITLE(32)

      COMMON /ALPVCT/DUM1(4885)
      COMMON /COMM/LEAVE,IM,DUM2(10)
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10),XIBAR(15,10), XIS(15,31,10), YBAR(15,10),          &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)

!!!      EQUIVALENCE ( NAME, OMDLIB, OMANRD, OMANAE, AEROPR, INTRPQ,       &
!!!     &              FLTTER )

!!!      DATA NAME /4HR2T /, ENDE /4HEND /
      INTEGER:: errCode
      CHARACTER(LEN=80):: fileName
      CHARACTER(LEN=3):: fend
!----------------------------------------------------------------------------
      WRITE(*,*) 'kernel'
      DO
        WRITE(*,*) 'Enter the name of the input file:'
        READ(*,'(A)') fileName
        IF (LEN_TRIM(fileName)==0) STOP
        OPEN(UNIT=5,FILE=fileName,STATUS='OLD',                          &
     &    IOSTAT=errCode, ACTION='READ')
        IF (errCode==0) EXIT
        OPEN(UNIT=5,FILE=TRIM(fileName)//'.inp',STATUS='OLD',            &
     &    IOSTAT=errCode, ACTION='READ')
        IF (errCode==0) EXIT
        WRITE(*,*) 'Unable to open this file. Try again.'
      END DO
      OPEN(UNIT=6,FILE='kernel.out',STATUS='REPLACE',ACTION='WRITE')

      CALL GSTART (MARSHA)
      FREQ = 0.0
      IJOB = 1
      LEAVE = 0
      CALL LIB()
      CALL FREEFD()
!
!!!      CALL OVERLAY ( OMDLIB,1,0 )
      CALL Omdlib()

      IF ( LEAVE .GT. 0 ) GO TO 10
      GO TO 100
!
   10 STOP    !!! CALL EXIT
!
  100 IJOB = IJOB + 1
      CALL PROB(leave)
!
!!!      CALL OVERLAY ( OMANRD,2,0 )
      CALL Omanrd()

      IF ( LEAVE .NE. 0 ) GO TO 10
!
!
!!!      IF ( IOP1 .LE. 0 ) CALL OVERLAY ( OMANAE,3,0 )
      IF (IOP1 .LE. 0) CALL Omanae()
!
      IF ( IOP1 .NE. 0 ) REWIND(UNIT=NUNIT)
!!!      CALL OVERLAY ( AEROPR,4,0 )
      CALL Aeropr()


      IF ( IOP1 .LT. 0 ) END FILE NUNIT
      IF ( IOP1 .NE. 0 ) REWIND NUNIT
!
 1000 CALL FORMFD()
      READ (5,'(A)',IOSTAT=errCode) FEND
      IF (errCode < 0) GO TO 200
!!!      IF ( EOF(5) .NE. 0 ) GO TO 200
!!!  150 CONTINUE
!!! 1001 FORMAT ( A4 )
      IF ( FEND .NE. 'END' ) GO TO 100
!
  200 WRITE (6,1002)
 1002 FORMAT (1H1,////,21X,28H " " " " " " " " " " " " "  ,/,           &
     &                 21X,28H "                       "  ,/,           &
     &                 21X,28H " THIS JOB IS COMPLETE  "  ,/,           &
     &                 21X,28H "                       "  ,/,           &
     &                 21X,28H " " " " " " " " " " " " "     )
!
      GO TO 10
      END Program Kernel

      SUBROUTINE PROB ( LEAVE )
!!! 'fixed' this routine by commenting the next line.
!!!      IF ( EOF(5) .EQ. 0 ) GO TO 20
   10 LEAVE = 9999
      WRITE ( 6,15 )
   15 FORMAT (1H1,//,33H     " " " " " " " " " " " " " "                &
     &           , /,33H     "                         "                &
     &           , /,33H     "  THIS JOB IS COMPLETE   "                &
     &           , /,33H     "        (PROB)           "                &
     &           , /,33H     " " " " " " " " " " " " " "         )
   20 RETURN
      END Subroutine Prob


      SUBROUTINE FREEFD
      RETURN
      END Subroutine FreeFd
      SUBROUTINE FORMFD
      RETURN
      END Subroutine FormFd
      SUBROUTINE LIB
      RETURN
      END Subroutine Lib

      SUBROUTINE GSTART (M)
      INTEGER,INTENT(OUT):: m
      M = 1
      WRITE (6,10)
   10 FORMAT(1H1,14X,"  GENERAL DYNAMICS-CONVAIR AEROSPACE DIVISION ",  &
     &        //,14X,"    STEADY AND UNSTEADY AERODYNAMIC PROGRAM   ",  &
     &        //,14X,"  FOR SUBSONIC, TRANSONIC, AND SUPERSONIC FLOW",  &
     &        //,14X,"           WITH INTERFERENCE                  "  )
      RETURN
      END Subroutine Gstart

      SUBROUTINE STATUS (IA)
      DIMENSION IA(20)
      DO 10 I=1,20
   10 IA(I) = 1
      RETURN
      END Subroutine Status

      FUNCTION FETE(I1,I2,X,BRPT)
!
!  FUNCTION TO PERFORM A LINEAR INTERPOLATION BETWEEN BREAK CHORD POINTS
!
      DIMENSION BRPT(2,40)
      FETE=(BRPT(1,I1)-BRPT(1,I2))*(X-BRPT(2,I2))                       &
     &       /(BRPT(2,I1)-BRPT(2,I2))+BRPT(1,I2)
!
      RETURN
      END Function Fete

      SUBROUTINE CALMDS ( X, Y, BREF, NX, NY, ALP, DNWSH, DEF,          &
     &                    IOPT, NST, A, NXRO )
!
!     SUBROUTINE TO CALCULATE  ALPHAS ------- IOPT = 1
!                              DOWNWASH ----- IOPT = 2
!                              DEFLECTIONS -- IOPT = 3
!
      COMMON /ALPVCT/ COEF(20,200), XF(200), YF(200), NFS(10), NFS1(10),&
     &   GMASS(20,20), NSTRS, NMODES, DH, DW1, DW2, SPAN2(10),          &
     &   TAN12(10), TAN22(10), CAVG2(10), XR(10), YR(10)
!
      DIMENSION X(NXRO,NY),   ALP(100,20), DNWSH(2,70,20), DEF(100,20), &
     &          A(100), Y(NY)
!
!     NST IS THE STRUCTURAL SURFACE ASSIGNED TO
!     THE CURRENT AERODYNAMIC SURFACE.
!
      NF3 = 1
      DO 5 IS=1,NST
      IF ( IS .LT. NST ) NF3 = NF3 + NFS(IS) + 3
    5 END DO
!
      NF1 = NFS1(NST)
      NF2 = NFS(NST)
      SPAN = SPAN2(NST)
      TAN1 = TAN12(NST)
      TAN2 = TAN22(NST)
      CAVG = CAVG2(NST)
      IF ( IOPT-2 ) 10, 50, 100
!
   10 IR = 0
      DO 20 J=1,NY
!
      Y2 = 2.0*( Y(J)*BREF - YR(NST) )/SPAN - 1.0
      CETA = CAVG - ( TAN1 - TAN2 )*Y2
      XLE  = XR(NST) + Y2*TAN1
!
      DO 20 I=1,NX
      IR = IR + 1
!
      X2   = 2.0*( X(I,J)*BREF - XLE )/CETA - 1.0
!
      CALL MAGN ( X2, Y2, NF2, NMODES, XF(NF1), YF(NF1), DEF(1,1),      &
     &            DNWSH(1,1,1), COEF(1,NF3), 20 )
!
      X2   = X2 + 0.001
!
      CALL MAGN ( X2, Y2, NF2, NMODES, XF(NF1), YF(NF1), A,             &
     &            DNWSH(1,1,1), COEF(1,NF3), 20 )
!
      DELT = 2000.0/CETA
      DO 20 K=1,NMODES
   20 ALP(IR,K) = DELT*( A(K) - DEF(K,1) )*DW1
      GO TO 1000
!
   50 IR = 0
      DO 70 J=1,NY
!
      Y2 = 2.0*( Y(J)*BREF - YR(NST) )/SPAN - 1.0
!
      CETA = CAVG - ( TAN1 - TAN2 )*Y2
      XLE  = XR(NST) + Y2*TAN1
!
      DO 70 I=1,NX
      IR = IR + 1
!
      X2 = 2.0*( X(I,J)*BREF - XLE )/CETA - 1.0
!
      CALL MAGN ( X2, Y2, NF2, NMODES, XF(NF1), YF(NF1), DEF(1,1),      &
     &            A, COEF(1,NF3), 20 )
!
      X2 = X2 + 0.001
!
      CALL MAGN ( X2, Y2, NF2, NMODES, XF(NF1), YF(NF1), ALP(1,1),      &
     &            A, COEF(1,NF3), 20 )
!
      DELT = 2000.0/CETA
      DO 70 K=1,NMODES
      DNWSH(1,IR,K) = DELT*( ALP(K,1) - DEF(K,1) )*DW1
   70 DNWSH(2,IR,K) = DEF(K,1)*DW2
      GO TO 1000
!
  100 IR = 0
      DO 120 J=1,NY
!
      Y2 = 2.0*( Y(J)*BREF - YR(NST) )/SPAN - 1.0
!
      CETA = CAVG - ( TAN1 - TAN2 )*Y2
      XLE  = XR(NST) + Y2*TAN1
!
      DO 120 I=1,NX
      IR = IR + 1
      X2 = 2.0*( X(I,J)*BREF - XLE )/CETA - 1.0
!
      CALL MAGN ( X2, Y2, NF2, NMODES, XF(NF1), YF(NF1), A,             &
     &            ALP(1,1), COEF(1,NF3), 20 )
!
      DO 120 K=1,NMODES
  120 DEF(IR,K) = A(K)*DH
!
 1000 RETURN
      END Subroutine Calmds


      SUBROUTINE SPLINE ( XF, YF, NF, C, D, IA, IB, NROWS, IND )
!
!     SUBROUTINE TO CALCULATE AND INVERT THE INTERPOLATION MATRIX FOR
!     SPLINE FIT OF HARDER AND DESMARAIS.
!
      DIMENSION XF(NROWS), YF(NROWS), C(NROWS,NROWS), D(NROWS),         &
     &          IA(NROWS), IB(NROWS)
!
      N  = NF
      N1 = NF+1
      N2 = NF+2
      N3 = NF+3
!
      DO 40 J=1,N
      C(J,J) = 0.0
      J1 = J + 1
      IF ( J1 .GT. N ) GO TO 15
      DO 10 I=J1,N
      RIJ = ( XF(I) - XF(J) )**2 + ( YF(I) - YF(J) )**2
      IF ( RIJ .EQ. 0.0 ) GO TO 100
   10 C(I,J) = RIJ*ALOG( RIJ )
   15 J1 = J - 1
      IF ( J1 .LT. 1 ) GO TO 30
      DO 20 I=1,J1
   20 C(I,J) = C(J,I)
!
   30 C(N1,J) = 1.0
      C(N2,J) = XF(J)
      C(N3,J) = YF(J)
      C(J,N1) = 1.0
      C(J,N2) = XF(J)
      C(J,N3) = YF(J)
!
   40 END DO
!
      DO 50 J=N1,N3
      DO 50 I=N1,N3
   50 C(I,J) = 0.0
!
      CALL INVRT ( C, D, IA, IB, N3, NROWS, +0, IND )
!
      RETURN
!
  100 WRITE (6,101) I, J
  101 FORMAT (1H1,20H  STRUCTURAL POINTS , 2I4,26H HAVE IDENTICAL LOCATI&
     &ONS  )
      STOP   !!! CALL EXIT
!
      END Subroutine Spline

!----------------------------------------------------------------------
      SUBROUTINE MAGN ( X, Y, LFS, NMDS, XF, YF, VECT, D, C, NROWS )
!
      DIMENSION XF(NROWS), YF(NROWS), C(NROWS,NROWS), D(NROWS),         &
     &          VECT(NMDS)
!
      LFS3 = LFS + 3
      DO 10 K=1,LFS
      RK = ( X - XF(K) )**2 + ( Y - YF(K) )**2 + 1.0E-8
   10 D(K) = RK*ALOG(RK)
      D(LFS+1) = 1.0
      D(LFS+2) = X
      D(LFS+3) = Y
!
      DO 30 J=1,NMDS
      VECT(J) = 0.0
      DO 20 K=1,LFS3
   20 VECT(J) = VECT(J) + C(J,K)*D(K)
   30 END DO
!
      RETURN
      END Subroutine Magn

      SUBROUTINE CALCM ( XI, ETA, NC, NS, BMACH, TMACH, TN, A )
!
!     SUBROUTINE TO EVALUATE THE MACH NUMBER, TMACH, AT POINT XI, ETA,
!     FROM THE TSCHEBYCHEV POLYNOMIAL COEFICIENTS IN BMACH.
!
      DIMENSION  BMACH(NS,NC), TN(NS), A(NC)
!
      TN(1) = 1.0
      TN(2) = ETA
      ETA2  = 2*ETA
      DO 10 J=3,NS
   10 TN(J) = ETA2*TN(J-1) - TN(J-2)
!
      DO 20 I=1,NC
      A(I) = 0.0
      DO 20 J=1,NS
   20 A(I) = A(I) + TN(J)*BMACH(J,I)
!
      TN(2) = XI
      XI2   = 2*XI
      DO 30 I=3,NC
   30 TN(I) = XI2*TN(I-1) - TN(I-2)
!
      TMACH = 0.0
      DO 40 I=1,NC
   40 TMACH = TMACH + TN(I)*A(I)
!
      RETURN
      END Subroutine Calcm

      SUBROUTINE XMBY ( XM, BY, YY, BRPT, NEND )
!
      DIMENSION BRPT(2,40)
!
      Y = ABS( YY )
      DO 10 I=1,NEND
      NUM=2*I
      I1=NUM+1
      IF (Y .LE. BRPT(2,I1)) GO TO 11
   10 END DO
      GO TO 12
   11 I2=NUM-1
      XLE=FETE(I1,I2,Y,BRPT)
!
      I1=NUM+2
      I2=NUM
      XTE=FETE(I1,I2,Y,BRPT)
!
      BY = 0.5*(XTE - XLE)
      XM = 0.5*(XTE + XLE)
!
      RETURN
!
   12 WRITE (6,20) Y, ( (BRPT(I,I2), I=1,2 ), I2=1,I1 )
   20 FORMAT (1H1, 19H    A VALUE FOR Y =,E15.7,20H  HAS BEEN SPECIFIED &
     &       ,//, 52H     WHICH EXCEEDS THE INPUT GEOMETRY SPECIFICATION&
     &S , ///,    35H           X(I)                Y(I),/, (/,2E20.7) )
!
      STOP   !!! CALL EXIT
      END Subroutine Xmby

      SUBROUTINE TNXI ( TN, XI, JJ, MBAR, ICH )
!
      DIMENSION TN(10,22), XI(22)
!
      IF ( ICH .GT. 4 ) GO TO 200
      DO 100 J=1,JJ
      XI2 = 2.0*XI(J)
      GO TO ( 1, 2, 3, 4 ), ICH
!
    1 F = 1.0
      GO TO 10
    2 F = SQRT ( 1.0 + XI(J) )
      GO TO 10
    3 F = 1.0/SQRT( 1.0 - XI(J) )
      GO TO 10
    4 F = SQRT ( (1.0+XI(J))/(1.0-XI(J)) )
!
   10 TN(1,J) = 1.0*F
      TN(2,J) = ( XI2 + 1.0 )*F
      TN(3,J) = ( (XI2 + 1.0)*XI2 - 1.0 )*F
      IF ( MBAR .LT. 4 ) GO TO 100
!
      DO 20 I=4,MBAR
   20 TN(I,J) = XI2*TN(I-1,J) - TN(I-2,J)
!
  100 END DO
      GO TO 1000
!
  200 DO 220 J=1,JJ
!
      XI2 = 2*XI(J)
      TN(1,J) = 1.0
      TN(2,J) = XI2 + 1.0
!
      DO 210 I=3,MBAR
  210 TN(I,J) = XI2*TN(I-1,J) - TN(I-2,J)
!
  220 TN(1,J) = 1.0
!
 1000 RETURN
      END Subroutine Tnxi

      SUBROUTINE UNETA ( UN, ETA, NS, NR, LS, ISTYPE, LSPAN )
!
      DIMENSION UN(15,31), ETA(61)
!
      IF ( LSPAN .GT. 4 ) GO TO 300
!
      IF ( ISTYPE .GT. 1 ) GO TO 110
!
      IF ( LS .NE. 0 ) GO TO 50
!
      DO 20 I=1,NS
      UN(1,I) = 1.0
      UM1     = 2.0*ETA(I)
      UN(2,I) = 2.0*ETA(I)*UM1 - 1.0
      DO 20 J=3,NR
      UM1     = 2.0*ETA(I)*UN(J-1,I) - UM1
   20 UN(J,I) = 2.0*ETA(I)*UM1       - UN(J-1,I)
      GO TO 1000
!
   50 DO 60 I=1,NS
      UN(1,I) = 2.0*ETA(I)
      UM1     = 2.0*ETA(I) - 1.0
      UN(2,I) = 2.0*ETA(I)*UM1 - UN(1,I)
      DO 60 J=3,NR
      UM1     = 2.0*ETA(I)*UN(J-1,I) - UM1
   60 UN(J,I) = 2.0*ETA(I)*UM1       - UN(J-1,I)
      GO TO 1000
!
  110 NS2 = NS
      A1 = 1.0
      A2 = 1.0
      IF ( LSPAN .EQ. 3 ) A2 = -1.0
      IF ( ISTYPE .EQ. 3 ) NS2 = NS/2
      DO 230 I=1,NS2
      IF ( LSPAN .GT. 2 ) A1 = SQRT( 1.0 + A2*ETA(I) )
      UN(1,I) = A1
      UN(2,I) = A1*2.0*ETA(I)
      DO 230 J=3,NR
  230 UN(J,I) = 2.0*ETA(I)*UN(J-1,I) - UN(J-2,I)
!
      IF ( NS .EQ. 1 ) GO TO 1000
      IF ( ISTYPE .NE. 3 ) GO TO 1000
!
      NS3 = NS2+ 1
      A3 = 1.0
      IF ( LS .NE. 0 ) A3 = -1.0
      A1 = A1*A3
!
      DO 240 I=NS3,NS
      IF ( LSPAN .GT. 2 ) A1 = SQRT( 1.0 + A2*ETA(I) )*A3
      UN(1,I) = A1
      UN(2,I) = A1*2.0*ETA(I)
      DO 240 J=3,NR
  240 UN(J,I) = 2.0*ETA(I)*UN(J-1,I) - UN(J-2,I)
!
      GO TO 1000
!
  300 NS2 = NS
      IF ( ISTYPE .EQ. 3 .AND. NS .NE. 1 ) NS2 = NS/2
      A1 = 1.0
      N1 = 1
      N2 = NS2
  310 DO 320 I=N1,N2
      UN(1,I) = A1
      ETA2    = 2.0*ETA(I)
      UN(2,I) = A1*ETA2
      DO 320 J=3,NR
  320 UN(J,I) = ETA2*UN(J-1,I) - UN(J-2,I)
      IF ( N2 .EQ. NS ) GO TO 1000
      N1 = N2+1
      N2 = NS
      IF ( LS .NE. 0 ) A1 = -1.0
      GO TO 310
!
 1000 RETURN
      END Subroutine Uneta

      SUBROUTINE CHDTSS ( W, XIB, XIBAR, JJ, ICH, LSPAN, IND, KINT,     &
     &                    XMACH, BRPT, ETAS, JSUROP )
!
!     SUBROUTINE TO COMPUTE THE CHORDWISE INTEGRATION WEIGHTING
!     FUNCTION FOR SUPERSONIC FLOW AND COMPUTE THE SUPERSONIC
!     PRESSURE BASE FUNCTION FROM CONICAL FLOW THEORY.
!
      DIMENSION XIB(JJ), XIBAR(JJ), W(JJ), BRPT(2,4)
!
      COMMON /ETC/ XIS(67), TRM(2,3), XSC(4), YSC(4)
!
      IF ( JSUROP .EQ. 0 ) GO TO 300
      IF ( IND .NE. 0 ) GO TO 100
      IND  = 1
      BETA = SQRT( XMACH*XMACH - 1.0 )
      XSC(1) = BRPT(1,1)
      XSC(2) = BRPT(1,3)
      XSC(3) = BRPT(1,2)
      XSC(4) = BRPT(1,4)
      YSC(1) = BRPT(2,1)
      YSC(2) = BRPT(2,3)
      YSC(3) = BRPT(2,2)
      YSC(4) = BRPT(2,4)
!
      A1 = XSC(2) - XSC(1)
      A2 = YSC(2) - YSC(1)
      A1 = A1/A2
      A2 = ( XSC(4) - XSC(3) )/A2
      TRM(1,3) = A1
      TRM(2,3) = A2
!
      IF ( A1 .EQ. 0.0 ) GO TO 25
      IST = ICH
      TRM(1,1) = BETA/A1
!
      IF ( TRM(1,1) .GT. 1.0 ) GO TO 20
!
      TRM(1,2) = 7.0/( ( 1.75 + 1.0/TRM(1,1) )*BETA )
      GO TO 30
!
   20 TRM(1,2) = 4.0*TRM(1,1)/( BETA*SQRT( TRM(1,1)**2 - 1.0 ) )
!
      TRM(2,2) = 1.0 - 7.0/( BETA*TRM(1,2)*( 1.75 + 1.0/TRM(1,1) ) )
!
      GO TO 30
!
   25 TRM(1,1) = 1.0E+7
      IST = 0
!
   30 TRM(2,1) = 1.0E+7
      IF ( A2 .NE. 0.0 ) TRM(2,1) = BETA/A2
!
  100 XLE = XSC(1) + ETAS*TRM(1,3)
      XTE = XSC(3) + ETAS*TRM(2,3)
      XM  = 0.5*( XTE + XLE )
      BY  = 0.5*( XTE - XLE )
!
      DO 110 J=1,JJ
  110 XIS(J) = XM + BY*XIB(J)
!
      CALL SSFUNC( XSC, YSC, XIS, ETAS, JJ, BETA, TRM, IST, W )
!
      IF ( KINT .NE. 0 ) GO TO 200
!
      DO 120 J=1,JJ
  120 W(J) = W(J)*SQRT( 1.0 - XIBAR(J)**2 )
!
  200 GO TO 400
!
  300 IF ( KINT .NE. 0 ) GO TO 340
      IF ( XIB(JJ) .EQ. XIBAR(JJ) ) GO TO 320
!
      DO 310 J=1,JJ
      GO TO ( 304, 303, 302, 301 ), ICH
!
  301 F = 1.0
      GO TO 310
  302 F = 1.0/SQRT( 1.0 + XIB(J) )
      GO TO 310
  303 F = SQRT( 1.0 - XIB(J) )
      GO TO 310
  304 F = SQRT( ( 1.0 - XIB(J) )/( 1.0 + XIB(J) ) )
!
  310 W(J) = F*SQRT( 1.0 - XIBAR(J)**2 )
      GO TO 400
!
  320 DO 330 J=1,JJ
!
      GO TO ( 324, 323, 322, 321 ), ICH
!
  321 W(J) = SQRT( 1.0 - XIB(J)**2 )
      GO TO 330
!
  322 W(J) = SQRT( 1.0 - XIB(J) )
      GO TO 330
!
  323 W(J) = ( 1.0 - XIB(J) )*SQRT( 1.0 + XIB(J) )
      GO TO 330
!
  324 W(J) = ( 1.0 - XIB(J) )
!
  330 END DO
      GO TO 400
!
  340 DO 350 J=1,JJ
      GO TO ( 344, 343, 342, 341 ), ICH
!
  341 W(J) = 1.0
      GO TO 350
  342 W(J) = 1.0/SQRT( 1.0 + XIB(J) )
      GO TO 350
  343 W(J) = SQRT( 1.0 - XIB(J) )
      GO TO 350
  344 W(J) = SQRT( ( 1.0 - XIB(J) )/(  1.0 + XIB(J) ) )
!
  350 END DO
!
  400 RETURN
      END Subroutine Chdtss

      SUBROUTINE SSFUNC ( XSC, YSC, XIS, ETAS, JJ, BETA, TRM, IST, P )
!
!     SUBROUTINE TO COMPUTE THE SUPERSONIC WEIGHTING FUNCTION FOR
!     THE ASSUMED PRESSURE MODES OVER TRAPEZOIDAL SURFACES.
!
      DIMENSION XSC(4), YSC(4), XIS(JJ), TRM(2,3), P(JJ)
!
      IF ( IST .NE. 0 ) GO TO 100
!
      DO 50 J=1,JJ
      A1 = BETA*( YSC(2) - ETAS )/( XIS(J) - XSC(2) )
      A2 = BETA*( YSC(2) + ETAS )/( XIS(J) - XSC(2) )
      P1 = 0.0
      P2 = 0.0
!
      P(J) = 4.0/BETA
      IF ( A1 .GE. 1.0 ) GO TO 10
!
      P1 = 1.0 - 0.6366198* ASIN( SQRT( A1 ) )
!
   10 IF ( A2 .GE. 1.0 ) GO TO 50
!
      P2 = 1.0 - 0.6366198* ASIN( SQRT( A2 ) )
!
   50 P(J) = P(J)*( 1.0 - P1 - P2 )
!
      GO TO 1000
!
  100 IF ( TRM(1,1) .LE. 1.0 ) GO TO 200
!
      DO 150 J=1,JJ
      A1 = 2.0
      B1 = XIS(J) - XSC(1)
      IF ( B1 .LE. 0.0 ) GO TO 104
      A1 = BETA*( ETAS - YSC(1) )/B1
!
  104 A2 = 2.0
      B1 = XIS(J) - XSC(2)
      IF ( B1 .LE. 0.0 ) GO TO 106
      A2 = BETA*( YSC(2) - ETAS )/B1
!
  106 P1 = 0.0
      P2 = 0.0
      P(J) = TRM(1,2)
!
      IF ( A1 .GT. 1.0 ) GO TO 110
!
      P1 = TRM(2,2)*SQRT( 1.0 - A1*A1 )
!
  110 IF ( A2 .GT. 1.0 ) GO TO 150
!
      P2 = 1.0 - 0.6366198* ASIN( SQRT( A2 ) )
!
  150 P(J) = P(J)*( 1.0 - P1 - P2 )
!
      GO TO 1000
!
  200 B1 = 0.5/( TRM(1,1)*( 1.0 + TRM(1,1) ) )
!
      DO 250 J=1,JJ
!
      A1 = 2.0
      B2 = XIS(J) - XSC(1)
      IF ( B2 .LE. 0.0 ) GO TO 204
      A1 = BETA*( ETAS - YSC(1) )/B2
!
  204 A2 = 2.0
      B2 = XIS(J) - XSC(2)
      IF ( B2 .LE. 0.0 ) GO TO 206
      A2 = BETA*( YSC(2) - ETAS )/B2
!
  206 IF ( A2 .LT. 1.0 ) GO TO 210
!
      P(J) = TRM(1,2)/SQRT( 1.0 - ( A1/TRM(1,1) )**2 )
      GO TO 220
!
  210 A3 = BETA*( ETAS - YSC(1) )/( XSC(2) + (YSC(2)-ETAS)*BETA-XSC(1) )
!
      P(J) = ( TRM(1,2)/SQRT( 1.0 - ( A3/TRM(1,1) )**2 ) )              &
     &      *( 1.0 - SQRT( (1.0+A3)*(TRM(1,1)+A3)*B1 ) )
!
  220 IF ( ABS(TRM(2,1)) .GE. 1.0 .OR. IST .GT. 2 ) GO TO 250
      IF ( TRM(2,1) .LT. 0.0 ) GO TO 230
!
      B2 = XIS(J) - XSC(3)
      ETA2 = ETAS - YSC(3)
      GO TO 240
!
  230 B2 = XIS(J) - XSC(4)
      ETA2 = YSC(4) - ETAS
!
  240 IF ( B2 .LE. 0.0 .OR. ETA2 .EQ. 0.0 ) GO TO 250
!
      A3 = ABS(TRM(2,1))
      A3 = ( ETA2 - B2*A3/BETA )/( ETA2*( 1.0 - A3 ) )
      IF ( A3 .GT. 0.999999 ) GO TO 250
!
      P(J) = P(J)*0.6366198* ASIN( SQRT( A3 ) )
!
  250 END DO
!
 1000 RETURN
      END Subroutine SSfunc

      SUBROUTINE INVRT (C1,C2,N1,N2,N,NROWS,NTIME,IND)
      DIMENSION  C1(NROWS,NROWS), C2(NROWS), N1(NROWS), N2(NROWS)
!
!     SUBROUTINE TO PERFORM DOUBLE PRECISION INVERSION OF A MATRIX C
!     USING SIMPLE ELIMINATION AND FULL SEARCH BEFORE ELIMINATION.
      IND = 0
      IF (NTIME) 2,3,2
    2 CALL STATUS (N1)
      IT1 = N1(8)
    3 D1 = 0.0
      IF ( N .EQ. 1 ) GO TO 151
      II3 = N
      II2 = N - 1
      DO 20 J=1,N
      N1(J) = J
      N2(J) = J
      DO 10 I=1,N
      IF (D1 -  ABS(C1(I,J))) 5,10,10
    5 D1 =  ABS(C1(I,J))
      I1 = I
      J1 = J
   10 END DO
   20 END DO
      DO 150 K6=2,N
      IF (C1(I1,J1)) 50,30,50
   30 K5 = K6- 1
      GO TO 1000
   50 D1 = 1.0/C1(I1,J1)
      D2 = C1(I1,II3)
      D3 = C1(II3,J1)
      D4 = C1(II3,II3)
      DO 60 I=1,II2
      C2(I) = C1(I,J1)
      C1(I,J1)   = C1(I,II3)
      C1(I,II3)  = -C2(I)*D1
      D5         = -C1(I1,I)*D1
      C1(I1,I)   = C1(II3,I)
      C1(II3,I)  = D5
   60 END DO
      C2(I1) = D3
      C1(I1,J1)  = D4
      C1(II3,J1) = -D2*D1
      C1(I1,II3) = -D3*D1
      C1(II3,II3)= D1
      IF (II3-N) 70,110,110
   70 II4 = II3 + 1
      DO 80 I=II4,N
      C2(I) = C1(I,J1)
      C1(I,J1)  = C1(I,II3)
      C1(I,II3) = C2(I)
      D6   = C1(I1,I)*D1
      C1(I1,I)  = C1(II3,I)
      C1(II3,I) = D6
   80 END DO
      DO 100 J=II4,N
      DO 90  I=1,II2
      C1(I,J) = C1(I,J) - C2(I)*C1(II3,J)
   90 END DO
  100 END DO
  110 I = N1(I1)
      N1(I1) = N1(II3)
      N1(II3) = I
      I = N2(J1)
      N2(J1) = N2(II3)
      N2(II3) = I
      D7 = 0.0
      DO 140 J=1,II2
      DO 130 I=1,II2
      C1(I,J) = C1(I,J) + C2(I)*C1(II3,J)
      D8 =  ABS(C1(I,J))
      IF (D7 - D8) 120,130,130
  120 D7 = D8
      I1 = I
      J1 = J
  130 END DO
  140 END DO
      II3 = II3 - 1
      II2 = II2 - 1
  150 END DO
  151 IF (C1(1,1)) 160,155,160
  155 K5 = N
      GO TO 1000
  160 C1(1,1) = 1.0/C1(1,1)
      IF ( N .EQ. 1 ) GO TO 249
      DO 170 J=2,N
      C1(1,J) = C1(1,J)*C1(1,1)
  170 END DO
!     NOW THE FIRST ROW SOLUTION HAS BEEN OBTAINED.  THE REVERSE PRO-
!     CEDURE WILL BE STARTED.
!
      DO 210  K=2,N
      KM1 = K - 1
      DO 180 J=1,KM1
      C2(J) = C1(K,J)
      C1(K,J) = 0.0
  180 END DO
      DO 200 J=1,N
      DO 190 I=1,KM1
      C1(K,J) = C1(K,J) + C1(I,J)*C2(I)
  190 END DO
  200 END DO
!     THIS COMPLETES THE SOLUTION FOR THE K'TH ROW.
!
  210 END DO
!
!     NOW THE INVERSE HAS BEEN COMPUTED. THE NEXT STEPS WILL RE-ARRANGE
!     IT BACK INTO ITS ORIGINAL ORDER.
      DO 240 J=1,N
      J1 = N1(J)
      N1(J) = J
  215 IF (J1 - J) 220,240,220
  220 DO 230 I=1,N
      D1 = C1(I,J)
      C1(I,J) = C1(I,J1)
      C1(I,J1) = D1
  230 END DO
      K = N1(J1)
      N1(J1) = J1
      J1 = K
      GO TO 215
  240 END DO
      DO 248 I=1,N
      I1 = N2(I)
      N2(I) = I
  242 IF (I1 - I) 244,248,244
  244 DO 246 J=1,N
      D1 = C1(I,J)
      C1(I,J) = C1(I1,J)
      C1(I1,J) = D1
  246 END DO
      K = N2(I1)
      N2(I1) = I1
      I1 = K
      GO TO 242
  248 END DO
!     THE MATRIX HAS BEEN RE-ARRANGED IN ITS ORIGINAL ORDER.
!     IF NTIME N.E. ZERO, THEN THE INVERSION TIME IS PRINTED OUT.
  249 IF (NTIME) 250,270,250
  250 CALL STATUS (N1)
      TIME=(N1(8)-IT1)*0.01
      WRITE (6,260) N,TIME
  260 FORMAT (1H1,////,42H   THE TOTAL TIME FOR INVERTING THE MATRIX,// &
     &                 12H   OF ORDER ,I3,6H , IS ,E12.5,8H SECONDS)
      GO TO 270
 1000 WRITE (6,1001) K5
 1001 FORMAT (1H1, 4X,57H THE REDUCED MATRIX WAS FOUND TO BE SINGULAR ON&
     & ITERATION,I4 )
      IND = 1
  270 RETURN
      END Subroutine Invrt

!!!      OVERLAY(R2T,1,0)
      SUBROUTINE Omdlib()
!
      COMMON /COMM/ LEAVE, DUM1(11)

      CHARACTER(LEN=132):: title
      COMMON /FLTLE/ title   !!! TITLE(32)
!!!      REAL LI, INT
!!!      DATA LI/4HLIBR/, INT/4HINTQ/, FLT /4HFLTR/
!
      DUM1(11) = 0.0
!
      CALL FORMFD
      READ (5,'(A)') title(1:4)   !!!TITLE(1)
!!!   10 FORMAT (A4)
!
      IF ( title(1:4) .EQ. 'LIBR' ) GO TO 20
      IF ( title(1:4) .EQ. 'INTQ' ) LEAVE = -2
      IF ( title(1:4) .EQ. 'FLTR' ) LEAVE = -1
      BACKSPACE 5
      GO TO 100
!
   20 CALL RDLIB (LEAVE)
!
      IF ( LEAVE .EQ. 0 ) CALL TRNSMD()
!
      DUM1(11) = 1.0
!
  100 CONTINUE
      END Subroutine Omdlib


      SUBROUTINE RDLIB (LEAVE)
!
!     SUBROUTINE TO READ THE BY7 MODES AND THE DATA FOR SURFACE
!     SPLINE INTERPOLATION.
!
      COMMON /ALPVCT/ COEF(20,200), XF(200), YF(200), NFS(10), NFS1(10),&
     &   GMASS(20,20), NSURF, NMODES, DH, DW1, DW2, SPAN2(10),          &
     &   TAN12(10), TAN22(10), CAVG2(10), XR(10), YR(10)
!
      COMMON /MANE2/  H(200,20), IF(150,12), AMASS(200), NF, NUNIT,     &
     &   NMDTP, L1, L2, L3, L4, XMODE, XMASS, BREF, RHO, XIS(2,10),     &
     &   XOS(2,10), YCS(2,10), C(53,53), D(200), XFP(200), YFP(1120)

!!!      INTEGER KFM(8)
      CHARACTER(LEN=40):: kfm
!----------------------------------------------------------------------------
      CALL FREEFD
!
      READ (5,10) NF, NSURF, NMODES, NUNIT, NMDTP, L1, L2, L3, L4
      READ (5,11) XMODE, XMASS, DH, DW1, DW2, BREF, RHO
   10 FORMAT ( 6I10 )
   11 FORMAT ( 6F10.0 )
!
      IF ( NF .LT. 201 ) GO TO 20
!
      WRITE ( 6,15 ) NF
   15 FORMAT (1H1,"   THE STRUCTURE SIZE IS GREATER THAN 200, NF=",I10)
      GO TO 999
!
   20 IF ( NSURF .LT. 11 ) GO TO 30
!
      WRITE ( 6,25 ) NSURF
   25 FORMAT (1H1,"  THE NUMBER OF SURFACES EXCEED 10, NSURF=",I10 )
      GO TO 999
!
   30 IF ( NMODES .LT. 21 ) GO TO 40
!
      WRITE ( 6,35 ) NMODES
   35 FORMAT (1H1,"  THE NUMBER OF MODES EXCEED 20, NMODES=", I10 )
      GO TO 999
!
   40 IF ( NUNIT .EQ. 0 ) GO TO 50
      IF ( NUNIT .LT. 14 .AND. NUNIT .GT. 8 ) GO TO 50
!
      WRITE ( 6,45 ) NUNIT
   45 FORMAT (1H1,"  THE BY7 MODE TAPE UNIT IS NOT 9 TO 13, NUNIT=",I10)
      GO TO 999
!
   50 WRITE ( 6,55 ) NF, NSURF, NMODES
   55 FORMAT (1H1,16X," LIBRARY DATA FOR THE MODES AND MASSES ",//,     &
     &            16X,"    NUMBER OF STRUCTURAL POINTS = ",I3,//,       &
     &            16X,"    NUMBER OF LIFTING SURFACES  = ",I3,//,       &
     &            16X,"    NUMBER OF NATURAL MODES     = ",I3,//  )
!
      IF ( L1 .NE. 0 ) GO TO 70
      WRITE ( 6,60 )
   60 FORMAT (1H ,16X," THE MODES AND MASSES ARE OUTPUT FROM BY7", // )
      GO TO 80
!
   70 WRITE ( 6,75 )
   75 FORMAT (1H ,16X," THE MODES AND MASSES ARE OUTPUT FROM C28", // )
!
   80 WRITE ( 6,85 ) XMODE, XMASS, DH, DW1, DW2
   85 FORMAT (1H ,16X," THE SCALING CONSTANTS ARE AS FOLLOWS ", //,     &
     &            18X," XMODE = ",E14.7,     " XMASS = ",E14.7, //,     &
     &            18X," DH    = ",E14.7,     " DW1   = ",E14.7, //,     &
     &            18X," DW2   = ",E14.7, // )
!
      IF ( L4 .EQ. 0 ) GO TO 100
      WRITE ( 6,86 )
   86 FORMAT (1H ,10X,"  THE PITCH AND ROLL MODES WILL BE COMPUTED ",/  &
     &           ,10X,"  AS MODES (NMODES+1) AND (NMODES+2). ",//    )
!
  100 READ (5,'(A)') KFM
      READ (5,KFM)  ( XF(I), YF(I), I=1,NF )
!!!    1 FORMAT (8A10)
      WRITE( 6,110) ( XF(I), YF(I), I=1,NF )
  110 FORMAT (1H1,//,20X," THE STRUCTURAL POINT LOCATIONS",///, 4X,     &
     &  3 ("   XS(I)       YS(I)    " ), //,( 4X, 6E12.4 ) )
      WRITE ( 6,11 )
      DO 120 I=1,NSURF
      READ ( 5,10 ) NFF, ( IF(J,I), J=1,NFF )
      NFS(I) = NFF
      READ ( 5,11 )  XIS(1,I), XIS(2,I), YCS(1,I),                      &
     &               XOS(1,I), XOS(2,I), YCS(2,I)
  112 FORMAT (1H ,//, 5X," WING SURFACE " )
      WRITE ( 6,115 ) I, NFF, ( IF(J,I), J=1,NFF )
  115 FORMAT (1H ,//,10X," SURFACE ",I2," HAS ",I3," STRUCTURAL POINTS",&
     &            //,10X," THE POINTS ARE ",//, ( 5X,14I5 ) )
  120 WRITE ( 6,111 ) YCS(1,I), YCS(2,I), XIS(1,I), XOS(1,I),           &
     &                XIS(2,I), XOS(2,I)
  111 FORMAT (1H ,///, 17X,                                             &
     &   "       INBOARD                 OUTBOARD ",//,12X," Y     = ", &
     &   E11.4, 10X, E15.4, //,12X," X(LE) = ",E11.4, 10X,E15.4,        &
     &                      //,12X," X(TE) = ",E11.4, 10X,E15.4  )
!
!
      NF2 = NF
      IF ( L3 .EQ. 0 ) GO TO 170
      IF ( L3 .LT. 201  .AND. L3 .GT. 0 ) GO TO 140
      WRITE ( 6,130 ) L3
  130 FORMAT (1H1,"   THE NUMBER OF ZERO COORDINATES, L3 = ",I4,        &
     &         //,"   IS OUT OF BOUNDS. " )
      GO TO 999
  140 IF ( L3 .LE. NF ) GO TO 150
      WRITE ( 6,145 ) L3, NF
  145 FORMAT (1H1,"   THE NUMBER OF MODIFIED COORDINATES, L3= ",I4,     &
     &         //,"   EXCEEDS THE TOTAL NUMBER, NF= ",I4 )
      GO TO 999
!
  150 READ ( 5,10 ) NF2, ( IF(I,11), I=1,L3 )
      IF ( NF2 .LT. 201 .AND. NF2 .GT. 0 ) GO TO 160
      WRITE ( 6,155 ) NF2
  155 FORMAT (1H1,"   THE NUMBER OF INPUT MODE COORDINATES IS ",        &
     &         //,"   OUT OF BOUNDS, NF2 = ",I4 )
      GO TO 999
!
  160 WRITE ( 6,165 ) NF2, ( IF(I,11), I=1,L3 )
  165 FORMAT (1H ,//,10X," THE NUMBER OF INPUT MODE COORDINATES = ",I4, &
     &            //,10X,"   THE COORDINATES TO BE MODIFIED ARE ",      &
     &            //, (5X,14I5)  )
!
  170 READ (5,'(A)') KFM
      IF ( L2 .NE. 0 ) READ (5,KFM)( AMASS(I), I=1,NF2 )
      DO 180 N=1,NMODES
      IF ( L1 .EQ. 0 ) READ (5,11)
      READ (5,KFM) ZN
      READ (5,KFM)     ( H(I,N), I=1,NF2 )
      IF ( ZN .EQ. 0. ) ZN = 1.0
      RZN = 1.0/ZN
      DO 180 I=1,NF2
  180 H(I,N) = RZN*H(I,N)
!
      IF ( L2 .NE. 0 .AND. L3 .NE. 0 ) CALL ZEROS ( AMASS, H, IF(1,11), &
     &                                      L3, NF, NF2, NMODES )
!
      IF ( L2 .EQ. 0 .AND. L3 .NE. 0 ) CALL ZEROS ( H(1,1), H(1,2),     &
     &                                 IF(1,11), L3, NF, NF2, NMODES-1 )
!
      IF ( NF .LT. 0 ) GO TO 999
!
  200 IF ( L2 .EQ. 0 ) GO TO 205
!
      WRITE ( 6,210 ) ( AMASS(I), I=1,NF )
  210 FORMAT (1H1,26X,//," THE MASS DATA ",//, ( 2X, 5E15.7 ) )
  205 WRITE ( 6,215 )
  215 FORMAT (1H1,22X,//," THE MODE DEFLECTIONS" )
!
      DO 240 N=1,NMODES
  240 WRITE ( 6,220 ) N, ( H(I,N), I=1,NF )
  220 FORMAT (1H ,26X,//," MODE SHAPE ",I2,//, ( 2X, 5E15.7 ) )
      GO TO 1000
!
  999 LEAVE = 1
      GO TO 2000
!
 1000 IF ( L4 .EQ. 0 ) GO TO 2000
!
      N1 = NMODES + 1
      N2 = NMODES + 2
      DO 1010 I=1,NF
      H(I,N1) = XF(I)
 1010 H(I,N2) = YF(I)
!
 2000 RETURN
      END Subroutine RdLib

      SUBROUTINE ZEROS ( AM, H, IZ, L3, NF, NF2, NMODES )
!
      DIMENSION  AM(200), H(200,20), IZ(150)
!
      I1 = NF2 + 1
      I2 = 0
      I3 = 201
      NFK = NF + 1
      DO 50 I=1,NF
      I3 = I3 - 1
      NFK = NFK - 1
      DO 10 J=1,L3
      NZ = IABS( IZ(J) )
      IF ( NZ .NE. NFK ) GO TO 10
!
      IF ( IZ(J) ) 1, 2, 2
    1 IF ( I1 .EQ. 0 ) GO TO 999
      I1 = I1 - 1
      GO TO 30
    2 I2 = I2 + 1
      GO TO 30
!
   10 END DO
!
      IF ( I1 .EQ. 0 ) GO TO 999
      I1 = I1 - 1
      AM(I3) = AM(I1)
      DO 20 N=1,NMODES
   20 H(I3,N) = H(I1,N)
      GO TO 50
!
   30 AM(I3) = 1.0
      DO 40 N=1,NMODES
   40 H(I3,N) = 0.0
!
   50 END DO
!
      IF ( I1 .NE. 1 ) GO TO 999
!
      I3 = 200 - NF
      DO 100 I=1,NF
      I3 = I3 + 1
      AM(I) = AM(I3)
      DO 90 N=1,NMODES
   90 H(I,N) = H(I3,N)
  100 END DO
      GO TO 1000
!
  999 WRITE ( 6,101 )  NF2, NF, L3, I1, I2
  101 FORMAT (1H1,//," AN ERROR HAS OCCURED IN MODIFYING THE MODES ",   &
     &            //,"  TOTAL COORDINATES IN INPUT MODES = ", I4,       &
     &            //,"  TOTAL FINAL COORDINATES DESIRED  = ", I4,       &
     &            //,"  TOTAL COORDINATES TO BE MODIFIED = ", I4,       &
     &            //,"  NUMBER OF INPUT MODE COORDINATES   ", I4,       &
     &             /,"                     NOT USED      = ", I4,       &
     &            //,"  NUMBER OF COORDINATES SET TO ZERO= ", I4  )
!
      NF = -1
!
 1000 RETURN
      END Subroutine Zeros

      SUBROUTINE ARRANG ( A, X, NFS, NSURF, IZ )
      DIMENSION  A(200), X(200), NFS(10), IZ(150,10)
      NF2 = 0
      DO 10 IS=1,NSURF
      NF0 = NF2
      NF1 = NF2 + 1
      NF2 = NF2 + NFS(IS)
      DO 10 I=NF1,NF2
      I2 = IZ(I-NF0,IS)
   10 X(I) = A(I2)
      DO 20 I=1,NF2
   20 A(I) = X(I)
      RETURN
      END Subroutine Arrang

      SUBROUTINE TRNSMD
!
!     SUBROUTINE TO TRANSFORM THE MODES INTO SURFACE SPLINE COEFFICIENTS
!     AND TO REDUCE THE SET TO MINIMUM SIZE FOR FITTING INTO 'ALPVCT'.
!
      COMMON /ALPVCT/ COEF(20,200), XF(200), YF(200), NFS(10), NFS1(10),&
     &   GMASS(20,20), NSURF, NMODES, DH, DW1, DW2, SPAN2(10),          &
     &   TAN12(10), TAN22(10), CAVG2(10), XR(10), YR(10)
!
      COMMON /MANE2/  H(200,20), IF(150,12), AMASS(200), NF, NUNIT,     &
     &   NMDTP, L1, L2, L3, L4, XMODE, XMASS, BREF, RHO, XIS(2,10),     &
     &   XOS(2,10), YCS(2,10), E(53,53), D(200), XFP(200), YFP(1120)
      COMMON /MATR/ C(103,103)
      IF ( L2  .EQ.  0  ) GO TO 100
      IF ( RHO .EQ. 0.0 ) RHO = 0.0023769
!
! ***** SEA LEVEL STANDARD DENSITY IS THE DEFAULT VALUE. ****
!
      XMASS = XMASS*432.0/(32.2*RHO*(BREF**3))
!
! ***** RHO IS SLUGS/FT**3 AND BREF IS INCHES. ****
!
      DO 10 J=1,NF
      D(J) = AMASS(J)*XMASS
      DO 10 I=1,NMODES
   10 H(J,I) = H(J,I)*XMODE
!
!     THE GENERALIZED MASSES ARE NOW COMPUTED.
!
      DO 30 KS=1,NMODES
      DO 30 KR=1,NMODES
      GMASS(KR,KS) = 0.0
      DO 20 I=1,NF
   20 GMASS(KR,KS) = GMASS(KR,KS) + D(I)*H(I,KR)*H(I,KS)
   30 CONTINUE
!
      WRITE (6,40) BREF, RHO
   40 FORMAT (1H1,//,23X," THE GENERALIZED MASSES",//,                  &
     &               20X," REFERENCE LENGTH = ",E12.5,//,               &
     &               20X," DENSITY          = ",E12.5,//  )
!
      WRITE (6,50) (( GMASS(I,J), I=1,NMODES ), J=1,NMODES )
   50 FORMAT ( 5E15.7 )
!
!
  100 NF2 = 0
!
!     THE STRUCTURAL POINTS ARE TRANSFORMED INTO THE SQUARE PLANE
!     FOR EACH STRUCTURAL SURFACE.
!
      DO 600 IS=1,NSURF
!
      SPAN = YCS(2,IS) - YCS(1,IS)
      TAN1 = ( XOS(1,IS) - XIS(1,IS) )*0.5
      TAN2 = ( XOS(2,IS) - XIS(2,IS) )*0.5
      CAVG = ( XIS(2,IS) - XIS(1,IS) + XOS(2,IS) - XOS(1,IS) )*0.5
!
      SPAN2(IS) = SPAN
      TAN12(IS) = TAN1
      TAN22(IS) = TAN2
      CAVG2(IS) = CAVG
      XR(IS)    = ( XIS(1,IS) + XOS(1,IS) )*0.5
      YR(IS)    = YCS(1,IS)
!
      NF0 = NF2
      NF1 = NF2 + 1
      NF2 = NF2 + NFS(IS)
!
      DO 520 I=NF1,NF2
      IF2 = IF(I-NF0,IS)
      YFP(I) = 2.0*( YF(IF2) - YCS(1,IS) )/SPAN - 1.0
      CETA   = CAVG - ( TAN1 - TAN2 )*YFP(I)
      XLE    = XR(IS) + YFP(I)*TAN1
  520 XFP(I) = 2.0*( XF(IF2) - XLE )/CETA - 1.0
!
  600 END DO
!
      WRITE (6,530)
  530 FORMAT (1H1)
!
      NF2 = 0
      DO 700 IS=1,NSURF
      WRITE (6,535) IS
  535 FORMAT (1H ,//,16X," TRANSFORMED STRUCTURAL POINTS FOR SURFACE ", &
     &  I2 )
      NF1 = NF2 + 1
      NF2 = NF2 + NFS(IS)
      DO 550 I=NF1,NF2
      XF(I) = XFP(I)
  550 YF(I) = YFP(I)
!
      WRITE (6,560) ( XF(I), YF(I), I=NF1,NF2 )
  560 FORMAT (1H ,//,2(10X,"   XF          YF     "),//,(2(8X,2E12.4)) )
!
  700 END DO
!
      IF ( L2 .NE. 0 ) CALL ARRANG ( AMASS, XFP, NFS, NSURF, IF )
!
      NM2 = NMODES
      IF ( L4 .NE. 0 ) NM2 = NM2 + 2
      DO 710 J=1,NM2
  710 CALL ARRANG ( H(1,J), XFP, NFS, NSURF, IF )
!
      IR = 0
      NF2 = 0
      DO 750 IS=1,NSURF
      NF1 = NF2 + 1
      NF2 = NF2 + NFS(IS)
      NFS1(IS) = NF1
!
      CALL SPLINE( XF(NF1), YF(NF1), NFS(IS), C,D, XFP, YFP, 103, IND )
!
      LFS3 = NFS(IS) + 3
      DO 740 I=1,LFS3
      IR = IR + 1
!
      DO 730 J=1,NMODES
      D(J) = 0.0
      K = 1
      DO 730 K2=NF1,NF2
      D(J) = D(J) + H(K2,J)*C(K,I)
      K = K + 1
  730 CONTINUE
!
      DO 740 J=1,NMODES
  740 COEF(J,IR) = D(J)
      NF3 = IR - LFS3 + 1
!
  750 END DO
      RETURN
      END Subroutine Trnsmd

!!!      OVERLAY(R2T,2,0)
      SUBROUTINE Omanrd()
      CALL MAINRD
      END Subroutine OmanRd


      SUBROUTINE MAINRD
!
!     SUBROUTINE TO CALL APPROPRIATE READ SUBROUTINE ACCORDING TO OP1.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), HMACH
      COMMON /OPTION/ IOP1, IOP2, IOP3, IOP4, NUNIT, IQT, IOPLU, IND

      CHARACTER(LEN=128):: tle
!!!      COMMON /TITLE/ VEHC(16), SURF(16)

      COMMON /COMM/LEAVE,DUM1(11)

!!!      EQUIVALENCE ( NAME, OREDPR, OGEOTY )
!!!      DATA NAME /6HREADTY/, ENDE /4HEND /


      READ (5,'(A)') tle(1:64) !!!( VEHC(I), I=1,16 )
   10 FORMAT ( 16A4 )
      READ (5,'(A)') tle(65:128) !!! ( SURF(I), I=1,16 )
!
      WRITE (6,'(///4X,A)') tle(1:64)   !!! ( VEHC(I), I=1,16)
   20 FORMAT (1H1, ///, 4X, 16A4 )
      WRITE (6,'(/4X,A)') tle(65:128)   !!! ( SURF(I), I=1,16)
   21 FORMAT ( / , 4X, 16A4, // )
!
      CALL FREEFD()
      READ (5,9001) IOP1, IOP2, IOP3, IOP4, IQT, IOPLU
 9001 FORMAT ( 6I10 )
!
      IF ( IOP1 .GT. 0 ) GO TO 40
!
      WRITE (6,35)
   35 FORMAT (54H     OPTION 1   THE AERODYNAMIC MATRICES WILL BE       &
     & ,//,   45H                COMPUTED FOR ALL FREQUENCIES. , //    )
      GO TO 50
!
   40 WRITE (6,45)  IOP1
   45 FORMAT (54H     OPTION 1   THE AERODYNAMIC MATRICES WILL BE       &
     & ,//,   31H                READ FROM UNIT , I2, // )
!
   50 IF ( IOP2 .GE. 0 ) GO TO 60
!
      WRITE (6,55)
   55 FORMAT (54H     OPTION 2   THE AERODYNAMIC PRESSURE DISTRIBUTIONS &
     & ,//,   54H                AND INTEGRATED AERODYNAMIC CHARACTER-  &
     & ,//,   54H                ISTICS WILL BE CALCULATED FOR ALL      &
     & ,//,   28H                FREQUENCIES. , // )
      GO TO 92
!
   60 IF ( IOP2 .NE. 0 ) GO TO 70
!
      WRITE (6,65)
   65 FORMAT (54H     OPTION 2   THE GENERALIZED FORCES WILL BE CALCU-  &
     & ,//,   54H                LATED FOR ALL FREQUENCIES AS A SET OF  &
     & ,//,   39H                NMODESXNMODES MATRICES. , // )
      GO TO 80
!
   70 WRITE (6,75) IOP2
   75 FORMAT (54H     OPTION 2   THE GENERALIZED FORCES WILL BE CALCU-  &
     & ,//,   45H                LATED FOR ALL FREQUENCIES AS  ,I2,     &
     & ,//,   54H                SETS OF VECTORS WITH NMODES ELEMENTS.  &
     & ,// )
!
   80 IF ( IQT .NE. 0 ) WRITE (6,85) IQT
   85 FORMAT (54H     IQT        A Q-TERM MATRIX TAPE FOR THE INPUT     &
     & ,//,   52H                FREQUENCIES WILL BE PRODUCED ON UNIT,  &
     &  I3, // )
!
      IF ( IOP3 .NE. 0 ) WRITE (6,90) IOP3
   90 FORMAT (54H     OPTION 3   THE GENERALIZED FORCES WILL BE INTER-  &
     & ,//,   54H                POLATED AND THE EXPANDED SET WRITTEN   &
     & ,//,   24H                ON UNIT , I2, // )
!
      IF ( IOP4 .EQ. 0 ) WRITE (6,91)
   91 FORMAT (54H     OPTION 4   A FLUTTER SOLUTION WILL BE OBTAINED    &
     &        ,// )
!
   92 IF ( IOP1 .EQ. 0 ) GO TO 500
      NUNIT = IABS( IOP1 )
      REWIND NUNIT
!
      WRITE (6,95)
   95 FORMAT (54H     OPTIONS 1 AND LU                                  &
     & ,//,   51H                THE AERODYNAMIC MATRICES ARE IN THE,/ )
!
      IF ( IOPLU .NE. 0 ) GO TO 120
!
      IF ( IOP1 ) 100, 500, 110
!
  100 WRITE (6,105) NUNIT
  105 FORMAT (50H                INVERTED FORM ON OUTPUT TAPE UNIT ,    &
     &  I2, // )
      GO TO 500
!
  110 WRITE (6,115) NUNIT
  115 FORMAT (49H                INVERTED FORM ON INPUT TAPE UNIT  ,    &
     &  I2, // )
      GO TO 500
!
  120 IF ( IOP1 ) 130, 500, 140
!
  130 WRITE (6,135) NUNIT
  135 FORMAT (52H                UNINVERTED FORM ON OUTPUT TAPE UNIT ,  &
     &  I2, // )
      GO TO 500
!
  140 WRITE (6,145) NUNIT
  145 FORMAT (51H                UNINVERTED FORM ON INPUT TAPE UNIT ,   &
     &  I2, // )
!
!!!  500 CALL OVERLAY ( OREDPR,2,1 )
  500 CALL Oredpr()

      CALL Ogeoty()
!
 1000 CONTINUE
      RETURN
      END Subroutine MainRd


!!!      OVERLAY(READTY,2,1)
      SUBROUTINE Oredpr()
      CALL READPR
      END Subroutine Oredpr


      SUBROUTINE READPR
!
!     SUBROUTINE TO READ THE PROBLEM INPUT DATA.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), HMACH
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10),XIBAR(15,10), XIS(15,31,10), YBAR(15,10),          &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      DIMENSION XLE(63), YLE(63), XTE(63), YTE(63), C(63,63), D(63),    &
     &          GMACH(11,16)
      COMMON /COMM/LEAVE,DUM1(11)
      COMMON /GAMA/ GAMMA
      COMPLEX ALPRS   !!! never referenced ???
!
      EQUIVALENCE (XLE(1),ETAS(1,1)), (YLE(1),ETAS(1,4)),               &
     &            (XTE(1),ETAS(1,7)), (YTE(1),ETABAR(1,1)),             &
     &            (D(1),ETABAR(1,4)), (GMACH(1,1),ETABAR(1,7)),         &
     &            (C(1,1),XIS(1,1,1))
!
      CALL FREEFD
!
      READ (5,1) NSURF, LS, NALP, NK, NW, IDUMP, IDNWSH, JSUROP
      READ (5,3) XMACH, GAMMA, BREF, ( REFREQ(I), I=1,NK )
!
      IF ( GAMMA .EQ. 0.0 ) GAMMA = 1.4
!
    1 FORMAT ( 6I10 )
    3 FORMAT ( 6F10.0 )
!
      DO 2 I=1,NK
    2 RK(I) = REFREQ(I)
!
      BETA2 = ABS( 1.0 - XMACH**2 )
      BETA  = SQRT(BETA2)
      PI = 3.141593
!
      IF ( LS .EQ. 0 ) GO TO 20
!
      WRITE (6,10)
   10 FORMAT (1H1,///,17X,36H ANTI-SYMMETRIC AERODYNAMIC SOLUTION ,// )
      GO TO 40
   20 WRITE (6,30)
   30 FORMAT (1H1,///,17X,36H   SYMMETRIC AERODYNAMIC SOLUTION    ,// )
!
      IF ( JSUROP .NE. 0 ) WRITE (6,35)
   35 FORMAT (17X,36H WITH SUPERSONIC WEIGHTING FUNCTION  ,// )
   40 BRINV = 1.0/BREF
!
      WRITE (6,50) NSURF, NALP, NW, XMACH, GAMMA, BREF, NK,             &
     &             ( REFREQ(I), I=1,NK )
   50 FORMAT (15X,39H TOTAL NUMBER OF SURFACES            = ,I3   ,// , &
     &        15X,39H TOTAL SETS OF ALPHA VECTORS         = ,I3   ,// , &
     &        15X,39H TOTAL NUMBER OF CONTROL POINTS      = ,I3   ,// , &
     &        15X,39H MACH NUMBER                         = ,E11.4,// , &
     &        15X,39H GAMMA                               = ,E11.4,// , &
     &        15X,39H REFERENCE LENGTH                    = ,E11.4,// , &
     &        15X,39H TOTAL NUMBER OF REDUCED FREQUENCIES = ,I3   ,///, &
     &        15X,39H    THE REDUCED FREQUENCIES ARE              ,// , &
     &                 (3X, 5E15.7 )  )
!
      DO 500  IS=1,NSURF
!
      WRITE (6,48) IS
   48 FORMAT (1H1,20X,30H DATA FOR LIFTING SURFACE NO.  ,I2,///,        &
     &            22X,30H       GEOMETRIC DATA          )
!
      READ (5,1) NC(IS), NS(IS), NJ(IS), NSI(IS), ITRANS(IS), NLE, NTE, &
     & ICHORD(IS), IXI(IS), LSPAN(IS), LSYM(IS), ISTYPE(IS),            &
     & ( KSURF(IS,KS), KS=1,NSURF )
      READ (5,3) ALPO(IS), XV(IS), YV(IS), ZV(IS), THETA(IS)
!
      IST = ISTYPE(IS)
!
      IF ( ISTYPE(IS) .GT. 1 ) ISTYPE(IS) = ISTYPE(IS) + 1
!
      THETA(IS) = THETA(IS)*0.0174533
!
!     THETA IS CONVERTED TO RADIANS.
!
      IF ( JSUROP .EQ. 0 ) GO TO 65
      IF ( ISTYPE(IS) .NE. 1 ) GO TO 65
!
      WRITE (6,60)
   60 FORMAT (5( /,10X,20(2H ") ),//,                                   &
     &  59H   ISTYPE CANNOT BE (1) FOR SUPERSONIC WEIGHTING FUNCTION, , &
     &/,40H   IT MUST BE (2) OR (3) ---- TRY AGAIN.  )
!
      LEAVE = 1
      GO TO 2000
!
   65 NW2(IS) = NC(IS)*NS(IS)
      IF ( NJ(IS) .NE. 0 ) GO TO 51
!
      IF ( XMACH .GE. 1.0 ) GO TO 66
!
      NJ(IS) = NC(IS)
      IF ( IXI(IS) .EQ. 2 )  NJ(IS) = NC(IS) + 1
      IF ( IXI(IS) .EQ. 3 )  NJ(IS) = NC(IS) - 1
      GO TO 51
!
   66 NJ(IS) = 2*NC(IS)
!
   51 IF ( NSI(IS) .NE. 0 ) GO TO 52
!
      NSI(IS) = 2*NS(IS) + 1
      IF ( ISTYPE(IS) .EQ. 3 )  NSI(IS) = NSI(IS) + 1
      IF ( ISTYPE(IS) .EQ. 4 )  NSI(IS) = NS(IS)  + 1
      IF ( ISTYPE(IS) .GT. 1 )  NSI(IS) = 2*NSI(IS)
!
   52 READ (5,3) ( XLE(I), YLE(I), I=1,NLE )
      READ (5,3) ( XTE(I), YTE(I), I=1,NTE )
!
      WRITE (6,216) ( XLE(I), YLE(I), I=1,NLE )
      WRITE (6,217) ( XTE(I), YTE(I), I=1,NTE )
!
      CALL BreakPoint ( BRPT(1,1,IS), NBRPT(IS), XLE, YLE, XTE, YTE,         &
     &             NLE, NTE, LEAVE )
!
      IF ( LEAVE .NE. 0 ) GO TO 2000
!
      BO(IS) = 0.5*( BRPT(1,2,IS) - BRPT(1,1,IS) )
      I = NBRPT(IS)
      SO(IS) = BRPT(2,I,IS) - BRPT(2,1,IS)
!
      NB2      = NBRPT(IS)
      NEND(IS) = NBRPT(IS)/2 - 1
      BOINV(IS)= 1.0/BO(IS)
!
      BO2 = BO(IS) + BRPT(1,1,IS)
      Y   = BRPT(2,1,IS)
      DO 70 IX2=1,NB2
      BRPT(1,IX2,IS) = BRPT(1,IX2,IS) - BO2
   70 BRPT(2,IX2,IS) = BRPT(2,IX2,IS) - Y
!
      XV(IS) = ( XV(IS) + BO(IS) )*BRINV
      YV(IS) = YV(IS)*BRINV
      ZV(IS) = ZV(IS)*BRINV
!
!     CALCULATION OF THE ASPECT RATIO AND AREA OF THE SEMI-SPAN.
!
      AREA2 = 0.0
      I = 1
   80 I = I + 2
      AREA2 = AREA2 + ( BRPT(2,I,IS) - BRPT(2,I-2,IS) )                 &
     &           *0.5*( BRPT(1,I-1,IS) - BRPT(1,I-2,IS)                 &
     &                + BRPT(1,I+1,IS) - BRPT(1,I,IS) )
!
      IF ( I .LT. NB2-2 ) GO TO 80
!
      AR(IS) = 2.0*( SO(IS)**2 )/AREA2
      AREA(IS) = AREA2
!
      IF ( DUM1(11) .EQ. 0.0 ) GO TO 89
      READ (5,1) NST(IS)
      WRITE (6,85) IS, NST(IS)
   85 FORMAT (//,"  AERO SURFACE ", I2," WILL USE SLOPES AND ",         &
     &        //,"  DEFLECTIONS FROM STRUCTURAL SURFACE ",I2  )
!
   89 IF ( ITRANS(IS) .EQ. 0 ) GO TO 200
!
!     THE MACH NUMBER DISTRIBUTION IS NOW READ.
!
      NCP1 = NC(IS) + 1
      NSP1 = NS(IS) + 1
!
      IF ( ITRANS(IS) .LT. 0 ) GO TO 100
!
      READ (5,3) ( GMACH(I,1), I=1,NCP1 )
!
      DO 90 J=2,NSP1
      DO 90 I=1,NCP1
   90 GMACH(I,J) = GMACH(I,1)
      GO TO 110
!
  100 READ (5,3) ( ( GMACH(I,J), I=1,NCP1 ), J=1,NSP1 )
!
!
  110 DELT = 2.0/NC(IS)
      XTE(1) = -1.0
      NCP = NC(IS)
      NSP = NS(IS)
      DO 112 I=2,NCP
  112 XTE(I) = XTE(I-1) + DELT
      XTE(NCP1) = 1.0
!
      DELT = 2.0/NS(IS)
      YTE(1) = -1.0
      DO 114 J=2,NSP
  114 YTE(J) = YTE(J-1) + DELT
      YTE(NSP1) = 1.0
!
      K = 0
      DO 120 J=1,NSP1
      DO 120 I=1,NCP1
      K = K + 1
      XLE(K) = XTE(I)
  120 YLE(K) = YTE(J)
!
      N2 = K
!
      CALL SPLINE ( XLE, YLE, N2, C, D, IA, IB, 63, LEAVE )
      IF ( LEAVE .NE. 0 ) GO TO 2000
!
      N3 = N2 + 3
      DO 140 K=1,N3
      D(1) = 0.0
      L = 0
      DO 130 J=1,NSP1
      DO 130 I=1,NCP1
      L = L + 1
      BMACH(L,IS) = GMACH(I,J)
  130 D(1) = D(1) + C(L,K)*GMACH(I,J)
  140 C(1,K) = D(1)
!
      DELT = PI/(2.0*NCP)
      DO 145 I=1,NCP
  145 XTE(I) = -COS( DELT*(2*I-1) )
      DELT = PI/(2.0*NSP)
      DO 150 J=1,NSP
  150 YTE(J) = -COS( DELT*(2*J-1) )
!
      DO 160 J=1,NSP
      DO 160 I=1,NCP
!
  160 CALL MAGN ( XTE(I), YTE(J), N2, 1, XLE, YLE, GMACH(I,J), D, C, 63)
!
      DO 165 I=1,NCP
      XBAR(1,I) = 1.0
      XBAR(2,I) = XTE(I)
      DO 165 K=3,NCP
  165 XBAR(K,I) = 2*XTE(I)*XBAR(K-1,I) - XBAR(K-2,I)
!
      DO 170 J=1,NSP
      YBAR(1,J) = 1.0
      YBAR(2,J) = YTE(J)
      DO 170 K=3,NSP
  170 YBAR(K,J) = 2*YTE(J)*YBAR(K-1,J) - YBAR(K-2,J)
!
      DELT = 2.0/NCP
      DO 185 J=1,NSP
      DO 175 K=1,NCP
      D(K) = 0.0
      DO 175 I=1,NCP
  175 D(K) = D(K) + XBAR(K,I)*GMACH(I,J)
      DO 180 K=2,NCP
  180 GMACH(K,J) = D(K)*DELT
  185 GMACH(1,J) = D(1)/NCP
!
      DELT = 2.0/NSP
      DO 198 I=1,NCP
      DO 190 K=1,NSP
      D(K) = 0.0
      DO 190 J=1,NSP
  190 D(K) = D(K) + YBAR(K,J)*GMACH(I,J)
      DO 195 K=2,NSP
  195 GMACH(I,K) = D(K)*DELT
  198 GMACH(I,1) = D(1)/NSP
!
  200 WRITE (6,210) NC(IS), NS(IS), NJ(IS), NSI(IS), SO(IS), BO(IS),    &
     & AREA(IS), AR(IS)
!
  210 FORMAT( //, 22X,30H      AERODYNAMIC DATA         , //, 5X,       &
     &  "   NC = ",I2,"   NS = ",I2,"   NJ = ",I2,"   NSI = ",I2,"   SO &
     &= ",E11.4,//,5X,"   BO = ",E11.4,   "  AREA = ",E11.4,"    AR = ",&
     &  E11.4    )
!
      WRITE (6,215) XV(IS), YV(IS), ZV(IS), THETA(IS), ICHORD(IS),      &
     & IXI(IS),IST,        LSYM(IS), LSPAN(IS)
!
  215 FORMAT (/,5X,"   XV = ",E11.4,"    YV = ",E11.4,"    ZV = ",E11.4,&
     & //,      5X,"THETA = ",E11.4," ICORD = ",I2,9X,"   IXI = ",I2,9X,&
     & //,      5X,"ISTYP = ",I2,9X,"  LSYM = ",I2,9X," LSPAN = ",I2 )
!
  216 FORMAT (//,24X," LEADING EDGE POINTS", //,                        &
     &           24X,"XLE               YLE",//, (20X,E12.5,6X,E12.5) )
  217 FORMAT (//,24X," TRAILING EDGE POINTS",//,                        &
     &           24X,"XTE               YTE",//, (20X,E12.5,6X,E12.5) )
!
      WRITE (6,218) IS
  218 FORMAT (1H1,20X,30H DATA FOR LIFTING SURFACE NO.  ,I2,/ )
!
      IF ( ITRANS(IS) .EQ. 0 ) GO TO 400
!
      WRITE (6,220) ( YLE(J), J=1,N2,NCP1 )
  220 FORMAT (//,22X," MACH NUMBER DISTRIBUTION",//,9X,"YBAR",/,        &
     &  (7X,16F7.3),/,"  XBAR" )
      WRITE (6,225)
  225 FORMAT ( "  XBAR" )
!
      DO 230 I=1,NCP1
  230 WRITE (6,240) XLE(I), ( BMACH(J,IS), J=I,N2,NCP1 )
  240 FORMAT ( (17F7.3,/) )
!
      K = 0
      DO 250 I=1,NCP
      DO 250 J=1,NSP
      K = K + 1
  250 BMACH(K,IS) = GMACH(I,J)
!
      WRITE (6,260) ( YLE(J), J=1,N2,NCP1 )
  260 FORMAT (//,15X," RECALCULATED MACH NUMBER DISTRIBUTION ",//,      &
     &  9X, "YBAR",/, ( 7X,16F7.3 ) )
      WRITE (6,225)
!
      DO 280 I=1,NCP1
      K = 0
      DO 270 J=1,N2,NCP1
      K = K + 1
  270 CALL CALCM ( XLE(I), YLE(J), NCP, NSP, BMACH(1,IS),               &
     &             GMACH(K,1), IA, IB )
!
  280 WRITE (6,240) XLE(I), ( GMACH(K,1), K=1,NSP1 )
!
  400 CALL PROPT ( KSURF, NSURF, IS, ITRANS(IS), LEAVE )
!
  500 END DO
      LEAVE = 0
!
 2000 RETURN
      END Subroutine ReadPr

      SUBROUTINE PROPT ( KSURF, NSURF, IS, ITRANS, LEAVE )
!
!     SUBROUTINE TO PRINT THE AERODYNAMIC MATRIX OPTION DATA.
!
      DIMENSION KSURF(10,10)
!
      WRITE (6,10) IS, ITRANS
   10 FORMAT (///,23X," AERODYNAMIC MATRIX DATA ",//,                   &
     &           ,23X,"     ITRANS(",I2,") = ",I2     )
!
      WRITE (6,15) IS, IS, KSURF(IS,IS)
   15 FORMAT (//,"        KSURF(",I2,",",I2,") = ",I3  )
!
      IF ( KSURF(IS,IS) ) 20, 30, 40
!
   20 WRITE (6,25) IS
   25 FORMAT (/, "     SURFACE ",I2," IS TRANSONIC WITH A SHOCK AT THE",&
     &        /, "     LEADING EDGE AND ALL SUBSONIC FLOW. " )
      GO TO 50
!
   30 WRITE (6,35) IS
   35 FORMAT (/, "     SURFACE ",I2," IS NOT A TRANSONIC SURFACE." )
      GO TO 50
!
   40 WRITE (6,45) IS
   45 FORMAT (/, "     SURFACE ",I2," IS TRANSONIC WITH A SHOCK AT THE",&
     &        /, "     TRAILING EDGE AND ALL SUPERSONIC OR SUBSONIC TO",&
     &        /, "     SUPERSONIC ACCELERATING FLOW. " )
!
   50 IF ( NSURF .EQ. 1 ) GO TO 200
      DO 100 KS=1,NSURF
      IF ( IS .EQ. KS ) GO TO 100
      WRITE (6,15) IS, KS, KSURF(IS,KS)
      IF ( KSURF(IS,KS) ) 60, 70, 80
!
   60 WRITE (6,65) KS, IS
   65 FORMAT (/, "     SURFACE ",I2," WILL NOT HAVE ANY INTERFERENCE ", &
     &        /, "     EFFECTS ON SURFACE ",I2,"." )
      GO TO 100
!
   70 WRITE (6,75) KS, IS
   75 FORMAT (/, "     SURFACE ",I2," WILL HAVE INTERFERENCE EFFECTS ", &
     &        /, "     ON SURFACE ",I2,", BUT CORRECTION TERMS WILL ",  &
     &        /, "     NOT BE CALCULATED." )
      GO TO 100
!
   80 WRITE (6,85) KS, IS
   85 FORMAT (/, "     SURFACE ",I2," WILL HAVE INTERFERENCE EFFECTS ", &
     &        /, "     ON SURFACE ",I2,", AND CORRECTION TERMS WILL ",  &
     &        /, "     BE CALCULATED." )
!
  100 END DO
!
  200 RETURN
      END Subroutine Propt

      SUBROUTINE BreakPoint( BRPT, NBRPT, XLE, YLE, XTE, YTE, NLE, NTE, IND )
!
!     SUBROUTINE TO CONVERT LEADING AND TRAILING EDGE BREAKPOINT
!     ARRAYS INTO BREAKCHORD ARRAYS.
!
      DIMENSION BRPT(2,40), XLE(20), YLE(20), XTE(20), YTE(20)
!
      IF ( YLE(NLE) .EQ. YTE(NTE) ) GO TO 70
!
      WRITE (6,65) YLE(NLE), YTE(NTE)
   65 FORMAT (1H1," THE LEADING AND TRAILING EDGE END COORDINATES", //, &
     &            " DO NOT MATCH, YLE(NLE)=",E14.7,"  YTE(NTE)=",E14.7 )
      GO TO 999
!
   70 IF ( NLE .LT. 21 .AND. NTE .LT. 21 ) GO TO 80
      WRITE (6,75) NLE, NTE
   75 FORMAT (1H1," THE PLANFORM BOUNDARY POINTS EXCEED 20 ", //,       &
     &            " NLE = ", I10,"        NTE = ", I10 )
      GO TO 999
!
   80 BRPT(1,1) = XLE(1)
      BRPT(2,1) = YLE(1)
      BRPT(1,2) = XTE(1)
      BRPT(2,2) = YTE(1)
      L = 2
      I = 2
      K = 3
!
   85 IF ( I .GT. NLE .OR. L .GT. NTE ) GO TO 110
      IF ( YLE(I) - YTE(L) ) 87,90,100
!
   87 BRPT(1,K) = XLE(I)
      BRPT(2,K) = YLE(I)
      BRPT(1,K+1) = XTE(L-1) + (XTE(L)-XTE(L-1))*(YLE(I)-YTE(L-1))      &
     &                                          /(YTE(L)-YTE(L-1))
      BRPT(2,K+1) = YLE(I)
      I = I+1
      K = K+2
      GO TO 85
!
   90 BRPT(1,K) = XLE(I)
      BRPT(2,K) = YLE(I)
      BRPT(1,K+1) = XTE(L)
      BRPT(2,K+1) = YTE(L)
      I = I+1
      L = L+1
      K = K+2
      GO TO 85
!
  100 BRPT(1,K) = XLE(I-1)+(XLE(I)-XLE(I-1))*(YTE(L)-YLE(I-1))          &
     &                                      /(YLE(I)-YLE(I-1))
      BRPT(2,K) = YTE(L)
      BRPT(1,K+1) = XTE(L)
      BRPT(2,K+1) = YTE(L)
      L = L+1
      K = K+2
      GO TO 85
!
  110 NBRPT = K-1
      IND = 0
      GO TO 1000
!
  999 IND = 1
 1000 RETURN
      END Subroutine BreakPoint

!!!      OVERLAY(READTY,2,2)
      SUBROUTINE Ogeoty()
      CALL GEOMTY
      END Subroutine Ogeoty


      SUBROUTINE GEOMTY
!
!     SUBROUTINE TO COMPUTE THE SUBSONIC OR SUPERSONIC GEOMETRY DATA
!     ACCORDING TO THE VALUE OF QMACH(1).
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), HMACH
!
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10),XIBAR(15,10), XIS(15,31,10), YBAR(15,10),          &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
      COMMON /OPTION/ IOP1, IOP2, IOP3, IOP4, NUNIT, IQT, IOPLU, IND
      COMMON /COMM/ LEAVE, DUM1(11)
!
      DO 1000 IS=1,NSURF
!
      J1 = 1
      MBAR = NC(IS)
      I  = IXI(IS)
      GO TO ( 10, 20, 30, 50 ),I
!
   10 B1 = PI/(2*MBAR+1)
!     REGULAR HSU POINTS ARE COMPUTED.
      JJ = 0
      GO TO 40
!
   20 B1 = PI/(2*MBAR+2)
      JJ = 0
!     TSCHEBYCHEV ROOTS CORRESPONDING TO FIRST KIND POLY. USED.
      GO TO 40
!
   30 B1 = PI/(2*MBAR-1)
      XBAR(1,IS) = -0.999
      JJ = 2
!     HSU POINTS USED WITH A LEADING EDGE POINT.
      J1 = 2
!
   40 DO 45 JI=J1,MBAR
      JP = 2*JI - JJ
   45 XBAR(JI,IS) = -COS(JP*B1)
      GO TO 60
!
   50 READ (5,1) ( XBAR(JI,IS), JI=1,MBAR )
    1 FORMAT ( 6F10.0 )
!
   60 NR = IABS(NS(IS))
      B1 = PI/(2*NR+1)
      IF ( ISTYPE(IS) .GT. 2 ) B1 = PI/(NR+1)
      COSTH = COS(THETA(IS))*BRINV
      SINTH = SIN(THETA(IS))*BRINV
      DO 70 JR=1,NR
      YBAR(JR,IS) = COS(JR*B1)
      Y = YBAR(JR,IS)*SO(IS)
      IF ( ISTYPE(IS) .GT. 2 ) Y = ( Y + SO(IS) )*0.5
      YS(JR,IS) = Y*COSTH + YV(IS)
      SIGY(JR,IS) = THETA(IS)
      ZY(JR,IS) = Y*SINTH + ZV(IS)
      CALL XMBY ( XM, BY, Y, BRPT(1,1,IS), NEND(IS) )
      BYBO (JR,IS) = BY*BRINV
      DO 70 JI=1,MBAR
   70 XS(JI,JR,IS) = ( BY*XBAR(JI,IS) + XM )*BRINV + XV(IS)
!
!
      J = NJ(IS)
      I = IXI(IS)
      GO TO ( 72, 74, 72, 85 ), I
!
   72 B1 = PI/(2*J+1)
      GO TO 80
!
   74 B1 = PI/(2*J)
!
   80 DO 82 J1=1,J
   82 XIBAR(J1,IS) = -COS((2*J1-1)*B1)
      GO TO 86
!
   85 READ (5,1) ( XIBAR(J1,IS), J1=1,J )
!
   86 IF ( ISTYPE(IS) .GT. 2 ) GO TO 130
!
!
      NS2 = NSI(IS)
      B1  = PI/(2*NS2)
      DO 90 JS=1,NS2
      ETABAR(JS,IS) = -COS( (2*JS-1)*B1 )
      ETA = ETABAR(JS,IS)*SO(IS)
      SIGETA(JS,IS) = THETA(IS)
      IF ( ETA .LT. 0 ) SIGETA(JS,IS) = -THETA(IS)
      ETAS(JS,IS) = ETA*BRINV*COS(SIGETA(JS,IS))
      ZETA(JS,IS) = ETA*BRINV*SIN(SIGETA(JS,IS)) + ZV(IS)
      CALL XMBY ( XM, BY, ETA, BRPT(1,1,IS), NEND(IS) )
      BETABO(JS,IS) = BY*BRINV
      XISLTE(1,JS,IS) = -BETABO(JS,IS) + XM*BRINV + XV(IS)
      XISLTE(2,JS,IS) =  BETABO(JS,IS) + XM*BRINV + XV(IS)
      DO 90 JJ=1,J
   90 XIS(JJ,JS,IS) = ( BY*XIBAR(JJ,IS) + XM )*BRINV + XV(IS)
      GO TO 200
!
!
  130 IF ( ISTYPE(IS) .EQ. 4 ) GO TO 160
!
      NS2 = NSI(IS)/2
      B1  = PI/(2*NS2)
      DO 140 JS=1,NS2
      ETABAR(JS,IS) = -COS( (2*JS-1)*B1 )
      ETABAR(JS+NS2,IS) = ETABAR(JS,IS)
      ETA = ( ETABAR(JS,IS) + 1.0 )*SO(IS)*0.5
      SIGETA(JS,IS) = THETA(IS)
      SIGETA(JS+NS2,IS) = -THETA(IS)
      ETAS(JS,IS) = ETA*BRINV*COS(SIGETA(JS,IS)) + YV(IS)
      ETAS(JS+NS2,IS) = -ETAS(JS,IS)
      ZETA(JS,IS) = ETA*BRINV*SIN(SIGETA(JS,IS)) + ZV(IS)
      ZETA(JS+NS2,IS) = ZETA(JS,IS)
!
      CALL XMBY ( XM, BY, ETA, BRPT(1,1,IS), NEND(IS) )
      BETABO (JS,IS) = BY*BRINV
      BETABO (JS+NS2,IS) = BETABO(JS,IS)
      XISLTE (1,JS,IS) = -BETABO(JS,IS) + XM*BRINV + XV(IS)
      XISLTE (2,JS,IS) =  BETABO(JS,IS) + XM*BRINV + XV(IS)
      XISLTE (1,JS+NS2,IS)= XISLTE(1,JS,IS)
      XISLTE (2,JS+NS2,IS)= XISLTE(2,JS,IS)
!
      DO 140 JJ=1,J
      XIS(JJ,JS,IS) = ( BY*XIBAR(JJ,IS) + XM )*BRINV + XV(IS)
  140 XIS(JJ,JS+NS2,IS) = XIS(JJ,JS,IS)
!
      GO TO 200
!
  160 NS2 =  NSI(IS)
      B1  = PI/(2*NS2)
      DO 170 JS=1,NS2
      ETABAR(JS,IS) = -COS( (2*JS-1)*B1 )
      ETA = ( ETABAR(JS,IS) + 1.0 )*SO(IS)*0.5
      SIGETA(JS,IS) = THETA(IS)
      ETAS(JS,IS) = ETA*BRINV*COS(SIGETA(JS,IS)) + YV(IS)
      ZETA(JS,IS) = ETA*BRINV*SIN(SIGETA(JS,IS)) + ZV(IS)
      CALL XMBY ( XM, BY, ETA, BRPT(1,1,IS), NEND(IS) )
      BETABO(JS,IS) = BY*BRINV
      XISLTE(1,JS,IS) = -BETABO(JS,IS) + XM*BRINV + XV(IS)
      XISLTE(2,JS,IS) =  BETABO(JS,IS) + XM*BRINV + XV(IS)
!
      DO 170 JJ=1,J
  170 XIS(JJ,JS,IS) = ( BY*XIBAR(JJ,IS) + XM )*BRINV + XV(IS)
!
!
  200 MBAR = NC(IS)
      NR   = IABS(NS(IS))
!
!
      J = NJ (IS)
      NS2= IABS(NSI(IS))
!
  300 IF ( IDUMP .EQ. 0 ) GO TO 1000
!
      WRITE (6,350) IS
  350 FORMAT (1H1,17X,37H DOWNWASH AND CONTROL POINT LOCATIONS          &
     &   ,//,     17X,27H   FOR LIFTING SURFACE NO. ,I2,// )
!
      WRITE (6,360)
  360 FORMAT (1H ,19X," THE DOWNWASH CHORD LOCATIONS ",//,              &
     &   3(25H   YBAR(I)       YS(I)   ), / )
      WRITE (6,380) ( YBAR(I,IS), YS(I,IS), I=1,NR )
  380 FORMAT (3(1X,2E12.4) )
!
      WRITE (6,382)
  382 FORMAT( /, 3(25H   SIGY(I)       ZY(I)   ), / )
      WRITE (6,380) ( SIGY(I,IS), ZY(I,IS), I=1,NR )
!
      WRITE (6,390)
  390 FORMAT (1H ,//,17X," THE CHORDWISE DOWNWASH LOCATIONS",/ )
      DO 399 JI=1,MBAR
      WRITE (6,396) JI, XBAR(JI,IS)
  396 FORMAT ( /," XBAR(",I2,")=",E12.4 )
      WRITE (6,397) JI
  397 FORMAT (" THE VALUES OF XS(JI,JR) FOR JI = ", I2, " ARE ", / )
      WRITE (6,398) ( XS(JI,JR,IS), JR=1,NR )
  398 FORMAT (1X, 6E12.4 )
!
  399 END DO
!

      WRITE (6,400)
  400 FORMAT (1H1,//,18X," THE INTEGRATION CHORD LOCATIONS ", //,       &
     &  3(25H  ETABAR(I)     ETAS(I)  ),/ )
      WRITE (6,380) ( ETABAR(I,IS), ETAS(I,IS), I=1,NS2 )
!
      WRITE (6,402)
  402 FORMAT( /, 3( 25H  SIGETA(I)     ZETA(I)  ), / )
      WRITE (6,380) ( SIGETA(I,IS), ZETA(I,IS), I=1,NS2 )
!
      WRITE (6,410)
  410 FORMAT (1H ,//,16X," THE CHORDWISE INTEGRATION LOCATIONS", / )
!
      DO 420 JJ=1,J
      WRITE (6,416) JJ, XIBAR(JJ,IS)
  416 FORMAT ( /, " XIBAR(",I2,")=", E11.4 )
      WRITE (6,417) JJ
  417 FORMAT (" THE VALUES OF XIS(JJ,JS) FOR JJ=",I2," ARE " )
      WRITE (6,398) ( XIS(JJ,JS,IS), JS=1,NS2 )
  420 END DO
!
!
 1000 END DO
!
 2000 LEAVE = 0
      RETURN
      END Subroutine Geomty


!!!      OVERLAY(R2T,3,0)
      SUBROUTINE Omanae()
      CALL MAINAE()
      END Subroutine Omanae

      SUBROUTINE MAINAE
!
!     SUBROUTINE TO CALL THE SUBROUTINES FOR COMPUTING THE SUBSONIC
!     OR SUPERSONIC, STEADY OR UNSTEADY AERODYNAMIC MATRICES.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), HMACH
!
      COMMON /OPTION/ IOP1, IOP2, IOP3, IOP4, NUNIT, IQT, IOPLU, IND
      COMMON /COMM/ LEAVE,IM,K,DUM(9)


!!!      EQUIVALENCE (NAME,TRANST,TRANUN)
!!!      DATA NAME/6HCOMPTY/
!
      REWIND 1
      IF ( RK(1) .GT. FREQ ) GO TO 20
!
!!!      CALL OVERLAY ( TRANST,3,7 )
      CALL Transt()

      IF ( NK .EQ. 1 ) GO TO 30
!
!!!   20 CALL OVERLAY ( TRANUN,3,4 )
   20 CALL Tranun()
!
   30 REWIND 1
      RETURN
      END Subroutine MainAE


      SUBROUTINE PHIX ( PHIX1, PHIX2, XMACH )
!
      COMMON /GAMA/ GAMMA
      G1    = 2.0/( GAMMA*XMACH*XMACH )
      G2    = GAMMA - 1.0
      G3    = GAMMA/G2
      G4    = G2*XMACH*XMACH
!
      PHIX1 = G1*( ( (2.0+G4)/(2.0+G2*PHIX1*PHIX1) )**G3 - 1.0 )
!
      PHIX2 = G1*( ( (2.0+G4)/(2.0+G2*PHIX2*PHIX2) )**G3 - 1.0 )
!
      RETURN
      END Subroutine Phix

      SUBROUTINE LOGSNG ( TL1, TL2, Y, Z, ETA, NS3, SY2, THETA, RK2 )
!
!     SUBROUTINE TO COMPUTE THE LOG SINGULARITY CORRECTION TERM
!     FOR NON-COPLANAR SURFACES IN UNSTEADY FLOW.
!
      DIMENSION  ETA(61), SY2(61)
!
      TL1 = 0.0
      TL2 = 0.0
      P1  = 3.141563/NS3
      Z2 = Z*Z
      DO 30 JS=1,NS3
      S1 = ETA(JS) - Y
      S2 = ALOG( ABS( S1*S1 + Z2 ) )*SY2(JS)
      TL1 = TL1 + S2
   30 TL2 = TL2 + S2*S1
!
      S1 = 1 + Y
      S2 = 1 - Y
      S3 = S1*S1 + Z2
      S4 = S2*S2 + Z2
      S5 = ALOG( S3 )
      S6 = ALOG( S4 )
!
      TL1 = -TL1*P1 + S1*S5 + S2*S6 + 2.0*Z*(ATAN2(S1,Z) +ATAN2(S2,Z) ) &
     &      - 4.0
!
      TL2 = -TL2*P1 - 0.5*( S3*S5 - S4*S6 - 4.0*Y )
!
      S1 = 0.25*( RK2*RK2 )*COS( THETA )
!
      TL1 = S1*TL1
      TL2 = S1*TL2
!
      RETURN
      END Subroutine LogSng

      FUNCTION DKERN ( Y, Z, THETA )
!
      IF ( Z .NE. 0.0 ) GO TO 10
!
      DKERN = 1.0/( Y*Y )
      GO TO 30
!
   10 IF ( THETA .NE. 0.0 ) GO TO 20
!
      R = 1.0/( Y*Y + Z*Z )
      DKERN = R - 2*((R*Z)**2)
      GO TO 30
!
   20 T = Z*( Z*COS(THETA) - Y*SIN(THETA) )
      R = 1.0/( Y*Y + Z*Z )
      DKERN = R*COS(THETA) - T*2.0*R*R
!
   30 RETURN
      END Function Dkern

      SUBROUTINE TRANSF ( THETA, YV, ZV, YR, ZR, YPR, ZPR )
!
!     SUBROUTINE TO ROTATE THE DOWNWASH POINT LOCATION THROUGH THETA
!     TO THE YPR AND ZPR COORDINATE SYSTEM.
!
      ATHETA = ABS(THETA)
      IF ( ATHETA .LT. 0.01 ) GO TO 10
      IF ( ATHETA .LT. 1.56 ) GO TO 20
      IF ( ATHETA .GT. 1.58 ) GO TO 20
!
      YPR = ( ZR - ZV )*THETA/ATHETA
      ZPR = ( YV - YR )*THETA/ATHETA
      GO TO 100
!
   10 YPR = ( YR - YV )
      ZPR = ( ZR - ZV )
      GO TO 100
!
   20 COST = COS( THETA )
      SINT = SIN( THETA )
      TANT = SINT/COST
      TAN2 = TANT*TANT
      TAN3 = 1.0/( 1.0 + TAN2 )
      DY   = YV - YR
      DZ   = ZV - ZR
!
      Y2   = ( YR + YV*TAN2 - DZ*TANT )*TAN3
      Z2   = ( ZV + ZR*TAN2 - DY*TANT )*TAN3
!
      DY   = Y2 - YV
      Y2   = Y2 - YR
      DZ   = Z2 - ZV
      Z2   = Z2 - ZR
!
      YPR  = DY*COST + DZ*SINT
      ZPR  = Y2*SINT - Z2*COST
!
  100 RETURN
      END Subroutine TranSF

      SUBROUTINE WEIGHT ( SY2, ETA, NS3, T, TN )
!
      DIMENSION SY2(NS3), ETA(NS3), T(NS3), TN(NS3)
!
      T(1) = 4.0/3.141593
      DO 10 J=3,NS3,2
   10 T(J) = -T(1)/(J*(J-2.0) )
      T(1) = T(1)*0.5
!
      TN(1) = 1.0
      DO 40 JS=1,NS3
      TN(2) = ETA(JS)
      ETA2  = 2*ETA(JS)
      DO 20 J=3,NS3
   20 TN(J) = ETA2*TN(J-1) - TN(J-2)
!
      SY2(JS) = T(1)
      DO 40 J=3,NS3,2
   40 SY2(JS) = SY2(JS) + TN(J)*T(J)
!
      RETURN
      END Subroutine Weight

      SUBROUTINE PLNTRM ( T1, T2, T3, Y, Z, THETA, ETA, LSPAN, NS3,     &
     &                    SY2, ETAA, ETAB )
!
      DIMENSION  ETA(61), SY2(61)
!
      P2  = 3.14159/NS3
      IF ( Z .NE. 0.0 ) GO TO 20
!
      F0 = 0.0
      F1 = 0.0
      F2 = 0.0
      DO 10 JS=1,NS3
      IF ( ETA(JS) .LT. ETAA .OR. ETA(JS) .GT. ETAB ) GO TO 10
      OY1 = 1.0/( ETA(JS) - Y )
      F0  = F0 + SY2(JS)*OY1*OY1
      F1  = F1 + SY2(JS)*OY1
      F2  = F2 + SY2(JS)
   10 END DO
!
      T1 = 1.0/(ETAA-Y) - 1.0/(ETAB-Y) - P2*F0
      T2 = ALOG( ABS( (ETAB-Y)/(ETAA-Y) ) ) - P2*F1
      T3 = ETAB - ETAA - P2*F2
      GO TO 100
!
   20 IF ( THETA .NE. 0.0 ) GO TO 50
      Z2 = Z*Z
      F0 = 0.0
      F1 = 0.0
      DO 30 JS=1,NS3
      IF ( ETA(JS) .LT. ETAA  .OR. ETA(JS).GT. ETAB ) GO TO 30
      R = 1.0/( (Y-ETA(JS))**2 + Z2 )
      G = R - 2*( (R*Z)**2 )
      G1= G*( ETA(JS) - Y )
!
      F0 = F0 + SY2(JS)*G
      F1 = F1 + SY2(JS)*G1
   30 END DO
!
      YP = Y - ETAA
      YM = Y - ETAB
      RP = 1.0/( YP*YP + Z2 )
      RM = 1.0/( YM*YM + Z2 )
!
      Q0 = ( YM*RM - YP*RP )
      Q1 = ( ALOG( ABS(RP/RM) )*0.5 + Z2*( RM - RP ) )
!
      T1 = ( Q0 - P2*F0 )
      T2 = ( Q1 - P2*F1 )
!
      GO TO 100
!
   50 Z2 = Z*Z
      F0 = 0.0
      F1 = 0.0
      F2 = 0.0
      DO 60 JS=1,NS3
      IF ( ETA(JS) .LT. ETAA  .OR. ETA(JS).GT. ETAB ) GO TO 60
      Y0 = Y - ETA(JS)
      R  = 1.0/( Y0*Y0 + Z2 )
      R2 = R*R
      G  = R - 2.0*R2*Z2
      G1 = -G*Y0
      H  = R2*Y0*Z
!
      F0 = F0 + SY2(JS)*G
      F1 = F1 + SY2(JS)*G1
      F2 = F2 + SY2(JS)*H
   60 END DO
!
      YP = Y - ETAA
      YM = Y - ETAB
      RP = 1.0/( YP*YP + Z2 )
      RM = 1.0/( YM*YM + Z2 )
!
      Q0 = ( YM*RM - YP*RP )
      Q1 = ( ALOG( ABS(RP/RM) )*0.5 + Z2*( RM - RP ) )
      Q2 = Z*( RM - RP )
!
      T1 = ( Q0 - P2*F0 )
      T2 = ( Q1 - P2*F1 )
      T3 = ( Q2 - P2*F2 )
!
      COSP = COS( THETA )
      SINP = SIN( THETA )
!
      T1 = T1*COSP + T3*SINP
      T2 = T2*COSP
!
  100 RETURN
      END Subroutine PlnTrm

      SUBROUTINE ETALIM ( ETAA, ETAB, SLE, BETA, XXLE, Y, Z )
!
      ETAA = -1.0
      ETAB = +1.0
!
      IF ( Z .NE. 0.0 ) GO TO 10
!
      IF ( BETA .LE. SLE ) GO TO 5
      ETAA = Y - XXLE/( BETA - SLE )
!
    5 IF ( BETA .LE. -SLE ) GO TO 100
      ETAB = Y + XXLE/( BETA + SLE )
      GO TO 100
!
   10 B2 = BETA**2
      B1 = B2 - SLE**2
      Z2 = Z*Z
      X2 = XXLE**2
!
      IF ( B1 .NE. 0.0 ) GO TO 20
!
      IF ( BETA .EQ. -SLE ) GO TO 15
      ETAB = Y + ( X2 - Z2*B2 )/( 2.0*BETA*XXLE )
      GO TO 100
!
   15 ETAA = Y - ( X2 - Z2*B2 )/( 2.0*BETA*XXLE )
      GO TO 100
!
   20 Y3 = X2 - Z2*B1
      IF ( Y3 .GT. 0.0 ) GO TO 25
      ETAB = -1.0
      GO TO 1000
!
   25 Y2 = SLE*XXLE/B1
      Y3 = SQRT( Y3 )*BETA/B1
!
      ETAA = Y - ( Y2 + Y3 )
      ETAB = Y - ( Y2 - Y3 )
!
      IF ( ETAA .LT. ETAB ) GO TO 100
!
      IF ( ETAA .GT. Y ) ETAA = -1.0
      IF ( ETAB .LT. Y ) ETAB = +1.0
!
  100 IF ( ABS(ETAA) .GT. 1.0 ) ETAA = SIGN( 1.0, ETAA )
      IF ( ABS(ETAB) .GT. 1.0 ) ETAB = SIGN( 1.0, ETAB )
!
 1000 RETURN
      END Subroutine EtaLim


!!!      OVERLAY(COMPTY,3,7)
      SUBROUTINE TRANST()
!
!     SUBROUTINE TO CALCULATE THE AERODYNAMIC MATRIX FOR MULTIPLE
!     ARBITRARY SURFACES IN STEADY SUBSONIC, MIXED TRANSONIC OR
!     SUPERSONIC FLOW.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), GMACH
!
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10), XIBAR(15,10), XIS(15,31,10), YBAR(15,10),         &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      COMMON /FUNCT/ GETA(61), GYY, HPSI(10,31), TN(10,67), DUNY(30),   &
     & SY2(61), TNX(10,2), UN(15,31), UNY(30), CK(6), YQ(6),QC(6),IC(6),&
     & ETAB2(62), ETAS2(62), SIGS(62), ZETAS(62), XIB3(22), XIS3(22),   &
     & XIB2(67), XIS2(67), W2(22), W1(67), XICT(67), NJ3, S6, Y, Z, XMC
!
      COMMON /GAMA/ GAMMA
      COMMON /KERN/ GJMU(15,31), CKERNL(67)
      DIMENSION ACOEF(140), D(20)
      DIMENSION EK(31), GMACHP(31), GMACHM(31)
!
      EQUIVALENCE ( ACOEF(1), A(1) )
      EQUIVALENCE ( EK(1), B(23) ), ( GMACHP(1), B(65) ),               &
     &                              ( GMACHM(1), B(96) )
!
      KR = 0
      DO 4000 IS=1,NSURF
!
      IF ( IS .EQ. 1 ) GO TO 1
!
      KR = KR + NW2(IS-1)
      IF ( KSURF(IS-1,IS-1) .LT. 0 ) KR = KR + NS(IS-1)
!
    1 MBAR = NC(IS)
      NR   = NS(IS)
!
      DO 4000 KS=1,NSURF
!
      IF ( JSUROP .NE. 0 ) GO TO 5000
!
      ISTP = ISTYPE(KS)
      LSPA = LSPAN(KS)
      ICHD = ICHORD(KS)
      GO TO 5001
!
 5000 LSPA = 5
      ISTP = 3
      IF ( ISTYPE(KS) .EQ. 4 ) ISTP = 4
      ICHD = 1
!
 5001 IF ( IDUMP .NE. 0 ) WRITE (6,5) IS, KS
    5 FORMAT (1H1,13X," THE DOWNWASH MATRIX FOR THE DOWNWASH INDUCED "  &
     &  ,/,       20X," ON SURFACE ", I2,", BY SURFACE ",I2,// )
      IF ( IS .EQ. KS ) GO TO 6
!
      IF ( KSURF(IS,KS) .GE. 0 ) GO TO 6
      MUNU = NW2(KS)
      IF ( KSURF(KS,KS) .LT. 0 ) MUNU = MUNU + NS(KS)
      N = NW2(IS) + KR
      IF ( KSURF(IS,IS) .LT. 0 ) N = N + NS(IS)
!
      KRP1 = KR + 1
      DO 3 IR=KRP1,N
      DO 2 I =1,MUNU
    2 ACOEF(I) = 0.0
      WRITE (1) ( ACOEF(I), I=1,MUNU )
      IF ( IDUMP .NE. 0 ) WRITE (6,2700) IR, (ACOEF(I),I=1,MUNU)
    3 END DO
      GO TO 4000
!
    6 IR = KR
      SOBR = SO(KS)*BRINV
      IF ( ISTYPE(KS) .GT. 2 ) SOBR = 0.5*SOBR
      NS2  = NSI(KS)
      JJ   = NJ(KS)
      JJMIN = 4
      NR2  = NS(KS)
      MBAR2= NC(KS)
      NS3  = NS2
      IF ( ISTYPE(KS) .EQ. 3 ) NS3 = NS2/2
!
      COSKS = COS( THETA(KS) )
      SINKS = SIN( THETA(KS) )
      IF ( ABS( COSKS ) .LT. 1.0E-5 ) COSKS = 0.0
      IF ( ABS( SINKS ) .LT. 1.0E-5 ) SINKS = 0.0
!
      P2 = 2.0*PI*PI/( (2*JJ+1)*NS3 )
      IF ( IXI(KS) .EQ. 2 ) P2 = PI*PI/( JJ*NS3 )
      P3 = P2*NS3/PI
!
      IF ( IS .EQ. KS ) GO TO 20
!
      CALL UNETA ( UN, ETABAR(1,KS), NS2, NR2, LSYM(KS), ISTP, LSPA )
!
      IF ( LSPA .NE. 1 ) GO TO 12
!
      DO 10 JS=1,NS2
      SETA    = 1.0 - ETABAR(JS,KS)**2
      SY2(JS) = SQRT( SETA )
      DO 10 JNU=1,NR2
   10 UN(JNU,JS) = UN(JNU,JS)*SETA
      GO TO 20
!
   12 CALL WEIGHT ( SY2(1), ETABAR(1,KS), NS3, A, A(71) )
!
      IF ( NS3 .NE. NS2 ) CALL WEIGHT ( SY2(NS3+1), ETABAR(NS3+1,KS),   &
     &                                  NS3, A, A(71) )
!
      DO 14 JS=1,NS2
      DO 14 JNU=1,NR2
   14 UN(JNU,JS) = UN(JNU,JS)*SY2(JS)
!
   20 ICH = 1
      JS1 = NS3
      JS2 = 1
!
      YV2 = YV(KS)
      ZV2 = ZV(KS)
      IF ( ISTYPE(KS) .LT. 3 ) GO TO 22
      YV2 = YV2 + SOBR*COSKS
      ZV2 = ZV2 + SOBR*SINKS
!
   22 DO 3000 JR=1,NR
!
      JUMP = 1
      IF ( IS .NE. KS ) GO TO 24
!
      Y = YBAR(JR,IS)
      Z = 0.0
      DTHETA = 0.0
      GO TO 25
!
   24 CALL TRANSF ( THETA(KS), YV2, ZV2, YS(JR,IS), ZY(JR,IS), Y, Z )
!
      Y = Y/SOBR
      Z = Z/SOBR
      Z2 = 0.01*( 1.0 - ABS(Y) )
      IF ( SO(IS) .EQ. SO(KS) .AND. NS(IS) .EQ. NS(KS) ) GO TO 23
      IF ( ABS(Z ) .LT. Z2 ) Z = SIGN( Z2, Z )
   23 DTHETA = THETA(IS) - THETA(KS)
!
      IF ( KSURF(IS,KS) .EQ. 0 ) GO TO 30
      IF ( ABS(Y) .GT. 0.999   ) GO TO 30
!
      JUMP = 0
      Y3 = Y + 0.001
      CALL UNETA ( DUNY, Y3, 1, NR2, LSYM(KS), ISTP, LSPA )
!
   25 CALL UNETA ( UNY,  Y,  1, NR2, LSYM(KS), ISTP, LSPA )
!
   30 CALL GEOMTS ( JR, IS, KS, NS2, NS3, JS1, JS2, YV2, ZV2,           &
     &  COSKS*SOBR, SINKS*SOBR, LSPA, ISTP )
!
      IF ( ABS(Y) .LT. 0.999 ) Y2 = SQRT( 1.0 - Y*Y )
!
      IF ( JUMP .NE. 0 ) GO TO 60
!
      DO 55 I=1,NR2
   55 DUNY(I) = 1000.0*( DUNY(I) - UNY(I) )
!
   60 SIGQ = SIGS(1) + 1.0
      COSP = COS( SIGY(JR,IS) )
      SINP = SIN( SIGY(JR,IS) )
      IF ( ABS(COSP) .LT. 1.0E-5 ) COSP = 0.0
      IF ( ABS(SINP) .LT. 1.0E-5 ) SINP = 0.0
!
      DO 3000 JI=1,MBAR
!
      Y5 = YBAR(JR,IS)
      IF ( ISTYPE(IS) .LT. 3 ) Y5 = 2*Y5 - 1.0
!
      GMACH = XMACH
      IF ( ITRANS(IS) .NE. 0 )                                          &
     &CALL CALCM ( XBAR(JI,IS), Y5, MBAR, NR, BMACH(1,IS), GMACH,       &
     &             B(1), B(12) )
!
      JJ = NJ(KS)
      IR = IR + 1
      BETA2 = 1.0 - GMACH**2
      BETA  = SQRT( ABS(BETA2) )
!
      IF ( JUMP .NE. 0 .AND. IS .NE. KS ) GO TO 250
      ETAA = -1.0
      ETAB =  1.0
      IF ( GMACH .LT. 1.0 ) GO TO 240
!
      SLE = ( BRPT(1,3,KS) - BRPT(1,1,KS) )/SO(KS)
      YP  = Y*SOBR
      IF ( ISTYPE(KS) .GT. 1 ) YP = YP + SOBR
      XXLE = (XS(JI,JR,IS) - ( XV(KS) - BO(KS)*BRINV + YP*SLE ) )/SOBR
!
      CALL ETALIM ( ETAA, ETAB, SLE, BETA, XXLE, Y, Z )
!
!
  240 CALL PLNTRM ( T1, T2, T3, Y, Z, DTHETA, ETAB2, LSPA, NS3, SY2,    &
     &              ETAA, ETAB )
!
  250 JS3 = 0
      DELT = 1.0 - ABS( YBAR(1,IS) )
  201 DO 200 JS=1,NS2
      MULT = 1
      IF ( ABS( ETAB2(JS) - Y ) .LT. DELT .OR. JS3 .EQ. 62 ) MULT = 3
!
      IF ( GMACH .GE. 1.0 ) GO TO 70
!
      IF ( MULT .NE. 1 ) GO TO 62
!
      JJ = NJ(KS)
      DO 61 J=1,JJ
   61 XICT(J) = XIBAR(J,KS)
      GO TO 64
!
   62 JJ = NJ(KS)*MULT
      IF ( IXI(KS) .NE. 2 ) JJ = JJ + 1
      P2 = PI/(2*JJ+1)
      IF ( IXI(KS) .EQ. 2 ) P2 = PI/(2*JJ)
      DO 63 J=1,JJ
   63 XICT(J) = -COS( (2*J-1)*P2 )
!
   64 P3 = PI/JJ
      IF ( IXI(KS) .NE. 2 ) P3 = 2.0*PI/( 2*JJ+1.0 )
!
      CALL TNXI ( TN, XICT, JJ, MBAR2, ICHD )
      IF ( JSUROP .NE. 0 ) GO TO 100
!
      DO 65 J=1,JJ
      Z3 = 1.0 - XICT(J)
      DO 65 JMU=1,MBAR2
      W1(J) = 1.0
   65 TN(JMU,J) = TN(JMU,J)*Z3
      GO TO 100
!
   70 IF ( IS .NE. KS ) GO TO 100
!
      JJ = ( JJMIN + JI )*MULT
      NJ3 = JJ
      P2  = PI/( 2*NJ3 )
      DO 80 J=1,NJ3
   80 XICT(J) = -COS( (2*J-1)*P2 )
!
  100 IF ( JS3 .EQ. 62 ) GO TO 202
!
!
      CALL TSTGEO ( JR, JI, JJ, JS, P3, NS3, IS, KS, COSKS, COSKS2, IR )
!
      IF ( P3 .NE. 0.0 ) GO TO 110
!
      DO 105 JMU=1,MBAR2
  105 GJMU(JMU,JS) = 0.0
      GO TO 200
!
  110 CALL TKERNS ( CKERNL, GMACH, JJ, 1,                               &
     &     XS(JI,JR,IS), YS(JR,IS), ZY(JR,IS), COSP, SINP, SIGS(JS),    &
     &     XIS2, ETAS2(JS), ZETAS(JS), COSKS , SINKS )
!
      IF ( GMACH .GE. 1.0 ) CALL TNXI ( TN, XIB2, JJ, MBAR2, +1 )
!
  120 CALL TSTPSI ( JR, JI, IS, KS, JS, P3, JJ, MBAR2, NS3, IR,         &
     &              COSP, SINP, COSKS, SINKS )
      DO 130 JMU=1,MBAR2
      GETA(JMU) =  HPSI(JMU,1)
      DO 130 J=1,JJ
  130 GETA(JMU) = GETA(JMU) + ( TN(JMU,J)*W1(J) )*CKERNL(J)
!
      DO 140 JMU=1,MBAR2
  140 GJMU(JMU,JS) = GETA(JMU)*P3
!
  200 END DO
!
      IF ( JUMP .NE. 0 .AND. IS .NE. KS ) GO TO 210
!
      JS3 = 62
      GO TO 201
!
  202 CALL TSTGEO ( JR, JI, JJ, 62, P3, NS3, IS, KS, COSKS, COSKS2, IR )
!
      CALL TSTPSI ( JR, JI, IS, KS, 62, P3, JJ, MBAR2, NS3, IR,         &
     &              COSP, SINP, COSKS, SINKS )
!
  210 MUNU = 0
      Y4 = Y*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
!
      DO 2000 JMU=1,MBAR2
      DO 2000 JNU=1,NR2
      MUNU = MUNU + 1
!
!     CALCULATION OF THE DOUBLE INTEGRALS WITH THE CHORDWISE
!     INTEGRALS FROM AB0VE CONTAINED IN ARRAY GJMU(JMU,JS).
!
      D(1) = 0.0
      DO 300 JS=1,NS2
      GETA(JS) = GJMU(JMU,JS)*UN(JNU,JS)
  300 D(1) = D(1) + GETA(JS)
      D(1) = D(1)*(SOBR**2)*PI/NS3
!
      IF ( IS .NE. KS ) GO TO 400
      IF ( JNU .GT. 1 ) GO TO 310
!
!     WRITE (6,2223) D(1), ( GETA(JS), JS=1,NS2 )
!2223 FORMAT (1H ,/,* D(1) = *,E14.7,* GETA = *,//, (10E12.4) )
      IF ( GMACH .GE. 1.0 ) GO TO 301
      D(2) = HPSI(JMU,1)*2.0*P3
!
!     D(2) IS THE CORRECTED DOWNWASH CHORD INTEGRAL COMPUTED
!     AS HPSI(JMU,62) IN TSTPSI.
!
      IF ( LSPAN(IS) .EQ. 1 ) D(2) = D(2)*Y2
  301 M1 = NS3
      M2 = 1
!
      ETA1 = ETAB2(M1) - YBAR(JR,IS)
      ETA2 = ETAB2(M2) - YBAR(JR,IS)
      S1   = ( (ETA1*SOBR)**2 )/SY2(M1)
      S2   = ( (ETA2*SOBR)**2 )/SY2(M2)
      IF ( GMACH .LT. 1.0 ) GO TO 310
!
      ETA4  = ETAB2(M2+1) - YBAR(JR,IS)
      Y1Y2  = ( ETA4/ETA1 )**2
      QC(1) = ( (ETA4*SOBR)**2 )/SY2(M1-1)
      QC(2) = S1
      QC(3) = S2
      QC(4) = ( (ETA4*SOBR)**2 )/SY2(M2+1)
!
      D(2) = -( GETA(M1-1)*QC(1) + GETA(M2+1)*QC(4) - ( GETA(M1)*QC(2)  &
     &        + GETA(M2)*QC(3) )*Y1Y2 )/( 2*UNY(1)*( 1 - Y1Y2 ) )
!
  310 GYY   = D(2)*UNY(JNU)
      D(11) = -GETA(M1)*S1
      D(13) = -GETA(M2)*S2
      D(12) = GYY
!     WRITE (6,2222) ( D(I), I=11,13 )
!2222 FORMAT (1H ,/,* D(JR-1) = *,E14.7,*  D(JR) = *,E14.7,*  D(JR+1) =
!    1*,E14.7 )
!
      D(20) = ( ( D(11) - D(12) ) - ( D(13) - D(12) )*ETA1/ETA2 )/      &
     &          ( ETA1**2 - ETA1*ETA2 )
      D(19) = ( D(11) - D(12) )/ETA1 - D(20)*ETA1
      D(18) = D(12)
!
      D(1)  = D(1) - ( D(18)*T1 + D(19)*T2 + D(20)*T3 )
!
      GO TO 500
!
  400 IF ( JUMP .NE. 0 ) GO TO 500
!
      IF ( JNU .GT. 1 ) GO TO 410
!
      D(2) = HPSI(JMU,1)*2.0*P3
!
      Q1 = T1
      Q2 = T2
!
      IF ( LSPAN(KS) .NE. 1 ) GO TO 410
!
      Q1 = T1*Y2 - T2*Y/Y2
      Q2 = T2*Y2
!
  410 D(1) = D(1) - D(2)*( Q1*UNY(JNU) + Q2*DUNY(JNU) )
!
  500 IF ( ISTYPE(KS) .GT. 2 ) D(1) = 2.0*D(1)
!
 2000 ACOEF(MUNU) = D(1)*XMACH/GMACH
!
      IF ( KSURF(KS,KS) .GE. 0 ) GO TO 2600
!
!     FOR THE KS'TH SURFACE TO BE THE SUBSONIC SURFACE OF A TRANSONIC
!     PAIR, THE DOWNWASH DUE TO THE SH0CK LINE DOUBLET MUST BE
!     CALCULATED ON THE IS'TH SURFACE. THIS IS PERFORMED IF
!     KSURF(KS,KS) IS LT. ZERO.
!
      DO 2200 JS=1,NS2
!
      IF ( JS .LE. NS3 .AND. IS .EQ. KS ) GO TO 2015
      BY  = BETABO(JS,KS)
      XSH = XISLTE(1,JS,KS)
      GO TO 2020
!
 2015 Y5 = ETAB2(JS)*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y5 = 0.5*( Y5 + SO(KS) )
!
      CALL XMBY ( XM, BY, Y5, BRPT(1,1,KS), NEND(KS) )
!
      XSH = ( XM - BY )*BRINV + XV(KS)
!
 2020 IF ( GMACH .LT. 1.0 ) GO TO 2100
!
      R = SQRT( ( YS(JR,IS) - ETAS2(JS) )**2                            &
     &        + ( ZY(JR,IS) - ZETAS(JS) )**2 )*BETA
      XMC = XS(JI,JR,IS) - R
      IF ( XMC .GT. XSH ) GO TO 2100
!
      GJMU(1,JS) = 0.0
      GO TO 2200
!
 2100 CALL TKERNS ( CKER, GMACH, 1, 1,                                  &
     &     XS(JI,JR,IS), YS(JR,IS), ZY(JR,IS), COSP, SINP, SIGS(JS),    &
     &     XSH, ETAS2(JS), ZETAS(JS), COSKS, SINKS )
!
      GJMU(1,JS) = CKER
!
 2200 END DO
!
      D(2) = 2.0
      IF ( IS .EQ. KS .AND. LSPAN(IS) .EQ. 1 ) D(2) = 2*Y2
!
      DO 2500 JNU=1,NR2
      MUNU = MUNU + 1
      D(1) = 0.0
      DO 2250 JS=1,NS2
      GETA(JS) = GJMU(1,JS)*UN(JNU,JS)
 2250 D(1) = D(1) + GETA(JS)
      D(1) = D(1)*(SOBR**2)*AR(KS)*PI/NS3
!
      IF ( IS .NE. KS ) GO TO 2300
!
      D(11) = -GETA(M1)*S1
      D(12) = D(2)*UNY(JNU)
      D(13) = -GETA(M2)*S2
!     WRITE (6,2224) D(1), T1, T2, T3
!2224 FORMAT (/,*    D(1) = *,E14.7, *  T1,T2,T3 = *, 3E14.7    )
!     WRITE (6,2222) ( D(I), I=11,13 )
!
      D(20) = ( ( D(11) - D(12) ) - ( D(13) - D(12) )*ETA1/ETA2 )/      &
     &          ( ETA1**2 - ETA1*ETA2 )
      D(19) = ( D(11) - D(12) )/ETA1 - D(20)*ETA1
      D(18) = D(12)
!
      D(1)  = D(1) - ( D(18)*T1 + D(19)*T2 + D(20)*T3 )*AR(KS)
!
      GO TO 2400
!
 2300 IF ( JUMP .NE. 0 ) GO TO 2400
!
      D(1) = D(1) - D(2)*( Q1*UNY(JNU) + Q2*DUNY(JNU) )*AR(KS)
!
 2400 IF ( ISTYPE(KS) .GT. 2 ) D(1) = 2.0*D(1)
!
 2500 ACOEF(MUNU) = D(1)*XMACH/( GMACH*SOBR )
!
 2600 WRITE (1) ( ACOEF(I), I=1,MUNU )
!
      IF ( IDUMP .NE. 0 ) WRITE (6,2700) IR, ( ACOEF(I), I=1,MUNU )
!
 2700 FORMAT (/,9X,"A COEFFICIENTS FOR ROW ",I3, //, ( 3X, 5E15.7 ) )
!
 3000 CONTINUE
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 4000
!
!     THE PRESSURE AND POTENTIAL VALUES ACROSS THE SHOCK ARE CALCULATED
!     AT THE YBAR STATIONS FOR SATISFYING THE NORMAL SHOCK B.C.'S
!      (PHIX+) - K*(PHI+)  + M(BAR)*(PHIX-) + K*(PHI-) = 0
!
!     THIS SECTION IS CURRENTLY RESTRICTED TO COPLANAR SURFACES FOR
!     WHICH THE SPAN OF DOWNSTREAM SURFACES IS .LE. TO THOSE UPSTREAM.
!
      ZMBAR = ( GAMMA-1.0 + 2.0/((  XMACH)**2) )/( GAMMA+1.0 )
!
      DO 3500 JR=1,NR
!
      IF ( IS .EQ. KS ) GO TO 3010
!
      DO 3005 I=1,MUNU
 3005 ACOEF(I) = 0.0
      IF( IABS( KSURF(KS,KS) ) .NE. IS ) GO TO 3450
!
!
      CALL TRANSF ( THETA(KS), YV2, ZV2, YS(JR,IS), ZY(JR,IS), Y5, Z )
!
      Y  = Y5/SOBR
      Z  =  Z/SOBR
      Z2 = 0.01*( 1.0 - ABS(Y) )
      IF ( ABS(Z) .LT. Z2 ) Z = SIGN( Z2,Z )
      DTHETA = THETA(IS) - THETA(KS)
      Y5 = Y5*BREF
      IF ( ISTYPE(KS) .GT. 2 ) Y5 = 0.5*( Y5 + SO(KS) )
!
      CALL XMBY ( XM, BY, Y5, BRPT(1,1,KS), NEND(KS) )
      XSH = ( XM + BY )*BRINV + XV(KS)
      GO TO 3020
!
 3010 Y  = YBAR(JR,IS)
      Z  = 0.0
      BY = BYBO(JR,IS)
      XSH = XS(1,JR,IS) - ( XBAR(1,IS) - 1.0 )*BY
!
 3020 CALL UNETA ( UNY, Y, 1, NR2, LSYM(KS), ISTYPE(KS), LSPAN(KS) )
!
      Y2 = -SQRT( 1.0 - Y*Y )
!
!
      YB = YBAR(JR,IS)
      IF ( ISTYPE(IS) .LT. 3 ) YB = 2*YB - 1.0
!
      CALL CALCM ( -1.0, YB, MBAR, NR, BMACH(1,IS), PHIP, B(1), B(12) )
      CALL CALCM ( -.99, YB, MBAR, NR, BMACH(1,IS), PHIPX,B(1), B(12) )
      GMACHP(JR) =  PHIP
      CALL PHIX ( PHIP, PHIPX, XMACH )
      PHIPX = 100.0*( PHIPX - PHIP )/BYBO(JR,IS)
!
      KT = -KSURF(IS,IS)
!
      CALL CALCM ( +1.0, YB, MBAR2, NR2, BMACH(1,KT), PHIM, B(1),B(12))
      CALL CALCM ( +.99, YB, MBAR2, NR2, BMACH(1,KT), PHIMX,B(1),B(12))
      GMACHM(JR) = PHIM
      CALL PHIX ( PHIM, PHIMX, XMACH )
      PHIMX = 100.0*( PHIM - PHIMX )/BYBO(JR,KT)
!
      EK(JR) = ( PHIPX + ZMBAR*PHIMX )/( PHIP - PHIM )
!
 3030 IF ( IS .EQ. KS ) GO TO 3100
!
!
 3050 MUNU = 0
      DO 3060 JMU=1,MBAR2
      DO 3060 JNU=1,NR2
      MUNU = MUNU + 1
!
 3060 ACOEF(MUNU) = 0.0
!
      IF ( KSURF(KS,KS) .NE. IS ) GO TO 3080
!
      X3 = 0.999
      Q3 = SQRT ( 0.001/1.999 )
!
      CALL TNXI ( TNX, X3, 1, MBAR2, ICHORD(KS) )
!
      CONST = 2*SOBR*ZMBAR/BYBO(JR,KS)
      IF ( ISTYPE(KS) .GT. 2 ) CONST = 2.0*CONST
!
      MUNU = 0
      DO 3070 JMU=1,MBAR2
      DO 3070 JNU=1,NR2
      MUNU = MUNU + 1
!
 3070 ACOEF(MUNU) = ACOEF(MUNU) + TNX(JMU,1)*UNY(JNU)*CONST*Q3
 3080 IF ( KSURF(KS,KS) ) 3120, 3400, 3400
!
 3100 X3 = -0.999
      Q3 = SQRT ( 1.999/0.001 )
!
      CALL TNXI ( TNX, X3, 1, MBAR2, ICHORD(KS) )
!
      CONST = 2*SOBR/BYBO(JR,KS)
      IF ( ISTYPE(KS) .GT. 2 ) CONST = 2.0*CONST
!
      MUNU = 0
      DO 3110 JMU=1,MBAR2
      DO 3110 JNU=1,NR2
      MUNU = MUNU + 1
!
 3110 ACOEF(MUNU) = TNX(JMU,1)*UNY(JNU)*CONST*Q3
!
 3120 CONST = -2*EK(JR)*AR(KS)
!
      DO 3130 JNU=1,NR2
      MUNU = MUNU + 1
!
 3130 ACOEF(MUNU) =  UNY(JNU)*CONST
!
 3400 IF ( LSPAN(KS) .NE. 1 ) Y2 = -1.0
!
      DO 3440 I=1,MUNU
 3440 ACOEF(I) = Y2*ACOEF(I)
!
 3450 WRITE (1) ( ACOEF(I), I=1,MUNU )
!
 3405 IF ( IDUMP .NE. 0 ) WRITE (6,3410) JR, ( ACOEF(I), I=1,MUNU )
!
 3410 FORMAT (/,9X,"SHOCK BC.S FOR DOWNWASH ROW ",I2,//,( 3X,5E15.7 ) )
!
 3500 END DO
!
 4000 CONTINUE
      END Subroutine TranST

      SUBROUTINE TSTPSI ( JR, JI, IS, KS, JS, P3, JJ, MBAR2, NS3, IR,   &
     &                    COSP, SINP, COSKS, SINKS )
!     SUBROUTINE TO CALCULATE THE CHORDWISE INTEGRAL CORRECTION TERM
!     AND THE DOWNWASH CHORD INTEGRAL IN TRANSONIC FLOW.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), GMACH
!
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10), XIBAR(15,10), XIS(15,31,10), YBAR(15,10),         &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      COMMON /FUNCT/ GETA(61), GYY, HPSI(10,31), TN(10,67), DUNY(30),   &
     & SY2(61), TNX(10,2), UN(15,31), UNY(30), CK(6), YQ(6),QC(6),IC(6),&
     & ETAB2(62), ETAS2(62), SIGS(62), ZETAS(62), XIB3(22), XIS3(22),   &
     & XIB2(67), XIS2(67), W2(22), W1(67), XICT(67), NJ3, S6, Y, Z, XMC
!
      SOBR = SO(KS)*BRINV
      IF ( ISTYPE(KS) .GT. 2 ) SOBR = SOBR*0.5
      Z1 = Z*SOBR
!
      IF ( GMACH .GE. 1.0 ) GO TO 30
!
      IF ( JS .NE. 62 ) GO TO 200
!
      IF ( Z .NE. 0.0 ) GO TO 50
!
      DO 20 JMU=1,MBAR2
      HPSI(JMU,1) = 0.0
      DO 10 J=1,JJ
      IF ( XIS2(J) .GT. XS(JI,JR,IS)) GO TO 20
   10 HPSI(JMU,1) = HPSI(JMU,1) + W1(J)*TN(JMU,J)
   20 END DO
      GO TO 1000
!
   30 IF ( JS .NE. 62 .AND. Z .EQ. 0.0 ) GO TO 200
!
      CALL TNXI ( TN, XIB2, JJ, MBAR2, +1 )
!
      IF ( Z .NE. 0.0 ) GO TO 50
!
      DO 40 JMU=1,MBAR2
      HPSI(JMU,1) = 0.0
      DO 40 J=1,JJ
   40 HPSI(JMU,1) = HPSI(JMU,1) + W1(J)*TN(JMU,J)
      GO TO 1000
!
   50 COSP2 = COSP*COSKS + SINP*SINKS
      SINP2 = SINP*COSKS - COSP*SINKS
!
      IF ( JS .NE. 62 ) GO TO 100
!
      CKINF = Z1*Z1*0.5/COSP2
!
      CALL TKERNS ( HPSI(1,2), GMACH, JJ, 1, XS(JI,JR,IS),              &
     &     0.0,Z1, COSP2, SINP2, (SIGETA(1,KS)-SIGY(JR,IS)),            &
     &     XIS2, 0.0, 0.0, 1.0, 0.0 )
!
      DO 60 J=1,JJ
   60 HPSI(J,2) = HPSI(J,2)*CKINF*W1(J)
!
      DO 70 JMU=1,MBAR2
      HPSI(JMU,1) = 0.0
      DO 70 J=1,JJ
   70 HPSI(JMU,1) = HPSI(JMU,1) + TN(JMU,J)*HPSI(J,2)
!
      IF ( GMACH .LT. 1.0 ) GO TO 1000
      GO TO 120
!
  100 DO 110 JMU=1,MBAR2
  110 HPSI(JMU,1) = 0.0
      Y0 = YS(JR,IS) - ETAS2(JS)
      Z0 = ZY(JR,IS) - ZETAS(JS)
      R2 = Y0**2 + Z0**2
      T2 =-(Z0*COSP-Y0*SINP)*(Z0*COSKS-Y0*SINKS)*2.0*BETA2/R2
      GO TO 130
!
  120 R2 = Z1*Z1
      T2 = -BETA2*0.5/COSP2
!
  130 XLE = YQ(1)
      BY = YQ(2)
      XISH = (XMC-XLE)*(XIB2(JJ)-XIB2(1))/(XIS2(JJ)-XIS2(1)) - 1.0
      XSH = XMC
      IF ( XISH .GT. 1.0 ) GO TO 1000
!
  135 CALL TNXI ( TNX(1,1), XISH, 1 , MBAR2, 1 )
      IND = 0
      CALL CHDTSS ( PSH, XISH, XISH, 1, ICHORD(KS), LSPAN(KS), IND, 1,  &
     &              XMACH, BRPT(1,1,KS), YQ(3), JSUROP )
!
      C1 = 1.0/( P3*BY*SQRT( (XLE-XS(JI,JR,IS))**2 + BETA2*R2 ) )
      IF ( JSUROP .EQ. 0 ) GO TO 139
      IF ( ISTYPE(KS) .GT. 2 ) C1 = C1*2
      C1 = C1*BREF/SO(KS)
!
  139 DO 140 J=1,JJ
      X0 = XS(JI,JR,IS) - XIS2(J)
  140 C1 = C1 + X0*SQRT( (1-XICT(J)**2)/( ( X0**2 + BETA2*R2 )**3 ) )
      C1 = C1*T2*PSH
      DO 150 JMU=1,MBAR2
  150 HPSI(JMU,1) = HPSI(JMU,1) + C1*TNX(JMU,1)
      GO TO 1000
!
  200 DO 210 JMU=1,MBAR2
  210 HPSI(JMU,1) = 0.0
!
 1000 RETURN
      END Subroutine TstPsi

      SUBROUTINE TSTGEO ( JR, JI, JJ, JS, P3, NS3, IS, KS,              &
     &                    COSKS, COSK2, IR )
!
!     SUBROUTINE TO CALCULATE THE CHORDWISE INTEGRATION POINT ARRAY
!     AT THE JS STATION AND OTHER GEOMETRIC DATA.
!
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), GMACH
!
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10), XIBAR(15,10), XIS(15,31,10), YBAR(15,10),         &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      COMMON /FUNCT/ GETA(61), GYY, HPSI(10,31), TN(10,67), DUNY(30),   &
     & SY2(61), TNX(10,2), UN(15,31), UNY(30), CK(6), YQ(6),QC(6),IC(6),&
     & ETAB2(62), ETAS2(62), SIGS(62), ZETAS(62), XIB3(22), XIS3(22),   &
     & XIB2(67), XIS2(67), W2(22), W1(67), XICT(67), NJ3, S6, Y, Z, XMC
!
      IF ( JS .EQ. 1 ) IND = 0
      BRAT = 0.0
      IF ( JSUROP .EQ. 0 .AND. GMACH .LT. 1.0 ) GO TO 500
      IF ( JS .NE. 62 ) GO TO 10
      IF ( IS .NE. KS ) GO TO 5
      Y4  = Y*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
      BY  = BYBO(JR,IS)
      R   = 0
      XMC = XS(JI,JR,IS)
      XLE = XMC - BYBO(JR,IS)*( XBAR(JI,IS) + 1.0 )
      GO TO 100
!
    5 Z4  = Z*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Z4 = 0.5*Z4
      R   = ABS( Z4*BRINV )*BETA
      XMC = XS(JI,JR,IS) - R
      Y4  = Y*SO(KS)
      GO TO 15
!
   10 R   = SQRT( (YS(JR,IS)-ETAS2(JS))**2 + (ZY(JR,IS)-ZETAS(JS))**2 ) &
     &      *BETA
      XMC = XS(JI,JR,IS) - R
      Y4  = ETAB2(JS)*SO(KS)
!
      IF ( JS .LE. NS3 .AND. IS .EQ. KS ) GO TO 15
      BY  = BETABO(JS,KS)
      XLE = XISLTE(1,JS,KS)
      XTE = XISLTE(2,JS,KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
      GO TO 20
!
   15 IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
!
      CALL XMBY ( XM, BY, Y4, BRPT(1,1,KS), NEND(KS) )
!
      XLE = ( XM - BY )*BRINV + XV(KS)
      XTE = ( XM + BY )*BRINV + XV(KS)
      BY  = BY*BRINV
!
   20 IF ( XMC .GT. XLE ) GO TO 30
      P3 = 0.0
      GO TO 1000
!
!     P3 IS SET TO ZERO IF THE MACH CONE IS FORWARD OF THE LEADING EDGE.
!
   30 IF ( XMC .LT. XTE ) GO TO 100
      IF ( JJ .NE. NJ(KS) ) GO TO 101
      IF ( IS .EQ. KS .OR. JS .EQ. 62 ) GO TO 101
!
      DO 40 J=1,JJ
      XICT(J) = XIBAR(J,KS)
      XIB2(J) = XIBAR(J,KS)
   40 XIS2(J) = XIS(J,JS,KS)
!
      CALL CHDTSS ( W1, XIB2, XIB2, JJ, ICHORD(KS), LSPAN(KS), IND, 0,  &
     &              XMACH, BRPT(1,1,KS), Y4, JSUROP )
!
   60 P3 = PI/JJ
      IF ( IXI(KS) .NE. 2 ) P3 = 2.0*PI/( 2.0*JJ + 1.0 )
      GO TO 1000
!
  100 XTE = XMC
  101 XM2 = 0.5*( XTE + XLE )
      BY2 = 0.5*( XTE - XLE )
      BRAT = BY2/BY
!
      IF ( IS .EQ. KS .OR. JJ .NE. NJ(KS) ) GO TO 200
!
  105 DO 110 J=1,JJ
      XICT(J) = XIBAR(J,KS)
      XIB2(J) = BRAT*( 1.0 + XIBAR(J,KS) ) - 1.0
  110 XIS2(J) = XM2 + BY2*XIBAR(J,KS)
!
      CALL CHDTSS ( W1, XIB2, XIBAR(1,KS), JJ, ICHORD(KS), LSPAN(KS),   &
     &              IND, 0, XMACH, BRPT(1,1,KS), Y4, JSUROP )
!
      GO TO 60
!
  200 P3 = PI/JJ
!
      DO 210 J=1,JJ
      XIB2(J) = BRAT*( 1.0 + XICT(J) ) - 1.0
  210 XIS2(J) = XM2 + BY2*XICT(J)
!
      CALL CHDTSS ( W1, XIB2, XICT, JJ, ICHORD(KS), LSPAN(KS), IND, 0,  &
     &              XMACH, BRPT(1,1,KS), Y4, JSUROP )
!
      GO TO 1000
!
  500 Y4  = Y*SO(KS)
      IF ( JS .NE. 62 ) Y4 = ETAB2(JS)*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
!
      CALL XMBY ( XM, BY, Y4, BRPT(1,1,KS), NEND(KS) )
!
      BY  = BY*BRINV
      XM  = XM*BRINV + XV(KS)
!
      DO 510 J=1,JJ
  510 XIS2(J) = XM + XICT(J)*BY
!
 1000 IF ( BRAT .NE. 0.0 ) P3 = P3*BRAT
      IF ( JSUROP .NE. 0 ) P3 = P3*BREF*BY/SO(KS)
      YQ(1) = XLE
      YQ(2) = BY
      YQ(3) = Y4
      RETURN
      END Subroutine TstGeo

      SUBROUTINE TKERNS ( CK, XMACH, JJ, NS2, X, Y, Z, COSP, SINP, SIG, &
     &                    XI, ETA, ZETA, COSQ, SINQ )
!
!     SUBROUTINE TO COMPUTE THE PLANAR OR NONPLANAR SUBSONIC OR
!     SUPERSONIC KERNEL FUNCTION FOR STEADY TRANSONIC FLOW.
!
      DIMENSION CK(15,31), XI(15,31), ETA(61), ZETA(61), SIG(61)
!
      BETA2 = 1.0 - XMACH*XMACH
!
      IF ( ABS(Z-ZETA(1)) .GT. ABS(Y-ETA(1))*0.01 ) GO TO 100
      IF ( ABS(SINP-SINQ) .GT. ABS(COSP+COSQ)*0.005 ) GO TO 100
!
      IF ( XMACH .LT. 1.0 ) GO TO 20
!
      Y0 = Y - ETA(1)
      Y2 = Y0*Y0
      DO 10 J=1,JJ
      X0 = X - XI(J,1)
   10 CK(J,1) = -2.*X0/( Y2*SQRT( X0*X0 + BETA2*Y2 ) )
      GO TO 300
!
   20 DO 30 I=1,NS2
      Y0 = Y - ETA(I)
      Y2 = Y0*Y0
      DO 30 J=1,JJ
      X0 = X - XI(J,I)
   30 CK(J,I) = - ( 1.0 + X0/SQRT( X0*X0 + BETA2*Y2 ) )/Y2
      GO TO 300
!
  100 T1 = COSP*COSQ + SINP*SINQ
!
      IF ( XMACH .LT. 1.0 ) GO TO 150
!
      Y0 = Y - ETA(1)
      Z0 = Z - ZETA(1)
      T2 = ( Z0*COSP - Y0*SINP )*( Z0*COSQ - Y0*SINQ )
      R2 = Y0*Y0 + Z0*Z0
      OR2 = 1.0/R2
!
      DO 110 J=1,JJ
      X0 = X - XI(J,1)
      RR2= 1.0/( X0*X0 + BETA2*R2 )
      RR = SQRT( RR2 )
      CK1= 2.0*X0*RR
!
  110 CK(J,1) = -OR2*CK1*( T1 - T2*( 2.0*OR2 + BETA2*RR2 ) )
      GO TO 300
!
  150 NS  = NS2
      NS1 = 1
      NS3 = NS2/2
      COSQ2 = COSQ
      SINQ2 = SINQ
  160 DO 170 JS=NS1,NS
      Y0  = Y - ETA(JS)
      Z0  = Z - ZETA(JS)
      T2  = ( Z0*COSP - Y0*SINP )*( Z0*COSQ2 - Y0*SINQ2 )
      R2  = Y0*Y0 + Z0*Z0
      OR2 = 1.0/R2
!
      DO 170 J=1,JJ
      X0  = X - XI(J,JS)
      RR2 = 1.0/( X0*X0 + BETA2*R2 )
      RR  = SQRT( RR2 )
      CK1 = 1.0 + X0*RR
!
  170 CK(J,JS) = -OR2*( T1*CK1 - T2*( 2*OR2*CK1 + BETA2*X0*RR*RR2 ) )
      IF ( NS .EQ. NS2 ) GO TO 300
!
      NS1 = NS + 1
      NS  = NS2
      COSQ2 = -COSQ
      T1    = -COSP*COSQ + SINP*SINQ
      GO TO 160
!
  300 RETURN
      END Subroutine Tkerns

      SUBROUTINE GEOMTS ( JR, IS, KS, NS2, NS3, JS1, JS2, YV2, ZV2,     &
     &                       CKSS, SKSS, LSPA, ISTP )
!
!     SUBROUTINE TO CALCULATE THE TRANSFORMED SPANWISE INTEGRATION
!     POINT ARRAY FOR IS .EQ. KS.
!     IF IS .NE. KS, THE ARRAYS ARE SET UP BUT NO TRANSFORMATION IS
!     PERFORMED.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), GMACH
!
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10), XIBAR(15,10), XIS(15,31,10), YBAR(15,10),         &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      COMMON /FUNCT/ GETA(61), GYY, HPSI(10,31), TN(10,67), DUNY(30),   &
     & SY2(61), TNX(10,2), UN(15,31), UNY(30), CK(6), YQ(6),QC(6),IC(6),&
     & ETAB2(62), ETAS2(62), SIGS(62), ZETAS(62), XIB3(22), XIS3(22),   &
     & XIB2(67), XIS2(67), W2(22), W1(67), XICT(67), NJ3, S6, Y, Z, XMC
!
      NR2  = NS(KS)
      IF ( IS .NE. KS ) GO TO 45
!
      DO 10 JS=1,NS3
      IF ( Y .LT. ETABAR(JS,KS) ) GO TO 20
   10 END DO
!
   20 JS1 = NS3 + 1 - JS
      JS2 = JS1 + 1
!
      DO 25 JS=1,JS1
      ETAB2(JS) = Y + ETABAR(JS,KS) + 1.0
      ETAS2(JS) = YV2 + ETAB2(JS)*CKSS
      ZETAS(JS) = ZV2 + ETAB2(JS)*SKSS
   25 SIGS(JS)  = SIGETA(JS,KS)
!
      DO 30 JS=JS2,NS3
      ETAB2(JS) = Y + ETABAR(JS,KS) - 1.0
      ETAS2(JS) = YV2 + ETAB2(JS)*CKSS
      ZETAS(JS) = ZV2 + ETAB2(JS)*SKSS
   30 SIGS(JS)  = SIGETA(JS,KS)
!
      IF ( ISTYPE(KS) .NE. 3 ) GO TO 50
!
      NS3P1 = NS3 + 1
   35 DO 40 JS=NS3P1,NS2
      ETAB2(JS) = ETABAR(JS,KS)
      ETAS2(JS) = ETAS(JS,KS)
      ZETAS(JS) = ZETA(JS,KS)
   40 SIGS(JS)  = SIGETA(JS,KS)
      IF ( IS - KS ) 100, 50, 100
!
   45 NS3P1 = 1
      GO TO 35
!
   50 CALL UNETA ( UN, ETAB2, NS2, NR2, LSYM(KS), ISTP, LSPA )
!
      IF ( LSPA .NE. 1 ) GO TO 60
!
      DO 55 JS=1,NS2
      SY2(JS) = SQRT( 1.0 - ETABAR(JS,KS)**2 )
      SETA    = SY2(JS)*SQRT( 1.0 - ETAB2(JS)**2 )
      DO 55 JNU=1,NR2
   55 UN(JNU,JS) = UN(JNU,JS)*SETA
      GO TO 100
!
   60 CALL WEIGHT ( SY2(1), ETABAR(1,KS), NS3, A, A(71) )
!
      IF ( NS3 .NE. NS2 )                                               &
     &CALL WEIGHT ( SY2(NS3+1), ETABAR(NS3+1,KS), NS3, A, A(71) )
!
      DO 65 JS=1,NS2
      DO 65 JNU=1,NR2
   65 UN(JNU,JS) = UN(JNU,JS)*SY2(JS)
!
  100 RETURN
      END Subroutine GeomTs

!!!      OVERLAY(COMPTY,3,4)
      SUBROUTINE Tranun()
!
!     SUBROUTINE TO CALCULATE THE AERODYNAMIC MATRIX FOR MULTIPLE
!     ARBITRARY SURFACES OSCILLATING IN SUBSONIC, MIXED TRANSONIC
!     OR SUPERSONIC FLOW.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), GMACH
!
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10), XIBAR(15,10), XIS(15,31,10), YBAR(15,10),         &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      COMMON /UNFUN/ GETA(61), GYY, HPSI(10,31), TN(10,67), DUNY(30),   &
     & SY2(61), TNX(10,2), UN(15,31), UNY(30), CK(6), YQ(6),QC(6),IC(6),&
     & ETAB2(62), ETAS2(62), SIGS(62), ZETAS(62), XIB3(22), XIS3(22),   &
     & XIB2(67), XIS2(67), W2(22), W1(67), XICT(67), NJ3, S6, Y, Z, XMC
!
      COMMON /GAMA/ GAMMA
      COMMON /KERNL/  GJMU(15,31), CKERNL(67)
      COMPLEX GETA, GYY, HPSI, CKERNL, ACOEF, D, GJMU, RK3
      COMPLEX EK, CONST, CKER
      DIMENSION ACOEF(70), D(20)
      DIMENSION EK(31), GMACHP(31), GMACHM(31)
!
      EQUIVALENCE ( ACOEF(1), A(1) )
      EQUIVALENCE ( EK(1), B( 1) ), ( GMACHP(1), B(65) ),               &
     &                              ( GMACHM(1), B(96) )
!
      NK1 = 1
      IF ( RK(1) .LE. FREQ ) NK1 = 2
      DO 5000 IK=NK1,NK
      KR = 0
      DO 4000 IS=1,NSURF
!
      IF ( IS .EQ. 1 ) GO TO 1
!
      KR = KR + NW2(IS-1)
      IF ( KSURF(IS-1,IS-1) .LT. 0 ) KR = KR + NS(IS-1)
!
    1 MBAR = NC(IS)
      NR   = NS(IS)
!
      DO 4000 KS=1,NSURF
!
      IF ( JSUROP .NE. 0 ) GO TO 6000
!
      ISTP = ISTYPE(KS)
      LSPA = LSPAN(KS)
      ICHD = ICHORD(KS)
      GO TO 6001
!
 6000 LSPA = 5
      ISTP = 3
      IF ( ISTYPE(KS) .EQ. 4 ) ISTP = 4
      ICHD = 1
!
 6001 IF ( IDUMP .NE. 0 ) WRITE (6,5) IS, KS
    5 FORMAT (1H1,13X," THE DOWNWASH MATRIX FOR THE DOWNWASH INDUCED "  &
     &  ,/,       20X," ON SURFACE ", I2,", BY SURFACE ",I2,// )
      IF ( IS .EQ. KS ) GO TO 6
!
      IF ( KSURF(IS,KS) .GE. 0 ) GO TO 6
      MUNU = NW2(KS)
      IF ( KSURF(KS,KS) .LT. 0 ) MUNU = MUNU + NS(KS)
      N = NW2(IS) + KR
      IF ( KSURF(IS,IS) .LT. 0 ) N = N + NS(IS)
!
      KRP1 = KR + 1
      DO 3 IR=KRP1,N
      DO 2 I =1,MUNU
    2 ACOEF(I) = 0.0
      WRITE (1) ( ACOEF(I), I=1,MUNU )
      IF ( IDUMP .NE. 0 ) WRITE (6,2700) IR, (ACOEF(I),I=1,MUNU)
    3 END DO
      GO TO 4000
!
    6 IR = KR
      SOBR = SO(KS)*BRINV
      IF ( ISTYPE(KS) .GT. 2 ) SOBR = 0.5*SOBR
      NS2  = NSI(KS)
      JJ   = NJ(KS)
      JJMIN= 4
      NR2  = NS(KS)
      MBAR2= NC(KS)
      NS3  = NS2
      IF ( ISTYPE(KS) .EQ. 3 ) NS3 = NS2/2
!
      COSKS = COS( THETA(KS) )
      SINKS = SIN( THETA(KS) )
      IF ( ABS( COSKS ) .LT. 1.0E-5 ) COSKS = 0.0
      IF ( ABS( SINKS ) .LT. 1.0E-5 ) SINKS = 0.0
!
      P2 = 2.0*PI*PI/( (2*JJ+1)*NS3 )
      IF ( IXI(KS) .EQ. 2 ) P2 = PI*PI/( JJ*NS3 )
      P3 = P2*NS3/PI
!
      IF ( IS .EQ. KS ) GO TO 20
!
      CALL UNETA ( UN, ETABAR(1,KS), NS2, NR2, LSYM(KS), ISTP, LSPA )
!
      IF ( LSPA .NE. 1 ) GO TO 12
!
      DO 10 JS=1,NS2
      SETA    = 1.0 - ETABAR(JS,KS)**2
      SY2(JS) = SQRT( SETA )
      DO 10 JNU=1,NR2
   10 UN(JNU,JS) = UN(JNU,JS)*SETA
      GO TO 20
!
   12 CALL WEIGHT ( SY2(1), ETABAR(1,KS), NS3, A, A(71) )
!
      IF ( NS3 .NE. NS2 ) CALL WEIGHT ( SY2(NS3+1), ETABAR(NS3+1,KS),   &
     &                                  NS3, A, A(71) )
!
      DO 14 JS=1,NS2
      DO 14 JNU=1,NR2
   14 UN(JNU,JS) = UN(JNU,JS)*SY2(JS)
!
   20 ICH = 1
      JS1 = NS3
      JS2 = 1
!
      YV2 = YV(KS)
      ZV2 = ZV(KS)
      IF ( ISTYPE(KS) .LT. 3 ) GO TO 22
      YV2 = YV2 + SOBR*COSKS
      ZV2 = ZV2 + SOBR*SINKS
!
   22 DO 3000 JR=1,NR
!
      JUMP = 1
      IF ( IS .NE. KS ) GO TO 24
!
      Y = YBAR(JR,IS)
      Z = 0.0
      DTHETA = 0.0
      GO TO 25
!
   24 CALL TRANSF ( THETA(KS), YV2, ZV2, YS(JR,IS), ZY(JR,IS), Y, Z )
!
      Y = Y/SOBR
      Z = Z/SOBR
      IF ( SO(IS) .EQ. SO(KS) .AND. NS(IS) .EQ. NS(KS) ) GO TO 23
      Z2 = 0.01*( 1.0 - ABS(Y) )
      IF ( ABS(Z ) .LT. Z2 ) Z = SIGN( Z2, Z )
   23 DTHETA = THETA(IS) - THETA(KS)
!
      IF ( KSURF(IS,KS) .EQ. 0 ) GO TO 30
      IF ( ABS(Y) .GT. 0.999   ) GO TO 30
!
      JUMP = 0
      Y3 = Y + 0.001
      CALL UNETA ( DUNY, Y3, 1, NR2, LSYM(KS), ISTP, LSPA )
!
      CALL LOGSNG ( TL1, TL2, Y, Z, ETAB2, NS3, SY2, DTHETA, 1.0 )
!
   25 CALL UNETA ( UNY,  Y,  1, NR2, LSYM(KS), ISTP, LSPA )
!
   30 CALL GEOMTU ( JR, IS, KS, NS2, NS3, JS1, JS2, YV2, ZV2,           &
     &              COSKS*SOBR, SINKS*SOBR, LSPA, ISTP )
!
      IF ( ABS(Y) .LT. 0.999 ) Y2 = SQRT( 1.0 - Y*Y )
!
      IF ( JUMP .NE. 0 ) GO TO 60
!
      DO 55 I=1,NR2
   55 DUNY(I) = 1000.0*( DUNY(I) - UNY(I) )
!
   60 SIGQ = SIGS(1) + 1.0
      COSP = COS( SIGY(JR,IS) )
      SINP = SIN( SIGY(JR,IS) )
      IF ( ABS(COSP) .LT. 1.0E-5 ) COSP = 0.0
      IF ( ABS(SINP) .LT. 1.0E-5 ) SINP = 0.0
!
      DO 3000 JI=1,MBAR
!
      Y5 = YBAR(JR,IS)
      IF ( ISTYPE(IS) .LT. 3 ) Y5 = 2*Y5 - 1.0
      GMACH = XMACH
!
      IF ( ITRANS(IS) .NE. 0 )                                          &
     &CALL CALCM ( XBAR(JI,IS), Y5, MBAR, NR, BMACH(1,IS), GMACH,       &
     &             D(1), D(12) )
!
      JJ = NJ(KS)
      IR = IR + 1
      BETA2 = 1.0 - GMACH**2
      BETA  = SQRT( ABS( BETA2 ) )
      RK2 = RK(IK)*XMACH/GMACH
!
      IF ( JUMP .NE. 0 .AND. IS .NE. KS ) GO TO 250
      ETAA = -1.0
      ETAB =  1.0
      IF ( GMACH .LT. 1.0 ) GO TO 240
!
      SLE = ( BRPT(1,3,KS) - BRPT(1,1,KS) )/SO(KS)
      YP  = Y*SOBR
      IF ( ISTYPE(KS) .GT. 1 ) YP = YP + SOBR
      XXLE = (XS(JI,JR,IS) - ( XV(KS) - BO(KS)*BRINV + YP*SLE ) )/SOBR
!
      CALL ETALIM ( ETAA, ETAB, SLE, BETA, XXLE, Y, Z )
!
  240 CALL PLNTRM ( T1, T2, T3, Y, Z, DTHETA, ETAB2, LSPA, NS3, SY2,    &
     &              ETAA, ETAB )
!
  250 JS3 = 0
      DELT = 1.0 - ABS( YBAR(1,IS) )
  201 DO 200 JS=1,NS2
      MULT = 1
      IF ( ABS( ETAB2(JS) - Y ) .LT. DELT .OR. JS3 .EQ. 62 ) MULT = 3
!
      IF ( GMACH .GE. 1.0 ) GO TO 70
!
      IF ( MULT .NE. 1 ) GO TO 62
!
      JJ = NJ(KS)
      DO 61 J=1,JJ
   61 XICT(J) = XIBAR(J,KS)
      GO TO 64
!
   62 JJ = NJ(KS)*MULT
      IF ( IXI(KS) .NE. 2 ) JJ = JJ + 1
      P2 = PI/(2*JJ+1)
      IF ( IXI(KS) .EQ. 2 ) P2 = PI/(2*JJ)
      DO 63 J=1,JJ
   63 XICT(J) = -COS( (2*J-1)*P2 )
!
   64 P3 = PI/JJ
      IF ( IXI(KS) .NE. 2 ) P3 = 2.0*PI/( 2*JJ+1.0 )
!
      CALL TNXI ( TN, XICT, JJ, MBAR2, ICHD )
      IF ( JSUROP .NE. 0 ) GO TO 100
!
      DO 65 J=1,JJ
      Z3 = 1.0 - XICT(J)
      DO 65 JMU=1,MBAR2
      W1(J) = 1.0
   65 TN(JMU,J) = TN(JMU,J)*Z3
      GO TO 100
!
   70 IF ( IS .NE. KS ) GO TO 100
!
      JJ = ( JJMIN + JI )*MULT
      NJ3 = JJ
      P2  = PI/( 2*NJ3 )
      DO 80 J=1,NJ3
   80 XICT(J) = -COS( (2*J-1)*P2 )
!
  100 IF ( JS3 .EQ. 62 ) GO TO 202
!
!
      CALL TUNGEO ( JR, JI, JJ, JS, P3, NS3, IS, KS, COSKS, COSKS2, IR )
!
      IF ( P3 .NE. 0.0 ) GO TO 110
!
      DO 105 JMU=1,MBAR2
  105 GJMU(JMU,JS) = 0.0
      GO TO 200
!
  110 CALL TKERNU ( CKERNL,       GMACH, JJ, 1, RK2,                    &
     &     XS(JI,JR,IS), YS(JR,IS), ZY(JR,IS), COSP, SINP, SIGS(JS),    &
     &     XIS2, ETAS2(JS), ZETAS(JS), COSKS , SINKS )
!
      IF ( GMACH .GE. 1.0 ) CALL TNXI ( TN, XIB2, JJ, MBAR2, +1 )
!
  120 CALL TUNPSI ( JR, JI, IS, KS, JS, P3, JJ, MBAR2, NS3, IR,         &
     &              COSP, SINP, COSKS, SINKS, RK2 )
      DO 130 JMU=1,MBAR2
      GETA(JMU) =  HPSI(JMU,1)
      DO 130 J=1,JJ
  130 GETA(JMU) = GETA(JMU) + ( TN(JMU,J)*W1(J) )*CKERNL(J)
!
      DO 140 KMU=1,MBAR2
  140 GJMU(KMU,JS) = GETA(KMU)*P3
!
  200 END DO
!
      IF ( JUMP .NE. 0 .AND. IS .NE. KS ) GO TO 210
!
      JS3 = 62
      GO TO 201
!
  202 CALL TUNGEO ( JR, JI, JJ, 62, P3, NS3, IS, KS,                    &
     &              COSKS, COSKS2, IR )
!
      IF ( P3 .EQ. 0 ) GO TO 210
!
      CALL TUNPSI ( JR, JI, IS, KS, 62, P3, JJ, MBAR2, NS3, IR,         &
     &              COSP, SINP, COSKS, SINKS, RK2 )
!
  210 MUNU = 0
      Y4 = Y*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
!
      DO 2000 JMU=1,MBAR2
      DO 2000 JNU=1,NR2
      MUNU = MUNU + 1
!
!     CALCULATION OF THE DOUBLE INTEGRALS WITH THE CHORDWISE
!     INTEGRALS FROM AB0VE CONTAINED IN ARRAY GJMU(JMU,JS).
!
      D(1) = 0.0
      DO 300 JS=1,NS2
      GETA(JS) = GJMU(JMU,JS)*UN(JNU,JS)
  300 D(1) = D(1) + GETA(JS)
      D(1) = D(1)*(SOBR**2)*PI/NS3
!
      IF ( IS .NE. KS ) GO TO 400
      IF ( JNU .GT. 1 ) GO TO 310
!
!     WRITE (6,2223) D(1), ( GETA(JS), JS=1,NS2 )
 2223 FORMAT (1H ,/," D(1)= ",2E14.7," GETA = ",//, (10E12.4) )
      IF ( GMACH .GE. 1.0 ) GO TO 301
      D(2) = HPSI(JMU,1)*2.0*P3
!
!     D(2) IS THE CORRECTED DOWNWASH CHORD INTEGRAL COMPUTED
!     AS HPSI(JMU,62) IN TSTPSI.
!
      IF ( LSPAN(IS) .EQ. 1 ) D(2) = D(2)*Y2
  301 M1 = NS3
      M2 = 1
      RK3 = CMPLX( 0.0, RK2*BYBO(JR,IS) )
      ETA1 = ETAB2(M1) - YBAR(JR,IS)
      ETA2 = ETAB2(M2) - YBAR(JR,IS)
      S1   = ( (ETA1*SOBR)**2 )/SY2(M1)
      S2   = ( (ETA2*SOBR)**2 )/SY2(M2)
      IF ( GMACH .LT. 1.0 ) GO TO 310
!
      ETA4  = ETAB2(M2+1) - YBAR(JR,IS)
      Y1Y2  = ( ETA4/ETA1 )**2
      QC(1) = ( (ETA4*SOBR)**2 )/SY2(M1-1)
      QC(2) = S1
      QC(3) = S2
      QC(4) = ( (ETA4*SOBR)**2 )/SY2(M2+1)
!
      D(2) = -( GETA(M1-1)*QC(1) + GETA(M2+1)*QC(4) - ( GETA(M1)*QC(2)  &
     &        + GETA(M2)*QC(3) )*Y1Y2 )/( 2*UNY(1)*( 1 - Y1Y2 ) )
!
      DELT1 = 0.0
      DELT2 = 0.0
  310 GYY   = D(2)*UNY(JNU)
      D(11) = -GETA(M1)*S1
      D(13) = -GETA(M2)*S2
      D(12) = GYY
!     WRITE (6,2222) ( D(I), I=11,13 )
 2222 FORMAT (1H ,/," D(JR-1,JR,JR+1) =", 3(2X,2E14.6) )
!
      D(20) = ( ( D(11) - D(12) ) - ( D(13) - D(12) )*ETA1/ETA2 )/      &
     &          ( ETA1**2 - ETA1*ETA2 )
      D(19) = ( D(11) - D(12) )/ETA1 - D(20)*ETA1
      D(18) = D(12)
!
      D(1)  = D(1) - ( D(18)*T1 + D(19)*T2 + D(20)*T3 )
!
      GO TO 500
!
  400 IF ( JUMP .NE. 0 ) GO TO 500
!
      IF ( JNU .GT. 1 ) GO TO 410
!
      D(2) = HPSI(JMU,1)*2.0*P3
!
      IF ( JMU .GT. 1 ) GO TO 410
!
      Q1 = T1 + TL1*RK2*RK2
      Q2 = T2 + TL2*RK2*RK2
!
      IF ( LSPAN(KS) .NE. 1 ) GO TO 410
!
      Q1 = Q1*Y2 - Q2*Y/Y2
      Q2 = Q2*Y2
!
  410 D(1) = D(1) - D(2)*( Q1*UNY(JNU) + Q2*DUNY(JNU) )
!
  500 IF ( ISTYPE(KS) .GT. 2 ) D(1) = 2.0*D(1)
!
 2000 ACOEF(MUNU) = D(1)*XMACH/GMACH
!
      IF ( KSURF(KS,KS) .GE. 0 ) GO TO 2600
!
!     FOR THE KS'TH SURFACE TO BE THE SUBSONIC SURFACE OF A TRANSONIC
!     PAIR, THE DOWNWASH DUE TO THE SH0CK LINE DOUBLET MUST BE
!     CALCULATED ON THE IS'TH SURFACE. THIS IS PERFORMED IF
!     KSURF(KS,KS) IS LT. ZERO.
!
      DO 2200 JS=1,NS2
!
      IF ( JS .LE. NS3 .AND. IS .EQ. KS ) GO TO 2015
      BY  = BETABO(JS,KS)
      XSH = XISLTE(1,JS,KS)
      GO TO 2020
!
 2015 Y5 = ETAB2(JS)*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y5 = 0.5*( Y5 + SO(KS) )
!
      CALL XMBY ( XM, BY, Y5, BRPT(1,1,KS), NEND(KS) )
!
      XSH = ( XM - BY )*BRINV + XV(KS)
!
 2020 IF ( GMACH .LT. 1.0 ) GO TO 2100
!
      R = SQRT( ( YS(JR,IS) - ETAS2(JS) )**2                            &
     &        + ( ZY(JR,IS) - ZETAS(JS) )**2 )*BETA
      XMC = XS(JI,JR,IS) - R
      IF ( XMC .GT. XSH ) GO TO 2100
!
      GJMU(1,JS) = 0.0
      GO TO 2200
!
 2100 CALL TKERNU ( CKER, GMACH, 1, 1, RK2,                             &
     &     XS(JI,JR,IS), YS(JR,IS), ZY(JR,IS), COSP, SINP, SIGS(JS),    &
     &     XSH, ETAS2(JS), ZETAS(JS), COSKS, SINKS )
!
      GJMU(1,JS) = CKER
!
 2200 END DO
!
      IF ( JUMP .NE. 0 .AND. IS .NE. KS ) GO TO 2210
      CALL XMBY ( XM, BY, Y4, BRPT(1,1,KS), NEND(KS) )
      XSH = ( XM - BY )*BRINV + XV(KS)
      D(2) = 2.0*CEXP( CMPLX( 0.0, RK2*(XSH-XS(JI,JR,IS)) ) )
      IF ( IS .EQ. KS .AND. LSPAN(IS) .EQ. 1 ) D(2) = D(2)*Y2
!
 2210 DO 2500 JNU=1,NR2
      MUNU = MUNU + 1
      D(1) = 0.0
      DO 2250 JS=1,NS2
      GETA(JS) = GJMU(1,JS)*UN(JNU,JS)
 2250 D(1) = D(1) + GETA(JS)
      D(1) = D(1)*(SOBR**2)*AR(KS)*PI/NS3
!
      IF ( IS .NE. KS ) GO TO 2300
!
      D(11) = -GETA(M1)*S1
      D(12) = D(2)*UNY(JNU)
      D(13) = -GETA(M2)*S2
!     WRITE (6,2224) D(1), T1, T2, T3
!     WRITE (6,2222) ( D(I), I=11,13 )
 2224 FORMAT (/,"    D(1) = ",E14.7, "  T1,T2,T3 = ", 3E14.7    )
!
      D(20) = ( ( D(11) - D(12) ) - ( D(13) - D(12) )*ETA1/ETA2 )/      &
     &          ( ETA1**2 - ETA1*ETA2 )
      D(19) = ( D(11) - D(12) )/ETA1 - D(20)*ETA1
      D(18) = D(12)
!
      D(1)  = D(1) - ( D(18)*T1 + D(19)*T2 + D(20)*T3 )*AR(KS)
!
      GO TO 2400
!
 2300 IF ( JUMP .NE. 0 ) GO TO 2400
!
      D(1) = D(1) - D(2)*( Q1*UNY(JNU) + Q2*DUNY(JNU) )*AR(KS)
!
 2400 IF ( ISTYPE(KS) .GT. 2 ) D(1) = 2.0*D(1)
!
 2500 ACOEF(MUNU) = D(1)*XMACH/( GMACH*SOBR )
!
 2600 WRITE (1) ( ACOEF(I), I=1,MUNU )
!
      IF ( IDUMP .NE. 0 ) WRITE (6,2700) IR, ( ACOEF(I), I=1,MUNU )
!
 2700 FORMAT (/,9X,"A COEFFICIENTS FOR ROW ",I3, //, ( 3X, 5E15.7 ) )
!
 3000 CONTINUE
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 4000
!
!     THE PRESSURE AND POTENTIAL VALUES ACROSS THE SHOCK ARE CALCULATED
!     AT THE YBAR STATIONS FOR SATISFYING THE NORMAL SHOCK B.C.'S
!      (PHIX+) - K*(PHI+)  + M(BAR)*(PHIX-) + K*(PHI-) = 0
!
!     THIS SECTION IS CURRENTLY RESTRICTED TO COPLANAR SURFACES FOR
!     WHICH THE SPAN OF DOWNSTREAM SURFACES IS .LE. TO THOSE UPSTREAM.
!
      ZMBAR = ( GAMMA-1.0 + 2.0/((  XMACH)**2) )/( GAMMA+1.0 )
!
      RK3   = CMPLX( 0.0, RK(IK)*( 2.0*(GAMMA-1.0)/(GAMMA+1.0) ) )
!
! ***** NU = RK3 = I(RK2*(2(GAMMA-1)/(GAMMA+1))),GAMMA = 1.4
!
      DO 3500 JR=1,NR
!
      IF ( IS .EQ. KS ) GO TO 3010
!
      DO 3005 I=1,MUNU
 3005 ACOEF(I) = 0.0
      IF( IABS( KSURF(KS,KS) ) .NE. IS ) GO TO 3450
!
!
      CALL TRANSF ( THETA(KS), YV2, ZV2, YS(JR,IS), ZY(JR,IS), Y5, Z )
!
      Y  = Y5/SOBR
      Z  =  Z/SOBR
      Z2 = 0.01*( 1.0 - ABS(Y) )
      IF ( ABS(Z) .LT. Z2 ) Z = SIGN( Z2,Z )
      DTHETA = THETA(IS) - THETA(KS)
      Y5 = Y5*BREF
      IF ( ISTYPE(KS) .GT. 2 ) Y5 = 0.5*( Y5 + SO(KS) )
!
      CALL XMBY ( XM, BY, Y5, BRPT(1,1,KS), NEND(KS) )
      XSH = ( XM + BY )*BRINV + XV(KS)
      GO TO 3020
!
 3010 Y  = YBAR(JR,IS)
      Z  = 0.0
      BY = BYBO(JR,IS)
      XSH = XS(1,JR,IS) - ( XBAR(1,IS) - 1.0 )*BY
!
 3020 CALL UNETA ( UNY, Y, 1, NR2, LSYM(KS), ISTYPE(KS), LSPAN(KS) )
!
      Y2 = -SQRT( 1.0 - Y*Y )
!
!
      YB = YBAR(JR,IS)
      IF ( ISTYPE(IS) .LT. 3 ) YB = 2*YB - 1.0
!
      CALL CALCM ( -1.0, YB, MBAR, NR, BMACH(1,IS), PHIP, D(1), D(12) )
      CALL CALCM ( -.99, YB, MBAR, NR, BMACH(1,IS), PHIPX,D(1), D(12) )
      GMACHP(JR) =  PHIP
      CALL PHIX ( PHIP, PHIPX, XMACH )
      PHIPX = 100.0*( PHIPX - PHIP )/BYBO(JR,IS)
!
      KT = -KSURF(IS,IS)
!
      CALL CALCM ( +1.0, YB, MBAR2, NR2, BMACH(1,KT), PHIM, D(1),D(12))
      CALL CALCM ( +.99, YB, MBAR2, NR2, BMACH(1,KT), PHIMX,D(1),D(12))
      GMACHM(JR) = PHIM
      CALL PHIX ( PHIM, PHIMX, XMACH )
      PHIMX = 100.0*( PHIM - PHIMX )/BYBO(JR,KT)
!
      EK(JR) = ( PHIPX + ZMBAR*PHIMX + RK3*( PHIM - 2.0/(GAMMA-1.0)))/  &
     &                                     ( PHIP - PHIM )
!
 3030 IF ( IS .EQ. KS ) GO TO 3100
!
      CONST = 2.0*SOBR*RK3
      IF ( ISTYPE(KS) .GT. 2 ) CONST = 2*CONST
      IF ( JR .GT. 1 ) GO TO 3050
!
      JJ = NJ(KS)
      CALL TNXI ( TN, XIBAR(1,KS), JJ, MBAR2, ICHORD(KS) )
      P3 = 2.0*PI/( 2*JJ+1.0 )
      DO 3040 J=1,JJ
 3040 W1(J) = P3*( 1.0 - XIBAR(J,KS) )
 3050 DO 3045 JMU=1,MBAR2
      GETA(JMU) = 0.0
      DO 3045 J=1,JJ
 3045 GETA(JMU) = GETA(JMU) + W1(J)*TN(JMU,J)*CEXP( CMPLX( 0.0,         &
     &            RK(IK)*( XIBAR(J,KS) - 1.0 )*BYBO(JR,KS) ) )
!
      MUNU = 0
      DO 3060 JMU=1,MBAR2
      DO 3060 JNU=1,NR2
      MUNU = MUNU + 1
!
 3060 ACOEF(MUNU) = GETA(JMU)*UNY(JNU)*CONST
!
      IF ( KSURF(KS,KS) .NE. IS ) GO TO 3080
!
      X3 = 0.999
      Q3 = SQRT ( 0.001/1.999 )
!
      CALL TNXI ( TNX, X3, 1, MBAR2, ICHORD(KS) )
!
      CONST = 2*SOBR*ZMBAR/BYBO(JR,KS)
      IF ( ISTYPE(KS) .GT. 2 ) CONST = 2.0*CONST
!
      MUNU = 0
      DO 3070 JMU=1,MBAR2
      DO 3070 JNU=1,NR2
      MUNU = MUNU + 1
!
 3070 ACOEF(MUNU) = ACOEF(MUNU) + TNX(JMU,1)*UNY(JNU)*CONST*Q3
 3080 IF ( KSURF(KS,KS) ) 3120, 3400, 3400
!
 3100 X3 = -0.999
      Q3 = SQRT ( 1.999/0.001 )
!
      CALL TNXI ( TNX, X3, 1, MBAR2, ICHORD(KS) )
!
      CONST = 2*SOBR/BYBO(JR,KS)
      IF ( ISTYPE(KS) .GT. 2 ) CONST = 2.0*CONST
!
      MUNU = 0
      DO 3110 JMU=1,MBAR2
      DO 3110 JNU=1,NR2
      MUNU = MUNU + 1
!
 3110 ACOEF(MUNU) = TNX(JMU,1)*UNY(JNU)*CONST*Q3
!
 3120 CONST = -2*EK(JR)*AR(KS)
!
      DO 3130 JNU=1,NR2
      MUNU = MUNU + 1
!
 3130 ACOEF(MUNU) =  UNY(JNU)*CONST
!
 3400 IF ( LSPAN(KS) .NE. 1 ) Y2 = -1.0
!
      DO 3440 I=1,MUNU
 3440 ACOEF(I) = Y2*ACOEF(I)
!
 3450 WRITE (1) ( ACOEF(I), I=1,MUNU )
!
 3405 IF ( IDUMP .NE. 0 ) WRITE (6,3410) JR, ( ACOEF(I), I=1,MUNU )
!
 3410 FORMAT (/,9X,"SHOCK BC.S FOR DOWNWASH ROW ",I2,//,( 3X,5E15.7 ) )
!
 3500 END DO
!
 4000 CONTINUE
!
 5000 END DO
      END Subroutine TranUn

      SUBROUTINE TUNPSI ( JR, JI, IS, KS, JS, P3, JJ, MBAR2, NS3, IR,   &
     &                    COSP, SINP, COSKS, SINKS, RK2 )
!     SUBROUTINE TO CALCULATE THE CHORDWISE INTEGRAL CORRECTION TERM
!     AND THE DOWNWASH CHORD INTEGRAL IN TRANSONIC FLOW.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), GMACH
!
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10), XIBAR(15,10), XIS(15,31,10), YBAR(15,10),         &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      COMMON /UNFUN/ GETA(61), GYY, HPSI(10,31), TN(10,67), DUNY(30),   &
     & SY2(61), TNX(10,2), UN(15,31), UNY(30), CK(6), YQ(6),QC(6),IC(6),&
     & ETAB2(62), ETAS2(62), SIGS(62), ZETAS(62), XIB3(22), XIS3(22),   &
     & XIB2(67), XIS2(67), W2(22), W1(67), XICT(67), NJ3, S6, Y, Z, XMC
      COMPLEX GETA, GYY, HPSI, RK3, C2
!
      RK3 = CMPLX( 0.0, RK2 )
      SOBR = SO(KS)*BRINV
      IF ( ISTYPE(KS) .GT. 2 ) SOBR = SOBR*0.5
      Z1 =  Z*SOBR
!
      IF ( GMACH .GE. 1.0 ) GO TO 30
!
      IF ( JS .NE. 62 ) GO TO 200
!
      IF ( Z .NE. 0.0 ) GO TO 50
!
      DO 20 JMU=1,MBAR2
      HPSI(JMU,1) = 0.0
      DO 10 J=1,JJ
      IF ( XIS2(J) .GT. XS(JI,JR,IS)) GO TO 20
   10 HPSI(JMU,1) = HPSI(JMU,1) + W1(J)*TN(JMU,J)                       &
     &                          *CEXP( RK3*(XIS2(J)-XS(JI,JR,IS)) )
   20 END DO
      GO TO 1000
!
   30 IF ( JS .NE. 62 .AND. Z .EQ. 0.0 ) GO TO 200
!
      CALL TNXI ( TN, XIB2, JJ, MBAR2, +1 )
!
      IF ( Z .NE. 0.0 ) GO TO 50
!
      DO 40 JMU=1,MBAR2
      HPSI(JMU,1) = 0.0
      DO 40 J=1,JJ
   40 HPSI(JMU,1) = HPSI(JMU,1) + W1(J)*TN(JMU,J)                       &
     &                          *CEXP( RK3*(XIS2(J)-XS(JI,JR,IS)) )
      GO TO 1000
!
   50 COSP2 = COSP*COSKS + SINP*SINKS
      SINP2 = SINP*COSKS - COSP*SINKS
!
      IF ( JS .NE. 62 ) GO TO 100
!
      CKINF = Z1*Z1*0.5/COSP2
!
      CALL TKERNU ( HPSI(1,2), GMACH, JJ, 1, RK2, XS(JI,JR,IS),         &
     &     0.0,Z1, COSP2, SINP2, (SIGETA(1,KS)-SIGY(JR,IS)),            &
     &     XIS2, 0.0, 0.0, 1.0, 0.0 )
!
      DO 60 J=1,JJ
   60 HPSI(J,2) = HPSI(J,2)*CKINF*W1(J)
!
      DO 70 JMU=1,MBAR2
      HPSI(JMU,1) = 0.0
      DO 70 J=1,JJ
   70 HPSI(JMU,1) = HPSI(JMU,1) + TN(JMU,J)*HPSI(J,2)
!
      IF ( GMACH .LT. 1.0 ) GO TO 1000
      GO TO 120
!
  100 DO 110 JMU=1,MBAR2
  110 HPSI(JMU,1) = 0.0
      Y0 = YS(JR,IS) - ETAS2(JS)
      Z0 = ZY(JR,IS) - ZETAS(JS)
      R2 = Y0**2 + Z0**2
      T2 =-(Z0*COSP-Y0*SINP)*(Z0*COSKS-Y0*SINKS)*2.0*BETA2/R2
      GO TO 130
!
  120 R2 = Z1*Z1
      T2 = -BETA2*0.5/COSP2
!
  130 XLE = YQ(1)
      BY = YQ(2)
      XISH = (XMC-XLE)*(XIB2(JJ)-XIB2(1))/(XIS2(JJ)-XIS2(1)) - 1.0
      XSH = XMC
      IF ( XISH .GT. 1.0 ) GO TO 1000
!
  135 CALL TNXI ( TNX(1,1), XISH, 1 , MBAR2, 1 )
      IND = 0
      CALL CHDTSS ( PSH, XISH, XISH, 1, ICHORD(KS), LSPAN(KS), IND, 1,  &
     &              XMACH, BRPT(1,1,KS), YQ(3), JSUROP )
!
      C1 = 1.0/( P3*BY*SQRT( (XLE-XS(JI,JR,IS))**2 + BETA2*R2 ) )
      IF ( JSUROP .NE. 0 ) C1 = C1*BREF/SO(KS)
!
      DO 140 J=1,JJ
      X0 = XS(JI,JR,IS) - XIS2(J)
  140 C1 = C1 + X0*SQRT( (1-XICT(J)**2)/( ( X0**2 + BETA2*R2 )**3 ) )
      C2 = C1*T2*PSH*CEXP( RK3*( XSH - XS(JI,JR,IS) ) )
      DO 150 JMU=1,MBAR2
  150 HPSI(JMU,1) = HPSI(JMU,1) + C2*TNX(JMU,1)
      GO TO 1000
!
  200 DO 210 JMU=1,MBAR2
  210 HPSI(JMU,1) = 0.0
!
 1000 RETURN
      END Subroutine TunPsi

      SUBROUTINE TUNGEO ( JR, JI, JJ, JS, P3, NS3, IS, KS,              &
     &                    COSKS, COSK2, IR )
!
!     SUBROUTINE TO CALCULATE THE CHORDWISE INTEGRATION POINT ARRAY
!     AT THE JS STATION AND OTHER GEOMETRIC DATA.
!
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), GMACH
!
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10), XIBAR(15,10), XIS(15,31,10), YBAR(15,10),         &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      COMMON /UNFUN/ GETA(61), GYY, HPSI(10,31), TN(10,67), DUNY(30),   &
     & SY2(61), TNX(10,2), UN(15,31), UNY(30), CK(6), YQ(6),QC(6),IC(6),&
     & ETAB2(62), ETAS2(62), SIGS(62), ZETAS(62), XIB3(22), XIS3(22),   &
     & XIB2(67), XIS2(67), W2(22), W1(67), XICT(67), NJ3, S6, Y, Z, XMC
!
      COMPLEX GETA, GYY, HPSI
!
      IF ( JS .EQ. 1 ) IND = 0
      BRAT = 0.0
      IF ( JSUROP .EQ. 0 .AND. GMACH .LT. 1.0 ) GO TO 500
      IF ( JS .NE. 62 ) GO TO 10
      IF ( IS .NE. KS ) GO TO 5
      Y4  = Y*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
      BY  = BYBO(JR,IS)
      R   = 0
      XMC = XS(JI,JR,IS)
      XLE = XMC - BYBO(JR,IS)*( XBAR(JI,IS) + 1.0 )
      GO TO 100
!
    5 Z4  = Z*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Z4 = 0.5*Z4
      R   = ABS( Z4*BRINV )*BETA
      XMC = XS(JI,JR,IS) - R
      Y4  = Y*SO(KS)
      GO TO 15
!
   10 R   = SQRT( (YS(JR,IS)-ETAS2(JS))**2 + (ZY(JR,IS)-ZETAS(JS))**2 ) &
     &      *BETA
      XMC = XS(JI,JR,IS) - R
      Y4  = ETAB2(JS)*SO(KS)
!
      IF ( JS .LE. NS3 .AND. IS .EQ. KS ) GO TO 15
      BY  = BETABO(JS,KS)
      XLE = XISLTE(1,JS,KS)
      XTE = XISLTE(2,JS,KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
      GO TO 20
!
   15 IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
!
      CALL XMBY ( XM, BY, Y4, BRPT(1,1,KS), NEND(KS) )
!
      XLE = ( XM - BY )*BRINV + XV(KS)
      XTE = ( XM + BY )*BRINV + XV(KS)
      BY  = BY*BRINV
!
   20 IF ( XMC .GT. XLE ) GO TO 30
      P3 = 0.0
      GO TO 1000
!
!     P3 IS SET TO ZERO IF THE MACH CONE IS FORWARD OF THE LEADING EDGE.
!
   30 IF ( XMC .LT. XTE ) GO TO 100
      IF ( JJ .NE. NJ(KS) ) GO TO 101
      IF ( IS .EQ. KS .OR. JS .EQ. 62 ) GO TO 101
!
      DO 40 J=1,JJ
      XICT(J) = XIBAR(J,KS)
      XIB2(J) = XIBAR(J,KS)
   40 XIS2(J) = XIS(J,JS,KS)
!
      CALL CHDTSS ( W1, XIB2, XIB2, JJ, ICHORD(KS), LSPAN(KS), IND, 0,  &
     &              XMACH, BRPT(1,1,KS), Y4, JSUROP )
!
   60 P3 = PI/JJ
      IF ( IXI(KS) .NE. 2 ) P3 = 2.0*PI/( 2.0*JJ + 1.0 )
      GO TO 1000
!
  100 XTE = XMC
  101 XM2 = 0.5*( XTE + XLE )
      BY2 = 0.5*( XTE - XLE )
      BRAT = BY2/BY
!
      IF ( IS .EQ. KS .OR. JJ .NE. NJ(KS) ) GO TO 200
!
  105 DO 110 J=1,JJ
      XICT(J) = XIBAR(J,KS)
      XIB2(J) = BRAT*( 1.0 + XIBAR(J,KS) ) - 1.0
  110 XIS2(J) = XM2 + BY2*XIBAR(J,KS)
!
      CALL CHDTSS ( W1, XIB2, XIBAR(1,KS), JJ, ICHORD(KS), LSPAN(KS),   &
     &              IND, 0, XMACH, BRPT(1,1,KS), Y4, JSUROP )
!
      GO TO 60
!
  200 P3 = PI/JJ
!
      DO 210 J=1,JJ
      XIB2(J) = BRAT*( 1.0 + XICT(J) ) - 1.0
  210 XIS2(J) = XM2 + BY2*XICT(J)
!
      CALL CHDTSS ( W1, XIB2, XICT, JJ, ICHORD(KS), LSPAN(KS), IND, 0,  &
     &              XMACH, BRPT(1,1,KS), Y4, JSUROP )
!
      GO TO 1000
!
  500 Y4  = Y*SO(KS)
      IF ( JS .NE. 62 ) Y4 = ETAB2(JS)*SO(KS)
      IF ( ISTYPE(KS) .GT. 2 ) Y4 = 0.5*( Y4 + SO(KS) )
!
      CALL XMBY ( XM, BY, Y4, BRPT(1,1,KS), NEND(KS) )
!
      BY  = BY*BRINV
      XM  = XM*BRINV + XV(KS)
!
      DO 510 J=1,JJ
  510 XIS2(J) = XM + XICT(J)*BY
!
 1000 IF ( BRAT .NE. 0.0 ) P3 = P3*BRAT
      IF ( JSUROP .NE. 0 ) P3 = P3*BREF*BY/SO(KS)
      YQ(1) = XLE
      YQ(2) = BY
      YQ(3) = Y4
      RETURN
      END Subroutine TunGeo

      SUBROUTINE GEOMTU ( JR, IS, KS, NS2, NS3, JS1, JS2, YV2, ZV2,     &
     &                    CKSS, SKSS, LSPA, ISTP )
!
!     SUBROUTINE TO CALCULATE THE TRANSFORMED SPANWISE INTEGRATION
!     POINT ARRAY FOR IS .EQ. KS.
!     IF IS .NE. KS, THE ARRAYS ARE SET UP BUT NO TRANSFORMATION IS
!     PERFORMED.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), GMACH
!
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), ITRANS(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10), XIBAR(15,10), XIS(15,31,10), YBAR(15,10),         &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      COMMON /UNFUN/ GETA(61), GYY, HPSI(10,31), TN(10,67), DUNY(30),   &
     & SY2(61), TNX(10,2), UN(15,31), UNY(30), CK(6), YQ(6),QC(6),IC(6),&
     & ETAB2(62), ETAS2(62), SIGS(62), ZETAS(62), XIB3(22), XIS3(22),   &
     & XIB2(67), XIS2(67), W2(22), W1(67), XICT(67), NJ3, S6, Y, Z, XMC
      COMPLEX GETA, GYY, HPSI
!
      NR2  = NS(KS)
      IF ( IS .NE. KS ) GO TO 45
!
      DO 10 JS=1,NS3
      IF ( Y .LT. ETABAR(JS,KS) ) GO TO 20
   10 END DO
!
   20 JS1 = NS3 + 1 - JS
      JS2 = JS1 + 1
!
      DO 25 JS=1,JS1
      ETAB2(JS) = Y + ETABAR(JS,KS) + 1.0
      ETAS2(JS) = YV2 + ETAB2(JS)*CKSS
      ZETAS(JS) = ZV2 + ETAB2(JS)*SKSS
   25 SIGS(JS)  = SIGETA(JS,KS)
!
      DO 30 JS=JS2,NS3
      ETAB2(JS) = Y + ETABAR(JS,KS) - 1.0
      ETAS2(JS) = YV2 + ETAB2(JS)*CKSS
      ZETAS(JS) = ZV2 + ETAB2(JS)*SKSS
   30 SIGS(JS)  = SIGETA(JS,KS)
!
      IF ( ISTYPE(KS) .NE. 3 ) GO TO 50
!
      NS3P1 = NS3 + 1
   35 DO 40 JS=NS3P1,NS2
      ETAB2(JS) = ETABAR(JS,KS)
      ETAS2(JS) = ETAS(JS,KS)
      ZETAS(JS) = ZETA(JS,KS)
   40 SIGS(JS)  = SIGETA(JS,KS)
      IF ( IS - KS ) 100, 50, 100
!
   45 NS3P1 = 1
      GO TO 35
!
   50 CALL UNETA ( UN, ETAB2, NS2, NR2, LSYM(KS), ISTP, LSPA )
!
      IF ( LSPA .NE. 1 ) GO TO 60
!
      DO 55 JS=1,NS2
!
      SETA = SQRT( 1.0 - ETAB2(JS)**2 )
      DO 55 JNU=1,NR2
   55 UN(JNU,JS) = UN(JNU,JS)*SETA
!
!
   60 CALL WEIGHT ( SY2(1), ETABAR(1,KS), NS3, A, A(71) )
!
      IF ( NS3 .NE. NS2 )                                               &
     &CALL WEIGHT ( SY2(NS3+1), ETABAR(NS3+1,KS), NS3, A, A(71) )
!
      DO 65 JS=1,NS2
      DO 65 JNU=1,NR2
   65 UN(JNU,JS) = UN(JNU,JS)*SY2(JS)
!
  100 RETURN
      END Subroutine GeomTU

      SUBROUTINE TKERNU ( CK, XMACH, JJ, NS2, FK, X, Y, Z, COSP, SINP,  &
     &                    SIG, XI, ETA, ZETA, COSQ, SINQ )
!
!     SUBROUTINE TO COMPUTE THE PLANAR OR NONPLANAR SUBSONIC OR
!     SUPERSONIC KERNEL FUNCTION FOR UNSTEADY TRANSONIC FLOW.
!
      DIMENSION CK(JJ), XI(JJ)
      COMPLEX   CK, CFK, CK1, CK2, CK3, CK4, FK1, EK1, EK2
!
      BETA2= 1.0 - XMACH*XMACH
      CFK  = CMPLX( 0.0, -FK )
      IF ( ABS(Z-ZETA) .GT. ABS(Y-ETA)*0.01 ) GO TO 100
      IF ( ABS(SINP-SINQ) .GT. ABS(COSP+COSQ)*0.005 ) GO TO 100
!
      Y0 = Y - ETA
      R2 = Y0*Y0
      R  = ABS( Y0 )
      FK1 = CFK*R
      IF ( XMACH .LT. 1.0 ) GO TO 20
!
      DO 10 J=1,JJ
      X0  = X - XI(J)
      RR2 = 1.0/( X0*X0 + BETA2*R2 )
      RR  = SQRT( RR2 )
      U1  = ( XMACH/RR - X0 )/( BETA2*R )
      U2  = (-XMACH/RR - X0 )/( BETA2*R )
!
      X1  = X0*RR
      EK1 = CEXP( FK1*U1 )
      EK2 = CEXP( FK1*U2 )
!
      CALL SPRI1 ( U1, U2, -FK1, CK1, CK2, EK1, EK2 )
!
      CK1 = -CK1 + (X1+1.0)*EK1 + CK2 + (X1-1.0)*EK2
!
   10 CK(J)   = -CK1*CEXP( CFK*X0 )/R2
      GO TO 1000
!
   20 DO 30 J=1,JJ
      X0  = X - XI(J)
      RR2 = 1.0/( X0*X0 + BETA2*R2 )
      RR  = SQRT( RR2 )
      U1  = ( XMACH/RR - X0 )/( BETA2*R )
!
      X1  = 1.0 + X0*RR
      EK1 = CEXP( FK1*U1 )
!
      CALL SUBI1 ( U1, FK*R, CK1 )
!
   30 CK(J)   = -( X1*EK1 - CK1 )*CEXP( CFK*X0 )/R2
      GO TO 1000
!
  100 T1 = COSP*COSQ + SINP*SINQ
!
      Y0 = Y - ETA
      Z0 = Z - ZETA
      T2 = ( Z0*COSP - Y0*SINP )*( Z0*COSQ - Y0*SINQ )
      R2 = Y0*Y0 + Z0*Z0
      OR2 = 1.0/R2
      R   = SQRT( R2 )
      FK1 = CFK*R
!
      IF ( XMACH .LT. 1.0 ) GO TO 150
!
      DO 110 J=1,JJ
      X0  = X - XI(J)
      RR2 = 1.0/( X0*X0 + BETA2*R2 )
      RR  = SQRT( RR2 )
      U1  = ( XMACH/RR - X0 )/( BETA2*R )
      U2  = (-XMACH/RR - X0 )/( BETA2*R )
!
      X1  = X0*RR
      EK1 = CEXP( FK1*U1 )
      EK2 = CEXP( FK1*U2 )
!
      CALL SPRI2 ( U1, U2, -FK1, CK1, CK2, CK3, CK4, EK1, EK2 )
!
      CK3 = CK3 - ( 2*(X1+1.0) + R2*RR2*( BETA2*X1 + FK1*(1.0-BETA2)/   &
     &                           SQRT( 1.0 + U1*U1 ) ) )*EK1            &
     &    - CK4 - ( 2*(X1-1.0) + R2*RR2*( BETA2*X1 - FK1*(1.0-BETA2)/   &
     &                           SQRT( 1.0 + U2*U2 ) ) )*EK2
!
      CK3 = CK3*OR2
!
      CK1 = -CK1 + (X1+1.0)*EK1 + CK2 + (X1-1.0)*EK2
!
  110 CK(J) = -OR2*( CK1*T1 + CK3*T2 )*CEXP( CFK*X0 )
      GO TO 1000
!
  150 DO 160 J=1,JJ

      X0  = X - XI(J)
      RR2 = 1.0/( X0*X0 + BETA2*R2 )
      RR  = SQRT( RR2 )
      U1  = ( XMACH/RR - X0 )/( BETA2*R )
!
      X1  = 1.0 + X0*RR
      EK1 = CEXP( FK1*U1 )
!
      CALL SUBI2 ( U1, FK*R, CK1, CK2 )
!
      CK2 = CK2*OR2 - ( 2*X1*OR2 + RR2*( BETA2*X0*RR                    &
     &                - FK1*U1*XMACH*XMACH/SQRT( 1.0 + U1*U1 ) ) )*EK1
!
      CK1 = X1*EK1 - CK1
!
  160 CK(J) = -OR2*( T1*CK1 + T2*CK2 )*CEXP( CFK*X0 )
!
 1000 RETURN
      END Subroutine Tkernu

      SUBROUTINE SUBI1 ( U1, RK1, CK1 )
!
!     SUBROUTINE TO COMPUTE THE I1' INTEGRAL BY THE USE OF
!     LASCHKA'S EXPONENTIAL APPROXIMATION.
!
      COMPLEX  CK1, C2, C3
      COMMON / COEF / AN(11), BN(11)
!
      C1  = -0.372*U1
      C2  = CMPLX( 0.0, RK1 )
      C3  = CEXP ( CMPLX( 0.0,-RK1*U1 ) )
      C4 = 0.0
      A1 = ABS(C1)
      IF ( A1 .LT. 50.0 ) C4 = EXP( -A1 )
!
      IF ( U1 .LT. 0.0 ) GO TO 20
      CK1 = 0.0
      C5  = 1.0
      DO 10 I=1,11
      C5  = C5*C4
   10 CK1 = CK1 + AN(I)*C5*( C2/CMPLX( I*0.372, RK1 ) - 1.0 )
!
      CK1 = ( CK1 + 1.0 - U1/SQRT( 1.0 + U1*U1 ) )*C3
      GO TO 40
!
   20 CK1 = -2.0 + ( 1.0 - U1/SQRT( 1.0 + U1*U1 ) )*C3
      RK1S = RK1*RK1
      C3  = C3
      DO 30 I=1,11
      C3  = C3*C4
   30 CK1 = CK1 + AN(I)*( C3*(  C2/CMPLX(I*0.372,-RK1) + 1.0  )         &
     &                  + 2.0*RK1S/( (I*0.372)**2 + RK1S ) )
!
   40 RETURN
      END Subroutine SubI1

      SUBROUTINE SUBI2  ( U1, RK1, CK1, CK2 )
!
!     SUBROUTINE TO COMPUTE BOTH THE I1' AND I22' INTEGRALS USING
!     EXPONENTIAL APPROXIMATIONS.
!
      COMPLEX CK1, CK2, C2, C3
      COMMON / COEF / AN(11), BN(11)
!
      C1  = -0.372*U1
      C2  = CMPLX( 0.0, RK1 )
      C3  = CEXP ( CMPLX( 0.0,-RK1*U1 ) )
      C4 = 0.0
      A1 = ABS(C1)
      IF ( A1 .LT. 50.0 ) C4 = EXP( -A1 )
      C6  = C4*C4
!
      IF ( U1 .LT. 0.0 ) GO TO 20
      CK1 = 0.0
      CK2 = 0.0
      C5  = 1.0
      C7  = 1.0
      DO 10 I=1,11
      C5  = C5*C4
      C7  = C7*C6
      CK1 = CK1 + ( AN(I)*C5 )/CMPLX( I*0.372, RK1 )
   10 CK2 = CK2 + ( BN(I)*C7 )/CMPLX( I*0.744, RK1 )
!
      CK1 = CK1*C3*C2
      CK2 = CK2*C3*C2
      GO TO 40
!
   20 CK1 = 2.0*( C3 - 1.0 )
      CK2 = CK1
      RK1S = RK1*RK1
      C2  = C3*C2
      C3  = C2
      DO 30 I=1,11
      C2  = C2*C4
      C3  = C3*C6
      CK1 = CK1 + AN(I)*( C2/CMPLX( I*0.372,-RK1 )                      &
     &                  + 2.0*RK1S/( (I*0.372)**2 + RK1S ) )
   30 CK2 = CK2 + BN(I)*( C3/CMPLX( I*0.744,-RK1 )                      &
     &                  + 2.0*RK1S/( (I*0.744)**2 + RK1S ) )
!
   40 CK2 = 3*CK1 - CK2
      RETURN
      END Subroutine SubI2

      SUBROUTINE SPRI1 ( U1, U2, FK1, CK1, CK2, EK1, EK2 )
!
!     SUBROUTINE TO COMPUTE THE I1' INTEGRAL FOR SUPERSONIC
!     FLOW USING LASCHKA'S APPROXIMATION.
!
      COMPLEX FK1, CK1, CK2, EK1, EK2, C7
!
      COMMON /COEF/ AN(11), BN(11)
!
      C11 = -0.372*ABS( U1 )
      C12 = -0.372*ABS( U2 )
      C2  = AIMAG( FK1 )
      C41 = 0.0
      IF ( C11 .GT. -50.0 ) C41 = EXP( C11 )
      C42 = 0.0
      IF ( C12 .GT. -50.0 ) C42 = EXP( C12 )
!
      IF ( U1 .LT. 0.0 ) GO TO 20
!
      CK1 = 0.0
      CK2 = 0.0
      C5  = 1.0
      C6  = 1.0
      DO 10 I=1,11
      C5  = C5*C41
      C6  = C6*C42
      C7  = AN(I)/CMPLX( I*0.372, C2 )
!
      CK1 = CK1 + C5*C7
   10 CK2 = CK2 + C6*C7
!
      CK1 = CK1*EK1*FK1
      CK2 = CK2*EK2*FK1
      GO TO 40
!
   20 CK1 = 0.0
      CK2 = 0.0
      C3  = C2*C2
      C5  = 2.0*C3
      C6  = 1.0
      C7  = EK1*FK1
      DO 30 I=1,11
      C7  = C7*C41
      C6  = C6*C42
!
      CK1 = CK1 + AN(I)*( C7/CMPLX( I*0.372,-C2 )                       &
     &                  + C5/( (I*0.372)**2 + C3 ) )
!
   30 CK2 = CK2 + AN(I)*( C6/CMPLX( I*0.372, C2 ) )
!
!
      CK1 = CK1 + 2.0*( EK1 - 1.0 )
      CK2 = CK2*EK2*FK1
!
   40 RETURN
      END Subroutine SprI1

      SUBROUTINE SPRI2 ( U1, U2, FK1, CK1, CK2, CK3, CK4, EK1, EK2 )
!
!     SUBROUTINE TO COMPUTE THE I1' AND I2' INTEGRALS FOR SUPERSONIC
!     FLOW USING LASCHKA'S AND CUNNINGHAM'S APPROXIMATIONS RESP.
!
      COMPLEX FK1, CK1, CK2, CK3, CK4, EK1, EK2, C7, C8
!
      COMMON /COEF/ AN(11), BN(11)
!
      C11 = -0.372*ABS( U1 )
      C12 = -0.372*ABS( U2 )
      C2  = AIMAG( FK1 )
      C41 = 0.0
      IF ( C11 .GT. -50.0 ) C41 = EXP( C11 )
      C42 = 0.0
      IF ( C12 .GT. -50.0 ) C42 = EXP( C12 )
!
      IF ( U1 .LT. 0.0 ) GO TO 20
!
      CK1 = 0.0
      CK2 = 0.0
      CK3 = 0.0
      CK4 = 0.0
      C5  = 1.0
      C6  = 1.0
      DO 10 I=1,11
      C5  = C5*C41
      C6  = C6*C42
      C7  = AN(I)/CMPLX( I*0.372, C2 )
      C8  = BN(I)/CMPLX( I*0.744, C2 )
!
      CK1 = CK1 + C5*C7
      CK2 = CK2 + C6*C7
!
      CK3 = CK3 + C5*C5*C8
   10 CK4 = CK4 + C6*C6*C8
!
      CK1 = CK1*EK1*FK1
      CK2 = CK2*EK2*FK1
      CK3 = CK3*EK1*FK1
      CK4 = CK4*EK2*FK1
      GO TO 40
!
   20 CK1 = 0.0
      CK2 = 0.0
      CK3 = 0.0
      CK4 = 0.0
      C3  = C2*C2
      C5  = 2.0*C3
      C6  = 1.0
      C7  = EK1*FK1
      C8  = C7
      DO 30 I=1,11
      C6  = C6*C42
      C7  = C7*C41
      C8  = C8*C41*C41
!
      CK1 = CK1 + AN(I)*( C7/CMPLX( I*0.372,-C2 )                       &
     &                  + C5/( (I*0.372)**2 + C3 ) )
!
      CK2 = CK2 + AN(I)*( C6/CMPLX( I*0.372,  C2 ) )
!
      CK3 = CK3 + BN(I)*( C8/CMPLX( I*0.744,-C2 )                       &
     &                  + C5/( (I*0.744)**2 + C3 ) )
!
   30 CK4 = CK4 + BN(I)*( C6*C6/CMPLX( I*0.744, C2 ) )
!
      CK1 = CK1 + 2.0*( EK1 - 1.0 )
      CK2 = CK2*EK2*FK1
      CK3 = CK3 + 2.0*( EK1 - 1.0 )
      CK4 = CK4*EK2*FK1
!
   40 CK3 = 3*CK1 - CK3
      CK4 = 3*CK2 - CK4
      RETURN
      END Subroutine SprI2

      BLOCK DATA
!
      COMMON / COEF / AN(11), BN(11)
!
      DATA AN / +.2418620 , -2.791803 , +24.99108 , -111.5920 ,         &
     &          +271.4355 , -305.7529 , -41.18363 , +545.9854 ,         &
     &          -644.7816 , +328.7276 , -64.27951 /
!
      DATA BN / +3.509407, -57.17120,+624.7548,-3830.151,+14538.51,     &
     &          -35718.32, +57824.14,-61303.92,+40969.58,-15660.04,     &
     &          +2610.093  /
!
      END Block Data


!!!      OVERLAY(R2T,4,0)
      SUBROUTINE Aeropr()
      COMMON /MANE/ DUM1(8), FREQ, DUM2(50), RK(50), NK, DUM3(1128)
      COMMON /OPTION/ IOP1, IOP2, IOP3, IOP4, NUNIT, IQT, IOPLU, IND

!!!      EQUIVALENCE ( NAME, STYCOF, UNSYCO, QGENF )
!!!      DATA NAME /6HCOMPAE/
!
      REWIND 1
      IF ( RK(1) .GT. FREQ ) GO TO 10
!!!      CALL OVERLAY ( STYCOF, 4, 1 )
      CALL Stycof()
      IF ( NK .EQ. 1 ) GO TO 30
   10 IF ( IOP2 .GE. 0 ) GO TO 20
!!!      CALL OVERLAY ( UNSYCO, 4, 2 )
      CALL Unsyco()
      GO TO 30
!!!   20 CALL OVERLAY ( QGENF,  4, 3 )
   20 CALL Qgenf()
   30 REWIND 1
      END Subroutine AeroPr

      SUBROUTINE MATINV ( A, N, NROWS )
      DIMENSION A(NROWS,NROWS)
      COMPLEX PIVOT, A
      N1=N+1
      DO 1 I=1,N
      DO 1 J=1,N
      J1=N+2-J
      J2=J1-1
    1 A(I,J1)=A(I,J2)
      J=1
      I=0
    2 I=I+1
      J1=J+1
      ID=0
      NN=0
    3 IF(REAL(A(I,J1)))5,4,5
    4 JJ=J1+101
      IF(AIMAG(A(I,J1)))5,6,5
    5 IF (ID) 6,10,6
    6 IF (I-N) 8,7,8
    7 I=0
    8 I=I+1
      NN=NN+1
      IF (N-NN) 20,20,3
   10 ID=J
      DO 101 L=1,N
  101 A(L,J)=CMPLX(0.0,0.0)
      A(I,J)=CMPLX(1.0,0.0)
      PIVOT=A(I,J1)
      DO 11 M=1,N1
   11 A(I,M)=A(I,M)/PIVOT
      DO 17 L=1,N
      IF (I-L) 12,12,15
   12 KK=L+1
      IF (N-KK) 18,13,13
   13 PIVOT=A(KK,J1)
      DO 14 M=1,N1
   14 A(KK,M)=A(KK,M)-PIVOT*A(I,M)
      GO TO 17
   15 PIVOT=A(L,J1)
      DO 16 M=1,N1
   16 A(L,M)=A(L,M)-PIVOT*A(I,M)
   17 END DO
   18 IF (N-J) 20,20,19
   19 J=J1
      GO TO 2
   20 RETURN
      END Subroutine Matinv

      SUBROUTINE AECOEF ( AC, D, NW2, NS, NSURF, KSURF, NMODES, RK,     &
     & BREF )
!
!     SUBROUTINE TO CONVERT THE DOWNWASH VECTORS STORED IN ARRAY 'AC'
!     TO PRESSURE SERIES COEFICIENT VECTORS IN THE SAME ARRAY.
!       THE INCOMING AERO MATRICES ON NUNIT MAY BE EITHER INVERTED
!     OR UNINVERTED ACCORDING TO THE OPTIONS IOP1 AND IOPLU.
!
      COMMON /OPTION/ IOP1, IOP2, IOP3, IOP4, NUNIT, IQT, IOPLU, IND
!
      COMMON /MATRIK/ C(70,70), ALPRS(2,70,20)
!
      DIMENSION AC(70,20), D(100), NW2(10), NS(10), KSURF(10,10)
!
      COMPLEX C, AC, D
      INTEGER:: errCode
!----------------------------------------------------------------------------
      IND = 0
!
      IF ( IOP1 .GT. 0 ) GO TO 50
      IF ( IOPLU .NE. 0 ) WRITE (NUNIT) RK, D
      I2 = 0
      DO 20 IS=1,NSURF
      I1 = I2 + 1
      I2 = I2 + NW2(IS)
      IF ( KSURF(IS,IS) .LT. 0 ) I2 = I2 + NS(IS)
!
      M2 = 0
      DO 10 KS=1,NSURF
      M1 = M2 + 1
      M2 = M2 + NW2(KS)
      IF ( KSURF(KS,KS) .LT. 0 ) M2 = M2 + NS(KS)
      DO 10 IT=I1,I2
   10 READ  (1) ( C(I,IT), I=M1,M2 )
!
      IF ( IOPLU .EQ. 0 ) GO TO 20
!
      DO 15 IT=I1,I2
   15 WRITE (NUNIT) ( C(I,IT), I=1,M2 )
!
   20 END DO
      GO TO 100
!
   50 M2 = 0
      DO 60 IS=1,NSURF
      M2 = M2 + NW2(IS)
      IF ( KSURF(IS,IS) .LT. 0 ) M2 = M2 + NS(IS)
   60 END DO
!
   70 READ (NUNIT,IOSTAT=errCode) RK3, D
!!!   71 IF ( EOF(NUNIT) .NE. 0 ) GO TO 999
   71 IF (errCode < 0) GO TO 999

   72 CONTINUE
      IF ( ABS( (RK3-RK)/RK ) .LT. 0.001 ) GO TO 80
      DO 75 IT=1,M2
   75 READ (NUNIT)
      GO TO 70
!
   80 DO 90 IT=1,M2
   90 READ (NUNIT) ( C(I,IT), I=1,M2 )
!
      IF ( IOPLU .EQ. 0 ) GO TO 150
!
  100 CALL MATINV ( C, M2, 70 )
!
      IF ( IOP1 .GE. 0 ) GO TO 150
      IF ( IOPLU .NE. 0 ) GO TO 150
!
      WRITE (NUNIT) RK, D
      DO 110 IT=1,M2
  110 WRITE (NUNIT) ( C(I,IT), I=1,M2 )
!
  150 DO 180 JA=1,NMODES
      DO 160 IT=1,M2
      D(IT) = ( 0.0, 0.0 )
      DO 160 I=1,M2
  160 D(IT) = D(IT) + C(I,IT)*CMPLX( ALPRS(1,I,JA), RK*ALPRS(2,I,JA)/   &
     &        BREF )
      DO 180 IT=1,M2
  180 AC(IT,JA) = D(IT)
!
      GO TO 1000
!
  999 IND = 1
      WRITE (6,900) NUNIT, RK
  900 FORMAT (///,"  AN UNEXPECTED END OF FILE WAS ENCOUNTERED ON ",    &
     &         //,"  THE AERODYNAMIC MATRIX TAPE ON UNIT ",I2,          &
     &         //,"  WHILE SEARCHING FOR THE MATRIX FOR RK = ",E12.4 )
!
 1000 RETURN
      END Subroutine AeCoef

!!!      OVERLAY(COMPAE,4,1)
      SUBROUTINE Stycof()
!
!     SUBROUTINE TO COMPUTE THE STEADY AERODYNAMIC CHARACTERISTICS.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), HMACH
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), NBREAK(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10),XIBAR(15,10), XIS(15,31,10), YBAR(15,10),          &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
!
      COMMON /MATRIK/ C(100,100), ALPRS(100,20), CP(50), XP(50), G(50), &
     &                E(50), D(100), F(300)
!
      COMMON /ALPVCT/ COEF(20,200), XF(200), YF(200), NFS(10), NFS1(10),&
     &   GMASS(20,20), NSTRS, NMODES, DH, DW1, DW2, SPAN2(10),          &
     &   TAN12(10), TAN22(10), CAVG2(10), XR(10), YR(10)
!
      COMMON /COMM/ DUM(12)
      DIMENSION AC(100,20)
!
      IF ( DUM(12) .NE. 0.0 ) NALP = NMODES
      REWIND 1
!
      I2 = 0
      DO 50 IS=1,NSURF
      I1  = I2 + 1
      I2  = I2 + NW2(IS)
      IF ( DUM(12) .EQ. 0.0 ) GO TO 10
!
      CALL CALMDS ( XS(1,1,IS), YS(1,IS), BREF, NC(IS), NS(IS),         &
     &              ALPRS(I1,1), A, B, +1, NST(IS), D, +10 )
      GO TO 20
!
   10 IF ( IDNWSH .NE. 0 ) GO TO 16
!
      DO 14 K=1,NALP
   14 READ (5,5) ( ALPRS(I,K), I=I1,I2 )
    5 FORMAT ( 6F10.0 )
      GO TO 20
!
   16 DO 18 K=1,NALP
      DO 18 I=I1,I2
   18 ALPRS(I,K) = 0.0
!
   20 DO 22 K=1,NALP
      DO 22 I=I1,I2
   22 ALPRS(I,K) = ALPRS(I,K) + ALPO(IS)
!
      WRITE (6,24) IS
   24 FORMAT (1H1,17X," ALPHA DISTRIBUTIONS FOR SURFACE ", I2, // )
!
      DO 26 K=1,NALP
   26 WRITE (6,28) K, ( ALPRS(I,K), I=I1,I2 )
   28 FORMAT (/,"    ALPHA VECTOR NO. ",I2, //, ( 5E15.4 ) )
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 40
!
      I3 = I2 + 1
      I2 = I2 + NS(IS)
      DO 30 K=1,NALP
      DO 30 I=I3,I2
   30 ALPRS(I,K) = 0.0
!
   40 CONTINUE
   50 END DO
      M2 = I2
!
      CALL AECSTD ( NW2, NS, NSURF, KSURF, NALP, M2, RK(1), IA, IB, AC)
!
      DO 1500 JA=1,NALP
!
      IF ( IDUMP .NE. 0 ) WRITE (6,65) JA, (AC(I,JA), I=1,M2 )
   65 FORMAT (1H1,19X,32H THE PRESSURE COEFICIENTS FOR                  &
     &       ,//, 19X,20H       ALPHA VECTOR , I2, // (3X,5E15.7)  )
!
   70 FORMAT (1H1,19X,32H THE PRESSURE DISTRIBUTIONS FOR                &
     &       ,//, 19X,20H       ALPHA VECTOR , I2        )
      DO 75 I=11,24
   75 D(I) = 0.0
!
      I2 = 0
      DO 1000 IS=1,NSURF
      WRITE (6,70) JA
      IF ( JSUROP .NE. 0 ) GO TO 76
!
      ICH = ICHORD(IS)
      LSP = LSPAN(IS)
      IST = ISTYPE(IS)
      GO TO 77
!
   76 ICH = 1
      LSP = 5
      IST = 4
      IND = 0
!
   77 I1 = I2 + 1
      I2 = I2 + NW2(IS)
      IF ( KSURF(IS,IS) .LT. 0 ) I2 = I2 + NS(IS)
      JR = 0
      DO 80 I=I1,I2
      JR = JR + 1
   80 A(JR) = AC(I,JA)*PI
!
      WRITE (6,90) IS
   90 FORMAT (//, 4X,21H LIFTING SURFACE NO. ,I2,/  )
      XP(1) = -0.95
      DO 100 I=2,20
  100 XP(I) = -1.0 + 0.1*(I-1)
      XP(21) = +0.999
!
      NC2 = NC(IS)
      NS3 = IABS( NS(IS) )
      DO 170 IY=1,21
!
      B2 = SO(IS)
      IF ( ISTYPE(IS) .GT. 2 ) GO TO 101
!
      YP = 0.05*(IY-1)
      YPSO = SO(IS)*YP
      IF ( LSPAN(IS) .EQ. 1 ) B2 = B2*SQRT(1.0 - YP*YP)
      GO TO 102
!
  101 YP = 0.1*(IY-1) - 1.0
      YPSO = SO(IS)*0.05*(IY-1)
      IF ( LSPAN(IS) .EQ. 1 ) B2 = B2*SQRT(1.0 - YP*YP)
!
  102 CALL UNETA ( G, YP, 1, NS3, LSYM(IS), IST, LSP )
      IF ( ISTYPE(IS) .GT. 2 ) YP = 0.05*(IY-1)
      CALL XMBY ( XTE, BYP, YPSO, BRPT(1,1,IS), NEND(IS) )
      IF ( BYP .NE. 0 ) GO TO 110
!
      WRITE (6,105) YP
  105 FORMAT ( //, "  THE CHORD AT YP =",E15.7," IS ZERO",// )
      GO TO 170
!
  110 B2 = B2/BYP
      YPSO = 0.9999*YPSO
!
      IF ( JSUROP .NE. 0 ) CALL CHDTSS ( B, XP, XP, 21, ICHORD(IS),     &
     & LSPAN(IS), IND, 1, XMACH, BRPT(1,1,IS), YPSO, JSUROP )
!
      DO 150 IX=1,21
      MUNU = 0
      CALL TNXI ( E, XP(IX), 1, NC(IS), ICH )
      Q1   = SQRT( (1.0 - XP(IX))/(1.0 + XP(IX)) )
!
  130 D(1) = 0.0
      DO 140 JI=1,NC2
      D(2) = 0.0
      DO 135 JR=1,NS3
      MUNU = MUNU + 1
  135 D(2) = D(2) + G(JR)*A(MUNU)
  140 D(1) = D(1) + D(2)*E(JI)
!
      IF ( JSUROP .NE. 0 ) GO TO 145
!
      CP(IX) = 8.0*B2*D(1)*Q1
      GO TO 150
!
  145 CP(IX) = 8.0*D(1)*B(IX)
!
  150 END DO
!
      WRITE (6,160) YP
  160 FORMAT (1H ,/,20X," SPAN STATION, YP=",F8.5,/,                    &
     &            2(31H         XP       CP(XP,YP)    ),/ )
      WRITE (6,165) ( XP(IX), CP(IX), IX=1,21 )
  165 FORMAT ( 2(7X,F8.5,E16.7) )
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 170
!
!
      JR = NW2(IS)
      D(1) = 0.0
      DO 180 J=1,NS3
      JR = JR + 1
  180 D(1) = D(1) + G(J)*A(JR)
!
      D(1) = D(1)*8.0*(AR(IS)**2)*BREF/SO(IS)
      IF ( LSPAN(IS) .EQ. 1 ) D(1) = D(1)*SQRT( 1.0 - YP*YP )
      WRITE (6,190)  D(1)
  190 FORMAT (/,"  SHOCK LOAD = ", E16.7, / )
!
  170 END DO
!
!
!     CALCULATION OF THE INTEGRATED CHARACTERISTICS.
!
      D(1) = 0.0
      D(2) = 0.0
      D(3) = 0.0
      NS2  = IABS(NSI(IS))
      NC3  = IABS(NJ (IS))
      B2   = 2.0*PI/(2.0*NC3+1.0)
      IND  = 0
      KS = ISTYPE(IS)
      GO TO ( 300, 300, 320, 330 ), KS
  300 JS2 = 1 + NS2/2
      JS3 = NS2
      GO TO 380
  320 JS2 = 1
      JS3 = NS2/2
      GO TO 380
  330 JS2 = 1
      JS3 = NS2
!
  380 KS = 0
      DO 500 JS=JS2,JS3
      KS = KS + 1
      S1 = SQRT(1.0-ETABAR(JS,IS)**2)
      S2 = XISLTE(1,JS,IS) - ( XV(IS) - BO(IS)*BRINV )
      CALL UNETA ( G, ETABAR(JS,IS), 1, NS3, LSYM(IS), IST, LSP )
!
      YPSO = SO(IS)*ETABAR(JS,IS)
      IF ( ISTYPE(IS) .GT. 1 ) YPSO = 0.5*( YPSO + SO(IS) )
      IF ( JSUROP .NE. 0 ) CALL CHDTSS ( F, XIBAR(1,IS), XIBAR(1,IS),   &
     &  NC3, ICHORD(IS), LSPAN(IS), IND, 0, XMACH, BRPT(1,1,IS),        &
     &   YPSO, JSUROP )
!
      D(4) = 0.0
      D(5) = 0.0
!
      DO 450 JJ=1,NC3
      CALL TNXI ( E, XIBAR(JJ,IS), 1, NC(IS), ICH )
      Q1 = F(JJ)
      IF ( JSUROP .EQ. 0 ) Q1 = 1.0 - XIBAR(JJ,IS)
      MUNU = 0
      D(7) = 0.0
      DO 400 JI=1,NC2
      D(6) = 0.0
!
      DO 390 JR=1,NS3
      MUNU = MUNU + 1
  390 D(6) = D(6) + G(JR)*A(MUNU)
  400 D(7) = D(7) + D(6)*E(JI)
      D(7) = D(7)*Q1
      D(4) = D(4) + D(7)
  450 D(5) = D(5) + D(7)*XIBAR(JJ,IS)
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 490
!
      D(100) = 0.0
      JR = NW2(IS)
      DO 460 J=1,NS3
      JR = JR + 1
  460 D(100) = D(100) + G(J)*A(JR)
      D(100) = D(100)*2*AR(IS)*BREF/(B2*SO(IS))
      D(4) = D(4) + D(100)
      D(5) = D(5) - D(100)
!
  490 D(8) = D(4)*B2
      IF ( LSP .EQ. 1 ) D(8) = D(8)*S1
      IF ( JSUROP .NE. 0 ) D(8) = D(8)*BETABO(JS,IS)/SO(IS)
      B(KS) = D(8)*4.0*AR(IS)
      XP(KS)= ( D(5)/D(4) + 1.0)*0.5
      D(8)  = B(KS)*S1
      IF ( JS .EQ. JS2 .AND. ISTYPE(IS) .LT. 3 ) D(8) = D(8)*0.5
      D(1)  = D(1) + D(8)
      D(2) = D(2) + D(8)*ETABAR(JS,IS)
      D(3)  = D(3) + ( XP(KS)*2.0*BETABO(JS,IS) + S2 )*D(8)
!
  500 END DO
!
  510 CP(1) = D(1)*PI/JS3
!
      IF ( ISTYPE(IS) .EQ. 3 ) CP(1) = CP(1)*0.5
!
!
  515 CP(2) = D(2)/D(1)
      IF ( ISTYPE(IS) .GT. 2 ) CP(2) = 0.5*( CP(2) + 1.0 )
      CP(3) = D(3)/D(1)
      K0 = KS - 1
      IF ( ISTYPE(IS) .NE. 1 ) K0 = 0
!
!     CP(1) =  LIFT COEFICIENT.
!     CP(2) =  SPANWISE CENTER OF PRESSURE, FRACTION OF SO(IS).
!     CP(3) =  CHORDWISE CENTER OF PRESSURE, FRACTION OF BO(IS) ABOUT
!              LEADING EDGE APEX.
!     B(JS) =  CLC/CAVG
!     XP(JS)=  LOCAL XCP IN THE TRANSFORMED COORDINATE SYSTEM.
!
      WRITE (6,520) IS
  520 FORMAT (1H1,18X," INTEGRATED AERODYNAMIC RESULTS", //             &
     &           ,18X,"       FOR SURFACE NO. ",I2,      //   )
      IF (LS .EQ. 1) GO TO 531
      WRITE (6,530) CP(1), CP(3), CP(2)
  530 FORMAT(14X,"SYMMETRIC LIFT COEFFICIENT=",ES12.5,//,                &
     &    11X, "XCPAVG = ",ES12.5, "  YCPAVG = ",ES12.5, //,              &
     &       13X,47H  ETA             CLC/CAVG             XCP(ETA),/)
      GO TO 533
  531 WRITE (6,532) CP(1), CP(3), CP(2)
  532 FORMAT (14X,"ANTI-SYMMETRIC LIFT COEFFICIENT = ", ES12.5, //,      &
     &    11X, "XCPAVG = ",ES12.5, "  YCPAVG = ",ES12.5, //,              &
     &       13X,47H  ETA             CLC/CAVG             XCP(ETA),/)
  533 WRITE (6,540) ( ETABAR(JS+K0,IS), B(JS), XP(JS), JS=1,KS )
  540 FORMAT(11X,F8.5, ES21.6,ES21.6)
!
      S2 = AREA(IS)
      IF ( ISTYPE(IS) .EQ. 4 ) S2 = 0.5*S2
      COSTH = COS( THETA(IS) )
      SINTH = SIN( THETA(IS) )
!
      IF ( ABS( COSTH ) .LT. 0.001 ) COSTH = 0.0
      IF ( ABS( SINTH ) .LT. 0.001 ) SINTH = 0.0
!
      D(11) = D(11) + S2*COSTH
      D(21) = D(21) + S2*SINTH
      D(1)  = D(1)*S2*PI/JS3
      IF ( ISTYPE(IS) .GT. 2 ) D(1) = D(1)*0.5
      D(12) = D(12) + D(1)*COSTH
      D(22) = D(22) + D(1)*SINTH
      CP(2) = SO(IS)*BRINV*D(1)*CP(2)
      D(13) = D(13) + CP(2)*COSTH + YV(IS)*D(1)
      D(23) = D(23) + CP(2)*SINTH + ZV(IS)*D(1)
      CP(3) = ( XV(IS) - BO(IS)*BRINV + CP(3) )*D(1)
      D(14) = D(14) + CP(3)*COSTH
      D(24) = D(24) + CP(3)*SINTH
!
 1000 END DO
!
      CP(1) = D(12)/D(11)
      CP(2) = D(13)/D(12)
      CP(3) = D(14)/D(12)
!
      WRITE (6,1010)  CP(1), CP(3), CP(2)
 1010 FORMAT (1H1,//,14X,"THE INTEGRATED CONFIGURATION CHARACTERISTICS",&
     &            //,14X,"    (A) VERTICAL COMPONENT, ALPHA " ,         &
     &            //,14X,"           CL  = ", ES14.6,                    &
     &            //,14X,"           XCP = ", ES14.6,                    &
     &            //,14X,"           YCP = ", ES14.6  )
!
      IF ( D(21) .EQ. 0.0 .OR. D(22) .EQ. 0.0 ) GO TO 1500
      CP(4) = D(22)/D(21)
      CP(5) = D(23)/D(22)
      CP(6) = D(24)/D(22)
!
!
      WRITE (6,1020) CP(4), CP(6), CP(5)
 1020 FORMAT (1H ,//,14X,"    (B) LATERAL COMPONENT, BETA ",            &
     &            //,14X,"           CL  = ",ES14.6,                     &
     &            //,14X,"           XCP = ",ES14.6,                     &
     &            //,14X,"           ZCP = ",ES14.6  )
!
 1500 END DO
!
 2000 CONTINUE
      END Subroutine StyCof

      SUBROUTINE AECSTD ( NW2, NS, NSURF, KSURF, NMODES, M2, RK, IA, IB,&
     &                    AC )
!
!
!     SUBROUTINE TO CONVERT THE DOWNWASH VECTORS STORED IN ARRAY ALPRS
!     TO PRESSURE SERIES COEFICIENT VECTORS IN THE ARRAY AC.
!     *** THE INCOMING AERO MATRICES ON NUNIT MAY BE EITHER INVERTED
!     OR UNINVERTED ACCORDING TO THE OPTIONS IOP1 AND IOPLU.
!
      COMMON /OPTION/ IOP1, IOP2, IOP3, IOP4, NUNIT, IQT, IOPLU, IND
!
      COMMON /MATRIK/ C(100,100), ALPRS(100,20), CP(50), XP(50), G(50), &
     &                E(50), D(100), F(300)
!
      DIMENSION NW2(10), NS(10), KSURF(10,10), IA(100), IB(100)
      DIMENSION AC(100,20)
      INTEGER:: errCode
!----------------------------------------------------------------------------
      IND = 0
!
      IF ( IOP1 .GT. 0 ) GO TO 70
      IF ( IOPLU .NE. 0 ) WRITE (NUNIT) RK, D
      I2 = 0
      DO 20 IS=1,NSURF
      I1 = I2 + 1
      I2 = I2 + NW2(IS)
      IF ( KSURF(IS,IS) .LT. 0 ) I2 = I2 + NS(IS)
!
      M2 = 0
      DO 10 KS=1,NSURF
      M1 = M2 + 1
      M2 = M2 + NW2(KS)
      IF ( KSURF(KS,KS) .LT. 0 ) M2 = M2 + NS(KS)
      DO 10 IT=I1,I2
   10 READ (1) ( C(I,IT), I=M1,M2 )
!
      IF ( IOPLU .EQ. 0 ) GO TO 20
!
      DO 15 IT=I1,I2
   15 WRITE (NUNIT) ( C(I,IT), I=1,M2 )
!
   20 END DO
      GO TO 100
!
   70 READ (NUNIT,IOSTAT=errCode) RK3, D
!!!   71 IF ( EOF(NUNIT) .NE. 0 ) GO TO 999
   71 IF (errCode < 0) GO TO 999

   72 CONTINUE
      IF ( ABS(RK3-RK) .LT. 0.01 ) GO TO 80
!
      WRITE (6,50) NUNIT, RK3, RK
   50 FORMAT (///,"  THE FIRST MATRIX ON THE INPUT TAPE ON UNIT ",I2,   &
     &         //,"  IS NOT A STEADY FLOW MATRIX. THE K VALUE ON"       &
     &         //,"  THE TAPE IS",E12.4,"  THE VALUE REQUESTED IS",E12.4&
     &  )
      IND = 1
      GO TO 1000
!
   80 DO 90 IT=1,M2
   90 READ (NUNIT) ( C(I,IT), I=1,M2 )
!
      IF ( IOPLU .EQ. 0 ) GO TO 150
!
  100 CALL INVRT ( C, D, IA, IB, M2, 100, 0 )
!
      IF ( IOP1 .GE. 0 ) GO TO 150
      IF ( IOPLU .NE. 0 ) GO TO 150
!
      WRITE (NUNIT) RK, D
      DO 110 IT=1,M2
  110 WRITE (NUNIT) ( C(I,IT), I=1,M2 )
!
  150 DO 180 JA=1,NMODES
      DO 160 IT=1,M2
      D(IT) = 0.0
      DO 160 I=1,M2
  160 D(IT) = D(IT) + C(I,IT)*ALPRS(I,JA)
      DO 180 IT=1,M2
  180 AC(IT,JA) = D(IT)
!
      GO TO 1000
!
  999 IND = 1
      WRITE (6,900) NUNIT, RK
  900 FORMAT (///,"  AN UNEXPECTED END OF FILE WAS ENCOUNTERED ON ",    &
     &         //,"  THE AERODYNAMIC MATRIX TAPE ON UNIT ",I2,          &
     &         //,"  WHILE SEARCHING FOR THE MATRIX FOR RK = ",E12.4 )
!
 1000 RETURN
      END Subroutine AecStd

!!!      OVERLAY(COMPAE,4,2)
      SUBROUTINE Unsyco()
!
!     SUBROUTINE TO COMPUTE THE STEADY AERODYNAMIC CHARACTERISTICS.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), HMACH
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), NBREAK(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10),XIBAR(15,10), XIS(15,31,10), YBAR(15,10),          &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
      COMMON /OPTION/ IOP1, IOP2, IOP3, IOP4, NUNIT, IQT, IOPLU, IND
!
      COMMON /AEROC/  AC(70,20), D(100), F(100), CP(21), G(30), E(10),  &
     &                XP(21)
      COMMON /ALPVCT/ COEF(20,200), XF(200), YF(200), NFS(10), NFS1(10),&
     &   GMASS(20,20), NSTRS, NMODES, DH, DW1, DW2, SPAN2(10),          &
     &   TAN12(10), TAN22(10), CAVG2(10), XR(10), YR(10)
!
      COMMON /MATRIK/ C(70,70), ALPRS(2,70,20)
!
      COMMON /COMMON/ DUM1(280)
      COMMON /COMM/ DUM(12)
      COMPLEX C, AC, D, CP, F
      COMPLEX XCP(70), CLC(70)
      EQUIVALENCE ( XCP(1), A(1) ), ( CLC(1), IA(1) )
!
      IF ( DUM(12) .NE. 0.0 ) NALP = NMODES
      I2 = 0
      DO 50 IS=1,NSURF
      I1 = I2 + 1
      I2 = I2 + NW2(IS)
      IF ( DUM(12) .EQ. 0.0 ) GO TO 10
!
      CALL CALMDS ( XS(1,1,IS), YS(1,IS), BREF, NC(IS), NS(IS),         &
     &              A, ALPRS(1,I1,1),B,+2, NST(IS), D, +10 )
      GO TO 20
!
   10 IF ( IDNWSH .NE. 0 ) GO TO 16
!
      DO 14 K=1,NALP
      READ (5,5) ( ALPRS(1,I,K), I=I1,I2 )
   14 READ (5,5) ( ALPRS(2,I,K), I=I1,I2 )
    5 FORMAT (6F10.0)
      GO TO 20
!
   16 DO 18 K=1,NALP
      DO 18 I=I1,I2
      ALPRS(1,I,K) = 0.0
   18 ALPRS(2,I,K) = 0.0
!
   20 DO 22 K=1,NALP
      DO 22 I=I1,I2
   22 ALPRS(1,I,K) = ALPRS(1,I,K) + ALPO(IS)
!
      WRITE (6,24) IS
   24 FORMAT (1H1,17X," DOWNWASH VECTORS FOR SURFACE ", I2, // )
!
      DO 26 K=1,NALP
      WRITE (6,27) K, ( ALPRS(1,I,K), I=I1,I2 )
   26 WRITE (6,28)    ( ALPRS(2,I,K), I=I1,I2 )
   27 FORMAT (/,"    SLOPES FOR MODE NO. ",I2, //, ( 5E15.4 ) )
   28 FORMAT (/,"    DEFLECTIONS", //, ( 5E15.4 ) )
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 40
!
      I3 = I2 + 1
      I2 = I2 + NS(IS)
      DO 30 K=1,NALP
      DO 30 I=I3,I2
      ALPRS(1,I,K) = 0.0
   30 ALPRS(2,I,K) = 0.0
!
   40 CONTINUE
   50 END DO
      M2 = I2
!
      REWIND 1
      DO 2000 IK=1,NK
      RK2 = RK(IK)
!
       CALL AECOEF ( AC, D, NW2, NS, NSURF, KSURF, NALP, RK2, BREF )
!
      IF ( IND .NE. 0 ) GO TO 3000
!
!
      DO 1500 JA=1,NALP
!
      IF ( IDUMP .NE. 0 ) WRITE (6,65) JA, (AC(I,JA), I=1,M2 )
   65 FORMAT (1H1,19X,32H THE PRESSURE COEFICIENTS FOR                  &
     &       ,//, 19X,20H       ALPHA VECTOR , I2, // (3X,5E15.7)  )
!
   70 FORMAT (1H1,19X,32H THE PRESSURE DISTRIBUTIONS FOR                &
     &       ,//, 19X,20H       ALPHA VECTOR , I2        )
   66 FORMAT (1H ,/,19X,17H     FREQUENCY = ,E11.4 )
      D(11) = 0.0
      D(12) = 0.0
      D(13) = 0.0
      D(14) = 0.0
      D(21) = 0.0
      D(22) = 0.0
      D(23) = 0.0
      D(24) = 0.0
!
      I2 = 0
      DO 1000 IS=1,NSURF
      WRITE (6,70) JA
      WRITE (6,66) RK(IK)
      IF ( JSUROP .NE. 0 ) GO TO 76
!
      ICH = ICHORD(IS)
      LSP = LSPAN(IS)
      IST = ISTYPE(IS)
      GO TO 77
!
   76 ICH = 1
      LSP = 5
      IST = 4
      IND = 0
!
   77 I1 = I2 + 1
      I2 = I2 + NW2(IS)
      IF ( KSURF(IS,IS) .LT. 0 ) I2 = I2 + NS(IS)
      JR = 0
      DO 80 I=I1,I2
      JR = JR + 1
   80 F(JR) = AC(I,JA)*PI
!
      IF ( IOP2 .LT. -1 ) WRITE (6,67)
   67 FORMAT (/,11X,48H *** PRESSURE DATA IN AMPLITUDE/PHASE FORMAT ***)
!
      WRITE (6,90) IS
   90 FORMAT (//, 4X,21H LIFTING SURFACE NO. ,I2,/  )
      XP(1) = -0.95
      DO 100 I=2,21
  100 XP(I) = -1.0 + 0.1*(I-1)
      XP(21) = 0.999
!
      NC2 = NC(IS)
      NS3 = IABS( NS(IS) )
      DO 170 IY=1,21
!
      IF ( ISTYPE(IS) .GT. 2 ) GO TO 101
!
      YP = 0.05*(IY-1)
      YPSO = SO(IS)*YP
      GO TO 102
!
  101 YP = 0.1*(IY-1) - 1.0
      YPSO = SO(IS)*0.05*(IY-1)
!
  102 SYP = SQRT(1.0 - YP*YP)
!
      CALL UNETA ( G, YP, 1, NS3, LSYM(IS), IST, LSP )
      IF ( ISTYPE(IS) .GT. 2 ) YP = 0.05*(IY-1)
      CALL XMBY ( XTE, BYP, YPSO, BRPT(1,1,IS), NEND(IS) )
      IF ( BYP .NE. 0 ) GO TO 110
!
      WRITE (6,105) YP
  105 FORMAT ( //, "  THE CHORD AT YP =",E15.7," IS ZERO",// )
      GO TO 170
!
  110 B2 = SO(IS)/BYP
      IF (LSPAN(IS) .EQ. 1 ) B2 = B2*SYP
!
      YPSO = 0.9999*YPSO
!
      IF ( JSUROP .NE. 0 ) CALL CHDTSS ( B, XP, XP, 21, ICHORD(IS),     &
     &   LSPAN(IS), IND, 1, XMACH, BRPT(1,1,IS), YPSO, JSUROP )
!
      DO 150 IX=1,21
      MUNU = 0
      CALL TNXI ( E, XP(IX), 1, NC(IS), ICH )
      Q1   = SQRT( (1.0 - XP(IX))/(1.0 + XP(IX)) )
!
  130 D(1) = 0.0
      DO 140 JI=1,NC2
      D(2) = 0.0
      DO 135 JR=1,NS3
      MUNU = MUNU + 1
  135 D(2) = D(2) + G(JR)*F(MUNU)
  140 D(1) = D(1) + D(2)*E(JI)
!
      IF ( JSUROP .NE. 0 ) GO TO 145
!
      CP(IX) = 8.0*B2*D(1)*Q1
      GO TO 150
!
  145 CP(IX) = 8.0*D(1)*B(IX)
!
  150 END DO
!
      WRITE (6,160) YP
  160 FORMAT (1H ,/,25X,"SPAN STATION, YP=", F8.5, /,                   &
     & 2("      XP             CP(XP,YP)          "), / )
!
      IF ( IOP2 .LT. -1 ) CALL POLAR ( CP, 21 )
!
      WRITE (6,165) ( XP(IX), CP(IX), IX=1,21 )
  165 FORMAT (2(3X,F8.5,2E14.7) )
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 170
!
      JR = NW2(IS)
      D(1) = 0.0
      DO 180 J=1,NS3
      JR = JR + 1
  180 D(1) = D(1) + G(J)*F(JR)
      D(1) = D(1)*8.0*(AR(IS)**2)*BREF/SO(IS)
      IF ( LSPAN(IS) .EQ. 1 ) D(1) = D(1)*SQRT( 1.0 - YP*YP )
!
      IF ( IOP2 .LT. -1 ) CALL POLAR ( D(1), 1 )
!
      WRITE (6,190) D(1)
  190 FORMAT (/,"  SHOCK LOAD = ", 2E16.7, / )
!
  170 END DO
!
!
!     CALCULATION OF THE INTEGRATED CHARACTERISTICS
!
      D(1) = 0.0
      D(2) = 0.0
      D(3) = 0.0
      NS2  = IABS(NSI(IS))
      NC3  = IABS(NJ (IS))
      B2   = 2.0*PI/(2.0*NC3+1.0)
      IND = 0
      KS = ISTYPE(IS)
      GO TO ( 300, 300, 320, 330 ), KS
  300 JS2 = 1 + NS2/2
      JS3 = NS2
      GO TO 380
  320 JS2 = 1
      JS3 = NS2/2
      GO TO 380
  330 JS2 = 1
      JS3 = NS2
!
  380 KS = 0
      DO 500 JS=JS2,JS3
      KS = KS + 1
      S1 = SQRT(1.0-ETABAR(JS,IS)**2)
      S2 = XISLTE(1,JS,IS) - ( XV(IS) - BO(IS)*BRINV )
      CALL UNETA ( G, ETABAR(JS,IS), 1, NS3, LSYM(IS), IST, LSP )
!
      YPSO = SO(IS)*ETABAR(JS,IS)
      IF ( ISTYPE(IS) .GT. 1 ) YPSO = 0.5*( YPSO + SO(IS) )
      IF ( JSUROP .NE. 0 ) CALL CHDTSS ( B, XIBAR(1,IS), XIBAR(1,IS),   &
     &     NC3, ICHORD(IS), LSPAN(IS), IND, 0, XMACH, BRPT(1,1,IS),     &
     &     YPSO, JSUROP )
!
      D(4) = 0.0
      D(5) = 0.0
!
      DO 450 JJ=1,NC3
      CALL TNXI ( E, XIBAR(JJ,IS), 1, NC(IS), ICH )
      Q1 = B(JJ)
      IF ( JSUROP .EQ. 0 ) Q1 = 1.0 - XIBAR(JJ,IS)
      MUNU = 0
      D(7) = 0.0
      DO 400 JI=1,NC2
      D(6) = 0.0
!
      DO 390 JR=1,NS3
      MUNU = MUNU + 1
  390 D(6) = D(6) + G(JR)*F(MUNU)
  400 D(7) = D(7) + D(6)*E(JI)
      D(7) = D(7)*Q1
      D(4) = D(4) + D(7)
  450 D(5) = D(5) + D(7)*XIBAR(JJ,IS)
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 490
!
      D(100) = 0.0
      JR = NW2(IS)
      DO 460 J=1,NS3
      JR = JR + 1
  460 D(100) = D(100) + G(J)*F(JR)
      D(100) = D(100)*2*AR(IS)*BREF/(B2*SO(IS))
!
      D(4)   = D(4) + D(100)
      D(5)   = D(5) - D(100)
!
  490 D(8) = D(4)*B2
      IF ( LSP .EQ. 1 ) D(8) = D(8)*S1
      IF ( JSUROP .NE. 0 ) D(8) = D(8)*BETABO(JS,IS)/SO(IS)
      CLC(KS) = D(8)*4.0*AR(IS)
      IF ( CABS(D(4)) .NE. 0.0 ) GO TO 495
      XCP(KS) = 0.0
      GO TO 496
  495 XCP(KS) = 0.5*CMPLX ( ( REAL(D(5))/REAL(D(4)) + 1.0 ),            &
     &                      ( AIMAG(D(5))/AIMAG(D(4)) + 1.0 ) )
  496 D(8)    = CLC(KS)*S1
      IF ( JS .EQ. JS2 .AND. ISTYPE(IS) .LT. 3 ) D(8) = D(8)*0.5
      D(1)  =  D(1) + D(8)
      D(2) = D(2) + D(8)*ETABAR(JS,IS)
      XCPR = REAL(XCP(KS))*2.0*BETABO(JS,IS) + S2
      XCPI = AIMAG(XCP(KS))*2.0*BETABO(JS,IS) + S2
!
      D(3) = D(3) + CMPLX ( XCPR*REAL(D(8)), XCPI*AIMAG(D(8)) )
!
  500 END DO
!
  510 CP(1) = D(1)*PI/JS3
!
      IF ( ISTYPE(IS) .EQ. 3 ) CP(1) = CP(1)*0.5
!
!
  515 FRR = REAL(D(1))
      FRI = AIMAG(D(1))
!
      IF ( CABS(D(1)) .NE. 0.0 ) GO TO 516
      CP(2) = 0.0
      CP(3) = 0.0
      GO TO 517
  516 CP(2) = CMPLX ( REAL(D(2))/FRR, AIMAG(D(2))/FRI )
      IF ( ISTYPE(IS) .GT. 2 ) CP(2) = 0.5*( CP(2) + CMPLX(1.0,1.0) )
      CP(3) = CMPLX ( REAL(D(3))/FRR, AIMAG(D(3))/FRI )
  517 K0 = KS - 1
      IF ( ISTYPE(IS) .NE. 1 ) K0 = 0
!
!     CP(1) = LIFT COEFICIENT.
!     CP(2) = SPANWISE CENTER OF PRESSURE, FRACTION OF SO(IS)
!     CP(3) = CHORDWISE CENTER OF PRESSURE, FRACTION OF BO(IS) ABOUT
!             LEADING EDGE APEX.
!
      WRITE (6,520) IS, RK(IK)
  520 FORMAT (1H1,18X," INTEGRATED AERODYNAMIC RESULTS ",//             &
     &           ,18X,"       FOR SURFACE NO. ", I2,     //             &
     &           ,18X," REDUCED FREQUENCY = ", E11.4,    // )
      IF (LS .EQ. 1) GO TO 531
      WRITE (6,530) CP(1), CP(3), CP(2)
  530 FORMAT(14X,"SYMMETRIC LIFT COEFFICIENT=",2E14.6,/,                &
     &       14X,"XCPAVG=",2E14.6,/,14X,"YCPAVG=",2E14.6,//,            &
     &       13X,47H  ETA             CLC/CAVG             XCP(ETA),/)
      GO TO 533
  531 WRITE (6,532) CP(1), CP(3), CP(2)
  532 FORMAT(14X,"ANTI-SYMMETRIC LIFT COEFFICIENT=",2E14.6,/,           &
     &       14X,"XCPAVG=",2E14.6,/,14X,"YCPAVG=",2E14.6,//,            &
     &       13X,47H  ETA             CLC/CAVG             XCP(ETA),/)
  533 WRITE (6,540) ( ETABAR(JS+K0,IS), CLC(JS), XCP(JS), JS=1,KS )
  540 FORMAT(11X,F8.5,3E14.6,E14.6)
!
      S2 = AREA(IS)
      S1 = XV(IS) - BO(IS)*BRINV
      IF ( ISTYPE(IS) .EQ. 4 ) S2 = 0.5*S2
      COSTH = COS( THETA(IS) )
      SINTH = SIN( THETA(IS) )
!
      IF ( ABS( COSTH ) .LT. 0.001 ) COSTH = 0.0
      IF ( ABS( SINTH ) .LT. 0.001 ) SINTH = 0.0
!
      D(11) = D(11) + S2*COSTH
      D(21) = D(21) + S2*SINTH
      D(1)  = D(1)*S2*PI/JS3
      IF ( ISTYPE(IS) .GT. 2 ) D(1) = D(1)*0.5
      D(12) = D(12) + D(1)*COSTH
      D(22) = D(22) + D(1)*SINTH
!
      FRR = REAL(D(1))
      FRI = AIMAG(D(1))
!
      CP(2) = SO(IS)*BRINV*CMPLX( (  FRR*REAL(CP(2)) ),                 &
     &                            ( FRI*AIMAG(CP(2)) )  )
!
      D(13) = D(13) + CP(2)*COSTH + YV(IS)*D(1)
      D(23) = D(23) + CP(2)*SINTH + ZV(IS)*D(1)
!
      CP(3) = CMPLX ( ( S1 +  REAL(CP(3)) )*FRR,                        &
     &                ( S1 + AIMAG(CP(3)) )*FRI  )
!
      D(14) = D(14) + CP(3)*COSTH
      D(24) = D(24) + CP(3)*SINTH
!
 1000 END DO
!
      IF ( CABS(D(11)) .EQ. 0.0 .OR. CABS(D(12)) .EQ. 0.0 ) GO TO 1005
!
      CP(1) = D(12)/D(11)
!
      FRR = REAL(D(12))
      FRI = AIMAG(D(12))
!
      CP(2) = CMPLX ( REAL(D(13))/FRR, AIMAG(D(13))/FRI )
      CP(3) = CMPLX ( REAL(D(14))/FRR, AIMAG(D(14))/FRI )
      GO TO 1006
!
 1005 CP(1) = 0.0
      CP(2) = 0.0
      CP(3) = 0.0
!
 1006 WRITE (6,1010) CP(1), CP(3), CP(2)
 1010 FORMAT (1H1,//,14X,"THE INTEGRATED CONFIGURATION CHARACTERISTICS",&
     &            //,14X,"    (A) VERTICAL COMPONENT, ALPHA ",          &
     &            //,14X,"      CL  =", 2E15.7,                         &
     &            //,14X,"      XCP =", 2E15.7,                         &
     &            //,14X,"      YCP =", 2E15.7  )
!
      IF (  CABS(D(21)) .EQ. 0.0 .OR.  CABS(D(22)) .EQ. 0.0) GO TO 1500
!
      CP(4) = D(22)/D(21)
!
      FRR = REAL(D(22))
      FRI = AIMAG(D(22))
!
      CP(5) = CMPLX ( REAL(D(23))/FRR, AIMAG(D(23))/FRI )
      CP(6) = CMPLX ( REAL(D(24))/FRR, AIMAG(D(24))/FRI )
!
      WRITE (6,1020) CP(4), CP(6), CP(5)
 1020 FORMAT (1H ,//,14X,"    (B) LATERAL COMPONENT, BETA ",            &
     &            //,14X,"      CL  =", 2E15.7,                         &
     &            //,14X,"      XCP =", 2E15.7,                         &
     &            //,14X,"      YCP =", 2E15.7   )
!
 1500 END DO
!
 2000 END DO
!
!
 3000 CONTINUE
      END Subroutine Unsyco

      SUBROUTINE POLAR ( C, N )
      COMPLEX C(N)
      RADDEG = 45.0/ATAN(1.0)
!
      DO 10 I=1,N
      AMP =  CABS( C(I) )
      R = REAL( C(I) )
      A = AIMAG( C(I) )
      IF ( R .NE. 0.0 .AND. A .NE. 0.0 ) GO TO 5
      IF ( AMP .EQ. 0.0 ) PHASE = 0.0
      IF ( A .NE. 0.0 ) PHASE = SIGN( 90.0, A )
      IF ( R .NE. 0.0 ) PHASE = 90.0 - SIGN( 90.0, R )
      GO TO 10
    5 PHASE = RADDEG*ATAN2( A, R )
   10 C(I) = CMPLX( AMP, PHASE )
!
      RETURN
      END Subroutine Polar


!!!      OVERLAY(COMPAE,4,3)
      SUBROUTINE Qgenf()
!
!     SUBROUTINE TO CALCULATE GENERALIZED FORCES.
!
      COMMON /MANE/ BREF, XMACH, BETA, BETA2, BRINV, LS, NW, PI, FREQ,  &
     & REFREQ(50), RK(50), NK, MARSHA, NALP, IA(140), IB(140), NSURF,   &
     & TAREA, A(140), B(140), JSUROP, IDUMP, IDNWSH, ICHORD(10), IXI(10)&
     & , ISTYPE(10), LSYM(10), LSPAN(10), BMACH(50,10), NST(10), HMACH
      COMMON /MANE2/ BO(10), SO(10), XV(10), YV(10), ZV(10), THETA(10), &
     & ALPO(10), NBREAK(10), NEND(10), NBRPT(10), AR(10), ETAS(31,10),  &
     & ETABAR(31,10),XIBAR(15,10), XIS(15,31,10), YBAR(15,10),          &
     & YS(15,10), XBAR(10,10), XS(10,15,10), BYBO(15,10), NW2(10),      &
     & NC(10), NS(10), NJ(10), NSI(10), AREA(10), BOINV(10), KSURF(10,  &
     & 10), BRPT(2,40,10), XISLTE(2,31,10), BETABO(31,10), SIGY(15,10), &
     & SIGETA(31,10), ZY(15,10), ZETA(31,10)
      COMMON /OPTION/ IOP1, IOP2, IOP3, IOP4, NUNIT, KQT, IOPLU, IND
!
      COMMON /AEROC/  AC(70,20), D(100), F(100), CP(21), G(30), E(10),  &
     &                XP(21)
      COMMON /ALPVCT/ COEF(20,200), XF(200), YF(200), NFS(10), NFS1(10),&
     &   GMASS(20,20), NSTRS, NMODES, DH, DW1, DW2, SPAN2(10),          &
     &   TAN12(10), TAN22(10), CAVG2(10), XR(10), YR(10)
!
      COMMON /MATRIK/ QRS(20,20), HRS(100,20), DUM2(7000),              &
     &                ALPRS(2,70,20)
!
      COMMON /COMMON/ DUM1(280)
      COMMON /COMM  / DUM(12)
!
      CHARACTER(LEN=128):: tle
      COMMON /TITLE / TLE
!
      COMPLEX AC, D, CP, F, QRS
!
      IQT = IABS(KQT)
      IF ( IQT .EQ. 0 ) IQT = 4
      REWIND IQT
      NQ2 = NMODES
      NQ = NQ2
      IF ( IOP2 .GT. 0 ) NQ  = IOP2
!
      IF ( IDNWSH .GE. 0 ) GO TO 2
!
      NQ = 1
      READ (5,5) UGUST, VFS, CONST
      WRITE (6,3) UGUST, VFS, CONST
    3 FORMAT (1H1,//,17X,35H GUST QRF TERMS WILL BE CALCULATED  ,       &
     &            //,17X,28H      GUST VELOCITY (FPS) = , E14.7 ,       &
     &             /,17X,28H      FREE STREAM V (FPS) = , E14.7 ,       &
     &             /,17X,28H      MULT. CONST         = , E14.7   )
!
    2 IF (  KQT .EQ. 0 ) GO TO 1
!
      CALL STATUS (IA)
  950 REWIND IQT
      WRITE (IQT) IA(1), IA(3), ( IA(I), I=7,18 )
      WRITE (IQT) TLE
      WRITE (IQT) NK, NQ, NQ2, BREF, XMACH
      WRITE (IQT) ( IA(I), I=1,10 )
      WRITE (IQT) ( RK(I), I=1,NK )
!
!     NEXT THE DOWNWASH VECTORS WILL BE CONSTRUCTED.
!
    1 DO 2000 IK=1,NK
!
      IF ( IK .NE. 1 .AND. IDNWSH .GE. 0 ) GO TO 56
      I2 = 0
      DO 50 IS=1,NSURF
      I1 = I2 + 1
      I2 = I2 + NW2(IS)
!
      IF ( IDNWSH .GE. 0 ) GO TO 9
!
      C2 = CONST*UGUST/VFS
      NX = NC(IS)
      NY = NS(IS)
      I = I1 - 1
      DO 8 K=1,NY
      DO 8 J=1,NX
      I = I+1
      ALPRS(1,I,1) = C2*COS( XS(J,K,IS)*RK(IK) )
    8 ALPRS(2,I,1) =-C2*SIN( XS(J,K,IS)*RK(IK) )*BREF/RK(IK)
      GO TO 23
!
    9 IF ( DUM(12) .EQ. 0.0 ) GO TO 10
!
      IF ( THETA(IS) .EQ. 0.0 ) GO TO 13
!
      NY = NS(IS)
      DO 11 J=1,NY
   11 YS(J,IS) = ZY(J,IS)
      NS2 = NSI(IS)
      DO 12 J=1,NS2
   12 ETAS(J,IS) = ZETA(J,IS)
!
   13 IF ( IOP2 .GT. 0 ) GO TO 10
!
      CALL CALMDS ( XS(1,1,IS), YS(1,IS), BREF, NC(IS), NS(IS),         &
     &              A, ALPRS(1,I1,1), B, 2, NST(IS), D, +10   )
      GO TO 20
!
   10 IF ( IDNWSH .NE. 0 ) GO TO 16
!
      DO 14 K=1,NQ
      READ (5,5) ( ALPRS(1,I,K), I=I1,I2 )
   14 READ (5,5) ( ALPRS(2,I,K), I=I1,I2 )
    5 FORMAT ( 6F10.0 )
      GO TO 20
!
   16 DO 18 K=1,NQ
      DO 18 I=I1,I2
      ALPRS(1,I,K) = 0.0
   18 ALPRS(2,I,K) = 0.0
!
   20 DO 22 K=1,NQ
      DO 22 I=I1,I2
   22 ALPRS(1,I,K) = ALPRS(1,I,K) + ALPO(IS)
!
   23 WRITE (6,24) IS
   24 FORMAT (1H1,17X," DOWNWASH VECTORS FOR SURFACE ", I2, // )
!
      DO 26 K=1,NQ
      WRITE (6,27) K, ( ALPRS(1,I,K), I=I1,I2 )
   26 WRITE (6,28)    ( ALPRS(2,I,K), I=I1,I2 )
   27 FORMAT (/,"    SLOPES FOR MODE NO. ",I2, //, ( 8E15.4 ) )
   28 FORMAT (/,"    DEFLECTIONS", //,             ( 8E15.4 ) )
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 50
!
      I3 = I2 + 1
      I2 = I2 + NS(IS)
      DO 30 K=1,NQ
      DO 30 I=I3,I2
      ALPRS(1,I,K) = 0.0
   30 ALPRS(2,I,K) = 0.0
!
      WRITE (6,35) I3, I2
   35 FORMAT (/,"    ROWS ",I3," THROUGH ", I3," ARE SET TO ZERO FOR ", &
     &        /,"    NORMAL SHOCK BOUNDARY CONDITIONS. "  )
!
   50 END DO
!
      WRITE (6,55)
   55 FORMAT (1H1)
   56 RK2 = RK(IK)
!
!     NEXT, THE PRESSURE SERIES COEFICIENTS WILL BE CALCULATED
!     FOR THE NQ DOWNWASH VECTORS AT THE BEGINNING OF THE FREQUENCY
!     LOOP FOR NK FREQUENCIES.
!
      CALL AECOEF ( AC, D, NW2, NS, NSURF, KSURF, NQ, RK2, BREF )
      IF ( IND .NE. 0 ) GO TO 3000
!
!     THE GENERALIZED FORCES WILL BE CALCULATED IN ARRAY QRS AND
!     WRITTEN ON TAPE IQT.
!
      IF ( KQT .LT.0 ) GO TO 70
!
      DO 60 J=1,NQ
      DO 60 I=1,NQ2
   60 QRS(I,J) = ( 0.0, 0.0 )
      GO TO 90
!
   70 READ ( 4 ) ( ( QRS(I,J), I=1,NQ2 ), J=1,NQ )
!
   90 I2 = 0
      DO 1000 IS=1,NSURF
      I1 = I2 + 1
      I2 = I2 + NW2(IS)
      IF ( KSURF(IS,IS) .LT. 0 ) I2 = I2 + NS(IS)
      IF ( JSUROP .NE. 0 ) GO TO 94
!
      ICH = ICHORD(IS)
      LSP = LSPAN(IS)
      IST = ISTYPE(IS)
      GO TO 95
!
   94 ICH = 1
      LSP = 5
      IST = 4
!
   95 IND = 0
      MBAR = NC(IS)
      NR   = NS(IS)
      JJ   = NJ(IS)
      NS2  = NSI(IS)
      KS   = ISTYPE(IS)
      GO TO ( 100, 100, 120, 130 ), KS
  100 JS2 = 1 + NS2/2
      JS3 = NS2
      GO TO 140
  120 JS2 = 1
      JS3 = NS2/2
      GO TO 140
  130 JS2 = 1
      JS3 = NS2
!
  140 B2 = 8.0*( ( PI*SO(IS))**2 )/(JS3*BREF)
      IF ( ISTYPE(IS) .GT. 2 ) B2 = B2*0.5
      IF ( IDNWSH     .LT. 0 ) B2 = B2/12.0
!
      B1 = B2*PI*2.0/( 2.0*JJ+1.0 )
      B4 = B2*2.0*BREF*AR(IS)/SO(IS)
!
      IF ( ISTYPE(IS) .NE. 4 ) GO TO 190
      B1 = 0.5*B1
      B2 = 0.5*B2
!
  190 DO 500 JS=JS2,JS3
      S1 = SQRT( 1.0 - ETABAR(JS,IS)**2 )
      IF ( LSP .EQ. 1 ) S1 = S1*S1
      IF ( JSUROP .NE. 0 ) S1 =     ( S1 )*BETABO(JS,IS)/SO(IS)
      B3 = S1*B1
!
      CALL UNETA ( G, ETABAR(JS,IS), 1, NR, LSYM(IS), IST, LSP )
!
      IF ( JSUROP .EQ. 0 ) GO TO 385
      YPSO = SO(IS)*ETABAR(JS,IS)
      IF ( ISTYPE(IS) .GT. 1 ) YPSO = 0.5*( YPSO + SO(IS) )
      CALL CHDTSS ( B, XIBAR(1,IS), XIBAR(1,IS), NC3, ICHORD(IS),       &
     &  LSPAN(IS), IND, 0, XMACH, BRPT(1,1,IS), YPSO, JSUROP )
!
  385 CALL CALMDS ( XIS(1,JS,IS), ETAS(JS,IS), BREF, JJ, 1,             &
     &              A, B, HRS, +3, NST(IS), D, +15 )
!
      DO 400 J=1,JJ
      CALL TNXI ( E, XIBAR(J,IS), 1, MBAR, ICH )
      Q1 = B3*( 1.0 - XIBAR(J,IS) )
      IF ( JSUROP .NE. 0 ) Q1 = B3*B(JJ)
!
      IF ( JS .EQ. JS2 .AND. ISTYPE(IS) .LT. 3 ) Q1 = Q1*0.5
!
      DO 300 JQ=1,NQ
!
      D(1) = (0.0,0.0)
      MUNU = I1 - 1
      DO 210 JI=1,MBAR
      D(2) = (0.0,0.0)
      DO 200 JR=1,NR
      MUNU = MUNU + 1
  200 D(2) = D(2) + G(JR)*AC(MUNU,JQ)
  210 D(1) = D(1) + D(2)*E(JI)
      D(1) = D(1)*Q1
      DO 300 IQ=1,NQ2
  300 QRS(IQ,JQ) = QRS(IQ,JQ) + D(1)*HRS(J,IQ)
!
  400 END DO
!
      IF ( KSURF(IS,IS) .GE. 0 ) GO TO 500
!
      CALL CALMDS ( XISLTE(1,JS,IS), ETAS(JS,IS), BREF, 1, 1,           &
     &              A, B, HRS,  +3, NST(IS), D, +2 )
!
      Q1 = S1*B4
!
      IF ( JS .EQ. JS2 .AND. ISTYPE(IS) .LT. 3 ) Q1 = Q1*0.5
!
      DO 450 JQ=1,NQ
!
      D(1) = (0.0,0.0)
      I4 = MUNU
      DO 410 JR=1,NR
      I4 = I4 + 1
  410 D(1) = D(1) + G(JR)*AC(I4,JQ)
      D(1) = D(1)*Q1
!
      DO 450 IQ=1,NQ2
  450 QRS(IQ,JQ) = QRS(IQ,JQ) + D(1)*HRS(1,IQ)
!
  500 END DO
!
 1000 END DO
!
      WRITE (IQT) ( ( QRS(I,J), I=1,NQ2 ), J=1,NQ )
!
!
 1009 WRITE ( 6,1010 ) RK2
 1010 FORMAT (//,37H     THE GENERALIZED FORCES FOR K =  ,E11.4, / )
      DO 1020 J=1,NQ
 1020 WRITE ( 6,1030 ) J, ( QRS(I,J), I=1,NQ2 )
 1030 FORMAT (/, 25H         PRESSURE MODE =  , I2, //, ( 8E15.6 ) )
!
 2000 END DO
!
      END FILE IQT
      REWIND IQT
!
 3000 CONTINUE
      END Subroutine Qgenf
