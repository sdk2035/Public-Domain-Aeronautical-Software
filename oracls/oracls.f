C ORACLS LIBRARY
C modifications by PDAS:
C  21 Dec 2005 changed common block /LINES/ to be
C      INTEGER:: nlp,lin
C      CHARACTER(LEN=80):: title
C      COMMON/LINES/ nlp,lin,title
C  21 Dec 2005 changed common block /FORM/ to be
C      INTEGER:: nepr
C      CHARACTER(LEN=16):: fmt1,fmt2
C      COMMON/FORM/ nepr,fmt1,fmt2
C  15 Dec 2009 changed common block /CONV/ 
C   from
C      COMMON/CONV/SUMCV,MAXSUM,RICTCV,SERCV
C   to
C      COMMON/CONV/SUMCV,RICTCV,SERCV,MAXSUM
C   to avoid messages that padding is needed (for alignment of 8-byte words)

C   15 Dec 2009 changed many (all?) dimension of (1) to (*) in dummy arguments

C   16 Dec 2009 changed all calls to PRNT that had integer zero as third
C     argument to have '    ' as third argument

      SUBROUTINE ADD (A,NA,B,NB,C,NC)
C 
C   PURPOSE:
C      Perform matrix addition C = A + B for given matrices A and B.
C 
C   Subroutines employed by ADD: LNCNT
C   Subroutines employing ADD: ASYREG, CNTREG, DISREG, DSTAB, EXMDFL,
C      EXPINT, EXPSER, IMMDFL, RICNWT, SAMPL, SUM, TRNSIT
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),C(*),NA(2),NB(2),NC(2)
      IF( (NA(1) .NE. NB(1)) .OR. (NA(2) .NE. NB(2)) )  GO TO 999
      NC(1)=NA(1)
      NC(2)=NA(2)
      L=NA(1)*NA(2)
      IF( NA(1) .LT. 1  .OR.  L .LT. 1 )  GO TO 999
      DO 300 I=1,L
  300 C(I)=A(I)+B(I)
      GO TO 1000
  999 CALL LNCNT (1)
      WRITE(6,50) NA,NB
   50 FORMAT  (' DIMENSION ERROR IN ADD     NA=',2I6,5X,'NB=',2I6)
 1000 RETURN
      END
      SUBROUTINE ASYFIL(A,NA,G,NG,H,NH,Q,NQ,R,NR,F,NF,P,NP,IDENT,DISC,N
     1EWT,STABLE,FNULL,ALPHA,IOP,DUMMY)
C 
C   PURPOSE:
C      Solve either the continuous or discrete time-invariant asymptotic
C      optimal Kalman-Bucy filter problem.  Computation of both the
C      discrete and continuous versions of the optimal filter problems
C      is performed using duality theory.  No computations involve the
C      matrix W, therefore, no data for W are required.
C 
C   REFERENCES:
C      Kwakernaak, Huibert; and Sivan, Raphael: Linear Optimal Control
C        Systems.  John Wiley & Sons, Inc., c.1972.
C 
C   Subroutines employed by ASYFIL: ASYREG, EQUATE, LNCNT, PRNT, TRANP
C   Subroutines employing ASYFIL: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),G(*),H(*),Q(*),R(*),F(*),P(*),DUMMY(*)
      DIMENSION NA(2),NG(2),NH(2),NQ(2),NR(2),NF(2),NP(2),IOPT(5)
      DIMENSION NDUM1(2),IOP(*)
      LOGICAL  IDENT,DISC,NEWT,STABLE,FNULL
C 
      IF( IOP(1) .EQ. 0 ) GO TO 100
      CALL LNCNT(4)
      IF(DISC)  WRITE(6,15)
      IF( .NOT. DISC )  WRITE(6,25)
   15 FORMAT(//' PROGRAM TO SOLVE THE DISCRETE INFINITE-DURATION OPTIMAL
     1 FILTER PROBLEM'/)
   25 FORMAT(//' PROGRAM TO SOLVE THE CONTINUOUS INFINITE-DURATION OPTIM
     1AL FILTER PROBLEM'/)
      CALL PRNT(A,NA,' A  ',1)
      IF( .NOT.  IDENT )  GO TO 35
      CALL LNCNT(3)
      WRITE(6,30)
   30 FORMAT(/' G IS AN IDENTITY MATRIX'/)
      GO TO 40
   35 CONTINUE
      CALL PRNT(G,NG,' G  ',1)
   40 CONTINUE
      CALL PRNT(H,NH,' H  ',1)
      CALL LNCNT(3)
      WRITE(6,45)
   45 FORMAT(/' INTENSITY MATRIX FOR COVARIANCE OF MEASUREMENT NOISE'/)
      CALL PRNT(R,NR,' R  ',1)
C 
      IF( .NOT. IDENT ) GO TO 65
      CALL LNCNT(3)
      WRITE(6,55)
   55 FORMAT(/' INTENSITY MATRIX FOR COVARIANCE OF PROCESS NOISE'/)
C 
   65 CONTINUE
      CALL PRNT(Q,NQ,' Q  ',1)
C 
  100 CONTINUE
      IOPT(1)=IOP(2)
      IOPT(2)=IOP(3)
      IOPT(3)=IOP(4)
      IOPT(4)=IOP(5)
      IOPT(5)=0
      K = 0
C 
  200 CONTINUE
      CALL TRANP(A,NA,DUMMY,NA)
      CALL EQUATE(DUMMY,NA,A,NA)
      CALL TRANP(H,NH,DUMMY,NDUM1)
      CALL EQUATE(DUMMY,NDUM1,H,NH)
      IF( IDENT )  GO TO 250
      CALL TRANP(G,NG,DUMMY,NDUM1)
      CALL EQUATE(DUMMY,NDUM1,G,NG)
  250 CONTINUE
      IF ( K .EQ. 1 ) RETURN
C 
      K = K+1
      CALL ASYREG(A,NA,H,NH,G,NG,Q,NQ,R,NR,F,NF,P,NP,IDENT,DISC,NEWT,ST
     1ABLE,FNULL,ALPHA,IOPT,DUMMY)
C 
      N1=(NA(1)**2)+3*NA(1)+1
      CALL TRANP(F,NF,DUMMY(N1),NDUM1)
      CALL EQUATE(DUMMY(N1),NDUM1,F,NF)
C 
      IF( IOP(1) .EQ. 0 ) GO TO 200
C 
      IF(IDENT) GO TO 300
      CALL LNCNT(3)
      WRITE(6,55)
      CALL PRNT(Q,NQ,'GQGT',1)
C 
  300 CONTINUE
      CALL LNCNT(3)
      WRITE(6,325)
  325 FORMAT(/' FILTER GAIN'/)
      CALL PRNT(F,NF,' F  ',1)
      CALL LNCNT(3)
      WRITE(6,350)
  350 FORMAT(/' STEADY-STATE VARIANCE MATRIX OF RECONSTRUCTION ERROR'/)
      CALL PRNT(P,NP,' P  ',1)
      NDUM1(1)=NP(1)
      NDUM1(2)=1
      CALL LNCNT(3)
      WRITE(6,375)
  375 FORMAT(/' EIGENVALUES OF P '/)
      CALL PRNT(DUMMY,NDUM1,'EVLP',1)
      N1 = NP(1) + 1
      N = NA(1)**2
      N2 = N1 + N + 2*NA(1)
      CALL TRANP(DUMMY(N1),NA,DUMMY(N2),NA)
      CALL PRNT(DUMMY(N2),NA,'A-FH',1)
      N2 = N1 + N
      CALL LNCNT(3)
      WRITE(6,385)
  385 FORMAT(/' EIGENVALUES OF A-FH MATRIX'/)
      NDUM1(1) = NA(1)
      NDUM1(2) =  2
      CALL PRNT(DUMMY(N2),NDUM1,'    ', 3)
C 
      GO TO 200
C 
      END
      SUBROUTINE ASYREG(A,NA,B,NB,H,NH,Q,NQ,R,NR,F,NF,P,NP,IDENT,DISC,N
     1EWT,STABLE,FNULL,ALPHA,IOP,DUMMY)
C 
C   PURPOSE:
C      Solve either the continuous or discrete time-invariant asymp-
C      totic linear optimal output regulator problem with noise-free
C      measurements.  ASYREG does not evaluate the optimal values of
C      the performance criteria.  Therefore, no V data are input.  The
C      option of solving the appropriate steady-state Riccati equation
C      using either of the subroutines DISREG, CNTREG, or RICNWT is
C      provided.
C 
C   REFERENCES:
C      Kwakernaak, Huibert; and Sivan, Raphael: Linear Optimal Control
C        Systems.  John Wiley & Sons, Inc., c. 1972.
C 
C   Subroutines employed by ASYREG: ADD, CNTREG, CSTAB, DISREG, DSTAB,
C      EIGEN, EQUATE, JUXTC, LNCNT, MULT, PRNT, RICNWT, SCALE, SUBT,
C      TESTST, TRANP
C   Subroutines employing ASYREG: ASYFIL, EXMDFL, IMMDFL
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),H(*),Q(*),R(*),F(*),P(*),DUMMY(*)
      DIMENSION NA(2),NB(2),NH(2),NQ(2),NR(2),NF(2),NP(2),IOP(5),IOPT(3)
     1,NDUM1(2),NDUM2(2),NDUM3(2)
      LOGICAL IDENT,DISC,NEWT,STABLE,FNULL,SING
C 
      N = NA(1)**2
      N1= N+1
      IOPTT=0
C 
      IF ( .NOT. NEWT ) GO TO 600
      IF( STABLE )  GO TO 500
      IF ( FNULL ) GO TO 100
      CALL MULT(B,NB,F,NF,DUMMY,NA)
      CALL SUBT(A,NA,DUMMY,NA,DUMMY,NA)
      CALL TESTST(DUMMY,NA,ALPHA,DISC,STABLE,IOPTT,DUMMY(N1))
      GO TO 200
C 
  100 CONTINUE
      CALL TESTST(A,NA,ALPHA,DISC,STABLE,IOPTT,DUMMY)
C 
  200 CONTINUE
      IF( STABLE ) GO TO 500
      IF( DISC ) GO TO 230
      J = -NA(1)
      NAX = NA(1)
      DO 210 I =1,NAX
      J = J + NAX +1
      A(J) = A(J)-ALPHA
  210 CONTINUE
      SCLE = 3.
      IOPT(1)=IOP(1)
      IOPT(2) = 1
      IOPT(3)=1
      CALL CSTAB(A,NA,B,NB,F,NF,IOPT,SCLE,DUMMY)
      J = -NA(1)
      DO 220 I=1,NAX
      J = J + NAX + 1
      A(J) = A(J) + ALPHA
  220 CONTINUE
  225 CONTINUE
      CALL MULT(B,NB,F,NF,DUMMY,NA)
      CALL SUBT(A,NA,DUMMY,NA,DUMMY,NA)
      CALL TESTST(DUMMY,NA,ALPHA,DISC,STABLE,IOPTT,DUMMY(N1))
      GO TO 300
C 
  230 CONTINUE
      J = 2*NA(1) + 1
      IF( .NOT. FNULL )  J = J + N
      SING = .FALSE.
      IF( DUMMY(J) .EQ. 0.0 )  SING = .TRUE.
      IOPT(1) = IOP(1)
      IOPT(2) = 1
      DSCLE = 0.5
      ALPHAT = 1./ALPHA
      CALL SCALE(A,NA,A,NA,ALPHAT)
      CALL SCALE(B,NB,B,NB,ALPHAT)
      CALL DSTAB(A,NA,B,NB,F,NF,SING,IOPT,DSCLE,DUMMY)
      CALL SCALE(A,NA,A,NA,ALPHA)
      CALL SCALE(B,NB,B,NB,ALPHA)
      GO TO 225
C 
  300 CONTINUE
      IF( STABLE) GO TO 400
      CALL LNCNT(5)
      IF( DISC ) GO TO 330
      WRITE(6,310) ALPHA
  310 FORMAT(//' IN ASYREG, CSTAB HAS FAILED TO FIND A STABILIZING GAIN
     1 MATRIX (F) RELATIVE TO '/' ALPHA = ',D16.8/)
      RETURN
  330 CONTINUE
      WRITE(6,340) ALPHA
  340 FORMAT(//' IN ASYREG, DSTAB HAS FAILED TO FIND A STABILIZING GAIN
     1 MATRIX (F) RELATIVE TO '/' ALPHA = ',D16.8/)
      RETURN
C 
  400 CONTINUE
      FNULL = .FALSE.
C 
  500 CONTINUE
      CALL RICNWT(A,NA,B,NB,H,NH,Q,NQ,R,NR,F,NF,P,NP,IOP,IDENT,DISC,FNU
     1LL,DUMMY)
      GO TO 750
C 
  600 CONTINUE
      IF( DISC ) GO TO 700
      NW = 4*N + 1
      NLAM = NW + 4*N
      NDUM = NLAM + N
      IOP(3) = 1
      CALL CNTREG(A,NA,B,NB,H,NH,Q,NQ,R,NR,DUMMY,DUMMY(NW),DUMMY(NLAM),
     1S,F,NF,P,NP,T,IOP,IDENT,DUMMY(NDUM))
      GO TO 750
  700 CONTINUE
      CALL DISREG(A,NA,B,NB,H,NH,Q,NQ,R,NR,F,NF,P,NP,IOP,IDENT,DUMMY)
C 
  750 CONTINUE
C 
      IF( IOP(4) .EQ. 0 ) GO TO 1100
C 
      N2= N1 + N
      N3= N2 + N
C 
      IF( DISC ) GO TO 800
      CALL MULT(P,NP,B,NB,DUMMY,NB)
      CALL MULT(DUMMY,NB,F,NF,DUMMY(N1),NP)
      CALL TRANP(DUMMY(N1),NP,DUMMY,NP)
      CALL ADD(DUMMY,NP,DUMMY(N1),NP,DUMMY,NP)
      CALL SCALE(DUMMY,NP,DUMMY,NP,0.5D0)
      CALL SUBT(Q,NQ,DUMMY,NP,DUMMY,NP)
      CALL MULT(P,NP,A,NA,DUMMY(N1),NP)
      CALL ADD(DUMMY,NP,DUMMY(N1),NP,DUMMY,NP)
      CALL TRANP(DUMMY(N1),NP,DUMMY(N2),NP)
      CALL ADD(DUMMY,NP,DUMMY(N2),NP,DUMMY,NP)
      GO TO 900
C 
  800 CONTINUE
      CALL MULT(R,NR,F,NF,DUMMY,NF)
      CALL TRANP(F,NF,DUMMY(N1),NB)
      CALL MULT(DUMMY(N1),NB,DUMMY,NF,DUMMY(N2),NA)
      CALL ADD(DUMMY(N2),NA,Q,NQ,DUMMY,NA)
      CALL MULT(B,NB,F,NF,DUMMY(N1),NA)
      CALL SUBT(A,NA,DUMMY(N1),NA,DUMMY(N1),NA)
      CALL MULT(P,NP,DUMMY(N1),NA,DUMMY(N2),NA)
      CALL TRANP(DUMMY(N1),NA,DUMMY(N3),NA)
      CALL MULT(DUMMY(N3),NA,DUMMY(N2),NA,DUMMY(N1),NA)
      CALL ADD(DUMMY,NA,DUMMY(N1),NA,DUMMY,NA)
      CALL SUBT(P,NP,DUMMY,NA,DUMMY,NA)
C 
  900 CONTINUE
      CALL LNCNT(4)
      WRITE(6,1000)
 1000 FORMAT(//' RESIDUAL ERROR IN RICCATI EQUATION '/)
      CALL PRNT(DUMMY,NP,'EROR',1)
C 
 1100 CONTINUE
      N2= N1+NA(1)
      N3= N2+NA(1)
      ISV = 0
      CALL EQUATE(P,NP,DUMMY,NP)
      CALL EIGEN(NA(1),NA(1),DUMMY,DUMMY(N1),DUMMY(N2),ISV,ISV,V,DUMMY(N
     13),IERR)
      NEVL = NA(1)
      IF( IERR .EQ. 0) GO TO 1300
      NEVL=NA(1)-IERR
      CALL LNCNT(4)
      WRITE(6,1200) IERR
 1200 FORMAT(//' IN ASYREG, THE ',I5 ,' EIGENVALUE OF P  HAS NOT BEEN C
     1OMPUTED AFTER 30 ITERATIONS '/)
C 
 1300 CONTINUE
      NDUM1(1) = NEVL
      NDUM1(2) = 1
      CALL EQUATE(DUMMY(N1),NDUM1,DUMMY,NDUM1)
      N1 = NDUM1(1) +1
      CALL MULT(B,NB,F,NF,DUMMY(N1),NA)
      CALL SUBT(A,NA,DUMMY(N1),NA,DUMMY(N1),NA)
      N2 = N1+N
      CALL EQUATE(DUMMY(N1),NA,DUMMY(N2),NA)
      N3=N2+N
      N4=N3+NA(1)
      N5=N4+NA(1)
      CALL EIGEN(NA(1),NA(1),DUMMY(N2),DUMMY(N3),DUMMY(N4),ISV,ISV,V,DUM
     1MY(N5),IERR)
      NEVL = NA(1)
      IF( IERR .EQ. 0 ) GO TO 1500
      NEVL=NA(1)-IERR
      CALL LNCNT(4)
      WRITE(6,1400) IERR
 1400 FORMAT(//' IN ASYREG, THE ',I5,' EIGENVALUE OF A-BF HAS NOT BEEN
     1COMPUTED AFTER 30 ITERATIONS'/)
C 
 1500 CONTINUE
      NDUM2(1) = NEVL
      NDUM2(2) = 1
      CALL JUXTC(DUMMY(N3),NDUM2,DUMMY(N4),NDUM2,DUMMY(N2),NDUM3)
C 
      IF ( IOP(5) .EQ. 0 ) RETURN
C 
      CALL LNCNT(4)
      WRITE(6,1600)
 1600 FORMAT(//' EIGENVALUES OF P '/)
      CALL PRNT(DUMMY,NDUM1,'EVLP',1)
      CALL LNCNT(4)
      WRITE(6,1700)
 1700 FORMAT(//' CLOSED-LOOP RESPONSE MATRIX A-BF '/)
      CALL PRNT(DUMMY(N1),NA,'A-BF',1)
      CALL LNCNT(3)
      WRITE(6,1800)
 1800 FORMAT(//' EIGENVALUES OF A-BF')
      CALL PRNT(DUMMY(N2),NDUM3,'    ',3)
C 
      RETURN
      END
      SUBROUTINE ATXPXA(A,U,C,N,NA,NU,NC,EPS,FAIL)
C 
C   PURPOSE:
C      Solve the real matrix equation A'X + XA = C, where A and C are
C      constant matrices of dimension n x n with C=C'.  The matrix A
C      is transformed into upper Schur form and the transformed system
C      is solved by back substitution.  The option is provided to input
C      the Schur form directly and bypass the Schur decomposition.
C   REFERENCES:
C      Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
C        the Matrix Equation AX + XB = C.  Commun. ACM, vol. 15, no. 9,
C        Sept. 1972, pp. 820-826.
C 
C   Subroutines employed by ATXPXA: BCKMLT, HSHLDR, SCHUR, SYMSLV
C   Subroutines employing ATXPXA: BARSTW
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NA,1),U(NU,1),C(NC,1),EPS
      INTEGER N,NA,NU,NC,FAIL,N1,NM1,I,J,K
      N1 = N+1
      NM1 = N-1
C 
C IF REQUIRED, REDUCE A TO LOWER REAL SCHUR FORM.
C 
      IF(EPS .LT. 0.) GO TO 15
      CALL HSHLDR(A,N,NA)
      CALL BCKMLT(A,U,N,NA,NU)
      DO 10 I=1,NM1
        A(I+1,I) = A(I,N1)
   10 CONTINUE
      CALL SCHUR(A,U,N,NA,NU,EPS,FAIL)
      IF(FAIL .NE. 0) RETURN
C 
C TRANSFORM C.
C 
   15 DO 20 I=1,N
          C(I,I)=C(I,I)/2.
   20 CONTINUE
      DO 40 I=1,N
        DO 30 J=1,N
          A(N1,J) = 0.
          DO 30 K=I,N
            A(N1,J) = A(N1,J) + C(I,K)*U(K,J)
   30   CONTINUE
          DO 40 J=1,N
          C(I,J) = A(N1,J)
   40 CONTINUE
      DO 60 J=1,N
        DO 50 I=1,N
          A(I,N1) = 0.
          DO 50 K=1,N
            A(I,N1) = A(I,N1) + U(K,I)*C(K,J)
   50   CONTINUE
        DO 60 I=1,N
          C(I,J) = A(I,N1)
   60 CONTINUE
      DO 70 I=1,N
        DO 70 J=I,N
          C(I,J) = C(I,J) + C(J,I)
          C(J,I) = C(I,J)
   70 CONTINUE
C 
C SOLVE THE TRANSFORMED SYSTEM.
C 
      CALL SYMSLV(A,C,N,NA,NC)
C 
C TRANSFORM C BACK TO THE SOLUTION.
C 
      DO 80 I=1,N
        C(I,I) = C(I,I)/2.
   80 CONTINUE
      DO 100 I=1,N
        DO 90 J=1,N
          A(N1,J) = 0.
          DO 90 K=I,N
            A(N1,J) = A(N1,J) + C(I,K)*U(J,K)
   90   CONTINUE
        DO 100 J=1,N
          C(I,J) = A(N1,J)
  100 CONTINUE
      DO 120 J=1,N
        DO 110 I=1,N
          A(I,N1) = 0.
          DO 110 K=1,N
            A(I,N1) = A(I,N1) + U(I,K)*C(K,J)
  110   CONTINUE
        DO 120 I=1,N
          C(I,J) = A(I,N1)
  120 CONTINUE
      DO 130 I=1,N
        DO 130 J=I,N
          C(I,J) = C(I,J) + C(J,I)
          C(J,I) = C(I,J)
  130 CONTINUE
      RETURN
      END
      SUBROUTINE AXPXB(A,U,M,NA,NU,B,V,N,NB,NV,C,NC,EPSA,
     1EPSB,FAIL)
C 
C   PURPOSE:
C      Solve the real matrix equation AX + XB = C, where A, B, and C
C      are constant matrices of order m x n, n x n, and m x n.  The ma-
C      trices are transformed into real lower and upper Schur form, and
C      the transformed system is solved by back substitution.  The op-
C      tion is provided to input the Schur forms directly and bypass
C      the Schur decomposition.
C 
C   REFERENCES:
C      Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
C        the Matrix Equation AX + XB = C.  Commun. ACM, vol. 15, no. 9,
C        Sept. 1972, p. 820-826.
C 
C   Subroutines employed by AXPXB: BCKMLT, HSHLDR, SCHUR, SHRSLV
C   Subroutines employing AXPXB: BARSTW
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NA,1),U(NU,1),B(NB,1),V(NV,1),C(NC,1),EPSA,EPSB,TEMP
      INTEGER M,NA,NU,N,NB,NV,NC,FAIL,M1,MM1,N1,NM1,I,J,K
      M1 = M+1
      MM1 = M-1
      N1 = N+1
      NM1 = N-1
C 
C IF REQUIRED, REDUCE A TO UPPER REAL SCHUR FORM.
C 
      IF(EPSA .LT. 0.) GO TO 35
      DO 10 I=1,M
        DO 10 J=I,M
          TEMP = A(I,J)
          A(I,J) = A(J,I)
          A(J,I) = TEMP
   10 CONTINUE
      CALL HSHLDR(A,M,NA)
      CALL BCKMLT(A,U,M,NA,NU)
      IF(MM1 .EQ. 0) GO TO 25
      DO 20 I=1,MM1
        A(I+1,I) = A(I,M1)
   20 CONTINUE
      CALL SCHUR(A,U,M,NA,NU,EPSA,FAIL)
      IF(FAIL .NE. 0) RETURN
   25 DO 30 I=1,M
        DO 30 J=I,M
          TEMP = A(I,J)
          A(I,J) = A(J,I)
          A(J,I) = TEMP
   30 CONTINUE
C 
C IF REQUIRED, REDUCE B TO UPPER REAL SCHUR FORM.
C 
   35 IF(EPSB .LT. 0.) GO TO 45
      CALL HSHLDR(B,N,NB)
      CALL BCKMLT(B,V,N,NB,NV)
      IF(NM1 .EQ. 0) GO TO 45
      DO 40 I=1,NM1
        B(I+1,I) = B(I,N1)
   40 CONTINUE
      CALL SCHUR(B,V,N,NB,NV,EPSB,FAIL)
      FAIL = -FAIL
      IF(FAIL .NE. 0) RETURN
C 
C TRANSFORM C.
C 
   45 DO 60 J=1,N
        DO 50 I=1,M
          A(I,M1) = 0.
          DO 50 K=1,M
            A(I,M1) = A(I,M1) + U(K,I)*C(K,J)
   50 CONTINUE
      DO 60 I=1,M
        C(I,J) = A(I,M1)
   60 CONTINUE
      DO 80 I=1,M
        DO 70 J=1,N
          B(N1,J) = 0.
          DO 70 K=1,N
            B(N1,J) = B(N1,J) + C(I,K)*V(K,J)
   70 CONTINUE
      DO 80 J=1,N
        C(I,J) = B(N1,J)
   80 CONTINUE
C 
C SOLVE THE TRANSFORMED SYSTEM.
C 
      CALL SHRSLV(A,B,C,M,N,NA,NB,NC)
C 
C TRANSFORM C BACK TO THE SOLUTION.
C 
      DO 100 J=1,N
        DO 90 I=1,M
          A(I,M1) = 0.
          DO 90 K=1,M
            A(I,M1) = A(I,M1) + U(I,K)*C(K,J)
   90 CONTINUE
      DO 100 I=1,M
        C(I,J) = A(I,M1)
  100 CONTINUE
      DO 120 I=1,M
        DO 110 J=1,N
          B(N1,J) = 0.
          DO 110 K=1,N
            B(N1,J) = B(N1,J) + C(I,K)*V(J,K)
  110    CONTINUE
         DO 120 J=1,N
           C(I,J) = B(N1,J)
  120  CONTINUE
       RETURN
       END
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
C 
C   PURPOSE:
C      Balance a real square matrix for the calculation of eigenvalues
C      and eigenvectors and isolate the eigenvalues whenever possible.
C      The computation follows that of Parlett and Reinsch.
C 
C   REFERENCES:
C      Parlett, B.N.; and Reinsch, C.: Balancing a Matrix for Calcula-
C        tion of Eigenvalues and Eigenvectors. Numer. Math., Bd. 13,
C        Heft 4, 1969, pp.293-304.
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Computa-
C        tion.  Volume II - Linear Algebra.  Springer-Verlag, 1971.
C 
C   Subroutines employed by BALANC: None
C   Subroutines employing BALANC: EIGEN
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL*8 A(NM,N),SCALE(N)
      REAL*8 C,F,G,R,S,B2,RADIX
C     REAL*8 DABS
      LOGICAL NOCONV
C 
C 
C     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.
C 
C 
      RADIX = 2.
C 
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     ********** IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE **********
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C 
      DO 30 I = 1, L
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
C 
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
C 
   50 GO TO (80,130), IEXC
C     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN **********
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
C     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C 
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (A(J,I) .NE. 0.0) GO TO 120
  110    CONTINUE
C 
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C 
      GO TO 140
C     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT **********
  130 K = K + 1
C 
  140 DO 170 J = K, L
C 
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (A(I,J) .NE. 0.0) GO TO 170
  150    CONTINUE
C 
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********
      DO 180 I = K, L
  180 SCALE(I) = 1.0
C     ********** ITERATIVE LOOP FOR NORM REDUCTION **********
  190 NOCONV = .FALSE.
C 
      DO 270 I = K, L
         C = 0.0
         R = 0.0
C 
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + DABS(A(J,I))
            R = R + DABS(A(I,J))
  200    CONTINUE
C     ********** GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW **********
         IF (C .EQ. 0.0 .OR. R .EQ. 0.0) GO TO 270
         G = R / RADIX
         F = 1.0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     ********** NOW BALANCE **********
  240    IF ((C + R) / F .GE. 0.95 * S) GO TO 270
         G = 1.0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C 
         DO 250 J = K, N
  250    A(I,J) = A(I,J) * G
C 
         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
C 
  270 CONTINUE
C 
      IF (NOCONV) GO TO 190
C 
  280 LOW = K
      IGH = L
      RETURN
C     ********** LAST CARD OF BALANC **********
      END
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
C 
C   PURPOSE:
C      Form the eigenvectors of a real square matrix by back transform-
C      ing those of the corresponding balanced matrix determined by
C      subroutine BALANC.
C 
C   REFERENCES:
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Computa-
C        tion.  Volume II - Linear Algebra.  Springer-Verlag, 1971.
C 
C   Subroutines employed by BALBAK: None
C   Subroutines employing BALBAK: EIGEN
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL*8 SCALE(N),Z(NM,M)
      REAL*8 S
C 
C 
C 
      IF (M .EQ. 0) GO TO 200
      IF (IGH .EQ. LOW) GO TO 120
C 
      DO 110 I = LOW, IGH
         S = SCALE(I)
C     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0/SCALE(I). **********
         DO 100 J = 1, M
  100    Z(I,J) = Z(I,J) * S
C 
  110 CONTINUE
C     ********- FOR I=LOW-1 STEP -1 UNTIL 1,
C               IGH+1 STEP 1 UNTIL N DO -- **********
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C 
         DO 130 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
  130    CONTINUE
C 
  140 CONTINUE
C 
  200 RETURN
C     ********** LAST CARD OF BALBAK **********
      END
      SUBROUTINE BARSTW(A,NA,B,NB,C,NC,IOP,SYM,EPSA,EPSB,DUMMY)
C 
C   PURPOSE:
C      Solve the matrix equation AX + XB = C, where A and B are real
C      constant matrices of dimension n x n and m x m.  The matrix C is
C      real constant and of dimension n x m.  Assumed that
C              A(i) + B(i) # 0  (i=1,2,...,n;  j=1,2,...,m)
C      where A(i) and B(i) are eigenvalues of A and B.
C 
C   REFERENCES:
C      Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
C        the Matrix Equation AX + XB = C. Commun. ACM, vol. 15, no. 9,
C        Sept. 1972, pp. 820-826.
C 
C   Subroutines employed by BARSTW: ATXPXA, AXPXB, EQUATE, JUXTR, LNCNT,
C      NULL, PRNT, TRANP
C   Subroutines employing BARSTW: CSTAB, DSTAB, EXMDFL, RICNWT, VARANC
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NC(2),NDUM1(2),NDUM2(2),NDUM3(2),NDUM4(2)
      LOGICAL  SYM
C 
      IF ( IOP .EQ. 0 )  GO TO 250
      IF(SYM) GO TO 100
      CALL LNCNT(3)
      WRITE(6,50)
   50 FORMAT(//' LINEAR EQUATION SOLVER     AX + XB = C ')
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      GO TO 200
  100 CONTINUE
      CALL LNCNT(3)
      WRITE(6,150)
  150 FORMAT(//' LINEAR EQUATION SOLVER  ( B TRANSPOSE )X + XB = C')
      CALL TRANP(A,NA,DUMMY,NDUM1)
      CALL PRNT(DUMMY,NDUM1,' B  ',1)
  200 CONTINUE
      CALL PRNT(C,NC,' C  ',1)
C 
  250 CONTINUE
      CALL EQUATE(A,NA,DUMMY,NDUM1)
      N1=(NA(1)**2)+1
      N2=N1+NA(1)-1
      DO 300I=N1,N2
      DUMMY(I)=0.0
  300 CONTINUE
C 
      NDUM1(2)=NDUM1(2)+1
      NDUM2(1)=1
      NDUM2(2)=NDUM1(2)
      N1=NDUM1(1)*NDUM1(2)+1
      CALL NULL(DUMMY(N1),NDUM2)
      LU=(NA(1)+1)**2 + 1
      CALL JUXTR(DUMMY,NDUM1,DUMMY(N1),NDUM2,DUMMY(LU),NDUM3)
      CALL EQUATE(DUMMY(LU),NDUM3,DUMMY,NDUM1)
      N=NA(1)+1
C 
      IF(SYM ) GO TO 500
C 
      CALL EQUATE(B,NB,DUMMY(LU),NDUM2)
      M1=LU+NB(1)**2
      M2=M1+NB(1)-1
      DO400I=M1,M2
      DUMMY(I)=0.0
  400 CONTINUE
C 
      NDUM2(2)=NDUM2(2)+1
      NDUM3(1)=1
      NDUM3(2)=NDUM2(2)
      M1=NDUM2(1)*NDUM2(2)+LU
      CALL NULL(DUMMY(M1),NDUM3)
      M2=LU+(NB(1)+1)**2
      CALL JUXTR(DUMMY(LU),NDUM2,DUMMY(M1),NDUM3,DUMMY(M2),NDUM4)
      CALL EQUATE(DUMMY(M2),NDUM4,DUMMY(LU),NDUM2)
      M=NB(1)+ 1
      LNB = LU
      LU = LU + (NB(1)+1)**2
      LV = LU +  NA(1)**2
      CALL AXPXB(DUMMY,DUMMY(LU),NA(1),N,NA(1),DUMMY(LNB),DUMMY(LV),NB(1
     1),M,NB(1),C,NC(1),EPSA,EPSB,NFAIL)
      GO TO 600
C 
  500 CONTINUE
      CALL TRANP(DUMMY,NDUM1,DUMMY(LU),NDUM2)
      CALL EQUATE(DUMMY(LU),NDUM2,DUMMY,NDUM1)
      CALL ATXPXA(DUMMY,DUMMY(LU),C,NA(1),N,NA(1),NC(1),EPSA,NFAIL)
C 
  600 CONTINUE
      IF(NFAIL .EQ. 0 ) GO TO 700
      CALL LNCNT(3)
      WRITE(6,650)
  650 FORMAT(//' IN BARSTW, EITHER THE SUBROUTINE AXPXB  OR  ATXPXA  WAS
     1 UNABLE TO REDUCE A OR B TO SCHUR FORM ')
      RETURN
C 
  700 CONTINUE
C 
      IF( IOP .NE. 0 )  CALL PRNT(C,NC,'  X ',1)
      RETURN
      END
      SUBROUTINE BCKMLT(A,U,N,NA,NU)
C 
C   PURPOSE:
C      Compute the orthogonal matrix that reduces the output matrix A
C      from subroutine HSHLDR, to upper Hessenberg form.
C 
C   REFERENCES:
C      Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
C        the Matrix Equation AX + XB = C.  Commun. ACM, vol. 15, no. 9,
C        Sept. 1972, pp. 820-826.
C 
C   Subroutines employed by BCKMLT: None
C   Subroutines employing BCKMLT: ATXPXA, AXPXB
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8
     1A(NA,1),U(NU,1),SUM,P
      INTEGER
     1N,NA,N1,NM1,NM2,LL,L,L1,I,J
      N1 = N+1
      NM1 = N-1
      NM2 = N-2
      U(N,N) = 1.
      IF(NM1 .EQ. 0) RETURN
      U(NM1,N) = 0.
      U(N,NM1) = 0.
      U(NM1,NM1) = 1.
      IF(NM2 .EQ. 0) RETURN
      DO 40 LL=1,NM2
        L = NM2-LL+1
        L1 = L+1
        IF(A(N1,L) .EQ. 0.) GO TO 25
        DO 20 J=L1,N
          SUM = 0.
          DO 10 I=L1,N
            SUM = SUM + A(I,L)*U(I,J)
   10     CONTINUE
          P = SUM/A(N1,L)
          DO 20 I=L1,N
            U(I,J) = U(I,J) - A(I,L)*P
   20   CONTINUE
   25   DO 30 I=L1,N
          U(I,L) = 0.
          U(L,I) = 0.
   30   CONTINUE
        U(L,L) = 1.
   40 CONTINUE
      RETURN
      END
      SUBROUTINE BILIN(A,NA,B,NB,C,NC,IOP,BETA,SYM,DUMMY)
C 
C   PURPOSE:
C      Solve the matrix equation, AX + XB = C, where A and B are real
C      constant matrices of dimension n x n and m x m. The matrix C is
C      real constant and of dimension n x m. Assumed that all eigenva-
C      lues of A and B have strictly negative real parts.
C 
C   REFERENCES:
C      Smith, R.A.: Matrix Equation XA + BX = C. SIAM J. Appl. Math.,
C        vol. 16, no. 2, Mar. 1968, pp. 198-201.
C 
C   Subroutines employed by BILIN: DGECO, DGESL, EIGEN, EQUATE, LNCNT,
C      MULT, NORMS, PRNT, SCALE, SUM, TRANP, UNITY
C   Subroutines employing BILIN: CSTAB, RICNWT, VARANC
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NC(2),NDUM(2)
      DIMENSION IOP(2)
      LOGICAL SYM
C 
      IF( IOP(1) .EQ. 0 )  GO TO 300
      IF(SYM) GO TO 100
      CALL LNCNT(3)
      WRITE(6,50)
   50 FORMAT(//' LINEAR EQUATION SOLVER  AX + XB = C ')
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      GO TO 200
  100 CONTINUE
      CALL LNCNT(3)
      WRITE(6,150)
  150 FORMAT(//' LINEAR EQUATION SOLVER  ( B TRANSPOSE )X + XB = C ')
      CALL TRANP(A,NA,DUMMY,NDUM)
      CALL PRNT(DUMMY,NDUM,' B  ',1)
  200 CONTINUE
      CALL PRNT(C,NC,' C  ',1)
  300 CONTINUE
C 
      IOPTT = 0
      N=NA(1)**2
      M=NB(1)**2
C 
      IF( IOP(2) .EQ. 0 )  GO TO 500
C 
      N1 = N + 1
      CALL EQUATE(A,NA,DUMMY,NA)
      N2 = N1 + NA(1)
      N3 = N2 + NA(1)
      ISV = 0
      ILV = 0
      NEVL = NA(1)
      CALL EIGEN(NA(1),NA(1),DUMMY,DUMMY(N1),DUMMY(N2),ISV,ILV,V,DUMMY(N
     13),IERR)
      IF (IERR .EQ. 0) GO TO 350
      CALL LNCNT(3)
      WRITE(6,325) IERR
  325 FORMAT(//' IN BILIN, THE ',I4,' EIGENVALUE OF A HAS NOT BEEN  DETE
     1RMINED AFTER 30 ITERATIONS')
      IERR=1
      CALL NORMS(NEVL,NEVL,NEVL,A,IERR,BETA)
      BETA=2.*BETA
      GO TO 385
  350 CONTINUE
      J= N1 + NEVL -1
      K = N2 + NEVL -1
      CO = DSQRT(DUMMY(N1)**2 + DUMMY(N2)**2)
      CN = DSQRT(DUMMY(J)**2 + DUMMY(K)**2)
      CD = DUMMY(J)-DUMMY(N1)
      IF(CD .EQ. 0.0)  GO TO 365
      BETA = (DUMMY(N1)*CN-DUMMY(J)*CO)/CD
      IF(BETA .LE. 0.0)  GO TO 365
      BETA = DSQRT(BETA)
      GO TO 385
C 
  365 CONTINUE
C 
      BETA = 0.0
      DO 375 I = 1,NEVL
      J = N1 + I -1
      K = N2 + I -1
      IF(DUMMY(J) .GE. 0.0)  GO TO 375
      BETA = BETA + DSQRT(DUMMY(J)**2 + DUMMY(K)**2)
  375 CONTINUE
      BETA = BETA/NEVL
C 
  385 CONTINUE
C 
      IF( SYM ) GO TO 500
      CALL EQUATE(B,NB,DUMMY,NB)
      N1=M+1
      N2 = N1 +NB(1)
      N3 = N2 +NB(1)
      NEVL = NB(1)
      CALL EIGEN(NB(1),NB(1),DUMMY,DUMMY(N1),DUMMY(N2),ISV,ILV,V,DUMMY(N
     13),IERR)
      IF(IERR .EQ. 0) GO TO 450
      CALL LNCNT(3)
      WRITE(6,400) IERR
  400 FORMAT(//' IN BILIN, THE ',I4,' EIGENVALUE OF B HAS NOT BEEN FOUND
     1 AFTER 30 ITERATIONS')
      IERR=1
      CALL NORMS(NEVL,NEVL,NEVL,B,IERR,BETA1)
      BETA1=2.*BETA1
      GO TO 485
  450 CONTINUE
      J = N1 + NEVL -1
      K = N2 + NEVL -1
      CO = DSQRT(DUMMY(N1)**2 + DUMMY(N2)**2)
      CN = DSQRT(DUMMY(J)**2 + DUMMY(K)**2)
      CD = DUMMY(J)-DUMMY(N1)
      IF(CD .EQ. 0.0)  GO TO 465
      BETA1 = (DUMMY(N1)*CN - DUMMY(J)*CO)/CD
      IF(BETA1 .LE. 0.0)  GO TO 465
      BETA1 = DSQRT(BETA1)
      GO TO 485
C 
  465 CONTINUE
C 
      BETA1 = 0.0
      DO 475 I= 1,NEVL
      J = N1 + I -1
      K = N2 + I -1
      IF(DUMMY(J) .GE. 0.0)  GO TO 475
      BETA1 = BETA1 + DSQRT(DUMMY(J)**2 + DUMMY(K)**2)
  475 CONTINUE
      BETA1 = BETA1/NEVL
C 
  485 CONTINUE
      BETA = (BETA + BETA1)/2.
C 
  500 CONTINUE
C 
C 
      IF( IOP(1) .EQ. 0 )  GO TO 520
      CALL LNCNT(4)
      WRITE(6,515) BETA
  515 FORMAT(//' BETA = ',D16.8/)
  520 CONTINUE
C 
      N1 = N+1
      CALL EQUATE(A,NA,DUMMY,NA)
      CALL EQUATE(A,NA,DUMMY(N1),NA)
      CALL SCALE(DUMMY,NA,DUMMY,NA,-1.0D0)
      L = -NA(1)
      NAX = NA(1)
      DO 525 I=1,NAX
      L = L + NAX +1
      M1 = L + N
      DUMMY(L) = BETA - A(L)
      DUMMY(M1)= BETA + A(L)
  525 CONTINUE
      N2 = N1 + N
      CALL EQUATE(C,NC,DUMMY(N2),NDUM)
      NDUM(2)= NDUM(2) + NA(1)
      N3 = N2 + NC(1)*NC(2)
      GAM = -2.*BETA
C 
      IF( .NOT. SYM ) GO  TO 600
C 
      CALL UNITY(DUMMY(N3),NA)
      N4 = N3 + N
      NDUM(2) = NDUM(2) + NA(1)
      N5 = N4 + NA(1)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY,NA(1),NA(1),DUMMY(N4),RCOND,DUMMY(N5))
      IF ((1.0 + RCOND) .EQ. 1.0) WRITE(6,625) RCOND
      NT = N1
      DO 530 M2 = 1,NDUM(2)
         CALL DGESL(DUMMY,NA(1),NA(1),DUMMY(N4),DUMMY(NT),0)
         NT = NT + NA(1)
  530 CONTINUE
      CALL EQUATE(DUMMY(N1),NA,DUMMY,NA)
      CALL EQUATE(DUMMY(N2),NC,C,NC)
      CALL TRANP(DUMMY,NA,DUMMY(N1),NA)
      CALL TRANP(DUMMY(N3),NA,DUMMY(N2),NA)
      CALL MULT(C,NC,DUMMY(N2),NA,DUMMY(N3),NA)
      CALL SCALE(DUMMY(N3),NC,C,NC,GAM)
C 
C 
      CALL SUM(DUMMY,NA,C,NC,DUMMY(N1),NA,IOPTT,SYM,DUMMY(N2))
      GO TO 700
  600 CONTINUE
      N4 = N3 +NA(1)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY,NA(1),NA(1),DUMMY(N3),RCOND,DUMMY(N4))
      IF ((1.0 + RCOND) .EQ. 1.0) WRITE(6,625) RCOND
      NT = N1
      DO 620 M2 = 1,NDUM(2)
         CALL DGESL(DUMMY,NA(1),NA(1),DUMMY(N3),DUMMY(NT),0)
         NT = NT + NA(1)
  620 CONTINUE
  625 FORMAT(//' IN BILIN, THE MATRIX  (BETA)I - A IS SINGULAR, INCREASE
     1 BETA, RCOND = ',D16.8)
      CALL EQUATE(DUMMY(N1),NA,DUMMY,NA)
      CALL EQUATE(DUMMY(N2),NC,C,NC)
      N2 = M + N1
      CALL EQUATE(B,NB,DUMMY(N1),NB)
      CALL EQUATE(B,NB,DUMMY(N2),NB)
      CALL SCALE(DUMMY(N1),NB,DUMMY(N1),NB,-1.0D0)
      L=-NB(1)
      NAX=NB(1)
      DO650I =1,NAX
      L=L + NAX +1
      L1 = L + N
      M1 = L + N2-1
      DUMMY(L1)= BETA- B(L)
      DUMMY(M1)= BETA + B(L)
  650 CONTINUE
C 
      N3 = N2 + M
      CALL TRANP(DUMMY(N1),NB,DUMMY(N3),NB)
      CALL EQUATE(DUMMY(N3),NB,DUMMY(N1),NB)
      CALL TRANP(DUMMY(N2),NB,DUMMY(N3),NB)
      CALL EQUATE(DUMMY(N3),NB,DUMMY(N2),NB)
      CALL TRANP(C,NC,DUMMY(N3),NDUM)
      NSDUM = NDUM(2)
      NDUM(2)= NDUM(2) + NB(2)
      N4=N3+NC(1)*NC(2)
      N5=N4+NB(1)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY(N1),NB(1),NB(1),DUMMY(N4),RCOND,DUMMY(N5))
      IF ((1.0 + RCOND) .EQ. 1.0) WRITE(6,680) RCOND
      NT = N2
      DO 660 M2 = 1,NDUM(2)
         CALL DGESL(DUMMY(N1),NB(1),NB(1),DUMMY(N4),DUMMY(NT),0)
         NT = NT + NB(1)
  660 CONTINUE
  680 FORMAT(//' IN BILIN, THE MATRIX (BETA)I - B IS SINGULAR, INCREASE
     1 BETA, RCOND = ',D16.8)
      CALL TRANP(DUMMY(N2),NB,DUMMY(N1),NB)
      NDUM(2)= NSDUM
      CALL TRANP(DUMMY(N3),NDUM,C,NC)
      CALL SCALE(C,NC,C,NC,GAM)
      N2 = N + M + 1
      CALL SUM(DUMMY,NA,C,NC,DUMMY(N1),NB,IOPTT,SYM,DUMMY(N2))
C 
  700 CONTINUE
      IF( IOP (1).EQ. 0 ) RETURN
      CALL PRNT(C,NC,' X  ',1)
      RETURN
      END
      SUBROUTINE CDIV (X,Y,R,S,T1,T2)
C 
C   PURPOSE:
C      Perform a double precision complex divide and returns the result
C      in T1 and T2.
C 
C   Subroutines employed by CDIV: None
C   Subroutines employing CDIV: HQR2, INVIT
C 
C       T1 + J*T2= (X + J*Y)/(R + J*S)
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 R,S,T1,T2,TEMP,X,Y
C 
      TEMP= R*R + S*S
      T1= (X*R + Y*S)/TEMP
      T2= (Y*R - X*S)/TEMP
      RETURN
      END
      SUBROUTINE CNTREG(A,NA,B,NB,H,NH,Q,NQ,R,NR,Z,W,LAMBDA,S,F,NF,P,NP
     1,T,IOP,IDENT,DUMMY)
C 
C   PURPOSE:
C      Solve the time-invariant continuous-time linear optimal output
C      regulator problem with noise-free measurements.
C 
C   REFERENCES:
C      Kwakernaak, Huibert; and Sivan Raphael: Linear Optimal Control
C        Systems. John Wiley & Sons, Inc., c. 1972.
C      Vaughan, David R.: A Negative Exponential Solution for the Matrix
C        Riccati Equation. IEEE Trans. Autom. Control, vol. AC-14, no.
C        1, Feb. 1969, pp. 72-75.
C 
C   Subroutines employed by CNTREG: ADD, DGECO, DGESL, DPOCO, DPOSL,
C      EIGEN, EQUATE, EXPSER, JUXTR, LNCNT, MAXEL, MULT, NULL, PRNT,
C      SCALE, SUBT, TRANP
C   Subroutines employing CNTREG: ASYREG
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),H(*),Q(*),R(*),Z(*),W(*),LAMBDA(*),S(*)
      DIMENSION F(*),P(*),T(*),DUMMY(*)
      DIMENSION NA(2),NB(2),NH(2),NQ(2),NR(2),NF(2),NP(2),IOP(3)
      DIMENSION NDUM1(2),NDUM2(2)
      LOGICAL IDENT
      REAL*8 LAMBDA
CCC      COMMON/CONV/SUMCV,MAXSUM,RICTCV,SERCV
      COMMON/CONV/SUMCV,RICTCV,SERCV,MAXSUM
C 
      IF( IOP(1). EQ. 0 ) GO TO 65
      CALL LNCNT(5)
      IF( IOP(3) .EQ. 0 ) WRITE(6,25)
   25 FORMAT(//' PROGRAM TO SOLVE THE TIME-INVARIANT FINITE-DURATION CON
     1TINUOUS OPTIMAL'/' REGULATOR PROBLEM WITH NOISE-FREE MEASUREMENTS'
     2)
      IF( IOP(3) .NE. 0 ) WRITE(6,30)
   30 FORMAT(//' PROGRAM TO SOLVE THE TIME-INVARIANT INFINITE-DURATION C
     1ONTINUOUS OPTIMAL'/' REGULATOR PROBLEM  WITH NOISE-FREE MEASUREMEN
     2TS')
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      CALL PRNT(Q,NQ,' Q  ',1)
      IF( .NOT. IDENT ) GO TO 45
      CALL LNCNT(3)
      WRITE(6,35)
   35 FORMAT(/' H IS AN IDENTITY MATRIX'/)
      GO TO 55
C 
   45 CONTINUE
      CALL PRNT(H,NH,' H  ',1)
      CALL MULT(Q,NQ,H,NH,DUMMY,NH)
      N1= NH(1)*NH(2)+1
      CALL TRANP(H,NH,DUMMY(N1),NDUM1)
      CALL MULT(DUMMY(N1),NDUM1,DUMMY,NH,Q,NQ)
      CALL LNCNT(3)
      WRITE(6,50)
   50 FORMAT(//' MATRIX (H TRANSPOSE)QH')
      CALL PRNT(Q,NQ,'    ',3)
   55 CONTINUE
      CALL PRNT(R,NR,' R  ',1)
C 
      IF( IOP(3) .NE. 0 ) GO TO 65
      CALL LNCNT(4)
      WRITE(6,60)
   60 FORMAT(//' WEIGHTING ON TERMINAL VALUE OF STATE VECTOR'/)
      CALL PRNT(P,NP,' P  ',1)
C 
   65 CONTINUE
      CALL EQUATE(R,NR,DUMMY,NR)
      N = NA(1)**2
      N1 = NB(1)*NB(2)+1
      CALL TRANP(B,NB,DUMMY(N1),NDUM1)
      N2 = N1 + N
      L = NR(1)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DPOCO(DUMMY,L,L,RCOND,DUMMY(N2),IERR)
      IF( IERR .EQ. 0 ) GO TO 100
      CALL LNCNT(4)
      WRITE(6,75)
   75 FORMAT(//' IN CNTREG, THE SUBROUTINE  DPOCO HAS FOUND THE MATRIX
     1 R NOT POSITIVE DEFINITE'/)
      RETURN
C 
  100 CONTINUE
      NT = N1
      DO 150 M1 = 1,NB(1)
         CALL DPOSL(DUMMY,L,L,DUMMY(NT))
         NT = NT + L
  150 CONTINUE
      CALL EQUATE(DUMMY(N1),NDUM1,DUMMY,NDUM1)
      CALL MULT(B,NB,DUMMY(N1),NDUM1,DUMMY(N2),NA)
      CALL SCALE(DUMMY(N2),NA,DUMMY(N1),NA,-1.0D0)
      N3 = N2 + N
      IF( IDENT .OR. (IOP(1) .NE. 0) ) GO TO 200
      CALL MULT(Q,NQ,H,NH,DUMMY(N2),NH)
      CALL TRANP(H,NH,DUMMY(N3),NDUM1)
      CALL MULT(DUMMY(N3),NDUM1,DUMMY(N2),NH,Q,NQ)
C 
  200 CONTINUE
      CALL SCALE(Q,NQ,Q,NQ,-1.0D0)
      CALL JUXTR(A,NA,Q,NQ,Z,NDUM1)
      CALL TRANP(A,NA,DUMMY(N2),NA)
      CALL SCALE( DUMMY(N2),NA,DUMMY(N2),NA,-1.0D0)
      L = 2*N + 1
      CALL JUXTR(DUMMY(N1),NA,DUMMY(N2),NA,Z(L),NDUM1)
      CALL SCALE(Q,NQ,Q,NQ,-1.0D0)
      NDUM2(1) = 2*NA(1)
      NDUM2(2) = NDUM2(1)
      IF( IOP(1) .NE. 0 )  CALL PRNT(Z,NDUM2,' Z  ',1)
      CALL EQUATE(Z,NDUM2,DUMMY(N1),NDUM2)
      M = 4*N
      N2 = M + N1
      L = 2*NA(1)
      N3 = N2 + L
      N4 = N3 + L
      ISV = L
      ILV = 0
      CALL EIGEN(L,L,DUMMY(N1),DUMMY(N2),DUMMY(N3),ISV,ILV,W,DUMMY(N4),I
     1ERR)
      IF( IERR .EQ. 0 ) GO TO 300
      CALL LNCNT(4)
      IF( IERR .GT. 0 ) GO TO 250
      WRITE(6,225) IERR
  225 FORMAT(//' IN CNTREG, EIGEN FAILED TO COMPUTE THE ' ,I6 ,' EIGENV
     1ECTOR OF Z '/)
      RETURN
  250 CONTINUE
      WRITE(6,275) IERR
  275 FORMAT(//' IN CNTREG, THE ',I6 ,' EIGENVALUE OF Z HAS NOT BEEN FO
     1UND AFTER 30 ITERATIONS IN EIGEN'/)
      RETURN
C 
  300 CONTINUE
      IF( IOP(1) .EQ. 0 ) GO TO 400
      CALL LNCNT(3)
      WRITE(6,325)
  325 FORMAT(//' EIGENVALUES OF Z')
      NDUM1(1) = L
      NDUM1(2) = 2
      CALL PRNT(DUMMY(N2),NDUM1,'    ',3)
      CALL LNCNT(3)
      WRITE(6,350)
  350 FORMAT(//' CORRESPONDING EIGENVECTORS')
      CALL PRNT(W,NDUM2,'    ',3)
C 
  400 CONTINUE
      CALL EQUATE(W,NDUM2,DUMMY(N1),NDUM2)
      J1 = 1
      J2 = 1
      M = 2*N
      NDUM1(1) = L
      NDUM1(2) = 1
      K4 = N4
C 
      I=1
  415 CONTINUE
      IF( I .GT. L )  GO TO 515
      K1 = N2+I-1
      K2 = N1+(I-1)*L
      K3 = N3+I-1
      IF(DUMMY(K1) .GT. 0.0 ) GO TO 425
      J = (J1-1)*L+M+1
      J1 = J1+1
      IF(DUMMY(K3).NE. 0.0) J1=J1+1
      GO TO 450
  425 CONTINUE
      DUMMY(K4)=I
      K4 = K4+1
      J = (J2-1)*L+1
      J2 = J2+1
      IF( DUMMY(K3) .NE. 0.0 )  J2 = J2 + 1
  450 CONTINUE
      CALL EQUATE(DUMMY(K2),NDUM1,W(J),NDUM1)
      IF(DUMMY(K3) .EQ. 0.0) GO TO 500
      I = I+1
      K2 = K2+L
      J = J+L
      CALL EQUATE(DUMMY(K2),NDUM1,W(J),NDUM1)
  500 CONTINUE
      I=I+1
      GO TO 415
  515 CONTINUE
C 
      CALL NULL(LAMBDA,NA)
      K0 = -1
      J = -NA(1)
      NAX = NA(1)
      I=1
  520 CONTINUE
      IF( I .GT. NAX )  GO TO 530
      J = NAX + J + 1
      K0 = K0 + 1
      K1 = N4 + K0
      K2 = DUMMY(K1)
      K = N2+K2-1
      LAMBDA(J) = DUMMY(K)
      K3 = N3+K2-1
      IF( DUMMY(K3) .EQ. 0.0 ) GO TO 525
      K4 = J+1
      LAMBDA(K4) = -DUMMY(K3)
      K4 = K4+NAX
      LAMBDA(K4) = DUMMY(K)
      K4 = K4-1
      LAMBDA(K4) = DUMMY(K3)
      K5 = M + (I-1)*L + 1
      K6 = K5 + L
      CALL EQUATE(W(K5),NDUM1,DUMMY(N1),NDUM1)
      CALL EQUATE(W(K6),NDUM1,W(K5),NDUM1)
      CALL EQUATE(DUMMY(N1),NDUM1,W(K6),NDUM1)
      I = I+1
      J = NAX + J +1
  525 CONTINUE
      I=I+1
      GO TO 520
  530 CONTINUE
C 
      IF( IOP(1) .EQ. 0 ) GO TO 700
      CALL LNCNT(3)
      WRITE(6,535)
  535 FORMAT(//' REORDERED EIGENVECTORS')
      CALL PRNT(W,NDUM2,'    ',3)
      CALL LNCNT(4)
      WRITE(6,545)
  545 FORMAT(//' LAMBDA MATRIX OF EIGENVALUES OF Z WITH POSITIVE REAL PA
     1RTS'/)
      CALL PRNT(LAMBDA,NA,'    ',3)
C 
      CALL MULT(Z,NDUM2,W,NDUM2,DUMMY(N1),NDUM2)
      L = NDUM2(1)
      M = L**2
      N2 = N1+M
      CALL EQUATE(W,NDUM2,DUMMY(N2),NDUM2)
      N3 = N2+M
      N4 = N3+L
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY(N2),L,L,DUMMY(N3),RCOND,DUMMY(N4))
      IF ((1.0 + RCOND) .NE. 1.0) GO TO 600
      CALL LNCNT(4)
      WRITE(6,550) RCOND
  550 FORMAT(//' IN CNTREG, DGECO HAS FOUND THE REORDERED MATRIX W TO B
     1E SINGULAR, RCOND = ',D16.8,/)
  600 CONTINUE
      NT = N1
      DO 650 M1 = 1,L
         CALL DGESL(DUMMY(N2),L,L,DUMMY(N3),DUMMY(NT),0)
         NT = NT + L
  650 CONTINUE
      CALL PRNT(DUMMY(N1),NDUM2,'WIZW',1)
C 
  700 CONTINUE
      NDUM1(1) = 2*NA(1)
      NDUM1(2) = NA(1)
      N2 = 2*N + N1
      CALL TRANP(W,NDUM1,DUMMY(N2),NDUM2)
      NW11 = N1
      NDUM1(1) = NA(1)
      CALL TRANP(DUMMY(N2),NDUM1,DUMMY(NW11),NDUM1)
      L = N2+N
      NW21 = NW11+N
      CALL TRANP(DUMMY(L),NDUM1,DUMMY(NW21),NDUM1)
      L = 2*N+1
      NDUM1(1)=2*NA(1)
      N3 = N2 + 2*N
      CALL TRANP(W(L),NDUM1,DUMMY(N3),NDUM2)
      NDUM1(1) = NA(1)
      NW12 = NW21+N
      CALL TRANP(DUMMY(N3),NDUM1,DUMMY(NW12),NDUM1)
      L = N3 + N
      NW22 = NW12 + N
      CALL TRANP(DUMMY(L),NDUM1,DUMMY(NW22),NDUM1)
C 
      IF( IOP(1) .EQ. 0 ) GO TO 800
      CALL PRNT(DUMMY(NW11),NA,'W11 ',1)
      CALL PRNT(DUMMY(NW21),NA,'W21 ',1)
      CALL PRNT(DUMMY(NW12),NA,'W12 ',1)
      CALL PRNT(DUMMY(NW22),NA,'W22 ',1)
C 
  800 CONTINUE
      IF( IOP(3) .NE. 0 ) GO TO 900
      N2 = N1+4*N
      CALL MULT(P,NP,DUMMY(NW12),NA,S,NA)
      CALL MULT(P,NP,DUMMY(NW11),NA,DUMMY(N2),NA)
      CALL SUBT(S,NA,DUMMY(NW22),NA,S,NA)
      CALL SUBT(DUMMY(NW21),NA,DUMMY(N2),NA,DUMMY(N2),NA)
      N3 = N2+N
      L = NA(1)
      N4 = N3+NA(1)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY(N2),L,L,DUMMY(N3),RCOND,DUMMY(N4))
      IF ((1.0 + RCOND) .NE. 1.0) GO TO 850
      CALL LNCNT(4)
      WRITE(6,825) RCOND
  825 FORMAT(//' IN CNTREG, DGECO HAS FOUND THE MATRIX  W21 - P1XW11 TO
     1 BE SINGULAR, RCOND = ',D16.8,/)
      RETURN
C 
  850 CONTINUE
      NT = 1
      DO 860 M1 = 1,L
         CALL DGESL(DUMMY(N2),L,L,DUMMY(N3),S(NT),0)
         NT = NT + L
  860 CONTINUE
      IF( IOP(1) .EQ. 0 ) GO TO 1000
      CALL PRNT(S,NA,' S  ',1)
      NDUM1(1) = NR(1)
      NDUM1(2) = NA(1)
      CALL LNCNT(3)
      WRITE(6,875)
  875 FORMAT(//' MATRIX (R INVERSE)X(B TRANSPOSE)')
      CALL PRNT(DUMMY,NDUM1,'    ',3)
      GO TO 1000
C 
  900 CONTINUE
      N2 = N1+4*N
      CALL TRANP(DUMMY(NW12),NA,DUMMY(N2),NA)
      CALL TRANP(DUMMY(NW22),NA,P,NP)
      N3 = N2+N
      L = NA(1)
      N4 = N3 + NA(1)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY(N2),L,L,DUMMY(N3),RCOND,DUMMY(N4))
      IF ((1.0 + RCOND) .NE. 1.0) GO TO 950
      CALL LNCNT(4)
      WRITE(6,925) RCOND
  925 FORMAT(//' IN CNTREG, DGECO HAS FOUND THE MATRIX W12 TO BE SINGUL
     1AR, RCOND = ',D16.8,/)
      RETURN
  950 CONTINUE
      NT = 1
      DO 975 M1 = 1,L
         CALL DGESL(DUMMY(N2),L,L,DUMMY(N3),P(NT),0)
         NT = NT + L
  975 CONTINUE
      NDUM1(1) = NR(1)
      NDUM1(2) = NA(1)
      CALL MULT(DUMMY,NDUM1,P,NP,F,NF)
      IF( IOP(1) .EQ. 0 ) RETURN
      CALL PRNT(P,NP,' P  ',1)
      CALL PRNT(F,NF,' F  ',1)
      RETURN
C 
 1000 CONTINUE
      NMAX = T(1)/T(2)
      I = NMAX
      CALL EQUATE(LAMBDA,NA,DUMMY(N2),NA)
      TT = -T(2)
      N4 = N3+N
      N5 = N4+N
      N6 = N5+N
      N7 = N6+NA(1)
      KSS = 0
      NDUM1(1) = NR(1)
      NDUM1(2) = NA(1)
      CALL EXPSER(DUMMY(N2),NA,DUMMY(N3),NA,TT,KSS,DUMMY(N4))
      CALL EQUATE(DUMMY(N3),NA,DUMMY(N2),NA)
      IF( IOP(1) .EQ. 0 )  GO TO 1075
      CALL LNCNT(3)
      WRITE(6,1050) T(2)
 1050 FORMAT(//' EXP(-LAMBDA X ',D16.8 ,')')
      CALL PRNT(DUMMY(N2),NA,'    ',3)
 1075 CONTINUE
      IF( NMAX .LE. 0 )  RETURN
      CALL EQUATE(S,NA,DUMMY(N3),NA)
 1100 CONTINUE
      TIME = I*T(2)
      IF( I .NE. NMAX )  CALL EQUATE(DUMMY(N5),NA,P,NP)
      CALL MULT(DUMMY(N3),NA,DUMMY(N2),NA,DUMMY(N4),NA)
      CALL MULT(DUMMY(N2),NA,DUMMY(N4),NA,DUMMY(N3),NA)
      CALL MULT(DUMMY(NW11),NA,DUMMY(N3),NA,DUMMY(N4),NA)
      CALL ADD(DUMMY(NW12),NA,DUMMY(N4),NA,DUMMY(N4),NA)
      CALL TRANP(DUMMY(N4),NA,DUMMY(N5),NA)
      CALL EQUATE(DUMMY(N5),NA,DUMMY(N4),NA)
      CALL MULT(DUMMY(NW21),NA,DUMMY(N3),NA,DUMMY(N5),NA)
      CALL ADD(DUMMY(NW22),NA,DUMMY(N5),NA,DUMMY(N5),NA)
      CALL TRANP(DUMMY(N5),NA,DUMMY(N6),NA)
      CALL EQUATE(DUMMY(N6),NA,DUMMY(N5),NA)
      L = NA(1)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY(N4),L,L,DUMMY(N6),RCOND,DUMMY(N7))
      IF ((1.0 + RCOND) .NE. 1.0) GO TO 1200
      CALL LNCNT(3)
      WRITE(6,1150) TIME,RCOND
 1150 FORMAT(//' IN CNTREG AT TIME ',D16.8,'  P CANNOT BE COMPUTED DUE
     1 TO MATRIX SINGULARITY IN DGECO, RCOND = ',D16.8)
      RETURN
C 
 1200 CONTINUE
      NT = N5
      DO 1220 M1 = 1,L
         CALL DGESL(DUMMY(N4),L,L,DUMMY(N6),DUMMY(NT),0)
         NT = NT + L
 1220 CONTINUE
      CALL MAXEL(P,NP,ANORM1)
      CALL SUBT(DUMMY(N5),NA,P,NP,DUMMY(N4),NA)
      CALL MAXEL(DUMMY(N4),NA,ANORM2)
      IF( ANORM1 .NE. 0.0 ) GO TO 1225
      GO TO 1300
C 
 1225 CONTINUE
      IF(ANORM1 .GT. 1.0 ) GO TO 1250
      IF( ANORM2/ANORM1 .LT. RICTCV ) KSS=1
      GO TO 1300
 1250 CONTINUE
      IF( ANORM2 .LT. RICTCV ) KSS=1
C 
 1300 CONTINUE
      CALL MULT(DUMMY,NDUM1,P,NP,F,NF)
      IF( IOP(2) .EQ. 0 ) GO TO 1400
      CALL LNCNT(5)
      WRITE(6,1350) TIME
 1350 FORMAT(///' TIME = ',D16.8/)
      CALL PRNT(P,NP,' P  ',1)
      IF( I .NE. NMAX )  CALL PRNT(F,NF,' F  ',1)
C 
 1400 CONTINUE
      IF( KSS .EQ. 1 ) GO TO 1500
      I = I-1
      IF( I .GE. 0 ) GO TO 1100
      GO TO 1600
 1500 CONTINUE
      CALL LNCNT(4)
      WRITE(6,1550)
 1550 FORMAT(//' STEADY-STATE SOLUTION HAS BEEN REACHED IN CNTREG'/)
C 
 1600 CONTINUE
      IF( IOP(2) .NE. 0 ) RETURN
      IF( IOP(1) .EQ. 0 ) RETURN
      CALL LNCNT(5)
      WRITE(6,1350) TIME
      CALL PRNT(P,NP,' P  ',1)
      CALL PRNT(F,NF,' F  ',1)
C 
      RETURN
      END
      SUBROUTINE CSTAB(A,NA,B,NB,F,NF,IOP,SCLE,DUMMY)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      DIMENSION A(1),B(1),F(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NF(2),IOP(3),NDUM(2)
      DIMENSION IOPT(2)
      LOGICAL SYM
      COMMON/TOL/EPSAM,EPSBM,IACM
C
      N = NA(1)**2
      N1=N+1
C
      IF(IOP(2) .EQ. 0 ) GO TO 100
      CALL EQUATE(A,NA,DUMMY,NA)
      N2=N1+NA(1)
      N3=N2+NA(1)
      ISV =0
      ILV =0
      CALL EIGEN(NA(1),NA(1),DUMMY,DUMMY(N1),DUMMY(N2),ISV,ILV,V,DUMMY(N
     13),IERR)
C
      M=NA(1)
      IF(IERR .EQ. 0) GO TO 50
      CALL LNCNT(3)
      PRINT 25,IERR
   25 FORMAT(//' IN CSTAB, THE SUBROUTINE EIGEN FAILED TO DETERMINE THE
     1',I4,' EIGENVALUE FOR THE MATRIX A  AFTER 30 ITERATIONS')
      IERR=1
      CALL NORMS(M,M,M,A,IERR,BETA)
      BETA=2.*BETA
      GO TO 200
   50 CONTINUE
C
      BETA = 0.0
      DO 75 I = 1,M
      J = N1 + I - 1
      BETA1 = ABS(DUMMY(J))
      IF(BETA1 .GT. BETA)  BETA = BETA1
   75 CONTINUE
      BETA = SCLE*(BETA + .001)
      GO TO 200
C
  100 CONTINUE
      BETA = SCLE
  200 CONTINUE
C
      CALL TRANP(B,NB,DUMMY,NDUM)
      CALL MULT(B,NB,DUMMY,NDUM,DUMMY(N1),NA)
      CALL SCALE(DUMMY(N1),NA,DUMMY,NA,-2.0)
      CALL SCALE(A,NA,DUMMY(N1),NA,-1.0)
      J = -NA(1)
      NAX = NA(1)
      DO 225 I=1,NAX
      J = J+NAX+1
      K = N1+J-1
      DUMMY(K)=DUMMY(K)-BETA
  225 CONTINUE
      N2 = N1 + N
      SYM = .TRUE.
      IOPT(1)=0
C
      IF( IOP(3) .NE. 0 )  GO TO 300
      EPSA=EPSAM
      CALL BARSTW(DUMMY(N1),NA,A,NA,DUMMY,NA,IOPT,SYM,EPSA,EPSA,DUMMY(N2
     1))
      GO TO 350
  300 CONTINUE
      IOPT(2) = 1
      CALL BILIN(DUMMY(N1),NA,A,NA,DUMMY,NA,IOPT,ASCLE,SYM,DUMMY(N2))
  350 CONTINUE
C
      CALL EQUATE(B,NB,DUMMY(N1),NB)
      IOPT (1)= 3
      IAC =IACM
      N3 = N2 + NA(1)
      CALL SNVDEC(IOPT,NA(1),NA(1),NA(1),NA(1),DUMMY,NB(2),DUMMY(N1),IAC
     1,ZTEST,DUMMY(N2),DUMMY(N3),IRANK,APLUS,IERR)
      IF(IERR .EQ. 0 ) GO TO 400
      CALL LNCNT(5)
      IF(IERR .GT. 0 ) PRINT 360,IERR
      IF(IERR .EQ. -1) PRINT 370,ZTEST,IRANK
  360 FORMAT(//' IN CSTAB, SNVDEC HAS FAILED TO CONVERGE TO THE ',I4,' S
     1INGULARVALUE AFTER 30 ITERATIONS'//)
  370 FORMAT(//' IN CSTAB, THE MATRIX SUBMITTED TO SNVDEC USING ZTEST =
     1',E16.8,' IS CLOSE TO A  MATRIX OF LOWER  RANK '/' IF THE ACCURACY
     2 IAC IS REDUCED THE RANK MAY ALSO BE REDUCED'/' CURRENT RANK =',I4
     2)
      IF( IERR .GT. 0 ) RETURN
      NDUM(1) = NA(1)
      NDUM(2) =1
      CALL PRNT(DUMMY(N2),NDUM,'SGVL',1)
  400 CONTINUE
C
      CALL TRANP(DUMMY(N1),NB,F,NF)
      IF ( IOP(1) .EQ. 0 )   RETURN
      CALL LNCNT(4)
      PRINT 500
  500 FORMAT(//' COMPUTATION OF F MATRIX SUCH THAT A-BF IS ASYMPTOTICALL
     1Y STABLE IN THE CONTINUOUS SENSE '/)
      CALL PRNT(A,NA,' A  ',1)
      CALL LNCNT(4)
      PRINT 550,BETA
  550 FORMAT(//' BETA = ',E16.8/)
      CALL PRNT(B,NB,' B  ',1)
      CALL PRNT(F,NF,' F  ',1)
      CALL MULT(B,NB,F,NF,DUMMY,NA)
      CALL SUBT(A,NA,DUMMY,NA,DUMMY,NA)
      CALL PRNT(DUMMY,NA,'A-BF',1)
      N2 = N1+NA(1)
      N3 = N2+NA(1)
      ISV = 0
      ILV = 0
      CALL EIGEN(NA(1),NA(1),DUMMY,DUMMY(N1),DUMMY(N2),ISV,ILV,V,DUMMY(N
     13),IERR)
      M = NA(1)
      IF( IERR .EQ. 0 ) GO TO 600
      M = NA(1)-IERR
      CALL LNCNT(3)
      PRINT 25,IERR
  600 CONTINUE
      CALL LNCNT(4)
      PRINT 650
  650 FORMAT(//' EIGENVALUES OF A-BF'/)
  675 FORMAT(10X,2E16.8)
      CALL LNCNT(M)
      DO 700 I=1,M
      J = N1+I-1
      K = N2+I-1
      PRINT 675,DUMMY(J),DUMMY(K)
  700 CONTINUE
C
      RETURN
      END
      SUBROUTINE CTROL(A,NA,B,NB,C,NC,IOP,IAC,IRANK,DUMMY)
C 
C   PURPOSE:
C      Evaluate the controllability matrix C=[B,AB,...,(A**(n-1))B]
C      for a real constant (A,B) pair.  The matrix A is n x n and B is
C      n x r with r <= n.  Options are provided to compute both the
C      rank and singular values of C along with the controllabilty can-
C      onical form for the (A,B) pair.
C 
C   REFERENCES:
C      Kwakernaak, Huibert; and Sivan, Raphael: Linear Optimal Control
C        Systems.  John Wiley & Sons, Inc., c.1972.
C 
C   Subroutines employed by CTROL: EQUATE, JUXTC, LNCNT, MULT, PRNT,
C      SNVDEC, TRANP
C   Subroutines employing CTROL: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NC(2),NV(2),IOP(5)
C 
      N = NA(1)*NB(2)
      N1 = N+1
      N2 = N1+N
      K = NA(1)-1
      J = 1
C 
      CALL EQUATE(B,NB,DUMMY(N2),NV)
      CALL EQUATE(B,NB,DUMMY,NB)
100   CONTINUE
      CALL MULT(A,NA,DUMMY,NB,DUMMY(N1),NB)
      CALL JUXTC(DUMMY(N2),NV,DUMMY(N1),NB,C,NC)
C 
      IF( J .EQ. K ) GO TO 200
C 
      CALL EQUATE(DUMMY(N1),NB,DUMMY,NB)
      CALL EQUATE(C,NC,DUMMY(N2),NV)
      J = J + 1
      GO TO 100
C 
  200 CONTINUE
C 
      IF(IOP(1) .EQ. 0 ) GO TO 300
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      CALL LNCNT(4)
      WRITE(6,250)
  250 FORMAT(//' THE MATRIX C IS THE CONTROLLABILITY MATRIX FOR THE  A/B
     1 PAIR'/)
      CALL PRNT(C,NC,' C  ',1)
C 
  300  IF( IOP(2) .EQ. 0 ) RETURN
      NOS = 0
      IOPT = 2
      K = NC(2)
      NC(2) = NB(2)*(NA(2)-NB(2)+1)
      N = NC(1)*NC(2)
      CALL TRANP(C,NC,DUMMY,NV)
      NC(2) = K
      N1 = N + 1
      N2 = N1 + NV(2)
      CALL SNVDEC(IOPT,NV(1),NV(2),NV(1),NV(2),DUMMY,NOS,B,IAC,ZTEST,DUM
     1MY(N1),DUMMY(N2),IRANK,A,IERR)
      IF( IERR .EQ. 0 ) GO TO 340
      CALL LNCNT(5)
      IF( IERR .GT. 0 ) WRITE(6,310) IERR
      IF( IERR .EQ. -1 ) WRITE(6,320) ZTEST,IRANK
  310 FORMAT(//' IN CTROL, SNVDEC HAS FAILED TO CONVERGE TO THE ',I4 ,'
     1SINGULAR VALUE AFTER 30 ITERATIONS ')
  320 FORMAT(//' IN CTROL, THE MATRIX SUBMITTED TO SNVDEC USING ZTEST =
     1',D16.8,' IS CLOSE TO A MATRIX WHICH IS OF LOWER RANK'/' IF THE AC
     2CURACACY IS REDUCED THE RANK MAY ALSO BE REDUCED'/' CURRENT RANK =
     3',I4)
      IF( IERR .GT. 0 ) RETURN
C 
  340 CONTINUE
      IF( IOP(3) .EQ. 0 ) GO TO 400
      CALL LNCNT(6)
      WRITE(6,350) ZTEST,IRANK
  350 FORMAT(//' BASED ON THE ZERO-TEST ',D16.8,' THE RANK OF THE CONTRO
     1LLABILITY MATRIX IS ',I4/' THE SINGULAR VALUES ARE '/)
      IOPT = 0
      NV(1)= NV(2)
      NV(2)= 1
      CALL PRNT(DUMMY(N1),NV,'    ',3)
C 
  400 IF( IOP(4) .EQ. 0 ) RETURN
      N = NA(1)**2
      CALL EQUATE(DUMMY(N2),NA,DUMMY,NA)
      N1 = N + 1
      N2 = N1 + N
      CALL MULT(A,NA,DUMMY,NA,DUMMY(N1),NA)
      CALL TRANP(DUMMY,NA,DUMMY(N2),NA)
      CALL EQUATE(DUMMY(N2),NA,DUMMY,NA)
      CALL MULT(DUMMY,NA,DUMMY(N1),NA,DUMMY(N2),NA)
      CALL MULT(DUMMY,NA,B,NB,DUMMY(N1),NB)
C 
      IF( IOP(5) .EQ. 0 ) RETURN
      CALL LNCNT(5)
      WRITE(6,500)
  500 FORMAT(//' CONTROLLABILITY CANONICAL FORM '/ ' (V TRANSPOSE) A V')
      CALL PRNT(DUMMY(N2),NA,'    ',3)
      CALL LNCNT(2)
      WRITE(6,510)
  510 FORMAT(/' (V TRANSPOSE ) B ')
      CALL PRNT(DUMMY(N1),NB,'    ',3)
      CALL LNCNT(2)
      WRITE(6,520)
  520 FORMAT(/' V TRANSPOSE')
      CALL PRNT(DUMMY,NA,'    ',3)
C 
      RETURN
      END

      DOUBLE PRECISION FUNCTION DAMCON(JOB)
      INTEGER JOB
C
C     ADAPTED FROM BLAS SUBROUTINE SMACH/DMACH FOR ORACLS
C
C     SMACH COMPUTES MACHINE PARAMETERS OF FLOATING POINT
C     ARITHMETIC FOR USE IN TESTING ONLY.  NOT REQUIRED BY
C     LINPACK PROPER.
C
C     IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES,
C     THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS.
C     ASSUME THE COMPUTER HAS
C
C        B = BASE OF ARITHMETIC
C        T = NUMBER OF BASE  B  DIGITS
C        L = SMALLEST POSSIBLE EXPONENT
C        U = LARGEST POSSIBLE EXPONENT
C
C     THEN
C
C        EPS = B**(1-T)
C        TINY = 100.0*B**(-L+T)
C        HUGE = 0.01*B**(U-T)
C
C     DMACH SAME AS SMACH EXCEPT T, L, U APPLY TO
C     DOUBLE PRECISION.
C
C     CMACH SAME AS SMACH EXCEPT IF COMPLEX DIVISION
C     IS DONE BY
C
C        1/(X+I*Y) = (X-I*Y)/(X**2+Y**2)
C
C     THEN
C
C        TINY = SQRT(TINY)
C        HUGE = SQRT(HUGE)
C
C **** ORACLS REQUIRES NUMBERING REVERSED FROM THAT USED BY LINPACK
C
C     JOB IS 3, 2 OR 1 FOR EPSILON, TINY AND HUGE, RESPECTIVELY.
C
      DOUBLE PRECISION EPS,TINY,HUGE,S
C
      EPS = 1.0D0
   10 EPS = EPS/2.0D0
      S = 1.0D0 + EPS
      IF (S .GT. 1.0D0) GO TO 10
      EPS = 2.0D0*EPS
C
      S = 1.0D0
   20 TINY = S
      S = S/16.0D0
      IF (S*1.0 .NE. 0.0D0) GO TO 20
      TINY = (TINY/EPS)*100.0
      HUGE = 1.0D0/TINY
C
      IF (JOB .EQ. 3) DAMCON = EPS
      IF (JOB .EQ. 2) DAMCON = TINY
      IF (JOB .EQ. 1) DAMCON = HUGE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DTEMP
      INTEGER I,INCX,M,MP1,N,NINCX
C
      DASUM = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DTEMP = DTEMP + DABS(DX(I))
   10 CONTINUE
      DASUM = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I + 1)) + DABS(DX(I + 2))
     *  + DABS(DX(I + 3)) + DABS(DX(I + 4)) + DABS(DX(I + 5))
   50 CONTINUE
   60 DASUM = DTEMP
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IXIY,M,MP1,N   ! ixiy never referenced
C
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
      SUBROUTINE DETFAC(NMAX,N,A,IPIVOT,IDET,DETERM,ISCALE,WK,IERR)
C 
C   PURPOSE:
C      Factor a real square matrix A as PA=LU, where P is a permutation
C      matrix representing row pivotal strategy, L is a unit lower tri-
C      angular matrix, and U is an upper triangular matrix.  Options
C      are provided to compute the determinant of A with and without
C      A input in factored form.
C 
C   Subroutines employed by DETFAC: None
C   Subroutines employing DETFAC: GELIM
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NMAX,1),IPIVOT(1),WK(1)
C 
      DATA R1,R2/1,10/
C 
      ISCALE=0
      NM1=N-1
      IERR=0
C 
C     DETERMINANT CALCULATION TEST
C 
      IF(IDET.EQ.1)GO TO 230
C 
C     TEST FOR A SCALAR MATRIX
C 
      IF(NM1.GT.0)GO TO 20
      DETERM=A(1,1)
      RETURN
C 
C     COMPUTE SCALING FACTORS
C 
   20 DO 60 I=1,N
      P=0.0
      DO 30 J=1,N
      Q=DMAX1(P,DABS(A(I,J)))
      IF(Q.GT.P)P=Q
   30 CONTINUE
      IF(P)60,40,60
   60 WK(I)=P
C 
      DO 210 M=1,NM1
C 
C     PIVOTAL LOGIC SETUP
C 
      P=0.0
      DO 110 I=M,N
      Q=DABS(A(I,M)/WK(I))
      IF(Q-P)110,110,100
  100 P=Q
      IP=I
  110 CONTINUE
C 
      IPIVOT(M)=IP
C 
      IF(P.EQ.0.)GO TO 40
      IF(M.EQ.IP)GO TO 155
C 
C     PIVOT THE M-TH ROW OF THE A MATRIX
C 
      DO 150 I=1,N
      P=A(IP,I)
      A(IP,I)=A(M,I)
  150 A(M,I)=P
C 
      P=WK(IP)
      WK(IP)=WK(M)
      WK(M)=P
C 
  155 MP1=M+1
C 
C      L/U FACTORIZATION LOGIC
C 
      P=A(M,M)
      DO 180 I=MP1,N
      A(I,M)=A(I,M)/P
      Q=A(I,M)
      DO 180 K=MP1,N
  180 A(I,K)=A(I,K)-Q*A(M,K)
C 
  210 CONTINUE
C 
      IPIVOT(N)=N
      IF (A(N,N) .EQ. 0.0) GO TO 40
C 
C     CALCULATION OF THE DETERMINANT OF A
C 
      IF(IDET.EQ.0)RETURN
C 
  230 CONTINUE
      DETERM=1.0
C 
C     ADJUST SIGN OF DETERMINANT DUE TO PIVOTAL STRATEGY
C 
      DO 250 I=1,NM1
      IF(I-IPIVOT(I))240,250,240
  240 DETERM = - DETERM
  250 CONTINUE
C 
      DO 340 I=1,N
C 
C 
  290 DETERM = DETERM * A(I,I)
C 
  300 CONTINUE
      IF(R1.GT.DABS(DETERM))GO TO 320
      DETERM=DETERM*R2
      ISCALE=ISCALE+1
      GO TO 300
C 
  320 CONTINUE
      IF(R2.LT.DABS(DETERM))GO TO 340
      DETERM=DETERM*R1
      ISCALE=ISCALE-1
      GO TO 320
C 
  340 CONTINUE
C 
      RETURN
   40 DETERM=0.0
      IERR=1
      RETURN
      END
      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(*)
      DOUBLE PRECISION A(LDA,1),Z(*)
      DOUBLE PRECISION RCOND
C
C     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DGEFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,DSIGN
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,1)
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,1),B(1)
C
C     DGESL SOLVES THE DOUBLE PRECISION SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0
C        OR DGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DISREG(A,NA,B,NB,H,NH,Q,NQ,R,NR,F,NF,P,NP,IOP,IDENT,DU
     1MMY)
C
C   PURPOSE:
C      Solve the time-invariant discrete-time linear optimal output
C      regulator problem with noise-free measurements.
C
C   REFERENCES:
C      Astrom, Karl J.: Introduction to Stochastic Control Theory.
C        Academic Press, Inc., 1970.
C
C   Subroutines employed by DISREG: ADD, DSICO, DSISL, EQUATE, LNCNT,
C      MAXEL, MULT, PRNT, SUBT, TRANP
C   Subroutines employing DISREG: ASYREG
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),Q(*),R(*),F(*),P(*),DUMMY(*)
      DIMENSION NA(2),NB(2),NQ(2),NR(2),NF(2),NP(2)
      DIMENSION IOP(3)
      DIMENSION H(*),NH(2),NDUM(2)   ! ndum never referenced
      LOGICAL  IDENT
      COMMON/TOL/EPSAM,EPSBM,IACM
C      COMMON/CONV/SUMCV,MAXSUM,RICTCV,SERCV
      COMMON/CONV/SUMCV,RICTCV,SERCV,MAXSUM
C
      N = NA(1)**2
      N1= N +1
      N2= N1+N
      N3= N2+N
C
      KSS = 0
      I=IOP(3)
C
      IF(IOP(1) .EQ. 0)  GO TO 85
      CALL LNCNT(5)
      WRITE(6,25)
   25 FORMAT(//' PROGRAM TO SOLVE THE TIME-INVARIANT FINITE-DURATION OPT
     1IMAL'/' DIGITAL REGULATOR PROBLEM WITH NOISE-FREE MEASUREMENTS'/)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      CALL PRNT(Q,NQ,' Q  ',1)
      IF( .NOT. IDENT )  GO TO 45
      CALL LNCNT(3)
      WRITE(6,35)
   35 FORMAT(/' H IS AN IDENTITY MATRIX'/)
      GO TO 65
   45 CONTINUE
      CALL PRNT(H,NH,' H  ',1)
      CALL MULT(Q,NQ,H,NH,DUMMY,NH)
      CALL TRANP(H,NH,DUMMY(N1),NF)
      CALL MULT(DUMMY(N1),NF,DUMMY,NH,Q,NQ)
      CALL LNCNT(3)
      WRITE(6,55)
   55 FORMAT(/' MATRIX ( H TRANSPOSE )QH'/)
      CALL PRNT(Q,NQ,'HTQH',1)
   65 CONTINUE
      CALL PRNT(R,NR,' R  ',1)
      CALL LNCNT(4)
      WRITE(6,75)
   75 FORMAT(//' WEIGHTING ON TERMINAL VALUE OF STATE VECTOR'/)
      CALL PRNT(P,NP,' P  ',1)
C
   85 CONTINUE
      IF((IOP(1) .NE. 0)  .OR. IDENT)  GO TO 100
      CALL MULT(Q,NQ,H,NH,DUMMY,NH)
      CALL TRANP(H,NH,DUMMY(N1),NF)
      CALL MULT(DUMMY(N1),NF,DUMMY,NH,Q,NQ)
C
  100 CONTINUE
      I=I-1
      CALL EQUATE(P,NP,DUMMY,NP)
      CALL MULT(P,NP,A,NA,DUMMY(N1),NA)
      CALL TRANP(B,NB,DUMMY(N2),NF)
      CALL MULT(DUMMY(N2),NF,DUMMY(N1),NA,F,NF)
      CALL MULT(P,NP,B,NB,DUMMY(N1),NB)
      CALL MULT(DUMMY(N2),NF,DUMMY(N1),NB,DUMMY(N3),NR)
      CALL ADD(R,NR,DUMMY(N3),NR,DUMMY(N1),NR)
      IOPT = 3
      IAC=IACM
      MF = NR(1)
C
C   * * *  CALL TO MATHLIB FUNCTIONS * * *
      CALL DSICO(DUMMY(N1),MF,MF,DUMMY(N2),RCOND,DUMMY(N3))
      IF ((1.0 + RCOND) .NE. 1.0) GO TO 300
      CALL LNCNT(5)
      WRITE(6,200) I,RCOND
  200 FORMAT(//' IN DISREG, THE MATRIX R + B(TRANS)*P(',I2,')*B IS SIN',
     1 'GULAR TO WORKING PRECISION, RCOND = ',D16.8)
      RETURN
C
  300 CONTINUE
      NT = 1
      DO 325 M1 = 1,NF(2)
         CALL DSISL(DUMMY(N1),MF,MF,DUMMY(N2),F(NT))
         NT = NT + MF
  325 CONTINUE
      CALL MULT(R,NR,F,NF,DUMMY(N1),NF)
      CALL TRANP(F,NF,DUMMY(N2),NB)
      CALL MULT(DUMMY(N2),NB,DUMMY(N1),NF,P,NP)
      CALL ADD(Q,NQ,P,NP,P,NP)
      CALL MULT(B,NB,F,NF,DUMMY(N1),NA)
      CALL SUBT(A,NA,DUMMY(N1),NA,DUMMY(N1),NA)
      CALL MULT(DUMMY,NA,DUMMY(N1),NA,DUMMY(N2),NA)
      CALL TRANP(DUMMY(N1),NA,DUMMY(N3),NA)
      CALL MULT(DUMMY(N3),NA,DUMMY(N2),NA,DUMMY(N1),NA)
      CALL ADD(P,NP,DUMMY(N1),NA,P,NP)
C
      IF( IOP(2) .EQ. 0 )  GO TO 400
      CALL LNCNT(5)
      WRITE(6,350) I
  350 FORMAT(///' STAGE ',I5,/)
      CALL PRNT(F,NF,' F  ',1)
      CALL PRNT(P,NP,' P  ',1)
C
  400 CONTINUE
      IF( I .EQ. 0 )  GO TO 600
      CALL MAXEL(DUMMY,NP,ANORM1)
      CALL SUBT(DUMMY,NP,P,NP,DUMMY(N2),NP)
      CALL MAXEL(DUMMY(N2),NP,ANORM2)
      IF( ANORM1 .NE. 0.0 )  GO TO 500
      GO TO 100
C
  500 CONTINUE
      IF(ANORM1 .GT. 1.0 ) GO TO 550
      IF( ANORM2/ANORM1 .LT. RICTCV ) KSS = 1
      GO TO 575
  550 CONTINUE
      IF( ANORM2 .LT. RICTCV ) KSS=1
  575 CONTINUE
      IF( KSS .EQ. 1) GO TO 600
      GO TO 100
C
  600 CONTINUE
      K = IOP(1) + IOP(2)
      IF( K .EQ. 0 ) RETURN
      IF( KSS .EQ. 0) GO TO 700
      CALL LNCNT(4)
      WRITE(6,650)
  650 FORMAT(//' STEADY-STATE SOLUTION HAS BEEN REACHED IN  DISREG'/)
C
  700 CONTINUE
      IF( IOP(2) .NE. 0 )  RETURN
      IF( IOP(1) .EQ. 0 )  RETURN
      CALL LNCNT(3)
      I = IOP(3)-I
      WRITE(6,800) I
  800 FORMAT(/' F AND P AFTER ',I5 ,' STEPS'/)
      CALL PRNT(F,NF,' F  ',1)
      CALL PRNT(P,NP,' P  ',1)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER          NEXT
      DOUBLE PRECISION   DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
      SUBROUTINE DPOCO(A,LDA,N,RCOND,Z,INFO)
      INTEGER LDA,N,INFO
      DOUBLE PRECISION A(LDA,1),Z(1)
      DOUBLE PRECISION RCOND
C
C     DPOCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
C     MATRIX AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DPOFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DPOCO BY DPOSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DPOCO BY DPOSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DPOCO BY DPODI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DPOCO BY DPODI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
C                DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
C                WHERE  TRANS(R)  IS THE TRANSPOSE.
C                THE STRICT LOWER TRIANGLE IS UNALTERED.
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C                IF  INFO .NE. 0 , Z  IS UNCHANGED.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DPOFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,DREAL,DSIGN
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER I,J,JM1,K,KB,KP1
C
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
      DO 30 J = 1, N
         Z(J) = DASUM(J,A(1,J),1)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + DABS(A(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
C
C     FACTOR
C
      CALL DPOFA(A,LDA,N,INFO)
      IF (INFO .NE. 0) GO TO 180
C
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C        SOLVE TRANS(R)*W = E
C
         EK = 1.0D0
         DO 50 J = 1, N
            Z(J) = 0.0D0
   50    CONTINUE
         DO 110 K = 1, N
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
            IF (DABS(EK-Z(K)) .LE. A(K,K)) GO TO 60
               S = A(K,K)/DABS(EK-Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = DABS(WK)
            SM = DABS(WKM)
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
            KP1 = K + 1
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + DABS(Z(J)+WKM*A(K,J))
                  Z(J) = Z(J) + WK*A(K,J)
                  S = S + DABS(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*A(K,J)
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
C
C        SOLVE R*Y = W
C
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 120
               S = A(K,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  130    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
C
         YNORM = 1.0D0
C
C        SOLVE TRANS(R)*V = Y
C
         DO 150 K = 1, N
            Z(K) = Z(K) - DDOT(K-1,A(1,K),1,Z(1),1)
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 140
               S = A(K,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/A(K,K)
  150    CONTINUE
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
C        SOLVE R*Z = V
C
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. A(K,K)) GO TO 160
               S = A(K,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  170    CONTINUE
C        MAKE ZNORM = 1.0
         S = 1.0D0/DASUM(N,Z,1)
         CALL DSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
  180 CONTINUE
      RETURN
      END
      SUBROUTINE DPOFA(A,LDA,N,INFO)
      INTEGER LDA,N,INFO
      DOUBLE PRECISION A(LDA,1)
C
C     DPOFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
C     MATRIX.
C
C     DPOFA IS USUALLY CALLED BY DPOCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DPOCO) = (1 + 18/N)*(TIME FOR DPOFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
C                DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
C                WHERE  TRANS(R)  IS THE TRANSPOSE.
C                THE STRICT LOWER TRIANGLE IS UNALTERED.
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DDOT
C     FORTRAN DSQRT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION S
      INTEGER J,JM1,K
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            S = A(J,J) - S
C     ......EXIT
            IF (S .LE. 0.0D0) GO TO 40
            A(J,J) = DSQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
      SUBROUTINE DPOSL(A,LDA,N,B)
      INTEGER LDA,N
      DOUBLE PRECISION A(LDA,1),B(1)
C
C     DPOSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
C     SYSTEM A * X = B
C     USING THE FACTORS COMPUTED BY DPOCO OR DPOFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DPOCO OR DPOFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
C        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
C        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
C        CORRECTLY AND  INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DPOCO(A,LDA,N,RCOND,Z,INFO)
C           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DPOSL(A,LDA,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB
C
C     SOLVE TRANS(R)*Y = B
C
      DO 10 K = 1, N
         T = DDOT(K-1,A(1,K),1,B(1),1)
         B(K) = (B(K) - T)/A(K,K)
   10 CONTINUE
C
C     SOLVE R*X = Y
C
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   20 CONTINUE
      RETURN
      END

      SUBROUTINE  DROT (N,DX,INCX,DY,INCY,C,S)
C
C     APPLIES A PLANE ROTATION.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP,C,S
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = C*DX(IX) + S*DY(IY)
        DY(IY) = C*DY(IY) - S*DX(IX)
        DX(IX) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        DTEMP = C*DX(I) + S*DY(I)
        DY(I) = C*DY(I) - S*DX(I)
        DX(I) = DTEMP
   30 CONTINUE
      RETURN
      END
      SUBROUTINE DROTG(DA,DB,C,S)
C
C     CONSTRUCT GIVENS PLANE ROTATION.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DB,C,S,ROE,SCALE,R,Z
C
      ROE = DB
      IF( DABS(DA) .GT. DABS(DB) ) ROE = DA
      SCALE = DABS(DA) + DABS(DB)
      IF( SCALE .NE. 0.0D0 ) GO TO 10
         C = 1.0D0
         S = 0.0D0
         R = 0.0D0
         GO TO 20
   10 R = SCALE*DSQRT((DA/SCALE)**2 + (DB/SCALE)**2)
      R = DSIGN(1.0D0,ROE)*R
      C = DA/R
      S = DB/R
   20 Z = 1.0D0
      IF( DABS(DA) .GT. DABS(DB) ) Z = S
      IF( DABS(DB) .GE. DABS(DA) .AND. C .NE. 0.0D0 ) Z = 1.0D0/C
      DA = R
      DB = Z
      RETURN
      END
      SUBROUTINE  DSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DSICO(A,LDA,N,KPVT,RCOND,Z)
      INTEGER LDA,N,KPVT(1)
      DOUBLE PRECISION A(LDA,1),Z(1)
      DOUBLE PRECISION RCOND
C
C     DSICO FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION
C     WITH SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE
C     MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DSIFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DSICO BY DSISL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSICO BY DSISL.
C     TO COMPUTE  INVERSE(A) , FOLLOW DSICO BY DSIDI.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSICO BY DSIDI.
C     TO COMPUTE  INERTIA(A), FOLLOW DSICO BY DSIDI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     OUTPUT
C
C        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DSIFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,IABS,DSIGN
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,EK,T
      DOUBLE PRECISION ANORM,S,DASUM,YNORM
      INTEGER I,INFO,J,JM1,K,KP,KPS,KS
C
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
      DO 30 J = 1, N
         Z(J) = DASUM(J,A(1,J),1)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + DABS(A(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
C
C     FACTOR
C
      CALL DSIFA(A,LDA,N,KPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE U*D*W = E
C
      EK = 1.0D0
      DO 50 J = 1, N
         Z(J) = 0.0D0
   50 CONTINUE
      K = N
   60 IF (K .EQ. 0) GO TO 120
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = IABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 90
               S = DABS(A(K,K))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   90       CONTINUE
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 110
  100    CONTINUE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
      GO TO 60
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(U)*Y = W
C
      K = 1
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE U*D*V = Y
C
      K = N
  170 IF (K .EQ. 0) GO TO 230
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = IABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL DAXPY(K-KS,Z(K),A(1,K),1,Z(1),1)
            IF (KS .EQ. 2) CALL DAXPY(K-KS,Z(K-1),A(1,K-1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 200
               S = DABS(A(K,K))/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
            IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         GO TO 220
  210    CONTINUE
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = Z(K)/A(K-1,K)
            BKM1 = Z(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
      GO TO 170
  230 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE TRANS(U)*Z = V
C
      K = 1
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + DDOT(K-1,A(1,K),1,Z(1),1)
            IF (KS .EQ. 2)
     *         Z(K+1) = Z(K+1) + DDOT(K-1,A(1,K+1),1,Z(1),1)
            KP = IABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         K = K + KS
      GO TO 240
  270 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
      SUBROUTINE DSIFA(A,LDA,N,KPVT,INFO)
      INTEGER LDA,N,KPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
C
C     DSIFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION
C     WITH SYMMETRIC PIVOTING.
C
C     TO SOLVE  A*X = B , FOLLOW DSIFA BY DSISL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DSIFA BY DSISL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DSIFA BY DSIDI.
C     TO COMPUTE  INERTIA(A) , FOLLOW DSIFA BY DSIDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DSIFA BY DSIDI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA,N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.
C                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH
C                WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        KPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
C                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE,
C                     BUT IT DOES INDICATE THAT DSISL OR DSIDI MAY
C                     DIVIDE BY ZERO IF CALLED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSWAP,IDAMAX
C     FORTRAN DABS,DMAX1,DSQRT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP,IDAMAX
      LOGICAL SWAP
C
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
      ALPHA = (1.0D0 + DSQRT(17.0D0))/8.0D0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (A(1,1) .EQ. 0.0D0) INFO = 1
C     ......EXIT
            GO TO 200
   20    CONTINUE
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         ABSAKK = DABS(A(K,K))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = IDAMAX(K-1,A(1,K),1)
         COLMAX = DABS(A(IMAX,K))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0D0
            IMAXP1 = IMAX + 1
            DO 40 J = IMAXP1, K
               ROWMAX = DMAX1(ROWMAX,DABS(A(IMAX,J)))
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = IDAMAX(IMAX-1,A(1,IMAX),1)
               ROWMAX = DMAX1(ROWMAX,DABS(A(JMAX,IMAX)))
   50       CONTINUE
            IF (DABS(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (DMAX1(ABSAKK,COLMAX) .NE. 0.0D0) GO TO 100
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 X 1 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 120
C
C              PERFORM AN INTERCHANGE.
C
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  T = A(J,K)
                  A(J,K) = A(IMAX,J)
                  A(IMAX,J) = T
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            DO 130 JJ = 1, KM1
               J = K - JJ
               MULK = -A(J,K)/A(K,K)
               T = MULK
               CALL DAXPY(J,T,A(1,K),1,A(1,J),1)
               A(J,K) = MULK
  130       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTERCHANGE.
C
               CALL DSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  T = A(J,K-1)
                  A(J,K-1) = A(IMAX,J)
                  A(IMAX,J) = T
  150          CONTINUE
               T = A(K-1,K)
               A(K-1,K) = A(IMAX,K)
               A(IMAX,K) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = A(K,K)/A(K-1,K)
               AKM1 = A(K-1,K-1)/A(K-1,K)
               DENOM = 1.0D0 - AK*AKM1
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  BK = A(J,K)/A(K-1,K)
                  BKM1 = A(J,K-1)/A(K-1,K)
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = MULK
                  CALL DAXPY(J,T,A(1,K),1,A(1,J),1)
                  T = MULKM1
                  CALL DAXPY(J,T,A(1,K-1),1,A(1,J),1)
                  A(J,K) = MULK
                  A(J,K-1) = MULKM1
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
      SUBROUTINE DSISL(A,LDA,N,KPVT,B)
      INTEGER LDA,N,KPVT(1)
      DOUBLE PRECISION A(LDA,1),B(1)
C
C     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM
C     A * X = B
C     USING THE FACTORS COMPUTED BY DSIFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA,N)
C                THE OUTPUT FROM DSIFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        KPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DSIFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO MAY OCCUR IF  DSICO  HAS SET RCOND .EQ. 0.0
C        OR  DSIFA  HAS SET INFO .NE. 0  .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DSIFA(A,LDA,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DSISL(A,LDA,N,KPVT,C(1,J))
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C     FORTRAN IABS
C
C     INTERNAL VARIABLES.
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP
      INTEGER K,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
      K = N
   10 IF (K .EQ. 0) GO TO 80
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY(K-1,B(K),A(1,K),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/A(K,K)
            K = K - 1
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 2) GO TO 60
               KP = IABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY(K-2,B(K),A(1,K),1,B(1),1)
               CALL DAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = B(K)/A(K-1,K)
            BKM1 = B(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)
               B(K+1) = B(K+1) + DDOT(K-1,A(1,K+1),1,B(1),1)
               KP = IABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
      SUBROUTINE DSTAB(A,NA,B,NB,F,NF,SING,IOP,SCLE,DUMMY)
C 
C   PURPOSE:
C      Primary use in ORACLS is to generate a stabilizing gain matrix
C      for initializing the quasilinearization method for solving the
C      discrete steady-state Riccati Equation.
C 
C   REFERENCES:
C      Hewer, Gary A.: An Iterative Technique for the Computation of the
C        Steady State Gains for the Discrete Optimal Regulator. IEEE
C        Trans. Autom. Control, vol. AC-16, no. 4, Aug. 1971, pp. 382-38
C      Armstrong, Ernest S.; and Rublein, George T.: A Discrete Analog
C        of the Extended Bass Algorithm for Stabilizing Constant Linear
C        Systems.  Proceedings of the 1976 IEEE Conference on Decision
C        and Control Including the 15th Symposium on Adaptive Processes,
C        76CH1150-2CS, Dec. 1976, pp. 1129-1131.
C      Armstrong, Ernest S.; and Rublein, George T.: A Stabilization
C        Algorithm for Linear Discrete Constant Systems. IEEE Trans.
C        Autom. Control, vol. AC-21, no. 4, Aug. 1976, pp. 629-631.
C 
C   Subroutines employed by DSTAB: ADD, BARSTW, CSTAB, DGECO, DGESL,
C      DLSQHTM, EIGEN, EQUATE, JUXTC, LNCNT, MULT, PRNT, SCALE, SUBT,
C      TRANP
C   Subroutines employing DSTAB: ASYREG
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),F(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NF(2),NDUM(2),IOP(2),IOPT(3),NDUM1(2)
      LOGICAL SING,SYM
      COMMON/TOL/EPSAM,EPSBM,IACM
C 
      N = NA(1)**2
      N1 = N + 1
      N2 = N1 + N
      IF( .NOT. SING ) GO TO 100
      IOPT(1)=IOP(1)
      IOPT(2) = 1
      IOPT(3) = 0
      CSCLE=1.05
      CALL CSTAB(A,NA,B,NB,F,NF,IOPT,CSCLE,DUMMY)
      CALL MULT(B,NB,F,NF,DUMMY,NA)
      CALL SUBT(A,NA,DUMMY,NA,DUMMY,NA)
      CALL EQUATE(DUMMY,NA,DUMMY(N1),NA)
      GO TO 200
C 
  100 CONTINUE
      CALL EQUATE(A,NA,DUMMY,NA)
      CALL EQUATE(A,NA,DUMMY(N1),NA)
C 
  200 CONTINUE
      IF( IOP(2) .EQ. 0 ) GO TO 300
      N3 = N2 + NA(1)
      N4 = N3 + NA(1)
      ISV = 0
      CALL EIGEN(NA(1),NA(1),DUMMY(N1),DUMMY(N2),DUMMY(N3),ISV,ISV,V,DUM
     1MY(N4),IERR)
      CALL EQUATE(DUMMY,NA,DUMMY(N1),NA)
      M = NA(1)
      IF( IERR .EQ. 0 ) GO TO 250
      CALL LNCNT(3)
      WRITE(6,225)
  225 FORMAT(//' IN DSTAB, THE PROGRAM EIGEN FAILED TO DETERMINE THE ',I
     15,' EIGENVALUE FOR THE MATRIX A-BG AFTER 30 ITERATIONS')
      CALL PRNT(DUMMY,NA,'A-BG',1)
      IF( SING ) CALL PRNT(F,NF,' G  ',1)
      RETURN
C 
  250 CONTINUE
      ALPHA = 1.0
      DO 275 I =1,M
      I1 = N2 + I -1
      I2 = N3 + I -1
      ALPHA1 = DSQRT(DUMMY(I1)**2 + DUMMY(I2)**2)
      IF( ALPHA1 .LT. ALPHA .AND. ALPHA1 .NE. 0 ) ALPHA = ALPHA1
  275 CONTINUE
      ALPHA = SCLE*ALPHA
      GO TO 400
C 
  300 CONTINUE
      ALPHA = SCLE
C 
  400 CONTINUE
      J = -NA(1)
      NAX = NA(1)
      DO 425 I = 1,NAX
      J = J + NAX + 1
      K = N1 + J -1
      DUMMY(K) = DUMMY(J) - ALPHA
      DUMMY(J) = DUMMY(J) + ALPHA
  425 CONTINUE
      CALL EQUATE(B,NB,DUMMY(N2),NB)
      N3 = N2 + NA(1)*NB(2)
      NRHS = NA(1)+NB(2)
      N4 = N3 + NA(1)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY,NA(1),NA(1),DUMMY(N3),RCOND,DUMMY(N4))
      IF ((1.0 + RCOND) .NE. 1.0) GO TO 500
      CALL LNCNT(3)
      IF( .NOT. SING ) GO TO 445
      WRITE(6,435) RCOND
  435 FORMAT(//' IN DSTAB, DGECO HAS FOUND THE MATRIX ( A-BG) + (ALPHA)I
     1 SINGULAR, RCOND = ',D16.8)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(F,NF,' G  ',1)
      GO TO 465
  445 CONTINUE
      CALL LNCNT(3)
      WRITE(6,455) RCOND
  455 FORMAT(//' IN DSTAB, DGECO HAS FOUND THE MATRIX A + (ALPHA)I SINGU
     1LAR, RCOND = ',D16.8)
      CALL PRNT(A,NA,' A  ',1)
  465 CONTINUE
      CALL LNCNT(3)
      WRITE(6,475) ALPHA
  475 FORMAT(//' ALPHA = ',D16.8)
      RETURN
C 
  500 CONTINUE
      NT = N1
      DO 520 M1 = 1,NRHS
         CALL DGESL(DUMMY,NA(1),NA(1),DUMMY(N3),DUMMY(NT),0)
         NT = NT + NA(1)
  520 CONTINUE
      CALL EQUATE(DUMMY(N1),NA,DUMMY,NA)
      CALL TRANP(DUMMY(N2),NB,DUMMY(N1),NDUM)
      N3 = N2 + N
      CALL MULT(DUMMY(N2),NB,DUMMY(N1),NDUM,DUMMY(N3),NA)
      CALL SCALE(DUMMY(N3),NA,DUMMY(N1),NA,4.0D0)
      SYM = .TRUE.
      IOPT(1) = 0
      EPSA=EPSAM
      CALL BARSTW(DUMMY,NA,B,NB,DUMMY(N1),NA,IOPT,SYM,EPSA,EPSA,DUMMY(N2
     1))
      CALL EQUATE(DUMMY(N1),NA,DUMMY,NA)
      CALL TRANP(B,NB,DUMMY(N1),NDUM)
      CALL MULT(B,NB,DUMMY(N1),NDUM,DUMMY(N2),NA)
      CALL ADD(DUMMY,NA,DUMMY(N2),NA,DUMMY,NA)
      CALL EQUATE(A,NA,DUMMY(N1),NA)
      IF( .NOT. SING ) GO TO 600
      CALL MULT(B,NB,F,NF,DUMMY(N1),NA)
      CALL SUBT(A,NA,DUMMY(N1),NA,DUMMY(N1),NA)
C 
  600 CONTINUE
      IOPT(1) = 3
      M = NA(1)
      IAC=IACM
      CALL SNVDEC(IOPT,M,M,M,M,DUMMY,M,DUMMY(N1),IAC,ZTEST,DUMMY(N2),DUM
     1MY(N3),IRANK,APLUS,IERR)
      IF( IERR  .EQ. 0 ) GO TO 700
      CALL LNCNT(5)
      IF( IERR .GT. 0 ) PRINT 625,IERR
      IF( IERR .EQ. -1) PRINT 650,ZTEST,IRANK
  625 FORMAT(//' IN DSTAB, SNVDEC HAS FAILED TO CONVERGE TO THE ',I5,' S
     1INGULAR VALUE AFTER 30 ITERATIONS')
  650 FORMAT(//' IN DSTAB, THE MATRIX SUBMITTED TO SNVDEC, USING ZTEST =
     1 ',E16.8,' , IS CLOSE TO A MATRIX OF LOWER RANK'/' IF THE ACCURACY
     2 IAC IS REDUCED THE RANK MAY ALSO BE REDUCED'/' CURRENT RANK = ',I
     24)
      IF( IERR  .GT. 0 ) RETURN
      NDUM(1)= NA(1)
      NDUM(2)= 1
      CALL PRNT(DUMMY(N2),NDUM,'SGVL',1)
C
  700 CONTINUE
      CALL TRANP(B,NB,DUMMY(N2),NDUM)
      CALL MULT(DUMMY(N2),NDUM,DUMMY(N1),NA,DUMMY,NF)
      IF( .NOT. SING ) GO TO 800
      CALL ADD(F,NF,DUMMY,NF,F,NF)
      GO TO 900
C 
  800 CONTINUE
      CALL EQUATE(DUMMY,NF,F,NF)
C 
  900 CONTINUE
      IF( IOP(1) .EQ. 0 ) RETURN
      CALL LNCNT(4)
      WRITE(6,1000)
 1000 FORMAT(//' COMPUTATION OF F SUCH THAT A-BF IS ASYMPTOTICALLY STABL
     1E IN THE DISCRETE SENSE'/)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      CALL LNCNT(4)
      WRITE(6,1100) ALPHA
 1100 FORMAT(//' ALPHA = ',D16.8/)
      CALL PRNT(F,NF,' F  ',1)
      CALL MULT(B,NB,F,NF,DUMMY,NA)
      CALL SUBT(A,NA,DUMMY,NA,DUMMY,NA)
      CALL PRNT(DUMMY,NA,'A-BF',1)
      CALL LNCNT(3)
      WRITE(6,1200)
 1200 FORMAT(//' EIGENVALUES OF A-BF')
      NDUM(1) = NA(1)
      NDUM(2) = 1
      N2 = N1 + NA(1)
      N3 = N2 + NA(1)
      ISV = 0
      CALL EIGEN(NA(1),NA(1),DUMMY,DUMMY(N1),DUMMY(N2),ISV,ISV,V,DUMMY(N
     13),IERR)
      IF( IERR .EQ. 0 ) GO TO 1300
      CALL LNCNT(3)
      WRITE(6,1250)
 1250 FORMAT(//' IN DSTAB, THE PROGRAM EIGEN FAILED TO DETERMINE THE ',I
     15,' EIGENVALUE FOR THE A-BF MATRIX AFTER 30 ITERATIONS')
      NDUM(1)=NA(1)-IERR
C 
 1300 CONTINUE
      CALL JUXTC(DUMMY(N1),NDUM,DUMMY(N2),NDUM,DUMMY,NDUM1)
      CALL PRNT(DUMMY,NDUM1,'EIGN',1)
      CALL LNCNT(4)
      WRITE(6,1400)
 1400 FORMAT(//' MODULI OF EIGENVALUES OF A-BF'/)
      M =NDUM(1)
      DO 1500 I = 1,M
      J = N1 + I - 1
      K = N2 + I - 1
      DUMMY(I)=DSQRT(DUMMY(J)**2 + DUMMY(K)**2)
 1500 CONTINUE
      CALL PRNT(DUMMY,NDUM,'MOD ',1)
C 
      RETURN
      END
      SUBROUTINE DSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      DOUBLE PRECISION X(LDX,1),S(1),E(1),U(LDU,1),V(LDV,1),WORK(1)
C
C
C     DSVDC IS A SUBROUTINE TO REDUCE A DOUBLE PRECISION NXP MATRIX X
C     BY ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
C     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
C     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
C     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
C
C     ON ENTRY
C
C         X         DOUBLE PRECISION(LDX,P), WHERE LDX.GE.N.
C                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
C                   DECOMPOSITION IS TO BE COMPUTED.  X IS
C                   DESTROYED BY DSVDC.
C
C         LDX       INTEGER.
C                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C         N         INTEGER.
C                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X.
C
C         P         INTEGER.
C                   P IS THE NUMBER OF ROWS OF THE MATRIX X.
C
C         LDU       INTEGER.
C                   LDU IS THE LEADING DIMENSION OF THE ARRAY U.
C                   (SEE BELOW).
C
C         LDV       INTEGER.
C                   LDV IS THE LEADING DIMENSION OF THE ARRAY V.
C                   (SEE BELOW).
C
C         WORK      DOUBLE PRECISION(N).
C                   WORK IS A SCRATCH ARRAY.
C
C         JOB       INTEGER.
C                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR
C                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB
C                   WITH THE FOLLOWING MEANING
C
C                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR
C                                  VECTORS.
C                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS
C                                  IN U.
C                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR
C                                  VECTORS IN U.
C                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
C                                  VECTORS.
C                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
C                                  IN V.
C
C     ON RETURN
C
C         S         DOUBLE PRECISION(MM), WHERE MM=MIN(N+1,P).
C                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
C                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
C                   ORDER OF MAGNITUDE.
C
C         E         DOUBLE PRECISION(P).
C                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
C                   DISCUSSION OF INFO FOR EXCEPTIONS.
C
C         U         DOUBLE PRECISION(LDU,K), WHERE LDU.GE.N.  IF
C                                   JOBA.EQ.1 THEN K.EQ.N, IF JOBA.GE.2
C                                   THEN K.EQ.MIN(N,P).
C                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
C                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X
C                   IN THE SUBROUTINE CALL.
C
C         V         DOUBLE PRECISION(LDV,P), WHERE LDV.GE.P.
C                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N,
C                   THEN V MAY BE IDENTIFIED WITH X IN THE
C                   SUBROUTINE CALL.
C
C         INFO      INTEGER.
C                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
C                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
C                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
C                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
C                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
C                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX
C                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
C                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U)
C                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR
C                   VALUES OF X AND B ARE THE SAME.
C
C     LINPACK. THIS VERSION DATED 03/19/79 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     DSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     EXTERNAL DROT
C     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2,DROTG
C     FORTRAN DABS,DMAX1,MAX0,MIN0,MOD,DSQRT
C
C     INTERNAL VARIABLES
C
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,
     *        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      DOUBLE PRECISION DDOT,T,R   ! R is never referenced
      DOUBLE PRECISION B,C,CS,EL,EMM1,F,G,DNRM2,SCALE,SHIFT,SL,SM,SN,
     *                 SMM1,T1,TEST,ZTEST
      LOGICAL WANTU,WANTV
C
C
C     SET THE MAXIMUM NUMBER OF ITERATIONS.
C
      MAXIT = 30
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN0(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
C
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
C
      INFO = 0
      NCT = MIN0(N-1,P)
      NRT = MAX0(0,MIN0(P-2,N))
      LU = MAX0(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
C
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
C           PLACE THE L-TH DIAGONAL IN S(L).
C
            S(L) = DNRM2(N-L+1,X(L,L),1)
            IF (S(L) .EQ. 0.0D0) GO TO 10
               IF (X(L,L) .NE. 0.0D0) S(L) = DSIGN(S(L),X(L,L))
               CALL DSCAL(N-L+1,1.0D0/S(L),X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .EQ. 0.0D0) GO TO 30
C
C              APPLY THE TRANSFORMATION.
C
               T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
            E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
C
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
C           MULTIPLICATION.
C
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
C
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
C           L-TH SUPER-DIAGONAL IN E(L).
C
            E(L) = DNRM2(P-L,E(LP1),1)
            IF (E(L) .EQ. 0.0D0) GO TO 80
               IF (E(LP1) .NE. 0.0D0) E(L) = DSIGN(E(L),E(LP1))
               CALL DSCAL(P-L,1.0D0/E(L),E(LP1),1)
               E(LP1) = 1.0D0 + E(LP1)
   80       CONTINUE
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0D0) GO TO 120
C
C              APPLY THE TRANSFORMATION.
C
               DO 90 I = LP1, N
                  WORK(I) = 0.0D0
   90          CONTINUE
               DO 100 J = LP1, P
                  CALL DAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL DAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
C
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
C              BACK MULTIPLICATION.
C
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
C
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
C
      M = MIN0(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0D0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0D0
C
C     IF REQUIRED, GENERATE U.
C
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = 0.0D0
  180       CONTINUE
            U(J,J) = 1.0D0
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0D0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
               DO 210 J = LP1, NCU
                  T = -DDOT(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL DAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL DSCAL(N-L+1,-1.0D0,U(L,L),1)
               U(L,L) = 1.0D0 + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = 0.0D0
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = 0.0D0
  260          CONTINUE
               U(L,L) = 1.0D0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C     IF IT IS REQUIRED, GENERATE V.
C
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0D0) GO TO 320
               DO 310 J = LP1, P
                  T = -DDOT(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL DAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = 0.0D0
  330       CONTINUE
            V(L,L) = 1.0D0
  340    CONTINUE
  350 CONTINUE
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
      MM = M
      ITER = 0
  360 CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
C     ...EXIT
         IF (M .EQ. 0) GO TO 620
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
         IF (ITER .LT. MAXIT) GO TO 370
            INFO = M
C     ......EXIT
            GO TO 620
  370    CONTINUE
C
C        THIS SECTION OF THE PROGRAM INSPECTS FOR
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
C
C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
C
         DO 390 LL = 1, M
            L = M - LL
C        ...EXIT
            IF (L .EQ. 0) GO TO 400
            TEST = DABS(S(L)) + DABS(S(L+1))
            ZTEST = TEST + DABS(E(L))
            IF (ZTEST .NE. TEST) GO TO 380
               E(L) = 0.0D0
C        ......EXIT
               GO TO 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
         IF (L .NE. M - 1) GO TO 410
            KASE = 4
         GO TO 480
  410    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 430 LLS = LP1, MP1
               LS = M - LLS + LP1
C           ...EXIT
               IF (LS .EQ. L) GO TO 440
               TEST = 0.0D0
               IF (LS .NE. M) TEST = TEST + DABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + DABS(E(LS-1))
               ZTEST = TEST + DABS(S(LS))
               IF (ZTEST .NE. TEST) GO TO 420
                  S(LS) = 0.0D0
C           ......EXIT
                  GO TO 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
            IF (LS .NE. L) GO TO 450
               KASE = 3
            GO TO 470
  450       CONTINUE
            IF (LS .NE. M) GO TO 460
               KASE = 1
            GO TO 470
  460       CONTINUE
               KASE = 2
               L = LS
  470       CONTINUE
  480    CONTINUE
         L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C
         GO TO (490,520,540,570), KASE
C
C        DEFLATE NEGLIGIBLE S(M).
C
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0D0
            DO 510 KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               IF (K .EQ. L) GO TO 500
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
  500          CONTINUE
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
         GO TO 610
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0D0
            DO 530 K = L, M
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL DROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
         GO TO 610
C
C        PERFORM ONE QR STEP.
C
  540    CONTINUE
C
C           CALCULATE THE SHIFT.
C
            SCALE = DMAX1(DABS(S(M)),DABS(S(M-1)),DABS(E(M-1)),
     *                    DABS(S(L)),DABS(E(L)))
            SM = S(M)/SCALE
            SMM1 = S(M-1)/SCALE
            EMM1 = E(M-1)/SCALE
            SL = S(L)/SCALE
            EL = E(L)/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0
            C = (SM*EMM1)**2
            SHIFT = 0.0D0
            IF (B .EQ. 0.0D0 .AND. C .EQ. 0.0D0) GO TO 550
               SHIFT = DSQRT(B**2+C)
               IF (B .LT. 0.0D0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  550       CONTINUE
            F = (SL + SM)*(SL - SM) - SHIFT
            G = SL*EL
C
C           CHASE ZEROS.
C
            MM1 = M - 1
            DO 560 K = L, MM1
               CALL DROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL DROTG(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)
     *            CALL DROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
C
C        CONVERGENCE.
C
  570    CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE.
C
            IF (S(L) .GE. 0.0D0) GO TO 580
               S(L) = -S(L)
               IF (WANTV) CALL DSCAL(P,-1.0D0,V(1,L),1)
  580       CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
  590       IF (L .EQ. MM) GO TO 600
C           ...EXIT
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)
     *            CALL DSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)
     *            CALL DSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            ITER = 0
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
      END
      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
      END
      SUBROUTINE EIGEN(MAX, N, A, ER, EI, ISV, ILV, V, WK, IERR)
C 
C   PURPOSE:
C      Compute all the eigenvalues and eigenvectors of a real n x n
C      matrix A stored as a variable-dimensioned two-dimensional array.
C 
C   REFERENCES:
C      Parlett, B.N.; and Reinsch, C.: Balancing a Matrix for Calcula-
C        tion of Eigenvalues and Eigenvectors. Numer. Math., Bd 13,
C        Heft 4, 1969, pp. 293-304.
C      Martin, R.S.; and Wilkinson, J.H.: Similarity Reduction of a
C        General Matrix to Hessenberg Form. Numer. Math., Bd. 12, Heft
C        5, 1968, pp. 349-368.
C      Martin, R.S.; Peters, G.; and Wilkinson, J.H.: The QR Algorithm
C        for Real Hesenberg Matrices. Numer. Math., Bd. 14, Heft 3,
C        1970, pp. 219-231.
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Compu-
C        tation. Volume II - Linear Algebra. Springer-Verlag, 1971.
C 
C   Subroutines employed by EIGEN: BALANC, BALBAK, CDIV, DAMCON, ELMHES,
C      ELTRAN, HQR, HQR2
C   Subroutines employing EIGEN: ASYREG, BILIN, CNTREG, CSTAB, DSTAB,
C      TESTST, POLE
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(MAX,N),ER(N),EI(N),V(MAX,1),WK(N,*)
C 
      IF (((ISV .NE. 0) .AND. (ISV .NE. N))
     1   .OR.  (ILV .NE. 0))   GO TO 999
C 
C   *** PRELIMINARY REDUCTION ***
      CALL BALANC (MAX,N,A,LOW,IGH,WK)
      CALL ELMHES (MAX,N,LOW,IGH,A,WK(1,2))
      IV = ISV + ILV
C 
      IF (IV .EQ. 0) THEN
C        *** COMPUTE ONLY EIGENVALUES ***
         CALL HQR (MAX,N,LOW,IGH,A,ER,EI,IERR)
      ELSE
C        *** COMPUTE EIGENVALUES AND EIGENVECTORS ***
         CALL ELTRAN (MAX,N,LOW,IGH,A,WK(1,2),V)
         CALL HQR2 (MAX,N,LOW,IGH,A,ER,EI,V,IERR)
      END IF
C 
      IF (IERR .NE. 0)   GO TO 9999
C 
      IF (IV .NE. 0) CALL BALBAK (MAX,N,LOW,IGH,WK,N,V)
C 
      DO I = 1,N
         WK(I,1) = ER(I)
         WK(I,2) = EI(I)
         WK(I,3) = ER(I)**2 + EI(I)**2
      END DO
C 
C   ***** NORMALIZE EIGENVECTORS *****
      J = 0
      DO WHILE (J .LT. IV)
         J = J + 1
         M = 1
         IF (EI(J) .EQ. 0) THEN
C           *** NORMALIZE REAL EIGENVECTOR ***
            DO I = 2,N
               IF (DABS(V(I,J)) .GT. DABS(V(M,J)))  M = I
            END DO
            IF (V(M,J) .NE. 0.0D0)  THEN
               DO I = 1,N
                  WK(I,J+3) = V(I,J) / V(M,J)
               END DO
            END IF
         ELSE
C           *** NORMALIZE COMPLEX EIGENVECTOR ***
            JP1 = J + 1
            DO I = 2,N
               IF ((DABS(V(I,J)) + DABS(V(I,JP1))) .GT.
     1             (DABS(V(M,J)) + DABS(V(M,JP1))))  M = I
            END DO
            IF ((V(M,J) .NE. 0.0D0) .OR. (V(M,JP1) .NE. 0.0D0)) THEN
               DO I = 1,N
                  CALL CDIV (V(I,J),V(I,JP1),V(M,J),V(M,JP1),
     1                       WK(I,J+3),WK(I,J+4))
               END DO
            END IF
            J = JP1
         END IF
      END DO
C 
C   ***** ORDER EIGENVALUES AND EIGENVECTORS FOR OUTPUT *****
C 
      DO I = 1,N
         P = DAMCON(1)
         DO J = 1,N
            IF (WK(J,3) .LT. P) THEN
               K = J
               P = WK(J,3)
            END IF
         END DO
         ER(I) = WK(K,1)
         EI(I) = WK(K,2)
         IF (IV .NE. 0) THEN
C           *** MOVE EIGENVECTOR ***
            DO J = 1,N
               V(J,I) = WK(J,K+3)
            END DO
         END IF
         WK(K,3) = DAMCON(1)
      END DO
      RETURN
C   *** ERROR MESSAGES ***
  999 WRITE(6,1001) ISV,ILV,N
 1001 FORMAT (' ILLEGAL USE OF EIGEN: ISV = ',I4,', ILV = ',I4,/,
     1' ISV MUST EQUAL ',I4,' AND ILV MUST EQUAL 0')
 9999 RETURN
      END
      SUBROUTINE ELMBAK(NM,LOW,IGH,A,INT,M,Z)
C 
C   PURPOSE:
C      Form the eigenvectors of a real square matrix A by transforming
C      those of the corresponding upper Hessenberg matrix determined by
C      subroutine ELMHES.
C 
C   REFERENCES:
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Computa-
C        tion.  Volume II - Linear Algebra.  Springer-Verlag, 1971.
C 
C   Subroutines employed by ELMBAK: None
C   Subroutines employing ELMBAK: EIGEN
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL*8 A(NM,IGH),Z(NM,M)
      REAL*8 X
      INTEGER INT(IGH)
C 
C 
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     ********** FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- **********
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         MP1 = MP + 1
C 
         DO 110 I = MP1, IGH
            X = A(I,MP-1)
            IF (X .EQ. 0.0) GO TO 110
C 
            DO 100 J = 1, M
  100       Z(I,J) = Z(I,J) + X * Z(MP,J)
C 
  110    CONTINUE
C 
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C 
         DO 130 J = 1, M
            X = Z(I,J)
            Z(I,J) = Z(MP,J)
            Z(MP,J) = X
  130    CONTINUE
C 
  140 CONTINUE
C 
  200 RETURN
C     ********** LAST CARD OF ELMBAK **********
      END
      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
C 
C   PURPOSE:
C      Reduce a submatrix situated in rows and columns LOW through
C      IGH to upper Hesenberg form by stabilized elementary similarity
C      transformations.  The computational method follows that of Martin
C      and Wilkinson.
C 
C   REFERENCES:
C      Martin, R.S.; and Wilkinson, J.H.: Similarity Reduction of a
C        General Matrix to Hessenberg Form.  Numer. Math., Bd. 12, Heft
C        5, 1968, pp. 349-368.
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Computa-
C        tion.  Volume II - Linear Algebra.  Springer-Verlag, 1971.
C 
C   Subroutines employed by ELMHES: None
C   Subroutines employing ELMHES: EIGEN
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL*8 A(NM,N)
      REAL*8 X,Y
C     REAL*8 DABS
      INTEGER INT(IGH)
C 
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C 
      DO 180 M = KP1, LA
         MM1 = M - 1
         X = 0.0
         I = M
C 
         DO 100 J = M, IGH
            IF (DABS(A(J,MM1)) .LE. DABS(X)) GO TO 100
            X = A(J,MM1)
            I = J
  100    CONTINUE
C 
         INT(M) = I
         IF (I .EQ. M) GO TO 130
C    ********** INTERCHANGE ROWS AND COLUMNS OF A **********
         DO 110 J = MM1, N
            Y = A(I,J)
            A(I,J) = A(M,J)
            A(M,J) = Y
  110    CONTINUE
C 
         DO 120 J = 1, IGH
            Y = A(J,I)
            A(J,I) = A(J,M)
            A(J,M) = Y
  120    CONTINUE
C    ********** END INTERCHANGE **********
  130    IF (X .EQ. 0.0) GO TO 180
         MP1 = M + 1
C 
         DO 160 I = MP1, IGH
            Y = A(I,MM1)
            IF (Y .EQ. 0.0) GO TO 160
            Y = Y / X
            A(I,MM1) = Y
C 
            DO 140 J = M, N
  140       A(I,J) = A(I,J) - Y * A(M,J)
C 
            DO 150 J = 1, IGH
  150       A(J,M) = A(J,M) + Y * A(J,I)
C 
  160    CONTINUE
C 
  180 CONTINUE
C 
  200 RETURN
C    ********** LAST CARD OF ELMHES **********
      END
      SUBROUTINE ELTRAN(NM,N,LOW,IGH,A,INT,Z)
C 
C   PURPOSE:
C      Accumulates the stabilized elementary similarity transformations
C      used in the reduction of a real general matrix to upper Hessen-
C      berg form by ELMHES.
C 
C   References:
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Computa-
C        tion.  Volume II - Linear Algebra.  Springer-Verlag, 1971.
C 
C   Subroutines employed by ELTRAN: None
C   Subroutines employing ELTRAN: EIGEN
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL*8 A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
C 
C     ********** INITIALIZE Z TO IDENTITY MATRIX **********
C 
      DO 80 I = 1,N
C 
         DO 60 J = 1,N
   60    Z(I,J) = 0.0
C 
         Z(I,I) = 1.0
   80 CONTINUE
C 
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
C     ********** FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- **********
      DO 140 MM = 1,KL
         MP = IGH - MM
         MP1 = MP + 1
C 
         DO 100 I = MP1,IGH
  100    Z(I,MP) = A(I,MP-1)
C 
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C 
         DO 130 J = MP,IGH
            Z(MP,J) = Z(I,J)
            Z(I,J) = 0.0
  130    CONTINUE
C 
         Z(I,MP) = 1.0
  140 CONTINUE
C 
  200 RETURN
      END
      SUBROUTINE EQUATE(A,NA,B,NB)
C 
C   PURPOSE:
C      Store a matrix in an alternate computer location.
C 
C   Subroutines employed by EQUATE: LNCNT
C   Subroutines employing EQUATE: ASYFIL, ASYREG, BARSTW, BILIN, CNTREG,
C      CSTAB, CTROL, DISREG, DSTAB, EXMDFL, EXPINT, EXPSER, FACTOR,
C      IMMDFL, PREFIL, RICNWT, SAMPL, SUM, TESTST, TRNSIT, VARANC
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),NA(2),NB(2)
      NB(1) = NA(1)
      NB(2) = NA(2)
      L=NA(1)*NA(2)
      IF( NA(1) .LT. 1 .OR. L .LT. 1 ) GO TO 999
      DO 300 I=1,L
  300 B(I)=A(I)
 1000 RETURN
  999 CALL LNCNT (1)
      WRITE  (6,50)  NA
   50 FORMAT  (' DIMENSION ERROR IN EQUATE  NA=',2I6)
      RETURN
      END
      SUBROUTINE EXMDFL (A,NA,B,NB,H,NH,AM,NAM,HM,NHM,Q,NQ,R,NR,F,NF,P,
     1NP,HIDENT,HMIDEN,DISC,NEWT,STABLE,FNULL,ALPHA,IOP,DUMMY)
C 
C   PURPOSE:
C      Solve either the continuous or discrete time-invariant asymptotic
C      explicit (model-in-the-system) model-following problem.
C 
C   REFERENCES:
C      Anderson, Brian D.O.; and Moore, John B.: Linear Optimal Con-
C        trol.  Prentice-Hall, Inc., c.1971.
C      Armstrong, E.S.: Digital Explicit Model Following With Unstable
C        Model Dynamics.  AIAA Paper No. 74-888, AIAA Mechanics and
C        Control of Flight Conference, Aug. 1974.
C 
C   Subroutines employed by EXMDFL: ADD, ASYREG, BARSTW, DPOCO, DPOSL,
C      EQUATE, LNCNT, MULT, PRNT, SCALE, SUBT, SUM, TRANP
C   Subroutines employing EXMDFL: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),H(*),AM(*),HM(*),Q(*),R(*),F(*),P(*),DUMMY(*)
      DIMENSION NA(2),NB(2),NH(2),NAM(2),NHM(2),NQ(2),NR(2),NF(2),NP(2)
      DIMENSION IOP(*),IOPT(5),NDUM1(2),NDUM2(2),NDUM3(2)
      LOGICAL  HIDENT,HMIDEN,DISC,NEWT,STABLE,FNULL,SYM
      COMMON/TOL/EPSAM,EPSBM,IACM
C 
      IF( IOP(1) .EQ. 0 ) GO TO 300
      CALL LNCNT(6)
      IF( DISC ) WRITE(6,25)
      IF( .NOT. DISC ) WRITE(6,50)
   25 FORMAT(/' PROGRAM TO SOLVE ASYMPTOTIC DISCRETE EXPLICIT MODEL-FOLL
     1OWING PROBLEM'//' PLANT DYNAMICS'/)
   50 FORMAT(/' PROGRAM TO SOLVE ASYMPTOTIC CONTINUOUS EXPLICIT MODEL-FO
     1LLOWING PROBLEM'//' PLANT DYNAMICS'/)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      IF( HIDENT ) GO TO 75
      CALL PRNT(H,NH,' H  ',1)
      GO TO 100
   75 CONTINUE
      CALL LNCNT(3)
      WRITE(6,85)
   85 FORMAT(/' H IS AN IDENTITY MATRIX'/)
C 
  100 CONTINUE
      CALL LNCNT(4)
      WRITE(6,125)
  125 FORMAT(//' MODEL DYNAMICS'/)
      CALL PRNT(AM,NAM,' AM ',1)
      IF( HMIDEN ) GO TO 175
      CALL PRNT(HM,NHM,' HM ',1)
      GO TO 200
  175 CONTINUE
      CALL LNCNT(3)
      WRITE(6,185)
  185 FORMAT(/' HM IS AN IDENTITY MATRIX '/)
C 
  200 CONTINUE
      CALL LNCNT(4)
      WRITE(6,225)
  225 FORMAT(//' WEIGHTING MATRICES '/)
      CALL PRNT(Q,NQ,' Q  ',1)
      CALL PRNT(R,NR,' R  ',1)
C 
  300 CONTINUE
      IF( IOP(2) .EQ. 0 ) GO TO 400
      NF(1) = NB(2)
      NF(2) = NA(1)
      NP(1) = NA(1)
      NP(2) = NA(1)
      IOPT(1) = IOP(3)
      IOPT(2) = IOP(4)
      IOPT(3) = IOP(5)
      IOPT(4) =  0
      IOPT(5) =  0
      N1 = NA(1)*NA(2) + 1
      CALL EQUATE(Q,NQ,DUMMY,NDUM1)
      CALL ASYREG(A,NA,B,NB,H,NH,DUMMY,NDUM1,R,NR,F,NF,P,NP,HIDENT,DISC
     1,NEWT,STABLE,FNULL,ALPHA,IOPT,DUMMY(N1))
C 
  400 CONTINUE
      IF( IOP(1) .EQ. 0 ) GO TO 600
      CALL LNCNT(4)
      WRITE(6,425)
  425 FORMAT(//' CONTROL LAW U = -F( COL.(X,XM) ),  F = (F11,F12)'/)
      CALL LNCNT(3)
      WRITE(6,450)
  450 FORMAT(/' PART OF F MULTIPLYING  X '/)
      CALL PRNT(F,NF,' F11',1)
      IF( .NOT. DISC  .AND. IOP(2) .EQ. 0 ) GO  TO 600
      CALL PRNT(P,NP,' P11',1)
      IF(  IOP(2) .EQ. 0 ) GO TO 600
      CALL LNCNT(2)
      WRITE(6,475)
  475 FORMAT(/' EIGENVALUES OF P11')
      NDUM1(1) = NA(1)
      NDUM1(2) = 1
      CALL PRNT(DUMMY(N1),NDUM1,'    ',3)
      N1 = N1 + NDUM1(1)
      NDUM1(2) = NA(1)
      CALL LNCNT(2)
      WRITE(6,500)
  500 FORMAT(/' PLANT CLOSED-LOOP RESPONSE MATRIX A - BF11')
      CALL PRNT(DUMMY(N1),NDUM1,'    ',3)
      CALL LNCNT(2)
      WRITE(6,525)
  525 FORMAT(/' EIGENVALUES OF CLOSED-LOOP RESPONSE MATRIX')
      N1 = N1 + NDUM1(1)*NDUM1(2)
      NDUM1(2) = 2
      CALL PRNT(DUMMY(N1),NDUM1,'    ',3)
C 
  600 CONTINUE
      NF(1)= NB(2)
      NF(2)= NA(1)
      CALL MULT(B,NB,F,NF,DUMMY,NA)
      CALL SUBT(A,NA,DUMMY,NA,DUMMY,NA)
      IF(  IOP(1).EQ. 0  .OR.  IOP(2) .NE. 0 ) GO TO 700
      CALL LNCNT(2)
      WRITE(6,500)
      CALL PRNT(DUMMY,NA,'    ',3)
C 
  700 CONTINUE
      N1 =  NA(1)**2 +1
      CALL TRANP(DUMMY,NA,DUMMY(N1),NA)
      CALL EQUATE(DUMMY(N1),NA,DUMMY,NA)
      NF(2) = NA(1) + NAM(1)
      NP(2) = NF(2)
      IF( .NOT. DISC .AND. IOP(2).EQ. 0 ) NP(2) = NAM(2)
      IOPTT = 0
      SYM = .FALSE.
      CALL EQUATE( Q,NQ,DUMMY(N1),NDUM2)
      IF( HMIDEN ) GO TO 725
      CALL MULT(Q,NQ,HM,NHM,DUMMY(N1),NDUM2)
  725 CONTINUE
      IF( HIDENT ) GO TO 750
      N2 = N1 + NQ(1)*NHM(2)
      CALL TRANP(H,NH,DUMMY(N2),NDUM1)
      N3 = N2 + NH(1)*NH(2)
      CALL MULT(DUMMY(N2),NDUM1,DUMMY(N1),NHM,DUMMY(N3),NDUM2)
      CALL EQUATE(DUMMY(N3),NDUM2,DUMMY(N1),NDUM2)
  750 CONTINUE
      N2 = NA(1)**2 + NA(1)*NHM(2) + 1
      N3 = NA(1)**2 + 1
      IF( .NOT. DISC .AND. IOP(2) .EQ. 0 ) N3 = 1
      CALL EQUATE(DUMMY(N1),NDUM2,P(N3),NDUM2)
      IF( DISC ) GO TO 800
      EPSA = EPSAM
      CALL BARSTW(DUMMY,NA,AM,NAM,P(N3),NDUM2,IOPTT,SYM ,EPSA,EPSA,DUMMY
     1(N2))
      GO TO 900
C 
  800 CONTINUE
      CALL SCALE(P(N3),NDUM2,P(N3),NDUM2,-1.0D0)
      N4 = N2 +NAM(1)**2
      CALL EQUATE(AM,NAM,DUMMY(N2),NAM)
      CALL SUM(DUMMY,NA,P(N3),NDUM2,DUMMY(N2),NAM,IOPTT,SYM,DUMMY(N4))
C 
  900 CONTINUE
      N2 = NB(2)*NA(1) + 1
      CALL TRANP(B,NB,DUMMY,NDUM1)
      CALL MULT(DUMMY,NDUM1,P(N3),NDUM2,F(N2),NDUM3)
      IF( .NOT. DISC ) GO TO 1000
      N1 = NB(1)*NB(2) + 1
      CALL MULT(DUMMY,NDUM1,P,NA,DUMMY(N1),NDUM2)
      CALL MULT(DUMMY(N1),NDUM2,B,NB,DUMMY,NR)
      CALL ADD(R,NR,DUMMY,NR,DUMMY,NR)
      GO TO 1100
C 
 1000 CONTINUE
      CALL EQUATE(R,NR,DUMMY,NR)
C 
 1100 CONTINUE
      N1 = NR(1)**2 + 1
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DPOCO(DUMMY,NR(1),NR(1),RCOND,DUMMY(N1),IERR)
      IF( IERR .EQ. 0 ) GO TO 1200
      CALL LNCNT(3)
      WRITE(6,1150)
 1150 FORMAT(/' IN EXMDFL, THE COEFFICIENT MATRIX FOR DPOCO IS NOT POSIT
     1IVE DEFINITE'/)
      RETURN
C 
 1200 CONTINUE
      NT = N2
      DO 1250 M1 = 1,NHM(2)
         CALL DPOSL(DUMMY,NR(1),NR(1),F(NT))
         NT = NT + NR(1)
 1250 CONTINUE
      IF( .NOT. DISC ) GO TO 1300
      CALL MULT(F(N2),NDUM3,AM,NAM,DUMMY,NDUM1)
      CALL EQUATE(DUMMY,NDUM1,F(N2),NDUM1)
C 
 1300 CONTINUE
      IF( IOP(1) .EQ. 0 ) RETURN
      CALL LNCNT(3)
      WRITE(6,1325)
 1325 FORMAT(/' PART OF F MULTIPLYING XM '/)
      CALL PRNT(F(N2),NDUM3,' F12',1)
      NDUM1(1) = NA(1)
      NDUM1(2) = NAM(1)
      CALL PRNT(P(N3),NDUM1,' P12',1)
      RETURN
      END
      SUBROUTINE EXPADE (MAX, N, A, EA, IDIG, WK, IERR)
C 
C   PURPOSE:
C      Compute the matrix exponential e**A, where A is a real square
C      matrix stored as a variable-dimensioned two-dimensional array.
C      Computation is by the Pade approximation method.  The exponen-
C      tial is computed using the approximation given by the ninth di-
C      agonal term in the Pade table for exponential approximations.
C 
C   REFERENCES:
C      Ward, R.C.: Numerical Computation of the Matrix Exponential with
C        Accuracy Estimate.  UCCND-CSD-24,Nov. 1975.
C 
C   Subroutines employed by EXPADE: GAUSEL
C   Subroutines employing EXPADE: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(MAX,N),EA(MAX,N),WK(N,1),C(9)
      IERR = 0
C ****
C     CALCULATE NORM OF A
C ****
      ANORM = 0.
      DO 10 I=1,N
         S = 0.
         DO 5 J=1,N
            S = S + DABS(A(I,J))
 5       CONTINUE
         IF (S .GT. ANORM)  ANORM = S
 10   CONTINUE
C ****
C     CALCULATE ACCURACY ESTIMATE
C ****
      DIGC = 24.*DFLOAT(N)
      IF (ANORM .GT. 1.) DIGC = DIGC*ANORM
      IDIG = 15 - IDNINT(DLOG10(DIGC))
C ****
C     DETERMINE POWER OF TWO AND NORMALIZATION FACTOR
C ****
      M = 0
      IF (ANORM .LE. 1.)  GO TO 27
      FACTOR =2.
      DO 15 M=1,46
         IF (ANORM .LE. FACTOR) GO TO 20
         FACTOR = FACTOR*2.
 15   CONTINUE
      GO TO 125
 20   CONTINUE
C ****
C     NORMALIZE MATRIX
C ****
      DO 25 I=1,N
         DO 25 J=1,N
            A(I,J) = A(I,J)/FACTOR
 25   CONTINUE
 27   CONTINUE
C ****
C     SET COEFFICIENTS FOR (9,9) PADE TABLE ENTRY
C ****
      C(1) = .5
      C(2) = 1.1764705882352D-01
      C(3) = 1.7156862745098D-02
      C(4) = 1.7156862745098D-03
      C(5) = 1.2254901960784D-04
      C(6) = 6.2845651080945D-06
      C(7) = 2.2444875386051D-07
      C(8) = 5.1011080422845D-09
      C(9) = 5.6678978247605D-11
C ****
C     CALCULATE PADE NUMERATOR AND DENOMINATOR BY COLUMNS
C ****
      NP1 = N+1
      NP7 = N+7
      DO 95 J=1,N
C ****
C        COMPUTE JTH COLUMN OF FIRST NINE POWERS OF A
C ****
         DO 35 I=1,N
            S = 0.
            DO 30 L=1,N
               S = S + A(I,L)*A(L,J)
 30         CONTINUE
            WK(I,NP1) = S
 35      CONTINUE
         DO 45 K=NP1,NP7
            KP1 = K+1
            DO 45 I=1,N
               S = 0.
               DO 40 L=1,N
                  S = S + A(I,L)*WK(L,K)
 40            CONTINUE
               WK(I,KP1) = S
 45      CONTINUE
C ****
C        COLLECT TERMS FOR JTH COLUMN OF NUMERATOR AND DENOMINATOR
C ****
         DO 85 I=1,N
            S = 0.
            U = 0.
            DO 65 L=1,8
               K = N+9-L
               KN1 = K-N+1
               P = C(KN1)*WK(I,K)
               S = S + P
              IEO = MOD(KN1,2)
              IF (IEO.EQ.0) GO TO 55
               U = U - P
               GO TO 65
 55            CONTINUE
               U = U + P
 65         CONTINUE
            P = C(1)*A(I,J)
            S = S + P
            U = U - P
            IF (I .NE. J) GO TO 80
            S = S + 1.
            U = U + 1.
 80         CONTINUE
            EA(I,J) = S
            WK(I,J) = U
 85      CONTINUE
 95   CONTINUE
C ****
C     CALCULATE NORMALIZED EXP(A) BY  WK * EXP(A) = EA
C ****
      CALL GAUSEL (MAX,N,WK,N,EA,IERR)
      IF (IERR .NE. 0) GO TO 130
      IF (M .EQ. 0)  GO TO 130
C ****
C     TAKE OUT EFFECT OF NORMALIZATION ON EXP(A)
C ****
      DO 120 K=1,M
         DO 110 I=1,N
            DO 110 J=1,N
               S = 0.
               DO 105 L=1,N
                  S = S + EA(I,L)*EA(L,J)
 105           CONTINUE
               WK(I,J) = S
 110     CONTINUE
         DO 115 I=1,N
            DO 115 J=1,N
               EA(I,J) = WK(I,J)
 115     CONTINUE
 120  CONTINUE
C ****
C     UN-NORMALIZE A
C ****
      DO 122 I=1,N
         DO 122 J=1,N
            A(I,J) = A(I,J)*FACTOR
 122  CONTINUE
      GO TO 130
C ****
C     NORM OF A IS EXCESSIVE
C ****
 125  CONTINUE
      IERR = 1
C ****
C     EXIT ROUTINE
C ****
 130  CONTINUE
      RETURN
      END
      SUBROUTINE EXPINT(A,NA,B,NB,C,NC,T,IOP,DUMMY)
C 
C   PURPOSE:
C      Compute both the matrix exponential e**AT and the integral
C      (0 to T) e**As ds, for a square real matrix A and scalar T.
C      Computation is based on the finite-series algorithm described
C      by Kallstom.
C 
C   REFERENCES:
C      Kallstrom, Claes: Computing Exp(A) and integral(Exp(A)ds).  Rep.
C        7309, Lund Inst. Technol. (Sweden), Mar. 1973.
C 
C   Subroutines employed by EXPINT: ADD, EQUATE, LNCNT, MAXEL, MULT,
C      NORMS, PRNT, SCALE, UNITY
C   Subroutines employing EXPINT: SAMPL, TRNSIT
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NC(2)
C      COMMON/CONV/SUMCV,MAXSUM,RICTCV,SERCV
      COMMON/CONV/SUMCV,RICTCV,SERCV,MAXSUM
C 
      N = NA(1)
      L = (N**2)+1
      NC(1) = NA(1)
      NC(2) = NA(2)
      NB(1) = NA(1)
      NB(2) = NA(2)
      TT = T
C 
      IOPT = 1
      CALL NORMS(N,N,N,A,IOPT,COL)
      IOPT = 3
      CALL NORMS(N,N,N,A,IOPT,ROW)
      ANAA = COL
      IF( ANAA .GT. ROW )  ANAA = ROW
      TMAX = 1./ANAA
      K = 0
  100 CONTINUE
      IF( TMAX - TT ) 125,150,150
  125 CONTINUE
      K = K + 1
      TT = T/(2**K)
      IF( K - 1000 )100,600,600
C 
  150 CONTINUE
      SC = TT
      CALL SCALE(A,NA,A,NA,TT)
      CALL UNITY(B,NB)
      CALL SCALE(B,NB,DUMMY,NB,TT)
      S =  TT/2.
      CALL SCALE(A,NA,DUMMY(L),NA,S)
      II = 2
      CALL ADD(DUMMY,NA,DUMMY(L),NA,DUMMY(L),NA)
      CALL ADD(A,NA,B,NB,DUMMY,NA)
      CALL EQUATE(A,NA,C,NC)
  200 CONTINUE
      CALL MULT(A,NA,C,NC,B,NB)
      S = 1./II
      CALL SCALE(B,NB,C,NC,S)
      CALL MAXEL(DUMMY,NA,TOT)
      CALL MAXEL(C,NC,DELT)
      IF( TOT .GT. 1.0 ) GO TO 300
      IF( DELT/TOT .LT. SERCV )  GO TO 400
      GO TO 350
  300 CONTINUE
      IF( DELT .LT. SERCV )  GO TO 400
  350 CONTINUE
      S = TT/(II + 1)
      CALL SCALE(C,NC,B,NB,S)
      CALL ADD(B,NB,DUMMY(L),NB,DUMMY(L),NB)
      CALL ADD(C,NC,DUMMY,NC,DUMMY,NC)
      II = II + 1
      GO TO 200
C 
  400 CONTINUE
      CALL EQUATE(DUMMY,NB,B,NB)
      IF( K ) 425,500,450
  425 CONTINUE
      CALL LNCNT(1)
      WRITE(6,435)
  435 FORMAT('  ERROR IN EXPINT, K IS NEGATIVE')
      RETURN
C 
  450 CONTINUE
      DO 475 J = 1,K
      TT = 2*TT
      CALL EQUATE(B,NB,DUMMY,NB)
      CALL MULT(DUMMY,NA,DUMMY(L),NA,C,NC)
      CALL ADD(DUMMY(L),NC,C,NC,DUMMY(L),NC)
      CALL MULT(DUMMY,NB,DUMMY,NB,B,NB)
  475 CONTINUE
      T = TT
C 
  500 CONTINUE
      CALL EQUATE(DUMMY(L),NC,C,NC)
      S = 1./SC
      CALL SCALE(A,NA,A,NA,S)
C 
      IF( IOP .EQ. 0 ) RETURN
      CALL LNCNT(5)
      WRITE(6,550)
  550 FORMAT(//' COMPUTATION OF THE MATRIX EXPONENTIAL  EXP(A T)'/' AND
     1ITS INTEGRAL OVER  (0,T) BY THE SERIES METHOD '/)
      CALL PRNT(A,NA,' A  ',1)
      CALL LNCNT(3)
      WRITE(6,575) T
  575 FORMAT(/'  T = ' ,D16.8/)
      CALL  PRNT(B,NB,'EXPA',1)
      CALL PRNT(C,NC,'INT ',1)
      RETURN
C 
  600 CONTINUE
      CALL LNCNT(1)
      WRITE(6,650)
  650 FORMAT( ' ERROR IN EXPINT, K = 1000 ')
      RETURN
C 
      END
      SUBROUTINE EXPSER(A,NA,EXPA,NEXPA,T,IOP,DUMMY)
C 
C   PURPOSE:
C      Evaluate the matrix exponential e**AT for a real square matrix
C      A and scalar T.  Computation is based on the finite-series
C      algorithm described by Kallstrom.
C 
C   REFERENCES:
C      Kallstrom, Claes: Computing Exp(A) and integral(Exp(A)ds). Rep.
C        7309, Lund Inst. Technol. (Sweden), Mar. 1973.
C 
C   Subroutines employed by EXPSER: ADD, DAMCON, EQUATE, LNCNT, MAXEL,
C      MULT, NORMS, PRNT, SCALE, TRCE, UNITY
C   Subroutines employing EXPSER: CNTREG, SAMPL, TRNSIT
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),EXPA(1),DUMMY(1)
      DIMENSION NA(2),NEXPA(2)
C      COMMON/CONV/SUMCV,MAXSUM,RICTCV,SERCV
      COMMON/CONV/SUMCV,RICTCV,SERCV,MAXSUM
!!!      REAL*8:: ZER0 = 0.0    ! something weird going on here
C 
      N = NA(1)
      L = (N**2) + 1
      TT = T
      NEXPA(1)=NA(1)
      NEXPA(2)=NA(2)
C 
      CALL MAXEL(A,NA,ANAA)
      ANAA = ANAA*TT
      IF ( DABS(ANAA) .GT. DAMCON(3) ) GO TO 100
      CALL UNITY(EXPA,NEXPA)
      GO TO 800
C 
  100 CONTINUE
      IOPT=2
      CALL NORMS(N,N,N,A,IOPT,ZERO)
      ZERO = ZERO * DAMCON(3)   ! changed by RLC  guessing that ZER0 was a bug
      CALL TRCE(A,NA,TR)
      TR = TR/N
      DO 200 I =1,N
      M =I+N*(I-1)
      A(M) = A(M) - TR
  200 CONTINUE
C 
      IOPT = 1
      CALL NORMS(N,N,N,A,IOPT,COL)
      IOPT = 3
      CALL NORMS(N,N,N,A,IOPT,ROW)
      ANORM = ROW
      IF( ANORM .GT. COL )  ANORM = COL
      TMAX = 1./ANORM
      K= 0
  300 CONTINUE
      IF( TMAX - TT ) 325,350,350
  325 CONTINUE
      K=K+1
      TT = T/(2**K)
      IF( K - 1000 ) 300,700,700
  350 CONTINUE
      SC = TT
      CALL SCALE(A,NA,A,NA,TT)
      CALL UNITY(EXPA,NEXPA)
      II = 2
      CALL ADD(A,NA,EXPA,NEXPA,DUMMY,NA)
      CALL EQUATE(A,NA,DUMMY(L),NA)
  400 CONTINUE
      CALL MULT(A,NA,DUMMY(L),NA,EXPA,NEXPA)
      S = 1./II
      CALL SCALE(EXPA,NEXPA,DUMMY(L),NA,S)
      CALL ADD(DUMMY(L),NA,DUMMY,NA,EXPA,NEXPA)
      CALL MAXEL(DUMMY,NA,TOT)
      CALL MAXEL(DUMMY(L),NA,DELT)
      IF( TOT .GT. 1.0 ) GO TO 500
      IF( DELT/TOT .LT. SERCV )  GO TO 600
      GO TO 550
  500 CONTINUE
      IF( DELT .LT. SERCV )  GO TO 600
  550 CONTINUE
      CALL EQUATE(EXPA,NEXPA,DUMMY,NA)
      II = II + 1
      GO TO 400
C 
  600 CONTINUE
      IF( K ) 625,675,650
  625 CONTINUE
      CALL LNCNT(1)
      WRITE(6,635)
  635 FORMAT( '   ERROR IN EXPSER,  K IS NEGATIVE ' )
      RETURN
C 
  650 CONTINUE
      DO 660 I =1,K
      TT = 2*TT
      CALL EQUATE(EXPA,NEXPA,DUMMY,NA)
      CALL EQUATE(DUMMY,NA,DUMMY(L),NA)
      CALL MULT(DUMMY(L),NA,DUMMY,NA,EXPA,NEXPA)
  660 CONTINUE
      T = TT
  675 CONTINUE
      S = 1./SC
      CALL SCALE(A,NA,A,NA,S)
      DO 685 I = 1,N
      M = I + N*(I-1)
      A(M) = A(M) + TR
      IF( DABS(A(M)) .LE. ZERO )  A(M) = 0.0
  685 CONTINUE
C 
      TR=TR*T
      S =  DEXP(TR)
      CALL SCALE(EXPA,NEXPA,EXPA,NEXPA,S)
      GO TO 800
C 
  700 CONTINUE
      CALL LNCNT(1)
      WRITE(6,750)
  750 FORMAT('  ERROR IN EXPSER,  K = 1000 ')
      RETURN
C 
  800 CONTINUE
      IF( IOP .EQ. 0 ) RETURN
      CALL LNCNT(4)
      WRITE(6,825)
  825 FORMAT(//' COMPUTATION OF THE MATRIX EXPONENTIAL EXP(A T) BY THE S
     1ERIES METHOD '/)
      CALL PRNT(A,NA,' A  ',1)
      CALL LNCNT(3)
      WRITE(6,850) T
  850 FORMAT(/' T = ' ,D16.8/)
      CALL PRNT(EXPA,NEXPA,'EXPA',1)
      RETURN
      END
      SUBROUTINE FACTOR(Q,NQ,D,ND,IOP,IAC,DUMMY)
C 
C   PURPOSE:
C      Compute a real m x n (m <= n) matrix D of rank m such that a
C      real n x n nonnegative definite matrix Q can be factored as
C                Q = D'D
C 
C   Subroutines employed by FACTOR: EQUATE, LNCNT, MULT, PRNT, SCALE,
C      SNVDEC, TRANP
C   Subroutine employing FACTOR: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(1),D(1),DUMMY(1)
      DIMENSION NQ(2),ND(2),NDUM(2)
C 
      IOPT = 2
      N = NQ(1)
      M = N**2
      N1 = M + 1
      N2 = N1 + N
C 
      CALL EQUATE(Q,NQ,DUMMY,NQ)
      CALL SNVDEC(IOPT,N,N,N,N,DUMMY,NOS,B,IAC,ZTEST,DUMMY(N1),D,IRANK,A
     1PLUS,IERR)
      IF( IERR .EQ. 0 ) GO TO 200
      CALL LNCNT(5)
      IF( IERR .GT. 0 ) WRITE(6,100) IERR
      IF( IERR .EQ. -1) WRITE(6,150) ZTEST,IRANK
  100 FORMAT(BZ,//' IN FACTOR , SNVDEC HAS FAILED TO CONVERGE TO THE ',I
     14,' SINGULARVALUE AFTER 30 ITERATIONS')
  150 FORMAT(//' IN FACTOR, THE MATRIX Q SUBMITTED TO SNVDEC IS CLOSE TO
     1 A MATRIX OF LOWER RANK USING ZTEST = ',D16.8/' IF THE ACCURACY IS
     2 REDUCED  THE RANK MAY ALSO BE REDUCED'/' CURRENT RANK =',I4)
      NDUM(1)=N
      NDUM(2)=1
      IF(IERR .EQ. -1)  CALL PRNT(DUMMY(N1),NDUM,'SNVL',1)
      IF( IERR .GT. 0 ) RETURN
C 
  200 CONTINUE
      NDUM(1) = N
C 
      DO 250 J =1,N
      M1 = (J-1)*N + 1
      M2 = J*N
      DO 250 I =M1,M2
      K = N2+I-1
      L = N1+J-1
      IF( DUMMY(L) .EQ. 0.0) GO TO 300
      DUMMY(K) = DSQRT(DUMMY(L))*DUMMY(I)
  250 CONTINUE
      NDUM(2)=N
      GO TO 350
C 
  300 NDUM(2) = J - 1
  350 CONTINUE
      IF( DUMMY(N2) .LT. 0.0 ) CALL SCALE(DUMMY(N2),NDUM,DUMMY(N2),NDUM,
     1-1.0D0)
      CALL TRANP(DUMMY(N2),NDUM,D,ND)
C 
      IF( IOP .EQ. 0 ) RETURN
      CALL LNCNT(4)
      WRITE(6,400)
  400 FORMAT(//' FACTOR Q AS (D TRANSPOSE)XD '/)
      CALL PRNT(Q,NQ,' Q  ',1)
      CALL PRNT(D,ND,' D  ',1)
      CALL MULT(DUMMY(N2),NDUM,D,ND,DUMMY,NQ)
      CALL PRNT(DUMMY,NQ,'DTXD',1)
C 
      RETURN
      END
      SUBROUTINE GAUSEL (MAX, N, A, NR, B, IERR)
C 
C   PURPOSE:
C      Solve a set of linear equations, AX=B, by the method of Gaussian
C      elimination.  The constant matrices A and B are of dimension
C      n x n and n x r.  No information is returned on the pivotal
C      strategy or the value of the determinant of A.
C 
C   Subroutines employed by GAUSEL: None
C   Subroutines employing GAUSEL: EXPADE
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),B(MAX,NR)
      NM1 = N-1
      IF (NM1 .EQ. 0) GO TO 140
C 
C     FIND LARGEST REMAINING ELEMENT IN I-TH COLUMN FOR PIVOT
C 
      DO 100 I=1,NM1
         BIG = 0.
         DO 20 K=I,N
            TERM = DABS(A(K,I))
            IF (TERM - BIG) 20,20,10
  10        BIG = TERM
            L = K
  20     CONTINUE
         IF (BIG) 40,30,40
  40     IF (I-L) 50,80,50
C 
C     PIVOT ROWS OF A AND B
C 
  50     CONTINUE
         DO 60 J=1,N
            TEMP = A(I,J)
            A(I,J) = A(L,J)
            A(L,J) = TEMP
  60     CONTINUE
         DO 70 J=1,NR
            TEMP = B(I,J)
            B(I,J) = B(L,J)
            B(L,J) = TEMP
  70     CONTINUE
  80     CONTINUE
C 
C     STORE PIVOT AND PERFORM COLUMN OPERATIONS ON A AND B
C 
         IP1 = I+1
         DO 100 II=IP1,N
            A(II,I) = A(II,I)/A(I,I)
            X3 = A(II,I)
            DO 90 K=IP1,N
               A(II,K) = A(II,K) - X3*A(I,K)
  90        CONTINUE
            DO 100 K=1,NR
               B(II,K) = B(II,K) - X3*B(I,K)
 100  CONTINUE
C 
C     PERFORM BACK SUBSTITUTION
C 
      DO 110 IC=1,NR
         B(N,IC) = B(N,IC)/A(N,N)
 110  CONTINUE
      DO 130 KK=1,NM1
         I = N-KK
         IP1 = I+1
         DO 130 J=1,NR
            SUM = B(I,J)
            DO 120 K=IP1,N
               SUM = SUM - A(I,K)*B(K,J)
 120        CONTINUE
            B(I,J) = SUM/A(I,I)
 130  CONTINUE
      RETURN
 140  CONTINUE
      IF (A(1,1) .EQ. 0.) GO TO 30
      DO 150 J=1,NR
         B(1,J) = B(1,J)/A(1,1)
 150  CONTINUE
      RETURN
  30     IERR = 2
      RETURN
      END
      SUBROUTINE GELIM(NMAX,N,A,NRHS,B,IPIVOT,IFAC,WK,IERR)
C 
C   PURPOSE:
C      Solve the real matrix equation, AX=B, where A is required to be
C      square and nonsingular and B is a matrix of constant vectors.
C      Solution is by Gaussian elimination or LU factorization.
C 
C   REFERENCES:
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Compu-
C        tation. Volume II - Linear Algebra. Springer-Verlag, 1971.
C 
C   Subroutines employed by GELIM: DETFAC
C   Subroutines employing GELIM: RICNWT
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NMAX,1),B(NMAX,1),IPIVOT(1),WK(1)
C 
      IERR=0
C 
C     TEST FOR L/U FACTORIZATION
C 
      IF(IFAC.EQ.1)GO TO 10
      CALL DETFAC(NMAX,N,A,IPIVOT,0,DETERM,ISCALE,WK,IERR)
      IF(IERR.GT.0)RETURN
   10 NM1=N-1
C 
C     TEST FOR SCALAR A MATRIX
C 
      IF(NM1.GT.0)GO TO 40
      IF(A(1,1).EQ.0.)GO TO 30
      DO 20 I=1,NRHS
   20 B(1,I)=B(1,I)/A(1,1)
      RETURN
   30 IERR=1
      RETURN
C 
   40 DO 100 M=1,NRHS
C 
C     PIVOT THE M-TH COLUMN OF B MATRIX
C 
      DO 50 I=1,NM1
      KI=IPIVOT(I)
      P=B(KI,M)
      B(KI,M)=B(I,M)
   50 B(I,M)=P
C 
C     FORWARD SUBSTITUTION
C 
      WK(1)=B(1,M)
C 
      DO 70 I=2,N
      IM1=I-1
      P=0.0
      DO 60 K=1,IM1
   60 P=P+A(I,K)*WK(K)
   70 WK(I)=B(I,M)-P
C 
C     BACK SUBSTITUTION
C 
      B(N,M)=WK(N)/A(N,N)
C 
      DO 90 J=1,NM1
      I=N-J
      IP1=I+1
      P=WK(I)
      DO 80 K=IP1,N
   80 P=P-A(I,K)*B(K,M)
   90 B(I,M)=P/A(I,I)
C 
  100 CONTINUE
      RETURN
      END
      SUBROUTINE HQR(NM,N,LOW,IGH,H,WR,WI,IERR)
C 
C   PURPOSE:
C      Find the eigenvalues of a real square upper Hessenberg matrix H
C      by the QR algorithm.  The computational method follows that of
C      Martin, Peters, and Wilkinson.
C 
C   REFERENCES:
C      Martin, R.S.; Peters, G.; and Wilkinson, J.H.: The QR Algorithm
C        for Real Hessenberg Matrices.  Numer. Math., Bd. 14, Heft 3,
C        1970, pp. 219-231.
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Computa-
C        tion.  Volume II - Linear Algebra.  Springer-Verlag, 1971.
C 
C   Subroutines employed by HQR: DAMCON
C   Subroutines employing HQR: EIGEN
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITS,LOW,MP2,ENM2,IERR
      REAL*8 H(NM,N),WR(N),WI(N)
      REAL*8 P,Q,R,S,T,W,X,Y,ZZ,NORM,MACHEP
C     REAL*8 DSQRT,DABS,DSIGN
C     INTEGER MIN0
      LOGICAL NOTLAS
C 
C 
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C 
C 
      MACHEP = DAMCON(3)
C 
      IERR = 0
      NORM = 0.0
      K = 1
C     ********** STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM **********
      DO 50 I = 1, N
C 
         DO 40 J = K, N
   40    NORM = NORM + DABS(H(I,J))
C 
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0
   50 CONTINUE
C 
      EN = IGH
      T = 0.0
C     ********** SEARCH FOR NEXT EIGENVALUES **********
   60 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- **********
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = DABS(H(L-1,L-1)) + DABS(H(L,L))
         IF (S .EQ. 0.0) S = NORM
         IF (DABS(H(L,L-1)) .LE. MACHEP * S) GO TO 100
   80 CONTINUE
C     ********** FORM SHIFT **********
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITS .EQ. 30) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     ********** FORM EXCEPTIONAL SHIFT **********
      T = T + X
C 
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C 
      S = DABS(H(EN,NA)) + DABS(H(NA,ENM2))
      X = 0.75 * S
      Y = X
      W = -0.4375 * S * S
  130 ITS = ITS + 1
C     ********** LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- **********
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = DABS(P) + DABS(Q) + DABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         IF (DABS(H(M,M-1)) * (DABS(Q) + DABS(R)) .LE. MACHEP * DABS(P)
     X    * (DABS(H(M-1,M-1)) + DABS(ZZ) + DABS(H(M+1,M+1)))) GO TO 150
  140 CONTINUE
C 
  150 MP2 = M + 2
C 
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0
  160 CONTINUE
C     ********** DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN **********
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0
         IF (NOTLAS) R = H(K+2,K-1)
         X = DABS(P) + DABS(Q) + DABS(R)
         IF (X .EQ. 0.0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = DSIGN(DSQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
C     ********** ROW MODIFICATION **********
         DO 210 J = K, EN
            P = H(K,J) + Q * H(K+1,J)
            IF (.NOT. NOTLAS) GO TO 200
            P = P + R * H(K+2,J)
            H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
            H(K,J) = H(K,J) - P * X
  210    CONTINUE
C 
         J = MIN0(EN,K+3)
C     ********** COLUMN MODIFICATION **********
         DO 230 I = L, J
            P = X * H(I,K) + Y * H(I,K+1)
            IF (.NOT. NOTLAS) GO TO 220
            P = P + ZZ * H(I,K+2)
            H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
            H(I,K) = H(I,K) - P
  230    CONTINUE
C 
  260 CONTINUE
C 
      GO TO 70
C     ********** ONE ROOT FOUND **********
  270 WR(EN) = X + T
      WI(EN) = 0.0
      EN = NA
      GO TO 60
C     ********** TWO ROOTS FOUND **********
  280 P = (Y - X) / 2.0
      Q = P * P + W
      ZZ = DSQRT(DABS(Q))
      X = X + T
      IF (Q .LT. 0.0) GO TO 320
C     ********** REAL PAIR **********
      ZZ = P + DSIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0
      WI(EN) = 0.0
      GO TO 330
C     ********** COMPLEX PAIR **********
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = EN
 1001 RETURN
C     ********** LAST CARD OF HQR **********
      END
      SUBROUTINE HQR2 (NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C 
C   PURPOSE:
C      Compute the eigenvalues and eigenvectors of a real upper Hes-
C      senberg matrix using the QR method.
C 
C   REFERENCES:
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Computa-
C        tion.  Volume II - Linear Algebra.  Springer-Verlag, 1971.
C 
C   Subroutines employed by HQR2: CDIV, DAMCON
C   Subroutines employing HQR2: EIGEN
C 
C  THIS SUBROUTINE IS A DOUBLE PRECISION FORM AND THE COMPLEX DIVIDES
C  HAVE BEEN REPLACED BY CALLS TO SUBROUTINE CDIV.
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 H(NM,N),T3(2),WR(N),WI(N),Z(NM,N)
      REAL*8 P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,IGH,ITS,LOW,MP2,
     1        ENM2,IERR
      LOGICAL NOTLAS
C 
C  ************ MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFICATION
C               THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C 
C     ***************************************************************
C 
      MACHEP= DAMCON(3)
C 
      IERR= 0
C  ****************** STORE ROOTS ISOLATED BY BALANC
      DO 50 I= 1,N
        IF (I.GE.LOW .AND. I.LE.IGH) GO TO 50
        WR(I)= H(I,I)
        WI(I)= 0.0
   50 CONTINUE
C 
      EN= IGH
      T= 0.0D0
C  ***************** SEARCH FOR NEXT EIGENVALUES
   60 IF (EN.LT.LOW) GO TO 340
      ITS= 0
      NA= EN-1
      ENM2= NA-1
C  ***************** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                    FOR L= EN STEP -1 UNTIL LOW DO ---  **************
   70 DO 80 LL= LOW,EN
        L= EN + LOW - LL
        IF (L.EQ.LOW) GO TO 100
        IF (DABS(H(L,L-1)).LE.MACHEP*(DABS(H(L-1,L-1)) +
     1                               DABS(H(L,L)))) GO TO 100
   80 CONTINUE
C  ***************** FORM SHIFT ******************************
  100 X= H(EN,EN)
      IF (L.EQ.EN) GO TO 270
      Y= H(NA,NA)
      W= H(EN,NA)*H(NA,EN)
      IF (L.EQ.NA) GO TO 280
      IF (ITS.EQ.30) GO TO 1000
      IF (ITS.NE.10 .AND. ITS.NE.20) GO TO 130
C  ***************** FORM EXCEPTIONAL SHIFT ******************
      T= T + X
C 
      DO 120 I= LOW,EN
        H(I,I)= H(I,I) - X
  120 CONTINUE
C 
      S= DABS(H(EN,NA)) + DABS(H(NA,ENM2))
      X= .75D0*S
      Y= X
      W= -.4375D0*S*S
  130 ITS= ITS+1
C  ******************* LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEME
C                      FOR M= EM-2 STEP -1 UNTIL L DO ---  *************
      DO 140 MM= L,ENM2
        M= ENM2 + L - MM
        ZZ= H(M,M)
        R= X - ZZ
        S= Y - ZZ
        P= (R*S - W)/H(M+1,M) + H(M,M+1)
        Q= H(M+1,M+1) - ZZ - R - S
        R= H(M+2,M+1)
        S= DABS(P) + DABS(Q) + DABS(R)
        P= P/S
        Q= Q/S
        R= R/S
        IF (M.EQ.L) GO TO 150
        IF (DABS(H(M,M-1))*(DABS(Q) + DABS(R)).LE.MACHEP*DABS(P)*
     1     (DABS(H(M-1,M-1)) + DABS(ZZ) + DABS(H(M+1,M+1)))) GO TO 150
  140 CONTINUE
C 
  150 MP2= M+2
C 
      DO 160 I= MP2,EN
        H(I,I-2)= 0.0D0
        IF (I.EQ.MP2) GO TO 160
          H(I,I-3)= 0.0D0
  160 CONTINUE
C  ****************** DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                     COLUMNS M TO EN
      DO 260 K= M,NA
        NOTLAS= K.NE.NA
        IF (K.EQ.M) GO TO 170
          P= H(K,K-1)
          Q= H(K+1,K-1)
          R= 0.0D0
          IF (NOTLAS) R= H(K+2,K-1)
          X= DABS(P) + DABS(Q) + DABS(R)
          IF (X.EQ.0.0D0) GO TO 260
          P= P/X
          Q= Q/X
          R= R/X
  170   S= DSIGN (DSQRT(P*P + Q*Q + R*R),P)
        IF (K.EQ.M) GO TO 180
          H(K,K-1)= -S*X
          GO TO 190
  180   IF (L.NE.M) H(K,K-1)= -H(K,K-1)
  190   P= P + S
        X= P/S
        Y= Q/S
        ZZ= R/S
        Q= Q/P
        R= R/P
C  ****************** ROW MODIFICATION ********************************
        DO 210 J= K,N
          P= H(K,J) + Q*H(K+1,J)
          IF (.NOT.NOTLAS) GO TO 200
            P= P + R*H(K+2,J)
            H(K+2,J)= H(K+2,J) - P*ZZ
  200     H(K+1,J)= H(K+1,J) - P*Y
          H(K,J)= H(K,J) - P*X
  210   CONTINUE
C 
        J= MIN0(EN,K+3)
C  ********************* COLUMN MODIFICATION *************************
        DO 230 I= 1,J
          P= X*H(I,K) + Y*H(I,K+1)
          IF (.NOT.NOTLAS) GO TO 220
            P= P + ZZ*H(I,K+2)
            H(I,K+2)= H(I,K+2) - P*R
  220     H(I,K+1)= H(I,K+1) - P*Q
          H(I,K)= H(I,K) - P
  230   CONTINUE
C  ******************* ACCUMULATE TRANSFORMATIONS ********************
        DO 250 I= LOW,IGH
          P= X*Z(I,K) + Y*Z(I,K+1)
          IF (.NOT.NOTLAS) GO TO 240
            P= P + ZZ*Z(I,K+2)
            Z(I,K+2)= Z(I,K+2) - P*R
  240     Z(I,K+1)= Z(I,K+1) - P*Q
          Z(I,K)= Z(I,K) - P
  250   CONTINUE
C 
  260 CONTINUE
C 
      GO TO 70
C  ****************** ONE ROOT FOUND **********************************
  270 H(EN,EN)= X + T
      WR(EN)= H(EN,EN)
      WI(EN)= 0.0D0
      EN= NA
      GO TO 60
C  ****************** TWO ROOTS FOUND *********************************
  280 P= (Y - X)/2.0D0
      Q= P*P + W
      ZZ= DSQRT(DABS(Q))
      H(EN,EN)= X + T
      X= H(EN,EN)
      H(NA,NA)= Y + T
      IF (Q.LT.0.0D0) GO TO 320
C  ******************** REAL PAIR **************************************
      ZZ= P + DSIGN(ZZ,P)
      WR(NA)= X + ZZ
      WR(EN)= WR(NA)
      IF (ZZ.NE.0.0D0) WR(EN)= X - W/ZZ
      WI(NA)= 0.0D0
      WI(EN)= 0.0D0
      X= H(EN,NA)
      R= DSQRT (X*X + ZZ*ZZ)
      P= X/R
      Q= ZZ/R
C  ******************* ROW MODIFICATION ********************************
      DO 290 J= NA,N
        ZZ= H(NA,J)
        H(NA,J)= Q*ZZ + P*H(EN,J)
        H(EN,J)= Q*H(EN,J) - P*ZZ
  290 CONTINUE
C  ******************* COLUMN MODIFICATION *****************************
      DO 300 I= 1,EN
        ZZ= H(I,NA)
        H(I,NA)= Q*ZZ + P*H(I,EN)
        H(I,EN)= Q*H(I,EN) - P*ZZ
  300 CONTINUE
C  ******************* ACCUMULATE TRANSFORMATIONS **********************
      DO 310 I= LOW,IGH
        ZZ= Z(I,NA)
        Z(I,NA)= Q*ZZ + P*Z(I,EN)
        Z(I,EN)= Q*Z(I,EN) - P*ZZ
  310 CONTINUE
C 
      GO TO 330
C  ********************** COMPLEX PAIR *********************************
  320 WR(NA)= X + P
      WR(EN)= X + P
      WI(NA)= ZZ
      WI(EN)= -ZZ
  330 EN= ENM2
      GO TO 60
C  ************************* ALL ROOTS FOUND. BACKSUBSTITUTE TO FIND
C                            VECTORS OF UPPER TRIANGULAR FORM  *********
  340 NORM= 0.0D0
      K= 1
C 
      DO 360 I= 1,N
C 
        DO 350 J= K,N
          NORM= NORM + DABS(H(I,J))
  350   CONTINUE
C 
        K= I
  360 CONTINUE
C 
      IF (NORM.EQ.0.0D0) GO TO 1001
C  **************** FOR EN= N STEP -1 UNTIL 1 DO --- *******************
      DO 800 NN= 1,N
        EN= N+1-NN
        P= WR(EN)
        Q= WI(EN)
        NA= EN-1
        IF (Q) 710,600,800
C  **************************** REAL VECTOR ****************************
  600   M= EN
        H(EN,EN)= 1.0D0
        IF (NA.EQ.0) GO TO 800
C  ******************* FOR I= EN-1 STEP -1 UNTIL 1 DO --- **************
        DO 700 II= 1,NA
          I= EN-II
          W= H(I,I) - P
          R= H(I,EN)
          IF (M.GT.NA) GO TO 620
C 
            DO 610 J= M,NA
  610         R= R + H(I,J)*H(J,EN)
C 
  620     IF (WI(I).GE.0.0D0) GO TO 630
            ZZ= W
            S= R
            GO TO 700
  630     M= I
          IF (WI(I).NE.0.0D0) GO TO 640
            T= W
            IF (W.EQ.0.0D0) T= MACHEP*NORM
            H(I,EN)= -R/T
            GO TO 700
C  ********************** SOLVE REAL EQUATIONS *************************
  640     X= H(I,I+1)
          Y= H(I+1,I)
          Q= (WR(I) - P)*(WR(I) - P) + WI(I)*WI(I)
          T= (X*S - ZZ*R)/Q
          H(I,EN)= T
          IF (DABS(X).LE.DABS(ZZ)) GO TO 650
            H(I+1,EN)= (-R - W*T)/X
            GO TO 700
  650     H(I+1,EN)= (-S - Y*T)/ZZ
  700   CONTINUE
C  *********************** END REAL VECTOR *****************************
        GO TO 800
C  *********************** COMPLEX VECTOR ******************************
  710   M= NA
C  *********************** LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THA
C                          EIGENVECTOR MATRIX IS TRIANGULAR **********
        IF (DABS(H(EN,NA)).LE.DABS(H(NA,EN))) GO TO 720
          H(NA,NA)= Q/H(EN,NA)
          H(NA,EN)= -(H(EN,EN) - P)/H(EN,NA)
          GO TO 730
  720   CALL CDIV (0.0,-H(NA,EN),H(NA,NA)-P,Q,T3(1),T3(2))
        H(NA,NA)= T3(1)
        H(NA,EN)= T3(2)
  730   H(EN,NA)= 0.0D0
        H(EN,EN)= 1.0D0
        ENM2= NA-1
        IF (ENM2.EQ.0) GO TO 800
C 
        DO 790 II= 1,ENM2
          I= NA-II
          W= H(I,I) - P
          RA= 0.0D0
          SA= H(I,EN)
C 
          DO 760 J= M,NA
            RA= RA + H(I,J)*H(J,NA)
            SA= SA + H(I,J)*H(J,EN)
  760     CONTINUE
C 
          IF (WI(I).GE.0.0D0) GO TO 770
            ZZ= W
            R= RA
            S= SA
            GO TO 790
  770     M= I
          IF (WI(I).NE.0.0D0) GO TO 780
            CALL CDIV (-RA,-SA,W,Q,T3(1),T3(2))
            H(I,NA)= T3(1)
            H(I,EN)= T3(2)
            GO TO 790
C  ************************** SOLVE COMPLEX EQUATIONS ******************
  780     X= H(I,I+1)
          Y= H(I+1,I)
          VR= (WR(I) - P)*(WR(I) - P) + WI(I)*WI(I) - Q*Q
          VI= (WR(I) - P)*2.0D0*Q
          IF (VR.EQ.0.0D0 .AND. VI.EQ.0.0D0)
     1     VR= MACHEP*NORM*(DABS(W)+DABS(Q)+DABS(X)+DABS(Y)+DABS(ZZ))
          CALL CDIV (X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI,T3(1),T3(2))
          H(I,NA)= T3(1)
          H(I,EN)= T3(2)
          IF (DABS(X).LE.DABS(ZZ) + DABS(Q)) GO TO 785
            H(I+1,NA)= (-RA - W*H(I,NA) + Q*H(I,EN))/X
            H(I+1,EN)= (-SA - W*H(I,EN) - Q*H(I,NA))/X
            GO TO 790
  785     CALL CDIV (-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q,T3(1),T3(2))
          H(I+1,NA)= T3(1)
          H(I+1,EN)= T3(2)
  790   CONTINUE
C  *********************** END COMPLEX VECTOR **************************
  800 CONTINUE
C  *********************** END BACK SUBSTITUTION
C                          GET VECTORS OF ISOLATED ROOTS ***************
      DO 840 I= 1,N
        IF (I.GE.LOW .AND. I.LE.IGH) GO TO 840
C 
        DO 820 J= I,N
  820     Z(I,J)= H(I,J)
C 
  840 CONTINUE
C  ********************* MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                        VECTORS OF ORIGINAL FULL MATRIX
C                        FOR J= N STEP -1 UNTIL LOW DO ---  ************
      DO 880 JJ= LOW,N
        J= N+LOW-JJ
        M= MIN0(J,IGH)
C 
        DO 880 I= LOW,IGH
          ZZ= 0.0D0
          DO 860 K= LOW,M
  860       ZZ= ZZ + Z(I,K)*H(K,J)
          Z(I,J)= ZZ
  880 CONTINUE
C 
      GO TO 1001
C  ********************* SET ERROR  --  NO CONVERGENCE TO AN
C                        EIGENVALUE AFTER 30 ITERATIONS ****************
 1000 IERR= EN
 1001 RETURN
      END
      SUBROUTINE HSHLDR(A,N,NA)
C 
C   PURPOSE:
C      Reduce a real n x n matrix A to upper Hessenberg form by House-
C      holder's method of elementary Hermitian transformations.
C   REFERENCES:
C      Wilkinson, J.H.: The Algebraic Eigenvalue Problem.  Clarendon
C        Press (Oxford), 1965.
C      Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
C        the Matrix Equation AX + XB = C.  Commun. ACM, vol. 15, no. 9,
C        Sept. 1972, pp. 820-826.
C 
C   Subroutines employed by HSHLDR: None
C   Subroutines employing HSHLDR: ATXPXA, AXPXB
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8
     1A(NA,1),MAX,SUM,S,P
      INTEGER
     1N,NA,NM2,N1,L,L1,I,J
      NM2 = N-2
      N1 = N+1
      IF(N .EQ. 1) RETURN
      IF(N .GT. 2) GO TO 5
      A(1,N1) = A(2,1)
      RETURN
    5 DO 80 L=1,NM2
        L1 = L+1
        MAX = 0.
        DO 10 I=L1,N
          MAX = DMAX1(MAX,DABS(A(I,L)))
   10   CONTINUE
        IF(MAX .NE. 0.) GO TO 20
        A(L,N1) = 0.
        A(N1,L) = 0.
        GO TO 80
   20   SUM = 0.
        DO 30 I=L1,N
          A(I,L) = A(I,L)/MAX
          SUM = SUM + A(I,L)**2
   30   CONTINUE
        S = DSIGN(DSQRT(SUM),A(L1,L))
        A(L,N1) = -MAX*S
        A(L1,L) = S + A(L1,L)
        A(N1,L) = S*A(L1,L)
        DO 50 J=L1,N
          SUM = 0.
          DO 40 I=L1,N
            SUM = SUM + A(I,L)*A(I,J)
   40     CONTINUE
          P = SUM/A(N1,L)
          DO 50 I=L1,N
            A(I,J) = A(I,J) - A(I,L)*P
   50   CONTINUE
        DO 70 I=1,N
          SUM = 0.
          DO 60 J=L1,N
            SUM = SUM + A(I,J)*A(J,L)
   60     CONTINUE
          P = SUM/A(N1,L)
          DO 70 J=L1,N
            A(I,J) = A(I,J) - P*A(J,L)
   70   CONTINUE
   80 CONTINUE
      A(N-1,N1) = A(N,N-1)
      RETURN
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N
C
      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
      SUBROUTINE IMMDFL(A,NA,B,NB,H,NH,AM,NAM,BM,NBM,Q,NQ,R,NR,F,NF,P,N
     1P,IDENT,DISC,NEWT,STABLE,FNULL,ALPHA,IOP,DUMMY)
C 
C   PURPOSE:
C      Solve either the continuous or discrete time-invariant asymptotic
C      implicit (model-in-the-performance-index) model-following
C      problem.
C 
C   Subroutines employed by IMMDFL: ADD, ASYREG, EQUATE, LNCNT, MULT,
C      PREFIL, PRNT, SCALE, SUBT, TRANP
C   Subroutines employing IMMDFL: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),H(*),AM(*),BM(*),Q(*),R(*),F(*),P(*),DUMMY(*)
      DIMENSION NA(2),NB(2),NH(2),NAM(2),NBM(2),NQ(2),NR(2),NF(2),NP(2)
      DIMENSION IOP(*),IOPT(5),NDUM1(2)
      LOGICAL IDENT,DISC,NEWT,STABLE,FNULL,HIDENT
C 
      IF( IOP(1) .EQ. 0 ) GO TO 200
      CALL LNCNT(6)
      IF( DISC ) WRITE(6,25)
      IF( .NOT. DISC ) WRITE(6,50)
   25 FORMAT(/' PROGRAM TO SOLVE ASYMPTOTIC DISCRETE IMPLICIT MODEL-FOLL
     1OWING PROBLEM'//' PLANT DYNAMICS '/)
   50 FORMAT(/' PROGRAM TO SOLVE ASYMPTOTIC CONTINUOUS IMPLICIT MODEL-FO
     1LLOWING PROBLEM'//' PLANT DYNAMICS'/)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      IF( IDENT ) GO TO 75
      CALL PRNT(H,NH,' H  ',1)
      GO TO 100
   75 CONTINUE
      CALL LNCNT(3)
      WRITE(6,85)
   85 FORMAT(/' H IS AN IDENTITY MATRIX'/)
C 
  100 CONTINUE
      CALL LNCNT(4)
      WRITE(6,125)
  125 FORMAT(//' MODEL DYNAMICS'/)
      CALL PRNT(AM,NAM,' AM ',1)
      CALL PRNT(BM,NBM,' BM ',1)
      CALL LNCNT(4)
      WRITE(6,150)
  150 FORMAT(//' WEIGHTING MATRICES'/)
      CALL PRNT(Q,NQ,' Q  ',1)
      CALL PRNT(R,NR,' R  ',1)
C 
  200 CONTINUE
      N = NA(1)**2
      N1 = N + 1
      IF( .NOT. IDENT ) GO TO 300
      CALL SUBT(A,NA,AM,NAM,DUMMY,NH)
      CALL SUBT(B,NB,BM,NBM,DUMMY(N1),NB)
      GO TO 400
C 
  300 CONTINUE
      CALL MULT(H,NH,A,NA,DUMMY,NH)
      CALL MULT(AM,NAM,H,NH,DUMMY(N1),NH)
      CALL SUBT(DUMMY,NH,DUMMY(N1),NH,DUMMY,NH)
      CALL MULT(H,NH,B,NB,DUMMY(N1),NBM)
      CALL SUBT(DUMMY(N1),NBM,BM,NBM,DUMMY(N1),NBM)
C 
  400 CONTINUE
      IF( IOP(1) .EQ. 0 ) GO TO 500
      CALL LNCNT(3)
      WRITE(6,450)
  450 FORMAT(//' MATRIX HA - AMH')
      CALL PRNT(DUMMY,NH,'    ',3)
      CALL LNCNT(3)
      WRITE(6,475)
  475 FORMAT(//' MATRIX HB - BM')
      CALL PRNT(DUMMY(N1),NBM,'    ',3)
C 
  500 CONTINUE
      N2 = N1 + N
      N3 = N2 + N
      N4 = N3 + N
      CALL MULT(Q,NQ,DUMMY,NH,DUMMY(N2),NH)
      CALL MULT(Q,NQ,DUMMY(N1),NBM,DUMMY(N3),NBM)
      CALL TRANP(DUMMY,NH,DUMMY(N4),NDUM1)
      CALL MULT(DUMMY(N4),NDUM1,DUMMY(N2),NH,DUMMY,NA)
      CALL MULT(DUMMY(N4),NDUM1,DUMMY(N3),NBM,DUMMY(N2),NB)
      CALL TRANP(DUMMY(N1),NBM,DUMMY(N4),NDUM1)
      CALL SCALE(DUMMY(N2),NB,DUMMY(N1),NB,2.0D0)
      CALL MULT(DUMMY(N4),NDUM1,DUMMY(N3),NBM,DUMMY(N2),NR)
      CALL ADD(DUMMY(N2),NR,R,NR,DUMMY(N2),NR)
      IF( IOP(1) .EQ. 0 ) GO TO 600
      CALL LNCNT(3)
      WRITE(6,525)
  525 FORMAT(//' MATRIX ( HA - AMH TRANSPOSE)Q( HA - AMH)')
      CALL PRNT(DUMMY,NA,'    ',3)
      CALL LNCNT(3)
      WRITE(6,550)
  550 FORMAT(//' MATRIX 2( HA - AMH  TRANSPOSE)Q( HB - BM)')
      CALL PRNT(DUMMY(N1),NB,'    ',3)
      CALL LNCNT(3)
      WRITE(6,575)
  575 FORMAT(//' MATRIX ( HB - BM  TRANSPOSE)Q( HB - BM ) + R')
      CALL PRNT(DUMMY(N2),NR,'    ',3)
C 
  600 CONTINUE
      IOPT(1)= 0
      IOPT(2)= 1
      IOPT(3)= 1
      N5 = N4 + N
      CALL EQUATE(A,NA,DUMMY(N3),NA)
      CALL PREFIL(DUMMY(N3),NA,B,NB,DUMMY,NA,DUMMY(N1),NB,DUMMY(N2),NR,D
     1UMMY(N4),NF,IOPT,DUMMY(N5))
      IF(IOP(1) .EQ. 0 ) GO TO 700
      CALL LNCNT(3)
      WRITE(6,625)
  625 FORMAT(//' PREFILTER GAIN')
      CALL PRNT(DUMMY(N4),NF,'    ',3)
      CALL LNCNT(3)
      WRITE(6,650)
  650 FORMAT(//' MATRIX A - B(PREFILTER)')
      CALL PRNT(DUMMY(N3),NA,'    ',3)
      CALL LNCNT(3)
      WRITE(6,675)
  675 FORMAT(//' MODIFIED STATE VECTOR WEIGHTING MATRIX')
      CALL PRNT(DUMMY,NA,'    ',3)
C 
  700 CONTINUE
      CALL EQUATE(DUMMY(N4),NF,DUMMY(N1),NF)
C 
      IF( IOP(2) .EQ. -1000 ) RETURN
C 
      NF(1)=NB(2)
      NF(2)=NA(1)
      NP(1)=NA(1)
      NP(2)=NA(1)
      IOPT(1) = IOP(2)
      IOPT(2) = IOP(3)
      IOPT(3) = IOP(4)
      IOPT(4) = 0
      IOPT(5) = 0
      HIDENT = .TRUE.
      CALL ASYREG(DUMMY(N3),NA,B,NB,H,NH,DUMMY,NA,DUMMY(N2),NR,F,NF,P,N
     1P,HIDENT,DISC,NEWT,STABLE,FNULL,ALPHA,IOPT,DUMMY(N5))
      IF( IOP(1) .EQ. 0 ) GO TO 800
      CALL LNCNT(3)
      WRITE(6,725)
  725 FORMAT(//' GAIN FROM ASYREG')
      CALL PRNT(F,NF,'    ',3)
      CALL LNCNT(3)
      WRITE(6,750)
  750 FORMAT(//' SOLUTION OF ASSOCIATED STEADY-STATE RICCATI EQUATION')
      CALL PRNT(P,NP,'    ',3)
      CALL LNCNT(3)
      WRITE(6,775)
  775 FORMAT(//' EIGENVALUES OF P')
      NDUM1(1)= NA(1)
      NDUM1(2)= 1
      CALL PRNT(DUMMY(N5),NDUM1,'    ',3)
C 
  800 CONTINUE
      CALL ADD(F,NF,DUMMY(N1),NF,F,NF)
      IF( IOP(1) .EQ. 0 ) RETURN
      CALL LNCNT(4)
      WRITE(6,825)
  825 FORMAT(//' GAIN FOR MODEL-FOLLOWING CONTROL LAW, U = - F X  , F =
     1(PREFILTER) + (ASYREG)'/)
      CALL PRNT(F,NF,' F  ',1)
      N6 = N5 + NA(1)
      CALL PRNT(DUMMY(N6),NA,'A-BF',1)
      NDUM1(2) = 2
      N6 = N6 + N
      CALL LNCNT(3)
      WRITE(6,850)
  850 FORMAT(//' EIGENVALUES OF A-BF')
      CALL PRNT(DUMMY(N6),NDUM1,'    ',3)
C 
      RETURN
      END
      SUBROUTINE INVIT(NM,N,A,WR,WI,SELECT,MM,M,Z,IERR,RM1,RV1,RV2)
C 
C   PURPOSE:
C      Find those eigenvectors of a real square upper Hessenberg matrix
C      corresponding to specified eigenvalues using inverse iteration.
C 
C   REFERENCES:
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Computa-
C        tion.  Volume II - Linear Algebra.  Springer-Verlag, 1971.
C 
C   Subroutines employed by INVIT: CDIV, DAMCON
C   Subroutines employing INVIT: EIGEN
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,K,L,M,N,S,II,IP,MM,MP,NM,NS,N1,UK,IP1,ITS,KM1,IERR
      REAL*8 A(NM,N),WR(N),WI(N),Z(NM,MM),RM1(N,N),RV1(N),RV2(N)
      REAL*8 T,W,X,Y,EPS3,NORM,NORMV,GROWTO,ILAMBD,MACHEP,RLAMBD,UKROOT
C     REAL*8 DSQRT,DABS,DFLOAT
C     INTEGER IABS
      LOGICAL SELECT(N)
C 
C 
      MACHEP = DAMCON(3)
C 
      IERR = 0
      UK = 0
      S = 1
C     ********** IP = 0, REAL EIGENVALUE
C                     1, FIRST OF CONJUGATE COMPLEX PAIR
C                    -1, SECOND OF CONJUGATE COMPLEX PAIR **********
      IP = 0
      N1 = N - 1
C 
      DO 980 K = 1, N
         IF (WI(K) .EQ. 0.0 .OR. IP .LT. 0) GO TO 100
         IP = 1
         IF (SELECT(K) .AND. SELECT(K+1)) SELECT(K+1) = .FALSE.
  100    IF (.NOT. SELECT(K)) GO TO 960
         IF (WI(K) .NE. 0.0) S = S + 1
         IF (S .GT. MM) GO TO 1000
         IF (UK .GE. K) GO TO 200
C     ********** CHECK FOR POSSIBLE SPLITTING **********
         DO 120 UK = K, N
            IF (UK .EQ. N) GO TO 140
            IF (A(UK+1,UK) .EQ. 0.0) GO TO 140
  120    CONTINUE
C     ********** COMPUTE INFINITY NORM OF LEADING UK BY UK
C                (HESSENBERG) MATRIX **********
  140    NORM = 0.0
         MP = 1
C 
         DO 180 I = 1, UK
            X = 0.0
C 
            DO 160 J = MP, UK
  160       X = X + DABS(A(I,J))
C 
            IF (X .GT. NORM) NORM = X
            MP = I
  180    CONTINUE
C     ********** EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
C                AND CLOSE ROOTS ARE MODIFIED BY EPS3 **********
         IF (NORM .EQ. 0.0) NORM = 1.0
         EPS3 = MACHEP * NORM
C     ********** GROWTO IS THE CRITERION FOR THE GROWTH **********
         UKROOT = DSQRT(DFLOAT(UK))
         GROWTO = 1.0E-1 / UKROOT
  200    RLAMBD = WR(K)
         ILAMBD = WI(K)
         IF (K .EQ. 1) GO TO 280
         KM1 = K - 1
         GO TO 240
C     ********** PERTURB EIGENVALUE IF IT IS CLOSE
C                TO ANY PREVIOUS EIGENVALUE **********
  220    RLAMBD = RLAMBD + EPS3
C     ********** FOR I=K-1 STEP -1 UNTIL 1 DO -- **********
  240    DO 260 II = 1, KM1
            I = K - II
            IF (SELECT(I) .AND. DABS(WR(I)-RLAMBD) .LT. EPS3 .AND.
     X         DABS(WI(I)-ILAMBD) .LT. EPS3) GO TO 220
  260    CONTINUE
C 
         WR(K) = RLAMBD
C     ********** PERTURB CONJUGATE EIGENVALUE TO MATCH **********
         IP1 = K + IP
         WR(IP1) = RLAMBD
C     ********** FORM UPPER HESSENBERG A-RLAMBD*I (TRANSPOSED)
C                AND INITIAL REAL VECTOR **********
  280    MP = 1
C 
         DO 320 I = 1, UK
C 
            DO 300 J = MP, UK
  300       RM1(J,I) = A(I,J)
C 
            RM1(I,I) = RM1(I,I) - RLAMBD
            MP = I
            RV1(I) = EPS3
  320    CONTINUE
C 
         ITS = 0
         IF (ILAMBD .NE. 0.0) GO TO 520
C     ********** REAL EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3 **********
         IF (UK .EQ. 1) GO TO 420
C 
         DO 400 I = 2, UK
            MP = I - 1
            IF (DABS(RM1(MP,I)) .LE. DABS(RM1(MP,MP))) GO TO 360
C 
            DO 340 J = MP, UK
               Y = RM1(J,I)
               RM1(J,I) = RM1(J,MP)
               RM1(J,MP) = Y
  340       CONTINUE
C 
  360       IF (RM1(MP,MP) .EQ. 0.0) RM1(MP,MP) = EPS3
            X = RM1(MP,I) / RM1(MP,MP)
            IF (X .EQ. 0.0) GO TO 400
C 
            DO 380 J = I, UK
  380       RM1(J,I) = RM1(J,I) - X * RM1(J,MP)
C 
  400    CONTINUE
C 
  420    IF (RM1(UK,UK) .EQ. 0.0) RM1(UK,UK) = EPS3
C     ********** BACK SUBSTITUTION FOR REAL VECTOR
C                FOR I=UK STEP -1 UNTIL 1 DO -- **********
  440    DO 500 II = 1, UK
            I = UK + 1 - II
            Y = RV1(I)
            IF (I .EQ. UK) GO TO 480
            IP1 = I + 1
C 
            DO 460 J = IP1, UK
  460       Y = Y - RM1(J,I) * RV1(J)
C 
  480       RV1(I) = Y / RM1(I,I)
  500    CONTINUE
C 
         GO TO 740
C     ********** COMPLEX EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3.  STORE IMAGINARY
C                PARTS IN UPPER TRIANGLE STARTING AT (1,3) **********
  520    NS = N - S
         Z(1,S-1) = -ILAMBD
         Z(1,S) = 0.0
         IF (N .EQ. 2) GO TO 550
         RM1(1,3) = -ILAMBD
         Z(1,S-1) = 0.0
         IF (N .EQ. 3) GO TO 550
C 
         DO 540 I = 4, N
  540    RM1(1,I) = 0.0
C 
  550    DO 640 I = 2, UK
            MP = I - 1
            W = RM1(MP,I)
            IF (I .LT. N) T = RM1(MP,I+1)
            IF (I .EQ. N) T = Z(MP,S-1)
            X = RM1(MP,MP) * RM1(MP,MP) + T * T
            IF (W * W .LE. X) GO TO 580
            X = RM1(MP,MP) / W
            Y = T / W
            RM1(MP,MP) = W
            IF (I .LT. N) RM1(MP,I+1) = 0.0
            IF (I .EQ. N) Z(MP,S-1) = 0.0
C 
            DO 560 J = I, UK
               W = RM1(J,I)
               RM1(J,I) = RM1(J,MP) - X * W
               RM1(J,MP) = W
               IF (J .LT. N1) GO TO 555
               L = J - NS
               Z(I,L) = Z(MP,L) - Y * W
               Z(MP,L) = 0.0
               GO TO 560
  555          RM1(I,J+2) = RM1(MP,J+2) - Y * W
               RM1(MP,J+2) = 0.0
  560       CONTINUE
C 
            RM1(I,I) = RM1(I,I) - Y * ILAMBD
            IF (I .LT. N1) GO TO 570
            L = I - NS
            Z(MP,L) = -ILAMBD
            Z(I,L) = Z(I,L) + X * ILAMBD
            GO TO 640
  570       RM1(MP,I+2) = -ILAMBD
            RM1(I,I+2) = RM1(I,I+2) + X * ILAMBD
            GO TO 640
  580       IF (X .NE. 0.0) GO TO 600
            RM1(MP,MP) = EPS3
            IF (I .LT. N) RM1(MP,I+1) = 0.0
            IF (I .EQ. N) Z(MP,S-1) = 0.0
            T = 0.0
            X = EPS3 * EPS3
  600       W = W / X
            X = RM1(MP,MP) * W
            Y = -T * W
C 
            DO 620 J = I, UK
               IF (J .LT. N1) GO TO 610
               L = J - NS
               T = Z(MP,L)
               Z(I,L) = -X * T - Y * RM1(J,MP)
               GO TO 615
  610          T = RM1(MP,J+2)
               RM1(I,J+2) = -X * T - Y * RM1(J,MP)
  615          RM1(J,I) = RM1(J,I) - X * RM1(J,MP) + Y * T
  620       CONTINUE
C 
            IF (I .LT. N1) GO TO 630
            L = I - NS
            Z(I,L) = Z(I,L) - ILAMBD
            GO TO 640
  630       RM1(I,I+2) = RM1(I,I+2) - ILAMBD
  640    CONTINUE
C 
         IF (UK .LT. N1) GO TO 650
         L = UK - NS
         T = Z(UK,L)
         GO TO 655
  650    T = RM1(UK,UK+2)
  655    IF (RM1(UK,UK) .EQ. 0.0 .AND. T .EQ. 0.0) RM1(UK,UK) = EPS3
C     ********** BACK SUBSTITUTION FOR COMPLEX VECTOR
C                FOR I=UK STEP -1 UNTIL 1 DO -- **********
  660    DO 720 II = 1, UK
            I = UK + 1 - II
            X = RV1(I)
            Y = 0.0
            IF (I .EQ. UK) GO TO 700
            IP1 = I + 1
C 
            DO 680 J = IP1, UK
               IF (J .LT. N1) GO TO 670
               L = J - NS
               T = Z(I,L)
               GO TO 675
  670          T = RM1(I,J+2)
  675          X = X - RM1(J,I) * RV1(J) + T * RV2(J)
               Y = Y - RM1(J,I) * RV2(J) - T * RV1(J)
  680       CONTINUE
C 
  700       IF (I .LT. N1) GO TO 710
            L = I - NS
            T = Z(I,L)
            GO TO 715
  710       T = RM1(I,I+2)
  715       CALL CDIV(X,Y,RM1(I,I),T,RV1(I),RV2(I))
  720    CONTINUE
C     ********** ACCEPTANCE TEST FOR REAL OR COMPLEX
C                EIGENVECTOR AND NORMALIZATION **********
  740    ITS = ITS + 1
         NORM = 0.0
         NORMV = 0.0
C 
         DO 780 I = 1, UK
            IF (ILAMBD .EQ. 0.0) X = DABS(RV1(I))
            IF (ILAMBD .NE. 0.0) X = DSQRT(RV1(I)**2 + RV2(I)**2)
            IF (NORMV .GE. X) GO TO 760
            NORMV = X
            J = I
  760       NORM = NORM + X
  780    CONTINUE
C 
         IF (NORM .LT. GROWTO) GO TO 840
C     ********** ACCEPT VECTOR **********
         X = RV1(J)
         IF (ILAMBD .EQ. 0.0) X = 1.0 / X
         IF (ILAMBD .NE. 0.0) Y = RV2(J)
C 
         DO 820 I = 1, UK
            IF (ILAMBD .NE. 0.0) GO TO 800
            Z(I,S) = RV1(I) * X
            GO TO 820
  800       CALL CDIV(RV1(I),RV2(I),X,Y,Z(I,S-1),Z(I,S))
  820    CONTINUE
C 
         IF (UK .EQ. N) GO TO 940
         J = UK + 1
         GO TO 900
C     ********** IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR **********
  840    IF (ITS .GE. UK) GO TO 880
         X = UKROOT
         Y = EPS3 / (X + 1.0)
         RV1(1) = EPS3
C 
         DO 860 I = 2, UK
  860    RV1(I) = Y
C 
         J = UK - ITS + 1
         RV1(J) = RV1(J) - EPS3 * X
         IF (ILAMBD .EQ. 0.0) GO TO 440
         GO TO 660
C     ********** SET ERROR -- UNACCEPTED EIGENVECTOR **********
  880    J = 1
         IERR = -K
C     ********** SET REMAINING VECTOR COMPONENTS TO ZERO **********
  900    DO 920 I = J, N
            Z(I,S) = 0.0
            IF (ILAMBD .NE. 0.0) Z(I,S-1) = 0.0
  920    CONTINUE
C 
  940    S = S + 1
  960    IF (IP .EQ. (-1)) IP = 0
         IF (IP .EQ. 1) IP = -1
  980 CONTINUE
C 
      GO TO 1001
C     ********** SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
C                SPACE REQUIRED **********
 1000 IF (IERR .NE. 0) IERR = IERR - N
      IF (IERR .EQ. 0) IERR = -(2 * N + 1)
 1001 M = S - 1 - IABS(IP)
      RETURN
C     ********** LAST CARD OF INVIT **********
      END
      SUBROUTINE JUXTC(A,NA,B,NB,C,NC)
C 
C   PURPOSE:
C      Construct a matrix [A,B] from given matrices A and B.
C 
C   Subroutines employed by JUXTC: LNCNT
C   Subroutines employing JUXTC: ASYREG, CTROL, DSTAB, TESTST
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),NA(2),NB(2),NC(2)
      IF (NA(1).NE.NB(1)) GO TO 600
      NC(1)=NA(1)
      NC(2)=NA(2)+NB(2)
      L=NA(1)*NA(2)
      NNC=NC(1)*NC(2)
      IF( NA(1) .LT. 1 .OR. L .LT. 1 ) GO TO 600
      IF( NC(2) .LT. 1  )  GO TO 600
      MS=NA(1)*NA(2)
      DO 10 I=1,MS
   10 C(I)=A(I)
      MBS=NA(1)*NB(2)
      DO 20 I=1,MBS
      J=MS+I
   20 C(J)=B(I)
      RETURN
  600 CALL LNCNT(1)
      WRITE (6,1600) NA,NB
 1600 FORMAT (' DIMENSION ERROR IN JUXTC,  NA=',2I6,5X,'NB=',2I6)
      RETURN
      END
      SUBROUTINE JUXTR(A,NA,B,NB,C,NC)
C 
C   PURPOSE:
C                            A
C      Construct a matrix  [ B ] from given matrices A and B.
C 
C   Subroutines employed by JUXTR: LNCNT
C   Subroutines employing JUXTR: BARSTW, CNTREG
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),NA(2),NB(2),NC(2)
      IF(NA(2).NE.NB(2))GO TO 600
      NC(2)=NA(2)
      NC(1)=NA(1)+NB(1)
      L=NA(1)*NA(2)
      IF( NA(1) .LT. 1 .OR. L .LT. 1 ) GO TO 600
      IF( NC(2) .LT. 1 ) GO TO 600
      MCA=NA(2)
      MRA=NA(1)
      MRB=NB(1)
      MRC=NC(1)
      DO 10 I=1,MCA
      DO 10 J=1,MRA
      K=J+MRA*(I-1)
      L=J+MRC*(I-1)
   10 C(L)=A(K)
      DO 20 I=1,MCA
      DO 20 J=1,MRB
      K=J+MRB*(I-1)
      L=MRA+J+MRC*(I-1)
   20 C(L)=B(K)
      RETURN
  600 CALL LNCNT(1)
      WRITE(6,1600) NA,NB
 1600 FORMAT(' DIMENSION ERROR IN JUXTR,  NA=',2I6,5X,'NB=',2I6)
      RETURN
      END
      SUBROUTINE LEVIER(A,NA,B,NB,H,NH,C,NC,HCB,NHCB,D,ND,BIDENT,HIDENT,
     1IOP,DUMMY)
C 
C   PURPOSE:
C                                           -1
C      Evaluate the transfer matrix H(sI - A)  B  for the linear
C                             .
C      time-invariant system  x(t) = Ax(t) + Bu(t) with output
C      y(t) = Hx(t).
C 
C   Subroutines employed by LEVIER: ADD, EQUATE, LNCNT, MULT, PRNT,
C      SCALE, TRCE, UNITY
C   Subroutines employing LEVIER: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),H(1),C(1),HCB(1),D(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NH(2),NC(2),NHCB(2),ND(2),NDUM1(2),NDUM2(2)
      LOGICAL  BIDENT,HIDENT
C 
      IF( IOP .EQ. 0 ) GO TO 100
      CALL LNCNT(6)
      WRITE(6,10)
   10 FORMAT(//' THE TRANSFER MATRIX  H((SI-A)INVERSE)B IS COMPUTED'/' F
     1OR THE  (A,B,H) SYSTEM'/)
      CALL PRNT(A,NA,' A  ',1)
      IF( BIDENT .AND. HIDENT ) GO TO 40
      IF( .NOT. BIDENT ) CALL PRNT(B,NB,' B  ',1)
      IF( .NOT. HIDENT ) CALL PRNT(H,NH,' H  ',1)
      IF( (.NOT. BIDENT) .AND. (.NOT. HIDENT) )  GO TO 50
      CALL LNCNT(3)
      IF( BIDENT ) WRITE(6,20)
      IF( HIDENT ) WRITE(6,30)
   20 FORMAT(/' B IS AN IDENTITY MATRIX'/)
   30 FORMAT(/' H IS AN IDENTITY MATRIX'/)
      GO TO 50
   40 CONTINUE
      CALL LNCNT(6)
      WRITE(6,20)
      WRITE(6,30)
   50 CONTINUE
      CALL LNCNT(8)
      WRITE(6,60)
   60 FORMAT(//22X,'N-1',6X,'N-2')
      WRITE(6,65)
   65 FORMAT(6X,'(SI-A)INVERSE=(S   C(0)+S   C(1)+...+SC(N-2)+C(N-1))/D(
     1S)'/)
      WRITE(6,70)
   70 FORMAT(14X,'N  N-1',6X,'N-2')
      WRITE(6,75)
   75 FORMAT(8X,'D(S)=S +S   D(1)+S   D(2)+...+SD(N-1)+D(N)'/)
C 
  100 CONTINUE
      N = NA(1)
      M = N**2
      N1= M+1
C 
      CALL UNITY(C,NA)
      CALL EQUATE(C,NA,DUMMY,NA)
      K = N + 1
C 
      DO 600 I = 1,K
      L = I-1
      IF( I .LT. K ) GO TO 200
C 
      IF( IOP .EQ. 0 )  GO TO 700
      CALL LNCNT(3)
      WRITE(6,150)
  150 FORMAT(/' ERROR IN SATISFYING CAYLEY-HAMILTON THEOREM '/)
      CALL PRNT(DUMMY,NA,'EROR',1)
      GO TO 700
C 
  200 CONTINUE
      IF( IOP .EQ. 0 ) GO TO 300
      CALL LNCNT(4)
      WRITE(6,250) L
  250 FORMAT(//' MATRIX C(',I3,' ) IN (SI-A)INVERSE EQUATION '/)
      IM = L*M + 1
      CALL PRNT(C(IM),NA,'    ',3)
  300 CONTINUE
      IF( BIDENT .AND. HIDENT ) GO TO 500
      IF( I .GT. 1) GO TO 400
      IF( BIDENT )  CALL EQUATE(H,NH,HCB,NHCB)
      IF ( HIDENT )  CALL EQUATE(B,NB,HCB,NHCB)
      IF( (.NOT. BIDENT) .AND. (.NOT. HIDENT) ) CALL MULT(H,NH,B,NB,HCB,
     1NHCB)
      IF( IOP .EQ. 0 ) GO TO 500
      CALL LNCNT(4)
      WRITE(6,350) L
  350 FORMAT(//' MATRIX HC(',I3,' )B'/)
      CALL PRNT(HCB,NHCB,'    ',3)
      GO TO 500
  400 CONTINUE
      IF( BIDENT )  IHB = L*NH(1)*NA(2) + 1
      IF( HIDENT )  IHB = L*NA(1)*NB(2) + 1
      IF( (.NOT. BIDENT) .AND. (.NOT. HIDENT) ) IHB = L*NH(1)*NB(2) + 1
      IF( BIDENT ) CALL EQUATE(C(IM),NA,DUMMY(N1),NDUM1)
      IF( .NOT. BIDENT )  CALL MULT(C(IM),NA,B,NB,DUMMY(N1),NDUM1)
      IF( HIDENT ) CALL EQUATE(DUMMY(N1),NDUM1,HCB(IHB),NDUM2)
      IF( .NOT. HIDENT )  CALL MULT(H,NH,DUMMY(N1),NDUM1,HCB(IHB),NDUM2)
      IF( IOP .EQ. 0 ) GO TO 500
      CALL LNCNT(4)
      WRITE(6,350) L
      CALL PRNT(HCB(IHB),NDUM2,'    ',3)
C 
  500 CONTINUE
      IM = I*M + 1
      IF( I .EQ. 1) CALL EQUATE(A,NA,C(IM),NA)
      IF( I .GT. 1) CALL MULT(DUMMY,NA,A,NA,C(IM),NA)
      CALL TRCE(C(IM),NA,TR)
      D(I) = -TR/I
      S = D(I)
      CALL UNITY(DUMMY,NA)
      CALL SCALE(DUMMY,NA,DUMMY,NA,S)
      CALL ADD(C(IM),NA,DUMMY,NA,C(IM),NA)
      CALL EQUATE(C(IM),NA,DUMMY,NA)
  600 CONTINUE
C 
  700 CONTINUE
      NC(1) = N
      NC(2) = N*(N+1)
      ND(1)   = N
      ND(2)   = 1
      IF( HIDENT .AND. BIDENT )  GO TO 750
      IF( HIDENT )  NHCB(1) = NA(1)
      IF( (.NOT. HIDENT) )  NHCB(1) = NH(1)
      IF( BIDENT )  NHCB(2) = N*NA(1)
      IF( (.NOT. BIDENT) )  NHCB(2) = N*NB(2)
  750 CONTINUE
C 
      IF( IOP .EQ. 0 ) RETURN
      CALL PRNT(D,ND,' D  ',1)
C 
      RETURN
      END
      SUBROUTINE LNCNT  (N)
C 
C   PURPOSE:
C      Keeps track of the number of lines printed and automatically pag-
C      inates the output.  Page length is controlled by the variable
C      NLP set in subroutine RDTITL.  LNCNT prints the contents of
C      TITLE and TIL at the top of each page.
C 
C   Subroutines employed by LNCNT: None
C   Subroutines employing LNCNT: ADD, ASYFIL, ASYREG, BARSTW, BILIN,
C      CNTREG, CSTAB, CTROL, DISREG, DSTAB, EQUATE, EXMDFL, EXPINT,
C      EXPSER, FACTOR, IMMDFL, JUXTC, JUXTR, MULT, PREFIL, PRNT, RDTITL,
C      READ1, RICNWT, SAMPL, SCALE, SUBT, SUM, TESTST, TRANP, TRCE,
C      TRNSIT, UNITY, VARANC
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER:: nlp,lin
      CHARACTER(LEN=80):: title
      COMMON/LINES/NLP,LIN,TITLE    !!!(10),TIL(2)
      LIN=LIN+N
      IF  (LIN.LE.NLP)   GO TO 20
      WRITE(6,*) ACHAR(12)   ! new page
      WRITE(6,*) TITLE
!!! 1010 FORMAT(1H1,10A8,2A8/)
      LIN=2+N
      IF  (N.GT.NLP)  LIN=2
   20 RETURN
      END
      SUBROUTINE MAXEL(A,NA,ELMAX)
C 
C   PURPOSE:
C      Compute the maximum of the absolute values of the elements of a
C      real matrix.
C 
C   Subroutines employed by MAXEL: None
C   Subroutines employing MAXEL: CNTREG, DISREG, EXPINT, EXPSER, RICNWT,
C      SAMPL, SUM
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),NA(2)
C 
      N = NA(1)*NA(2)
C 
      ELMAX = DABS( A(1))
      DO 100 I = 2,N
      ELMAXI = DABS( A(I) )
      IF( ELMAXI .GT. ELMAX )  ELMAX = ELMAXI
  100 CONTINUE
C 
      RETURN
      END
      SUBROUTINE MULT(A,NA,B,NB,C,NC)
C 
C   PURPOSE:
C      Perform matrix multiplication C = AB for given matrices A and B.
C 
C   Subroutines employed by MULT: LNCNT
C   Subroutines employing MULT: ASYREG, BILIN, CNTREG, CSTAB, CTROL,
C      DISREG, DSTAB, EXMDFL, EXPINT, EXPSER, FACTOR, IMMDFL, PREFIL,
C      RICNWT, SAMPL, SUM, TRNSIT, VARANC
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),NA(2),NB(2),NC(2)
C     DOUBLE PRECISION V1,V2,V3,V4
      NC(1) = NA(1)
      NC(2) = NB(2)
      IF(NA(2).NE.NB(1)) GO TO 999
      NAR = NA(1)
      NAC = NA(2)
      NBC = NB(2)
      NAA=NAR*NAC
      NBB=NAR*NBC
      IF ( NAR .LT. 1 .OR. NAA .LT. 1 .OR. NBB .LT. 1 ) GO TO 999
      IR = 0
      IK=-NAC
      DO 350 K=1,NBC
      IK = IK + NAC
      DO 350 J=1,NAR
      IR=IR+1
      IB=IK
      JI=J-NAR
      V1=0.0
      DO 300 I=1,NAC
      JI = JI + NAR
      IB=IB+1
      V3=A(JI)
      V4=B(IB)
      V2=V3*V4
      V1=V1+V2
  300 CONTINUE
      C(IR)=V1
  350 CONTINUE
      GO TO 1000
  999 CALL LNCNT (1)
      WRITE(6,500) (NA(I),I=1,2),(NB(I),I=1,2)
  500 FORMAT  (' DIMENSION ERROR IN MULT    NA=',2I6,5X,'NB=',2I6)
 1000 RETURN
      END
      SUBROUTINE NORMS(MAXROW,M,N,A,IOPT,RLNORM)
C 
C   PURPOSE:
C      Compute either the l(1), l(2) (Euclidean), or l(infinite) matrix
C      norms for a real m x n matrix A stored as a variable-dimensioned
C      two-dimensional array.
C 
C   Subroutines employed by NORMS: None
C   Subroutines employing NORMS: BILIN, CSTAB, EXPINT, EXPSER, SAMPL
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1)
C 
C  INITIALIZATION
C 
      SUM=0.
      RLNORM=0.
      I=-MAXROW
C 
C  TRANSFER TO APPROPRIATE LOOP TO COMPUTE THE DESIRED NORM
C 
      IF(IOPT-2)5,20,30
C 
C  THIS LOOP COMPUTES THE ONE-NORM
C 
    5 DO 15 K=1,N
      I=I+MAXROW
      DO 10 J=1,M
      L=I+J
   10 SUM=DABS(A(L))+SUM
      IF(SUM.GT.RLNORM)RLNORM=SUM
   15 SUM=0.
      RETURN
C 
C  THIS LOOP COMPUTES THE EUCLIDEAN NORM
C 
   20 DO 25 K=1,N
      I=I+MAXROW
      DO 25 J=1,M
      L=I+J
      SUM=A(L)
   25 RLNORM=SUM*SUM+RLNORM
      RLNORM=DSQRT(RLNORM)
      RETURN
C 
C  THIS LOOP COMPUTES THE INFINITY-NORM
C 
   30 DO 40 J=1,M
      L=I+J
      DO 35 K=1,N
      L=L+MAXROW
   35 SUM=DABS(A(L))+SUM
      IF(SUM.GT.RLNORM)RLNORM=SUM
   40 SUM=0.
      RETURN
      END
      SUBROUTINE NULL(A,NA)
C 
C   PURPOSE:
C      Generate a null matrix.
C 
C   Subroutines employed by NULL: LNCNT
C   Subroutines employing NULL: BARSTW, CNTREG
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1)
      DIMENSION NA(2)
      N=NA(1)*NA(2)
      IF( NA(1) .LT. 1 .OR.  N .LT. 1 )  GO TO 999
      DO 10I=1,N
   10 A(I) = 0.0
      RETURN
C 
  999 CONTINUE
      WRITE (6,50) NA
   50 FORMAT(' DIMENSION ERROR IN NULL  NA =',2I6)
      RETURN
      END
      SUBROUTINE POLE(A,NA,B,NB,EVAL,NUMR,F,NF,IOP,DUMMY)
C 
C   PURPOSE:
C      Calculate an (1 x n) matrix F which causes the eigenvalues
C      of the matrix A-BF to assume n specified values.
C 
C   Subroutines employed by POLE: DGECO, DGESL, EIGEN, EQUATE, JUXTC,
C      LNCNT, MULT, NULL, PRNT, SCALE, SUBT, TRANP, TRCE
C   Subroutines employing POLE: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),EVAL(1),F(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NF(2),NDUM1(2),NDUM2(2)
C 
C 
      IF( IOP .EQ. 0 ) GO TO 100
      CALL LNCNT(6)
      WRITE(6,25)
   25 FORMAT('  GIVEN THE COMPLETELY CONTROLLABLE PAIR (A,B), THE MATRIX
     1 F IS COMPUTED'/' SUCH THAT A-BF HAS PRESCRIBED EIGENVALUES'///)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      CALL LNCNT(6)
      WRITE(6,50) NUMR
   50 FORMAT(///' THE PRESCRIBED EIGENVALUES ARE STORED IN EVAL'/' PAST
     1ENTRY ',I4,' NUMBERS CORRESPOND TO REAL AND IMAGINARY PARTS'/)
      CALL PRNT(EVAL,NB,'EVAL',1)
  100 CONTINUE
C 
      N = NA(1)
      M = N**2
      N1 = M + 1
      N2 = N1 + M
      N3 = N2 + M
      N4 = N3 + M
      J = 1
      CALL EQUATE(B,NB,DUMMY,NDUM2)
      CALL EQUATE(B,NB,DUMMY(N1),NB)
  200 CONTINUE
      CALL MULT(A,NA,DUMMY(N1),NB,DUMMY(N2),NB)
      IF( J .EQ. N ) GO TO 300
      CALL JUXTC(DUMMY,NDUM2,DUMMY(N2),NB,DUMMY(N3),NDUM1)
      CALL EQUATE(DUMMY(N2),NB,DUMMY(N1),NB)
      CALL EQUATE(DUMMY(N3),NDUM1,DUMMY,NDUM2)
      J = J + 1
      GO TO 200
  300 CONTINUE
C 
      CALL SCALE(DUMMY(N2),NB,F,NF,-1.0D0)
      CALL EQUATE(DUMMY,NA,DUMMY(N4),NA)
      NRHS = 1
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY,N,N,DUMMY(N1),RCOND,DUMMY(N2))
      IF ((1.0 + RCOND) .NE. 1.0) GO TO 400
      CALL LNCNT(1)
      WRITE(6,350) RCOND
  350 FORMAT(' IN POLE, THE CONTROLLABILITY MATRIX, C, HAS BEEN DECLARED
     1 SINGULAR BY DGECO, RCOND = ',D16.8)
      RETURN
  400 CONTINUE
      NT = 1
      DO 420 M1 = 1,NRHS
         CALL DGESL(DUMMY,N,N,DUMMY(N1),F(NT),0)
         NT = NT + N
  420 CONTINUE
C 
      IF( IOP .EQ. 0 ) GO TO 500
      CALL LNCNT(10)
      WRITE(6,425)
  425 FORMAT(//' DENOTE CHARACTERISTIC POLYNOMIAL OF A BY D(S)')
      WRITE(6,450)
  450 FORMAT(14X,'N  N-1',6X,'N-2')
      WRITE(6,475)
  475 FORMAT(8X,'D(S)=S +S   D(1)+S   D(2)+...+SD(N-1)+D(N)'//)
      WRITE(6,485)
  485 FORMAT(//'  AL = COL.(D(N),D(N-1),...,D(1))')
      CALL PRNT(F,NF,' AL ',1)
C 
  500 CONTINUE
      NDUM1(1) = 2
      NDUM1(2) = 2
      CALL NULL(DUMMY,NDUM1)
      IF( NUMR .EQ. 0 )   GO TO 535
C 
      DO 525 I =1,N
      DO 525 J =1,NUMR
      DUMMY(I) = DUMMY(I) + EVAL(J)**I
  525 CONTINUE
  535 CONTINUE
      N11 = N1
      N21 = N11 + 1
      N12 = N21 + 1
      N22 = N12 + 1
      IF( NUMR .EQ. N ) GO TO 600
      L = NUMR + 1
      DO 550 I= L,N,2
      DUMMY(N11) =  EVAL(I)
      DUMMY(N12) =  EVAL(I+1)
      DUMMY(N21) = -EVAL(I+1)
      DUMMY(N22) =  EVAL(I)
      CALL EQUATE(DUMMY(N1),NDUM1,DUMMY(N2),NDUM1)
      DO 550 K =1,N
      CALL TRCE(DUMMY(N2),NDUM1,TR)
      DUMMY(K) = DUMMY(K) + TR
      IF( K .EQ. N ) GO TO 550
      CALL MULT(DUMMY(N1),NDUM1,DUMMY(N2),NDUM1,DUMMY(N3),NDUM1)
      CALL EQUATE(DUMMY(N3),NDUM1,DUMMY(N2),NDUM1)
  550 CONTINUE
  600 CONTINUE
      NDUM1(1) = N
      NDUM1(2) = 1
      DUMMY(N1) = - DUMMY(1)
      DO 650 K = 2,N
      NK = N1+K-1
      DUMMY(NK) = DUMMY(K)
      KM1 = K-1
      DO 625 J = 1,KM1
      NJ = N1+J-1
      NKJ = NK-NJ
      DUMMY(NK) = DUMMY(NK)+ DUMMY(NJ)*DUMMY(NKJ)
  625 CONTINUE
      DUMMY(NK) = - DUMMY(NK)/K
  650 CONTINUE
      DO 675 J = 1,N
      I = N+1-J
      DUMMY(J) = DUMMY(N1-1+I)
  675 CONTINUE
C 
      IF( IOP .EQ. 0 ) GO TO 800
      CALL LNCNT(7)
      WRITE(6,700) N
  700 FORMAT(//'  FOR THE ',I4,' EIGENVALUES STORED IN EVAL THE CHARACTE
     1RISTIC POLYNOMIAL IS'/14X,'N  N-1      N-2'/8X,'D(S)=S +S   D(1)+S
     2   D(2)+...+SD(N-1)+D(N)'//)
      WRITE(6,750)
  750 FORMAT(//'  ALTL = COL.(D(N),D(N-1),...,D(1))')
      CALL PRNT(DUMMY,NB,'ALTL',1)
C 
  800 CONTINUE
      CALL NULL(DUMMY(N1),NA)
      DO 850 I =1,N
      N1I = I*N-I+N1
      DUMMY(N1I) = 1.0
  850 CONTINUE
      NM1 = N-1
      DO 875 J = 1,NM1
      K = J+1
      DO 875 I = K,N
      NJ = N1 + (J-1)*N + I-K
      DUMMY(NJ) = F(I)
  875 CONTINUE
      CALL MULT(DUMMY(N4),NA,DUMMY(N1),NA,DUMMY(N2),NA)
      CALL SUBT(DUMMY,NB,F,NF,F,NF)
      CALL TRANP(DUMMY(N2),NA,DUMMY,NA)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(DUMMY,N,N,DUMMY(N1),RCOND,DUMMY(N2))
      IF ((1.0 + RCOND) .NE. 1.0) GO TO 900
      CALL LNCNT(1)
      WRITE(6,885) RCOND
  885 FORMAT('  IN POLE, THE MATRIX V HAS BEEN DECLARED SINGULAR BY DGEC
     1O, RCOND = ',D16.8)
      RETURN
C 
  900 CONTINUE
      NT = 1
      DO 925 M1 = 1,NRHS
         CALL DGESL(DUMMY,N,N,DUMMY(N1),F(NT),0)
         NT = NT + N
  925 CONTINUE
      NF(1) = 1
      NF(2) = N
      CALL MULT(B,NB,F,NF,DUMMY,NDUM1)
      CALL SUBT(A,NA,DUMMY,NDUM1,DUMMY,NDUM1)
      CALL EQUATE(DUMMY,NDUM1,DUMMY(N1),NDUM1)
      ISV = N
      ILV = 0
      N21 = N2
      N22 = N2+N
      N23 = N22+N
      CALL EIGEN(N,N,DUMMY(N1),DUMMY(N21),DUMMY(N22),ISV,ILV,DUMMY(N3),D
     1UMMY(N4),IERR)
      IF( IERR .EQ. 0 ) GO TO 1000
      CALL LNCNT(1)
      WRITE(6,950)
  950 FORMAT('  IN POLE, AFTER CALL TO EIGEN WITH  A-BF, IERR= ',I4)
C 
 1000 CONTINUE
      CALL EQUATE(DUMMY(N3),NA,DUMMY(N23),NA)
      IF( IOP .EQ. 0 ) RETURN
      CALL PRNT(F,NF,' F  ',1)
      CALL PRNT(DUMMY,NA,'A-BF',1)
      CALL LNCNT(3)
      WRITE(6,1100)
 1100 FORMAT(//'  EIGENVALUES(VAL MATRIX) AND EIGENVECTORS(VECT MATRIX)
     1OF A-BF')
      NDUM1(1) = N
      NDUM1(2) = 2
      CALL PRNT(DUMMY(N2),NDUM1,' VAL',1)
      CALL PRNT(DUMMY(N23),NA,'VECT',1)
C 
      RETURN
      END
      SUBROUTINE PREFIL(A,NA,B,NB,Q,NQ,W,NW,R,NR,F,NF,IOP,DUMMY)
C 
C   PURPOSE:
C      Compute an r x n (r <= n) matrix F which, when used in the vector
C      equation u = -Fx + v eliminates the cross-product term in the
C      quadratic scalar function, x'Qx + x'Wu + u'Ru, where Q = Q'>=0,
C      W, and R=R'>0 are constant matrices.  Specifically, F=(R**-1)
C      (W/2).
C 
C   Subroutines employed by PREFIL: DPOCO, DPOSL, EQUATE, LNCNT, MULT,
C      PRNT, SCALE, SUBT, TRANP
C   Subroutines employing PREFIL: IMMDFL
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),Q(1),W(1),R(1),F(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NQ(2),NW(2),NR(2),NF(2),IOP(3)
C 
      IF( IOP(1) .EQ. 0 ) GO TO 100
      CALL LNCNT(5)
      WRITE(6,25)
   25 FORMAT(// ' PROGRAM TO COMPUTE PREFILTER GAIN F TO ELIMINATE  CROS
     1S-PRODUCT TERM '/' IN QUADRATIC PERFORMANCE INDEX '/)
      IF( IOP(3) .EQ. 0 ) GO TO 50
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
   50 CONTINUE
      CALL PRNT(Q,NQ,' Q  ',1)
      CALL PRNT(W,NW,' W  ',1)
      CALL PRNT(R,NR,' R  ',1)
C 
  100 CONTINUE
      CALL TRANP(W,NW,F,NF)
      CALL SCALE(F,NF,F,NF,0.5D0)
      CALL EQUATE(R,NR,DUMMY,NR)
      N1 = NR(1)**2 + 1
      M = NR(1)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DPOCO(DUMMY,M,M,RCOND,DUMMY(N1),IERR)
      IF( IERR .EQ. 0 ) GO TO 200
      CALL LNCNT(4)
      WRITE(6,150)
  150 FORMAT(//' IN PREFIL, THE MATRIX R IS NOT POSITIVE DEFINITE'/)
      RETURN
C 
  200 CONTINUE
      NT = 1
      DO 250 M1 = 1,NF(2)
         CALL DPOSL(DUMMY,M,M,F(NT))
         NT = NT + M
  250 CONTINUE
      IF( IOP(2) .EQ. 0 ) GO TO 300
      CALL MULT(W,NW,F,NF,DUMMY,NQ)
      CALL SCALE(DUMMY,NQ,DUMMY,NQ,0.5D0)
      CALL SUBT(Q,NQ,DUMMY,NQ,Q,NQ)
C 
  300 CONTINUE
      IF( IOP(3) .EQ. 0 ) GO TO 400
      CALL MULT(B,NB,F,NF,DUMMY,NA)
      CALL SUBT(A,NA,DUMMY,NA,A,NA)
C 
  400 CONTINUE
      IF( IOP(1) .EQ. 0 ) RETURN
      CALL PRNT(F,NF,' F  ',1)
      IF( IOP(2) .EQ. 0 ) GO TO 500
      CALL LNCNT(3)
      WRITE(6,450)
  450 FORMAT(/ ' MATRIX  Q - (W/2)F '/)
      CALL PRNT(Q,NQ,'NEWQ',1)
C 
  500 CONTINUE
      IF( IOP(3) .EQ. 0 ) RETURN
      CALL PRNT(A,NA,'NEWA',1)
      RETURN
      END

        SUBROUTINE PRNT(A,NA,NAM,IOP)
C
C   PURPOSE:
C      Print a single matrix with or without a descriptive heading
C      either on the same page or on a new page.  Format for the print-
C      ing is stored in the COMMON block FORM.
C
C   Subroutines employed by PRNT: LNCNT
C   Subroutines employing PRNT: ASYFIL, ASYREG, BARSTW, BILIN, CNTREG,
C      CSTAB, CTROL, DISREG, DSTAB, EXMDFL, EXPINT, EXPSER, FACTOR,
C      IMMDFL, PREFIL, READ1, RICNWT, SAMPL, SUM, TESTST, TRNSIT,
C      VARANC
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(*),NA(2)
        CHARACTER(LEN=4),INTENT(IN):: NAM
        INTEGER:: nepr
        CHARACTER(LEN=16):: fmt1,fmt2
        COMMON/FORM/nepr,fmt1,fmt2
!!!     COMMON/FORM/NEPR,FMT1(2),FMT2(2)
        INTEGER:: nlp,lin
        CHARACTER(LEN=80):: title
        COMMON/LINES/NLP,LIN,TITLE   !!!(10),TIL(2)
C- NOTE NLP NO. LINES/PAGE VARIES WITH THE INSTALLATION.

C       changed 15 Dec 2009 by RLC
CCC     DATA KZ,KW,KB /1H0,1H1,1H /
        CHARACTER(LEN=1),PARAMETER:: KZ='0', KW='1', KB=' ' 
        CHARACTER(LEN=4):: name

        NAME = NAM
        II = IOP
        NR = NA(1)
        NC = NA(2)
        NLST = NR * NC
        IF( NLST .LT. 1 .OR. NR .LT. 1 ) GO TO 999

CCC     IF(NAME  .EQ. 0) NAME = KB

C- SKIP HEADLINE IF REQUESTED.
        GO TO (20,10,40,30),      II
   10   CALL LNCNT(100)
   20   CALL LNCNT(2)
        WRITE(6,100) NAME,NR,NC
  100   FORMAT(/5X,A4,'  MATRIX',5X,I3,' ROWS',5X,I3,' COLUMNS')
        GO TO 50

   30   CALL LNCNT(100)
        GO TO 50

   40   CALL LNCNT(2)
        WRITE (6,150)
  150   FORMAT (/)

   50   CONTINUE
        DO I = 1,NC,NEPR
           J = I + (NEPR-1)
           IF (J .GT. NC) J=NC

           IF (NC .GT. NEPR) THEN
              CALL LNCNT(2)
              IF (I .EQ. J) THEN
                WRITE(6,200) I
  200           FORMAT(/'    COL ',I4)
              ELSE
                WRITE(6,250) I,J
  250           FORMAT(/'               COLS ',I4,' TO ',I4)
              END IF
           END IF

           DO K = 1,NR
              CALL LNCNT(1)
              WRITE(6,FMT1) (A(NR*(L-1)+K),L=I,J)
           END DO

        END DO
        RETURN

C   * * * ERROR MESSAGE * * *
  999   CALL LNCNT(1)
        WRITE (6,300) NAM,NA
  300   FORMAT  (' ERROR IN PRNT  MATRIX ',A4,' HAS NA=',2I6)
        RETURN
        END
      SUBROUTINE RDTITL
C 
C   PURPOSE:
C      Read a single card of Hollerith input which is loaded into the
C      array TITLE of COMMON block LINES of RDTITL and automatically
C      printed at the top of each page of output through the subroutine
C      LNCNT.  Subroutine RDTITL also initializes several COMMON blocks
C      used in ORACLS.  RDTITL must be in every executive program.
C 
C   Subroutines employed by RDTITL: LNCNT
C   Subroutines employing RDTITL: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER:: nlp,lin
      CHARACTER(LEN=80):: title
      COMMON/LINES/NLP,LIN,TITLE   !!! (10),TIL(2)
      INTEGER:: nepr
      CHARACTER(LEN=16):: fmt1,fmt2
      COMMON/FORM/nepr,fmt1,fmt2
!!!      COMMON/FORM/NEPR,FMT1(2),FMT2(2)
      COMMON/TOL/EPSAM,EPSBM,IACM
C      COMMON/CONV/SUMCV,MAXSUM,RICTCV,SERCV
      COMMON/CONV/SUMCV,RICTCV,SERCV,MAXSUM
C     NLP = NO. LINES/PAGE VARIES WITH THE INSTALLATION
!!!      DATA LIN,NLP/1,60/
!!!      DATA NEPR,FMT1/7,8H(1P7D16.,8H7)      /
!!!      DATA TIL/8H ORACLS ,8H PROGRAM/
!!!      DATA FMT2/8H(3X,1P7D,8H16.7)   /
      NEPR=7
      FMT1 = '(7D16.7)'
      FMT2 = '(3X,7D16.7)'
      LIN=1
      NLP=60
      EPSAM=1.D-10
      EPSBM=1.D-10
      IACM=12
      SUMCV=1.D-8
      RICTCV=1.D-8
      SERCV=1.D-8
      MAXSUM=50
      READ(5,'(A)',END=90) TITLE
   91 CONTINUE
  100 FORMAT(10A8)
      CALL LNCNT(100)
      RETURN
   90 CONTINUE
      STOP 1
      END
      SUBROUTINE READ(I,A,NA,B,NB,C,NC,D,ND,E,NE)
C 
C   PURPOSE:
C      Read from one to five matrices along with their names and dimen-
C      sions and print the same information.
C 
C   Subroutines employed by READ: READ1
C   Subroutines employing READ: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),C(*),D(*),E(*)
      DIMENSION NA(2),NB(2),NC(2),ND(2),NE(2),NZ(2)
      CHARACTER(LEN=4):: LAB
      READ(5,100,END=10000) LAB,            NZ(1), NZ(2)
10000 CALL READ1(A, NA,NZ,  LAB)
      IF(I .EQ. 1) GO TO 999
      READ(5,100,END=10001) LAB,            NZ(1), NZ(2)
10001 CALL READ1(B, NB,NZ,  LAB)
      IF(I .EQ. 2) GO TO 999
      READ(5,100,END=10002) LAB,            NZ(1), NZ(2)
10002 CALL READ1(C,NC,NZ,LAB)
      IF(I .EQ. 3) GO TO 999
      READ(5,100,END=10003) LAB,            NZ(1), NZ(2)
10003 CALL READ1(D, ND,NZ,  LAB)
      IF(I .EQ. 4) GO TO 999
      READ(5,100,END=10004) LAB,            NZ(1), NZ(2)
10004 CALL READ1(E, NE,NZ,  LAB)
  100 FORMAT(BZ,A4,4X,2I4)
  999 RETURN
      END
      SUBROUTINE READ1 (A,NA,NZ,NAM)
C 
C   PURPOSE:
C      Read in a single matrix and print the matrix using subroutine
C      PRNT. Each row of the matrix starts on a new line.
C 
C   Subroutines employed by READ1: LNCNT, PRNT
C   Subroutines employing READ1: READ
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*)
      DIMENSION NA(2),NZ(2)
      CHARACTER(LEN=4),INTENT(IN):: nam
      IF  (NZ(1).EQ.0)  GO TO 410
      NR=NZ(1)
      NC=NZ(2)
      NLST=NR*NC
      IF( NLST .LT. 1 .OR. NR .LT. 1 ) GO TO 16
      DO 400 I = 1, NR
      READ (5,101,END=  400) (A(  J), J = I,NLST,NR)
  400 CONTINUE
      NA(1)=NR
      NA(2)=NC
  410 CALL  PRNT (A,NA,NAM,1)
  101 FORMAT(BZ,8E10.2)
      RETURN
   16 CALL LNCNT(1)
      WRITE  (6,916)  NAM,NR,NC
  916 FORMAT  (' ERROR IN READ1   MATRIX ',A4,' HAS NA=',2I6)
      RETURN
      END
      SUBROUTINE RICNWT(A,NA,B,NB,H,NH,Q,NQ,R,NR,F,NF,P,NP,IOP,IDENT,DI
     1SC,FNULL,DUMMY)
C 
C   PURPOSE:
C      Solve either the continuous or discrete steady-state Riccati
C      equation by the Newton algorithms described by Kleinman and
C      Hewer.
C 
C   REFERENCES:
C      Kleinman, David L.: On an Iterative Technique for Riccati Equa-
C        tion Computations. IEEE Trans. Autom. Control, vol. AC-13, no.
C        1, Feb. 1968, pp. 114-115.
C      Hewer, Gary A.: An Iterative Technique for the Computation of
C        the Steady State Gains for the Discrete Optimal Regulator.
C        IEEE Trans. Autom. Control, vol. AC-16, no. 4, Aug. 1971,
C        pp. 382-383.
C      Kantorovich, L.V.; and Akilov, G.P. (D. E. Brown, trans.): Func-
C        tional Analysis in Normed Spaces. A.P. Robertson, ed., Mac-
C        Millan Co., 1964.
C      Sandell, Nils R., Jr.: On Newton's Method for Riccati Equation
C        Solution. IEEE Trans. Autom. Control, vol. AC-19, no. 3, June
C        1974, pp. 254-255.
C      Kwakernaak, Huibert; and Sivan, Raphael: Linear Optimal Control
C        Systems. John Wiley & Sons, Inc., c. 1972.
C 
C   Subroutines employed by RICNWT: ADD, BARSTW, BILIN, DPOCO, DPOSL,
C      EQUATE, GELIM, LNCNT, MAXEL, MULT, PRNT, SCALE, SUBT, SUM,
C      SYMPDS, TRANP, UNITY
C   Subroutines employing RICNWT: ASYREG
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),Q(1),R(1),F(1),P(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NQ(2),NR(2),NF(2),NP(2),IOP(3)
      DIMENSION H(1),NH(2),IOPT(2)
      LOGICAL  IDENT,DISC,FNULL,SYM
      COMMON/TOL/EPSAM,EPSBM,IACM
C      COMMON/CONV/SUMCV,MAXSUM,RICTCV,SERCV
      COMMON/CONV/SUMCV,RICTCV,SERCV,MAXSUM
C 
      I=1
      IOPT(1)=0
      SYM = .TRUE.
C 
      N = NA(1)**2
      N1 = N +1
      IF( .NOT. DISC) N1 = NA(1)*NR(1) + 1
      N2= N1+N
      N3= N2+N
      N4 = N3+N
C 
      IF( IOP(1) .EQ. 0 )  GO TO 210
      CALL LNCNT(4)
      IF(.NOT. DISC)WRITE(6,100)
      IF( DISC )WRITE(6,150)
  100 FORMAT(BZ,//' PROGRAM TO SOLVE CONTINUOUS STEADY-STATE RICCATI EQU
     1ATION BY THE NEWTON ALGORITHM'/)
  150 FORMAT(//' PROGRAM TO SOLVE DISCRETE STEADY-STATE RICCATI EQUATION
     1 BY THE NEWTON ALGORITHM'/)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      CALL PRNT(Q,NQ,' Q  ',1)
      IF( .NOT. IDENT )GO TO 185
      CALL LNCNT(3)
      WRITE(6,180)
  180 FORMAT(/' H IS AN IDENTITY MATRIX'/)
      GO TO 200
  185 CONTINUE
      CALL PRNT(H,NH,' H  ',1)
      CALL MULT(Q,NQ,H,NH,DUMMY,NH)
      CALL TRANP(H,NH,DUMMY(N2),NP)
      CALL MULT(DUMMY(N2),NP,DUMMY,NH,Q,NQ)
      CALL LNCNT(3)
      WRITE(6,195)
  195 FORMAT(/' MATRIX (H TRANSPOSE)QH '/)
      CALL PRNT(Q,NQ,'HTQH',1)
  200 CONTINUE
      CALL PRNT(R,NR,' R  ',1)
      IF( FNULL )  GO TO 210
      CALL LNCNT(3)
      WRITE(6,205)
  205 FORMAT(/' INITIAL F MATRIX'/)
      CALL PRNT(F,NF,' F  ',1)
C 
  210 CONTINUE
      IF((IOP(1) .NE. 0)  .OR. IDENT)  GO TO 220
      CALL MULT(Q,NQ,H,NH,DUMMY,NH)
      CALL TRANP(H,NH,DUMMY(N2),NP)
      CALL MULT(DUMMY(N2),NP,DUMMY,NH,Q,NQ)
  220 CONTINUE
C 
      IF (DISC) GO TO 900
C 
      CALL TRANP(B,NB,P,NP)
      CALL EQUATE(R,NR,DUMMY,NR)
      CALL SYMPDS(NR(1),NR(1),DUMMY,NR(2),P,IOPT,IOPT,DET,ISCALE,DUMMY(N
     11),IERR)
      IF(IERR .NE. 0) GO TO 240
C      GENERATE R-INVERSE
      CALL UNITY(P,NR)
      CALL GELIM(NR(1),NR(1),DUMMY,NR(2),P,DUMMY(N2),0,DUMMY(N1),IERR)
      IF(IERR.EQ.0) GO TO 250
  240 CALL LNCNT(3)
      WRITE(6,225)
  225 FORMAT(/' IN RICNWT, A MATRIX WHICH IS  NOT SYMMETRIC POSITIVE DE
     1FINITE HAS BEEN SUBMITTED TO  SYMPDS'/)
      RETURN
C 
  250 CONTINUE
C     CALCULATE (R-INVERSE)*(B-TRANSPOSE) AND STORE IN DUMMY(N1)
      CALL TRANP(B,NB,DUMMY(N2),NF)
      CALL EQUATE(P,NR,DUMMY,NR)
      CALL MULT(DUMMY,NR,DUMMY(N2),NF,DUMMY(N1),NF)
      CALL PRNT(DUMMY,NR,'RINV',1)
C 
      IF(FNULL) GO TO 300
C 
C     GENERATE THE RECURSIVE EQUATION FOR THE INITIAL RUN
      CALL MULT(B,NB,F,NF,DUMMY(N2),NA)
      CALL SUBT(A,NA,DUMMY(N2),NA,DUMMY(N2),NA)
      CALL TRANP(DUMMY(N2),NA,DUMMY(N3),NA)
      CALL EQUATE(DUMMY(N3),NA,DUMMY(N2),NA)
      CALL MULT(R,NR,F,NF,DUMMY(N3),NF)
      CALL TRANP(F,NF,P,NP)
      CALL MULT(P,NP,DUMMY(N3),NF,DUMMY(N4),NA)
      CALL TRANP(DUMMY(N4),NA,DUMMY(N3),NA)
      CALL ADD(DUMMY(N4),NA,DUMMY(N3),NA,DUMMY(N3),NA)
      CALL SCALE(DUMMY(N3),NA,DUMMY(N3),NA,0.5D0)
      CALL ADD(DUMMY(N3),NA,Q,NQ,P,NP)
      CALL SCALE(P,NP,P,NP,-1.0D0)
      GO TO 350
C 
  300 CONTINUE
      CALL TRANP(A,NA,DUMMY(N2),NA)
      CALL SCALE(Q,NQ,P,NP,-1.0D0)
C 
  350 CONTINUE
      IF(IOP(3) .NE. 0) GO TO 400
      EPSA= EPSAM
      CALL BARSTW(DUMMY(N2),NA,B,NB,P,NP,IOPT,SYM,EPSA,EPSA,DUMMY(N3))
      GO TO 450
C 
  400 CONTINUE
      IOPT(2)=1
      CALL BILIN(DUMMY(N2),NA,B,NB,P,NP,IOPT,SCLE,SYM,DUMMY(N3))
C 
  450 CONTINUE
      CALL EQUATE(P,NP,DUMMY(N2),NP)
      IF(IOP(2).EQ. 0) GO TO 550
      CALL LNCNT(3)
      WRITE(6,500) I
  500 FORMAT(/' ITERATION  ',I5/)
      CALL PRNT(P,NP,' P  ',1)
C 
  550 CONTINUE
C     GENERATE THE RECURSIVE EQUATION ON THE LOOP
      CALL MULT(DUMMY(N1),NF,P,NP,DUMMY(N3),NF)
      CALL EQUATE(DUMMY(N3),NF,F,NF)
      CALL MULT(R,NR,F,NF,DUMMY(N3),NF)
      CALL TRANP(F,NF,P,NP)
      CALL MULT(P,NP,DUMMY(N3),NF,DUMMY(N4),NA)
      CALL TRANP(DUMMY(N4),NA,DUMMY(N3),NA)
      CALL ADD(DUMMY(N4),NA,DUMMY(N3),NA,DUMMY(N3),NA)
      CALL SCALE(DUMMY(N3),NA,DUMMY(N3),NA,0.5D0)
      CALL ADD(DUMMY(N3),NA,Q,NQ,P,NP)
      CALL SCALE(P,NP,P,NP,-1.0D0)
      CALL MULT(B,NB,F,NF,DUMMY(N3),NA)
      CALL SUBT(A,NA,DUMMY(N3),NA,DUMMY(N4),NA)
      CALL TRANP(DUMMY(N4),NA,DUMMY(N3),NA)
C 
      IF(IOP(3) .NE. 0 ) GO TO 650
      CALL BARSTW(DUMMY(N3),NA,B,NB,P,NP,IOPT,SYM,EPSA,EPSA,DUMMY(N4))
      GO TO 675
C 
  650 CONTINUE
      CALL BILIN(DUMMY(N3),NA,B,NB,P,NP,IOPT,SCLE,SYM,DUMMY(N4))
C 
  675 CONTINUE
      I=I+1
      CALL MAXEL(DUMMY(N2),NA,ANORM1)
      CALL SUBT(P,NP,DUMMY(N2),NA,DUMMY(N3),NA)
      CALL MAXEL(DUMMY(N3),NA,ANORM2)
      IF(ANORM1 .GT. 1.0) GO TO 700
      IF( ANORM2/ANORM1 .LT. RICTCV ) GO TO 800
      GO TO 750
C 
  700 CONTINUE
      IF( ANORM2 .LT. RICTCV ) GO TO 800
C 
  750 CONTINUE
      IF( I .LE. 101) GO TO 450
      CALL LNCNT(3)
      WRITE(6,775)
  775 FORMAT(/' THE SUBROUTINE RICNWT HAS EXCEEDED 100 ITERATIONS WITHO
     1UT CONVERGENCE'/)
      IOP(1) = 1
C 
  800 CONTINUE
      CALL MULT(DUMMY(N1),NF,P,NP,F,NF)
      GO TO 1300
C 
  900 CONTINUE
      IF( .NOT. FNULL ) GO TO 950
C 
      CALL EQUATE(Q,NQ,P,NP)
      CALL EQUATE(A,NA,DUMMY(N1),NA)
      CALL TRANP(A,NA,DUMMY(N2),NA)
      GO TO 1000
  925 CONTINUE
C 
      I=I+1
      CALL EQUATE(P,NP,DUMMY,NP)
  950 CONTINUE
C 
      CALL MULT(R,NR,F,NF,DUMMY(N1),NF)
      CALL TRANP(F,NF,P,NP)
      CALL MULT(P,NP,DUMMY(N1),NF,DUMMY(N2),NA)
      CALL TRANP(DUMMY(N2),NA,DUMMY(N1),NA)
      CALL ADD(DUMMY(N1),NA,DUMMY(N2),NA,DUMMY(N1),NA)
      CALL SCALE(DUMMY(N1),NA,DUMMY(N1),NA,0.5D0)
      CALL ADD(Q,NQ,DUMMY(N1),NA,P,NP)
      CALL MULT(B,NB,F,NF,DUMMY(N1),NA)
      CALL SUBT(A,NA,DUMMY(N1),NA,DUMMY(N1),NA)
      CALL TRANP(DUMMY(N1),NA,DUMMY(N2),NA)
C 
 1000 CONTINUE
      CALL SUM(DUMMY(N2),NA,P,NP,DUMMY(N1),NA,IOPT,SYM,DUMMY(N3))
      IF(IOP(2) .EQ. 0) GO TO 1100
      CALL LNCNT(3)
      WRITE(6,500) I
      CALL PRNT(P,NP,' P  ',1)
C 
 1100 CONTINUE
      CALL MULT(P,NP,A,NA,DUMMY(N1),NA)
      CALL MULT(P,NP,B,NB,DUMMY(N2),NB)
      CALL TRANP(B,NB,DUMMY(N3),NF)
      CALL MULT(DUMMY(N3),NF,DUMMY(N1),NA,F,NF)
      CALL MULT(DUMMY(N3),NF,DUMMY(N2),NB,DUMMY(N1),NR)
      CALL TRANP(DUMMY(N1),NR,DUMMY(N2),NR)
      CALL ADD(DUMMY(N1),NR,DUMMY(N2),NR,DUMMY(N1),NR)
      CALL SCALE(DUMMY(N1),NR,DUMMY(N1),NR,0.5D0)
      CALL ADD(R,NR,DUMMY(N1),NR,DUMMY(N1),NR)
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DPOCO(DUMMY(N1),NR(1),NR(1),RCOND,DUMMY(N2),IERR)
      IF(IERR .EQ. 0) GO TO 1150
      CALL LNCNT(3)
      WRITE(6,1125)
 1125 FORMAT(/' IN RICNWT, A MATRIX WHICH IS NOT POSITIVE DEFINITE HAS B
     1EEN SUBMITTED TO DPOCO'/)
      RETURN
C 
 1150 CONTINUE
      NT = 1
      DO 1175 M1 = 1,NA(1)
         CALL DPOSL(DUMMY(N1),NR(1),NR(1),F(NT))
         NT = NT + NR(1)
 1175 CONTINUE
      IF( I .EQ. 1) GO TO 925
      CALL MAXEL(DUMMY,NA,ANORM1)
      CALL SUBT(P,NP,DUMMY,NA,DUMMY(N1),NA)
      CALL MAXEL(DUMMY(N1),NA,ANORM2)
      IF( ANORM1 .GT. 1) GO TO 1200
      IF( ANORM2/ANORM1 .LT. RICTCV ) GO TO 1300
      GO TO 1250
 1200 CONTINUE
      IF( ANORM2 .LT. RICTCV ) GO TO 1300
C 
 1250 CONTINUE
      IF( I .LE. 101) GO TO  925
      CALL LNCNT(3)
      WRITE(6,775)
      IOP(1) = 1
C 
 1300 CONTINUE
      IF(IOP(1) .EQ. 0 ) RETURN
      CALL LNCNT(4)
      WRITE(6,1350) I
 1350 FORMAT(//' FINAL VALUES OF P AND F AFTER',I5,' ITERATIONS TO CONVE
     1RGE'/)
      CALL PRNT(P,NP,' P  ',1)
      CALL PRNT(F,NF,' F  ',1)
C 
      RETURN
      END
      SUBROUTINE SAMPL(A,NA,B,NB,Q,NQ,R,NR,W,NW,T,IOP,DUMMY)
C 
C   REFERENCES:
C      Kwakernaak, Huibert; and Sivan, Raphael: Linear Optimal Control
C        Systems. John Wiley & Sons, Inc., c. 1972.
C      Dorato, Peter; and Levis, Alexander H.: Optimal Linear Regula-
C        tors:  The Discrete-Time Case. IEEE Trans. Autom. Control,
C        vol. AC-16, no. 6, Dec. 1971, pp. 613-620.
C      Levis, Alexander H.; Schlueter, Robert A.; and Athans, Michael:
C        On the Behaviour of Optimal Linear Sampled-Data Regulators.
C        Int. J. Control, vol. 13, no. 2, Feb. 1971, pp. 343-361.
C      Armstrong, Ernest S.; and Caglayan, Alper K.: An Algorithm for
C        the Weighting Matrices in the Sampled-Data Optimal Linear
C        Regulator Problem. NASA TN D-8372, 1976.
C 
C   Subroutines employed by SAMPL: ADD, EQUATE, EXPINT, EXPSER, LNCNT,
C      MAXEL, MULT, NORMS, PRNT, SCALE, TRANP
C   Subroutines employing SAMPL: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),Q(1),R(1),W(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NQ(2),NR(2),NW(2),IOP(2),NDUM(2)
C      COMMON/CONV/SUMCV,MAXSUM,RICTCV,SERCV
      COMMON/CONV/SUMCV,RICTCV,SERCV,MAXSUM
C 
      IF(  IOP(1) .EQ. 0 ) GO TO 100
C 
      IF(  IOP(2) .EQ. 0 ) GO TO 50
C 
      CALL LNCNT(5)
      WRITE(6,25)
   25 FORMAT(//' COMPUTATION OF WEIGHTING MATRICES FOR THE OPTIMAL SAMPL
     1ED-DATA REGULATOR PROBLEM'//)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      CALL LNCNT(3)
      WRITE(6,35)
   35 FORMAT(/' CONTINUOUS PERFORMANCE INDEX WEIGHTING MATRICES'/)
      CALL PRNT(Q,NQ,' Q  ',1)
      CALL PRNT(R,NR,' R  ',1)
      CALL LNCNT(3)
      WRITE(6,45) T
   45 FORMAT(/' SAMPLE TIME = ',D16.8/)
C 
      GO TO 100
C 
   50 CONTINUE
      CALL LNCNT(8)
      WRITE(6,75)
   75 FORMAT(//' COMPUTATION OF THE RECONSTRUCTIBILITY GRAMIAN'/' FOR TH
     1E (A,H) SYSTEM OVER THE INTERVAL (0,T) '/' THE MATRIX Q IS  ( H TR
     2ANSPOSE ) X H'//)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(Q,NQ,' Q  ',1)
      CALL LNCNT(3)
      WRITE(6,85) T
   85 FORMAT(/' T = ',D16.8/)
C 
  100 CONTINUE
C 
      N = NA(1)
      L = ( N**2)
      N1 = L + 1
      N2 = N1 + L
      TT  = T
C 
      IOPT = 1
      CALL NORMS(N,N,N,A,IOPT,ANORM)
      IOPT = 3
      CALL NORMS(N,N,N,A,IOPT,ROWA)
      IF( ANORM .GT. ROWA ) ANORM = ROWA
C 
      IF( ANORM .LE. 1.E-15 ) GO TO 900
C 
      TMAX = 1.0/ANORM
      K = 0
C 
  125 CONTINUE
      IF( TMAX - TT ) 150,150,200
  150 CONTINUE
      K = K + 1
      TT = T/( 2**K)
      IF( K - 1000 ) 125,800,800
C 
  200 CONTINUE
C 
      I = 0
      SC = TT
      CALL SCALE(A,NA,A,NA,TT)
      CALL SCALE(Q,NQ,Q,NQ,TT)
      CALL EQUATE(Q,NQ,DUMMY,NQ)
C 
      IF( IOP(2) .NE. 0 ) GO TO 500
C 
  225 CONTINUE
      II = I + 2
      I = I + 1
      F = 1.0/II
      CALL SCALE(A,NA,DUMMY(N1),NA,F)
      CALL MULT(DUMMY,NA,DUMMY(N1),NA,DUMMY(N2),NA)
      CALL TRANP(DUMMY(N2),NA,DUMMY(N1),NA)
      CALL ADD(DUMMY(N1),NA,DUMMY(N2),NA,DUMMY,NA)
C 
      CALL MAXEL(Q,NQ,TOT)
      CALL MAXEL(DUMMY,NA,DELT)
      IF( TOT .GT. 1.0 ) GO TO 250
      IF( DELT/TOT .LT. SERCV )  GO TO 300
      GO TO 275
  250 CONTINUE
      IF( DELT .LT. SERCV )  GO TO 300
  275 CONTINUE
      CALL ADD(Q,NQ,DUMMY,NA,Q,NQ)
      GO TO 225
C 
  300 CONTINUE
C 
      IF( K .EQ. 0 ) GO TO 400
      N3 = N2 + L
      G = 1.0
      IOPT = 0
      CALL EXPSER(A,NA,DUMMY,NA,G,IOPT,DUMMY(N1))
C 
  350 CONTINUE
      IF( K .EQ. 0 ) GO TO 400
      K = K-1
C 
      CALL TRANP(DUMMY,NA,DUMMY(N1),NA)
      CALL MULT(Q,NQ,DUMMY,NA,DUMMY(N2),NA)
      CALL MULT(DUMMY(N1),NA,DUMMY(N2),NA,DUMMY(N3),NA)
      CALL ADD(Q,NQ,DUMMY(N3),NA,Q,NQ)
      CALL MULT(DUMMY,NA,DUMMY,NA,DUMMY(N1),NA)
      CALL EQUATE(DUMMY(N1),NA,DUMMY,NA)
C 
      GO TO 350
C 
  400 CONTINUE
      S =  1.0/SC
      CALL SCALE(A,NA,A,NA,S)
C 
      IF( IOP(1) .EQ. 0 ) RETURN
      CALL PRNT(Q,NQ,'GRAM',1)
      RETURN
C 
  500 CONTINUE
      CALL SCALE(B,NB,B,NB,TT)
      N3 = N2 + L
      N4 = N3 + L
      N5 = N4 + L
      N6 = N5 + L
C 
  525 CONTINUE
      II = I + 2
      I = I + 1
      F = 1.0/II
      CALL SCALE(A,NA,DUMMY(N1),NA,F)
      CALL TRANP(DUMMY(N1),NA,DUMMY(N2),NA)
      CALL MULT(DUMMY,NA,DUMMY(N1),NA,DUMMY(N3),NA)
      CALL TRANP(DUMMY(N3),NA,DUMMY(N1),NA)
      CALL MULT(DUMMY,NA,B,NB,DUMMY(N5),NW)
      CALL ADD(DUMMY(N1),NA,DUMMY(N3),NA,DUMMY,NA)
      CALL SCALE(DUMMY(N5),NW,DUMMY(N1),NW,F)
      IF( I .NE. 1 ) GO TO 550
      CALL EQUATE(DUMMY(N1),NW,W,NW)
      CALL EQUATE(DUMMY(N1),NW,DUMMY(N6),NW)
      CALL ADD(Q,NQ,DUMMY,NQ,Q,NQ)
      GO TO 525
C 
  550 CONTINUE
      CALL MULT(DUMMY(N2),NA,DUMMY(N6),NW,DUMMY(N5),NW)
      CALL ADD(DUMMY(N5),NW,DUMMY(N1),NW,DUMMY(N1),NW)
      CALL TRANP(B,NB,DUMMY(N2),NDUM)
      CALL SCALE(DUMMY(N2),NDUM,DUMMY(N2),NDUM,F)
      CALL MULT(DUMMY(N2),NDUM,DUMMY(N6),NW,DUMMY(N3),NR)
      CALL TRANP(DUMMY(N3),NR,DUMMY(N5),NR)
      CALL ADD(DUMMY(N3),NR,DUMMY(N5),NR,DUMMY(N3),NR)
      CALL EQUATE(DUMMY(N1),NW,DUMMY(N6),NW)
      IF(  I .NE. 2 ) GO TO 575
      CALL ADD(Q,NQ,DUMMY,NQ,Q,NQ)
      CALL ADD(W,NW,DUMMY(N1),NW,W,NW)
      CALL EQUATE(DUMMY(N3),NR,DUMMY(N4),NR)
      GO TO 525
C 
  575 CONTINUE
      CALL MAXEL(Q,NQ,TOT)
      CALL MAXEL(DUMMY,NQ,DELT)
      IF( TOT .GT. 1.0 )  GO TO 580
      IF( DELT/TOT .LT. SERCV )  GO TO 585
      GO TO 595
C 
  580 CONTINUE
      IF( DELT .LT. SERCV )  GO TO 585
      GO TO 595
C 
  585 CONTINUE
      CALL MAXEL(DUMMY(N4),NR,TOT)
      CALL MAXEL(DUMMY(N3),NR,DELT)
      IF( TOT .GT. 1.0 )  GO TO 590
      IF( DELT/TOT .LT. SERCV )  GO TO 600
      GO TO 595
C 
  590 CONTINUE
      IF( DELT .LT. SERCV )  GO TO 600
C 
  595 CONTINUE
      CALL ADD (Q,NQ,DUMMY,NQ,Q,NQ)
      CALL ADD(W,NW,DUMMY(N1),NW,W,NW)
      CALL ADD(DUMMY(N4),NR,DUMMY(N3),NR,DUMMY(N4),NR)
      GO TO 525
C 
  600 CONTINUE
      IF( K .EQ. 0 ) GO TO 700
      G = 1.0
      IOPT = 0
      CALL EXPINT(A,NA,DUMMY,NA,DUMMY(N1),NA,G,IOPT,DUMMY(N2))
      CALL MULT(DUMMY(N1),NA,B,NB,DUMMY(N2),NB)
      CALL EQUATE(DUMMY(N2),NB,DUMMY(N1),NB)
C 
  650 CONTINUE
      IF( K .EQ. 0 ) GO TO 700
      K = K - 1
      CALL MULT(Q,NQ,DUMMY,NA,DUMMY(N2),NA)
      CALL TRANP(DUMMY,NA,DUMMY(N3),NA)
      CALL MULT(DUMMY(N3),NA,DUMMY(N2),NA,DUMMY(N5),NA)
      CALL MULT(Q,NQ,DUMMY(N1),NB,DUMMY(N2),NB)
      CALL ADD(Q,NQ,DUMMY(N5),NA,Q,NQ)
      CALL MULT(DUMMY(N3),NA,DUMMY(N2),NB,DUMMY(N5),NB)
      CALL MULT(DUMMY(N3),NA,W,NW,DUMMY(N6),NW)
      CALL ADD(DUMMY(N5),NW,DUMMY(N6),NW,DUMMY(N5),NW)
      CALL TRANP(DUMMY(N1),NB,DUMMY(N6),NDUM)
      CALL MULT(DUMMY(N6),NDUM,W,NW,DUMMY(N3),NR)
      CALL ADD(W,NW,DUMMY(N5),NW,W,NW)
      CALL MULT(DUMMY(N6),NDUM,DUMMY(N2),NB,DUMMY(N5),NR)
      CALL ADD(DUMMY(N5),NR,DUMMY(N3),NR,DUMMY(N5),NR)
      CALL TRANP(DUMMY(N3),NR,DUMMY(N6),NR)
      CALL ADD(DUMMY(N5),NR,DUMMY(N6),NR,DUMMY(N6),NR)
      CALL SCALE(DUMMY(N4),NR,DUMMY(N4),NR,2.0D0)
      CALL ADD(DUMMY(N6),NR,DUMMY(N4),NR,DUMMY(N4),NR)
      CALL MULT(DUMMY,NA,DUMMY(N1),NB,DUMMY(N3),NB)
      CALL ADD(DUMMY(N3),NB,DUMMY(N1),NB,DUMMY(N1),NB)
      CALL MULT(DUMMY,NA,DUMMY,NA,DUMMY(N3),NA)
      CALL EQUATE(DUMMY(N3),NA,DUMMY,NA)
      GO TO 650
C 
  700 CONTINUE
      CALL SCALE(R,NR,R,NR,T)
      CALL ADD(R,NR,DUMMY(N4),NR,R,NR)
      CALL SCALE(W,NW,W,NW,2.0D0)
      S = 1.0/SC
      CALL SCALE(A,NA,A,NA,S)
      CALL SCALE(B,NB,B,NB,S)
      IF( IOP(1) .EQ. 0 ) RETURN
C 
      CALL LNCNT(3)
      WRITE(6,750)
  750 FORMAT(/' DISCRETE PERFORMANCE INDEX WEIGHTING MATRICES'/)
      CALL PRNT(Q,NQ,' Q  ',1)
      CALL PRNT(W,NW,' W  ',1)
      CALL PRNT(R,NR,' R  ',1)
      RETURN
C 
  800 CONTINUE
      CALL LNCNT(1)
      WRITE(6,850)
  850 FORMAT(' ERROR IN SAMPL , K = 1000')
      RETURN
C 
  900 CONTINUE
      CALL SCALE(Q,NQ,Q,NQ,T)
      IF( IOP(2) .NE. 0 ) GO TO 925
      IF( IOP(1) .NE. 0 ) CALL PRNT(Q,NQ,'GRAM',1)
      RETURN
C 
  925 CONTINUE
      CALL MULT(Q,NQ,B,NB,W,NW)
      CALL SCALE(W,NW,W,NW,T)
      CALL TRANP(B,NB,DUMMY,NDUM)
      CALL MULT(DUMMY,NDUM,W,NW,DUMMY(N1),NR)
      TT = T/3.
      CALL SCALE(DUMMY(N1),NR,DUMMY,NR,TT)
      CALL SCALE(R,NR,R,NR,T)
      CALL ADD(R,NR,DUMMY,NR,R,NR)
      IF( IOP(1) .EQ. 0 ) RETURN
      CALL LNCNT(3)
      WRITE(6,750)
      CALL PRNT(Q,NQ,' Q  ',1)
      CALL PRNT(W,NW,' W  ',1)
      CALL PRNT(R,NR,' R  ',1)
      RETURN
C 
      END
      SUBROUTINE SCALE (A, NA, B, NB, S)
C 
C   PURPOSE:
C      Perform scalar multiplication on a given matrix.
C 
C   Subroutines employed by SCALE: LNCNT
C   Subroutines employing SCALE: ASYREG, BILIN, CNTREG, CSTAB, DSTAB,
C      EXMDFL, EXPINT, EXPSER, FACTOR, IMMDFL, PREFIL, RICNWT, SAMPL,
C      TRNSIT, VARANC
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),NA(2),NB(2)
      NB(1) = NA(1)
      NB(2) =NA(2)
      L = NA(1)*NA(2)
      IF( NA(1) .LT. 1 .OR. L .LT. 1 ) GO TO 999
      DO 300 I=1,L
  300 B(I)=A(I)*S
 1000 RETURN
  999 CALL LNCNT(1)
      WRITE  (6,50) NA
   50 FORMAT  (' DIMENSION ERROR IN SCALE   NA=',2I6)
      RETURN
      END
      SUBROUTINE SCHUR(H,U,NN,NH,NU,EPS,FAIL)
C 
C   PURPOSE:
C      Reduce an n x n upper Hessenberg matrix H to real Schur form.
C      Computation is by the QR algorithm with implicit origin shifts.
C      The product of the transformations used in the reduction is
C      accumulated.
C 
C   REFERENCES:
C      Martin, R.S.; Peters, G.; and Wilkinson, J.H.: The QR Algorithm
C        for Real Hessenberg Marices.  Numer. Math., Bd. 14, Heft 3,
C        1970, pp. 219-231.
C      Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
C        the Matrix Equation AX + XB = C.  Commun. ACM, vol. 15, no. 9,
C        Sept. 1972, pp. 820-826.
C 
C   Subroutines employed by SCHUR: None
C   Subroutines employing SCHUR: ATXPXA, AXPXB
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8
     1H(NH,1),U(NU,1),EPS,HN,RSUM,TEST,P,Q,R,S,W,X,Y,Z
      INTEGER
     1NN,NA,NH,FAIL,I,ITS,J,JL,K,L,LL,M,MM,M2,M3,N
      LOGICAL
     1LAST
      N = NN
      HN = 0.
      DO 20 I=1,N
        JL = MAX0(1,I-1)
        RSUM = 0.
      DO 10 J=JL,N
          RSUM = RSUM + DABS(H(I,J))
   10   CONTINUE
        HN = DMAX1(HN,RSUM)
   20 CONTINUE
      TEST = EPS*HN
      IF(HN .EQ. 0.) GO TO 230
   30 IF(N .LE. 1) GO TO 230
      ITS = 0
      NA = N-1
      NM2 = N-2
   40 DO 50 LL=2,N
      L = N-LL+2
        IF(DABS(H(L,L-1)) .LE. TEST) GO TO 60
   50 CONTINUE
      L = 1
      GO TO 70
   60 H(L,L-1) = 0.
   70 IF(L .LT. NA) GO TO 72
      N = L-1
      GO TO 30
   72 X = H(N,N)/HN
      Y = H(NA,NA)/HN
      R = (H(N,NA)/HN)*(H(NA,N)/HN)
      IF(ITS .LT. 30) GO TO 75
      FAIL = N
      RETURN
   75 IF(ITS.EQ.10 .OR. ITS.EQ.20) GO TO 80
      S = X + Y
      Y = X*Y - R
      GO TO 90
   80 Y = (DABS(H(N,NA)) + DABS(H(NA,NM2)))/HN
      S = 1.5*Y
      Y = Y**2
   90 ITS = ITS + 1
      DO 100 MM=L,NM2
        M = NM2-MM+L
        X = H(M,M)/HN
        R = H(M+1,M)/HN
        Z = H(M+1,M+1)/HN
        P = X*(X-S) + Y + R*(H(M,M+1)/HN)
        Q = R*(X+Z-S)
        R = R*(H(M+2,M+1)/HN)
        W = DABS(P) + DABS(Q) + DABS(R)
        P = P/W
        Q = Q/W
        R = R/W
        IF(M .EQ. L) GO TO 110
      IF(DABS(H(M,M-1))*(DABS(Q)+DABS(R)) .LE. DABS(P)*TEST)
     1GO TO 110
  100 CONTINUE
  110 M2 = M+2
      M3 = M+3
      DO 120 I=M2,N
        H(I,I-2) = 0.
  120 CONTINUE
      IF(M3 .GT. N) GO TO 140
      DO 130 I=M3,N
        H(I,I-3) = 0.
  130 CONTINUE
  140 DO 220 K=M,NA
        LAST = K.EQ.NA
        IF(K .EQ. M) GO TO 150
        P = H(K,K-1)
        Q = H(K+1,K-1)
        R = 0.
        IF(.NOT.LAST) R = H(K+2,K-1)
        X = DABS(P) + DABS(Q) + DABS(R)
        IF(X .EQ. 0.) GO TO 220
        P = P/X
        Q = Q/X
        R = R/X
  150   S = DSQRT(P**2 + Q**2 + R**2)
        IF(P .LT. 0.) S = -S
        IF(K .NE. M) H(K,K-1) = -S*X
        IF(K.EQ.M .AND. L.NE.M) H(K,K-1) = -H(K,K-1)
        P = P + S
        X = P/S
        Y = Q/S
        Z = R/S
        Q = Q/P
        R = R/P
        DO 170 J=K,NN
          P = H(K,J) + Q*H(K+1,J)
          IF(LAST) GO TO 160
          P = P + R*H(K+2,J)
          H(K+2,J) = H(K+2,J) - P*Z
  160     H(K+1,J) = H(K+1,J) - P*Y
          H(K,J) = H(K,J) - P*X
  170   CONTINUE
        J = MIN0(K+3,N)
        DO 190 I=1,J
          P = X*H(I,K) + Y*H(I,K+1)
          IF(LAST) GO TO 180
          P = P + Z*H(I,K+2)
          H(I,K+2) = H(I,K+2) - P*R
  180     H(I,K+1) = H(I,K+1) - P*Q
          H(I,K) = H(I,K) - P
  190   CONTINUE
        DO 210 I=1,NN
          P = X*U(I,K) + Y*U(I,K+1)
          IF(LAST) GO TO 200
          P = P + Z*U(I,K+2)
          U(I,K+2) = U(I,K+2) - P*R
  200     U(I,K+1) = U(I,K+1) - P*Q
          U(I,K) = U(I,K) - P
  210   CONTINUE
  220 CONTINUE
      GO TO 40
  230 FAIL = 0
      RETURN
      END
       SUBROUTINE SHRSLV(A,B,C,M,N,NA,NB,NC)
C 
C   PURPOSE:
C      Solve the real matrix equation AX + XB = C, where A is an m x m
C      matrix in lower real Schur form and B is an n x n matrix in
C      upper real Schur form.
C 
C   REFERENCES:
C      Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
C        the Matrix Equation AX + XB = C. Commun. ACM, vol. 15, no. 9,
C        Sept. 1972, pp. 820-826.
C 
C   Subroutines employed by SHRSLV: SYSSLV
C   Subroutines employing SHRSLV: AXPXB
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8
     1A(NA,1),B(NB,1),C(NC,1),T,P
      INTEGER
     1M,N,NA,NB,NC,K,KM1,DK,KK,L,LM1,DL,LL,I,IB,J,JA,NSYS
      COMMON/SLVBLK/T(5,5),P(5),NSYS
      L = 1
   10   LM1 = L-1
        DL = 1
        IF(L .EQ. N) GO TO 15
        IF(B(L+1,L) .NE. 0.) DL = 2
   15   LL = L+DL-1
        IF(L .EQ. 1) GO TO 30
        DO 20 J=L,LL
          DO 20 I=1,M
            DO 20 IB=1,LM1
              C(I,J) = C(I,J) - C(I,IB)*B(IB,J)
   20   CONTINUE
   30   K = 1
   40     KM1 = K-1
          DK = 1
          IF(K .EQ. M) GO TO 45
          IF(A(K,K+1) .NE. 0.) DK = 2
   45     KK = K+DK-1
          IF(K .EQ. 1) GO TO 60
          DO 50 I=K,KK
            DO 50 J=L,LL
              DO 50 JA=1,KM1
                C(I,J) = C(I,J) - A(I,JA)*C(JA,J)
  50      CONTINUE
  60      IF(DL .EQ. 2) GO TO 80
          IF(DK .EQ. 2) GO TO 70
          T(1,1) = A(K,K) + B(L,L)
          IF(T(1,1) .EQ. 0.) STOP
          C(K,L) = C(K,L)/T(1,1)
          GO TO 100
  70      T(1,1) = A(K,K) + B(L,L)
          T(1,2) = A(K,KK)
          T(2,1) = A(KK,K)
          T(2,2) = A(KK,KK) + B(L,L)
          P(1) = C(K,L)
          P(2) = C(KK,L)
          NSYS = 2
          CALL SYSSLV
          C(K,L) = P(1)
          C(KK,L) = P(2)
          GO TO 100
  80      IF(DK .EQ. 2) GO TO 90
          T(1,1) = A(K,K) + B(L,L)
          T(1,2) = B(LL,L)
          T(2,1) = B(L,LL)
          T(2,2) = A(K,K) + B(LL,LL)
          P(1) = C(K,L)
          P(2) = C(K,LL)
          NSYS = 2
          CALL SYSSLV
          C(K,L) = P(1)
          C(K,LL) = P(2)
          GO TO 100
  90      T(1,1) = A(K,K) + B(L,L)
          T(1,2) = A(K,KK)
          T(1,3) = B(LL,L)
          T(1,4) = 0.
          T(2,1) = A(KK,K)
          T(2,2) = A(KK,KK) + B(L,L)
          T(2,3) = 0.
          T(2,4) = T(1,3)
          T(3,1) = B(L,LL)
          T(3,2) = 0.
          T(3,3) = A(K,K) + B(LL,LL)
          T(3,4) = T(1,2)
          T(4,1) = 0.
          T(4,2) = T(3,1)
          T(4,3) = T(2,1)
          T(4,4) = A(KK,KK) + B(LL,LL)
          P(1) = C(K,L)
          P(2) = C(KK,L)
          P(3) = C(K,LL)
          P(4) = C(KK,LL)
          NSYS = 4
          CALL SYSSLV
          C(K,L) = P(1)
          C(KK,L) = P(2)
          C(K,LL) = P(3)
          C(KK,LL) = P(4)
  100   K = K + DK
        IF(K .LE. M) GO TO 40
      L = L + DL
      IF(L .LE. N) GO TO 10
      RETURN
      END
      SUBROUTINE SNVDEC(IOP,MD,ND,M,N,A,NOS,B,IAC,ZTEST,Q,V,IRANK,APLUS,
     1 IERR)
C
      INTEGER IAC,IERR,IOP,IRANK,M,MD,N,ND,NOS
      DOUBLE PRECISION A(MD,N),APLUS(ND,N),B(MD,NOS),Q(N),V(ND,N)
      DOUBLE PRECISION ZTEST
C
C   PURPOSE:
C     Compute the singular-value decomposition of a real m x n matrix A.
C
C   REFERENCES:
C     Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Compu-
C        tation. Volume II - Linear Algebra. Springer-Verlag, 1971.
C        See pp.135-151 for computational procedure.
C
C   Subroutines employed by SNVDEC: DSVDC
C   Fortran-77 functions employed by SNVDEC: SQRT,MIN,POWER [X**Y]
C   Subroutines employing SNVDEC: CTROL, CSTAB, DSTAB, FACTOR
C
C   Modified to LINPACK form by: John L. Tietze, Boeing Aerospace Co.
C                Latest Version: November 10, 1981
C
C
C   ON ENTRY:
C
C     IOP      INTEGER
C              OPTION CODE WITH FOLLOWING MEANINGS:
C                = 1  RETURN RANK AND SINGULAR VALUES OF A
C                = 2  IOP= 1 INFORMATION PLUS U AND V MATRICES
C                = 3  IOP= 2 INFORMATION PLUS LEAST SQUARES
C                     SOLUTION OF A*X= B
C                = 4  IOP= 2 INFORMATION PLUS PSEUDOINVERSE OF A
C                = 5  IOP= 4 INFORMATION PLUS LEAST SQUARES
C                     SOLUTION OF A*X= B
C
C     MD       INTEGER
C              LEADING DIMENSION OF A AND B MATRICES
C
C     ND       INTEGER
C              LEADING DIMENSION OF V AND APLUS MATRICES
C
C     M        INTEGER
C              NUMBER OF ROWS IN A MATRIX
C
C     N        INTEGER
C              NUMBER OF COLUMNS IN A MATRIX
C
C     A        DOUBLE PRECISION (MD,N)
C              INPUT MATRIX, DESTROYED BY COMPUTATION
C
C     NOS      INTEGER
C              NUMBER OF COLUMNS IN B MATRIX
C
C     B        DOUBLE PRECISION (MD,NOS)
C              B CONTAINS THE NOS RIGHT HAND SIDES TO BE
C              SOLVED FOR IF IOP= 3 OR IOP= 5. B MAY BE
C              A DUMMY VARIABLE IF NOT REFERENCED.
C
C     IAC      INTEGER
C              NUMBER OF DECIMAL DIGITS OF ACCURACY IN THE
C              ELEMENTS OF THE A MATRIX. THIS PARAMETER
C              IS USED IN THE RANK DETERMINATION TEST.
C
C  ON RETURN:
C
C     A        DOUBLE PRECISION (MD,N)
C              A CONTAINS THE ORTHOGONAL MATRIX U WHEN IOP>1
C
C     B        DOUBLE PRECISION (MD,NOS)
C              B CONTAINS THE LEAST SQUARES SOLUTION WHEN
C              IOP= 3 OR IOP= 5
C
C     ZTEST    DOUBLE PRECISION
C              ZERO CRITERION FOR RANK TEST, COMPUTED AS
C              NORM(A)*10**-IAC. NORM(A) IS THE EUCLIDEAN
C              MATRIX NORM. NORM(A)= 1 IF N= 1.
C
C     Q        DOUBLE PRECISION (MIN(M+1,N))
C              SINGULAR VALUES OF A MATRIX IN DESCENDING ORDER.
C
C     V        DOUBLE PRECISION (ND,N)
C              LEFT SINGULAR VECTORS OF A MATRIX WHEN IOP>1.
C
C     IRANK    INTEGER
C              RANK OF A MATRIX DETERMINED FROM THE SINGULAR
C              VALUES AND ZTEST.
C
C     APLUS    DOUBLE PRECISION (ND,N)
C              PSEUDOINVERSE OF A MATRIX IF IOP= 4 OR IOP= 5.
C              APLUS MAY BE A DUMMY VARIABLE IF NOT REFERENCED.
C
C     IERR     INTEGER
C               ERROR INDICATOR:
C                  = 0     A NORMAL RETURN
C                  = K > 0 SINGULAR VALUES FROM K+1 ON HAVE NOT BEEN
C                          FOUND BECAUSE THE QR ITERATION TOOK MORE
C                          THAN 30 PASSES.
C                  = -1    USING THE GIVEN IAC, THE A MATRIX IS CLOSE
C                          TO A MATRIX OF LOWER RANK.
C                  = -2    MATRIX HAS ZERO RANK
C                  = -3    HAVE EXCEEDED THE INTERNAL STORAGE CAPACITY
C                          OF THE E VECTOR.
C
C  LOCAL VARIABLES:
      DOUBLE PRECISION E(150)
      DOUBLE PRECISION SUM
      INTEGER I,J,JOB,K,MINMN
C
C  TEST FOR  SCALAR OR VECTOR A
C
      IF( N .GE. 2 ) GO TO 3000
C
          IERR = 0
          ZTEST = 10.**(-IAC)
          SUM = 0.0
          DO 1000 I=1,M
              SUM = SUM + A(I,1)*A(I,1)
 1000     CONTINUE
          SUM = SQRT(SUM)
          IRANK = 0
          IF( SUM .GT. ZTEST ) IRANK = 1
          Q(1) = SUM
C
          IF( IOP .EQ. 1) GO TO 3500
          V(1,1) = 1.0
          IF( IRANK .EQ. 0 ) GO TO 1200
              DO 1100 I =1,M
                  A(I,1) = A(I,1)/SUM
 1100         CONTINUE
              GO TO 1300
 1200     CONTINUE
          A(1,1) = 1.0
 1300     CONTINUE
C
          IF( IOP .EQ. 2 ) GO TO 3500
          IF( IOP .EQ. 4 ) GO TO 1850
              IF( IRANK .EQ. 0 ) GO TO 1600
                  DO 1500 J = 1,NOS
                      SUM = 0.0
                      DO 1400 I = 1,M
                          SUM = SUM + A(I,1)*B(I,J)/SUM
 1400                 CONTINUE
                      B(1,J) = SUM
 1500             CONTINUE
                  GO TO 1800
 1600         CONTINUE
              DO 1700 J =1,NOS
                  B(1,J) = 0.0
 1700         CONTINUE
 1800         CONTINUE
C
              IF( IOP .EQ. 3 ) GO TO 3500
 1850     CONTINUE
          IF( IRANK .EQ. 0 ) GO TO 2000
              DO 1900 I =1,M
                  APLUS(1,I) = A(I,1)/SUM
 1900         CONTINUE
              GO TO 3500
 2000     CONTINUE
          DO 2100 I=1,M
              APLUS(1,I) = 0.0
 2100     CONTINUE
          GO TO 3500
C
C
 3000 CONTINUE
C
C  COMPUTE THE E-NORM OF MATRIX A AS ZERO TEST FOR SINGULAR VALUES
C
      SUM=0.0D0
      DO 3010 I=1,M
        DO 3010 J=1,N
 3010     SUM = SUM + A(I,J)**2
      ZTEST = SQRT(SUM)*10.0**(-IAC)
C
C   EXECUTIVE CALL TO LINPACK SINGULAR VALUE SUBROUTINE HERE
C   CHECK FOR PROBLEM BIGGER THAN LOCAL ARRAY FIRST
C
      MINMN= MIN (M,N)
      IF (M + MINMN.GT.150) GO TO 4000
C
C  MAY WANT TO MODIFIY THIS FOR LARGE SYSTEMS WHERE M>>N, SEE LINPACK
C  WRITE UP ON SINGULAR VALUE FACTORIZATION.
      JOB= 0
      IF (IOP.GT.1) JOB= 11
C****
      CALL DSVDC (A,MD,M,N,Q,E(M+1),A,MD,V,ND,E,JOB,IERR)
C****
      IF (IERR.GT.0) GO TO 3500
C
C  DETERMINE RANK
      DO 3020 I= 1,MINMN
          IRANK= MINMN+1-I
          IF (Q(IRANK).GT.ZTEST) GO TO 3030
 3020 CONTINUE
      IRANK= 0
      IERR= -2
      GO TO 3500
 3030 CONTINUE
      IF (ZTEST/Q(IRANK).GT..0625D0) IERR= -1
      IF (IOP.LT.3) GO TO 3500
      IF (IOP.EQ.4) GO TO 3100
C
C  GET LEAST-SQUARES SOLUTION, E IS TEMPORARY STORAGE AND
C  A CONTAINS THE LEFT SINGULAR VECTORS
      DO 3080  K=1,NOS
         DO 3050  J=1,IRANK
            SUM=0.0
            DO 3040  I=1,M
 3040          SUM =SUM + A(I,J)*B(I,K)
 3050    E(J)= SUM/Q(J)
C
         DO 3070  J=1,N
            SUM=0.0
            DO 3060  I=1,IRANK
 3060         SUM =SUM  + V(J,I)*E(I)
 3070    B(J,K)=SUM
 3080 CONTINUE
      IF (IOP.LT.4) GO TO 3500
C
C  COMPUTE THE PSEUDOINVERSE  -  IS THIS REALLY NEEDED ???
C
 3100 DO 3130  J=1,M
         DO 3120  I=1,N
            SUM=0.0
            DO 3110  K=1,IRANK
 3110          SUM =SUM + V(I,K)*A(J,K)/Q(K)
 3120       APLUS(I,J)= SUM
 3130 CONTINUE
C
C------------------------------------------- END -----------------------
C
 3500 CONTINUE
      RETURN
C
C  ERROR EXIT
 4000 J= M + MIN(M,N)
      WRITE (6,4010) J
 4010 FORMAT (//15X,'**** M + MIN(M,N) =',I4,' IN SNVDEC.'
     1 /15X,'THIS EXCEEDS THE INTERNAL SPECIFICATION ON',
     2 ' THE E VECTOR OF 150.'/15X,'SNVDEC WILL BE BYPASSED.')
      IERR= -2
      RETURN
C
      END
      SUBROUTINE SUBT(A,NA,B,NB,C,NC)
C 
C   PURPOSE:
C      Perform matrix subtraction C = A - B for given matrices A and B.
C 
C   Subroutines employed by SUBT: LNCNT
C   Subroutines employing SUBT: ASYREG, CNTREG, CSTAB, DISREG, DSTAB,
C      EXMDFL, IMMDFL, PREFIL, RICNWT, TRNSIT
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),NA(2),NB(2),NC(2)
      IF((NA(1).NE.NB(1)).OR.(NA(2).NE.NB(2))) GO TO 999
      NC(1)=NA(1)
      NC(2)=NA(2)
      L=NA(1)*NA(2)
      IF( NA(1) .LT. 1 .OR. L .LT. 1 ) GO TO 999
      DO 300 I=1,L
  300 C(I)=A(I)-B(I)
      GO TO 1000
  999 CALL LNCNT (1)
      WRITE(6,50) NA,NB
   50 FORMAT  (' DIMENSION ERROR IN SUBT    NA=',2I6,5X,'NB=',2I6)
 1000 RETURN
      END
      SUBROUTINE SUM(A,NA,B,NB,C,NC,IOP,SYM,DUMMY)
C 
C   PURPOSE:
C      Evaluate until convergence the matrix series
C         X = summation(i=0,infinity) (A**i)B(C**i)
C      where A and C are n x n and m x m real constant matrices.  The
C      matrix B is real constant and n x m.
C 
C   REFERENCES:
C      Armstrong, E.S.: Digital Explicit Model Following with Unstable
C        Model Dynamics. AIAA Paper No. 74-888, AIAA Mechanics and Con-
C        trol of Flight Conference, Aug. 1974.
C 
C   Subroutines employed by SUM: ADD, EQUATE, LNCNT, MAXEL, MULT, PRNT,
C      TRANP
C   Subroutines employing SUM: BILIN, EXMDFL, RICNWT, VARANC
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),DUMMY(1)
      DIMENSION NA(2),NB(2),NC(2)
      LOGICAL SYM
C      COMMON/CONV/SUMCV,MAXSUM,RICTCV,SERCV
      COMMON/CONV/SUMCV,RICTCV,SERCV,MAXSUM
C 
      IF( IOP  .EQ. 0  ) GO TO 100
      WRITE(6,50)
   50 FORMAT(//' LINEAR EQUATION SOLVER    X = AXC + B ')
      CALL PRNT(A,NA,' A  ',1)
      IF( SYM ) GO TO 75
      CALL PRNT(C,NC,' C  ',1)
      GO TO 85
   75 CONTINUE
      WRITE(6,80)
   80 FORMAT(/ ' C = A TRANSPOSE '/)
   85 CONTINUE
      CALL PRNT(B,NB,' B  ',1)
C 
  100 CONTINUE
      N1 = 1 + NA(1)*NC(1)
      I=1
  200 CONTINUE
      CALL MULT(A,NA,B,NB,DUMMY,NB)
      CALL MULT(DUMMY,NB,C,NC,DUMMY(N1),NB)
      CALL MAXEL(B,NB,WNS)
      CALL MAXEL(DUMMY(N1),NB,WNDX)
      IF(WNS .GE. 1.) GO TO 225
      IF( WNDX/WNS .LT. SUMCV ) GO TO 300
      GO TO 235
  225 IF( WNDX .LT. SUMCV ) GO TO 300
  235 CONTINUE
      CALL ADD(B,NB,DUMMY(N1),NB,B,NB)
      CALL MULT(A,NA,A,NA,DUMMY,NA)
      CALL EQUATE(DUMMY,NA,A,NA)
      IF( SYM ) GO TO 245
      CALL MULT(C,NC,C,NC,DUMMY,NC)
      CALL EQUATE(DUMMY,NC,C,NC)
      GO TO  250
  245 CONTINUE
      CALL TRANP(A,NA,C,NC)
  250 CONTINUE
      I=I+1
      IF( I .LE. MAXSUM ) GO TO 200
      CALL LNCNT(3)
      WRITE(6,275) MAXSUM
  275 FORMAT(//' IN SUM, THE SEQUENCE OF PARTIAL SUMS HAS EXCEEDED STAGE
     1 ',I5,' WITHOUT CONVERGENCE')
  300 CONTINUE
      IF(IOP .EQ. 0) RETURN
      CALL PRNT(B,NB,' X  ',1)
      RETURN
      END
      SUBROUTINE SYMPDS (MAXN,N,A,NRHS,B,IOPT,IFAC,DETERM,ISCALE,P,IERR)
C 
C   PURPOSE:
C      Solve the matrix equation, AX=B, where A is a symmetric positive
C      definite matrix and B is a matrix of constant vectors.
C 
C   REFERENCES:
C      Wilkinson, J.H.; and Reinsch, C.: Handbook for Automatic Compu-
C        tation. Volume II - Linear Algebra. Springer-Verlag, 1971.
C 
C   Subroutines employed by SYMPDS: None
C   Subroutines employing SYMPDS: RICNWT
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(MAXN,N),B(MAXN,NRHS),P(N)
C 
      DATA R1,R2/1,10/
C 
C             TEST FOR A SCALAR MATRIX (IF COEFFICIENT MATRIX IS A
C             SCALAR-- SOLVE  AND COMPUTE DETERMINANT IF DESIRED)
      IERR = 0
      NM1 = N-1
      IF (NM1 .GT. 0) GO TO 20
C 
      IF( A(1,1) .LE. 0.0 )  IERR = 1
      ISCALE = 0
      DETERM = A(1,1)
      P(1) = 1.0/A(1,1)
      DO 10 J=1,NRHS
      B(1,J) = B(1,J)/DETERM
   10 CONTINUE
      RETURN
C 
C             TEST TO DETERMINE IF CHOLESKY DECOMPOSITION OF COEFFICIENT
C             MATRIX IS DESIRED
   20 IF (IFAC .EQ. 1) GO TO 160
C 
C             INITIALIZE DETERMINANT EVALUATION PARAMETERS
      DETERM=1.0
      ISCALE=0
C 
C              "LOOP" TO PERFORM CHOLESKY DECOMPOSTION ON THE COEF-
C             FICIENT MATRIX A (I.E. MATRIX A WILL BE DECOMPOSED INTO
C             THE PRODUCT OF A UNIT LOWER TRIANGULAR MATRIX (L), A
C             DIAGONAL MATRIX (D), AND THE TRANSPOSE OF L (LTRANSPOSE).)
C 
   30 DO 150 I=1,N
      IM1 = I-1
C 
      DO 150 J=1,I
      X = A(J,I)
C 
C             DETERMINE IF ELEMENTS ARE ABOVE OR BELOW THE DIAGONAL
      IF (I .GT. J) GO TO 110
C 
C             USING THE DIAGONAL ELEMENTS OF MATRIX A, THIS SECTION
C             COMPUTES DIAGONAL MATRIX AND DETERMINES IF MATRIX A IS
C             SYMMETRIC POSITIVE DEFINITE
      IF (IM1 .EQ. 0) GO TO 50
      DO 40 K=1,IM1
      Y = A(I,K)
      A(I,K) = Y*P(K)
      X = X - Y*A(I,K)
   40 CONTINUE
C 
C             TEST IF COEFFICIENT MATRIX IS POSITIVE DEFINITE
   50 IF( X .LE. 0.0 )  IERR = 1
C 
C             COMPUTE INVERSE OF DIAGONAL MATRIX D**-1 = 1/P
      P(I) = 1.0 / X
C 
C             TEST TO SEE IF DETERMINANT IS TO BE EVALUATED
      IF (IOPT .EQ. 0) GO TO 150
C 
C 
C             SCALE THE DETERMINANT (COMPUTE THE DETERMINANT EVALUATION
C             PARAMETERS DETERM AND ISCALE)
      DETERM = DETERM * X
   60 IF(DABS(DETERM).LT.R1) GO TO 70
      DETERM = DETERM*R2
      ISCALE = ISCALE+1
      GO TO 60
   70 IF(DABS(DETERM).GT.R2) GO TO 150
      DETERM = DETERM*R1
      ISCALE = ISCALE-1
      GO TO 70
C 
C 
C             USING THE LOWER TRIANGULAR ELEMENTS OF MATRIX A, THIS
C             SECTION COMPUTES THE UNIT LOWER TRIANGULAR MATRIX
  110 JM1 = J-1
      IF (JM1 .EQ. 0) GO TO 140
      DO 120 K=1,JM1
      X = X - A(I,K)*A(J,K)
  120 CONTINUE
C 
  140 A(I,J) = X
C 
  150 CONTINUE
C 
C             SECTION TO APPLY BACK SUBSTITUTION TO SOLVE L*Y = B FOR
C             UNIT LOWER TRIANGULAR MATRIX AND CONSTANT COLUMN VECTOR B
C 
  160 IF( IFAC .EQ. 2 )  RETURN
      DO 180 I=2,N
      IM1=I-1
C 
      DO 180 J=1,NRHS
      X = B(I,J)
C 
      DO 170 K=1,IM1
      X = X - A(I,K)*B(K,J)
  170 CONTINUE
C 
      B(I,J) = X
  180 CONTINUE
C 
C             SECTION TO SOLVE (LTRANSROSE)*X = (D**-1)*Y FOR TRANSPOSE
C             OF UNIT LOWER TRIANGULAR MATRIX AND INVERSE OF DIAGONAL
C             MATRIX
C 
      Y = P(N)
      DO 190 J=1,NRHS
      B(N,J) = B(N,J)*Y
  190 CONTINUE
C 
  200 I = NM1+1
      Y = P(NM1)
C 
      DO 220 J=1,NRHS
      X = B(NM1,J)*Y
C 
      DO 210 K=I,N
      X = X - A(K,NM1)*B(K,J)
  210 CONTINUE
C 
      B(NM1,J) = X
  220 CONTINUE
C 
C 
C             TEST TO DETERMINE IF SOLUTIONS HAVE BEEN DETERMINED FOR
C             ALL COLUMN VECTORS
      NM1 = NM1-1
      IF (NM1 .GT. 0) GO TO 200
C 
      RETURN
      END
      SUBROUTINE SYMSLV(A,C,N,NA,NC)
C 
C   PURPOSE:
C      Solve the real matrix equation A'X + XA = C, where C=C' and A
C      is n x n and in upper real Schur form.
C 
C   REFERENCES:
C      Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
C        the Matrix Equation AX + XB = C.  Commun. ACM, vol. 15, no. 9,
C        Sept. 1972, pp. 820-826.
C 
C   Subroutines employed by SYMSLV: SYSSLV
C   Subroutines employing SYMSLV: ATXPXA
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8
     1A(NA,1),C(NC,1),T,P
      INTEGER
     1N,NA,NC,K,KK,DK,KM1,L,LL,DL,LDL,I,IA,J,NSYS
      COMMON/SLVBLK/T(5,5),P(5),NSYS
      L = 1
   10   DL = 1
        IF(L .EQ. N) GO TO 20
        IF(A(L+1,L) .NE. 0.) DL = 2
   20   LL = L+DL-1
        K = L
   30     KM1 = K-1
          DK = 1
          IF(K .EQ. N) GO TO 35
          IF(A(K+1,K) .NE. 0.) DK = 2
   35     KK = K+DK-1
          IF(K .EQ. L) GO TO 45
          DO 40 I=K,KK
            DO 40 J=L,LL
              DO 40 IA=L,KM1
                C(I,J) = C(I,J) - A(IA,I)*C(IA,J)
   40     CONTINUE
   45     IF(DL .EQ. 2) GO TO 60
          IF(DK .EQ. 2 ) GO TO 50
          T(1,1) = A(K,K) + A(L,L)
          IF(T(1,1) .EQ. 0.) STOP
          C(K,L) = C(K,L)/T(1,1)
          GO TO 90
  50      T(1,1) = A(K,K) + A(L,L)
          T(1,2) = A(KK,K)
          T(2,1) = A(K,KK)
          T(2,2) = A(KK,KK) + A(L,L)
          P(1) = C(K,L)
          P(2) = C(KK,L)
          NSYS = 2
          CALL SYSSLV
          C(K,L) = P(1)
        C(KK,L) = P(2)
          GO TO 90
  60      IF(DK .EQ. 2) GO TO 70
          T(1,1) = A(K,K) + A(L,L)
          T(1,2) = A(LL,L)
          T(2,1) = A(L,LL)
          T(2,2) = A(K,K) + A(LL,LL)
          P(1) = C(K,L)
          P(2) = C(K,LL)
          NSYS = 2
          CALL SYSSLV
          C(K,L) = P(1)
          C(K,LL) = P(2)
          GO TO 90
  70      IF(K .NE. L) GO TO 80
          T(1,1) = A(L,L)
          T(1,2) = A(LL,L)
          T(1,3) = 0.
          T(2,1) = A(L,LL)
          T(2,2) = A(L,L) + A(LL,LL)
          T(2,3) = T(1,2)
          T(3,1) = 0.
          T(3,2) = T(2,1)
          T(3,3) = A(LL,LL)
          P(1) = C(L,L)/2.
          P(2) = C(LL,L)
          P(3) = C(LL,LL)/2.
          NSYS = 3
          CALL SYSSLV
          C(L,L) = P(1)
          C(LL,L) = P(2)
          C(L,LL) = P(2)
          C(LL,LL) = P(3)
          GO TO 90
  80      T(1,1) = A(K,K) + A(L,L)
          T(1,2) = A(KK,K)
          T(1,3) = A(LL,L)
          T(1,4) = 0.
          T(2,1) = A(K,KK)
          T(2,2) = A(KK,KK) + A(L,L)
          T(2,3) = 0.
          T(2,4) = T(1,3)
          T(3,1) = A(L,LL)
          T(3,2) = 0.
          T(3,3) = A(K,K) + A(LL,LL)
          T(3,4) = T(1,2)
          T(4,1) = 0.
          T(4,2) = T(3,1)
          T(4,3) = T(2,1)
          T(4,4) = A(KK,KK) + A(LL,LL)
          P(1) = C(K,L)
          P(2) = C(KK,L)
          P(3) = C(K,LL)
          P(4) = C(KK,LL)
          NSYS = 4
          CALL SYSSLV
          C(K,L) = P(1)
          C(KK,L) = P(2)
          C(K,LL) = P(3)
          C(KK,LL) = P(4)
  90    K = K + DK
        IF(K .LE. N) GO TO 30
        LDL = L + DL
        IF(LDL .GT. N) RETURN
        DO 120 J=LDL,N
          DO 100 I=L,LL
            C(I,J) = C(J,I)
  100     CONTINUE
          DO 120 I=J,N
            DO 110 K=L,LL
              C(I,J) = C(I,J) - C(I,K)*A(K,J) - A(K,I)*C(K,J)
  110     CONTINUE
          C(J,I) = C(I,J)
  120 CONTINUE
      L = LDL
      GO TO 10
      END
      SUBROUTINE SYSSLV
C 
C   PURPOSE:
C      Solve the linear system Ax = b, where A is an n x n (n <= 5)
C      matrix and b is n-dimensional vector.  Solution is by Crout re-
C      duction.  The matrix A, the vector b, and order n are contained
C      in the arrays A, B, and the variable N of the COMMON block
C      SLVBLK.
C 
C   REFERENCES:
C      Bartels, R.H.; and Stewart, G.W.: Algorithm 432 - Solution of
C        the Matrix Equation AX + XB = C.  Commun. ACM, vol. 15, no. 9,
C        Sept. 1972, pp. 820-826.
C 
C   Subroutines employed by SYSSLV: None
C   Subroutines employing SYSSLV: SHRSLV, SYMSLV
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SLVBLK/A(5,5),B(5),N
      REAL*8 MAX
    1 NM1 = N - 1
      N1 = N+1
C 
C COMPUTE THE LU FACTORIZATION OF A.
C 
      DO 80 K=1,N
        KM1 = K-1
        IF(K.EQ.1) GO TO 20
        DO 10 I=K,N
          DO 10 J=1,KM1
            A(I,K) = A(I,K) - A(I,J)*A(J,K)
   10   CONTINUE
   20   IF(K.EQ.N) GO TO 100
        KP1 = K+1
      MAX = DABS(A(K,K))
        INTR = K
        DO 30 I=KP1,N
          AA = DABS(A(I,K))
          IF(AA .LE. MAX) GO TO 30
          MAX = AA
          INTR = I
   30   CONTINUE
        IF(MAX .EQ. 0.) STOP
        A(N1,K) = INTR
        IF(INTR .EQ. K) GO TO 50
        DO 40 J=1,N
          TEMP = A(K,J)
          A(K,J) = A(INTR,J)
          A(INTR,J) = TEMP
   40   CONTINUE
   50   DO 80 J=KP1,N
          IF(K.EQ.1) GO TO 70
          DO 60 I=1,KM1
            A(K,J) = A(K,J) - A(K,I)*A(I,J)
   60     CONTINUE
   70     A(K,J) = A(K,J)/A(K,K)
   80 CONTINUE
C 
C INTERCHANGE THE COMPONENTS OF B.
C 
  100 DO 110 J=1,NM1
        INTR = A(N1,J)
        IF(INTR .EQ. J) GO TO 110
        TEMP = B(J)
        B(J) = B(INTR)
        B(INTR) = TEMP
  110 CONTINUE
C 
C SOLVE LX = B.
C 
  200 B(1) = B(1)/A(1,1)
      DO 220 I=2,N
        IM1 = I-1
        DO 210 J=1,IM1
          B(I) = B(I) - A(I,J)*B(J)
  210   CONTINUE
        B(I) = B(I)/A(I,I)
  220 CONTINUE
C 
C SOLVE UX = B.
C 
  300 DO 310 II=1,NM1
          I = NM1-II+1
        I1 = I+1
        DO 310 J=I1,N
          B(I) = B(I) - A(I,J)*B(J)
  310 CONTINUE
      RETURN
      END
      SUBROUTINE TESTST(A,NA,ALPHA,DISC,STABLE,IOP,DUMMY)
C 
C   PURPOSE:
C      Compute and test the eigenvalues of a constant real matrix A for
C      stability relative to a parameter 'alpha' in either the continu-
C      ous or digital sense.
C 
C   Subroutines employed by TESTST: EIGEN, EQUATE, JUXTC, LNCNT, PRNT
C   Subroutines employing TESTST: ASYREG
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),DUMMY(1)
      DIMENSION NA(2),NDUM1(2),NDUM2(2)
      LOGICAL  DISC,STABLE
C 
      STABLE = .FALSE.
C 
      CALL EQUATE(A,NA,DUMMY,NA)
      N1= NA(1)**2 + 1
      N2= N1+NA(1)
      N3= N2+NA(1)
      ISV = 0
      CALL EIGEN(NA(1),NA(1),DUMMY,DUMMY(N1),DUMMY(N2),ISV,ISV,V,DUMMY(N
     13),IERR)
      NEVL = NA(1)
      IF( IERR .EQ. 0 ) GO TO 200
      CALL LNCNT(4)
      WRITE(6,100) IERR
  100 FORMAT(BZ,//' IN TESTST, THE ',I5,' EIGENVALUE OF A HAS NOT BEEN
     1FOUND AFTER 30 ITERATIONS'/)
      RETURN
C 
  200 CONTINUE
      NDUM1(1) = NEVL
      NDUM1(2) = 1
      CALL JUXTC(DUMMY(N1),NDUM1,DUMMY(N2),NDUM1,DUMMY,NDUM2)
C 
      IF( DISC ) GO TO 400
      DO 300 I=1,NEVL
      IF( DUMMY(I) .GE. ALPHA ) GO TO 600
  300 CONTINUE
      GO TO 550
  400 CONTINUE
      N = NDUM2(1)*NDUM2(2)+1
      DO 500 I =1,NEVL
      K = I + NEVL
      L=N +I -1
      DUMMY(L) = DSQRT((DUMMY(I)**2)+(DUMMY(K)**2))
  500 CONTINUE
C 
      IF( DUMMY(L) .GE. ALPHA ) GO TO 600
C 
  550 CONTINUE
      STABLE =.TRUE.
  600 CONTINUE
      IF( IOP .EQ. 0 ) RETURN
      CALL LNCNT(4)
      WRITE(6,700)
  700 FORMAT(//' PROGRAM TO TEST THE RELATIVE STABILITY OF THE MATRIX A
     1'/)
      CALL PRNT(A,NA,' A  ',1)
      CALL LNCNT(4)
  750 FORMAT(//' EIGENVALUES OF A '/)
      CALL PRNT(DUMMY,NDUM2,'EVLA',1)
      IF(  .NOT. DISC ) GO TO 850
      CALL LNCNT(4)
      WRITE(6,800)
  800 FORMAT(//' MODULI OF EIGENVALUES OF A'/)
      CALL PRNT(DUMMY(N),NDUM1,'MODA',1)
C 
  850 CONTINUE
      CALL LNCNT(4)
      IF(STABLE) WRITE(6,900) ALPHA
      IF( .NOT. STABLE) WRITE(6,950) ALPHA
  900 FORMAT(//' MATRIX A  IS STABLE RELATIVE TO ',D16.8,/)
  950 FORMAT(//' MATRIX A  IS UNSTABLE RELATIVE TO ' ,D16.8,/)
C 
      RETURN
      END
      SUBROUTINE TRANP(A,NA,B,NB)
C 
C   PURPOSE:
C      Compute the transpose A' of a given matrix A.
C 
C   Subroutines employed by TRANP: LNCNT
C   Subroutines employing TRANP: ASYFIL, ASYREG, BARSTW, BILIN, CNTREG,
C      CSTAB, CTROL, DISREG, DSTAB, EXMDFL, FACTOR, IMMDFL, PREFIL,
C      RICNWT, SAMPL, SUM, TRNSIT, VARANC
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),NA(2),NB(2)
      NB(1)=NA(2)
      NB(2)=NA(1)
      NR=NA(1)
      NC=NA(2)
      L=NR*NC
      IF( NR .LT. 1 .OR. L .LT. 1  )  GO TO 999
      IR=0
      DO 300 I=1,NR
      IJ=I-NR
      DO 300 J=1,NC
      IJ=IJ+NR
      IR=IR+1
  300 B(IR)=A(IJ)
      RETURN
  999 CALL LNCNT(1)
      WRITE  (6,50)  NA
   50 FORMAT  (' DIMENSION ERROR IN TRANP   NA=',2I6)
      RETURN
      END
      SUBROUTINE TRCE  (A,NA,TR)
C 
C   PURPOSE:
C      Compute the trace of a square matrix.
C 
C   Subroutines employed by TRCE: LNCNT
C   Subroutines employing TRCE: EXPSER
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1)
      DIMENSION NA(2)
      IF (NA(1).NE.NA(2)) GO TO 600
      N=NA(1)
      TR=0.0
      IF( N .LT. 1 ) GO TO 600
      DO 10 I=1,N
      M=I+N*(I-1)
   10 TR=TR+A(M)
      RETURN
  600 CALL LNCNT(1)
      WRITE (6,1600) NA
 1600 FORMAT (' TRACE REQUIRES SQUARE MATRIX    NA=',2I6)
      RETURN
      END

      SUBROUTINE TRNSIT(A,NA,B,NB,H,NH,G,NG,F,NF,V,NV,T,X,NX,DISC,
     1  STABLE,IOP,DUMMY)
C 
C   PURPOSE:
C      Compute and print the transient response of the time-invariant
C      continuous or discrete system.
C 
C   Subroutines employed by TRNSIT: ADD, DGECO, DGESL, EQUATE, EXPINT,
C      EXPSER, LNCNT, MULT, PRNT, SCALE, SUBT, TRANP, UNITY
C   Subroutines employing TRNSIT: none
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),B(*),H(*),G(*),F(*),V(*),X(*),DUMMY(:)
      DIMENSION NA(2),NB(2),NH(2),NG(2),NF(2),NV(2),NX(2),T(2),IOP(4)
      DIMENSION NDUM1(2),NDUM2(2)
      LOGICAL  DISC,STABLE
C 
      N = NA(1)*NA(2)
      N1 = N + 1
      N2 = N + N1
      N3 = N + N2
      N4 = N + N3
      N5 = N + N4
      N6 = N + N5
C 
      CALL LNCNT(4)
      IF(DISC) WRITE(6,100)
      IF( .NOT. DISC ) WRITE(6,120)
  100 FORMAT(BZ,//
     & ' COMPUTATION OF TRANSIENT RESPONSE FOR THE DIGITAL SYSTEM '/)
  120 FORMAT(//
     & ' COMPUTATION OF TRANSIENT RESPONSE FOR THE CONTINUOUS SYSTEM'/)
      CALL PRNT(A,NA,' A  ',1)
      CALL PRNT(B,NB,' B  ',1)
      IF( (IOP(1) .NE. 1) .AND. (IOP(1) .NE. 0) )  GO TO 180
      CALL LNCNT(3)
      IF( IOP(1) .EQ. 0 ) WRITE(6,140)
      IF( IOP(1) .EQ. 1 ) WRITE(6,160)
  140 FORMAT(//' H IS A NULL MATRIX ')
  160 FORMAT(//' H IS AN IDENTITY MATRIX ')
      GO TO 200
  180 CONTINUE
      CALL PRNT(H,NH,' H  ',1)
  200 CONTINUE
      IF( (IOP(2) .NE. 1) .AND. (IOP(2) .NE. 0) )  GO TO 260
      CALL LNCNT(3)
      IF( IOP(2) .EQ. 0 ) WRITE(6,220)
      IF( IOP(2) .EQ. 1 ) WRITE(6,240)
  220 FORMAT(//' G IS A NULL MATRIX')
  240 FORMAT(//' G IS AN IDENTITY MATRIX')
      GO TO 280
  260 CONTINUE
      CALL PRNT(G,NG,' G  ',1)
  280 CONTINUE
      CALL PRNT(F,NF,' F  ',1)
      IF( (IOP(3) .NE. 0) .AND. (IOP(3) .NE. 1) )  GO TO 295
      CALL LNCNT(3)
      IF(IOP(3).EQ.0) WRITE(6,285)
      IF(IOP(3).EQ.1) WRITE(6,290)
  285 FORMAT(//' V IS A NULL MATRIX')
  290 FORMAT(//' V IS AN IDENTITY MATRIX')
      GO TO 300
  295 CONTINUE
      CALL PRNT(V,NV,' V  ',1)
C 
  300 CONTINUE
      CALL EQUATE(A,NA,DUMMY(N6),NA)
      CALL MULT(B,NB,F,NF,DUMMY,NA)
      CALL SUBT(A,NA,DUMMY,NA,A,NA)
C 
      IF(DISC) GO TO 350
      NMAX = T(1)/T(2)
      IOPT = 1
      TT = T(2)
      IF( IOP(3) .NE. 0 )  GO TO 315
      CALL EXPSER(A,NA,DUMMY,NA,TT,IOPT,DUMMY(N1))
      GO TO 400
  315 CONTINUE
      CALL EXPINT(A,NA,DUMMY,NA,DUMMY(N1),NA,TT,IOPT,DUMMY(N2))
      CALL MULT(DUMMY(N1),NA,B,NB,DUMMY(N2),NB)
      IF( IOP(3) .NE. 1 ) GO TO 325
      CALL EQUATE(DUMMY(N2),NB,DUMMY(N1),NX)
      GO TO 400
  325 CONTINUE
      CALL MULT(DUMMY(N2),NB,V,NV,DUMMY(N1),NX)
      GO TO 400
  350 CONTINUE
      NMAX = IOP(4)
      CALL EQUATE(A,NA,DUMMY,NA)
      IF( IOP(3) .EQ. 0 ) GO TO 400
      CALL MULT(B,NB,V,NV,DUMMY(N1),NX)
C 
  400 CONTINUE
      CALL LNCNT(4)
      WRITE(6,420)
  420 FORMAT(//' STRUCTURE OF PRINTING TO FOLLOW'/)
      CALL LNCNT(6)
      WRITE(6,440)
  440 FORMAT('   TIME OR STAGE '/'  STATE - X TRANSPOSE - FROM DX = AX +
     1 BU'/'  OUTPUT - Y TRANSPOSE - FROM Y = HX + GU   IF DIFFERENT FRO
     2M X'/'  CONTROL - U TRANSPOSE - FROM U = -FX + V'//)
C 
      K = 0
      L = 0
      CALL SCALE(F,NF,F,NF,-1.0D0)
C 
  450 CONTINUE
      IF( K .GT. NMAX ) GO TO 800
      CALL MULT(F,NF,X,NX,DUMMY(N2),NV)
      IF( IOP(3) .NE. 0 ) CALL ADD(DUMMY(N2),NV,V,NV,DUMMY(N2),NV)
      CALL MULT(DUMMY,NA,X,NX,DUMMY(N3),NX)
      IF( IOP(3) .EQ. 0 ) GO TO 475
      CALL ADD(DUMMY(N1),NX,DUMMY(N3),NX,DUMMY(N3),NX)
  475 CONTINUE
      IF( IOP(2) .EQ. 0 ) GO TO 525
      IF( IOP(2) .EQ. 1 ) GO TO 500
      CALL MULT(G,NG,DUMMY(N2),NV,DUMMY(N4),NDUM1)
      GO TO 525
  500 CONTINUE
      CALL EQUATE(DUMMY(N2),NV,DUMMY(N4),NDUM1)
  525 CONTINUE
      IF( IOP(1) .EQ. 0 ) GO TO 575
      IF( IOP(1) .EQ. 1 ) GO TO 550
      CALL MULT(H,NH,X,NX,DUMMY(N5),NDUM1)
      GO TO 575
  550 CONTINUE
      CALL EQUATE(X,NX,DUMMY(N5),NDUM1)
  575 CONTINUE
      IF( IOP(2) .EQ. 0 ) GO TO 600
      IF( IOP(1) .EQ. 0 ) GO TO 700
      CALL ADD(DUMMY(N4),NDUM1,DUMMY(N5),NDUM1,DUMMY(N4),NDUM1)
      GO TO 700
  600 CONTINUE
      IF( IOP(1) .NE. 0 ) CALL EQUATE(DUMMY(N5),NDUM1,DUMMY(N4),NDUM1)
C 
  700 CONTINUE
      CALL LNCNT(5)
      IF( .NOT. DISC ) GO TO 720
      WRITE(6,710) K
  710 FORMAT(////,I5)
      GO TO 740
  720 CONTINUE
      TIME=K*T(2)
      WRITE(6,730) TIME
  730 FORMAT(////,D16.7)
  740 CONTINUE
      CALL TRANP(X,NX,DUMMY(N5),NDUM2)
      CALL PRNT(DUMMY(N5),NDUM2,'    ',3)
      IF( (IOP(2) .EQ. 0) .AND. ( (IOP(1) .EQ. 0) .OR. (IOP(1) .EQ. 1) )
     1) GO TO 750
      CALL TRANP(DUMMY(N4),NDUM1,DUMMY(N5),NDUM2)
      CALL PRNT(DUMMY(N5),NDUM2,'    ',3)
  750 CONTINUE
      CALL TRANP(DUMMY(N2),NV,DUMMY(N5),NDUM2)
      CALL PRNT(DUMMY(N5),NDUM2,'    ',3)
C 
      CALL EQUATE(DUMMY(N3),NX,X,NX)
      K = K + 1
      GO TO 450
C 
C 
  800 CONTINUE
      CALL SCALE(F,NF,F,NF,-1.0D0)
      IF( .NOT. STABLE  .OR.  IOP(3) .EQ. 0  )   GO TO 900
      IF( IOP(3) .EQ. 1 )  GO TO 820
      CALL MULT(B,NB,V,NV,DUMMY,NX)
      GO TO 840
  820 CONTINUE
      CALL EQUATE(B,NB,DUMMY,NX)
  840 CONTINUE
      IF( .NOT. DISC )  GO TO 860
      CALL UNITY(DUMMY(N1),NA)
      CALL SUBT(DUMMY(N1),NA,A,NA,A,NA)
  860 CONTINUE
C 
C   * * * CALL TO MATHLIB FUNCTIONS * * *
      CALL DGECO(A,NA(1),NA(1),DUMMY(N1:),RCOND,DUMMY(N2:))   ! changed RLC
      IF ((1.0 + RCOND) .NE. 1.0) GO TO 880
      CALL LNCNT(3)
      IF( .NOT. DISC )  WRITE(6,865) RCOND
      IF( DISC )  WRITE(6,870) RCOND
  865 FORMAT(//' IN TRNSIT, THE MATRIX A-BF SUBMITTED TO DGECO IS SINGU
     1LAR, RCOND = ',D16.8)
  870 FORMAT(//' IN TRNSIT, THE MATRIX  I - (A-BF) SUBMITTED TO DGECO I
     1S SINGULAR, RCOND = ',D16.8)
      GO TO 900
  880 CONTINUE
      NT = 1
      DO 885 M1 = 1,NX(2)
         CALL DGESL(A,NA(1),NA(1),DUMMY(N1:),DUMMY(NT:),0)   ! changed RLC
         NT = NT + NA(1)
  885 CONTINUE
      IF( .NOT. DISC )  CALL SCALE(DUMMY,NX,DUMMY,NX,-1.0D0)
      CALL LNCNT(5)
      WRITE(6,890)
  890 FORMAT(////' STEADY-STATE VALUE OF  X TRANSPOSE')
      CALL TRANP(DUMMY,NX,DUMMY(N5),NDUM2)
      CALL PRNT(DUMMY(N5),NDUM2,'    ',3)
C 
  900 CONTINUE
      CALL EQUATE(DUMMY(N6),NA,A,NA)
C 
      RETURN
      END
      SUBROUTINE UNITY(A,NA)
C 
C   PURPOSE:
C      Generate an identity matrix.
C 
C   Subroutines employed by UNITY: LNCNT
C   Subroutines employing UNITY: BILIN, EXPINT, EXPSER, TRNSIT
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),NA(2)
      IF(NA(1).NE.NA(2)) GO TO 999
      L=NA(1)*NA(2)
      DO 100 IT=1,L
  100 A(IT)=0.0
      J = - NA(1)
      NAX = NA(1)
      DO 300 I=1,NAX
      J=NAX +J+1
  300 A(J)=1.
      GO TO 1000
  999 CALL LNCNT (1)
      WRITE(6, 50)(NA(I),I=1,2)
   50 FORMAT  (' DIMENSION ERROR IN UNITY   NA=',2I6)
 1000 RETURN
      END
      SUBROUTINE VARANC(A,NA,G,NG,Q,NQ,W,NW,IDENT,DISC,IOP,DUMMY)
C 
C   PURPOSE:
C      Compute the steady-state variance matrix of the state of the con-
C      tinuous or discrete linear time-invariant system.
C 
C   REFERENCES:
C      Kwakernaak, Huibert; and Sivan, Raphael: Linear Optimal Control
C        Systems.  John Wiley & Sons, Inc., c.1972.
C 
C   Subroutines employed by VARANC: BARSTW, BILIN, EQUATE, LNCNT, MULT,
C      PRNT, SCALE, SUM, TRANP
C   Subroutines employing VARANC: None
C 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*),G(*),Q(*),W(*),DUMMY(:)
      DIMENSION NA(2),NG(2),NQ(2),NW(2),NDUM1(2),IOP(3),IOPT(2)
      LOGICAL IDENT,DISC,SYM
      COMMON/TOL/EPSAM,EPSBM,IACM
      INTEGER:: ndummy
C 
      ndummy=SIZE(dummy)
      IF( IOP(1) .EQ. 0 ) GO TO 100
      CALL LNCNT(5)
      IF( DISC ) WRITE(6,25)
      IF( .NOT. DISC ) WRITE(6,35)
   25 FORMAT(//' PROGRAM TO SOLVE FOR THE STEADY-STATE VARIANCE MATRIX'/
     /' FOR A LINEAR DISCRETE SYSTEM'/)
   35 FORMAT(//' PROGRAM TO SOLVE FOR THE STEADY-STATE VARIANCE MATRIX'/
     1' FOR A LINEAR CONTINUOUS SYSTEM'/)
      CALL PRNT(A,NA,' A  ',1)
      IF( .NOT. IDENT ) GO TO 55
      CALL LNCNT(3)
      WRITE(6,45)
   45 FORMAT(/' G IS AN IDENTITY MATRIX '/)
      GO TO 65
   55 CONTINUE
      CALL PRNT(G,NG,' G  ',1)
   65 CONTINUE
      IF ( .NOT. IDENT ) GO TO 85
      CALL LNCNT(3)
      WRITE(6,75)
   75 FORMAT(/' INTENSITY MATRIX FOR COVARIANCE OF PROCESS NOISE '/)
C 
   85 CONTINUE
      CALL PRNT(Q,NQ,' Q  ',1)
C 
  100 CONTINUE
      IF( IDENT ) GO TO 200
      CALL MULT(G,NG,Q,NQ,DUMMY,NG)
      N1 = NG(1)*NG(2) + 1
      CALL TRANP(G,NG,DUMMY(N1),NDUM1)
      CALL MULT(DUMMY,NG,DUMMY(N1),NDUM1,Q,NQ)
C 
      IF( IOP(1) .EQ. 0 ) GO TO 200
      CALL LNCNT(3)
      WRITE(6,75)
      CALL PRNT(Q,NQ,'GQGT',1)
C 
  200 CONTINUE
      CALL EQUATE(Q,NQ,W,NW)
      IF(.NOT. DISC) CALL SCALE(W,NW,W,NW,-1.0D0)
      IOPT(1) = IOP(2)
      IOPT(2) = 1
      SYM = .TRUE.
      IF( DISC ) GO TO 300
      IF( IOP(3) .EQ. 0 ) GO TO 250
      CALL BILIN(A,NA,A,NA,W,NW,IOPT,BETA,SYM,DUMMY)
      GO TO 400
C 
  250 CONTINUE
      EPSA=EPSAM
      CALL BARSTW(A,NA,A,NA,W,NW,IOPT(1),SYM,EPSA,EPSA,DUMMY)   ! RLC
      GO TO 400
C 
  300 CONTINUE
      CALL EQUATE(A,NA,DUMMY,NA)
      N = NA(1)**2
      N1 = N + 1
      CALL TRANP(A,NA,DUMMY(N1),NA)
      N2 = N1 + N
      CALL SUM(DUMMY,NA,W,NW,DUMMY(N1:ndummy),NA,                          &
     & IOPT(1),SYM,DUMMY(N2:ndummy))   ! RLC
C 
  400 CONTINUE
      IF( IOP(1) .EQ. 0 ) RETURN
      CALL LNCNT(3)
      WRITE(6,450)
  450 FORMAT(/ ' VARIANCE MATRIX '/)
      CALL PRNT(W,NW,' W  ',1)
C 
      RETURN
      END
