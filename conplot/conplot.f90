!+
MODULE Smooth
! ---------------------------------------------------------------------------
! PURPOSE -

IMPLICIT NONE

  INTEGER,PARAMETER,PRIVATE:: SP=KIND(1.0), DP=KIND(1.0D0)
  PRIVATE:: Hypot
  PUBLIC:: SmoothSurface
  PRIVATE:: SolveSVD,SVD,SVDbackSubstitution

CONTAINS

!+
      SUBROUTINE SmoothSurface(x, y, z, znew, iexp, jexp,  &
        ncoef, cnorm, ipowr, jpowr)
! ---------------------------------------------------------------------------
! PURPOSE - Performs the optional smoothing of data before triangulation
!  of the plane is initiated.  The surface defined by z = f(x,y) is
!  smoothed via a polynomial curve fit defined by a least squares criteria.

IMPLICIT NONE

  REAL,INTENT(IN),DIMENSION(:):: x,y,z   ! arrays of values defining the
                                         ! known surface (points in space
                                         ! for the function z=f(x,y))
  INTEGER,INTENT(IN):: iexp,jexp         ! exponents for the smoothing
                                         ! polynomial as selected by the user.
  REAL,INTENT(OUT),DIMENSION(:):: znew   ! array of smoothed values for the
                                         ! function (znew will contain the
                                         ! original z data on return if the
                                         ! smoothing operation fails, in
                                         ! which case ncoef will be set to -1).
  INTEGER,INTENT(OUT):: ncoef            ! the number of terms in the
                                         ! polynomial resulting from the
                                         ! values of i and j.
  REAL,INTENT(OUT),DIMENSION(:):: cnorm  ! array of ncoef computed coefficients
  INTEGER,INTENT(OUT),DIMENSION(:):: ipowr ! array of i exponents for each term
  INTEGER,INTENT(OUT),DIMENSION(:):: jpowr ! array of j exponents for each term
!                  (each element of c, ipowr and jpowr is associated
!                  with the ncoef terms of the polynomial, in order)
!

  REAL,ALLOCATABLE,DIMENSION(:,:):: am
  REAL,ALLOCATABLE,DIMENSION(:):: b
  REAL,ALLOCATABLE,DIMENSION(:):: c,ave
!  REAL,DIMENSION(23):: h
  INTEGER:: i,j
  INTEGER:: i1,j1
!  INTEGER,PARAMETER:: IA=500
  INTEGER:: ier
  INTEGER:: iex,jex
  INTEGER:: ii
!  INTEGER,DIMENSION(23):: ip
  INTEGER:: k,l
!  INTEGER:: m
  INTEGER:: jj
!  INTEGER:: kbasis
  INTEGER:: kcol,krow
  INTEGER:: ki1
  INTEGER:: l1
  INTEGER:: maxCoef
  INTEGER:: n
  REAL:: realn
!  REAL:: tol
  REAL,ALLOCATABLE,DIMENSION(:):: xx
  REAL:: x2,y2
  REAL:: xp,yp

!----------------------------------------------------------------------------
! (A) initialize local variables and range check
!
      n=MIN(SIZE(x),SIZE(y),SIZE(z),SIZE(znew))
      realn=REAL(n)
      maxCoef=MIN(SIZE(ipowr),SIZE(jpowr))
      i=MAX(iexp,1)
      j=MAX(jexp,1)
      i1=i+1
      j1=j+1
      ncoef=0
 
 
! (B) Determine the x and y exponents to be used and save them in
! arrays ipowr and jpowr

      ncoef=0
      k=MAX(i1,j1)
!      IF (k.eq.0) GO TO 110
      DO ii=1,i1
        ki1=k-ii+1
        l=MIN(ki1,j1)
        DO jj=1,l
          ncoef=ncoef+1
          IF (ncoef > maxCoef) THEN
            znew(1:n)=z(1:n)           ! error, no room to store coeff.
            ncoef=-1
            RETURN
          END IF
          ipowr(ncoef)=ii-1
          jpowr(ncoef)=jj-1
        END DO  
      END DO   
      WRITE(2,'(3I5)') (k,ipowr(k),jpowr(k),k=1,ncoef) 
 
! (C) Using the exponent lists from above and the known xy data points,
! construct the matrix am

      ALLOCATE(am(n,ncoef), ave(ncoef), b(n), c(ncoef), xx(ncoef) )
      DO kcol=1,ncoef
        iex=ipowr(kcol)
        jex=jpowr(kcol)
        DO krow=1,n
          x2=x(krow)
          IF (x2.eq.0.0) x2=1.0
          xp=x2**iex
          y2=y(krow)
          IF (y2.eq.0.0) y2=1.0
          yp=y2**jex
          am(krow,kcol)=xp*yp
        END DO  
      END DO   
      krow=ncoef
 
 
! (D) Normalize each value in each column of am by the column average
 
      ave(1)=1.0
      DO l1=2,ncoef
        ave(l1)=SUM(ABS(am(1:n,l1)) )  ! column l1
!        ave(l1)=0.0
!        DO l2=1,n
!          ave(l1)=ave(l1)+abs(am(l2,l1))
!        END DO
        ave(l1)=ave(l1)/realn
        IF (ave(l1).eq.0.) ave(l1)=1.0
        am(1:n,l1)=am(1:n,l1)/ave(l1)
!        DO l2=1,n
!          am(l2,l1)=am(l2,l1)/ave(l1)
!        END DO
      END DO   
 
 
! (E,F,G) Use imsl routine llsqf to solve (via least-squares)
!        the system  am*c = z  for matrix c

!      m=n
      ier=0
!      kbasis=ncoef
!      tol=0.0
      b(1:n)=z(1:n)

      CALL SolveSVD(am,b,xx,ier)
!      CALL LLSQF (AM,IA,M,NCOEF,B,TOL,KBASIS,XX,H,IP,IER)
!      xx=0.0    !!!!!!!!!!!!!!!!!!!!!!!!!!!kludge
      DEALLOCATE(am,b)

      IF (ier/=0) THEN
        znew(1:n)=z(1:n)  ! error return
        ncoef=-1
        RETURN
      END IF

 
 
! (H) Divide out the scale factor from the solution matrix and
! establish the coefficients
      c(1:ncoef)=xx(1:ncoef)
      cnorm(1:ncoef)=c(1:ncoef)/ave(1:ncoef)
  write(2,'(3I5,2F15.6)') (k,ipowr(k),jpowr(k),c(k),cnorm(k),k=1,ncoef)
      DEALLOCATE(ave,c,xx)

!      DO 90 l3=1,ncoef
!        c(l3)=xx(l3)
!        cnorm(l3)=c(l3)/ave(l3)
!   90 CONTINUE
 
 
! (I) Establish the new z values by evaluating the polynomial for each
! known x-y pair

      DO k=1,n
        znew(k)=-1.0*Polyx2(0.0, x(k),y(k), cnorm,ipowr,jpowr,ncoef)
      END DO   
  WRITE(2,'(I4,4F12.6)') (k,x(k),y(k),z(k),znew(k),k=1,n)

  RETURN
END Subroutine SmoothSurface   ! --------------------------------------------

!+
FUNCTION Polyx2 (z, x, y, c, ipowr, jpowr, ncoef) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Polyx2 is the polynomial evaluation function used when the
!  smoothing option has been invoked.  x and y are the known values of the
!  independent variables.  c is the list of coefficients for each term.
!  ipowr and jpowr are the exponents for each term and ncoef is the number
!  of terms in the polynomial. z is an offset term when evaluating for a
!  constant x value.
! 
IMPLICIT NONE
  REAL,INTENT(IN):: z,x,y
  REAL,INTENT(IN),DIMENSION(:):: c
  INTEGER,INTENT(IN),DIMENSION(:):: ipowr,jpowr
  INTEGER,INTENT(IN):: ncoef

  REAL:: f
  INTEGER:: k
  REAL:: term1,term2

!----------------------------------------------------------------------------
      IF ((x.eq.0.0) .OR. (y.eq.0.0)) GO TO 30
      f=c(1)
      DO k=2,ncoef
        f=f+(x**ipowr(k))*(y**jpowr(k))*c(k)
      END DO   
      f=z-f
      RETURN


   30 f=c(1)
      DO k=2,ncoef
        term1=1.0
        IF (ipowr(k).NE.0) term1=x**ipowr(k)
        term2=1.0
        IF (jpowr(k).NE.0) term2=y**jpowr(k)
        f=f+term1*term2*c(k)
      END DO   
      f=z-f
      RETURN
END Function PolyX2   ! -----------------------------------------------------

!+
SUBROUTINE SolveSVD(a,b,x,ier)
! ---------------------------------------------------------------------------
! PURPOSE - Solve the system Ax=b in least-squares sense using SVD
IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:,:):: a
  REAL,INTENT(IN),DIMENSION(:):: b
  REAL,INTENT(OUT),DIMENSION(:):: x
  INTEGER,INTENT(OUT):: ier

  INTEGER:: m,n
  REAL,ALLOCATABLE,DIMENSION(:,:):: u,v
  REAL,ALLOCATABLE,DIMENSION(:):: w

!----------------------------------------------------------------------------
  m=SIZE(a,1)
  n=SIZE(a,2)
  ier=0
  IF (SIZE(b)<n) ier=-10
  IF (SIZE(x)<n) ier=-11
  IF (ier < 0) RETURN
  ALLOCATE(u(m,n), w(n), v(n,n))

  CALL SVD(a,u,w,v,ier)
  IF (ier /= 0) RETURN

  CALL SVDbackSubstitution(u,w,v,b,x)
  DEALLOCATE(v,w,u)

  RETURN
END Subroutine SolveSVD   ! -------------------------------------------------

!+
SUBROUTINE SVD(a,u,w,v,ierr)
!SUBROUTINE Svd(nm, m, n, a, w, matu, u, matv, v, ierr, rv1)
!----------------------------------------------------------------------------
! PURPOSE - Perform the singular value decomposition of a rectangular matrix.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SVD-S)
!***AUTHOR  (UNKNOWN)

IMPLICIT NONE

!!! nm, m, n, a, w, matu, u, matv, v, ierr, rv1)

  REAL,INTENT(IN),DIMENSION(:,:):: a   ! matrix to be decomposed
  REAL,INTENT(OUT),DIMENSION(:,:):: u  ! the matrix U 
  REAL,INTENT(OUT),DIMENSION(:):: w
  REAL,INTENT(OUT),DIMENSION(:,:):: v
  INTEGER,INTENT(OUT):: ierr

!***DESCRIPTION
! 
!     This subroutine is a translation of the ALGOL procedure SVD,
!     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
!     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
! 
!     This subroutine determines the singular value decomposition
!          T
!     A=USV  of a REAL M by N rectangular matrix.  Householder
!     bidiagonalization and a variant of the QR algorithm are used.
! 
!     On Input
! 
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A, U and V, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!          Note that NM must be at least as large as the maximum
!          of M and N.
! 
!        M is the number of rows of A and U.
! 
!        N is the number of columns of A and U and the order of V.
! 
!        A contains the rectangular input matrix to be decomposed.  A is
!          a two-dimensional REAL array, dimensioned A(NM,N).
! 
!        MATU should be set to .TRUE. if the U matrix in the
!          decomposition is desired, and to .FALSE. otherwise.
!          MATU is a LOGICAL variable.
! 
!        MATV should be set to .TRUE. if the V matrix in the
!          decomposition is desired, and to .FALSE. otherwise.
!          MATV is a LOGICAL variable.
! 
!     On Output
! 
!        A is unaltered (unless overwritten by U or V).
! 
!        W contains the N (non-negative) singular values of A (the
!          diagonal elements of S).  They are unordered.  If an
!          error exit is made, the singular values should be correct
!          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
!          REAL array, dimensioned W(N).
! 
!        U contains the matrix U (orthogonal column vectors) of the
!          decomposition if MATU has been set to .TRUE.  Otherwise,
!          U is used as a temporary array.  U may coincide with A.
!          If an error exit is made, the columns of U corresponding
!          to indices of correct singular values should be correct.
!          U is a two-dimensional REAL array, dimensioned U(NM,N).
! 
!        V contains the matrix V (orthogonal) of the decomposition if
!          MATV has been set to .TRUE.  Otherwise, V is not referenced.
!          V may also coincide with A if U does not.  If an error
!          exit is made, the columns of V corresponding to indices of
!          correct singular values should be correct.  V is a two-
!          dimensional REAL array, dimensioned V(NM,N).
! 
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          K          if the K-th singular value has not been
!                     determined after 30 iterations.
! 
!        RV1 is a one-dimensional REAL array used for temporary storage,
!          dimensioned RV1(N).
! 
!     CALLS Hypot(A,B) for sqrt(A**2 + B**2).
! 
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
! 
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  Hypot
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   991125  Revised the subroutine arguments
!***END PROLOGUE  SVD
! 
      INTEGER:: i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its
!      REAL a(nm,*),w(*),u(nm,*),v(nm,*),rv1(*)
      REAL(SP):: c
      REAL(SP):: f,g,h
      REAL(SP):: s,s1
      REAL(SP):: scale
      REAL(SP):: x,y,z
      LOGICAL:: matu,matv
  REAL(SP),PARAMETER:: ZERO=0.0, HALF=0.5, ONE=1.0

  REAL(SP),ALLOCATABLE,DIMENSION(:):: rv1
!----------------------------------------------------------------------------
  m=SIZE(a,1)
  n=SIZE(a,2)

  ierr=0
  IF (SIZE(u,1) < m) ierr=-1
  IF (SIZE(u,2) < n) ierr=-2
  IF (SIZE(w) < n) ierr=-3
  IF (SIZE(v,1) < n) ierr=-4
  IF (SIZE(v,2) < n) ierr=-5
  IF (ierr < 0) RETURN

  ALLOCATE(rv1(n))

  u(1:m,1:n)=a
!      DO 10 i=1,m
!        DO 10 j=1,n
!        u(i,j)=a(i,j)
!   10 CONTINUE

!     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
      g=ZERO
      scale=ZERO
      s1=ZERO
 
      DO 170 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=ZERO
        s=ZERO
        scale=ZERO
        IF (i.gt.m) GO TO 80
 
        DO k=i,m
          scale=scale+ABS(u(k,i))                  ! scale=SUM(ABS(u(i:m,i)))
        END DO
 
        IF (scale==ZERO) GO TO 80
 
        DO k=i,m
          u(k,i)=u(k,i)/scale                    ! u(i:m,i)=u(i:m,i)/scale
          s=s+u(k,i)**2                          ! s=SUM(u(i:m,i)**2)
        END DO   
 
        f=u(i,i)
        g=-SIGN(SQRT(s),f)
        h=f*g-s
        u(i,i)=f-g
        IF (i==n) GO TO 60
 
        DO j=l,n
          s=ZERO 
          DO k=i,m            ! DOT_PRODUCT(u(i:m,i),u(i:m,j))
            s=s+u(k,i)*u(k,j)
          END DO
 
          f=s/h
 
          DO k=i,m
            u(k,j)=u(k,j)+f*u(k,i)   ! u(i:m,j)=u(i:m,j)+f*u(i:m,i)
          END DO
        END DO   
 
   60   DO k=i,m
          u(k,i)=scale*u(k,i)        ! u(i:m,i)=u(i:m,i)*scale
        END DO
 
   80   w(i)=scale*g
        g=ZERO
        s=ZERO
        scale=ZERO
        IF (i.gt.m .OR. i.eq.n) GO TO 160
 
        DO k=l,n
          scale=scale+ABS(u(i,k))                  ! scale=SUM(ABS(u(i,1:n)))
        END DO
 
        IF (scale==ZERO) GO TO 160
 
        DO k=l,n
          u(i,k)=u(i,k)/scale
          s=s+u(i,k)**2
        END DO  
 
        f=u(i,l)
        g=-SIGN(SQRT(s),f)
        h=f*g-s
        u(i,l)=f-g
 
        DO k=l,n
          rv1(k)=u(i,k)/h                               ! rv1(L:n)=u(i,L:n)/h
        END DO
 
        IF (i.eq.m) GO TO 140
 
        DO j=l,m
          s=ZERO
 
          DO k=l,n
            s=s+u(j,k)*u(i,k)
          END DO
 
          DO k=l,n
            u(j,k)=u(j,k)+s*rv1(k)
          END DO
        END DO  
! 
  140   DO k=l,n
          u(i,k)=scale*u(i,k)
        END DO
! 
  160   s1=MAX(s1,ABS(w(i))+ABS(rv1(i)))
  170 END DO


!     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS ..........
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 240 ii=1,n
        i=n+1-ii
        IF (i.eq.n) GO TO 230
        IF (g==ZERO) GO TO 210
 
        DO j=l,n
          v(j,i)=(u(i,j)/u(i,l))/g  !  DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
        END DO
 
        DO j=l,n
          s=ZERO 
          DO k=l,n
            s=s+u(i,k)*v(k,j)
          END DO
 
          DO k=l,n
            v(k,j)=v(k,j)+s*v(k,i)
          END DO
        END DO   
 
  210   DO j=l,n
          v(i,j)=ZERO
          v(j,i)=ZERO
        END DO  
 
  230   v(i,i)=ONE
        g=rv1(i)
        l=i
  240 END DO

!     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ..........

!     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- ..........
      mn=MIN(m,n)
 
      DO 350 ii=1,mn
        i=mn+1-ii
        l=i+1
        g=w(i)
        IF (i==n) GO TO 270
 
        DO j=l,n
          u(i,j)=ZERO
        END DO
 
  270   IF (g==ZERO) GO TO 320
        IF (i==mn) GO TO 300
 
        DO j=l,n
          s=ZERO
          DO k=l,m
            s=s+u(k,i)*u(k,j)
          END DO

          f=(s/u(i,i))/g   ! DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW 
 
          DO k=i,m
            u(k,j)=u(k,j)+f*u(k,i)
          END DO
        END DO  
 
  300   DO j=i,m
          u(j,i)=u(j,i)/g
        END DO
 
        GO TO 340
 
  320   DO j=i,m
          u(j,i)=ZERO
        END DO
 
  340   u(i,i)=u(i,i)+ONE
  350 END DO

!     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
!!!  360 CONTINUE

!     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
      DO 500 kk=1,n
        k1=n-kk
        k=k1+1
        its=0
!     .......... TEST FOR SPLITTING.
!                FOR L=K STEP -1 UNTIL 1 DO -- ..........
  370   DO ll=1,k
          l1=k-ll
          l=l1+1
          IF (s1+ABS(rv1(l))==s1) GO TO 420
!     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
          IF (s1+ABS(w(l1))==s1) GO TO 390
        END DO  
!     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 ..........
  390   c=ZERO
        s=ONE
 
        DO i=l,k
          f=s*rv1(i)
          rv1(i)=c*rv1(i)
          IF (s1+ABS(f)==s1) GO TO 420
          g=w(i)
          h=Hypot(f,g)
          w(i)=h
          c=g/h
          s=-f/h
 
          DO j=1,m
            y=u(j,l1)
            z=u(j,i)
            u(j,l1)=y*c+z*s
            u(j,i)=-y*s+z*c
          END DO  
        END DO  

!     .......... TEST FOR CONVERGENCE ..........
  420   z=w(k)
        IF (l==k) GO TO 480
!     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
        IF (its==30) GO TO 510
        its=its+1
        x=w(l)
        y=w(k1)
        g=rv1(k1)
        h=rv1(k)
        f=HALF*(((g+z)/h)*((g-z)/y)+y/h-h/y)
        g=Hypot(f,ONE)
        f=x-(z/x)*z+(h/x)*(y/(f+sign(g,f))-h)

!     .......... NEXT QR TRANSFORMATION ..........
        c=ONE
        s=ONE
 
        DO 470 i1=l,k1
          i=i1+1
          g=rv1(i)
          y=w(i)
          h=s*g
          g=c*g
          z=Hypot(f,h)
          rv1(i1)=z
          c=f/z
          s=h/z
          f=x*c+g*s
          g=-x*s+g*c
          h=y*s
          y=y*c
 
          DO j=1,n
            x=v(j,i1)
            z=v(j,i)
            v(j,i1)=x*c+z*s
            v(j,i)=-x*s+z*c
          END DO  
! 
          z=Hypot(f,h)
          w(i1)=z

!     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
          IF (z==ZERO) GO TO 450
          c=f/z
          s=h/z
  450     f=c*g+s*y
          x=-s*g+c*y
 
          DO j=1,m
            y=u(j,i1)
            z=u(j,i)
            u(j,i1)=y*c+z*s
            u(j,i)=-y*s+z*c
          END DO  
 
   470  END DO
 
        rv1(l)=ZERO
        rv1(k)=f
        w(k)=x
        GO TO 370

!     .......... CONVERGENCE ..........
  480   IF (z < ZERO) THEN
          w(k)=-z   !     .......... w(k) is made non-negative ..........
          DO j=1,n
            v(j,k)=-v(j,k)
          END DO
        END IF
  500 END DO
 
      GO TO 520
!     .......... SET ERROR -- NO CONVERGENCE TO A
!                SINGULAR VALUE AFTER 30 ITERATIONS ..........
  510 ierr=k

  520 DEALLOCATE(rv1)
  RETURN
END Subroutine Svd   ! ------------------------------------------------------

!+
SUBROUTINE SVDbackSubstitution(u,w,v,b,x)
! ---------------------------------------------------------------------------
! PURPOSE -
IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:,:):: u
  REAL,INTENT(IN),DIMENSION(:)::   w
  REAL,INTENT(IN),DIMENSION(:,:):: v
  REAL,INTENT(IN),DIMENSION(:)::   b
  REAL,INTENT(OUT),DIMENSION(:)::  x

  REAL(SP),DIMENSION(SIZE(x)):: tmp
  REAL(SP),PARAMETER:: ZERO=0.0
!----------------------------------------------------------------------------
  WHERE (w/=ZERO)
    tmp=MATMUL(b,u)/w
  ELSEWHERE
    tmp=ZERO
  END WHERE

  x=MATMUL(v,tmp)
  RETURN
END Subroutine SVDbackSubstitution   ! ----------------------------------------

  
!+
PURE FUNCTION Hypot(x,y) RESULT(z)   ! a PURE function
! ---------------------------------------------------------------------------
! PURPOSE - Solve for SQRT(x*x+y*y) carefully to avoid overflow
!  (will be an intrinsic in Fortran 2008)
  REAL,INTENT(IN):: x,y
  REAL:: z

  REAL:: a,b
!----------------------------------------------------------------------------
  a=ABS(x)
  b=ABS(y)

  IF (a > b) THEN
    z=a*SQRT(1.0+(b/a)**2)
  ELSE IF (b==0) THEN
    z=0.0
  ELSE
    z=b*SQRT(1.0+(a/b)**2)
  END IF

  RETURN
END Function Hypot   ! -----------------------------------------------------

END Module Smooth   ! =======================================================

!+
MODULE ContourPlot
! ---------------------------------------------------------------------------
! PURPOSE
! AUTHORS -

USE Smooth
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: CONPLOT_VERSION= "1.0 (5Oct99)"
  INTEGER,PARAMETER,PRIVATE:: SP=KIND(1.0), DP=KIND(1.0D0)

  PUBLIC:: ContourLines
  PRIVATE:: Interpolate,Middle,Triangulate

INTERFACE
  SUBROUTINE CNTCRV (xx,yy,npoint,zcon)
    IMPLICIT NONE
    REAL,INTENT(IN),DIMENSION(:):: xx,yy
    INTEGER,INTENT(IN):: npoint
    REAL,INTENT(IN):: zcon
  END Subroutine CntCrv
END INTERFACE

CONTAINS
  
!+
SUBROUTINE ContourLines(x,y,z,ismopt,iexp,jexp,clist,epslon,ierr)
! ---------------------------------------------------------------------------
! PURPOSE - Driver program for computing and drawing contour lines of
!   constant z for the function z = f(x,y).

IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: x  ! input list of x values
  REAL,INTENT(IN),DIMENSION(:):: y  ! input list of y values
  REAL,INTENT(IN),DIMENSION(:):: z  ! input list of z values
!  INTEGER,INTENT(IN):: n  !   number of values in x,y and z
  INTEGER,INTENT(IN):: ismopt  ! smoothing option flag (0=no/off, 1=yes/on)
  INTEGER,INTENT(IN):: iexp !   = i exponent value for smoothing
  INTEGER,INTENT(IN):: jexp !   =j exponent value for smoothing
!  INTEGER,INTENT(IN):: ncntrs  !  number of contour lines to be drawn
!                  (self computing if ncntrs.le.0)
  REAL,INTENT(IN),DIMENSION(:):: clist  ! = list of constant contour values
  REAL,INTENT(OUT):: epslon  !  error function (normalized value) returned to
!                  caller if ismopt is non-zero
  INTEGER,INTENT(OUT):: ierr !   = return error flag
!                  = 0 for normal return
!                  = 1 for invalid value for n
!                  = 2 for number of ismopt coefficients greater
!                      than 'maxcof' or n
!

  REAL(SP),DIMENSION(23):: coef
  INTEGER:: j,k, i1,j1
  INTEGER,DIMENSION(1494):: ibe,lambda
  INTEGER,DIMENSION(1494,2):: ie
  INTEGER,DIMENSION(1494,4):: ite
  INTEGER,DIMENSION(23):: ipowr,jpowr
  INTEGER:: ledges
  INTEGER,PARAMETER:: MAXCOF=23, MAXPTS=500
  INTEGER:: n
  INTEGER:: nc
  INTEGER:: ncoef=0
!  INTEGER:: nmax,nmin
  REAL(SP),DIMENSION(1494):: xi,eta
  REAL(SP):: zcon, zmin,zmax
  REAL(SP),DIMENSION(500):: znew

!-----------------------------------------------------------------------
!!! allocate  ibe,ie,ite,lambda,xi,eta,znew

! (A) initialize local variables and check inputs for errors! 
 
      ierr=0
      epslon=0.0
      n=MIN(SIZE(x), SIZE(y), SIZE(z))
      IF (n<3 .or. n>maxpts) THEN
        ierr=1
        RETURN
      END IF
 
 
! (B) triangulate x-y data points
 
      CALL Triangulate(x, y, n, ledges, ie, ibe, ite)

! (C) smooth data if requested
 
      IF (ismopt==0) THEN
        znew(1:n)=z(1:n)
      ELSE 
        i1=iexp+1   ! (D)  check requested exponent values for errors
        j1=jexp+1
!        nmin=MIN(i1,j1)
!        nmax=MAX(i1,j1)
        IF (jexp < iexp) THEN
          nc=(jexp+1)*(iexp+1-jexp/2)
        ELSE
          nc=(iexp+1)*(jexp+1-iexp/2)
        END IF
        
        IF (nc>n .or. nc>maxcof) THEN
          znew(1:n)=z(1:n)
          ierr=2
          RETURN
        END IF

        ipowr(:)=0
        jpowr(:)=0 
        CALL SmoothSurface(x,y,z, znew(1:n), iexp,jexp, &
           ncoef,coef, ipowr,jpowr)
        IF (ncoef.lt.0) GO TO 50
        DO k=1,n
          epslon=epslon+(z(k)-znew(k))**2
        END DO  
        epslon=sqrt(epslon)/REAL(n)
      END IF
! 

! (F) determine the range of the z data under consideration

   50 zmin=MINVAL(z(1:n))
      zmax=MAXVAL(z(1:n))
 
! (J) interpolate for contour line data points
! (K,L) any data points found? . .
!   call subroutine cntour to sort the interpolated points
!   on the contour line and draw it

  DO k=1,SIZE(clist)
    zcon=clist(k)
    IF (zcon < zmin) CYCLE
    IF (zcon > zmax) CYCLE
    CALL Interpolate(x, y, znew, zcon, ledges, ie, ismopt, lambda, &
      xi, eta, j, coef, ipowr, jpowr, ncoef)
      IF (j>0) CALL cntour (zcon, xi, eta, lambda, j, ibe, ite)
  END DO
 
  RETURN
END Subroutine ContourLines   ! ---------------------------------------------

!+
      SUBROUTINE Triangulate(xd, yd, n, l, e, be, te) 
! ------------------------------------------------------------------
! PURPOSE - A set of n data points are known (x(i),y(i),i=1,n)  They are to
!   be connected by lines to form a set of triangles (for n.LE. maxpts).
!   the final triangulation establishes a convex polygon defined by linked
!   lists of edge numbers, end points and boundary edges.
!
 IMPLICIT NONE
!      IMPLICIT INTEGER (P,B)   ! oh dear
!
!
!     SUBROUTINE INPUT
  REAL,INTENT(IN),DIMENSION(:):: xd !  = array of abscissas
  REAL,INTENT(IN),DIMENSION(:):: yd !  = array of ordinates
  INTEGER,INTENT(IN):: n  ! = number of points in x and y
!
!     SUBROUTINE OUTPUT
  INTEGER,INTENT(OUT):: l  ! = number of edges listed in e, be and te
  INTEGER,INTENT(OUT),DIMENSION(:,:):: e  ! = list of indices of each traingle edge
  INTEGER,INTENT(OUT),DIMENSION(:):: be  ! = 1 if i of e is a boundary edge
  INTEGER,INTENT(OUT),DIMENSION(:,:):: te  != indices of neighboring edges for each traingle
!
!     LOCAL VARIABLES

  INTEGER,DIMENSION(500):: b   ! index of points on the boundary ..  inorder
  INTEGER:: b1,b2
  INTEGER:: bflag
  INTEGER:: bj
  REAL(DP):: d,d1
!  REAL:: dlxinv,dlyinv
  REAL:: dmin
  REAL:: dst
  INTEGER:: i1,i2
  INTEGER:: j   ! no. of values in list p
  INTEGER:: j1
  INTEGER:: jcnt
  INTEGER:: k   ! no. of points listed in array b
  INTEGER:: k0,k1,k2,k3
  INTEGER:: kcnt
  INTEGER:: kk
  INTEGER:: kknt
  INTEGER:: kl
  INTEGER:: km
  INTEGER:: kn
  INTEGER:: kp1
  INTEGER:: kt
  INTEGER:: l1,l2,l3
  INTEGER:: lambda
  INTEGER:: lcnt
  INTEGER:: lj,ll
  INTEGER:: m   ! no. of rows used in array t
  INTEGER:: mcnt
  INTEGER:: nq
  INTEGER,DIMENSION(500):: p   ! indices of points outside the boundary
  INTEGER:: p1,p2
  INTEGER:: pj
  INTEGER:: pk3,pk0,pkl
  REAL(DP):: s
  INTEGER,DIMENSION(995,3):: t   ! indices of adjacent triangle edges
  REAL(DP):: tc
  REAL(DP):: term

  REAL,DIMENSION(500):: x,y   ! arrays of scaled data
  REAL(DP):: x12,y12
  REAL:: x1b,y1b
  REAL:: xmax,xmin, ymax,ymin
  REAL:: xqb,yqb

!----------------------------------------------------------------------------

! (A) The procedure begins with no boundary, no edges, and all x-y data
!   points under consideration. Scale the x,y data and initialize local
!   variables.
 
      j=n
      k=0
      l=0
      m=0
      kknt=0
      DO jcnt=1,SIZE(x)
        p(jcnt)=jcnt
      END DO
      xmax=MAXVAL(xd(1:n))
      xmin=MINVAL(xd(1:n))
      ymax=MAXVAL(yd(1:n))
      ymin=MINVAL(yd(1:n))
      x(1:n)=xd(1:n)/(xmax-xmin)
      y(1:n)=yd(1:n)/(ymax-ymin)

!      dlxinv=1.0/(xmax-xmin)
!      dlyinv=1.0/(ymax-ymin)
!      DO 30 k=1,n
!        x(k)=xd(k)*dlxinv
!        y(k)=yd(k)*dlyinv
!   30 CONTINUE
! 
! 
! (B) Begin by taking the last pair of points (x(j),y(j)) in the
!   list to be the first boundary point.
 
      b(1)=j
      j=j-1
! 
! 
! (C) from the remaining points, find the point nearest the first
  
      i2=1
      i1=b(1)
      dmin=(x(i1)-x(1))**2+(y(i1)-y(1))**2
      DO j1=2,j
        dst=(x(i1)-x(j1))**2+(y(i1)-y(j1))**2
        IF (dst.ge.dmin) CYCLE
        i2=j1
        dmin=dst
      END DO   
! 
! (D) now b(1) to b(i2) is the first edge.
!   there is one edge and two boundary points.

! 
      j=j-1
      DO jcnt=i2,j
        p(jcnt)=p(jcnt+1)
      END DO  
      k=2
      b(2)=i2
      l=1
      e(1,1)=MIN(b(1),b(2))
      e(1,2)=MAX(b(1),b(2))
! 
! 
! (E) now begin circling around the boundary of the polygon, considering,
!   in order, each boundary edge.  maintain the following indices -
!     k1 = b array index of the current edge - point 1
!     k2 = b array index of the current edge - point 2
!     b1,b2 = indices of boundary point coordinates
 
      k1=0
      kt=0
   70 k1=k1+1
      IF (k1.gt.k) k1=1
   80 k2=k1+1
      IF (k2.gt.k) k2=1
      b1=b(k1)
      b2=b(k2)
      kt=kt+1
! 
! (F) consider the boundary edge from b1 to b2.  for all points not yet
!     triangulated (the j points remaining in p), find the point that,
!     when traingulated with b1,b2, minimizes the length of the two new
!     edges to be drawn.


      d1=0.
      j1=0
      bflag=0
      
      DO lj=1,j
        pj=p(lj)
        term=(y(pj)-y(b1))*(x(b2)-x(b1))-(x(pj)-x(b1))*(y(b2)-y(b1))
        IF (term.le.0.) CYCLE   
        d=sqrt((x(pj)-x(b1))**2+(y(pj)-y(b1))**2)+sqrt((x(pj)-x(b2))**2+ &
          (y(pj)-y(b2))**2)
        IF (j1.ne.0.and.d1.lt.d) CYCLE   
        j1=lj
        d1=d
      END DO  
! 
!        (G)
!        if less than three edges exist (no triangle defined yet),
!        then there are no adjacent boundary points to be considered.
!        so go to section J.
! 
      IF (k.le.3) GO TO 120
! 
! 
!        (H)
!        consider the adjacent boundary point of the next edge of the
!        polygon.  call its index number k3 and see if its closer to
!        the current edge than p(j1).
! 
      k3=k2+1   ! was labelled 100
      IF (k3.gt.k) k3=1
      pk3=b(k3)
      term=(y(pk3)-y(b1))*(x(b2)-x(b1))-(x(pk3)-x(b1))*(y(b2)-y(b1))
      IF (term.le.0.) GO TO 110
      d=sqrt((x(pk3)-x(b1))**2+(y(pk3)-y(b1))**2) + &
        sqrt((x(pk3)-x(b2))**2+(y(pk3)-y(b2))**2)
      IF (j1.ne.0.and.d1.lt.d) GO TO 110
      j1=k3
      d1=d
      bflag=1
! 
!        (I)
!        consider the adjacent boundary point of the previous edge of
!        the polygon.  call its index number k0 and see if its closer
!        to the current edge than p(j1) and b(k3).
! 
  110 k0=k1-1
      IF (k0.lt.1) k0=k
      pk0=b(k0)
      term=(y(pk0)-y(b1))*(x(b2)-x(b1))-(x(pk0)-x(b1))*(y(b2)-y(b1))
      IF (term.le.0.) GO TO 120
      d=sqrt((x(pk0)-x(b1))**2+(y(pk0)-y(b1))**2) + &
        sqrt((x(pk0)-x(b2))**2+(y(pk0)-y(b2))**2)
      IF (j1.ne.0.and.d1.lt.d) GO TO 120
      j1=k0
      d1=d
      bflag=-1

!        (J)
!        skip the next section if j1 is still zero, since a candidate
!        point for triangulation with edge b1,b2 was not found.
 
  120 IF (j1.eq.0) GO TO 130
  
! (K,L) If the search for a candidate point has already considered each
! boundary edge at least once (kt.gt.k) or if the boundary is being
!  checked for concave edges (j=0), then the next section
!        (section m) can be ommitted.
! 
      IF (kt.gt.k.or.j.eq.0) GO TO 130
! 
!        (M)
!        at this point the user may insert any additional constraint
!        on the triangle to be formed by the point pj1.  If the
!        candidate triangle fails the test, it is deleted from
!        consideration by setting the variable j1 to zero.
! 
!!!  130 CONTINUE
! 
! 
!        (N,O)
!        the next procedure checks all boundary edges of the polygon
!        for intersection with the candidate triangle.  If any existing
!        boundary edge intersects any of the edges to be formed by the
!        candidate triangle, then the candidate point is rejected.  If
!        bflag is not zero, then the edge defined by j1=k0 or j1=k3 is
!        exempt form this test.
! 
!        If there are three or less existing boundary edges or if
!        j1 has been set to zero, this test is ommitted.
! 
  130 IF (k.le.3.or.j1.eq.0) GO TO 160
      IF (bflag.eq.0) nq=p(j1)
      IF (bflag.eq.1) nq=b(k3)
      IF (bflag.eq.-1) nq=b(k0)
      DO kcnt=1,k
        IF (kcnt.eq.k1) GO TO 150
        kn=kcnt+1
        IF (kcnt.eq.k) kn=1
        IF (bflag.eq.-1.and.(kcnt.eq.k0.or.kn.eq.k0)) GO TO 150
        IF (bflag.eq.1.and.(kcnt.eq.k3.or.kn.eq.k3)) GO TO 150
        p1=b(kcnt)
        p2=b(kn)
        DO jcnt=1,2
          IF (jcnt.eq.1.and.(bflag.eq.0.or.bflag.eq.1).and.kcnt.eq.k0) &
     &     GO TO 150
          IF (jcnt.eq.2.and.(bflag.eq.0.or.bflag.eq.-1).and.kcnt.eq.k2) &
     &     GO TO 150
          bj=b1
          IF (jcnt.eq.2) bj=b2
          xqb=x(nq)-x(bj)
          yqb=y(nq)-y(bj)
          x12=x(p1)-x(p2)
          y12=y(p1)-y(p2)
          d=xqb*y12-yqb*x12
          IF (d.eq.0.) CYCLE
          x1b=x(p1)-x(bj)
          y1b=y(p1)-y(bj)
          s=(x1b*y12-y1b*x12)/d
          IF (s.lt.0..or.s.gt.1.) CYCLE
          tc=(xqb*y1b-yqb*x1b)/d
          IF (tc.lt.0.or.tc.gt.1.) CYCLE
          j1=0
          GO TO 160
    END DO       ! was 140
    150 b1=b1
  END DO         ! was 150
! 
! 
!        (P,Q)
!        If j1 is zero, then the candidate point did not pass the above
!        tests or no point was found.  If bflag is not zero, then a
!        point on the boundary was found.
! 
 160  IF (j1.eq.0) GO TO 270
!!!      IF (bflag) 240,170,210
         IF (bflag > 0) GO TO 210
         IF (bflag==0) GO TO 170
         IF (bflag < 0) GO TO 240
! 
! 
!        The triangulated point is outside the boundary.  Establish two
!        new edges, a new boundary point and delete one point from
!        outside the boundary.
! 
! 
  170 e(l+1,1)=MIN(p(j1),b(k1))
      e(l+1,2)=MAX(p(j1),b(k1))
      e(l+2,1)=MIN(p(j1),b(k2))
      e(l+2,2)=MAX(p(j1),b(k2))
      kt=0
      l=l+2
      m=m+1
      t(m,1)=MIN(p(j1),b(k1),b(k2))
      t(m,2)=middle(p(j1),b(k1),b(k2))
      t(m,3)=MAX(p(j1),b(k1),b(k2))
      IF (k1.eq.k) GO TO 190
      km=k
      kp1=k1+1

  180 b(km+1)=b(km)
      km=km-1
      IF (km.ge.kp1) GO TO 180

  190 b(k1+1)=p(j1)
      k=k+1
      j=j-1
      IF (j1.gt.j) GO TO 270
      DO jcnt=j1,j
        p(jcnt)=p(jcnt+1)
      END DO
      GO TO 270
! 
! 
!        (S)
!        The triangulated point is the next point on the boundary.
!        Establish one new edge (from b(k1) to b(k3)), one new
!        triangle (from b(k1) to b(k2) to b(k3)), and delete one point
!        from the boundary (b(k2)).
! 
! 
  210 e(l+1,1)=MIN(b(k3),b(k1))
      e(l+1,2)=MAX(b(k3),b(k1))
      kk=0
      kknt=0
      kt=0
      l=l+1
      k=k-1
      m=m+1
      t(m,1)=MIN(b(k1),b(k2),b(k3))
      t(m,2)=middle(b(k1),b(k2),b(k3))
      t(m,3)=MAX(b(k1),b(k2),b(k3))

      DO kcnt=k2,k
        b(kcnt)=b(kcnt+1)
      END DO

  230 IF (k2.eq.1) k1=k1-1
      GO TO 270
! 
!        (R)
!        The triangulated point is the previous point on the boundary.
!        Establish a new edge (from b(k0) to b(k2)), one new triangle
!        (from b(k0) to b(k1) to b(k2)), and delete one point from the
!        boundary (b(k1))
! 
  240 e(l+1,1)=MIN(b(k0),b(k2))
      e(l+1,2)=MAX(b(k0),b(k2))
      kk=0
      kknt=0
      kt=0
      l=l+1
      k=k-1
      m=m+1
      t(m,1)=MIN(b(k0),b(k1),b(k2))
      t(m,2)=Middle(b(k0),b(k1),b(k2))
      t(m,3)=MAX(b(k0),b(k1),b(k2))

      DO kcnt=k1,k
        b(kcnt)=b(kcnt+1)
      END DO
  260 k1=k1-1
      IF (k1.lt.1) k1=k

! (T) If there are any points remaining outside the boundary, then
! repeat the procedure for the next edge.
 
  270 IF (j.gt.0.and.j1.ne.0) GO TO 80
      IF (j.gt.0) GO TO 70
 
 
! (U,V,W) All points have been triangulated.  Check that all boundary
! edges form a concave polygon.

      IF (kk.ne.0) GO TO 280
      kk=1
      kl=0
  280 kknt=kknt+1
      IF (kknt.gt.n) GO TO 300
  290 kl=kl+1
      k2=kl+1
      IF (k2.gt.k) k2=1
      k1=kl-1
      IF (k1.lt.1) k1=k
      pkl=b(kl)
      b1=b(k1)
      b2=b(k2)
      term=(y(pkl)-y(b1))*(x(b2)-x(b1))-(x(pkl)-x(b1))*(y(b2)-y(b1))
      IF (term.lt.0.) GO TO 70
      IF (kl.lt.k) GO TO 290
 
! (X) The triangulation is complete and has been checked for a
! concave boundary.  Now identify the boundary edges.

  300 DO 340 lcnt=1,l
        be(lcnt)=0
        kl=0
  310   kl=kl+1
        IF (e(lcnt,1).ne.b(kl)) GO TO 330
        k1=kl+1
        IF (k1.gt.k) k1=1
        IF (e(lcnt,2).ne.b(k1)) GO TO 320
        be(lcnt)=1
        CYCLE   !!! GO TO 340
  320   k1=kl-1
        IF (k1.lt.1) k1=k
        IF (e(lcnt,2).ne.b(k1)) GO TO 330
        be(lcnt)=1
        CYCLE   !!! GO TO 340
  330   IF (kl.lt.k) GO TO 310
 340  END DO   
 
! (Y) Finally, establish the indices of adjacent edges for each edge
! in the triangulation.  Each boundary edge will have two adjacent edges -
! each interior edge will have four.

       te(1:L,1:4)=0
!      DO 350 ll=1,4
!        DO 350 lcnt=1,l
!  350   te(lcnt,ll)=0

      DO mcnt=1,m
        DO ll=1,l
          IF (e(ll,1).eq.t(mcnt,1).and.e(ll,2).eq.t(mcnt,2)) l1=ll
          IF (e(ll,1).eq.t(mcnt,2).and.e(ll,2).eq.t(mcnt,3)) l2=ll
          IF (e(ll,1).eq.t(mcnt,1).and.e(ll,2).eq.t(mcnt,3)) l3=ll
        END DO  
        lambda=0
        IF (te(l1,1).ne.0) lambda=2
        te(l1,lambda+1)=l2
        te(l1,lambda+2)=l3
        lambda=0
        IF (te(l2,1).ne.0) lambda=2
        te(l2,lambda+1)=l1
        te(l2,lambda+2)=l3
        lambda=0
        IF (te(l3,1).ne.0) lambda=2
        te(l3,lambda+1)=l1
        te(l3,lambda+2)=l2
      END DO   
! 
      RETURN
      END Subroutine Triangulate



!+
      FUNCTION Middle (i, j, k) RESULT(m)
! 
! 
!        This function subprogram is used by the triangulation algorithm
!        to find the middle value of the three integer arguments (the
!        value which is neither a minimum or a maximum).  i, j and k are
!        are assumed to be discrete values with no two equal.
!
  INTEGER,INTENT(IN):: i,j,k

  INTEGER:: m
!----------------------------------------------------------------------------
      IF (j.lt.i.and.i.lt.k) GO TO 10
      IF (k.lt.i.and.i.lt.j) GO TO 10
      IF (i.lt.j.and.j.lt.k) GO TO 20
      IF (k.lt.j.and.j.lt.i) GO TO 20
      m=k
      RETURN
   10 m=i
      RETURN
   20 m=j
      RETURN
      END Function Middle

!+
      SUBROUTINE Interpolate (x, y, u, zcon, ledges, ie, ismopt,   &
        lambda, xi, eta, j, c, ipowr, jpowr, ncoef)
! ---------------------------------------------------------------------------
! PURPOSE - Subroutine Interpolate is given a constant u value (bigu)
!  for which the contour line is to be drawn.  Check all given triangle
!  edges, (array ie) and check the values of u at the endpoints.
!  Interpolate for all possible values on the triangle edges.
!  If ismopt = 0, then use a linear interpolation, if ismopt not zero,
!  then evaluate for a non-linear surface using the coefficients
!  from smsrf and function subroutine Polyx2.
!
  REAL,INTENT(IN),DIMENSION(:):: x,y,u   ! dependent and independent values
                                         ! for the relation u=f(x,y)

!  INTEGER,INTENT(IN):: n
  REAL,INTENT(IN):: zcon ! CONSTANT VALUE OF Z FOR WHICH INTERPOLATION
                         ! IS REQUIRED

  INTEGER,INTENT(IN):: ledges   ! NO. OF EDGES IN THE TRIANGULATION
  INTEGER,INTENT(IN),DIMENSION(:,:):: ie   ! edge endpoint indices from
                                           ! triangulation
  INTEGER,INTENT(IN):: ismopt   ! smoothing option flag,  0=off, 1=on
  INTEGER,INTENT(OUT),DIMENSION(:):: lambda   ! index of edges for interpolated points
  REAL,INTENT(OUT),DIMENSION(:):: xi   ! x-coordinates of interpolated points
  REAL,INTENT(OUT),DIMENSION(:):: eta  ! y-coordinates of interpolated points
  INTEGER,INTENT(OUT):: j  ! length of xi and eta
  REAL,INTENT(IN),DIMENSION(:):: c   ! coefficients of each term of the equation.
  INTEGER,INTENT(IN),DIMENSION(:):: ipowr,jpowr
  INTEGER,INTENT(IN):: ncoef


!        ipowr,jpowr are the list of exponents for each term of
!        the polynomial used to smooth the data  (input).
!        ncoef     = number of terms in the polynomial
!                    (ipowr,jpowr,c, and ncoef are input)
! 
! 
!     ------------------------------------------------------------------
! 
! 
!      DIMENSION x(n), y(n), u(n)
!      DIMENSION ie(1494,2), xi(1494), eta(1494), lambda(1494)
!      DIMENSION ipowr(23), jpowr(23), c(23)

!
  REAL:: f1,fn
  INTEGER:: i1,i2
  INTEGER:: k
  INTEGER:: lcnt
  INTEGER:: n
  LOGICAL:: smoothed
  REAL:: t1,t2
  REAL:: temp
  REAL:: u1,u2
  REAL:: x1,x2
  REAL:: xn,yn
  REAL:: y1,y2
!----------------------------------------------------------------------------
      n=MIN(SIZE(x), SIZE(y), SIZE(u))
      smoothed=ismopt/=0
      IF (ncoef < 1) smoothed=.FALSE.

!!!      IF (ncoef.lt.1) ismopt=0   ! deal with this
      j=0
! 
      DO 70 lcnt=1,ledges
! 
!        (A)
!        Determine x,y,z for the endpoints of the nextedge - order them
! 
        i1=ie(lcnt,1)
        i2=ie(lcnt,2)
        x1=x(i1)
        x2=x(i2)
        y1=y(i1)
        y2=y(i2)
        u1=u(i1)
        u2=u(i2)
! 
!        (B)
!        Function values equal at endpoints or
!        constant zc not between them? . .
! 
        IF (u1.eq.u2) CYCLE   !!! GO TO 70
        IF (u1.lt.u2) GO TO 10
        temp=u2
        u2=u1
        u1=temp
        temp=x2
        x2=x1
        x1=temp
        temp=y2
        y2=y1
        y1=temp

   10   IF (zcon.lt.u1.or.u2.lt.zcon) CYCLE   !!! GO TO 70
        IF (u2.eq.zcon) u2=1.000001*zcon
        j=j+1
! 
!        (C)
!        Has data been smoothed? . .
!        If not, goto section e (statement label 101)
! 
        IF (smoothed) THEN
!         (D,F) Non-linear interpolation is required on this edge
          f1=Polyx2(zcon,x1,y1,c,ipowr,jpowr,ncoef)
          DO k=1,10
            xn=0.5*(x1+x2)
            yn=0.5*(y1+y2)
            fn=Polyx2(zcon,xn,yn,c,ipowr,jpowr,ncoef)
            IF (fn==0.0) EXIT
            IF (fn.lt.0..and.f1.lt.0.) THEN
              x1=xn
              y1=yn
            ELSE IF (fn.gt.0..and.f1.gt.0.) THEN
              x1=xn
              y1=yn
            ELSE
              x2=xn
              y2=yn
            END IF
          END DO 
          xi(j)=(x1+x2)*0.5
          eta(j)=(y1+y2)*0.5
        ELSE
!         (E,F) Linear interpolation is required for this edge
          t1=(u2-zcon)/(u2-u1)
          t2=1.0-t1
          xi(j)=t1*x1+t2*x2
          eta(j)=t1*y1+t2*y2
        END IF
        lambda(j)=lcnt
   70 END DO

      RETURN
END Subroutine Interpolate   ! ----------------------------------------------

!+
SUBROUTINE Cntour (zcon, xi, eta, lambda, j, ibe, ite)
! ---------------------------------------------------------------------------
! PURPOSE - 
!     A set of j interpolated points for z=zcon  (xi(i),eta(i) on edge
!     lambda(i) for i=1,j), the contour lines must now be drawn.  There
!     may be several lines, either open or closed contours.  This
!     algorithm will use the triangulation relationships to sort out
!     each line in order.  As each contour line is established, user
!     supplied program CNTCRV is called to output it to the graphics
!     device being used.
! 
! 
!     ARGUMENTS (ALL ARE INPUTS) -
!        ZCON    = CONSTANT VALUE OF Z UNDER CONSIDERATION
!        XI(J)   = ARRAY OF X COORDIANTES OF INTERPOLATED POINTS
!        ETA(J)  = ARRAY OF Y COORDIANTES OF INTERPOLATED POINTS
!        LAMBDA(J) = ARRAY OF EDGE NUMBERS FOR J-TH INTERPOLATED POINT
!        J       = NUMBER OF POINTS IN THE LIST OF INTERPOLATED POINTS
!        IBE     = THE LIST OF BOUNDARY EDGES TAKEN FROM THE TRIANGULATI
!        ITE     = LINKED LIST OF INDICES OF ADJACENT EDGES PROVIDED
!                  BY THE TRIANGULATION PROCEDURE.
 

  REAL,INTENT(IN):: zcon                         ! z-value of the contour line
  REAL,INTENT(IN OUT),DIMENSION(:):: xi          ! x coordiantes of
                                                 ! interpolated points
  REAL,INTENT(IN OUT),DIMENSION(:):: eta         ! array of y coordiantes of
                                                 ! interpolated points
  INTEGER,INTENT(IN OUT),DIMENSION(:):: lambda   ! array of edge numbers for
                                                 ! j-th interpolated point
  INTEGER,INTENT(IN OUT):: j                     ! number of points in the
                                                 ! list of interpolated points
  INTEGER,INTENT(IN),DIMENSION(:):: ibe          ! the list of boundary edges
                                                 ! from the triangulation
  INTEGER,INTENT(IN),DIMENSION(:,:):: ite        ! linked list of indices of
                                                 ! adjacent edges provided
                                                 ! by the triangulation
                                                 ! procedure.

  INTEGER:: i
  INTEGER:: j1,l1
  INTEGER:: j1big
  INTEGER:: jcnt, lcnt
  INTEGER:: jj,kk
  INTEGER:: knt
  INTEGER:: npoint
  REAL,DIMENSION(1494):: xx,yy

!!!      DIMENSION  XI(1494),ETA(1494),LAMBDA(1494),IBE(1494),XX(1494),    &
!!!     &           YY(1494),ITE(1494,4)


!----------------------------------------------------------------------------
 
! (A) Initialize local variables
 
      IF (j.eq.0) RETURN
   10 j1=0
 
! (B,C) Search the list of edges for a boundary edge (be(i)=1)

   20 j1=j1+1
      l1=lambda(j1)
      IF (ibe(l1).eq.1) GO TO 30
      IF (j1.lt.j) GO TO 20
      GO TO 120

! (D) Search for a boundary edge and put it at the top of the list.
!     Put this interpolated point at the top of the list for this contour,
!     set j1

   30 IF (j1.eq.j) GO TO 50
      xi(j+1)=xi(j1)
      eta(j+1)=eta(j1)
      lambda(j+1)=lambda(j1)
      DO 40 jcnt=j1,j
        xi(jcnt)=xi(jcnt+1)
        eta(jcnt)=eta(jcnt+1)
   40   lambda(jcnt)=lambda(jcnt+1)
 
! (E) Search the remaining points for an adjacent (common) edge

   50 j1big=j
      lcnt=l1
   60 j1big=j1big-1
      j1=0
   70 j1=j1+1
      l1=lambda(j1)
      DO i=1,4
        IF (l1==ite(lcnt,i)) GO TO 90
      END DO   

! (F) ERROR - THERE IS NO NEXT POINT.
      IF (j1.lt.j1big) GO TO 70
      GO TO 200

! (G) Put this point at the top of the list.  Continue if it is
! not a boundary edge.

   90 xi(j+1)=xi(j1)
      eta(j+1)=eta(j1)
      lambda(j+1)=lambda(j1)
      DO jcnt=j1,j
        xi(jcnt)=xi(jcnt+1)
        eta(jcnt)=eta(jcnt+1)
        lambda(jcnt)=lambda(jcnt+1)
      END DO
      lcnt=l1
      IF (ibe(l1).ne.1) GO TO 60
 
! (H) Draw the open contour line through the points
!  xi(j1),eta(j1) ...... xi(j1+1),eta(j1+1) ...... xi(j),eta(j)
!  Then reset j and continue

      npoint=j-j1big+1
      IF (npoint.le.1) GO TO 110
      WRITE (*,*) "calling cntcrv (open) with",npoint,"points"
      CALL cntcrv (xi(j1big:j), eta(j1big:j), npoint, zcon)
 
  110 j=j1big-1
! 
! (I)  are there any more points left? . .
! 
!      IF (j) 200,200,10
       IF (j > 0) THEN
         GO TO 10
       ELSE
         GO TO 200
       END IF

! (J) Now draw internal lines (closed contours that do not start or stop
! at boundary edges).  The point at j1big=j in the list is chosen to
! start the contour.
! 
  120 j1big=j+1
      lcnt=lambda(j)

! (K,M,P) Find the next point for this contour (on an edge with a common
!  end point).  Put it at the top of the list, and repeat until no more
!  common edges remain for this line.

  130 j1big=j1big-1
      j1=0
      IF (j1big.gt.j) j1=1
  140 j1=j1+1
      l1=lambda(j1)
      DO i=1,4
        IF (l1.eq.ite(lcnt,i)) GO TO 160
      END DO  
      IF (j1.lt.j1big) GO TO 140

! (L) Otherwise, no adjacent edge was found.
! This contour line is complete, go draw it.
      GO TO 180

  160 xi(j+1)=xi(j1)
      eta(j+1)=eta(j1)
      lambda(j+1)=lambda(j1)
      DO jcnt=j1,j
        xi(jcnt)=xi(jcnt+1)
        eta(jcnt)=eta(jcnt+1)
        lambda(jcnt)=lambda(jcnt+1)
      END DO
      lcnt=l1
      IF (j1big.ne.1) GO TO 130

! (O) Draw the closed contour line, the interpolated line through
!  xi(j1),eta(j1) ...... xi(j),eta(j) ...... xi(j1),eta(j1)
 
  180 jj=j1big
      IF (j1big.ne.1) jj=j1big+1
      knt=0
      DO kk=jj,j
        knt=knt+1
        xx(knt)=xi(kk)
        yy(knt)=eta(kk)
      END DO

      xx(knt+1)=xx(1)
      yy(knt+1)=yy(1)
      npoint=knt+1
      WRITE (*,*) "calling cntcrv (closed) with",npoint," points"
      CALL cntcrv (xx(1:npoint), yy(1:npoint), npoint, zcon)
 
! (P) Reset j.  Establish the next contour line for remaining points
!  or quit the procedure if no more points remain.
 
      j=j1big-1
      IF (j.gt.0) j=j+1
!      IF (j) 200,200,120
      IF (j > 0) GO TO 120
  200 RETURN
END Subroutine Cntour   ! ---------------------------------------------------


END Module ContourPlot   ! ==================================================
