
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 02Dec01  0.5   RLC   Collected from other directories
! 18Dec01  0.6   RLC   Added FillArray and GetTimeDateStr
! 06Jan02  1.0   RLC   Final cleanup for release of PDAS v7
! 12Mar19  1.1   RLC   Added output statements   



!+
      MODULE Helpers
! ------------------------------------------------------------------------------
! PURPOSE - A collection of useful procedures for this and other programs.
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
! NOTE - The function HypotPDAS should never be needed. It was originally named
!   HYPOT and was used before the intrinsic function was added to modern
!   Fortran in 2008. It might be useful if you have a very old Fortran compiler.


!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 02Dec01  0.5   RLC   Collected from other directories
! 18Dec01  0.6   RLC   Added FillArray and GetTimeDateStr
! 06Jan02  1.0   RLC   Final cleanup for release of PDAS v7

! NOTE - The source is maintained in compatible fixed/free format.

      IMPLICIT NONE

      CHARACTER(LEN=*),PARAMETER:: HELPERS_VERSION =              &
     &  "1.0 (6 January 2002)"


      PUBLIC:: FillArray
      PUBLIC:: GetDateTimeStr
      PUBLIC:: HypotPDAS
      PUBLIC:: PrintArrays
      PUBLIC:: PrintArraysNumbered
      PUBLIC:: PrintMatrix
      PUBLIC:: Swap
      
!-------------------------------------------------------------------------------

      CONTAINS
!+
      SUBROUTINE FillArray(start,end, array, spacingCode)
! ------------------------------------------------------------------------------
! PURPOSE - fill an array from start to end. The intermediate points are
!    computed according to various spacing rules.


      REAL,INTENT(IN):: start,end
      REAL,INTENT(OUT),DIMENSION(:):: array
      INTEGER,INTENT(IN),OPTIONAL:: spacingCode
                                ! =2 full cosine
                                ! =3 half cosine
                                ! =4 half sine
                                ! anything else (or nothing)= uniform spacing
      INTEGER:: k,n
      REAL,PARAMETER:: PI = 3.14159265358979323846264338, HALFPI=PI/2
      REAL,ALLOCATABLE,DIMENSION(:):: temp
!-------------------------------------------------------------------------------
      n=SIZE(array)
      IF (n <= 0) RETURN

      array(n)=end
      array(1)=start
      IF (n <= 2) RETURN

      ALLOCATE(temp(n-2))
      temp= (/ (REAL(k), k=1,n-2) /) / REAL(n-1)

      IF (Present(spacingCode)) THEN
        SELECT CASE(spacingCode)
          CASE (2)
            temp=0.5*(1.0-COS(PI*temp))   ! full cosine, dense near both ends
          CASE (3)
            temp=1.0-COS(HALFPI*temp)         ! half cosine, dense near start
          CASE (4)
            temp=SIN(HALFPI*temp)                 ! half sine, dense near end
        END SELECT
      END IF

      array(2:n-1)=start + (end-start)*temp

      DEALLOCATE(temp)
      RETURN
      END Subroutine FillArray   ! ---------------------------------------------

!+
      FUNCTION GetDateTimeStr() RESULT(s)
! ------------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time
!   This is now standard Fortran. It should work on ALL compilers!
!   You can change the first I2.2 below to I2 if you don't like a
!   leading zero on early morning times, i.e., 6:45 instead of 06:45
      CHARACTER(LEN=*),PARAMETER:: MONTH=                               &
     & 'JanFebMarAprMayJunJulAugSepOctNovDec'
      CHARACTER(LEN=*),PARAMETER:: FMT = '(I2.2,A1,I2.2,I3,A3,I4)'
      CHARACTER(LEN=15):: s
      INTEGER,DIMENSION(8):: v
!-------------------------------------------------------------------------------
      CALL DATE_AND_TIME(VALUES=v)

      WRITE(s,FMT) v(5), ':', v(6), v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
      RETURN
      END FUNCTION GetDateTimeStr   ! ------------------------------------------

!+
      FUNCTION HypotPDAS(x,y) RESULT(z)
! ------------------------------------------------------------------------------
! PURPOSE - Compute SQRT(x**2 + y**2) with care for overflow.
!   Compatible with functions in C and Pascal
!   Identical to the intrinsic function HYPOT added to modern Fortran in 2008.
      REAL,INTENT(IN):: x,y
      REAL:: z
      REAL:: xx,yy
!-------------------------------------------------------------------------------
      xx=ABS(x)
      yy=ABS(y)
      IF (xx > yy) THEN
        z=xx*SQRT(1.0+(yy/xx)**2)
      ELSE
        IF (yy==0.0) THEN
          z=0.0
        ELSE
          z=yy*SQRT(1.0+(xx/yy)**2)
        END IF
      END IF
      RETURN
      END Function HypotPDAS   ! -----------------------------------------------

!+
      SUBROUTINE PrintMatrix(efu,a,s)
! ------------------------------------------------------------------------------
! PURPOSE - Print the contents of a matrix by rows

      INTEGER,INTENT(IN):: efu  ! file number
      REAL,INTENT(IN),DIMENSION(:,:):: a
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: s

      INTEGER:: k
!-------------------------------------------------------------------------------
      IF (Present(s)) THEN
        WRITE(efu,*) "Printing contents of matrix "//Trim(s)
      ELSE
        WRITE(efu,*) "Printing a matrix"
      END IF

      DO k=1,SIZE(a,1)
        WRITE(efu,*) 'Row', k
        WRITE(efu,'(5ES15.6)') a(k,:)
      END DO

      RETURN
      END Subroutine PrintMatrix   ! -------------------------------------------

!+
      SUBROUTINE PrintArrays(efu, x,y, a,b,c,d)
! ------------------------------------------------------------------------------
! PURPOSE - Print arrays, neatly formatted. X and Y are always printed.
!   The next 4 optional variables are printed when present.

      INTEGER,INTENT(IN):: efu   ! unit number of the external file
      REAL,INTENT(IN),DIMENSION(:):: x,y
      REAL,INTENT(IN),DIMENSION(:),OPTIONAL:: a,b,c,d

      CHARACTER(LEN=90):: buffer
      CHARACTER(LEN=*),PARAMETER:: FMT = '(2ES15.7)'
      INTEGER:: k
!-------------------------------------------------------------------------------
      buffer=""
      DO k=1,SIZE(x)
        WRITE(buffer(1:30),FMT) x(k),y(k)
        IF (Present(a)) WRITE(buffer(31:45),FMT) a(k)
        IF (Present(b)) WRITE(buffer(46:60),FMT) b(k)
        IF (Present(c)) WRITE(buffer(61:75),FMT) c(k)
        IF (Present(d)) WRITE(buffer(76:90),FMT) d(k)
        WRITE(efu,*) Trim(buffer)
      END DO
      RETURN
      END Subroutine PrintArrays   ! -------------------------------------------

!+
      SUBROUTINE PrintArraysNumbered(efu, x,y, a,b,c,d)
! ------------------------------------------------------------------------------
! PURPOSE - Print arrays, neatly formatted. X and Y are always printed.
!   The next 4 optional variables are printed when present.
!   The format FMT was chosen to agree with NASA TN D-7397. 
      INTEGER,INTENT(IN):: efu   ! unit number of the external file
      REAL,INTENT(IN),DIMENSION(:):: x,y
      REAL,INTENT(IN),DIMENSION(:),OPTIONAL:: a,b,c,d

      CHARACTER(LEN=95):: buffer
      CHARACTER(LEN=*),PARAMETER:: FMT = '(2ES15.7)'
      INTEGER:: k
!----------------------------------------------------------------------------
      buffer=""
      DO k=1,SIZE(x)
        WRITE(buffer(1:5),'(I5)') k
        WRITE(buffer(6:35),FMT) x(k),y(k)
        IF (Present(a)) WRITE(buffer(36:50),FMT) a(k)
        IF (Present(b)) WRITE(buffer(51:65),FMT) b(k)
        IF (Present(c)) WRITE(buffer(66:80),FMT) c(k)
        IF (Present(d)) WRITE(buffer(81:95),FMT) d(k)
        WRITE(efu,*) Trim(buffer)
      END DO
      RETURN
      END Subroutine PrintArraysNumbered   ! --------------------------------

!+
      SUBROUTINE Swap(x,y)
! ---------------------------------------------------------------------------
! PURPOSE - Interchange the values of x and y
      REAL,INTENT(IN OUT):: x,y
      REAL::xx
!----------------------------------------------------------------------------
      xx=x
      x=y
      y=xx
      RETURN
      END Subroutine Swap   !------------------------------------------------

      END Module Helpers   !=================================================



!+
MODULE SmoothingSpline
! ---------------------------------------------------------------------------
! PURPOSE - Define a cubic spline that is an approximation to a set of
!   data points.
! AUTHORS - Robert E. Smith, Jr., NASA Langley Research Center
!           Joseph M. Price, NASA Langley Research Center
!           Lona M. Howser, NASA Langley Research Center
!           Ralph L. Carmichael, Public Domain Aeronautical Software

IMPLICIT NONE

!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   Feb74  1.0 RES,JMP,LMH  Publication of NASA TN D-7397
! 12Aug01  1.1   RLC   Began conversion to modern Fortran
! 13Aug01  1.11  RLC   Made CountDataPoints
! 14Aug01  1.12  RLC   Started subroutine DefineSpline
! 16Aug01  1.13  RLC   Added PrintMatrix to help debugging
! 02Dec01  1.14  RLC   Renamed as fair
! 18Dec01  1.15  RLC   Incorporated PClookup2 and Decomp and Solve

!----------------------------------------------------------------------------

  INTEGER,PARAMETER,PRIVATE:: IDBG=3  ! only for debugging use
  CHARACTER(LEN=*),PARAMETER:: SMOOTHINGSPLINE_VERSION = "1.15 (18 Dec 01)"

  PUBLIC:: CountDataPoints
  PRIVATE:: Decomp,Solve
  PUBLIC:: DefineSpline
  PRIVATE:: Lookup

  INTERFACE PClookup2
    MODULE PROCEDURE PClookup2scalar,PClookup2vector
  END INTERFACE

CONTAINS

!+
SUBROUTINE CountDataPoints(x,r,n)
! ---------------------------------------------------------------------------
! PURPOSE - Count the number of data points in each curve region.
! NOTE - The coding in this subroutine does not require that x be sorted,
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(IN),DIMENSION(:):: r
  INTEGER,INTENT(OUT),DIMENSION(:):: n

  INTEGER:: j,k
  INTEGER:: nCurves
  INTEGER:: nData
!----------------------------------------------------------------------------
  nData=SIZE(x)
  nCurves=SIZE(r) - 1

  n(:)=0

  DO k=1,nData
    IF (x(k) <= r(2)) n(1)=n(1)+1
  END DO

  DO j=2,nCurves
    DO k=1,nData
      IF (x(k)>r(j) .AND. x(k)<=r(j+1)) n(j)=n(j)+1
    END DO
  END DO

  RETURN
END Subroutine CountDataPoints   ! ------------------------------------------

!+
SUBROUTINE DefineSpline(x,y,w,r,closed,s,spp)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the second derivatives and functional values at the
!   endpoints and at the junction points. Based on the original subroutine
!   named CUSPFIT.
USE Helpers,ONLY: PrintMatrix
IMPLICIT NONE

  REAL,INTENT(IN),DIMENSION(:):: x,y,w
  REAL,INTENT(IN),DIMENSION(:):: r
  LOGICAL,INTENT(IN):: closed
  REAL,INTENT(OUT),DIMENSION(:):: s    ! value of spline at the joints
  REAL,INTENT(OUT),DIMENSION(:):: spp  ! 2nd deriv of spline at the joints

  REAL,ALLOCATABLE,DIMENSION(:,:):: a
  REAL,ALLOCATABLE,DIMENSION(:,:):: c
  REAL:: conditionNumber
  INTEGER,ALLOCATABLE,DIMENSION(:):: dataPoints
!  INTEGER:: errCode
  REAL,ALLOCATABLE,DIMENSION(:):: h
  INTEGER,ALLOCATABLE,DIMENSION(:):: ipivot
  INTEGER:: I,II,IK,IL,IM,IN
  INTEGER:: J
  INTEGER:: k,kk
  INTEGER:: LL,LR
  INTEGER:: LLJ
  INTEGER:: m,mm
  INTEGER:: nCurves
  INTEGER:: nData
  INTEGER:: nJoints
  INTEGER:: nRowsC
!  INTEGER:: NKR1

  REAL:: XI,XIM1

!----------------------------------------------------------------------------
  nJoints=SIZE(r)
  nCurves=nJoints-1
  nData=SIZE(x)

  ALLOCATE(a(nData,4))
  ALLOCATE(h(nJoints))
  ALLOCATE(dataPoints(nCurves+1))

  nRowsC=1 + 3*nCurves   ! number of rows in the C-matrix
  IF (closed) nRowsC=nRowsC+2
!  LA=LO
  LR=nRowsC+1            ! number of columns in the C-matrix
  LL=LR-nCurves
  IF (closed) LL=LL-2

!  NKR1=nCurves+1   ! maybe not needed
  WRITE(*,*) "nRowsC,LL=", nRowsC,LL
  WRITE(IDBG,*) "nRowsC,LL=", nRowsC,LL

  CALL CountDataPoints(x,r,dataPoints)

  ALLOCATE(c(nRowsC,nRowsC+1))   ! extra column for right-hand-side
  c(:,:)=0.0
  WRITE(IDBG,*) "c zeroed"

  h(2:nJoints)=r(2:nJoints)-r(1:nCurves)    ! h(1) is never used; not defined
!  COMPUTES COEFFICIENT MATRIX(A)
  IN=0
  IM=1
  II=0

  bigLoop: DO J=2,nCurves+1   ! start of big loop
    ik=dataPoints(j-1)
    IF (ik==0) GO TO 19
!    L(J)=R(J)-R(J-1)
    in=ik+in  ! updated end count of loop
    DO i=im,in
      k=i-im+1    ! starts at 1
      xi=x(i)-r(j)
      xim1=x(i)-r(j-1)
      a(k,1)= -xi**3/(6.0*h(j)) + h(j)*xim1/6.0 - h(j)**2/6.0
      a(k,2)=1.0-xim1/h(j)
      a(k,3)=xim1**3/(6.0*h(j))-h(j)*xim1/6.0
      a(k,4)=xim1/h(j)
    END DO

!!!!    CALL PrintMatrix(IDBG,a(1:in-im+1,1:4), "A-matrix")  

!      WEIGHT POINTS W * A * X = W * Y
    DO i=im,in
      k=i-im+1
      a(k,1:4)=a(k,1:4)*w(i)
!      A(K,1)=A(K,1)*W(I)
!      A(K,2)=A(K,2)*W(I)
!      A(K,3)=A(K,3)*W(I)
!      A(K,4)=A(K,4)*W(I)
    END DO

    kk=0
    il=1
    DO  k=1,4   !  computes a transpose  a  ; a is in-im+1 by 4
      kk=k+ii
      DO m=il,4
        mm=m+ii
        DO i=1,ik
          c(kk,mm)=a(i,k)*a(i,m)+c(kk,mm)
        END DO
        c(mm,kk)=c(kk,mm)
      END DO
      il=il+1
    END DO

    DO m=1,4   !      computes a transpose * y
      mm=m+ii
      DO i=im,in
        k=i-im+1
        c(mm,lr)=a(k,m)*y(i)+c(mm,lr)
      END DO
    END DO

      IF (j.EQ.2) GO TO 19

!!!      J=J-1
      llj=ll+j-2   ! was llj=ll+j-1
      k=ii-2      ! computes conditions forcing the first derivative to be continuous
      c(llj,1+k)=h(j-1)/6.0   ! was j
      c(llj,2+k)=-1.0/h(j-1)   ! was j
      c(llj,3+k)=(h(j-1)+h(j))/3.0   ! was j and j+1
      c(llj,4+k)=(h(j-1)+h(j))/(h(j)*h(j-1))  ! was j j+1 j+1 j
      c(llj,5+k)=h(j)/6.0   ! was j+1
      c(llj,6+k)=-1.0/h(j)   ! was j+1
!      c(1+k,llj)=c(llj,1+k)
!      c(2+k,llj)=c(llj,2+k)
!      c(3+k,llj)=c(llj,3+k)
!      c(4+k,llj)=c(llj,4+k)
!      c(5+k,llj)=c(llj,5+k)
!      c(6+k,llj)=c(llj,6+k)
      c(k+1:k+6,llj)=c(llj,k+1:k+6)

!      IF (.NOT.closed) GO TO 18
!      IF (J-1 .NE. NKR1-1) GO TO 18   ! was j
    IF (closed .AND. j==nCurves+1) THEN

!  FORCE EQUALITY OF THE FIRST DERIVATIVES AT THE ENDPOINTS
      llj=llj+1
      c(llj,3+k)=-h(j)/6.0   ! was j+1
      c(llj,4+k)=1.0/h(j)   ! was j+1
      c(llj,5+k)=-h(j)/3.   ! was j+1
      c(llj,6+k)=-1.0/h(j)   ! was j+1
      c(llj,1)=-h(2)/3.0
      c(llj,2)=-1.0/h(2)
      c(llj,3)=-h(2)/6.0
      c(llj,4)=1.0/h(2)
!      c(3+k,llj)=c(llj,3+k)
!      c(4+k,llj)=c(llj,4+k)
!      c(5+k,llj)=c(llj,5+k)
!      c(6+k,llj)=c(llj,6+k)
      c(k+3:k+6,llj)=c(llj,k+3:k+6)
!      c(1,llj)=c(llj,1)
!      c(2,llj)=c(llj,2)
!      c(3,llj)=c(llj,3)
!      c(4,llj)=c(llj,4)
      c(1:4,llj)=c(llj,1:4)

!  force equality of the functional values at the endpoints
      llj=llj+1
      c(llj,2)=1.0
      c(2,llj)=1.0
      c(llj,6+k)=-1.0
      c(6+k,llj)=-1.0
    END IF

!!!   18 continue       ! was J=J+1
   19 ii=ii+2
      IM=IN+1

  END DO bigLoop

  WRITE(IDBG,*) "Ready to solve ", nRowsC, " equations."
  CALL PrintMatrix(IDBG,c(:,:),"C-matrix, before")
  ALLOCATE(ipivot(nRowsC))
  CALL Decomp(c(1:nRowsC,1:nRowsC), conditionNumber, ipivot)
  CALL Solve(c(1:nRowsC,1:nRowsC), c(1:nRowsC,nRowsC+1), ipivot)
  CALL PrintMatrix(IDBG, c(:,:), "C-matrix, after")
  WRITE(IDBG,*) "Matrix computed and solved"
  WRITE(IDBG,*) "condition number=",conditionNumber
  WRITE(IDBG,*) nJoints,nCurves,2*nCurves+2,SIZE(s),SIZE(spp),SIZE(c,1),SIZE(c,2)
  s(1:nCurves+1)=c(2:2*nCurves+2:2,LR)
  spp(1:nCurves+1)=c(1:2*nCurves+1:2,LR)

  DEALLOCATE(ipivot,c,a,dataPoints)
  WRITE(*,*) "DefineSpline completed"

  RETURN
END Subroutine DefineSpline   ! ---------------------------------------------

!+
FUNCTION Lookup(xtab,x) RESULT (i)
! ---------------------------------------------------------------------------
! PURPOSE - Search a sorted (increasing) array to find the interval
!  bounding a given number. If n is the size of the array a,
!  return 0 if number x is < a(1)
!  return n if x > a(n).
!  return i if  a(i) <= x < a(i+1).
!  If x is exactly equal to a(n), return n-1.

  REAL,INTENT(IN),DIMENSION(:)::  xtab    ! input array                        
  REAL,INTENT(IN):: x  ! input number              

  INTEGER:: i  ! index of interval such that xtab(i) <= x < xtab(i+1)
               ! =0 if x < xtab(1) and =n if x > xtab(i)
  INTEGER:: j,k,n
!----------------------------------------------------------------------------
  n=SIZE(xtab)
  IF (n <= 0) THEN
    i=-1
    RETURN
  END IF

  IF (x < xtab(1)) THEN
    i=0
    RETURN
  END IF

  IF (x > xtab(n)) THEN
    i=n
    RETURN
  END IF

  i=1 
  j=SIZE(xtab)
  DO
    IF (j <= i+1) EXIT
    k=(i+j)/2                                     ! integer division
    IF (x < xtab(k)) THEN
      j=k
    ELSE
      i=k
    END IF      
  END DO

  RETURN
END Function Lookup   ! -----------------------------------------------------

!+
SUBROUTINE PClookup2scalar(x,y,ypp, u, f)
! ---------------------------------------------------------------------------
! PURPOSE - Interpolate in a cubic spline at one point 
  REAL,INTENT(IN),DIMENSION(:):: x,y,ypp   ! defines the cubic spline
  REAL,INTENT(IN):: u   ! point where spline is to be evaluated
  REAL,INTENT(OUT):: f   ! f(u)

  REAL:: a,fa,fppa, b,fb,fppb
  REAL:: aa,bb,h
!  REAL:: z,zp,zpp,zppp
  INTEGER:: k
!----------------------------------------------------------------------------
  k=LookUp(x,u)
  k=MAX(1, MIN(SIZE(x)-1,k) )
  a=x(k)
  fa=y(k)
  fppa=ypp(k)
  b=x(k+1)
  fb=y(k+1)
  fppb=ypp(k+1)

  h=b-a
  aa=(b-u)/h
  bb=(u-a)/h   ! = 1-aa
  f=aa*fa + bb*fb + ((aa**3-aa)*fppa + (bb**3-bb)*fppb)*h*h/6.0
  RETURN
END Subroutine PClookup2scalar   ! ---------------------------------------------

!+
SUBROUTINE PClookup2vector(x,y,ypp, u, f)
! ---------------------------------------------------------------------------
! PURPOSE - Interpolate in a cubic spline at one point 
  REAL,INTENT(IN),DIMENSION(:):: x,y,ypp   ! defines the cubic spline
  REAL,INTENT(IN),DIMENSION(:):: u   ! points where spline is to be evaluated
  REAL,INTENT(OUT),DIMENSION(:):: f   ! f(u)

  REAL:: a,fa,fppa, b,fb,fppb
  REAL:: aa,bb,h
  INTEGER:: j,k
!----------------------------------------------------------------------------
  DO j=1,MIN(SIZE(u),SIZE(f))
    k=LookUp(x,u(j))
    k=MAX(1, MIN(SIZE(x)-1,k) )
    a=x(k)
    fa=y(k)
    fppa=ypp(k)
    b=x(k+1)
    fb=y(k+1)
    fppb=ypp(k+1)

    h=b-a
    aa=(b-u(j))/h
    bb=1.0-aa
    f(j)=aa*fa + bb*fb + ((aa**3-aa)*fppa + (bb**3-bb)*fppb)*h*h/6.0
  END DO
  RETURN
END Subroutine PClookup2vector   ! -------------------------------------------------

!+
SUBROUTINE Solve(a, b, ipvt)
!   --------------------------------------------------------------------
! PURPOSE - Solve the linear system a*x=b

  IMPLICIT NONE

  REAL,INTENT(IN),DIMENSION(:,:):: a
  REAL,INTENT(IN OUT),DIMENSION(:):: b
  INTEGER,INTENT(IN),DIMENSION(:):: ipvt

  INTEGER:: i,k,m
  INTEGER:: n
  REAL:: t
!----------------------------------------------------------------------------
  n=SIZE(a,1)
  DO k=1,n-1   ! forward elimination
    m=ipvt(k)
    t=b(m)
    b(m)=b(k)
    b(k)=t
    DO i=k+1,n
      b(i)=b(i)+a(i,k)*t
    END DO
  END DO       

  DO k=n,1,-1      ! back substitution
    b(k)=b(k)/a(k,k)
    t=-b(k)
    DO i=1,k-1
      b(i)=b(i)+a(i,k)*t
    END DO       
  END DO
  RETURN
END Subroutine Solve   ! ----------------------------------------------------

!+
SUBROUTINE Decomp(a, cond, ipvt)
! --------------------------------------------------------------------
! FUNCTION- Matrix triangularization by Gaussian elimination and
!    estimation of the condition of the matrix.

  IMPLICIT NONE

  REAL,INTENT(IN OUT),DIMENSION(:,:):: a   ! matrix to be decomposed
  REAL,INTENT(OUT):: cond  ! condition number
  INTEGER,INTENT(OUT),DIMENSION(:):: ipvt   !    index of pivot rows

   INTEGER:: j,k,m
   INTEGER:: n
   REAL:: anorm
   REAL:: t
   REAL:: ek
   REAL:: ynorm,znorm

   REAL,    DIMENSION(SIZE(a,1)) :: work  ! scratch array
   INTEGER, DIMENSION(1) :: mloc  ! receives the result of MAXLOC
!-----------------------------------------------------------------------
   n=SIZE(a,1)
   WRITE(*,*) 'The size of the matrix is ', n
   IF (n <= 1) THEN
      IF (a(1,1) == 0.0) THEN
         cond=HUGE(1.0)  ! 1.0E32      ! really wanted this to be HUGE
      ELSE
         cond=1.0
      END IF
      RETURN        ! note abnormal RETURN
   END IF

  anorm=0.0   ! compute 1-norm of a
  DO j=1,n
    anorm=MAX(anorm, SUM(ABS(a(1:n,j))))
  END DO   

   WRITE(*,*) '1-norm is ', anorm

   DO k=1,n-1   ! Gaussian elimination with partial pivoting
      mloc= MAXLOC(ABS(a(k:n,k)))   ! pivot row
      m=k-1+mloc(1)
!!!   DO i=k+1,n
!!!      IF (ABS(a(i,k)) > ABS(a(m,k)) ) m=i
!!!   END DO
      ipvt(k)=m
      if(m /= k) ipvt(n)=-ipvt(n)
      t=a(m,k)
      a(m,k)=a(k,k)
      a(k,k)=t
      IF (t /= 0.0) THEN
         t=1.0/t
         a(k+1:n,k)=-t*a(k+1:n,k)
         DO j=k+1,n   ! interchange and eliminate by columns.......
            t=a(m,j)
            a(m,j)=a(k,j)
            a(k,j)=t
            IF (t /= 0.0) a(k+1:n,j)=a(k+1:n,j) + t*a(k+1:n,k)
         END DO
      END IF
   END DO

   DO k=1,n   ! solve (a-transpose)*y=e
      t=0.0
      IF (k > 1) t=DOT_PRODUCT(a(1:k-1,k), work(1:k-1))
      ek=1.0
      if (t .lt. 0.) ek=-1.0
      if (a(k,k) .eq. 0.) THEN    ! singular matrix
         cond=1.0E32   ! wanted HUGE
         RETURN                            ! note abnormal RETURN
      END IF
      work(k)=-(ek+t)/a(k,k)
   END DO

   DO k=n-1,1,-1 
      t=0.0
      t=work(k)*SUM(a(k+1:n,k))
      work(k)=t
      m=ipvt(k)
      if (m /= k) THEN
         t=work(m)
         work(m)=work(k)
         work(k)=t
      END IF
   END DO

   ynorm=SUM(ABS(work(1:n)))
   CALL Solve(a,work,ipvt)   ! solve a*z=y
   znorm=SUM(ABS(work(1:n)))
   cond=anorm*znorm/ynorm   ! estimate condition
   if (cond < 1.0) cond=1.0
  RETURN
END Subroutine Decomp   ! ---------------------------------------------------



END Module SmoothingSpline   ! ==============================================

!+
PROGRAM FairData
! ---------------------------------------------------------------------------
! PURPOSE - Define a cubic spline that is an approximation to a set of
!   data points.
! AUTHORS - Robert E. Smith, Jr., NASA Langley Research Center
!           Joseph M. Price, NASA Langley Research Center
!           Lona M. Howser, NASA Langley Research Center
!           Ralph L. Carmichael, Public Domain Aeronautical Software

USE SmoothingSpline
IMPLICIT NONE

!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   Feb74  0.1 RES,JMP,LMH  Publication of NASA TN D-7397
! 12Aug01  0.2   RLC   Began conversion to modern Fortran
! 13Aug01  0.3   RLC   Added extended graphics metafile coding
! 16Aug01  0.4   RLC   Added debugging code
! 02Dec01  0.5   RLC   Renamed program and files
! 15Dec01  0.6   RLC   Added gnuplot output; add'l debugging
! 19Dec01  0.7   RLC   Made Elf90 compatible
! 06Jan02  1.0   RLC   Final cleanup for release of PDAS v7

  INTEGER,PARAMETER:: IN=1, OUT=2, DBG=3, GNU=4
  CHARACTER(LEN=15):: dateTimeStr
  CHARACTER(LEN=*),PARAMETER:: VERSION = "1.0 (6 January 2002)"
!----------------------------------------------------------------------------
  CALL Welcome()
  CALL ReadAndProcessData()
  WRITE(*,*) "Normal termination of program fair"
  STOP

CONTAINS

!+
SUBROUTINE CountInputArrays(preset, xData,yData, joints, nData,nJoints)
! ---------------------------------------------------------------------------
! PURPOSE - Count the number of elements of an array that have been input
!   thru a NAMELIST call. This is to save the user from the error-prone
!   process of counting the number of entries in a list. The array is
!   preset to a very negative value that no one would ever use for real data.
!   This routine simply searches for the first entry in the array that still
!   still has this value.
  REAL,INTENT(IN):: preset   ! the value that was originally placed in each
  REAL,INTENT(IN),DIMENSION(:):: xData,yData,joints   ! arrays to ne counted
  INTEGER,INTENT(OUT):: nData,nJoints

  INTEGER:: k,nx,ny
!----------------------------------------------------------------------------
  DO k=1,SIZE(xData)
    IF (xData(k) <= preset) EXIT
  END DO
  nx=k-1   ! failed on k; k-1 was the last good one

  DO k=1,SIZE(yData)
    IF (yData(k) <= preset) EXIT
  END DO
  ny=k-1

  IF (nx==ny) THEN
    nData=nx
  ELSE
    WRITE(*,*) "WARNING: xData and yData have different lengths!"
    WRITE(*,*) "x-length=",nx, "   y-length=",ny
    nData=MIN(nx,ny)
    WRITE(*,*) "Using ", nData
  END IF

  DO k=1,SIZE(joints)
    IF (joints(k) <= preset) EXIT
  END DO
  nJoints=k-1

  RETURN
END Subroutine CountInputArrays   ! -----------------------------------------

!+
SUBROUTINE ReadAndProcessData()
! ---------------------------------------------------------------------------
! PURPOSE -
!USE Helpers,ONLY: PrintArraysNumbered

  LOGICAL:: closed = .FALSE.
  INTEGER:: errCode
  REAL,DIMENSION(500):: joints
  INTEGER:: nData
  INTEGER:: nJoints
  LOGICAL:: parametric = .FALSE.
  REAL,PARAMETER:: PRESET=-1E20   ! smaller that anybody's data
  REAL,DIMENSION(1000):: xData,yData,weights

  NAMELIST /NAM1/ xData,yData,weights, joints, parametric, closed

!----------------------------------------------------------------------------
  weights=1.0
  xData(:)=PRESET
  yData(:)=PRESET
  joints(:)=PRESET

  READ(UNIT=IN, NML=NAM1, IOSTAT=errCode)
  IF (errCode==0) THEN
    WRITE(DBG,*) "Data read from namelist without error."
  ELSE
    WRITE(*,*) "Error reading data in namelist NAM1"
    WRITE(DBG,*) "Error reading data in namelist NAM1"
  END IF

  CALL CountInputArrays(PRESET, xData,yData,joints, nData,nJoints)

! WRITE(DBG,*) "DATA POINTS AS ENTERED"
! CALL PrintArraysNumbered(DBG,xData(1:nData),yData(1:nData),weights(1:nData))
  
  WRITE(DBG,*) "In ReadAndProcessData, nData=",nData
  IF (parametric) THEN
    CALL ProcessParametric(xData(1:nData),yData(1:nData),               &
      weights(1:nData), joints(1:nJoints),closed)
  ELSE
    CALL ProcessNonParametric(xData(1:nData),yData(1:nData),            &
      weights(1:nData), joints(1:nJoints),closed)    
  END IF

  RETURN
END Subroutine ReadAndProcessData   ! ---------------------------------------

!+
SUBROUTINE FillArrayOfJoints(interiorJoints,start,end,joints,nJoints)
! ---------------------------------------------------------------------------
! PURPOSE - With the array interiorJoints as input, make a new array called
!   joints that includes the interiorJoints with start placed at the beginning
!   and end at the end. BUT, if the user happened to put these values in the
!   array of interiorJoints, don't add them.
!   Also, sort the array joints, so that it is monotone increasing and has
!   no duplicate values. This seems like a lot of work when we could simply
!   put the burden on the user to input "correct" data, but the goal is to
!   make a program that is robust and forgiving and easy to use.
USE Helpers,ONLY: Swap
  REAL,INTENT(IN),DIMENSION(:):: interiorJoints
  REAL,INTENT(IN):: start,end
  REAL,INTENT(OUT),DIMENSION(:):: joints
  INTEGER,INTENT(OUT):: nJoints

  CHARACTER(LEN=*),PARAMETER:: FMT1 = "(A/(5F15.5))"
  INTEGER:: j
  INTEGER,DIMENSION(1):: kk
  INTEGER:: n
  REAL,ALLOCATABLE,DIMENSION(:):: work

!----------------------------------------------------------------------------
  n=SIZE(interiorJoints)
  ALLOCATE(work(n+2))   ! two extra for start and end

!... add start and end to work
  work(1:n)=interiorJoints                       ! never use interiorJoints again
  work(n+1)=start
  work(n+2)=end
  n=n+2
  WRITE(DBG,FMT1) "Unsorted, unculled joints", work

!... first, sort work into increasing order
  DO j=1,n
    kk=MINLOC(work(j:n))              ! location in work(j:n), not in work(1:n)
    CALL Swap(work(j),work(kk(1)+j-1) )
    WRITE(DBG,*) "Sorting, j=",j, kk(1), work
  END DO
  WRITE(DBG,FMT1) "Sorted, unculled joints", work

!... Check for duplicate entries and squeeze them out if they exist
  j=1
  DO
    IF (work(j) >= work(j+1)) THEN
      work(j:n-1)=work(j+1:n)
      n=n-1
    ELSE
      j=j+1
    END IF
    IF (j==n) EXIT
  END DO
  WRITE(DBG,FMT1) "Sorted, culled joints", work(1:n)

!... return the results
  njoints=MIN(n,SIZE(joints))                ! we certainly hope it will be n
  joints(1:nJoints)=work(1:nJoints)
  DEALLOCATE(work)
  RETURN
END Subroutine FillArrayOfJoints   ! ----------------------------------------

!+
SUBROUTINE ProcessNonParametric(xData,yData,weights, interiorJoints,closed)
! ------------------------------------------------------------------------------
! PURPOSE - Perform the spline fit for non-parametric data
USE Helpers,ONLY: PrintArraysNumbered

  REAL,INTENT(IN),DIMENSION(:):: xData,yData,weights
  REAL,INTENT(IN),DIMENSION(:):: interiorJoints
  LOGICAL,INTENT(IN):: closed

  REAL,ALLOCATABLE,DIMENSION(:):: joints
  INTEGER:: nData,nJoints
  REAL,ALLOCATABLE,DIMENSION(:):: s,spp, ycomputed
!-------------------------------------------------------------------------------
  nData=SIZE(xData)
  ALLOCATE(joints(2+SIZE(interiorJoints)) )
  CALL FillArrayOfJoints(interiorJoints,xData(1),xData(nData),joints,nJoints)

  ALLOCATE(s(nJoints),spp(nJoints),ycomputed(nData))

  CALL DefineSpline(xData(1:nData),yData(1:nData),weights(1:nData),     &
    joints(1:nJoints),closed,s,spp)
  write(DBG,*) "spline"
  CALL PrintArraysNumbered(DBG,joints(1:nJoints),s,spp)


  CALL PrintDataNonParametric(xData(1:nData),yData(1:nData), &
     joints(1:nJoints),s,spp)
  

  CALL PlotDataNonParametric(xData(1:nData),yData(1:nData),             &
    joints(1:nJoints),s,spp)

  DEALLOCATE(spp,s)
  RETURN
END Subroutine ProcessNonParametric   ! ----------------------------------------

!+
SUBROUTINE ProcessParametric(xData,yData,weights, interiorJoints,closed)
! ------------------------------------------------------------------------------
! PURPOSE - Perform the fit for parametric data
USE Helpers,ONLY: PrintArraysNumbered

  REAL,INTENT(IN),DIMENSION(:):: xData,yData,weights
  REAL,INTENT(IN),DIMENSION(:):: interiorJoints
  LOGICAL,INTENT(IN):: closed

  REAL,ALLOCATABLE,DIMENSION(:):: joints
  INTEGER:: k
  INTEGER:: nData,nJoints
  REAL,ALLOCATABLE,DIMENSION(:):: sx,sxpp,sy,sypp
  REAL,ALLOCATABLE,DIMENSION(:):: t
!----------------------------------------------------------------------------
  nData=SIZE(xData)
  ALLOCATE(t(nData))
  t(1)=0.0
  DO k=2,nData
    t(k)=t(k-1) + Hypot(xData(k)-xData(k-1), yData(k)-yData(k-1))
  END DO
  WRITE(DBG,*) "DATA FOR PARAMETRIC CURVE FITTING"
  WRITE(DBG,*) "  #     t         x          y"
  CALL PrintArraysNumbered(DBG,t,xData,yData)

  ALLOCATE(joints(2+SIZE(interiorJoints)))
  CALL FillArrayOfJoints(interiorJoints,0.0,1.0,joints,nJoints)
  joints(1:nJoints)=t(nData)*joints(1:nJoints)
  WRITE(DBG,*) "JOINTS FOR PARAMETRIC CURVE FITTING"
  WRITE(DBG,'(5F15.6)') joints(1:nJoints)

  ALLOCATE(sx(nJoints),sxpp(nJoints), sy(nJoints),sypp(nJoints))
  CALL DefineSpline(t,xData,weights, joints(1:njoints),closed,sx,sxpp)
  CALL DefineSpline(t,yData,weights, joints(1:njoints),closed,sy,sypp)

  write(DBG,*) "spline"
  CALL PrintArraysNumbered(DBG,joints(1:nJoints),sx,sxpp,sy,sypp)

  CALL PrintDataParametric(xData(1:nData),yData(1:nData), &
    joints(1:nJoints), sx,sxpp, sy,sypp)

  CALL PlotDataParametric(xData(1:nData),yData(1:nData), &
    joints(1:nJoints), sx,sxpp, sy,sypp)

  DEALLOCATE(sypp,sy,sxpp,sx)

  RETURN
END Subroutine ProcessParametric   ! -------------------------------------------

!+
SUBROUTINE PrintDataNonParametric(x,y, joints,s,spp)
! ------------------------------------------------------------------------------
! PURPOSE - Add output to file fair.out for non-parametric
USE Helpers,ONLY: PrintArraysNumbered
USE SmoothingSpline,ONLY: PClookup2 
  REAL,INTENT(IN),DIMENSION(:):: x,y   ! the data points
  REAL,INTENT(IN),DIMENSION(:):: joints,s,spp
  REAL,ALLOCATABLE,DIMENSION(:):: ycomputed 
  REAL:: sigma
!-------------------------------------------------------------------------------
  ALLOCATE(ycomputed(SIZE(x)))

  CALL PClookup2vector(joints,s,spp,x,ycomputed)
  WRITE(OUT,'(10X,A)') "****************","DATA FOR X VS. Y","****************"
  WRITE(OUT,'(A,I0,5X,A,I0)') "NO. OF CURVES = ", SIZE(joints)-1,               &
     "NO. OF DATA POINTS = ", SIZE(x)
  WRITE(OUT,*) "ENDPOINTS AND JUNCTION POINTS"
  WRITE(OUT,'(5ES15.7)') joints
  WRITE(OUT,'(//14X,A,14X,A,9X,A,9X,A)') "X","Y","Y-COMPUTED","RESIDUAL"
  CALL PrintArraysNumbered(OUT,x,y,ycomputed,ycomputed-y)   
  sigma=SQRT(SUM((ycomputed-y)**2))/SIZE(ycomputed)  ! UN-weighted std dev
  WRITE(OUT,'(A,ES15.7)') " UNWEIGHTED STD DEVIATION = ", sigma
  
  DEALLOCATE(ycomputed)
  
  RETURN 
END Subroutine PrintDataNonParametric   ! --------------------------------------  


!+
SUBROUTINE PrintDataParametric(x,y, joints,sx,sxpp, sy,sypp)
! ------------------------------------------------------------------------------
! PURPOSE - Add output to file fair.out for non-parametric
USE Helpers,ONLY: PrintArraysNumbered
USE SmoothingSpline,ONLY: PClookup2 
  REAL,INTENT(IN),DIMENSION(:):: x,y   ! the data points
  REAL,INTENT(IN),DIMENSION(:):: joints,sx,sxpp, sy,sypp
  REAL,ALLOCATABLE,DIMENSION(:):: t
!-------------------------------------------------------------------------------
  ALLOCATE(t(SIZE(x)))

!  CALL PClookup2vector(joints,s,spp,x,ycomputed)
  WRITE(OUT,'(10X,A)') "***************"," DATA FOR X VS. T","***************"
  WRITE(OUT,'(//A,I0,5X,A,I0)') "NO. OF CURVES = ", SIZE(joints)-1,               &
     "NO. OF DATA POINTS = ", SIZE(x)
  WRITE(OUT,*) "ENDPOINTS AND JUNCTION POINTS"
  WRITE(OUT,'(//14X,A,16X,A,12X,A,13X,A)') "T","X","X-COMPUTED","RESIDUAL"
!  CALL PrintArraysNumbered(OUT,x,y,ycomputed,ycomputed-y)   

  WRITE(OUT,'(10X,A)') "***************"," DATA FOR Y VS. T","***************"
  WRITE(OUT,'(//A,I0,5X,A,I0)') "NO. OF CURVES = ", SIZE(joints)-1,               &
     "NO. OF DATA POINTS = ", SIZE(x)
  WRITE(OUT,*) "ENDPOINTS AND JUNCTION POINTS"
  WRITE(OUT,'(//14X,A,16X,A,12X,A,12X,A)') "T","Y","Y-COMPUTED","RESIDUAL"
  
  DEALLOCATE(t)
  
  RETURN 
END Subroutine PrintDataParametric   ! --------------------------------------  



!+
SUBROUTINE PlotDataNonParametric(x,y, joints,s,spp)
! ---------------------------------------------------------------------------
! PURPOSE - Make files fair.dat, fair.crv, and fair.jnt for non-parametric
USE Helpers,ONLY: FillArray,PrintArrays
USE SmoothingSpline,ONLY: PClookup2
  REAL,INTENT(IN),DIMENSION(:):: x,y   ! the data points
  REAL,INTENT(IN),DIMENSION(:):: joints,s,spp

  INTEGER:: k,k1,k2
  INTEGER,PARAMETER:: MPLOT = 9
  INTEGER:: nJoints,nPlot
!  REAL,DIMENSION(SIZE(x)):: yFaired
  REAL,ALLOCATABLE,DIMENSION(:):: xPlot,yPlot
!----------------------------------------------------------------------------
  OPEN(UNIT=GNU,FILE='fair.dat',STATUS='REPLACE',ACTION='WRITE')
  CALL PrintArrays(GNU,x,y)
  CLOSE(UNIT=GNU)

  nJoints=SIZE(joints)
  OPEN(UNIT=GNU,FILE='fair.crv',STATUS='REPLACE',ACTION='WRITE')
  nPlot=nJoints + MPLOT*(nJoints-1)
  WRITE(DBG,*) 'nplot=', nplot
  ALLOCATE(xPlot(nPlot),yPlot(nPlot))
  k2=1
  DO k=1,nJoints-1
    k1=k2
    k2=k1+MPLOT+1
    CALL FillArray(joints(k),joints(k+1),xPlot(k1:k2))
  END DO
  CALL PClookup2(joints,s,spp, xPlot,yPlot)
  CALL PrintArrays(GNU,xPlot,yPlot)
  CLOSE(UNIT=GNU)

  OPEN(UNIT=GNU,FILE='fair.jnt',STATUS='REPLACE',ACTION='WRITE')
  CALL PrintArrays(GNU,joints,s)
  CLOSE(UNIT=GNU)

  RETURN
END Subroutine PlotDataNonParametric   ! ------------------------------------

!+
SUBROUTINE PlotDataParametric(x,y, joints,sx,sxpp, sy,sypp)
! ---------------------------------------------------------------------------
! PURPOSE - Make files fair.dat, fair.crv, and fair.jnt for non-parametric
USE Helpers,ONLY: FillArray,PrintArrays,PrintArraysNumbered
USE SmoothingSpline,ONLY: PClookup2
  REAL,INTENT(IN),DIMENSION(:):: x,y   ! the data points
  REAL,INTENT(IN),DIMENSION(:):: joints,sx,sxpp, sy,sypp

  INTEGER:: k,k1,k2
  INTEGER,PARAMETER:: MPLOT = 9
  INTEGER:: nJoints,nPlot
  REAL,ALLOCATABLE,DIMENSION(:):: tPlot,xPlot,yPlot
!----------------------------------------------------------------------------
  OPEN(UNIT=GNU,FILE='fair.dat',STATUS='REPLACE',ACTION='WRITE')
  CALL PrintArrays(GNU,x,y)
  CLOSE(UNIT=GNU)

  nJoints=SIZE(joints)
  OPEN(UNIT=GNU,FILE='fair.crv',STATUS='REPLACE',ACTION='WRITE')
  nPlot=nJoints + MPLOT*(nJoints-1)
  WRITE(DBG,*) 'nplot=', nplot
  ALLOCATE(tPlot(nPlot),xPlot(nPlot),yPlot(nPlot))
  k2=1
  DO k=1,nJoints-1
    k1=k2
    k2=k1+MPLOT+1
    CALL FillArray(joints(k),joints(k+1),tPlot(k1:k2))
  END DO

  CALL PClookup2(joints,sx,sxpp, tPlot,xPlot)
  CALL PClookup2(joints,sy,sypp, tPlot,yPlot)
  CALL PrintArrays(GNU,xPlot,yPlot)
  CLOSE(UNIT=GNU)
  CALL PrintArraysNumbered(DBG,tPlot,xPlot,yPlot)

  OPEN(UNIT=GNU,FILE='fair.jnt',STATUS='REPLACE',ACTION='WRITE')
  CALL PrintArrays(GNU,sx,sy)
  CLOSE(UNIT=GNU)

  RETURN
END Subroutine PlotDataParametric   ! ---------------------------------------

!+
SUBROUTINE Welcome()
! ---------------------------------------------------------------------------
! PURPOSE - Get the input file name and open it. Open all other files.
!   Get the current date and time.
USE ISO_FORTRAN_ENV
USE Helpers,ONLY: GetDateTimeStr,HELPERS_VERSION
  INTEGER:: errCode
  CHARACTER(LEN=132):: fileName
  CHARACTER(LEN=*),PARAMETER:: GREETING = "fair - smoothing algorithm"
  INTEGER:: paramCount,paramLength
!----------------------------------------------------------------------------
  WRITE(*,'(3A/A)') 'This program was compiled by ', compiler_version(),   &
                    ' using the options ',           compiler_options()

  WRITE(*,*) GREETING
  WRITE(*,*) "PDAS version " // VERSION
  dateTimeStr=GetDateTimeStr()
  
  paramCount=COMMAND_ARGUMENT_COUNT()
  IF (paramCount > 0) THEN
     CALL GET_COMMAND_ARGUMENT(1, fileName, paramLength)
  END IF
 
  OPEN(UNIT=IN, FILE=fileName, STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
 
  IF (errCode /= 0) THEN 
     DO
        WRITE(*,*) "Enter the name of the input file: "
        READ(*,'(A)') fileName
        IF (Len_Trim(fileName)==0) STOP
        OPEN(UNIT=IN, FILE=fileName, STATUS='OLD', &
           IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
        IF (errCode==0) EXIT
        OPEN(UNIT=IN, FILE=Trim(fileName)//'.nml', STATUS='OLD', &
           IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
        IF (errCode==0) EXIT
        WRITE(*,*) "Unable to open this file. Try again."
     END DO
  END IF 
     
  INQUIRE(UNIT=IN,NAME=fileName)
  WRITE(*,*) "Reading from "//Trim(fileName)

  OPEN(UNIT=DBG,FILE='fair.dbg',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=OUT,FILE='fair.out',STATUS='REPLACE',ACTION='WRITE')
  WRITE(DBG,*) "Program fair, version "//VERSION//"   "//dateTimeStr
  WRITE(DBG,*) "Reading from "//Trim(fileName)
  WRITE(DBG,*) "Use SmoothingSpline, version "//SMOOTHINGSPLINE_VERSION
  WRITE(DBG,*) "Use Helpers, version "//HELPERS_VERSION

  RETURN
END Subroutine Welcome   ! --------------------------------------------------

END Program FairData   ! ====================================================


