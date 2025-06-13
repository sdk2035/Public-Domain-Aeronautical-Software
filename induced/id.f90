!+
MODULE InducedDrag
! PURPOSE - Compute the induced drag of symmetrical and asymmetrical
!    span loadings.
! ------------------------------------------------------------------------------
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 23Jul96  0.5   RLC   Translation to Fortran 90
! 24Dec96  0.6   RLC   Restructured the modules
! 01Jan05  0.8   RLC   Changed order of arguments in FillArray to agree
!                         with other programs
!                                                                         
! NOTES-Based on note by J.L.Lundry, J.Aircraft, March 1977, p.309

IMPLICIT NONE

  REAL,PARAMETER,PRIVATE:: PI=3.14159265
  REAL,PARAMETER,PRIVATE:: HALFPI = 0.5*PI

  PUBLIC:: AsymmetricLoadingInducedDrag
  PUBLIC:: SymmetricLoadingInducedDrag
  PUBLIC:: DragFromCoefficients
  PUBLIC:: ComputeFourierCoefficients
  PUBLIC:: DefineLoadingFromFourierCoeff
  PUBLIC:: FillArray
  PRIVATE:: LUdecompose
  PRIVATE:: Solve

!-------------------------------------------------------------------------------
CONTAINS

!+
FUNCTION AsymmetricLoadingInducedDrag(y,ccl,span,nTerms) RESULT(drag)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN),DIMENSION(:):: y,ccl
  REAL,INTENT(IN):: span
  INTEGER,INTENT(IN):: nTerms

  REAL:: drag   ! drag/q,  units of area
  REAL,ALLOCATABLE,DIMENSION(:):: coeff
!----------------------------------------------------------------------------
  ALLOCATE(coeff(nTerms))
  CALL ComputeFourierCoefficients(y,ccl,span,.FALSE.,coeff)
  drag=DragFromCoefficients(.FALSE.,coeff,span)
  DEALLOCATE(coeff)
  RETURN
END Function AsymmetricLoadingInducedDrag   ! -------------------------------


!+
FUNCTION SymmetricLoadingInducedDrag(y,ccl,span,nTerms) RESULT(drag)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN),DIMENSION(:):: y,ccl
  REAL,INTENT(IN):: span
  INTEGER,INTENT(IN):: nTerms

  REAL:: drag   ! drag/q,  units of area
  REAL,ALLOCATABLE,DIMENSION(:):: coeff
!----------------------------------------------------------------------------
  ALLOCATE(coeff(nTerms))
  CALL ComputeFourierCoefficients(y,ccl,span,.TRUE.,coeff)
  drag=DragFromCoefficients(.TRUE.,coeff,span)
  DEALLOCATE(coeff)
  RETURN
END Function SymmetricLoadingInducedDrag   ! --------------------------------

!+
FUNCTION DragFromCoefficients(symLoading,coeff,span) RESULT(drag)
! ---------------------------------------------------------------------------
! NOTE - This is a PURE function

  LOGICAL,INTENT(IN):: symLoading
  REAL,INTENT(IN),DIMENSION(:):: coeff
  REAL,INTENT(IN):: span

  REAL:: drag           ! this will be draq/q, units of area
  REAL,ALLOCATABLE,DIMENSION(:):: b

  INTEGER::  i ,nTerms
!----------------------------------------------------------------------------
  nTerms=SIZE(coeff)
  ALLOCATE(b(nTerms))
  IF (symLoading) THEN
    DO i=1,nTerms
      b(i)=REAL(i)
    END DO
  ELSE
    DO i=1,nTerms
      b(i)=REAL(i+i-1)
    END DO
  END IF
  drag=PI*span*span*DOT_PRODUCT(b,coeff**2)

  DEALLOCATE(b)
  RETURN
END Function DragFromCoefficients   ! ---------------------------------------

!+
SUBROUTINE ComputeFourierCoefficients(y,ccl,span,symLoading,coeff)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN),DIMENSION(:):: y,ccl
  REAL,INTENT(IN):: span   ! FULL span, not semispan!!!
  LOGICAL,INTENT(IN):: symLoading
  REAL,INTENT(OUT),DIMENSION(:):: coeff

  INTEGER:: i,j
  INTEGER:: nLoad, nTerms, nTotal
  REAL:: cond 
  REAL,ALLOCATABLE,DIMENSION(:,:):: bMatrix
  REAL,ALLOCATABLE,DIMENSION(:):: rhs,theta
  INTEGER,ALLOCATABLE,DIMENSION(:):: ipvt

!----------------------------------------------------------------------------
  nLoad=SIZE(y)
  nTerms=SIZE(coeff)
  nTotal=nTerms+nLoad
  ALLOCATE( bMatrix(nTotal,nTotal), ipvt(nTotal) )
  ALLOCATE(theta(nLoad), rhs(nTotal) )
  rhs=0
  theta=ACOS(y*(2.0/span))
  rhs(1:nLoad)=ccl/(4.0*span)
 
  bMatrix = 0.0   ! whole matrix

  IF(symLoading)THEN
    DO j=1,nTerms
      DO i=1,nLoad
        bMatrix(i,j+nLoad)=Sin((j+j-1)*theta(i))
      END DO
    END DO
  ELSE
    DO j=1,nTerms
      DO i=1,nLoad 
        bMatrix(i,j+nLoad)=Sin(j*theta(i))
      END DO
    END DO
  END IF

  bMatrix(nLoad+1:nTotal,1:nload)=Transpose(bmatrix(1:nLoad,nLoad+1:nTotal))

  IF (symLoading) THEN
    DO i=1,nTerms 
      bMatrix(i+nLoad,i+nLoad)=(PI/2)*(i+i-1)**4   ! diagonal
    END DO
  ELSE
    DO i=1,nTerms 
      bMatrix(i+nLoad,i+nLoad)=PI*i**4
    END DO
  END IF

  CALL LUdecompose(bMatrix, ipvt, cond)
  CALL Solve(bMatrix, ipvt, rhs)
  coeff=rhs(nLoad+1:nTotal)

!!!  WRITE(*,*) "Coeff of Fourier sine series"
!!!  WRITE(*,'(I3,F12.7)' ) (i, coeff(i),i=1,nTerms)
  WRITE(*,*) "Condition number=", cond
  DEALLOCATE(bMatrix, ipvt, rhs, theta)
  RETURN
END Subroutine ComputeFourierCoefficients   ! -------------------------------

!+
SUBROUTINE DefineLoadingFromFourierCoeff(coeff,span,symLoading,y,ccl)
! ---------------------------------------------------------------------------
! PURPOSE - Define the loading (y,ccl) defined by a set of Fourier sine
!    series coefficients, coeff. Usually this will be a loading at a much
!    denser set of points than the original loading data used to compute
!    the coefficients.
!
  REAL,INTENT(IN),DIMENSION(:):: coeff
  REAL,INTENT(IN):: span
  LOGICAL,INTENT(IN):: symLoading
  REAL,INTENT(OUT),DIMENSION(:):: y,ccl

  REAL,PARAMETER:: PI = 3.14159265

  INTEGER:: i,j
  INTEGER:: nrich   ! number of points to be computed
  INTEGER:: nTerms  ! number of coefficients
  REAL,ALLOCATABLE,DIMENSION(:):: theta
  REAL,ALLOCATABLE,DIMENSION(:):: sinTerms
!----------------------------------------------------------------------------
  nrich=SIZE(y)
  nTerms=SIZE(coeff)
  ALLOCATE(theta(nrich))
  ALLOCATE(sinTerms(nterms) )

  IF(symLoading) THEN
    CALL FillArray(0.0, PI/2, theta, 1)    
  ELSE
    CALL FillArray(0.0, PI, theta, 1)
  END IF
  y=0.5*span*COS(theta)

  DO j=1,nrich
    IF(symLoading) THEN
      DO i=1,nTerms
        sinTerms(i)=SIN((i+i-1)*theta(j))
      END DO
    ELSE
      DO i=1,nterms
        sinTerms(i)=SIN((i)*theta(j))
      END DO
    END IF
    ccl(j)=DOT_PRODUCT(coeff, sinTerms)
  END DO

  ccl=4.0*span*ccl
  DEALLOCATE(theta)
  DEALLOCATE(sinTerms)
  RETURN
END Subroutine DefineLoadingFromFourierCoeff   ! ----------------------------

!+
SUBROUTINE FillArray(startValue,endValue, array, spacingCode)
! ---------------------------------------------------------------------------
! PURPOSE - fill an array from start to end. The intermediate points are
!    computed according to various spacing rules.

  REAL,INTENT(IN):: startValue,endValue
  REAL,INTENT(OUT),DIMENSION(:):: array
  INTEGER,INTENT(IN):: spacingCode  ! =2 full cosine
                                    ! =3 half cosine
                                    ! =4 half sine
                                    ! anything else = uniform spacing

  INTEGER:: i,n
  REAL,ALLOCATABLE,DIMENSION(:):: temp
!----------------------------------------------------------------------------
  n=SIZE(array)
  IF (n <= 0) RETURN

  array(n)=endValue
  array(1)=startValue
  IF (n <= 2) RETURN

  ALLOCATE(temp(n-2))
  DO i=1,n-2
    temp(i)=REAL(i)
  END DO
  temp=temp*(1.0/REAL(n-1))

  SELECT CASE(spacingCode)
    CASE (2)
      temp=0.5*(1.0-COS(PI*temp))         ! full cosine, dense near both ends
    CASE (3)
      temp=1.0-COS(HALFPI*temp)               ! half cosine, dense near start
    CASE (4)
      temp=SIN(HALFPI*temp)                       ! half sine, dense near end
  END SELECT

  array(2:n-1)=startValue + (endValue-startValue)*temp

  DEALLOCATE(temp)

  RETURN
END Subroutine FillArray   ! ------------------------------------------------

!+
SUBROUTINE LuDecompose(a,ip,conditionNumber)
!   -------------------------------------------------------------------------
! PURPOSE LU decomposition of matrix for Gaussian elimination
!     ACM Algorithm 423, April 1972, by Cleve Moler
!     Thanks to the Office of Naval Research, Contract NR 044-377.
IMPLICIT NONE
!----------------------------------------------------------------------------
  REAL,INTENT(IN OUT),DIMENSION(:,:):: a
  INTEGER,INTENT(OUT),DIMENSION(:)::ip
  REAL,INTENT(OUT),optional:: conditionNumber

  REAL:: anorm,ynorm,znorm,ek
  INTEGER:: n
  INTEGER:: i,j,k,m
  REAL:: t
  REAL,ALLOCATABLE,DIMENSION(:):: work

  n=SIZE(a(:,1))
!  CALL Assert(SIZE(ip) >= n, 'ip too small')
!  CALL Assert(n==SIZE(a(1,:)), 'matrix not square')

  anorm=0.0    ! compute the a-norm of a
  DO j=1,n
    t=SUM(ABS(a(:,j)))
    anorm=MAX(anorm,t)
  END DO

  IP(N)=1
  DO k=1,n-1
    m=k
    DO i=k+1,n
      IF (ABS(a(i,k)) .GT. ABS(a(m,k)) ) m=i
    END DO

    ip(k)=m
    IF (m.NE.k) ip(n)=-ip(n)
    t=a(m,k)
    a(m,k)=a(k,k)
    a(k,k)=t
    IF (t .EQ. 0.) CYCLE
    a(k+1:n,k) = (-1.0/t)*a(k+1:n,k)
    DO j=k+1,n
      t=a(m,j)
      a(m,j)=a(k,j)
      a(k,j)=t
     IF (t /= 0.0) a(k+1:n,j)=a(k+1:n,j) + t*a(k+1:n,k)
    END DO
  END DO

!....Solve for (a-transpose)*y=e
  ALLOCATE(work(n))
  DO k=1,n
    IF (a(k,k)==0.0) THEN
      conditionNumber = 1E32 ! something HUGE
      RETURN
    END IF
    t=DOT_PRODUCT(a(1:k-1,k), work(1:k-1))
    IF (t < 0.0) THEN
      ek=-1.0
    ELSE
      ek=1.0
    END IF
    work(k)=(ek+t)/a(k,k)
  END DO
  
  DO k=n-1,1,-1
    work(k)=DOT_PRODUCT(a(k+1:n,k), work(k+1:n))
    m=ip(k)
    IF (m /= k) THEN
      t=work(m)
      work(m)=work(k)
      work(k)=t
    END IF
  END DO

  ynorm=SUM(ABS(work))
  CALL Solve(a,ip,work)
  znorm=SUM(ABS(work))
  IF (Present(conditionNumber)) conditionNumber=anorm*znorm/ynorm

  DEALLOCATE(work)
  RETURN
END Subroutine LuDecompose  ! -----------------------------------------------

SUBROUTINE Solve(a,ip,b)
  IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:,:)::a
  INTEGER,INTENT(IN),DIMENSION(:):: ip
  REAL,INTENT(IN OUT),DIMENSION(:):: b

  INTEGER:: n
  INTEGER:: i,k,m
  REAL:: t

  n=SIZE(a(:,1))
      
  DO k=1,n-1
    m=ip(k)
    t=b(m)
    b(m)=b(k)
    b(k)=t
    b(k+1:n)=b(k+1:n) + t*a(k+1:n,k)
!!!    DO 7 I=KP1,N
!!!          B(I)=B(I)+A(I,K)*T
!!!    7   CONTINUE
  END DO

  DO k=n,1,-1
    b(k)=b(k)/a(k,k)
    t=-b(k)
!    b(1:k-1) = b(1:k-1) + t*a(1:k-1,k)
    DO i=1,k-1
      b(i)=b(i)+a(i,k)*t
    END DO
  END DO

  RETURN
END Subroutine Solve   ! ----------------------------------------------------
 

END Module InducedDrag   ! ==================================================
