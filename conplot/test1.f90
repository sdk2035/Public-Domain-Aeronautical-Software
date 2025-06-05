INCLUDE 'conplot.f90'

!+
PROGRAM TestCase
! ---------------------------------------------------------------------------
USE ContourPlot
IMPLICIT NONE

  INTEGER,PARAMETER:: NPTS=121
  REAL,PARAMETER,DIMENSION(8):: c = &
    (/ 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0 /)
  
  REAL:: eps
  INTEGER:: errCode
  INTEGER,PARAMETER:: GNU=1, DBG=2
  INTEGER:: ierr
!  INTEGER:: iexp=0, jexp=0, ism=0
  INTEGER:: iexp=2, jexp=3, ism=1
  INTEGER:: i,j,k
  INTEGER:: nc
  REAL,DIMENSION(NPTS):: x,y,z
!----------------------------------------------------------------------------
  OPEN(UNIT=DBG, FILE='test.log', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')

  OPEN(UNIT=GNU, FILE='test.gnu', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')

!  DO k=1,NPTS
!    CALL Random_Number(x(k))
!    CALL Random_Number(y(k))
!    z(k)=(x(k)-0.5)**2 + (y(k)-0.5)**2
!  END DO


  k=0
  DO j=0,10
    DO i=0,10
      k=k+1
      x(k)=0.1*REAL(i)
      y(k)=0.1*REAL(j)
    END DO
  END DO

  z=(x-0.5)**2 + (y-0.5)**2
  z(:)=16.0*z(:)
  WRITE(DBG,'(I4,3F12.6)') (k,x(k),y(k),z(k),k=1,NPTS)

  CALL ContourLines(x,y,z, ism,iexp,jexp, c, eps,ierr)

  WRITE(*,*) "File test.gnu added to your directory"
  STOP

END   ! =====================================================================

!+
SUBROUTINE CNTCRV(x,y,n,z)
! ---------------------------------------------------------------------------
! PURPOSE
IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: x,y
  INTEGER,INTENT(IN):: n
  REAL,INTENT(IN):: z

  INTEGER:: GNU=1
  INTEGER:: k
!----------------------------------------------------------------------------
  WRITE(GNU,'(2F12.5)') (x(k),y(k),k=1,n)   ! replace n with SIZE(x)
  WRITE(GNU,'(A)') " "

  RETURN
END Subroutine CntCrv   ! ---------------------------------------------------

