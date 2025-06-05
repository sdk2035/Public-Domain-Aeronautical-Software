INCLUDE 'conplot.f90'

!+
PROGRAM TestCase2
! ---------------------------------------------------------------------------
! PURPOSE - Test the contour plotting algorithm on a set of data points
!   read from an external file

! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY                                                        
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 21Dec99  1.0   RLC   Original coding


USE ContourPlot
IMPLICIT NONE

  CHARACTER(LEN=132):: buffer
  REAL:: eps
  INTEGER:: errCode
  CHARACTER(LEN=132):: fileName
  INTEGER,PARAMETER:: IN=1, DBG=2, GNU=3
  INTEGER:: ierr
  INTEGER:: iexp, jexp, ism
  INTEGER:: i,j,k
  INTEGER:: n
  INTEGER:: nLines
  REAL,ALLOCATABLE,DIMENSION(:):: x,y,z, zContour
  CHARACTER:: yorn
!----------------------------------------------------------------------------
  OPEN(UNIT=DBG, FILE='test2.log', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')

  DO
    WRITE(*,*) 'Enter the name of the data file:'
    READ(*,'(A)') fileName
    IF (fileName== " ") STOP
    OPEN(UNIT=IN, FILE=fileName, STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) THEN
      WRITE(DBG,*) "Reading data from "//Trim(fileName)
      EXIT
    ELSE
      WRITE(*,*) "Unable to open this file. Try again."
    END IF
  END DO

  k=0
  DO
    READ(IN,'(A)',IOSTAT=errCode) buffer
    IF (errCode < 0) THEN
      n=k
      EXIT
    END IF
    IF (Len_Trim(buffer)==0) CYCLE
    k=k+1
  END DO
  WRITE(DBG,*) n, " data points read from "//Trim(fileName)

  ALLOCATE (x(n),y(n),z(n),zContour(n))
  REWIND(UNIT=IN)
  DO k=1,n
    READ(IN,*) x(k),y(k),z(k)
  END DO
  WRITE(DBG,*) "DATA from "//Trim(fileName)
  WRITE(DBG,'(I5,3F15.6)') (k,x(k),y(k),z(k),k=1,n)
  CLOSE(UNIT=IN)

  DO
    WRITE(*,*) 'Enter the name of the file of contour levels:'
    READ(*,'(A)') fileName
    IF (fileName== " ") STOP
    OPEN(UNIT=IN, FILE=fileName, STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) THEN
      WRITE(DBG,*) "Reading data from "//Trim(fileName)
      EXIT
    ELSE
      WRITE(*,*) "Unable to open this file. Try again."
    END IF
  END DO

  DO k=1,n
    READ(IN,*,IOSTAT=errCode) zContour(k)
    IF (errCode < 0) THEN
      nLines=k-1
      EXIT
    END IF
  END DO
  WRITE(DBG,*) "Contour levels from "//Trim(fileName)
  WRITE(DBG,'(I5,F15.6)') (k,zContour(k),k=1,nLines)
  CLOSE(UNIT=IN)

  WRITE(*,*) "Do you want the data smoothed before plotting? "
  READ(*,'(A)') yorn
  IF (yorn=='Y' .OR. yorn=='y') THEN
    ism=1
    WRITE(*,*) "What is iexp? "
    READ(*,*) iexp
    WRITE(*,*) "What is jexp?"
    READ(*,*) jexp
    WRITE(DBG,*) "Data will be smoothed with", iexp,jexp
  ELSE
    ism=0
    iexp=0
    jexp=0
    WRITE(DBG,*) "No data smoothing"
  END IF





  OPEN(UNIT=GNU, FILE='test.gnu', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')

  CALL ContourLines(x,y,z, ism,iexp,jexp, zContour(1:nLines), eps,ierr)

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

  INTEGER:: GNU=3
  INTEGER:: k
!----------------------------------------------------------------------------
  WRITE(GNU,'(2F12.5)') (x(k),y(k),k=1,n)   ! replace n with SIZE(x)
  WRITE(GNU,'(A)') " "
  RETURN
END Subroutine CntCrv   ! ---------------------------------------------------

