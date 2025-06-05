!+
!PROGRAM MakeTablesOfCelestialTransformations

! AUTHORS - Nadine A. Bicket and Gilmer A. Gary, Space Sciences Laboratory
!           Ralph L. Carmichael, Public Domain Aeronautical Software
!
! REVISION HISTORY
!   DATE  VERS RERSON  STATEMENT OF CHANGES
! 16Sep69  0.1 NAB&GAC Release of NASA TM 53943
! 03May08  0.2   RLC   Scanned document for source code
! 25Nov09  0.3   RLC   Conversion to Fortran 90 module
! 12Dec09  0.4   RLC   Created this program to match tables in NASA document
! 14Dec09  0.5   RLC   Output to files; added header

! NOTE - This program creates a subset of the tables in the reference
!  NASA document (TM 53943). If you want an extensive set of tables as given
!  in the report, it should be a simple reprogramming task.
! You may notice that the results differ after the seventh or eighth digit.
! The PDAS results used an Intel CPU with double precision.
INCLUDE 'celest.f90'

!+
PROGRAM MakeTablesOfCelestialTransformations
! ------------------------------------------------------------------------------
! PURPOSE - Reproduce a subset of the tables in NASA TM 53943
USE CelestialTransformations
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(/T13,A,T43,A)'
  CHARACTER(LEN=*),PARAMETER:: FMT2 = '(/T9,A,T24,A,T39,A,T54,A)'
  INTEGER:: j,k
  INTEGER,PARAMETER:: OUT = 2
  REAL(WP):: alpha,delta,lambda,beta, lambdaGal,betaGal
!-------------------------------------------------------------------------------

! change these from 0,360,30 to 0,360,5 and
!              from -90,90,30 to -90,90,5 to make the extensive
! tables shown in the report. You will get 36 times as much output.

  OPEN(UNIT=OUT,FILE='table1.txt',STATUS='REPLACE',ACTION='WRITE')
  WRITE(OUT,*) 'EQUATORIAL TO ECLIPTIC  (angles in degrees)'
  WRITE(OUT,FMT1) 'EQUATORIAL', 'ECLIPTIC'
  WRITE(OUT,FMT2) 'alpha','delta','lambda','beta'
  DO j=0,360,30
    alpha=REAL(j)
    DO k=-90,90,30
      delta=REAL(k)
      CALL EquatorialToEcliptic(alpha,delta, lambda,beta)
      WRITE(OUT,'(4F15.7)') alpha,delta, lambda,beta
    END DO
  END DO
  CLOSE(UNIT=OUT)   ! **************************************************

  OPEN(UNIT=OUT,FILE='table2.txt',STATUS='REPLACE',ACTION='WRITE')
  WRITE(OUT,*) 'ECLIPTIC TO EQUATORIAL   (angles in degrees)'
  WRITE(OUT,FMT1) 'ECLIPTIC', 'EQUATORIAL'
  WRITE(OUT,FMT2) 'lambda','beta', 'alpha','delta'

  DO j=0,360,30
    lambda=REAL(j)
    DO k=-90,90,30
      beta=REAL(k)
      CALL EclipticToEquatorial(lambda,beta, alpha,delta)
      WRITE(OUT,'(4F15.7)') lambda,beta, alpha,delta
    END DO
  END DO
  CLOSE(UNIT=OUT)   ! **************************************************

  OPEN(UNIT=OUT,FILE='table3.txt',STATUS='REPLACE',ACTION='WRITE')
  WRITE(OUT,*) 'ECLIPTIC TO GALACTIC   (angles in degrees)'
  WRITE(OUT,FMT1) 'ECLIPTIC', 'GALACTIC'
  WRITE(OUT,FMT2) 'lambda','beta', 'L','B'

  DO j=0,360,30
    lambda=REAL(j)
    DO k=-90,90,30
      beta=REAL(k)
      CALL EclipticToEquatorial(lambda,beta, alpha,delta)
      CALL EquatorialToGalactic(alpha,delta,lambdaGal,betaGal)
      WRITE(OUT,'(4F15.7)') lambda,beta, lambdaGal,betaGal
    END DO
  END DO
  CLOSE(UNIT=OUT)   ! **************************************************

  OPEN(UNIT=OUT,FILE='table4.txt',STATUS='REPLACE',ACTION='WRITE')
  WRITE(OUT,*) 'EQUATORIAL TO GALACTIC   (angles in degrees)'
  WRITE(OUT,FMT1) 'EQUATORIAL', 'GALACTIC'
  WRITE(OUT,FMT2) 'lambda','beta', 'L','B'

  DO j=0,360,30
    alpha=REAL(j)
    DO k=-90,90,30
      delta=REAL(k)
      CALL EquatorialToGalactic(alpha,delta,lambdaGal,betaGal)
      WRITE(OUT,'(4F15.7)') alpha,delta,lambdaGal,betaGal
    END DO
  END DO
  CLOSE(UNIT=OUT)   ! **************************************************

  OPEN(UNIT=OUT,FILE='table5.txt',STATUS='REPLACE',ACTION='WRITE')
  WRITE(OUT,*) 'GALACTIC TO EQUATORIAL   (angles in degrees)'
  WRITE(OUT,FMT1) 'GALACTIC', 'EQUATORIAL'
  WRITE(OUT,FMT2) 'L','B', 'alpha','delta'

  DO j=0,360,30
    lambdaGal=REAL(j)
    DO k=-90,90,30
      betaGal=REAL(k)
      CALL GalacticToEquatorial(lambdaGal,betaGal, alpha,delta)   !!!,dummy1,dummy2)
      WRITE(OUT,'(4F15.7)') lambdaGal,betaGal, alpha,delta
    END DO
  END DO
  CLOSE(UNIT=OUT)   ! **************************************************

  STOP
END PROGRAM MakeTablesOfCelestialTransformations   ! ===========================

 