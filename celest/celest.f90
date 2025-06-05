!+
MODULE CelestialTransformations
! ------------------------------------------------------------------------------
! PURPOSE - Perform transformations among the equatorial, 
!  ecliptic, and galactic coordinate systems

! AUTHORS - Nadine A. Bicket and Gilmer A. Gary, Space Sciences Laboratory
!           Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS RERSON  STATEMENT OF CHANGES
! 16Sep69  0.1 NAB&GAC Release of NASA TM 53943
! 03May08  0.2   RLC   Scanned document for source code
! 25Nov09  0.3   RLC   Conversion to Fortran 90 module
! 03Dec09  0.35  RLC   Added INTENT, other style changes
! 12Dec09  0.4   RLC   Made EclipticToEquatorial
! 13Dec09  0.5   RLC   Debugged and made checks against tables
! 03Jan13  0.6   RLC   Corrected error in GalacticToEquatorial
!                        (thanks to Nathaniel Cook for spotting this)

! NOTE - The code in this module is taken from the NASA reference.
!  A number of changes have been made to bring it up to modern Fortran 
!  standards. Procedures are provided for conversion from equatorial to
!  ecliptic and the inverse and also from equatorial to galactic and
!  the inverse. We could add code for ecliptic to galactic, but it is easy
!  to combine two transformations. The table that corresponds to Table 3 in
!  the reference document was made this way.
!  The code is written for double precision. Change to 
!     SELECTED_REAL_KIND(6) in the code below, if you want single precision.
IMPLICIT NONE
  INTEGER,PARAMETER:: WP = SELECTED_REAL_KIND(14)   ! 14 decimal digits
  CHARACTER(LEN=*),PARAMETER:: &
    CELESTIAL_TRANSFORMATIONS_VERSION = "0.6 (3 January 2013)"

  REAL(WP),PARAMETER:: ZERO = 0.0_WP, HALF = 0.5_WP, ONE=1.0_WP, TWO=2.0_WP
  REAL(WP),PARAMETER:: A90 = 90, A180 = 180, A270 = 270, A360 = 360
  REAL(WP),PARAMETER:: PI = 3.141592653589793238462643383279502884197_WP
  REAL(WP),PARAMETER:: HALFPI=PI/2, TWOPI=PI+PI
  REAL(WP),PARAMETER:: DEG2RAD = PI/180,   RAD2DEG = 180/PI

  REAL(WP),PARAMETER:: EPSILON_DEG = 23.45227222_WP
  REAL(WP),PARAMETER:: EPSILON_RAD = .4093184946_WP
  REAL(WP),PARAMETER:: SINEPS=0.397985013_WP, COSEPS=0.917391916_WP

  REAL(WP),PARAMETER:: GALACTIC_OBLIQUITY_DEG = 62.6_WP
  REAL(WP),PARAMETER:: GALACTIC_OBLIQUITY_RAD = 1.09257611_WP
  REAL(WP),PARAMETER:: SINGAL=0.887815_WP, COSGAL=0.460199784_WP
!-------------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE EquatorialToEcliptic(alpha,delta, lambda,beta)
! ------------------------------------------------------------------------------
! PURPOSE - Convert equatorial coordinates to ecliptic coordinates
  REAL(WP),INTENT(IN):: alpha,delta  ! degrees
  REAL(WP),INTENT(OUT):: lambda,beta  ! degrees

  REAL(WP):: alphaLocal,deltaLocal
!  REAL(WP):: array,drray
  REAL(WP):: sinalp,cosalp
  REAL(WP):: sindel,cosdel
  REAL(WP):: sinbet,cosbet
  REAL(WP):: sinlam,coslam
!-------------------------------------------------------------------------------
  alphaLocal=alpha
  deltaLocal=delta

! The next steps take care of the trouble spots where the
!     spherical triangle becomes an arc

  IF (alphaLocal==A90) THEN
    beta=deltaLocal-EPSILON_DEG
    lambda=alphaLocal
    IF (lambda < ZERO) lambda=lambda+A360
    RETURN
  END IF

  IF (alphaLocal==A270) THEN
    beta=deltaLocal+EPSILON_DEG
    lambda=alphaLocal
    IF (lambda < ZERO) lambda=lambda+A360
    RETURN
  END IF

!      if(alp-90.)5,9,5
!    5 if(alp-270.)10,8,10

!    8 beta=del+e
!      flam=alp
!      go to 11

!    9 beta=del-e
!      flam=alp
!      go to 11

!   10 continue 

  sinalp=SIN(alphaLocal*DEG2RAD)
  cosalp=COS(alphaLocal*DEG2RAD)
  sindel=SIN(deltaLocal*DEG2RAD)
  cosdel=COS(deltaLocal*DEG2RAD)

  sinbet=COSEPS*sindel - SINEPS*cosdel*sinalp
  cosbet=SQRT(ONE-sinbet**2)
  IF (cosbet==ZERO) THEN
    beta=RAD2DEG*ATAN2(sinbet,cosbet)
    lambda=ZERO
    WRITE(*,*) 'cosbet==0 in EquatorialToEcliptic'
    RETURN
  END IF

  sinlam=(SINEPS*sindel+COSEPS*cosdel*sinalp)/cosbet
  coslam=(cosdel*cosalp)/cosbet
  beta=RAD2DEG*ATAN2(sinbet,cosbet)
  lambda=RAD2DEG*ATAN2(sinlam,coslam)
  IF (lambda < ZERO) lambda=lambda+A360

  RETURN
END Subroutine EquatorialToEcliptic   ! ----------------------------------------


!+
SUBROUTINE EclipticToEquatorial(lambda,beta, alpha,delta)
! ------------------------------------------------------------------------------
! PURPOSE - Convert ecliptic coordinates to equatorial coordinates
  REAL(WP),INTENT(IN):: lambda,beta  ! degrees
  REAL(WP),INTENT(OUT):: alpha,delta  ! degrees

  REAL(WP):: lambdaLocal,betaLocal
  REAL(WP):: array,drray
  REAL(WP):: sinalp,cosalp
  REAL(WP):: sindel,cosdel
  REAL(WP):: sinbet,cosbet
  REAL(WP):: sinlam,coslam
!-------------------------------------------------------------------------------
  lambdaLocal=lambda
  betaLocal=beta

  sinlam=SIN(lambdaLocal*DEG2RAD)
  coslam=COS(lambdaLocal*DEG2RAD)
  sinbet=SIN(betaLocal*DEG2RAD)
  cosbet=COS(betaLocal*DEG2RAD)

  sindel=COSEPS*sinbet + SINEPS*cosbet*sinlam
  cosdel=SQRT(ONE-sindel**2)
  IF (cosdel==ZERO) THEN
    delta=RAD2DEG*ATAN2(sindel,cosdel)
    alpha=ZERO
    WRITE(*,*) 'cosdel==0 in EclipticToEquatorial'
    RETURN
  END IF

  sinalp=(-SINEPS*sinbet + COSEPS*cosbet*sinlam)/cosdel
  cosalp=(cosbet*coslam)/cosdel
  delta=RAD2DEG*ATAN2(sindel,cosdel)
  alpha=RAD2DEG*ATAN2(sinalp,cosalp)
  IF (alpha < ZERO) alpha=alpha+A360

  RETURN
END Subroutine EclipticToEquatorial   ! ----------------------------------------

! PROGRAM LISTING FOR TABLES 1, 2, AND 3 (Continued)
!+
SUBROUTINE EquatorialToGalactic(alpha,delta, lambdaGal,betaGal)
! ------------------------------------------------------------------------------
! PURPOSE - Convert Equatorial to Galactic coordinates

  REAL(WP),INTENT(IN):: alpha,delta   ! equatorial coordinates, degrees
  REAL(WP),INTENT(OUT):: lambdaGal,betaGal   ! galactic coordinates, degrees
  REAL(WP):: alphaLocal,deltaLocal
!  REAL(WP):: array,drray
  REAL(WP):: sinbet,cosbet
  REAL(WP):: sinlam,coslam
  REAL(WP):: sinalp,cosalp
  REAL(WP):: sindel,cosdel

  REAL(WP):: a
!-------------------------------------------------------------------------------

  alphaLocal=alpha
  deltaLocal=delta
  a=alpha-282.25


!     the following {six} steps take into account the displacement of the
!     intersection of the great equatorial circles of the two systems
!     from the first point of Aries

  IF (a < ZERO) THEN
    alphaLocal=A360 + a
  ELSE
    alphaLocal=a
  END IF

!     the next {6} steps take care of the trouble spots where the
!     spherical triangle becomes an arc
  IF ( (alphaLocal==A90) .OR. (alphaLocal==A270) ) THEN
!!!      if(alp-90.)21,22,21
!!!   21 if(alp-270.)23,22,23
!!!   22 flam=90.
    lambdaGal=A90
    betaGal=delta-GALACTIC_OBLIQUITY_DEG
  ELSE
    sinalp=SIN(DEG2RAD*alphaLocal)
    cosalp=COS(DEG2RAD*alphaLocal)
    sindel=SIN(DEG2RAD*deltaLocal)
    cosdel=COS(DEG2RAD*deltaLocal)

    sinbet=COSGAL*sindel-SINGAL*cosdel*sinalp 
    cosbet=SQRT(ONE-sinbet**2)
    IF (cosbet==ZERO) THEN
      betaGal=RAD2DEG*ATAN2(sinbet,cosbet)
      lambdaGal=ZERO
      WRITE(*,*) 'cosbet==0 in EquatorialToGalactic'
      RETURN
    END IF

    sinlam=(SINGAL*sindel+COSGAL*cosdel*sinalp)/cosbet 
    coslam=(cosdel*cosalp)/cosbet
    betaGal=RAD2DEG*ATAN2(sinbet,cosbet)
    lambdaGal=RAD2DEG*ATAN2(sinlam,coslam)
  END IF

  IF (lambdaGal < ZERO) lambdaGal=lambdaGal + A360

  RETURN
END Subroutine EquatorialToGalactic   ! ----------------------------------------

!+
SUBROUTINE GalacticToEquatorial(lambdaGal,betaGal, alpha,delta) !,frray,brray)
! ------------------------------------------------------------------------------
! PURPOSE - Transform galactic coordinates to equatorial

  REAL(WP),INTENT(IN OUT):: lambdaGal,betaGal  ! degrees
  REAL(WP),INTENT(OUT):: alpha,delta   ! equatorial coordinates, degrees

  REAL(WP):: lambdaLocal,betaLocal
  REAL(WP):: sinalp,cosalp
  REAL(WP):: sindel,cosdel
  REAL(WP):: sinlam,coslam
  REAL(WP):: sinbet,cosbet
!-------------------------------------------------------------------------------
  lambdaLocal=lambdaGal
  betaLocal=betaGal

  IF (lambdaLocal==A90) THEN
    IF (betaLocal-27.4_WP <= ZERO) THEN
      alpha=lambdaLocal
      delta=betaLocal-GALACTIC_OBLIQUITY_DEG
    ELSE
      alpha=A180+lambdaLocal
      delta=A180-betaLocal+GALACTIC_OBLIQUITY_DEG
    END IF
  ELSE IF (lambdaLocal==A270) THEN
    IF (betaLocal+27.4_WP <= ZERO) THEN
      alpha=A180-lambdaLocal
      delta=-(A180-GALACTIC_OBLIQUITY_DEG+betaLocal)
    ELSE
      alpha=lambdaLocal
      delta=betaLocal+GALACTIC_OBLIQUITY_DEG
    END IF
  ELSE
!!!    lambdaLocal=lambdaLocal*DEG2RAD
!!!    betaLocal=betaLocal*DEG2RAD
    sinlam=SIN(DEG2RAD*lambdaLocal)
    coslam=COS(DEG2RAD*lambdaLocal)
    sinbet=SIN(DEG2RAD*betaLocal)
    cosbet=COS(DEG2RAD*betaLocal)

    sindel=COSGAL*sinbet - SINGAL*cosbet*sinlam
    cosdel=SQRT(ONE-sindel**2)
    IF (cosdel==ZERO) THEN
      delta=RAD2DEG*ATAN2(sindel,cosdel)
      alpha=ZERO
      WRITE(*,*) 'cosdel==0 in GalacticToEquatorial'
    ELSE
      sinalp=(SINGAL*sinbet + COSGAL*cosbet*sinlam)/cosdel 
      cosalp=(cosbet*coslam)/cosdel
      delta=RAD2DEG*ATAN2(sindel,cosdel)
      alpha=RAD2DEG*ATAN2(sinalp,cosalp)
    END IF
  END IF

  IF (alpha >= 77.75_WP) THEN
    alpha=alpha-77.75_WP
  ELSE
    alpha=alpha+282.25_WP
  END IF

  RETURN
END Subroutine GalacticToEquatorial   ! ----------------------------------------

END Module CelestialTransformations   ! ========================================
