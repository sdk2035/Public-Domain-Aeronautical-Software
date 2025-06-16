!+
! PROGRAM MeanAerodynamicChord
! ------------------------------------------------------------------------------
! PURPOSE - Compute the mean aerodynamic chord of a wing.
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 28Feb92  0.5   RLC   Original coding (Fortran 77) (from old fragments)
! 22Oct00  0.6   RLC   Recoded
! 28Dec02  0.7   RLC   Write output to a file
! 16Oct09  0.8   RLC   Write xte in output
! 21Oct09  0.81  RLC   Minor adjustments to output
! 03Nov09  0.85  RLC   All reals to double precision
! 01Mar21  1.0   RLC   Prints date/time in ISO 8601

! NOTE - This program was included in the PDAS collection for several years with
!  the name getmac. Name changed in 2021 to avoid conflict with Microsoft
!  Windows system command with name getmac.

!+
MODULE MeanAerodynamicChordProcedures
! ------------------------------------------------------------------------------
! PURPOSE -Collect the procedures used by the MeanAerodynamicChord program.
  IMPLICIT NONE 

  INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(14)  ! 14 decimal digits precision
  
! A wing is defined by a set of chords. Each chord has the following components:
  TYPE::WingChord
    REAL(DP):: y,xle,c
  END TYPE WingChord
  
  PUBLIC:: GetDateTimeStr
  PUBLIC:: MACofWing
  PRIVATE:: MACofOneSegment
!-------------------------------------------------------------------------------

CONTAINS 


!+
FUNCTION GetDateTimeStr() RESULT (s)
!   ---------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time (ISO 8601) 
  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: MONTH='JanFebMarAprMayJunJulAugSepOctNovDec'
  CHARACTER(LEN=*),PARAMETER:: FMT='(I4,A3,I2.2,1X,I2.2,A1,I2.2)'
  CHARACTER(LEN=15):: s
  INTEGER,DIMENSION(8):: v
!------------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)
  WRITE(s,FMT) v(1), MONTH(3*v(2)-2:3*v(2)), v(3), v(5), ':', v(6)
  RETURN
END FUNCTION GetDateTimeStr   ! ------------------------------------------------

!+
SUBROUTINE MACofOneSegment(c1,c2,area,mac,ymac,xlemac)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the mean aerodynamic chord of one trapezoidal wing 
!  segment defined by two chords.
  TYPE(WingChord),INTENT(IN):: c1,c2  ! the two chords
  REAL(DP),INTENT(OUT):: area,mac,ymac,xlemac

  REAL(DP):: frac
  REAL(DP):: span
  REAL(DP):: taper,taper1
!-------------------------------------------------------------------------------
  span=c2%y - c1%y
  area=0.5*(c1%c + c2%c)*span
  taper=c2%c/c1%c
  taper1=taper+1.0
  frac=(taper+taper1)/(3.0*taper1)
  mac=c1%c*(taper*taper + taper + 1.0)/(1.5*taper1)
  ymac=c1%y + frac*span
  xlemac=c1%xle + frac*(c2%xle-c1%xle)
  RETURN
END Subroutine MACofOneSegment   ! ---------------------------------------------

!+
SUBROUTINE MACofWing(c,area,mac,ymac,xlemac)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the mean aerodynamic chord (MAC) of a wing made defined 
!  by an array of chords.
  TYPE(WingChord),INTENT(IN),DIMENSION(:):: c  ! the chords defining the wing
  REAL(DP),INTENT(OUT):: area    ! area of the wing
  REAL(DP),INTENT(OUT):: mac     ! length of the MAC
  REAL(DP),INTENT(OUT):: ymac    ! y-coordinate of the MAC
  REAL(DP),INTENT(OUT):: xlemac  ! x-coordinate of the leading edge of the MAC

  REAL(DP):: aseg   ! area of one segment
  REAL(DP):: cseg   ! length of MAC of one segment
  REAL(DP):: xseg,yseg  ! x-le and y of MAC of one segment
  INTEGER:: k
!-------------------------------------------------------------------------------
  area=0.0
  mac=0.0
  ymac=0.0
  xlemac=0.0

  DO k=2,SIZE(c)
    CALL MACofOneSegment(c(k-1),c(k),aseg,cseg,yseg,xseg)
    area=area+aseg
    mac=mac+cseg*aseg
    ymac=ymac+yseg*aseg
    xlemac=xlemac+xseg*aseg
  END DO

  IF (area <= 0.0) RETURN
  
  mac=mac/area
  ymac=ymac/area
  xlemac=xlemac/area
  RETURN
END Subroutine MACofWing   ! ---------------------------------------------------


END Module MeanAerodynamicChordProcedures   ! ==================================

!+
PROGRAM MeanAerodynamicChord
! ------------------------------------------------------------------------------
USE ISO_FORTRAN_ENV
USE MeanAerodynamicChordProcedures

IMPLICIT NONE


  CHARACTER(LEN=*),PARAMETER:: GREETING = " Mean Aerodynamic Chord Calculator"
  CHARACTER(LEN=*),PARAMETER:: AUTHOR = &
    " Ralph Carmichael, Public Domain Aeronautical Software"
  CHARACTER(LEN=*),PARAMETER:: VERSION = " Version 1.0 (2021 November 16)"

  INTEGER,PARAMETER:: MAXCHORDS = 200
  TYPE(WINGCHORD),DIMENSION(MAXCHORDS):: ch
  CHARACTER(LEN=15):: dateTimeString
  INTEGER:: errCode = -1
  CHARACTER(LEN=80):: fileName
  CHARACTER(LEN=*),PARAMETER:: FMT = '(T3,A,T20,F15.4)'
  INTEGER,PARAMETER:: IDAT=1, OUT=2
  INTEGER:: k
  INTEGER:: n ! the number of chords read successfully from the input file
  REAL(DP):: area,mac,ymac,xlemac
  REAL(DP):: y,xle,c
  REAL(DP):: time1,time2

  NAMELIST /CHORD/ y,xle,c
  INTRINSIC:: TRIM
!-------------------------------------------------------------------------------

  CALL CPU_TIME(time1)
  WRITE(*,'(A/A/A)') GREETING,AUTHOR,VERSION
  WRITE(*,'(3A/A)') 'This program was compiled by ', compiler_version(),   &
                    ' using the options ',           compiler_options()
  
  k=COMMAND_ARGUMENT_COUNT()
  IF (k > 0) THEN   ! get fileName from the command line 
    CALL GET_COMMAND_ARGUMENT(1,fileName) 
    OPEN(UNIT=IDAT, FILE=fileName, IOSTAT=errCode, &
      STATUS='OLD', ACTION='READ', POSITION='REWIND')
  END IF
  
  IF (errCode /= 0) THEN    
    DO
      WRITE(*,*) "Enter the name of the input file: "
      READ(*,'(A)') fileName
      IF (fileName == ' ') STOP
      OPEN(UNIT=IDAT, FILE=fileName, IOSTAT=errCode, &
        STATUS='OLD', ACTION='READ', POSITION='REWIND')
      IF (errCode==0) EXIT
      OPEN(UNIT=IDAT, FILE=Trim(fileName)//'.nml', IOSTAT=errCode, &
        STATUS='OLD', ACTION='READ', POSITION='REWIND')
      IF (errCode==0) EXIT  
      OPEN(UNIT=IDAT, FILE=Trim(fileName)//'.inp', IOSTAT=errCode, &
        STATUS='OLD', ACTION='READ', POSITION='REWIND')
      IF (errCode==0) EXIT
      WRITE(*,*) "Unable to open this file. Try again."
    END DO
    
  END IF 
  
  INQUIRE(UNIT=IDAT, NAME=fileName)
  WRITE(*,*) "Reading from "//Trim(fileName)
  
  n=0
  DO k=1,MAXCHORDS
    READ(IDAT,CHORD,IOSTAT=errCode)
    IF (errCode /= 0) EXIT
    n=n+1
    ch(n)%y=y
    ch(n)%xle=xle
    ch(n)%c=c
  END DO
  
  IF (n < 2) THEN
    WRITE(*,*) 'At least two chords must be defined. ABORT.'
    STOP
  END IF

  OPEN(UNIT=OUT, FILE='mac.out', STATUS='REPLACE', ACTION='WRITE')
  WRITE(OUT,'(3A/A)') 'This program was compiled by ', compiler_version(),   &
                    ' using the options ',           compiler_options()
                    
  dateTimeString = GetDateTimeStr()
  WRITE(OUT,*) "Date and time of execution " // dateTimeString
  WRITE(OUT,'(T2,A,A)') 'MEAN AERODYNAMIC CHORD DEFINED BY ', TRIM(fileName)
  CALL MACofWing(ch(1:n), area,mac,ymac,xlemac)
  WRITE(OUT,'(T3,A)') "#        y        xle          c          xte"
  DO k=1,n
    WRITE(OUT,'(I3,4F12.5)')  k,ch(k),ch(k)%xle+ch(k)%c
  END DO
  WRITE(OUT,*) ' '
  WRITE(OUT,FMT) 'area of wing=   ', area
  WRITE(OUT,FMT) 'length of MAC=  ', mac
  WRITE(OUT,FMT) 'y of MAC=       ', ymac
  WRITE(OUT,FMT) 'xle of MAC=     ', xlemac
  WRITE(OUT,FMT) 'xte of MAC=     ', xlemac+mac
  WRITE(OUT,FMT) 'x of c/4 of MAC=', xlemac+0.25*mac
  WRITE(*,FMT) 'length of MAC=    ', mac
  WRITE(*,FMT) 'x of c/4 of MAC=  ', xlemac+0.25*mac
  WRITE(*,*) 'File mac.out added to your directory'
  WRITE(OUT,*) 'Values are for right wing. For symmetrical configuration, ' // &
    'double the area (but NOT the length)'
  CALL CPU_TIME(time2)
  WRITE(*,'(A,I0,A)') "Program normal termination, elapsed time = ", &
    INT(1000.0_DP*(time2-time1)), " milliseconds."
  STOP

END Program MeanAerodynamicChord   ! ===========================================
