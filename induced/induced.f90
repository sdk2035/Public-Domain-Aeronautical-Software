INCLUDE 'id.f90'
!+
!PROGRAM InducedDragFromSpanLoading                     ! \induced\induced.f90
! ---------------------------------------------------------------------------
! PURPOSE - Compute lift and induced drag from sparse span load data.
! AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   Oct77  0.1   RLC   Original coding (FORTRAN)
!   May81  0.2   RLC   Extensively revised
! 09Jan86  0.3   RLC   New graphics
! 08Feb93  0.4   RLC   Translation to Pascal                              
! 23Jul96  0.5   RLC   Translation to Fortran 90
! 24Dec96  0.6   RLC   Restructured the modules
! 28Dec96  0.7   RLC   Cosmetic cleanup. No graphics
! 01Jan05  0.8   RLC   Added INCLUDE 'id.f90' as first line
! 28Jan09  0.9   RLC   Labelled outer DO loop in ProcessCommands
!                      Created ProgramProcedures 
!                                                                         
! NOTES-Based on note by J.L.Lundry, J.Aircraft, March 1977, p.309
!


!+
MODULE InducedKeywords
! ---------------------------------------------------------------------------
! PURPOSE - Define the dictionary used to parse the input file.
!    Provide a public function called TestKeyWord that will compare a
!    character variable to each string in the dictionary and return the
!    index of the matching string, if any.
!    Provide a public function called ParseLine that will identify the
!    tokens (maximal strings with no embedded blanks) in a character
!    variable ( similar to strtok in C or C++).
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

IMPLICIT NONE

  INTEGER,PARAMETER,PRIVATE::NDICT=6, WORDLENGTH=7
  CHARACTER(LEN=WORDLENGTH),DIMENSION(NDICT),PRIVATE:: dict = &
    (/ "TITLE  ", "SREF   ", "SPAN   ", "TERMS  ", "LOADING", "ASYM   " /)
!----------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE ParseLine (string, nfound, startToken, endToken)
! ---------------------------------------------------------------------------
! PURPOSE-Parse a string to delimit the tokens
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!          recoded from a note by Sergio Gelato
! NOTES-If entire string is blank, then nfound =0
!
  CHARACTER(LEN=*),INTENT(IN):: string   ! character variable to be analyzed
  INTEGER,INTENT(OUT):: nfound           ! number of tokens found
  INTEGER,INTENT(OUT):: startToken(:)    ! starting position of each token
  INTEGER,INTENT(OUT):: endToken(:)      ! ending position of each token.
!                          i.e. token i is string(startToken(i):endToken(i))
  INTEGER:: i
  LOGICAL:: inword
!------------------------------------------------------------------------------
  inword=.FALSE.
  nfound=0
  DO i=1,Len_Trim(string)
    IF (string(i:i)==" ") THEN
      IF (inword) THEN
        endToken(nfound)=i-1
        inword=.FALSE.
      END IF
    ELSE 
      IF (.NOT.inword) THEN
        nfound=nfound+1
        IF (nfound > SIZE(startToken)) RETURN   ! don't take any more 
        IF (nfound > SIZE(endToken)  ) RETURN   ! than you can eat
        startToken(nfound)=i                 
        inword=.TRUE.
      END IF
    END IF
  END DO
  IF (inword) endToken(nfound)=Len_Trim(string)


  RETURN
END Subroutine ParseLine  ! -------------------------------------------------

!+
SUBROUTINE UpCase(a)
! ---------------------------------------------------------------------------
! PURPOSE - Make all elements of a character variable upper case
!
  CHARACTER(LEN=*), INTENT(IN OUT):: a
  INTEGER:: i,j
!----------------------------------------------------------------------------
  DO i=1,LEN(a)
    j=IACHAR(a(i:i))
    IF (j>96 .AND. j<123) a(i:i)=ACHAR(j-32)  ! assumes the ASCII char. set
  END DO
  RETURN
END Subroutine UpCase  ! ----------------------------------------------------

!+
FUNCTION TestKeyWord (test) RESULT (k)
! ---------------------------------------------------------------------------
! PURPOSE - Compare the variable test to each word in the dictionary. If a
!    match is found, return the index. If not, return 0.
  CHARACTER (LEN=*),INTENT(IN):: test
  INTEGER:: k

  CHARACTER (LEN=WORDLENGTH):: testWord
  INTEGER:: i
!----------------------------------------------------------------------------
  k=0
  testWord=test         ! makes all the tests between strings of equal length
  CALL UpCase(testWord)     ! make it upper case for the test
  DO i=1,NDICT
    IF (testWord == dict(i)) THEN
      k=i
      RETURN
    END IF
  END DO
  RETURN
END Function TestKeyWord   ! ------------------------------------------------

END Module InducedKeywords   ! ==============================================

!+
MODULE ProgramProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Collect the various procedures for InducedDragFromSpanLoading

IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: VERSION ="0.9 (28 January 2009)"

  PUBLIC:: ErrorHandler
  PUBLIC:: GetDateTimeStr
  PUBLIC:: PrintCoefficients
  PUBLIC:: Welcome
!-------------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE ErrorHandler(n,k)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: n,k
!----------------------------------------------------------------------------
  WRITE(*,*) "Error", n, k
  RETURN
END Subroutine ErrorHandler   ! ---------------------------------------------

!+
FUNCTION GetDateTimeStr() RESULT(s)
! ---------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time
!   This is now standard Fortran. It should work on ALL compilers!
!   You can change the first I2.2 below to I2 if you don't like a leading
!   zero on early morning times, i.e., 6:45 instead of 06:45
  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: MONTH='JanFebMarAprMayJunJulAugSepOctNovDec'
  CHARACTER(LEN=*),PARAMETER:: FMT = '(I2.2,A1,I2.2,I3,A3,I4)'
  CHARACTER(LEN=15):: s
  INTEGER,DIMENSION(8):: v
!----------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)

  WRITE(s,FMT) v(5), ':', v(6), v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
  RETURN
END Function GetDateTimeStr   ! ---------------------------------------------

!+
SUBROUTINE PrintCoefficients(efu,sym,a)
! ------------------------------------------------------------------------------
! PURPOSE - 
  INTEGER,INTENT(IN):: efu   ! external file unit for output
  LOGICAL,INTENT(IN):: sym
  REAL,INTENT(IN),DIMENSION(:):: a

  INTEGER:: i,j
!-------------------------------------------------------------------------------
  WRITE(efu,*) "Fourier Coefficients of Sine Series"
  DO i=1,SIZE(a)
    IF (sym) THEN
      j=i+i-1
    ELSE
      j=i
    END IF
    WRITE(efu, '(I4,F12.7)' ) j, a(i)
  END DO
  RETURN
END Subroutine PrintCoefficients   ! -------------------------------------------


!+
SUBROUTINE Welcome(efu1,efu2)
! ------------------------------------------------------------------------------
! PURPOSE -
  INTEGER,INTENT(IN):: efu1,efu2
  CHARACTER(LEN=*),PARAMETER:: GREETING = &
     "induced - A program to compute lift and induced drag."
  CHARACTER(LEN=*),PARAMETER:: AUTHOR = &
     " Ralph L. Carmichael, Public Domain Aeronautical Software"

  CHARACTER(LEN=15):: dateTimeStr
  INTEGER:: errCode
  CHARACTER(LEN=80):: fileName
!----------------------------------------------------------------------------
  dateTimeStr=GetDateTimeStr()
  WRITE(*,*) "It is now "//dateTimeStr
  WRITE(*,*) GREETING
  WRITE(*,*) "Version "//VERSION
  WRITE(*,*) AUTHOR

  DO
    WRITE(*,*) "Enter the name of the input file:"
    READ(*, '(A)') fileName
    IF (fileName .EQ. "") STOP
    OPEN(UNIT=efu1, FILE=fileName, &
       IOSTAT=errCode, STATUS='OLD', ACTION='READ', POSITION='REWIND')
    IF (errCode ==0) EXIT
    OPEN(UNIT=efu1, FILE=Trim(fileName)//".ccl", &
       IOSTAT=errCode, STATUS='OLD', ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) "Unable to open this file. Try again."
  END DO
  INQUIRE(UNIT=efu1, NAME=fileName)   ! gets the .ccl if you omitted it
  WRITE(*,*) 'Reading from '//TRIM(fileName)

  OPEN(UNIT=efu2, FILE="induced.out", &
       IOSTAT=errCode, STATUS='REPLACE', ACTION='WRITE', POSITION='REWIND')
  IF (errCode==0) THEN
    WRITE(efu2,*) "INDUCED DRAG FROM SPAN LOADING"
    WRITE(efu2,*) "   Program version "//VERSION
    WRITE(efu2,*) "   Input file: ",Trim(fileName)
    WRITE(efu2,*) "   Computed on "//dateTimeStr
    WRITE(efu2,*) " "
  ELSE
    WRITE(*,*) "Unable to open output file"
    STOP
  END IF

  RETURN
END Subroutine Welcome   ! -----------------------------------------------------

END Module ProgramProcedures   ! ===============================================

!+
PROGRAM InducedDragFromSpanLoading                     ! \induced\induced.f90
! ---------------------------------------------------------------------------
! PURPOSE - Compute lift and induced drag from sparse span load data.
! AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   Oct77  0.1   RLC   Original coding (FORTRAN)
!   May81  0.2   RLC   Extensively revised
! 09Jan86  0.3   RLC   New graphics
! 08Feb93  0.4   RLC   Translation to Pascal                              
! 23Jul96  0.5   RLC   Translation to Fortran 90
! 24Dec96  0.6   RLC   Restructured the modules
! 28Dec96  0.7   RLC   Cosmetic cleanup. No graphics
! 01Jan05  0.8   RLC   Added INCLUDE 'id.f90' as first line
!                                                                         
! NOTES-Based on note by J.L.Lundry, J.Aircraft, March 1977, p.309
!
USE InducedKeywords
USE InducedDrag
USE ProgramProcedures

IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: FAREWELL = &
     "File induced.out added to your directory."
  INTEGER,PARAMETER:: IDAT=1, IPRT=2
  REAL,PARAMETER:: PI=3.14159265
  REAL,ALLOCATABLE,DIMENSION(:):: y,ccl,coeff

  REAL:: cl,cd


  LOGICAL:: symLoading=.TRUE.
  REAL:: span
  REAL:: sref
  INTEGER:: nLoad
  INTEGER:: nTerms
  CHARACTER(LEN=132):: title

!------------------------------------------------------------------------------
  CALL Welcome(IDAT,IPRT)
  CALL ProcessCommands(IDAT,IPRT)
! first, compute the drag by getting the coefficients and then
! computing the drag from the coefficients. You have to do it this
! way if you also want the lift { from coeff(1) }
  ALLOCATE(coeff(nTerms))
  CALL ComputeFourierCoefficients(y,ccl,span,symLoading,coeff)
  CALL PrintCoefficients(IPRT, symLoading,coeff)
  cl=PI*span*span*coeff(1)/sref
  cd=DragFromCoefficients(symLoading,coeff,span)/sref
  WRITE(Iprt,'(A,F9.5,A,F9.5)') '  CL=',cl,'   CD=',cd
  DEALLOCATE(coeff)
! now check that the one-step function gives the same result ...............
  IF (symLoading) THEN
    cd=SymmetricLoadingInducedDrag(y,ccl,span,nTerms)/sref
  ELSE
    cd=AsymmetricLoadingInducedDrag(y,ccl,span,nTerms)/sref
  END IF
  WRITE(IPRT, '(A,F9.5)' ) 'CD from the simple function=', cd
  WRITE(*,*) FAREWELL
  WRITE(*,*) "InducedDrag has terminated successfully."
  STOP


CONTAINS

!+
SUBROUTINE ProcessCommands(efu1,efu2)
! ---------------------------------------------------------------------------
! PURPOSE - Read the records in the input file one by one.
!    Look at the first token on each record. If it matches a keyword
!    in the dictionary, then take the appropriate action.
!    Otherwise, assume it is a comment and proceed. This program requires
!    that the LOADING keyword be the last keyword, because it is followed
!    by a set of loading pairs (unspecified number) which are read until
!    an end of file condition is encountered. I think it can tolerate
!    blank lines and comments at the end of loading data, but not any
!    more valid keywords.
  INTEGER,INTENT(IN):: efu1,efu2

  INTEGER,PARAMETER:: MAX_TOKENS = 2

  CHARACTER(LEN=132):: buffer, token, dataString
  INTEGER,DIMENSION(MAX_TOKENS):: startToken,endToken
  INTEGER:: i, nFound
  INTEGER:: k
  INTEGER:: errCode
  REAL,ALLOCATABLE,DIMENSION(:):: temp1,temp2
!-------------------------------------------------------------------------------
  span=1.0
  sref=1.0
  title=" "
  nTerms=0
  ALLOCATE(temp1(1000), temp2(1000) )                ! most I can imagine
  outer: DO
    READ(efu1, '(A)', IOSTAT=errCode) buffer
    IF (errCode < 0) EXIT                             ! End of file test

    IF (Len_Trim(buffer)==0) THEN
      k=0
    ELSE
      CALL ParseLine(buffer, nfound, startToken, endToken)
      IF (nfound==0) THEN
        k=0
      ELSE
        token=buffer(startToken(1):endToken(1))
        k=TestKeyWord(token)
        IF (nfound > 1) dataString=buffer(startToken(2):endToken(2))
      END IF
    END IF
 
    SELECT CASE (k)
      CASE (1)                                                           ! TITLE
        title = buffer(endToken(1)+1:)
      CASE (2)                                                           !  SREF
        READ(dataString, *, IOSTAT=errCode) sref
      CASE (3)                                                            ! SPAN
        READ(dataString, *, IOSTAT=errCode) span
      CASE (4)                                                           ! TERMS
        READ(dataString, *, IOSTAT=errCode) nTerms
      CASE (5)                                                         ! LOADING
        nLoad=0           ! LOADING must come last because you read until EOF
        DO
          READ(efu1, '(A)', IOSTAT=errCode) buffer
          IF (errCode < 0) EXIT outer  ! jump out of both DO-loops
          IF (errCode/=0) CYCLE
          CALL ParseLine(buffer, nfound, startToken, endToken)
          IF (nfound < 2) CYCLE
          dataString=buffer(startToken(1):endToken(1))
          READ(dataString, *, IOSTAT=errCode) temp1(nLoad+1)
          IF (errCode /= 0) CYCLE
          dataString=buffer(startToken(2):endToken(2))
          READ(dataString, *, IOSTAT=errCode) temp2(nLoad+1)
          IF (errCode /= 0) CYCLE
          nLoad=nLoad+1           ! don't update it until 2 fields read OK
        END DO
      CASE (6)                                                            ! ASYM
        symLoading = .FALSE.
    END SELECT
    dataString=" "
  END DO outer

  Close(efu1)
  WRITE(*,*) 'Input data has been read.'
  IF (nTerms==0) nTerms=nLoad+nLoad   ! probably as good as anything
  ALLOCATE(y(nLoad), ccl(nLoad))

  y=temp1(1:nLoad)
  ccl=temp2(1:nLoad)
  DEALLOCATE(temp1,temp2)  ! get rid of these big guys before you
                           ! allocate memory for the matrix
 
  WRITE(efu2,*) "Summary of Input Data"
  WRITE(efu2,*) title
  WRITE(efu2,*) "span=", span, "   area=",sref, "   AR=", span*span/sref
  WRITE(efu2,*) "            y           ccl          eta"
  DO i=1,nLoad
    WRITE(efu2, '(I4,3F13.5)' ) i, y(i), ccl(i), y(i)/span
  END DO

  RETURN
END Subroutine ProcessCommands   ! ---------------------------------------------

END Program InducedDragFromSpanLoading   ! =====================================

