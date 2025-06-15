!+
MODULE Parse
! ---------------------------------------------------------------------------
! PURPOSE - Several procedures that enable parsing of input commands

! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS RERSON  STATEMENT OF CHANGES
! 25Apr98  0.1   RLC   Extracted from other programs

  IMPLICIT NONE
  CHARACTER(LEN=*),PUBLIC,PARAMETER:: PARSE_VERSION = "0.2 (25 April 1998)"

  PUBLIC:: GetToken
  PUBLIC:: ParseLine
  PUBLIC:: ScanString
  PUBLIC:: TestKeyWord
  PUBLIC:: UpCase

CONTAINS

!+
SUBROUTINE GetToken(s,token)
! ---------------------------------------------------------------------------
! PURPOSE - Return the first token in string. This will be the characters
!    from the first non-blank character up to but not including the next
!    blank character on the line. If s contains nothing but blanks, then
!    token will also be blank.
  CHARACTER(LEN=*),INTENT(IN):: s
  CHARACTER(LEN=*),INTENT(OUT):: token

  INTEGER:: nfound
  INTEGER,DIMENSION(1):: startToken,endToken

! NOTE - If the length of the character variable token is not enough to
!   hold the token from s, it is simply truncated. This is not considered
!   to be an error and no indication is returned to the calling program.
!   There are a number of programs that use this feature; for example,
!   they wish to see the first four characters of the first token of the
!   command line. You will have to make an alternate version of this
!   subroutine if you want to return a code that flags this condition.
!----------------------------------------------------------------------------

  CALL ParseLine(s, nfound, startToken,endToken)
  IF (nfound==0) THEN
    token=" "
  ELSE
    token=s(startToken(1):endToken(1))
  END IF
  RETURN
END Subroutine GetToken   ! -------------------------------------------------

!+
SUBROUTINE ParseLine(string, nfound, startToken,endToken)
! ---------------------------------------------------------------------------
! PURPOSE - Parse a line to delimit the tokens
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!          from Sergio Gelato
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 19Apr82  0.1   RLC   Original coding
!  6May83  1.0   RLC   Replaced SCAN_STRING with SCAN
!  4Dec84  1.1   RLC   Added test for INC>NFOUND
! 28Jul86  1.2   RLC   IMPLICIT NONE
!  2Dec91  1.3   RLC   END DO
!  4Feb95  1.4   RLC   Renamed Scan to ScanString
!  1Jun95  1.5   RLC   Added LenString (N2 parameter in ScanString)
!                         SCAN is an intrinsic in Fortran 90
! 10Jul96  1.6   RLC   Recoded, inspired by a code fragment posted to
!                         comp.lang.fortran by Sergio Gelato
!                         gelato@astrosun.tn.cornell.edu
!                      Also, omitted nmax from args. Use length of start,end

  CHARACTER(LEN=*),INTENT(IN):: string   ! character variable to be analyzed
  INTEGER,INTENT(OUT):: nfound           ! number of tokens found in string
  INTEGER,INTENT(OUT),DIMENSION(:):: startToken ! starting position 
                                                ! (in string) of each token
  INTEGER,INTENT(OUT),DIMENSION(:):: endToken ! ending position (in string)
                                              ! of each token
                               ! token i is string(startToken(i):endToken(i))

  INTEGER:: k,nmax
  LOGICAL:: inword
!----------------------------------------------------------------------------
  inword=.FALSE.
  nfound=0
  nmax=MIN(SIZE(startToken), SIZE(endToken))
  IF (nmax <= 0) RETURN

  DO k=1,LEN_TRIM(string)
    IF (string(k:k) == " ") THEN
      IF (inword) THEN
        endToken(nfound)=k-1
        IF (nfound == nmax) RETURN
        inword=.FALSE.
      END IF
    ELSE
      IF (.NOT.inword) THEN
        IF (nfound+1 > nmax) RETURN
        nfound=nfound+1
        startToken(nfound)=k
        inword=.TRUE.
      END IF
    END IF
  END DO
  IF (inword) endToken(nFound)=LEN_TRIM(string)

  RETURN
END Subroutine ParseLine   ! ------------------------------------------------

!+
SUBROUTINE ScanString(s, n1,n2,n3)
! ---------------------------------------------------------------------------
! PURPOSE - Scan the characters in a string, and denote the first
!    non-blank character, the last non-blank character and the first
!    non-blank following the first blank that follows a non-blank.
  CHARACTER(LEN=*),INTENT(IN):: s
  INTEGER,INTENT(OUT):: n1,n2,n3

  INTEGER:: nfound
  INTEGER,DIMENSION(2):: startToken,endToken

!----------------------------------------------------------------------------
  n2=LEN_TRIM(s)
  CALL ParseLine(s, nfound, startToken,endToken)
  SELECT CASE(nfound)
    CASE(0)
      n1=0
      n3=0
    CASE(1)   ! there is no second token
      n1=startToken(1)
      n3=0
    CASE(2)
      n1=startToken(1)
      n3=endToken(1)
  END SELECT
  RETURN
END Subroutine ScanString   ! -----------------------------------------------

!+
FUNCTION TestKeyWord(test,dict) RESULT(k)
! ---------------------------------------------------------------------------
  CHARACTER(LEN=*),INTENT(IN):: test
  CHARACTER(LEN=*),INTENT(IN),DIMENSION(:):: dict
  INTEGER:: k

!  CHARACTER(LEN=WORDLENGTH):: testWord   ! some compilers won't buy it
  CHARACTER(LEN=LEN(dict(1))):: testWord
  INTEGER:: j
!----------------------------------------------------------------------------
  k=0
  testWord=test         ! makes all the tests between strings of equal length
  CALL UpCase(testWord) ! assumes dict entries are all upper case

  DO j=1,SIZE(dict)
    IF (testWord == dict(j)) THEN
      k=j
      RETURN
    END IF
  END DO
  RETURN
END Function TestKeyWord   ! --------------------------------------------------

!+
SUBROUTINE UpCase(a)
! ---------------------------------------------------------------------------
! PURPOSE - Convert lower case elements of a character variable to upper
  CHARACTER(LEN=*),INTENT(IN OUT):: a
  INTEGER:: j,k
!----------------------------------------------------------------------------
  DO k=1,LEN(a)
    j=IACHAR(a(k:k))
    IF (j>96 .AND. j<123) a(k:k)=ACHAR(j-32)    ! assumes ASCII character set
  END DO
  RETURN
END Subroutine UpCase   ! ---------------------------------------------------

END Module Parse   ! ========================================================




!+
MODULE MakeWgsUtilities
! ---------------------------------------------------------------------------
! PURPOSE - Various constants, global variables, and procedures for
!   the MakeWgs program
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

IMPLICIT NONE

!   M O D U L E    C O N S T A N T S
  INTEGER,PARAMETER:: IDAT=1
  INTEGER,PARAMETER:: WGS=2
  INTEGER,PARAMETER:: DBG=3
  REAL,PARAMETER :: PI=3.14159265, TWOPI=PI+PI, HALFPI=0.5*PI, RAD=180/PI
  CHARACTER(LEN=1),PARAMETER:: APOSTROPHE = "'"

!   M O D U L E    V A R I A B L E S  
  CHARACTER(LEN=132):: fileName
  CHARACTER(LEN=80):: name

  INTEGER:: nets=0
  INTEGER:: ngrid=0
!
  INTEGER:: nxbody = 21
  REAL:: xstart  = 0.0
  REAL:: cpfract = 0.5

  REAL:: lNose   = 1.0
  REAL:: lBody   = 2.0
  REAL:: lTail   = 0.0
  REAL:: rNose   = 0.0
  REAL:: radius  = 1.0
  REAL:: rbase   = 0.0
  REAL:: xNose   = 0.0
  REAL:: yNose   = 0.0
  REAL:: zNose   = 0.0
  REAL:: yBase   = 0.0
  REAL:: zBase   = 0.0
  REAL:: yAxis   = 0.0
  REAL:: zAxis   = 0.0
  INTEGER:: bodyShape = 1
  INTEGER:: nfs       = 11
  INTEGER:: ntheta    = 5
  LOGICAL:: makeBase = .FALSE.
  LOGICAL:: asym     = .FALSE.

  REAL:: rootle = 0.0
  REAL:: rootte = 1.0
  REAL:: rooty  = 0.0
  REAL:: rootz  = 0.0
  REAL:: tiple  = 0.0
  REAL:: tipte  = 1.0
  REAL:: tipy   = 1.0
  REAL:: tipz   = 0.0
  REAL:: tcroot = 0.0
  REAL:: tctip  = 0.0
  INTEGER:: rows     = 4
  INTEGER:: cols     = 4
  INTEGER:: iroot    = 0
  INTEGER:: itip     = 0
  INTEGER:: ispan    = 0
  INTEGER:: sect     = 1
  LOGICAL:: lefttip  = .FALSE.
  LOGICAL:: righttip = .FALSE.
!----------------------------------------------------------------------------


CONTAINS


!+
PURE FUNCTION BodyRadius(bodyShape, x) RESULT(r)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: bodyShape  ! 0=cone; 1=parabolic; 2=Sears-Haack
                                  ! 3=vonKarman ogive; 4=ellipsoid
  REAL,INTENT(IN):: x  ! between 0 and 1
  REAL:: r
  REAL:: area
!----------------------------------------------------------------------------

  SELECT CASE(bodyShape) 
    CASE(1)                                                  ! parabolic body
      r=x*(2.0-x)  
    CASE(2)                                               ! Sears-Haack (L-V)
      r=(x*(2.0-x))**0.75  
    CASE(3)                                                ! Von Karman ogive
      area=2.0*(ASIN(SQRT(x)) - (1.0-x-x)*SQRT(x*(1.0-x)))
      r=SQRT(area/PI)
    CASE(4)                                                       ! ellipsoid
      r=SQRT(x*(2.0-x))    
    CASE DEFAULT
      r=x                                                              ! cone  
  END SELECT

  RETURN
END FUNCTION BodyRadius   ! -------------------------------------------------


!+
SUBROUTINE FillArray(start,end, array, spacingCode)
! ---------------------------------------------------------------------------
! PURPOSE - fill an array from start to end. The intermediate points are
!    computed according to various spacing rules.

  REAL,INTENT(IN):: start,end
  REAL,INTENT(OUT),DIMENSION(:):: array
  INTEGER,INTENT(IN),OPTIONAL:: spacingCode  ! =2 full cosine
                                             ! =3 half cosine
                                             ! =4 half sine
                                             ! anything else = uniform spacing
                                             ! if omitted, uniform spacing

  INTEGER:: i,n
  REAL,ALLOCATABLE,DIMENSION(:):: temp
!----------------------------------------------------------------------------
  n=SIZE(array)
  IF (n <= 0) RETURN

  array(n)=end
  array(1)=start
  IF (n <= 2) RETURN

  ALLOCATE(temp(n-2))
  temp(:)= (/ (REAL(i),i=1,n-2) /)   !  1. 2. 3. ...
  temp=temp*(1.0/REAL(n-1))

  IF (Present(spacingCode)) THEN
    SELECT CASE(spacingCode)
      CASE (2)
        temp=0.5*(1.0-COS(PI*temp))       ! full cosine, dense near both ends
      CASE (3)
        temp=1.0-COS(HALFPI*temp)             ! half cosine, dense near start
      CASE (4)
        temp=SIN(HALFPI*temp)                     ! half sine, dense near end
    END SELECT
  END IF

  array(2:n-1)=start + (end-start)*temp

  DEALLOCATE(temp)

  RETURN
END Subroutine FillArray   ! ------------------------------------------------


!+
SUBROUTINE GetBodyRadiusAndZ(x,r,z)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: x
  REAL,INTENT(OUT):: r,z
!----------------------------------------------------------------------------
  IF (x < 0.0) THEN
    r=rNose                                                   ! ahead of nose
    z=zNose
  ELSE IF (x > lBody) THEN
    r=rbase                                                     ! behind tail
    z=zBase
  ELSE IF (x < lNose) THEN
    r=rNose+(radius-rNose)*BodyRadius(bodyShape,x/lNose)               ! nose
    z=zNose*(1-x/lNose)*(1-x/lNose)
  ELSE IF (x <= (lBody-lTail)) THEN
    r=radius                                             ! cylindrical region
    z=0.0
  ELSE
    r=rbase-(rbase-radius)*BodyRadius(bodyShape,(lBody-x)/lTail)  ! afterbody
    z=zBase*(1-(lBody-x)/lTail)**2
  END IF

  RETURN
END Subroutine GetBodyRadiusAndZ  ! =========================================

!+
PURE FUNCTION GetBodyRadius(x) RESULT(r)
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: x
  REAL:: r
!----------------------------------------------------------------------------
  IF (x < 0.0) THEN
    r=rNose                                                   ! ahead of nose
  ELSE IF (x > lBody) THEN
    r=rbase                                                     ! behind tail
  ELSE IF (x < lNose) THEN
    r=rNose+(radius-rNose)*BodyRadius(bodyShape,x/lNose)               ! nose
  ELSE IF (x <= (lBody-lTail)) THEN
    r=radius                                             ! cylindrical region
  ELSE
    r=rbase-(rbase-radius)*BodyRadius(bodyShape,(lBody-x)/lTail)  ! afterbody
  END IF

  RETURN
END Function GetBodyRadius  ! ==================================================

!+
PURE FUNCTION GetBodyY(x) RESULT(y)
! ------------------------------------------------------------------------------
  REAL,INTENT(IN):: x
  REAL:: y
!-------------------------------------------------------------------------------
  IF (x < 0.0) THEN
    y=yAxis                                                      ! ahead of nose
  ELSE IF (x > lBody) THEN
    y=yAxis                                                        ! behind tail
  ELSE IF (x < lNose) THEN
    y=yaxis + yNose*(1-x/lNose)*(1-x/lNose)                       ! nose
  ELSE IF (x <= (lBody-lTail)) THEN
    y=yaxis                                                  ! cylindrical region
  ELSE
    y=yaxis + yBase*(1-(lBody-x)/lTail)**2                   ! afterbody
  END IF

  RETURN
END Function GetBodyY  ! =======================================================

!+
PURE FUNCTION GetBodyZ(x) RESULT(z)
! ------------------------------------------------------------------------------
  REAL,INTENT(IN):: x
  REAL:: z
!-------------------------------------------------------------------------------
  IF (x < 0.0) THEN
    z=zAxis + zNose                                              ! ahead of nose
  ELSE IF (x > lBody) THEN
    z=zAxis + zBase                                                ! behind tail
  ELSE IF (x < lNose) THEN
    z=zAxis + zNose*(1-x/lNose)*(1-x/lNose)                       ! nose
  ELSE IF (x <= (lBody-lTail)) THEN
    z=zAxis                                                 ! cylindrical region
  ELSE
    z=zAxis + zBase*(1-(lBody-x)/lTail)**2                   ! afterbody
  END IF

  RETURN
END Function GetBodyZ  ! =======================================================

!+
FUNCTION Hypot(x,y) RESULT(z)
! ---------------------------------------------------------------------------
! PURPOSE - An alternative to SQRT(x*x+y*y) that is less likely to overflow.

  REAL,INTENT(IN):: x,y
  REAL:: z

  REAL:: xx,yy
!----------------------------------------------------------------------------
  xx=ABS(x)
  yy=ABS(y)
  IF (xx==0.0 .AND. yy==0.0) THEN
    z=0.0
    RETURN
  END IF

  IF (xx > yy) THEN
    z=xx*SQRT(1.0+(yy/xx)**2)
  ELSE
    z=yy*SQRT(1.0+(xx/yy)**2)
  END IF
  RETURN
END Function Hypot   ! ======================================================

!+
FUNCTION WingThickness(sectionCode, x) RESULT(z)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: sectionCode   ! 1=parabolic arc; 2=double wedge
                                     ! 3=30-70 hex; 4=wedge(blunt base)
                                     ! 5=NACA 4 digit sharp t.e.
  REAL,INTENT(IN):: x
  REAL:: z

  REAL,PARAMETER:: a0=1.4845, a1= -0.63, a2= -1.758, a3=1.4215, a4= -0.518
!----------------------------------------------------------------------------
  SELECT CASE (sectionCode) 
    CASE(2)
      IF (x <= 0.5) THEN                          !   /*  double wedge  */
        z=x
      ELSE
        z=1.0-x
      END IF
    CASE(3)
      IF (x < 0.3) THEN                          !        /*  30-70 hex  */
        z=1.6*x
      ELSE IF (x <= 0.7) THEN
        z=0.5
      ELSE
        z=1.6*(1.0-x)
      END IF
    CASE(4)
      z=0.5*x   ! ; break;                      /*  wedge (thick base)  */
    CASE(5)            !  a4 changed from -0.5075 to make sharp t.e.  */
      z=a0*SQRT(x)+x*(a1+x*(a2+x*(a3+x*a4)))    !          /*  NACA 000x  */
    CASE DEFAULT
      z=2.0*x*(1.0-x)   ! break;                              ! parabolic arc 
  END SELECT
  RETURN
END Function WingThickness   !===============================================


!+
SUBROUTINE LoadWingThicknessArray(sect, x, tc, t)
! ---------------------------------------------------------------------------
  INTEGER,INTENT(IN):: sect
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(IN):: tc
  REAL,INTENT(OUT),DIMENSION(:):: t

  LOGICAL:: zeroThickness
  REAL:: chord
  INTEGER:: i,n
!----------------------------------------------------------------------------
  n=SIZE(x)
  IF (n <= 0) RETURN
  chord=x(n)-x(1)
  zeroThickness=(sect <= 0) .OR. (tc == 0.0) .OR. (chord <= 0.0)
  IF (zeroThickness) THEN
    t=0
  ELSE
    DO i=1,n
      t(i)=WingThickness(sect, (x(i)-x(1))/chord)
    END DO
    t=t*tc*chord
  END IF

  RETURN
END Subroutine LoadWingThicknessArray ! =====================================

!+
SUBROUTINE CreateBody()
! ---------------------------------------------------------------------------
  CHARACTER(LEN=*),PARAMETER:: transform = " 0    0 0 0   0 0 0   1 1 1    1"  
  REAL,ALLOCATABLE,DIMENSION(:):: fs
  REAL,ALLOCATABLE,DIMENSION(:):: theta
  REAL,ALLOCATABLE,DIMENSION(:):: stheta, ctheta
  REAL,ALLOCATABLE,DIMENSION(:,:):: coor
  REAL:: r, yLocal, zLocal, z
  INTEGER:: j
!----------------------------------------------------------------------------
  WRITE(*,*) "Creating body network named "//Trim(name)
  ALLOCATE(fs(nfs), theta(ntheta+1), stheta(ntheta+1), ctheta(ntheta+1) )
  CALL FillArray(xNose, xNose+lBody, fs)
  IF (asym) THEN
    CALL FillArray(0.0, TWOPI, theta)
    stheta=SIN(theta)
    ctheta=COS(theta)
    stheta(ntheta+1)=0.0
    ctheta(1)=1.0
    ctheta(ntheta+1)= 1.0
  ELSE
    CALL FillArray(0.0, PI, theta)
    stheta=SIN(theta)
    ctheta=COS(theta)
    stheta(ntheta+1)=0.0
    ctheta(1)=1.0
    ctheta(ntheta+1)= -1.0
  END IF
  WRITE(Dbg,*) "theta"
  WRITE(Dbg,*) theta
  DEALLOCATE(theta)
  WRITE(Dbg,*) "stheta", stheta
  WRITE(Dbg,*) "ctheta", ctheta

  ALLOCATE(coor(3,ntheta+1) )
  nets=nets+1
  ngrid=ngrid+nfs*(ntheta+1)
  WRITE(Wgs,'(A)') APOSTROPHE//Trim(name)//APOSTROPHE
  WRITE(Wgs, '(3I6,A)' ) nets, nfs, ntheta+1, transform
  DO j=1,nfs
!!!    CALL GetBodyRadiusAndZ(fs(j)-xNose, r, z)
    r=GetBodyRadius(fs(j)-xnose)
    yLocal=GetBodyY(fs(j)-xnose)
    zLocal=GetBodyZ(fs(j)-xnose)
    coor(1,:)=fs(j)
    coor(2,:)=yLocal+r*stheta
    coor(3,:)=zLocal+r*ctheta
    WRITE(Wgs,'(6ES13.5)' ) coor 
  END DO
  
  IF (makeBase) THEN
    nets=nets+1
    ngrid=ngrid+2*(ntheta+1)
    WRITE(Wgs,'(A)') APOSTROPHE//Trim(name)//'-BASE'//APOSTROPHE
    WRITE(WGS,'(3I6,A)' ) nets, 2, ntheta+1, transform
    WRITE(WGS,'(6ES13.5)' ) coor
    coor(2,:) = yAxis+yBase
    coor(3,:) = zAxis+zBase
    WRITE(WGS,'(6ES13.5)' ) coor
  END IF

  DEALLOCATE(coor, ctheta,stheta,fs)
  RETURN
END SUBROUTINE CreateBody  ! ================================================


!+
SUBROUTINE CreateWing()
! ---------------------------------------------------------------------------
  CHARACTER(LEN=*),PARAMETER:: transform = " 0    0 0 0   0 0 0   1 1 1    1"
  REAL,ALLOCATABLE,DIMENSION(:):: xroot,xtip
  REAL,ALLOCATABLE,DIMENSION(:):: yspan,zspan
  REAL,ALLOCATABLE,DIMENSION(:):: cmroot,cmtip,camber
  REAL,ALLOCATABLE,DIMENSION(:):: troot,ttip,thick
  REAL,ALLOCATABLE,DIMENSION(:):: eta
  REAL,ALLOCATABLE,DIMENSION(:,:):: coor
!!!  REAL,ALLOCATABLE,DIMENSION(:,:):: coor2
  REAL:: dy,dz,ds, dyds,dzds
  LOGICAL:: thinWing, centerlineFin
  INTEGER:: i,j,k
!----------------------------------------------------------------------------
  IF ( (rows <=0) .OR. (cols <= 0) ) RETURN
  WRITE(*,*) "Creating wing network named "//Trim(name)

  ALLOCATE(xroot(rows+1), xtip(rows+1) )
  ALLOCATE(cmroot(rows+1), cmtip(rows+1), camber(rows+1) )
  ALLOCATE(troot(rows+1), ttip(rows+1), thick(rows+1) )
  ALLOCATE(eta(cols+1), yspan(cols+1), zspan(cols+1) )
  ALLOCATE(coor(3,rows+1))
!!!  ALLOCATE(coor2(3,rows+1) )
  CALL FillArray(rootle, rootte, xroot, iroot)
  CALL FillArray(tiple,  tipte,  xtip, itip)
  CALL FillArray(0.0, 1.0, eta, ispan)
  dy=tipy-rooty
  dz=tipz-rootz
  ds=Hypot(dy,dz)
  dyds=dy/ds
  dzds=dz/ds
  yspan=rooty+eta*dy
  zspan=rootz+eta*dz

  thinWing=(tcroot == 0) .AND. (tctip == 0)
  centerlineFin = (rooty==0.0) .AND. (tipy==0.0)
  IF (thinWing) THEN
    troot=0
    ttip=0
  ELSE
    CALL LoadWingThicknessArray(sect,xroot,tcroot,troot)
    CALL LoadWingThicknessArray(sect,xtip,tctip,ttip)
  END IF

  cmroot=0   ! equations later
  cmtip=0

  
!..... First construct the wing upper surface (t.e. to l.e.) ................
  IF (.NOT. thinWing .AND. .NOT. centerlineFin) THEN
    nets=nets+1
    WRITE(Wgs,'(A)') APOSTROPHE//TRIM(name)//"-UPPER"//APOSTROPHE
    WRITE(Wgs,'(3I6,A)') nets,cols+1,rows+1, transform
    DO j=1,cols+1
      thick= ttip*eta(j) + troot*(1.0-eta(j))
      camber=cmtip*eta(j) + cmroot*(1.0-eta(j))
      coor(1,:)= xtip*eta(j) + xroot*(1.0-eta(j))
      coor(2,:)= yspan(j) - (thick+camber)*dzds
      coor(3,:)= zspan(j) + (thick+camber)*dyds
!!!      coor2(1:3,1:rows+1)=coor(1:3,rows+1:1)  ! reverse order
      WRITE(Wgs, '(6ES13.5)' ) ((coor(k,i),k=1,3),i=rows+1,1,-1)
    END DO
  END IF


!..... Write the lower surface (le to te) ...................................
    nets=nets+1
    WRITE(Wgs,'(A)') APOSTROPHE//TRIM(name)//"-LOWER"//APOSTROPHE
    WRITE(Wgs,'(3I6,A)') nets,cols+1,rows+1, transform
    DO j=1,cols+1
      thick= ttip*eta(j) + troot*(1.0-eta(j))
      camber=cmtip*eta(j) + cmroot*(1.0-eta(j))
      coor(1,:)= xtip*eta(j) + xroot*(1.0-eta(j))
      coor(2,:)= yspan(j) + (thick+camber)*dzds
      coor(3,:)= zspan(j) - (thick+camber)*dyds
      WRITE(Wgs, '(6ES13.5)' ) coor
    END DO

!
!..... Next, the left tip (if applicable) ...................................
!
  IF ( lefttip .AND. (.NOT.thinWing) .AND. (tiple < tipte)) THEN
    nets=nets+1
    WRITE(Wgs,'(A)') APOSTROPHE//TRIM(name)//"-TIP(L)"//APOSTROPHE
    WRITE(Wgs,'(3I6,A)') nets, 3,rows+1, transform
    coor(1,:)=xroot
    coor(2,:)=rooty
    coor(3,:)=rootz+troot+cmroot
    WRITE(Wgs, '(6ES13.5)' ) coor
    coor(3,:)=rootz+cmroot
    WRITE(Wgs, '(6ES13.5)' ) coor
    coor(3,:)=rootz-troot+cmroot
    WRITE(Wgs, '(6ES13.5)' ) coor
  END IF

!..... And the right tip (if applicable)
  IF (righttip .AND. (.NOT.thinWing) .AND. (tiple < tipte) ) THEN
    nets=nets+1
    WRITE(Wgs,'(A)') APOSTROPHE//TRIM(name)//"-TIP(R)"//APOSTROPHE
    WRITE(Wgs,'(3I6,A)') nets, 3,rows+1, transform
    coor(1,:)=xtip
    coor(2,:)=tipy
    coor(3,:)=tipz-ttip+cmtip
    WRITE(Wgs, '(6ES13.5)' ) coor
    coor(3,:)=tipz+cmtip
    WRITE(Wgs, '(6ES13.5)' ) coor
    coor(3,:)=tipz+ttip+cmtip
    WRITE(Wgs, '(6ES13.5)' ) coor

  END IF

  DEALLOCATE(xroot, xtip, yspan, zspan)
  DEALLOCATE(cmroot, cmtip, camber)
  DEALLOCATE(troot, ttip, thick)
  DEALLOCATE(eta)
!!!  DEALLOCATE(coor2)
  DEALLOCATE(coor)

  RETURN
END Subroutine CreateWing   ! ================================================

END MODULE MakeWgsUtilities   ! =============================================


!+
PROGRAM MakeWgs
!   ---------------------------------------------------------------------------
! PURPOSE - Make simple wireframe objects
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!
! REVISION HISTORY
!   DATE VERSI PERSON STATEMENT OF CHANGES                              
!  5Jul88  0.1   RLC  Original coding (based on PANGEOM, GEOMWBODY      
!  2Nov88  0.5   RLC  Several fixes                                     
! 27Jun89  0.6   RLC  SECTION was being ignored; fixedit                
!  8Jan90  0.7   RLC  Translated to C. Several updates                  
!  9Jan90  0.8   RLC  Removed the old ROOT, TIP commands                
! 26Jan90  0.9   RLC  Added SHIFT, time-of-day. Better bombproofing     
!  3Aug90  0.91  RLC  Use strcmpi in KeyWord; fixed bug in WingThickness
!  6Aug90  0.92  RLC  Added LEFTTIP, RIGHTTIP                           
! 20Aug90  0.93  RLC  changed treatment of point ahead of nose (..tail) 
! 20Jun93  0.95  RLC  added yNose and full 2pi if > 0                   
! 20Aug93  0.97  RLC  increased precision to 6 figures in output
! 10Jul96  0.98  RLC  converted to Fortran90
! 21Jun98  0.99  RLC  replaced MakeWgsKeywords with Parse; moved DICT
! 28Dec02  1.0   RLC  use module Parse  
! 28Jul03  1.1   RLC  added BASENET 
! 13Feb07  1.2   RLC  fixed bug where YNOSE was ignored

USE MakeWgsUtilities
USE Parse

  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: VERSION  = "1.2 (13 February 2007)"

!------------------------------------------------------------------------------
  CALL Welcome()
  CALL ProcessScript()
  CLOSE(UNIT=Idat)
  REWIND(Wgs)
  WRITE(*,*) "File make.wgs added to your directory"
  STOP


CONTAINS

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
END FUNCTION GetDateTimeStr   ! ---------------------------------------------


!+
SUBROUTINE ProcessScript()
! ---------------------------------------------------------------------------
! PURPOSE -

  INTEGER,PARAMETER:: MAX_TOKENS = 2
  INTEGER,PARAMETER::NDICT=70, WORDLENGTH=8
  CHARACTER(LEN=WORDLENGTH),DIMENSION(NDICT),PARAMETER:: DICT = &
    (/                                                             &
       "BODY    ", "WING    ", "NAME    ", "ASYM    ", "SHIFT   ", &
       "RTABLE  ", "        ", "        ", "        ", "        ", & ! 10
       "NXBODY  ", "SHAPE   ", "NFS     ", "NTHETA  ", "        ", &
       "        ", "        ", "        ", "        ", "        ", & ! 20
       "ROWS    ", "COLS    ", "IROOT   ", "ITIP    ", "ISPAN   ", &
       "SECTION ", "SECT    ", "LEFTTIP ", "RIGHTTIP", "YNOSE   ", & ! 30
       "LNOSE   ", "LBODY   ", "LTAIL   ", "RNOSE   ", "RADIUS  ", &
       "RBASE   ", "XNOSE   ", "ZNOSE   ", "ZBASE   ", "XSTART  ", & ! 40
       "CPFRACT ", "        ", "BASENET ", "YAXIS   ", "ZAXIS   ", &
       "YNOSE   ", "YBASE   ", "        ", "        ", "        ", & ! 50
       "ROOTLE  ", "LEROOT  ", "ROOTTE  ", "TEROOT  ", "ROOTY   ", &
       "YROOT   ", "ROOTZ   ", "ZROOT   ", "TIPLE   ", "LETIP   ", & ! 60
       "TIPTE   ", "TETIP   ", "TIPY    ", "YTIP    ", "TIPZ    ", &
       "ZTIP    ", "TCROOT  ", "TCTIP   ", "XCAMBER ", "ZCAMBER "  /)


  CHARACTER(LEN=255):: buffer,dataString   !
  INTEGER:: errCode
  INTEGER:: nfound 
  INTEGER:: k
  CHARACTER(LEN=255):: token    !

  INTEGER,DIMENSION(MAX_TOKENS):: start,end

!----------------------------------------------------------------------------
  DO
    READ(Idat,'(A)',IOSTAT=errCode) buffer
    IF (errCode < 0) EXIT                        ! negative means end-of-file

    WRITE(Dbg,*) "Processing:"//Trim(buffer)
    IF (Len_Trim(buffer)==0) THEN
      k=0
    ELSE
      CALL ParseLine(buffer, nfound, start,end)
      token=buffer(start(1):end(1))
      k=TestKeyWord(token,DICT)
      IF (nfound > 1) dataString=buffer(start(2):end(2))
      WRITE(Dbg,*) k, " token number for:"//Trim(AdjustL(token))
    END IF

    SELECT CASE (k)
      CASE(1)                                                          ! BODY
        CALL CreateBody()                           
      CASE(2)                                                          ! WING
        CALL CreateWing()                        
      CASE(3)                                                          ! NAME
        name=ADJUSTL(buffer(END(1)+1:LEN_TRIM(buffer)))
      CASE(4)                                                          ! ASYM
        asym = .TRUE.
      CASE(5)                                                         ! SHIFT
        rootle=tiple
        rootte=tipte
        rooty=tipy
        rootz=tipz
      CASE(6)           ! RTABLE ??
      CASE(7)
      CASE(8)
      CASE(9)
      CASE(10)
      CASE(11)
        READ(dataString, *, IOSTAT=errCode) nxbody                   ! NXBODY
      CASE(12)
        READ(dataString, *, IOSTAT=errCode) bodyShape                 ! SHAPE
      CASE(13)
        READ(dataString, *, IOSTAT=errCode) nfs                         ! NFS
      CASE(14)
        READ(dataString, *, IOSTAT=errCode) ntheta                   ! NTHETA
      CASE(15)
      CASE(16)
      CASE(17)
      CASE(18)
      CASE(19)

      CASE(20)
      CASE(21)
        READ(dataString, *, IOSTAT=errCode) rows                       ! ROWS
      CASE(22)
        READ(dataString, *, IOSTAT=errCode) cols                       ! COLS
      CASE(23)
        READ(dataString, *, IOSTAT=errCode) iroot                     ! IROOT
      CASE(24)
        READ(dataString, *, IOSTAT=errCode) itip                       ! ITIP
      CASE(25)
        READ(dataString, *, IOSTAT=errCode) ispan                     ! ISPAN
      CASE(26,27)
        READ(dataString, *, IOSTAT=errCode) sect               ! SECT,SECTION
      CASE(28)                                             
        lefttip = .TRUE.                                            ! LEFTTIP
      CASE(29)
        righttip = .TRUE.                                          ! RIGHTTIP
      CASE(30)
        READ(dataString, *, IOSTAT=errCode) yNose                     ! YNOSE
      CASE(31)
        READ(dataString, *, IOSTAT=errCode) lnose                     ! LNOSE
      CASE(32)
        READ(dataString, *, IOSTAT=errCode) lBody                     ! LBODY
      CASE(33)
        READ(dataString, *, IOSTAT=errCode) lTail                     ! LTAIL
      CASE(34)
        READ(dataString, *, IOSTAT=errCode) rNose                     ! RNOSE
      CASE(35)
        READ(dataString, *, IOSTAT=errCode) radius                   ! RADIUS
      CASE(36)
        READ(dataString, *, IOSTAT=errCode) rbase                     ! RBASE
      CASE(37)
        READ(dataString, *, IOSTAT=errCode) xnose                     ! XNOSE
      CASE(38)
        READ(dataString, *, IOSTAT=errCode) zNose                     ! ZNOSE
      CASE(39)
        READ(dataString, *, IOSTAT=errCode) zBase                     ! ZBASE
      CASE(40)
        READ(dataString, *, IOSTAT=errCode) xStart                   ! XSTART
      CASE(41)
        READ(dataString, *, IOSTAT=errCode) cpFract                 ! CPFRACT
      CASE(42)
                                                                        ! xxx
      CASE(43)                                                      ! BASENET
        makeBase=.TRUE.                                         
      CASE(44)                                                        ! YAXIS
        READ(dataString, *, IOSTAT=errCode) yaxis
      CASE(45)                                                        ! ZAXIS
        READ(dataString, *, IOSTAT=errCode) zaxis
      CASE(46)                                                        ! YNOSE
        READ(dataString, *, IOSTAT=errCode) ynose
      CASE(47)                                                        ! YBASE
        READ(dataString, *, IOSTAT=errCode) yBase       
      CASE(48)
        READ(dataString, *, IOSTAT=errCode) zBase                    ! NXBODY
      CASE(49)
        READ(dataString, *, IOSTAT=errCode) zBase                    ! NXBODY
      CASE(50)
        READ(dataString, *, IOSTAT=errCode) zBase                    ! NXBODY
      CASE(51,52)
        READ(dataString, *, IOSTAT=errCode) rootle            ! ROOTLE,LEROOT
      CASE(53,54)
        READ(dataString, *, IOSTAT=errCode) rootte            ! ROOTTE,TEROOT
      CASE(55,56)
        READ(dataString, *, IOSTAT=errCode) rooty               ! ROOTY,YROOT
      CASE(57,58)
        READ(dataString, *, IOSTAT=errCode) rootz               ! ROOTZ,ZROOT
      CASE(59,60)
        READ(dataString, *, IOSTAT=errCode) tiple               ! TIPLE,LETIP
      CASE(61,62)
        READ(dataString, *, IOSTAT=errCode) tipte               ! TIPTE,TETIP
      CASE(63,64)
        READ(dataString, *, IOSTAT=errCode) tipY                  ! TIPY,YTIP
      CASE(65,66)
        READ(dataString, *, IOSTAT=errCode) tipZ                 ! TIPZ, ZTIP
      CASE(67)
        READ(dataString, *, IOSTAT=errCode) tcRoot                   ! TCROOT
      CASE(68)
        READ(dataString, *, IOSTAT=errCode) tcTip                     ! TCTIP
      CASE DEFAULT

    END SELECT
    IF (errCode /= 0) THEN
      WRITE(*,*) "Error"   ! need more than this
      errCode=0
    END IF
  END DO

  RETURN
END SUBROUTINE ProcessScript   ! --------------------------------------------

!+
SUBROUTINE Welcome()
! ---------------------------------------------------------------------------
  CHARACTER(LEN=*),PARAMETER:: GREETING = &
     "makewgs - A program to make WGS files."
  CHARACTER(LEN=*),PARAMETER:: AUTHOR  = &
     "Ralph L. Carmichael, Public Domain Aeronautical Software"

  CHARACTER(LEN=15):: dateTime
  INTEGER:: errCode
!----------------------------------------------------------------------------
  dateTime=GetDateTimeStr()
  WRITE(*,*) GREETING
  WRITE(*,*) AUTHOR
  WRITE(*,*) "Version ", VERSION

  DO
    WRITE(*,*) "Enter the name of the input file: "
    READ(*,'(A)') fileName
    IF (Len_Trim(fileName) == 0) STOP
    OPEN(UNIT=Idat, FILE=fileName, &
      IOSTAT=errCode, STATUS='OLD', ACTION='READ', POSITION='REWIND')
    IF (errCode == 0) EXIT
    OPEN(UNIT=Idat, FILE=TRIM(fileName)//".dat", &
      IOSTAT=errCode, STATUS='OLD', ACTION='READ', POSITION='REWIND')
    IF (errCode == 0) EXIT
    OPEN(UNIT=Idat, FILE=TRIM(fileName)//".mak", &
      IOSTAT=errCode, STATUS='OLD', ACTION='READ', POSITION='REWIND')
    IF (errCode == 0) EXIT
    WRITE(*,*) "Unable to open this file.  Try again"
  END DO
  INQUIRE(UNIT=IDAT,NAME=fileName)
  WRITE(*,*) "Reading from "//Trim(fileName)
    
  OPEN(UNIT=Wgs, FILE='make.wgs', &
    IOSTAT=errCode, STATUS='REPLACE', ACTION='READWRITE', POSITION='REWIND')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open output WGS file"
    STOP
  END IF
  WRITE(Wgs,*) "Created by makewgs from "//Trim(fileName)//" on "//dateTime

  OPEN(UNIT=Dbg, FILE='make.dbg', &
    IOSTAT=errCode, STATUS='REPLACE', ACTION='WRITE', POSITION='REWIND')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open debugging file"
    STOP
  END IF
  WRITE(Dbg,*) "Created by makewgs from "//Trim(fileName)//" on "//dateTime

  RETURN
END Subroutine Welcome   ! --------------------------------------------------

END Program MakeWgs   ! =====================================================

