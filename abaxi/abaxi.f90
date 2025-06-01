!+
! PROGRAM AbAxi (NASA Langley Program D2430, LAR-11049 )
! ------------------------------------------------------------------------------
! PURPOSE - Analyze the transient response of an ablating axisymmetric body, 
!  including the effect of shape change. 
! AUTHORS - Lona M. Howser, NASA Langley Research Center
!           Stephen S. Tompkins, NASA Langley Research Center
!           Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1971   1.0   LMH   Original coding at NASA Langley
!   1995   1.1   RLC   Acquisition of source code from COSMIC
!   1998   1.2   RLC   Conversion to F90 free form source
!   1999   1.21  RLC   Compilation with Fortran 90; minor changes
! 21Jul08  1.3   RLC   Reordered to put all procedures in modules
! 05Sep08  1.4   RLC   Removed COMMON statements; INTRINSIC NONE
! 06Sep08  1.5   RLC   Compiled and ran with -chk option under LF95.
                        ! corrected some problems with undefined vars.
! 08Sep08  1.6   RLC   Created subroutine InitializeInputVariables
! 25Sep08  1.7   RLC   Made 1 the default of all interpolation orders.
! 26Sep08  1.71  RLC   Added the SAVE attribute to module AbaxiProcedures.
! 02Feb09  1.75  RLC   Changed output file names back to abaxi

! NOTES -
!   2 Feb 09 - still problems with test case2

!  Transient Response of Ablating Axisymmetric Bodies Including the 
!   Effects of Shape Change 
! ( NASA Langley Research Center Program D2430, LAR-11049)
!
! A computer program has been developed to analyze the transient 
! response of an ablating axisymmetric body, including the effect of 
! shape change. The governing differential equation, the boundary 
! conditions for the analysis on which the computer program is based, 
! and the method of solution of the resulting finite-difference 
! equations are discussed in the documentation.
 
! Some of the features of the analysis and the associated program are:
!  (1) the ablation material is considered to be orthotropic with 
! temperature-dependent thermal properties; 
!  (2) the thermal response of the entire body is considered 
! simultaneously; 
!  (3) the heat transfer and pressure distribution over the body are 
! adjusted to the new geometry as ablation occurs; 
!  (4) the governing equations and several boundary-condition options 
! are formulated in terms of generalized orthogonal coordinates for 
! fixed points in a moving coordinate system; 
!  (5) the finite-difference equations are solved implicitly; and 
!  (6) other instantaneous body shapes can be displayed with a
! user-supplied plotting routine.
 
! The physical problem to be modeled with the analysis is described
! by FORTRAN input variables. For example, the external body geometry is 
! described in the W, Z coordinates; material density is given; and the 
! stagnation cold-wall heating rate is given in a time-dependent array. 
! Other input variables are required which control the solution, specify 
! boundary conditions, and determine output from the program. The 
! equations have been programmed so that either the International System 
! of Units or the U. S. Customary Units may be used.
 
! REFERENCES
! 1. Howser, Lona M.; and Tompkins, Stephen S.: Computer Program for the 
! Transient Response of Ablating Axisymmetric Bodies Including the Effects of
! Shape Change. NASA TM X-2375, 1971.
! 2. Tompkins, Stephen S.; Moss, James N.; Pittman, Claud M.; and Howser, 
! Lona M.: Numerical Analysis of the Transient Response of Ablating Axisymmetric 
! Bodies Including the Effects of Shape Change. NASA TN D-6220, 1971.
! 3. Gavril, Bruce D.; and Lane, Frank: Finite Difference Equation and Their 
! Solution for the Transient Temperature Distribution in Composite, Anisotropic,
! Generalized Bodies of Revolution. Tech. Rep. No. 230 (Contract No. NOrd 18053),
! Gen. Appl. Sci, Lab., Inc., May 26, 1961.
! 4. Hovanessian, Shahen A.; and Pipes, Louis A.: Digital Computer Methods in 
! Engineering. McGraw-Hill Book Co., Inc., c.1969.

!-------------------------------------------------------------------------------
!!!! include 'tableLookup.f90'

!+
MODULE AbaxiConstants
! ------------------------------------------------------------------------------
! PURPOSE -  Collect the global constants
  INTEGER,PARAMETER:: DBG = 3, GNU1 = 12, GNU2=13, GNU3=14
!  INTEGER,PARAMETER:: WP = SELECTED_REAL_KIND(6,20)
  INTEGER,PARAMETER:: WP = SELECTED_REAL_KIND(10,50)
  CHARACTER(LEN=*),PARAMETER:: VERSION = "1.75 (2 February 2009)"
END Module AbaxiConstants   ! ==================================================

!+
MODULE InterpolationMethods
! ------------------------------------------------------------------------------
! PURPOSE - Collect the interpolation procedures used in AbAxi
USE AbaxiConstants, ONLY: WP
IMPLICIT NONE
PUBLIC:: Ftlup
PRIVATE:: Uns
PRIVATE:: Lagran
  
CONTAINS

!+
SUBROUTINE Ftlup (x,y,m,n,vari,vard)
! ------------------------------------------------------------------------------
!     PURPOSE - Interpolation

!     NOTES-
!     ***DOCUMENT DATE 09-12-69    SUBROUTINE REVISED 07-07-69 *********
!        MODIFICATION OF LIBRARY INTERPOLATION SUBROUTINE  FTLUP


      IMPLICIT NONE
!!!      IMPLICIT REAL*8(A-H,O-Z)

  REAL(WP),INTENT(IN):: x   ! interpolation point
  REAL(WP),INTENT(OUT):: y  ! interpolation value (the output)
  INTEGER,INTENT(IN) :: m   ! order of interpolation
  INTEGER,INTENT(IN) :: n   ! size of arrays vari and vard
  REAL(WP),INTENT(IN),DIMENSION(:):: vari  ! independent variable data
  REAL(WP),INTENT(IN),DIMENSION(:):: vard  ! dependent variable data

  INTEGER:: i,j,k
  INTEGER:: in
  INTEGER:: l
  INTEGER:: li
  INTEGER:: ma
  REAL(WP),PARAMETER:: ONE=1.0
  REAL(WP)::sk
  REAL(WP),DIMENSION(3):: v
  REAL(WP),DIMENSION(2):: yy
  INTEGER,DIMENSION(43):: ii

!-------------------------------------------------------------------------------

  IF (m==0) THEN
    WRITE(*,*) 'ftlup entered with m=0'
  END IF

!      INITIALIZE ALL INTERVAL POINTERS TO -1.0   FOR MONOTONICITY CHECK
!!!DATA (ii(j),j=1,43)/43*-1/
  ii=-1
  ma=ABS(m)

!            ASSIGN INTERVAL POINTER FOR GIVEN VARI TABLE
!      THE SAME POINTER WILL BE USED ON A GIVEN VARI TABLE EVERY TIME
  li=1
  i=ii(li)
  IF (i >= 0) GO TO 6
  IF (n < 2) GO TO 6

!     MONOTONICITY CHECK
  IF (vari(2)-vari(1) > 0.0) THEN
    GO TO 4
  ELSE
    GO TO 2
  END IF

!     ERROR IN MONOTONICITY
1     k=j   ! LOCF(VARI(1))
  PRINT 17, j,k,(vari(j),j=1,n),(vard(j),j=1,n)
  STOP

!     MONOTONIC DECREASING
2     DO 3  j=2,n
  IF (vari(j)-vari(j-1) < 0.0) THEN
    CYCLE
  ELSE
    GO TO 1
  END IF
3 END DO
  GO TO 6
!     MONOTONIC INCREASING
4     DO 5 j=2,n
  IF (vari(j)-vari(j-1) > 0.0) THEN
    CYCLE
  ELSE
    GO TO     1
  END IF
5 END DO

!     INTERPOLATION
6 IF (i <= 0) i=1
  IF (i >= n) i=n-1
  IF (n <= 1) GO TO 7
  IF (ma /= 0) GO TO 8
!     ZERO ORDER
7     y=vard(1)
GO TO 16

!     LOCATE I INTERVAL (X(I).LE.X.LT.X(I+1))
8     IF ((vari(i)-x)*(vari(i+1)-x) > 0.0) THEN
  GO TO     9
ELSE
  GO TO    11
END IF
!     IN GIVES DIRECTION FOR SEARCH OF INTERVALS
9     in=INT(SIGN(one, (vari(i+1)-vari(i))*(x-vari(i))))  ! little fix by RLC
!     IF X OUTSIDE ENDPOINTS, EXTRAPOLATE FROM END INTERVAL
10    IF ((i+in) <= 0) GO TO 11
  IF ((i+in) >= n) GO TO 11
  i=i+in
  IF ((vari(i)-x)*(vari(i+1)-x) > 0.0) THEN
    GO TO    10
  END IF
11    IF (ma == 2) GO TO 12

!     FIRST ORDER
y=(vard(i)*(vari(i+1)-x)-vard(i+1)*(vari(i)-x))/(vari(i+1)-vari(i) )
GO TO 16

!     SECOND ORDER
12    IF (n == 2) GO TO 1
IF (i == (n-1)) GO TO 14
IF (i == 1) GO TO 13
!     PICK THIRD POINT
sk=vari(i+1)-vari(i)
IF ((sk*(x-vari(i-1))) < (sk*(vari(i+2)-x))) GO TO 14
13    l=i
GO TO 15
14    l=i-1
15    v(1)=vari(l)-x
v(2)=vari(l+1)-x
v(3)=vari(l+2)-x
yy(1)=(vard(l)*v(2)-vard(l+1)*v(1))/(vari(l+1)-vari(l))
yy(2)=(vard(l+1)*v(3)-vard(l+2)*v(2))/(vari(l+2)-vari(l+1))
y=(yy(1)*v(3)-yy(2)*v(1))/(vari(l+2)-vari(l))
16    ii(li)=i
RETURN


17    FORMAT (' Table below out of order for ftlup  at position ',  &
    I5,/' x table is stored in location ', I6//(8G15.8))
END Subroutine Ftlup   ! ----------------------------------------------------



!+
SUBROUTINE DISCOT (XA,ZA,TABX,TABY,TABZ,NC,NY,NZ,ANS)
! ------------------------------------------------------------------------------
! PURPOSE - Interpolation for discontinuous functions
!******** DOCUMENT DATE 08-01-68   SUBROUTINE REVISED 08-01-68 *********
!     THE DIMENSIONS IN THIS SUBROUTINE ARE ONLY DUMMY DIMENSIONS.      
IMPLICIT NONE
  REAL(WP),INTENT(IN):: xa,za
  REAL(WP),INTENT(IN),DIMENSION(:):: tabx,taby,tabz  !!    DIMENSION TABX(2),TABY(2),TABZ(2)
  INTEGER,INTENT(IN):: nc   ! some sort of code
  INTEGER,INTENT(IN):: ny,nz
  REAL(WP),INTENT(OUT):: ans

  INTEGER:: i
  INTEGER:: ia
  INTEGER:: idx,idz
  INTEGER:: ims
  INTEGER:: ip1x,ip1z
  INTEGER:: is
  INTEGER:: jj
  INTEGER:: ll
  INTEGER:: nloc, nlocy
  INTEGER:: nn, nnn
  INTEGER:: npx(8),npy(8)
  INTEGER:: npz, npzl
  INTEGER:: nx
  REAL(WP):: yy(8) 
  REAL(WP):: zarg
!     DIMENSION TABX(2),TABY(2),TABZ(2),NPX(8),NPY(8),YY(8)     
!-------------------------------------------------------------------------------        
      CALL UNS (NC,IA,IDX,IDZ,IMS) 
      IF (NZ-1)   5,5,10 
    5 CALL DISSER (XA,TABX(1:),1,NY,IDX,NN) 
      NNN=IDX+1 
      CALL LAGRAN (XA,TABX(NN:),TABY(NN:),NNN,ANS) ! modified by RLC
      GOTO 70 
   10 ZARG=ZA 
      IP1X=IDX+1 
      IP1Z=IDZ+1 
      IF (IA)   15,25,15 
   15 IF (ZARG-TABZ(NZ))   25,25,20 
   20 ZARG=TABZ(NZ) 
   25 CALL DISSER(ZARG,TABZ(1:),1,NZ,IDZ,NPZ) 
      NX=NY/NZ 
      NPZL=NPZ+IDZ 
      I=1 
      IF (IMS)   30,30,40 
   30 CALL DISSER (XA,TABX(1:),1,NX,IDX,NPX(1)) 
      DO 35 JJ=NPZ,NPZL 
      NPY(I)=(JJ-1)*NX+NPX(1) 
      NPX(I)=NPX(1) 
   35 I=I+1 
      GOTO 50 
   40 DO 45 JJ=NPZ,NPZL 
      IS=(JJ-1)*NX+1 
      CALL DISSER (XA,TABX(1:),IS,NX,IDX,NPX(I)) 
      NPY(I)=NPX(I) 
   45 I=I+1 
   50 DO 55 LL=1,IP1Z 
      NLOC=NPX(LL) 
      NLOCY=NPY(LL) 
   55 CALL LAGRAN(  XA, TABX(NLOC:), TABY(NLOCY:), IP1X, YY(LL)) 
      CALL LAGRAN(ZARG,  TABZ(NPZ:),       YY(1:), IP1Z,    ANS) 
   70 RETURN 
END Subroutine Discot   ! ------------------------------------------------------

!+
SUBROUTINE UNS (IC,IA,IDX,IDZ,IMS) 
! ------------------------------------------------------------------------------
! PURPOSE - Interpret the control number for subroutine Discot
IMPLICIT NONE
      INTEGER,INTENT(IN):: ic    ! the control number
      INTEGER,INTENT(OUT):: ia   ! =1 if hundreds digit is 1; =0 if 0
      INTEGER,INTENT(OUT):: idx  ! tens digit of ic
      INTEGER,INTENT(OUT):: idz  ! units digit of ic
      INTEGER,INTENT(OUT):: ims  ! =0 if ic > 0; =1 otherwise

      INTEGER:: nc
!-------------------------------------------------------------------------------
      IF (IC)   5,5,10 
    5 IMS=1 
      NC=-IC 
      GOTO 15 
   10 IMS=0 
      NC=IC 


   15 IF (NC-100)   20,25,25 
   20 IA=0 
      GOTO 30 
   25 IA=1 
      NC=NC-100 
   30 IDX=NC/10 
      IDZ=NC-IDX*10 
      RETURN 
END Subroutine Uns   ! ---------------------------------------------------------

!+
SUBROUTINE DISSER (XA,TAB,I,NX,ID,NPX) 
! ------------------------------------------------------------------------------
! PURPOSE -
IMPLICIT NONE
      REAL(WP),INTENT(IN):: xa
      REAL(WP),INTENT(IN),DIMENSION(:):: tab
      INTEGER,INTENT(IN):: i
      INTEGER,INTENT(IN):: nx
      INTEGER,INTENT(IN):: id
      INTEGER,INTENT(OUT):: npx

      INTEGER:: ii,jj
      INTEGER:: ndis
      INTEGER:: nloc
      INTEGER:: nlow
      INTEGER:: nl,nu
      INTEGER:: npb,npt,npu
      INTEGER:: nupp

!     DIMENSION TAB(2)                   
!-------------------------------------------------------------------------------
      NPT=ID+1 
      NPB=NPT/2 
      NPU=NPT-NPB 
      IF (NX-NPT)   10,5,10 
    5 NPX=I 
      RETURN 
   10 NLOW=I+NPB 
      NUPP=I+NX-(NPU+1) 

  DO II=NLOW,NUPP   ! was 15
    NLOC=II 
    IF (TAB(II)>=XA)   GO TO 20
  END DO 

  NPX=NUPP-NPB+1 
  RETURN 

20 CONTINUE
  NL=NLOC-NPB 
  NU=NL+ID 
  DO 25 JJ=NL,NU 
    NDIS=JJ 
    IF (TAB(JJ)==TAB(JJ+1)) GO TO 30
25 END DO 

  NPX=NL 
  RETURN 

30 IF (TAB(NDIS)-XA)   40,35,35 
   35 NPX=NDIS-ID 
      RETURN 
40 NPX=NDIS+1 
      RETURN 
END Subroutine Disser   ! ------------------------------------------------------

!+
SUBROUTINE LAGRAN (xa,x,y,n,ans)
! ------------------------------------------------------------------------------
! PURPOSE - Evaluate the unique polynomial thru {(x(i),y(i)) : i=1,n} at xa.
! NOTE - Why isn't this a function?
IMPLICIT NONE
      REAL(WP),INTENT(IN):: xa
      REAL(WP),INTENT(IN),DIMENSION(:):: x,y
      INTEGER,INTENT(IN):: n   ! length of arrays x and y
      REAL(WP),INTENT(OUT):: ans

      REAL(WP):: a,b
      INTEGER:: i,j
      REAL(WP):: sumprod,prod   ! original sum changed to sumprod by RLC
!     DIMENSION X(2),Y(2) 
!     DIMENSION X(2),Y(2)
!-------------------------------------------------------------------------------
      sumprod=0.0 
      DO  i=1,n 
        prod=y(i) 
        DO j=1,n 
          a=x(i)-x(j) 
          IF (a==0.0) CYCLE
          b=(xa-x(j))/a 
          prod=prod*b 
        END DO 
        sumprod=sumprod+prod 
      END DO
      ans=sumprod
      RETURN 
END Subroutine Lagran   ! ------------------------------------------------------


END Module InterpolationMethods   ! ============================================


!+
MODULE PickCommon
! ------------------------------------------------------------------------------
! PURPOSE - Define the variables originally held in COMMON /PICK/
USE AbaxiConstants, ONLY: WP
  REAL(WP):: A(10,20)     ! Elements in coefficient matrix for the column solution
  REAL(WP):: AA(20)       ! partial of delta w.r.t. x
  REAL(WP):: AB(10,20)    ! Elements in coefficient matrix for the row solution
  REAL(WP):: ALPHA(20)    ! alpha
  REAL(WP):: B(20)        ! Major diagonal elements in coefficient matrix
  REAL(WP):: BS1(10,20)   ! Major diagonal elements in coefficient matrix for the
                           ! column solution minus partial of T w.r.t. tau term
  REAL(WP):: BS1B(10,20)  ! Major diagonal elements in coefficient matrix for the
                           ! row solution minus partial of T w.r.t. tau term
  REAL(WP):: C(10,20)     ! Elements in coefficient matrix for the column solution
  REAL(WP):: CB(10,20)    ! Elements in coefficient matrix for the row solution
  REAL(WP):: CK(10)       ! Temporary storage used to define the thermal
                           ! conductivity at a half station
  REAL(WP):: CKETA(10,20) ! k-sub-nu at the station
  REAL(WP):: CKXI(10,20)  ! k-sub-xi at the station
  REAL(WP):: COST(20)     ! cos theta
  REAL(WP):: CP(10,20)    ! Cp
  REAL(WP):: D(10,20)     ! h2*h3*k-sub-xi/h1
  REAL(WP):: DC(20)       ! Right-hand side of the matrix solution
  REAL(WP):: DELESQ       ! (delta-eta)**2
  REAL(WP):: DELETA       ! delta-eta
  REAL(WP):: DELTA(20)    ! delta
  REAL(WP):: DELXI        ! delta-sub-xi
  REAL(WP):: DELXISQ      ! (delta-sub-xi)**2
  REAL(WP):: E(10,20)     ! h1*h3*k-sub-eta/h2
  REAL(WP):: EIGHT3       ! Constant, 8.0/3,0
  REAL(WP):: ELAM(20)     ! lambda
  REAL(WP):: ETA(10)      ! eta
  REAL(WP):: EXPG         ! Computed constant used in computing new heating
                           ! distribution
  REAL(WP):: F(10,20)     ! fraction, see TM X-2375
  REAL(WP):: GG           ! Computed constant used in computing new heating
                           ! distribution
  REAL(WP):: GIMACH       ! Computed constant used in computing new pressure
                           ! distribution
  REAL(WP):: H1(10,20)    ! h1
  REAL(WP):: H2(10,20)    ! h2
  REAL(WP):: H3(10,20)    ! h3
  REAL(WP):: HC(20)       ! delta-H-sub-S
  REAL(WP):: HCOMB(20)    ! delta-H-sub-c
  REAL(WP):: HE           ! delta-H-sub-e
  REAL(WP):: HW(20)       ! H-sub-w
  INTEGER:: IFIRST        ! 0 for first time step in calculation, 
                          ! 1 for any time after first time step
  INTEGER:: IROCOL        ! 1 for column solution, 2 for row solution
  INTEGER:: ITC           ! Number of iterations during the column solution
  INTEGER:: ITR           ! Number of iterations during the row solution
  INTEGER:: ITT           ! Number of iterations during a solution
  INTEGER:: ITTO          ! Total number of iterations from the initial time
  INTEGER:: LM1           ! Computed constant (L-1)
  INTEGER:: LM2           ! Computed constant (L-2)
  REAL(WP):: MCDOT(20)    ! m-dot-sub-c
  REAL(WP):: MDOT(22)     ! m-dot
  REAL(WP):: MSDOT(20)    ! m-dot-sub-S
  REAL(WP):: PID2         ! Constant 1.5707963268 (half-PI)
  REAL(WP):: PRELOC(20)   ! Local wall pressure
  REAL(WP):: QC(20)       ! Adjusted convective heating rate
  REAL(WP):: QC1          ! q-sub-C
  REAL(WP):: QCNET(20)    ! q-sub-C,net
  REAL(WP):: QCOMB(20)    ! Heat due to combustion for oxidation
  REAL(WP):: QR(20)       ! Adjusted radiant heating rate
  REAL(WP):: QR1          ! q-sub-r
  REAL(WP):: QS(20)       ! Net heat input
  REAL(WP):: RNS          ! Nose radius
  REAL(WP):: RODPC        ! t'' * rho'' * c-sub-p'' / delta-tau
  REAL(WP):: ROPCPP       ! t' * rho' * c-sub-p' / delta-tau
  REAL(WP):: RSS(22)      ! Coordinate used to define body geometry, w
  REAL(WP):: RSTO2        ! Computed constant, ratio of molecular weight of free
                           ! stream to molecular weight of diatomic oxygen used 
                           ! in oxidation equation
  REAL(WP):: SIG          ! Computed constant sigma*epsilon
  REAL(WP):: SIGDP        ! Computed constant sigma*epsilon''
  REAL(WP):: SIGMA        ! sigma
  REAL(WP):: SIGP         ! Computed constant sigma-epsilon'
  REAL(WP):: SINT(20)     ! sin theta
  INTEGER:: SM1           ! Computed constant (S-1)
  INTEGER:: SM2           ! Computed constant (S-2)
  REAL(WP):: TAU          ! Time at which calculation is being made
  REAL(WP):: TB           ! Temperature to which back surfaces radiate
  REAL(WP):: TT(10,20)    ! Estimated temperatures at T
  REAL(WP):: TTF(10,20)   ! 
  REAL(WP):: TWDELXI      ! Computed constant 2*delta-xi
  REAL(WP):: TWOGI        ! Computed constant used in computing new heating
                           ! distribution
  REAL(WP):: V(20)        ! Elements in coefficient matrix for column solution
  REAL(WP):: VB(10)       ! Elements in coefficient matrix for row solution
  REAL(WP):: X(22)        ! Curvilinear coordinate
  REAL(WP):: XDXISQ       ! Computed constant (x-sub-b*delta-xi)**2
  REAL(WP):: XODXI        ! Computed constant x-sub-b*delta-xi
  REAL(WP):: Y(10,20)     ! Curvilinear coordinate
  REAL(WP):: Z(20)        ! Elements in coefficient matrix for the column solution
  REAL(WP):: ZB(10)       ! Elements in coefficient matrix for the row solution

END Module PickCommon

MODULE InputsCommon
USE AbaxiConstants,ONLY: WP

  REAL(WP):: AEXP   ! Coefficient of the exponential term when the Arrhenius
                     ! expression is used for calculating mc-dot
  REAL(WP):: ALCTAB(10)   ! Aerodynamic-blocking coefficient for heat and
                     ! mass transfer associated with MCDOT, 
                     ! a function of time (TTALC)
  REAL(WP):: ALPHAT(10)  ! Absorptance of surface, a function of temperature
                     ! (TALPHA)
  REAL(WP):: ALSTAB(10)  ! Aerodynamic-blocking coefficient for heat and mass 
     ! transfer associated with MSDOT, a function of time (TTALS)
  REAL(WP):: ASEXP  !     Coefficient in the expression for calculating MSDOT.
  REAL(WP):: BETA   !     Determines whether ablation or transpiration theory
    ! will be used for effect of mass transfer on heat transfer; 
    ! for ablation theory, BETA=1; for transpiration theory,, BETA=0
  REAL(WP):: BEXP  !      Power of the exponential term in the Arrhenius
            !  expression for calculating MCDOT
  REAL(WP):: BSEXP !      Power of the exponential term in the expression
         !   for calculating MSDOT
  REAL(WP):: CE     !     Oxygen concentration, by mass, at edge of boundary layer
  REAL(WP):: CKETATB(50)  ! Thermal conductivity in n-direction, a function of
        ! eta (ETATAB) and temperature (TTCKETA) 
  REAL(WP):: CKXITAB(50)   ! Thermal conductivity in -direction, a function of
       ! xi (XITAB) and temperature (TTCIOCI)
  REAL(WP):: CORDSY  ! Trigger to indicate coordinate system; 
                     ! if curvilinear coordinates, CORDSY=0; 
                     ! if Cartesian coordinates, CORDSY=1
  REAL(WP):: CPDP   ! c-sub-p-double prime, Specific heat of layer along y=0
  REAL(WP):: CPP    ! c-sub-p-prime     Specific heat of layer along x=L
  REAL(WP):: CPTAB(10) ! c-sub-p  Specific heat, a function of temperature (TTABCP)
  REAL(WP):: DELTAO(20) ! Initial material thickness, must have L values
  REAL(WP):: DELTAU  ! Initial computing time interval
  REAL(WP):: DELTMIN  ! Minimum value allowed for DELTA
  REAL(WP):: DTMAX   ! Maximum DELTAU which can be used; if no value
       ! is given, DTMAX=2.0
  REAL(WP):: ELAMTB(28) ! lambda, the ratio of mass of material removed per unit 
                        ! mass of oxygen that reaches the surface, a function of 
                        ! pressure (PELAM) and temperature (TTELAM)
  REAL(WP):: ENDTAU   ! Time at which calculation stops
  REAL(WP):: EPSONE   ! epsilon, the emittance of front surface
  REAL(WP):: EPSONEP  ! epsilon-prime, the emittance of layer along x=L
  REAL(WP):: EPSONPP  ! epsilon-double-prime, the emittance of layer along y=0
  REAL(WP):: ERRORT   ! Acceptable relative error in temperature
  REAL(WP):: ETATAB(5) ! ETA table for CKETATB
  REAL(WP):: GAMBAR  ! Mean ratio of specific heats behind bow shock wave, used
                 ! only in computation of heating-rate distribution around body
  REAL(WP):: GAMINF  ! Ratio of specific heats in free stream, used only in 
                     ! computing heating-rate distribution around body
  REAL(WP):: HCOMBTB(28) ! delta-H-sub-c, the heat of combustion, a function of
                         ! pressure (PHCOMB) and temperature (TTHCOMB)
  REAL(WP):: HCTAB(28) ! delta-H-sub-s, the heat of sublimation, a function of 
                       ! pressure (PHC) and temperature (TTABHC)
  REAL(WP):: HETAB(10)  ! H-sub-e, the total free-stream enthalpy, a function 
                         ! of time(TTABHE)
  REAL(WP):: HWTAB(15)  ! H-sub-w, the enthalpy of gas at the wall temperature,
                        ! a function of temperature (TTABHW)

  INTEGER:: IADJUST ! Trigger for adjusting heating-rate and pressure 
                    ! distributions to shape change; if IADJUST=O, 
                    ! QRAT and PRAT are not adjusted; if IADJUST/=0, 
                    ! QRAT and PRAT will be adjusted to shape change
  INTEGER:: IPLOT   ! Trigger for plotting routine; if IPLOT=O, no plots;
                    ! if IPLOT /= 0, the following plots will be made: 
                    ! 1) RSS versus ZS at times indicated in PLTIME table; 
                    ! 2) MDOT versus x at each PRFREQ time; and
                    ! 3) T(M,N) versus x indicated in NTP array at each PREREQ
  INTEGER:: L       ! Number of stations in the x-direction
  REAL(WP):: MACHNO  ! Free-stream Mach number
  INTEGER:: MALPHA  ! Order of interpolation for ALPHAT
  INTEGER:: MALPHC  ! Order of interpolation for ALCTAB
  INTEGER:: MALPHS  ! Order of interpolation for ALSTAB
  INTEGER:: MAXITT  ! Maximum iteration count; when iteration count exceeds 
                    ! this number, DELTAU will be halved until DELTAU is less 
                    ! than 1.0E-6, then the program will stop and a message will
                    ! be printed
  INTEGER:: MCP     ! Order of interpolation for CPTAB
  REAL(WP):: MDMAX   ! Maximum expected MDOT; this must be given to get a 
                    ! reasonable scale for plots; not needed if IPLOT=0
  REAL(WP):: MDOTO(20) ! Initial mass loss rate at surface, must have L values
  INTEGER:: MHE    ! Order of interpolation for HETAB
  INTEGER:: MHW    ! Order of interpolation for HWTAB
  INTEGER:: MPSTAG ! Order of interpolation for PSTAGTB
  INTEGER:: MQC    ! Order of interpolation for QCTAB
  INTEGER:: MQR    ! Order of interpolation for QRTAB
  INTEGER:: MTB    ! Order of interpolation for TBTAB
  REAL(WP):: MWO2   ! Molecular weight of diatomic oxygen used in oxidation equation
  REAL(WP):: MWSTR  ! Molecular weight of free stream used in oxidation equation
  INTEGER:: NALPHA ! Number of entries in ALPHAT
  INTEGER:: NALPHC ! Number of entries in ALCTAB
  INTEGER:: NALPHS ! Number of entries in ALSTAB
  INTEGER:: NCKETA ! Number of entries in CKETATB
  INTEGER:: NCKXI  ! Number of entries in CKXITAB
  INTEGER:: NCP    ! Number of entries in CPTAB
  INTEGER:: NELAM  ! Number of entries in ELAMTB
  INTEGER:: NETA   ! Number of entries in ETATAB
  INTEGER:: NHC    ! Number of entries in HCTAB
  INTEGER:: NHCOMB ! Number of entries in HCOMBTB
  INTEGER:: NHE    ! Number of entries in HETAB
  INTEGER:: NHW    ! Number of entries in HWTAB
  INTEGER:: NPELAM ! Number of entries in PELAM
  INTEGER:: NPHC   ! Number of entries in PHC
  INTEGER:: NPHCOMB! Number of entries in PHCOMB
  INTEGER:: NPSTAG ! Number of entries in PSTAGTB
  INTEGER:: NQC    ! Number of entries in QCTAB
  INTEGER:: NQR    ! Number of entries in QRTAB
  INTEGER:: NTB    ! Number of entries in TBTAB
  INTEGER:: NTP(7) ! Array of seven values which specify the temperatures to 
                   ! be plotted; NTP(1) = the number of temperature rows to be 
                   ! plotted (may be six or less); NTP(2) through NTP(7), 
                   ! the row number of the temperatures to be plotted. 
                   ! For example, NTP(1)=3, NTP(2)=1, NTP(4)=10, specifies that
                   ! three (3) rows of temperature will be plotted and these 
                   ! rows are 1, 5, and 10
  INTEGER:: NXI    ! Number of entries in XITAB

  REAL(WP):: PELAM(4)   !     Pressure table for ELAMTB
  REAL(WP):: PHC(4)    ! Pressure table for HCTAB
  REAL(WP):: PHCOMB(4)  !  Pressure table for HCOMBTB
  REAL(WP):: PLTIME(15)  ! Times at which RSS versus ZS, that is, the body
                         ! shape, will be plotted; not needed if IPLOT=0
  REAL(WP):: PRAT(20)  ! Initial ratio of local to stagnation pressure, must
                       ! have L values, not needed if IADJUSWO
  REAL(WP):: PRFREQ  ! Printing time frequency for output data
  REAL(WP):: PSEXP  ! p, the exponent of pressure term in sublimation equation
  REAL(WP):: PSTAGTB(10) ! Stagnation pressure, a function of time (TTPSTAG)
  REAL(WP):: PTMAX   ! Maximum expected value of T, used to get reasonable scale 
    ! in plotting, not needed if IPLOT=0
  REAL(WP):: PTMIN   ! Minimum expected value of T, used to get reasonable scale 
    ! in plotting; not needed if IPLOT=0
  REAL(WP):: QCTAB(10) ! q-sub-C, the cold-wall convective heating rate, 
                       ! a function of time (TTABQC)
  REAL(WP):: QRAT(20)  ! Initial convective heating-rate distribution. Must
                       ! have L values, not needed if IADJUSTL0
  REAL(WP):: QRRAT(20) ! Radiant heating-rate distribution over body, must have L values
  REAL(WP):: QRTAB(10) ! q-sub-r, the radiant heating-rate tables, 
                       ! a function of time (TTABQR)
  REAL(WP):: R(20)   ! R, the radius of curvature of base curve at node points,
                     !  must have L values
  REAL(WP):: RIEXP   ! r, the exponent of nose-radius term in MSDOT equation
  REAL(WP):: RNSI    ! Initial nose radius
  REAL(WP):: RO      ! rho, the material density
  REAL(WP):: RODP    ! rho-double-prime, the density of layer along y=0
  REAL(WP):: ROP     ! rho-prime, the density of layer along x=L
  REAL(WP):: RS(20)  ! Rcyl    Cylindrical radius from body axis of symmetry to
                     ! node points on the base curve, must have L values
  REAL(WP):: RSSMAX  ! Maximum expected value of RSS, used to get a reasonable 
                     ! scale for plots, not needed if IPLOT=0
  INTEGER:: S             ! Number of stations in y-direction
  REAL(WP):: STEBOL       ! sigma, the Stefan-Boltzmann constant for radiation
  REAL(WP):: T(10,20)     ! Initial temperature, must have S*L values
  REAL(WP):: TALPHA(10)   ! Temperature table for ALPHAT
  REAL(WP):: TAUO         ! Initial time
  REAL(WP):: TBTAB(10)    ! Temperature to which back surface is radiating,
                          ! a function of time (TTABTB)
  REAL(WP):: TDPRIME      ! Thickness of layer along y=0
  REAL(WP):: THETA(20)    ! theta, the angle (in degrees) less than or equal 
                          ! to 90° between RS and R, must have L values
!  REAL(WP):: TMIN  ! Minimum temperature value; if TMIN/=0 and a computed 
!                   ! temperature goes below TMIN, the temperature will be set
!                   ! equal to TMIN; if TMIN=0, no restraint will be made on 
!                   ! the computed temperatures
  REAL(WP):: TPRIME       ! Thickness of layer along x=L
  REAL(WP):: TTABCP(10)   ! Temperature table for CPTAB
  REAL(WP):: TTABHC(7)    ! Temperature table for HCTAB
  REAL(WP):: TTABHE(10)   ! Time table for HETAB
  REAL(WP):: TTABHW(15)   ! Temperature table for HWTAB
  REAL(WP):: TTABQC (10)  ! Time table for QCTAB
  REAL(WP):: TTABQR(10)   ! Time table for QRTAB
  REAL(WP):: TTABTB(10)   ! Time table for TBTAB
  REAL(WP):: TTALC(10)    ! Time table for ALCTAB
  REAL(WP):: TTALS(10)    ! Time table for ALSTAB
  REAL(WP):: TTCKETA(10)  ! Temperature table for CKETATB
  REAL(WP):: TTCKXI(10)   ! Temperature table for CKXITAB
  REAL(WP):: TTELAM(7)    ! Temperature table for ELAMTB
  REAL(WP):: TTHCOMB(7)   ! Temperature table for HCOMBTB
  REAL(WP):: TTPSTAG(10)  ! Time table for PSTAGTB
  REAL(WP):: XITAB(5)     ! xi, the table of values of CKXITAB
  REAL(WP):: XO     ! x-sub-b, the length of base curve
  REAL(WP):: XORDER ! Order of oxidation
  REAL(WP):: ZS(20) ! Initial distance from the initial stagnation point to RSS
                    ! along body axis of symmetry, must have L values
  REAL(WP):: ZSMAX  ! Maximum expected value of ZS, used to get 
                    ! reasonable scale for plotting RSS versus ZS, 
                    ! not needed if IPLOT=O




CONTAINS

!+
SUBROUTINE InitializeInputQuantities()
! ------------------------------------------------------------------------------
! PURPOSE - Give each input quantity a default value

  REAL(WP),PARAMETER:: ZERO = 0
!-------------------------------------------------------------------------------


   aexp = ZERO
   alctab(:) = ZERO
   alphat(:) = ZERO
   alstab(:) = ZERO
   asexp = ZERO !     coefficient in the expression for calculating msdot.
   beta = ZERO    !     determines whether ablation or transpiration theory
    ! will be used for effect of mass transfer on heat transfer; 
    ! for ablation theory, beta=1; for transpiration theory,, beta=0
   bexp = ZERO   !      power of the exponential term in the arrhenius
            !  expression for calculating mcdot
   bsexp = ZERO  !      power of the exponential term in the expression
         !   for calculating msdot
   ce = ZERO      !     oxygen concentration, by mass, at edge of boundary layer
   cketatb(:) = ZERO   ! thermal conductivity in n-direction, a function of
        ! eta (etatab) and temperature (ttcketa) 
   ckxitab(:) = ZERO    ! thermal conductivity in -direction, a function of
       ! xi (xitab) and temperature (ttcioci)
   cordsy = ZERO   ! trigger to indicate coordinate system; 
                     ! if curvilinear coordinates, cordsy=0; 
                     ! if cartesian coordinates, cordsy=1
   cpdp = ZERO    ! c-sub-p-double prime, specific heat of layer along y=0
   cpp = ZERO     ! c-sub-p-prime     specific heat of layer along x=l
   cptab(:) = ZERO  ! c-sub-p  specific heat, a function of temperature (ttabcp)
   deltao(:) = ZERO  ! initial material thickness, must have l values
   deltau = ZERO   ! initial computing time interval
   deltmin = ZERO   ! minimum value allowed for delta
   dtmax = ZERO    ! maximum deltau which can be used; if no value
       ! is given, dtmax=2.0
   elamtb(:) = ZERO  ! lambda, the ratio of mass of material removed per unit 
                        ! mass of oxygen that reaches the surface, a function of 
                        ! pressure (pelam) and temperature (ttelam)
   endtau = ZERO    ! time at which calculation stops
   epsone = ZERO    ! epsilon, the emittance of front surface
   epsonep = ZERO   ! epsilon-prime, the emittance of layer along x=l
   epsonpp = ZERO   ! epsilon-double-prime, the emittance of layer along y=0
   errort = ZERO    ! acceptable relative error in temperature
   etatab(5) = ZERO  ! eta table for cketatb
   gambar = ZERO   ! mean ratio of specific heats behind bow shock wave, used
                 ! only in computation of heating-rate distribution around body
   gaminf = ZERO   ! ratio of specific heats in free stream, used only in 
                     ! computing heating-rate distribution around body
   hcombtb(:) = ZERO  ! delta-h-sub-c, the heat of combustion, a function of
                         ! pressure (phcomb) and temperature (tthcomb)
   hctab(:) = ZERO  ! delta-h-sub-s, the heat of sublimation, a function of 
                       ! pressure (phc) and temperature (ttabhc)
   hetab(:) = ZERO   ! h-sub-e, the total free-stream enthalpy, a function 
                         ! of time(ttabhe)
   hwtab(:) = ZERO   ! h-sub-w, the enthalpy of gas at the wall temperature,
                        ! a function of temperature (ttabhw)

   iadjust = 0 ! trigger for adjusting heating-rate and pressure 
                    ! distributions to shape change; if iadjust=o, 
                    ! qrat and prat are not adjusted; if iadjust/=0, 
                    ! qrat and prat will be adjusted to shape change
   iplot = 0    ! trigger for plotting routine; if iplot=o, no plots;
                    ! if iplot /= 0, the following plots will be made: 
                    ! 1) rss versus zs at times indicated in pltime table; 
                    ! 2) mdot versus x at each prfreq time; and
                    ! 3) t(m,n) versus x indicated in ntp array at each prereq
   l = 0        ! number of stations in the x-direction
   machno = ZERO   ! free-stream mach number
   malpha = 1   ! order of interpolation for alphat
   malphc = 1   ! order of interpolation for alctab
   malphs = 1   ! order of interpolation for alstab
   maxitt = 0   ! maximum iteration count; when iteration count exceeds 
                    ! this number, deltau will be halved until deltau is less 
                    ! than 1.0e-6, then the program will stop and a message will
                    ! be printed
   mcp = 1      ! order of interpolation for cptab
   mdmax = ZERO    ! maximum expected mdot; this must be given to get a 
                    ! reasonable scale for plots; not needed if iplot=0
   mdoto(20) = ZERO  ! initial mass loss rate at surface, must have l values
   mhe = 1     ! order of interpolation for hetab
   mhw = 1     ! order of interpolation for hwtab
   mpstag = 1  ! order of interpolation for pstagtb
   mqc = 1     ! order of interpolation for qctab
   mqr = 1     ! order of interpolation for qrtab
   mtb = 1     ! order of interpolation for tbtab
   mwo2 = ZERO    ! molecular weight of diatomic oxygen used in oxidation equation
   mwstr = ZERO   ! molecular weight of free stream used in oxidation equation
   nalpha = 0  ! number of entries in alphat
   nalphc = 0  ! number of entries in alctab
   nalphs = 0  ! number of entries in alstab
   ncketa = 0  ! number of entries in cketatb
   nckxi = 0   ! number of entries in ckxitab
   ncp = 0     ! number of entries in cptab
   nelam = 0   ! number of entries in elamtb
   neta = 0    ! number of entries in etatab
   nhc = 0     ! number of entries in hctab
   nhcomb = 0  ! number of entries in hcombtb
   nhe = 0     ! number of entries in hetab
   nhw = 0     ! number of entries in hwtab
   npelam = 0  ! number of entries in pelam
   nphc = 0    ! number of entries in phc
   nphcomb = 0 ! number of entries in phcomb
   npstag = 0  ! number of entries in pstagtb
   nqc = 0     ! number of entries in qctab
   nqr = 0     ! number of entries in qrtab
   ntb = 0     ! number of entries in tbtab
   ntp(7) = 0  ! array of seven values which specify the temperatures to 
                   ! be plotted; ntp(1) = the number of temperature rows to be 
                   ! plotted (may be six or less); ntp(2) through ntp(7), 
                   ! the row number of the temperatures to be plotted. 
                   ! for example, ntp(1)=3, ntp(2)=1, ntp(4)=10, specifies that
                   ! three (3) rows of temperature will be plotted and these 
                   ! rows are 1, 5, and 10
   nxi = 0     ! number of entries in xitab

   pelam(:) = ZERO    !     pressure table for elamtb
   phc(:) = ZERO     ! pressure table for hctab
   phcomb(:) = ZERO   !  pressure table for hcombtb
   pltime(:) = ZERO   ! times at which rss versus zs, that is, the body
                         ! shape, will be plotted; not needed if iplot=0
   prat(:) = ZERO   ! initial ratio of local to stagnation pressure, must
                       ! have l values, not needed if iadjuswo
   prfreq = ZERO   ! printing time frequency for output data
   psexp = ZERO   ! p, the exponent of pressure term in sublimation equation
   pstagtb(:) = ZERO  ! stagnation pressure, a function of time (ttpstag)
   ptmax = ZERO    ! maximum expected value of t, used to get reasonable scale 
    ! in plotting, not needed if iplot=0
   ptmin = ZERO    ! minimum expected value of t, used to get reasonable scale 
    ! in plotting; not needed if iplot=0
   qctab(:) = ZERO  ! q-sub-c, the cold-wall convective heating rate, 
                       ! a function of time (ttabqc)
   qrat(:) = ZERO   ! initial convective heating-rate distribution. must
                       ! have l values, not needed if iadjustl0
   qrrat(:) = ZERO  ! radiant heating-rate distribution over body, must have l values
   qrtab(:) = ZERO  ! q-sub-r, the radiant heating-rate tables, 
                       ! a function of time (ttabqr)
   r(:) = ZERO    ! r, the radius of curvature of base curve at node points,
                     !  must have l values
   riexp = ZERO    ! r, the exponent of nose-radius term in msdot equation
   rnsi = ZERO     ! initial nose radius
   ro = ZERO       ! rho, the material density
   rodp = ZERO     ! rho-double-prime, the density of layer along y=0
   rop = ZERO      ! rho-prime, the density of layer along x=l
   rs(:) = ZERO   ! rcyl    cylindrical radius from body axis of symmetry to
                     ! node points on the base curve, must have l values
   rssmax = ZERO   ! maximum expected value of rss, used to get a reasonable 
                     ! scale for plots, not needed if iplot=0
   s = ZERO              ! number of stations in y-direction
   stebol = 5.67E-8 ! sigma, the stefan-boltzmann constant for radiation
                    ! 5.67E-8 in SI units 
   t(:,:) = ZERO      ! initial temperature, must have s*l values
   talpha(:) = ZERO    ! temperature table for alphat
   tauo = ZERO          ! initial time
   tbtab(:) = ZERO     ! temperature to which back surface is radiating,
                          ! a function of time (ttabtb)
   tdprime = ZERO       ! thickness of layer along y=0
   theta(:) = ZERO     ! theta, the angle (in degrees) less than or equal 
                          ! to 90° between rs and r, must have l values
!!!   tmin = ZERO   ! minimum temperature value; if tmin/=0 and a computed 
                   ! temperature goes below tmin, the temperature will be set
                   ! equal to tmin; if tmin=0, no restraint will be made on 
                   ! the computed temperatures
   tprime = ZERO        ! thickness of layer along x=l
   ttabcp(:) = ZERO    ! temperature table for cptab
   ttabhc(:) = ZERO     ! temperature table for hctab
   ttabhe(:) = ZERO    ! time table for hetab
   ttabhw(:) = ZERO    ! temperature table for hwtab
   ttabqc (:) = ZERO   ! time table for qctab
   ttabqr(:) = ZERO    ! time table for qrtab
   ttabtb(:) = ZERO    ! time table for tbtab
   ttalc(:) = ZERO     ! time table for alctab
   ttals(:) = ZERO     ! time table for alstab
   ttcketa(:) = ZERO   ! temperature table for cketatb
   ttckxi(:) = ZERO    ! temperature table for ckxitab
   ttelam(:) = ZERO     ! temperature table for elamtb
   tthcomb(:) = ZERO    ! temperature table for hcombtb
   ttpstag(:) = ZERO   ! time table for pstagtb
   xitab(:) = ZERO      ! xi, the table of values of ckxitab
   xo = ZERO      ! x-sub-b, the length of base curve
   xorder = ZERO  ! order of oxidation
   zs(:) = ZERO  ! initial distance from the initial stagnation point to rss
                    ! along body axis of symmetry, must have l values
   zsmax = ZERO   ! maximum expected value of zs, used to get 
                    ! reasonable scale for plotting rss versus zs, 
                    ! not needed if iplot=o

  RETURN
END Subroutine InitializeInputQuantities   ! -----------------------------------

END Module InputsCommon   ! ====================================================

!+
MODULE HoldCommon
! ------------------------------------------------------------------------------
! PURPOSE -
USE AbaxiConstants,ONLY: WP

  REAL(WP):: tmin = 0.0_WP
!!!      COMMON /HOLD/ TMIN 
END Module HoldCommon   ! ======================================================


MODULE AbAxiProcedures
USE AbaxiConstants
IMPLICIT NONE
SAVE



CONTAINS

!+
SUBROUTINE Column()
! ------------------------------------------------------------------------------
! PURPOSE - Solve the matrix column by column for one iteration.!!!
!  Solves  m (no. of rows)  sets of simultaneous equations  n (no. of columns)
!  times, then  returns to main program to test for convergence.
USE PickCommon
USE InputsCommon
IMPLICIT NONE

!!!      REAL MDOTO,MDOT,MCDOT,MSDOT,MWSTR,MWO2,MACHNO,MDMAX 
!!!      INTEGER S,SM1,SM2 
  INTEGER:: n,n1,n2 
!-------------------------------------------------------------------------------
! COMPUTE  COLUMN 1                                                     
      N1 =2 
      N2 =SM1 
      CALL COLXO (N1,N2) 
      CALL SOLMAT (A(1:,1),B,C(1:,1),Z(1),V(1),DC,TTF(1:,1),S) 

! COMPUTE  COLUMN 2 THRU L-1  (LM1)                                            
      DO 300 N=2,LM1 
      CALL COLMN (N1,N2,N) 
      CALL SOLMAT (A(1:,N),B,C(1:,N),Z(N),V(N),DC,TTF(1:,N),S) 
  300 END DO 

! Compute  column  L                                                    
  CALL Colxl(n1,n2) 
  CALL Solmat (a(1:,l),b,c(1:,l),z(l),v(l),dc,ttf(1:,l),s) 
  RETURN 
END Subroutine Column   ! ------------------------------------------------------

!+
SUBROUTINE Row()
! ------------------------------------------------------------------------------
! PURPOSE - Solves the matrix row by row for one iteration.
!  Solves  n (no. of columns) sets of simultaneous eqs. m(no.of rows) times
!    then returns to main program to check for convergence.
USE PickCommon
USE InputsCommon
IMPLICIT NONE

!!!      REAL MDOTO,MDOT,MCDOT,MSDOT,MWSTR,MWO2,MACHNO,MDMAX 
!!!      INTEGER SM1 ,S 
!!!      DIMENSION ANS(20), ATEMP(20), CTEMP(20) 


  REAL(WP),DIMENSION(20):: ans, atemp, ctemp
  INTEGER:: m,n
  INTEGER:: n1,n2

!-------------------------------------------------------------------------------
! COMPUTE  ROW 1                                                        
      N1 =2 
      N2 =LM1 
      CALL COLXO (N1,N2) 
      DO  300 N=2,LM1 
      CALL COLMN (N1,N2,N) 
  300 END DO 
      CALL COLXL(N1,N2) 
      DO 320 N =1,L 
        ATEMP(N) = AB(1,N) 
        CTEMP(N) = CB(1,N) 
  320 END DO
      CALL SOLMAT (ATEMP,B,CTEMP,ZB(1),VB(1),DC,ANS(1:),L) 
      DO 400 N=1,L 
        TTF(1,N)=ANS(N) 
  400 END DO

! COMPUTE  ROW  2  THRU SM1                                             
      DO  600  M=2,SM1 
      N1 =M 
      N2 =M 
      CALL COLXOM (N1,N2) 
      DO 500 N=2,LM1 
      CALL COLMNMN(N1,N2,N) 
  500 END DO 
      CALL COLXLM (N1,N2) 
      DO 510 N=1,L 
      ATEMP(N) = AB(M,N) 
  510 CTEMP(N) = CB(M,N) 
      CALL SOLMAT (ATEMP,B,CTEMP,ZB(M),VB(M),DC,ANS(1:),L) 
      DO 590 N=1,L 
  590 TTF(M,N)=ANS(N) 
  600 END DO 

! COMPUTE ROW S                                                         
      CALL  COLXO1(N1,N2) 
      DO  800 N=2,LM1 
      CALL COLMN1 (N1,N2,N) 
  800 END DO 
      CALL COLXL1(N1,N2) 
      DO 810 N=1,L 
      ATEMP(N) = AB(S,N) 
      CTEMP(N) = CB(S,N) 
  810 END DO
      CALL SOLMAT  (ATEMP,B,CTEMP,ZB(S),VB(S),DC,ANS(1:),L) 
      DO 890 N=1,L 
        TTF(S,N)=ANS(N) 
  890 CONTINUE
  RETURN 
END Subroutine Row   ! ---------------------------------------------------------

!+
SUBROUTINE COLXO(N1,N2) 
! ------------------------------------------------------------------------------
! PURPOSE - COMPUTE COEF. FOR XI=0,   COLUMN IMPLICIT                             
! IROCOL = 1       COLUMN IMPLICIT                                      
! IROCOL = 2       ROW  IMPLICIT                                        
USE PickCommon
USE InputsCommon
IMPLICIT NONE
  INTEGER,INTENT(IN):: n1,n2

!!!      REAL MDOTO,MDOT,MCDOT,MSDOT,MWSTR,MWO2,MACHNO,MDMAX 
!!!      INTEGER S,SM1,SM2 

  REAL(WP):: bsave
  REAL(WP):: cord
  REAL(WP):: dd
  REAL(WP):: ddqs
  REAL(WP):: delde
  REAL(WP):: dr
  REAL(WP):: ept4, eptb
  REAL(WP):: ff
  REAL(WP):: fp
  REAL(WP):: g
  REAL(WP):: h
  REAL(WP):: h1r

  INTEGER:: i
  INTEGER:: m, mm1, mp1

  REAL(WP):: p
  REAL(WP),SAVE:: p817   ! SAVE added by RLC 6 Sept 2008
  REAL(WP):: part1,part2
  REAL(WP):: sc
  REAL(WP):: u
  REAL(WP):: xx
!-------------------------------------------------------------------------------
! STATION (1,1)    XI=0  ,  ETA=0                                       
!                                                                       
      DO 60 I=1,SM1 
        CK(I)= (CKETA(I,1)+ CKETA(I+1,1))/2.0 
   60 END DO
      DELDE = DELTA(1)* DELETA 
      PART2= H1(1,1) **2  * XDXISQ 
      PART1=RODPC 
      H1R = H1(1,1) * R(1) 
      FF=CKXI(1,1)*(2.0-CORDSY)/(2.0*PART2) 
      G=RO*CP(1,1)/DELTAU-2.0*PART1/H1R+8.0*PART1/(3.0*DELDE) 
      H = 1.0/( H2(1,1)**2 * DELTA(1)**2) 
      SC= H /(3.0* DELESQ) 
      EPT4=SIGDP* (2.0/(H1R*H2(1,1)**2) - EIGHT3/DELDE) 
      EPTB= EPT4 *TB 
      EPT4= EPT4 *T(1,1)**3 
      BSAVE  = G 
      GO TO  (70, 80), IROCOL 
   70 CONTINUE 
      A(1,1) = 0.0 
      BS1(1,1) = -SC*9.0 *CK(1) 
      C(1,1)=  SC * (9.0 *CK(1) + CK(2) ) 
      Z(1) = -SC * CK(2) 
      B(1)= BS1(1,1) - BSAVE + EPT4 
      IF (IFIRST.EQ.0 )  GO TO 80 
   78 DC(1)=(-BSAVE-BS1B(1,1))*T(1,1) - CB(1,1)*T(1,2)- ZB(1)*T(1,3) + EPTB
      GO TO 99 
   80 FP=FF 
      BS1B(1,1)= -7.0* FP 
      CB(1,1)= 8.0 *FP 
      ZB(1) = -FP 
      IF (IFIRST.EQ.0 )  GO TO 78 
   86 B(1) = BS1B(1,1)- BSAVE + EPT4 
      DC(1)=(-BSAVE -BS1(1,1))*T(1,1) -C(1,1)*T(2,1) - Z(1)*T(3,1) + EPTB
   99 GO TO (101,600),IROCOL 
!                                                                       
! STATION(M,1)   , XI=0   ,  ETA LESS THAN 1 , GREATER THAN 0           
!                                                                       
      ENTRY COLXOM(n1,n2)
  101 DO 200 M=N1,N2 
      DELDE=DELTA(1)*DELETA 
      MP1=M+1 
      MM1=M-1 
      P817=8.0*DELTA(2) - DELTA(3) - 7.0*DELTA(1) 
      PART2 = H1(M,1)**2  * XDXISQ 
      CORD= (2.0-CORDSY)/2.0 
      FF= CKXI(M,1)*CORD/PART2 
      G = RO *CP(M,1)/DELTAU 
      SC = 1.0 /(H2(M,1)* DELDE **2) 
      H = FF* P817/(2.0* DELDE)  *ETA(M) 
      P =  CKETA(M,1)/(H2(M,1)**2 *H1(M,1)* R(1) * DELDE) 
      BSAVE =G 
      GO TO (170,180), IROCOL 
  170 CONTINUE 
      U= ETA(M)*MDOT(1) * CP(M,1)/(2.0*DELTA(1) * DELETA) 
      A(M,1)=  H -P + SC* CK(MM1) +U 
      BS1(M,1) = SC * (-CK(MM1) - CK(M)) 
      C(M,1)= -H + P + SC* CK(M) -U 
      B(M) = BS1(M,1) - BSAVE 
      IF (IFIRST.EQ.0 )  GO TO 180 
  178 DC(M)=(-BSAVE -BS1B(M,1))* T(M,1)-ZB(M)*T(M,3)-CB(M,1)*T(M,2) 
      CYCLE   !!! GO TO 200 
  180 ZB(M) =  -FF 
      CB(M,1)= 8.0 * FF 
      BS1B(M,1)= -7.0*FF 
      IF (IFIRST.EQ.0 )  GO TO 178 
  190 B(1) = BS1B(M,1) - BSAVE 
      DC(1)= (-BSAVE - BS1(M,1))*T(M,1)-A(M,1)*T(MM1,1)-C(M,1)*T(MP1,1) 
  200 END DO 
      GO TO (202,600),IROCOL 
!                                                                       
! STATION (S,1)   ,XI=0   ,  ETA=1                                      
!                                                                       
      ENTRY COLXO1(n1,n2)
  202 CORD=(2.0-CORDSY)/2.0 
      FF=CKXI(S,1)*CORD/H1(S,1)**2 
      DELDE=DELTA(1)*DELETA 
      P =FF/XDXISQ 
      H = 1.0/(H2(S,1)**2 *3.0* DELDE**2) 
      G = RO*  CP(S,1)/ DELTAU 
      SC = -9.0 * CK(SM1) * H 
      BSAVE = G 
      GO TO (270,280) ,IROCOL 
  270 CONTINUE 
      XX=CP(S,1)*MDOT(1)/(2.0*DELTA(1)*DELETA) 
      V(1)= -CK(SM2)*H - XX 
      A(S,1) = -SC + CK(SM2)*H + 4.0*XX 
      DR=P*P817/CKETA(S,1)*H2(S,1) 
      DD = DR - 2.0/(H1(S,1)*R(1)*H2(S,1))-EIGHT3/(H2(S,1)*DELDE)                                                  
      DDQS=DD*QS(1) 
      BS1(S,1)=DD*SIG*T(S,1)**3 + SC - 3.0*XX 
      B(S) = BS1(S,1) -BSAVE 
      IF (IFIRST.EQ.0 )  GO TO 280 

  278 CONTINUE
      DR=P*P817/CKETA(S,1)*H2(S,1) ! copied here by RLC  6 Sept 2008
      DD = DR - 2.0/(H1(S,1)*R(1)*H2(S,1))-EIGHT3/(H2(S,1)*DELDE) ! copied by RLC
      DDQS=DD*QS(1)   ! copied here by RLC
      DC(S) = DDQS + (-BSAVE - BS1B(S,1))*T(S,1)- CB(S,1)*T(S,2) - ZB(S)*T(S,3)
      GO TO 600 

  280 CB(S,1)=8.0*P 
      ZB(S) = -P 
      BS1B(S,1) = -7.0*P 
      IF (IFIRST.EQ.0 )  GO TO 278 
  290 B(1) = BS1B(S,1) - BSAVE 
      DR=P*P817/CKETA(S,1)*H2(S,1)   ! copied here by RLC
      DD = DR - 2.0/(H1(S,1)*R(1)*H2(S,1))-EIGHT3/(H2(S,1)*DELDE) ! copied by RLC
      DDQS=DD*QS(1)  ! copied here by RLC
      DC(1)=(-BSAVE-BS1(S,1))*T(S,1) - V(1)*T(SM2,1) - A(S,1)*T(SM1,1)  +DDQS                              
  600 RETURN 
  201 FORMAT    (7E18.7) 
END Subroutine Colxo   ! -------------------------------------------------------

!+
SUBROUTINE Colmn(N1,N2,N) 
! PURPOSE - 
! IROCOL = 1       COLUMN IMPLICIT                                      
! IROCOL = 2       ROW  IMPLICIT                                        
USE PickCommon
USE InputsCommon
IMPLICIT NONE
  INTEGER,INTENT(IN):: n1,n2,n   ! ??

!!!      REAL MDOTO,MDOT,MCDOT,MSDOT,MWSTR,MWO2,MACHNO,MDMAX 
!!!      DIMENSION DDQS(20),DDQSR(20) 
!!!      INTEGER S,SM1,SM2 

  REAL(WP):: bsave
  REAL(WP):: d1,d2
  REAL(WP):: d1nm1, d1np1
  REAL(WP):: dd
  REAL(WP),DIMENSION(20),SAVE:: ddqs, ddqsr   ! SAVE added by RLC 6 Sept 2008
  REAL(WP):: denom
  REAL(WP):: dmm12n, dmp12n
  REAL(WP):: dmnm12, dmnp12
  REAL(WP):: dsm12n, dsm32n
  REAL(WP):: dsnm12, dsnp12
  REAL(WP):: ds
  REAL(WP):: e32n, e52n
  REAL(WP):: emm12n, emp12n
  REAL(WP):: esm12n, esm32n
  REAL(WP):: ept4, eptb
  REAL(WP):: fmm12n, fmp12n
  REAL(WP):: fmnm12, fmnp12
  REAL(WP):: fsm12n, fsm32n
  REAL(WP):: fsnm12, fsnp12
  REAL(WP):: g
  REAL(WP):: g1n
  REAL(WP):: h1h3
  REAL(WP):: hmn
  INTEGER:: m, mm1, mp1
  INTEGER:: nm1, np1
  REAL(WP):: p1nm1, p1np1
  REAL(WP):: part
  REAL(WP):: qsn, qsnm1, qsnp1
  REAL(WP):: sst
  REAL(WP):: u
  REAL(WP):: vv  
  REAL(WP):: w
  REAL(WP):: xx
  REAL(WP):: yy
!-------------------------------------------------------------------------------
! STATION (1,N)    XI GREATER THAN 0,LESS THAN 1     ETA=0 

!!!  201 FORMAT    (7E18.7) 
      NM1 = N-1 
      NP1 = N+1 
      E32N=(H1(2,N)+H1(1,N))*(H3(2,N)+H3(1,N))*   &
        (CKETA(2,N)+CKETA(1,N))/(4.*(H2(2,N)+H2(1,N)))
      E52N=(H1(3,N)+H1(2,N))*(H3(3,N)+H3(2,N))*   &
        (CKETA(3,N)+CKETA(2,N))/(4.*(H2(3,N)+H2(2,N)))
      VV=1.0/ (3.0* DELTA(N)**2  * DELESQ ) 
      P1NP1=(H3(1,NP1)+H3(1,N))/(H1(1,NP1)+H1(1,N)) !!! never used???
      P1NM1=(H3(1,NM1)+H3(1,N))/(H1(1,NM1)+H1(1,N)) !!! never used???
      W=  H1(1,N)*H3(1,N)* DELETA *DELTA(N) *8.0 
      G1N = H1(1,N)* H2(1,N) * H3(1,N) * RO *CP(1,N) 
      YY=(-VV*W*RODPC-G1N/DELTAU) 
      EPT4= -VV *W *SIGDP 
      EPTB=  EPT4 * TB 
      EPT4 = EPT4 * T(1,N)**3 
      BSAVE = YY 
      GO TO (170,180),IROCOL 
  170 CONTINUE 
      BS1(1,N) = -VV* 9.0* E32N 
      C(1,N)= VV *(9.0* E32N + E52N) 
      Z(N) = -VV * E52N 
      B(1)= BS1(1,N) + BSAVE +  EPT4 
      IF (IFIRST.EQ.0 )  GO TO 180 
  178 DC(1) = (BSAVE - BS1B(1,N))*T(1,N) -AB(1,N)*T(1,NM1)-CB(1,N)*     &
     & T(1,NP1) + EPTB                                                  
      GO TO 200 
  180 D1NP1=(H2(1,NP1)+H2(1,N))*(H3(1,NP1)+H3(1,N))*(CKXI(1,NP1)+CKXI(1,&
     &N))/(4.*XDXISQ*(H1(1,NP1)+H1(1,N)))                               
      D1NM1=(H2(1,NM1)+H2(1,N))*(H3(1,NM1)+H3(1,N))*(CKXI(1,NM1)+CKXI(1,&
     &N))/(4.*XDXISQ*(H1(1,NM1)+H1(1,N)))                               
      AB(1,N)=D1NM1 
      BS1B(1,N)=-D1NP1- D1NM1 
      CB(1,N)=D1NP1 
      IF (IFIRST.EQ.0 )  GO TO 178 
  190 B(N)= BS1B(1,N) + BSAVE  + EPT4 
      DC(N) = (BSAVE -BS1(1,N))*T(1,N) -C(1,N)*T(2,N) - Z(N)*T(3,N)     &
     & + EPTB                                                           
  200 CONTINUE 
      GO TO (202,800), IROCOL 
!                                                                       
! STATION (M,N)   XI GREATER THAN 0, LESS THAN 1                        
!                ETA GREATER THAN 0, LESS THAN 1                        
!                                                                       
      ENTRY COLMNMN(n1,n2,n) 
      NP1=N+1 
      NM1=N-1 
  202 DO 400 M=N1,N2 
      MM1 = M-1 
      MP1 = M+1 
      VV= 1.0 /(DELTA(N)**2 * DELESQ) 
      XX  =  ETA(M)*AA(N)/(DELTA(N)*  DELESQ) 
      G   =  H1(M,N)* H2(M,N) * H3(M,N) *RO * CP(M,N) 
      EMM12N=(H1(MM1,N)+H1(M,N))*(H3(MM1,N)+H3(M,N))*(CKETA(MM1,N)+CKETA&
     &(M,N))/(4.*(H2(MM1,N)+H2(M,N)))*VV                                
      EMP12N=(H1(MP1,N)+H1(M,N))*(H3(MP1,N)+H3(M,N))*(CKETA(MP1,N)+CKETA&
     &(M,N))/(4.*(H2(MP1,N)+H2(M,N)))*VV                                
      DMM12N=(H2(MM1,N)+H2(M,N))*(H3(MM1,N)+H3(M,N))*(CKXI(MM1,N)+CKXI(M&
     &M1,N))/(4.*(H1(MM1,N)+H1(M,N)))                                   
      DMP12N=(H2(MP1,N)+H2(M,N))*(H3(MP1,N)+H3(M,N))*(CKXI(MP1,N)+CKXI(M&
     &P1,N))/(4.*(H1(MP1,N)+H1(M,N)))                                   
      FMM12N=XX*DMM12N*AA(N)*(ETA(MM1)+ETA(M))/(DELTA(N)*2.) 
      FMP12N=XX*DMP12N*AA(N)*(ETA(MP1)+ETA(M))/(DELTA(N)*2.) 
      W = 4.0 * XO * DELXI * DELETA 
      DENOM= 4.0* (DELTA(  NM1) + DELTA (  N)) *( H1(M,NM1) + H1 (M,N)) 
      FMNM12= (H3(M,NM1)+H3(M,N)) *(H2(M,NM1)+H2(M,N))* (CKXI(M,NM1)    &
     &+ CKXI(M,N))* (AA(  NM1) + AA(  N))* ETA(M)/DENOM                 
      DENOM= 4.0* (DELTA(  NP1)+DELTA(  N))*(H1(M,NP1)+H1(M,N)) 
      FMNP12= (H3(M,NP1)+H3(M,N))*(H2(M,NP1) +H2(M,N))*(CKXI(M,NP1)     &
     & +CKXI(M,N))*(AA(  NP1)+AA(  N))*ETA(M)/DENOM                     
      D1 = (FMNP12*(T(MP1,NP1)-T(MM1,NP1)+T(MP1,N)-T(MM1,N))-FMNM12*    &
     & (T(MP1,N)-T(MM1,N)+T(MP1,NM1)-T(MM1,NM1)))/W                     
      D2 = ETA(M) *AA(N)* CKXI(M,N)* (H2(MP1,N)* H3(MP1,N)* (T(MP1,NP1) &
     & - T(MP1,NM1))/H1(MP1,N)   - H2(MM1,N) * H3(MM1,N)* (T(MM1,NP1)   &
     & - T(MM1,NM1))/H1(MM1,N) ) /(DELTA(N) * W)                        
      DS = D1  + D2  - G *T(M,N )/ DELTAU 
      BSAVE = G/DELTAU 
      GO TO  (370,380),IROCOL 
  370 CONTINUE 
      HMN = ETA(M) * MDOT(N)/(DELTA(N) * RO) 
      YY=  G * HMN/(2.0* DELETA) 
      A(M,N) = EMM12N + FMM12N + YY 
      BS1(M,N) = -EMM12N - EMP12N -FMP12N - FMM12N 
      C(M,N) = EMP12N + FMP12N - YY 
      B(M) = BS1(M,N) - BSAVE 
      IF (IFIRST.EQ.0 )  GO TO 380 
  378 DC(M) = DS- BS1B(M,N)* T (M,N) -AB(M,N)*T(M,NM1)-CB(M,N)*T(M,NP1) 
      CYCLE   !!! GO TO 400 
  380 DMNM12=(H2(M,NM1)+H2(M,N))*(H3(M,NM1)+H3(M,N))*   &
        (CKXI(M,NM1)+CKXI(M,N))/(4.*(H1(M,NM1)+H1(M,N)))
      AB(M,N)=DMNM12/XDXISQ 
      DMNP12=(H2(M,NP1)+H2(M,N))*(H3(M,NP1)+H3(M,N))*   &
        (CKXI(M,NP1)+CKXI(M,N))/(4.*(H1(M,NP1)+H1(M,N)))
      CB(M,N)=DMNP12/XDXISQ 
      BS1B(M,N)= -AB(M,N) - CB(M,N) 
      IF (IFIRST.EQ.0 )  GO TO 378 
  390 B(N) = BS1B(M,N) - BSAVE 
      DC(N) = DS - BS1(M,N)*T(M,N) - A(M,N)*T(MM1,N)- C(M,N)*T(MP1,N) 
  400 END DO 
      GO TO  (401,800), IROCOL 
!                                                                       
! STATION (S,N)   XI GREATER THAN 0, LESS THAN 1 ,     ETA =1           
!                                                                       
      ENTRY  COLMN1(n1,n2,n) 
      NP1=N+1 
      NM1=N-1 
  401 H1H3 = H1(S,N)* H3(S,N) 
      XX= 3.0 * DELTA(N)**2 * DELESQ 
      U = AA(N)/ (3.0                            *DELESQ* DELTA(N)   ) 
      G =            H1H3  *H2(S,N) * RO *CP(S,N) 
      PART=AA(N)/(DELTA(N)*4.0*DELETA*DELXI*XO) 
      SST=     H3(S,N)*CKXI(S,N)*3./H1(S,N) 
      DS=PART*(SST*T(S,NP1)-SST*T(S,NM1)                                &
     &    -4.0*H2(SM1,N)*H3(SM1,N)*CKXI(SM1,N)*(T(SM1,NP1)-T(SM1,NM1))/ &
     &H1(SM1,N)+H2(SM2,N)*H3(SM2,N)*CKXI(SM2,N)* (T(SM2,NP1)-T(SM2,NM1))&
     &/H1(SM2,N))                                                       
      ESM32N=(H1(SM1,N)+H1(SM2,N))*(H3(SM1,N)+H3(SM2,N))*(CKETA(SM1,N)+ &
     &CKETA(SM2,N))/(4.*(H2(SM1,N)+H2(SM2,N))*XX)                       
      ESM12N=(H1(SM1,N)+H1(S,N))*(H3(SM1,N)+H3(S,N))*(CKETA(SM1,N)+CKETA&
     &(S,N))/(4.*(H2(SM1,N)+H2(S,N))*XX)*9.                             
      DSM12N=(H2(SM1,N)+H2(S,N))*(H3(SM1,N)+H3(S,N))*(CKXI(SM1,N)+CKXI(S&
     & ,N))/ (4.0*(H1(SM1,N)+ H1(S,N)))                                 
      FSM12N=DSM12N*AA(N)*(ETA(SM1)+ETA(S))/(DELTA(N)*2.)*9.*U 
      DSM32N=(H2(SM2,N)+H2(SM1,N))*(H3(SM2,N)+H3(SM1,N))*(CKXI(SM2,N)+  &
     & CKXI(SM1,N))/(4.*(H1(SM2,N)+H1(SM1,N)))                          
      FSM32N=DSM32N*AA(N)*(ETA(SM2)+ETA(SM1))/(DELTA(N)*2.)*U 
      BSAVE = G/DELTAU 
      GO TO  (570,580),IROCOL 
  570 CONTINUE 
      YY=G*MDOT(N)/ (RO*2.0*DELTA(N)*DELETA) 
      V(N)= -ESM32N - FSM32N -YY 
      A(S,N) = ESM12N + ESM32N + FSM12N + FSM32N + 4.0*YY 
      DD=8.*H1H3*DELTA(N)*DELETA/XX  + 8.*U*                            &
     &H2(S,N)*DELTA(N)*F(S,N)*DELETA/CKETA(S,N)                         
      DDQS(N)=DD*QS(N) 
      BS1(S,N)=-DD*SIG*T(S,N)**3-ESM12N-FSM12N-3.0*YY 
      B(S) = BS1(S,N) - BSAVE 
      IF (IFIRST.EQ.0 )  GO TO 580 
  578 DC(S) =-DDQS(N)+ DS+(-BSAVE-BS1B(S,N))*T(S,N)-AB(S,N)*T(S,NM1)  &
     & -CB(S,N)*T(S,NP1)+ DDQSR(N)                                      
      GO TO 600 

  580 DSNM12=(H2(S,NM1)+H2(S,N))*(H3(S,NM1)+H3(S,N))*(CKXI(S,NM1)+CKXI(S&
     &,N))/(4.*(H1(S,NM1)+H1(S,N)) *XDXISQ)                             
      DSNP12=(H2(S,NP1)+H2(S,N))*(H3(S,NP1)+H3(S,N))*(CKXI(S,NP1)+CKXI(S&
     &,N))/(4.*(H1(S,NP1)+H1(S,N)) *XDXISQ)                             
      DENOM=4.0*(DELTA(  NM1)+DELTA(  N))*(H2(S,NM1)+H1(S,N)) 
      FSNM12= (H3(S,NM1)+H3(S,N))*(H2(S,NM1)+H2(S,N))* (CKXI(S,NM1)     &
     & +CKXI(S,N))*(AA(  NM1)+AA(  N))/DENOM                            
      DENOM=4.0*(DELTA(  NP1)+DELTA(  N))*(H1(S,NP1)+H1(S,N)) 
      FSNP12= (H3(S,NP1)+H3(S,N))*(H2(S,NP1)+H2(S,N))*(CKXI(S,NP1)      &
     & +CKXI(S,N))*(AA(  NP1)+AA(  N))/DENOM                            
      DENOM=2.0*XO*DELXI 
      QSN= DELTA(N)*H2(S,N)/(CKETA(S,N)*DENOM) 
      QSNP1= DELTA(NP1)* H2(S,NP1)/(CKETA(S,NP1)* DENOM) 
      QSNM1= DELTA(NM1)* H2(S,NM1)/(CKETA(S,NM1)* DENOM) 
      DDQSR(N)= FSNP12* (QSNP1*QS(NP1)+ QSN*QS(N))-FSNM12*              &
     &(QSN*QS(N)+ QSNM1*QS(NM1))                                        
      AB(S,N)=DSNM12-FSNM12*SIG*QSNM1*T(S,NM1)**3 
      CB(S,N)=DSNP12+FSNP12*QSNP1*SIG*T(S,NP1)**3 
      BS1B(S,N)=-DSNP12-DSNM12+SIG*T(S,N)**3*QSN *(FSNP12-FSNM12) 
      IF (IFIRST.EQ.0 )  GO TO 578 
  590 B(N)=BS1B(S,N)-BSAVE 
      DC(N) = -DDQS(N) +(-BSAVE-BS1(S,N))*T(S,N)+DS                    &
     &-A(S,N)*T(SM1,N)-V(N)*T(SM2,N)+ DDQSR(N)                          
  600 CONTINUE 
  800 RETURN
END Subroutine Colmn   ! -------------------------------------------------------

!+
SUBROUTINE COLXL(N1,N2)
! PURPOSE - computes coef. for  xi=1 ( x=l)   column implicit
! IROCOL = 1       COLUMN IMPLICIT                                      
! IROCOL = 2       ROW  IMPLICIT                                        
USE PickCommon
USE InputsCommon
IMPLICIT NONE
!!! SAVE

!!!      REAL MDOTO,MDOT,MCDOT,MSDOT,MWSTR,MWO2,MACHNO,MDMAX 
!!!      DIMENSION AL(10) 
!!!      INTEGER S,SM1,SM2 

  REAL(WP):: add,add1,add2
  REAL(WP),DIMENSION(10):: al
  REAL(WP):: aj,am,an
  REAL(WP):: bsave
  REAL(WP):: d1,d2
  REAL(WP):: d1lm12, d1lm32
  REAL(WP):: ddqsr
  REAL(WP):: ddsl
  REAL(WP):: dedeta
  REAL(WP):: dhk,dhk1,dhk2
  REAL(WP):: dmlm12,dmlm32
  REAL(WP):: dn,dn1
  REAL(WP):: dslm1,dslm3
  REAL(WP):: dsm12l,dsm32l
  REAL(WP):: e32l, e52l
  REAL(WP):: emm12, emp12
  REAL(WP):: ept4, eptb
  REAL(WP):: esm12l, esm32l
  REAL(WP):: ff
  REAL(WP):: g
  REAL(WP):: gh
  REAL(WP):: gsl
  REAL(WP):: h
  INTEGER:: m
  INTEGER:: mm1,mp1
  INTEGER:: n1,n2
  REAL(WP):: part, part1, part2
  REAL(WP):: partd1, parte1
  REAL(WP):: partd3, parte3
  REAL(WP):: partw
  REAL(WP):: qsave   ! SAVE added by RLC 
  REAL(WP):: qsl, qslm1, qslm2
  REAL(WP):: sc
  REAL(WP):: sp
  REAL(WP):: twdel
  REAL(WP):: u
  REAL(WP):: u1,u2
  REAL(WP):: w
  REAL(WP):: wsq
  REAL(WP):: wxodxi
  REAL(WP):: xx
  REAL(WP):: xx1
  REAL(WP):: xy
  REAL(WP):: yy
  REAL(WP):: zz2

  201 FORMAT    (7E18.7) 
!-------------------------------------------------------------------------------

! STATION  (1,L)     X= L  ,   ETA =0,
      W= 3.0* XO**2 * DELXI 
      U= 8.0*H2(1,L)*H3(1,L) *XO 
      XX = 8.0 * H1(1,L) * H3(1,L)* DELTA(L) 
      SC=  3.0*  DELTA(L)**2  * DELETA 
      G= -U*ROPCPP/W-XX*RODPC/SC - H1(1,L)*H2(1,L)*H3(1,L)*RO*CP(1,L)/DELTAU
      PART1 =   SC * DELETA 
      E32L=(H1(2,L)+H1(1,L))*(H3(2,L)+H3(1,L))*(CKETA(2,L)+CKETA(1,L))/ &
     &(4.*(H2(2,L)+H2(1,L)))*9.                                         
      E52L=(H1(3,L)+H1(2,L))*(H3(3,L)+H3(2,L))*(CKETA(3,L)+CKETA(2,L))/ &
     &(4.*(H2(3,L)+H2(2,L)))                                            
      D1LM32=(H2(1,LM1)+H2(1,LM2))*(H3(1,LM1)+H3(1,LM2))*(CKXI(1,LM1)+  &
     &CKXI(1,LM2))/(4.*(H1(1,LM1)+H1(1,LM2)))                           
      D1LM12=(H2(1,LM1)+H2(1,L))*(H3(1,LM1)+H3(1,L))*(CKXI(1,LM1)+CKXI(1&
     &,L))/(4.*(H1(1,LM1)+H1(1,L)))                                     
      EPT4=  (-U*SIGP/W  -XX *SIGDP/SC) 
      EPTB= EPT4* TB 
      EPT4= EPT4* T(1,L) **3 
      BSAVE = G 

      GO TO ( 150,180),IROCOL 
  150 CONTINUE 
      BS1(1,L)= -E32L/PART1 
      C(1,L)= (E52L + E32L)/PART1 
      Z(L)= -E52L/PART1 
      B(1)= BS1(1,L) + BSAVE  + EPT4 
      IF (IFIRST.EQ.0 )  GO TO 180 
  178 DC(1)=(BSAVE - BS1B(1,L))*T(1,L) - VB(1)*T(1,LM2) - AB(1,L)*T(1,LM1) + EPTB
      GO TO 198 

  180 CONTINUE 
      VB(1)=- D1LM32/(W*DELXI) 
      AB(1,L)= (D1LM32+ 9.0*D1LM12)/(W*DELXI) 
      BS1B(1,L)=-9.0*D1LM12/ (W*DELXI) 
      IF (IFIRST.EQ.0 )  GO TO 178 
  190 B(L) = BS1B(1,L) + BSAVE  + EPT4 
      DC(L)=(BSAVE -BS1(1,L))*T(1,L)- C(1,L)*T(2,L) - Z(L)*T(3,L) + EPTB
  198 CONTINUE 
      GO TO (202,800),IROCOL 

! STATION (M,L)    X=L      ETA GREATER THAN 0, LESS THAN 1 

      ENTRY  COLXLM(n1,n2) 
  202 DO 210 M=1,S 
  210 AL(M) =  H2(M,L)*H3(M,L)/H1(M,L) 
      W= 3.0 * XO * DELXI 
      YY = DELTA(L) **2  * DELESQ 
      DO 300 M=N1,N2 
      MM1 = M-1 
      MP1 = M+1 
      XX = ETA(M) *(AA(L)  + AA(LM1))/(4.0* (DELTA(L)+ DELTA(LM1))*DELETA) 
      XX1= ETA(M) *(AA(LM1)+ AA(LM2))/(4.0*(DELTA(LM1)+DELTA(LM2))*DELETA)
      XY = 8.0* D(M,L) * H1(M,L)/ CKXI(M,L) 
      AN = ETA(M) *AA(L)* CKXI(M,L)/ DELTA(L) 
      AM = AN/(DELTA(L) * DELESQ) 
      G = H1(M,L)* H2(M,L)* H3(M,L) * RO * CP(M,L) 
      AJ =  AN /(4.0* DELETA * XO  * DELXI ) 
      U1= (H2(MP1,L)+H2(M,L))* (H3(MP1,L)+H3(M,L)) *(ETA(MP1)+ETA(M))   &
     & /(4.0* (H1(MP1,L)+H1(M,L)))*AA(L)                                
      U2= (H2(MM1,L)+H2(M,L)) * (H3(MM1,L)+H3(M,L)) *(ETA(MM1)+ETA(M))  &
     &/ (4.0* (H1(MM1,L) +H1(M,L)))*AA(L)                               
      DMLM32=(H2(M,LM1)+H2(M,LM2))*(H3(M,LM1)+H3(M,LM2))*(CKXI(M,LM1)+  &
     &CKXI(M,LM2))/(4.*(H1(M,LM1)+H1(M,LM2)))                           
      DMLM12=(H2(M,LM1)+H2(M,L))*(H3(M,LM1)+H3(M,L))*   &
       (CKXI(M,LM1)+CKXI(M,L))/(4.*(H1(M,LM1)+H1(M,L)))                                     
      D1=-9.0*DMLM12*XX*(T(MP1,L)-T(MM1,L) + T(MP1,LM1) - T(MM1,LM1))
      D2=DMLM32*(-XX1)*(T(MP1,LM1)-T(MM1,LM1) + T(MP1,LM2)- T(MM1,LM2))
      DN =- (D1 +D2)/W 
      DN1 = AJ *( AL(MP1)* (3.0*T(MP1,L)-4.0*T(MP1,LM1)+T(MP1,LM2))     &
     &- AL(MM1  ) *(3.0*T(MM1,L)- 4.0*T(MM1,LM1) + T(MM1,LM2)) )        
      BSAVE = -ROPCPP * XY/ W  - G/DELTAU 
      EPT4=  -SIGP * XY /W 
      EPTB= EPT4 * TB 
      EPT4= EPT4 * T(M,L) **3 

      GO TO ( 240, 280),IROCOL 
  240 CONTINUE 
      EMM12 = (H1(M,L)+H1(MM1,L))*(H3(M,L)+H3(MM1,L))*   &
        (CKETA(M,L)+CKETA(MM1,L))/ (4.0*(H2(M,L)+H2(MM1,L)))                       
      EMP12=(H1(M,L)+H1(MP1,L))*(H3(M,L)+H3(MP1,L))*   &
       (CKETA(M,L)+ CKETA(MP1,L))/ (4.0*(H2(M,L)+H2(MP1,L)))                      
      GH =G * ETA(M)* MDOT(L)/(DELTA(L)*RO *2.0* DELETA) 
      A(M,L)= AM*U2 + EMM12/YY +GH 
      C(M,L)= AM*U1 + EMP12/YY -GH 
      BS1 (M,L)= AM* (-U1-U2) + (-EMM12 -EMP12)/YY 
      B(M) = BS1(M,L) + BSAVE  + EPT4 
      IF (IFIRST.EQ.0 )  GO TO 280 
  278 DC(M) = DN + DN1 + BSAVE*T(M,L) - VB(M  )*T(M,LM2) -AB(M,L)*      &
     & T(M,LM1)- BS1B(M,L)* T(M,L) + EPTB                               
      CYCLE   !!! GO TO 300 
  280 CONTINUE 
      PART =W * XO * DELXI 
      PART2=DMLM32/PART 
      PART1= 9.0 * DMLM12/PART 
      VB(M)    = - PART2 
      AB(M,L)  = PART1  + PART2 
      BS1B(M,L) =- PART1 
      IF (IFIRST.EQ.0 )  GO TO 278 
  290 B(L)=  BS1B(M,L) + BSAVE + EPT4 
      DC(L)=DN+DN1+(BSAVE-BS1(M,L))*  T(M,L)- A(M,L)*T(MM1,L)- C(M,L)   &
     & * T(MP1,L) + EPTB                                                
  300 END DO 
      GO TO  (301,800),IROCOL 

! STATION (S,L)    XI =1, (X=L)  ,   ETA=1, 

      ENTRY COLXL1(n1,n2)
  301 CONTINUE 
      W = 3.0 * XODXI 
      WSQ = 3.0* XDXISQ   !!! never used???
      DEDETA =  DELTA(L) * DELETA 
      TWDEL =  2.0* DELTA(L) !!! never used???
      U1=(AA(L)+AA(LM1))/(2.*(DELTA(L)+DELTA(LM1))) 
      U2=(AA(LM1)+AA(LM2))/(2.*(DELTA(LM1)+DELTA(LM2))) 
      SP=(H1(S,L)* XODXI+ 2.0*TPRIME)/(H1(S,L)*XODXI) 
      DHK = DELTA(L) * H2(S,L)/ CKETA(S,L)   *SP 
      DHK1=  DELTA(LM1)* H2(S,LM1)/ CKETA(S,LM1) 
      DHK2=  DELTA(LM2)* H2(S,LM2)/ CKETA(S,LM2) 
      ZZ2 =8.0* DELETA * E(S,L)* DELTA(L) * H2(S,L)* SP/CKETA(S,L) 
      FF=1.0/(3.0*DEDETA**2) 
      H = 8.0 * H1(S,L) * D(S,L)/CKXI(S,L) 
      PART = AA(L) /DEDETA 
      ADD = PART/3.0 
      ADD1 = (1.0 + ETA(SM1)) * PART/2.0 
      ADD2 =   (ETA(SM1) + ETA(SM2)) *PART/2.0 
      PART =       3.0* T(SM1,L)-4.0*T(SM1,LM1)+ T(SM1,LM2) 
      DSM32L=(H2(SM2,L)+H2(SM1,L))*(H3(SM2,L)+H3(SM1,L))*(CKXI(SM2,L)+  &
     & CKXI(SM1,L))/(4.*(H1(SM2,L)+H1(SM1,L)))                          
      PART2=DSM32L*(3.*T(SM2,L)-4.*T(SM2,LM1)+T(SM2,LM2)+PART) 
      DSM12L= (H2(SM1,L)+H2(S,L))* (H3(SM1,L)+H3(S,L))* (CKXI(SM1,L)    &
     & +CKXI(S,L))/(4.0* (H1(SM1,L)+H1(S,L)))                           
      PART1= -9.0*DSM12L*(3.0*T(S,L)-4.0*T(S,LM1)+T(S,LM2)+PART) 
      GSL =  H1(S,L)* H2(S,L)*H3(S,L)* RO *CP(S,L) 
      PARTW = -1.0/W + ADD 
      EPT4= SIGP * H*PARTW 
      EPTB = EPT4 * TB 
      EPT4 = EPT4 * T(S,L) **3 
      DN = ADD * (PART1 + PART2)/(4.0* XODXI)  + EPTB 
      BSAVE=H*ROPCPP*(PARTW)-GSL/DELTAU 
      GO TO (550,650),IROCOL 
  550 CONTINUE 
      AJ=GSL *MDOT(L)/ (RO*2.0*DELTA(L)*DELETA) 
      DDSL= -FF*ZZ2 
      QSAVE= DDSL* QS(L) 
      ESM32L=(H1(SM2,L)+H1(SM1,L))* (H3(SM2,L)+H3(SM1,L))* (CKETA(SM2,L)&
     &+ CKETA(SM1,L))/ (4.0*(H2(SM2,L)+H2(SM1,L)))                      
      PARTE3=FF*ESM32L 
      PARTD3= ADD*ADD2*DSM32L 
      V(L)= -PARTD3- PARTE3- AJ 
      ESM12L= (H1(SM1,L)+H1(S,L))*(H3(SM1,L)+H3(S,L))*(CKETA(SM1,L)     &
     & +CKETA(S,L))/(4.0* (H2(SM1,L)+H2(S,L)))                          
      PARTE1 = FF*9.0*ESM12L 
      PARTD1= ADD*ADD1*9.0*DSM12L 
      A(S,L)= PARTD1 + PARTD3 + PARTE3 + PARTE1 + 4.0*AJ 
      BS1(S,L)= DDSL*SIG*T(S,L)**3 - PARTD1 - PARTE1 -3.0*AJ 
      B(S)=  BS1(S,L) +  BSAVE + EPT4 
      IF (IFIRST.EQ.0 )  GO TO 650 
  648 DC( S) =  DN - VB(S  ) *T(S,LM2) -AB(S,L)*T(S,LM1)- (BS1B(S,L)    &
     & -BSAVE) * T(S,L)+ QSAVE + DDQSR                                  
      GO TO 800 
  650 CONTINUE 
      WXODXI = W* XODXI 
      DSLM3=(H2(S,LM1)+H2(S,LM2))*(H3(S,LM1)+H3(S,LM2))* (CKXI(S,LM1)   &
     &+CKXI(S,LM2))/(4.0* (H1(S,LM1)+H1(S,LM2)))                        
      DSLM1=9.0* (H2(S,L)+H2(S,LM1)) *(H3(S,L)+H3(S,LM1)) *(CKXI(S,L)+  &
     & CKXI (S,LM1))/ (4.0*(H1(S,L)+H1(S,LM1)))                         
      QSLM1 = (-DSLM1 *U1 + DSLM3* U2) *DHK1/W 
      QSLM2 = DSLM3*U2 *DHK2/W 
      QSL=-U1* DHK* DSLM1/W 
      DDQSR= QSLM1* QS(LM1) +QSLM2 *QS(LM2) + QSL*QS(L) 
      VB(S)=-DSLM3/WXODXI+QSLM2*SIG*T(S,LM2)**3 
      AB(S,L)=(DSLM1+DSLM3)/WXODXI+QSLM1*SIG*T(S,LM1)**3 
      BS1B(S,L)=-DSLM1/WXODXI+QSL*SIG*T(S,L)**3 
      IF (IFIRST.EQ.0 )  GO TO 648 
  690 B(L) = BS1B(S,L) + BSAVE + EPT4 
      DC(L)=DN+QSAVE + DDQSR                                            &
     & - V(  L) *T(SM2,L) - A(S,L) *T(SM1,L) - (BS1(S,L)-BSAVE)*T(S,L)  
  800 RETURN 
END Subroutine ColXL   ! -------------------------------------------------------

!+
SUBROUTINE SQAERO()
! ------------------------------------------------------------------------------
! PURPOSE - This routine computes the heating rates and the mass loss rates.
USE InterpolationMethods,ONLY: Ftlup,Discot
USE PickCommon
USE InputsCommon
IMPLICIT NONE
!!!SAVE

!!!      REAL MDOTO,MDOT,MCDOT,MSDOT,MWSTR,MWO2,MACHNO,MDMAX 
!!!      INTEGER S,SM1,SM2 
  INTEGER:: n
  REAL(WP):: abc
  REAL(WP):: alphac
  REAL(WP):: alphas
  REAL(WP):: bat
  REAL(WP):: block
  REAL(WP):: cat
  REAL(WP):: cell
  REAL(WP):: coll
  REAL(WP):: part
  REAL(WP):: pstag
  REAL(WP):: test

!-------------------------------------------------------------------------------
! LOOK UP  CP,  CPBAR,  CKN  ,ETC.  AS  FUNCTIONS  OF  TEMPERATURE      
      DO 11 N=1,L 
        CALL Ftlup(tt(s,n),alpha(n),malpha,nalpha,talpha,alphat) 
        CALL Ftlup(tt(s,n),hw(n),mhw,nhw,ttabhw,hwtab) 
   11 END DO

      IF (ITT.NE.1) GO TO 100 
! LOOK UP  FUNCTIONS  OF  TIME                                          
      CALL FTLUP (TAU,ALPHAC,MALPHC,NALPHC,TTALC,ALCTAB) 
      CALL FTLUP (TAU,ALPHAS,MALPHS,NALPHS,TTALS,ALSTAB) 
      CALL FTLUP (TAU,HE,MHE,NHE,TTABHE,HETAB) 
      CALL FTLUP (TAU,PSTAG,MPSTAG,NPSTAG,TTPSTAG,PSTAGTB) 
      CALL FTLUP (TAU,QC1  ,MQC,NQC,TTABQC,QCTAB) 
      CALL FTLUP (TAU,QR1,MQR,NQR,TTABQR,QRTAB) 
      CALL FTLUP (TAU,TB,MTB,NTB,TTABTB,TBTAB) 
      TB =TB**4 

! ADJUST CONVECTIVE AND RADIANT HEATING RATES AND THE PRESSURE AND      
! HEATING DISTRIBUTION TO SHAPE CHANGE  (ADJUST  QC1,QR1,PRAT,QRAT )    

      IF (CORDSY.NE.0) GO TO 20 
      CALL ADJUST()
   20 DO 30 N=1,L 
      DELTAO(N)=DELTA(N) 
      QR(N) = QR1 * QRRAT(N) 
      QC(N)= QC1 *QRAT(N) 
      PRELOC(N) = PSTAG * PRAT(N) 
      CALL Discot(tt(s,n),preloc(n),ttabhc,hctab,phc,11,28,4,hc(n)) 
      CALL Discot(tt(s,n),preloc(n),ttelam,elamtb,pelam,11,28,4,elam(n) ) 
      CALL Discot(tt(s,n),preloc(n),tthcomb,hcombtb,phcomb,11,28,4,hcomb(n))
  30 END DO                                                     
! COMPUTE QS ACROSS FRONT SURFACE                                       
      BAT = 1.0 - BETA 
  100 DO 200 N=1,L 
      CELL =HE /QC(N) 
      CAT = QC(N) * (1.0 - HW(N)/HE) 
      BLOCK=(ALPHAC *MCDOT(N) + ALPHAS *MSDOT(N))* CELL 
      QCNET(N) =  CAT *(1.0 - BAT *(0.6* BLOCK - 0.084 * BLOCK**2) - BETA * BLOCK) 
      QCOMB(N)= MCDOT(N) * HCOMB(N) 
      QS(N)= QCNET(N) + ALPHA(1)* QR(N)- MSDOT(N)*HC(N)+ QCOMB(N) ! subscript 1 on alpha added by RLC 7Sep05
  200 END DO 
      RETURN 

! THIS PART OF ROUTINE  COMPUTES  MDOTS  

      ENTRY  SQAEROM() 
      DO 1000 N=1,L 

! COMPUTE  MSDOT--- MASS LOSS RATE DUE TO SUBLIMATION                   

      IF (ASEXP ) 310,305,310 
  305 MSDOT(N)=0.0 
      GO TO 330 
  310 BLOCK =-BSEXP/TTF(S,N) 
      MSDOT(N)=   ASEXP *  PRELOC(N) **PSEXP * EXP(BLOCK)*R(1)**RIEXP 
  330 COLL = (HE-HW(N))/(QCNET(N)*ELAM(N)) 

! COMPUTE  MCDOT--- MASS LOSS RATE DUE TO OXIDATION                     

! HALF ORDER OXIDATION                                                  

  380 IF (AEXP) 390,385,390 
  385 MCDOT(N) =0.0 
      GO TO 900 
  390 MCDOT(N) = AEXP * EXP(-BEXP/TTF(S,N)) 
      IF (XORDER-0.5) 900,400,600 
  400 ABC = 4.0* MCDOT(N)**2 * PRELOC(N) * CE * RSTO2 
      PART = COLL * MCDOT(N)**2 * PRELOC(N) * RSTO2 
      TEST = ABC/ PART**2 
      IF (TEST.LT.7.E-12)GO TO  420 
      MCDOT(N) =.5*((-PART) + SQRT (PART**2 + ABC)) 
      GO TO 900 
  420 MCDOT(N) = CE /COLL 
      GO TO 900 

! First order oxidation

  600 MCDOT(N) = MCDOT(N)* PRELOC(N)* RSTO2 * CE/(1.0 + MCDOT(N)*PRELOC &
     & (N)* COLL*RSTO2)                                                 

! mdot is equal to the larger of msdot and mcdot

  900 IF (MCDOT(N).LT.MSDOT(N)) GO TO 950 
      MDOT(N)= MCDOT(N) 
      MSDOT(N)= 0.0 
      CYCLE   !!! GO TO 1000 
  950 MDOT(N)= MSDOT(N) 
      MCDOT(N)=0.0 
 1000 END DO 
      RETURN 
END Subroutine SqAero   ! ------------------------------------------------------

!+
SUBROUTINE ADJUST() 
! PURPOSE - This routine adjusts the convective and radiant heating rates,the pres
! and heating distribution to shape change (adjust qc1,qr1,prat,qrat )  
USE PickCommon
USE InputsCommon
IMPLICIT NONE

!!!      REAL MDOTO,MDOT,MCDOT,MSDOT,MWSTR,MWO2,MACHNO,MDMAX 
!!!      INTEGER S,SM1,SM2 

  REAL(WP),DIMENSION(20):: aint,al     ! aint conflicts with intrinsic
  REAL(WP):: ainto
  REAL(WP):: anum
  REAL(WP):: coef1,coef2,coef3
  INTEGER:: i,k
  INTEGER:: m,n
  INTEGER:: nm1,nm2,nmk,np1
  REAL(WP):: p0x0,p1x1,p2x2
  REAL(WP):: phi
  REAL(WP),DIMENSION(20):: psi
  REAL(WP):: sqrns
  REAL(WP):: sumh1
  REAL(WP):: tanphi
  REAL(WP),DIMENSION(20):: ueui
  REAL(WP),DIMENSION(3):: yy
!!!      DIMENSION  UEUI(20), AL(20),AINT(20),YY(3) 
!-------------------------------------------------------------------------------
!!!      NSP1 = NSTEP + 1    ! removed by RLC  7Sep05
      DO 50 N=1,L 
      RSS(N) =  RS(N) + DELTA(N)*COST(N) 
      ZS(N) = ZS(N) + (DELTAO(N) - DELTA(N))* SINT(N) 
   50 END DO
      IF (IADJUST.EQ.0) RETURN 
      RNS=(ZS(2)**2 + RSS(2)**2 -2.0*ZS(2)*ZS(1) + ZS(1)**2)/           &
     &(2.0*(ZS(2)-ZS(1)))                                               
      SQRNS =  SQRT (RNS) 

! Adjust rate to shape change
      QC1 = QC1 * SQRT ( RNSI/RNS ) 
      QR1 = QR1 * RNS/ RNSI 
      PSI(1)=0. 
      M=1  !!! never used???

  100 DO 200 N=2,L 
      NP1 = N+1 
      NM1 = N-1 
      IF (N.EQ. L) GO TO 130 
      TANPHI = (RSS(NP1) -RSS(NM1))/(ZS(NP1)- ZS(NM1)) 
      GO TO 150 
  130 TANPHI=  (RSS(L)-RSS(LM1))/(ZS(L)-ZS(LM1)) 
  150 PHI = ATAN (TANPHI) 
      PSI(N)=PID2-PHI 
  200 END DO 
! NEW PRESSURE  DISTRIBUTION                                            
      DO 250 N=1,L 
      PRAT(N) =(1.0 - GIMACH) *COS(PSI(N))**2   + GIMACH 
      UEUI(N) = SQRT((1.0+TWOGI) *(1.0-PRAT(N)**EXPG) ) 
  250 END DO 
! OBTAIN  NEW  HEAT DISTRIBUTION                                        
!                                                                       
! EVALUATE INTEGRAL AT L =0                                             
      AL(1)=0.0 
      AINTO=PRAT(1)*UEUI(1)* RSS(1)**2 
  270 CONTINUE     !!! no jumps to here???

      QRAT(1)=1.0 
      DO 600 N=2,L 
      NM1=N-1 
      NM2 =N-2 
      AINT=AINTO 
      SUMH1=0. 
      IF (N.EQ.2) GO TO 310    ! unnecessary...
      DO 300 I=2,NM1 
  300 SUMH1=SUMH1+H1(S,I) 
  310 AL(N)= X(2) *(SUMH1 + (H1(S,1)+ H1(S,N))/2.0) 

! Evaluate integral............

      IF (N.EQ. 2)  GO TO 500 
! Evaluate y(0),y(1),y(3) ...
      DO 400 K= 1,3 
      NMK = N-(3-K)    ! n-2, n-1, n
  400 YY(K)= PRAT(NMK)*UEUI(NMK)*(RSS(NMK)**2) 

  COEF2= AL(NM2)-AL(N) 
  P0X0= (AL(NM2)-AL(NM1))* COEF2 
  P1X1=(AL(NM1)-AL(NM2))* (AL(NM1)- AL(N)) 
  P2X2= (AL(N)-AL(NM2)) *  (AL(N)-AL(NM1)) 
  COEF1= (3.0*AL(NM1)-2.0* AL(NM2) - AL(N))/P0X0 
  COEF3=(2.0*AL(N) + AL(NM2)- 3.0* AL(NM1))/ P2X2 
  AINT(N)=((AL(N)- AL(NM2))**2/6.0)*(YY(1)*COEF1 + YY(2)*COEF2/P1X1 + YY(3)*COEF3)
  IF (N.GT.3)  AINT (N) = AINT (NM2) + AINT(N) 
  GO TO 590 

! N= 2                                                                  
  500 YY(2)=(PRAT(1)+PRAT(2))*(UEUI(1)+UEUI(2))*((RSS(1)+RSS(2))/2.0)**2/4.0
      YY(3)=PRAT(2)*UEUI(2)*(RSS(2)**2) 
      AINT(N)=AL(2)*(4.0*YY(2)+YY(3))/6.0 
  590 ANUM=PRAT(N)*UEUI(N)*RSS(N)*SQRNS 
      QRAT(N) = ANUM / (SQRT(AINT(N))*GG) 
  600 END DO 
      RETURN 
END Subroutine Adjust   ! ------------------------------------------------------

!+
SUBROUTINE ZPRINT() 
! ------------------------------------------------------------------------------
! PURPOSE - Print the results at each time step chosen for output.
USE PickCommon
USE InputsCommon
IMPLICIT NONE

!!!      REAL MDOTO,MDOT,MCDOT,MSDOT,MWSTR,MWO2,MACHNO,MDMAX 
!!!      INTEGER S,SM1,SM2 

  INTEGER:: m,mm,n
  REAL(WP),DIMENSION(20):: qrr

! following code removed by RLC 5 Sep 2008. 
!!!      EQUIVALENCE  (QRR(1),H1(1,1)) 
!-------------------------------------------------------------------------------
  WRITE(DBG,*) 'Entering Zprint. stebol,sig=', stebol,sig
!      DO 10 N=1,L 
!   10 QRR(N)= SIG * TTF(S,N)**4 
      qrr(1:L)=sig * ttf(s,1:L)**4

!!!      WRITE (6, 98) 
!!!   98 FORMAT ( '0') 
      WRITE (6,100) TAU,DELTAU 
  100 FORMAT (/,' TAU=',F10.4,14X,'DELTAU=',F9.6) 

      WRITE (6,101)  QC1, QR1, HE 
  101 FORMAT (/,14X,'QC=',ES11.4,5X,'QR=',ES11.4,5X,'HE=',ES11.4) 

      WRITE (6,102) T(S,1) 
  102 FORMAT (15X,'T(S,1)=',E11.4) 

      WRITE (6,105) 
  105 FORMAT (/,14X,'TEMPERATURE (M,N)') 

      WRITE (6,110)(X (N),N=1,L) 
  110 FORMAT (' ETA',6X,'X=',15F8.5/(12X,15F8.5)) 

      DO 115 M=1,S 
      MM= S- (M-1) 
  115 WRITE  (6,120) ETA(MM),(TTF(MM,N),N=1,L) 
  120 FORMAT (F6.3,6X,15F8.1/(12X,15F8.1)) 

  140 FORMAT (' ETA',6X,'X=',10(F9.5,3X)/(12X,10(F9.5,3X))) 
  150 FORMAT (F6.3,6X,10ES12.4/(12X,10ES12.4)) 
                                                                       
      WRITE (6,155) 
  155 FORMAT (/,14X,'MDOT(N)--SURFACE MASS LOSS RATE') 

      WRITE (6,140) (X (N),N=1,L) 
      WRITE (6,150) ETA(S),(MDOT(N),N=1,L) 
                                                                       
      WRITE (6,165) 
  165 FORMAT (/,14X,'MCDOT(N)--SURFACE MASS LOSS RATE DUE TO OXIDATION')

      WRITE (6,150)  ETA(S),(MCDOT(N),N=1,L) 
                                                                       
      WRITE (6,170) 
  170 FORMAT (/,14X,'DELTA(N)--MATERIAL THICKNESS') 

      WRITE (6,150)  ETA(S), (DELTA(N),N=1,L) 
                                                                       
      WRITE (6,175) 
  175 FORMAT (/,14X,'QRAT(N)--RATIO OF LOCAL HEATING TO STAGNATION HEATING')
      WRITE (6,150) ETA(S),( QRAT(N),N=1,L) 
                                                                       
      WRITE(6,176) 
  176 FORMAT(/,14X,'PRAT(N)--RATIO OF LOCAL PRESS TO STAG PRESS') 
      WRITE(6,150) ETA(S),(PRAT(N),N=1,L) 
                                                                       
      WRITE (6,180) 
  180 FORMAT (/,14X,'QS(N)--NET HEAT INPUT') 
      WRITE (6,150) ETA(S), (QS    (N),N=1,L) 
                                                                       
      WRITE (6,190) 
  190 FORMAT (/,14X,'QRR(N)--RERADIATION')
      WRITE (6,150)  ETA(S), (QRR(N),N=1,L) 
!                                                                       
      WRITE (6,200) 
  200 FORMAT (/,14X,'QCOMB(N)--HEAT DUE TO COMBUSTION FOR OXIDATION') 
      WRITE (6,150) ETA(S), (QCOMB(N),N=1,L) 
!                                                                       
      WRITE (6,400) ITC,ITR,ITTO ,IROCOL 
  400 FORMAT (/,'     NO. ITER. COL.=',I4,5X,'NO. ITER. ROW=',I4,5X, &
       'TOTAL NO. ITER.=',I8,5X,'IROCOL=',I3)                                      
      RETURN
END Subroutine Zprint   ! ------------------------------------------------------

!+
SUBROUTINE Solmat(a,b,c,z,v,d,t,n)
! PURPOSE - This routine solves the tridiagonal (except two elements) matrix.
!  Method of solution is equivalent to Gaussian elimination.
USE HoldCommon
IMPLICIT NONE 
      REAL(WP),INTENT(IN),DIMENSION(:):: a,b,c,d
      REAL(WP),INTENT(IN):: z,v
      REAL(WP),INTENT(OUT),DIMENSION(:):: t
      INTEGER,INTENT(IN):: n

      REAL(WP),DIMENSION(20):: g,sv,w
      INTEGER:: k,kk,km1
      INTEGER:: nm1,nm2
      REAL(WP):: x

!      DIMENSION W(20),SV(20),G(20),T(20),A(20),B(20),C(20),D(20) 
!!!      REAL:: tmin
!!!      COMMON /HOLD/ TMIN 
!----------------------------------------------------------------------------
      w(1)=b(1) 
      sv(1)= c(1) / b(1) 
      x= z/b(1) 
      g(1)= d(1)/w(1) 
      nm1=n-1 
      nm2=n-2 

      DO 100 k=2,n 
      km1 = k-1 
      IF (k.EQ.n) GO TO 20 
      w(k) = b(k) - a(k)*sv(km1) 
      IF (k.EQ.2) GO TO 10 
      sv(k)= c(k)/w(k)    ! was label 4
    5 g(k) = (d(k)- a(k)*g(km1))/w(k) 
      CYCLE   !!! GO TO 100 

   10 sv(2) =(c(2)-x*a(2))/w(2)   ! only if k==2
      GO TO 5 
   20 w(n)= b(n)- (a(n)- v*sv(nm2))*sv(nm1)  ! only if k==n
      g(n)=(d(n)-a(n)*g(km1)-v*g(nm2)+v*sv(nm2)*g(km1))/w(n)  ! was label 30
  100 END DO

      t(n)=g(n) 
      DO 200 k=1,nm2 
      kk= n-k 
      t(kk)= g(kk)- sv(kk)*t(kk+1) 
  200 END DO 

      t(1)= g(1)- sv(1)*t(2)- x*t(3) 
      IF (tmin.EQ.0.) RETURN 

      DO 300 k=1,n 
      IF(t(k) .LT. tmin)  t(k)=tmin 
  300 END DO 
      RETURN 
END Subroutine Solmat   ! ------------------------------------------------------
!+
SUBROUTINE DumpVector(efu,a)
! ------------------------------------------------------------------------------
! PURPOSE - Write a real vector to the debug file
IMPLICIT NONE
  INTEGER,INTENT(IN):: efu
  REAL(WP),INTENT(IN),DIMENSION(:):: a
!-------------------------------------------------------------------------------
  WRITE(efu,'(5E15.5)') a(:)
  RETURN
END Subroutine DumpVector   ! --------------------------------------------------
 
!+
SUBROUTINE DumpMatrix(efu,a)
! ------------------------------------------------------------------------------
! PURPOSE - Write a real matrix to the debug file
IMPLICIT NONE
  INTEGER,INTENT(IN):: efu
  REAL(WP),INTENT(IN),DIMENSION(:,:):: a
  INTEGER:: j,rows,cols
!-------------------------------------------------------------------------------
  rows=SIZE(a,1)
  cols=SIZE(a,2)
  DO j=1,cols
    WRITE(efu,'(A,I0)') ' COLUMN ', j
    WRITE(efu,'(5E15.5)') a(1:rows,j)
  END DO
  RETURN
END Subroutine DumpMatrix   ! --------------------------------------------------
 



END Module AbAxiProcedures   ! =================================================

!!LAR-11049                                                             
!!PROGRAM FOR THE TRANSIENT RESPONSE OF ABLATING AXISYMMETRIC BODIES INC
!!EFFECTS OF SHAPE CHANGE                                               
!!JOB,1,0700,75000.           D2430      32736T                   BIN204
!!USER.HOWSER, LONA M                 000425325  11160                  
!!RUN(S)                                                                
!!SETINDF.                                                              
!!LGO.                                                                  
!!      PROGRAM D2430 (INPUT,OUTPUT,TAPE5=INPUT,TAPE6=OUTPUT,TAPE7=201, 
!!     1TAPE8=201,TAPE9=201)        

!+                                    
PROGRAM LAR2430 
!                                                                       
! AXISYMMETRIC ABLATION PROGRAM                                         
! TWO-DIMENSIONAL ABLATION ANALYSIS FOR AXIALLY SYMMETRIC BODIES OF REVOLUTION
! AT HIGH HEATING RATES, CONSIDERING SHAPE CHANGE

! THIS IS THE MAIN PROGRAM - IT CONTROLS THE GENERAL FLOW OF PROGRAM    
! ------------------------------------------------------------------------------
USE AbaxiConstants
USE InterpolationMethods,ONLY: Ftlup,Discot
USE AbAxiProcedures
USE PickCommon
USE InputsCommon
USE HoldCommon
IMPLICIT NONE

!      DATA XLABEL,YLABEL,X2L,Y2L,Y3L/ 2HZB,3HRSS,1HX,4HMDOT,12HTEMPERATU&
!     &RES/                                                              
  CHARACTER(LEN=*),PARAMETER:: XLABEL = 'ZB'   ! never used???
  CHARACTER(LEN=*),PARAMETER:: YLABEL = 'RSS'  ! never used???
  CHARACTER(LEN=*),PARAMETER:: X2L = 'X'       ! never used???
  CHARACTER(LEN=*),PARAMETER:: Y2L = 'MDOT'    ! never used???
  CHARACTER(LEN=*),PARAMETER:: Y3L = 'TEMPERATURES'  ! never used???
  CHARACTER(LEN=80):: genTitle,fileName
  INTEGER:: errCode
  INTEGER,PARAMETER:: FIG=11
                                       
  REAL(WP):: abstt,absttf   ! convergence tests
  REAL(WP):: alm1   ! L-1 converted to REAL 
  REAL(WP):: am     !
  REAL(WP),DIMENSION(10,20):: delt
  REAL(WP):: delt1  ! something to do with convergence
  REAL(WP):: deltn
  REAL(WP):: delx
  REAL(WP):: dtau1
  INTEGER:: i,m,mm,n

  INTEGER:: iec   ! something to do with plotting...
  INTEGER:: iplt,ipltk   ! something to do with plotting...
  INTEGER:: isym  ! plotting symbol???

  INTEGER:: inop
  INTEGER:: irow
  INTEGER:: nntp
  REAL(WP):: tauoo
  REAL(WP):: test   ! convergence test quantity
  REAL(WP),DIMENSION(22):: zz                             
!!!      REAL MDOTO,MDOT,MCDOT,MSDOT,MWSTR,MWO2,MACHNO,MDMAX 
!!!      INTEGER S,SM1,SM2 
!!!      COMMON /HOLD/ TMIN 

!      NAMELIST /D2430/       AEXP,ALCTAB,TTALC,MALPHC,NALPHC,ALPHAT,    &
!     & TALPHA,MALPHA,NALPHA,ALSTAB,TTALS,MALPHS,NALPHS,ASEXP,           &
!     & BETA,BEXP,BSEXP,CE,CKETATB,ETATAB,TTCKETA,NCKETA,NETA,           &
!     & CKXITAB,XITAB,TTCKXI,NCKXI,NXI,CORDSY,CPDP,CPP,CPTAB,TTABCP,MCP, &
!     & NCP,DELTAO,DELTAU,DELTMIN,DTMAX,ELAMTB,TTELAM,PELAM,NELAM,NPELAM,&
!     & ENDTAU,EPSONE,EPSONEP,EPSONPP,ERRORT,GAMBAR,GAMINF,HCOMBTB,      &
!     & TTHCOMB,PHCOMB,NHCOMB,NPHCOMB,HCTAB,TTABHC,PHC,NHC,NPHC,HETAB,   &
!     & TTABHE,MHE,NHE,HWTAB,TTABHW,MHW,NHW,IADJUST,IPLOT,L,MACHNO,      &
!     & MAXITT,MDMAX,MDOTO,MWO2,MWSTR,NTP,PLTIME,PRAT,PRFREQ,PSEXP,      &
!     & PSTAGTB,TTPSTAG,MPSTAG,NPSTAG,PTMAX,PTMIN,QCTAB,TTABQC,MQC,      &
!     & NQC,QRAT,QRRAT,QRTAB,TTABQR,MQR,NQR,R,RIEXP,RNSI,RO,RODP,ROP,RS, &
!     & RSSMAX,S,STEBOL,T,TAUO,TBTAB,TTABTB,MTB,NTB,TDPRIME,THETA,       &
!     & TMIN,TPRIME,XO,XORDER,ZS,ZSMAX          

      NAMELIST /D2430/ AEXP,ALCTAB,ALPHAT,ALSTAB,ASEXP,BETA,BEXP,        &
     & BSEXP,CE,CKETATB,CKXITAB,CORDSY,CPDP,CPP,CPTAB,DELTAO,DELTAU,     &
     & DELTMIN,DTMAX,ELAMTB,ENDTAU,EPSONE,EPSONEP,EPSONPP,ERRORT,        &
     & ETATAB,GAMBAR,GAMINF,HCOMBTB,HCTAB,HETAB,HWTAB,IADJUST,IPLOT,     &
     & L,MACHNO,MALPHA,MALPHC,MALPHS,MAXITT,MCP,MDMAX,MDOTO,MHE,MHW,     &
     & MPSTAG,MQC,MQR,MTB,MWO2,MWSTR,NALPHA,NALPHC,NALPHS,NCKETA,        &
     & NCKXI,NCP,NELAM,NETA,NHC,NHCOMB,NHE,NHW,NPELAM,NPHC,NPHCOMB,      &
     & NPSTAG,NQC,NQR,NTB,NTP,NXI,PELAM,PHC,PHCOMB,PLTIME,PRAT,PRFREQ,   &
     & PSEXP,PSTAGTB,PTMAX,PTMIN,QCTAB,QRAT,QRRAT,QRTAB,R,RIEXP,RNSI,    &
     & RO,RODP,ROP,RS,RSSMAX,S,STEBOL,T,TALPHA,TAUO,TBTAB,TDPRIME,       &
     & THETA,TMIN,TPRIME,TTABCP,TTABHC,TTABHE,TTABHW,TTABQC,TTABQR,      &
     & TTABTB,TTALC,TTALS,TTCKETA,TTCKXI,TTELAM,TTHCOMB,TTPSTAG,XITAB,   &
     & XO,XORDER,ZS,ZSMAX
!-------------------------------------------------------------------------------                         

  DO
    WRITE(*,*) 'Enter the name of the input file:'
    READ(*,'(A)',IOSTAT=errCode) fileName
    IF (LEN_TRIM(fileName)==0) STOP
    OPEN(UNIT=5,FILE=fileName,STATUS='OLD',IOSTAT=errCode, &
      ACTION='READ',POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) 'Unable to open this file. Try again.'
  END DO
  INQUIRE(UNIT=5,NAME=fileName)
  WRITE(*,*) 'Reading from ', TRIM(fileName)
  OPEN(UNIT=DBG,FILE='abaxi.dbg',STATUS='REPLACE',ACTION='WRITE')
  WRITE(DBG,*) 'Reading from '//Trim(fileName)
  OPEN(UNIT=6,FILE='abaxi.out',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=7,FILE='abaxi.f7',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='READWRITE')
  OPEN(UNIT=8,FILE='abaxi.f8',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='READWRITE')
  OPEN(UNIT=9,FILE='abaxi.f9',FORM='UNFORMATTED',STATUS='REPLACE',ACTION='READWRITE')
  OPEN(UNIT=FIG,FILE='abaxi.fig',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=GNU1,FILE='abaxi.gnu1',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=GNU2,FILE='abaxi.gnu2',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=GNU3,FILE='abaxi.gnu3',STATUS='REPLACE',ACTION='WRITE')
                       
! zero all of the input quantities                                                 
      TMIN = 0
  CALL InitializeInputQuantities()
!!!      DO I=1,934 
!!!        DUMMY(I)=0.0 
!!!      END DO
!!!  CALL Junk()

      DTMAX=2. 
!(RLC) Program input. Reads one title line and then the namelist input. That is all.
! keep coming back here for more cases...
    1 WRITE(DBG,*) 'Starting a new case at statement #1'
      READ(5,'(A)',IOSTAT=errCode) genTitle
      IF (errCode < 0) STOP
      READ(5,D2430,IOSTAT=errCode)
      IF (errCode < 0) THEN
        WRITE(*,*) 'End of file reading namelist D2430'
        STOP
      END IF

!    1 READ (5,100) 
!  100 FORMAT (80H                                                       &
!     &                         )                                        
!      IF (EOF,5) 2,3 
!    2 STOP 
!    3 READ (5,D2430) 

      WRITE (DBG,D2430) 
      WRITE (6,*) genTitle
!                                                                       
! SET INITIAL VALUES                                                    
!                                                                       
      NNTP= NTP(1) 
      PID2 = 1.5707963268 
      TWOGI  = 2.0 /((GAMINF - 1.0) * MACHNO **2) 
      EXPG =(GAMBAR - 1.0)/ GAMBAR 
      GIMACH=    1./(GAMINF * MACHNO **2) 
      GG= SQRT( EXPG * (1.0 + TWOGI) * (1.0 - GIMACH)) 
      GG= SQRT (GG) * 2.0 
      INOP=0 
      IROW=0 
!!!      IDT=1          !!! never used???
!!!      DTAUO=1.0      !!! never used???
      DTAU1=DELTAU 
      IROCOL =1 
! WILL PRINT ONLY AFTER A COL. AND ROW COMPUTATION HAS BEEN MADE        
      TAUOO = TAUO + PRFREQ   
      ITTO=0 

      DO 11 M=1,S 
      DO 11 N=1,L 
      DELT(M,N)=1000. 
   11 TT(M,N)= T(M,N) 

      DELTAU=DELTAU/2. 
      TAU=TAUO+DELTAU 
  WRITE(DBG,'(A,5F12.4)') 'Initial tauZero,deltau,tau', tauo,deltau,tau
      IFIRST=0 
      ITT=1 
      LM1 = L-1 
      ALM1 = LM1 
      LM2 = L- 2 
      SM1 = S- 1 
      SM2 = S- 2 
      DELXI =1./ALM1 
      DELX =XO/ALM1 
      RSTO2 = MWSTR/MWO2 
      X(1) =0. 

      DO  N=2,L 
        X(N) = X(N-1) + DELX 
      END DO

      DELETA = 1./SM1 
      DELXISQ = DELXI**2 
      DELESQ= DELETA**2 
      TWDELXI = 2.0*DELXI 
      EIGHT3=8.0/3.0 

      DO M=1,S 
        AM=M-1 
        ETA(M)=DELETA*AM 
      END DO

      SIGMA=STEBOL 
      SIG = SIGMA* EPSONE 
      SIGP = SIGMA * EPSONEP 
      SIGDP= SIGMA * EPSONPP 
      XODXI = XO * DELXI 
      RODPC  = TDPRIME*RODP * CPDP / DELTAU 
      ROPCPP = TPRIME * ROP * CPP/ DELTAU 
!!!      RODT= RO/DELTAU                        !!! never used???
      XDXISQ = XO**2 * DELXISQ 

      DO N=1,L 
        MDOT(N)=MDOTO(N) 
        MCDOT(N)=MDOTO(N) 
        MSDOT(N)=MDOTO(N) 
        DELTA(N)= DELTAO(N) 
        THETA(N)=.0174532925*THETA(N) 
        SINT(N) = SIN(THETA(N)) 
        ZZ(N)= ZS(N)+DELTAO(N)*SINT(N) 
        COST(N)= COS(THETA(N)) 
      END DO

      IF (IPLOT.EQ.0) GO TO 23 
! PLOT BASE CURVE IF PLOTTING IS CALLED FOR                             
      REWIND 7 
      REWIND 8 
      REWIND 9 
!!!      CALL CALCOMP 
      IPLT=1 
      IPLTK=0 
      IF (CORDSY==0.0) THEN
        WRITE(7)(ZZ(N),RS(N),N=1,L) 
      ELSE
        WRITE(7)(ZS(N),DELTA(N),N=1,L) 
      END IF

! COMPUTE  H-S                                                          
   23 DO 25 M=1,S 
      DO 26 N=1,L 
      Y(M,N)=ETA(M)*DELTA(N) 
      H1(M,N)= 1.0 + ETA(M)* DELTA(N)/R(N) 
      H2(M,N)= 1.0 
      H3(M,N)= RS(N) + Y(M,N)*COST(N) 
   26 END DO
   25 END DO
      WRITE(DBG,*) 'Scale factors h1,h2,h3 have been computed.'
!  WRITE(DBG,*) 'h1'
!  CALL DumpMatrix(DBG,h1(1:s,1:l))
!  WRITE(DBG,*) 'h2'
!  CALL DumpMatrix(DBG,h2(1:s,1:l))
!  WRITE(DBG,*) 'h3'
!  CALL DumpMatrix(DBG,h3(1:s,1:l))

   95 DO 101 M=1,S 
      DO 102 N=1,L 
      CALL FTLUP (TT(M,N),CP(M,N),MCP,NCP,TTABCP,CPTAB)
      CALL DISCOT(TT(M,N),  X(N), TTCKXI,CKXITAB, XITAB,11, NCKXI, NXI, CKXI(M,N))
      CALL DISCOT(TT(M,N),Y(M,N),TTCKETA,CKETATB,ETATAB,11,NCKETA,NETA,CKETA(M,N))
  102 END DO
  101 END DO   
  WRITE(DBG,*) 'Finished DISCOT at tau=', tau
                                                 
      AA(1)=0.0 
      DO 109 N=2,LM1 
        AA(N)= (DELTA(N+1)-DELTA(N-1))/(TWDELXI*XO) 
  109 END DO

      AA(L)=(3.0*DELTA(L)-4.0*DELTA(LM1)+DELTA(LM2))/(TWDELXI*XO) 
      DO 110 N=1,L 
      DO 111 M=1,S 
      D(M,N)=H2(M,N)*H3(M,N)*CKXI(M,N)/H1(M,N) 
      E(M,N)=H1(M,N)*H3(M,N)*CKETA(M,N)/H2(M,N) 
      F(M,N)=D(M,N)*ETA(M)*AA(N)/DELTA(N) 
  111 CONTINUE
  110 CONTINUE

      CALL SQAERO()
      WRITE(DBG,*) 'Finished SqAero. IROCOL=', irocol

      GO TO (310,320), IROCOL 
  310 CALL COLUMN()
      ITC=ITT 
      IFIRST=1 
      GO TO 350 
  320 CALL ROW()
      ITR=ITT 
      IF (IROW.EQ.0) IROW=2 
  350 CONTINUE 
      WRITE(DBG,*) 'Finished column/row. Dump TTF'
      CALL DumpMatrix(DBG,ttf(1:s,1:l))

! IF ANY TEMPERATURES ARE NEGATIVE  STOP  CALCULATIONS                  
      DO 360 N=1,L 
      DO 361 M=1,S 
!      IF (TTF(M,N).LE.0) GO TO 411 
      IF (TTF(M,N).LE. 0.0) THEN
        CALL ZPRINT()
        WRITE(DBG,*) 'Program halted because of negative temperatures'
        STOP
      END IF
  361 END DO
  360 END DO

! TEST TO SEE IF TEMPERATURES HAVE CONVERGED

      ITTO=ITTO+1 
      DO 400 N=1,L 
      DO 401 M=1,S 
      ABSTT=ABS(TT(M,N)) 
      ABSTTF=ABS(TTF(M,N)) 
      TEST=ABS(ABSTTF-ABSTT)/ABSTT 
      IF (TEST > ERRORT) GO TO 700
  401 END DO
  400 END DO
  WRITE(DBG,*) 'Loop 400 did not jump to 700'

! COMPUTE MDOT
      CALL SQAEROM()

! COMPUTE DELTA
      DO 410 N=1,L 
      DELTA(N)=DELTAO(N)-(MDOTO(N)+MDOT(N))*DELTAU/(2.0*RO) 
! RESET DELTAO AND MDOTO
      MDOTO(N)=MDOT(N) 
  410 END DO

! IF DELTA BECOMES LESS THAN DELTMIN (SOME MINIMUM DELTA INPUT) STOP    
! THE CALCULATIONS                                                      

      DO N=1,L 
!      IF (DELTA(N).GT. DELTMIN) GO TO 412 
!  411 CALL ZPRINT 
        IF (DELTA(N) .GT. DELTMIN) CYCLE
        CALL ZPRINT()
        WRITE(DBG,*) 'Program halted because delta less than deltmin'
        STOP 
      END DO 


      IF (INOP.EQ.1) GO TO 418 
      IF (TAU.LT.TAUOO) GO TO 420 
      IF (IROCOL.EQ.1) GO TO 418 
      INOP=1 
      GO TO 420 
  418 INOP =0 
      TAUOO=TAUOO+ PRFREQ 

      CALL ZPRINT() 

      IF (IPLOT.EQ.0) GO TO 420 
      IPLTK= IPLTK + 1 
      WRITE(8) (MDOT(N), N=1,L) 
      IF (NNTP.EQ.0) GO TO 420 
      DO 419 M=1,NNTP 
        I= NTP(M+1) 
        WRITE (9) (TTF(I,N),N=1,L) 
  419 END DO

  420 IF (IROW-1) 540,490,484 
  484 DELTAU=DELTAU*2.0 
      IROW=1 
!!!      KFRE=KFRE+1    ! this statement seems to have no purpose ????

! OBTAIN DELTAU AS A FUNCTION OF ITERATION OF PREVIOUS TIME STEP
  490 DTAU1 = DELTAU 
      IF (IROCOL.EQ.1) GO TO 540 
      IF (ITT-2) 495,540,530 
  495 DELTAU=2.0*DTAU1 
      IF (DELTAU.GT.DTMAX) DELTAU=DTMAX 
      GO TO 540 
  530 DELTAU=DTAU1/2. 
      IF (DELTAU.LT.1.E-6) GO TO 900 
  540 TAUO = TAU 

! CHECK TO SEE IF IT IS TIME TO PLOT                                    
      IF (IPLOT.EQ.0) GO TO 543 
      IF (TAU.LT.PLTIME(IPLT)) GO TO 543 
      IPLT=IPLT+1 
      IF (CORDSY.NE.0) GO TO 542 
      WRITE (7) (ZS(N),RSS(N),N=1,L) 
      GO TO 543 
  542 WRITE (7) (ZS(N),DELTA(N),N=1,L) 

! INCREMENT TIME AND REPEAT CYCLE ALTERNATING ROW AND COLUMN SOLUTION
  543 TAU=TAU+DELTAU 
      RODPC  = TDPRIME*RODP * CPDP / DELTAU 
      ROPCPP = TPRIME * ROP * CPP/ DELTAU 
!!!      RODT= RO/DELTAU   ! never used ??
  WRITE(DBG,'(A,5F12.4)') ' Update tau,deltau,rodpc,ropcpp', tau,deltau,rodpc,ropcpp
      IF (TAU.GT.ENDTAU)  GO TO 950 ! normal finish

! EXTRAPOLATE TO GET NEW GUESS TEMP(TT)

      DO 446 M=1,S 
      DO 447 N=1,L 
      DELT(M,N)=1000. 
      DELTN=TTF(M,N)-T(M,N) 
      T(M,N)=TTF(M,N) 
      TT(M,N)=TTF(M,N)+(DELTAU/DTAU1)*DELTN 
  447 END DO
  446 END DO

      GO TO (550,650),IROCOL 
  550 IROCOL = 2 
      ITT=1 
      GO TO 23 
  650 IROCOL = 1 
      ITT=1 
      GO TO 23    ! these are the only jumps back to 23

! TEMP. DOES NOT MEET ERROR CRITERIA,  MUST ITERATE AGAIN 
! NEW GUESS IS TEMP. OF  PREVIOUS ITERATION   TT =TTF 

  700 ITT =ITT +1 
      IF (ITT > MAXITT) GO TO 800
      DO N=1,L 
        DO M=1,S 
          DELT1 = ABS(TTF(M,N)- TT(M,N)) 
          IF (DELT1 .LT. 10.) GO TO 718 
          IF (DELT1-DELT(M,N)) 718,750,750 
  718     DELT(M,N)=DELT1 
        END DO
      END DO

  WRITE(DBG,*) 'temp does not meet error criteria. Jump back to 95'
      DO M=1,S 
        DO N=1,L 
          TT(M,N)= TTF(M,N)   ! use previous iteration
        END DO
      END DO
      GO TO 95 

  750 IF (ITT.LT.3) GO TO 700   !! was 718 changed by RLC

! PROGRAMED STOPS

      WRITE (6,752) 
  752 FORMAT (/,'TEMPERATURE IS DIVERGING ----- WHY') 
  758 WRITE (6,759) 
  759 FORMAT (/' TT(M,N)') 
      DO M=1,S 
        MM=S-(M-1) 
        WRITE (6,766) ETA(MM),(TT(MM,N),N=1,L) 
      END DO
  766 FORMAT (F6.3,6X,15F8.1/(12X,15F8.1)) 
      WRITE (6,767) IROCOL 
  767 FORMAT (/'IROCOL=',I3) 
      CALL ZPRINT()
      WRITE(DBG,*) 'Program halted because of diverging temperature'
      STOP 

  800 IF (IROCOL.EQ.1) GO TO 803 
      WRITE (6,801) 
  801 FORMAT (/,'THIS IS A ROW SOLUTION, DELTAU CANNOT CHANGE') 
      GO TO 758 
!                                                                       
  803 DTAU1= DELTAU 
      DELTAU = DELTAU/2.0 
      WRITE (6,805) DELTAU ,TAU 
  805 FORMAT (/,' I DID IT--      DELTAU=',E14.5,'TAU=',E14.5) 
      IF (DELTAU .LT. 1.E-6)  GO TO 900 
      TAU = TAU - DELTAU 
      DO 810 M=1,S 
      DO 811 N=1,L 
      DELT(M,N)=1000. 
      TT(M,N) = T(M,N) 
  811 END DO
  810 END DO

      ITT  = 1 
      GO TO 95 
  900 WRITE (6,901) 
  901 FORMAT  (/,' TEMPERATURE ITERATION DOES NOT CONVERGE') 
      GO TO 758 
!                                                                       
! PLOT  ZS  VS. RSS ,   X  VS MDOT ,  X  VS BACK SURFACE TEMPERATURE    
!                                                                       
  950 CALL ZPRINT()
      IF (IPLOT.EQ.0) GO TO 1 
      END FILE 7 
      END FILE 8 
      END FILE 9 
      REWIND 7 
      REWIND 8 
      REWIND 9 
      IEC = 0    ! never used???

      DO M=1,IPLT 
        READ (7)  (ZZ(N), RSS(N),N=1,L) 
        IF (M.EQ.IPLT )  IEC =1 
!  960 CALL INFOPLT (IEC,L,ZZ,1,RSS,1,0.,ZSMAX,0.,RSSMAX,1.,10,XLABEL,10,&
!     & YLABEL,0)
        WRITE(FIG,*) 'DATA ', L
        WRITE(FIG, '(2ES12.4)' ) (zz(n),rss(n),n=1,L)
        WRITE(GNU1, '(2ES12.4)' ) (zz(n),rss(n),n=1,L)
        WRITE(GNU1,*) ' '
      END DO
  WRITE(FIG,*) 'ENDFRAME'
  
  WRITE(FIG,*) '--- TAPE 8 plots...'                                                  
      IEC =0 
      DO M=1,IPLTK 
        READ(8) (MDOT(N),N=1,L) 
        IF (M.EQ.IPLTK)  IEC= 1 
!  970 CALL INFOPLT (IEC,L,X,1,MDOT,1,0.,0.,0.,MDMAX,1.,10,X2L,10,Y2L,0) 
        WRITE(FIG,*) 'DATA ', L
        WRITE(FIG, '(2ES12.4)' ) (x(n),mdot(n),n=1,L)
        WRITE(GNU2, '(2ES12.4)' ) (x(n),mdot(n),n=1,L)
        WRITE(GNU2,*) ' '
      END DO
  WRITE(FIG,*) 'ENDFRAME ------------------------------------------------------'

  WRITE(FIG,*) '--- TAPE 9 plots...'
      IEC =0 
      IF (NNTP.EQ.0) GO TO 1 
      DO M=1,IPLTK 
        ISYM=10 
        DO I=1,NNTP 
        READ (9)   (ZZ(N),N=1,L) 
        IF (M.EQ.IPLTK .AND. I.EQ.NNTP)  IEC =1 
        ISYM= ISYM + 1 
!  980 CALL INFOPLT (IEC,L,X,1,ZZ,1,0.,0.,PTMIN,PTMAX,1.,10,X2L,20,Y3L,  &
!     & ISYM)                                        
        WRITE(FIG,*) 'DATA ', L
        WRITE(FIG, '(2ES12.4)' ) (x(n),zz(n),n=1,L)                    
        WRITE(GNU3, '(2ES12.4)' ) (x(n),zz(n),n=1,L)                    
        WRITE(GNU3, *) ' '
        END DO
      END DO 
  WRITE(FIG,*) 'ENDFRAME ------------------------------------------------------'

      GO TO 1 
END PROGRAM LAR2430   ! ==================================================
                                           
