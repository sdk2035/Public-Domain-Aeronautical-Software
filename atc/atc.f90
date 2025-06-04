!+
! PROGRAM AircraftMotions
! ------------------------------------------------------------------------------
! PURPOSE - Analyze the ....
! AUTHORS - unknown, NASA Ames Research Center
!           Ralph L. Carmichael, Public Domain Aeronautical Software
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1971   1.0   ???   Original coding at NASA Ames
!   1995   1.1   RLC   Acquisition of source code from COSMIC
!   1998   1.2   RLC   Conversion to F90 free form source
!   1999   1.21  RLC   Compilation with Fortran 90; minor changes
! 21Jul08  1.3   RLC   xxx

! NOTES - 
!  Tasks for this program
!  1. Make the output fit letter paper, portrait orientation
!  2. Replace TAINT with a simple linear interpolation scheme
!  3. Remove unreferenced FORMATs
!  4. Write a description of the equations used for calculation
!  5. Don't use SAVE as a variable name


! This program was developed by Ames Research Center, in cooperation with 
! the National Transportation Safety Board, as a technique for deriving 
! time histories of an aircraft's motion from Air Traffic Control (ATC)
! radar records. This technique uses the radar range and azimuth data, 
! along with the downlinked altitude data, to derive an expanded set of 
! data which includes airspeed, lift, attitude angles (pitch, roll, and 
! heading), etc.
! This technique should prove useful as a source of data in the investigation
! of commercial airline accidents and in the analysis of accidents involving
! aircraft which do not have onboard data recorders (e.g., military, 
! short-haul, and general aviation).

! The technique used to determine the aircraft motions involves smoothing 
! of raw radar data. These smoothed results, in combination with other
! available information (wind profiles and aircraft performance data), are
! used to derive the expanded set of data. This program uses a cubic 
! least-square fit to smooth the raw data. This moving-arc procedure provides 
! a smoothed time history of the aircraft position, the inertial velocities,
! and accelerations. Using known winds, these inertial data are transformed
! to aircraft stability axes to provide true airspeed, thrust-drag, lift, and
! roll angle. Further derivation, based on aircraft dependent performance
! data, can determine the aircraft angle of attack, pitch, and heading angle.
! Results of experimental tests indicate that values derived from ATC radar
! records using this technique agree favorably with airborne measurements.

!+
MODULE AircraftMotionsProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Collect all of the procedures used in the program AircraftMotions

CONTAINS

!+
SUBROUTINE Smooth (x, dt, n, xs, xdot, xaccel)
! ---------------------------------------------------------------------------
! PURPOSE:
IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: x
  REAL,INTENT(IN):: dt
  INTEGER,INTENT(IN):: n
  REAL,INTENT(OUT):: xs
  REAL,INTENT(OUT):: xdot
  REAL,INTENT(OUT):: xaccel
  
  REAL:: a2
  REAL:: c11,c13,c22,c31,c33
  INTEGER:: i
  INTEGER:: p
  REAL:: pn
  REAL:: q1,q2
  REAL:: s1,s2,s3
!----------------------------------------------------------------------------
      pn=n
      p=n-((n-1)/2)
      q1=3.*(3.*pn**2-7.)
      q2=4.*pn*(pn**2-4.)
      c11=q1/q2
      c13=-15./(pn*(pn**2-4.))
      c22=12./(pn*(pn**2-1.))
      c31=c13
      c33=c22*(15./(pn**2-4.))
      s1=0.0
      s2=0.0
      s3=0.0
      DO i=1,n
        s1=s1+x(i)
        s2=s2+x(i)*(i-p)
        s3=s3+x(i)*(i-p)**2
      END DO
      xs=c11*s1+c13*s3
      xdot=(c22*s2)/dt
      a2=c31*s1+c33*s3
      xaccel=2.*a2/(dt**2)
  RETURN
END Subroutine Smooth   ! ---------------------------------------------------

!+
SUBROUTINE Taint (xtab, ftab, x, fx, n, k, ner, mon)
! ---------------------------------------------------------------------------
! PURPOSE: Table interpolation
IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: xtab,ftab
  REAL,INTENT(IN):: x
  REAL,INTENT(OUT):: fx
  INTEGER,INTENT(IN):: n
  INTEGER,INTENT(IN):: k
  INTEGER,INTENT(OUT):: ner
  INTEGER,INTENT(IN OUT):: mon

  REAL,DIMENSION(10):: c,t
  INTEGER:: i,j,jsave,kp1,l,m,nm1
!----------------------------------------------------------------------------
!!!      DIMENSION xtab(*), ftab(*), t(10), c(10)
!!!PS0400  TAINT SUBROUTINE- IN FORTRAN II.
      IF (n-k) 10,10,20
   10 ner=2
      RETURN
   20 IF (k-9) 30,30,10
   30 IF (mon) 50,50,40
   40 IF (mon-2) 140,110,50
   50 j=0
      nm1=n-1
      DO 90 i=1,nm1
        IF (xtab(i)-xtab(i+1)) 70,60,80
   60   ner=3
        RETURN
   70   j=j-1
        CYCLE   !!!GO TO 90
   80   j=j+1
   90 END DO
      mon=1
      IF (j) 100,140,140
  100 mon=2
  110 DO 130 i=1,n
!!!        IF (x-xtab(i)) 120,120,130
        IF (x > xtab(i)) CYCLE
  120   j=i
        GO TO 180
  130 END DO
      GO TO 170

  140 DO 160 i=1,n
!!!        IF (x-xtab(i)) 160,150,150
        IF (x < xtab(i)) CYCLE
  150   j=i
        GO TO 180
  160 END DO
  170 j=n
  180 j=j-(k+1)/2
      IF (j) 190,190,200
  190 j=1
  200 m=j+k
      IF (m-n) 220,220,210
  210 j=j-1
      GO TO 200
  220 kp1=k+1
      jsave=j
      DO 230 l=1,kp1
        c(l)=x-xtab(j)
        t(l)=ftab(j)
  230   j=j+1
      DO 250 j=1,k
        i=j+1
  240   t(i)=(c(j)*t(i)-c(i)*t(j))/(c(j)-c(i))
        i=i+1
!!!        IF (i-kp1) 240,240,250
        IF (i <= kp1) GO TO 240
  250 END DO
      fx=t(kp1)
      ner=1
  RETURN
END Subroutine Taint   ! ----------------------------------------------------

!+
FUNCTION SingleLookup(xtab,x) RESULT (i)
! ---------------------------------------------------------------------------
! PURPOSE - Search a sorted (increasing) array to find the interval
!  bounding a given number. If n is the size of the array a,
!  return 0 if number x is < a(1)
!  return n if x > a(n).
!  return i if  a(i) <= x < a(i+1).
!  If x is exactly equal to a(n), return n-1.
  INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6)
  REAL(SP),INTENT(IN),DIMENSION(:)::  xtab    ! input array                        
  REAL(SP),INTENT(IN):: x  ! input number              

  INTEGER:: i  ! index of interval such that xtab(i) <= x < xtab(i+1)
               ! =0 if x < xtab(1) and =n if x > xtab(i)
  INTEGER:: j,k,n
  INTEGER,SAVE:: isave = 1
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

  IF (isave>0 .AND. isave<n) THEN
    IF ((x >= xtab(isave)) .AND. (x < xtab(isave+1)) ) THEN
      i=isave
      RETURN
    END IF
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
  isave=i

  RETURN
END Function SingleLookup   ! ===============================================



!+
SUBROUTINE Xplot(efu1,efu2)
! ---------------------------------------------------------------------------
! PURPOSE:
IMPLICIT NONE
  INTEGER,INTENT(IN):: efu1,efu2

!!!      REAL line
!!!      DIMENSION nwy(30), nwx(30)
  CHARACTER(LEN=30):: nwx,nwy 
  CHARACTER(LEN=1),DIMENSION(102):: line
  CHARACTER(LEN=1),PARAMETER:: DOT='*',BLANK=' ',PLUS='+'


  INTEGER:: errCode
  CHARACTER(LEN=80):: fileName
  INTEGER:: i,j,m,n
  INTEGER:: mm

  INTEGER:: nx,ny
  REAL:: start,stop
  REAL:: startx,stopx
  REAL:: t1,t2,t3,tm,tx,ty
  REAL:: value

  REAL:: dt
  INTEGER:: nb,ne
  REAL,DIMENSION(30,300):: save
  REAL,DIMENSION(300):: time
  COMMON /dsave/ save,nb,ne,dt,time
!----------------------------------------------------------------------------
  DO
    WRITE(*,*) 'Enter name of the plot file:'
    READ(*,'(A)') fileName
    IF (LEN_TRIM(fileName)==0) STOP
    OPEN(UNIT=efu1,FILE=fileName,STATUS='OLD', &
      IOSTAT=errCode,ACTION='READ',POSITION='REWIND')
    IF (errCode==0) EXIT
    OPEN(UNIT=efu1,FILE=TRIM(fileName)//'.inp',STATUS='OLD', &
      IOSTAT=errCode,ACTION='READ',POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) 'Unable to open this file. Try again.'
  END DO
  INQUIRE(UNIT=efu1,NAME=fileName)
  WRITE(*,*) 'Reading from ', Trim(fileName)


   10 FORMAT (I10,A30,2F10.2)
   20 READ (efu1,10,IOSTAT=errCode) ny,nwy,start,stop
!      IF (eof(5).ne.0) GO TO 140
      IF (errCode < 0) GO TO 140
!   30 FORMAT ('1')
      WRITE (efu2,'(A)') ACHAR(12)
      tm=stop-start
      t1=start+tm*0.25
      t2=start+tm*0.5
      t3=start+tm*0.75
      ty=100.0/tm
      READ (efu1,10) nx,nwx,startx,stopx
      tx=40.0/(stopx-startx)
      DO 120 m=1,41
        line(:)=blank
!        DO 40 j=1,102
!   40     line(j)=blank

        line(1)=dot
        line(101)=dot
        IF (m.eq.1) GO TO 60
        IF (m.eq.41) GO TO 60
        IF (m.eq.2) GO TO 50
        IF (m.eq.40) GO TO 50
        GO TO 80
   50   CONTINUE
        line(26)=dot
        line(51)=dot
        line(76)=dot
        GO TO 80
   60   CONTINUE
        DO 70 j=1,101
   70     line(j)=dot
   80   CONTINUE
        i=42-m

        DO 100 n=nb,ne
          mm=tx*(save(nx,n)-startx)+1.5   ! converts real to integer
          IF (mm.ne.i) GO TO 100
          j=ty*(save(ny,n)-start)+1.5
          IF (j.lt.1) GO TO 90
          IF (j.gt.102) GO TO 90
          line(j)=plus
   90     CONTINUE
  100   CONTINUE

        value=startx+(i-1)/tx
        WRITE (efu2,110) value,line
  110   FORMAT (' ',f15.2,5x,102A1)
  120 END DO

      WRITE (efu2,130) start,t1,t2,t3,stop,nwx,nwy
  130 FORMAT (14x,f10.2,4(15x,f10.2)/1x,A/55x,A)
      GO TO 20
  140 CONTINUE
      CLOSE(UNIT=efu1)

  RETURN
END Subroutine Xplot   ! -------------------------------------------------------

END Module AircraftMotionsProcedures




!+
PROGRAM AircraftMotions
! ---------------------------------------------------------------------------
! PURPOSE:
USE AircraftMotionsProcedures

!IMPLICIT NONE
  INTEGER,PARAMETER:: IN=1, OUT=2, DBG=3
  INTEGER,PARAMETER:: NMAX=300
  REAL,DIMENSION(30,NMAX):: save
  REAL,DIMENSION(NMAX):: time, xraw,yraw, hraw,draw,vraw,araw
  INTEGER,DIMENSION(NMAX):: hour,minute,second
  REAL,DIMENSION(20):: atab,vtab,dtab,ttab
  REAL,DIMENSION(60):: b,xsum,ysum,hsum
!!!      COMMON /dsave/ save(30,300),nb,ne,dt,time(300)
!!!      DIMENSION atab(20), vtab(20), dtab(20), ttab(20), b(60)
!!!      DIMENSION xraw(300), yraw(300), hraw(300), draw(300), vraw(300)
  CHARACTER(LEN=80):: fileName
!!!      DIMENSION xsum(60), ysum(60), hsum(60), araw(300)
  COMMON /DSAVE/ save,nb,ne,dt,time
!!!      DATA dot/'*'/,blank/' '/

  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(A,F10.2)'
  CHARACTER(LEN=*),PARAMETER:: FMT2 = '(A,F10.4)'
  CHARACTER(LEN=*),PARAMETER:: FMT3 = '(A,I4)'
  CHARACTER(LEN=*),PARAMETER:: FMT20 = '(A,F10.2)'
!  INTEGER:: hour,minute,second
  INTEGER:: nd
  INTEGER:: nh
  INTEGER:: nr
  INTEGER:: nv
  CHARACTER(LEN=20):: aircraftID
  CHARACTER(LEN=1),PARAMETER:: DOT='*', BLANK=' '
  CHARACTER(LEN=1):: p   ! gets set to DOT or BLANK
  INTEGER:: errCode, monoCode
!----------------------------------------------------------------------------

  DO   ! read in constants
    WRITE(*,*) 'Enter name of the case file:'
    READ(*,'(A)') fileName
    IF (LEN_TRIM(fileName)==0) STOP
    OPEN(UNIT=IN,FILE=fileName,STATUS='OLD', &
      IOSTAT=errCode,ACTION='READ',POSITION='REWIND')
    IF (errCode==0) EXIT
    OPEN(UNIT=IN,FILE=TRIM(fileName)//'.inp',STATUS='OLD', &
      IOSTAT=errCode,ACTION='READ',POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) 'Unable to open this file. Try again.'
  END DO
  INQUIRE(UNIT=IN,NAME=fileName)
  WRITE(*,*) 'Reading from ', Trim(fileName)
  READ(IN,'(A)') aircraftID
  READ(IN,*) dt,ws,clao,alphao,var
  READ(IN,*) nr,nv,nd,nh
  CLOSE(UNIT=IN)

  OPEN(UNIT=OUT,FILE='atc.out',STATUS='REPLACE',ACTION='WRITE')
  WRITE(OUT,*) 'PROGRAM TO PROCESS ATC RADAR DATA'
  WRITE(OUT,*) TRIM(aircraftID)
  WRITE(OUT,FMT1) 'time difference, seconds', dt
  WRITE(OUT,FMT1) 'wing loading (weight/wingArea), psf', ws
  WRITE(OUT,FMT2) 'lift curve slope, per degree', clao
  WRITE(OUT,FMT1) 'angle of attack for zero lift', alphao
  WRITE(OUT,FMT1) 'magnetic variation, degrees', var
  WRITE(OUT,FMT3) 'points for pre-smoothing', nr
  WRITE(OUT,FMT3) 'points for ground speed smoothing', nv
  WRITE(OUT,FMT3) 'points for track angle smoothing', nd
  WRITE(OUT,FMT3) 'points for altitude smoothing', nh
!!  WRITE (6,150) dt,ws,clao,alphao,var,nr,nv,nd,nh
  WRITE (OUT,190)
!!!  dd=0.0

  DO   ! read in winds and temperature .........................................
    WRITE(*,*) 'Enter name of the wind and temperature file:'
    READ(*,'(A)') fileName
    IF (LEN_TRIM(fileName)==0) STOP
    OPEN(UNIT=IN,FILE=fileName,STATUS='OLD', &
      IOSTAT=errCode,ACTION='READ',POSITION='REWIND')
    IF (errCode==0) EXIT
    OPEN(UNIT=IN,FILE=TRIM(fileName)//'.inp',STATUS='OLD', &
      IOSTAT=errCode,ACTION='READ',POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) 'Unable to open this file. Try again.'
  END DO
  INQUIRE(UNIT=IN,NAME=fileName)
  WRITE(*,*) 'Reading from ', Trim(fileName)

   20 FORMAT (f10.3,3f10.2)

  DO nw=1,SIZE(atab)
    READ(IN,*,IOSTAT=errCode) h,d,v,t
    IF (errCode < 0) EXIT
    WRITE (OUT,20) h,d,v,t
    atab(nw)=h*1000.
    vtab(nw)=v*1.69
    ttab(nw)=t
    dtab(nw)=d
  END DO
  CLOSE(UNIT=IN)

  DO   ! read in radar data.....................................................
    WRITE(*,*) 'Enter name of the radar data file:'
    READ(*,'(A)') fileName
    IF (LEN_TRIM(fileName)==0) STOP
    OPEN(UNIT=IN,FILE=fileName,STATUS='OLD', &
      IOSTAT=errCode,ACTION='READ',POSITION='REWIND')
    IF (errCode==0) EXIT
    OPEN(UNIT=IN,FILE=TRIM(fileName)//'.inp',STATUS='OLD', &
      IOSTAT=errCode,ACTION='READ',POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) 'Unable to open this file. Try again.'
  END DO
  INQUIRE(UNIT=IN,NAME=fileName)
  WRITE(*,*) 'Reading from ', Trim(fileName)


  WRITE(OUT,*) "                     PRINTOUT OF INPUT DATA"
  WRITE(OUT,*) " POINT                      EAST     NORTH"
  WRITE(OUT,*) " NO   HOUR  MIN   SEC      RANGE     RANGE     ALTITUDE"
  WRITE(OUT,*) "                           N. MI     N. MI       FEET"

!!!   50 FORMAT(2I5,3F10.4,F10.1)
  ne=0
  DO j=1,SIZE(time)   ! r
    READ(IN,*,IOSTAT=errCode) hour(j),minute(j),second(j),y,x,h
    IF (errCode < 0) EXIT
    ne=ne+1   ! counts the lines of good data
    WRITE (OUT,'(4I5,3F12.2)') j,hour(j),minute(j),second(j),y,x,h
!!!  160 FORMAT(3I5,F9.3,2F10.3,F10.0)
    time(j)=REAL(hour(j)*3600 + minute(j)*60 + second(j))
    xraw(j)=x*6080.3
    yraw(j)=y*6080.3
    hraw(j)=h
    araw(j)=0.0
    vraw(j)=0.0
    draw(j)=0.0
  END DO
  CLOSE(UNIT=IN)

  nb=1

!....................... PRE-SMOOTHING OF RAW DATA  ....................
  mmin=(nr-1)/2+1
  mmax=ne-mmin+1
  DO j=mmin,mmax
    DO l=1,nr
      ll=j+l-1-(nr-1)/2
      xsum(l)=xraw(ll)
      ysum(l)=yraw(ll)
      hsum(l)=hraw(ll)
    END DO
    CALL Smooth (xsum, dt, nr, xsmoth, xrate, xaccl)
    CALL Smooth (ysum, dt, nr, ysmoth, yrate, yaccl)
    CALL Smooth (hsum, dt, nr, hsmoth, hrate, haccl)
    vraw(j)=SQRT(xrate**2+yrate**2)
    draw(j)=ATAN2(yrate,xrate)
    araw(j)=hrate
  END DO

!!!  WRITE (OUT,210)
  WRITE(OUT,*) "                    PRINTOUT  OF  OUTPUT  DATA"
  WRITE(OUT,*)
  WRITE(OUT,*) " POINT             GROUND  TRACK  VERT.  FLIGHT     FORCES     ......ANGLES......... AIRSPEED "
  WRITE(OUT,*) " NO          ALT.   SPEED  ANGLE   VEL.   PATH   LIFT    T-D   ROLL  PITCH HEADING  TRUE   IND.  MACH"
  WRITE(OUT,*) "     TIME    FEET   KNOTS   DEG    FPS     DEG    G,S    G'S   DEG    DEG    DEG   KNOTS  KNOTS"

  DO 140 j=nb,ne
    p=blank
!....................... FINAL SMOOTHING FOR INERTIAL DATA  ............
    DO l=1,nd
      ll=j+l-1-(nd-1)/2
      IF (ll.lt.mmin) p=dot
      IF (ll.gt.mmax) p=dot
      IF (ll.lt.mmin) ll=mmin
      IF (ll.gt.mmax) ll=mmax
      b(l)=draw(ll)
      IF (l.eq.1) GO TO 100
      bx=b(l)-b(l-1)
      IF (bx.gt.3.141) b(l)=b(l)-6.283185
      IF (bx.lt.-3.141) b(l)=b(l)+6.283185
  100     CONTINUE
      ysum(l)=b(l)
    END DO 

    DO l=1,nv
      ll=j+l-1-(nv-1)/2
      IF (ll.lt.mmin) p=dot
      IF (ll.gt.mmax) p=dot
      IF (ll.gt.mmax) ll=mmax
      IF (ll.lt.mmin) ll=mmin
      xsum(l)=vraw(ll)
    END DO

    DO l=1,nh
      ll=j+l-1-(nh-1)/2
      IF (ll.lt.mmin) p=dot
      IF (ll.gt.mmax) p=dot
      IF (ll.gt.mmax) ll=mmax
      IF (ll.lt.mmin) ll=mmin
      hsum(l)=araw(ll)
    END DO

    CALL Smooth (ysum, dt, nd, drate, daccl, dxx)
    CALL Smooth (xsum, dt, nv, vrate, vaccl, vxx)
    CALL Smooth (hsum, dt, nh, hrate, haccl, hxx)
    xrate=vrate*cos(drate)
    yrate=vrate*sin(drate)
    xaccl=vaccl*cos(drate)-daccl*sin(drate)*vrate
    yaccl=vaccl*sin(drate)+daccl*cos(drate)*vrate

!....................... AIR DATA  .....................................
    hw=hraw(j)
    CALL Taint (atab, vtab, hw, vw, nw, 1, ner, monoCode)
    CALL Taint (atab, dtab, hw, dw, nw, 1, ner, monoCode)
    CALL Taint (atab, ttab, hw, tp, nw, 1, ner, monoCode)
    delta=(1.-(0.003566*hw)/518.688)**5.2561
    thbar=(tp+459.4)/518.688
    f=0.457072*(10.**(-6.))
    xwind=vw*cos(dw/57.3)
    ywind=vw*sin(dw/57.3)
    vx=xrate+xwind
    vy=yrate+ywind
    vh=hrate
    vv=SQRT(vx**2+vy**2)
    vt=SQRT(vx**2+vy**2+vh**2)
    v1=vt/1.69
    vc=(((((((((v1*v1)*(f/thbar)+1.)**3.5)-1.)*delta)+1.)**0.285714)-1.)/f)**0.5
    q=0.00119*vt**2*(1.0-(hw*.000001)*6.87535)**4.2561
    xm=vt/(49.04*sqrt(tp+459.4))
    cla=clao*(3.8/(1.8+2.0*((1.0-xm**2)**0.5)))

!....................... FORCES  .......................................
    at=SQRT(yaccl**2+xaccl**2+(haccl+32.2)**2)/32.2
    av=haccl/32.2+1.0
    ad=(yaccl*vy+xaccl*vx+(haccl+32.2)*vh)/(32.2*vt)
    al=SQRT(at**2-ad**2)

!....................... ANGLES  .......................................
    roll=-57.3*asin((xaccl*vy-yaccl*vx)/(vv*at*32.2))
    alpha=alphao+(al*ws)/(q*cla)
    theta=ATAN2(vh,vv)*57.3+alpha*COS(roll/57.3)
    hdg=ATAN2(vy,vx)*57.3+alpha*SIN(roll/57.3)+var

!....................... SAVE DATA FOR PRINT AND PLOTS  ........ .......
        SAVE  (1,j)=hraw(j)
        SAVE  (2,j)=vrate/1.69
        SAVE  (3,j)=drate*57.3
        SAVE  (4,j)=hrate
        SAVE  (5,j)=atan2(hrate,vrate)*57.3
        SAVE  (6,j)=al
        SAVE  (7,j)=ad
        SAVE  (8,j)=roll
        SAVE  (9,j)=theta
        SAVE  (10,j)=hdg
        SAVE  (11,j)=vt/1.69
        SAVE  (12,j)=vc
        SAVE  (13,j)=xm
        SAVE  (14,j)=araw(j)
        SAVE  (15,j)=vraw(j)/1.69
        SAVE  (16,j)=draw(j)*57.3
        SAVE  (17,j)=atan2(vh,vv)*57.3
        SAVE  (18,j)=av
        SAVE  (19,j)=at
        SAVE  (20,j)=alpha
        SAVE  (21,j)=q
        SAVE  (22,j)=xraw(j)/6080.3
        SAVE  (23,j)=yraw(j)/6080.3
        SAVE  (24,j)=dw
        SAVE  (25,j)=vw/1.69
        SAVE  (26,j)=time(j)
        SAVE  (27,j)=(j-1)*dt
!!!        hour=INT(time(j)/3600.)
!!!        minute=INT((time(j)-hour*3600.)/60.)
!!!        second=INT(time(j)-hour*3600.-minute*60.)
        WRITE(OUT,170) p,j,hour(j),minute(j),second(j),INT(save(1,j)), (save(n,j),n=2,13)
  140 END DO



      WRITE (OUT,220)
      CALL Xplot(IN,OUT)
!!!  150 FORMAT(A20,2F10.2,F10.4,2F10.2,4I10)

  170 FORMAT(1X,A1,2I3,':',I2.2,':',I2.2,I6,F6.1,11F7.2)
!  180 FORMAT('1',' PROGRAM TO PROCESS ATC RADAR DATA                '// &
!     &//23x,'.................CONSTANTS.......................   ',     &
!     & '.....NO. OF POINTS FOR SMOOTHING....... ',//                    &
!     & ' AIRCRAFT  ID       ',3x,                                       &
!     & '  TIME     WEIGHT/    LIFT      ALPHA      MAG.   ','    PRE-    &
!     & GROUND     TRACK     ALTITUDE                  '/23x,'  DIF.     &
!     & AREA      CURVE     ZERO       VAR.   ','   SMOOTH    SPEED      &
!     &ANGLE                               '/23x,'  SEC      LB/SQFT   PE&
!     &R DEG    DEG        DEG.   '/)
  190 FORMAT (/////4x,'ALTITUDE        WINDS          TEMPERATURE'/4x,' &
     &1000 FT    DEG      KNOTS     DEG F      '/)
  220 FORMAT (///' * SMOOTHED VALUES ARE APPROXIMATE NEAR END POINTS')

  WRITE(*,*) "Normal termination of program AircraftMotions"
  STOP

END Program AircraftMotions   ! ================================================
