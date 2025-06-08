INCLUDE 'g.for'
!+
PROGRAM Problem2
! ---------------------------------------------------------------------------
! PURPOSE - Calculate the end point (exhaust conditions) of fluid
!   propertiesfor a methane expander with inlet conditions of
!   100 atmospheres and 200 degrees K and exhaust conditions at the
!   critical pressure.   Assume the expansion is isentropic.
    
    
!     DO NOT USE K AS AN INDEX IT MEANS THERMAL CONDUCTIVITY
!     DO NOT FORGET YOUR COMMON STATEMENT
!     DO NOT FORGET TO CALL SETUP

       REAL mu,mul,muv,k,kl,kv
       COMMON/PROPTY/ku,dl,dv,hl,hv,s,sl,sv,cv,cvl,cvv,cp,cpl,cpv,gamma, &
         gammal,gammav,c,cl,cvp,mu,mul,muv,k,kl,kv,sigma,excesk,excl,excv
    
!!!         DATA PIN,TIN/100.,200./    

      REAL,PARAMETER:: PIN = 100.0
      REAL,PARAMETER:: TIN=200.0
      CHARACTER(LEN=3),PARAMETER:: nam = "CH4"
!----------------------------------------------------------------------------
      CALL Setup(nam)

      ku=2   !   select units of calculation - set ku=2
      kr=0  ! assume region of calculation unknown. set kr=0
      ks=1
      kp=31  ! everything but surface tension

      CALL Gasp(ks,kp,tin,pin,din,hin,kr)
      sin=s
      cpin=cp
      cvin=cv
      sonin=c
      visin=mu
      conin=k
!!!         uin=hin-pin*0.101325/din
      akrin=kr

      WRITE(6,*) "tin ", tin
      WRITE(6,*) "pin ", pin
      WRITE(6,*) "din ", din
      WRITE(6,*) "hin ", hin
      WRITE(6,*) "sin ", sin
      WRITE(6,*) "cpin ", cpin
      WRITE(6,*) "cvin ", cvin
      WRITE(6,*) "sonin ", sonin
      WRITE(6,*) "visin ", visin
      WRITE(6,*) "conin ", conin
      WRITE(6,*) "akrin ", akrin

!   isentropic expansion to outlet  sout=sin; p=pout; KS=5
!   leave S unchanged
      ks=5
      p=45.66   ! critical pressure of methane (atm.)

      CALL Gasp(ks,kp,tout,p,dout,hout,kr)
      sout = s
      cpout=cp
      cvout=cv
      sonout=c
      visout=mu
      conout=k
!!!         uout=hout-p*0.101325/din
      akrout=kr

      WRITE(6,*) "tout ", tout
      WRITE(6,*) "p ", p
      WRITE(6,*) "dout ", dout
      WRITE(6,*) "hout ", hout
      WRITE(6,*) "sout ", sout
      WRITE(6,*) "cpout ", cpout
      WRITE(6,*) "cvout ", cvout
      WRITE(6,*) "sonout ", sonout
      WRITE(6,*) "visout ", visout
      WRITE(6,*) "conout ", conout
      WRITE(6,*) "akrout ", akrout
  WRITE(*,*) "End of Problem 2"
  STOP 
END Program Problem2   ! ====================================================

