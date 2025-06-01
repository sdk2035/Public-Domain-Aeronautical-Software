      PROGRAM aipp


!     ******************************************************************
!     **  ATMOSPHERIC INTERACTION PLUME PROGRAM        5/1/74         **
!     **          AEROCHEM RESEARCH LABORATORIES, INC.                **
!     **                    PRINCETON, NEW JERSEY   O8540             **
!     ******************************************************************
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      DIMENSION id(20)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      REAL massmd,l1,l2,lbar
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      COMMON /mainsh/ kount,kint
      COMMON /ibugs1/ ibugsh
      COMMON /chem1/ zid(5),irr(30),irt(30),rc(30,3),irrr(30,5),av,     &
     &cm(21,21),fi(21),wp(21),wm(21),wdot(21,30),csi(21),qx(21)
      COMMON /intl/ xlx,iki,rex
      COMMON /inpux/ itd,x0,y0,frac,kfl
      COMMON /shoput/ xshw,rshw,psiw,pshw,tshw,ushw,cshw(20)
      COMMON /inshop/ poo,too,uoo,phoo,coo(20)
      COMMON /xlm/ xlmd
      COMMON /dmd/ kmnf
      COMMON /sub/ xdis,l1,l2,phio,rstar,rmd,massmd,rnoz,thtp,thru
      COMMON /dsl/ rdsl,phdsl,xinit
      COMMON /tripleBlock/ del1,pp1,tt1,uu1,th3,pp3,tt3,uu3,slang,shang
      COMMON /turple/ del2,cpshw,gamshw,zmwsh1
      COMMON /utpp/ psh11,tsh11,ush11,flang,csh11(20)
      INTEGER thru
      DIMENSION error(3)
      EQUIVALENCE (move,iextra(15))
      EQUIVALENCE (kalarm,iextra(20))
      DATA error/'MAIN','BDRY','STEP'/
!     SETTING UP THE PROBLEM
!     REVISED SUBSCRIPTS (ZEROS DELETED)
      WRITE (6,10)

   10 FORMAT ('1',//////////40x,'ATMOSPHERIC INTERACTION PLUME PROGRAM',&
     &///40x,'AEROCHEM RESEARCH LABORATORIES, INC.',/45x,'PRINCETON, NEW&
     & JERSEY 08540',/////45x,'REVISED - JUNE 7,1974')
      idata=0
      jtest=0
      READ (5,70) ndata,egsf1,egsf2
   20 READ (5,30) (id(i),i=1,20)
      WRITE (6,40) (id(i),i=1,20)
   30 FORMAT (20a4)
   40 FORMAT ('1',20a4)
      WRITE (6,50) ndata,egsf1,egsf2
   50 FORMAT ('0',/10x,'NO. OF DATA SETS ',i5,6x,'EGSF1=',f6.2,6x,'EGSF2&
     &=',f6.2)
      DO 60 i=1,7300
   60   a(i,1)=0.0
      betap=1.013e06
      hj=4.182e07
      rcon=1.9871
      kalarm=0
      kount=0
      kint=0
      ll=-1
      READ (5,70) nsec,dxlss,xlmax,xlmdd,imd,nnt
      WRITE (6,80) nsec,dxlss,xlmax,xlmdd,imd,nnt
   70 FORMAT (i10,3e12.5,2i2)
   80 FORMAT (' ',4x,'NSEC=',i5,4x,'DXLSS=',1pe10.3,4x,'XLMAX=',e10.3,  &
     &4x,'XLMDD=',e10.3,4x,'IMD=',i5,4x,'NNT=',i5///)
      xlmmx=dxlss
      fsec=nsec
      told=0.
      xlmd=1.e10
      IF (isb.ne.-1) GO TO 90
      END FILE 12
   90 CONTINUE
      CALL putin
      IF (ipart.ne.0) CALL pputin
      IF (ikine.ne.0) CALL cputin (ident, ikine, nds)
      IF (kmn.gt.1) GO TO 220
      IF (ikind.lt.3.or.iki.ne.0) GO TO 170
      cpb=0.0
      zmb=0.0
      DO 100 i=1,nds
        zmb=zmb+cstrem(i)/mw(i)
  100   cpb=cpb+cstrem(i)*a(i,3)
      zmb=1.0/zmb
      emf=ustrem*sqrt((cpb*zmb-rcon)/(cpb*tstrem*ru))
      pe=0.0
      em=0.0
      DO 110 i=2,kmax
        em=em+ma(i)
  110   pe=pe+p(i)
      em=em/(kmax-1)
      pe=pe/(kmax-1)
      IF (ipart.eq.0) GO TO 130
      rhg=0.0
      pgm=0.0
      DO 120 i=2,kmax
        rhg=rhg+rho(i)
        DO 120 j=1,npg
        n=i-1
  120   pgm=pgm+rhp2(n,j)
      pgm=pgm/rhg
      GO TO 140
  130 pgm=0.0
  140 xlmd=rex*1.38*em*sqrt(gama*pe/pstrem)/(emf**0.25)
!     CARD OMITTED
  150 FORMAT (e12.5)
      IF (ikind.lt.3) xlmd=1.e10
      IF (imd.ne.0) xlmd=xlmdd
      WRITE (9,150) xlmd
      WRITE (6,160) xlmd
  160 FORMAT (//3x,'THE MACH DISC IS LOCATED AT X=',e11.3,1x,'CM'//)
      GO TO 210
  170 IF (ishock.ne.1) GO TO 220
      penf=pstrem
      gamm=cpb/(cpb-rcon/zmb)
      emex=ma(kmax)
      semx2=sqrt(emex**2-1.0)
      emf2=emf**2
      thn=phi(kmax)
      pex=p(kmax)
      emb=emex*(pex/penf)**((gama-1.)/gama*0.5)
      tl=1.e-4*emb
      c1=sqrt((gama+1.)/(gama-1.))
      c3=gama/(gama-1.)
      thc=c1*atan(semx2/c1)-atan(semx2)-thn
      DO 180 it=1,100
        emb2=emb**2
        sqemb=sqrt(emb2-1.0)
        zta=c1*atan(sqemb/c1)-atan(sqemb)-thc
        pp=pex*(1.+(gama-1.)/2.*emex**2)**c3/((1.+(gama-1.)/2.*emb2)**  &
     &   c3)
        bta=asin(sqrt(1.+(gamm+1.)/(2.*gamm)*(pp/penf-1.))/emf)
        sinbta=sin(bta)
        snbta2=sinbta**2
        cosbta=cos(bta)
        e=((emf*sinbta)**2-1.)/tan(bta)/((gamm+1.)/2.*emf2-(emf*sinbta)*&
     &   *2+1.)-tan(zta)
        dnom=(gamm+1.)/2.*emf2-(emf*sinbta)**2+1.0
        dzmb=emb/(sqemb*(1.+(emb2-1.)/c1**2))-1./(emb*sqemb)
        dpb=-pp*emb/(1+(gama-1.)/2.*emb2)*gama
        dbta=(gamm+1.)/(4.*gamm*emf2*sinbta*cosbta)*dpb
        dhdmb=-dzmb/((cos(zta))**2)+((-((emf*sinbta)**2-1.)*dbta/       &
     &   (snbta2)+2.*(emf*cosbta)**2*dbta)/dnom+((emf*sinbta)**2-1.)*2.*&
     &   (emf*cosbta)**2*dbta/dnom**2)
        IF (abs(e).lt.tl) GO TO 200
        emb=emb-e/dhdmb
  180 END DO
      WRITE (6,190) it,e,dhdmb,emb,zta,bta,pp
  190 FORMAT (8x,'ITERATION DOES NOT CONVERGE',i4/2x,6e12.5)
  200 CONTINUE
      psi=bta
      psh=pp
      rsh=r(kmax)
      xsh=x(kmax)
      gm1=(gamm+1.)/(gamm-1.)
      tsh=tstrem*(pp/penf+gm1)/(gm1+penf/pp)
      un=(pp/penf+gm1)/(gm1*pp/penf+1.)*ustrem*sinbta
      ush=sqrt((ustrem*cosbta)**2+un**2)
  210 CONTINUE
  220 kmx2=kmx+2
      icont=0
      intger=alphap
      nucase=0
      ierror=0
      ll=0
      s=0.0
      xlx=1.0e5
      spdy=0.0
      ikk=0
      lint=0
      sled=(1.0e-8)*rex
!     ESTABLISHING STABLE STEP SIZE
  230 CALL second (time)
      tleft=fsec-time
      tstep=time-told
      IF (ll.le.0) GO TO 240
      IF (tstep.gt.(0.5*tleft)) GO TO 420
  240 told=time
      CALL stable
      icount=0
      move=0
      IF ((ll.gt.lplane).or.(x(2).gt.xlmax)) GO TO 720
      IF ((x2(1).ge.xlmd).and.(ishock.ne.4).and.(kmn.eq.0)              &
     &.and.(isb.ne.1)) GO TO 470
      IF ((x2(1).lt.xlmd).or.(kmn.ge.1).or.(isb.eq.1)) GO TO 750
      WRITE (6,250) xshw,xlmd
  250 FORMAT (/5x,e10.3,2x,e10.3)
      WRITE (6,260) del2,p2(2),t2(2),u2(2),rho2(2),zmwk(2)
  260 FORMAT (/5x,6e10.3)
!
!.....COMPUTING PROPERTIES IN SUBSONIC AND SUPERSONIC REGIONS
!     DOWNSTREAM OF TRIPLE POINT AND WRITING ON TAPE4 FOR INPUT IN
!     SS REGION CALCULATIONS.
      iconst=1
      CALL triple
      WRITE (4,270) (csh11(i),i=1,nds)
      WRITE (4,280) x2(1),r2(1),p2(2),t2(2),u2(2),phi2(2)
      WRITE (4,270) (c2(i,2),i=1,nds)
      WRITE (4,290) flang,slang,pp1,tt1,uu1,shang,pp3,tt3,uu3
  270 FORMAT (e10.3)
  280 FORMAT (6e12.5)
  290 FORMAT (9e12.5)
      itype=1
      ishock=1
      mmax=20
      xst=2.*y(kmax)/tan(asin(1./ma(2)))
      DO 300 m=2,19
        xw(m)=x2(1)+abs(float(1-m)/float(20-m)*r2(1)/tan(phi(1)))
  300   rw(m)=r2(1)*(float(20-m))/19.
      xw(1)=x2(1)
      rw(1)=r2(1)
      xw(20)=1.e8
      rw(20)=0.0
      WRITE (6,310)
  310 FORMAT (/10x,'###FICTITIOUS WALL COORDINATES FOR CONTINUING SHOCK &
     &LAYER SOLUTION DOWNSTREAM OF MACH DISC###'/)
      WRITE (6,400)
      sw(1)=s
      swa=sw(1)
      phiw(1)=atan((rw(2)-rw(1))/(xw(2)-xw(1)))
      m=1
  320 m=m+1
      IF (m-mmax) 330,340,340
  330 phiw(m)=atan((rw(m+1)-rw(m-1))/(xw(m+1)-xw(m-1)))
      GO TO 350
  340 phiw(m)=2.0*atan((rw(m)-rw(m-1))/(xw(m)-xw(m-1)))
      phiw(m)=phiw(m)-phiw(m-1)
  350 fkw=(sin(phiw(m))-sin(phiw(m-1)))/(xw(m)-xw(m-1))
      IF (abs(fkw/(xw(m)-xw(m-1)))-1.e-06) 370,370,360
  360 delswa=(phiw(m)-phiw(m-1))/fkw
      GO TO 380
  370 delswa=sqrt((xw(m)-xw(m-1))**2+(rw(m)-rw(m-1))**2)
  380 swa=swa+delswa
      sw(m)=swa
      WRITE (6,390) m,xw(m),rw(m),phiw(m),pw(m),sw(m)
      IF (m-mmax) 320,410,410
  390 FORMAT (' ',1x,i2,1x,6(1pe10.3,5x))
  400 FORMAT (' ',/2x,'M',2x,'XW',13x,'RW',13x,'PHIW',11x,'PW',13x,'SW')
  410 xlmd=1.e20
      iki=0
      xlx=x2(2)
      CALL stable
      GO TO 750
  420 k=1
      IF (iflag.eq.0) GO TO 720
      WRITE (7,430) x(k),r(k),phi(k)
  430 FORMAT (6e12.5)
      DO 450 k=2,kmax
        WRITE (7,430) x(k),r(k),phi(k),p(k),t(k),u(k)
        WRITE (7,440) (c(i,k),i=1,nds)
  440   FORMAT (8e10.3)
  450 END DO
      x11=1.e8
      x22=1.e-4
      WRITE (7,460) x(1),r(1),x11,r(1),x(kmax),r(kmax),x22,r(kmax)
  460 FORMAT (2e12.5/2e12.5/2e12.5/2e12.5)
      IF (ikind.lt.4) GO TO 720
      WRITE (7,430) x(kmax),r(kmax),psi,psh,tsh,ush
      WRITE (7,440) (cstrem(i),i=1,nds)
      WRITE (7,430) xstrem,rstrem,phstrm,pstrem,tstrem,ustrem
      WRITE (7,440) (cstrem(i),i=1,nds)
      IF (itype.ne.3) GO TO 720
      WRITE (7,430) xshw,rshw,psiw,pshw,tshw,ushw
      WRITE (7,440) (cshw(i),i=1,nds)
      WRITE (7,430) xoo,roo,phoo,poo,too,uoo
      WRITE (7,440) (coo(i),i=1,nds)
      GO TO 720
!     ..... BEHIND THE MACH DISC CONDITIONS .....
  470 IF (iki.eq.0) GO TO 740
      WRITE (6,490)
      DO 480 k=2,kmax
  480   WRITE (6,500) x2(k),r2(k),rho2(k),p2(k),t2(k),u2(k),phi2(k),    &
     &   ma(k),gama,zmwk(k)
  490 FORMAT (//5x,'PROPERTIES UPSTREAM OF THE MACH DISC ARE '//5x,'X(CM&
     &)',7x,'R(CM)',7x,'RHO(G/CM3)',2x,'P(ATM)     T(K)        U(CM/SEC &
     &  PHI(RAD)    MA          GAMMA     MOL.WT.')
  500 FORMAT (5x,10(1pe12.5))
      ka=k
      DO 630 k=2,kmax
        tl=t2(k)*1.e-4
        rhoxy=.15
!     CARD OMITTED
        rhoy=rho2(k)/rhoxy
        DO 540 it=1,100
          py=p2(k)+rho2(k)*u2(k)**2*(1.-rhoxy)/betap
          hy=h(k)+u2(k)**2/(2.*hj)*(1.-rhoxy**2)
          ty=t2(k)*py/p2(k)*rhoxy
          DO 520 itt=1,100
            dhyt=0.0
            hyt=0.0
            DO 510 i=1,nds
              DO 510 jj=1,nn
              hyt=hyt+a(i,jj)*ty**(jj-2)*c2(i,k)
  510         dhyt=dhyt+a(i,jj)*float(jj-2)*ty**(jj-3)*c2(i,k)
            d=hyt-hy
            IF (abs(d).lt.tl) GO TO 530
            ty=ty-d/dhyt
  520     CONTINUE
          WRITE (6,550) it,k,hy,hyt,dhyt
          STOP
  530     rhoy1=py*zmwk(k)/rv/ty
          IF ((abs(rhoy1/rhoy)-1.).lt..001) GO TO 560
          rhoy=(rhoy1+rhoy)*.5
          rhoxy=rho2(k)/rhoy
  540   CONTINUE
        WRITE (6,550) it,k,rhoy,rhoy1,rhoxy
  550   FORMAT (2x,'ITERATIONS DO NOT CONVERGE',/2x,'IT=',i3,2x,'K=',i3,&
     &   2x,'HY=',e12.5,2x,'H(K)=',e12.5,2x,'DHYT=',e12.5)
  560   uy=u2(k)*rho2(k)/rhoy1
!     CARD OMITTED
!     CARD OMITTED
!     CARD OMITTED
        IF (k.gt.2) GO TO 580
        WRITE (6,570)
  570   FORMAT (//10x,'THE MACH DISC HAS BEEN REACHED- POST-NORMAL SHOCK&
     &   PROPERTIES ALONG NEXT DOWNSTREAM SURFACE ARE-'//,2x,'X(CM)',8x,&
     &   'R(CM)',2x,'RHO(GM/CM3)',3x,'P(ATM)',8x,'T(K)',4x,'U(CM/SEC)'/)
  580   WRITE (6,590) x2(k),r2(k),rhoy1,py,ty,uy
  590   FORMAT (6e12.5)
        IF (k.gt.2) GO TO 630
        WRITE (6,600)
  600   FORMAT (//10x,'SPECIES MASS FRACTIONS ON THE AXIS AT THE MACH DI&
     &SC ARE'/)
        DO 610 i=1,nds
  610     WRITE (6,620) i,c2(i,2)
  620   FORMAT (10x,'C(',i2,')=',e14.4)
  630 END DO
      IF (ipart.eq.0) GO TO 710
      WRITE (6,640)
  640 FORMAT (//15x,'PARTICLE PROPERTIES AT THE MACH DISC'/)
      DO 650 k=1,kmax
        WRITE (6,660) k
        WRITE (6,670) (w(k,j),j=1,npg)
        WRITE (6,680) (v(k,j),j=1,npg)
        WRITE (6,690) (tp(k,j),j=1,npg)
        WRITE (6,700) (rhp(k,j),j=1,npg)
  650 END DO
  660 FORMAT (//' ','K=',i5/)
  670 FORMAT (' ','DOWNSTREAM VELOCITY',2x,1p8e12.3)
  680 FORMAT (' ','CROSS-STREAM VELOCITY',1p8e12.3)
  690 FORMAT (' ','PARTICLE TEMPERATURE',1x,1p8e12.3)
  700 FORMAT (' ','PARTICLE DENSITY',5x,1p8e12.3)
  710 k=ka
  720 CONTINUE
      IF (x2(1).ge.xlmd.and.ishock.ne.4) GO TO 730
      idata=idata+1
      IF (idata.lt.ndata) GO TO 20
      END FILE 13
           !callexit
      STOP
  730 lplane=ll+30
  740 xlmd=1.e20
  750 iflag=0
      IF (par(2).eq.0.0) GO TO 770
      ipar2=par(2)
      IF (mod(ll,ipar2).eq.0) WRITE (6,760) ll,x2(2),dels
  760 FORMAT (35x,'LL=',i5,4x,'X=',f8.1,4x,'STEP SIZE=',e11.4)
  770 CONTINUE
      IF (dels.lt.sled) GO TO 1290
      l=1
      IF (icount.eq.0) GO TO 780
      s=s-delso
      delso=0.5*delso
      IF (kmn.le.1) GO TO 800
      IF (karm.eq.11) delso=0.1*delso
      karm=0
      GO TO 800
  780 ksm=iextra(7)
      fcare=cos(phi(ksm)-phi(l))
      IF (fcare.gt.0.25) GO TO 790
      delso=dels*0.25
      GO TO 800
  790 delso=dels*fcare
!     INNER BOUDARY OR WALL CONDITIONS
  800 s=s+delso
      dls(1)=delso
  810 nbound=0
      l=1
      CALL bndry
      extra(1)=y(l)
      extra(2)=r(l)
      IF (mu.ne.0.0) um(k-1)=mu
      IF (k.ge.kmax) GO TO 1100
      IF (nucase.ne.0) GO TO 1240
!     CALCULATION OF STREAMLINE CURVATURE
  820 sumpv=0.0
      sump=0.0
      sumpe=0.0
      IF (((k-1).gt.kmx2).or.(ipart.eq.0)) GO TO 840
      CALL drag
      DO 830 j=1,npg
        sump=sump+xmomp(j)*rhp(k-1,j)
        sumpv=sumpv+xmomv(j)*rhp(k-1,j)
        sumpe=sumpe+rhp(k-1,j)*(enep1(j)+enep2(j))
  830 END DO
  840 IF ((kmn.gt.1).and.(k.le.kmn)) GO TO 990
      IF (iflag.ne.0) GO TO 860
  850 dpdy(k)=2.0*(p(k+1)-p(k))/(y(k+1)-y(k-1))
      dphids(k)=-2.0*(p(k+1)-p(k))/(rho(k+1)*u(k+1)**2*dely(k+1)+rho(k)*&
     &u(k)**2*dely(k))*betap
      dphids(k)=dphids(k)+sumpv*2.0/(rho(k+1)*u(k+1)**2+rho(k)*u(k)**2)
      GO TO 880
  860 dpdy(k)=2.0*(p2(k+1)-p2(k))/(y2(k+1)-y2(k-1))
      dphods=-2.0*(p2(k+1)-p2(k))/(rho2(k+1)*u2(k+1)**2*dely(k+1)+      &
     &rho2(k)*u2(k)**2*dely(k))*betap
      dphods=dphods+sumpv*2.0/(rho2(k+1)*u2(k+1)**2+rho2(k)*u2(k)**2)
      dphids(k)=.45*dphids(k)+.55*dphods
!     LOCATING NEW GRID POINTS
!     IK FIRST PT.,IKI ONE INITIAL POINT
      IF ((ikind.eq.2).or.(nnt.ne.0).or.(k.gt.kfl).or.(ikind.eq.1)      &
     &.or.(r(k).le.rex)) GO TO 880
      IF (ikk.ne.0) GO TO 870
      IF (dpdy(k).le.0.0) GO TO 880
      IF (ikk.eq.0) ik=k
      ikk=1
  870 dpy=p2(k+1)-p2(k)
      spdy=dpy+spdy
      p2p=spdy/p2(k+1)
      oat=1.-1./egsf1
      IF ((p2p.le.oat).or.(iki.ne.0)) GO TO 880
      xlx=x2(2)
      lint=lint+1
      IF (lint.eq.1) GO TO 880
      IF (kmn.le.1) CALL ental (ik)
      ikind=3
      ishock=0
      ustrem=1.e-6
      tstu=t2(kfl)
      IF (kmax.gt.kfl) tstu=t2(kfl+1)
      pstrem=p2(kfl)*(250./tstu)**(gama/(gama-1.))
      pstrem=egsf2*pstrem
      IF (kmn.gt.1) pstrem=1.0e-6
  880 phi1=phi(k)+dphids(k)*dls(k-1)
      theta=pi/2.0-0.5*(phi2(k-1)+phi1)
      ttheta=tan(theta)
      phibar=.5*(phi(k)+phi1)
      tphbr=tan(phibar)
      IF (abs(theta).lt.(pi/4.).or.abs(theta).gt.(3.*pi/4.)) GO TO 900
      IF (abs(theta-pi/2.).gt.1.e-04) GO TO 890
      x2(k)=x2(k-1)
      r2(k)=r(k)+(x2(k)-x(k))*tan(.5*(phi(k)+phi1))
      GO TO 920
  890 r2(k)=((x2(k-1)-x(k))*ttheta*tphbr+r(k)*ttheta+r2(k-1)*tphbr)/    &
     &(ttheta+tphbr)
      x2(k)=x2(k-1)-(r2(k)-r2(k-1))/ttheta
      GO TO 920
  900 IF (abs(theta).gt.1.e-04) GO TO 910
      r2(k)=r2(k-1)
      x2(k)=x(k)+(r2(k)-r(k))/tphbr
      GO TO 920
  910 x2(k)=(r2(k-1)-r(k)+x2(k-1)*ttheta+x(k)*tan(.5*(phi(k)+phi1)))/   &
     &(tan(.5*(phi(k)+phi1))+ttheta)
      r2(k)=r2(k-1)-(x2(k)-x2(k-1))*ttheta
  920 d=sqrt((x2(k)-x(k))**2+(r2(k)-r(k))**2)
      IF (abs(phi1-phi(k))-(1.0e-6)) 930,940,940
  930 dls(k)=d
      GO TO 950
  940 dls(k)=d*(phi1-phi(k))/(2.0*sin(.5*(phi1-phi(k))))
  950 phi2(k)=phi(k)+dls(k)*dphids(k)
!     CALCULATION OF STREAMTUBE PROPERTIES
      b=sqrt((x2(k)-x2(k-1))**2+(r2(k)-r2(k-1))**2)
      IF (abs(phi2(k)-phi2(k-1))-(10.**(-4))) 960,960,970
  960 dely(k)=b
      GO TO 980
  970 dely(k)=(phi2(k)-phi2(k-1))*b/(2.0*sin(.5*(phi2(k)-phi2(k-1))))
  980 aa(k)=pi*(r2(k)+r2(k-1))**delta*dely(k)
      y2(k)=y2(k-1)+dely(k)
      GO TO 1060
  990 tanph=tan(phio)
      xbar=x2(1)-xinit
      lbar=l2-l1
      thru=0
!
!     CALCULATING DIVIDING STREAMLINE RADIUS AND SLOPE FOR SS REGION
!
      IF (x2(1).gt.xinit+l1) GO TO 1000
      rdsl=rmd+tanph*xbar*(1.0-(xbar/(2.0*l1)))
      rdsl=par(6)*xbar**3+rdsl
      rmax=rmd+tanph*l1/2.0+par(6)*l1**3
      phdsl=atan(tanph*(1.0-(xbar/l1)))
      GO TO 1020
 1000 IF (x2(1).gt.xinit+l2) GO TO 1010
      rdsl=rmax+(rmax-rstar)*(xbar-l1)**2/lbar**2*(-3.0+2.0*(xbar-l1)/  &
     &lbar)
      rdsl=rdsl+par(7)*(xbar-l1)**4/lbar**4
      phdsl=atan((6.0*(rmax-rstar)*(xbar-l1)/lbar**2)*((xbar-l1)/lbar-  &
     &1.0))
      GO TO 1020
 1010 CONTINUE
      u(2)=u(2)*1.01/ma(2)
      ma(2)=1.01
      GO TO 850
      thru=1
!
!     CALCULATING NEW SUBSONIC STREAMTUBE COORDINATES AND AREA
!
 1020 x2(k)=x2(1)
      r2(k)=rdsl
      IF ((x2(1).lt.xinit).or.(x2(1).gt.xinit+l2).or.(par(8).eq.0.0)    &
     &.or.(iflag.eq.0)) GO TO 1040
      ipar8=par(8)
      IF (mod(ll,ipar8).eq.0) WRITE (6,1030) x2(1),rdsl
 1030 FORMAT (35x,'DSL COORDINATES ARE',2x,'X=',f8.1,4x,'R=',f8.1)
 1040 CONTINUE
      IF ((mod(ll,20).ne.0).or.(x2(1).gt.xinit+l2).or.(iflag.eq.0)) GO  &
     &TO 1050
      aaa=x2(1)/rnoz
      bbb=rdsl/rnoz
      WRITE (11) aaa,bbb
 1050 CONTINUE
      phi2(k)=phdsl
      dely(k)=r2(k)-r2(k-1)
      aa(k)=pi*(r2(k)+r2(k-1))**delta*dely(k)
      y2(k)=y2(k-1)+dely(k)
 1060 CALL step
      IF ((iflag.ne.0).and.((k-1).le.kmx2).and.(ipart.ne.0)) CALL part
      IF (kalarm.le.0) GO TO 1080
      move=0
 1070 kalarm=0
      IF (icount-25) 750,1230,1230
 1080 IF (nucase.ne.0) GO TO 1240
      IF (move.ne.0) GO TO 1100
      IF (k.gt.kmn) GO TO 1090
      dls(k)=dls(1)
      rbar(k)=r(k)
 1090 k=k+1
      IF (k-kmax) 820,1100,1100
!     OUTER BOUNDARY CONDITIONS
 1100 nbound=1
      ikk=0
      spdy=0.0
      CALL bndry
      IF (move.ne.0) GO TO 880
      move=0
      IF (kalarm) 1110,1110,1070
 1110 IF (nucase) 1240,1120,1240
 1120 y2(k)=y2(k-1)+dely(k)
      IF (iflag) 1130,1220,1130
 1130 max=kmax+1
      DO 1180 l=2,max
        j=l-1
        r(j)=r2(j)
        x(j)=x2(j)
        phi(j)=phi2(j)
        IF (j-1) 1180,1180,1140
 1140   u(j)=u2(j)
        p(j)=p2(j)
        t(j)=t2(j)
        rho(j)=rho2(j)
        zmw=0.0
        DO 1150 i=1,nds
 1150     zmw=zmw+c2(i,j)/mw(i)
        zmw=1.0/zmw
        zmwk(j)=zmw
        DO 1160 i=1,nds
          c(i,j)=c2(i,j)
 1160     xs(i,j)=c(i,j)*zmw/mw(i)
        h(j)=0.0
        DO 1170 i=1,nds
          DO 1170 jj=1,nn
 1170     h(j)=h(j)+a(i,jj)*t(j)**(jj-2)*c(i,j)
        y(j)=y(j-1)+dely(j)
        sx(j)=sx(j)+dls(j)
 1180 END DO
      IF (ikind.gt.3) CALL shocke
      icont=icont+1
      IF ((iflag.ne.0).and.(ipart.ne.0)) CALL pbndy
      kmx2=kmx+2
      ll=ll+1
      icount=0
      IF ((isb.ne.1).or.(x2(1).le.xinit+l2).or.(jtest.ne.0)) GO TO 1190
      jtest=1
      aaa=(xinit+l2)/rnoz
      bbb=rstar/rnoz
      WRITE (11) aaa,bbb
      END FILE 11
 1190 CONTINUE
      IF (mod(icont,intger)) 1200,1210,1200
 1200 IF (x(2).lt.xlmmx) GO TO 230
      xlmmx=xlmmx+dxlss
 1210 k=kmax
      CALL putout (1)
      IF ((par(1).eq.0.0).and.(iflag.ne.0).and.(ipart.ne.0)) CALL       &
     &pptout (1)
      GO TO 230
 1220 IF (ikind.ge.4) CALL shocke
      iflag=1
      GO TO 810
!     ERROR MESSAGES
 1230 nucase=876
      ierror=1
 1240 IF (ierror.eq.0) GO TO 1260
      WRITE (6,1310) error(ierror),nucase
      WRITE (6,1270) iflag,icont,icount,move,kalarm
!     CALL PDUMP (A(1,1),N,5)
      DO 1250 i=1,nds
        IF (c2(i,k).lt.0.0) GO TO 1290
 1250 END DO
      GO TO 1320
 1260 WRITE (6,1280) ierror
      WRITE (6,1270) iflag,icont,icount,move,kalarm
 1270 FORMAT ('0',10x,'PROGRAM EXIT CONDITIONS',/15x,'IFLAG=',i3,4x,'ICO&
     &NT=',i3,4x,'ICOUNT=',i3,4x,'MOVE=',i3,4x,'KALARM=',i3)
 1280 FORMAT ('0',/10x,'IERROR=',i3/)
      GO TO 1320
 1290 WRITE (6,1300)
 1300 FORMAT (2x,'NEGATIVE C2')
 1310 FORMAT ('1ERROR ORIGIN IN SUBROUTINE',1x,a6,',  STATEMENT NUMBER',&
     &i6,'.')
 1320 CONTINUE
      STOP   ! CALL exit
      END Program Aipp


      SUBROUTINE Triple()
! *** THIS ROUTINE OBTAINS THE TRIPLE-POINT SOLUTION FOR A CONFIGURATION
! ***        CONSISTING OF A WEAK INCIDENT SHOCK, A WEAK REFLECTED SHOCK
! ***        AND A STRONG SHOCK (MACH DISC).  PROPERTIES BEHIND THE
! ***        INCIDENT SHOCK AND IN THE UNSHOCKED FLOW ARE GIVEN.  THE
! ***        ROUTINE ITERATES ON THE WEAK AND STRONG SHOCKS SEPARATELY
! ***        AND USES A SECANT PROCEDURE BASED ON THE PRESSURES BEHIND
! ***        THE TWO SHOCKS AT THE SAME DEFLECTION ANGLE (TOTAL) TO
! ***        OBTAIN THE TRIPLE-POINT SOLUTION.
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /tripleBlock/ del1,pp1,tt1,uu1,th3,pp3,tt3,uu3,slang,shang
      COMMON /turple/ del2,cpshw,gamshw,zmwsh1
      COMMON /utpp/ psh11,tsh11,ush11,flang,csh11(20)
      DIMENSION cshw(20), csh5(20), chh(20)
      DIMENSION g(4), dgdx(4,4), z1(4), z2(4), z3(4,4), dum(4)
      REAL mw,mshw,msh1
      WRITE (6,10)
   10 FORMAT ('0',' TRIPLE IS CALLED')
      ido=0
      pshw=p2(2)
      tshw=t2(2)
      ushw=u2(2)
      rhoshw=rho2(2)
      zmwshw=zmwk(2)
      delphi=del2
      psh5=psh11
      tsh5=tsh11
      ush5=ush11
      DO 20 i=1,nds
   20   csh5(i)=csh11(i)
      WRITE (6,30) del2,cpshw,gamshw,zmwsh1,pshw,tshw,ushw,rhoshw,      &
     &zmwshw
   30 FORMAT (5x,9e10.3)
      DO 40 i=1,nds
        cshw(i)=c2(i,2)
   40 END DO
      mshw=ushw/sqrt(cpshw/(cpshw-rcon/zmwshw)*hj*rcon*tshw/zmwshw)
      pi=3.141592654
      twopi=2.*pi
      cpsh1=0
      DO 60 i=1,nds
        DO 50 j=1,nn
   50     cpsh1=cpsh1+float(j-2)*a(i,j)*tsh5**(j-3)*csh5(i)
   60 END DO
      msh1=ush5/sqrt(cpsh1/(cpsh1-rcon/zmwsh1)*hj*rcon*tsh5/zmwsh1)
      gamsh1=cpsh1/(cpsh1-rcon/zmwsh1)
      rhosh1=zmwsh1*psh5/(rv*tsh5)
      angm=asin(1./mshw)
      WRITE (6,70)
   70 FORMAT ('0',' PROPERTIES UPSTREAM OF STRONG SHOCK ARE')
      WRITE (6,150) psh5,tsh5,ush5,cpsh1,msh1,zmwsh1,gamsh1
      DO 80 i=1,nds
   80   WRITE (6,90) csh5(i)
   90 FORMAT (10x,' MASS FRACTIONS ARE',e10.3)
      WRITE (6,100) flang
  100 FORMAT (10x,'PHISH1=',e10.3)
      WRITE (6,110)
  110 FORMAT ('0',' PROPERTIES UPSTREAM OF WEAK SHOCK ARE')
      WRITE (6,150) pshw,tshw,ushw,cpshw,mshw,zmwshw,gamshw
      DO 120 i=1,nds
  120   WRITE (6,130) cshw(i)
  130 FORMAT (10x,' MASS FRACTIONS ARE',e10.3)
      WRITE (6,140) delphi
  140 FORMAT (10x,'DELPHI=',e10.3)
  150 FORMAT (10x,' P=',e10.3,' T=',e10.3,' U=',e10.3,' CP=',e10.3,' MA=&
     &',e10.3,' MOL. WEIGHT=',e10.3,' GAMMA=',e10.3)
! *** NORMAL SHOCK SOLUTION
! *** INITIAL GUESSES
      th1=pi/2.
      del1=0.
      tt1=3000.
      pp1=rhosh1*ush5*ush5/betap
      uu1=5.e+04
      fmass=rhosh1*ush5
      fforce=psh5+rhosh1*ush5*ush5/betap
      ht=ush5*ush5/(2.*hj)
      DO 170 i=1,nds
        DO 160 j=1,nn
  160     ht=ht+a(i,j)*tsh5**(j-2)*csh5(i)
  170 END DO
      hti=ht
      a1=zmwsh1/rv
! *** NEWTON-RAPHSON ITERATION PROCEDURE ON P,T, AND U
      iter=1
  180 a2=pp1*uu1/tt1
      g(1)=a1*a2-fmass
      g(2)=pp1*(1.+a1*uu1*uu1/(tt1*betap))-fforce
      g(3)=uu1*uu1/(2.*hj)-ht
      cpsh1=0.
      DO 200 i=1,nds
        DO 190 j=1,nn
          g(3)=g(3)+a(i,j)*tt1**(j-2)*csh5(i)
  190     cpsh1=cpsh1+float(j-2)*a(i,j)*tt1**(j-3)*csh5(i)
  200 END DO
      dgdx(1,1)=a1*uu1/tt1
      dgdx(1,2)=-a1*a2/tt1
      dgdx(1,3)=a1*pp1/tt1
      dgdx(2,1)=1.+a1*uu1*uu1/(tt1*betap)
      dgdx(2,2)=-a1*a2/(tt1*betap)
      dgdx(2,3)=2.*a1*a2/betap
      dgdx(3,1)=0.
      dgdx(3,2)=cpsh1
      dgdx(3,3)=uu1/hj
      DO 210 i=1,3
  210   g(i)=-g(i)
      yyy=0.
      DO 230 i=1,3
        z2(i)=g(i)
        DO 220 j=1,3
  220     z3(i,3)=dgdx(i,j)
  230 END DO
      z1(1)=pp1
      z1(2)=tt1
      z1(3)=uu1
      CALL mate7 (4, 3, 1, dgdx, g, yyy, dum, j)
! *** J NOT =1  MATE7 COULD NOT FIND A SOLUTION
      IF (j.eq.1) GO TO 260
      WRITE (6,240)
  240 FORMAT ('O','NORMAL SHOCK MATRIX FAILURE IN TRIPLE POINT')
  250 CALL pdump (z1(1), z2(4), 1)
      GO TO 790
! *** J=1  SOLUTION OBTAINED.  CHECK FOR CONVERGENCE
  260 DO 270 i=1,3
  270   dum(i)=dgdx(i,1)
      pp1=pp1+dum(1)
      tt1=tt1+dum(2)
      uu1=uu1+dum(3)
      IF (abs(dum(1)/pp1).gt.tol.or.abs(dum(2)/tt1).gt.tol.or.abs(dum(3)&
     &/uu1).gt.tol) GO TO 290
      DO 280 i=1,3
        IF (abs(g(i)).gt.tol) GO TO 290
  280 END DO
! *** ITERATIONS HAVE CONVERGED.  PRINT RESULTS AND BEGIN WEAK SHOCK
! ***        ITERATIONS
      GO TO 320
  290 IF (iter.le.niter) GO TO 310
! *** ITERATION LIMIT REACHED.  ADMIT FAILURE
      WRITE (6,300)
  300 FORMAT ('0','NORMAL SHOCK ITERATION FAILURE IN TRIPLE POINT')
      GO TO 250
! *** ITERATIONS HAVE NOT CONVERGED.  TRY AGAIN
  310 iter=iter+1
      GO TO 180
! *** INITIAL WEAK SHOCK SOLUTION CORRESPONDS TO NORMAL STRONG SHOCK
  320 del=abs(delphi)
      iway=0
! *** GENERAL SHOCK SOLUTION
! ***        IWAY=0 INDICATES WEAK SHOCK, SOLVE FOR THETA AND T
! ***        IWAY=1 INDICATES STRONG SHOCK, SOLVE FOR DEL AND T
  330 iter=0
! *** INITIAL GUESSES
      IF (iway.eq.1) GO TO 360
      th=del+0.1
      IF (th.lt.angm) th=angm
      emsq=(mshw*sin(th))**2
      tt=tshw*((1.+.5*(gamshw-1.)*emsq)*(2.*gamshw/(gamshw-1.)*emsq-1.)*&
     &2.*(gamshw-1.)/((gamshw+1.)**2*emsq))
      fmass=rhoshw*ushw
      fforce=rhoshw*ushw**2/betap
      pold=pshw
      uold=ushw
      told=tshw
      ht=uold*uold/(2.*hj)
      DO 350 i=1,nds
        DO 340 j=1,nn
  340     ht=ht+a(i,j)*told**(j-2)*cshw(i)
        chh(i)=cshw(i)
  350 END DO
      GO TO 380
  360 tt=tt1
      pold=psh5
      uold=ush5
      told=tsh5
      fmass=rhosh1*ush5
      pp=rhosh1*ush5*ush5/betap
      fforce=pp
      ht=hti
      th=th1
      DO 370 i=1,nds
  370   chh(i)=csh5(i)
! *** NEWTON-RAPHSON ITERATION PROCEDURE ON T AND EITHER THETA OR DEL
  380 CONTINUE
      sit=sin(th)
      cot=cos(th)
      sid=sin(th-del)
      cod=cos(th-del)
      g(1)=-ht+.5*(uold*cot/cod)**2/hj
      cp=0.
      hh=0.
      DO 400 i=1,nds
        DO 390 j=1,nn
          hh=hh+a(i,j)*tt**(j-2)*chh(i)
  390     cp=cp+float(j-2)*a(i,j)*tt**(j-3)*chh(i)
  400 END DO
      g(1)=g(1)+hh
      g(2)=pold*told*cot/sit+told*fforce*sit*cot*(1.-sid*cot/(cod*sit))-&
     &pold*tt*cod/sid
      dgdx(1,1)=cp
      dgdx(1,2)=-(cot*cot*sid/cod**3)*uold*uold/hj
      IF (iway.eq.0) dgdx(1,2)=-dgdx(1,2)-sit*cot/(cod*cod)*uold*uold/  &
     &hj
      dgdx(2,1)=-pold*cod/sid
      dgdx(2,2)=told*fforce*cot*cot/(cod*cod)-pold*tt/(sid*sid)
      IF (iway.eq.0) dgdx(2,2)=-dgdx(2,2)-pold*told/(sit*sit)+told*     &
     &fforce*((1.-sid*cot/(sit*cod))*(cot*cot-sit*sit)+sid*cot/(sit*cod)&
     &)
      DO 410 i=1,2
  410 END DO
      yyy=0.
      DO 420 i=1,2
        z2(i)=-g(i)
        g(i)=-g(i)
        DO 420 j=1,2
  420   z3(i,j)=dgdx(i,j)
      z1(1)=tt
      z1(2)=th
      IF (iway.ne.0) z1(2)=del
      CALL mate7 (4, 2, 1, dgdx, g, yyy, dum, j)
      IF (j.eq.1) GO TO 470
! *** J NOT =1  MATE7 COULD NOT FIND A SOLUTION
      IF (iway.ne.1) GO TO 440
      WRITE (6,430)
  430 FORMAT ('0','STRONG SHOCK MATRIX FAILURE IN TRIPLE POINT')
      GO TO 460
  440 WRITE (6,450)
  450 FORMAT ('0','WEAK SHOCK MATRIX FAILURE IN TRIPLE POINT')
  460 CALL pdump (z1(1), z3(4,4), 1)
      GO TO 790
! *** J=1  SOLUTION OBTAINED.  CHECK FOR CONVERGENCE
  470 DO 480 i=1,2
  480   dum(i)=dgdx(i,1)
      tt=tt+dum(1)
      IF (iway.ne.0) GO TO 490
      th=th+dum(2)
      IF (th.gt.twopi) th=th-twopi
      IF (th.lt.angm) th=angm+del
      IF (abs(dum(1)/tt).gt.tol.or.abs(dum(2)/th).gt.tol) GO TO 520
      GO TO 500
  490 del=del+dum(2)
      IF (abs(dum(1)/tt).gt.tol.or.abs(dum(2)/del).gt.tol) GO TO 520
  500 DO 510 i=1,2
        IF (abs(g(i)).gt.tol) GO TO 520
  510 END DO
      GO TO 570
  520 IF (iter.le.niter) GO TO 560
! *** ITERATION LIMIT REACHED.  ADMIT FAILURE
      IF (iway.ne.1) GO TO 540
      WRITE (6,530)
  530 FORMAT ('0','STRONG SHOCK ITERATION FAILURE IN TRIPLE POINT')
      GO TO 460
  540 WRITE (6,550)
  550 FORMAT ('0','WEAK SHOCK ITERATION FAILURE IN TRIPLE POINT')
      GO TO 460
! *** ITERATIONS HAVE NOT CONVERGED.  TRY AGAIN
  560 iter=iter+1
      GO TO 380
! *** ITERATIONS HAVE CONVERGED
  570 IF (iway.ne.0) GO TO 680
! *** WEAK SHOCK SOLUTION
      pp=pold*(1.+a1/(betap*told)*uold*uold*sit*sit*(1.-sid/cod*cot/sit)&
     &)
      uu=sqrt(2.*hj*(ht-hh))
! *** SECANT PROCEDURE FOR BALANCING WEAK AND STRONG SHOCK RESULTS
      th3=th
      ido=ido+1
      pp3=pp
      dp=pp1-pp3
      tt3=tt
      uu3=uu
      IF (ido.ne.1) GO TO 590
! *** FIRST TIME THROUGH.  CHANGE TH1, GO BACK AND GET A SECOND POINT
      dp1=dp
      tha=th1
      igo=0
      th1=tha+sign(1.,dp)*.017453295
  580 iway=1
      tandel=((gamsh1+1.)*msh1*msh1/(2.*((msh1*sin(th1))**2-1.))-1.)*   &
     &tan(th1)
      del=atan(1./tandel)
      GO TO 650
  590 IF (ido.ne.2) GO TO 610
! *** SECOND TIME THROUGH.  CHECK TO SEE IF TRIPLE-POINT SOLUTION HAS
! ***        BEEN BRACKETED
      dp2=dp
      thb=th1
      IF (sign(1.,dp2)*sign(1.,dp1).lt.0.) GO TO 600
! *** NOT BRACKETED.  CHANGE TH1 IN SAME DIRECTION AND TRY AGAIN
      th1=thb+sign(1.,dp)*.017453295
      GO TO 580
! *** BRACKETED.  SET IGO=1 AND USE SECANT PROCEDURE TO CHANGE TH1
  600 igo=1
      th1=thb-dp2/(dp2-dp1)*(thb-tha)
      GO TO 580
! *** THIRD OR GREATER TIME THROUGH.  IGO=1 INDICATES SOLUTION HAS
! ***        PREVIOUSLY BEEN BRACKETED
  610 IF (igo.eq.1) GO TO 630
! *** SOLUTION NOT PREVIOUSLY BRACKETED.  CHECK TO SEE ON WHICH SIDE OF
! ***        SOLUTION THE CURRENT VALUE OF DP FALLS, UPDATE DP1 OR DP2
! ***        AS REQUIRED, AND CHANGE TH1.  TRY AGAIN
      IF (sign(1.,dp)*sign(1.,dp2).lt.0) GO TO 620
      dp1=dp2
      tha=thb
      dp2=dp
      thb=th1
      th1=thb+sign(1.,dp)*.017453295
      GO TO 580
  620 dp1=dp2
      tha=thb
      dp2=dp
      thb=th1
      th1=thb-dp2/(dp2-dp1)*(thb-tha)
      igo=1
      GO TO 580
  630 IF (sign(1.,dp)*sign(1.,dp2).lt.0) GO TO 640
      thb=th1
      dp2=dp
      th1=thb-dp2/(dp2-dp1)*(thb-tha)
      GO TO 580
  640 tha=thb
      dp1=dp2
      thb=th1
      dp2=dp
      th1=thb-dp2/(dp2-dp1)*(thb-tha)
      GO TO 580
! *** HAS SECANT ITERATION LIMIT BEEN EXCEEDED
  650 IF (ido.lt.niter) GO TO 690
! *** ITERATION LIMIT REACHED.  ADMIT FAILURE
      WRITE (6,660)
  660 FORMAT ('0',' TRIPLE POINT ITERATIONS FAILED')
      WRITE (6,670) pp1,pp3,del1,del3,th1,th3
  670 FORMAT ('  PP1=',e10.3,' PP3=',e10.3,' DEL1=',e10.3,'DEL2=',e10.3,&
     &' TH=',e10.3,' TH3=',e10.3)
      GO TO 460
! *** STRONG SHOCK SOLUTION OBTAINED.  GO BACK AND GET A CORRESPONDING
! ***        WEAK SHOCK SOLUTION
  680 iway=0
      pp=pold*(1.+a1/(betap*told)*uold*uold*sit*sit*(1.-sid/cod*cot/sit)&
     &)
      uu=sqrt(2.*hj*(ht-hh))
      del1=del
      th1=th
      pp1=pp
      tt1=tt
      uu1=uu
      del3=abs(delphi)-del1
      del=del3
      GO TO 330
! *** CHECK FOR CONVERGENCE OF SECANT ITERATIONS
  690 IF (abs((pp1-pp3)/(pp1+pp3)).lt.0.001) GO TO 700
      GO TO 330
! *** ITERATIONS HAVE CONVERGED.   PRINT OUT RESULTS AND RETURN
  700 CONTINUE
      fm3=uu3/sqrt(cp/(cp-rcon/zmwshw)*tt3*rcon*hj/zmwshw)
      WRITE (6,710)
  710 FORMAT ('0',' PROPERTIES DOWNSTREAM OF WEAK SHOCK ARE')
      WRITE (6,720) pp3,tt3,uu3,th3,del3,fm3
  720 FORMAT (10x,' P=',e10.3,' T=',e10.3,' U=',e10.3,' TH3=',e10.3,' DE&
     &L3=',e10.3,' MA=',e10.3)
      cp=0.
      DO 740 i=1,nds
        DO 730 j=1,nn
  730     cp=cp+a(i,j)*float(j-2)*tt1**(j-3)*csh5(i)
  740 END DO
      fm1=uu1/sqrt(cp/(cp-rcon/zmwsh1)*tt1*rcon*hj/zmwsh1)
      WRITE (6,750)
  750 FORMAT ('0',' PROPERTIES DOWNSTREAM OF STRONG SHOCK ARE')
      WRITE (6,760) pp1,tt1,uu1,th1,del1,fm1
  760 FORMAT (10x,' P=',e10.3,' T=',e10.3,' U=',e10.3,' TH1=',e10.3,' DE&
     &L1=',e10.3,' MA=',e10.3)
      thtp=flang
      bao=abs(th1)-flang
      slang=thtp+abs(del1)
      shang=abs(th3)-abs(del2-thtp)
      WRITE (6,770) shang,bao,slang
  770 FORMAT ('0,'' WEAK SHOCK ANGLE=',e10.3,' STRONG SHOCK ANGLE=',    &
     &e10.3,' SLIP LINE ANGLE=',e10.3)
      WRITE (6,780) ido
  780 FORMAT (5x,' NUMBER OF ITERATION=',i5)
  790 RETURN
      END Subroutine Triple

      SUBROUTINE bndry
!     THIS SUBROUTINE HANDLES MOST BOUNDARY CONDITIONS
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      EQUIVALENCE (move,iextra(15))
      EQUIVALENCE (pwsh,extra(11))
      IF (move) 640,10,640
   10 IF (nbound) 310,20,310
!     WALL CONDITIONS
   20 k=1
      IF (itype-2) 140,30,40
   30 pwsh=pw(1)+(pw(2)-pw(1))*s/sw(2)
      GO TO 40
!     SHOCK BOUNDARY
   40 IF (iflag) 60,50,60
   50 dpdy(k)=2.0*(p(k+1)-pwsh)/(y(k+1)-y(k))
      dphids(k)=-dpdy(k)*betap/(rho(k+1)*u(k+1)**2)
      GO TO 70
   60 dpdy(k)=2.0*(p2(k+1)-pwsh)/(y2(k+1)-y(k))
      dphods=-dpdy(k)*betap/(rho2(k+1)*u2(k+1)**2)
      dphids(k)=.45*dphids(k)+.55*dphods
   70 phi2(k)=phi(k)+dphids(k)*dls(k)
      IF (abs(phi2(k)-phi(k))-10.0**(-4)) 80,80,90
   80 x2(k)=x(k)+dls(k)*cos(phi(k))
      r2(k)=r(k)+dls(k)*sin(phi(k))
      GO TO 100
   90 radw=dls(k)/(phi(k)-phi2(k))
      xo=x(k)+radw*sin(phi(k))
      ro=r(k)-radw*cos(phi(k))
      x2(k)=xo-radw*sin(phi2(k))
      r2(k)=ro+radw*cos(phi2(k))
  100 y2(k)=0.0
      y(k)=0.0
      max=kmax
      IF (itype-3) 200,110,110
  110 CALL shocke
      IF (kmax-max) 120,120,130
  120 k=1
      GO TO 200
  130 k=2
      GO TO 200
  140 m=1
      max=kmax
!     FIXED WALL
  150 IF (s-sw(m+1)) 170,170,160
  160 m=m+1
      IF (m-mmax) 150,150,790
  170 k=1
      y(k)=0.0
      y2(k)=0.0
      phi2(k)=phiw(m)+(phiw(m+1)-phiw(m))*(s-sw(m))/(sw(m+1)-sw(m))
      IF (abs(phiw(m+1)-phiw(m))-10.**(-4)) 180,180,190
  180 x2(k)=xw(m)+(xw(m+1)-xw(m))*(s-sw(m))/(sw(m+1)-sw(m))
      r2(k)=rw(m)+(rw(m+1)-rw(m))*(s-sw(m))/(sw(m+1)-sw(m))
      GO TO 200
  190 radw=delso/(phi(k)-phi2(k))
      xo=x(k)+radw*sin(phi(k))
      ro=r(k)-radw*cos(phi(k))
      x2(k)=xo-radw*sin(phi2(k))
      r2(k)=radw*cos(phi2(k))+ro
  200 IF (iflag) 210,240,210
  210 rbar(k)=.5*(r(k)+r2(k))
      rbar(k+1)=.5*(r(k+1)+r2(k+1))
      yabar=.5*(y(k+1)+y2(k+1))
      ybbar=0.0
      ubar(k+1)=sqrt(.5*(u(k+1)**2+u2(k+1)**2))
      rhabar=.5*(rho(k+1)+rho2(k+1))
      tabar=.5*(t(k+1)+t2(k+1))
      pabar=.5*(p(k+1)+p2(k+1))
      mw2=0.0
      DO 220 i=1,nds
        cabar(i)=.5*(c(i,k+1)+c2(i,k+1))
  220   mw2=mw2+cabar(i)/mw(i)
      mw2=1.0/mw2
      DO 230 i=1,nds
        xabar(i)=cabar(i)*mw2/mw(i)
        h2pm(i)=0.0
        hpm(i)=0.0
        DO 230 j=1,nn
  230   hpm(i)=hpm(i)+a(i,j)*tabar**(j-2)
      GO TO 260
  240 rbar(k)=r(k)
      rbar(k+1)=r(k+1)
      yabar=y(k+1)
      ybbar=0.0
      ubar(k+1)=u(k+1)
      rhabar=rho(k+1)
      tabar=t(k+1)
      pabar=p(k+1)
      DO 250 i=1,nds
        cabar(i)=c(i,k+1)
        xabar(i)=xs(i,k+1)
        h2pm(i)=0.0
        hpm(i)=0.0
        DO 250 j=1,nn
  250   hpm(i)=hpm(i)+a(i,j)*t(k+1)**(j-2)
  260 IF (iextra(1).eq.0) GO TO 270
      kk=k+1
      CALL transp (tabar, pabar, kk)
  270 taw(k)=0.0
      q(k)=0.0
      DO 280 i=1,nds
  280   zj(i,k)=0.0
      IF (kmax-max) 290,290,300
  290 k=2
      GO TO 820
  300 k=3
      GO TO 820
!     OUTER BOUNDARY CONDITIONS
  310 n=1
      IF (ll) 320,320,350
  320 d=sqrt((x(k)-xb(n))**2+(r(k)-rb(n))**2)
      IF (abs(phi(k)-phib(n))-1.e-04) 330,330,340
  330 sx(k)=d
      GO TO 350
  340 sx(k)=d*(phi(k)-phib(n))/(2.*sin(0.5*(phi(k)-phib(n))))
  350 GO TO (360,650,660,750),ikind
!     FIXED BOUNDARY
  360 n=1
      ssx=sx(k)+dls(k-1)
  370 IF (ssx-sb(n+1)) 390,390,380
  380 n=n+1
      IF (n-nmax) 370,800,800
  390 IF (abs(phib(n+1)-phib(n))-1.e-04) 400,400,440
  400 phibar=0.5*(phi2(k-1)+phib(n))
      tphbr=tan(phibar)
      IF (abs(phibar)-1.e-04) 410,410,420
  410 x2(k)=x2(k-1)
      GO TO 430
  420 x2(k)=(r2(k-1)-r(k)+tan(phib(n))*x(k)+x2(k-1)/tphbr)/(tan(phib(n))&
     &+1./tphbr)
  430 r2(k)=r(k)+tan(phib(n))*(x2(k)-x(k))
      phi2(k)=phib(n)
      GO TO 570
  440 radb=(sb(n+1)-sb(n))/(phib(n)-phib(n+1))
      ell1=sqrt((xb(n+1)-xb(n))**2+(rb(n+1)-rb(n))**2)
      ell2=sqrt(radb**2-ell1**2/4.)
      fiota=atan((rb(n+1)-rb(n))/(xb(n+1)-xb(n)))
      IF (radb) 460,460,450
  450 xx=xb(n)+.5*(xb(n+1)-xb(n))+ell2*sin(fiota)
      rr=rb(n)+.5*(rb(n+1)-rb(n))-ell2*cos(fiota)
      GO TO 470
  460 xx=xb(n)+.5*(xb(n+1)-xb(n))-ell2*sin(fiota)
      rr=rb(n)+.5*(rb(n+1)-rb(n))+ell2*cos(fiota)
  470 l=1
      eks=x(k)+dls(k-1)*cos(phi(k))
      IF (radb) 480,490,490
  480 ahr=rr-sqrt(radb**2-(eks-xx)**2)
      GO TO 500
  490 ahr=rr+sqrt(radb**2-(eks-xx)**2)
  500 phibar=.5*(phi2(k-1)+atan((eks-xx)/(rr-ahr)))
      tphbr=tan(phibar)
      gee=(r2(k-1)-rr-(eks-x2(k-1))/tphbr)**2-radb**2+(eks-xx)**2
      drdx=(eks-xx)/sqrt(radb**2-(eks-xx)**2)
      IF (radb) 520,520,510
  510 drdx=-drdx
  520 dphdx=-.5*((ahr-rr)-(eks-xx)*drdx)/radb**2
      dtph=dphdx/(cos(phibar))**2
      dgeedx=2.*((eks-xx)-(r2(k-1)-rr-(eks-x2(k-1))/tphbr)*(1./tphbr-   &
     &(eks-x2(k-1))*dtph/(tan(phibar))**2))
      dx=-gee/dgeedx
      l=l+1
      eks=eks+dx
      IF (radb) 530,540,540
  530 ahr=rr-sqrt(radb**2-(eks-xx)**2)
      GO TO 550
  540 ahr=rr+sqrt(radb**2-(eks-xx)**2)
  550 IF (l-4) 500,560,560
  560 x2(k)=eks
      r2(k)=ahr
      phi2(k)=-atan((eks-xx)/(ahr-rr))
  570 d=sqrt((x2(k)-x(k))**2+(r2(k)-r(k))**2)
      IF (abs(phi2(k)-phi(k))-1.e-04) 580,580,590
  580 dls(k)=d
      GO TO 600
  590 dls(k)=d*(phi2(k)-phi(k))/(2.*sin(0.5*(phi2(k)-phi(k))))
  600 b=sqrt((x2(k)-x2(k-1))**2+(r2(k)-r2(k-1))**2)
      IF (abs(phi2(k)-phi2(k-1))-1.e-04) 610,610,620
  610 dely(k)=b
      GO TO 630
  620 dely(k)=b*(phi2(k)-phi2(k-1))/(2.0*sin(.5*(phi2(k)-phi2(k-1))))
  630 aa(k)=pi*(r2(k)+r2(k-1))**delta*dely(k)
      y2(k)=y2(k-1)+dely(k)
      CALL step
  640 move=0
      GO TO 820
!     FREE BOUNDARY
  650 n=1
      pbs=pb(n)+(pb(n+1)-pb(n))*(s2x-sb(n))/(sb(n+1)-sb(n))
      GO TO 750
!     NEWTONIAN BOUNDARY
  660 fmstrm=0.0
      DO 670 i=1,nds
  670   fmstrm=fmstrm+cstrem(i)/mw(i)
      fmstrm=1.0/fmstrm
      qstrm=pstrem*fmstrm*(ustrem**2)/(rv*tstrem*betap)
      IF (iflag) 690,680,690
  680 sigma=phi(k)-phstrm
      GO TO 700
  690 sigma=.5*(phi(k)+phi2(k))-phstrm
  700 IF (sigma) 710,710,720
  710 pbs=pstrem
      GO TO 750
  720 IF (sigma-0.5*pi) 740,740,730
  730 pbs=pstrem+qstrm
      GO TO 750
  740 pbs=pstrem+qstrm*(sin(sigma))**2
  750 IF (iflag) 770,760,770
  760 dpdy(k)=2.0*(pbs-p(k))/(y(k)-y(k-1))
      dphids(k)=-dpdy(k)*betap/(rho(k)*u(k)**2)
      GO TO 780
  770 dpdy(k)=2.0*(pbs-p2(k))/(y2(k)-y2(k-1))
      dphods=-dpdy(k)*betap/(rho2(k)*u2(k)**2)
      dphids(k)=.45*dphids(k)+.55*dphods
  780 move=1
      RETURN
!     ERROR MESSAGES
      nucase=100
      GO TO 810
  790 nucase=400
      GO TO 810
  800 nucase=3200
  810 ierror=2
  820 RETURN
      END Subroutine Bndry

      SUBROUTINE chem (ikine)
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /chem1/ zid(5),irr(30),irt(30),rc(30,3),irrr(30,5),av,     &
     &cm(21,21),fi(21),wp(21),wm(21),wdot(21,30),csi(21),qx(21)
      COMMON /ibugs1/ ibugsh
      DIMENSION g(20), rp(20), rm(20)
!   GAS CONSTANT IN 1.987 CAL/G-MOLE-K
      rrt=1.987*t(k)
!   GAS CONSTANT IN 82.06 CM3-ATM/G-MOLE-K
      rrt1=82.06*t(k)
      xllgt=alog(t(k)/1000.0)
      DO 20 ir=1,nds
        fi(ir)=c(ir,k)/mw(ir)
        qx(ir)=0.0
        wp(ir)=0.0
        wm(ir)=0.0
        wdot(ir,k)=0.0
        DO 10 jr=1,nds
   10     cm(ir,jr)=0.0
   20 END DO
!     CALCULATING GIBBS ENERGY
      DO 40 ir=1,nds
        g(ir)=a(ir,1)*.5/t(k)+(a(ir,3)*(1.-xllgt)-csi(ir))*t(k)+a(ir,2)
        IF (nn.lt.4) GO TO 30
        g(ir)=g(ir)-a(ir,4)*t(k)**2-.5*a(ir,5)*t(k)**3-a(ir,6)/3.*t(k)**&
     &   4
   30   g(ir)=g(ir)*mw(ir)
!  CALCULATE GIBBS ENERGY IN CAL/MOLE
   40 END DO
      DO 370 ir=1,ikine
        kirt=irt(ir)
        GO TO (60,70,80,90,100,110,120,50),kirt
   50   rkf=rc(ir,1)*(t(k)**rc(ir,2))*exp(rc(ir,3)/rrt)
        GO TO 130
   60   rkf=rc(ir,1)
        GO TO 130
   70   rkf=rc(ir,1)/t(k)
        GO TO 130
   80   rkf=rc(ir,1)/t(k)/t(k)
        GO TO 130
   90   rkf=rc(ir,1)/sqrt(t(k))
        GO TO 130
  100   rkf=rc(ir,1)*exp(rc(ir,3)/rrt)
        GO TO 130
  110   rkf=rc(ir,1)*exp(rc(ir,3)/rrt)/t(k)
        GO TO 130
  120   rkf=rc(ir,1)/t(k)/sqrt(t(k))
  130   CONTINUE
        kirr=irr(ir)
        GO TO (140,160,180,200,220,240,260,280,300,320),kirr
  140   j1=irrr(ir,1)
        j2=irrr(ir,2)
        j3=irrr(ir,3)
        j4=irrr(ir,4)
        e=(g(j1)+g(j2)-g(j3)-g(j4))/rrt
        IF (e.lt.-40.0) e=-40.0
        IF (e.gt.+40.0) e=+40.0
        e=exp(e)
        crr=rkf*rho(k)
        rp(ir)=crr*fi(j1)*fi(j2)
        rm(ir)=crr*fi(j3)*fi(j4)/e
        DO 150 j=1,4
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr*fi(j2)
          cm(irow,j2)=cm(irow,j2)+sign*crr*fi(j1)
          cm(irow,j3)=cm(irow,j3)-sign*crr*fi(j4)/e
          cm(irow,j4)=cm(irow,j4)-sign*crr*fi(j3)/e
          qx(irow)=qx(irow)+sign*(rp(ir)-rm(ir))
  150   CONTINUE
        GO TO 350
  160   j1=irrr(ir,1)
        j2=irrr(ir,2)
        j3=irrr(ir,3)
        e=(g(j1)+g(j2)-g(j3))/rrt
        IF (e.lt.-40.0) e=-40.0
        IF (e.gt.+40.0) e=+40.0
        e=exp(e)
        crr=rkf*rho(k)/zmw*av
        rp(ir)=crr*rho(k)*fi(j1)*fi(j2)
        rm(ir)=crr*fi(j3)/(e*rrt1)
        DO 170 j=1,3
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr*rho(k)*fi(j2)
          cm(irow,j2)=cm(irow,j2)+sign*crr*rho(k)*fi(j1)
          cm(irow,j3)=cm(irow,j3)-sign*crr/(e*rrt1)
          qx(irow)=qx(irow)+sign*rp(ir)
  170   CONTINUE
        GO TO 360
  180   j1=irrr(ir,1)
        j2=irrr(ir,2)
        j3=irrr(ir,3)
        j4=irrr(ir,4)
        j5=irrr(ir,5)
        e=(g(j1)+g(j2)-g(j3)-g(j4)-g(j5))/rrt
        IF (e.lt.-40.0) e=-40.0
        IF (e.gt.+40.0) e=+40.0
        e=exp(e)
        crr=rkf*rho(k)
        rp(ir)=crr*fi(j1)*fi(j2)
        rm(ir)=crr*fi(j3)*fi(j4)*fi(j5)*rho(k)*rrt1/e
        DO 190 j=1,5
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr*fi(j2)
          cm(irow,j2)=cm(irow,j2)+sign*crr*fi(j1)
          cm(irow,j3)=cm(irow,j3)-sign*crr*fi(j4)*fi(j5)*rho(k)*rrt1/e
          cm(irow,j4)=cm(irow,j4)-sign*crr*fi(j3)*fi(j5)*rho(k)*rrt1/e
          cm(irow,j5)=cm(irow,j5)-sign*crr*fi(j3)*fi(j4)*rho(k)*rrt1/e
          qx(irow)=qx(irow)+sign*(rp(ir)-2.0*rm(ir))
  190   CONTINUE
        GO TO 340
  200   j1=irrr(ir,1)
        j2=irrr(ir,2)
        j3=irrr(ir,3)
        e=(g(j1)+g(j2)-g(j3))/rrt
        IF (e.lt.-40.0) e=-40.0
        IF (e.gt.+40.0) e=+40.0
        e=exp(e)
        crr=rkf*rho(k)
        rp(ir)=crr*fi(j1)*fi(j2)
        rm(ir)=crr*fi(j3)/e/rrt1*rho(k)
        DO 210 j=1,3
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr*rho(k)*fi(j2)
          cm(irow,j2)=cm(irow,j2)+sign*crr*rho(k)*fi(j1)
          cm(irow,j3)=cm(irow,j3)-sign*crr/(e*rrt1)*rho(k)
          qx(irow)=qx(irow)+sign*rp(ir)
  210   CONTINUE
        GO TO 360
  220   j1=irrr(ir,1)
        j2=nds+1
        j3=irrr(ir,3)
        j4=irrr(ir,4)
        e=(g(j1)-g(j3)-g(j4))/rrt
        IF (e.lt.-40.0) e=-40.0
        IF (e.gt.+40.0) e=+40.0
        e=exp(e)
        crr=rkf*rho(k)/zmw
        rp(ir)=crr*fi(j1)
        rm(ir)=crr*rrt1*rho(k)*fi(j3)*fi(j4)/e
        DO 230 j=1,4
          IF (j.eq.2) GO TO 230
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr
          cm(irow,j3)=cm(irow,j3)-sign*crr*rrt1*rho(k)*fi(j4)/e
          cm(irow,j4)=cm(irow,j4)-sign*crr*rrt1*rho(k)*fi(j3)/e
          qx(irow)=qx(irow)-sign*rm(ir)
  230   CONTINUE
        GO TO 350
  240   j1=irrr(ir,1)
        j2=irrr(ir,2)
        j3=irrr(ir,3)
        j4=irrr(ir,4)
        crr=rkf*rho(k)
        rp(ir)=crr*fi(j1)*fi(j2)
        rm(ir)=0.0
        DO 250 j=1,4
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr*fi(j2)
          cm(irow,j2)=cm(irow,j2)+sign*crr*fi(j1)
          qx(irow)=qx(irow)+sign*rp(ir)
  250   CONTINUE
        GO TO 350
  260   j1=irrr(ir,1)
        j2=irrr(ir,2)
        j3=irrr(ir,3)
        crr=rkf*rho(k)*av/zmw
        rp(ir)=crr*rho(k)*fi(j1)*fi(j2)
        rm(ir)=0.0
        DO 270 j=1,3
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr*rho(k)*fi(j2)
          cm(irow,j2)=cm(irow,j2)+sign*crr*rho(k)*fi(j1)
          qx(irow)=qx(irow)+sign*rp(ir)
  270   CONTINUE
        GO TO 360
  280   j1=irrr(ir,1)
        j2=irrr(ir,2)
        j3=irrr(ir,3)
        j4=irrr(ir,4)
        j5=irrr(ir,5)
        crr=rkf*rho(k)
        rp(ir)=crr*fi(j1)*fi(j2)
        rm(ir)=0.0
        DO 290 j=1,5
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr*fi(j2)
          cm(irow,j2)=cm(irow,j2)+sign*crr*fi(j1)
          qx(irow)=qx(irow)+sign*rp(ir)
  290   CONTINUE
        GO TO 340
  300   j1=irrr(ir,1)
        j2=irrr(ir,2)
        j3=irrr(ir,3)
        crr=rkf
        rp(ir)=crr*rho(k)*fi(j1)*fi(j2)
        rm(ir)=0.0
        DO 310 j=1,3
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr*rho(k)*fi(j2)
          cm(irow,j2)=cm(irow,j2)+sign*crr*rho(k)*fi(j1)
          qx(irow)=qx(irow)+sign*rp(ir)
  310   CONTINUE
        GO TO 360
  320   j1=irrr(ir,1)
        j2=nds+1
        j3=irrr(ir,3)
        j4=irrr(ir,4)
        crr=rkf*rho(k)/zmw
        rp(ir)=crr*fi(j1)
        rm(ir)=0.0
        DO 330 j=1,4
          IF (j.eq.2) GO TO 330
          sign=1.0
          IF (j.gt.2) sign=-1.0
          irow=irrr(ir,j)
          cm(irow,j1)=cm(irow,j1)+sign*crr
  330   CONTINUE
        GO TO 350
  340   wp(j5)=wp(j5)+rp(ir)
        wm(j5)=wm(j5)+rm(ir)
  350   wp(j4)=wp(j4)+rp(ir)
        wm(j4)=wm(j4)+rm(ir)
  360   wp(j3)=wp(j3)+rp(ir)
        wm(j3)=wm(j3)+rm(ir)
        wp(j2)=wp(j2)-rp(ir)
        wm(j2)=wm(j2)-rm(ir)
        wp(j1)=wp(j1)-rp(ir)
        wm(j1)=wm(j1)-rm(ir)
  370 END DO
      DO 380 j=1,nds
!      PRODUCTION RATE IN MOLE/GM-CM
        wdot(j,k)=(wp(j)-wm(j))/u(k)
        wdot(j,k)=wdot(j,k)*mw(j)
        IF (t(k).lt.350.) wdot(j,k)=0.
  380 END DO
      dsuk=dls(k)/u(k)
      DO 410 ir=1,nds
        dswuk=dsuk
        qx(ir)=c(ir,k)/mw(ir)+rhs(ir)/mw(ir)/mdot(k)+dswuk*qx(ir)
        DO 390 jr=1,nds
          cm(ir,jr)=cm(ir,jr)*dswuk
          IF (ir.eq.jr) cm(ir,jr)=1.0+cm(ir,jr)
  390   CONTINUE
        IF (ibugsh.ne.0) WRITE (6,400) k,ir,(cm(ir,jr),jr=1,nds),qx(ir)
  400   FORMAT (1x,2i5,1p8e12.3)
  410 END DO
      CALL sldp (qx, cm, nds)
      DO 420 ir=1,nds
        IF (ibugsh.ne.0) WRITE (6,400) k,ir,(cm(ir,jr),jr=1,nds),qx(ir)
        c2(ir,k)=qx(ir)*mw(ir)
!     CALCULATE NEW MASS FRACTION OF SPECIES
  420 END DO
      RETURN
      END Subroutine Chem

      SUBROUTINE combo (l)
!     THIS SUBROUTINE COMBINES TWO STREAMTUBES, CONSERVING MASS,
!     MOMENTIM, AND ENERGY, WHENEVER THE ALLOWED STEPPING DISTANCE
!     BECOMES TOO LOW.
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /ibugs1/ ibugsh
      DIMENSION ftest(3), ntest(3)
      IF ((ll.lt.10).and.(kmn.gt.1)) GO TO 690
      IF (par(10).eq.0.0) GO TO 10
      jazz=par(10)
      IF ((ll.lt.jazz).and.(kmn.gt.1)) GO TO 690
   10 CONTINUE
      sumdot=0.0
      DO 20 kk=2,kmax
   20   sumdot=sumdot+mdot(kk)
      IF ((isb.eq.1).and.(kmax.gt.28)) GO TO 30
      IF (mdot(l)/(.5*(r(l)+r(l-1)))**2-sumdot/(1.5*(.5*(r(2)+r(kmax)))*&
     &*2*fstep)) 30,690,690
!     SETUP OF TESTS TO DETERMINE IF GRADIENTS ARE TOO HIGH TO ALLOW
!     COMBINING TUBES
   30 ki=l-1
      IF ((l.le.3).and.(isb.eq.1)) GO TO 690
      kf=l+1
      DO 40 l2=1,3
   40   ntest(l2)=0
      ntry=0
      GO TO (50,70,90,110),iptuc
   50 DO 60 kj=ki,kf
        iarg=kj-ki+1
   60   ftest(iarg)=p(kj)
      GO TO 180
   70 DO 80 kj=ki,kf
        iarg=kj-ki+1
   80   ftest(iarg)=u(kj)
      GO TO 180
   90 DO 100 kj=ki,kf
        iarg=kj-ki+1
  100   ftest(iarg)=u(kj)
      GO TO 180
  110 DO 120 kj=ki,kf
        iarg=kj-ki+1
  120   ftest(iarg)=c(1,kj)
      IF (abs(ftest(1)-ftest(2))-1.e-05) 150,150,130
  130 IF (iextra(1)) 150,140,150
  140 ntest(1)=1
  150 IF (abs(ftest(3)-ftest(2))-1.e-05) 180,180,160
  160 IF (iextra(1)) 180,170,180
  170 ntest(3)=1
  180 IF (l-kmax) 190,210,210
  190 IF (l-2) 230,230,200
  200 IF (mdot(l+1)-mdot(l-1)) 230,230,210
  210 IF (l-2) 690,690,220
  220 k2=ki
      l2=1
      rup=r(l)
      rlo=r(ki-1)
      xup=x(l)
      xlo=x(ki-1)
      phiup=phi(l)
      philo=phi(ki-1)
      GO TO 250
  230 IF (l-kmax) 240,690,690
  240 k2=kf
      l2=3
      rup=r(kf)
      rlo=r(ki)
      xup=x(kf)
      xlo=x(ki)
      phiup=phi(kf)
      philo=phi(ki)
  250 IF (ftest(2)) 290,260,290
  260 IF (ftest(l2)) 270,280,270
  270 IF (abs(ftest(2)-ftest(l2))/ftest(l2)-grad) 300,310,310
  280 IF (grad) 310,300,310
  290 IF (abs(ftest(2)-ftest(l2))/ftest(2)-grad) 300,310,310
  300 IF (ntest(l2)) 310,330,310
  310 IF (ntry) 690,320,690
  320 ntry=1
      IF (l2-2) 230,690,210
!     COMBINING TUBES
  330 d=sqrt((xup-xlo)**2+(rup-rlo)**2)
      IF (abs(phiup-philo)-(10.**(-4))) 350,350,340
  340 d=d*(phiup-philo)/(2.0*sin(0.5*(phiup-philo)))
  350 dely(l)=d
      fmom=mdot(l)*u(l)+mdot(k2)*u(k2)+(aa(l)*p(l)+aa(k2)*p(k2))*betap
      aa(l)=pi*(rup+rlo)**delta*dely(l)
      fmdot=mdot(l)
      mdot(l)=fmdot+mdot(k2)
      fh=(fmdot*h(l)+mdot(k2)*h(k2)+(fmdot*u(l)**2+mdot(k2)*u(k2)**2)/  &
     &(2.0*hj))/mdot(l)
      szmwc=0.0
      DO 360 i=1,nds
        c(i,l)=(fmdot*c(i,l)+mdot(k2)*c(i,k2))/mdot(l)
        zmwc=c(i,l)/mw(i)
  360   szmwc=szmwc+zmwc
      zmw=1.0/szmwc
      iter=0
      cp=0.0
      cp1=0.0
      cp2=0.0
      ts=0.5*(t(l)+t(k2))
      us=0.5*(u(l)+u(k2))
      DO 370 i=1,nds
        xs(i,l)=c(i,l)*zmw/mw(i)
        DO 370 j=1,nn
        cp1=cp1+a(i,j)*float(j-2)*c(i,l)*t(k2)**(j-3)
        cp2=cp2+a(i,j)*float(j-2)*c(i,l)*t(k)**(j-3)
  370   cp=cp+a(i,j)*float(j-2)*c(i,l)*t(l)**(j-3)
      b1=fmom/(betap*aa(l))
      b2=mdot(l)/(betap*aa(l))
      b3=zmw*aa(l)/(rv*mdot(l))
      IF (iconst) 380,390,380
  380 IF (0.5*(ma(l)+ma(k2))-20.) 440,390,390
!     CONSTANT HEAT CAPACITY
  390 b4=2.0*hj*(cp1*t(k2)*mdot(k2)+cp2*t(l)*fmdot)/mdot(l)+(mdot(k2)*  &
     &u(k2)**2+fmdot*u(l)**2)/mdot(l)
      b5=hj*cp*b3*b1
      b6=2.0*hj*cp*b3*b2-1.0
      discr=b5**2-b6*b4
      IF (discr) 400,430,430
  400 WRITE (6,410)
  410 FORMAT ('1NEGATIVE DISCRIMINANT IN COMBO')
  420 CALL pdump (a(1,1), n, 5)
      STOP   ! CALL exit
  430 u(l)=(b5+sqrt(discr))/b6
      GO TO 500
!     VARIABLE HEAT CAPACITY
  440 iter=iter+1
      b4=2.0*hj*fh
      dum1=0.0
      dum2=0.0
      DO 450 i=1,nds
        DO 450 j=1,nn
        dum1=dum1+a(i,j)*c(i,l)*(b3*(b1-b2*us)*us)**(j-2)
  450   dum2=dum2+a(i,j)*c(i,l)*float(j-2)*(b1*us-b2*us**2)**(j-3)*b3** &
     &   (j-2)
      g=2.0*hj*dum1+us**2-b4
      dgdu=2.0*us+2.0*hj*dum2*(b1-2.0*b2*us)
      deltau=-g/dgdu
      us=us+deltau
      IF (abs(deltau/us)-tol) 490,460,460
  460 IF (iter-niter) 440,440,470
  470 WRITE (6,480)
  480 FORMAT ('1ITERATIONS FAILED IN COMBO')
      GO TO 420
  490 u(l)=us
  500 p(l)=b1-b2*u(l)
      t(l)=b3*p(l)*u(l)
      h(l)=0.0
      DO 510 i=1,nds
        DO 510 j=1,nn
  510   h(l)=h(l)+a(i,j)*c(i,l)*t(l)**(j-2)
      rho(l)=zmw*p(l)/(rv*t(l))
      fmax=1.0
      IF (iextra(1)) 530,520,530
  520 re(l)=1.e+20
      GO TO 580
  530 CALL transp (t(l), p(l), l)
      IF (kay-mu*cp) 550,550,540
  540 fmax=kay/(cp*mu)
  550 DO 570 i=1,nds
        max=i
        DO 570 j=1,max
        IF (dih(i,j)*rho(l)/mu-fmax) 570,570,560
  560   fmax=dih(i,j)*rho(l)/mu
  570 CONTINUE
      dely(l)=dely(l)+dely(k2)
      re(l)=rho(l)*u(l)*dely(l)/(mu*fmax)
  580 IF (kmn.gt.1) GO TO 610
      sonic=(cp*zmw-rcon)/(cp*t(l)*ru)
      IF (sonic.gt.0.0) GO TO 600
      WRITE (6,590)
  590 FORMAT ('NEG SONIC')
      GO TO 420
  600 ma(l)=u(l)*sqrt(sonic)
!     REINDEXING STREAMTUBE PROPERTIES
  610 IF (k2-k) 620,620,630
  620 l=l-1
      GO TO 640
  630 phi(l)=phi(l+1)
      r(l)=r(l+1)
      sx(l)=sx(l+1)
      x(l)=x(l+1)
      y(l)=y(l+1)
      l=l+1
  640 kmax=kmax-1
      IF (l-kmax) 650,650,680
  650 DO 670 l2=l,kmax
        aa(l2)=aa(l2+1)
        DO 660 i=1,nds
          c(i,l2)=c(i,l2+1)
  660     xs(i,l2)=xs(i,l2+1)
        dely(l2)=dely(l2+1)
        h(l2)=h(l2+1)
        mdot(l2)=mdot(l2+1)
        p(l2)=p(l2+1)
        phi(l2)=phi(l2+1)
        phi2(l2)=phi2(l2+1)
        r(l2)=r(l2+1)
        re(l2)=re(l2+1)
        rho(l2)=rho(l2+1)
        sx(l2)=sx(l2+1)
        t(l2)=t(l2+1)
        u(l2)=u(l2+1)
        x(l2)=x(l2+1)
        y(l2)=y(l2+1)
  670   ma(l2)=ma(l2+1)
  680 RETURN
  690 RETURN
      END  Subroutine Combo

      SUBROUTINE cputin (ident, ikine, nds)
!       ---------------------------------------------------------------
      DIMENSION ident(25), izd(5)
      COMMON /chem1/ zid(5),irr(30),irt(30),rc(30,3),irrr(30,5),av,     &
     &cm(21,21),fi(21),wp(21),wm(21),wdot(21,30),csi(21),qx(21)
      av=6.03e23
      DO 40 i=1,ikine
        READ (5,50) (izd(j),j=1,5),irr(i),irt(i),(rc(i,kr),kr=1,3)
        WRITE (6,60) i,(izd(j),j=1,5),irr(i),irt(i),(rc(i,kr),kr=1,3)
        DO 30 j=1,5
          irrr(i,j)=0
          DO 20 l=1,nds
            IF (izd(j)-ident(l)) 20,10,20
   10       irrr(i,j)=l
   20     CONTINUE
   30   CONTINUE
        rc(i,1)=rc(i,1)*av
   40 END DO
   50 FORMAT (a4,3x,a4,10x,a4,3x,a4,3x,a4,9x,i2,i1,e8.2,f4.1,f9.1)
   60 FORMAT (2x,i2,2x,a4,3x,a4,10x,a4,3x,a4,3x,a4,9x,i2,i1,e8.2,f4.1,  &
     &f9.1)
      RETURN
      END  Subroutine Cputin

      SUBROUTINE drag
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      COMMON /cff/ ipto
      kzz=n
      n=k-1
      gamp1=gama+1.
      sqgam=sqrt(gama)
      sqtpi=1.7724
      DO 50 j=1,npg
        IF (n.gt.nbl(j)) GO TO 40
        axx1=((u(k)-w(n,j))**2+v(n,j)**2)
        axx=sqrt(axx1)
        um(k)=0.0
        DO 10 i=1,nds
   10     um(k)=um(k)+xs(i,k)*mu0(i)*(t(k)/t0(i))**omega(i)
!        CALCULATE DRAG COEFFICIENT
        rep(n,j)=rho(k)*axx*rp(j)/um(k)
        map=ma(k)*axx/u(k)
        rmr=map/rep(n,j)
        cdo=24./rep(n,j)
        IF (rep(n,j).ge.60.) cdo=0.40
        targ=2.*alog(map)
        tahy=(exp(2.*targ)-1.)/(exp(2.*targ)+1.)
        cdi=.66+.26*tahy+.17*exp(-2.5*(alog(map)/1.4)**2)
        IF (map.lt..3) cdi=.40
        cdb=1.0
        cd=cdi
        xkn=1.26*sqgam*rmr
        IF (xkn.gt.4.) GO TO 20
        IF (xkn.lt.1.e-8) GO TO 30
        xx1=xkn**0.4
        xx2=exp(1.2*xkn**0.5)
        gkn=xx1*xx2/(1.0+xx1*xx2)
        xx3=xkn**0.6
        xx4=exp(xkn)
        dre=1.0-exp(-xx3*xx4*(cdo-.4)*rep(n,j)/8.)
        cdb=gkn*dre
   20   CONTINUE
        s1=sqrt(0.5*gama)*map
        s2=s1*s1
        s3=s2*s1
        s4=s3*s1
        s32=s2*0.5
        sq=s1*sqrt(t(k)/tp(n,j))
        cdfm=exp(-s32)*(1.0+2.0*s2)/sqtpi/s3+(4.0*(s4+s2)-1.0)/2.0/s4*  &
     &   erf(s1)+0.667*sqtpi/sq
        cd=cdb*(cdfm-cdi)+cdi
   30   fq=cd/24.*rep(n,j)*fff
!       CALCULATE NUSSELT NUMBER
        xuo=2.0+0.459*(rep(n,j)**0.55)*(pq**0.33)
        xx5=5.0*(gama**1.5)*rmr/gamp1/pq*xuo
        xukdp=xuo/(1.0+xx5)
        xup=xukdp+(gamp1/gama*rep(n,j)*pq*exp(-0.5/rmr))
        gp=0.5*xup*ffg
        xx7=4.5*um(k)*fq/(rhss*rp(j)*rp(j))
        xmomv(j)=xx7*v(k-1,j)
        xmomp(j)=xx7*(u(k)-w(k-1,j))
        enep1(j)=xx7*axx1
        rfm=gama/gamp1*(2.+.67*exp(-map**2./3.))
        rcex=0.0
        IF (rmr.ge..01) rcex=exp(-.5/rmr)
        rec=.9+(rfm-.9)*rcex
        dlt=rec*axx1/(2.*hj*cp)
        tr(k)=t(k)+dlt
        enep2(j)=xx7*hj*2./3.*gp*cp*(tp(k-1,j)-tr(k))/fq/pq
!     CARD OMITTED
!     CARD OMITTED
!     CARD OMITTED
        GO TO 50
   40   CONTINUE
        xmomv(j)=0.0
        xmomp(j)=0.
        enep1(j)=0.
        enep2(j)=0.
   50 END DO
      n=kzz
      RETURN
      END  Subroutine Drag

      SUBROUTINE ental (ik)
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      COMMON /mainsh/ kount,kint
      COMMON /chem1/ zid(5),irr(30),irt(30),rc(30,3),irrr(30,5),av,     &
     &cm(21,21),fi(21),wp(21),wm(21),wdot(21,30),csi(21),qx(21)
      COMMON /inshop/ poo,too,uoo,phoo,coo(20)
      COMMON /shoput/ xshw,rshw,psiw,pshw,tshw,ushw,cshw(20)
      COMMON /intl/ xlx,iki,rex
      COMMON /turbul/ tle,tpr,eddyk,iturb,delmix
      COMMON /inpux/ itd,x0,y0,frac,kfl
      sph2mx=sin(phi2(kmax))
      cph2mx=cos(phi2(kmax))
!     ESTABLISHING LAST POINT FOR WHICH DP>0
      kkx=kfl-1
      i=k-1
   10 i=i+1
      dpy=p2(i+1)-p2(i)
      IF (dpy.le.0.0) GO TO 20
      IF (i.lt.kkx) GO TO 10
   20 ksh=(i-ik)/2+ik
!     SHOCK LOCATION
      rshw=(r2(ksh)+r2(ksh-1))*.5
      xshw=(x2(ksh)+x2(ksh-1))*.5
      IF (ksh.eq.ik) rshw=r2(ksh)
      IF (ksh.eq.ik) xshw=x2(ksh)
      kfll=kmax-kfl
      IF (kfll.gt.4) kfll=4
      IF (ishock.eq.0) kfll=4
      kfll=4
      iki=kfll+5
      WRITE (9,30) iki,delmix
   30 FORMAT (i2,2e12.5)
      smt=0.0
      smp=0.0
      smu=0.0
      DO 40 j=i,kfl
        smp=smp+p2(j)
        smt=smt+t2(j)
   40   smu=smu+u2(j)
      pshw=smp/(kfl-i+1)
      tshw=smt/(kfl-i+1)
      ushw=smu/(kfl-i+1)
      p2p1=pshw/p2(ik)
      DO 60 l=1,nds
        smc=0.0
        DO 50 j=i,kfl
   50     smc=smc+c2(l,j)
   60   cshw(l)=smc/(kfl-i+1)
      IF (iconst) 80,70,80
!     FOR CONSTANT GAMA
   70 gamp1=gama+1.
      pww=asin(sqrt((gamp1*p2p1+gama-1.)/(2.*gama*ma(ik)**2)))
      psiw=phi2(ik)-pww
!     NON CONSTANT SPECIFIC HEAT
      GO TO 110
   80 hsh=0.0
      DO 90 l=1,nds
        DO 90 jj=1,nn
   90   hsh=hsh+cshw(l)*a(l,jj)*tshw**(jj-2)
      fct=(pshw-p2(ik))*betap/(rho2(ik)*u2(ik)**2)
      gct=(hsh-h(ik))/u2(ik)**2*hj*2.
      pssw=0
      IF (2.0*fct.le.gct) GO TO 100
      pssw=fct/sqrt(2.*fct-gct)
  100 psiw=phi2(ik)-pssw
  110 IF (ipart.ne.2) GO TO 130
      ndsp=nds+1
      csh(ndsp)=0.
      cshw(ndsp)=0.
      DO 120 j=ik,kmax
  120   c2(ndsp,j)=0.
  130 IF (ipart.ne.2) ndsp=nds
      zms=0.0
      DO 140 j=1,nds
  140   zms=zms+cshw(j)/mw(j)
      zms=1./zms
      smdt=0.0
      kssh=ksh+1
      DO 150 j=kssh,kfl
  150   smdt=smdt+mdot(j)
      psss=phi2(kfl)
      rshw=sqrt(r2(kfl)**2-cos(psss)*smdt*tshw*rv/(pshw*zms*ushw*       &
     &3.14159))
      poo=p2(ik)
      too=t2(ik)
      uoo=u2(ik)
      xoo=xshw
      roo=rshw
      phoo=phi2(ik)
      DO 160 l=1,ndsp
  160   coo(l)=c2(l,ik)
      delr=(r2(kfl)-rshw)*.25
      dlex=(x2(kfl)-x2(ksh))*.25
      dephi=(phi2(kfl)-phi2(ksh))*.25
      xsl=x2(ksh)
      rsl=rshw
      phsl=phi2(ksh)
      WRITE (9,260) xsl,rsl,phsl
      DO 170 ii=1,4
        xsl=xsl+dlex
        rsl=rsl+delr
        phsl=phsl+dephi
        WRITE (9,260) xsl,rsl,phsl,pshw,tshw,ushw
        WRITE (9,270) (cshw(i),i=1,ndsp)
  170 END DO
      DO 180 i=1,nds
  180   csh(i)=cstrem(i)
      IF (ishock.eq.1) GO TO 220
  190 psh=pshw
      pspr=psh/pstrem
      fstrm=0.0
      bpb=0.0
      DO 200 i=1,nds
        fstrm=fstrm+cstrem(i)/mw(i)
  200   bpb=bpb+csh(i)*a(i,3)
      fstrm=1./fstrm
      gamm=bpb/(bpb-rcon/fstrm)
      emf=ustrem*sqrt((bpb*fstrm-rcon)/(bpb*tstrem*ru))
      emfns=((pspr-1.0)*(gamm+1.)/(2.*gamm)+1.)
      tsh=tstrem*(1.+2.*(gamm-1.)/(gamm+1.)**2*(emfns-1.)*(gamm+1./     &
     &emfns))
      ush=sqrt(2.*(bpb*(tstrem-tsh)*hj+ustrem**2/2.))
      psi=asin(sqrt(emfns)/emf)
      ptut=psh*tstrem*ush/(tsh*pstrem*ustrem)
      rsh=sqrt((rex**2-ptut*r2(kmax)**2/cph2mx)/(1.-ptut/cph2mx))
      xsh=x2(kmax)-(rsh-r2(kmax))*sph2mx/cph2mx
      delr=(rsh-r2(kmax))*.25
      dlex=(xsh-x2(kmax))*.25
      zzx=(emf*sin(psi))**2
      z=(zzx-1.)/tan(psi)/((gamm+1.)/2.*emf**2-zzx+1.)
      zta=atan(z)
      dephi=(zta-phi2(kmax))*.25
      xsl=x2(kfl)
      rsl=r2(kfl)
      phsl=phi2(kfl)
      DO 210 ii=1,4
        xsl=xsl+dlex
        rsl=rsl+delr
        phsl=phsl+dephi
        WRITE (9,260) xsl,rsl,phsl,psh,tsh,ush
        WRITE (9,270) (csh(i),i=1,ndsp)
  210 END DO
      GO TO 250
  220 kfll=kmax-kfl
      IF (kfll.eq.0) GO TO 190
      kfll=4
      kfl1=kfl+1
      r2kmx=(x2(kmax)-x2sh)*tan(psi-phi2(kmax))*cph2mx+r2(kmax)
      x2kmx=x2(kmax)-(x2(kmax)-x2sh)*tan(psi-phi2(kmax))*sph2mx
      xsh=x2kmx
      rsh=r2kmx
      psh=p2(kmax)
      tsh=t2(kmax)
      ush=u2(kmax)
!     ..... A GROUP OF CARDS WAS TAKEN FROM HERE .....
      delr=(r2kmx-r2(kfl))*.25
      dlex=(x2kmx-x2(kfl))*.25
      dephi=(phi2(kmax)-phi2(kfl))*.25
      xsl=x2(kfl)
      rsl=r2(kfl)
      phsl=phi2(kfl)
      smt=0.0
      smp=0.0
      smu=0.0
      DO 230 j=kfl1,kmax
        smp=smp+p2(j)
        smt=smt+t2(j)
  230   smu=smu+u2(j)
      tsh=smt/(kmax-kfl)
      ush=smu/(kmax-kfl)
      psh=smp/(kmax-kfl)
      DO 240 ii=1,4
        xsl=xsl+dlex
        rsl=rsl+delr
        phsl=phsl+dephi
        WRITE (9,260) xsl,rsl,phsl,psh,tsh,ush
        WRITE (9,270) (csh(i),i=1,ndsp)
  240 END DO
  250 xww=xshw
      xwy=xww+cos(phi2(ksh))*1.e-2*xww
      rww=rshw
      rwy=rww+sin(phi2(ksh))*(1.e-2)*rww
      WRITE (9,260) xww,rww,xwy,rwy
      xbb=xsh
      xbc=xsh+cph2mx*(1.e-2)*xbb
      rbb=rsh
      rbc=r2(kmax)+sph2mx*(1.e-2)*rbb
      WRITE (9,260) xbb,rbb,xbc,rbc
      WRITE (9,260) xsh,rsh,psi,psh,tsh,ush
      WRITE (9,270) (csh(i),i=1,ndsp)
  260 FORMAT (6e12.4)
      WRITE (9,260) xshw,rshw,psiw,pshw,tshw,ushw
      WRITE (9,270) (cshw(l),l=1,ndsp)
      WRITE (9,260) xoo,roo,phoo,poo,too,uoo
      WRITE (9,270) (coo(l),l=1,ndsp)
  270 FORMAT (8e10.3)
      RETURN
      END Subroutine Ental

      FUNCTION erf (x)
!     ----- THIS CALCULATES THE ERROR FUNCTION AND ITS COMPLEMENT -----
      DATA epsiln/1.0e-10/,n/100/,const/1.1283791679550/
      kompl=+1
      GO TO 10
      ENTRY  erfc
      kompl=-1
   10 IF (x) 20,110,30
   20 z=-x
      GO TO 40
   30 z=x
   40 IF (z-5.) 50,50,90
   50 zsq=z*z
      term=const*z*exp(-zsq)
      erf=term
      i=0
      r=1.0
   60 i=i+1
      r=r+2.0
      term=2.0*zsq*term/r
      erf=erf+term
      IF (i-n) 70,80,80
   70 IF (term-epsiln*erf) 80,60,60
   80 IF (erf-1.0) 100,100,90
   90 erf=1.0
  100 IF (x) 120,110,130
  110 erf=0.0
      GO TO 130
  120 erf=-erf
  130 IF (kompl) 140,150,150
!     ----- COMPUTE COMPLEMENTARY ERROR FUNCTION ERF -----
  140 erf=1.0-erf
  150 RETURN
      END Function Erf

      SUBROUTINE flux
!     THIS SUBROUTINE CALCULATES THE FLUXES OF MASS, MOMENTUM, AND HEAT
!     ASSUMING TURBULENT TRANSPORT
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /inpux/ itd,x0,y0,frac,kfl
      DIMENSION cdenom(20), pdenom(20), tdenom(20)
      IF (nbound) 10,30,10
   10 q(k)=0.0
      taw(k)=0.0
      DO 20 i=1,nds
   20   zj(i,k)=0.0
      GO TO 150
!     ----- SETUP OF AVERAGED FLOW QUANTITIES NEEDED -----
   30 ybbar=yabar
      tbbar=tabar
      rhbbar=rhabar
      DO 40 i=1,nds
        cbbar(i)=cabar(i)
   40   xbbar(i)=xabar(i)
      pbbar=pabar
      IF (iflag) 50,80,50
   50 yabar=.5*(y(k+1)+y2(k+1))
      tabar=.5*(t(k+1)+t2(k+1))
      ku=k+1
      ubar(ku)=sqrt(.5*(u(k+1)**2+u2(k+1)**2))
      rhabar=.5*(rho(k+1)+rho2(k+1))
      rbvg=0.25*(r(k)+r(k-1)+r2(k)+r2(k-1))
      ravg=0.25*(r(k+1)+r(k)+r2(k+1)+r2(k))
      zmw2=0.0
      DO 60 i=1,nds
        cabar(i)=.5*(c(i,k+1)+c2(i,k+1))
   60   zmw2=zmw2+cabar(i)/mw(i)
      zmw2=1.0/zmw2
      DO 70 i=1,nds
   70   xabar(i)=cabar(i)*zmw2/mw(i)
      pabar=.5*(p(k+1)+p2(k+1))
      GO TO 100
   80 tabar=t(k+1)
      ku=k+1
      ubar(ku)=u(k+1)
      rhabar=rho(k+1)
      rbvg=0.5*(r(k)+r(k-1))
      ravg=0.5*(r(k+1)+r(k))
      DO 90 i=1,nds
        cabar(i)=c(i,k+1)
   90   xabar(i)=xs(i,k+1)
      pabar=p(k+1)
  100 CALL transp (tabar, pabar, k+1)
!     ----- CALCULATION OF HEAT FLUX AND SHEAR -----
      denomu=(mdot(k+1)+mdot(k))/(rhabar*ubar(k+1)*ravg**delta*mu+      &
     &rhbbar*ubar(k)*rbvg**delta*mu2)/(2.0**delta*pi)
      denomt=(mdot(k+1)+mdot(k))/(rhabar*ubar(k+1)*ravg**delta*kay+     &
     &rhbbar*ubar(k)*rbvg**delta*kay2)/(2.0**delta*pi)
      taw(k)=(ubar(k+1)-ubar(k))/denomu
      q(k)=(tabar-tbbar)/denomt
!     ----- CALCULATION OF MASS FLUX -----
      DO 110 i=1,nds
        cdenom(i)=(mdot(k)+mdot(k+1))/(rhbbar**2*d2ih(1,2)*ubar(k)*pi*  &
     &   (2.0*rbvg)**delta+rhabar**2*dih(1,2)*ubar(k+1)*pi*(2.0*ravg)** &
     &   delta)
  110   zj(i,k)=(cabar(i)-cbbar(i))/cdenom(i)
      IF (ishock.ne.1) GO TO 120
      IF (kmax.gt.kfl) GO TO 120
      x0=x(1)
  120 l=1
      dist=x(l)-x0
      IF (abs(dist)-frac*y0) 130,150,150
  130 flimit=abs(dist)/(frac*y0)
      q(k)=q(k)*flimit
      taw(k)=taw(k)*flimit
      DO 140 i=1,nds
  140   zj(i,k)=zj(i,k)*flimit
  150 RETURN
      END  Subroutine Flux

      SUBROUTINE locate (xb, yb, x1, y1, x2, y2, x3, y3, ifgood)
      DIMENSION xloc(3)
!     IFGOOD = 1  POINT IS BOUNDED
!            = 2  POINT IS NOT BOUNDED
      DO 60 i=1,3
        GO TO (10,20,30),i
   10   xx2=x1
        yy2=y1
        GO TO 40
   20   xx2=x2
        yy2=y2
        GO TO 40
   30   xx2=x3
        yy2=y3
   40   d=sqrt((xx2-xb)**2+(yy2-yb)**2)
        cost=(xx2-xb)/d
        IF (yy2.ge.yb) GO TO 50
        xloc(i)=3.+cost
        GO TO 60
   50   xloc(i)=1.-cost
   60 END DO
      IF (xloc(3).lt.xloc(2)) GO TO 90
      IF (xloc(1).ge.xloc(2).and.xloc(1).le.xloc(3)) GO TO 80
   70 ifgood=2
      GO TO 100
   80 ifgood=1
      GO TO 100
   90 IF (xloc(1).ge.xloc(2).or.xloc(1).le.xloc(3)) GO TO 80
      GO TO 70
  100 RETURN
      END  Subroutine Locate

      SUBROUTINE mate7 (idim, nn, mm, a, b, det, ind, nogo)
!SIMEQ           SOLVES THE EQUATION AX=B WHERE A,X, AND B ARE MATRICES
      DIMENSION a(idim,1), b(idim,1), ind(1)
      n=nn
      m=mm
!        SET INITIAL CONSTANTS
      nogo=1
!        SPECIAL CONSIDERATION WHEN THE ORDER OF MATRIX A IS 1.
      IF (n.ne.1) GO TO 20
      IF (a(1,1).eq.0.) GO TO 170
      hold=a(1,1)
      DO 10 i=1,m
   10   a(1,i)=b(1,i)/hold
      det=det*hold
      RETURN
   20 nm1=n-1
      np1=n+1
!        INITIALIZE DETERMINANT
      detm=det
!        INITIALIZE COLUMN INDICATORS
      DO 30 i=1,n
   30   ind(i)=i
!        BEGIN TRIANGULARIZATION TO GET UPPER TRIANGLE
      DO 130 k=1,nm1
        kp1=k+1
        kr=k
        kc=k
!        SEARCH FOR PIVOTAL ELEMENT
        biga=abs(a(k,k))
        DO 40 i=k,n
          DO 40 j=k,n
          IF (biga.ge.abs(a(i,j))) GO TO 40
          biga=abs(a(i,j))
          kr=i
          kc=j
   40   CONTINUE
!        TEST FOR SINGULAR MATRIX
        IF (biga.eq.0.) GO TO 170
!        UPDATE DETERMINANT
        detm=detm*a(kr,kc)
!        INTERCHANGE ROWS
        IF (kr.eq.k) GO TO 70
        DO 50 i=k,n
          hold=a(k,i)
          a(k,i)=a(kr,i)
   50     a(kr,i)=hold
!        INTERCHANGE ELEMENTS OF RIGHT HAND SIDES
        DO 60 l=1,m
          hold=b(k,l)
          b(k,l)=b(kr,l)
   60     b(kr,l)=hold
!        CHANGE SIGN OF DETERMINANT DUE TO ROW INTERCHANGE
        detm=-detm
!        INTERCHANGE COLUMNS
   70   IF (kc.eq.k) GO TO 90
        DO 80 j=1,n
          hold=a(j,k)
          a(j,k)=a(j,kc)
   80     a(j,kc)=hold
!        INTERCHANGE COLUMN INDICATORS
        i=ind(k)
        ind(k)=ind(kc)
        ind(kc)=i
!        CHANGE SIGN OF DETERMINANT DUE TO COLUMN INTERCHANGE
        detm=-detm
!        DIVIDE REDUCED EQUATION-ON, BY LEADING ELEMENT
   90   DO 100 i=kp1,n
  100     a(k,i)=a(k,i)/a(k,k)
        DO 110 l=1,m
  110     b(k,l)=b(k,l)/a(k,k)
!        REDUCE MATRIX AND RIGHT HAND SIDES
        DO 130 i=kp1,n
        DO 120 j=kp1,n
  120     a(i,j)=a(i,j)-a(i,k)*a(k,j)
        DO 130 l=1,m
  130   b(i,l)=b(i,l)-a(i,k)*b(k,l)
!        FINAL TEST FOR SINGULAR MATRIX
      IF (a(n,n).eq.0.) GO TO 170
!        COMPUTE FINAL DETERMINANT
      det=detm*a(n,n)
!        BACK SUBSTITUE TO OBTAIN SOLUTION VECTORS
      DO 150 l=1,m
        b(n,l)=b(n,l)/a(n,n)
        DO 150 i=1,nm1
        hold=0.
        j=n-i
        DO 140 kc=1,i
          k=np1-kc
  140     hold=hold+a(j,k)*b(k,l)
  150   b(j,l)=b(j,l)-hold
!     REARRANGE SOLUTION VECTORS TO ORIGINAL ORDER
      DO 160 i=1,n
        j=ind(i)
        DO 160 l=1,m
  160   a(j,l)=b(i,l)
      RETURN
!        SINGULAR MATRIX - REDUNDANT SET OF EQUATIONS
  170 nogo=3
      det=0.
      RETURN
      END  Subroutine Mate7

      SUBROUTINE orthog (side, lim1, lim2)
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      INTEGER side,dirctn,prop1,prop2
      COMMON /intrp1/ nstrm(2),tsfc(2),psfc(2),usfc(2),phisfc(2),csfc(2,&
     &25),dptsfc(2)
      COMMON /orshok/ xtry,rtry
      ax=fx(side,lim2)-fx(side,lim1)
      bx=fr(side,lim2)-fr(side,lim1)
      cx=ax*xtry+bx*rtry
      IF (ax.eq.0.) GO TO 10
      temp=-bx/ax*fx(side,lim1)+fr(side,lim1)
      temp1=ax+bx*bx/ax
      xint=(cx-bx*temp)/temp1
      rint=temp+bx/ax*xint
      GO TO 20
   10 xint=fx(side,lim1)
      rint=(cx-ax*xint)/bx
   20 dptsfc(side)=sqrt((xint-xtry)**2+(rint-rtry)**2)
      d1sq=(xint-fx(side,lim1))**2+(rint-fr(side,lim1))**2
      xmid1=(fx(side,lim1)+fx(side,lim2))/2.
      rmid1=(fr(side,lim1)+fr(side,lim2))/2.
      prop1=lim2
      d2sq=(xmid1-fx(side,lim1))**2+(rmid1-fr(side,lim1))**2
      dirctn=1
      IF (d2sq.lt.d1sq) dirctn=2
      GO TO (30,40),dirctn
   30 j=lim1
      IF (j.eq.1) GO TO 50
      j=j-1
      xmid2=(fx(side,lim1)+fx(side,j))/2.
      rmid2=(fr(side,lim1)+fr(side,j))/2.
      prop2=lim1
      GO TO 60
   40 j=lim2
      IF (j.eq.nstrm(side)) GO TO 50
      j=j+1
      xmid2=(fx(side,lim2)+fx(side,j))/2.
      rmid2=(fr(side,lim2)+fr(side,j))/2.
      prop2=j
      GO TO 60
   50 prop2=prop1
      IF (j.eq.1) GO TO 60
      GO TO 70
   60 dptmd1=sqrt((xint-xmid1)**2+(rint-rmid1)**2)
      dmd12=sqrt((xmid1-xmid2)**2+(rmid1-rmid2)**2)
      prprat=dptmd1/dmd12
      usfc(side)=fu(side,prop1)+prprat*(fu(side,prop2)-fu(side,prop1))
      psfc(side)=fp(side,prop1)+prprat*(fp(side,prop2)-fp(side,prop1))
   70 tsfc(side)=ft(side,prop1)+prprat*(ft(side,prop2)-ft(side,prop1))
      DO 80 i=1,nds
   80   csfc(side,i)=fc(side,prop1,i)+prprat*(fc(side,prop2,i)-fc(side, &
     &   prop1,i))
      d3sq=(fx(side,lim2)-fx(side,lim1))**2+(fr(side,lim2)-fr(side,lim1)&
     &)**2
      phisfc(side)=fphi(side,lim1)+(fphi(side,lim2)-fphi(side,lim1))*   &
     &sqrt(d1sq/d3sq)
      RETURN
      END  Subroutine Orthog

      SUBROUTINE part
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      k=k-1
      kzz=n
      n=k+1
      del=(dls(n)+dls(n-1))*.5
      DO 70 j=1,npg
        IF (k.gt.nbl(j)) GO TO 60
        xx6=del*v(k,j)/w(k,j)
        IF (k.eq.1) GO TO 10
        k1=k-1
        rx1=(.5*(r(k1+1)+r(k1)))**delta
        rx2=(.5*(r(k1+2)+r(k1+1)))**delta
        delly=sqrt((x(k1+2)-x(k1))**2+(r(k1+2)-r(k1))**2)*.5
        then=(phi(k1+2)-phi(k1))/delly*.5
        dwdn=(w(k1+1,j)-w(k1,j))/delly
        dvdn=(v(k1+1,j)-v(k1,j))/delly
        dtdn=(tp(k1+1,j)-tp(k1,j))/delly
        drrvdn=(rx2*rhp(k1+1,j)*v(k1+1,j)-rx1*rhp(k1,j)*v(k1,j))/delly
        GO TO 20
   10   dwdn=0.0
        delly=sqrt(r(2)**2+(x(2)-x(1))**2)*.5
        dvdn=v(1,j)/delly
        dtdn=0.0
        drrvdn=rhp(1,j)*v(1,j)/delly*(r(2)/2.)**delta
        then=phi(2)/(2.*delly)
   20   thes=(phi2(n)+phi2(n-1)-phi(n)-phi(n-1))/(dls(n)+dls(n-1))
        w2(k,j)=w(k,j)-xx6*dwdn+xx6*v(k,j)*then+del*v(k,j)*thes+xmomp(j)&
     &   *del/w(k,j)
        v2(k,j)=v(k,j)-xx6*dvdn-del*v(k,j)*then-w(k,j)*del*thes-xmomv(j)&
     &   *del/w(k,j)
!      ENERGY EQUATION
        cc=cl
        IF (tp(k,j).le.tps) cc=cs
        IF (icond(k,j).eq.0) GO TO 30
        xx8=wtp(j)*enep2(j)*del/w(k,j)/hj
        xx9=4.0*pi*rp(j)*rp(j)*ep*sig*(tps**4.0)/w(k,j)*del/hj
        deng(k,j)=deng(k,j)+xx8+xx9
        IF (abs(deng(k,j)).lt.hc(j)) GO TO 50
        icond(k,j)=0
        WRITE (6,80) rp(j),s,k,deng(k,j)
   30   CONTINUE
        xx8=enep2(j)/hj/w(k,j)/cc
        xx9=ep*sig*3.0/(w(k,j)*rp(j)*rhss*cc)*(tp(k,j)**4.0)/hj
        tp2(k,j)=tp(k,j)-xx6*dtdn-(xx8+xx9)*del
        IF (((tp2(k,j)-tps)*(tp(k,j)-tps)).ge.0.0) GO TO 50
        icond(k,j)=1
        tp2(k,j)=tps
        WRITE (6,40) rp(j),s,k
   40   FORMAT (' ',1pe10.3,' CM PARTICLE  GROUP START SOLIDIFYING AT S=&
     &',0pf10.4,' CM STREAMTUBE NO. ',i5/)
   50   CONTINUE
!        CONTINUITY EQUATION
        rx1=(.5*(r(k+1)+r(k)))**delta
        r2x1=(.5*(r2(k+1)+r2(k)))**delta
        xx10=rx1*rhp(k,j)*w(k,j)
        xx11=rx1*rhp(k,j)*v(k,j)
        rhp2(k,j)=1.0/(r2x1*w2(k,j))*(xx10*(1.0-del*then)-drrvdn*del+   &
     &   xx11*thes*del)
        GO TO 70
   60   v2(k,j)=0.0
        w2(k,j)=0.0
        tp2(k,j)=0.0
        rhp2(k,j)=0.0
        rep(k,j)=0.0
   70 END DO
   80 FORMAT (' ',1pe10.3,' CM SIZE PARTICLE GROUP SOLIDIFIED AT S=',   &
     &0pf10.4,' CM  STREAMTUBE NO.',i5,3x,1pe10.3,'CAL'/)
      k=k+1
      n=kzz
      RETURN
      END  Subroutine Part

      SUBROUTINE pbndy
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      kmx=0
      DO 60 j=1,npg
        nblj=nbl(j)
        nbl1=nblj+1
        dndsp(j)=dndsp(j)+dls(nbl1)*v2(nblj,j)/w2(nblj,j)
        trbdy(j)=y2(nbl1)+dndsp(j)
        DO 10 iem=1,kmax
          IF (trbdy(j).lt.y2(iem+1)) GO TO 20
   10   CONTINUE
        IF (par(1).eq.0.0) WRITE (6,80)
   20   CONTINUE
        IF (iem.ne.nbl(j)) dndsp(j)=trbdy(j)-y2(iem+1)
        IF (iem.gt.kmx) kmx=iem
        DO 50 lk=1,iem
          IF (lk.gt.nbl(j)) GO TO 30
          GO TO 40
   30     w2(lk,j)=w2(lk-1,j)
          v2(lk,j)=v2(lk-1,j)
          tp2(lk,j)=tp2(lk-1,j)
          rhp2(lk,j)=rhp2(lk-1,j)
   40     CONTINUE
          v(lk,j)=v2(lk,j)
          w(lk,j)=w2(lk,j)
          tp(lk,j)=tp2(lk,j)
          rhp(lk,j)=rhp2(lk,j)
   50   CONTINUE
        nbl(j)=iem
   60 END DO
      IF (par(1).eq.0.0) WRITE (6,70) s,(trbdy(j),j=1,npg)
   70 FORMAT (' ',' BOUNDARY OF PARTICLE PHASE AT S=',f11.5,1p8e11.3/)
   80 FORMAT ('  -------NO PARTICLES IN THIS FLOW FIELD '//)
      RETURN
      END  Subroutine Pbndy

!+
SUBROUTINE Pdump(a,n,nn)
! ---------------------------------------------------------------------------
! PURPOSE - ??
  REAL,INTENT(IN):: a
  INTEGER,INTENT(IN):: n,nn
!----------------------------------------------------------------------------
  WRITE(*,*) 'Pdump called'
  WRITE(*,*) 'a=', a
  WRITE(*,*) 'n=', n
  WRITE(*,*) 'nn=', nn
  RETURN
END Subroutine Pdump   ! ----------------------------------------------------


      SUBROUTINE pptout (kkkk)
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      DIMENSION ptp(8)
      IF (kkkk.eq.1) GO TO 10
      lki=k-1
      lkmax=k-1
      GO TO 20
   10 lki=1
      lkmax=kmx+2
   20 CONTINUE
      DO 40 lk=lki,lkmax
        lk1=lk+1
        ptpp=0.0
        DO 30 j=1,npg
          ptp(j)=rhp2(lk,j)*w2(lk,j)**2/(2.*betap)
          ptpp=ptpp+ptp(j)
   30   CONTINUE
        WRITE (6,60) lk,s,u2(lk1),t2(lk1),rho(lk1),re(lk1),y2(lk1),     &
     &   r2(lk1)
        WRITE (6,50)
        WRITE (6,120) (j,j=1,npg)
        WRITE (6,70) (w2(lk,j),j=1,npg)
        WRITE (6,80) (v2(lk,j),j=1,npg)
        WRITE (6,90) (rep(lk,j),j=1,npg)
        WRITE (6,100) (tp2(lk,j),j=1,npg)
        WRITE (6,110) (rhp2(lk,j),j=1,npg)
        WRITE (6,130) ptpp
   40 END DO
   50 FORMAT (/' ','----- PARTICLE PHASE PROPERTIES -----  '/)
   60 FORMAT (//' ','---STREAMTUBE= ',i5,'  S=',f12.5///,2x,'-- LOCAL GA&
     &S PROPERTIES'//,2x,'U=',1pe10.3,2x,'T=',e10.3,2x,'DENSITY=',e10.3,&
     &2x,'REYNOLDS NO.=',e10.3,2x,'Y2=',e10.3,2x,'R2=',e10.3)
   70 FORMAT (' ','DOWNSTREAM VELOCITY',2x,1p8e12.3)
   80 FORMAT (' ','CROSS-STREAM VELOCITY',1p8e12.3)
   90 FORMAT (' ','PARTICLE REYNOLDS NO.',1p8e12.3)
  100 FORMAT (' ','PARTICLE TEMPERATURE',1x,1p8e12.3)
  110 FORMAT (' ','PARTICLE DENSITY',5x,1p8e12.3)
  120 FORMAT (' ','PARTICLE GROUP',8(7x,i5))
  130 FORMAT (','' PARTICLE MOMENTUM FLUX ',2x,1p8e12.3)
      RETURN
      END  Subroutine PptOut

      SUBROUTINE pputin
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      DIMENSION vi(8), wi(8), tpi(8), rhpi(8), dengi(8)
      NAMELIST / nam1/fff,ffg,umi,cl,cs,tps,rhss,wt,htran,sig,ep,npg,nc,&
     &kmx,rp,wi,vi,tpi,rhpi
      READ (5,160) fff,ffg,cl,cs,htran,wt
      READ (5,160) rhss,sig,ep,tps
      READ (5,220) npg,nc
      READ (5,220) (nbl(j),j=1,npg)
      DO 10 j=1,npg
   10   nbl(j)=nbl(j)-1
      kmx=nbl(1)
      IF (npg.le.1) GO TO 30
      DO 20 j=2,npg
        IF (kmx.lt.nbl(j)) kmx=nbl(j)
   20 END DO
   30 CONTINUE
      READ (5,160) (rp(j),j=1,npg)
      IF (nc.eq.1) GO TO 70
      READ (5,160) (wi(j),j=1,npg)
      READ (5,160) (vi(j),j=1,npg)
      READ (5,160) (tpi(j),j=1,npg)
      READ (5,160) (rhpi(j),j=1,npg)
      READ (5,160) (dengi(j),j=1,npg)
      DO 50 j=1,npg
        kma=nbl(j)
        DO 40 k=1,kma
          w(k,j)=wi(j)
          v(k,j)=(k-1)*vi(j)/(kma-1)
          tp(k,j)=tpi(j)
          rhp(k,j)=rhpi(j)
          deng(k,j)=dengi(j)
          w2(k,j)=w(k,j)
          v2(k,j)=v(k,j)
          tp2(k,j)=tp(k,j)
          rhp2(k,j)=rhp(k,j)
   40   CONTINUE
        dndsp(j)=0.0
        trbdy(j)=y(kma+1)
   50 END DO
      DO 60 j=1,npg
        kma=nbl(j)
        DO 60 k=1,kma
        icond(k,j)=0
   60 CONTINUE
      GO TO 100
   70 CONTINUE
      DO 90 k=1,kmx
        READ (5,160) (w(k,j),j=1,npg)
        READ (5,160) (v(k,j),j=1,npg)
        READ (5,160) (tp(k,j),j=1,npg)
        READ (5,160) (rhp(k,j),j=1,npg)
        READ (5,220) (icond(k,j),j=1,npg)
        READ (5,160) (deng(k,j),j=1,npg)
        DO 80 j=1,npg
          w2(k,j)=w(k,j)
          v2(k,j)=v(k,j)
          tp2(k,j)=tp(k,j)
          rhp2(k,j)=rhp(k,j)
   80   CONTINUE
   90 END DO
      READ (5,160) (trbdy(j),j=1,npg)
      READ (5,160) (dndsp(j),j=1,npg)
  100 CONTINUE
      WRITE (6,nam1)
      DO 110 j=1,npg
        volp(j)=1.333333*pi*(rp(j)**3)
        wtp(j)=volp(j)*rhss
  110   hc(j)=wtp(j)*htran
      DO 140 j=1,npg
        IF (kmax-nbl(j)) 140,140,120
  120   kma1=nbl(j)+1
        DO 130 k=kma1,kmax
          w(k,j)=0.0
          v(k,j)=0.0
          tp(k,j)=0.0
          rhp(k,j)=0.0
          w2(k,j)=0.0
          v2(k,j)=0.0
          tp2(k,j)=0.0
          rhp2(k,j)=0.0
          deng(k,j)=0.0
          icond(k,j)=0
  130   CONTINUE
  140 END DO
      DO 150 k=1,kmax
        WRITE (6,170) k
        WRITE (6,180) (w(k,j),j=1,npg)
        WRITE (6,190) (v(k,j),j=1,npg)
        WRITE (6,200) (tp(k,j),j=1,npg)
        WRITE (6,210) (rhp(k,j),j=1,npg)
  150 END DO
  160 FORMAT (8e10.3)
  170 FORMAT (//' ','K=',i5/)
  180 FORMAT (' ','DOWNSTREAM VELOCITY',2x,1p8e12.3)
  190 FORMAT (' ','CROSS-STREAM VELOCITY',1p8e12.3)
  200 FORMAT (' ','PARTICLE TEMPERATURE',1x,1p8e12.3)
  210 FORMAT (' ','PARTICLE DENSITY',5x,1p8e12.3)
  220 FORMAT (8i5)
      RETURN
      END  Subroutine PputIn

      SUBROUTINE putin
!     THIS SUBROUTINE READS IN AND INITIALIZES ALL DATA EXCEPT THAT
!     HAVING TO DO WITH AN EXTERNAL NON-UNIFORM FLOW FIELD
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      REAL massmd,l1,l2
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      COMMON /ibugs1/ ibugsh
      COMMON /chem1/ zid(5),irr(30),irt(30),rc(30,3),irrr(30,5),av,     &
     &cm(21,21),fi(21),wp(21),wm(21),wdot(21,30),csi(21),qx(21)
      COMMON /inpux/ itd,x0,y0,frac,kfl
      COMMON /trnspt/ alftd(25),dexp(25)
      COMMON /inpro/ icsh,iin
      EQUIVALENCE (ivisc,iextra(1)),(ipch,iextra(5))
      COMMON /shoput/ xshw,rshw,psiw,pshw,tshw,ushw,cshw(20)
      COMMON /inshop/ poo,too,uoo,phoo,coo(20)
      COMMON /intl/ xlx,iki,rex
      COMMON /turbul/ tle,tpr,eddyk,iturb,delmix
      COMMON /cff/ ipto
      COMMON /xlm/ xlmd
      COMMON /pergff/ pergf
      COMMON /dmd/ kmnf
      COMMON /sub/ xdis,l1,l2,phio,rstar,rmd,massmd,rnoz,thtp,thru
      COMMON /dsl/ rdsl,phdsl,xinit
      INTEGER thru
      DIMENSION ab(7)
      EQUIVALENCE (pwsh,extra(11)),(itapes,iextra(2))
      iin=0
!     INTEGRAL PARAMETERS,LIMITS,AND SWITCHES
      card1=1.
      READ (5,10) itype,ikind,ipd,iptuc,ishock,ibugsh,itd,icsh,itapes,  &
     &ipart,iki,iturb,ipto,isb
      READ (5,10) kmax,nn,kp,mmax,nmax,nds,niter,lplane,ikine,ipch,kmn, &
     &kmnf
   10 FORMAT (16i5)
      IF (ipart.eq.2) irog=1
      WRITE (6,1140) kmax,nn,kp,ip,itype,ikind,mmax,nmax,nds,ipd,niter, &
     &idiff,lplane,iptuc,ipch,ishock,ibugsh,itd,icsh,itapes,ikine,ipart,&
     &iki,iturb,ipto,isb,kmn,kmnf
      IF (iki.eq.0) GO TO 60
      REWIND  9
      REWIND  10
      DO 20 i=1,nds
   20   READ (9,160) (a(i,j),j=1,nn),csi(i)
      DO 30 i=1,nds
   30   READ (9,270) ident(i),mu0(i),t0(i),omega(i),pr(i),sc(i),mw(i)
      READ (9,40) xlmd
   40 FORMAT (e12.5)
      IF (kmn.ge.1) kmaxx=kmax
      READ (9,50) kmax,delmix
      IF (kmn.ge.1) kmax=kmaxx
   50 FORMAT (i2,2e12.5)
      REWIND  3
      IF (iki.eq.1) REWIND  4
   60 extra(3)=kmax
      kfl=kmax
      alphap=kp
   70 FORMAT (6e12.4)
!     NON-INTEGRAL PARAMETERS
      READ (5,100) alphah,epslon,tol,delta,atol,fstep,grad,frac,fractn
      WRITE (6,1150) alphah,epslon,tol,delta,atol,fstep,grad,frac,      &
     &fractn
      READ (5,80) (par(i),i=1,10)
   80 FORMAT (3f5.1,4e10.2,f5.1,2e10.2)
      WRITE (6,90) (par(i),i=1,10)
   90 FORMAT (' ',/5x,'PAR(1)',6x,'PAR(2)',6x,'PAR(3)',6x,'PAR(4)',6x,'P&
     &AR(5)',6x,'PAR(6)',6x,'PAR(7)',6x,'PAR(8)',6x,'PAR(9)',6x,'PAR(10)&
     &',/5x,f5.1,7x,f5.1,7x,f5.1,4x,3(e10.4,2x),1x,e10.4,4x,f5.1,3x,    &
     &e10.4,4x,e10.4)
  100 FORMAT (8f9.4,1f8.4)
      extra(4)=fractn
      IF (itd) 110,130,110
  110 READ (5,360) (alftd(i),i=1,nds)
      READ (5,360) (dexp(i),i=1,nds)
      WRITE (6,1160)
      DO 120 i=1,nds
  120   WRITE (6,1170) i,alftd(i),dexp(i)
  130 pi=3.14159265
      ru=rcon*hj
      rv=ru/betap
!     ENTHALPY POLYNOMIAL CONSTANTS
      IF (iki.ne.0) GO TO 150
      DO 140 i=1,nds
        READ (5,160) (a(i,j),j=1,nn),csi(i)
  140   WRITE (9,160) (a(i,j),j=1,nn),csi(i)
  150 CONTINUE
  160 FORMAT (4e13.5/3e13.5)
      DO 170 jk=1,6
  170   ab(jk)=0.0
!     TRANSPORT PROPERTIES
      ivisc=0
      WRITE (6,1180)
      DO 200 i=1,nds
        IF (iki.ne.0) GO TO 180
        READ (5,270) ident(i),mu0(i),t0(i),omega(i),pr(i),sc(i),mw(i)
        WRITE (9,270) ident(i),mu0(i),t0(i),omega(i),pr(i),sc(i),mw(i)
  180   WRITE (6,1190) ident(i),mu0(i),t0(i),omega(i),pr(i),sc(i),mw(i)
        DO 190 jk=1,6
          ab(jk)=a(i,jk)
  190   CONTINUE
        a(i,1)=-ab(5)*1.0e+6/mw(i)
        a(i,2)=ab(6)/mw(i)*1000.0
        a(i,3)=ab(1)/mw(i)
        a(i,4)=ab(2)/(2.0e3*mw(i))
        a(i,5)=ab(3)/(3.0e6*mw(i))
        a(i,6)=ab(4)/(4.0e9*mw(i))
        csi(i)=csi(i)/mw(i)
  200 END DO
      IF (irog.eq.1.and.ipart.eq.0) GO TO 210
      GO TO 230
  210 nds=nds+1
      a(nds,1)=cs
      DO 220 i=2,nn
  220   a(nds,i)=0.
      csi(nds)=0.
      ident(nds)=alox
      mu0(nds)=0.
      t0(nds)=300.
      omega(nds)=0.75
      pr(nds)=1.4
      sc(nds)=1.4
      mw(nds)=1.e+10
  230 CONTINUE
      DO 260 j=4,nn
        DO 260 i=1,nds
        IF (a(i,1)) 250,240,250
  240   IF (a(i,j)) 250,260,250
  250   iconst=1
  260 CONTINUE
  270 FORMAT (a4,8x,5e12.4,1e8.2)
      delmix=0.0
      IF (iturb.eq.0) GO TO 280
      READ (5,840) tle,tpr,eddyk,pergf
      WRITE (6,1200) tle,tpr,eddyk,pergf
      ivisc=1
      GO TO 310
  280 DO 300 i=1,nds
        IF (mu0(i)) 290,300,290
  290   ivisc=1
  300 END DO
      IF (iturb.eq.0.and.ipart.ne.0) ivisc=0
  310 mu=10.0**(-20)
      mu2=10.0**(-20)
      IF (kmn.gt.1) GO TO 900
!     INNER STREAMLINE POSITION
      k=1
      IF (iki.eq.0) READ (5,70) x(k),r(k),phi(k)
      WRITE (6,1210)
      IF (iki.ne.0) READ (9,70) x(k),r(k),phi(k)
      WRITE (6,1220) k,x(k),r(k),phi(k)
      x0=x(k)
      y(k)=0.0
!     OTHER STREAMLINE POSITION AND STREAMTUBE PROPER
      k=2
  320 IF (iki.eq.0) READ (5,70) x(k),r(k),phi(k),p(k),t(k),u(k)
      IF (iki.ne.0) READ (9,70) x(k),r(k),phi(k),p(k),t(k),u(k)
      WRITE (6,1220) k,x(k),r(k),phi(k),p(k),t(k),u(k)
      b=sqrt((x(k)-x(k-1))**2+(r(k)-r(k-1))**2)
      IF (abs(phi(k)-phi(k-1))-1.0e-06) 330,330,340
  330 dely(k)=b
      GO TO 350
  340 dely(k)=(phi(k)-phi(k-1))*b/(2.0*sin(.5*(phi(k)-phi(k-1))))
  350 aa(k)=pi*(r(k)+r(k-1))**delta*dely(k)
      y(k)=y(k-1)+dely(k)
!     STREAMTUBE COMPOSITION (SPECIES MASS FRACTIONS)
      IF (iki.eq.0) READ (5,360) (c(i,k),i=1,nds)
      IF (iki.ne.0) READ (9,360) (c(i,k),i=1,nds)
  360 FORMAT (8e10.3)
!     SETUP CALCULATIONS
      zmw=0.0
      DO 370 i=1,nds
  370   zmw=zmw+c(i,k)/mw(i)
      zmw=1.0/zmw
      zmwk(k)=zmw
      DO 380 i=1,nds
  380   xs(i,k)=c(i,k)*zmw/mw(i)
      rho(k)=zmw*p(k)/(rv*t(k))
      rhabar=rho(k)
      mdot(k)=rho(k)*u(k)*aa(k)
      fmax=1.0
      kay=0.0
      DO 400 i=1,nds
        DO 390 j=1,nn
          dih(i,j)=0.0
  390   CONTINUE
  400 END DO
      IF (ivisc) 430,410,430
  410 DO 420 i=1,nds
        cps(i)=0.0
        DO 420 j=1,nn
  420   cps(i)=cps(i)+float(j-2)*a(i,j)*t(k)**(j-3)
      GO TO 440
  430 CALL transp (t(k), p(k), k)
  440 cp=0.0
      DO 450 i=1,nds
  450   cp=cp+cps(i)*c(i,k)
      IF (kay-mu*cp) 470,470,460
  460 fmax=kay/(cp*mu)
  470 DO 490 i=1,nds
        max=i
        DO 490 j=1,max
        IF (dih(i,j)*rho(k)/mu-fmax) 490,490,480
  480   fmax=dih(i,j)*rho(k)/mu
  490 CONTINUE
      re(k)=rho(k)*u(k)*dely(k)/(mu*fmax)
      ma(k)=u(k)*sqrt((cp*zmw-rcon)/(cp*t(k)*ru))
      gama=cp/(cp-rcon/zmw)
      pq=4.0*gama/(9.*gama-5.0)
      h(k)=0.0
      DO 500 i=1,nds
        DO 500 j=1,nn
  500   h(k)=h(k)+a(i,j)*c(i,k)*t(k)**(j-2)
      IF (k-kmax) 510,520,520
  510 k=k+1
      GO TO 320
  520 y0=y(kmax)
      rex=r(kmax)
      IF (iki.eq.0) rnoz=rex
!     WALL CONDITIONS
      IF (iki.ne.0) GO TO 540
      WRITE (6,1230)
      DO 530 i=1,mmax
        READ (5,70) xw(i),rw(i),phiw(i),pw(i),sw(i)
  530   WRITE (6,1220) i,xw(i),rw(i),phiw(i),pw(i),sw(i)
      GO TO 550
  540 READ (9,70) xw(1),rw(1),xw(2),rw(2)
  550 swa=0.0
      m=1
      pwsh=pw(m)
      iphiw=phiw(m)
      IF (iphiw) 570,560,570
  560 phiw(m)=phi(m)
  570 m=m+1
      iphiw=phiw(m)
      IF (iphiw) 610,580,610
  580 IF (m-mmax) 590,600,600
  590 phiw(m)=atan((rw(m+1)-rw(m-1))/(xw(m+1)-xw(m-1)))
      GO TO 610
  600 phiw(m)=2.0*atan((rw(m)-rw(m-1))/(xw(m)-xw(m-1)))
      phiw(m)=phiw(m)-phiw(m-1)
  610 isw=sw(m)
      IF (isw) 660,620,660
  620 fkw=(sin(phiw(m))-sin(phiw(m-1)))/(xw(m)-xw(m-1))
      IF (abs(fkw/(xw(m)-xw(m-1)))-1.e-06) 640,640,630
  630 delswa=(phiw(m)-phiw(m-1))/fkw
      GO TO 650
  640 delswa=sqrt((xw(m)-xw(m-1))**2+(rw(m)-rw(m-1))**2)
  650 swa=swa+delswa
      sw(m)=swa
  660 IF (m-mmax) 570,670,670
  670 IF (iki.ne.0) GO TO 690
      WRITE (6,1240)
      DO 680 n=1,nmax
        READ (5,70) xb(n),rb(n),phib(n),pb(n),sb(n)
  680   WRITE (6,1220) n,xb(n),rb(n),phib(n),pb(n),sb(n)
      GO TO 700
  690 READ (9,70) xb(1),rb(1),xb(2),rb(2)
  700 swa=0.0
!     OUTER BOUNDARY CONDITIONS
      n=1
      iphib=phib(n)
      IF (iphib) 720,710,720
  710 phib(n)=phi(kmax)
  720 n=n+1
      iphib=phib(n)
      IF (iphib) 760,730,760
  730 IF (n-nmax) 740,750,750
  740 phib(n)=atan((rb(n+1)-rb(n-1))/(xb(n+1)-xb(n-1)))
      GO TO 760
  750 phib(n)=2.0*atan((rb(n)-rb(n-1))/(xb(n)-xb(n-1)))
      phib(n)=phib(n)-phib(n-1)
  760 isbq=sb(n)
      IF (isbq) 810,770,810
  770 fkw=(sin(phib(n))-sin(phib(n-1)))/(xb(n)-xb(n-1))
      IF (abs(fkw/(xb(n)-xb(n-1)))-1.e-06) 790,790,780
  780 delswa=(phib(n)-phib(n-1))/fkw
      GO TO 800
  790 delswa=sqrt((xb(n)-xb(n-1))**2+(rb(n)-rb(n-1))**2)
  800 swa=swa+delswa
      sb(n)=swa
  810 IF (n-nmax) 720,820,820
  820 IF (ikind-3) 860,830,830
!     CONDITIONS BEHIND OUTER SHOCK
  830 CONTINUE
      IF (iki.eq.0) READ (5,70) xsh,rsh,psi,psh,tsh,ush
      IF (iki.eq.0) READ (5,360) (csh(i),i=1,nds)
      IF (iki.ne.0) READ (9,70) xsh,rsh,psi,psh,tsh,ush
      WRITE (6,1250) xsh,rsh,psi,psh,tsh,ush
      IF (iki.ne.0) READ (9,360) (csh(i),i=1,nds)
  840 FORMAT (8e10.3)
      pbs=psh
!     EXTERNAL FLOW PROPERTIES AT FIRST SHOCK POINT
      READ (5,70) xstrem,rstrem,phstrm,pstrem,tstrem,ustrem
      WRITE (6,1260) xstrem,rstrem,phstrm,pstrem,tstrem,ustrem
      IF (irog.eq.1.and.ipart.eq.0) GO TO 850
      READ (5,840) (cstrem(i),i=1,nds)
      GO TO 860
  850 ndd=nds-1
      READ (5,840) (cstrem(i),i=1,ndd)
      cstrem(nds)=0.
  860 IF (itype-3) 880,870,870
!     CONDITIONS BEHIND INNER SHOCK
  870 CONTINUE
      IF (iki.eq.0) READ (5,70) xshw,rshw,psiw,pshw,tshw,ushw
      IF (iki.eq.0) READ (5,360) (cshw(i),i=1,nds)
      IF (iki.ne.0) READ (9,70) xshw,rshw,psiw,pshw,tshw,ushw
      WRITE (6,1270) xshw,rshw,psiw,pshw,tshw,ushw
      IF (iki.ne.0) READ (9,360) (cshw(i),i=1,nds)
      extra(11)=pshw
  880 IF (itype-3) 1130,890,890
  890 CONTINUE
      IF (iki.eq.0) READ (5,70) xoo,roo,phoo,poo,too,uoo
      IF (iki.eq.0) READ (5,360) (coo(i),i=1,nds)
      IF (iki.ne.0) READ (9,70) xoo,roo,phoo,poo,too,uoo
      WRITE (6,1280) xoo,roo,phoo,poo,too,uoo
      IF (iki.ne.0) READ (9,360) (coo(i),i=1,nds)
      pbs=psh
      GO TO 1130
  900 WRITE (6,1210)
      fmax=1.0
      kay=0.0
      y(1)=0.0
      IF (isb.eq.2) GO TO 950
      READ (4,920) (c(i,2),i=1,nds)
      READ (4,910) x(1),r(2),pstrem,tstrem,ustrem,phstrm
      READ (4,920) (cstrem(i),i=1,nds)
      READ (4,930) thtp,phi(2),p(2),t(2),u(2),psi,psh,tsh,ush
  910 FORMAT (6e12.5)
  920 FORMAT (1e10.3)
  930 FORMAT (9e12.5)
      phi(1)=0.0
      r(1)=0.0
      x(2)=x(1)
      xstrem=x(1)
      rstrem=r(2)
      xsh=x(1)
      rsh=r(2)
      phio=phi(2)
      DO 940 i=1,nds
  940   csh(i)=cstrem(i)
      GO TO 960
  950 READ (5,70) x(1),r(1),phi(1)
      READ (5,70) x(2),r(2),phi(2),p(2),t(2),u(2)
      READ (5,360) (c(i,2),i=1,nds)
      READ (5,360) (cstrem(i),i=1,nds)
      READ (5,70) xstrem,rstrem,phstrm,pstrem,tstrem,ustrem
      READ (5,70) xsh,rsh,psi,psh,tsh,ush
      READ (5,360) (csh(i),i=1,nds)
  960 r(kmn)=rstrem
      rmd=rstrem
      massmd=0.0
      dr21=r(kmn)/(kmn-1)
      DO 990 k=2,kmn
        r(k)=r(k-1)+dr21
        x(k)=xstrem
        b=sqrt((x(k)-x(k-1))**2+(r(k)-r(k-1))**2)
        dely(k)=b
        aa(k)=pi*(r(k)+r(k-1))**delta*dely(k)
        y(k)=y(k-1)+dely(k)
        p(k)=p(2)
        t(k)=t(2)
        u(k)=u(2)
        WRITE (6,1220) k,x(k),r(k),phi(k),p(k),t(k),u(k)
        zmw=0.0
        DO 970 i=1,nds
          c(i,k)=c(i,2)
  970     zmw=zmw+c(i,k)/mw(i)
        zmw=1.0/zmw
        zmwk(k)=zmw
        rho(k)=zmw*p(k)/(rv*t(k))
        mdot(k)=rho(k)*u(k)*aa(k)
        massmd=massmd+mdot(k)
        h(k)=0.0
        DO 980 i=1,nds
          DO 980 j=1,nn
  980     h(k)=h(k)+a(i,j)*c(i,k)*t(k)**(j-2)
  990 END DO
      WRITE (6,1250) xsh,rsh,psi,psh,tsh,ush
      WRITE (6,1260) xstrem,rstrem,phstrm,pstrem,tstrem,ustrem
      pbs=psh
      kxx=kmn+1
      kxxx=kmn+2
      DO 1020 k=kxx,kxxx
        zmw=0.0
        DO 1000 i=1,nds
          c(i,k)=csh(i)
 1000     zmw=zmw+c(i,k)/mw(i)
        zmw=1.0/zmw
        zmwk(k)=zmw
        p(k)=psh
        t(k)=tsh
        u(k)=ush
        rho(k)=zmw*p(k)/(rv*t(k))
        h(k)=0.0
        DO 1010 i=1,nds
          DO 1010 j=1,nn
 1010     h(k)=h(k)+a(i,j)*c(i,k)*t(k)**(j-2)
 1020 END DO
!
!     CALCULATING SIGNIFICENT DIVIDING STREAMLINE DIMENSIONS-L1 AND L2
!
      l2=130.*rmd/(thtp*57.3)**1.55
      l2=par(4)*l2
      l1=0.264*l2/(rmd/rnoz)**0.37
      l1=par(5)*l1
      ddis=0.02*l1
      rdslp=rmd+tan(phio)*ddis*(1.0-(ddis/2./l1))
      phdslp=atan(tan(phio)*(1.0-ddis/l1))
      xinit=xsh
      xdis=xsh+ddis
      x(1)=xdis
      DO 1070 k=2,kmn
        r(k)=rdslp
        dely(k)=r(k)-r(k-1)
        x(k)=xdis
        y(k)=y(k-1)+dely(k)
        DO 1030 i=1,nds
 1030     xs(i,k)=c(i,k)*zmw/mw(i)
        DO 1040 i=1,nds
          cps(i)=0.0
          DO 1040 j=1,nn
 1040     cps(i)=cps(i)+float(j-2)*a(i,j)*t(k)**(j-3)
        cp=0.0
        DO 1050 i=1,nds
 1050     cp=cp+cps(i)*c(i,k)
        re(k)=rho(k)*u(k)*dely(k)/(mu*fmax)
        ma(k)=u(k)*sqrt((cp*zmw-rcon)/(cp*t(k)*ru))
        gama=cp/(cp-rcon/zmw)
        pq=4.0*gama/(9.*gama-5.0)
        phi(k)=phdslp
        IF (k.gt.2) GO TO 1070
!
!     CALCULATING TROAT RADIUS FROM MASS FLOW THRU MACH DISC-
!    1   ASSUMUNG CONSTANT GAMMA
!
        factor=(1.+(gama-1.)*0.5*ma(2)**2)
        pstg=p(2)*factor**(gama/(gama-1.0))
        tstg=t(2)*factor
        rug1=(2.0/(gama+1.0))**((gama+1.0)/(gama-1.0))
        rug2=sqrt(gama*zmwk(k)*rug1*betap/rv)
        astar=massmd*sqrt(tstg)/rug2/pstg
        rstar=sqrt(astar/pi)
        rstar=par(9)*rstar
        WRITE (6,1060) massmd,rstar,l2,l1,rnoz,rmd
 1060   FORMAT (/5x,'MASS FLOW THRU MACH DISC=',1pe10.3,' GM/SEC'/5x,'TH&
     &ROAT RADIUS =',1pe10.3,' CM',/5x,'DISTANCE FROM MACH DISC TO THROA&
     &T  =',1pe10.3,' CM'/5x,'DISTANCE TO DSL MAXIMUM =',1pe10.3,' CM'/ &
     &   5x,'NOZZLE EXIT RADIUS =',1pe10.3,' CM'/5x,'MACH DISC RADIUS=  &
     &  ',1pe10.3,' CM'/)
        zaaa=xinit/rnoz
        zbbb=rmd/rnoz
        WRITE (11) zaaa,zbbb
 1070 END DO
      aaa=tan(psi)
      bbb=rmd-xsh*tan(psi)
      ccc=-tan(pi/2.-phdslp)
      ddd=rdslp-xdis*ccc
      xsh=-(bbb-ddd)/(aaa-ccc)
      rsh=aaa*xsh+bbb
      nsuper=kxxx-kxx+1
      delr=(rsh-rdslp)/nsuper
      delx=(xdis-xsh)/nsuper
      DO 1110 k=kxx,kxxx
        r(k)=r(k-1)+delr
        x(k)=x(k-1)-delx
        dely(k)=sqrt((x(k)-x(k-1))**2+(r(k)-r(k-1))**2)
        y(k)=y(k-1)+dely(k)
        phi(k)=phdslp
        aa(k)=pi*(r(k)+r(k-1))**delta*dely(k)
        DO 1080 i=1,nds
 1080     xs(i,k)=c(i,k)*zmw/mw(i)
        mdot(k)=rho(k)*u(k)*aa(k)
        DO 1090 i=1,nds
          cps(i)=0.0
          DO 1090 j=1,nn
 1090     cps(i)=cps(i)+float(j-2)*a(i,j)*t(k)**(j-3)
        cp=0.0
        DO 1100 i=1,nds
 1100     cp=cp+cps(i)*c(i,k)
        re(k)=rho(k)*u(k)*dely(k)/(mu*fmax)
        ma(k)=u(k)*sqrt((cp*zmw-rcon)/(cp*t(k)*ru))
        gama=cp/(cp-rcon/zmw)
        pq=4.0*gama/(9.*gama-5.0)
 1110 END DO
      xw(1)=xdis
      xw(2)=1.e8
      rw(1)=0.0
      rw(2)=0.0
      rb(1)=r(kmax)
      rb(2)=r(kmax)
      xb(1)=x(kmax)
      xb(2)=x(kmax)+0.01
      y0=y(kmax)
      sw(1)=0.0
      sw(2)=1.e8
      phiw(1)=0.0
      phiw(2)=0.0
      sb(1)=0.0
      sb(2)=1.e-2
      phib(1)=phi(kmax)
!     CARD OMITTED
      k=kxxx
      xsh=x(k)
      rsh=r(k)
      xstrem=x(k)
      rstrem=r(k)
      WRITE (6,1120) xsh,rsh,xstrem,rstrem
 1120 FORMAT (/5x,'REVISED INITIAL SHOCK AND STREAM COORDINATES ARE-    &
     &  XSH=',1pe10.3,' RSH=',1pe10.3,' XSTREM=',1pe10.3,' RSTREM=',    &
     &1pe10.3,/)
 1130 CALL putout (1)
 1140 FORMAT (' ',/5x,'KMAX',4x,'NN',6x,'KP',6x,'IP',6x,'ITYPE',3x,'IKIN&
     &D',3x,'MMAX',4x,'NMAX',4x,'NDS',5x,'IPD',5x,'NITER',3x,'IDIFF',4x,&
     &'LPLANE',/5x,13(i5,3x),/5x,'IPTUC',3x,'IPCH',4x,'ISHOCK',2x,'IBUGS&
     &H',2x,'ITD',5x,'ICSH',4x,'ITAPES',2x,'IKINE',3x,'IPART',3x,'IKI', &
     &5x,'ITURB',3x,'IPTO',5x,'ISB',5x,'KMN',5x,'KMNF',/5x,15(i5,3x))
 1150 FORMAT (' ',/5x,'ALPHAH',4x,'EPSLON',4x,'TOL',7x,'DELTA',5x,'ATOL'&
     &,6x,'FSTEP',5x,'GRAD',6x,'FRAC',6x,'FRACTN',/3x,9(f8.4,2x))
 1160 FORMAT (','/5x,'I=',8x,'ALFTD=',9x,'DEXP=')
 1170 FORMAT (5x,i5,5x,1pe10.3,5x,e10.3)
 1180 FORMAT (' ',/5x,'IDENT',10x,'MUO',12x,'TO',12x,'OMEGA',10x,'PR',  &
     &13x,'SC',13x,'MW')
 1190 FORMAT (5x,a4,9x,6(1pe10.3,5x))
 1200 FORMAT (' ',/5x,'VISCOUS OPTION VALUES%',5x,'TLE=',1pe10.3,5x,'TPR&
     &=',e10.3,5x,'EDDYK=',e10.3,5x,'PERGF=',e10.3)
 1210 FORMAT (' ',/2x,'K',2x,'X',14x,'R',14x,'PHI',12x,'P',14x,'T',14x,'&
     &U')
 1220 FORMAT (' ',1x,i2,1x,6(1pe10.3,5x))
 1230 FORMAT (' ',/2x,'I',2x,'XW',13x,'RW',13x,'PHIW',11x,'PW',13x,'SW')
 1240 FORMAT (' ',/2x,'N',2x,'XB',13x,'RB',13x,'PHIB',11x,'PB',13x,'SB')
 1250 FORMAT (' ',/5x,'XSH=',1pe10.3,6x,'RSH=',e10.3,6x,'PSI=',e10.3,6x,&
     &'PSH=',e10.3,6x,'TSH=',e10.3,6x,'USH=',e10.3)
 1260 FORMAT (' ',/5x,'XSTRM=,'1pe10.3,3x,'RSTREM=',e10.3,3x,'PHSTRM=', &
     &e10.3,3x,'PSTREM=',e10.3,3x,'TSTREM=',e10.3,3x,'USTREM=',e10.3)
 1270 FORMAT (' ',/5x,'XSHW=',1pe10.3,5x,'RSHW=',e10.3,5x,'PSIW=',e10.3,&
     &5x,'PSHW=',e10.3,5x,'TSHW=',e10.3,5x,'USHW=',e10.3)
 1280 FORMAT (' ',/5x,'XOO=',1pe10.3,6x,'ROO=',e10.3,6x,'PHOO=',e10.3,  &
     &5x,'POO=',e10.3,6x,'TOO=',e10.3,6x,'UOO=',e10.3)
      RETURN
      END  Subroutine Putin

      SUBROUTINE putout (iout)
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      COMMON /chem1/ zid(5),irr(30),irt(30),rc(30,3),irrr(30,5),av,     &
     &cm(21,21),fi(21),wp(21),wm(21),wdot(21,30),csi(21),qx(21)
      COMMON /shoput/ xshw,rshw,psiw,pshw,tshw,ushw,cshw(20)
      COMMON /cff/ ipto
      COMMON /intl/ xlx,iki,rex
      COMMON /sub/ xdis,l1,l2,phio,rstar,rmd,massmd,rnoz,thtp,thru
      DIMENSION sumi(24), slmd(20)
      REAL lmd(20,30)
!      DATA iq1/'CO  '/
!      DATA iq2/'CO2 '/
!      DATA iq3/'H2O '/
!      DATA iq4/'HCL '/
!      DATA iq5/'N2  '/
!      DATA iq6/'H2  '/
!      DATA iq7/'E   '/

  CHARACTER(LEN=4),PARAMETER:: IQ1='CO  ',IQ2='CO2 ', IQ3='H2O '
  CHARACTER(LEN=4),PARAMETER:: IQ4='HCL ',IQ5='N2  ', IQ6='H2  '
  CHARACTER(LEN=4),PARAMETER:: IQ7='E   '


      GO TO (20,430,10),iout
   10 l=k
      GO TO 60
   20 l=1
      sumdot=0.0
      DO 30 i=1,nds
   30   sumi(i)=0.0
      WRITE (6,40) ll,x2(l),r2(l),phi2(l)
   40 FORMAT ('0',i5,8x,'XW=',e11.4,3x,'RW=',e11.4,3x,'PHIW=',e11.4)
      IF ((isb.ne.-1).or.(x2(1).lt.xlmd)) GO TO 50
      zzz=x2(1)/rnoz
      WRITE (12) zzz,k
   50 CONTINUE
      l=2
   60 IF (ll) 70,80,80
   70 WRITE (6,90) l,x(l),r(l),phi(l),t(l),p(l),u(l)
      xl=x(l)/rnoz
      rl=r(l)/rnoz
      GO TO 100
   80 WRITE (6,90) l,x2(l),r2(l),phi2(l),t2(l),p2(l),u2(l)
      xl=x2(l)/rnoz
      rl=r2(l)/rnoz
   90 FORMAT ('0',5x,i5,3x,'X=',1p1e11.4,4x,'R=',1p1e11.4,4x,'PHI=',    &
     &1p1e11.4,2x,'T=',1p1e11.4,4x,'P=',1p1e11.4,4x,'U=',1p1e11.4)
  100 htk=h(l)+u(l)**2/(2.0*hj)
      ptk=p(l)+rho(l)*u(l)**2/(2.0*betap)
      sumdot=sumdot+mdot(l)
      IF (ll) 110,130,130
  110 WRITE (6,120) ma(l),dely(l),h(l),htk,taw(l),q(l),ptk,rho(l),sx(l),&
     &sumdot,xl,rl
  120 FORMAT (14x,'MA=',e11.4,3x,'DELY=',e11.4,1x,'H=',e11.4,4x,'HT=',  &
     &e11.4,3x,'TAW=',e11.4,2x,'Q=',e11.4/14x,'PT=',e11.4,3x,'RHO=',    &
     &e11.4,2x,'SX=',e11.4,3x,'SUMDOT=',e11.4,2x,'X/REX=',f6.2,2x,'R/REX&
     &=',f6.2)
      GO TO 140
  130 WRITE (6,120) ma(l),dely(l),h(l),htk,taw(l),q(l),ptk,rho2(l),sx(l)&
     &,sumdot,xl,rl
  140 WRITE (6,150)
  150 FORMAT ('0')
      IF (nds-1) 270,270,160
  160 DO 250 i=1,nds
        IF (iextra(1)) 210,170,210
  170   IF (ll) 180,190,190
  180   lmd(i,l)=mdot(l)*c(i,l)/u(l)
        WRITE (6,200) ident(i),c(i,l),xs(i,l),wdot(i,l),lmd(i,l)
        GO TO 250
  190   lmd(i,l)=mdot(l)*c(i,l)/u2(l)
        WRITE (6,200) ident(i),c2(i,l),xs(i,l),wdot(i,l),lmd(i,l)
  200   FORMAT (14x,'SPECIE',2x,a4,7x,'C=',e11.4,4x,'X=',e11.4,4x,'WDOT=&
     &',e11.4,4x,'LMD=',e11.4)
        GO TO 250
  210   sumi(i)=sumi(i)+mdot(l)*c(i,l)
        IF (ll) 220,240,240
  220   lmd(i,l)=mdot(l)*c(i,l)/u(l)
        WRITE (6,230) ident(i),c(i,l),xs(i,l),wdot(i,l),sumi(i),zj(i,l),&
     &   lmd(i,l)
  230   FORMAT (14x,'SPECIE',1x,a4,2x,'C=',1p1e11.4,3x,'X=',1p1e11.4,3x,&
     &   'WDOT=',e11.4,3x,'F1=',1p1e11.4,3x,'ZJ=',1p1e11.4,3x,'LMD=',   &
     &   1p1e11.4)
        GO TO 250
  240   lmd(i,l)=mdot(l)*c(i,l)/u2(l)
        WRITE (6,230) ident(i),c2(i,l),xs(i,l),wdot(i,l),sumi(i),zj(i,l)&
     &   ,lmd(i,l)
  250 END DO
      IF ((isb.ne.-1).or.(x2(1).lt.xlmd)) GO TO 260
      zzz=r2(l)/rnoz
      WRITE (12) zzz,lmd(3,l),lmd(13,l)
  260 CONTINUE
  270 IF (ipto.eq.0) GO TO 340
      xl=x2(l)/rex
      rl=r2(l)/rex
      iz2=25
      iz3=25
      iz4=25
      iz5=25
      iz6=25
      iz7=25
      DO 280 j=1,nds
        IF (ident(j).eq.iq1) iz1=j
        IF (ident(j).eq.iq2) iz2=j
        IF (ident(j).eq.iq3) iz3=j
        IF (ident(j).eq.iq4) iz4=j
        IF (ident(j).eq.iq5) iz5=j
        IF (ident(j).eq.iq6) iz6=j
        IF (ident(j).eq.iq7) iz7=j
  280 END DO
      ene=0.733e22*p2(l)/t2(l)*c2(iz7,l)
      vee=6.21e5*sqrt(t(l))
      smq=c2(iz1,l)*(2.08e-23*vee+2.46e-16)+c2(iz2,l)*4.7e-8/vee+c2(iz3,&
     &l)*5.9/(vee**2)+c2(iz4,l)*1.85/(vee*vee)+c2(iz5,l)*3.29e-23*vee+  &
     &c2(iz6,l)*(1.45e-23*vee+8.9e-16)
      enue=4.57e27*p2(l)/sqrt(t2(l))*smq
      sgm=7.14549e-4*ene/enue
      IF (ipto.gt.4) GO TO 320
      WRITE (6,290) xl,rl,ene,enue,sgm
  290 FORMAT (2x,'X/REX=',1pe11.4,5x,'R/REX=',e11.4,5x,'NE=',e11.4,5x,'N&
     &UE=',e11.4,5x,'SIGMA=',e11.4)
      IF (ipto.eq.1.or.ipto.eq.3) WRITE (7,300) xl,rl,ene,enue,sgm
  300 FORMAT (4e19.12/3e19.12)
      IF (ipto.eq.2.or.ipto.eq.4) WRITE (8,310) xl,rl,ene,enue,sgm
  310 FORMAT (7e19.12)
  320 WRITE (6,330) xl,rl,t2(l),p2(l),xs(iz2,l),xs(iz3,l),xs(iz1,l)
  330 FORMAT (5x,'X/REX=',1pe11.4,5x,'R/REX=',e11.4,9x,'T=',e11.4,9x,'P=&
     &',e11.4/10x,'(CO2)=',e11.4,6x,'(H2O)=',e11.4,7x,'(CO)=',e11.4)
      IF (ipto.eq.5.or.ipto.eq.3) WRITE (7,300) xl,rl,t2(l),p2(l),      &
     &xs(iz2,l),xs(iz3,l),xs(iz1,l)
      IF (ipto.eq.6.or.ipto.eq.3) WRITE (8,310) xl,rl,t2(l),p2(l),      &
     &xs(iz2,l),xs(iz3,l),xs(iz1,l)
  340 IF (l-k) 350,360,360
  350 l=l+1
      GO TO 60
  360 IF (isb.ne.1) GO TO 480
      IF (nds-1) 480,480,370
  370 DO 380 i=1,nds
  380   slmd(i)=0.0
      DO 390 i=1,nds
        DO 390 l=2,k
  390   slmd(i)=slmd(i)+lmd(i,l)
      WRITE (6,400)
  400 FORMAT ('0',13x,'TOTAL LINEAR MASS DENSITIES')
      WRITE (6,150)
      DO 410 i=1,nds
  410   WRITE (6,420) ident(i),slmd(i)
  420 FORMAT (14x,'SPECIE',1x,a4,2x,'SUMLMD=',1p1e11.4)
      IF (ll.le.10) GO TO 480
      zzz=x2(1)/rnoz
      xxx=r2(k)/rnoz
      WRITE (13) zzz,xxx,slmd(3),slmd(13)
      GO TO 480
  430 IF (nbound) 460,440,460
  440 k=1
      delphi=phi(k)-phish1
      WRITE (6,470) ll,k,x2(k),r2(k),phi2(k),psiw,delphi
      IF (isb.ne.-1) GO TO 450
      aaa=x2(k)/rnoz
      bbb=r2(k)/rnoz
      WRITE (11) aaa,bbb
  450 CONTINUE
      k=2
      GO TO 480
  460 delphi=phi(kmax)-phish1
      WRITE (6,470) ll,kmax,x(kmax),r(kmax),phi(kmax),psi,delphi
  470 FORMAT ('0','SHOCK',2x,'LL=',i4,2x,'K=',i2,2x,'X=',1p1e11.4,2x,'R=&
     &',1p1e11.4,2x,'PHI=',1p1e11.4,2x,'PSI=',1p1e11.4,2x,'DELPHI=',    &
     &1p1e11.4)
  480 RETURN
      END  Subroutine PutOut

      SUBROUTINE shocke
!     THIS SUBROUTINE DOES THE SHOCK CALCULATION
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /shoput/ xshw,rshw,psiw,pshw,tshw,ushw,cshw(20)
      COMMON /shosho/ x3sh,r3sh
      COMMON /turple/ del2,cpshw,gamshw,zmwsh1
      COMMON /mainsh/ kount,kint
      EQUIVALENCE (pwsh,extra(11))
      DIMENSION sumasb(25), csh2(25), csh3(25), sumasw(25)
      DIMENSION g(2), dgd(2,2), xx(2)
      fl=ll
      fplane=lplane
      fmult=fl/fplane
      IF (abs(extra(4)).gt.1.0e-7) GO TO 10
      gstep=fstep
      GO TO 20
   10 gstep=fstep-(fstep-extra(3))*(1.-fmult)/extra(4)
   20 IF (ll) 30,30,250
!     DETERMINING UPSTREAM PROPERTIES AT INITIAL SHOCK POINT
   30 IF (iflag) 250,40,250
   40 zmb=0.0
      rold=rshw
      GO TO (70,50,70,70,50),ishock
   50 IF (nbound) 60,70,60
   60 nk=k
      k=kmax
      r2sh=rsh
      x2sh=xsh
      CALL shopro
      phstrm=phish1
      pstrem=psh1
      tstrem=tsh1
      ustrem=ush1
      k=nk
   70 IF (ikind-3) 150,150,80
   80 cpb=0.0
      DO 90 i=1,nds
        zmb=zmb+csh(i)/mw(i)
   90   cpb=cpb+csh(i)*a(i,3)
      zmb=1.0/zmb
      gbold=zmb*pstrem*ustrem/(rv*tstrem)
      bmasum=0.0
      DO 100 i=1,nds
        sumasb(i)=0.
  100   csh2(i)=csh(i)
      IF (iconst) 120,110,120
  110 fhbold=gbold*(cpb*tstrem+ustrem**2/(2.0*hj))
      GO TO 140
  120 fhbold=gbold*ustrem**2/(2.0*hj)
      DO 130 i=1,nds
        DO 130 j=1,nn
  130   fhbold=fhbold+gbold*a(i,j)*csh1(i)*tstrem**(j-2)
  140 phbold=phstrm
      zwm=0.0
      IF (nbound) 250,150,250
  150 IF (ishock-4) 170,160,170
  160 mk=k
      k=1
      r3sh=rshw
      x3sh=xshw
      CALL shopro
      phoo=phish1
      poo=psh1
      too=tsh1
      uoo=ush1
      k=mk
  170 IF (itype-2) 250,250,180
  180 cpw=0.0
      DO 190 i=1,nds
        zwm=zwm+cshw(i)/mw(i)
  190   cpw=cpw+cshw(i)*a(i,3)
      zwm=1.0/zwm
      gwold=zwm*poo*uoo/(rv*too)
      wmasum=0.0
      whsum=0.0
      DO 200 i=1,nds
        sumasw(i)=0.
  200   csh3(i)=cshw(i)
      IF (iconst) 220,210,220
  210 fhwold=gwold*(cpw*too+uoo**2/(2.0*hj))
      GO TO 240
  220 fhwold=gwold*uoo**2/(2.0*hj)
      DO 230 i=1,nds
        DO 230 j=1,nn
  230   fhwold=fhwold+gwold*a(i,j)*csh1(i)*too**(j-2)
  240 phwold=phoo
!     LOCATING NEW SHOCK POINT
  250 tph2k=tan(phi2(k))
      IF (nbound) 260,270,260
  260 tpsi=tan(psi)
      r2sh=(rsh-xsh*tpsi+r2(k)*tph2k*tpsi+x2(k)*tan(psi))/(1.0+tpsi*    &
     &tph2k)
      x2sh=(r2(k)-r2sh)*tph2k+x2(k)
      GO TO 280
  270 tpsiw=tan(psiw)
      r3sh=(rshw-xshw*tpsiw+r2(k)*tph2k*tpsiw+x2(k)*tpsiw)/(1.0+tpsiw*  &
     &tph2k)
      x3sh=(r2(k)-r3sh)*tph2k+x2(k)
  280 CALL shopro
      cpsh=0.0
      mwsh=0.0
      DO 300 i=1,nds
        DO 290 j=1,nn
  290     cpsh=cpsh+a(i,j)*float(j-2)*csh1(i)*tsh1**(j-3)
  300   mwsh=mwsh+csh1(i)/mw(i)
      mwsh=1./mwsh
      gamma=cpsh/(cpsh-rcon/mwsh)
      mash1=ush1/sqrt(gamma*ru*tsh1/mwsh)
      IF (iconst) 600,310,600
!     SHOCK CALCULATION FOR CONSTANT GAMMA
!     DETERMINATION OF SHOCK ANGLE
  310 delphi=phi2(k)-phish1
      b1=-(mash1**2+2.0)/(mash1**2)-gamma*sin(delphi)**2
      b2=(2.0*mash1**2+1.0)/(mash1**4)
      b2=b2+((gamma+1.0)**2/4.0+(gamma-1.0)/(mash1**2))*sin(delphi)**2
      b3=cos(delphi)**2/(mash1**4)*(-1.0)
      iter=0
      IF (nbound) 360,320,360
  320 IF (delphi) 370,330,330
  330 zeta=(1.0/mash1)**2
      us=ush1*(1.+abs(delphi)/sqrt(mash1**2-1.))
      ts=tsh1-(us**2-ush1**2)/(2.0*hj*cpsh)
      ps=psh1*(ts/tsh1)**(gamma/(gamma-1.))
      IF (nbound) 350,340,350
  340 pshw=ps
      tshw=ts
      ushw=us
      GO TO 540
  350 psh=ps
      tsh=ts
      ush=us
      GO TO 470
  360 IF (delphi) 330,330,380
  370 zeta=(1.0/mash1-sin(delphi))**2
      GO TO 390
  380 zeta=(1.0/mash1+sin(delphi))**2
  390 iter=iter+1
      gg=zeta**3+b1*zeta**2+b2*zeta+b3
      dgdz=3.0*zeta**2+2.0*b1*zeta+b2
      dlzeta=-gg/dgdz
      IF (abs(dlzeta/zeta)-tol) 450,400,400
  400 IF (iter-niter) 410,410,420
  410 zeta=zeta+dlzeta
      GO TO 390
  420 WRITE (6,430)
  430 FORMAT ('1ITERATIONS FOR ZETA FAILED')
      WRITE (6,440) iter,niter,k,ll,iflag,tol,zeta,gg,b1,b2,b3
  440 FORMAT (1x,'ITER=',i5,5x,'NITER=',i5,5x,'K=',i5,5x,'LL=',i5,5x,'IF&
     &LAG=',i5/' TOL=',1pe15.6,5x,'ZETA=',1pe15.6,5x,'GG=',1pe15.6/1x,'B&
     &1=',1pe15.6,5x,'B2=',1pe15.6,5x,'B3=',1pe15.6)
      CALL pdump (a(1,1), n, 5)
!     DETERMINATION OF SHOCKED GAS PROPERTIES
  450 con4=cpsh*tsh1+ush1**2/(2.0*hj)
      IF (nbound) 460,530,460
  460 psh=(2.0*gamma*mash1**2*zeta-gamma+1.0)/(gamma+1.0)*psh1
      ttsh=(psh/psh1)*((gamma-1.0)*mash1**2*zeta+2.0)
      tsh=ttsh/((gamma+1.0)*mash1**2*zeta)*tsh1
      ush=sqrt(2.0*hj*(con4-cpsh*tsh))
  470 mash=mash1*ush/ush1*sqrt(tsh1/tsh)
      DO 480 i=1,nds
        csh(i)=csh1(i)
        IF (isb.eq.1) csh(i)=abs(csh1(i))
  480 END DO
      gb=mwsh*psh1*ush1/(rv*tsh1)
      IF (iconst) 500,490,500
  490 fhb=gb*(cpsh*tsh1+ush1**2/(2.0*hj))
      GO TO 520
  500 fhb=gb*ush1**2/(2.0*hj)
      DO 510 i=1,nds
        DO 510 j=1,nn
  510   fhb=fhb+gb*a(i,j)*csh1(i)*tsh1**(j-2)
  520 angle=0.5*(phish1+phbold)
      GO TO 840
  530 pshw=(2.0*gamma*mash1**2*zeta-gamma+1.0)/(gamma+1.0)*psh1
      ttsh=(pshw/psh1)*((gamma-1.0)*mash1**2*zeta+2.0)
      tshw=ttsh/((gamma+1.0)*mash1**2*zeta)*tsh1
      ushw=sqrt(2.0*hj*(con4-cpsh*tshw))
  540 mash=mash1*ushw/ush1*sqrt(tsh1/tshw)
      DO 550 i=1,nds
  550   cshw(i)=csh1(i)
      gw=mwsh*psh1*ush1/(rv*tsh1)
      IF (iconst) 570,560,570
  560 fhw=gw*(cpsh*tsh1+ush1**2/(2.0*hj))
      GO TO 590
  570 fhw=gw*ush1**2/(2.0*hj)
      DO 580 i=1,nds
        DO 580 j=1,nn
  580   fhw=fhw+gw*a(i,j)*csh1(i)*tsh1**(j-2)
  590 angle=0.5*(phish1+phwold)
      GO TO 840
  600 h1=0.0
      DO 610 j=1,nn
        DO 610 i=1,nds
  610   h1=h1+a(i,j)*csh1(i)*tsh1**(j-2)
      DO 640 i=1,nds
        IF (nbound) 620,630,620
  620   csh(i)=csh1(i)
        GO TO 640
  630   cshw(i)=csh1(i)
  640 END DO
      iter=0
      IF (nbound) 650,660,650
  650 ts=tsh
      ths=psi-phish1
      delphi=phi(kmax)-phish1
      IF (delphi) 330,330,680
  660 ts=tshw
      ths=phish1-psiw
      delphi=phish1-phi(k)
      del2=abs(delphi)
      IF (delphi) 670,670,680
  670 delphi=-delphi
      GO TO 330
  680 e=mash1**2
      f=gamma+1.
      term1=gamma*e*delphi/sqrt(e-1.)
      term2=term1*(f*e**2-4.*(e-1.))*delphi/(4.*(e-1.)**(3./2.))
      IF (term2.gt.0.001*term1) GO TO 690
      term3=term1*(f**2*e**4/32.-(7.+12.*gamma-3.*gamma**2)*e**3+.75*f* &
     &e**2-e+2./3.)*delphi**2/(e-1.)**3
      xi=1.+term1+term2+term3
      zeta=(f*xi+f-2.)/(2.*gamma*e)
      GO TO 450
  690 iter=iter+1
      h2=0.0
      cp=0.0
      DO 700 i=1,nds
        DO 700 j=1,nn
        h2=h2+a(i,j)*csh1(i)*ts**(j-2)
  700   cp=cp+float(j-2)*a(i,j)*ts**(j-3)*csh1(i)
      sit=sin(ths)
      sid=sin(ths-delphi)
      coot=cos(ths)
      cod=cos(ths-delphi)
      tat=tan(ths)
      tad=tan(ths-delphi)
      g(1)=-(h2-h1-ush1**2*(1.-(coot/cod)**2)/(2.*hj))
      g(2)=-(psh1*(tsh1/tat-ts/tad)+psh1*mwsh*ush1**2/(rv*betap)*sit*   &
     &coot*(1.-tad/tat))
      dgd(1,1)=cp
      dgd(1,2)=ush1**2*coot*(coot*tad-sit)/(hj*cod**2)
      dgd(2,1)=-psh1/tad
      dgd(2,2)=psh1*(ts/sid**2-tsh1/sit**2)+psh1*mwsh*ush1**2/(rv*betap)&
     &*((coot**2-sit**2)*(1.-tad/tat)-coot**2/cod**2+tad/tat)
      yyy=0.0
      CALL mate7 (2, 2, 1, dgd, g, yyy, xx, j)
      GO TO (760,710,730),j
  710 WRITE (6,720)
  720 FORMAT ('1DETERMINANT OVERFLOW IN SUBROUTINE SHOCKE')
      GO TO 750
  730 WRITE (6,740)
  740 FORMAT ('1SINGULAR DETERMINANT IN SUBROUTINE SHOCKE')
  750 CALL pdump (a(1,1), n, 1)
      STOP   ! CALL exit
  760 DO 770 i=1,2
  770   xx(i)=dgd(i,1)
      IF (abs(xx(1)/ts)-tol) 800,780,780
  780 ts=ts+xx(1)
      ths=ths+xx(2)
      IF (iter-niter) 680,680,790
  790 WRITE (6,430)
      WRITE (6,440) iter,niter,k,ll,iflag,tol,zeta,gg,b1,b2,b3
      GO TO 750
  800 IF (abs(xx(2)/ths)-sqrt(tol)) 810,780,780
  810 IF (nbound) 830,820,830
  820 tshw=ts
      ushw=sqrt(2.*hj*(h1-h2)+ush1**2)
      pshw=psh1*(1.+mwsh/(rv*betap*tsh1)*ush1**2*sit**2*(1.-tad/tat))
      zeta=sin(ths)**2
      delphi=-delphi
      GO TO 540
  830 tsh=ts
      ush=sqrt(2.*hj*(h1-h2)+ush1**2)
      psh=psh1*(1.+mwsh/(rv*betap*tsh1)*ush1**2*sit**2*(1.-tad/tat))
      zeta=sin(ths)**2
      GO TO 470
  840 IF (nbound) 850,1120,850
  850 pbs=psh
      IF (delphi) 870,860,860
  860 psitry=asin(sqrt(zeta))+phish1
      GO TO 880
  870 psitry=asin(sqrt(zeta))+phi2(k)
  880 IF (iflag) 900,890,900
  890 psi=0.5*(psi+psitry)
      GO TO 1500
!     EXTERNAL SHOCK -- CALCULATION OF PROPERTIES IN STREAMTUBE
!     BETWEEN SHOCK AND OUTERMOST STREAMLINE CARRIED
  900 b=sqrt((r2sh-r(k))**2+(x2sh-x(k))**2)
      bfact=sqrt((xsh-x2sh)**2+(rsh-r2sh)**2)*abs(sin(psi-angle))
      proja=bfact*(rsh+r2sh)**delta*pi
      resh=re(kmax)*psh*ush*t(kmax)**(1.0+omega(1))/(dely(kmax)*p(kmax)*&
     &u(kmax)*tsh**(1.0+omega(1)))
      sumdot=0.0
      DO 910 l=2,kmax
  910   sumdot=sumdot+mdot(l)/(.5*(r2(l)+r2(l-1)))**2
      bmasum=bmasum+proja*0.5*(gb+gbold)
      bhsum=bhsum+proja*0.5*(fhb+fhbold)
      DO 920 i=1,nds
  920   sumasb(i)=sumasb(i)+.5*proja*(gb*csh(i)+gbold*csh2(i))
      IF (bmasum/(.5*(r2sh+r2(kmax-1)))**2-1.2*sumdot/gstep) 1450,930,  &
     &930
!     ADDING NEW STREAMTUBE TO SHOCKED FLOW
  930 kmax=kmax+1
      x(kmax)=x2sh
      r(kmax)=r2sh
      DO 940 i=1,nds
  940   csh(i)=sumasb(i)/bmasum
      IF (delphi) 960,960,950
  950 phi(kmax)=phi(k)
      GO TO 970
  960 phi(kmax)=0.5*(phi(k)+phish1)
  970 p(kmax)=psh
      aa(kmax)=pi*(r(kmax)+r(k))**delta*b
      IF (kount) 980,1080,980
  980 z=(bmasum*rv/(aa(kmax)*mwsh*p(kmax)))**2/(2.0*hj)
      IF (iconst) 1000,990,1000
  990 t(kmax)=(sqrt(cpsh**2+4.0*z*bhsum/bmasum)-cpsh)/(2.0*z)
      GO TO 1070
 1000 iter=0
      ts=tsh
 1010 iter=iter+1
      gg=z*ts**2-bhsum/bmasum
      dgg=2.0*z*ts
      DO 1020 i=1,nds
        DO 1020 j=1,nn
        gg=gg+a(i,j)*csh(i)*ts**(j-2)
 1020   dgg=dgg+float(j-2)*a(i,j)*csh(i)*ts**(j-3)
      delt=-gg/dgg
      ts=ts+delt
      IF (abs(delt/ts)-tol) 1060,1030,1030
 1030 IF (iter-niter) 1010,1040,1040
 1040 WRITE (6,1050)
 1050 FORMAT ('1TEMPERATURE ITERATIONS FAILED')
      GO TO 750
 1060 t(kmax)=ts
 1070 u(kmax)=bmasum*rv*t(kmax)/(aa(kmax)*mwsh*p(kmax))
      GO TO 1090
 1080 t(kmax)=tsh
      u(kmax)=ush
      kount=1
 1090 bmasum=0.0
      bhsum=0.0
      DO 1100 i=1,nds
        c(i,kmax)=csh(i)
        sumasb(i)=0.
 1100   xs(i,kmax)=csh(i)*mwsh/mw(i)
      rho(kmax)=mwsh*p(kmax)/(rv*t(kmax))
      dely(kmax)=b
      y(kmax)=y(k)+b
      mdot(kmax)=rho(kmax)*u(kmax)*aa(kmax)
      re(kmax)=resh*b
      ma(kmax)=u(kmax)/(sqrt(cpsh*t(kmax)*ru/(cpsh*mwsh-rcon)))
      h(kmax)=0.0
      DO 1110 i=1,nds
        DO 1110 j=1,nn
 1110   h(kmax)=h(kmax)+a(i,j)*csh(i)*t(kmax)**(j-2)
      j=kmax
      x2(j)=x(j)
      r2(j)=r(j)
      t2(j)=t(j)
      p2(j)=p(j)
      u2(j)=u(j)
      phi2(j)=phi(j)
      GO TO 1440
 1120 pwsh=pshw
      IF (delphi) 1130,1130,1140
 1130 psiwtr=phish1-asin(sqrt(zeta))
      GO TO 1150
 1140 psiwtr=phi2(k)-asin(sqrt(zeta))
 1150 IF (iflag) 1180,1160,1180
 1160 psiw=0.5*(psiw+psiwtr)
      DO 1170 i=1,nds
 1170   zj(i,kmax)=zj(i,kmax-1)
      GO TO 1500
!     INTERNAL SHOCK -- CALCULATION OF PROPERTIES IN STREAMTUBE
!     BETWEEN SHOCK AND INNERMOST STREAMLINE CARRIED
 1180 b=sqrt((r3sh-r2(k))**2+(x3sh-x2(k))**2)
      DO 1190 i=1,nds
 1190   zj(i,kmax)=zj(i,kmax-1)
      bfact=sqrt((xshw-x3sh)**2+(rshw-r3sh)**2)*abs(sin(psiw-angle))
      projaw=bfact*(rshw+r3sh)**delta*pi
      resh=re(2)*pshw*ushw*t(2)**(1.0+omega(1))/(dely(2)*p(2)*u(2)*tshw*&
     &*(1.0+omega(1)))
      sumdot=0.0
      DO 1200 l=2,kmax
 1200   sumdot=sumdot+mdot(l)/(.5*(r2(l)+r2(l-1)))**2
      wmasum=wmasum+projaw*0.5*(gw+gwold)
      whsum=whsum+projaw*0.5*(fhw+fhwold)
      DO 1210 i=1,nds
 1210   sumasw(i)=sumasw(i)+.5*projaw*(gw*cshw(i)+gwold*csh3(i))
      j=1
      IF (ikind.ne.3) GO TO 1220
      IF ((r3sh-rold)-(1.25*(rold-rolder))) 1480,1230,1230
 1220 IF (wmasum/(.5*(r3sh+r2(j)))**2-1.2*sumdot/gstep) 1480,1230,1230
 1230 kmax=kmax+1
!     REINDEXING STREAMLINES IN SHOCKED GAS
      j=kmax
 1240 aa(j)=aa(j-1)
      DO 1250 i=1,nds
        c(i,j)=c(i,j-1)
        c2(i,j)=c2(i,j-1)
 1250   xs(i,j)=xs(i,j-1)
      dely(j)=dely(j-1)
      h(j)=h(j-1)
      ma(j)=ma(j-1)
      mdot(j)=mdot(j-1)
      p(j)=p(j-1)
      phi(j)=phi(j-1)
      r(j)=r(j-1)
      re(j)=re(j-1)
      rho(j)=rho(j-1)
      sx(j)=sx(j-1)
      t(j)=t(j-1)
      u(j)=u(j-1)
      x(j)=x(j-1)
      y(j)=y(j-1)+b
      p2(j)=p2(j-1)
      phi2(j)=phi2(j-1)
      r2(j)=r2(j-1)
      rho2(j)=rho2(j-1)
      t2(j)=t2(j-1)
      u2(j)=u2(j-1)
      x2(j)=x2(j-1)
      y2(j)=y2(j-1)+b
      delss(j)=delss(j-1)
      iextra(7)=iextra(7)+1
      j=j-1
      IF (j-2) 1260,1260,1240
!     ADDING NEW STREAMTUBE TO SHOCKED FLOW
 1260 j=2
      DO 1270 i=1,nds
 1270   cshw(i)=sumasw(i)/wmasum
      p(j)=pshw
      p2(j)=pshw
      aa(j)=pi*(r3sh+r(j-1))**delta*b
      IF (kint) 1280,1360,1280
 1280 z=(wmasum*rv/(aa(j)*mwsh*p(j)))**2/(2.0*hj)
      IF (iconst) 1300,1290,1300
 1290 t(j)=(sqrt(cpsh**2+4.0*z*whsum/wmasum)-cpsh)/(2.0*z)
      t2(j)=t(j)
      GO TO 1350
 1300 iter=0
      ts=tshw
 1310 iter=iter+1
      gg=z*ts**2-whsum/wmasum
      dgg=2.0*z*ts
      DO 1320 i=1,nds
        DO 1320 l=1,nn
        gg=gg+a(i,l)*cshw(i)*ts**(l-2)
 1320   dgg=dgg+float(l-2)*a(i,l)*cshw(i)*ts**(l-3)
      delt=-gg/dgg
      ts=ts+delt
      IF (abs(delt/ts)-tol) 1340,1330,1330
 1330 IF (iter-niter) 1310,1040,1040
 1340 t(j)=ts
      t2(j)=ts
 1350 u(j)=wmasum*rv*t(j)/(aa(j)*mwsh*p(j))
      u2(j)=u(j)
      GO TO 1370
 1360 t(j)=tshw
      t2(j)=tshw
      u(j)=ushw
      u2(j)=ushw
      kint=1
 1370 wmasum=0.0
      whsum=0.0
      DO 1380 i=1,nds
        c(i,j)=cshw(i)
        c2(i,j)=c(i,j)
        sumasw(i)=0.
 1380   xs(i,j)=c(i,j)*mwsh/mw(i)
      rho(j)=mwsh*p(j)/(rv*t(j))
      rho2(j)=rho(j)
      x(j)=x(j-1)
      x2(j)=x2(j-1)
      r(j)=r(j-1)
      r2(j)=r2(j-1)
      phi2(j)=phi2(j-1)
      phi(j)=phi(j-1)
      dely(j)=b
      y(j)=b
      y2(j)=b
      mdot(j)=rho(j)*u(j)*aa(j)
      re(j)=resh*b
      ma(j)=u(j)/(sqrt(cpsh*t(j)*ru/(cpsh*mwsh-rcon)))
      h(j)=0.0
      DO 1390 i=1,nds
        DO 1390 jj=1,nn
 1390   h(j)=h(j)+a(i,jj)*c(i,j)*t(j)**(jj-2)
      j=1
      x2(j)=x3sh
      r2(j)=r3sh
      IF (delphi) 1400,1410,1410
 1400 phi2(j)=phi2(1)
      GO TO 1420
 1410 phi2(j)=0.5*(phi2(1)+phish1)
 1420 y2(j)=0.0
      DO 1430 j=2,kmax
 1430   y2(j)=y2(j-1)+dely(j)
      IF (ikind.ne.3) GO TO 1440
      rolder=rold
      rold=r3sh
      IF (ikind.ne.3.or.nbound.ne.0) GO TO 1450
      k=2
      CALL putout (3)
 1440 CALL putout (2)
!        PREPARING FOR NEXT STEP DOWNSTREAM
 1450 IF (nbound) 1460,1480,1460
 1460 xsh=x2sh
      rsh=r2sh
      psi=psitry
      gbold=gb
      fhbold=fhb
      phbold=phish1
      DO 1470 i=1,nds
 1470   csh2(i)=csh1(i)
      GO TO 1500
 1480 xshw=x3sh
      rshw=r3sh
      psiw=psiwtr
      gwold=gw
      fhwold=fhw
      phwold=phish1
      DO 1490 i=1,nds
 1490   csh3(i)=csh1(i)
      k=2
 1500 RETURN
      END  Subroutine Shocke

      SUBROUTINE shopro
!     THIS SUBROUTINE DETERMINES THE CONDITIONS UPSTREAM OF THE SHOCK
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /intrp1/ nstrm(2),tsfc(2),psfc(2),usfc(2),phisfc(2),csfc(2,&
     &25),dptsfc(2)
      COMMON /turple/ del2,cpshw,gamshw,zmwsh1
      COMMON /jnd/ jndl(2,120),jndr(2,120)
      COMMON /inshop/ poo,too,uoo,phoo,coo(20)
      COMMON /inpro/ icsh,iin
      COMMON /orshok/ xtry,rtry
      COMMON /shosho/ x3sh,r3sh
      COMMON /utpp/ psh11,tsh11,ush11,flang,csh11(20)
      INTEGER r1s2,r2s2,r1,r3,r2m1,r1p1
      INTEGER:: errCode
      DIMENSION locf(2,120), ntot(2)
      DIMENSION arho(6,6), athet(6,6), ap(6,6)
                        !!!, csh11(20)
      DIMENSION cdum(25)
      EQUIVALENCE (itapes,iextra(2))
      IF (ll) 20,10,20
   10 resh=re(kmax)
!     TESTS FOR UNIFORM OR NON-UNIFORM FLOW
   20 GO TO (30,110,50,80,100),ishock
   30 tsh1=tstrem
      psh1=pstrem
      ush1=ustrem
      DO 40 i=1,nds
   40   csh1(i)=cstrem(i)
      phish1=phstrm
      GO TO 760
   50 IF (k-1) 30,60,30
   60 tsh1=too
      psh1=poo
      ush1=uoo
      DO 70 i=1,nds
   70   csh1(i)=coo(i)
      phish1=phoo
      GO TO 760
   80 IF (k-1) 30,90,30
   90 xtry=x3sh
      rtry=r3sh
      GO TO 120
  100 IF (k-1) 110,60,110
  110 xtry=x2sh
      rtry=r2sh
  120 IF (icsh.ne.0) GO TO 590
!     NON-UNIFORM EXTERNAL FLOW WITH TABULAR DATA
  130 FORMAT (///' BEGIN LOGIC TO FIND PROPERTIES FOR POINT PRECEDING SH&
     &OCK'/' POINT IS-- X='e15.6,' R='e15.6)
      IF (iflag.ne.0) GO TO 500
      IF (card1.gt.1.1) GO TO 360
      card1=2.
      mxstrm=60
      mxtot=120
      ictsfc=0
      i9=1
      GO TO 150
  140 i9=2
  150 ictsfc=ictsfc+1
      j=0
      IF (ictsfc.gt.1) j=1
  160 FORMAT (/' ORTHOGONAL SURFACE'i4/' L',9x,'X',14x,'R',13x,'PHI',13x,'P',&
     &14x,'T',14x,'U',8x,'JNDL  JNDR')
      DO 320 l=1,mxtot
        j=j+1
        IF ((itapes.ne.0).or.(isb.eq.1)) GO TO 180

        READ (5,220,IOSTAT=errCode) fx(i9,j),fr(i9,j),fphi(i9,j),fp(i9,j),ft(i9,j),    &
     &   fu(i9,j),jndl(i9,l),jndr(i9,l)
        IF (errCode < 0) GO TO 830
!!!        IF (eof(5)) 830,170
  170   READ (5,230) (fc(i9,j,i),i=1,nds)
  180   IF (isb.eq.1) GO TO 200


        READ (10,220,IOSTAT=errCode) fx(i9,j),fr(i9,j),fphi(i9,j),fp(i9,j),ft(i9,j),   &
     &   fu(i9,j),jndl(i9,l),jndr(i9,l)
        IF (errCode < 0) GO TO 830
!!!        IF (eof(10)) 830,190
  190   READ (10,230) (fc(i9,j,i),i=1,nds)
        GO TO 240


  200   READ(3,220,IOSTAT=errCode) fx(i9,j),fr(i9,j),fphi(i9,j),fp(i9,j),ft(i9,j),    &
     &   fu(i9,j),jndl(i9,l),jndr(i9,l)
!!!        IF (eof(3)) 830,210
        IF (errCode < 0) GO TO 830
  210   READ (3,230) (fc(i9,j,i),i=1,nds)
  220   FORMAT (6e12.5,2i4)
  230   FORMAT (8e10.3)
  240   IF (jndr(i9,l)-1) 250,270,250
  250   IF (j.eq.2.and.i9.eq.2) fr(2,1)=dum2
        IF (j.eq.2.and.i9.eq.2) fx(2,1)=dum1
        rrr=.5*(fr(i9,j)+fr(i9,j-1))
        xxx=.5*(fx(i9,j)+fx(i9,j-1))
        smw=0.0
        fcp=0.0
        DO 260 i=1,nds
          smw=smw+fc(i9,j,i)/mw(i)
          fcp=fcp+fc(i9,j,i)*a(i,3)
  260   CONTINUE
        smw=1.0/smw
        densty=smw*fp(i9,j)/(rv*ft(i9,j))
        fp(i9,j)=densty*(rrr**2+xxx**2)
  270   IF (jndr(i9,l).ge.1) GO TO 280
        j=j-1
        GO TO 320
  280   IF (j.ne.1) GO TO 290
        IF (jndl(i9,l).ne.1) WRITE (6,800)
        IF (jndl(i9,l).ne.1) GO TO 770
        GO TO 300
  290   IF (jndl(i9,l).eq.1) GO TO 330
  300   CONTINUE
  310   FORMAT (i3,1p6e15.6,2i6)
  320   locf(i9,l)=j
      nstrm(i9)=j
      ntot(i9)=mxtot
      GO TO 350
  330 dum1=fx(i9,j)
      dum2=fr(i9,j)
      dum3=fphi(i9,j)
      dum4=fp(i9,j)
      dum5=ft(i9,j)
      dum6=fu(i9,j)
      idum1=jndl(i9,j)
      idum2=jndr(i9,j)
      DO 340 i=1,nds
        cdum(i)=fc(i9,j,i)
  340 END DO
      nstrm(i9)=j-1
      ntot(i9)=l-1
  350 IF (i9.ne.2) GO TO 140
  360 notube=0
      l=0
  370 l=l+1
      IF (jndr(2,l).ge.1) GO TO 380
      IF (l.lt.ntot(2)) GO TO 370
      WRITE (6,810)
      GO TO 770
  380 r1=l
  390 l=l+1
      IF (jndr(2,l).gt.1) GO TO 400
      IF (l.lt.ntot(2)) GO TO 390
      IF (notube.eq.0) WRITE (6,820)
      IF (notube.eq.0) GO TO 770
      GO TO 550
  400 r3=l
  410 FORMAT (' POINTS FOUND    ON RIGHT. R1='i3,' R3='i3)
      notube=1
!     TEST TO SEE IF POINT IS BOUNDED BY STREAMTUBE
      r1s2=locf(2,r1)
      r2s2=locf(2,r3)
      IF (isb.eq.1) r1s2=locf(2,r1)-1
      IF (isb.eq.1) r2s2=locf(2,r3)-1
      CALL locate (fx(2,r1s2), fr(2,r1s2), xtry, rtry, fx(2,r2s2), fr(2,&
     &r2s2), fx(1,r1), fr(1,r1), ifgood)
      GO TO (430,420),ifgood
  420 r1=r3
      GO TO 390
  430 CALL locate (fx(1,r3), fr(1,r3), xtry, rtry, fx(1,r1), fr(1,r1),  &
     &fx(2,r2s2), fr(2,r2s2), ifgood)
      GO TO (440,420),ifgood
  440 IF (r3-r1.ne.1) GO TO 450
      l1=r1
      l2=r3
      GO TO 480
  450 r2m1=r3-1
      r1p1=r1+1
      DO 460 j=r1p1,r2m1
        xxx=fx(2,r1s2)+(fx(2,r2s2)-fx(2,r1s2))*(fx(1,j)-fx(1,r1))/(fx(1,&
     &   r3)-fx(1,r1))
        rrr=fr(2,r1s2)+(fr(2,r2s2)-fr(2,r1s2))*(fr(1,j)-fr(1,r1))/(fr(1,&
     &   r3)-fr(1,r1))
        CALL locate (fx(1,j), fr(1,j), xtry, rtry, fx(1,j-1), fr(1,j-1),&
     &    xxx, rrr, ifgood)
        GO TO (470,460),ifgood
  460 END DO
      l1=r2m1
      l2=r3
      GO TO 480
  470 l1=j-1
      l2=j
!     PROJECT POINT BACK TO ORTHOGONAL SURFACES ALONG NORMAL FROM POINT
!     TO THE SURFACES
  480 CONTINUE
  490 FORMAT (' POINT BRACKETED BOTH LEFT AND RIGHT. R1='i3,' R3='i3,' L&
     &1='i3,' L2='i3)
  500 CALL orthog (2, r1s2, r2s2)
      CALL orthog (1, l1, l2)
      ratio1=dptsfc(2)/(dptsfc(2)+dptsfc(1))
      tsh1=tsfc(2)+ratio1*(tsfc(1)-tsfc(2))
      psh1=psfc(2)+ratio1*(psfc(1)-psfc(2))
      ush1=usfc(2)+ratio1*(usfc(1)-usfc(2))
      phish1=phisfc(2)+ratio1*(phisfc(1)-phisfc(2))
      DO 510 i=1,nds
  510   csh1(i)=csfc(2,i)+ratio1*(csfc(1,i)-csfc(2,i))
      densty=psh1/(rtry**2+xtry**2)
      smw=0.0
      fcp=0.0
      DO 520 i=1,nds
        smw=smw+csh1(i)/mw(i)
  520   fcp=fcp+csh1(i)*a(i,3)
      smw=1.0/smw
      zmwsh1=smw
      psh1=densty*rv*tsh1/smw
      psh11=psh1
      tsh11=tsh1
      ush11=ush1
      flang=phish1
      DO 530 i=1,nds
  530   csh11(i)=csh1(i)
  540 FORMAT (/' FINAL PROPERTIES FOR SHOCK POINT.'/' T='e15.6,' P='    &
     &e15.6,' U='e15.6,' PHI='e15.6)
      GO TO 760
  550 nstrm2=nstrm(i9)
      ntot2=ntot(i9)
      DO 560 j=1,nstrm2
        fx(1,j)=fx(2,j)
        fr(1,j)=fr(2,j)
        fphi(1,j)=fphi(2,j)
        fp(1,j)=fp(2,j)
        ft(1,j)=ft(2,j)
        fu(1,j)=fu(2,j)
        DO 560 i=1,nds
  560   fc(1,j,i)=fc(2,j,i)
      DO 570 l=1,ntot2
        jndl(1,l)=jndl(2,l)
  570   jndr(1,l)=jndr(2,l)
      fx(2,1)=dum1
      fr(2,1)=dum2
      fphi(2,1)=dum3
      fp(2,1)=dum4
      ft(2,1)=dum5
      fu(2,1)=dum6
      jndl(2,1)=idum1
      jndr(2,1)=idum2
      DO 580 i=1,nds
        fc(2,1,i)=cdum(i)
  580 END DO
      nstrm(1)=nstrm(2)
      ntot(1)=ntot(2)
      GO TO 140
!     NON-UNIFORM FLOW WITH ANALYTICAL COEFFICIENTS
  590 IF (iin) 670,600,670
  600 READ (5,610) thetoo,rhoe,poe,toe,rex,nend,nemd,nent,nemt,nenp,    &
     &nemp
  610 FORMAT (5e12.6,6i3)
      DO 620 i=1,nend
        READ (5,660) (arho(i,j),j=1,nemd)
  620 END DO
      DO 630 i=1,nent
        READ (5,660) (athet(i,j),j=1,nemt)
  630 END DO
      DO 640 i=1,nenp
        READ (5,660) (ap(i,j),j=1,nemp)
  640 END DO
      smw=0.0
      DO 650 i=1,nds
  650   smw=smw+mw(i)*coo(i)
      rhoo=smw*poe/(rv*toe)
      iin=1
      tht=thetoo-0.5*pi
      xoe=rex*tan(tht)
  660 FORMAT (6e12.6)
  670 eta=rex/sqrt(rtry**2+(xtry-xoe)**2)
      IF (xtry-xoe) 680,690,700
  680 theta=pi-atan(rtry/(xoe-xtry))
      GO TO 710
  690 theta=0.5*pi

  700 theta=atan(rtry/(xtry-xoe))
  710 cp=0.0
      DO 720 i=1,nds
        csh1(i)=coo(i)
  720   cp=cp+a(i,3)*coo(i)
      gam=cp/(cp-rcon/smw)
      sumrho=0.0
      sumphi=0.0
      sump=0.0
      DO 730 i=1,nend
        DO 730 j=1,nemd
  730   sumrho=sumrho+arho(i,j)*eta**(j-1)*cos(pi*(float(i)-0.5)*theta/ &
     &   thetoo)
      DO 740 i=1,nent
        DO 740 j=1,nemt
  740   sumphi=sumphi+athet(i,j)*eta**(j-1)*cos(.5*pi*(float(2)-.5)*    &
     &   theta/thetoo)
      DO 750 i=1,nenp
        DO 750 j=1,nemp
  750   sump=sump+ap(i,j)*eta**(j-1)*cos(.5*pi*(float(i)-.5)*theta/     &
     &   thetoo)
      rhosh=rhoe*eta**2*sumrho**(2./(gam-1.))
      phish1=theta*(1.+eta*sumphi)
      posh=poe*(1.+eta*sumphi)
      rhoosh=posh*smw/(rv*toe)
      psh1=(rhosh/rhoosh)**gam*posh
      tsh1=psh1*smw/(rv*rhosh)
      ush1=sqrt(2.*hj*cp*(toe-tsh1))
  760 RETURN
  770 WRITE (6,780)
      CALL pdump (a(1,1), n, 5)
  780 FORMAT (' ERROR IN SHOPRO')
  790 FORMAT (2x,4i4,'I9 L JNDL JNDR')
  800 FORMAT (2x,'STATE 157')
  810 FORMAT (2x,'STATE 175')
  820 FORMAT (2x,'STATE 179')
  830 WRITE (6,410) r1,r3
      WRITE (6,130) xtry,rtry
      WRITE (6,160) ictsfc
      WRITE (6,790) i9,l,jndl(i9,l),jndr(i9,l)
      WRITE (6,310) l,fx(i9,j),fr(i9,j),fphi(i9,j),fp(i9,j),ft(i9,j),   &
     &fu(i9,j),jndl(i9,l),jndr(i9,l)
      WRITE (6,490) r1,r3,l1,l2
      WRITE (6,540) tsh1,psh1,ush1,phish1
      END FILE 13
      STOP   ! CALL exit
      END  Subroutine Shopro

      SUBROUTINE stable
!     THIS SUBROUTINE DETERMINES STABLE STEPPING DISTANCE AND PUNCHES
!     OUTPUT DATA WHEN CALLED FOR
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      DIMENSION index2(120)
      COMMON /intl/ xlx,iki,rex
      DIMENSION ca(30)
      IF ((isb.eq.-1).and.(x2(2).eq.xlx)) GO TO 460
      IF (ll) 10,10,40
!     INDEXING
   10 iextra(6)=0
      DO 20 i=1,kmax
   20   index2(i)=i
      max=kmax+1
      DO 30 i=max,120
   30   index2(i)=0
   40 k=2
      i=2
      IF (kmn.gt.1) i=1
      dels=1.0e10
   50 IF (iextra(5)-lplane) 60,60,80
   60 IF (index2(i)) 70,80,80
   70 i=i+1
      GO TO 60
!     VISCOUS STABILITY CRITERION
   80 del1=dely(k)*re(k)/2.0
!     INERTIAL STABILITY CRITERION
      IF (kmn.le.1) GO TO 110
!     CARD OMITTED
      IF (k.gt.kmn) GO TO 110
      IF (ll.ne.0) GO TO 100
      r2(k)=r(k)
      y2(k)=y(k)
      x2(k)=x(k)
      phi2(k)=phi(k)
      u2(k)=u(k)
      p2(k)=p(k)
      t2(k)=t(k)
      rho2(k)=rho(k)
      r(1)=0.0
      r2(1)=0.0
      y(1)=0.0
      y2(1)=0.0
      DO 90 i=1,nds
   90   c2(i,k)=c(i,k)
  100 kmnn=kmn+1
      amak=ma(kmnn)
      IF (amak-1.0-epslon) 120,120,140
  110 IF (ma(k)-1.0-epslon) 150,150,170
  120 WRITE (6,130) amak
  130 FORMAT (/'AMAK=',e12.5)
      amak=1.0+epslon
  140 del2=0.5*dely(k)*(amak**2-1.0)**(0.5)
      GO TO 180
  150 WRITE (6,160) k,ma(k)
  160 FORMAT (/'K=',i2,5x,'M=',e12.5)
      ma(k)=1.0+epslon
  170 del2=.5*dely(k)*(ma(k)**2-1.0)**(.5)
  180 CONTINUE
      IF (del1.eq.0.0) GO TO 190
      delss(k)=alphah/(1.0/del1+1.0/del2)
      GO TO 200
  190 delss(k)=alphah*del2
  200 CONTINUE
      IF ((isb.eq.1).and.(k.le.3)) GO TO 290
!     COMBINING SMALL TUBES
      IF (ll) 230,210,210
  210 fl=ll
      fplane=lplane
      fmult=fl/fplane
      gstep=fstep
      IF (extra(4).lt.1.e-5) GO TO 220
      fstep=gstep-(gstep-extra(3))*(1.-fmult)/extra(4)
      GO TO 230
  220 fstep=gstep
  230 l=k
      CALL combo (l)
      fstep=gstep
      IF (l-k) 240,290,250
  240 k=l
      IF (iextra(5)-lplane) 260,260,80
  250 k=l-1
      i=i+1
      IF (iextra(5)-lplane) 260,260,80
  260 ij=i
  270 ij=ij-1
      IF (index2(ij)) 270,280,280
  280 index2(ij)=-1
      iextra(6)=iextra(6)+1
      GO TO 60
!     AREA CHANGE LIMITATION
  290 dlna=(sin(.5*(phi(k)+phi(k-1)))/(.5*(r(k)+r(k-1)))*delta+(phi(k)- &
     &phi(k-1))/dely(k))*delss(k)
      IF (abs(dlna)-atol) 310,310,300
  300 delss(k)=delss(k)*atol/abs(dlna)
  310 IF (delss(k)-dels) 320,330,330
  320 IF ((isb.eq.1).and.(k.le.3)) GO TO 330
      dels=delss(k)
      iextra(7)=k
  330 IF (iextra(5)-lplane) 340,340,350
  340 index2(i)=k
      iextra(8)=i
      i=i+1
  350 k=k+1
      IF (k-kmax) 50,50,360
  360 IF (kmax-28) 370,370,380
!     COMBINING SMALLEST TUBE WHEN NUMBER OF TUBES GETS TOO LARGE
  370 itube=fstep
      IF (kmax-itube) 450,450,380
  380 l=iextra(7)
      k=l
      store1=grad
      store2=fstep
      store3=y(kmax)
      grad=10.**10
      fstep=10.**(-10)
      y(kmax)=10.**20
      CALL combo (l)
      grad=store1
      fstep=store2
      y(kmax)=store3
      k=kmax+1
      IF (iextra(5)-lplane) 390,390,730
  390 IF (l-iextra(7)) 410,400,400
  400 i=iextra(8)
      GO TO 420
  410 i=iextra(8)-1
  420 iextra(6)=iextra(6)+1
      jmax=kmax+iextra(6)
      jmin=i+1
      j=jmax
  430 index2(j)=index2(j-1)
      j=j-1
      IF (j-jmin) 440,430,430
  440 index2(i)=-1
  450 IF (iextra(5).gt.lplane) GO TO 730
      IF (x2(2).lt.xlx) GO TO 730
      IF (x2(2).eq.xlx) GO TO 460
      IF (mod(ll,iextra(5))) 730,460,730
!     PUNCHING OUTPUT CARDS
  460 i=1
      p(1)=p(2)
      t(1)=t(2)
      u(1)=u(2)
      DO 470 n=1,nds
  470   c(n,1)=c(n,2)
      nddp=nds
      IF (ipart.eq.2) nddp=nds+1
      IF (ipart.eq.2) c(nddp,1)=0.
      IF (iki.eq.9) int=1
      WRITE (6,480) iki,mmax,ll,x(i)
  480 FORMAT (6x,'IKI=',i4,4x,'MMAX=',i4,4x,'LL=',i5,4x,'X(1)=',e10.4)
      IF (int.eq.0.or.iki.ne.0) WRITE (10,490) x(i),r(i),phi(i),p(i),   &
     &t(i),u(i),i,i
      IF (isb.lt.0) WRITE (3,490) x(i),r(i),phi(i),p(i),t(i),u(i),i,i
  490 FORMAT (6e12.5,2i4)
      IF (int.eq.0.or.iki.ne.0) WRITE (10,500) (c(j,i),j=1,nddp)
      IF (isb.lt.0) WRITE (3,500) (c(j,i),j=1,nddp)
  500 FORMAT (8e10.3)
      kmxx=kmax+iextra(6)
      IF (ipart.le.1) GO TO 630
      kbao=kmax-1
      DO 620 i=2,kbao
        xmv=0.0
        xm=0.0
        xmm=0.0
        xmvv=0.
        rhtt=0.0
        eu=0.0
        ee=0.0
        ii=i+1
        DO 510 j=1,npg
          xm=xm+rhp2(ii,j)*w2(ii,j)**2
          xmm=xmm+rhp2(ii,j)*(w2(ii,j)+v2(ii,j))
          xmv=xmv+rhp2(ii,j)*v2(ii,j)**2
          xmvv=xmvv+rhp2(ii,j)*v2(ii,j)
          eu=eu+rhp2(ii,j)*(w2(ii,j)**3/2.0+v2(ii,j)**3/2.0)
          rhtt=rhtt+rhp2(ii,j)
  510     ee=ee+rhp2(ii,j)*tp2(ii,j)*cs*(w2(ii,j)+v2(ii,j))
        cw=xm+xmv+rho(i)*u(i)**2
        aw=xmm-xmvv+rho(i)*u(i)
        bw=xmvv
        up=cw*aw/(aw**2+bw**2)
        vp=bw*cw/(aw**2+bw**2)
        rhop=(aw**2+bw**2)/cw
        hp=(rho(i)*u(i)*(h(i)+u(i)**2/(2.*hj))+eu/hj+ee)/(aw+bw)-.5/hj* &
     &   (cw**2/(aw**2+bw**2)-cw**2*aw*bw/(aw**2+bw**2)**2)
        phia=phi(i)+atan(vp/up)
        ua=sqrt(up**2+vp**2)
        rhtot=rho(i)+rhtt
        WRITE (6,520) up,vp,rhop,hp,rhtot,rhtt
  520   FORMAT (6e12.5)
        ndsp=nds+npg
        nddp=nds+1
        DO 530 j=nddp,ndsp
          jj=j-nds
  530     ca(j)=rhp2(ii,jj)/rhtot
        DO 540 j=1,nds
  540     ca(j)=c(j,i)*rho(i)/rhtot
        tl=abs(1.e-3*hp)
        ta=t(i)
        DO 570 it=1,100
          ha=0.0
          dgh=0.0
          DO 550 j=1,nds
            DO 550 jj=1,nn
            ha=ha+a(j,jj)*ta**(jj-2)*ca(j)
  550       dgh=dgh+a(j,jj)*float(jj-2)*ta**(jj-3)*ca(j)
          hpp=0.0
          dhp=0.0
          DO 560 j=nddp,ndsp
            hpp=hpp+cs*ta*ca(j)
  560       dhp=dhp+cs*ca(j)
          g=ha+hpp-hp
          IF (abs(g).lt.tl) GO TO 590
          dg=dhp+dgh
          ta=ta-g/dg
  570   CONTINUE
        WRITE (6,580) it,ha,hpp,dgh,dhp
  580   FORMAT (8x,'ITERATION DOES NOT CONVERGE TO CORRECT TEMPERATURE',&
     &   i4,2x,'HA= ',1pe10.3,5x,'HPP=',e10.3,5x,'DGH=',e10.3,5x,'DHP=',&
     &   e10.3)
  590   CONTINUE
        zmw=0.0
        DO 600 j=1,nds
  600     zmw=zmw+ca(j)/mw(j)
        pa=rhtot*rv*ta*zmw
        nddd=nddp+1
        DO 610 j=nddd,ndsp
  610     ca(nddp)=ca(nddp)+ca(j)
        IF (int.eq.0.or.iki.ne.0) WRITE (10,490) x(i),r(i),phia,pa,ta,  &
     &   ua,i,i
        IF (int.eq.0.or.iki.ne.0) WRITE (10,500) (ca(j),j=1,nddp)
  620 END DO
      GO TO 700
  630 kchu=kmxx-1
      IF (iextra(6).eq.0) kchu=kmxx
      DO 690 i=2,kchu
        IF (index2(i)) 640,640,670
  640   jj=i
  650   jj=jj+1
        ii=index2(jj)
        IF (ii) 660,660,680
  660   IF ((kmax.eq.28).and.(index2(jj+1).eq.0)) GO TO 700
        GO TO 650
  670   ii=index2(i)
  680   IF (int.eq.0.or.iki.ne.0) WRITE (10,490) x(ii),r(ii),phi(ii),   &
     &   p(ii),t(ii),u(ii),i,index2(i)
        IF (isb.lt.0) WRITE (3,490) x(ii),r(ii),phi(ii),p(ii),t(ii),    &
     &   u(ii),i,index2(i)
        IF (int.eq.0.or.iki.ne.0) WRITE (10,500) (c(j,ii),j=1,nds)
        IF (isb.lt.0) WRITE (3,500) (c(j,ii),j=1,nds)
  690 END DO
  700 iextra(6)=0
      DO 710 i=2,kmax
  710   index2(i)=i
      max=kmax+1
      DO 720 i=max,120
  720   index2(i)=0
  730 RETURN
      DO 750 k=2,kmax
        WRITE (6,740) delss(k)
  740   FORMAT (' DELSS=',e11.4)
  750 END DO
      RETURN
      END  Subroutine Stable

      SUBROUTINE sldp (x, a, n)
!     THIS PROGRAM FINDS THE SOLUTIONS TO A SET OF N SIMULTANEOUS LINEAR
!     EQUATIONS BY USING THE GAUSS-GORDAN REDUCTION ALGORITHM WITH THE
!      DIAGONAL PIVOT STRATEGY
      DIMENSION a(21,21), x(21)
      DO 70 k=1,n
        IF (abs(a(k,k)).gt.1.e-10) GO TO 40
        WRITE (6,30) k,a(k,k)
        DO 10 ir=1,n
          WRITE (6,20) (a(ir,jr),jr=1,n),x(ir)
   10   CONTINUE
   20   FORMAT (1x,10e12.3)
   30   FORMAT (' ERROR--- SMALL PIVOT ',i5,e12.5)
        STOP
   40   kp1=k+1
        DO 50 j=kp1,n
   50     a(k,j)=a(k,j)/a(k,k)
        x(k)=x(k)/a(k,k)
        a(k,k)=1.0
        DO 70 i=1,n
        IF (i.eq.k.or.a(i,k).eq.0.) GO TO 70
        DO 60 j=kp1,n
   60     a(i,j)=a(i,j)-a(i,k)*a(k,j)
        x(i)=x(i)-a(i,k)*x(k)
        a(i,k)=0.
   70 CONTINUE
      RETURN
      END  Subroutine Sldp

      SUBROUTINE step
!     THIS SUBROUTINE CALCULATES STATE PROPERTIES IN A STREAMTUBE
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      REAL massmd,l1,l2
      COMMON /pt1/ xmomv(8),xmomp(8),enep1(8),rep(30,8),v(30,8),w(30,8),&
     &v2(30,8),w2(30,8),trbdy(8),tp(30,8),tp2(30,8),rp(8),rhp(30,8),    &
     &rhp2(30,8),enep2(8),nbl(8),um(30),volp(8),wtp(8),hc(8),deng(30,8),&
     &icond(30,8)
      COMMON /pt2/ fff,ffg,gama,umi,pq,npg,kmx,dndsp(8)
      COMMON /pt3/ cl,cs,tps,rhss,wt,htran,sig,ep,ipart,ikine,tr(30),   &
     &sump,sumpv,sumpe
      COMMON /turple/ del2,cpshw,gamshw,zmwsh1
      COMMON /chem1/ zid(5),irr(30),irt(30),rc(30,3),irrr(30,5),av,     &
     &cm(21,21),fi(21),wp(21),wm(21),wdot(21,30),csi(21),qx(21)
      COMMON /dsl/ rdsl,phdsl
      COMMON /sub/ xdis,l1,l2,phio,rstar,rmd,massmd,rnoz,thtp,thru
      INTEGER thru
      DIMENSION zao(25)
      EQUIVALENCE (ivisc,iextra(1))
      EQUIVALENCE (kalarm,iextra(20))
!     EVALUATION OF TRANSPORT TERMS
      IF (ivisc) 10,10,50
   10 rhsmom=0.0
      rhsen=0.0
      DO 20 i=1,nds
   20   rhs(i)=0.0
      IF (iflag) 30,40,30
   30 ubar(k)=u2(k)
      GO TO 80
   40 ubar(k)=u(k)
      GO TO 80
   50 CALL visco
      IF (k.lt.kmax) GO TO 70
      tbbar=tabar
      pbbar=pabar
      DO 60 i=1,nds
   60   xbbar(i)=xabar(i)
   70 CONTINUE
   80 aold=pi*(r(k)+r(k-1))**delta*(y(k)-y(k-1))
      abar=.5*(aa(k)+aold)
!     CARD OMITTED
!     CARD OMITTED
!     CARD OMITTED
      IF (t2(k).lt.400.) GO TO 90
      IF (ikine.ne.0.and.iflag.eq.1) GO TO 110
      IF (ikine.ne.0.and.iflag.eq.1.and.isb.eq.1) GO TO 110
   90 DO 100 i=1,nds
        wdot(i,k)=0.
  100   c2(i,k)=c(i,k)+rhs(i)/mdot(k)
      GO TO 120
  110 CALL chem (ikine)
  120 CONTINUE
      DO 130 i=1,nds
        IF (c2(i,k)) 140,130,130
  130 END DO
      GO TO 150
  140 dels=dels*.5
      icount=icount+1
      kalarm=2
      GO TO 670
  150 zmw=0.0
      DO 160 i=1,nds
  160   zmw=zmw+c2(i,k)/mw(i)
      zmw=1.0/zmw
      zmwk(k)=zmw
      DO 170 j=1,nn
        za(j)=0.0
        zao(j)=0.
        DO 170 i=1,nds
        zao(j)=zao(j)+a(i,j)*c(i,k)
  170   za(j)=za(j)+a(i,j)*c2(i,k)
      iter=0
      us=ubar(k)
      cp=0.0
      cpold=0.
      DO 180 j=1,nn
        cpold=cpold+float(j-2)*zao(j)*t(k)**(j-3)
  180   cp=cp+za(j)*float(j-2)*t(k)**(j-3)
  190 b1=(rhsmom+mdot(k)*u(k))/(betap*abar)+p(k)
      b1=b1-dls(k)*sump/betap
      b2=mdot(k)/(betap*abar)
      b3=zmw*aa(k)/rv/mdot(k)
      IF ((ma(k).le.1.002).and.(k.gt.kmn)) GO TO 620
      IF (iconst.ne.0) GO TO 230
!     CARD OMITTED
!     CONSTANT GAMMA
      b5=hj*cp*b3*b1
      b6=2.0*hj*cp*b3*b2-1.0
      b4=u(k)**2+2.0*hj*(rhsen/mdot(k)+cpold*t(k)+zao(2)-za(2))
      b4=b4+sumpe*abar*dls(k)/mdot(k)
      discr=b5**2-b6*b4
      IF (discr) 200,220,220
  200 IF (icount-9) 210,210,630
  210 icount=icount+1
      dels=dels/2.0
      kalarm=2
      GO TO 670
  220 us=(b5+sqrt(discr))/b6
      GO TO 340
!     VARIABLE HEAT CAPACITY
  230 iter=iter+1
      b4=2.0*hj*rhsen/mdot(k)+2.0*hj*h(k)+u(k)**2
      b4=b4+sumpe*abar*dls(k)/mdot(k)
      dum1=0.0
      dum2=0.0
      DO 240 i=1,nn
        dum1=dum1+za(i)*(b3*(b1-b2*us)*us)**(i-2)
  240   dum2=dum2+za(i)*float(i-2)*(b1*us-b2*us**2)**(i-3)*b3**(i-2)
      g=2.0*hj*dum1+us**2-b4
      dgdu=2.0*us+2.0*hj*dum2*(b1-2.0*b2*us)
      deltau=-g/dgdu
      IF (abs(deltau/us)-tol) 340,250,250
  250 IF (iter-niter) 330,260,260
  260 IF (icount-9) 320,320,270
  270 IF ((k.gt.2).or.(ma(2).lt.0.9)) GO TO 640
  280 WRITE (6,290)
  290 FORMAT (/5x,'FLOW HAS CHOKED-GO TO FULL SUPERSONIC SOLUTION')
      phi2(2)=0.
      dphids(2)=0.
      kmn=0
      WRITE (6,300) us
  300 FORMAT (/5x,'USOLD=',e12.5)
      us=us*1.005/ma(2)
      ma(2)=1.005
      icount=0
      iter=0
      aa(k)=aold
      abar=aold
      r2(k)=r(k)
      phi(k)=0.
      y2(k)=y(k)
      WRITE (6,310) us,x2(1)
  310 FORMAT (5x,'US=',e12.5,6x,'X=',e12.5)
      GO TO 190
  320 icount=icount+1
      kalarm=2
      dels=dels/2.0
      GO TO 670
  330 us=us+deltau
      GO TO 230
  340 u2(k)=us
      p2(k)=b1-b2*us
      IF ((iflag.eq.1).and.(k.eq.kmn+1)) deltap=p2(k)-p(k)
      t2(k)=b3*p2(k)*us
      rho2(k)=zmw*p2(k)/(rv*t2(k))
      IF (p2(k)) 350,350,380
  350 IF (icount.lt.9) GO TO 320
      WRITE (6,360)
  360 FORMAT ('1 PRESSURE BECAME NEGATIVE')
  370 CALL pdump (a(1,1), n, 1)
      STOP       ! CALL exit
  380 IF (t2(k)) 390,390,410
  390 IF (icount.lt.9) GO TO 320
      WRITE (6,400)
  400 FORMAT ('1 TEMPERATURE BECAME NEGATIVE')
      GO TO 370
!     STREAMTUBE MACH AND REYNOLDS NUMBERS
  410 fmax=1.0
      IF (k.ne.kmax) GO TO 430
      mu2=mu
      kay2=kay
      DO 420 i=1,nds
        d2eff(i)=deff(i)
        DO 420 j=1,nds
        d2ih(i,j)=dih(i,j)
  420 CONTINUE
  430 cp=0.
      DO 440 i=1,nn
  440   cp=cp+float(i-2)*za(i)*t2(k)**(i-3)
      IF (k.eq.2.and.isb.eq.-1) cpshw=cp
      IF (ivisc) 450,560,450
  450 IF (kay2-mu2*cp) 470,470,460
  460 fmax=kay2/(cp*mu2)
  470 IF (mult.ne.0.and.xs(ibig,k+1).gt.0.75) GO TO 530
      imax=nds-1
      DO 520 i=1,imax
        j=i
  480   j=j+1
        IF (d2ih(i,j)*rho2(k)/mu2-fmax) 510,510,490
  490   IF (mult.ne.0.and.xs(ibig,k+1).gt.0.75) GO TO 500
        fmax=d2ih(i,j)*rho2(k)/mu2
        GO TO 510
  500   fmax=d2eff(j)*rho2(k)/mu2
  510   IF (j.lt.nds) GO TO 480
  520 END DO
      GO TO 550
  530 DO 540 i=1,nds
        IF (d2eff(i)*rho2(k)/mu2.gt.fmax) fmax=d2eff(i)*rho2(k)/mu2
  540 END DO
  550 CONTINUE
      re(k)=rho2(k)*u2(k)*dely(k)/(mu2*fmax)
  560 ma(k)=u2(k)*sqrt((cp*zmw-rcon)/(cp*t2(k)*ru))
      IF ((isb.eq.1).and.(k.eq.2)) WRITE (6,590) ma(k)
      IF ((isb.ne.1).or.(k.ne.2).or.(par(3).eq.0.0)) GO TO 580
      ipar3=par(3)
      IF (mod(ll,ipar3).eq.0) WRITE (6,570) x2(k),p2(k),u2(k),t2(k),cp, &
     &zmw
  570 FORMAT (4x,'X=',e11.4,4x,'P=',e11.4,4x,'U=',e11.4,4x,'T=',e11.4,  &
     &4x,'CP=',e11.4,4x,'ZMW=',e11.4)
  580 CONTINUE
  590 FORMAT (/5x,'MACH NUMBER =',e12.5)
      IF (kmn.eq.0.and.isb.eq.1.and.ma(2).lt.1.1) GO TO 600
      IF (kmn.eq.0.or.k.gt.2) GO TO 610
      IF (ma(2).le.0.98.or.iflag.eq.0) GO TO 610
      GO TO 280
  600 u2(2)=1.10*u2(2)/ma(2)
      ma(2)=1.1
      y2(k)=y(k)
      phi(k)=0.
      dphids(2)=0.
      phi2(2)=0.
      r2(k)=r(k)
      abar=aold
      aa(k)=aold
  610 CONTINUE
      IF ((isb.eq.1).and.(k.eq.2).and.(ma(2).ge.1.0)) kmn=0
      gama=cp/(cp-rcon/zmw)
      IF (k.eq.2.and.isb.eq.-1) gamshw=gama
      pq=4.0*gama/(9.*gama-5.0)
      GO TO 670
  620 CONTINUE
!     ERROR MESSAGES
      nucase=2250
      GO TO 660
  630 nucase=730
      WRITE (6,650) nucase
      GO TO 660
  640 nucase=1010
      WRITE (6,650) nucase
  650 FORMAT ('1',15x,' ------------- NUCASE=',i5//)
  660 ierror=3
  670 RETURN
      END  Subroutine Step

      SUBROUTINE transp (tt, pp, kkk)
!     -- THIS SUBROUTINE CALCULATES THE TURBULENT TRANSPORT PROPERTIES -
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      COMMON /turbul/ tle,tpr,eddyk,iturb,delmix
      COMMON /inpux/ itd,x0,y0,frac,kfl
      COMMON /pergff/ pergf
      kay2=kay
      mu2=mu
      d2ih(1,2)=dih(1,2)
      IF (kkk-2) 10,10,90
   10 IF (ll.lt.0) GO TO 80
      IF (iturb-2) 50,20,80
   20 pumax=rho(2)*u(2)
      pumin=rho(2)*u(2)
      DO 40 l=2,kmax
        test=rho(l)*u(l)
        IF (test.gt.pumax) GO TO 30
        IF (test.ge.pumin) GO TO 40
        pumin=test
        GO TO 40
   30   pumax=test
   40 END DO
      delmix=y(kmax)
      mu=(delmix/eddyk)*(pumax-pumin)+.047896
      WRITE (6,140) pumax,pumin,delmix,mu,rho(kmax),u(kmax),rhabar
      GO TO 90
   50 kmm=kmax-1
!     ----- DETERMINE DIVIDING STREAMLINE -----
      tawm=taw(2)
      kdsl=2
      DO 60 if=2,kmm
        IF (abs(taw(if+1)).le.abs(tawm)) GO TO 60
        kdsl=if+1
        tawm=taw(kdsl)
   60 END DO
      rdv=r(kdsl)
!     ----- COMPUTE DELMIX AND MU -----
      perg=pergf*abs(tawm)
      inn=2
      DO 70 if=3,kdsl
        IF (abs(taw(if)).lt.perg) inn=if
   70 END DO
      rinn=r(inn)
      delmix=rdv-rinn
      mu=0.5*delmix*rho(kkk)*abs(u(2)-u(kmax))/eddyk
      IF (delmix.le.0.0) mu=1.0e-5
!     CARD OMITTED
      GO TO 90
   80 mu=1.0e-5
   90 cp=0.0
      DO 120 i=1,nds
        cps(i)=0.0
        DO 100 j=1,nn
  100     cps(i)=cps(i)+float(j-2)*a(i,j)*tt**(j-3)
        IF (ll.lt.0) GO TO 110
        cp=cp+cps(i)*cabar(i)
        GO TO 120
  110   cp=cp+cps(i)*c(i,kkk)
  120 END DO
      kay=cp*mu/tpr
      IF (ibugsh.ne.0) WRITE (6,150) cp,kay,rhabar
      IF (ll.lt.0) GO TO 130
      dih(1,2)=kay*tle/(rhabar*cp)
      GO TO 160
  130 dih(1,2)=kay*tle/(rho(kkk)*cp)
  140 FORMAT (' ','PUMAX=',1pe9.2,2x,'PUMIN,'e9.2,2x,'DELMIX=',e9.2,2x,'&
     &MU=',e9.2,2x,'RHO=',e9.2,2x,'U=',e9.2,2x,'RHABAR=',e9.2)
  150 FORMAT (' ','CP=',1pe9.2,2x,'KAY=',e9.2,2x,'RHABAR=',e9.2)
  160 RETURN
      END  Subroutine Transp

      SUBROUTINE visco
!     THIS SUBROUTINE CALCULATES THE FLUX CONTRIBUTIONS TO THE
!     CONSERVATION EQUATIONS
      COMMON a(20,7),aa(30),alfa(20,20),alphah,alphap,atol,betap,bmix,  &
     &c11(20),c(20,30),c12(20),c2(20,30),cabar(20),cbbar(20),cp,cps(20),&
     &cpsh,csh(20),csh1(20),cstrem(20),mxstrm,d2ih(20,20),deff(25),     &
     &d2eff(20),delta,d11,delss(30),dels,delso,dls(30),dih(20,20),d12,  &
     &dely(30),d13,dpdy(30),d14,dphids(30),epcon,epslon,extra(20),fstep,&
     &fmax,grad,h11,h(30),hh,hj,hpm(20),h2pm(20),h3pm(20),iconst,icount,&
     &ident(20),ierror,iextra(20),iflag,ikind,isb,iptuc,ishock,itype,   &
     &ipd,idiff,k,kay,kays,kay2,klo,kmax,kmn,kup,kw,ll,lplane,ma,mash,  &
     &mdot,mmax,mu0,mu,mu2,mus,mw,mw2,mwsh,nbound,nds,niter,nmax,nn
      COMMON nucase,omega(20),p11,p(30),p12,p2(30),pb(25),pabar,pbbar,  &
     &pbs,p15,phi(40),phib(25),phish1,phstrm,phiw(25),p18,phi2(30),     &
     &phibs,pi,pr(20),psh,psh1,psi,pstrem,pw(25),q11,q(30),qxtr1,qxtr2, &
     &qw(25),r11,r(30),rb(25),rbs,r13,rbar(30),r14,r2(30),rcon,r15,     &
     &re(30),resh,r16,rho(30),r17,rho2(30),rhs(20),rhsen,rhsmom,rhabar, &
     &rhbbar,ru,rsh,rstrem,rv,r19,rw(25),s,sb(25),sc(20),sw(25),s13,    &
     &sx(30),t11,t(30),t12,t2(30),tabar,tbbar,t0(20),t14,taw(30),txtr1, &
     &txtr2,tol,tsh,tsh1,tstrem,tw(25),tws,ts,u11,u(30),u12,u2(30)
      COMMON uxtr1,uxtr2,u14,ubar(30),ush,ush1,ustrem,uw(25),uws,x11,   &
     &x(30),x12,x2(30),xb(25),xbs,xabar(20),xbbar(20),x15,xs(20,30),xsh,&
     &xstrem,xw(25),zmwk(30),y11,y(30),y12,y2(30),yabar,ybbar,za(20),   &
     &z11(20),zj(20,30),zmw,s2x,r2sh,x2sh,fx(2,30),fr(2,30),fphi(2,30), &
     &fp(2,30),ft(2,30),fu(2,30),indl(2,30),indr(2,30),fc(2,30,20),     &
     &card1,m,n,par(10)
      REAL kay,kays(20),kay2,kw,mdot(30),ma(30),mu,mus(20),mu0(20),mu2, &
     &mw(20),mw2,mash,mash1,mwsh,map
      IF (p(k)) 10,30,30
   10 WRITE (6,20)
   20 FORMAT ('1NEGATIVE PRESSURE IN VISCO.')
      CALL pdump (a(1,1), n, 1)
   30 CALL flux
      IF (iflag) 40,50,40
   40 rbar(k)=.5*(r(k)+r2(k))
      GO TO 60
   50 rbar(k)=r(k)
   60 IF (rbar(k-1)) 70,100,90
   70 WRITE (6,80) rbar(k-1),k,delta
      STOP
   80 FORMAT (' IN VISCO RBAR(K-1)=',e16.8,' K=',i5,' DELTA=',e16.8)
   90 dumm=rbar(k-1)**delta
      GO TO 120
  100 IF (delta) 90,110,90
  110 dumm=1.0
!     MOMENTUM TRANSFER
  120 rhsmom=(rbar(k)**delta*taw(k)*dls(k)-dumm*taw(k-1)*dls(k-1))*pi*  &
     &2.0**delta
!     KINETIC ENERGY TRANSFER
      IF (k-2) 130,130,140
  130 dssptn=rbar(k)**delta*taw(k)*(ubar(k)+ubar(k+1))/2.0*dls(k)-dumm* &
     &taw(k)*uws*dls(k-1)
      GO TO 150
  140 dssptn=rbar(k)**delta*taw(k)*(ubar(k)+ubar(k+1))/2.0*dls(k)-dumm* &
     &taw(k-1)*(ubar(k)+ubar(k-1))/2.0*dls(k-1)
!     ENERGY TRANSFER BY DIFFUSION
  150 dfflx=0.0
      DO 200 i=1,nds
        h3pm(i)=h2pm(i)
        h2pm(i)=hpm(i)
        hpm(i)=0.0
        IF (k+1-kmax) 160,190,190
  160   IF (iflag) 180,170,180
  170   hpm(i)=cps(i)*t(k+1)
        GO TO 190
  180   hpm(i)=cps(i)*.5*(t(k+1)+t2(k+1))
  190   dfflx=dfflx+.5*rbar(k)**delta*zj(i,k)*dls(k)*(hpm(i)+h2pm(i))-  &
     &   .5*dumm*zj(i,k-1)*dls(k-1)*(h2pm(i)+h3pm(i))
!   SPECIES TRANSFER
  200   rhs(i)=(rbar(k)**delta*zj(i,k)*dls(k)-dumm*zj(i,k-1)*dls(k-1))* &
     &   pi*2.0**delta
!   TOTAL ENERGY TRANSFER
      rhsen=(dfflx+dssptn/hj+rbar(k)**delta*q(k)*dls(k)-dumm*q(k-1)*    &
     &dls(k-1))*pi*2.0**delta
      RETURN
      END  Subroutine Visco

