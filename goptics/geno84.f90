!          GENOPTICS - A GENERAL OPTICAL SYSTEMS EVALUATION PROGRAM     
!                                                                       
!          THIS VERSION ADAPTED FOR THE DEC VAX 11-780 USING VMS        
!                                                                       
!          GRAPHICS MODIFIED TO USE THE ZETA PLOT PACKAGE               
!                                                                       
!          ROUTINE STRUCTURE OF GENOPTICS     DEC VAX 11-780 VERSION    
!                                                                       
!             ********** MAY 1984 VERSION **********                    
!                                                                       
!     MAIN   - CARD INPUT AND DECODING ROUTINE                          
!     BESJN  - CALCULATES BESSEL FUNCTIONS FOR MTF                      
!     CNTCOM - COUNTS THE NUMBER OF COMMAS IN A CARD IMAGE              
!     FINDE  - DECODES ONE FLOATING POINT NUMBER FROM CARD IMAGES       
!     FINRAY - CREATES RAY LATTICE WITH RAYS COVERING EQUAL SOLID ANGLES
!     FLOTIN - DECODES PART OF A CARD IMAGE INTO FLOATING POINT ARRAY   
!     FOCUS  - DETERMINES IMAGE POSITION FOR MINIMUM SPOT SIZE          
!     FREARA - FLOATING POINT DECODE MONITOR ROUTINE                    
!     HEADIN - LEADS PAGES WITH PAGE HEADING INFORMATION                
!     LENSCL - SCALES THE CURRENT LENS SYSTEM BY A USER-SPECIFIED FACTOR
!     MAVEC  - PERFORMS VECTOR-MATRIX MULTIPLICATION                    
!     MTF    - MTF (MULTIPLE TRANSFER FUNCTION)                         
!     PARAX  - PERFORMS PARAXIAL RAY TRACE                              
!     PREPRT - PRINTS A PRESCRIPTION MATRIX                             
!     RED    - COMPUTES   RED (RADIAL ENERGY DISTRIBUTION)              
!     ROTM   - CONSTRUCTS MATRIX FOR ROTATION TRANSFORMATIONS           
!     SKEW   - PERFORMS REAL RAY TRACE OF SPECIFIED RAYS                
!     SRFSAG - COMPUTES THE SAG OF A ROTATIONALLY SYMMETRIC SURFACE     
!     SURFNO - COMPUTES THE SURFACE NO. TO WHICH INPUT DATA CARD APPLIES
!                                                                       
!     SURTYP -                                                          
!       THIS ROUTINE IS CALLED FROM PARAX AND ADDS THE CAPABILITY       
!       OF HANDLING SURFACES ENTIRELY DESCRIBED BY A POLYNOMIAL         
!       EXPRESSION AND CONIC SURFACES WHICH HAVE DEFORMATION CONSTANTS  
!     GRAPH -                                                           
!                ROUTINE GRAPH PLOTS SPOT,RED, AND MTF DATA             
!     FNDSPT -                                                          
!                ROUTINE SPOT PERFORMS SPOT PLOT AND PRINT              
!                CALCULATIONS                                           
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!                                                                       
!         SOME FLAGS - USED TO CONTROL CALCULATION INTERNALLY           
!                        SET ACCORDING TO INPUT                         
!                                                                       
!             IFLAG  - STORAGE FLAG FOR  HXINIT, HXDEL, HYINIT, HYDEL   
!                    = 0 => STORED AS ANGLES IN DEGREES                 
!                    = 1 => STORED AS LINEAR OBJECT SIZES               
!                                                                       
!             OFLAG  = SYSTEM UNITS FLAG                                
!                    = 1 => MM                                          
!                    = 2 => CM                                          
!                    = 3 => INCHES                                      
!                                                                       
!  /PMATX/    UFLAG  - UNITS FLAG                                       
!                    = 1 => MM                                          
!                    = 2 => CM                                          
!                    = 3 => INCHES                                      
!                    = 4 => TANGENTS OF OUTPUT ANGLES WITH RESPECT      
!                           TO OPTICAL AXIS                             
!                    = 5 => ANGLES OF INCIDENCE WITH RESPECT TO NORMAL  
!                           DEGREES.                                    
!                                                                       
!                                                                       
!    /CSPOT/     IOPA - SPOT PLOT (UNITS) SWITCH                        
!    /CRED/      IOPB - RADIAL ENERGY DIST PLOT SW                      
!    /CMTF/      IOPC - MODULATION TRANSFER FUNCTION PLOT SW            
!                                                                       
!     OPTNA   OPTION A   PRINTS PRIME IMAGE COORD ONLY                  
!     OPTNB   OPTION B   PRINTS COORD AND COSINES IN EP AND IMAGE       
!     OPTNC   OPTION C   PRINTS COORD FOR SURFACES                      
!     OPTND   OPTION D   CAUSES PRESCRIPTION MATRIX PRINT               
!     OPTNE   OPTION E   PRINTS RED TABLES                              
!     OPTNF   OPTION F   PRINTS MTF TABLES                              
!     OPTNG   OPTION G   CAUSES ANALYSIS FOR EACH HEIGHT, COLOR         
!                                                                       
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!                                                                       
!                     COMMONS                                           
!                                                                       
!                                                                       
!   /PMATX/ -  most of these are set in MAIN as the                     
!              result of input cards.                                   
!                                                                       
!             XSMIN       x(min) for spot diagrams    PLOT SCALE        
!             XSMAX       x(max) for spot diagrams    PLOT SCALE        
!             YSMIN       y(min) for spot diagrams    PLOT SCALE        
!             YSMAX       y(max) for spot diagrams    PLOT SCALE        
!                                                                       
!?            FCODE       focusing code               FOCUS             
!             FXYJ        historical variable                           
!             FXY         =RHO, set in SKEW and PARAX                   
!             FXNU        =RHO/(S-D);   "                               
!             FBY         =0. always                                    
!             FBNU        =-HYMAX/(S-D)                                 
!                                                                       
!             S           object - pupil distance     TH                
!             D           location of 1st real        SAY               
!                         surface with respect to                       
!                         the entrance pupil                            
!                                                                       
!             RHO         radius of entrance pupil    SAY               
!             UFLAG  -    Units flag                                    
!                         = 1 => millimeters          UNITS MM          
!                         = 2 => centimeters          UNITS CM          
!                         = 3 => inches               UNITS INCHES      
!                         = 4 => output angle tangent MODE AFOCAL       
!                         = 5 => incident angles(rad) MODE ANGLES       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!             RNOBJ       Number of object points     SCX or SCY        
!             HYINIT      Initial Y-object height     SCY               
!             HYDEL       Y-object height increment   SCY               
!             HXINIT      Initial X object height     SCX               
!             HXDEL       X-object height increment   SCX               
!                 Related IFLAG = 0 => H?INIT, H?DEL                    
!                                      stored as ang-                   
!                                      les in degrees                   
!                               = 1 => linear measure                   
!                                      (system units)                   
!                                                                       
!             APSTOP      System aperture stop        ASTOP             
!                         surface number                                
!             SMAX        MTF maximum spatial         PLOT GOTF or      
!                         frequency.                  PRINT GOTF        
!             RSMAX       RED maximum radius          PLOT RED or       
!                                                     PRINT RED         
!             FOCL        System focal length         CFL               
!             OBJN(3)     Indices of refraction of    GLASS or          
!                         object medium (3 colors)    AIR               
!             DELIMP      separation between image    IMAGE SEP         
!                         surfaces                                      
!             FPLANE      position of first image     IMAGE FIRST       
!                         plane w.r.t. prime plane                      
!             FAKEA       Surface number for STAT     STAT              
!                         output.                                       
!                     REL ATED:  IPR20 = 0 => dev. 6                    
!                                      = 1 => dev.20                    
!                                                                       
!             C(40)       Vertex curvature            CV                
!             T(40)       V-V' thickness              TH                
!             R(40)       Vertex radius of curvature  RD                
!             CONIC(40)   Conic Constant              CC                
!             FN(40,3)    Indices of refraction       GLASS or          
!                                                     AIR               
!                                                                       
!             FMASK(40)   = 1 => clear aperture       CLAP              
!                         = -1 => obscuration         COBS              
!                         = R > 0 => clear circular   CLAP              
!                         = -R <0 => circular obscur  COBS              
!             FAKEC(40)   = 0 => circular             CLAP or COBS      
!                         = 1 => rectangular mask     CLAP RECT or      
!                                                     COBS RECT         
!                         =-1 => Elliptical Mask      CLAP ELIP or      
!                                                     COBS ELIP         
!                                                                       
!                                                                       
!?            FAKEB(40)                                                 
!             XMN(40)     center of circle OR         CLAP,             
!                         left most X                 CLAP RECT         
!                         center of ellipse           CLAP ELIP         
!             XMX(40)     not used                    CLAP,             
!                         right most X                CLAP RECT         
!                         x-radius                    CLAP ELIP         
!             YMN(40)     same as XMN but for Y                         
!             YMX(40)     same as XMX but for Y                         
!                                                                       
!                                                                       
!             XDISP(40)   X decenter                  DEC               
!             YDISP(40)   Y decenter                  DEC               
!                                                                       
!             TILTX(40)   X tilt                      TILT              
!             TILTY(40)   Y tilt                      TILT              
!             TILTZ(40)   Z tilt                      TILT              
!                                                                       
!             ORDN(40,3)  grating orders (3 colors)   GORD              
!             SIDE(40)    direction of incidence for  CONV              
!                         ray                         CONC              
!             RDSPAC(40)  grid spacings               GRATX             
!                                                     GRATY             
!                                                                       
!                                                                       
!?            Y0(40)                                                    
!             SXY(40)     height of axial ray in      PY                
!                         thickness solve                               
!             SXNU(40)    angle of axial ray in       PIY               
!                         curvature solve                               
!             COEF(40, 4) higher order surface        ASPH              
!                         coefficients                                  
!             RX(40)      toric radius of curvature   RDX               
!             CX(40)      toric curvature             CVX               
!             FREF(40)    Flag for REFlection                           
!                         = 1 => transmissive         GLASS or AIR      
!                         = -1 => reflective          REFL              
!             FREF0       same as FREF but for object                   
!                                                                       
!             WAVL(3)     Wavelengths                 WV                
!                                                                       
! *********** *************************************** *******           
!                                                                       
!   /COLLAT/                                                            
!                                                                       
!             CLTRA(300)  Lattice points                                
!             RADIMG      Image surf. radius of curv- IMAGE RD          
!                         ature                                         
!             CVIMG       Image surf. curvature       IMAGE CV          
!             CONIMG      Image surf. conic constant  IMAGE CC          
!             NPLANE      Number of image surfs.      IMAGE SURF        
!             LATYPE      Lattice type flag                             
!                         = 1 => single ray           SPD RAY           
!                         = 2 => polar lattice        SPD POL           
!                         = 3 => rectangular lattice  SPD RECT          
!                         = 4 => FINRAY lattice                         
!                         = 5 => one point                              
!                         = 6 => read lattice from    SPD OPFILE or     
!                                                     SPD FILE          
!                         = 7 => Rim lattice          SPD RIM           
!             ICOL(3)     color numbers               SPD               
!             NCOL        number colors                                 
!             NSURF       number surfaces                               
!             IMODE       mode flag                                     
!             IPRINT      print flag                                    
!             IPLTPR      plot flag                                     
!             IWVFLG(3)   wavelength set flag                           
!             IPR20       dev. 20 print flag                            
!             IREF        reference surface number                      
!             IJK         historical vbl. always 0                      
!             IALLPL      = 1 => multiple spot plot                     
!                         = 0 => single spot plot                       
!CCCCCCCCCCCC                                                           
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!                                                                       
!     LAST UPDATE 4/27/84 BY JOHN PARKER OF SSAI                        
!                                                                       
! 4/19/84        PUT STATEMENT 'IFLAG = 0' IN MAIN JUST BEFORE LINE     
!                LABELED 1732                                           
!                FORMATS CORRECTED; COMMAS PUT BETWEEN BLOCKS IN ACCORD 
!                WITH ANS STANDARD.                                     
!                                                                       
! 4/11/84        UPDATE--      SKEW LOGIC LEADING TO RED FIXED.         
!                                                                       
! 4/10/84        UPDATE-- SUBROUTINE PARAX ERROR IN ABERRATIONS CALC-   
!                         LATION BROUGHT IN LINE WITH USER'S MANUAL.    
!                                                                       
! 4/09/84        ABBERATION COEFS. DEFINITION TABLE NOW IN PARAX.       
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!                                                                       
!     MAIN   - CARD INPUT AND DECODING ROUTINE (APPROX. 1675 LINES)     
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      INTEGER*4 ALPHA,ARRAY1 
      INTEGER *4 OPTNA,OPTNB,OPTNC,OPTND,OPTNE,OPTNF,OPTNG,ANULI,SECTRS 
!                                                                       
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /HEAD/   LINES,IPAGE,NSYS,LAMDA,NAME(160) 
      COMMON /PUPIL/  ENPUPR,ENPUPL,EXPUPR,EXPUPL 
      COMMON /SAGPAR/ YMAX,YMIN,DELY,REFCRV,CONST,WAVENM,ISRF 
!                                                                       
!                                                                       
      DIMENSION INPUT(80),IWORD1(40),IWORD2(40),ARRAY1(80),ARRAY(20) 
      DIMENSION ALPHA(26),NUM(10),INPUT2(80),WAVE(3),DEFWV(3) 
!                                                                       
      DATA ALPHA /'A','B','C','D','E','F','G','H','I','J','K','L','M',  &
     &            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/  
      DATA NUM   /'1','2','3','4','5','6','7','8','9','0'/ 
      DATA IZERO/0/,IONE/1/,ITWO/2/,ITHREE/3/,IFOUR/4/,I81/81/ 
      DATA ZERO/0.D0/,ONE/1.D0/,TWO/2.D0/,THREE/3.D0/,FOUR/4.D0/,       &
     &     FIVE/5.D0/                                                   
      DATA IBLANK/1H /,ISLASH/'/'/,IVAL/1H,/ 
      DATA PI/3.1415926535897932384626433832795D0/ 
!                                                                       
      EQUIVALENCE (TRASH(1),XSMIN),(TRASH(2),XSMAX),(TRASH(3),YSMIN) 
      EQUIVALENCE (TRASH(4),YSMAX),(TRASH(5),FCODE),(TRASH(6),FXYJ) 
      EQUIVALENCE (TRASH(7),FXY),(TRASH(8),FXNU),(TRASH(9),FBY) 
      EQUIVALENCE (TRASH(10),FBNU) 
!                                                                       
   10 FORMAT (T1,80A1) 
   20 FORMAT (10X,80A1) 
   30 FORMAT (T2,'END OF PROGRAM'/) 
   40 FORMAT (//) 
   50 FORMAT (15X,'FIRST WORD MUST START IN COLUMNS 1 THROUGH 6') 
   60 FORMAT (15X,'NO BLANK OR COMMA FOLLOWING WORD 1'/) 
   70 FORMAT (15X,'NO BLANK OR COMMA FOLLOWING WORD 2'/) 
   80 FORMAT (15X,'UNIDENTIFIED CARD, = ',80A1/) 
   84 FORMAT (15X,'SURFACE NOT SPECIFIED FOR "GRATD" CARD '/) 
   85 FORMAT (15X,'SURFACE NOT SPECIFED FOR "CLAPD" OR "COBSD" CARDS'/) 
   86 FORMAT (15X,'SURFACE NOT SPECIFIED FOR "ASPHD" CARD '/) 
   90 FORMAT (15X,'CANNOT IDENTIFY COLORS'/) 
   95 FORMAT (15X,'CALLED "SPD GEN" BEFORE CALLING "LATT GEN"'/) 
  100 FORMAT (15X,'UNIDENTIFIED WORD 2 ON "SCY" CARD',20A1/) 
  110 FORMAT (15X,'UNIDENTIFIED WORD2 ON "UNITS" CARD',20A1/) 
  112 FORMAT (15X,'UNIDENTIFIED WORD2 ON "MODE" CARD',20A1/) 
  120 FORMAT (15X,'TRYING TO FOCUS BEFORE "SPD" HAS BEEN CALLED'/) 
  124 FORMAT (15X,'UNIDENTIFIED TYPE OF PLOT: ',20A1/) 
  126 FORMAT (15X,'UNIDENTIFIED TYPE OF PRINTOUT: ',20A1/) 
  130 FORMAT (15X,'CALLED TEL- DESIGN ROUTINES BEFORE CALLING "PXTY"'/) 
  140 FORMAT (15X,'UNKNOWN TYPE OF TEL-DESIGN ROUTINE: ',20A1/) 
  150 FORMAT (15X,'UNKNOWN RAY LATTICE TYPE: ',20A1/) 
  160 FORMAT (15X,'UNKNOWN TYPE OF RAY FAN: ',20A1/) 
  165 FORMAT (15X,'LENS CARD NOT FIRST INPUT LINE'/) 
  166 FORMAT (15X,'UNIDENTIFIED "ASTOP" QUALIFIER: ',20A1/) 
!                                                                       
!****************************************************************       
!                                                                       
!             OPEN VAX FILES                                            
!                                                                       
      OPEN (UNIT=5,  TYPE='OLD', NAME='DATAIN') 
      OPEN (UNIT=20, TYPE='NEW', FORM='FORMATTED') 
      OPEN (UNIT=34, TYPE='SCRATCH', FORM='UNFORMATTED') 
      OPEN (UNIT=40, TYPE='SCRATCH', FORM='UNFORMATTED') 
      OPEN (UNIT=50, TYPE='SCRATCH', FORM='UNFORMATTED') 
!                                                                       
!****************************************************************       
!                                                                       
      REWIND 34 
      CALL GRAPH(1,0) 
!                                                                       
      NSYS=IONE 
      ICASE=IZERO 
      ICARD=IZERO 
!             LOOP FOR READING INPUT DATA CARDS AND DETERMINING         
!             WHAT TO DO WITH THEM                                      
  170 DO 171 I=1,80 
      INPUT(I)=IBLANK 
  171 END DO 
      READ (5,10,END=430) INPUT 
      WRITE (6,20) INPUT 
      ICARD = ICARD + IONE 
      LAST=80 
  172 LET1=IZERO 
      DO 180 I=1,6 
        IF (INPUT(I).EQ.IBLANK) GO TO 180 
      LET1=I 
      GOTO 183 
  180 END DO 
!             FIRST LETTER NOT FOUND                                    
      IF (LET1.EQ.IZERO) WRITE(6,50) 
!             FIND LOGICAL CARRIAGE RETURN, IF ANY                      
  183 LSL=I81 
      DO 185 I=LET1,LAST 
      IF (INPUT(I).NE.ISLASH) GO TO 185 
      LSL=I 
      GOTO 186 
  185 END DO 
  186 IF (LSL.EQ.IONE) GO TO 170 
      LSLM1=LSL-IONE 
      DO 188 I=LET1,LSLM1 
      INPUT2(I)=INPUT(I) 
  188 END DO 
!             SEE IF INPUT IS A COMMENT STATEMENT                       
      IF (INPUT(LET1).EQ.ALPHA(3).AND.INPUT(LET1+1).EQ.IVAL) GO TO 170 
      IF (INPUT(LET1).EQ.ALPHA(3).AND.INPUT(LET1+1).EQ.IBLANK           &
     &     .AND.INPUT(LET1+2).EQ.IBLANK) GO TO 170                      
      IF (LSL.EQ.I81) GO TO 191 
      DO 190 I=LSL,80 
      INPUT2(I)=IBLANK 
  190 END DO 
!             NEED TO INITIALIZE FIRST                                  
  191 CONTINUE 
      DO 200 I=1,40 
        IWORD1(I)=IBLANK 
        IWORD2(I)=IBLANK 
  200 END DO 
      DO 210 I=1,20 
        ARRAY(I)=ZERO 
  210 END DO 
      DO 220 I=1,80 
        ARRAY1(I)=IBLANK 
  220 END DO 
!             FIND FIRST BLANK OR COMMA,                                
!             STORE POSITION AS "L1"                                    
      L1=IZERO 
      DO 230 I=LET1,LAST 
        IF ((INPUT2(I).NE.IBLANK).AND.(INPUT2(I).NE.IVAL)) GO TO 230 
        L1=I 
        GO TO 240 
  230 END DO 
!             TEST TO SEE IF A BLANK OR COMMA IS FOUND;                 
!             IF NOT, ERROR EXISTS                                      
  240 IF (L1.NE.IZERO) GO TO 250 
!             L1 = 0                                                    
      WRITE (6,60) 
      GO TO 2083 
!             L1A = END OF WORD 1                                       
!             L1B = BEGINNING OF WORD 2 OR NUMERIC FIELD                
  250 L1A=L1-IONE 
      L1B=L1+IONE 
!             FIND SECOND BLANK, STORE POSITION AS "L2"                 
      L2=IZERO 
      IF (INPUT2(L1).NE.IVAL) GO TO 255 
      L2=L1B 
      GO TO 280 
  255 DO 260 I=L1B,LAST 
        IF ((INPUT2(I).NE.IBLANK).AND.(INPUT2(I).NE.IVAL)) GO TO 260 
        L2=I 
        GO TO 270 
  260 END DO 
!             FOUND SECOND BLANK, TEST TO SEE IF VALID                  
  270 IF (L2.NE.IZERO) GO TO 280 
!             L2 = 0                                                    
      WRITE (6,70) 
      GO TO 2083 
!             L2A = END OF WORD 2                                       
!             L2B = BEGINNING OF NUMERIC FIELD                          
  280 L2A=L2-IONE 
      L2B=L2+IONE 
!             DETERMINE WORD1, FOUND IN COLS. LET1 - L1A                
  285 M=IZERO 
      DO 290 I=LET1,L1A 
        M=M+IONE 
        IWORD1(M)=INPUT2(I) 
  290 END DO 
!             IF 2 BLANKS IN A ROW, NUMERIC FIELD STARTS AFTER 2ND ONE  
!             BLANKS, IF FOUND, ARE IN COLS. L1 AND L1B                 
      IF (L2.EQ.L1B) GO TO 301 
!             WE HAVE A NON-BLANK WORD 2, GET IT                        
      J=IZERO 
      DO 300 I=L1B,L2A 
        J=J+IONE 
        IWORD2(J)=INPUT2(I) 
  300 END DO 
!             FIND END OF NUMERIC FIELD                                 
  301 L3=LAST 
      K=L2B 
      IF (INPUT2(L1).EQ.IVAL) K=L2 
      DO 304 I=K,LAST 
      IF (INPUT2(I).NE.IBLANK) GOTO 303 
      IF (L3.EQ.LAST) L3=I 
      GO TO 304 
  303 L3=LAST 
  304 END DO 
!             STORE NUMERIC FIELD IN "ARRAY1"                           
  310 IF (L2.EQ.L1B) L2B=L1B+IONE 
      J=IZERO 
      K=L2B 
      IF (INPUT2(L1).EQ.IVAL) K=L2 
!         KBEG STORES FIRST NON BLANK CHARACTER IN NUMERIC FIELD        
      KBEG=K 
      DO 315 I=K,LAST 
        IF (INPUT2(I).EQ.IBLANK) GO TO 315 
        KBEG=I 
        GO TO 316 
  315 END DO 
  316 J=0 
      LENGTH=L3-KBEG 
      DO 320 I=KBEG,L3 
        J=J+IONE 
        ARRAY1(J)=INPUT2(I) 
  320 END DO 
!             COMPUTE LENGTH OF WORD 2                                  
      LEN2=L2A-L1B+IONE 
      IF (LEN2.EQ.IZERO) LEN2=IONE 
!                                                                       
!                                                                       
!             CONVERT NUMERIC FIELDS TO ACTUAL NUMBERS,                 
!             IF CARD IS NOT AN "LI" CARD                               
      IF ( IWORD1(1).EQ.ALPHA(12) .AND. IWORD1(2).EQ.ALPHA(9))          &
     &     GO TO 330                                                    
!    5/21/84                                                            
!     FORCE FLOTIN TO BE CALLED WITH LENGTH > 0                         
      LENG = LENGTH 
      IF(LENG .LT. 1) LENG = 20 
      CALL FLOTIN (0,ARRAY1,ARRAY,LENG) 
!                                                                       
!             AT THIS POINT, HAVE WORD1, WORD2, AND NUMERIC VALUES      
!             USE ALPHA ARRAY TO DISTINGUISH COMMANDS                   
!                                                                       
!             THE 'LENS' CARD MUST BE FIRST LINE OF INPUT               
!                                                                       
!             GO TO 340 IF 'LENS' CARD                                  
!                                                                       
  330 IF ((IWORD1(1).NE.ALPHA(12) .OR. IWORD1(2).NE.ALPHA(5)            &
     &   .OR. IWORD1(3).NE.ALPHA(14) .OR. IWORD1(4).NE.ALPHA(19))       &
     &   .AND. ICARD.EQ.IONE) GO TO 333                                 
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(12) .AND. IWORD1(2).EQ.ALPHA(5)           &
     &   .AND. IWORD1(3).EQ.ALPHA(14).AND.IWORD1(4).EQ.ALPHA(19))       &
     &   GO TO 340                                                      
      GO TO 337 
!             LENS NOT FIRST CARD, STOP                                 
!                                                                       
  333 WRITE(6,165) 
      GO TO 430 
!                                                                       
!             GO TO 430 IF 'EXIT' CARD                                  
!                                                                       
  337 IF ( IWORD1(1).EQ.ALPHA(5) .AND. IWORD1(2).EQ.ALPHA(24)           &
     &     .AND. IWORD1(3).EQ.ALPHA(9).AND.IWORD1(4).EQ.ALPHA(20))      &
     &     GO TO 430                                                    
!                                                                       
!             GO TO 440 IF 'RD' OR 'CV' CARD                            
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(18) .AND. IWORD1(2).EQ.ALPHA(4)) GO TO 440 
      IF (IWORD1(1).EQ.ALPHA(3) .AND. IWORD1(2).EQ.ALPHA(22)) GO TO 440 
!                                                                       
!             GO TO 580 IF 'TH' CARD                                    
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(20) .AND. IWORD1(2).EQ.ALPHA(8)) GO TO 580 
!                                                                       
!             GO TO 590 IF 'GLASS' OR 'REFL' OR 'AIR' CARD              
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(1) .AND. IWORD1(2).EQ.ALPHA(9)            &
     &     .AND. IWORD1(3).EQ.ALPHA(18))  GO TO 590                     
      IF ( IWORD1(1).EQ.ALPHA(18) .AND. IWORD1(2).EQ.ALPHA(5)           &
     &     .AND. IWORD1(3).EQ.ALPHA(6).AND.IWORD1(4).EQ.ALPHA(12))      &
     &     GO TO 590                                                    
      IF ( IWORD1(1).EQ.ALPHA(7) .AND. IWORD1(2).EQ.ALPHA(12)           &
     &     .AND. IWORD1(3).EQ.ALPHA(1).AND.IWORD1(4).EQ.ALPHA(19)       &
     &     .AND.IWORD1(5).EQ.ALPHA(19)) GO TO 590                       
!                                                                       
!             GO TO 640 IF 'ASPH' CARD                                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(1) .AND. IWORD1(2).EQ.ALPHA(19)           &
     &     .AND. IWORD1(3).EQ.ALPHA(16).AND.IWORD1(4).EQ.ALPHA(8))      &
     &     GO TO 640                                                    
!                                                                       
!             GO TO 720 IF 'CC' CARD                                    
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(3) .AND. IWORD1(2).EQ.ALPHA(3)) GO TO 720 
!                                                                       
!             GO TO 723 IF 'CONC' CARD                                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(3) .AND. IWORD1(2).EQ.ALPHA(15) .AND.     &
     &     IWORD1(3).EQ.ALPHA(14) .AND. IWORD1(4).EQ.ALPHA(3) )         &
     &     GO TO 723                                                    
!                                                                       
!             GO TO 726 IF 'CONV' CARD                                  
      IF ( IWORD1(1).EQ.ALPHA(3) .AND. IWORD1(2).EQ.ALPHA(15) .AND.     &
     &     IWORD1(3).EQ.ALPHA(14) .AND. IWORD1(4).EQ.ALPHA(22) )        &
     &     GO TO 726                                                    
!                                                                       
!             GO TO 730 IF OR 'DEC' CARD                                
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(4) .AND. IWORD1(2).EQ.ALPHA(5)            &
     &     .AND. IWORD1(3).EQ.ALPHA(3)) GO TO 730                       
!                                                                       
!             GO TO 760 IF 'TILT' CARD                                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(20) .AND. IWORD1(2).EQ.ALPHA(9)           &
     &     .AND. IWORD1(3).EQ.ALPHA(12).AND.IWORD1(4).EQ.ALPHA(20))     &
     &     GO TO 760                                                    
!                                                                       
!             GO TO 810 IF 'RTILT' CARD                                 
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(18) .AND. IWORD1(2).EQ.ALPHA(20)          &
     &     .AND. IWORD1(3).EQ.ALPHA(9).AND.IWORD1(4).EQ.ALPHA(12)       &
     &     .AND.IWORD1(5).EQ.ALPHA(20)) GO TO 810                       
!                                                                       
!             GO TO 850 IF 'GRAT' CARD OR 'GRATD' CARD                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(7) .AND. IWORD1(2).EQ.ALPHA(18)           &
     &     .AND. IWORD1(3).EQ.ALPHA(1).AND.IWORD1(4).EQ.ALPHA(20))      &
     &     GO TO 850                                                    
!                                                                       
!             GO TO 860 IF "GORD" CARD                                  
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(7).AND.IWORD1(2).EQ.ALPHA(15)              &
     &     .AND.IWORD1(3).EQ.ALPHA(18).AND.IWORD1(4).EQ.ALPHA(4))       &
     &     GO TO 860                                                    
!                                                                       
!             GO TO 910 IF 'CLAP' OR 'CLAPD' CARD                       
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(3) .AND. IWORD1(2).EQ.ALPHA(12)           &
     &     .AND. IWORD1(3).EQ.ALPHA(1).AND.IWORD1(4).EQ.ALPHA(16))      &
     &     GO TO 910                                                    
!                                                                       
!             GO TO 950 IF 'COBS' CARD                                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(3) .AND. IWORD1(2).EQ.ALPHA(15)           &
     &     .AND. IWORD1(3).EQ.ALPHA(2).AND.IWORD1(4).EQ.ALPHA(19))      &
     &     GO TO 950                                                    
!                                                                       
!             GO TO 1060 IF 'PXTY' CARD                                 
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(16) .AND. IWORD1(2).EQ.ALPHA(24)          &
     &     .AND. IWORD1(3).EQ.ALPHA(20).AND.IWORD1(4).EQ.ALPHA(25))     &
     &     GO TO 1060                                                   
!                                                                       
!             GO TO 1110 IF 'SCY' OR 'SCX' CARD                         
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(19) .AND. IWORD1(2).EQ.ALPHA(3)           &
     &     .AND.(IWORD1(3).EQ.ALPHA(24).OR.IWORD1(3).EQ.ALPHA(25)))     &
     &     GO TO 1110                                                   
!                                                                       
!             GO TO 1240 IF 'IMAGE' CARD                                
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(9) .AND. IWORD1(2).EQ.ALPHA(13)           &
     &     .AND. IWORD1(3).EQ.ALPHA(1).AND.IWORD1(4).EQ.ALPHA(7)        &
     &     .AND.IWORD1(5).EQ.ALPHA(5)) GO TO 1240                       
!                                                                       
!             GO TO 1340 IF 'UNITS' CARD                                
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(21) .AND. IWORD1(2).EQ.ALPHA(14)          &
     &     .AND. IWORD1(3).EQ.ALPHA(9).AND.IWORD1(4).EQ.ALPHA(20)       &
     &     .AND.IWORD1(5).EQ.ALPHA(19)) GO TO 1340                      
!                                                                       
!             GO TO 1375 IF "MODE" CARD                                 
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(13).AND.IWORD1(2).EQ.ALPHA(15)             &
     &     .AND.IWORD1(3).EQ.ALPHA(4).AND.IWORD1(4).EQ.ALPHA(5))        &
     &     GO TO 1375                                                   
!                                                                       
!             GO TO 1400 IF 'FOCUS' CARD                                
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(6) .AND. IWORD1(2).EQ.ALPHA(15)           &
     &     .AND. IWORD1(3).EQ.ALPHA(3)) GO TO 1400                      
!                                                                       
!             GO TO 1420 IF 'CFL' CARD                                  
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(3) .AND. IWORD1(2).EQ.ALPHA(6)             &
     &     .AND.IWORD1(3).EQ.ALPHA(12)) GO TO 1420                      
!                                                                       
!             GO TO 1430 IF 'SAG' CARD                                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(19) .AND. IWORD1(2).EQ.ALPHA(1)           &
     &     .AND. IWORD1(3).EQ.ALPHA(7)) GO TO 1430                      
!                                                                       
!             GO TO 1440 IF 'LI' CARD                                   
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(12).AND.IWORD1(2).EQ.ALPHA(9))            &
     &     GO TO 1440                                                   
!                                                                       
!             GO TO 1490 IF 'PLOT' CARD                                 
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(16) .AND. IWORD1(2).EQ.ALPHA(12)          &
     &     .AND. IWORD1(3).EQ.ALPHA(15).AND.IWORD1(4).EQ.ALPHA(20))     &
     &     GO TO 1490                                                   
!                                                                       
!             GO TO 1492 IF "NOPLOT' CARD                               
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(14).AND.IWORD1(2).EQ.ALPHA(15)             &
     &     .AND.IWORD1(3).EQ.ALPHA(16).AND.IWORD1(4).EQ.ALPHA(12)       &
     &     .AND.IWORD1(5).EQ.ALPHA(15).AND.IWORD1(6).EQ.ALPHA(20))      &
     &     GO TO 1492                                                   
!                                                                       
!             GO TO 1494 IF "PRINT" CARD                                
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(16).AND.IWORD1(2).EQ.ALPHA(18)             &
     &     .AND.IWORD1(3).EQ.ALPHA(9).AND.IWORD1(4).EQ.ALPHA(14)        &
     &     .AND.IWORD1(5).EQ.ALPHA(20)) GO TO 1494                      
!                                                                       
!             GO TO 1496 IF "NOPRNT" CARD                               
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(14).AND.IWORD1(2).EQ.ALPHA(15)             &
     &     .AND.IWORD1(3).EQ.ALPHA(16).AND.IWORD1(4).EQ.ALPHA(18)       &
     &     .AND.IWORD1(5).EQ.ALPHA(14).AND.IWORD1(6).EQ.ALPHA(20))      &
     &     GO TO 1496                                                   
!                                                                       
!             GO TO 1500 IF 'SAY' CARD                                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(19) .AND. IWORD1(2).EQ.ALPHA(1)           &
     &     .AND. IWORD1(3).EQ.ALPHA(25)) GO TO 1500                     
!                                                                       
!             GO TO 1510 IF 'DES' CARD                                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(4) .AND. IWORD1(2).EQ.ALPHA(5)            &
     &     .AND. IWORD1(3).EQ.ALPHA(19)) GO TO 1510                     
!                                                                       
!             GO TO 1560 IF 'ASTOP' CARD                                
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(1) .AND. IWORD1(2).EQ.ALPHA(19)           &
     &     .AND. IWORD1(3).EQ.ALPHA(20).AND.IWORD1(4).EQ.ALPHA(15)      &
     &     .AND. IWORD1(5).EQ.ALPHA(16)) GOTO 1560                      
!                                                                       
!             GO TO 1566 IF 'LATT' CARD                                 
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(12).AND.IWORD1(2).EQ.ALPHA(1)              &
     &     .AND.IWORD1(3).EQ.ALPHA(20).AND.IWORD1(4).EQ.ALPHA(20))      &
     &     GO TO 1566                                                   
!                                                                       
!             GO TO 1570 IF 'SPD' CARD                                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(19) .AND. IWORD1(2).EQ.ALPHA(16)          &
     &     .AND. IWORD1(3).EQ.ALPHA(4)) GO TO 1570                      
!                                                                       
!             GO TO 1830 IF 'BIL' CARD                                  
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(2) .AND. IWORD1(2).EQ.ALPHA(9)            &
     &     .AND. IWORD1(3).EQ.ALPHA(12)) GO TO 1830                     
!                                                                       
!             GO TO 1880 IF 'LEPRT' CARD                                
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(12).AND.IWORD1(2).EQ.ALPHA(5)              &
     &     .AND.IWORD1(3).EQ.ALPHA(16).AND.IWORD1(4).EQ.ALPHA(18)       &
     &     .AND.IWORD1(5).EQ.ALPHA(20)) GO TO 1880                      
!                                                                       
!             GO TO 1900 IF 'SC' CARD                                   
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(19) .AND. IWORD1(2).EQ.ALPHA(3)           &
     &     .AND. IWORD1(3).EQ.IBLANK) GO TO 1900                        
!                                                                       
!             GO TO 1910 IF 'WV' CARD                                   
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(23) .AND. IWORD1(2).EQ.ALPHA(22)          &
     &     .AND. IWORD1(3).EQ.IBLANK) GO TO 1910                        
!                                                                       
!             GO TO 1920 IF 'INSERT' CARD                               
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(9) .AND. IWORD1(2).EQ.ALPHA(14)           &
     &     .AND. IWORD1(3).EQ.ALPHA(19)) GO TO 1920                     
!                                                                       
!             GO TO 1930 IF 'DELETE' CARD                               
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(4) .AND. IWORD1(2).EQ.ALPHA(5)            &
     &     .AND. IWORD1(3).EQ.ALPHA(12)) GO TO 1930                     
!                                                                       
!             GO TO 1940 IF 'REFS' CARD                                 
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(18).AND.IWORD1(2).EQ.ALPHA(5)             &
     &     .AND.IWORD1(3).EQ.ALPHA(6).AND.IWORD1(4).EQ.ALPHA(19))       &
     &     GO TO 1940                                                   
!                                                                       
!             GO TO 1970 IF 'STAT' CARD                                 
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(19) .AND. IWORD1(2).EQ.ALPHA(20)          &
     &     .AND. IWORD1(3).EQ.ALPHA(1)) GO TO 1970                      
!                                                                       
!             TEST FOR "PY", "PIY", OR "CAY" CARD                       
!                                                                       
      IF ( IWORD1(1).EQ.ALPHA(16) .AND. IWORD1(2).EQ.ALPHA(25))         &
     &     GO TO 2010                                                   
      IF (IWORD1(1).EQ.ALPHA(16).AND.IWORD1(2).EQ.ALPHA(9)              &
     &     .AND.IWORD1(3).EQ.ALPHA(25)) GO TO 1990                      
      IF (IWORD1(1).EQ.ALPHA(3).AND.IWORD1(2).EQ.ALPHA(1)               &
     &     .AND.IWORD1(3).EQ.ALPHA(25)) GO TO 2000                      
!                                                                       
!             GO TO 2020 IF "SLVD" CARD                                 
!                                                                       
      IF (IWORD1(1).EQ.ALPHA(19).AND.IWORD1(2).EQ.ALPHA(12)             &
     &     .AND.IWORD1(3).EQ.ALPHA(22).AND.IWORD1(4).EQ.ALPHA(4))       &
     &     GO TO 2020                                                   
!                                                                       
!             AT THIS POINT, CARD IS UNIDENTIFIED; PRINT MESSAGE        
  335 WRITE (6,80) INPUT2 
      GO TO 2083 
!                                                                       
!             HERE IF "LENS" CARD                                       
!                             BEGIN INITIALIZATION                      
!                                                                       
  340 DO 350 I=1,10 
        TRASH(I)=ZERO 
  350 END DO 
      S=ZERO 
      D=ZERO 
      RHO=ZERO 
      UFLAG=ZERO 
      OFLAG=ZERO 
      RNOBJ=ZERO 
      HYINIT=ZERO 
      HYDEL=ZERO 
      HXINIT=ZERO 
      HXDEL=ZERO 
      APSTOP=ZERO 
      SMAX=ZERO 
      RSMAX=ZERO 
      IREF=IZERO 
      FOCL=ZERO 
      DO 360 I=1,3 
        WAVL(I)=ZERO 
        OBJN(I)=ZERO 
        WAVE(I)=ZERO 
  360 END DO 
      DELIMP=ZERO 
      FPLANE=ZERO 
      FAKEA=ZERO 
      DO 390 I=1,40 
        C(I)=ZERO 
        T(I)=ZERO 
        R(I)=ZERO 
        CONIC(I)=ZERO 
        DO 370 J=1,3 
          FN(I,J)=ZERO 
          ORDN(I,J)=ZERO 
  370   CONTINUE 
        FMASK(I)=ZERO 
        FAKEC(I)=ZERO 
        FAKEB(I)=ZERO 
        XDISP(I)=ZERO 
        YDISP(I)=ZERO 
        TILTX(I)=ZERO 
        TILTY(I)=ZERO 
        TILTZ(I)=ZERO 
        SIDE(I)=ZERO 
        RDSPAC(I)=ZERO 
        Y0(I)=ZERO 
        SXY(I)=ZERO 
        SXNU(I)=ZERO 
        DO 380 J=1,4 
          COEF(I,J)=ZERO 
  380   CONTINUE 
        XMN(I)=ZERO 
        XMX(I)=ZERO 
        YMN(I)=ZERO 
        YMX(I)=ZERO 
        RX(I)=ZERO 
        CX(I)=ZERO 
        FREF(I)=ZERO 
  390 END DO 
!                                                                       
      FREF0 = ZERO 
!                                                                       
      ENPUPR=ZERO 
      ENPUPL=ZERO 
      EXPUPR=ZERO 
      EXPUPL=ZERO 
!                                                                       
      ISRF   = IZERO 
      IFLAG  = IZERO 
      YMAX   = ZERO 
      YMIN   = ZERO 
      DELY   = ZERO 
      REFCRV = ZERO 
      CONST  = ZERO 
      WAVENM = ZERO 
!                                                                       
      DO 400 I=1,300 
        CLTRA(I)=ZERO 
  400 END DO 
      RADIMG=ZERO 
      CVIMG=ZERO 
      CONIMG=ZERO 
      NPLANE=IONE 
      LATYPE=IZERO 
      DO 410 I=1,3 
        ICOL(I)=IZERO 
  410 END DO 
      NCOL=IZERO 
      NSURF=IZERO 
      IMODE=IZERO 
      IPRINT=IZERO 
      IPLTPR=IZERO 
!                                                                       
      LINES=ITWO 
      IPAGE=IZERO 
      LAMDA=IZERO 
      DO 420 I=1,160 
        NAME(I)=IBLANK 
  420 END DO 
!                                                                       
      IF (ICASE.NE.IZERO) NSYS=NSYS+IONE 
      ICASE=IONE 
      NFOC=IZERO 
      LFOC=IZERO 
      IFOC=IZERO 
      ISKEW=IZERO 
      INDEX=IONE 
      IPEN=IZERO 
      ISSTRA=IZERO 
      IAP=IZERO 
      LATFLG=IZERO 
      IPR20=IZERO 
      IALLPL=IZERO 
!                                                                       
!                                                                       
!                                                                       
!             ISURF REPRESENTS THE CURRENT SURFACE NUMBER               
!             FINAL NUMBER OF SURFACES (NSURF) IS DETERMINED            
!             BY USING THE 'IMAGE' CARD                                 
!                                                                       
      ISURF=IZERO 
!             END OF INITIALIZATION                                     
!                                                                       
      GO TO 2083 
!                                                                       
!             HERE IF "EXIT" OR END-OF-FILE                             
  430 WRITE (6,30) 
      REWIND 34 
      CALL GRAPH(3,0) 
      GO TO 2090 
!                                                                       
!             HERE IF "RD" OR "CV" CARD                                 
!             DETERMINE SURFACE NUMBER TO USE FOR THIS DATA CARD        
  440 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      IF (JSURF.GT.IZERO) GO TO 450 
      WRITE(6,445) 
  445 FORMAT(T6,' OBJECT CANNOT HAVE RADIUS OR CURVATURE SPECIFIED ') 
      GO TO 2083 
!             DETERMINE IF TORIC (XZ) DATA OR MERIDIONAL (YZ) DATA      
  450   IF ( IWORD1(3).EQ.ALPHA(24)) GO TO 460 
!                                                                       
!             FOUND MERIDIONAL DATA                                     
      IF (IWORD1(1).EQ.ALPHA(3)) C(JSURF)=ARRAY(1) 
      IF (IWORD1(1).EQ.ALPHA(18)) R(JSURF)=ARRAY(1) 
      IF (R(JSURF).NE.ZERO) C(JSURF)=ONE/R(JSURF) 
      IF (R(JSURF).EQ.ZERO.AND.C(JSURF).NE.ZERO) R(JSURF)=ONE/C(JSURF) 
      GO TO 2083 
!             FOUND TORIC DATA                                          
  460 IF (IWORD1(1).EQ.ALPHA(3)) CX(JSURF)=ARRAY(1) 
      IF (IWORD1(1).EQ.ALPHA(18)) RX(JSURF)=ARRAY(1) 
      IF (RX(JSURF).NE.ZERO) CX(JSURF)=ONE/RX(JSURF) 
      IF ( RX(JSURF).EQ.ZERO .AND. CX(JSURF).NE.ZERO )                  &
     &     RX(JSURF) = ONE/CX(JSURF)                                    
      Y0(JSURF)=ONE 
      GO TO 2083 
!                                                                       
!             HERE IF "TH" CARD                                         
!                                                                       
  580 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      IF (JSURF.EQ.IZERO) GO TO 585 
      T(JSURF)=ARRAY(1) 
      GO TO 2083 
  585 S=ARRAY(1) 
      GO TO 2083 
!                                                                       
!             HERE IF "GLASS","AIR", OR "REFL" CARD                     
!                                                                       
  590 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      IF (IWORD1(1).EQ.ALPHA(1).AND.IWORD1(2).EQ.ALPHA(9))  GO TO 600 
      IF (IWORD1(1).EQ.ALPHA(18).AND.IWORD1(2).EQ.ALPHA(5)) GO TO 610 
      IF (IWORD1(1).EQ.ALPHA(7).AND.IWORD1(2).EQ.ALPHA(12)) GO TO 620 
      GO TO 2083 
!             FOUND 'AIR' CARD                                          
  600 IF (JSURF.EQ.IZERO) GO TO 605 
      FN(JSURF,1) = ONE 
      FN(JSURF,2) = ONE 
      FN(JSURF,3) = ONE 
      FREF(JSURF) = ONE 
      GO TO 630 
  605 OBJN(1) = ONE 
      OBJN(2) = ONE 
      OBJN(3) = ONE 
      FREF0 = ONE 
      GO TO 630 
!             FOUND 'REFL' CARD                                         
  610 IF (JSURF.EQ.IZERO) GO TO 615 
      FREF(JSURF) = -ONE 
      FN(JSURF,1) = ONE 
      FN(JSURF,2) = ONE 
      FN(JSURF,3) = ONE 
      KK=JSURF 
      GO TO 630 
  615 FREF0 = -ONE 
      OBJN(1) = ONE 
      OBJN(2) = ONE 
      OBJN(3) = ONE 
      GO TO 630 
!             FOUND 'GLASS' CARD                                        
  620 IF (JSURF.EQ.IZERO) GO TO 626 
      FN(JSURF,1) = ARRAY(1) 
      FN(JSURF,2) = ARRAY(2) 
      FN(JSURF,3) = ARRAY(3) 
      FREF(JSURF) = ONE 
      GO TO 630 
  626 OBJN(1) = ARRAY(1) 
      OBJN(2) = ARRAY(2) 
      OBJN(3) = ARRAY(3) 
      FREF0 = ONE 
!                                                                       
  630 IF (JSURF.EQ.ISURF) ISURF=ISURF+IONE 
      GO TO 2083 
!                                                                       
!             HERE IF "ASPH" CARD                                       
  640 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      IF (IWORD1(5).EQ.IBLANK) GO TO 700 
      GO TO 335 
!             STORE ALL COEFFICIENTS                                    
  700 DO 710 I=1,4 
        COEF(JSURF,I)=ARRAY(I) 
  710 END DO 
      GO TO 2083 
!                                                                       
!             HERE IF "CC" CARD                                         
  720 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      CONIC(JSURF)=ARRAY(1) 
      GO TO 2083 
!                                                                       
!             HERE IF "CONC" CARD                                       
  723 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      SIDE(JSURF) = ONE 
      GO TO 2083 
!                                                                       
!             HERE IF "CONV" CARD                                       
  726 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      SIDE(JSURF) = -ONE 
      GO TO 2083 
!             HERE IF "DEC" CARD                                        
  730 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
!             FOUND BOTH X AND Y DISPLACEMENTS                          
  750 YDISP(JSURF)=ARRAY(1) 
      XDISP(JSURF)=ARRAY(2) 
      GO TO 2083 
!                                                                       
!             HERE IF "TILT" CARD                                       
  760 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
!             FOUND X,Y, AND Z TILTS                                    
  800 TILTX(JSURF)=ARRAY(1) 
      TILTY(JSURF)=ARRAY(2) 
      TILTZ(JSURF)=ARRAY(3) 
      GO TO 2083 
!                                                                       
!             HERE IF "RTILT" CARD                                      
  810 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      DO 815 I=1,LEN2 
        IF ( IWORD2(I).EQ.IBLANK ) GO TO 840 
        IF ( IWORD2(I).EQ.ALPHA(4) .AND. IWORD2(I+1).EQ.ALPHA(5)        &
     &     .AND.IWORD2(I+2).EQ.ALPHA(3).AND.IWORD2(I+3).EQ.ALPHA(19))   &
     &     GO TO 830                                                    
        IF ( IWORD2(I).EQ.ALPHA(20) .AND. IWORD2(I+1).EQ.ALPHA(9)       &
     &     .AND.IWORD2(I+2).EQ.ALPHA(12).AND.IWORD2(I+3).EQ.ALPHA(20)   &
     &     .AND.IWORD2(I+4).EQ.ALPHA(19) ) GO TO 820                    
  815 END DO 
      GO TO 840 
!             RESTORE TILTS                                             
  820 FAKEB(JSURF)=TWO 
      GO TO 2083 
!             RESTORE DISPLACEMENTS                                     
  830 FAKEB(JSURF)=ONE 
      GO TO 2083 
!             RESTORE DISPLACEMENTS AND TILTS                           
  840 FAKEB(JSURF)=THREE 
      GO TO 2083 
!                                                                       
!             HERE IF "GRAT" OR "GRATD" CARD                            
  850 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      TEMP=ONE 
        IF ( IWORD1(5).EQ.ALPHA(24) ) GO TO 880 
        IF ( IWORD1(5).EQ.ALPHA(25) ) GO TO 890 
        IF ( IWORD1(5).EQ.ALPHA(4) ) GO TO 870 
      GO TO 890 
!             HERE IF "GORD" CARD                                       
  860 CALL SURFNO(IWORD2,LEN2,ISURF,JSURF) 
      CALL CNTCOM (NCOM,INPUT2(KBEG),LENGTH) 
      IF (NCOM.GT.1) GO TO 865 
        ORDN(JSURF,1) = ARRAY(1) 
        ORDN(JSURF,2) = ARRAY(1) 
        ORDN(JSURF,3) = ARRAY(1) 
      GO TO 2083 
  865   ORDN(JSURF,1) = ARRAY(1) 
        ORDN(JSURF,2) = ARRAY(2) 
        ORDN(JSURF,3) = ARRAY(3) 
      GO TO 2083 
  870 CALL SURFNO(IWORD2,LEN2,ISURF,JSURF) 
      IF (JSURF.GT.ZERO) GO TO 875 
!             NO SURFACE SPECIFIED                                      
      WRITE(6,84) 
      GO TO 2083 
  875 RDSPAC(JSURF) = ZERO 
      ORDN(JSURF,1) = ZERO 
      ORDN(JSURF,2) = ZERO 
      ORDN(JSURF,3) = ZERO 
      GO TO 2083 
!                                                                       
!             RULINGS ARE X-RULINGS                                     
  880 TEMP=-ONE 
      GO TO 900 
!                                                                       
!             RULINGS ARE Y-RULINGS (DEFAULT)                           
  890 TEMP=ONE 
!                                                                       
!             STORE SPACING                                             
  900 RDSPAC(JSURF)=TEMP*DABS(ARRAY(1)) 
      GO TO 2083 
!                                                                       
!             HERE IF "CLAP" CARD                                       
  910 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      IF (IWORD1(5).EQ.ALPHA(4)) GO TO 990 
      DO 920 I=1,LEN2 
        IF ( IWORD2(I).EQ.IBLANK ) GO TO 930 
        IF ( IWORD2(I).EQ.ALPHA(18) .AND. IWORD2(I+1).EQ.ALPHA(5)       &
     &       .AND.IWORD2(I+2).EQ.ALPHA(3).AND.IWORD2(I+3).EQ.ALPHA(20)) &
     &       GO TO 970                                                  
        IF ( IWORD2(I).EQ.ALPHA(5).AND.IWORD2(I+1).EQ.ALPHA(12)         &
     &       .AND.IWORD2(I+2).EQ.ALPHA(9).AND.IWORD2(I+3).EQ.ALPHA(16)) &
     &       GO TO 983                                                  
  920 END DO 
      GO TO 930 
!             MASK IS A CIRCULAR CLEAR APERTURE                         
  930 FMASK(JSURF)=DABS(ARRAY(1)) 
      YMN(JSURF)=ARRAY(2) 
      XMN(JSURF)=ARRAY(3) 
      FAKEC(JSURF)=ZERO 
      GO TO 2083 
!             MASK IS A CIRCULAR OBSCURATION                            
  940 FMASK(JSURF)=-DABS(ARRAY(1)) 
      YMN(JSURF)=ARRAY(2) 
      XMN(JSURF)=ARRAY(3) 
      FAKEC(JSURF)=ZERO 
      GO TO 2083 
!                                                                       
!             HERE IF "COBS" CARD                                       
  950 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      IF (IWORD1(5).EQ.ALPHA(4)) GO TO 990 
      DO 960 I=1,LEN2 
        IF ( IWORD2(I).EQ.IBLANK ) GO TO 940 
        IF ( IWORD2(I).EQ.ALPHA(18) .AND. IWORD2(I+1).EQ.ALPHA(5)       &
     &       .AND.IWORD2(I+2).EQ.ALPHA(3).AND.IWORD2(I+3).EQ.ALPHA(20)) &
     &       GO TO 980                                                  
        IF ( IWORD2(I).EQ.ALPHA(5).AND.IWORD2(I+1).EQ.ALPHA(12)         &
     &       .AND.IWORD2(I+2).EQ.ALPHA(9).AND.IWORD2(I+3).EQ.ALPHA(16)) &
     &       GO TO 986                                                  
  960 END DO 
      GO TO 940 
!             MASK IS A RECTANGULAR CLEAR APERTURE                      
  970 YMN(JSURF)=ARRAY(1) 
      YMX(JSURF)=ARRAY(2) 
      XMN(JSURF)=ARRAY(3) 
      XMX(JSURF)=ARRAY(4) 
      FAKEC(JSURF)=ONE 
      FMASK(JSURF)=ONE 
      GO TO 2083 
!             MASK IS A RECTANGULAR OBSCURATION                         
  980 YMN(JSURF)=ARRAY(1) 
      YMX(JSURF)=ARRAY(2) 
      XMN(JSURF)=ARRAY(3) 
      XMX(JSURF)=ARRAY(4) 
      FAKEC(JSURF)=ONE 
      FMASK(JSURF)=-ONE 
      GO TO 2083 
!                                                                       
!              MAP IS ELLIPTICAL CLEAR APERTURE                         
!                                                                       
  983   YMX(JSURF)=ARRAY(1) 
        XMX(JSURF)=ARRAY(2) 
        YMN(JSURF)=ARRAY(3) 
        XMN(JSURF)=ARRAY(4) 
        FAKEC(JSURF)=-ONE 
        FMASK(JSURF)=ONE 
        GO TO 2083 
!                                                                       
!              MAP IS ELLIPTICAL OBSCURATION                            
!                                                                       
  986   YMX(JSURF)=ARRAY(1) 
        XMX(JSURF)=ARRAY(2) 
        YMN(JSURF)=ARRAY(3) 
        XMN(JSURF)=ARRAY(4) 
        FAKEC(JSURF)=-ONE 
        FMASK(JSURF)=-ONE 
        GO TO 2083 
!                                                                       
!             HERE IF "CLAPD" CARD OR "COBSD" CARD                      
!                                                                       
  990 CALL SURFNO(IWORD2,LEN2,ISURF,JSURF) 
      IF (JSURF.GT.ZERO) GO TO 1000 
!             SURFACE NOT SPECIFIED                                     
      WRITE(6,85) 
      GO TO 2083 
 1000 FMASK(JSURF)=ZERO 
      FAKEC(JSURF)=ZERO 
      XMN(JSURF)=ZERO 
      XMX(JSURF)=ZERO 
      YMN(JSURF)=ZERO 
      YMX(JSURF)=ZERO 
      GO TO 2083 
!                                                                       
!             HERE IF "PXTY" CARD                                       
!             DETERMINE MAXIMUM SURFACE NUMBER                          
!             THIS IS DETERMINED BY FIRST SURFACE WITH A ZERO INDEX     
 1060 DO 1070 I=1,40 
        IF ( FN(I,1).EQ.ZERO ) GO TO 1080 
 1070 END DO 
 1080 ISURF=I 
      NSURF=ISURF-IONE 
!             CHANGE LINEAR HEIGHTS TO ANGLES                           
 1084 IF (IFLAG.EQ.0) GO TO 1086 
      HXDEL=-DATAN(HXDEL/(S-D))*(180.D0/PI) 
      HYDEL=-DATAN(HYDEL/(S-D))*(180.D0/PI) 
      HXINIT=-DATAN(HXINIT/(S-D))*(180.D0/PI) 
      HYINIT=-DATAN(HYINIT/(S-D))*(180.D0/PI) 
 1086 IF ((IPRINT/2)*2.NE.IPRINT) GO TO 1111 
!             DETERMINE COLORS FROM CARD                                
      CALL CNTCOM (NCOM,INPUT2(KBEG),LENGTH) 
      IF (NCOM.GT.IZERO) GO TO 1090 
!             CAN'T IDENTIFY COLORS                                     
      WRITE (6,90) 
      GO TO 2083 
 1090 IF (NCOM.GT.ITHREE) NCOM=ITHREE 
      NCOL=NCOM 
      DO 1100 I=1,NCOL 
        ICOL(I)=ARRAY(I) 
        IF (ICOL(I).LT.1) ICOL(I)=1 
        IF (ICOL(I).GT.3) ICOL(I)=3 
 1100 END DO 
!             CONVERT WAVELENGTH FROM NANOMETERS TO SYSTEM UNITS        
      DO 1105 I=1,NCOL 
      IF (UFLAG.EQ.ONE) WAVL(I)=WAVE(ICOL(I))*1.D-6 
      IF (UFLAG.EQ.TWO) WAVL(I)=WAVE(ICOL(I))*1.D-7 
      IF (UFLAG.EQ.THREE) WAVL(I)=(WAVE(ICOL(I))*1.D-6)/25.4D0 
 1105 END DO 
 1111 WRITE (6,40) 
      CALL PARAX 
      WRITE (6,40) 
      GO TO 2083 
!                                                                       
!             HERE IF "SCY" CARD                                        
 1110 IF (IWORD2(1).EQ.IBLANK)    GO TO 1150 
      IF (IWORD2(1).EQ.ALPHA(6).AND.IWORD2(2).EQ.ALPHA(1)               &
     &     .AND.IWORD2(3).EQ.ALPHA(14).AND.IWORD2(4).EQ.ALPHA(7))       &
     &     GO TO 1140                                                   
!             CAN'T IDENTIFY WORD2 ON "SCY" CARD                        
      WRITE(6,100) IWORD2 
      GO TO 2083 
!             FOUND "FANG" AS 2ND WORD                                  
 1140 RNOBJ=ARRAY(1) 
      HINIT=ARRAY(2) 
      HDEL=ARRAY(3) 
      IFLAG = IZERO 
      GO TO 1160 
!             FOUND BLANK AS 2ND WORD                                   
 1150 RNOBJ=ARRAY(1) 
      IFLAG = IONE 
      HINIT=ARRAY(2) 
      HDEL=ARRAY(3) 
 1160 IF (IWORD1(3).EQ.ALPHA(24)) GO TO 1180 
      IF (IWORD1(3).EQ.ALPHA(25)) GO TO 1170 
!              CAN'T IDENTIFY WORD2 ON "SCY" OR "SCX" CARD              
      WRITE (6,100) IWORD2 
      GO TO 2083 
!             FOUND Y OBJECT                                            
 1170 HYINIT=HINIT 
      HYDEL=HDEL 
      GO TO 2083 
!             FOUND X OJECT                                             
 1180 HXINIT=HINIT 
      HXDEL=HDEL 
      GO TO 2083 
!                                                                       
!             HERE IF "IMAGE" CARD                                      
 1240 IF (IWORD2(1).EQ.IBLANK) GO TO 1310 
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(4))             &
     &     GO TO 1250                                                   
      IF (IWORD2(1).EQ.ALPHA(6).AND.IWORD2(2).EQ.ALPHA(9)               &
     &     .AND.IWORD2(3).EQ.ALPHA(18).AND.IWORD2(4).EQ.ALPHA(19)       &
     &     .AND.IWORD2(5).EQ.ALPHA(20)) GO TO 1260                      
      IF (IWORD2(1).EQ.ALPHA(19).AND.IWORD2(2).EQ.ALPHA(21)             &
     &     .AND.IWORD2(3).EQ.ALPHA(18).AND.IWORD2(4).EQ.ALPHA(6))       &
     &     GO TO 1270                                                   
      IF (IWORD2(1).EQ.ALPHA(19).AND.IWORD2(2).EQ.ALPHA(5)              &
     &     .AND.IWORD2(3).EQ.ALPHA(16)) GO TO 1280                      
      IF (IWORD2(1).EQ.ALPHA(3).AND.IWORD2(2).EQ.ALPHA(3))   GO TO 1290 
      IF (IWORD2(1).EQ.ALPHA(3).AND.IWORD2(2).EQ.ALPHA(22))  GO TO 1300 
      GO TO 335 
!                                                                       
!             FOUND RADIUS                                              
 1250 RADIMG=ARRAY(1) 
      GO TO 1320 
!             FOUND FIRST IMAGE SURFACE NUMBER                          
 1260 FPLANE=ARRAY(1) 
      GO TO 1320 
!             FOUND NUMBER OF IMAGE SURFACES                            
 1270 NPLANE=ARRAY(1) 
      GO TO 1320 
!             FOUND IMAGE SURFACE SEPARATION                            
 1280 DELIMP=ARRAY(1) 
      GO TO 1320 
!             FOUND IMAGE SURFACE CONIC CONSTANT                        
 1290 CONIMG=ARRAY(1) 
      GO TO 1320 
!             FOUND IMAGE CURVATURE                                     
 1300 CVIMG=ARRAY(1) 
      GO TO 1320 
!             FOUND BLANK, EVERYTHING IS ON ONE CARD                    
 1310 NPLANE=ARRAY(1) 
      RADIMG=ARRAY(2) 
      CVIMG=ARRAY(3) 
      CONIMG=ARRAY(4) 
      DELIMP=ARRAY(5) 
      FPLANE=ARRAY(6) 
!                                                                       
 1320 IF (RADIMG.NE.ZERO) CVIMG=ONE/RADIMG 
      IF (RADIMG.EQ.ZERO.AND.CVIMG.NE.ZERO) RADIMG=ONE/CVIMG 
      IF (FPLANE.NE.ZERO) GO TO 2083 
      M=IONE+NPLANE/ITWO 
      M=M-IONE 
      FPLANE=-M 
      GO TO 2083 
!                                                                       
!             HERE IF "UNITS" CARD                                      
 1340 IF (IWORD2(1).EQ.ALPHA(13).AND.IWORD2(2).EQ.ALPHA(13)) GO TO 1350 
      IF (IWORD2(1).EQ.ALPHA(3).AND.IWORD2(2).EQ.ALPHA(13))  GO TO 1360 
      IF (IWORD2(1).EQ.ALPHA(9).AND.IWORD2(2).EQ.ALPHA(14)              &
     &     .AND.IWORD2(3).EQ.ALPHA(3).AND.IWORD2(4).EQ.ALPHA(8)         &
     &     .AND.IWORD2(5).EQ.ALPHA(5).AND.IWORD2(6).EQ.ALPHA(19))       &
     &     GO TO 1370                                                   
!             AT THIS POINT, UNKNOWN UNITS                              
      WRITE (6,110) IWORD2 
      GO TO 2083 
!                                                                       
!             UNITS ARE MILLIMETERS                                     
 1350 IF (OFLAG.EQ.ZERO) UFLAG=ONE 
      IF (OFLAG.NE.ZERO) OFLAG=ONE 
      GO TO 2083 
!             UNITS ARE CENTIMETERS                                     
 1360 IF (OFLAG.EQ.ZERO) UFLAG=TWO 
      IF (OFLAG.NE.ZERO) OFLAG=TWO 
      GO TO 2083 
!             UNITS ARE INCHES                                          
 1370 IF (OFLAG.EQ.ZERO) UFLAG=THREE 
      IF (OFLAG.NE.ZERO) OFLAG=THREE 
      GO TO 2083 
!                                                                       
!             HERE IF "MODE" CARD                                       
!                                                                       
 1375 IF (IWORD2(1).EQ.ALPHA(1).AND.IWORD2(2).EQ.ALPHA(6)               &
     &     .AND.IWORD1(3).EQ.ALPHA(15).AND.IWORD1(4).EQ.ALPHA(3)        &
     &     .AND.IWORD1(5).EQ.ALPHA(1).AND.IWORD1(6).EQ.ALPHA(12))       &
     &     GO TO 1380                                                   
      IF (IWORD2(1).EQ.ALPHA(6).AND.IWORD2(2).EQ.ALPHA(15)              &
     &     .AND.IWORD2(3).EQ.ALPHA(3).AND.IWORD2(4).EQ.ALPHA(1)         &
     &     .AND.IWORD2(5).EQ.ALPHA(12)) GO TO 1395                      
      IF (IWORD1(1).EQ.ALPHA(1).AND.IWORD1(2).EQ.ALPHA(14)              &
     &     .AND.IWORD1(3).EQ.ALPHA(7).AND.IWORD1(4).EQ.ALPHA(12)        &
     &     .AND.IWORD1(5).EQ.ALPHA(5).AND.IWORD1(6).EQ.ALPHA(19))       &
     &     GO TO 1390                                                   
!             AT THIS POINT, UNKNOWN MODE                               
      WRITE(6,112) IWORD2 
      GO TO 2083 
!      OUTPUT UNITS ARE TANGENTS OF ANGLES (W/ RESPECT TO AXIS)         
 1380 IF (OFLAG.EQ.ZERO) OFLAG=UFLAG 
      UFLAG=FOUR 
      GO TO 2083 
!      OUTPUT UNITS ARE ANGLES OF INCIDENCE(W/ RESPECT TO NORMAL)       
 1390 IF (OFLAG.EQ.ZERO) OFLAG=UFLAG 
      UFLAG=FIVE 
      GO TO 2083 
!      RETURN TO FOCAL MODE                                             
 1395 IF (OFLAG.NE.ZERO) UFLAG=OFLAG 
      OFLAG=ZERO 
      GO TO 2083 
!                                                                       
!             HERE IF "FOCUS" CARD                                      
!             FIRST TEST TO SEE IF SKEW ALREADY CALLED                  
 1400 IF (ISKEW.NE.IZERO) GO TO 1410 
!             HAVEN'T CALLED SKEW YET, BEFORE FOCUS                     
      WRITE (6,120) 
      GO TO 2083 
 1410 OBJNO=ARRAY(1) 
      COLNO=ARRAY(2) 
      IF (OBJNO.EQ.ZERO) OBJNO=ONE 
      IF (COLNO.EQ.ZERO) COLNO=ONE 
      FCODE=THREE*(OBJNO-ONE)+COLNO 
      NFOC=OBJNO 
      LFOC=COLNO 
      IF (IWORD2(1).EQ.IBLANK.OR.IWORD2(1).EQ.ALPHA(19)) IFOC=IZERO 
      IF (IWORD2(1).EQ.ALPHA(24)) IFOC=-IONE 
      IF (IWORD2(1).EQ.ALPHA(25)) IFOC=IONE 
      GO TO 2083 
!                                                                       
!             HERE IF "CFL" CARD (FOCAL LENGTH)                         
 1420 FOCL=ARRAY(1) 
      GO TO 2083 
!                                                                       
!             HERE IF "SAG" CARD                                        
 1430 IF ( IWORD2(1).EQ.IBLANK ) GO TO 1437 
      IF ( IWORD2(1).EQ.ALPHA(19).AND.IWORD2(2).EQ.ALPHA(21)            &
     &     .AND.IWORD2(3).EQ.ALPHA(18).AND.IWORD2(4).EQ.ALPHA(6))       &
     &     GO TO 1431                                                   
      IF ( IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(1) .AND.       &
     &     IWORD2(3).EQ.ALPHA(14).AND.IWORD2(4).EQ.ALPHA(7)             &
     &     .AND.IWORD2(5).EQ.ALPHA(5)) GO TO 1432                       
      IF ( IWORD2(1).EQ.ALPHA(3).AND.IWORD2(2).EQ.ALPHA(15)             &
     &     .AND.IWORD2(3).EQ.ALPHA(14).AND.IWORD2(4).EQ.ALPHA(19)       &
     &     .AND.IWORD2(5).EQ.ALPHA(20)) GO TO 1433                      
      IF ( IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(4) ) GO TO 1434 
      IF ( IWORD2(1).EQ.ALPHA(3).AND.IWORD2(2).EQ.ALPHA(22) ) GO TO 1435 
      IF ( IWORD2(1).EQ.ALPHA(23).AND.IWORD2(2).EQ.ALPHA(22)) GO TO 1436 
      GO TO 335 
!             FOUND "SURF" AS WORD 2                                    
 1431 ISRF=ARRAY(1) 
      GO TO 2083 
!             FOUND "RANGE" AS WORD 2                                   
 1432 YMIN=ARRAY(1) 
      YMAX=ARRAY(2) 
      DELY=ARRAY(3) 
      GO TO 2083 
!             FOUND "CONST" AS WORD 2                                   
 1433 CONST=ARRAY(1) 
      GO TO 2083 
!             FOUND "RD" AS WORD 2                                      
 1434 REFRAD = ARRAY(1) 
      REFCRV = ZERO 
      IF ( REFRAD.NE.ZERO ) REFCRV = ONE/REFRAD 
      GO TO 2083 
!             FOUND "CV" AS WORD 2                                      
 1435 REFCRV = ARRAY(1) 
      GO TO 2083 
!             FOUND "WV" AS WORD 2                                      
 1436 WAVENM = ARRAY(1) 
      GO TO 2083 
!             FOUND BLANK, EVERYTHING IS ON ONE LINE                    
 1437 ISRF=ARRAY(1) 
      YMIN=ARRAY(2) 
      YMAX=ARRAY(3) 
      DELY=ARRAY(4) 
      REFCRV=ARRAY(5) 
      CONST=ARRAY(6) 
      WAVENM=ARRAY(7) 
      GO TO 2083 
!                                                                       
!             HERE IF "LI" OR "LIC" CARD                                
 1440 IF ( IWORD1(3).EQ.IBLANK )   GO TO 1470 
      IF ( IWORD1(3).EQ.ALPHA(3)) GO TO 1450 
      GO TO 335 
!             FOUND PLOT LABEL                                          
 1450 K=L2B 
      IF (INPUT2(L1).EQ.IVAL) K=L2 
      DO 1460 I=1,25 
        NAME(80+I)=INPUT2(K+I-IONE) 
 1460 END DO 
      GO TO 2083 
!             HERE IF "LI" CARD                                         
!             FOUND LABEL USED FOR PRINTED OUTPUT                       
 1470 J=IZERO 
      K=L2B 
      IF (INPUT2(L1).EQ.IVAL) K=L2 
      DO 1480 I=K,80 
        J=J+IONE 
        NAME(J)=INPUT2(I) 
 1480 END DO 
      GO TO 2083 
!                                                                       
!                                                                       
!             HERE IF "PLOT" CARD                                       
 1490 IF (IWORD2(1).EQ.ALPHA(19).AND.IWORD2(2).EQ.ALPHA(16)             &
     &     .AND.IWORD2(3).EQ.ALPHA(4)) GO TO 1950                       
      IF (IWORD2(1).EQ.ALPHA(19).AND.IWORD2(2).EQ.ALPHA(3)              &
     &     .AND.IWORD2(3).EQ.ALPHA(1).AND.IWORD2(4).EQ.ALPHA(12)        &
     &     .AND.IWORD2(5).EQ.ALPHA(5)) GO TO 1498                       
      IF (IWORD2(1).EQ.ALPHA(7).AND.IWORD2(2).EQ.ALPHA(15)              &
     &     .AND.IWORD2(3).EQ.ALPHA(20).AND.IWORD2(4).EQ.ALPHA(6))       &
     &     GO TO 1760                                                   
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(5)              &
     &     .AND.IWORD2(3).EQ.ALPHA(4)) GO TO 1810                       
      IF (IWORD2(1).EQ.ALPHA(15).AND.IWORD2(2).EQ.ALPHA(16)             &
     &     .AND.IWORD2(3).EQ.ALPHA(19).AND.IWORD2(4).EQ.ALPHA(16))      &
     &     GO TO 1945                                                   
!             AT THIS POINT, UNIDENTIFIED PLOT OPTION                   
      WRITE(6,124) IWORD2 
      GO TO 2083 
!                                                                       
!             HERE IF "NOPLOT" CARD                                     
!                                                                       
 1492 IF (IWORD2(1).EQ.ALPHA(19).AND.IWORD2(2).EQ.ALPHA(16)             &
     &     .AND.IWORD2(3).EQ.ALPHA(4)) GO TO 1960                       
      IF (IWORD2(1).EQ.ALPHA(15).AND.IWORD2(2).EQ.ALPHA(16)             &
     &     .AND.IWORD2(3).EQ.ALPHA(19).AND.IWORD2(4).EQ.ALPHA(16))      &
     &     GO TO 1960                                                   
      IF (IWORD2(1).EQ.ALPHA(7).AND.IWORD2(2).EQ.ALPHA(15)              &
     &     .AND.IWORD2(3).EQ.ALPHA(20).AND.IWORD2(4).EQ.ALPHA(6))       &
     &     GO TO 1770                                                   
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(5)              &
     &     .AND.IWORD2(3).EQ.ALPHA(4)) GO TO 1820                       
      IF (IWORD2(1).EQ.ALPHA(15).AND.IWORD2(2).EQ.ALPHA(16)             &
     &     .AND.IWORD2(3).EQ.ALPHA(19).AND.IWORD2(4).EQ.ALPHA(16))      &
     &     GO TO 1965                                                   
      GO TO 2083 
!                                                                       
!             HERE IF "PRINT" CARD                                      
!                                                                       
 1494 IF (IWORD2(1).EQ.ALPHA(7).AND.IWORD2(2).EQ.ALPHA(15)              &
     &     .AND.IWORD2(3).EQ.ALPHA(20).AND.IWORD2(4).EQ.ALPHA(6))       &
     &     GO TO 1740                                                   
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(5)              &
     &     .AND.IWORD2(3).EQ.ALPHA(4)) GO TO 1790                       
      IF (IWORD2(1).EQ.ALPHA(16).AND.IWORD2(2).EQ.ALPHA(18)             &
     &     .AND.IWORD2(3).EQ.ALPHA(24).AND.IWORD2(4).EQ.ALPHA(25)       &
     &     .AND.IWORD2(5).EQ.ALPHA(26)) GO TO 1850                      
      IF (IWORD2(1).EQ.ALPHA(19).AND.IWORD2(2).EQ.ALPHA(20)             &
     &     .AND.IWORD2(3).EQ.ALPHA(1).AND.IWORD2(4).EQ.ALPHA(20)        &
     &     .AND.IWORD2(5).EQ.ALPHA(16).AND.IWORD2(6).EQ.ALPHA(19))      &
     &     GO TO 2040                                                   
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(5)              &
     &     .AND.IWORD2(3).EQ.ALPHA(6)) GO TO 2070                       
!             AT THIS POINT, UNIDENTIFIED PRINT OPTION                  
      WRITE(6,126) IWORD2 
      GO TO 2083 
!                                                                       
!             HERE IF "NOPRNT" CARD                                     
!                                                                       
 1496 IF (IWORD2(1).EQ.ALPHA(7).AND.IWORD2(2).EQ.ALPHA(15)              &
     &     .AND.IWORD2(3).EQ.ALPHA(20).AND.IWORD2(4).EQ.ALPHA(6))       &
     &     GO TO 1750                                                   
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(5)              &
     &     .AND.IWORD2(3).EQ.ALPHA(4)) GO TO 1800                       
      IF (IWORD2(1).EQ.ALPHA(16).AND.IWORD2(2).EQ.ALPHA(18)             &
     &     .AND.IWORD2(3).EQ.ALPHA(24).AND.IWORD2(4).EQ.ALPHA(25)       &
     &     .AND.IWORD2(5).EQ.ALPHA(26)) GO TO 1860                      
      IF (IWORD2(1).EQ.ALPHA(19).AND.IWORD2(2).EQ.ALPHA(20)             &
     &     .AND.IWORD2(3).EQ.ALPHA(1).AND.IWORD2(4).EQ.ALPHA(20)        &
     &     .AND.IWORD2(5).EQ.ALPHA(16).AND.IWORD2(6).EQ.ALPHA(19))      &
     &     GO TO 2050                                                   
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(5)              &
     &     .AND.IWORD2(3).EQ.ALPHA(6)) GO TO 2080                       
      GO TO 2083 
!             FOUND "SCALE" AS SECOND WORD                              
 1498 XSMIN=ARRAY(1) 
      XSMAX=ARRAY(2) 
      YSMIN=ARRAY(3) 
      YSMAX=ARRAY(4) 
      GO TO 2083 
!             HERE IF "SAY" CARD                                        
 1500 IF ( IWORD2(1).NE.IBLANK )  GO TO 335 
 1508 RHO=ARRAY(1) 
      D=ARRAY(2) 
      GO TO 2083 
!                                                                       
!             HERE IF "DES" CARD                                        
!             TEST TO SEE IF PARAX CALLED YET                           
 1510 IF (IPARAX.NE.IZERO) GO TO 1520 
!             TRIED TELDESIGN BEFORE CALLING PARAX                      
      WRITE (6,130) 
      GO TO 2083 
 1520 IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(3)) GO TO 1530 
      IF (IWORD2(1).EQ.ALPHA(3).AND.IWORD2(2).EQ.ALPHA(1)) GO TO 1540 
      IF (IWORD2(1).EQ.ALPHA(4).AND.IWORD2(2).EQ.ALPHA(11)) GO TO 1550 
!             UNKNOWN TYPE OF DESIGN                                    
      WRITE (6,140) IWORD2 
      GO TO 2083 
!             RITCHEY-CHRETIEN DESIGN                                   
 1530 CALL RCDES 
      GO TO 2083 
!             CASSEGRAIN DESIGN                                         
 1540 CALL CADES 
      GO TO 2083 
!             DAHL-KIRKHAM DESIGN                                       
 1550 CALL DKDES 
      GO TO 2083 
!                                                                       
!             HERE IF "ASTOP" CARD                                      
 1560 CALL SURFNO(IWORD2,LEN2,ISURF,JSURF) 
      APSTOP=DFLOAT(JSURF) 
      IAP=APSTOP 
      APSTOP=APSTOP*30.D0 
      IF (APSTOP.LT.ZERO) APSTOP=ZERO 
      GO TO 2083 
!                                                                       
!             HERE IF "LATT" CARD                                       
 1566 LATFLG=IONE 
      DO 1567 I=1,10 
      IF (IWORD2(1).EQ.NUM(I)) GO TO 1568 
 1567 END DO 
      IF (IWORD2(1).EQ.ALPHA(7).AND.IWORD2(2).EQ.ALPHA(5)               &
     &     .AND.IWORD2(3).EQ.ALPHA(14)) GO TO 1568                      
      GO TO 335 
 1568 CALL CNTCOM(NCOM,INPUT2(KBEG),LENGTH) 
      IF (NCOM.LE.IZERO.OR.NCOM.GT.20) GO TO 2083 
      DO 1569 I=1,NCOM 
        CLTRA(INDEX)=ARRAY(I) 
        INDEX=INDEX+IONE 
 1569 END DO 
      LATYPE=ITHREE 
      GO TO 2083 
!                                                                       
!             HERE IF "SPD" CARD                                        
 1570 IMODE=IZERO 
      NNSW=IZERO 
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(5)              &
     &     .AND.IWORD2(3).EQ.ALPHA(3)) GO TO 1590                       
      IF (IWORD2(1).EQ.ALPHA(16).AND.IWORD2(2).EQ.ALPHA(15)             &
     &     .AND.IWORD2(3).EQ.ALPHA(12)) GO TO 1640                      
      IF (IWORD2(1).EQ.ALPHA(6).AND.IWORD2(2).EQ.ALPHA(9)               &
     &     .AND.IWORD2(3).EQ.ALPHA(14)) GO TO 1650                      
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(1)              &
     &     .AND.IWORD2(3).EQ.ALPHA(25)) GO TO 1660                      
      IF (IWORD2(1).EQ.ALPHA(7).AND.IWORD2(2).EQ.ALPHA(5)               &
     &     .AND.IWORD2(3).EQ.ALPHA(14)) GO TO 1670                      
      IF (IWORD2(2).EQ.ALPHA(6).AND.IWORD2(3).EQ.ALPHA(1)               &
     &     .AND.IWORD2(4).EQ.ALPHA(14)) GO TO 1690                      
      IF (IWORD2(1).EQ.ALPHA(6).AND.IWORD2(2).EQ.ALPHA(9)               &
     &     .AND.IWORD2(3).EQ.ALPHA(12).AND.IWORD2(4).EQ.ALPHA(5))       &
     &     GO TO 1585                                                   
      IF (IWORD2(1).EQ.ALPHA(15).AND.IWORD2(2).EQ.ALPHA(16)             &
     &     .AND.IWORD2(3).EQ.ALPHA(6).AND.IWORD2(4).EQ.ALPHA(9)         &
     &     .AND.IWORD2(5).EQ.ALPHA(12).AND.IWORD2(6).EQ.ALPHA(5))       &
     &     GO TO 1583                                                   
      IF (IWORD2(1).EQ.ALPHA(18).AND.IWORD2(2).EQ.ALPHA(9)              &
     &     .AND.IWORD2(3).EQ.ALPHA(13)) GO TO 1587                      
      DO 1580 I=1,10 
        IF (IWORD2(1).EQ.NUM(I)) GO TO 1670 
 1580 END DO 
!             UNKNOWN RAY LATTICE TYPE                                  
      WRITE (6,150) IWORD2 
      GO TO 2083 
!             LATTICE ALREADY STORED ON UNIT 20                         
!             WANT TO ANALYZE MANY OBJECT POINTS AT SAME TIME           
 1583 LATYPE=6 
      NNSW=-1 
      GO TO 1721 
!             LATTICE ALREADY STORED ON UNIT 20                         
!             WANT TO ANALYZE ONE OBJECT POINT AT A TIME                
 1585 LATYPE=5 
      NNSW=-1 
      GO TO 1721 
!             RIM LATTICE                                               
 1587 CLTRA(1)=ARRAY(1) 
      NUMPTS=CLTRA(1) 
      DTHETA=(2*PI/CLTRA(1)) 
      CLTRA(2)=DTHETA 
      LATYPE=7 
      GO TO 1721 
!             RECTANGULAR LATTICE                                       
 1590 CLTRA(1)=ARRAY(1) 
      NUMPTS=CLTRA(1) 
      DX=DSQRT(PI/CLTRA(1)) 
      DX2=DX/TWO 
      DY=DX 
      DY2=DX2 
      CLTRA(2)=DY 
!             INITIALIZE COUNTERS                                       
      J=IZERO 
      N=IZERO 
!             INITIALIZE FIRST COLUMN LOCATION IN X                     
      X=-DX2 
!             START (IMPLIED) COLUMN LOOP                               
 1600 X=X+DX 
      IF (X.GT.ONE) GO TO 1630 
!             K IS NUMBER OF POINTS IN THIS COLUMN;                     
!             J IS USED FOR GENERATING CLTRA ARRAY;                     
      K=IZERO 
      J=J+ITHREE 
!             INITIALIZE Y FOR THIS COLUMN                              
      M=(DSQRT(ONE-X*X)-DY2)/DY 
      M=M+IONE 
      Y=DY2-M*DY 
      IF (DSQRT(X*X+Y*Y).GT.ONE) GO TO 1600 
!             SET UP CLTRA ARRAY FOR THIS COLUMN                        
      CLTRA(J+IONE)=X 
      CLTRA(J+ITWO)=Y 
!             START ROW (IMPLIED) LOOP                                  
 1610 K=K+IONE 
      N=N+IONE 
      Y=Y+DY 
      IF (DSQRT(X*X+Y*Y).GT.ONE) GO TO 1620 
      GO TO 1610 
!             FINISHED COUNTING ROWS IN THIS COLUMN, GO TO NEXT COLUMN  
 1620 CLTRA(J)=K 
      GO TO 1600 
!             FINISHED GENERATING RECTANGULAR ARRAY                     
 1630 CLTRA(1)=N 
      LATYPE=ITHREE 
      GO TO 1721 
!                                                                       
!             POLAR LATTICE                                             
 1640 ANNULI   = DSQRT( ARRAY(1) ) 
      IANN     = ANNULI/TWO 
      ANNULI   = IANN * ITWO 
      ISECT    = ARRAY(1) / ANNULI 
      ISECT    = (ISECT/ITWO) * ITWO 
      SECTRS   = ISECT 
      CLTRA(1) = ANNULI 
      CLTRA(2) = SECTRS 
      LATYPE   = ITWO 
      GO TO 1721 
!             FINRAY LATTICE                                            
 1650 CLTRA(1)=ARRAY(1) 
      LATYPE=IFOUR 
      GO TO 1721 
!             SINGLE RAY                                                
 1660 CLTRA(1)=ARRAY(1) 
      CLTRA(2)=ARRAY(2) 
      NNSW=IONE 
      LATYPE=IONE 
      IMODE=ITWO 
      GO TO 1721 
!              GENERAL LATTICE                                          
 1670 IF (LATFLG.EQ.IONE) GO TO 1685 
!              LATT GEN NOT CALLED YET                                  
      WRITE(6,95) 
      GO TO 2083 
 1685 LATYPE=ITHREE 
      NNSW=-1 
      GO TO 1721 
!                                                                       
!             RAY FAN                                                   
 1690 CLTRA(1)=ARRAY(1) 
      DELTA=TWO/(ARRAY(1)-ONE) 
      LATYPE=ITHREE 
      IMODE=ITWO 
      IF (IWORD2(1).EQ.ALPHA(24)) GO TO 1710 
      IF (IWORD2(1).EQ.ALPHA(25)) GO TO 1700 
!             UNKNOWN TYPE OF RAY FAN                                   
      WRITE (6,160) IWORD2 
      GO TO 1721 
!             YFAN OF RAYS                                              
 1700 CLTRA(2)=DELTA 
      CLTRA(3)=ARRAY(1) 
      CLTRA(4)=ZERO 
      CLTRA(5)=-ONE 
      GO TO 1721 
!             XFAN OF RAYS                                              
 1710 CLTRA(2)=ZERO 
      NRAYS=CLTRA(1) 
      J=ITHREE 
      DO 1720 I=1,NRAYS 
        CLTRA(J)=ONE 
        CLTRA(J+IONE)=-ONE+(I-1)*DELTA 
        CLTRA(J+ITWO)=ZERO 
        J=J+ITHREE 
 1720 END DO 
!                                                                       
!             PERFORM REAL RAY TRACE                                    
!             DETERMINE MAXIMUM SURFACE NUMBER                          
!             THIS IS DETERMINED BY FIRST SURFACE WITH A ZERO INDEX     
 1721 DO 1722 I=1,40 
        IF ( FN(I,1).EQ.ZERO ) GO TO 1723 
 1722 END DO 
 1723 ISURF=I 
      NSURF=ISURF-IONE 
      ISKEW=IONE 
      IF ((IPRINT/2)*2.NE.IPRINT) GO TO 1726 
!             DETERMINE COLORS FROM CARD                                
      CALL CNTCOM (NCOM,INPUT2(KBEG),LENGTH) 
      IF (NCOM.GT.IONE+NNSW) GO TO 1724 
!             CAN'T IDENTIFY COLORS                                     
      WRITE (6,90) 
      GO TO 2083 
 1724 IF (NCOM.GT.4+NNSW) NCOM=4+NNSW 
      NCOL=NCOM-(1+NNSW) 
      DO 1725 I=1,NCOL 
        ICOL(I)=ARRAY(I+1+NNSW) 
        IF (ICOL(I).LT.1) ICOL(I)=1 
        IF (ICOL(I).GT.3) ICOL(I)=3 
 1725 END DO 
!             DETERMINE MODE                                            
 1726 NCOM=IZERO 
      IX=IZERO 
      DO 1727 I=1,NSURF 
        IF (XDISP(I).EQ.ZERO.AND.TILTX(I).EQ.ZERO) GO TO 1727 
        IX=IONE 
        GO TO 1728 
 1727 END DO 
 1728 IF ( ( (IX.EQ.IONE).OR.(HXINIT.NE.ZERO).OR.(HXDEL.NE.ZERO) )      &
     &     .AND. IMODE.NE.ITWO) IMODE = IONE                            
!             CONVERT WAVELENGTH FROM NANOMETERS TO SYSTEM UNITS        
      DEFWV(1)=632.8 
      DEFWV(2)=587.6 
      DEFWV(3)=486.1 
      NEND=NCOL 
      IF ((IPRINT/2)*2.NE.IPRINT) NEND=3 
      DO 1730 I=1,NCOL 
      IF ((IPRINT/2)*2.NE.IPRINT) GO TO 1729 
      WVTEMP=WAVE(ICOL(I)) 
      IWVFLG(ICOL(I))=0 
      IF (WAVE(ICOL(I)).EQ.0) IWVFLG(ICOL(I))=1 
      IF (WAVE(ICOL(I)).EQ.0) WVTEMP=DEFWV(ICOL(I)) 
      IF (UFLAG.EQ.ONE) WAVL(ICOL(I))=WVTEMP*1.D-6 
      IF (UFLAG.EQ.TWO) WAVL(ICOL(I))=WVTEMP*1.D-7 
      IF (UFLAG.EQ.THREE) WAVL(ICOL(I))=(WVTEMP*1.D-6)/25.4D0 
      GO TO 1730 
 1729 WVTEMP=WAVE(I) 
      IF (WAVE(I).EQ.0) WVTEMP=DEFWV(I) 
      IF (UFLAG.EQ.ONE) WAVL(I)=WVTEMP*1.D-6 
      IF (UFLAG.EQ.TWO) WAVL(I)=WVTEMP*1.D-7 
      IF (UFLAG.EQ.THREE) WAVL(I)=(WVTEMP*1.D-6)/25.4D0 
 1730 END DO 
!             CONVERT LINEAR OBJECT SIZES TO ANGLES                     
 1731 IF (IFLAG.EQ.IZERO) GO TO 1732 
      HYINIT = -DATAN(HYINIT/(S-D)) * (180.D0/PI) 
      HYDEL = -DATAN(HYDEL/(S-D)) * (180.D0/PI) 
      HXINIT = -DATAN(HXINIT/(S-D)) * (180.D0/PI) 
      HXDEL = -DATAN(HXDEL/(S-D)) * (180.D0/PI) 
      IFLAG = IZERO 
! 4/19/84 THE ABOVE LINE INSERTED TO CORRECT A PROBLEM WITH             
!         Y OBJECT HEIGHTS;  IFLAG = 1 => H?INIT AND H?DEL              
!         ARE LINEAR                                                    
!         IFLAG = 0 => H?INIT AND H?DEL ARE ANGLES IN DEGREES           
!                                                                       
 1732 WRITE (6,40) 
      CALL SKEW (NFOC,LFOC,IFOC) 
      WRITE (6,40) 
      GO TO 2083 
!                                                                       
!             TURN ON MTF PRINT OPTION                                  
 1740 SMAX=ARRAY(1) 
      IF ((((IPRINT/32)/2)*2).NE.(IPRINT/32)) GO TO 2083 
      IPRINT=IPRINT+32 
      GO TO 2083 
!             TURN OFF MTF PRINT                                        
 1750 SMAX=ARRAY(1) 
      IF ((((IPRINT/32)/2)*2).EQ.(IPRINT/32)) GO TO 2083 
      IPRINT=IPRINT-32 
      GO TO 2083 
!             TURN ON MTF PLOT OPTION                                   
 1760 SMAX=ARRAY(1) 
      IF ((((IPLTPR/8)/2)*2).NE.(IPLTPR/8)) GO TO 2083 
      IPLTPR=IPLTPR+8 
      GO TO 2083 
!              TURN OFF MTF PLOT                                        
 1770 SMAX=ARRAY(1) 
      IF ((((IPLTPR/8)/2)*2).EQ.(IPLTPR/8)) GO TO 2083 
      IPLTPR=IPLTPR-8 
      GO TO 2083 
!                                                                       
!             TURN ON RED PRINT OPTION                                  
 1790 RSMAX = ARRAY(1) 
      IF ((((IPRINT/16)/2)*2).NE.(IPRINT/16)) GO TO 2083 
      IPRINT=IPRINT+16 
      GO TO 2083 
!             TURN OFF RED PRINT                                        
 1800 RSMAX=ARRAY(1) 
      IF ((((IPRINT/16)/2)*2).EQ.(IPRINT/16)) GO TO 2083 
      IPRINT=IPRINT-16 
      GO TO 2083 
!             TURN ON RED PLOT OPTION                                   
 1810 RSMAX=ARRAY(1) 
      IF ((((IPLTPR/4)/2)*2).NE.(IPLTPR/4)) GO TO 2083 
      IPLTPR=IPLTPR+4 
      GO TO 2083 
!             TURN OFF RED PLOT                                         
 1820 RSMAX=ARRAY(1) 
      IF ((((IPLTPR/4)/2)*2).EQ.(IPLTPR/4)) GO TO 2083 
      IPLTPR=IPLTPR-4 
      GO TO 2083 
!                                                                       
!             HERE IF "BIL" CARD                                        
 1830 Y0(I)=THREE 
      GO TO 2083 
!                                                                       
!             TURN ON SURFACE-BY-SURFACE PRINT OPTION                   
 1850 IF ((((IPRINT/6)/2)*2).NE.(IPRINT/6)) GO TO 2083 
      IPRINT=IPRINT+6 
      GO TO 2083 
!             TURN OFF PRINT CONTROL                                    
 1860 IF ((((IPRINT/6)/2)*2).EQ.(IPRINT/6)) GO TO 2083 
      IPRINT=IPRINT-6 
      GO TO 2083 
!                                                                       
!             HERE IF TURNING ON LEPRT PRINT CODE                       
 1880 IF (IPRINT.NE.IZERO.AND.(IPRINT/ITWO)*ITWO.NE.IPRINT) GO TO 2083 
      IPRINT=IPRINT+IONE 
      GO TO 1721 
!                                                                       
!             HERE IF "SC" CARD                                         
 1900 FACTOR=ARRAY(1) 
      CALL LENSCL (FACTOR) 
      GO TO 2083 
!                                                                       
!             HERE IF "WV" CARD                                         
 1910 WAVE(1)=ARRAY(1) 
      WAVE(2)=ARRAY(2) 
      WAVE(3)=ARRAY(3) 
      GO TO 2083 
!                                                                       
!             HERE IF "INSERT" CARD                                     
 1920 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      GO TO 2083 
!                                                                       
!             HERE IF "DELETE" CARD                                     
 1930 CALL SURFNO (IWORD2,LEN2,ISURF,JSURF) 
      GO TO 2083 
!                                                                       
!             HERE IF "REFS" CARD                                       
 1940 CALL SURFNO(IWORD2,LEN2,ISURF,JSURF) 
      IREF=ISURF 
      GO TO 2083 
!                                                                       
!             WANT MULTIPLE SPOT PLOT                                   
 1945 IALLPL=1 
      GO TO 1955 
!                                                                       
!             WANT SPOT UNIT (OR RADII) PLOT                            
 1950 IALLPL=0 
 1955 IPEN = ARRAY(1) 
      IF ( IPEN.LE.IZERO .OR. IPEN.GT.IFOUR ) IPEN = 1 
      CALL NEWPEN(IPEN) 
      IF (IPLTPR.NE.IZERO.AND.(IPLTPR/ITWO)*ITWO.NE.IPLTPR) GO TO 2083 
      IPLTPR = IPLTPR + IONE 
      GO TO 2083 
!             DON'T WANT SPOT UNIT (OR RADII) PLOT                      
 1960 IF ( IPLTPR.NE.IZERO.AND.(IPLTPR/ITWO)*ITWO.NE.IPLTPR)            &
     &     IPLTPR = IPLTPR-IONE                                         
      IALLPL=0 
      GO TO 2083 
!                                                                       
!             DON'T WANT MULTIPLE SPOT PLOT                             
 1965 IALLPL=0 
      GO TO 2083 
!                                                                       
!             HERE IF "STATISTICS" CARD                                 
 1970 IF (IWORD2(1).EQ.IBLANK) GO TO 1975 
      IF (IWORD2(1).EQ.ALPHA(6).AND.IWORD2(2).EQ.ALPHA(9)               &
     &     .AND.IWORD2(3).EQ.ALPHA(12).AND.IWORD2(4).EQ.ALPHA(5))       &
     &     GO TO 1980                                                   
      GO TO 335 
 1975 FAKEA=ARRAY(1) 
      IPR20=IZERO 
      GO TO 2083 
 1980 FAKEA=ARRAY(1) 
      IPR20=IONE 
      GO TO 2083 
!             HERE IF CURVATURE SOLVE                                   
 1990 CALL SURFNO(IWORD2,LEN2,ISURF,JSURF) 
      SXNU(JSURF)=ARRAY(2) 
      ORDN(JSURF,1)=ORDN(JSURF,1)+ONE 
      ORDN(JSURF,2)=ORDN(JSURF,2)+ONE 
      ORDN(JSURF,3)=ORDN(JSURF,3)+ONE 
      GO TO 2083 
!             HERE IF CLEAR APERTURE SOLVE                              
 2000 CALL SURFNO(IWORD2,LEN2,ISURF,JSURF) 
      Y0(JSURF)=ARRAY(2) 
      RDSPAC(JSURF)=ARRAY(3) 
      ORDN(JSURF,1)=ORDN(JSURF,1)+FOUR 
      ORDN(JSURF,2)=ORDN(JSURF,2)+FOUR 
      ORDN(JSURF,3)=ORDN(JSURF,3)+FOUR 
      GO TO 2083 
!             HERE IF THICKNESS SOLVE                                   
 2010 CALL SURFNO(IWORD2,LEN2,ISURF,JSURF) 
      SXY(JSURF)=ARRAY(2) 
      ORDN(JSURF,1)=ORDN(JSURF,1)+TWO 
      ORDN(JSURF,2)=ORDN(JSURF,2)+TWO 
      ORDN(JSURF,3)=ORDN(JSURF,3)+TWO 
      GO TO 2083 
!             HERE IF REMOVING SOLVES                                   
 2020 ORDN(JSURF,1)=ZERO 
      ORDN(JSURF,2)=ZERO 
      ORDN(JSURF,3)=ZERO 
      GO TO 2083 
!                                                                       
!             TURN ON OPTION FOR PRINTING PRIME IMAGE COORDINATES       
 2040 IF ((((IPRINT/8)/2)*2).NE.(IPRINT/8)) GO TO 2083 
      IPRINT=IPRINT+8 
      GO TO 2083 
!             TURN OFF PRINT CONTROL                                    
 2050 IF ((((IPRINT/8)/2)*2).EQ.(IPRINT/8)) GO TO 2083 
      IPRINT=IPRINT-8 
      GO TO 2083 
!                                                                       
!             HERE IF "PRINT REF" CARD                                  
!                                                                       
!             TURN ON PRINT CONTROL                                     
 2070 IF ((((IPRINT/2)/2)*2).NE.(IPRINT/2)) GO TO 2083 
      IPRINT=IPRINT+2 
      GO TO 2083 
!             TURN OFF PRINT CONTROL                                    
 2080 IF ((((IPRINT/2)/2)*2).EQ.(IPRINT/2)) GO TO 2083 
      IPRINT=IPRINT-2 
 2083 CONTINUE 
      IF (LSL.EQ.I81) GO TO 170 
      IENDO=LAST-LSL 
      DO 2085 I=1,IENDO 
        INPUT(I)=INPUT(I+LSL) 
 2085 END DO 
      IBEGDO=LAST+1 
      DO 2086 I=IBEGDO,80 
        INPUT(I)=IBLANK 
 2086 END DO 
      GOTO 172 
!                                                                       
!****************************************************************       
!                                                                       
!             CLOSE VAX FILES                                           
!                                                                       
 2090 CLOSE (UNIT=5,  DISPOSE='KEEP') 
      CLOSE (UNIT=10, DISPOSE='KEEP') 
      CLOSE (UNIT=20, DISPOSE='KEEP') 
      CLOSE (UNIT=34, DISPOSE='DELETE') 
      CLOSE (UNIT=40, DISPOSE='DELETE') 
      CLOSE (UNIT=50, DISPOSE='DELETE') 
!                                                                       
!****************************************************************       
!                                                                       
      STOP 
      END                                           
!                                                                       
!******************************************************                 
      SUBROUTINE BESSEL(X,N,BJ,D,IER) 
!******************************************************                 
      IMPLICIT REAL*8(A-H,O-Z), INTEGER*4(I-N) 
      BJ=0.D0 
      IF(N.GE.0) GO TO 10 
          IER= 1 
          GO TO 90 
   10 IF (X .GT. 0.D0) GO TO 30 
         IF(X .LT. 0.D0) GO TO 20 
            IER=0 
            BJ=1.0 
            GO TO 90 
   20    IER=2 
         GO TO 90 
!                                                                       
   30 IF (X .GT. 15.D0) GO TO 40 
         MAX= 20.D0 + 10.D0 * X - X**2 / 3.D0 
         GO TO 50 
   40    MAX= 90.D0 + X / 2.D0 
!                                                                       
   50 IF (N .LT. MAX) GO TO 60 
         IER= 4 
         GO TO 90 
!                                                                       
   60 IER= 0 
      BPREV= 0. 
!                                                                       
!             COMPUTE STARTING VALUE OF M                               
!                                                                       
      MA= X + 6.D0 
      IF (X .GE. 5.D0) MA= 1.4D0 * X + 60.D0 / X 
      MB= N + IDINT(X) / 4 + 2 
      MZERO= MAX0(MA, MB) 
      DO 80 M= MZERO, MAX, 3 
!                                                                       
!             SET INITIAL VALUES FOR F(M) AND F(M-1)                    
!                                                                       
      FM= 0.D0 
      FM1= 1.0D-28 
      ALPHA= 0.D0 
      ST= -1.D0 
      IF (M .NE. (M / 2) * 2) ST= 1.D0 
      M2= M - 2 
      DO 70 K= 1, M2 
      MK= M - K 
!                                                                       
!             CALCULATE FM2 FROM RECURRENCE RELATIONSHIP                
!                                                                       
      FM2= 2.D0 * DFLOAT(MK) * FM1 / X - FM 
      FM= FM1 
      FM1= FM2 
      IF (MK .EQ. N + 1) BJ= FM2 
      ST= -ST 
   70 ALPHA= ALPHA + FM2 * (1.D0 + ST) 
      F0= 2.D0 * FM1 / X - FM 
      IF (N .EQ. 0) BJ= F0 
      ALPHA= ALPHA + F0 
      BJ= BJ / ALPHA 
      IF (DABS(BJ - BPREV) .LT. DABS(D * BJ)) GO TO 90 
   80 BPREV= BJ 
      IER= 3 
   90 RETURN 
      END                                           
!                                                                       
!*************************************************                      
      SUBROUTINE CNTCOM(NCOM,ICARD,LENGTH) 
!*************************************************                      
!                                                                       
!                FIND THE NUMBER OF VARIABLE VALUES ON THIS CARD BY     
!                COUNTING THE NUMBER OF COMMAS OR BLANKS                
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      DIMENSION ICARD(LENGTH) 
      DATA IVAL/1H,/,IBLANK/1H /,ISLASH/1H// 
      NCOM = 0 
      DO 40 I=1,LENGTH 
      IF (ICARD(I).EQ.IBLANK) GOTO 40 
      ISTART=I 
      GO TO 50 
   40 END DO 
      GOTO 110 
   50 J0=0 
      DO 100 I=ISTART,LENGTH 
        IF (ICARD(I).NE.IVAL.AND.ICARD(I).NE.IBLANK) GOTO 70 
        IF (J0.EQ.0) J0=I 
        IF (J0.EQ.I) NCOM=NCOM+1 
        GO TO 100 
   70 J0=0 
        IF (I.EQ.LENGTH.AND.ICARD(I).NE.IBLANK.AND.ICARD(I).NE.IVAL)    &
     &      NCOM=NCOM + 1                                               
  100 END DO 
  110 RETURN 
      END                                           
!                                                                       
!******************************************************                 
      SUBROUTINE FINDE(STOR,ICARD,ICNT,LENGTH) 
!******************************************************                 
!                                                                       
!                THIS SUBR DECODES ONE FLOATING POINT NUMBER FROM       
!                THE CARD IMAGE ICARD, USING ICNT NUMBER OF CHARACTERS, 
!                AND PUTS THE VALUE IN STOR. ICARD IS AN INTEGER ARRAY  
!                ORIGINALLY READ UNDER XXA1 FORMAT, IE ONE HOLLERITH    
!                CHARACTER PER WORD.                                    
!                                                                       
!                INITIALIZE ALL SWITCHES                                
!                SWITCH   VALUE   USE                                   
!                ISSW       1     NEXT SIGN DECODED IS SIGN OF MANTISSA 
!                           2     NEXT SIGN DECODED IS SIGN OF EXPONENT 
!                           3     BOTH SIGNS ALREADY EXTRACTED          
!                IPSW       1     DEC POINT NOT FOUND IN MANTISSA       
!                           2     DEC POINT FOUND, ACCUMULATE ADJUSTMENT
!                                 FACTOR DEC                            
!                IESW       1     MANTISSA IS BEING EXTRACTED           
!                           2     EXPONENT IS BEING EXTRACTED           
!                                                                       
!                IDENTIFIER USE                                         
!                DEC        DECREMENT OF MANTISSA                       
!                IPOS       POSITION INDEX OF INPUT IMAGE CHARACTER     
!                           BEING EXAMINED                              
!                VAR        VALUE OF DIGIT BEING EXTRACTED              
!                TVAR       VALUE OF UNDECREMENTED MANTISSA OR EXPO-    
!                           NENT NOW BEING EXTRACTED                    
!                SING       =1 IF SIGN OF MANT IS + OR NOT SPECIFIED    
!                           =-1 IF SIGN OF MANT IS SPECIFIED MINUS      
!                ESING      AS FOR SING, SIGN ADJUSTMENT MULTIPLIER     
!                           OF EXPONENT                                 
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      DIMENSION ICARD(LENGTH),ILIST(17) 
      DATA  ILIST           /1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1H+&
     &,1HE,1H-,1H.,1H,,1HD,1H /                                         
      ISSW=1 
      DEC=1. 
      IPOS=1 
      IPSW=1 
      IESW=1 
      VAR=0. 
      TVAR=0. 
      STOR=0. 
      SING=1. 
      ESING=1. 
!                                                                       
!                EXAMINE LIST FOR CHARACTER MATCH                       
!                                                                       
   10 DO 20 I=1,17 
        IF ( ICARD(IPOS).EQ.ILIST(I) ) GO TO 40 
   20 END DO 
!                                                                       
!                IF NO MATCH, TREAT AS ZERO                             
!                                                                       
   30 I=1 
      GO TO 50 
!                                                                       
!                MATCH FOUND, GO TO PROPER ROUTINE                      
!                                                                       
   40 J=I-9 
      GO TO (50,80,110,80,120,130,110,130),J 
!                                                                       
!                DECODE DIGIT, INCREMENT TVAR                           
!                                                                       
   50 VAR=I-1 
      TVAR=TVAR*10.+VAR 
!                                                                       
!                IF PERIOD FOUND, INCREASE DECREMENT VARIABLE - WE      
!                ARE USING BORROWED DECIMAL PLACES                      
!                                                                       
      IF ( IPSW.EQ.2)DEC=DEC*10. 
!                                                                       
!                IF FIRST DIGIT NOT SIGN, WE CALL MANTISSA POSITIVE     
!                                                                       
   60 IF ( IPOS.EQ.1 ) ISSW=2 
!                                                                       
!                INCREMENT POSITION GO TO NEXT CHAR                     
!                                                                       
   70 IPOS=IPOS+1 
      GO TO 10 
!                                                                       
!                WE HAVE FOUND A SIGN, WHICH ONE                        
!                                                                       
   80 GO TO (90,100,30),ISSW 
!                                                                       
!                MANTISSA SIGN                                          
!                                                                       
   90 ISSW=2 
      SING=12-I 
      GO TO 70 
!                                                                       
!                EXPONENT SIGN                                          
!                                                                       
  100 ISSW=3 
      ESING=12-I 
!                                                                       
!                HERE IF D OR E FOUND IN IMAGE, MAKE STOR TAKE VALUE OF 
!                MANTISSA, RESET TVAR, SET E SWITCH                     
!                                                                       
      GO TO (110,70),IESW 
  110 STOR=SING*TVAR/DEC 
      TVAR=0. 
      IESW=2 
      GO TO 70 
!                                                                       
!                PERIOD SWITCH SET IF PERIOD FOUND                      
!                                                                       
  120 IPSW=2 
      GO TO 60 
!                                                                       
!                COMMA OR BLANK FOUND, END OF DECODE                    
!                                                                       
  130 GO TO (140,160),IESW 
!                                                                       
!                SET STOR TO MANTISSA VALUE IF NO E FOUND               
!                                                                       
  140 STOR=SING*TVAR/DEC 
  150 IF (IPOS.EQ.LENGTH) GOTO 156 
      IPOS1=IPOS+1 
      DO 155 J=IPOS1,LENGTH 
      IF (ICARD(J).EQ.ILIST(15).OR.ICARD(J).EQ.ILIST(17)) GOTO 155 
      ICNT=J-1 
      GOTO 156 
  155 END DO 
      ICNT=J-1 
  156 CONTINUE 
      RETURN 
!                                                                       
!                IF E FOUND, STOR IS ALREADY VALUE OF MANTISSA          
!                                                                       
  160 IPOW=TVAR*ESING 
      STOR=STOR*(10.D0**IPOW) 
      GO TO 150 
      END                                           
!                                                                       
!****************************************************                   
      SUBROUTINE FINRAY(NFOC) 
!****************************************************                   
!      6-1-81 VERSION                                                   
      IMPLICIT REAL*8(A-H,O-Z) 
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
!                                                                       
       DESX=0 
       DESY=0 
       N=0 
       IFLAG=0 
       Z=0.0 
       PI=3.141592653589793D0 
       TEMP = S-D 
       HYMAX = HYINIT + (RNOBJ-1.)*HYDEL 
       HYMAX = -DTAN( HYMAX*PI/180. ) * TEMP 
       HXMAX = HXINIT + (RNOBJ-1.)*HXDEL 
       HXMAX = -DTAN( HXMAX*PI/180. ) * TEMP 
       AY=DMAX1(DABS(HYMAX),DABS(HYINIT)) 
       AX=DMAX1(DABS(HXMAX),DABS(HXINIT)) 
       RHONEW=RHO+DSQRT(AX*AX+AY*AY) 
       A=RHONEW/TEMP 
       THISHY = -DTAN( HYINIT*PI/180. ) * TEMP 
       THISHX = -DTAN( HXINIT*PI/180. ) * TEMP 
       OBJNUM = 1.0 
       IF (NFOC.EQ.0) GO TO 5 
       THISHY = HYINIT + (NFOC-1)*HYDEL 
       THISHX = HXINIT + (NFOC-1)*HXDEL 
       THISHY = -DTAN( THISHY*PI/180. ) * TEMP 
       THISHX = -DTAN( THISHX*PI/180. ) * TEMP 
    5  OMEGA=DSQRT((2.*PI/CLTRA(1))*(1.-1./DSQRT(1.+A*A))) 
       BETA=OMEGA/2. 
       ALPHA=BETA 
       ALPHAN=ALPHA 
   10  CSALPH=DCOS(ALPHA) 
       P=DSIN(BETA)/CSALPH 
       IF( DABS(P).GE.1.0D0) GO TO 20 
       THETA = DASIN(P) 
       COSGAM=CSALPH*DCOS(THETA) 
       TANGAM=DSQRT((1./COSGAM**2)-1.) 
       IF( TANGAM.GT.DABS(A) ) GO TO 20 
       N=N+1 
       QX=DSIN(BETA) 
       QY=DSIN(ALPHA) 
       QZ=COSGAM 
       TEMP1 = (S-D)/QZ 
       X=TEMP1*QX + THISHX 
       Y=TEMP1*QY + THISHY 
       WRITE(40)X,Y,Z,QX,QY,QZ,THISHY,THISHX,DESX,DESY 
       QYM=-QY 
       YM=TEMP1*QYM + THISHY 
       WRITE(40)X,YM,Z,QX,QYM,QZ,THISHY,THISHX,DESX,DESY 
       IF(IMODE.NE.1) GO TO 15 
       QXM=-QX 
       XM=TEMP1*QXM + THISHX 
       WRITE(40)XM,YM,Z,QXM,QYM,QZ,THISHY,THISHX,DESX,DESY 
       WRITE(40)XM,Y,Z,QXM,QY,QZ,THISHY,THISHX,DESX,DESY 
   15  IFLAG=0 
       ALPHA=ALPHA+OMEGA 
       GO TO 10 
   20  IF(IFLAG.EQ.1) GO TO 25 
       BETA=BETA+OMEGA 
       ALPHA=ALPHAN 
       IFLAG=1 
       GO TO 10 
   25  CONTINUE 
       IF(NFOC.NE.0) GO TO 30 
       THISHY = HYINIT + OBJNUM*HYDEL 
       THISHX = HXINIT + OBJNUM*HXDEL 
       THISHY = -DTAN( THISHY*PI/180. ) * TEMP 
       THISHX = -DTAN( THISHX*PI/180. ) * TEMP 
       OBJNUM = OBJNUM + 1.0 
       IF ( OBJNUM.GT.RNOBJ ) GO TO 30 
       GO TO 5 
   30  RETURN 
      END                                           
!                                                                       
!*****************************************************                  
      SUBROUTINE FLOTIN(IND,ICARD,BUFR,LENGTH) 
!*****************************************************                  
!                                                                       
!                DECODE A STRING OF VARIABLES FROM CARD IMAGE           
!                INTO ONE CONTIGUOUS ARRAY                              
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      DIMENSION ICARD(LENGTH),BUFR(LENGTH) 
      CALL CNTCOM(NCOM,ICARD,LENGTH) 
      CALL FREARA(NCOM,ICARD,BUFR(IND+1),LENGTH) 
      RETURN 
      END                                           
!                                                                       
!*****************************************************                  
      SUBROUTINE FOCUS(ZF,BLOB,N,FID) 
!*****************************************************                  
      IMPLICIT REAL *8 (A-H,O-Z) 
!                                                                       
!             FOCUS DETERMINES IMAGE DISTANCE FOR MINIMUM SPOT SIZE     
!                                                                       
      DIMENSION ZF(10),BLOB(10) 
      A1=(BLOB(1)-BLOB(2))/(ZF(1)-ZF(2)) 
      N1=N-1 
      A2=(BLOB(N1)-BLOB(N))/(ZF(N1)-ZF(N)) 
      B1=BLOB(2)-A1*ZF(2) 
      B2=BLOB(N)-A2*ZF(N) 
      FID=(B2-B1)/(A1-A2) 
      RETURN 
      END                                           
!                                                                       
!*****************************************************                  
      SUBROUTINE FREARA(NCOM,ICARD,BUFR,LENGTH) 
!*****************************************************                  
!                                                                       
!             FREE FORM FLOATING POINT ARRAY INPUT                      
!                                                                       
!             CALLING SEQUENCE                                          
!             VARIABLE  USE                                             
!                                                                       
!             NCOM      NUMBER OF VALUES ON CARD IMAGE                  
!             ICARD     FWA OF CARD IMAGE                               
!             BUFR      FWA OF ARRAY INTO WHICH TO DECODE               
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      DIMENSION ICARD(LENGTH),BUFR(LENGTH) 
      IND=1 
      INDB=1 
      IF(NCOM.EQ.0)RETURN 
      DO 10 I=1,NCOM 
!             CALL FINDE TO DECODE ONE FLOATING POINT NUMBER            
        CALL FINDE(BUFR(INDB),ICARD(IND),ICNT,LENGTH) 
        INDB=INDB+1 
        IND=IND+ICNT 
   10 END DO 
      RETURN 
      END                                           
!                                                                       
!*******************************************************                
      SUBROUTINE HEADIN(LAMDX) 
!*******************************************************                
!                                                                       
!             HEADIN CALLED WITH COLOR NUMBER                           
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      COMMON /HEAD/   LINES,IPAGE,NSYS,LAMDA,NAME(160) 
   10 FORMAT(//' SYSTEM NO. ',I3,7X,'GENOPTICS - A GENERAL OPTICAL',    &
     &      ' SYSTEMS EVALUATION PROGRAM')                              
   20 FORMAT(13X,80A1) 
   30 FORMAT('       COLOR',I3) 
   40 FORMAT(1H ) 
      WRITE(6,10) NSYS 
!                                                                       
!             IF LAMDX IS ZERO, PRINT NO COLOR NUMBER                   
!                                                                       
      IF ( LAMDX ) 60,70,50 
   50 WRITE(6,30) LAMDX 
      GO TO 70 
   60 WRITE(6,30) LAMDA 
   70 WRITE(6,20) NAME 
   80 WRITE(6,40) 
      RETURN 
      END                                           
!                                                                       
!*****************************************************                  
      SUBROUTINE LENSCL(FACTOR) 
!*****************************************************                  
!                                                                       
!               THIS ROUTINE ADDED 9-12-80                              
!               IT SCALES A PREVIOUSLY DEFINED LENS SYSTEM              
!               BY AN AMOUNT "FACTOR"; ONLY FIRST ORDER                 
!               PROPERTIES ARE AFFECTED.                                
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
!                                                                       
      EQUIVALENCE (TRASH(1),XSMIN),(TRASH(2),XSMAX),(TRASH(3),YSMIN) 
      EQUIVALENCE (TRASH(4),YSMAX),(TRASH(5),FCODE),(TRASH(6),FXYJ) 
      EQUIVALENCE (TRASH(7),FXY),(TRASH(8),FXNU),(TRASH(9),FBY) 
      EQUIVALENCE (TRASH(10),FBNU) 
!                                                                       
      DATA ZERO/0.D0/ 
!                                                                       
      DIMENSION FACT(4) 
!             IF FACTOR LE 0, PRINT ERROR MESSAGE AND RETURN            
      IF ( FACTOR.LE.ZERO ) GO TO 30 
!               FACTOR IS VALID, PERFORM THE SCALING                    
      DO 4 I=1,40 
        IF (FN(I,1).EQ.ZERO) GO TO 8 
    4 END DO 
    8 NUMSRF = I-1 
      IF (NUMSRF.EQ.0) GO TO 50 
      DO 20 I=1,NUMSRF 
!               SCALE BASIC SURFACE PARAMETERS                          
         T(I)      = T(I)*FACTOR 
         R(I)      = R(I)*FACTOR 
         C(I)      = C(I)/FACTOR 
         RX(I)     = RX(I)*FACTOR 
         CX(I)     = CX(I)/FACTOR 
!               DISPLACEMENTS ARE SCALED, TILTS ARE NOT                 
         XDISP(I)  = XDISP(I)*FACTOR 
         YDISP(I)  = YDISP(I)*FACTOR 
!               GRATING SPACING IS SCALED                               
         RDSPAC(I) = RDSPAC(I)/FACTOR 
!               ACONIC COEFFICIENTS ARE SCALED ACCORDING TO POWERS      
         FACTSQ    = FACTOR*FACTOR 
         FACT(1)   = FACTSQ*FACTOR 
         FACT(2)   = FACTSQ*FACT(1) 
         FACT(3)   = FACTSQ*FACT(2) 
         FACT(4)   = FACTSQ*FACT(3) 
         DO 10 J=1,4 
            COEF(I,J) = COEF(I,J)/FACT(J) 
   10    CONTINUE 
!               MASK PARAMETERS ARE SCALED                              
         FMASK(I)  = FMASK(I)*FACTOR 
         IF( FAKEC(I).EQ.0.D0 )GO TO 20 
         FMASK(I)  = FMASK(I)/FACTOR 
         XMN(I)    = XMN(I)*FACTOR 
         XMX(I)    = XMX(I)*FACTOR 
         YMN(I)    = YMN(I)*FACTOR 
         YMX(I)    = YMX(I)*FACTOR 
   20 END DO 
      XSMIN  = XSMIN*FACTOR 
      XSMAX  = XSMAX*FACTOR 
      YSMIN  = YSMIN*FACTOR 
      YSMAX  = YSMAX*FACTOR 
      IF (UFLAG.EQ.1.0 .AND. FACTOR.EQ. .1D0 ) UFLAG=2 
      IF (UFLAG.EQ.1.0 .AND. FACTOR.EQ. .03937007874D0 ) UFLAG=3 
      IF (UFLAG.EQ.2.0 .AND. FACTOR.EQ. 10.D0 ) UFLAG=1 
      IF (UFLAG.EQ.2.0 .AND. FACTOR.EQ. .3937007874D0 ) UFLAG=3 
      IF (UFLAG.EQ.3.0 .AND. FACTOR.EQ. 25.4D0 ) UFLAG=1 
      IF (UFLAG.EQ.3.0 .AND. FACTOR.EQ. 2.54D0 ) UFLAG=2 
!                                                                       
      FXYJ   = FXYJ*FACTOR 
      D      = D*FACTOR 
      RHO    = RHO*FACTOR 
      FOCL   = FOCL*FACTOR 
      DO 25 I=1,3 
      WAVL(I) = WAVL(I)*FACTOR 
   25 END DO 
      GO TO 50 
!             PRINT ERROR MESSAGE(FOR BAD SCALE FACTOR)                 
   30 WRITE(6,40) FACTOR 
   40 FORMAT(1X,10X,'BAD SCALE FACTOR GIVEN ( ',1PE15.7,' ) ; ',        &
     &       'NO SCALING PERFORMED'/)                                   
!               FINISHED SCALING, RETURN TO MAIN PROGRAM                
   50 RETURN 
      END                                           
!                                                                       
!********************************************************               
      SUBROUTINE MAVEC(PRODA,MULT) 
!********************************************************               
!                                                                       
!                 MAVEC PERFORMS MATRIX*VECTOR AND                      
!                VECTOR*MATRIX MULTIPLICATION                           
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      COMMON /ROT/ ROTA(9) 
      DIMENSION PRODA(3) 
      DIMENSION PROD(3) 
      REAL *8 MULT(3) 
      IEN=7 
      INC=3 
      JNC=1 
      GO TO 10 
      ENTRY VECMA(PRODA,MULT) 
      IEN=3 
      JNC=3 
      INC=1 
   10 IST=1 
      DO 30 INT=1,3 
        PROD(INT)=0. 
        DO 20 I=IST,IEN,INC 
          JNT=3-(IEN-I)/INC 
          PROD(INT)=PROD(INT)+MULT(JNT)*ROTA(I) 
   20   CONTINUE 
        IST=IST+JNC 
        IEN=IEN+JNC 
   30 END DO 
      DO 40 I=1,3 
   40 PRODA(I)=PROD(I) 
      RETURN 
      END                                           
!                                                                       
!********************************************************               
      SUBROUTINE PARAX 
!********************************************************               
!                                                                       
!     PURPOSE                                                           
!            PARAXIAL RAY TRACING SUBROUTINE                            
!                                                                       
!     LAST UPDATE 4/10/84 BY JOHN PARKER (LINE 321 OF PARAX)            
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC              
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
!                                                                       
!                                                                       
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /HEAD/   LINES,IPAGE,NSYS,LAMDA,NAME(160) 
      COMMON /PUPIL/  ENPUPR,ENPUPL,EXPUPR,EXPUPL 
      COMMON /ABCODE/ ISURFX,ISURXN,RSAVE,TM3,DUM(40),RR(40),RCON(40) 
      DIMENSION RNU(40),B(40),U(40),BU(40),F(40),CC(40),E(40),P(40) 
      DIMENSION ACH(40),BCH(40),BAS(40),FAS(40),CAS(40),EAS(40) 
      DIMENSION X(40),BY(40),XY(40) 
!                                                                       
      DATA EPS/1.D-9/,ZERO/0.D0/,ONE/1.D0/,TWO/2.D0/ 
!                                                                       
      EQUIVALENCE (TRASH(1),XSMIN),(TRASH(2),XSMAX),(TRASH(3),YSMIN) 
      EQUIVALENCE (TRASH(4),YSMAX),(TRASH(5),FCODE),(TRASH(6),FXYJ) 
      EQUIVALENCE (TRASH(7),FXY),(TRASH(8),FXNU),(TRASH(9),FBY) 
      EQUIVALENCE (TRASH(10),FBNU) 
      EQUIVALENCE (RNU(1),CAS(1)) 
!                                                                       
   10 FORMAT (/10X,10HINPUT RAYS,                                       &
     &       //13X,2HY=,1PE13.6,2X,2HU=,1PE13.6,                        &
     &       2X,3HCY=,1PE13.6,2X,3HCU=,1PE13.6,2X,                      &
     &       14HTARGET HEIGHT=,1PE13.6)                                 
   30 FORMAT (//T3,7HSURFACE,T15,7HAXIAL Y,T33,7HAXIAL U,T51,7HCHIEF Y, &
     &       T69,7HCHIEF U,/)                                           
   40 FORMAT (T6,I2,T10,1PE16.8,T28,1PE16.8,T46,1PE16.8,T64,1PE16.8) 
   60 FORMAT (/13X,25HIMAGE DISTANCE AXIAL RAY=,1PE16.8,6X,             &
     &       15HMAGNIFICATION= ,1PE16.8,/13X,                           &
     &       25HIMAGE DISTANCE CHIEF RAY=,1PE16.8,6X,                   &
     &       14HFOCAL LENGTH =,1PE17.8,6X,10HF NUMBER =,1PE16.8//)      
   70 FORMAT (25H ABERRATIONS    SPHERICAL,9X,4HCOMA,10X,11HASTIGMATISM,&
     &       5X,10HDISTORTION,7X,7HPETZVAL,15X,9HCHROMATIC/8H SURFACE,  &
     &       89X,7HA COEF ,9X,7HB COEF )                                
   80 FORMAT(4X,I2,6X,1P7E16.8) 
   90 FORMAT(10X,'APERTURE STOP SURFACE MUST BE GREATER THAN 1; = ',I2/) 
  100 FORMAT(4X,I2,4H ASP,2X,1P4E16.8) 
  110 FORMAT (/9H TOTAL SP,3X,1P7E16.8) 
  120 FORMAT (10H TOTAL ASP,2X,1P4E16.8) 
  130 FORMAT (6H TOTAL,6X,1P7E16.8) 
  140 FORMAT (/13X,39HVECTOR SUM OF THIRD ORDER ABERRATIONS =,1PE16.8) 
  150 FORMAT (/13X,39HVECTOR SUM OF FIRST ORDER ABERRATIONS =,1PE16.8/) 
  170 FORMAT (//32H CONIC CONSTANT OF PRIMARY   =  ,1PE16.9) 
  180 FORMAT (32H CONIC CONSTANT OF SECONDARY =  ,1PE16.9///) 
  220 FORMAT (18H FOCAL POINTS AT  ,1PE16.9,9H AND AT  ,1PE16.9/) 
  230 FORMAT (/13X,8HSURFACE ,I2,18H IS A TORIC.RX(I)=,1PE16.8) 
  240 FORMAT (/13X,8HSURFACE ,I2,33H IS A CYLINDER WITH AXIS ALONG X.,  &
     &       6HR(I)= ,1PE16.8)                                          
  250 FORMAT (/13X,8HSURFACE ,I2,33H IS A CYLINDER WITH AXIS ALONG Y.,  &
     &       7HRX(I)= ,1PE16.8)                                         
  260 FORMAT(10X,'APERTURE STOP SIZE NOT SPECIFIED, WILL NOT LOCATE ',  &
     &       'PARAXIAL PUPILS'/)                                        
  270 FORMAT(10X,'FOUND TILTED/DECENTERED SURFACE WHILE ',              &
     &       'LOCATING PUPILS'/)                                        
  271 FORMAT(10X,'T(',I2,') = ',1PE16.8,' FROM THICKNESS SOLVE'/) 
  272 FORMAT(10X,'T(',I2,') = ',1PE16.8,' FROM CLEAR APERTURE SOLVE'/) 
  273 FORMAT(10X,'C(',I2,') = ',1PE16.8,' FROM CURVATURE SOLVE'/) 
!                                                                       
      KK=0 
      HYMAX=HYINIT+(RNOBJ-1.)*HYDEL 
      HYMAX=-DTAN(HYMAX*3.14159265358979D0/180.)*(S-D) 
      HINIT=-DTAN(HYINIT*3.14159265358979D0/180.)*(S-D) 
      IF ( DABS(HINIT).GT.DABS(HYMAX) ) HYMAX = HINIT 
!                                                                       
!                BEGIN CALCULATING POSITION OF EP                       
!                I = SURF NUMBER OF APERTURE STOP                       
!                                                                       
      I=APSTOP/30 
!                                                                       
!                NONE SPECIFIED                                         
!                                                                       
      IF (I.LE.1.AND.(IPRINT/2)*2.EQ.IPRINT) WRITE(6,90) I 
      IF (I.LE.1) RETURN 
      J=I 
      ISTOP=J 
      KK=I 
      ISW=1 
      L=J-1 
      LAMDA = 2 
      IF ( NCOL.EQ.1 ) LAMDA = ICOL(1) 
      IF ( (DABS(FXY)+DABS(FXNU)+DABS(FBY)+DABS(FBNU)) .EQ. ZERO )      &
     &   GO TO 275                                                      
      FXNU=RHO/(S-D) 
      FXY=FXNU*(S-D) 
      FBY=ZERO 
      FBNU = -HYMAX/(S-D) 
  275 XY(1)=FXY 
      U(1)=FXNU 
      BY(1)=FBY 
      BU(1)=FBNU 
      GO TO 370 
!                  TRACE UP TO APERTURE STOP IN CASE THERE ARE          
!                  CURVATURE/THICKNESS SOLVES TO BE DETERMINED          
!                FIND SIZE AND LOCATION OF ENTRANCE PUPIL               
  280 IF ( FMASK(ISTOP).GT.ZERO ) GO TO 290 
      IF ((IPRINT/2)*2.EQ.IPRINT) WRITE(6,260) 
      GO TO 682 
  290 AB=ZERO 
      ABB=FMASK(ISTOP) 
      SGN1 = FREF0 
      L = ISTOP-1 
      DO 295 I=1,L 
        IF (FREF(I).EQ.-ONE) SGN1 = -SGN1 
  295 END DO 
      SGNAP = SGN1*FREF(ISTOP) 
      IF (ISTOP.EQ.2) GO TO 300 
      L=ISTOP-2 
      TM1=.1D0 
      TM1B=TM1 
      GO TO 310 
  300 TM1=.1D0 
      TM1B=TM1 
      SGN1 = FREF0 
      IF (FREF(1).EQ.-ONE) SGN1 = -SGN1 
      GO TO 330 
  310 DO 320 I=1,L 
        M=ISTOP-I 
        IF ((IPRINT/2)*2.NE.IPRINT) GO TO 315 
        IF (TILTX(I).NE.ZERO .OR. TILTY(I).NE.ZERO .OR. TILTZ(I).NE.ZERO&
     &     .OR. XDISP(I).NE.ZERO .OR. YDISP(I).NE.ZERO ) WRITE(6,270)   
  315   SGN2 = SGN1 
        IF (FREF(M).EQ.-ONE) SGN2 = -SGN1 
        RN1  = FN(M-1,LAMDA)*SGN2 
        RN2  = FN(M,LAMDA)*SGN1 
        GNU  = RN2/RN1 
        ABB  = ABB + T(M)*TM1B 
        AB   = AB + T(M)*TM1 
        TM1B = GNU*TM1B + ABB*C(M)*(ONE-GNU) 
        TM1  = GNU*TM1 + AB*C(M)*(ONE-GNU) 
        SGN1 = SGN2 
  320 END DO 
  330 TM7=-AB/TM1 
      RHO=DABS(ABB+TM7*TM1B) 
      T(1)=TM7 
      D = ZERO 
      ENPUPR = RHO 
      ENPUPL = D 
      FXY=RHO 
      FXNU=RHO/(S-D) 
      FBY=ZERO 
      FBNU = -HYMAX/(S-D) 
!                                                                       
!             LOCATE AND SIZE PARAXIAL EXIT PUPIL; PRINT INFORMATION    
      AB=ZERO 
      ABB=FMASK(ISTOP) 
      L=ISTOP+1 
      SGN1 = SGNAP 
      TM1=.1D0 
      TM1B=TM1 
      DO 340 I=L,NSURF 
        IF ((IPRINT/2)*2.NE.IPRINT) GO TO 335 
        IF (TILTX(I).NE.ZERO .OR. TILTY(I).NE.ZERO .OR. TILTZ(I).NE.ZERO&
     &      .OR. XDISP(I).NE.ZERO .OR. YDISP(I).NE.ZERO ) WRITE(6,270)  
  335   SGN2 = SGN1 
        IF (FREF(I).EQ.-ONE) SGN2 = -SGN1 
        RN1  = FN(I,LAMDA)*SGN2 
        RN2  = FN(I-1,LAMDA)*SGN1 
        GNU  = RN2/RN1 
        ABB  = ABB + T(I-1)*TM1B 
        AB   = AB + T(I-1)*TM1 
        TM1B = GNU*TM1B + ABB*C(I)*(GNU-ONE) 
        TM1  = GNU*TM1 + AB*C(I)*(GNU-ONE) 
        SGN1 = SGN2 
  340 END DO 
      TM7=-AB/TM1 
      TM8=DABS(ABB+TM7*TM1B) 
      EXPUPR = TM8 
      EXPUPL = TM7 
!                END OF EXIT PUPIL LOCATION                             
!                                                                       
!             IF FOLLOWING 'LEPRT' COMMAND, RETURN TO MAIN ROUTINE      
!                                                                       
  350 IF ((IPRINT/2)*2.EQ.IPRINT) GO TO 355 
      CALL PREPRT 
      RETURN 
!                                                                       
!                                                                       
!             SET SURFACE PRINT SWITCH                                  
!                                                                       
  355 IPRINA=0 
      IPRINB=0 
      IF ((IPRINT/4)*2.NE.IPRINT/2) IPRINA=1 
      LAMDA=ICOL(1) 
!                                                                       
      ISW=2 
!                                                                       
!                TRACE NSURF NUMBER OF SURFACES FOR NCOL COLORS         
!                                                                       
!                                                                       
!        COLORS LOOP - ENTER                                            
!                                                                       
  360 DO 680 K=1,NCOL 
        L=NSURF 
        XY(1)=FXY 
        BY(1)=FBY 
        U(1)=FXNU 
        BU(1) = FBNU 
!             SET COLOR NUMBER                                          
        LAMDA=ICOL(K) 
!                                                                       
!             SGN1, SGN2 ARE USED TO DETERMINE SIGNS OF INDICES         
  370   SGN1=FREF0 
!             START OF SURFACE-BY-SURFACE PARAXIAL RAY TRACE            
        DO 450 I=1,L 
          SGN2 = SGN1 
          IF (FREF(I).EQ.-ONE) SGN2 = -SGN1 
          TM1 = SGN2*FN(I,LAMDA) 
          IF (I.EQ.1) TM2 = SGN1*OBJN(LAMDA) 
          IF (I.GT.1) TM2 = SGN1*FN(I-1,LAMDA) 
          TM3 = C(I) 
          GNU = TM2/TM1 
          IORDN=ORDN(I,LAMDA)+.5 
          IORDA=0 
!                                                                       
          IF (I.EQ.1) TM1A = TM1 
          IF (I.EQ.2) TM1B = TM1 
          IF (I.EQ.3) TM1C = TM1 
          IF (I.EQ.2) GNU2 = GNU 
          IF (I.EQ.3) GNU3 = GNU 
!                IF NORDN(I) AND 1, GO TO CURVATURE SOLVE               
!                                                                       
  390     IF ((IORDN/2)*2.NE.IORDN) GO TO 400 
          CALL SURTYP (I) 
          X(I)=TM3*(GNU-ONE) 
!                                                                       
!                U IS ANGLE OF AXIAL RAY                                
!                                                                       
          U(I+1) = GNU*U(I) + XY(I)*X(I) 
          GO TO 410 
!                                                                       
!                CURVATURE SOLVE ROUTINE                                
!                                                                       
!                THE USER SPECIFIES THE ANGLE AT WHICH THE RAY LEAVES   
!                THE SURFACE (SXNU), AND THE NECESSARY CURVATURE        
!                OF THE SURFACE (C(I)) IS CALCULATED                    
!                                                                       
  400     C(I)=(SXNU(I+1)-GNU*U(I))/(XY(I)*(GNU-ONE)) 
          IF ((IPRINT/2)*2.EQ.IPRINT) WRITE(6,273) I,C(I) 
          TM3=C(I) 
          X(I)=TM3*(GNU-ONE) 
          R(I)=1./C(I) 
          CALL SURTYP (I) 
          U(I+1) = GNU*U(I) + XY(I)*X(I) 
!                                                                       
!                IF AND(NORDN(I),2), GO TO THICKNESS SOLVE ROUTINE      
!                                                                       
  410     IF ((IORDN/4)*2.EQ.IORDN/2) GO TO 420 
!                                                                       
!                THICKNESS SOLVE (OR HEIGHT SOLVE) ROUTINE              
!                THE USER SPECIFIES THE HEIGHT OF THE RAY ON THE NEXT   
!                SURFACE (SXY(I+1)),CALCULATE THE THICKNESS-- DIST      
!                FROM SURFACE I TO SURFACE I+1 (T(I))                   
!                                                                       
          T(I) = (SXY(I+1)-XY(I))/U(I+1) 
          IF ((IPRINT/2)*2.EQ.IPRINT) WRITE(6,271) I,T(I) 
          IORDA=1 
!                                                                       
!                IF AND(NORDN(I),4) GO TO CLEAR APERTURE ROUTINE        
!                                                                       
  420     IF ((IORDN/8)*2.EQ.IORDN/4) GO TO 430 
!                                                                       
!                CLEAR APERTURE ROUTINE                                 
!                THE USER SPECIFIES Y SUB ZERO AND THE TO-THE-EDGE      
!                THICKNESS, CALCULATE T(I)                              
!                THAT IS, ADD THE DIFFERENCE OF SAG AT SURF I AND       
!                SAG AT I+1 (AT REFERENCE HEIGHT Y0) TO THE MINIMUM     
!                CLEAR APERTURE SPECIFIED BY USER - THIS MEANS THAT     
!                THE DIST BETWEEN SURF I AND I+1 IS A MINIMUM OF        
!                RDSPAC(I) FROM HEIGHT ZERO (VERTICES) THRU Y SUB ZERO  
!                                                                       
          ZA=C(I)*Y0(I)*Y0(I)/(1.+DSQRT(1.-C(I)*C(I)*Y0(I)*Y0(I))) 
          ZB=C(I+1)*Y0(I)*Y0(I)/(1.+DSQRT(1.-C(I+1)*C(I+1)*Y0(I)*Y0(I))) 
          T(I)=RDSPAC(I) 
          IF (ZA.GE.ZB) T(I)=T(I)+ZA-ZB 
          IF ((IPRINT/2)*2.EQ.IPRINT) WRITE(6,272) I,T(I) 
!                                                                       
!                CALCULATE HEIGHT OF AXIAL RAY, XY                      
!                                                                       
  430     XY(I+1) = XY(I) + U(I+1)*T(I) 
!                                                                       
!                CALCULATE ANGLE OF CHIEF RAY, BU                       
!                                                                       
          BU(I+1) = GNU*BU(I) + BY(I)*X(I) 
!                                                                       
!                CALCULATE HEIGHT OF CHIEF RAY, BY                      
!                                                                       
          BY(I+1) = BY(I) + BU(I+1)*T(I) 
          SGN1 = SGN2 
  450   CONTINUE 
        GO TO (280,460), ISW 
  460   IF ( DABS(U(NSURF)).LT.EPS ) GO TO 470 
!                                                                       
!                FTN IS FINAL THICKNESS                                 
!                                                                       
        FTN = ( FXYJ-XY(NSURF-1) )/U(NSURF) 
!                                                                       
!                RMAG IS LATERAL MAGNIFICATION                          
!                                                                       
        RMAG=U(1)/U(NSURF) 
        GO TO 480 
  470   FTN=1.E20 
        RMAG=1.E20 
  480   IF ( DABS(BU(NSURF)).LT.EPS ) GO TO 490 
!                                                                       
!                FTNB IS IMAGE DIST PRINCIPAL RAY                       
!                                                                       
        FTNB = -BY(NSURF-1)/BU(NSURF) 
        GO TO 500 
  490   FTNB=1.E20 
!                                                                       
!                PHI IS OPTICAL INVARIANT                               
!                                                                       
  500   PHI=TM1A*(BY(1)*U(1)-XY(1)*BU(1)) 
!                                                                       
!                COMPUTE ABERRATION COEFFICIENTS                        
!                                                                       
!             B(I)   -   SPHERICAL ABERRATION                           
!             CC(I)  -   ASTIGMATISM                                    
!             F(I)   -   COMA                                           
!             E(I)   -   DISTORTION                                     
!             P(I)   -   PETZVAL                                        
!             ACH(I) -   AXIAL CHROMATIC ABERRATION                     
!             BCH(I) -   LATERAL CHROMATIC ABERRATION                   
!                                                                       
!                    OTHER                                              
!                                                                       
!             TM2    -   SGN1*N                                         
!             TM1    -   SGN2*N'                                        
!                                                                       
        SP    = ZERO 
        DN1   = ZERO 
        DN2   = ZERO 
        CO    = ZERO 
        AS    = ZERO 
        DIS   = ZERO 
        PET   = ZERO 
        TACH  = ZERO 
        TCH   = ZERO 
        CSPH  = ZERO 
        CCOM  = ZERO 
        CAST  = ZERO 
        CDIST = ZERO 
        L=NSURF-1 
        SGN1 = FREF0 
!                                                                       
        DO 590 I = 1, L 
          SGN2 = SGN1 
          IF (FREF(I).EQ.-ONE) SGN2 = -SGN1 
          TM1=SGN2*FN(I,LAMDA) 
          IF (I.EQ.1) TM2 = SGN1*OBJN(LAMDA) 
          IF (I.GT.1) TM2 = SGN1*FN(I-1,LAMDA) 
          IF ( RR(I).EQ.ZERO ) TM3 = ZERO 
          IF ( RR(I).NE.ZERO ) TM3=1./RR(I) 
          GNU = TM2/TM1 
          EN=XY(I)*TM3+U(I) 
          SS=XY(I)*TM2*(GNU-ONE)*(U(I+1)+EN) 
          B(I)=SS*EN*EN 
          BEN=BY(I)*TM3+BU(I) 
          F(I)=SS*EN*BEN 
          CC(I)=SS*BEN*BEN 
          P(I)=C(I)*(GNU-ONE)/TM2 
          BS=BY(I)*TM2*(GNU-ONE)*(BU(I+1)+BEN) 
          E(I)=BS*EN*BEN+PHI*(BU(I)*BU(I)-BU(I+1)*BU(I+1)) 
          DFCN = ZERO 
          IF (NCOL.EQ.1) GO TO 560 
              IA=I-1 
              IF (I.NE.1) GO TO 550 
                 DN1=SGN1*((OBJN(3)-OBJN(1))/OBJN(2)) 
              GO TO 555 
!                                                                       
  550         DN1=SGN1*((FN(IA,3)-FN(IA,1))/FN(IA,2)) 
  555         DN2=SGN2*((FN(I,3)-FN(I,1))/FN(I,2)) 
  560     DFCN = DN1-DN2 
!                                                                       
          ACH(I)=XY(I) * TM2 * EN*DFCN 
          BCH(I)=XY(I) * TM2 * BEN*DFCN 
!                                                                       
!                ACCUMULATE TOTAL ABERRATIONS                           
!                                                                       
          SP=B(I)+SP 
          CO=F(I)+CO 
          AS=CC(I)+AS 
          DIS=E(I)+DIS 
          PET=P(I)+PET 
          TACH=TACH-ACH(I) 
          TCH=TCH-BCH(I) 
          YSQ=XY(I)*XY(I) 
!                                                                       
!                FIND ASPHERICAL COEFFICIENTS                           
!                                                                       
          BAS(I)=8.0*TM1*(GNU-ONE)*YSQ*YSQ*DUM(I) 
          FAS(I)=BAS(I)*BY(I)/XY(I) 
          IF (BAS(I).EQ.ZERO) GO TO 570 
          CAS(I)=FAS(I)*FAS(I)/BAS(I) 
          GO TO 580 
  570     CAS(I)=ZERO 
  580     EAS(I)=CAS(I)*BY(I)/XY(I) 
          CSPH=BAS(I)+CSPH 
          CCOM=FAS(I)+CCOM 
          CAST=CAS(I)+CAST 
          CDIST=EAS(I)+CDIST 
          SGN1 = SGN2 
  590   CONTINUE 
!                                                                       
! NSURF = NUMBER OF LAST SURFACE WHICH IS A DUMMY SURFACE               
!            L = NUMBER OF LAST REAL SURFACE                            
!            LAMDA = # OF THE COLOR (SET NEAR LINE LABELED 355 OR 360)  
!                                                                       
        TACH = TACH/( SGN1 * U(L) * FN(L, LAMDA) ) 
        TCH  = TCH /( SGN1 * U(L) * FN(L, LAMDA) ) 
!                                                                       
!                TAKE SUMS OF ABERRATIONS                               
!                                                                       
        ABB=CSPH+SP 
        ABF=CCOM+CO 
        ABC=CAST+AS 
        ABE=CDIST+DIS 
        VAL=ABB*ABB+ABF*ABF+ABC*ABC+ABE*ABE 
        VALUE=DSQRT(VAL) 
        VAL1=DSQRT(PET*PET+TACH*TACH+TCH*TCH) 
!                                                                       
!             PRINT THE WHOLE MESS AT ONE WHACK                         
!                                                                       
        IF (IPRINB.EQ.1 ) GO TO 655 
        IF (IPRINA.EQ.0 ) GO TO 610 
!                                                                       
        CALL HEADIN(LAMDA) 
!                                                                       
        WRITE(6,10)  XY(1),U(1),BY(1),BU(1),FXYJ 
        WRITE(6,30) 
        L=NSURF 
        DO 600 I=1,L 
          IF ( I.NE.L ) WRITE(6,40) I,XY(I),U(I+1),BY(I),BU(I+1) 
          IF ( I.EQ.L ) WRITE(6,40) I,XY(I),U(I),BY(I),BU(I) 
  600   CONTINUE 
!                                                                       
!             PRINT SURFACE TYPE IF TORIC/CYLINDER                      
!                                                                       
        IF (ISURFX.EQ.2) WRITE(6,230) ISURXN,RSAVE 
        IF (ISURFX.EQ.3) WRITE(6,240) ISURXN,RSAVE 
        IF (ISURFX.EQ.4) WRITE(6,250) ISURXN,RSAVE 
  610   IF ( DABS(U(1)).LT.EPS ) GO TO 630 
        U(1)=ZERO 
        XY(1)=RHO 
        SGN1 = FREF0 
        DO 620 I=1,L 
          SGN2 = SGN1 
          IF (FREF(I).EQ.-ONE) SGN2 = -SGN1 
          TM1=SGN2*FN(I,LAMDA) 
          IF (I.EQ.1) TM2 = SGN1*OBJN(LAMDA) 
          IF (I.GT.1) TM2 = SGN1*FN(I-1,LAMDA) 
          TM3=C(I) 
          CALL SURTYP (I) 
          GNU = TM2/TM1 
          XI = TM3 * (GNU-ONE) 
          U(I+1)=GNU*U(I)+XY(I)*XI 
          XY(I+1)=XY(I)+U(I+1)*T(I) 
          SGN1 = SGN2 
  620   CONTINUE 
  630   FOCL=-XY(1)/U(NSURF) 
        FNUM=DABS(.5/U(NSURF)) 
        WRITE(6,60) FTN,RMAG,FTNB,FOCL,FNUM 
        L = NSURF-1 
  655   WRITE(6,70) 
        DO 660 I=1,L 
          WRITE(6,80) I,B(I),F(I),CC(I),E(I),P(I),ACH(I),BCH(I) 
          IF (RCON(I).EQ.ZERO .AND. COEF(I,1).EQ.ZERO) GO TO 660 
          WRITE(6,100) I,BAS(I),FAS(I),CAS(I),EAS(I) 
  660   CONTINUE 
        WRITE(6,110) SP,CO,AS,DIS,PET,TACH,TCH 
        WRITE(6,120) CSPH,CCOM,CAST,CDIST 
        WRITE(6,130) ABB,ABF,ABC,ABE,PET,TACH,TCH 
        WRITE(6,140) VALUE 
        WRITE(6,150) VAL1 
  680 END DO 
!                                                                       
!       COLORS LOOP - EXIT                                              
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC           
!       PRINCIPLE RETURN                                                
!                                                                       
  682 RETURN 
!                                                                       
!***********************************************************            
!                                                                       
!                                                                       
!                RITCHEY-CHRETIEN- TELESCOPE DESIGN ROUTINE             
!                                                                       
      ENTRY RCDES 
      FP = -R(2)/TWO 
      FSY = -XY(1)/U(NSURF) 
      BF  = T(2) + (FXYJ-XY(3))/U(NSURF) 
      A   = FSY/FP 
      V = TWO*(FP+BF)/(A*A*(FSY-BF)) 
      CONIC(2) = -(ONE+V) 
      V = (TWO*FSY*(A+ONE))/((A-ONE)*(A-ONE)*(A-ONE)*(FSY-BF)) 
      V   = V + (4.D0*A/((A-ONE)*(A-ONE))) 
      CONIC(3) = -(ONE+V) 
      WRITE(6,170) CONIC(2) 
      WRITE(6,180) CONIC(3) 
      IPRINB=1 
      GO TO 370 
!                                                                       
!                CASSEGRAIN TELESCOPE DESIGN ROUTINE                    
!                                                                       
      ENTRY CADES 
      CONIC(2)=-ONE 
      BF=(FXYJ-XY(3))/U(4) 
      FP=-R(2)/2. 
      V=(BF+FP+T(2))/(BF-FP-T(2)) 
      CONIC(3)=-V*V 
      T(3)=BF 
      WRITE(6,170) CONIC(2) 
      WRITE(6,180) CONIC(3) 
      IPRINB=1 
      GO TO 370 
!                                                                       
!                DAHL-KIRKHAM TELESCOPE DESIGN ROUTINE                  
!                                                                       
      ENTRY DKDES 
      CONIC(3)=ZERO 
      BF=(FXYJ-XY(3))/U(4) + T(2) 
      PF=-R(2)/2. 
      A=FOCL/PF 
      V=(PF+BF)/(A+1.) 
      PP=PF+BF-V 
      RS=-(2.*V*PP)/(PP-V) 
      CONIC(2)=-(1.-4.*V*V*(PP+V)*(PP+V)/(RS*R(2)*PP*PP)) 
      WRITE(6,170) CONIC(2) 
      WRITE(6,180) CONIC(3) 
      ECC = CONIC(2) 
      F1=R(2)/(1.-ECC) 
      F2=R(2)/(1.+ECC) 
      WRITE(6,220) F1,F2 
      IPRINB=1 
      GO TO 370 
      END                                           
!                                                                       
!*******************************************************                
      SUBROUTINE PREPRT 
!*******************************************************                
!             PRINTS OPTICAL PRESCRIPTION DATA                          
      IMPLICIT REAL*8 (A-H,O-Z) 
!                                                                       
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /PUPIL/  ENPUPR,ENPUPL,EXPUPR,EXPUPL 
      COMMON /HEAD/   LINES,IPAGE,NSYS,LAMDA,NAME(160) 
!                                                                       
      DIMENSION XT(300),YT(300),IFLAG(40),XTH(300),YTH(300),IORDER(3) 
!                                                                       
      DATA PI/3.141592653589793D0/,ONE/1.D0/,TWO/2.D0/,THREE/3.D0/ 
      DATA FOUR/4.D0/,FIVE/5.D0/,ZERO/0.D0/,CLEAR/'CLEAR'/,OBSC/'OBSC'/ 
      DATA TYPE/'        '/ 
      DATA IZERO/0/,IONE/1/,ITWO/2/ 
!                                                                       
!             PRINT HEADING                                             
      CALL HEADIN(0) 
!                                                                       
!             PRINT 'SURFACE DATA' CARD                                 
      WRITE (6,200) 
!                                                                       
!             PRINT OBJECT DATA                                         
      WRITE (6,210) 
      NOBJ=RNOBJ 
      DO 10 I=1,NOBJ 
        IF ( NOBJ.EQ.0 ) GO TO 10 
        RI=I 
        XT(I)=HXINIT+(RI-ONE)*HXDEL 
        YT(I)=HYINIT+(RI-ONE)*HYDEL 
        XTH(I)=-DTAN(XT(I)*PI/180.D0)*(S-D) 
        YTH(I)=-DTAN(YT(I)*PI/180.D0)*(S-D) 
   10 END DO 
      WRITE (6,260) S 
      IF ( OBJN(2).EQ.ONE .AND. FREF0.EQ.ONE ) WRITE(6,262) 
      IF ( OBJN(2).EQ.ONE .AND. FREF0.EQ.-ONE) WRITE(6,264) 
      IF ( OBJN(2).NE.ONE .AND. FREF0.EQ.ONE ) WRITE(6,266) OBJN(2) 
      JOBJ=NOBJ 
      KADD=0 
      NPR1=6 
      IF (NOBJ.LT.6) NPR1=NOBJ 
      IF ( NOBJ.NE.0 ) WRITE (6,220) (XT(J),J=1,NPR1) 
   11 JOBJ=JOBJ-6 
      KADD=KADD+6 
      IF ( JOBJ.LE.0 ) GO TO 12 
      NPR2=6 
      IF (JOBJ.LT.6) NPR2=JOBJ 
      WRITE (6,222) (XT(J+KADD),J=1,NPR2) 
      GO TO 11 
   12 JOBJ=NOBJ 
      KADD=0 
      IF ( NOBJ.NE.0 ) WRITE (6,230) (YT(J),J=1,NPR1) 
   13 JOBJ=JOBJ-6 
      KADD=KADD+6 
      IF (JOBJ.LE.0) GO TO 14 
      NPR2=6 
      IF (JOBJ.LT.6) NPR2=JOBJ 
      WRITE(6,222) (YT(J+KADD),J=1,NPR2) 
      GO TO 13 
   14 JOBJ=NOBJ 
      KADD=0 
      IF ( NOBJ.NE.0 ) WRITE (6,240) (XTH(J),J=1,NPR1) 
   15 JOBJ=JOBJ-6 
      KADD=KADD+6 
      IF (JOBJ.LE.0) GO TO 16 
      NPR2=6 
      IF (JOBJ.LT.6) NPR2=JOBJ 
      WRITE(6,222) (XTH(J+KADD),J=1,NPR2) 
      GO TO 15 
   16 JOBJ=NOBJ 
      KADD=0 
      IF ( NOBJ.NE.0 ) WRITE (6,250) (YTH(J),J=1,NPR1) 
   17 JOBJ=JOBJ-6 
      KADD=KADD+6 
      IF (JOBJ.LE.0) GO TO 18 
      NPR2=6 
      IF (JOBJ.LT.6) NPR2=JOBJ 
      WRITE(6,222) (YTH(J+KADD),J=1,NPR2) 
      GO TO 17 
   18 IF ( NOBJ.EQ.0 ) WRITE(6,220) 
      IF ( NOBJ.EQ.0 ) WRITE(6,230) 
      IF ( NOBJ.EQ.0 ) WRITE(6,240) 
      IF ( NOBJ.EQ.0 ) WRITE(6,250) 
!                                                                       
!                                                                       
!             WRITE IMAGE SURFACE DATA                                  
      IF (RADIMG.NE.ZERO) WRITE(6,270) NPLANE,DELIMP,CVIMG,RADIMG,CONIMG 
      IF (RADIMG.EQ.ZERO) WRITE(6,275) NPLANE,DELIMP,CVIMG,CONIMG 
!                                                                       
!             WRITE SYSTEM UNITS                                        
      WRITE(6,277) 
      IF (UFLAG.EQ.ONE)   WRITE (6,280) 
      IF (UFLAG.EQ.TWO)   WRITE (6,290) 
      IF (UFLAG.EQ.THREE) WRITE (6,300) 
      IF (UFLAG.EQ.FOUR)  WRITE (6,310) 
      IF (UFLAG.EQ.FIVE)  WRITE (6,320) 
!             IF UFLAG NONE OF ABOVE, UNKNOWN UNITS                     
      IF ( UFLAG.LT.ONE .OR. UFLAG.GT.FIVE ) WRITE (6,330) UFLAG 
!                                                                       
!             WRITE ENTRANCE PUPIL DATA                                 
      DIAM=TWO*ENPUPR 
      WRITE (6,340) DIAM,ENPUPL 
!                                                                       
!             WRITE EXIT PUPIL DATA                                     
      DIAM=TWO*EXPUPR 
      WRITE (6,350) DIAM,EXPUPL 
      IF ( APSTOP.NE.ZERO ) ISTOP = APSTOP/30.D0 
      IF ( APSTOP.NE.ZERO ) WRITE(6,353) ISTOP 
      IF ( APSTOP.EQ.ZERO ) WRITE(6,355) 
      IF ( IREF.LE.0 ) WRITE(6,357) 
      IF ( IREF.GT.0 ) WRITE(6,358) IREF 
!                                                                       
!             WRITE BASIC SURFACE DATA                                  
      WRITE (6,360) 
      DO 20 I=1,NSURF 
        IF (R(I).EQ.ZERO .AND. FREF(I).EQ.ONE .AND. FN(I,2).EQ.ONE )    &
     &       WRITE(6,370) I,C(I),CONIC(I),T(I)                          
        IF ( R(I).NE.ZERO .AND. FREF(I).EQ.ONE .AND. FN(I,2).EQ.ONE )   &
     &       WRITE(6,372) I,C(I),R(I),CONIC(I),T(I)                     
        IF ( R(I).EQ.ZERO .AND. FREF(I).EQ.-ONE .AND. FN(I,2).EQ.ONE )  &
     &       WRITE(6,374) I,C(I),CONIC(I),T(I)                          
        IF ( R(I).NE.ZERO .AND. FREF(I).EQ.-ONE .AND. FN(I,2).EQ.ONE )  &
     &       WRITE(6,376) I,C(I),R(I),CONIC(I),T(I)                     
        IF ( R(I).EQ.ZERO .AND. FN(I,2).NE.ONE )                        &
     &       WRITE(6,377) I,C(I),CONIC(I),T(I),(FN(I,J),J=1,3)          
        IF ( R(I).NE.ZERO .AND. FN(I,2).NE.ONE )                        &
     &       WRITE(6,378) I,C(I),R(I),CONIC(I),T(I),(FN(I,J),J=1,3)     
   20 END DO 
!                                                                       
!             WRITE ACONIC SURFACE DATA, IF ANY                         
      ICONIC=0 
      DO 40 I=1,NSURF 
        IFLAG(I)=0 
        DO 30 J=1,4 
          IF (COEF(I,J).EQ.ZERO) GO TO 30 
          IFLAG(I)=1 
          ICONIC=1 
          GO TO 40 
   30   CONTINUE 
   40 END DO 
!                                                                       
      IF (ICONIC.EQ.0) GO TO 60 
      WRITE (6,380) 
      DO 50 I=1,NSURF 
        IF (IFLAG(I).EQ.0) GO TO 50 
        WRITE (6,390) I,(COEF(I,J),J=1,4) 
   50 END DO 
!                                                                       
!             WRITE TILT/DISPLACEMENT DATA, IF ANY                      
   60 ITILT=0 
      DO 70 I=1,NSURF 
        IFLAG(I)=0 
        IF (XDISP(I).EQ.ZERO.AND.YDISP(I).EQ.ZERO.AND.TILTX(I).EQ.ZERO.A&
     &  ND.TILTY(I).EQ.ZERO.AND.TILTZ(I).EQ.ZERO) GO TO 70              
        IFLAG(I)=1 
        ITILT=1 
   70 END DO 
!                                                                       
      IF (ITILT.EQ.0) GO TO 90 
      WRITE (6,400) 
      DO 80 I=1,NSURF 
        IF (IFLAG(I).EQ.0) GO TO 80 
        WRITE (6,410) I,XDISP(I),YDISP(I),TILTX(I),TILTY(I),TILTZ(I) 
        IF (DABS(FAKEB(I)).EQ.ONE) WRITE (6,420) 
        IF (DABS(FAKEB(I)).EQ.TWO) WRITE (6,430) 
        IF (DABS(FAKEB(I)).EQ.THREE) WRITE (6,440) 
   80 END DO 
!                                                                       
!             WRITE GRATING DATA, IF ANY                                
   90 IGRAT=0 
      DO 100 I=1,NSURF 
        IFLAG(I)=0 
        IF (RDSPAC(I).EQ.ZERO) GO TO 100 
        IFLAG(I)=1 
        IGRAT=1 
  100 END DO 
!                                                                       
      IF (IGRAT.EQ.0) GO TO 120 
      WRITE (6,450) 
      DO 110 I=1,NSURF 
        IF (IFLAG(I).EQ.0) GO TO 110 
        SPACE=DABS(RDSPAC(I)) 
        DO 105 K=1,3 
        IORDER(K)=ORDN(I,K) 
        IF (RDSPAC(I).LT.ZERO)                                          &
     &     WRITE(6,460) I,WAVL(K),IORDER(K),SPACE                       
        IF (RDSPAC(I).GT.ZERO)                                          &
     &     WRITE(6,470) I,WAVL(K),IORDER(K),SPACE                       
  105 END DO 
  110 END DO 
!                                                                       
!             WRITE TORIC DATA, IF ANY                                  
  120 ITOR=0 
      DO 130 I=1,NSURF 
        IFLAG(I)=0 
        TOR=Y0(I) 
        IF (TOR.NE.ONE) GO TO 130 
        IFLAG(I)=1 
        ITOR=1 
  130 END DO 
!                                                                       
      IF (ITOR.EQ.0) GO TO 150 
      WRITE (6,480) 
      DO 140 I=1,NSURF 
        IF (IFLAG(I).EQ.0) GO TO 140 
        IF ( RX(I).NE.ZERO ) WRITE (6,490) I,CX(I),RX(I) 
         IF ( RX(I).EQ.ZERO ) WRITE (6,495) I,CX(I) 
  140 END DO 
!                                                                       
!             WRITE MASK DATA, IF ANY                                   
  150 IMASKC=0 
      IMASKR=0 
      IMASKE=0 
      DO 160 I=1,NSURF 
        IFLAG(I)=0 
        IF (FMASK(I).EQ.ZERO) GO TO 160 
        IFLAG(I)=1 
        IF (FAKEC(I).GT.0) IMASKR=1 
        IF (FAKEC(I).EQ.0) IMASKC=1 
        IF (FAKEC(I).LT.0) IMASKE=1 
        IMASK=1 
  160 END DO 
!                                                                       
      IF (IMASKC+IMASKE+IMASKR.EQ.0) GO TO 190 
      WRITE(6,499) 
      IF (IMASKC.EQ.0) GO TO 182 
      WRITE (6,500) 
      DO 180 I=1,NSURF 
        IF (IFLAG(I).EQ.0) GO TO 180 
        IF ( FAKEC(I).NE.0) GO TO 180 
        IF (FMASK(I).LT.ZERO) TYPE=OBSC 
        IF (FMASK(I).GT.ZERO) TYPE=CLEAR 
        RMASK=FMASK(I) 
        XCEN=XMN(I) 
        YCEN=YMN(I) 
        XMIN=ZERO 
        XMAX=ZERO 
        YMIN=ZERO 
        YMAX=ZERO 
  170   WRITE (6,510) I,TYPE,RMASK,XCEN,YCEN 
  180 END DO 
!                                                                       
  182 IF (IMASKR.EQ.0) GO TO 186 
      WRITE(6,501) 
      DO 185 I=1,NSURF 
        IF (IFLAG(I).EQ.0) GO TO 185 
        IF (FAKEC(I).NE.1) GO TO 185 
        IF (FMASK(I).LT.0) TYPE=OBSC 
        IF (FMASK(I).GT.0) TYPE=CLEAR 
        RMASK=0 
        XCEN=0 
        YCEN=0 
        XMIN=XMN(I) 
        XMAX=XMX(I) 
        YMIN=YMN(I) 
        YMAX=YMX(I) 
        WRITE(6,511) I,TYPE,XMIN,XMAX,YMIN,YMAX 
  185 END DO 
!                                                                       
  186 IF (IMASKE.EQ.0) GO TO 190 
      WRITE(6,502) 
      DO 188 I=1,NSURF 
         IF (IFLAG(I).EQ.0) GO TO 188 
         IF (FAKEC(I).NE.-1) GO TO 188 
         IF (FMASK(I).LT.0) TYPE=OBSC 
         IF (FMASK(I).GT.0) TYPE=CLEAR 
         RMASK=0 
         YMIN=0 
         XMIN=0 
         YMAX=YMX(I) 
         XMAX=XMX(I) 
         YCEN=YMN(I) 
         XCEN=XMN(I) 
         WRITE(6,511) I,TYPE,XMAX,YMAX,XCEN,YCEN 
  188 END DO 
!                                                                       
  190 WRITE (6,520) 
!                                                                       
!             TURN OFF LEPRT SWITCH                                     
  195 IF ( IPRINT.NE.IZERO.AND.(IPRINT/ITWO)*ITWO.NE.IPRINT)            &
     &     IPRINT = IPRINT-IONE                                         
!                                                                       
      RETURN 
  200 FORMAT (///T45,'OPTICAL PRESCRIPTION DATA') 
  210 FORMAT (//T37,'----------     OBJECT DATA      ----------'/) 
  214 FORMAT (T10,'WRITING ONLY FIRST 6 OBJECT POINTS'//) 
  220 FORMAT (T10,'X OBJECT ANGLES : ',1X,6(1PE13.6,1X)) 
  222 FORMAT (T29,6(1PE13.6,1X)) 
  230 FORMAT (T10,'Y OBJECT ANGLES : ',1X,6(1PE13.6,1X)) 
  240 FORMAT (/T10,'X OBJECT HEIGHTS: ',1X,6(1PE13.6,1X)) 
  250 FORMAT (T10,'Y OBJECT HEIGHTS: ',1X,6(1PE13.6,1X)) 
  260 FORMAT (T10,'OBJECT DISTANCE = ',1X,1PE16.9) 
  262 FORMAT (T10,'OBJECT INDEX    =   AIR'/) 
  264 FORMAT (T10,'OBJECT INDEX    =   REFLECT'/) 
  266 FORMAT (T10,'OBJECT INDEX    = ',1X,0PF9.6/) 
  270 FORMAT (//T37,'----------  IMAGE SURFACE DATA  ----------',       &
     &       //T10,'NUMBER OF SURFACES = ',                             &
     &       1X,I2,/T10,'SEPARATION         = ',1X,1PE16.9,             &
     &       /T10,'CURVATURE          = ',1X,1PE16.9,                   &
     &       5X,'( RADIUS = ',1PE16.9,1X,')',                           &
     &       /T10,'CONIC CONSTANT     = ',1X,0PF9.6)                    
  275 FORMAT (//T37,'----------  IMAGE SURFACE DATA  ----------',       &
     &       //T10,'NUMBER OF SURFACES = ',                             &
     &       1X,I2,/T10,'SEPARATION         = ',1X,1PE16.9,             &
     &       /T10,'CURVATURE          = ',1X,1PE16.9,                   &
     &       5X,'( INFINITE RADIUS )',                                  &
     &       /T10,'CONIC CONSTANT     = ',1X,0PF9.6)                    
  277 FORMAT(//T37,'----------     SYSTEM UNITS     ----------') 
  280 FORMAT (//T10,'SYSTEM UNITS ARE MILLIMETERS') 
  290 FORMAT (//T10,'SYSTEM UNITS ARE CENTIMETERS') 
  300 FORMAT (//T10,'SYSTEM UNITS ARE INCHES') 
  310 FORMAT (//T10,'SYSTEM UNITS ARE TANGENTS OF ANGLES') 
  320 FORMAT (//T10,'SYSTEM UNITS ARE ANGLES OF INCIDENCE') 
  330 FORMAT (//T10,'UNDEFINED SYSTEM UNITS; UFLAG = ',0PF5.1) 
  340 FORMAT (//T37,'---------- ENTRANCE PUPIL DATA  ----------',       &
     &       //T10,'DIAMETER                   = ',1PE16.8,             &
     &       /T10,'DISTANCE TO FIRST SURFACE  = ',1PE16.8)              
  350 FORMAT (//T37,'----------   EXIT PUPIL DATA    ----------',       &
     &       //T10,'DIAMETER                    = ',1PE16.8,            &
     &       /T10,'DISTANCE FROM LAST SURFACE  = ',1PE16.8)             
  353 FORMAT (//T9,' APERTURE STOP AT SURFACE ',I4) 
  355 FORMAT (//T9,' NO APERTURE STOP DESIGNATED') 
  357 FORMAT (//T9,' NO REFERENCE SURFACE DESIGNATED') 
  358 FORMAT (//T9,' REFERENCE SURFACE AT SURFACE ',I4) 
  360 FORMAT (//T37,'----------  BASIC SURFACE DATA  ----------',       &
     &       //T5,'SURF',T13,'CURVATURE',T32,                           &
     &       'RADIUS',T48,'CONIC',T63,'THICKNESS',T82,'INDICES'/)       
  370 FORMAT(T5,I4,T10,1PE15.8,T31,'INFINITE',T46,0PF9.6,T60,1PE15.8,   &
     &       T78,'AIR')                                                 
  372 FORMAT(T5,I4,T10,1PE15.8,T28,1PE15.8,T46,0PF9.6,T60,1PE15.8,      &
     &       T78,'AIR')                                                 
  374 FORMAT(T5,I4,T10,1PE15.8,T31,'INFINITE',T46,0PF9.6,T60,1PE15.8,   &
     &       T78,'REFLECT')                                             
  376 FORMAT(T5,I4,T10,1PE15.8,T28,1PE15.8,T46,0PF9.6,T60,1PE15.8,      &
     &       T78,'REFLECT')                                             
  377 FORMAT(T5,I4,T10,1PE15.8,T31,'INFINITE',T46,0PF9.6,T60,1PE15.8,   &
     &       T77,3(0PF9.6,1X) )                                         
  378 FORMAT(T5,I4,T10,1PE15.8,T28,1PE15.8,T46,0PF9.6,T60,1PE15.8,      &
     &       T77,3(0PF9.6,1X) )                                         
  380 FORMAT (//T37,'----------    ASPHERIC DATA     ----------',       &
     &       //T10,'SURF',T23,'4TH',T43,'6TH',                          &
     &       T63,'8TH',T82,'10TH'/)                                     
  390 FORMAT (T10,I4,T17,1PE15.8,T37,1PE15.8,T57,1PE15.8,T77,1PE15.8) 
  400 FORMAT (//T37,'----------TILT/DISPLACEMENT DATA----------',       &
     &       //T73,'(TILTS ARE IN DEGREES)',/T10,'SURF',T22,'X-DEC',    &
     &       T42,'Y-DEC',T62,'X-TILT',T80,'Y-TILT',T101,'Z-TILT'/)      
  410 FORMAT (T10,I4,T17,1PE15.8,T37,1PE15.8,T57,1PE15.8,T77,1PE15.8,   &
     &       T97,1PE15.8)                                               
  420 FORMAT (T23,'(DISPLACEMENTS ARE RESTORED)') 
  430 FORMAT (T70,'(TILTS ARE RESTORED)') 
  440 FORMAT (T23,'(DISPLACEMENTS ARE RESTORED)',T69,                   &
     &        T70,'(TILTS ARE RESTORED)')                               
  450 FORMAT (//T37,'----------     GRATING DATA     ----------',       &
     &       //T10,'SURF',T21,'WAVELENGTH',T41,                         &
     &       'ORDER',T56,'SPACING'/)                                    
  460 FORMAT (T10,I4,T19,1PE15.8,T42,I4,                                &
     &       T52,1PE15.8,T72,'(X RULINGS)')                             
  470 FORMAT (T10,I4,T19,1PE15.8,T42,I4,                                &
     &       T52,1PE15.8,T72,'(Y RULINGS)')                             
  480 FORMAT (//T37,'----------      TORIC DATA      ----------',       &
     &       //T10,'SURF',T20,'CURVATURE',T39,                          &
     &       'RADIUS'/)                                                 
  490 FORMAT (T10,I4,T17,1PE15.8,T35,1PE15.8) 
  495 FORMAT(T10,I4,T17,1PE15.8,T37,'INFINITE') 
  499 FORMAT (//T37,'----------      MASK DATA       ----------') 
  500 FORMAT (    //T32,'CIRCULAR',                                     &
     &       //T3,'SURF',T10,'TYPE',T20,'RADIUS',T32,'X CENTER',T45,    &
     &       'Y CENTER'/)                                               
  501 FORMAT (//T38,'RECTANGULAR',//T3,'SURF',T10,'TYPE',               &
     &       T21,'X MIN',T34,'X MAX',T47,'Y MIN',T60,'Y MAX'/)          
  502 FORMAT (//T38,'ELLIPTICAL',//T3,'SURF',T10,'TYPE',T21,'X MAX',    &
     &       T34,'Y MAX',T45,'X CENTER',                                &
     &       T58,'Y CENTER'/)                                           
  510 FORMAT (T3,I4,T10,A5,T17,1PE12.5,T30,1PE12.5,T43,1PE12.5) 
  511 FORMAT (T3,I4,T10,A5,T17,1PE12.5,T30,1PE12.5,T43,1PE12.5,         &
     &       T56,1PE12.5)                                               
  520 FORMAT (//T45,'END OF PRESCRIPTION'///) 
      END                                           
!                                                                       
!*******************************************************                
      SUBROUTINE ROTM(ALPHA,BETA,GAMMA) 
!*******************************************************                
!                                                                       
!                ROTM CONSTRUCTS THE ROTA ROTATION MATRIX               
!                                                                       
!                ALPHA = TILT ABOUT X-AXIS, IN DEGREES                  
!                BETA  = TILT ABOUT (NEW) Y-AXIS, IN DEGREES            
!                GAMMA = TILT ABOUT (NEW) Z-AXIS, IN DEGREES            
!                ARAD  = TILT ABOUT X-AXIS, IN RADIANS                  
!                BRAD  = TILT ABOUT (NEW) Y-AXIS, IN RADIANS            
!                GRAD  = TILT ABOUT (NEW) Z-AXIS, IN RADIANS            
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      COMMON /ROT/ ROTA(3,3) 
      CONST = 3.14159265358979323846D0/180. 
      ARAD = ALPHA*CONST 
      BRAD = BETA*CONST 
      GRAD = GAMMA*CONST 
      CA   = DCOS(ARAD) 
      CB   = DCOS(BRAD) 
      CG   = DCOS(GRAD) 
      SA   = DSIN(ARAD) 
      SB   = DSIN(BRAD) 
      SG   = DSIN(GRAD) 
!                                                                       
      ROTA(1,1) = CB*CG 
      ROTA(1,2) = (SA*SB*CG + CA*SG) 
      ROTA(1,3) = -(CA*SB*CG - SA*SG) 
      ROTA(2,1) = -CB*SG 
      ROTA(2,2) = -(SA*SB*SG - CA*CG) 
      ROTA(2,3) = (CA*SB*SG + SA*CG) 
      ROTA(3,1) = SB 
      ROTA(3,2) = -SA*CB 
      ROTA(3,3) = CA*CB 
      RETURN 
      END                                           
!                                                                       
!*******************************************************                
      SUBROUTINE SKEW (NFOC,LFOC,IFOC) 
!*******************************************************                
!                ROUTINE SKEW PERFORMS TRACING OF RAYS THROUGH          
!                A SYSTEM OF UP TO 40 SURFACES                          
      IMPLICIT REAL *8 (A-H,O-Z) 
      INTEGER *4 OPTNA,OPTNB,OPTNC,OPTND,OPTNE,OPTNF,OPTNG,ANULI,SECTRS 
!                                                                       
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /HEAD/   LINES,IPAGE,NSYS,LAMDA,NAME(160) 
      COMMON /PUPIL/  ENPUPR,ENPUPL,EXPUPR,EXPUPL 
      COMMON /SAGPAR/ YMAX,YMIN,DELY,REFCRV,CONST,WAVENM,ISRF 
      COMMON /CRED/   RPP(40),EPP(40),XW(812),HYSAVE,HXSAVE,            &
     &                TIMES,RD,OPTNE,ISW,IOPB                           
      COMMON /CMTF/   AJ, AN, AR(40), DE(40), DELTCC, DELTEN, DELTPL,   &
     &                EN, IOPC, OPTNF, PRINAN, RSC, SUM, TX(51)         
      COMMON /CSPOT/  XTMAX,XTMIN,YTMAX,YTMIN,AVGX,AVGY,RP(812),SPOTP,  &
     &                XK(812),YK(812),XKALL(7600),YKALL(7600),JSKIP,    &
     &                IOPA,NTHRU                                        
!                                                                       
!                                                                       
!                                                                       
      DIMENSION YW(812),TOR(40),XWN(812),YWN(812) 
      DIMENSION ZF(10),BLOB(10) 
      DIMENSION XSUM(10),YSUM(10),XSUMSQ(10),YSUMSQ(10) 
      DIMENSION QT(3),XX(3),XYZ(3),GIBRSH(8) 
!                                                                       
      EQUIVALENCE (QT(1),QX),(QT(2),QY),(QT(3),QZ) 
      EQUIVALENCE (XX(1),XT),(XX(2),YT),(XX(3),ZT) 
      EQUIVALENCE (XYZ(1),X),(XYZ(2),Y),(XYZ(3),Z) 
      EQUIVALENCE (TRASH(1),XSMIN),(TRASH(2),XSMAX),(TRASH(3),YSMIN) 
      EQUIVALENCE (TRASH(4),YSMAX),(TRASH(5),FCODE),(TRASH(6),FXYJ) 
      EQUIVALENCE (TRASH(7),FXY),(TRASH(8),FXNU),(TRASH(9),FBY) 
      EQUIVALENCE (TRASH(10),FBNU),(Y0(1),TOR(1)) 
!                                                                       
!                                                                       
      DATA IZERO/0/,IONE/1/ 
      DATA ZERO/0.D0/,ONE/1.D0/,TWO/2.D0/,THREE/3.D0/,FOUR/4.D0/ 
      DATA FIVE/5.D0/,SIX/6.D0/,TEN/10.D0/,HUNDRD/100.D0/,THOUS/1000.D0/ 
      DATA EPS1/1.D-9/,EPS2/1.D-14/,PI/3.141592653589793D0/ 
      DATA EPSLON/1.D-8/ 
!                                                                       
   10 FORMAT (//10H  RAY NO. ,I3,/,2X,4HSURF,7X,1HX,15X,1HY,15X,1HZ,14X &
     &      ,2HQX,14X,2HQY,14X,2HQZ,8X,19HOPTICAL PATH LENGTH/)         
   20 FORMAT (3X,I2,1P8E16.8) 
   30 FORMAT (//,6X,1HX,14X,1HY,14X,1HZ) 
   40 FORMAT (1P3E15.7) 
   50 FORMAT (3X,I2,10H RAY MISS ) 
   60 FORMAT (3X,I2,15H RAY REFLECTION) 
   70 FORMAT (3X,I2,13H RAY VIGNETTE) 
   80 FORMAT (13X,3HLOC,15X,2HX=,15X,2HY=,15X,2HZ=) 
   90 FORMAT (/4X,12HIMAGE PLANES,10X,6HQX/QZ=,1PE15.8,10X,6HQY/QZ=,    &
     &       1PE15.8)                                                   
  100 FORMAT (5X,6(2X,1PE15.8)) 
  105 FORMAT (/4X,10HCOLOR NO. ,I1,4X,23HMULTIPLE OBJECT HEIGHTS,       &
     &       /4X,11HRAYS THRU =,I7,1X,7H REFL =,I4,1X,7H MISS =,        &
     &       I4,1X,7H VIGN =,I4/)                                       
  110 FORMAT (/4X,10HCOLOR NO. ,I1,4X,10HY HEIGHT =,1PE15.8,4X,         &
     &       11HX HEIGHT = ,1PE15.8,/4X,11HRAYS THRU =,I4,1X,7H REFL =, &
     &        I4,1X,7H MISS =,I4,1X,7H VIGN =,I4/)                      
  120 FORMAT (13X,3HLOC,14X,4HXBAR,13X,4HYBAR,13X,4HSDVX,13X,4HSDVY,9X, &
     &       11HSPOT RADIUS)                                            
  130 FORMAT (3X,2HEP,1P6E16.8) 
  140 FORMAT (10X,'NO ENTRANCE PUPIL SPECIFIED FOR FINRAY, WILL NOT ',  &
     &       'PERFORM RAY TRACE'/)                                      
  201 FORMAT(10X,'APERTURE STOP SIZE NOT SPECIFIED, WILL NOT LOCATE ',  &
     &       'PARAXIAL PUPILS'/)                                        
  202 FORMAT(10X,'FOUND DECENTERED SURFACE WHILE LOCATING PUPILS'/) 
  203 FORMAT(10X,'NUMBER OF IMAGE SURFACES = 0 OR BAD NUMBERING'/) 
  204 FORMAT(10X,'FIRST PLANE NUMBER IS NOT AN INTEGER'/) 
  205 FORMAT(10X,'OVERFLOW OF LATTICE ARRAY'/) 
  206 FORMAT(10X,'NO APERTURE STOP SPECIFIED',                          &
     &      ', WILL NOT LOCATE PARAXIAL PUPILS'/)                       
  210 FORMAT(10X,'QZ = 0 AT IMAGE, NO ANALYSIS PERFORMED'/) 
  211 FORMAT(6D16.8) 
  212 FORMAT(8D16.8) 
  213 FORMAT(10X,'UNABLE TO TRACE CHIEF RAY FOR Y HEIGHT = ',1PE14.8,   &
     &      ' X HEIGHT = ',1PE14.8/)                                    
  218 FORMAT(10X,'NUMBER OF SPOT PLOT RAYS TOO LARGE'/) 
!                                                                       
      QZOLD = 0.0D0 
      JKSAVE=0 
      IF ((IPRINT/2)*2.NE.IPRINT) JKSAVE=1 
      IF (JKSAVE.EQ.1) FFXY=FXY 
      IF (JKSAVE.EQ.1) FFXNU=FXNU 
      IF (JKSAVE.EQ.1) FFBY=FBY 
      IF (JKSAVE.EQ.1) FFBNU=FBNU 
      ISTOP=APSTOP/30 
      IF (LATYPE.NE.4.OR.RHO.NE.0.OR.ISTOP.NE.0) GO TO 219 
      WRITE(6,140) 
      GO TO 1910 
  219 LSTAT=FAKEA 
      IF (LSTAT.GE.1) THCK=T(LSTAT) 
      DLMP=DELIMP 
      FPLN=FPLANE 
      NPLN=NPLANE 
      REWIND 20 
      DIV=TEN 
      LIMIT=800 
      MAXIT=6 
      ENPUPR = RHO 
      ENPUPL = D 
      NNSURF=NSURF 
      NITRAT=0 
!                                                                       
!                A-INITIALIZATION-                                      
!                ISTOP = SURFACE NUMBER OF APERTURE STOP                
      ISTOP = APSTOP/30 
!                   IREF:  0 = FINISHED REFERENCE SURFACE ITERATIONS    
!                         -1 = NO REFERENCE SURFACE ITERATIONS DESIRED  
!                          1+ = NOW ITERATING TO SURFACE IREF           
      IF (IREF.LE.IZERO) IREF=-1 
      IIREF=IREF 
      LAMDA = 2 
      IF ( NCOL.EQ.1 ) LAMDA = 1 
      IF (LATYPE.EQ.5.OR.LATYPE.EQ.6) IREF=-1 
      IF (LATYPE.EQ.5.OR.LATYPE.EQ.6) GO TO 300 
!                IF NONE SPECIFIED, LEAVE ROUTINE                       
      IF ( ISTOP.LT.1.AND.(IPRINT/2)*2.EQ.IPRINT) WRITE(6,206) 
      IF ( ISTOP.LT.1 ) GO TO 300 
!                B ENTRANCE PUPIL                                       
!                B LOCATOR                                              
      IF ( FMASK(ISTOP).GT.ZERO ) GO TO 220 
      IF ((IPRINT/2)*2.EQ.IPRINT) WRITE(6,201) 
      GO TO 300 
  220 AB=ZERO 
      ABB=FMASK(ISTOP) 
      IF (ISTOP.EQ.1) GO TO 230 
      SGN1 = FREF0 
      L=ISTOP-1 
      DO 225 I=1,L 
        IF (FREF(I).EQ.-ONE) SGN1 = -SGN1 
  225 END DO 
      SGNAP = SGN1*FREF(ISTOP) 
      TM1=.1D0 
      TM1B=TM1 
      GO TO 240 
  230 TM1=.1D0 
      TM1B=TM1 
      SGNAP = FREF0 * FREF(1) 
      SGN1 = FREF0 
      GO TO 260 
!                BACK TRACE (PARAXIALLY) A RAY THROUGH THE SYSTEM       
  240 DO 250 I=1,L 
        M=ISTOP-I 
        IF ((IPRINT/2)*2.NE.IPRINT) GO TO 245 
!                IF TILTED OR DECENTERED SURFACE BETWEEN APERTURE STOP  
!                AND ENTRANCE PUPIL, TELL USER                          
        IF (TILTX(M).NE.ZERO .OR. TILTY(M).NE.ZERO .OR. TILTZ(M).NE.ZERO&
     &     .OR. XDISP(M).NE.ZERO .OR. YDISP(M).NE.ZERO ) WRITE(6,202)   
  245   SGN2 = SGN1 
        IF (FREF(M).EQ.-ONE) SGN2 = -SGN1 
        IF (M.EQ.1) RN1 = OBJN(LAMDA)*SGN2 
        IF (M.GT.1) RN1 = FN(M-1,LAMDA)*SGN2 
        RN2  = FN(M,LAMDA)*SGN1 
        GNU  = RN2/RN1 
        ABB  = ABB + T(M)*TM1B 
        AB   = AB + T(M)*TM1 
        TM1B = TM1B + ABB*C(M)*(ONE-GNU) 
        TM1  = TM1 + AB*C(M)*(ONE-GNU) 
        SGN1 = SGN2 
  250 END DO 
!                D IS INTERSECTION OF RAY AND AXIS WITH RESPECT TO      
!                SURFACE 1 COORD SYSTEM (D = 0 IF ISTOP =1)             
  260 TM7=-AB/TM1 
      D = TM7 
      RHO=DABS(ABB+TM7*TM1B) 
      ENPUPR = RHO 
      ENPUPL = D 
!             LOCATE AND SIZE PARAXIAL EXIT PUPIL; PRINT INFORMATION    
      AB=ZERO 
      ABB=FMASK(ISTOP) 
      L=ISTOP+1 
      SGN1 = SGNAP 
      TM1=.1D0 
      TM1B=TM1 
      DO 290 I=L,NSURF 
        IF ((IPRINT/2)*2.NE.IPRINT) GO TO 285 
        IF (TILTX(I).NE.ZERO .OR. TILTY(I).NE.ZERO .OR. TILTZ(I).NE.ZERO&
     &     .OR. XDISP(I).NE.ZERO .OR. YDISP(I).NE.ZERO ) WRITE(6,202)   
  285   SGN2 = SGN1 
        IF (FREF(I).EQ.-ONE) SGN2 = -SGN1 
        RN1  = FN(I,LAMDA)*SGN1 
        RN2  = FN(I-1,LAMDA)*SGN1 
        GNU  = RN2/RN1 
        ABB  = ABB + T(I-1)*TM1B 
        AB   = AB + T(I-1)*TM1 
        TM1B = TM1B + ABB*C(I)*(GNU-ONE) 
        TM1  = TM1 + AB*C(I)*(GNU-ONE) 
        SGN1 = SGN2 
  290 END DO 
      TM7=-AB*FN(NSURF,LAMDA)/TM1 
      TM8=DABS(ABB+TM7*TM1B/FN(NSURF,LAMDA)) 
      EXPUPR = TM8 
      EXPUPL = TM7 
!              END OF EXIT PUPIL LOCATION                               
!                                                                       
  300 HYMAX=HYINIT+(RNOBJ-ONE)*HYDEL 
      HYMAX=-DTAN(HYMAX*PI/180.)*(S-D) 
      HXMAX=HXINIT+(RNOBJ-ONE)*HXDEL 
      HXMAX=-DTAN(HXMAX*PI/180.)*(S-D) 
      FXNU=RHO/(S-D) 
      FXY=FXNU*(S-D) 
      FBNU = -HYMAX/(S-D) 
      FBY=ZERO 
!                END EP LOCATOR                                         
!                INITIALIZE SYSTEM                                      
!                XSUM, YSUM ARE ACCUMULATORS FOR AVERAGES FOR           
!                EACH IMAGE PLANE COORD                                 
!                XSUMSQ, YSUMSQ ARE USED FOR RMS CALCULATIONS           
  305 ITEMP=FPLANE 
      IF ((IPRINT/2)*2.NE.IPRINT) GO TO 308 
!                CHECK IMAGE PLANE NUMBERING                            
      IF ( NPLANE+ITEMP.LT.0 .OR. ITEMP.GT.0 .OR. NPLANE.LE.0 )         &
     &     WRITE(6,203)                                                 
  308 IF (NPLANE.LE.0) GO TO 320 
      DO 310 I=1,NPLANE 
        XSUM(I)=ZERO 
        YSUM(I)=ZERO 
        XSUMSQ(I)=ZERO 
        YSUMSQ(I)=ZERO 
  310 END DO 
!                NIMG = NUMBER OF PRIME IMAGE PLANE                     
  320 NIMG=1-(FPLANE) 
      TEMP=1-NIMG 
      PRINAN=ZERO 
      IF ((IPRINT/2)*2.NE.IPRINT) GO TO 325 
!                ALARM IF FIRST PLANE NO. NOT INTEGER                   
      IF ( TEMP.NE.FPLANE ) WRITE(6,204) 
  325 NMISS=0 
      NREFL=0 
      NVIGN=0 
      HYLAST=-DTAN(HYINIT*PI/180.)*(S-D) 
      HXLAST=-DTAN(HXINIT*PI/180.)*(S-D) 
      LAMDA=ICOL(1) 
!             VARIABLE           OPTION                                 
!             OPTION A   PRINTS PRIME IMAGE COORD ONLY                  
!             OPTION B   PRINTS COORD AND COSINES IN EP AND IMAGE       
!             OPTION C   PRINTS COORD FOR SURFACES                      
!             OPTION D   CAUSES PRESCRIPTION MATRIX PRINT               
!             OPTION E   PRINTS RED TABLES                              
!             OPTION F   PRINTS MTF TABLES                              
!             OPTION G   CAUSES ANALYSIS FOR EACH HEIGHT, COLOR         
      OPTNA=(IPRINT/8)-(IPRINT/16)*2 
      OPTNB=(IPRINT/2)-(IPRINT/4)*2 
      OPTNC=(IPRINT/4)-(IPRINT/8)*2 
      OPTND=IPRINT-(IPRINT/2)*2 
      OPTNE=(IPRINT/16)-(IPRINT/32)*2 
      OPTNF=(IPRINT/32)-(IPRINT/64)*2 
      OPTNG=1 
!                                                                       
!                IF NO VALID CODE, BYPASS ROUTINE                       
      IF ( IPLTPR.GT.15 .OR. IPLTPR.LT.0 ) IPLTPR = 0 
!                EXTRACT SWITCHES                                       
!                IOPA - SPOT PLOT (UNITS) SWITCH                        
!                IOPB - RADIAL ENERGY DIST PLOT SW                      
!                IOPC - MODULATION TRANSFER FUNCTION PLOT SW            
      IOPA = IPLTPR - (IPLTPR/2)*2 
      IOPB = (IPLTPR/4) - (IPLTPR/8)*2 
      IOPC = (IPLTPR/8) - (IPLTPR/16)*2 
!                                                                       
!             PRINT PRESCRIPTION MATRIX                                 
      IF (OPTND.EQ.0) GO TO 327 
      CALL PREPRT 
      FXY=FFXY 
      FXNU=FFXNU 
      FBY=FFBY 
      FBNU=FFBNU 
      GO TO 1910 
!                SET OPTIONB IF OPTIONC                                 
  327 IF (OPTNC.NE.0) OPTNB=1 
!                RESET OPTIONB AND OPTIONC IF OPTIONA                   
      OPTNB=-OPTNB*(OPTNA-1) 
      OPTNC=-OPTNC*(OPTNA-1) 
!                                                                       
      IF ( ISRF.NE.IZERO ) CALL SRFSAG 
      IF (LATYPE.NE.5.AND.LATYPE.NE.6) GO TO 328 
        D=0 
        GO TO 525 
!                                                                       
!                B LATTICE                                              
!                B GENERATION                                           
!                LATYPE =1 IF SINGLE RAY                                
!                       =2 IF POLAR LATTICE                             
!                       =3 IF RECTANGULAR LATTICE                       
!                IF IREF>1 GENERATE RAY COORDS                          
!                   AT REFERENCE SURFACE AS WELL AS                     
!                   AT ENTRANCEPUPIL                                    
!                FAKEC(IREF)  =-1  IF ELLIPTICAL MASK                   
!                             = 0  IF CIRCULAR MASK                     
!                             = 1  IF RECTANGULAR MASK                  
  328 GO TO (330,340,400,470,525,525,392), LATYPE 
!                C ONE RAY                                              
!                C LATTICE                                              
!                SET UP SINGLE RAY, GO FIND COSINES                     
  330 XW(1)=RHO*CLTRA(1) 
      YW(1)=RHO*CLTRA(2) 
      IF (IREF.LE.IZERO) GO TO 333 
!           CIRCULAR MASK                                               
      XWN(1)=FMASK(IREF)*CLTRA(1) 
      YWN(1)=FMASK(IREF)*CLTRA(2) 
      IF (FAKEC(IREF).EQ.0) GO TO 333 
!           RECTANGULAR OR ELLIPTICAL MASK                              
      XWN(1)=XMX(IREF)*CLTRA(1) 
      YWN(1)=YMX(IREF)*CLTRA(2) 
  333 NUMPTS=1 
      GO TO 480 
!                CC END ONE RAY LATTICE                                 
!                C POLAR                                                
!                C LATTICE                                              
!                SET UP POLAR LATTICE                                   
!                ANULI = NUMBER OF ANNULI                               
!                SECTRS = NUMBER OF SECTORS                             
!                XW = X COORD OF POINT IN 1/2 CIRCLE OF POINTS          
!                YW = Y COORD OF POINT IN 1/2 CIRCLE OF POINTS          
  340 ANULI=CLTRA(1) 
      SECTRS=CLTRA(2) 
      ITATS=ANULI*SECTRS 
      IF ( ITATS.EQ.0 ) GO TO 370 
!                TA*TA/2 POINTS ARE CREATED                             
      LTATS=ITATS/2 
      IF (LTATS.LE.LIMIT) GO TO 360 
!                ALARM IF POINTS O-FLOW                                 
  350 WRITE(6,205) 
      GO TO 1890 
  360 ITS2=SECTRS/2 
      IF (ITS2*2.EQ.SECTRS) GO TO 380 
!                ALARM IF ASSYMMETRIC                                   
  370 WRITE(6,205) 
      GO TO 1890 
!                CALCULATE POINTS COORD                                 
  380 PIOVTS=PI/SECTRS 
      ATA=ANULI 
      TEMP=FOUR*DSIN(PIOVTS/TWO)/(THREE*PIOVTS) 
      K=0 
      DO 390 I=1,ANULI 
        AI=I 
        RI=DSQRT(AI/ATA) 
        RIM1=DSQRT((AI-ONE)/ATA) 
        RHOBAR=TEMP*(RI**3-RIM1**3)/(RI*RI-RIM1*RIM1) 
        IF (IREF.LE.0) GO TO 385 
        IF (FAKEC(IREF))  381,382,383 
!             SCALE RHOBAR TO REFERENCE SURFACE                         
!             STORE RESULT IN REFBAR                                    
!                                                                       
!             ELLIPTICAL MASK                                           
  381   IF (XMX(IREF).GE.YMX(IREF)) ELIPRD=XMX(IREF) 
        IF (XMX(IREF).LT.YMX(IREF)) ELIPRD=YMX(IREF) 
        REFBAR=ELIPRD*RHOBAR 
        GO TO 385 
!             CIRCULAR MASK                                             
  382   REFBAR=FMASK(IREF)*RHOBAR 
        GO TO 385 
!             RECTANGULAR MASK                                          
  383   REFBAR=DSQRT(XMX(IREF)*XMX(IREF)+YMX(IREF)*YMX(IREF))*RHOBAR 
!             SCALE RHOBAR TO ENTRANCEPUPIL                             
!             STORE RESULT IN RHOBAR                                    
  385   RHOBAR=RHO*RHOBAR 
        DO 390 J=1,ITS2 
          ARG=PIOVTS*(TWO*J-ONE)-PI/TWO 
          K=K+1 
          XW(K)=RHOBAR*DCOS(ARG) 
          YW(K)=RHOBAR*DSIN(ARG) 
          IF (IREF.LE.IZERO) GO TO 390 
          XWN(K)=REFBAR*DCOS(ARG) 
          YWN(K)=REFBAR*DSIN(ARG) 
  390 CONTINUE 
!                GO CALCULATE COSINES                                   
      NUMPTS=K 
      GO TO 480 
!                CC END POLAR LATTICE                                   
!                RIM LATTICE ROUTINE                                    
!                NUMPTS = NUMBER POINTS TOTAL                           
!                DTHETA = ANGULAR INCREMENT                             
!                RIMANG = CURRENT ANGLE                                 
!                J = NUMBER POINTS THUSFAR MADE                         
  392 NUMPTS=CLTRA(1) 
      DTHETA=CLTRA(2) 
      IF (NUMPTS.LE.0.OR.NUMPTS.GT.LIMIT) GO TO 350 
      RIMANG=-PI/2.D0 
      J=0 
  393 J=J+1 
      RIMANG=RIMANG+DTHETA 
      XW(J)=DCOS(J*RIMANG) 
      YW(J)=DSIN(J*RIMANG) 
      IF (J.GT.LIMIT) GO TO 350 
      IF (RIMANG+DTHETA.LE.PI/2.D0) GO TO 393 
      NUMPTS=J 
      DO 399 I=1,NUMPTS 
        XWIHLD=XW(I) 
        YWIHLD=YW(I) 
!                SCALE POINTS AT ENTRANCEPUPIL                          
        XW(I)=RHO*XW(I) 
        YW(I)=RHO*YW(I) 
        IF (IREF.LE.IZERO) GO TO 399 
!                SCALE POINTS AT REFERENCE SURFACE                      
        IF (FAKEC(IREF)) 396,397,398 
!                ELLIPTICAL MASK                                        
  396  IF (XMX(IREF).GE.YMX(IREF)) RDWN=XMX(IREF) 
       IF (XMX(IREF).LT.YMX(IREF)) RDWN=YMX(IREF) 
       XWN(IREF)=RDWN*XWIHLD 
       YWN(IREF)=RDWN*YWIHLD 
       GO TO 399 
!                CIRCULAR MASK                                          
  397  XWN(I)=FMASK(IREF)*XWIHLD 
       YWN(I)=FMASK(IREF)*YWIHLD 
       GO TO 399 
!                RECTANGULAR OR ELLIPTICAL MASK                         
  398  XWN(I)=XMX(IREF)*XWIHLD 
       YWN(I)=YMX(IREF)*YWIHLD 
  399 END DO 
      GO TO 480 
!                CC END RIM LATTICE                                     
!                RECTANGULAR LATTICE ROUTINE                            
!                NUMPTS = NUMBER POINTS TOTAL                           
!                DELY = Y INCREMENT                                     
!                NUMCOL = NUM POINTS THIS COL                           
!                IND = FW INDEX OF NEXT COLUMN SET                      
!                NSUM = NUMBER POINTS THUSFAR MADE                      
  400 NUMPTS=CLTRA(1) 
      DELY=CLTRA(2) 
      NUMCOL=CLTRA(3) 
      XW(1)=CLTRA(4) 
      YW(1)=CLTRA(5) 
      IND=6 
      K=2 
      NSUM=NUMCOL 
      IF (NUMPTS.LE.0) GO TO 350 
      IF (NUMPTS.GT.LIMIT) GO TO 350 
  410 IF (NSUM.GT.LIMIT) GO TO 350 
!                THIS DO LOOP CREATES ALL POINTS ABOVE FIRST FOR EACH   
!                COLUMN                                                 
      IF (NUMCOL.LT.2) GO TO 430 
      DO 420 I=2,NUMCOL 
        XW(K)=XW(K-1) 
        YW(K)=YW(K-1)+DELY 
        K=K+1 
  420 END DO 
!                IS LATTICE GENERATED                                   
  430 IF (NSUM-NUMPTS) 440,450,370 
!                NO, START NEXT COLUMN                                  
  440 NUMCOL=CLTRA(IND) 
      XW(K)=CLTRA(IND+1) 
      YW(K)=CLTRA(IND+2) 
      IND=IND+3 
      IF (IND.GT.298) GO TO 350 
      NSUM=NSUM+NUMCOL 
      K=K+1 
      GO TO 410 
!                THROUGH BUILDING RECT LATTICE, SCALE POINTS            
  450 DO 460 I=1,NUMPTS 
        XWIHLD=XW(I) 
        YWIHLD=YW(I) 
!                SCALE POINTS AT ENTRANCEPUPIL                          
        XW(I)=RHO*XW(I) 
        YW(I)=RHO*YW(I) 
        IF (IREF.LE.IZERO) GO TO 460 
!                SCALE POINTS AT REFERENCE SURFACE                      
        IF (FAKEC(IREF)) 458,457,458 
!                CIRCULAR MASK                                          
  457  XWN(I)=FMASK(IREF)*XWIHLD 
       YWN(I)=FMASK(IREF)*YWIHLD 
       GO TO 460 
!                RECTANGULAR OR ELLIPTICAL MASK                         
  458  XWN(I)=XMX(IREF)*XWIHLD 
       YWN(I)=YMX(IREF)*YWIHLD 
  460 END DO 
      GO TO 480 
!                CC END COLUMN LATTICE                                  
!                C RAY                                                  
!                C GENERATION                                           
  470 REWIND 40 
      REWIND 50 
      CALL FINRAY (NFOC) 
      GO TO 520 
!                REWIND LATTICE TAPES                                   
  480 REWIND 40 
      REWIND 50 
!                CALCULATE DIRECTIONAL COSINES FOR EACH POINT IN THE    
!                LATTICE AND WRITE THEM ON TAPE                         
      Z=ZERO 
      DZ=S-D 
      NOBJ=RNOBJ 
      IF (NUMPTS.LE.0) GO TO 350 
      IF (NOBJ.LE.0) GO TO 1890 
      IF (DABS(DZ).LT.EPS1) GO TO 1890 
!                                                                       
      I=0 
  490 I=I+1 
      IF (I.GT.NOBJ) GO TO 520 
      IF (NFOC.GT.0) I=NFOC 
!                HXANG,HYANG ARE OBJECT POINT IN DEGREES                
      HXANG=HXINIT+(I-1)*HXDEL 
      HYANG=HYINIT+(I-1)*HYDEL 
!                HX,HY ARE OBJECT POINT IN LINEAR DIMENSIONS            
      HX=-DTAN(HXANG*PI/180.)*DZ 
      HY=-DTAN(HYANG*PI/180.)*DZ 
      SGN=ONE 
  500 CONTINUE 
      DO 510 J=1,NUMPTS 
        XP=SGN*XW(J) 
        DX=XP-HX 
        DY=YW(J)-HY 
        DENOM=DSQRT(DX*DX+DY*DY+DZ*DZ) 
        QX=DX/DENOM 
        QY=DY/DENOM 
        QZ=DZ/DENOM 
!            DESX,DESY = DESIRED RAY COORDS AT REFERENCE SURFACE        
        DESX=IZERO 
        DESY=IZERO 
        IF (IREF.GT.IZERO) DESX=SGN*XWN(J) 
        IF (IREF.GT.IZERO) DESY=YWN(J) 
        WRITE (40) XP,YW(J),Z,QX,QY,QZ,HY,HX,DESX,DESY 
  510 END DO 
      SGN=-SGN 
      IF (IMODE.EQ.1 .AND. SGN.EQ.-ONE) GO TO 500 
      IF (NFOC.GT.0) GO TO 520 
      GO TO 490 
!                                                                       
  520 FSTP=-TWO 
      WRITE (40) (FSTP,I=1,10) 
      END FILE 40 
      IEOF=2 
!                CC END RAY GENERATION                                  
!                BB END LATTICE PROCESSING                              
!                AA END INITIALIZATION                                  
!                A COLOR PROCESSING                                     
!                EXECUTE COMPLETE TRACE FOR EACH COLOR                  
!                LAMDA = COLOR NUMBER THIS TRIP                         
!                NTIMES = NUMBER RAYS THROUGH SYSTEM                    
!                NRAY = NUMBER OF THIS RAY (SEQUENCE)                   
  525 IF (NCOL.EQ.0) GO TO 1890 
      INDX=1 
      INDY=1 
      NDCOL=NCOL 
      IF (IMODE.NE.0) GO TO 530 
      INDX=0 
      INDY=2 
  530 IF (LFOC.GT.0) GO TO 540 
      GO TO 550 
  540 LAMDA=ICOL(LFOC) 
      NDCOL=1 
  550 DO 1880 LAMDX=1,NDCOL 
!           X2DEC,Y2DEC = X,Y DISPLACEMENTS OF CHIEF RAY                
!                 AT ENTPUP WHEN ITERATED TO GO THROUGH                 
!                 ORIGIN OF REFERENCE SURFACE                           
!           KFLAG  = 0  SINGLE SPOT PLOT OR NO SPOT PLOT                
!                  = 1  MULTIPLE SPOT PLOT                              
!           IJK    =    NUMBER OF RAYS IN MULTIPLE SPOT PLOT            
        X2DEC=0 
        Y2DEC=0 
        IJK=0 
        KFLAG=0 
        IF (LFOC.GT.0) GO TO 560 
        LAMDA=ICOL(LAMDX) 
        IF (IREF.LE.1) CALL HEADIN(LAMDA) 
  560   NTHRU=0 
        NRAY=0 
        REWIND 40 
        REWIND 50 
!             HEAD FOR PRIME IMAGE ONLY PRINT                           
        IF (OPTNA.NE.0.AND.IEOF.NE.1) WRITE(6,30) 
!                B RAY PROCESSING                                       
!                INPUT ONE RAY                                          
  570   IF (LATYPE.NE.5) GO TO 571 
        READ(20,212) ORIGX,ORIGY,PZ,(QT(I),I=1,3),THISHY,THISHX 
        GO TO 573 
  571   IF (LATYPE.NE.6) GO TO 572 
        READ(20,212) ORIGX,ORIGY,PZ,(QT(I),I=1,3),THISHY,THISHX 
        THISHY=0 
        THISHX=0 
        GO TO 573 
  572   IF (IREF.GE.IONE.OR.IREF.EQ.-1)                                 &
     &  READ (40) ORIGX,ORIGY,PZ,(QT(I),I=1,3),THISHY,THISHX,DESX,DESY  
        IF (NRAY.GE.1.AND.IREF.GT.0) ORIGX=ORIGX+X2DEC 
        IF (NRAY.GE.1.AND.IREF.GT.0) ORIGY=ORIGY+Y2DEC 
        ZHOLD=PZ 
        QXHOLD=QX 
        QYHOLD=QY 
        QZHOLD=QZ 
        IF (IREF.EQ.IZERO)                                              &
     &  READ (50) ORIGX,ORIGY,PZ,(QT(I),I=1,3),THISHY,THISHX,DESX,DESY  
  573   PX=ORIGX 
        PY=ORIGY 
        NITRAT=0 
        ITREAD=IZERO 
        IWHICH=IZERO 
  575   IF (IREF.GT.IZERO) ITREAD=IONE 
        IF (IREF.GT.IZERO) IWHICH=IONE 
        IEOF=2 
!                CHECK FOR LAST RAY FLAG                                
        IF (QT(1).EQ.-TWO) IEOF=1 
        IF (IEOF.EQ.1.AND.KFLAG.EQ.1) GO TO 1890 
        IF (IEOF.EQ.1.AND.ITREAD.EQ.1) GO TO 1880 
!                END OF COLOR                                           
        IF (OPTNG.NE.0.AND.IEOF.EQ.1) GO TO 1350 
!                DUMMY BRANCH IF EOF ONLY                               
        IF (IEOF.EQ.1) GO TO 1880 
!                END OF RAYS FOR THIS HEIGHT                            
        IF ((THISHY.NE.HYLAST.OR.THISHX.NE.HXLAST).AND.KFLAG.EQ.1)      &
     &     GO TO 1360                                                   
        IF ((THISHY.NE.HYLAST.OR.THISHX.NE.HXLAST).AND.OPTNG.NE.0       &
     &     .AND.IREF.GT.IZERO) GO TO 1360                               
        IF ( (THISHY.NE.HYLAST.OR.THISHX.NE.HXLAST) .AND. OPTNG.NE.0 )  &
     &     GO TO 1350                                                   
        IF (KFLAG.EQ.1) GO TO 570 
!                                                                       
!                SAVE ALL COORDS, TRACE CHIEF RAY                       
  580   IF (NRAY.GT.0.OR.IREF.LE.0.OR.KFLAG.EQ.1)                       &
     &     GO TO 582                                                    
!             GIBRSH(1..8) STORE ALL DATA FROM RECORD 1                 
!                  OF UNIT 40 UNTIL CHIEF RAY ITERATION IS OVER         
        GIBRSH(1)=ORIGX 
        GIBRSH(2)=ORIGY 
        GIBRSH(3)=PZ 
        GIBRSH(4)=QX 
        GIBRSH(5)=QY 
        GIBRSH(6)=QZ 
        GIBRSH(7)=DESX 
        GIBRSH(8)=DESY 
        ORIGX=0 
        ORIGY=0 
        DESX=0 
        DESY=0 
        PZ=0 
        ZHOLD=PZ 
        QR=DSQRT(THISHX*THISHX+THISHY*THISHY+DZ*DZ) 
        QX=-THISHX/QR 
        QY=-THISHY/QR 
        QZ=DZ/QR 
        QXHOLD=QX 
        QYHOLD=QY 
        QZHOLD=QZ 
        PX=ORIGX 
        PY=ORIGY 
        NITRAT=0 
        ITREAD=1 
        IWHICH=1 
        IEOF=2 
        GO TO 585 
!                                                                       
!                REASSIGN COORDS AFTER REFERENCE SURF ITERATIONS        
  581   ORIGX=GIBRSH(1)+X2DEC 
        ORIGY=GIBRSH(2)+Y2DEC 
        PZ=GIBRSH(3) 
        QX=GIBRSH(4) 
        QY=GIBRSH(5) 
        QZ=GIBRSH(6) 
        DESX=GIBRSH(7) 
        DESY=GIBRSH(8) 
        PX=ORIGX 
        PY=ORIGY 
        ZHOLD=PZ 
        QXHOLD=QX 
        QYHOLD=QY 
        QZHOLD=QZ 
!                UPDATE RAY NUMBER                                      
  582   NRAY=NRAY+1 
!                PN = PREVIOUS SURFACE INDEX OF REF                     
!                PT = DIST FROM LAST SURF                               
  585   PN=OBJN(LAMDA) 
        PT=D 
        IF (ITREAD.GT.IZERO) GO TO 605 
        IF (OPTNB) 590,600,590 
!             HEAD SURFACE PRINTOUT AND PRINT EP INFO                   
  590   IF (IREF.GT.0) GO TO 600 
        WRITE(6,10) NRAY 
        WRITE(6,130) PX,PY,PZ,QT 
!                C TRACING SURFACE BY SURFACE                           
!                START TRACE LOOP FOR NSURF SURFACES                    
!             TEST IF SURF(I) IS IMAGE SURFACE                          
!             IF YES, PRINT/PLOT X AND Y                                
  600   IF (FAKEA.GT.ZERO) GO TO 610 
        GO TO 620 
  605   NSURF=IREF 
        GO TO 620 
  610   NSURF=FAKEA 
        TEMP=ZERO 
        T(NSURF)=ZERO 
        DELIMP=ZERO 
        FPLANE=ZERO 
        NPLANE=1 
        NIMG=1 
  620   DO 1200 I=1,NSURF 
          TM1=FN(I,LAMDA) 
          TM2=PN 
          TM3=C(I) 
          ISW=2 
          NTOR=TOR(I)+ONE 
!              NOW NTOR = 1 FOR CONIC OR POLYNOMIAL                     
!                  NTOR = 2 FOR TORIC                                   
!                  NTOR = 4 FOR BILATERAL SYMMETRIC SURFACE             
          GO TO (640,630,640,640), NTOR 
  630     ISW=3 
          GO TO 670 
  640     DO 650 IA=1,4 
            IF (COEF(I,IA)) 660,650,660 
  650     CONTINUE 
          GO TO 670 
!                                                                       
!                ISW = 1 IF POLYNOMIAL SURFACE                          
!                ISW = 2 IF CONIC SURFACE                               
!                ISW = 3 IF TORIC                                       
  660     ISW=1 
!                D TRANSFER EQUATIONS                                   
!                E TRANSFER COORDINATES TO TILTED/DECENTERED            
!                E TANGENT PLANE TO SURFACE                             
  670     DN = (PT-PZ)*QZ - PY*QY - PX*QX 
          DN1 = DN 
          XT = PX + DN*QX - XDISP(I) 
          YT = PY + DN*QY - YDISP(I) 
          ZT = PZ + DN*QZ - PT 
          DN = ZT 
!                                                                       
!                (XT,YT,DZ) IS INTERSECTION OF RAY WITH CLOSEST POINT TO
!                VERTEX OF SURFACE I, UNTILTED, DECENTERED              
!                                                                       
!                IF ALL TILTS = 0, SKIP TILTING ROUTINE                 
          IF ( DABS(TILTX(I)).EQ.ZERO .AND. DABS(TILTY(I)).EQ.ZERO .AND.&
     &         DABS(TILTZ(I)).EQ.ZERO ) GO TO 680                       
!                                                                       
!                CALL ROTM TO CONSTRUCT ROTATION MATRIX                 
          ALPHA = TILTX(I) 
          BETA  = TILTY(I) 
          GAMMA = TILTZ(I) 
!                CALL ROTM TO CONSTRUCT ROTATION MATRIX                 
          CALL ROTM (ALPHA,BETA,GAMMA) 
!                CALL MAVEC TO ROTATE COORD, COSINES (SEE EQUIVALENCES) 
          CALL MAVEC (XX,XX) 
          CALL MAVEC (QT,QT) 
!                AT THIS POINT, COSINES ARE ROTATED AND (XT,YT,ZT) IS   
!                THE POINT OF INTERSEC OF RAY WITH DECENTERED, TILTED   
!                CLOSEST-POINT PLANE (NEW XT, NEW YT, NEW ZT)           
!                                                                       
  680     IF (UFLAG.EQ.FIVE) QZOLD=QZ 
!                END OF CLOSEST-POINT TRANSFER                          
!                TRANSFER RAY TO SURFACE (IE FIND COORDINATES OF RAY    
!                ON SURFACE)                                            
          IF (NTOR.EQ.4) ISW = 1 
          GO TO (690,690,730), ISW 
!                CONIC, SPHERIC, PLANE TRANSFER EQUATIONS               
  690     IF ( C(I).EQ.ZERO .AND. QZ.EQ.ZERO ) GO TO 780 
          DN2 = -ZT/QZ 
          DN = (DN-ZT)/QZ 
          IF ( C(I).EQ.ZERO ) GO TO 720 
          FTRA = -DSIGN( ONE,QZ*R(I) ) 
          IF ( DABS(QZ).LT.EPS1 .AND. SIDE(I).EQ.ZERO ) FTRA = ONE 
          IF ( DABS(QZ).LT.EPS1 .AND. SIDE(I).NE.ZERO ) FTRA = SIDE(I) 
          IF ( SIDE(I).NE.ZERO ) FTRA = SIDE(I) 
          HALPHA = ONE + CONIC(I)*QZ*QZ 
          HBETA  = XT*QX + YT*QY + (ONE+CONIC(I))*ZT*QZ - R(I)*QZ 
          HGAMMA = XT*XT + YT*YT + (ONE+CONIC(I))*ZT*ZT - TWO*R(I)*ZT 
          TEMP1  = HBETA*HBETA - HALPHA*HGAMMA 
!                TEST FOR MISS IF NEG RADICAL                           
          IF ( DABS(TEMP1).LT.EPS1 ) TEMP1 = ZERO 
          IF ( TEMP1.LT.ZERO ) GO TO 780 
          TEMP = FTRA * DSQRT(TEMP1) 
          IF ( DABS(HGAMMA).LT.EPS1 ) GO TO 700 
          DENOM = HBETA + TEMP 
          IF ( DABS(DENOM).LT.EPS1 ) GO TO 710 
          DN2 = -HGAMMA/DENOM 
          GO TO 720 
  700     IF (DABS(HALPHA).LT.EPS1) GO TO 710 
          DN2 = (-HBETA + TEMP)/HALPHA 
          GO TO 720 
  710     DN2 = ZERO 
  720     IF ( DABS(DN2).LT.EPS2 ) DN2 = ZERO 
          X = XT + DN2*QX 
          Y = YT + DN2*QY 
          Z = ZT + DN2*QZ 
!                AT THIS POINT, (X,Y,Z) IS THE POINT OF INTERSECTION    
!                WITH THE CONIC (SPHERE, PLANE) IN DECENTERED, TILTED   
!                COORDINATES                                            
!                END OF SPHERICAL-CONIC TRANS EQUATIONS                 
          GO TO (800,880,730), ISW 
!                CYLINDER, TORIC TRANSFER EQUATIONS                     
  730     FTRA = -DSIGN( ONE,QZ*R(I) ) 
          IF ( DABS(QZ).LT.EPS1 .AND. SIDE(I).EQ.ZERO ) FTRA = ONE 
          IF ( DABS(QZ).LT.EPS1 .AND. SIDE(I).NE.ZERO ) FTRA = SIDE(I) 
          IF ( SIDE(I).NE.ZERO ) FTRA = SIDE(I) 
          HALPHA = (ONE+CONIC(I))*QZ*QZ + QY*QY 
          HBETA  = YT*QY + (ONE+CONIC(I))*ZT*QZ -R(I)*QZ 
          HGAMMA = YT*YT + (ONE+CONIC(I))*ZT*ZT - TWO*R(I)*ZT 
          TEMP1  = HBETA*HBETA - HALPHA*HGAMMA 
!                TEST FOR RAY MISS IF NEGATIVE RADICAL                  
          IF ( DABS(TEMP1).LT.EPS1 ) TEMP1 = ZERO 
          IF ( TEMP1.LT.ZERO ) GO TO 780 
          TEMP = FTRA*DSQRT(TEMP1) 
          IF ( DABS(HGAMMA).LT.EPS1 ) GO TO 740 
          DENOM = HBETA + TEMP 
          IF ( DABS(DENOM).LT.EPS1 ) GO TO 750 
          DN2 = -HGAMMA/DENOM 
          GO TO 760 
  740     IF (DABS(HALPHA).LT.EPS1) GO TO 750 
          DN2 = (-HBETA+TEMP)/HALPHA 
          GO TO 760 
  750     DN2 = ZERO 
  760     X = XT + DN2*QX 
          Y = YT + DN2*QY 
          Z = ZT + DN2*QZ 
!                AT THIS POINT, (X,Y,Z) IS THE INTERSECTION POINT       
!                OF THE RAY WITH A CYLINDER; NEED TO ITERATE TO         
!                COMPLETE TRANSFER TO TORIC                             
          DN = DN2 
          RNX = ZERO 
          RNY = FTRA*Y 
          RNZ = FTRA*( (ONE+CONIC(I))*Z - R(I) ) 
          DOTNN = RNY*RNY + RNZ*RNZ 
          DOTNQ = RNY*QY + RNZ*QZ 
          IF ( CX(I).EQ.ZERO ) GO TO 880 
          CALL SURTOR(I) 
          TOR34 = CX(I)/TWO 
          TOR32 = TM3*(ONE + CONIC(I)) 
!                ITERATE TO TORIC                                       
          DO 770 K=1,6 
            TOR1 = TM3*Y*Y 
            TOR2 = ONE - TOR32*TOR1 
            IF ( TOR2.LE.ZERO ) GO TO 780 
            TOR3  = DSQRT(TOR2) 
            TOR4  = TOR1/(ONE+TOR3) 
            TOR2  = TM3*Y/TOR3 
            TORF  = Z - TOR34*( X*X + Z*Z ) - TOR4*( ONE - TOR34*TOR4 ) 
            TKS   = -CX(I)*X 
            TLS   = TOR2*( CX(I)*TOR4- ONE ) 
            TMS   = ONE - CX(I)*Z 
            DN2   = DN2 - TORF/DOTNQ 
            X     = XT + QX*DN2 
            Y     = YT + QY*DN2 
            Z     = QZ*DN2 
            DOTNQ = TKS*QX + TLS*QY + TMS*QZ 
            DOTNN = TKS*TKS + TLS*TLS + TMS*TMS 
            IF ( DABS(TORF)-EPS1 ) 880,880,770 
  770     CONTINUE 
          GO TO 880 
!                FF END QUADRIC TRANSFERS                               
!                F RAY MISS ACCUMULATOR                                 
  780     IF (OPTNB.NE.0.AND.IREF.LE.1) WRITE(6,50) I 
!                RAY MISSED SURFACE IF HERE, COUNT AND GO TO  NEXT RAY  
          IF (NRAY.EQ.0.AND.IREF.GT.1) GO TO 1210 
          NMISS=NMISS+INDY 
          IF (IREF.GT.IONE) GO TO 1209 
  790     GO TO 570 
!                FF END MISS MESSAGE                                    
!                F ASPHERIC TRANSFER                                    
  800     NITER=0 
          XND=X 
          YND=Y 
          ZND=Z 
!                 2   2   2                 2   2   2   2               
!                S = X + Y  (ROT. SYMM) OR S = X + Y COS  (BILAT. SYMM) 
  810     IF (NTOR.EQ.1) SSQ=XND*XND+YND*YND 
          IF (NTOR.EQ.4 .AND. TILTY(I).NE.ZERO ) TILT=TILTY(I)*PI/180. 
          IF (NTOR.EQ.4 .AND. TILTX(I).NE.ZERO ) TILT=TILTX(I)*PI/180. 
          IF (NTOR.EQ.4) SSQ=XND*XND+YND*YND*(ONE-DCOS(TILT)*DCOS(TILT)) 
!                        2 2      1/2                                   
!                W = (1-C S (1-B))                                      
          W=ONE-TM3*TM3*SSQ*(ONE+CONIC(I)) 
          IF (W) 780,840,820 
  820     IF (W-EPS2) 840,830,830 
  830     DW=DSQRT(W) 
          W=DW 
  840     NITER=NITER+1 
          RTMP=TM3/(ONE+W) 
          DF=ZERO 
          DO 850 IA=1,4 
            IND=5-IA 
            DF=(DF+COEF(I,IND))*SSQ 
  850     CONTINUE 
          DF=(DF+RTMP)*SSQ 
!                                   2                                   
!                                 CS                4   6   8   10      
!                Z (NEW) = --------------------- +ES +FS +GS +HS        
!                                 2 2      1/2                          
!                          (1+(1-C S (1-B))   )                         
!                DELTA Z = Z(NEW) - Z(LAST)                             
          DF=DF-ZND 
          IF (NTOR.EQ.4) TILT=TILTX(I)*PI/180. 
          IF (NTOR.EQ.4) DF=DF/DSQRT(ONE-DCOS(TILT)*DCOS(TILT)) 
          DL=ZERO 
          DO 860 IA=1,4 
            IND=5-IA 
            HARG=12-2*IA 
            DL=(DL+HARG*COEF(I,IND))*SSQ 
  860     CONTINUE 
          DL=TM3+W*DL 
!                              2    4    6     8                        
!                U = -X(C+W(4ES +6FS +8GS +10HS ))                      
          U=-XND*DL 
!                              2    4    6     8                        
!                V = -Y(C+W(4ES +6FS +8GS +10HS ))                      
          V=-YND*DL 
!                        W*DELTA Z                                      
!                G  = -----------------                                 
!                 0    (Q U +Q V +Q W)                                  
!                        X    Y    Z                                    
          DELA=DF*W/(QX*U+QY*V+QZ*W) 
          XND=XND+DELA*QX 
          YND=YND+DELA*QY 
          ZND=ZND+DELA*QZ 
          IF (NITER-MAXIT) 810,870,870 
!                 2   2  2  2                                           
!                P = U +V +W                                            
  870     GSQ=U*U+V*V+W*W 
!                F = Q U+Q V+Q W                                        
!                     X   Y   Z                                         
          GNMONE=QX*U+QY*V+QZ*W 
          DOTNQ=GNMONE 
          X=XND 
          Y=YND 
          Z=ZND 
          DN2=DELA 
!                FF END ASPHERIC TRANSFER                               
!                DD END TRANSFER SECTION                                
!                D VIGNETTE TEST CODING                                 
!                AT THIS POINT (X,Y,Z) IS THE INTERSECTION OF THE       
!                RAY WITH THE ASPHERIC (ACONIC) SURFACE IN DECENTERED   
!                TILTED COORDINATES                                     
!                RAD = RADIUS (DIST FROM SURF VERTEX TO POINT OF        
!                INTERSECTION) SQUARED                                  
  880     IF (FMASK(I)) 890,960,920 
!                                                                       
  890     IF (FAKEC(I)) 895,910,900 
!                ELLIPTICAL OBSCURATION                                 
  895     XTEMP=(X-XMN(I))/XMX(I) 
          YTEMP=(Y-YMN(I))/YMX(I) 
          ELIP=(XTEMP*XTEMP)+(YTEMP*YTEMP) 
          IF (ELIP.LT.1) GO TO 950 
          GO TO 960 
!                RECTANGULAR OBSCURATION                                
  900     TEMP=ZERO 
          TEMPR=ZERO 
          IF (X.GT.XMN(I).AND.X.LT.XMX(I)) TEMP=ONE 
          IF (Y.GT.YMN(I).AND.Y.LT.YMX(I)) TEMPR=ONE 
          IF (IREF.GT.1) GO TO 960 
          IF (TEMP.EQ.ONE.AND.TEMPR.EQ.ONE) GO TO 950 
          GO TO 960 
!                CIRCULAR OBSCURATION                                   
!                MASK HAS CENTER AT (XMN(II),YMN(I)), RADIUS = FMASK(I) 
  910     RAD=(X-XMN(I))*(X-XMN(I))+(Y-YMN(I))*(Y-YMN(I)) 
          TEMP=FMASK(I)*FMASK(I) 
          IF (IREF.GT.1) GO TO 960 
          IF (RAD-TEMP) 950,960,960 
!                                                                       
  920     IF (FAKEC(I)) 925,940,930 
!                ELLIPTICAL CLEAR APERTURE                              
  925     XTEMP=(X-XMN(I))/XMX(I) 
          YTEMP=(Y-YMN(I))/YMX(I) 
          ELIP=(XTEMP*XTEMP)+(YTEMP*YTEMP) 
          IF (ELIP.GT.1) GO TO 950 
          GO TO 960 
!                RECTANGULAR CLEAR APERATURE                            
  930     TEMP=ZERO 
          TEMPR=ZERO 
          IF (X.LT.XMN(I).OR.X.GT.XMX(I)) TEMP=ONE 
          IF (Y.LT.YMN(I).OR.Y.GT.YMX(I)) TEMPR=ONE 
          IF (IREF.GT.1) GO TO 960 
          IF (TEMP.EQ.ONE.OR.TEMPR.EQ.ONE) GO TO 950 
          GO TO 960 
!                CIRCULAR CLEAR APERATURE                               
!                MASK HAS CENTER AT (XMN(II),YMN(I)), RADIUS = FMASK(I) 
  940     RAD=(X-XMN(I))*(X-XMN(I))+(Y-YMN(I))*(Y-YMN(I)) 
          TEMP=FMASK(I)*FMASK(I) 
          IF (IREF.GT.IONE) GO TO 960 
          IF (RAD-TEMP) 960,960,950 
!                EEEEE  END VIGNETTE ACCOUNTING                         
  950     IF (OPTNB.NE.0) WRITE(6,70) I 
!                COUNT VIGNETTE AND GO TO  NEXT RAY                     
          IF (IREF.GT.IONE) GO TO 960 
          IF (OPTNB.NE.0.OR.OPTNC.NE.0) WRITE(6,20) I,XYZ,QT 
          NVIGN=NVIGN+INDY 
          GO TO 790 
!                DD END VIGNETTE CODING                                 
!                D REFRACTION (MODIFICATION OF COSINES) EQUATIONS       
!                FIND ZN, Z COMPONENT OF NON-UNIT NORMAL VECTOR TO      
!                THE SURFACE AT POINT (X,Y,Z) - MEANINGLESS QUANTITY    
!                FOR PLANAR SURFACE                                     
!                GNU = RATIO OF OLD INDEX TO NEW INDEX                  
  960     GNU    = TM2/TM1 
          GO TO (980,970,980), ISW 
  970     IF ( C(I).EQ.ZERO .AND. RDSPAC(I).EQ.ZERO ) GO TO 1060 
          RNX = FTRA*X 
          RNY = FTRA*Y 
          RNZ = FTRA*((ONE + CONIC(I))*Z - R(I)) 
!                NORMAL DOT NORMAL                                      
          DOTNN = RNX*RNX + RNY*RNY + RNZ*RNZ 
!                NORMAL DOT COSINE VECTOR                               
          DOTNQ = QX*RNX + QY*RNY + QZ*RNZ 
!                DIFFRACTION GRATING SURFACE                            
  980     IF (RDSPAC(I)) 1090,990,1090 
!                ASPHERIC SURFACE                                       
  990     GO TO (1000,1020,1020), ISW 
!                E ASPHERIC REFRACTION                                  
!                ASPHERIC REFRACTION EQUATIONS                          
!                             2      2                                  
!                 /  ( 2     N      N   2)1/2                           
!                F = (P (1- ---) + --- F )                              
!                    (       2      2    )                              
!                           N      N                                    
!                            -1     -1                                  
 1000     GN = (GNMONE*GNMONE-GSQ)*GNU*GNU+GSQ 
!                TEST FOR TOTAL INTERNAL REFLECTION (BREWSTERS ANGLE)   
          IF (GN.LE.ZERO) GO TO 1010 
          RTGN = FREF(I)*DSQRT(GN) 
          GN=RTGN 
!                     1 ( /   N   )                                     
!                G = ---(F - --- F)                                     
!                     2 (    N    )                                     
!                    P        -1                                        
          DP = (GN-GNU*GNMONE)/GSQ 
!                          N                                            
!                NEW Q  = --- Q  + GU                                   
!                     X   N    X                                        
!                          -1                                           
          QX = GNU*QX + U*DP 
!                          N                                            
!                NEW Q  = --- Q  + GV                                   
!                     Y   N    Y                                        
!                          -1                                           
          QY = GNU*QY + V*DP 
!                          N                                            
!                NEW Q  = --- Q  +  GW                                  
!                     Z   N    Z                                        
!                          -1                                           
          QZ = GNU*QZ + W*DP 
!                HERE QX,QY,QZ ARE NEW DIRECTION COSINES OF RAY IN      
!                DECENTERED, TILTED COORD AFTER REFRACTION THRU         
!                ASPHERIC                                               
          GO TO 1170 
!                EE END ASPHERIC REFRACTION                             
!                E INTERNAL REFLECTION ACCOUNTING                       
 1010     IF (OPTNB.NE.0) WRITE(6,60) I 
          IF (NRAY.EQ.0.AND.IREF.GT.1) GO TO 1210 
          NREFL=NREFL+INDY 
          GO TO 790 
!                EE END INT. REFL. ACTG.                                
!                E QUADRIC (SPHERE, PLANE) TRANSFER EQ.S                
 1020     IF ( DOTNN.EQ.ZERO .OR. DOTNQ.EQ.ZERO ) GO TO 1060 
          TEMP=(DOTNN/(DOTNQ*DOTNQ))*(ONE-GNU*GNU)+GNU*GNU 
          GMU=ZERO 
!                TEST FOR TOTAL INTERNAL REFLECTION                     
          IF (TEMP) 1010,1030,1030 
 1030     GMU = FREF(I)*DSQRT(TEMP) 
          GMU = (GMU-GNU)*DOTNQ/DOTNN 
          GO TO (1050,1050,1040), ISW 
 1040     IF ( CX(I).EQ.ZERO ) GO TO 1045 
!             GENERAL TORIC REFRACTION EQUATIONS                        
          QX = GNU*QX + GMU*TKS 
          QY = GNU*QY + GMU*TLS 
          QZ = GNU*QZ + GMU*TMS 
          GO TO 1170 
!             CYLINDER REFRACTION EQUATIONS                             
 1045     QX = GNU*QX + GMU*RNX 
          QY = GNU*QY + GMU*RNY 
          QZ = GNU*QZ + GMU*RNZ 
          GO TO 1170 
!             CONIC,ACONIC REFRACTION EQUATIONS                         
 1050     QX = GNU*QX + GMU*RNX 
          QY = GNU*QY + GMU*RNY 
          QZ = GNU*QZ + GMU*RNZ 
          GO TO 1170 
!             PLANE REFRACTION EQUATIONS                                
 1060     TEMP = ONE - GNU*GNU*(ONE-QZ*QZ) 
          GMU = ZERO 
          IF (TEMP) 1010,1080,1070 
 1070     GMU = DSIGN(ONE,QZ) * FREF(I) * DSQRT(TEMP) 
 1080     GMU = GMU - GNU*QZ 
          QX = GNU*QX 
          QY = GNU*QY 
          QZ = GNU*QZ + GMU 
          GO TO 1170 
!                EE END QUADRIC REFRACTION                              
!                E DIFFRACTION GRATING REFRACTION CODING                
!                IS THIS A PLANE GRATING                                
 1090     IF ( C(I) ) 1110,1100,1110 
!                YES, SET NORMAL VECTOR AND PRODUCTS                    
 1100     XND=ZERO 
          YND=ZERO 
          ZND=ONE 
          DOTNN=ONE 
          DOTNQ=QZ 
          GO TO 1120 
!                NOT PLANE                                              
 1110     XND = RNX 
          YND = RNY 
          ZND = RNZ 
!                RULING DIRECTION - RDSPAC LT 0 IF X, GT 0 IF Y         
 1120     IF (RDSPAC(I)) 1130,990,1140 
!                X GRATING, FIND CAPITAL LAMDA                          
 1130     SSQ=XND*XND+ZND*ZND 
          V=-ONE/DSQRT(ONE+YND*YND/SSQ) 
          U=-XND*V*YND/SSQ 
          W=-ZND*V*YND/SSQ 
!                USER SHALL SPECIFY THE SPACING OF DIFFRACTION GRATINGS 
!                AS NOT FORESHORTENED FROM NORMAL TO THE GRATING        
          RDLAM = ORDN(I,LAMDA)*WAVL(LAMDA)*DABS(V)/DABS(TM1*RDSPAC(I)) 
          GO TO 1150 
!                Y GRATING, FIND CAPITAL LAMDA                          
 1140     SSQ = YND*YND+ZND*ZND 
          U   = ONE/DSQRT(ONE + XND*XND/SSQ) 
          V   = -YND*U*XND/SSQ 
          W   = -ZND*U*XND/SSQ 
          RDLAM = ORDN(I,LAMDA)*WAVL(LAMDA)*U/DABS(TM1*RDSPAC(I)) 
!                SOLVE QUADRATIC EQUATION FOR GAMMA                     
 1150     B1 = (GNU*GNU - ONE + RDLAM*RDLAM - TWO*GNU*RDLAM*            &
     &       (U*QX + V*QY + W*QZ) )/DOTNN                               
          A = GNU*DOTNQ/DOTNN 
          ASQ = A*A 
          IF (B1-ASQ) 1160,1160,1010 
 1160     RDGAM = FREF(I) * DSIGN(ONE,QZ) * DSQRT(ASQ-B1) 
          RDGAM = -(A-RDGAM) 
!                FIND NEW COSINES FOR RAY                               
          QX = GNU*QX - RDLAM*U + RDGAM*XND 
          QY = GNU*QY - RDLAM*V + RDGAM*YND 
          QZ = GNU*QZ - RDLAM*W + RDGAM*ZND 
!                EE END OF DIFFRACTION GRATING CODING                   
!                DD END OF REFRACTION EQUATIONS                         
!                D RE-TRANSFORMATION TO OPTICAL AXIS COORDINATES AND    
!                D DIRECTION                                            
!                CONVERT BACK THRU OPPOSITE ROTATION                    
!      FAKEB=0 NO DEROTATE NO DETRANSLATE                               
!      FAKEB=1 TRANSLATES ONLY BACK TO ORIGINAL COORDINATES             
!      FAKEB=2 DEROTATE ONLY BACK TO ORIGINAL COORDINATES               
!      FAKEB=3 DEROTATES AND TRANSLATES BACK TO ORIGINAL COORDINATES    
 1170     IF ( FAKEB(I).EQ.ZERO ) GO TO 1190 
          IF ( FAKEB(I).EQ.ONE  )  GO TO 1180 
          IF ( FAKEB(I).EQ.TWO .OR. FAKEB(I).EQ.THREE ) GO TO 1175 
          GO TO 1190 
 1175     CALL VECMA (XYZ,XYZ) 
          CALL VECMA (QT,QT) 
          IF ( FAKEB(I).EQ.TWO )  GO TO 1190 
!                RE-CENTER COORD                                        
 1180     X = X + XDISP(I) 
          Y = Y + YDISP(I) 
!                                                                       
!                SIZEB = TRAVEL FROM SURFACE I-1 TO I                   
!                      = OPTICAL PATH LENGTH                            
 1190     SIZEB = (DN1 - DN + DN2)*TM2 
!             PRINT INFO THIS SURF                                      
          IF (OPTNC.NE.0.AND.IREF.LE.1) WRITE(6,20) I,XYZ,QT,SIZEB 
!                SET PREVIOUS SURFACE PARAMS                            
          PX = X 
          PY = Y 
          PZ = Z 
          PN = TM1 
          PT = T(I) 
!                END OF TRACE LOOP                                      
 1200   CONTINUE 
!                DD END TRANSFORMATIONS                                 
!                CC END SURFACE CALCULATIONS                            
!                C IMAGE PLANE(S) ANALYSIS AND ACTG. FOR RAY            
!                                                                       
!                PERFORM REFERENCE SURFACE ITERATIONS                   
!                                                                       
!                VARIABLES:                                             
!                                                                       
!                   ORIGX,ORIGY: ORIGINAL X,Y COORDS                    
!                                 AT ENTRANCEPUPIL                      
!                   X2,Y2: NEW X,Y COORDS AT ENTPUP                     
!                   DLTREF: DELTA X AND Y USED TO CALCULATE             
!                           PARTIAL DERIVATIVES                         
!                   TRACX1,TRACY1: REF SURF X,Y FOR ENTPUP X,Y          
!                   TRACX2,TRACY2: REF SURF X,Y FOR ENTPUP X+DLTREF,Y   
!                   TRACX3,TRACY3: REF SURF X,Y FOR ENTPUP X,Y+DLTREF   
!                   DESX,DESY: DESIRED X,Y AT REFERENCE SURFACE         
!                   ITREAD:  0 = NO REF SURF ITERATIONS PERFORMED       
!                            1 = FINDING PARTIAL DERIVATIVES            
!                            2 = FINDING NEW X,Y COORDS AT REF SURF     
!                   IWHICH:  1 = FINDING NEW X,Y COORDS AT REF SURF     
!                            2 = FINDING PARTIALS W.R.T. X              
!                            3 = FINDING PARTIALS W.R.T. Y              
!                   PX,PY: AFTER LABEL 1200 = X,Y COORDS AT REF SURF    
!                                                                       
!                FORMULAS:                                              
!                                                                       
!                   F1(X,Y) = PX                                        
!                   F2(X,Y) = PY                                        
!                   G1(X,Y) = PX-DESX                                   
!                   G2(X,Y) = PY-DESY                                   
!                                                                       
!                           (G2(X,Y)*(DF1/DY))-(G1(X,Y)*(DF2/DY))       
!                   CHX =  --------------------------------------       
!                          ((DF1/DX)*(DF2/DY)-(DF2/DX)*(DF1/DY))        
!                                                                       
!                           (G2(X,Y)*(DF1/DX)-(G2(X,Y)*(DF2/DX))        
!                   CHY =  ----------------------------------------     
!                          ((DF1/DY)*(DF2/DX)-(DF2/DY)*(DF1/DX))        
!                                                                       
        IF (ITREAD.EQ.0) GO TO 1216 
        NITRAT=NITRAT+1 
        IF (NITRAT.GT.48) GO TO 1209 
        IF (ITREAD.EQ.2) GO TO 1204 
        IF (IWHICH.EQ.1) GO TO 1201 
        IF (IWHICH.EQ.2) GO TO 1202 
        IF (IWHICH.EQ.3) GO TO 1203 
!                                                                       
 1201   X2=ORIGX 
        Y2=ORIGY 
        TRACX1=PX 
        TRACY1=PY 
        IF (DABS(PX-DESX).LE.EPSLON.AND.DABS(PY-DESY).LE.EPSLON)        &
     &     GO TO 1209                                                   
        ONETHO=1.D+3 
        DLTREF=1.0/ONETHO 
        IF (FMASK(1).GT.0) DLTREF=FMASK(1)/ONETHO 
        IF (ORIGX.GE.0) X2=ORIGX-DLTREF 
        IF (ORIGX.LT.0) X2=ORIGX+DLTREF 
        PX=X2 
        PY=Y2 
        PZ=ZHOLD 
        IWHICH=2 
        DX=X2-THISHX 
        DY=Y2-THISHY 
        DENOM=DSQRT(DX*DX+DY*DY+DZ*DZ) 
        QX=DX/DENOM 
        QY=DY/DENOM 
        QZ=DZ/DENOM 
        GO TO 585 
!                                                                       
 1202   TRACX2=PX 
        TRACY2=PY 
        X21=X2 
        IF (DABS(PX-DESX).LE.EPSLON.AND.DABS(PY-DESY).LE.EPSLON)        &
     &     GO TO 1209                                                   
        IF (ORIGY.GE.0) Y2=ORIGY-DLTREF 
        IF (ORIGY.LT.0) Y2=ORIGY+DLTREF 
        X2=ORIGX 
        PX=X2 
        PY=Y2 
        PZ=ZHOLD 
        IWHICH=3 
        DX=X2-THISHX 
        DY=Y2-THISHY 
        DENOM=DSQRT(DX*DX+DY*DY+DZ*DZ) 
        QX=DX/DENOM 
        QY=DY/DENOM 
        QZ=DZ/DENOM 
        GO TO 585 
!                                                                       
 1203   TRACX3=PX 
        TRACY3=PY 
        Y21=Y2 
        IF (DABS(PX-DESX).LE.EPSLON.AND.DABS(PY-DESY).LE.EPSLON)        &
     &     GO TO 1209                                                   
        IF (ORIGX.GE.0) X2=ORIGX-DLTREF 
        IF (ORIGX.LT.0) X2=ORIGX+DLTREF 
        IF (ORIGY.GE.0) Y2=ORIGY-DLTREF 
        IF (ORIGY.LT.0) Y2=ORIGY+DLTREF 
        CHX=X2-ORIGX 
        CHY=Y2-ORIGY 
        PX=X2 
        PY=Y2 
        PZ=ZHOLD 
        ITREAD=2 
        IWHICH=1 
        DX=X2-THISHX 
        DY=Y2-THISHY 
        DENOM=DSQRT(DX*DX+DY*DY+DZ*DZ) 
        QX=DX/DENOM 
        QY=DY/DENOM 
        QZ=DZ/DENOM 
        QXHOLD=QX 
        QYHOLD=QY 
        QZHOLD=QZ 
        GO TO 585 
!                                                                       
 1204   IF (IWHICH.EQ.1) GO TO 1206 
        IF (IWHICH.EQ.3) GO TO 1208 
!                                                                       
 1206   IF (DABS(PX-DESX).LE.EPSLON.AND.DABS(PY-DESY).LE.EPSLON)        &
     &     GO TO 1209                                                   
        FPR1X=(TRACX2-TRACX1)/CHX 
        FPR2X=(TRACY2-TRACY1)/CHX 
        FPR1Y=(TRACX3-TRACX1)/CHY 
        FPR2Y=(TRACY3-TRACY1)/CHY 
        G1XY=PX-DESX 
        G2XY=PY-DESY 
        UUJ=((FPR1X*FPR2Y)-(FPR2X*FPR1Y)) 
        UUK=((FPR1Y*FPR2X)-(FPR2Y*FPR1X)) 
        IF (UUJ.GT.0.AND.UUJ.LT.1.D-37) UUJ=1.D-37 
        IF (UUJ.LE.0.AND.UUJ.GT.-1.D-37) UUJ=-1.D-37 
        IF (UUK.GT.0.AND.UUK.LT.1.D-37) UUK=1.D-37 
        IF (UUK.LE.0.AND.UUK.GT.-1.D-37) UUK=-1.D-37 
        CHX=((G2XY*FPR1Y)-(G1XY*FPR2Y))/UUJ 
        CHY=((G2XY*FPR1X)-(G1XY*FPR2X))/UUK 
!                                                                       
 1208   IF (DABS(PX-DESX).LE.EPSLON.AND.DABS(PY-DESY).LE.EPSLON)        &
     &     GO TO 1209                                                   
        X2=ORIGX+CHX 
        Y2=ORIGY+CHY 
        PX=X2 
        PY=Y2 
        PZ=ZHOLD 
        IWHICH=1 
        ORIGX=X2 
        ORIGY=Y2 
        ITREAD=1 
        DX=X2-THISHX 
        DY=Y2-THISHY 
        DENOM=DSQRT(DX*DX+DY*DY+DZ*DZ) 
        QX=DX/DENOM 
        QY=DY/DENOM 
        QZ=DZ/DENOM 
        QXHOLD=QX 
        QYHOLD=QY 
        QZHOLD=QZ 
        GO TO 585 
!                                                                       
!               HERE IF ITERATED COORDINATES FOUND                      
 1209   IF (NRAY.EQ.0.AND.NITRAT.GT.48) GO TO 1210 
        IF (NRAY.EQ.0) GO TO 1211 
        DX=X2-THISHX 
        DY=Y2-THISHY 
        DENOM=DSQRT(DX*DX+DY*DY+DZ*DZ) 
        QX=DX/DENOM 
        QY=DY/DENOM 
        QZ=DZ/DENOM 
        WRITE(50) X2,Y2,ZHOLD,QX,QY,QZ,THISHY,THISHX,DESX,DESY 
        ITREAD=IONE 
        IWHICH=IONE 
        GO TO 570 
!                                                                       
!               HERE IF CHIEF RAY CANNOT BE TRACED                      
 1210   WRITE(6,213) THISY, THISX 
        X2DEC=0 
        Y2DEC=0 
        IREF=-1 
        ITREAD=0 
        KFLAG=1 
        NSURF=NNSURF 
        GO TO 570 
!               CHIEF RAY CAN BE TRACED                                 
!               PERFORM ITERATIONS NOW                                  
 1211   NRAY=0 
        X2DEC=X2 
        Y2DEC=Y2 
        ITREAD=IONE 
        IWHICH=IONE 
        KFLAG=0 
        NSURF=NNSURF 
        GO TO 581 
!               WRITE RAY DATA TO UNIT 20                               
 1216   IF (IPR20.EQ.IZERO) GO TO 1218 
        WRITE(20,212) PX,PY,PZ,QX,QY,QZ,THISHY,THISHX 
!               FIND FINAL ANGLE TANGENTS                               
 1218   IF (DABS(QZ).GT.ZERO) GO TO 1219 
        WRITE(6,210) 
        GO TO 1900 
 1219   TM1=QX/QZ 
        TM2=QY/QZ 
        PRINAN=PRINAN+TM2 
        IF (UFLAG.LT.FOUR) GO TO 1220 
        FPLANE=ZERO 
        NPLANE=1 
        FOCL=ONE 
        NFOC=0 
 1220   IF (OPTNB) 1230,1240,1230 
!             PRINT TANGENTS, IMAGE PLANES BREAKDOWN HEADING            
 1230   WRITE(6,90) TM1,TM2 
        WRITE(6,80) 
!             TRANSFER TO IMAGE SURFACES VIA CLOSEST POINT METHOD       
 1240   PT = T(NSURF) + FPLANE*DELIMP 
        ZLOC = PT 
        RECX = ZERO 
        RECY = ZERO 
        IF ( FAKEB(NSURF).EQ.ONE .OR. FAKEB(NSURF).EQ.THREE )           &
     &       RECX = XDISP(NSURF)                                        
        IF ( FAKEB(NSURF).EQ.ONE .OR. FAKEB(NSURF).EQ.THREE )           &
     &       RECY = YDISP(NSURF)                                        
        IF ( FAKEA.GT.ZERO ) DN = ZERO 
        I = IZERO 
 1250   I = I + IONE 
        IF ( I.GT.NPLANE ) GO TO 1345 
        IF ( UFLAG.NE.FOUR ) GO TO 1260 
        X = QX/QZ 
        Y = QY/QZ 
        GO TO 1300 
 1260   IF ( UFLAG.NE.FIVE ) GO TO 1270 
        IF ( C(NSURF).EQ.ZERO ) X = DACOS(QZOLD) 
        IF ( C(NSURF).NE.ZERO ) X = DACOS(DOTNQ*DABS(C(NSURF))) 
        Y = ZERO 
        GO TO 1300 
 1270   DN = (PT-PZ)*QZ - PY*QY - PX*QX 
        DN1 = DN 
        XT = X + DN*QX - RECX 
        YT = Y + DN*QY - RECY 
        ZT = Z + DN*QZ - PT 
!                END OF CLOSEST-POINT TRANSFER                          
!                TRANSFER RAY TO IMAGE SURFACE                          
!                CONIC, SPHERIC, PLANE TRANSFER EQUATIONS               
        IF ( CVIMG.EQ.ZERO .AND. QZ.EQ.ZERO ) GO TO 780 
        DN2 = -ZT/QZ 
        IF ( CVIMG.EQ.ZERO ) GO TO 1290 
        FTRA = -DSIGN( ONE,QZ*RADIMG ) 
        HALPHA = ONE + CONIMG*QZ*QZ 
        HBETA  = XT*QX + YT*QY + (ONE+CONIMG)*ZT*QZ - RADIMG*QZ 
        HGAMMA = XT*XT + YT*YT + (ONE+CONIMG)*ZT*ZT - TWO*RADIMG*ZT 
        TEMP1  = HBETA*HBETA - HALPHA*HGAMMA 
!                TEST FOR MISS IF NEG RADICAL                           
        IF ( DABS(TEMP1).LT.EPS1 ) TEMP1 = ZERO 
        IF ( TEMP1.LT.ZERO ) GO TO 780 
        TEMP = FTRA * DSQRT(TEMP1) 
        IF ( DABS(HGAMMA).LT.EPS1 ) GO TO 1280 
        DENOM = HBETA + TEMP 
        IF ( DABS(DENOM).LT.EPS1 ) GO TO 1285 
        DN2 = -HGAMMA/DENOM 
        GO TO 1290 
 1280   DN2 = (HBETA-TEMP)/HALPHA 
        GO TO 1290 
 1285   DN2 = ZERO 
 1290   X = XT + DN2*QX 
        Y = YT + DN2*QY 
        Z = ZT + DN2*QZ 
 1300   XSUM(I) = XSUM(I) + X*INDX 
        YSUM(I) = YSUM(I) + Y*INDY 
        XSUMSQ(I) = XSUMSQ(I) + X*X*INDY 
        YSUMSQ(I) = YSUMSQ(I) + Y*Y*INDY 
        XUS = X + RECX 
        YUS = Y + RECY 
!             NO PRINT UNLESS PRIME IMAGE PLANE                         
        IF ( I.NE.NIMG ) GO TO 1340 
!                WRITE X,Y ON SPOTPLOT TAPE                             
!                IF FIRST RAY THROUGH, INITIALIZE PLOT MAX, MIN         
        IF (NTHRU.NE.0) GO TO 1310 
        IF (IALLPL.EQ.1) GO TO 1310 
        XMINS = XUS 
        XMAXS = XUS 
        YMINS = YUS 
        YMAXS = YUS 
 1310   FAKEX = XUS 
!                ACCUMULATE MAX AND MIN X AND Y OVER IMAGE POINTS       
        IF ( YUS.GT.YMAXS ) YMAXS = YUS 
        IF ( YUS.LT.YMINS ) YMINS = YUS 
 1320   IF ( FAKEX.GT.XMAXS ) XMAXS = FAKEX 
        IF ( FAKEX.LT.XMINS ) XMINS = FAKEX 
        FAKEX = -FAKEX 
!                IF MODE = 0 CONSIDER + AND -X                          
        IF ( FAKEX.EQ.-XUS .AND. XUS.NE.ZERO .AND. IMODE.EQ.0 )         &
     &     GO TO 1320                                                   
        IF ( IPLTPR.NE.0 .OR. OPTNE.NE.0 .OR. OPTNF.NE.0 )              &
     &     WRITE(34) XUS,YUS                                            
        IF ( IPLTPR.NE.0 .AND. IALLPL.EQ.1) GO TO 1325 
        GO TO 1327 
 1325   IJK=IJK+1 
        IF (IJK.LE.7600) GO TO 1326 
        WRITE(6,218) 
        IALLPL=0 
        GO TO 1327 
 1326   XKALL(IJK)=XUS 
        YKALL(IJK)=YUS 
 1327   IF ( OPTNA.EQ.0 ) GO TO 1330 
!             PRINT X,Y,Z                                               
        WRITE(6,40) XUS,YUS,Z 
!             IF MODE = 0 WRITE, PRINT -X,Y,Z                           
 1330   IF ( IMODE.NE.0 ) GO TO 1340 
        FAKEX = -XUS 
        IF ( IPLTPR.NE.0 .OR. OPTNE.NE.0 .OR. OPTNF.NE.0 )              &
     &     WRITE(34) FAKEX,YUS                                          
        IF ( IPLTPR.NE.0 .AND. IALLPL.EQ.1) GO TO 1335 
        GO TO 1337 
 1335   IJK=IJK+1 
        IF (IJK.LE.7600) GO TO 1336 
        WRITE(6,218) 
        IALLPL=0 
        GO TO 1337 
 1336   XKALL(IJK)=FAKEX 
        YKALL(IJK)=YUS 
 1337   IF ( OPTNA.EQ.0 ) GO TO 1340 
        WRITE(6,40) FAKEX,YUS,Z 
 1340   IF ( OPTNB.NE.0 ) WRITE(6,100) ZLOC,XUS,YUS,Z 
        PT = DELIMP 
        ZLOC = ZLOC + DELIMP 
        GO TO 1250 
!                UPDATE NUMBER OF RAYS THRU                             
 1345   NTHRU=NTHRU+INDY 
        GO TO 570 
!                CC END IMAGE ANAL. FOR RAY                             
!                BB END RAY PROC FOR THIS COLOR                         
!                B IMAGE ANALYSIS FOR ALL RAYS, THIS HEIGHT, THIS COLOR 
!             PRINT FINAL IMAGE PLANE ANALYSIS HERE                     
!                FLOAT NUMBER RAYS THRU                                 
 1350   TIMES=NTHRU 
        IF (NTHRU+NMISS+NREFL+NVIGN.EQ.0) GO TO 1360 
!             PRINT ACCOUNT                                             
        IF (LATYPE.EQ.6) WRITE(6,105) LAMDA,NTHRU,NREFL,NMISS,NVIGN 
        IF (LATYPE.EQ.6) GO TO 1360 
        WRITE(6,110) LAMDA,HYLAST,HXLAST,NTHRU,NREFL,NMISS,NVIGN 
!                CLOSE PLOT                                             
!                NEW HEIGHT TO HOLDER                                   
 1360   HYSAVE=HYLAST 
        HXSAVE=HXLAST 
        HYLAST=THISHY 
        HXLAST=THISHX 
        IF (IREF.GT.IONE.OR.KFLAG.EQ.1) GO TO 1870 
!                SKIP IF NO RAYS THRU                                   
        IF (NTHRU.LT.1) GO TO 1870 
        WRITE(6,120) 
!                POSITION OF PLANE                                      
        Z=T(NSURF)+FPLANE*DELIMP 
        DO 1390 I=1,NPLANE 
!                FIND AVERAGES                                          
          XAVG=XSUM(I)/TIMES 
          YAVG=YSUM(I)/TIMES 
          RMSX=ZERO 
          RMSY=ZERO 
          SPOT=ZERO 
!                IF ONE RAY ONLY RMS =0                                 
          IF (TIMES.EQ.ONE) GO TO 1370 
!                CALCULATE ROOT MEAN SQUARES                            
!                NOTE STD DEV OF X AND Y ARE MISNAMED RMSX, RMSY        
          RMSX=XSUMSQ(I)/TIMES-XAVG*XAVG 
          IF ( DABS(RMSX).LT.EPS1 ) RMSX = DABS(RMSX) 
          RMSX=DSQRT(RMSX) 
!                                                                       
          RMSY=YSUMSQ(I)/TIMES-YAVG*YAVG 
          IF ( DABS(RMSY).LT.EPS1 ) RMSY = DABS(RMSY) 
          RMSY=DSQRT(RMSY) 
!                SPOT SIZE = MEASURE OF IMAGE SIZE                      
          SPOT=DSQRT(RMSX*RMSX+RMSY*RMSY) 
          IF (IFOC.LT.0) BLOB(I)=RMSX 
          IF (IFOC.EQ.0) BLOB(I)=SPOT 
          IF (IFOC.GT.0) BLOB(I)=RMSY 
          ZF(I)=Z 
          XAVG=XAVG+RECX 
          YAVG=YAVG+RECY 
 1370     WRITE(6,100) Z,XAVG,YAVG,RMSX,RMSY,SPOT 
!                RESTORE ACCUMULATORS                                   
!                SAVE X, Y AVERAGES, SPOT RADIUS IF PRIME IMAGE PLANE   
          IF (I.NE.NIMG) GO TO 1380 
          AVGX=XAVG 
          AVGY=YAVG 
          SPOTP=SPOT 
          IF (IALLPL.EQ.0) GO TO 1380 
          IHITE=IHITE+1 
          AHITE=IHITE 
          AVXALL=(((IHITE-1)*AVXALL)+AVGX)/AHITE 
          AVYALL=(((IHITE-1)*AVYALL)+AVGY)/AHITE 
          XALLSQ=XALLSQ+XSUMSQ(I) 
          YALLSQ=YALLSQ+YSUMSQ(I) 
 1380     XSUM(I)=ZERO 
          YSUM(I)=ZERO 
          XSUMSQ(I)=ZERO 
          YSUMSQ(I)=ZERO 
          Z=Z+DELIMP 
 1390   CONTINUE 
        IF (NFOC.EQ.0) GO TO 1400 
        NUMP=NPLANE 
        CALL FOCUS (ZF,BLOB,NUMP,FID) 
        NFOC=0 
        LFOC=0 
        T(NSURF)=FID 
        REWIND 34 
        GO TO 300 
!                                                                       
!             REWIND 34 BEFORE RED OR MTF CALCULATION                   
!                                                                       
 1400   REWIND 34 
!                                                                       
!   COMPUTE RED IF NEEDED; AVOID RED CALC. IF POSSIBLE                  
!                                                                       
!   RED  NEEDED IF    IPLTPR IS ON,                                     
!         OR          IOPB    "  ",                                     
!         OR          IOPC    "  ",                                     
!         OR          OPTNE   "  ",                                     
!         OR          OPTNF   "  ".                                     
!                                                                       
        IF( IPLTPR .NE. 0) GO TO 1540 
        IF( IOPB   .NE. 0) GO TO 1540 
        IF( IOPC   .NE. 0) GO TO 1540 
        IF( OPTNE  .NE. 0) GO TO 1540 
        IF( OPTNF  .NE. 0) GO TO 1540 
!                                                                       
!   SKIP RED CALCULATION (NOT NEEDED)                                   
        GO TO 1870 
!                                                                       
!            RED CALCULATION                                            
!            GET DEVIATIONS OF RAYS FROM CENTROID.                      
!       RP(I) - RMS DEVIATION OF ITH RAY                                
!                                                                       
 1540   DO 1560 I=1,NTHRU 
        READ (34) XK(I),YK(I) 
        RP(I)=DSQRT((XK(I)-AVGX)*(XK(I)-AVGX)+(YK(I)-AVGY)*             &
     &     (YK(I)-AVGY))                                                
 1560   CONTINUE 
!                                                                       
!     STORE RADII IN ASCENDING ORDER (SEE USER'S GUIDE P. 4-39)         
!                                                                       
        DO 1620 I=1,NTHRU 
          TEMP=RP(I) 
          TEMPR=TEMP 
          ITEMP=I 
          DO 1610 J=I,NTHRU 
            IF (TEMP.LT.RP(J)) GO TO 1610 
            TEMP=RP(J) 
            ITEMP=J 
 1610     CONTINUE 
          RP(I)=TEMP 
          RP(ITEMP)=TEMPR 
 1620   CONTINUE 
        IF (IOPA.NE.0.AND.IALLPL.EQ.0) CALL FNDSPT 
!                                                                       
        CALL RED 
!                                                                       
!                IF NO MTF NEEDED, WE ARE THROUGH                       
!                                                                       
 1750   IF ( IOPC.EQ.0 .AND. OPTNF.EQ.0 ) GO TO 1860 
!                                                                       
        CALL MTF 
!                                                                       
 1860   REWIND 34 
!                                                                       
!                RESET VARIABLES, RETURN TO NEXT PASS                   
!                                                                       
 1870   NRAY=0 
        NTHRU=0 
        NREFL=0 
        NMISS=0 
        NVIGN=0 
        KFLAG=0 
        ORIGX=ORIGX-X2DEC 
        ORIGY=ORIGY-Y2DEC 
        IF (IIREF.GT.1.AND.IREF.EQ.1) IREF=IIREF 
        IF (OPTNA.NE.0.AND.IEOF.NE.1) WRITE(6,30) 
        GO TO (1875,580), IEOF 
 1875 IF (IOPA.EQ.0.OR.IALLPL.EQ.0) GO TO 1880 
!               HERE IF MULTIPLE SPOT PLOT DESIRED                      
      AVGX=AVXALL 
      AVGY=AVYALL 
      TIMES=IJK 
!                CALCULATE ROOT MEAN SQUARES                            
!                NOTE STD DEV OF X AND Y ARE MISNAMED RMSX, RMSY        
      RMSX=XALLSQ/TIMES-AVGX*AVGX 
      IF ( DABS(RMSX).LT.EPS1 ) RMSX = DABS(RMSX) 
      RMSX=DSQRT(RMSX) 
!                                                                       
      RMSY=YALLSQ/TIMES-AVGY*AVGY 
      IF ( DABS(RMSY).LT.EPS1 ) RMSY = DABS(RMSY) 
      RMSY=DSQRT(RMSY) 
!                SPOT SIZE = MEASURE OF IMAGE SIZE                      
      SPOTP=DSQRT(RMSX*RMSX+RMSY*RMSY) 
      CALL GRAPH(2,4) 
 1880 END DO 
      IF (IREF.LE.0) GO TO 1890 
!               CLOSE UNIT 50                                           
      IREF=IZERO 
      GARB=-2 
      WRITE(50) (GARB, I=1,10) 
      END FILE 50 
      REWIND 50 
      LAMDA=ICOL(1) 
      NRAY=0 
      NTHRU=0 
      NREFL=0 
      NMISS=0 
      NVIGN=0 
      NSURF=NNSURF 
      GO TO 550 
 1890 CONTINUE 
!                BB END MACRO IMAGE ANALYSIS                            
!                AA END PROCESSING THIS COLOR                           
!                RESET IREF TO REFERENCE SURFACE NUMBER                 
 1900 IREF=IIREF 
      IF (FAKEA.EQ.ZERO) GO TO 1910 
      FAKEA=ZERO 
      T(NSURF)=THCK 
      NSURF=NNSURF 
      DELIMP=DLMP 
      FPLANE=FPLN 
      NPLANE=NPLN 
      IPR20=0 
!                                                                       
      GARB=-2 
      WRITE(20,212) (GARB, I=1,8) 
      REWIND 20 
 1910 RETURN 
      END                                           
!                                                                       
!*******************************************************                
      SUBROUTINE SRFSAG 
!*******************************************************                
!                                                                       
!            THIS ROUTINE ADDED 9-11-80; IT COMPUTES A TABLE            
!            OF SAGS OF SURFACE "ISURF",FOR HEIGHTS BETWEEN             
!            "YMIN" AND "YMAX" (INCLUSIVE) WITH A STEP OF "DELY".       
!            A CONSTANT "CONST" CAN BE ADDED AS WELL.  THE USER CAN     
!            CHOOSE TO HAVE THE DIFFERENCES IN SAG BETWEEN THE          
!            SURFACE AND A SURFACE OF CURVATURE "REFCRV" (E.G. A SPHERE)
!            COMPUTED AS WELL, OR HAVE THE SAGS COMPUTED IN TERMS       
!            OF WAVELENGTHS OF LIGHT USING "WAVENM" AS THE REFERENCE    
!            WAVELENGTH; "WAVENM" IS THE WAVELENGTH IN NANOMETERS       
!            THE SURFACE MUST BE ROTATIONALLY SYMMETRIC                 
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /SAGPAR/ YMAX,YMIN,DELY,REFCRV,CONST,WAVENM,ISRF 
!                                                                       
      DATA ZERO/0.D0/,ONE/1.D0/,TWO/2.D0/,THREE/3.D0/,IZERO/0/ 
!                 TEST TO SEE IF ISRF IS A VALID SURFACE NUMBER         
!                                                                       
      IF ( ISRF.GT.NSURF .OR. ISRF.LE.IZERO ) GO TO 180 
      I = ISRF 
!                 ISRF IS A VALID SURFACE NUMBER IF HERE                
!                 WRITE HEADER LABELS                                   
!                 CONVERT WAVELENGTH TO SYSTEM UNITS FOR CALCULATIONS   
      UNITS = UFLAG 
      IF ( UNITS.EQ.ZERO .OR. UNITS.EQ.ONE ) WAVE = WAVENM*1.E-6 
      IF ( UNITS.EQ.TWO ) WAVE = WAVENM*1.E-7 
      IF ( UNITS.EQ.THREE ) WAVE = WAVENM*1.E-7/2.54 
      IF ( REFCRV.EQ.ZERO .AND.WAVENM.EQ.ZERO) WRITE(6,10) I 
      IF ( REFCRV.EQ.ZERO .AND.WAVENM.NE.ZERO) WRITE(6,20) I 
      IF ( REFCRV.NE.ZERO .AND.WAVENM.EQ.ZERO) WRITE(6,30) I 
      IF ( REFCRV.NE.ZERO .AND.WAVENM.NE.ZERO) WRITE(6,40) I 
   10 FORMAT(//T2,'SURF ',I3,T14,'Y',T28,'SAG'/) 
   20 FORMAT(//T2,'SURF ',I3,T14,'Y',T28,'SAG',                         &
     &      T46,'SAG DIFFERENCE IN WAVELENGTHS'/)                       
   30 FORMAT(//T2,'SURF ',I3,T14,'Y',T28,'SAG',                         &
     &      T46,'SAG - REFERENCE SAG'/)                                 
   40 FORMAT(//T2,'SURF ',I3,T14,'Y',T28,'SAG',                         &
     &      T46,'SAG-REFERENCE SAG',T70,                                &
     &      T70,'SAG DIFFERENCE IN WAVELENGTHS'/)                       
!                 STARTING WITH Y = YMIN, COMPUTE SAGS                  
!                                                                       
      IF ( DABS(DELY).LT.1.D-8 ) DELY = ( YMAX-YMIN )/25. 
      Y = YMIN 
   50 CONTINUE 
      YSQ = Y*Y 
!                 COMPUTE CONIC CONTRIBUTION TO SAG                     
!                 CONIC CONSTANT = -(ECCENTRICITY)**2                   
      TEMP = C(I)*YSQ 
      DENOM = 1.+DSQRT(1.-(1.+CONIC(I))*C(I)*TEMP) 
      SAG = TEMP/DENOM 
!                 COMPUTE ACONIC CONTRIBUTION TO SAG, IF ANY            
!                 THIS INCLUDES 4TH THRU 10TH POWER TERMS               
      A4  = COEF(I,1) 
      A6  = COEF(I,2) 
      A8  = COEF(I,3) 
      A10 = COEF(I,4) 
!                                                                       
      ACONIC = DABS(A4) + DABS(A6) + DABS(A8) + DABS(A10) 
      IF( ACONIC.EQ.ZERO ) GO TO 60 
      TEMP = YSQ*(A2+YSQ*(A4+YSQ*(A6+YSQ*(A8+A10*YSQ)))) 
      SAG = SAG + TEMP 
!                                                                       
   60 CONTINUE 
!                 ADD CONSTANT TO SAG                                   
      SAG = SAG + CONST 
!                 SET "DIFF" = SAG IF SAG IN WAVELENGTHS DESIRED        
      DIFF = SAG 
!                                                                       
!                 COMPUTE REFERENCE SAG IF DESIRED                      
      IF ( REFCRV.EQ.ZERO ) GO TO 70 
      TEMP = REFCRV * YSQ 
      DENOM = 1. + DSQRT(1.-REFCRV*TEMP) 
      REFSAG = TEMP/DENOM 
      DIFF = SAG - REFSAG 
   70 CONTINUE 
!                 COMPUTE DIFFERENCE IN SAG OF THE SURFACE AND REFERENCE
!                 SURFACE IN WAVELENGTHS OF LIGHT IF DESIRED            
      IF ( WAVENM.EQ.ZERO ) GO TO 80 
      DIFFWV = DIFF/WAVE 
   80 CONTINUE 
!                 WRITE Y,SAG, DIFF, AND WAVES OF SAG                   
      IF ( REFCRV.EQ.ZERO .AND.WAVENM.EQ.ZERO ) WRITE(6,90) Y,SAG 
      IF ( REFCRV.EQ.ZERO .AND.WAVENM.NE.ZERO ) WRITE(6,100)Y,SAG,DIFFWV 
      IF ( REFCRV.NE.ZERO .AND.WAVENM.EQ.ZERO ) WRITE(6,110) Y,SAG,DIFF 
      IF( REFCRV.NE.ZERO.AND.WAVENM.NE.ZERO )                           &
     &   WRITE(6,120) Y,SAG,DIFF,DIFFWV                                 
   90 FORMAT(1X,T9,F8.3,T22,1PE15.7) 
  100 FORMAT(1X,T9,F8.3,T22,1PE15.7,T46,1PE15.7) 
  110 FORMAT(1X,T9,F8.3,T22,1PE15.7,T46,1PE15.7) 
  120 FORMAT(1X,T9,F8.3,T22,1PE15.7,T46,1PE15.7,T72,1PE15.7) 
!                                                                       
      Y = Y+DELY 
      IF ( Y.GT.YMAX ) GO TO 140 
      GO TO 50 
  140 WRITE(6,150) CONST 
  150 FORMAT(/6X,'CONSTANT ADDED TO EACH SAG',T34,'=',T36,1PE15.7) 
      WRITE(6,160) REFCRV 
  160 FORMAT(6X,'REFERENCE CURVATURE',T34,'=',T36,1PE15.7) 
      WRITE(6,170) WAVENM 
  170 FORMAT(6X,'REFERENCE WAVELENGTH',T34,'=',T36,1PE15.7,' NM'//) 
      GO TO 200 
  180 WRITE(6,190) ISRF 
  190 FORMAT(1X,/5X,'INVALID SURFACE NUMBER GIVEN FOR SAG OPTION',      &
     &       ', SURFACE NUMBER = ',I3//)                                
!                                                                       
!             RESET PARAMETERS, RETURN                                  
!                                                                       
  200 ISRF   = IZERO 
      YMAX   = ZERO 
      YMIN   = ZERO 
      DELY   = ZERO 
      REFCRV = ZERO 
      CONST  = ZERO 
      WAVENM = ZERO 
      RETURN 
      END                                           
!                                                                       
!******************************************************                 
      SUBROUTINE SURFNO (ICHAR,LEN2,ISURF,JSURF) 
!******************************************************                 
!                                                                       
      IMPLICIT REAL*8 (A-H,O-Z) 
!                                                                       
      DIMENSION ICHAR(LEN2),NUMBER(10) 
!                                                                       
!     THE PURPOSE OF THIS ROUTINE IS TO DETERMINE A SURFACE NUMBER      
!     FOR STORING DATA READ IN THE MAIN PROGRAM                         
!                                                                       
!     ICHAR  = WORD 2 OF INPUT CARD (SEE MAIN PROGRAM)                  
!     ISURF = CURRENT SURFACE NUMBER, EQUALS MAXIMUM                    
!     JSURF  = SURFACE NUMBER FROM CARD OR CURRENT SURFACE NUMBER       
!                                                                       
      DATA NUMBER /'1','2','3','4','5','6','7','8','9','0'/ 
      DATA IBLANK /' '/,IZERO/0/,ITEN/10/,IONE/1/,IFORTY/40/,ICOMM/','/ 
!                                                                       
      IF ( ICHAR(1).NE.IBLANK.AND.ICHAR(1).NE.ICOMM ) GO TO 20 
!             SURFACE NUMBER FROM CARD IS BLANK, USE CURRENT SURFACE NO.
   10 JSURF = ISURF 
      GO TO 100 
!             FOUND CHARACTERS ON CARD, DETERMINE SURFACE NUMBER        
   20 INDEX1 = IZERO 
      INDEX2 = IZERO 
      IVAL1  = IZERO 
      IVAL2  = IZERO 
      DO 40 I=1,LEN2 
        DO 30 J=1,10 
          IF ( ICHAR(I).EQ.NUMBER(J) ) GO TO 50 
   30   CONTINUE 
   40 END DO 
!             NO NUMBER, USE CURRENT SURFACE NUMBER                     
      GO TO 10 
!                                                                       
!             FOUND FIRST NUMBER, DECODE AND LOOK FOR SECOND NUMBER     
   50 INDEX1 = ITEN 
      IVAL1  = J 
      IF ( IVAL1.EQ.ITEN ) IVAL1 = IZERO 
      IF ( ICHAR(I+1).EQ.IBLANK ) GO TO 70 
!             TRY TO MATCH SECOND NUMBER                                
      DO 60 J=1,10 
        IF ( ICHAR(I+1).EQ.NUMBER(J) ) GO TO 80 
   60 END DO 
!             ONLY ONE DIGIT                                            
   70 INDEX1 = IONE 
      INDEX2 = IZERO 
      IVAL2  = IZERO 
      GO TO 90 
!                                                                       
!             FOUND MATCH, DECODE                                       
   80 INDEX2 = IONE 
      IVAL2 = J 
      IF ( IVAL2.EQ.ITEN ) IVAL2 = IZERO 
!                                                                       
!             DETERMINE SURFACE NUMBER                                  
   90 JSURF = INDEX1*IVAL1 + INDEX2*IVAL2 
!                                                                       
  100 IF ( JSURF.LT.IZERO .OR. JSURF.GT.IFORTY ) JSURF = ISURF 
      RETURN 
      END                                           
!                                                                       
!*******************************************************                
      SUBROUTINE SURTYP(I) 
!*******************************************************                
!       THIS ROUTINE IS CALLED FROM PARAX AND ADDS THE CAPABILITY       
!       OF HANDLING SURFACES ENTIRELY DESCRIBED BY A POLYNOMIAL         
!       EXPRESSION AND CONIC SURFACES WHICH HAVE DEFORMATION CONSTANTS  
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /ABCODE/ ISURFX,ISURXN,RSAVE,TM3,DUM(40),RR(40),RCON(40) 
!                                                                       
      DIMENSION TOR(40) 
!                                                                       
      EQUIVALENCE (Y0(1),TOR(1)) 
!                                                                       
      DATA ZERO/0.D0/ 
!                                                                       
      DUM(I) = ZERO 
      RR(I) = R(I) 
      RCON(I) =  CONIC(I) 
      IF ( R(I).NE.ZERO ) GO TO 10 
!        SURFTYPE MAY BE A PLANE                                        
      IF ( COEF(I,1).NE.ZERO ) GO TO 40 
!        SURFACE IS A PLANE NO ABERRATION                               
      GO TO 50 
   10 IF ( RCON(I).NE.ZERO ) GO TO 20 
!        SURFACE IS A SPHERE                                            
      GO TO 30 
!        SURFACE IS A CONIC                                             
   20 DUM(I) = RCON(I)*TM3*TM3*TM3/8. 
   30 IF ( COEF(I,1) .EQ. ZERO ) GO TO 50 
   40 DUM(I) = DUM(I)+COEF(I,1) 
   50 IF ( TOR(I).EQ.ZERO ) GO TO 60 
      ENTRY SURTOR(I) 
      ISURXN = I 
      ISURFX = 1 
      IF ( R(I).NE.ZERO .AND. RX(I).NE.ZERO ) ISURFX = 2 
      IF ( R(I).NE.ZERO .AND. RX(I).EQ.ZERO ) ISURFX = 3 
      IF ( R(I).EQ.ZERO .AND. RX(I).NE.ZERO ) ISURFX = 4 
      RSAVE = RX(I) 
      IF ( ISURFX.EQ.3 ) RSAVE = R(I) 
   60 RETURN 
      END                                           
!                                                                       
!*******************************************************                
      SUBROUTINE RED 
!*******************************************************                
!                ROUTINE RED PERFORMS RADIAL ENERGY DISTRIBUTION        
!                CALCULATIONS                                           
      IMPLICIT REAL *8 (A-H,O-Z) 
      INTEGER *4 OPTNA,OPTNB,OPTNC,OPTND,OPTNE,OPTNF,OPTNG,ANULI,SECTRS 
!                                                                       
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /HEAD/   LINES,IPAGE,NSYS,LAMDA,NAME(160) 
      COMMON /PUPIL/  ENPUPR,ENPUPL,EXPUPR,EXPUPL 
      COMMON /SAGPAR/ YMAX,YMIN,DELY,REFCRV,CONST,WAVENM,ISRF 
      COMMON /CRED/ RPP(40),EPP(40),XW(812),HYSAVE,HXSAVE,              &
     &             TIMES,RD,OPTNE,ISW,IOPB                              
      COMMON /CMTF/   AJ, AN, AR(40), DE(40), DELTCC, DELTEN, DELTPL,   &
     &                EN, IOPC, OPTNF, PRINAN, RSC, SUM, TX(51)         
      COMMON /CSPOT/  XTMAX,XTMIN,YTMAX,YTMIN,AVGX,AVGY,RP(812),SPOTP,  &
     &                XK(812),YK(812),XKALL(7600),YKALL(7600),JSKIP,    &
     &                IOPA,NTHRU                                        
!                                                                       
!                                                                       
      DIMENSION YW(812),TOR(40) 
      DIMENSION ZF(10),BLOB(10) 
      DIMENSION XSUM(10),YSUM(10),XSUMSQ(10),YSUMSQ(10) 
      DIMENSION QT(3),XX(3),XYZ(3) 
!                                                                       
      EQUIVALENCE (QT(1),QX),(QT(2),QY),(QT(3),QZ) 
      EQUIVALENCE (XX(1),XT),(XX(2),YT),(XX(3),ZT) 
      EQUIVALENCE (XYZ(1),X),(XYZ(2),Y),(XYZ(3),Z) 
      EQUIVALENCE (TRASH(1),XSMIN),(TRASH(2),XSMAX),(TRASH(3),YSMIN) 
      EQUIVALENCE (TRASH(4),YSMAX),(TRASH(5),FCODE),(TRASH(6),FXYJ) 
      EQUIVALENCE (TRASH(7),FXY),(TRASH(8),FXNU),(TRASH(9),FBY) 
      EQUIVALENCE (TRASH(10),FBNU),(Y0(1),TOR(1)) 
!                                                                       
!                                                                       
      DATA IZERO/0/,IONE/1/ 
      DATA ZERO/0.D0/,ONE/1.D0/,TWO/2.D0/,THREE/3.D0/,FOUR/4.D0/ 
      DATA FIVE/5.D0/,SIX/6.D0/,TEN/10.D0/,HUNDRD/100.D0/,THOUS/1000.D0/ 
      DATA EPS1/1.D-9/,EPS2/1.D-14/,PI/3.141592653589793D0/ 
  140 FORMAT (11H0Y OBJ.PNT.,5X,1PE15.8,2X,10HX OBJ.PNT.,5X,1PE15.8) 
  150 FORMAT (9H0COLOR   ,5X,I6) 
  160 FORMAT (1X,'  PERCENT ENERGY',7X,'RADIUS') 
  170 FORMAT (1X,I8,8X,1PE15.8) 
!                CALCULATE PERCENT ENERGY AND ADJUSTED RADIUS ARRAYS    
 1710   DO 1720 I=1,40 
          IPL=I*.025*TIMES 
          EPP(I)=.025*I 
          RPP(I)=RP(IPL) 
 1720   CONTINUE 
        IF (ISW.EQ.0) RSMAX=0. 
!            PRINT RED TABLE?                                           
 1730   IF (OPTNE.EQ.0) GO TO 1750 
        WRITE(6,140) HYSAVE,HXSAVE 
        WRITE(6,150) LAMDA 
        WRITE(6,160) 
        DO 1740 I=2,40,2 
          II=I*2.5 
          WRITE(6,170) II,RPP(I) 
 1740   CONTINUE 
!            PLOT RED?                                                  
       IF (IOPB.EQ.0) GO TO 1750 
       CALL GRAPH(2,2) 
 1750  RETURN 
      END                                           
!                                                                       
!*******************************************************                
      SUBROUTINE MTF 
!                                                                       
!        PURPOSE                                                        
!               PERFORM MODULATION TRANSFER CALCULATIONS                
!                                                                       
!        PARAMETERS -NONE                                               
!                                                                       
!        REFERENCE - GENOPTICS USER GUIDE CHAP. 4                       
!                                                                       
!        ACTIVE VARIABLES                                               
!                                                                       
!           INPUT:                                                      
!   /PMATX/    FOCL           SYSTEM FOCAL LENGTH                       
!                               DEFAULT IF(FOCL = 0.0) FOCL = 100.0     
!                                                                       
!   /CMTF/     IOPC           MTF PLOT SWITCH                           
!   /COLLAT/   IWVFLG(I)      WAVELENGTH SET FLAG                       
!   /HEAD/     LAMDA          RAY COLOR INDEX                           
!   /CMTF/     OPTNF          MTF PRINT OPTION SWITCH(INTEGER*4)        
!   /CMTF/     PRINAN         SUM OF (QY/QZ) FOR EACH RAY               
!   /CSPOT/    RP(I)          SEE USER'S GUIDE                          
!   /PMATX/    SMAX           MAX OBJECT DISTANCE (?)                   
!                               SEE LINE LABLED 1765 BELOW              
!                             *** NOTE  AS ISW IS USED HERE, IF         
!                             SMAX = 0.0 ON ENTRY, THEN SMAX = 0.0      
!                             ON EXIT.  THE DEFAULT IS USED ONLY        
!                             INTERNALLY TO  MTF.                       
!                                                                       
!   /CRED/     TIMES          NUMBER OF RAYS(?)                         
!   /PMATX/    WAVL(LAMDA)    USER SPECIFIED WAVELENGTH/DEFAULTS        
!           OUTPUT:                                                     
!  /CMTF/          AN         ANGLE OF PRINCIPLE RAY(?)                 
!  /CMTF/          DE(I)      ?                                         
!                  DEFWV(3)   DEFAULT WAVELENGTHS (NANOMETERS)          
!                               DEFWV(1) = 632.8                        
!                               DEFWV(2) = 587.6                        
!                               DEFWV(3) = 486.1                        
!                                                                       
!  /CMTF/          DELTCC     ?                                         
!  /CMTF/          DELTEN     ?                                         
!  /CMTF/          DELTPL     ?                                         
!  /CMTF/          EN                                                   
!  /CMTF/          EPP                                                  
!                  IPL                                                  
!  /CRED/          ISW                                                  
!  /CRED/          RPP                                                  
!  /CMTF/          RSC                                                  
!  /CMTF/          SUM                                                  
!  /CMTF/          TX                                                   
!  /PMATX/         WAVL                                                 
!                  X                                                    
!                                                                       
!  LAST UPTDATE 4/20/84 BY JOHN C PARKER OF SSAI                        
!                                                                       
!*******************************************************                
!     SUBROUTINE MTF                                                    
!                ROUTINE MTF PERFORMS MODULATION TRANSFER               
!                FUNCTION CALCULATIONS                                  
      IMPLICIT REAL *8 (A-H,O-Z) 
      INTEGER *4 OPTNA,OPTNB,OPTNC,OPTND,OPTNE,OPTNF,OPTNG,ANULI,SECTRS 
!                                                                       
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /HEAD/   LINES,IPAGE,NSYS,LAMDA,NAME(160) 
      COMMON /PUPIL/  ENPUPR,ENPUPL,EXPUPR,EXPUPL 
      COMMON /SAGPAR/ YMAX,YMIN,DELY,REFCRV,CONST,WAVENM,ISRF 
      COMMON /CRED/   RPP(40),EPP(40),XW(812),HYSAVE,HXSAVE,            &
     &                TIMES,RD,OPTNE,ISW,IOPB                           
      COMMON /CMTF/   AJ, AN, AR(40), DE(40), DELTCC, DELTEN, DELTPL,   &
     &                EN, IOPC, OPTNF, PRINAN, RSC, SUM, TX(51)         
      COMMON /CSPOT/  XTMAX,XTMIN,YTMAX,YTMIN,AVGX,AVGY,RP(812),SPOTP,  &
     &                XK(812),YK(812),XKALL(7600),YKALL(7600),JSKIP,    &
     &                IOPA,NTHRU                                        
!                                                                       
!                                                                       
!                                                                       
      DIMENSION YW(812),TOR(40),DEFWV(3) 
      DIMENSION ZF(10),BLOB(10) 
      DIMENSION XSUM(10),YSUM(10),XSUMSQ(10),YSUMSQ(10) 
      DIMENSION QT(3),XX(3),XYZ(3) 
!                                                                       
      EQUIVALENCE (QT(1),QX),(QT(2),QY),(QT(3),QZ) 
      EQUIVALENCE (XX(1),XT),(XX(2),YT),(XX(3),ZT) 
      EQUIVALENCE (XYZ(1),X),(XYZ(2),Y),(XYZ(3),Z) 
      EQUIVALENCE (TRASH(1),XSMIN),(TRASH(2),XSMAX),(TRASH(3),YSMIN) 
      EQUIVALENCE (TRASH(4),YSMAX),(TRASH(5),FCODE),(TRASH(6),FXYJ) 
      EQUIVALENCE (TRASH(7),FXY),(TRASH(8),FXNU),(TRASH(9),FBY) 
      EQUIVALENCE (TRASH(10),FBNU),(Y0(1),TOR(1)) 
!                                                                       
!                                                                       
      DATA IZERO/0/,IONE/1/ 
      DATA ZERO/0.D0/,ONE/1.D0/,TWO/2.D0/,THREE/3.D0/,FOUR/4.D0/ 
      DATA FIVE/5.D0/,SIX/6.D0/,TEN/10.D0/,HUNDRD/100.D0/,THOUS/1000.D0/ 
      DATA EPS1/1.D-9/,EPS2/1.D-14/,PI/3.141592653589793D0/ 
!                                                                       
  190 FORMAT (///9X,5HFREQ.,11X,3HMTF/,6X,11H(CYCLES/MR)/) 
  200 FORMAT (1X,1PE15.8,2X,1PE15.8) 
  208 FORMAT(//10X,'NO WAVELENGTH INPUT FOR MTF; ',F5.1,                &
     &       ' NM BEING USED'/)                                         
  209 FORMAT(10X,'NO FOCAL LENGTH INPUT FOR MTF; F.L. SET = 100'/) 
!                                                                       
!                                                                       
 1710   DO 1720 I=1,40 
          IPL=I * IDINT(.025*TIMES) 
          EPP(I)=.025 * DFLOAT(I) 
          RPP(I)=RP(IPL) 
 1720   CONTINUE 
        DEFWV(1)=632.8 
        DEFWV(2)=587.6 
        DEFWV(3)=486.1 
        ISW=0 
        IF (SMAX.NE.ZERO) GO TO 1770 
!                IF NO MAX CYCLES/MR SPECIFIED, CALCULATE               
        ISW=1 
        IF (IWVFLG(LAMDA).EQ.ZERO) GO TO 1760 
        WRITE(6,208) DEFWV(LAMDA) 
        IF (UFLAG.EQ.ONE) WAVL(LAMDA) = DEFWV(LAMDA) * 1.D-6 
        IF (UFLAG.EQ.TWO) WAVL(LAMDA) = DEFWV(LAMDA) * 1.D-7 
        IF (UFLAG.EQ.THREE) WAVL(LAMDA) = DEFWV(LAMDA) * 1.D-7/2.54D0 
 1760   IF ( FOCL.NE.ZERO ) GO TO 1765 
        WRITE(6,209) 
        FOCL = HUNDRD 
!             COMPUTE MAX SPATIAL FREQUENCY, CYCLES/RAD                 
 1765   SMAX=ONE/(WAVL(LAMDA)*TWO*RHO/DABS(FOCL)) 
!             CONVERT TO CYCLES/MR                                      
        SMAX=SMAX/THOUS 
        JSMAX = IDINT(SMAX) 
        SMAX = DFLOAT(JSMAX) 
!                AN IS ANGLE OF PRINCIPAL RAY                           
 1770 CONTINUE 
        AN=DATAN(PRINAN/TIMES) 
!                SET UP MTF CALCULATION PARAMETERS                      
        RSC=ONE/DCOS(AN) 
        DELTEN=SMAX/50. 
        DELTPL=DELTEN*FIVE 
        DELTCC=ONE/DELTPL 
        IF ( OPTNF.EQ.0 ) GO TO 1800 
!             HEAD MTF PRINT                                            
 1780   WRITE(6,190) 
!                CALCULATE MTF                                          
!                MODULATION TRANSFER FUNCTION (MTF) CAN BE CALCULATED   
!                FROM RADIAL ENERGY DISTRIBUTION DATA BY THE FOLLOWING  
!                FORMULA FOR EACH FREQUENCY                             
!                MTF = SUM FROM 1 TO M (HERE 40) OF                     
!                     DELTA E SUB I * J SUB ZERO OF                     
!                     (2.*PI*NU*R BAR SUB I)                            
!                WHERE                                                  
!                DELTA E SUB I IS THE CHANGE IN ENERGY E(I+1)-E(I)      
!                OVER EACH INTERVAL                                     
!                J SUB ZERO IS THE ZERO ORDER BESSEL FUNCTION           
!                R BAR SUB I IS THE AVERAGE RADIUS OVER THE             
!                INTERVAL (R(I+1)+R(I))/2.                              
!                NU IS THE FREQUENCY FOR WHICH THE FUNCTION IS          
!                CALCULATED.NU=EN*1000/!FOCL!                           
!                NU IS IN CYCLES PER UNIT LENGTH                        
 1800   DO 1810 I=1,39 
          AR(I)=(RPP(I+1)+RPP(I))/2 
          DE(I)=EPP(I+1)-EPP(I) 
 1810   CONTINUE 
        DO 1850 K=1,51 
          EN=(K-1)*DELTEN 
          SUM=ZERO 
          DO 1830 I=1,39 
!                                                                       
! NEXT LINE CORRECTED 4/11/84 BY JCP                                    
!                                                                       
            X=TWO*PI*EN*AR(I)*THOUS/(DABS(FOCL) * RSC) 
            CALL BESJN (X,0,AJ,ONE,IER) 
            IF ( IER.EQ.0 ) GO TO 1825 
            WRITE(6,1820) X,IER 
 1820       FORMAT (25X,'X =',3X,1PE15.8,5X,'IER =',I5) 
            AJ=ONE 
 1825       SUM=SUM+DE(I)*AJ 
 1830     CONTINUE 
          TX(K)=SUM/EPP(39) 
          IF (OPTNF.EQ.0) GO TO 1850 
!             PRINT MTF VALUES                                          
 1840     WRITE(6,200) EN,TX(K) 
 1850   CONTINUE 
        IF (IOPC.EQ.0) GO TO 1860 
        CALL GRAPH(2,3) 
 1860   IF (ISW.EQ.1) SMAX=ZERO 
 1900 RETURN 
      END                                           
!                                                                       
!********************************************************               
      SUBROUTINE GRAPH(IGRSW,IKIND) 
!********************************************************               
!                ROUTINE GRAPH PLOTS SPOT,RED, AND MTF DATA             
!                                                                       
      IMPLICIT REAL *8 (A-H,O-Z) 
      INTEGER *4 OPTNA,OPTNB,OPTNC,OPTND,OPTNE,OPTNF,OPTNG,ANULI,SECTRS 
!                                                                       
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /HEAD/   LINES,IPAGE,NSYS,LAMDA,NAME(160) 
      COMMON /PUPIL/  ENPUPR,ENPUPL,EXPUPR,EXPUPL 
      COMMON /SAGPAR/ YMAX,YMIN,DELY,REFCRV,CONST,WAVENM,ISRF 
      COMMON /CRED/   RPP(40),EPP(40),XW(812),HYSAVE,HXSAVE,            &
     &                TIMES,RD,OPTNE,ISW,IOPB                           
      COMMON /CMTF/   AJ, AN, AR(40), DE(40), DELTCC, DELTEN, DELTPL,   &
     &                EN, IOPC, OPTNF, PRINAN, RSC, SUM, TX(51)         
      COMMON /CSPOT/  XTMAX,XTMIN,YTMAX,YTMIN,AVGX,AVGY,RP(812),SPOTP,  &
     &                XK(812),YK(812),XKALL(7600),YKALL(7600),JSKIP,    &
     &                IOPA,NTHRU                                        
!                                                                       
!                                                                       
      DIMENSION YW(812),TOR(40) 
      DIMENSION ZF(10),BLOB(10) 
      DIMENSION XSUM(10),YSUM(10),XSUMSQ(10),YSUMSQ(10) 
      DIMENSION QT(3),XX(3),XYZ(3) 
      DIMENSION MSGA(4),MSGB(4),MSGC(6),MSGD(6),MSGE(3),MSGF(2),MSGG(4) 
      DIMENSION MSGH(3),MSGI(3),MSGJ(4),MSGK(6),MSGL(4),MSGM(11),MSGN(2) 
      DIMENSION MSGO(5),MSGP(5),MSGQ(7),MSGR(1),MSGS(6),MSGT(6),MSGU(6) 
      DIMENSION MSGV(4),MSGW(4),MSGX(4),MSGY(5),MSGZ(5) 
      DIMENSION MSG1(4),MSG2(1),MSG3(5),MSG4(3) 
!                                                                       
      EQUIVALENCE (QT(1),QX),(QT(2),QY),(QT(3),QZ) 
      EQUIVALENCE (XX(1),XT),(XX(2),YT),(XX(3),ZT) 
      EQUIVALENCE (XYZ(1),X),(XYZ(2),Y),(XYZ(3),Z) 
      EQUIVALENCE (TRASH(1),XSMIN),(TRASH(2),XSMAX),(TRASH(3),YSMIN) 
      EQUIVALENCE (TRASH(4),YSMAX),(TRASH(5),FCODE),(TRASH(6),FXYJ) 
      EQUIVALENCE (TRASH(7),FXY),(TRASH(8),FXNU),(TRASH(9),FBY) 
      EQUIVALENCE (TRASH(10),FBNU),(Y0(1),TOR(1)) 
!                                                                       
      DATA MSGA/4HX-AX,4HIS I,4HN IN,4HCHES/ 
      DATA MSGB/4HY-AX,4HIS I,4HN IN,4HCHES/ 
      DATA MSGC/4HX-AX,4HIS I,4HN CE,4HNTIM,4HETER,4HS   / 
      DATA MSGD/4HY-AX,4HIS I,4HN CE,4HNTIM,4HETER,4HS   / 
      DATA MSGE/4HSCAL,4HE FA,4HCTOR/ 
      DATA MSGF/4HCOLO,4HR   / 
      DATA MSGG/4HNUMB,4HER O,4HF PO,4HINTS/ 
      DATA MSGH/4HX-AV,4HERAG,4HE = / 
      DATA MSGI/4HY-AV,4HERAG,4HE = / 
      DATA MSGJ/4HRADI,4HUS I,4HN IN,4HCHES/ 
      DATA MSGK/4HRADI,4HUS I,4HN CE,4HNTIM,4HETER,4HS   / 
      DATA MSGL/4HPERC,4HENT ,4HENER,4HGY  / 
      DATA MSGM/4HSPAT,4HIAL ,4HFREQ,4HUENC,4HY, C,4HYCLE,4HS PE,4HR MI,&
     &          4HLLIR,4HADIA,4HN   /                                   
      DATA MSGN/4HCONT,4HRAST/ 
      DATA MSGO/4HGEOM,4HETRI,4HCAL ,4HMTF ,4H    / 
      DATA MSGP/4HGEOM,4HETRI,4HCAL ,4HRED ,4H    / 
      DATA MSGQ/4HOBJE,4HCT H,4HEIGH,4HT (R,4HADIA,4HNS) ,4H X= / 
      DATA MSGR/4H Y= / 
      DATA MSGS/4HX-AX,4HIS I,4HN MI,4HLLIM,4HETER,4HS   / 
      DATA MSGT/4HY-AX,4HIS I,4HN MI,4HLLIM,4HETER,4HS   / 
      DATA MSGU/4HRADI,4HUS I,4HN MI,4HLLIM,4HETER,4HS   / 
      DATA MSGV/4HX-AX,4HIS I,4HN UN,4HITS / 
      DATA MSGW/4HY-AX,4HIS I,4HN UN,4HITS / 
      DATA MSGX/4HRADI,4HUS I,4HN UN,4HITS / 
      DATA MSGY/4HX-AX,4HIS I,4HN RA,4HDIAN,4HS   / 
      DATA MSGZ/4HY-AX,4HIS I,4HN RA,4HDIAN,4HS   / 
      DATA MSG1/4HTAN(,4HANGL,4HE)  ,4H    / 
      DATA MSG2/4H    / 
      DATA MSG3/4HANGL,4HE OF,4H INC,4HIDEN,4HCE  / 
      DATA MSG4/4HSPOT,4H RAD,4HIUS=/ 
!                                                                       
      DATA IZERO/0/,IONE/1/, IFOUR/4/ 
      DATA ZERO/0.D0/,ONE/1.D0/,TWO/2.D0/,THREE/3.D0/,FOUR/4.D0/ 
      DATA FIVE/5.D0/,SIX/6.D0/,TEN/10.D0/,HUNDRD/100.D0/,THOUS/1000.D0/ 
      DATA EPS1/1.D-9/,EPS2/1.D-14/,PI/3.141592653589793D0/ 
  207 FORMAT(10X,'POINT SKIPPED ON SPOT PLOT') 
!                                                                       
      IF (IGRSW.EQ.1) GO TO 2000 
      IF (IGRSW.EQ.2) GO TO 10 
      IF (IGRSW.EQ.3) GO TO 3000 
      WRITE(6,1) 
    1 FORMAT(T5,' ERROR IN IGRSW ') 
      GO TO 4000 
   10 COL = LAMDA 
      HTNX=DATAN2(HXSAVE,S-D) 
      HTNY=DATAN2(HYSAVE,S-D) 
      CALL PLOT(.5,0.,3) 
      DIV=TEN 
      DELTPL=DELTEN*FIVE 
      IF (IKIND.EQ.1) GO TO 1400 
      IF (IKIND.EQ.2) GO TO 1600 
      IF (IKIND.EQ.3) GO TO 1800 
      IF (IKIND.EQ.4) GO TO 1400 
      WRITE(6,20) 
   20 FORMAT(T5,' ERROR IN IKIND') 
      GO TO 4000 
!            HERE IF SPOT PLOT                                          
 1400   IF (IOPA.EQ.0) GO TO 1420 
        ISW=0 
!                XD, YD = AXIS INCREMENTS                               
        XD=(XSMAX-XSMIN)/SIX 
        YD=(YSMAX-YSMIN)/SIX 
!                IF NO MAX, MIN SPECIFIED FOR AXES, SCALE AXES          
        IF (XD*YD.EQ.ZERO) GO TO 1410 
!                XW(4)= ADJUSTED X INCREMENT                            
!                ISW IS ZERO RESTORE SWITCH                             
        ISW=1 
        GO TO 1420 
!                SET UP SCALE SUBR PARAMETERS                           
 1410   XW(1)=XMINS 
        XW(2)=XMAXS 
        XW(5)=YMINS 
        XW(6)=YMAXS 
        IF (XMINS.EQ.XW(2)) XW(2)=XMINS*TWO 
        IF (YMINS.EQ.XW(6)) XW(6)=YMINS*TWO 
        CALL SCALE(XW(1),6.,2,2) 
        CALL SCALE(XW(5),6.,2,2) 
!                ON RETURN, XW(3)= ADJUSTED X MIN                       
!                XW(7)= ADJUSTED Y MIN                                  
!                XW(8)= ADJUSTED Y INC                                  
        XSMIN=XW(3) 
        XSMAX=XW(4)*SIX+XSMIN 
        YSMIN=XW(7) 
        YSMAX=XW(8)*SIX+YSMIN 
        XD=XW(4) 
        YD=XW(8) 
 1420   CONTINUE 
        XT=.5 
        YT=ZERO 
        DO 1430 I=1,42 
          YT=YT+.25 
          XT=.75-XT 
 1430   CALL PLOT (XT,YT,2) 
        XT=.75-XT 
        CALL PLOT (XT,YT,2) 
        DO 1440 I=1,42 
          XT=.75-XT 
          YT=YT-.25 
 1440   CALL PLOT (XT,YT,2) 
        SF = ONE/XD 
        IF (IALLPL.EQ.1) GO TO 1445 
        CALL SYMBOL (1.2,10.0,.14,MSGQ,0.,28) 
        CALL NUMBER (5.1,10.0,.14,HTNX,0.,6) 
        CALL SYMBOL (6.3,10.0,.14,MSGR,0.,3) 
        CALL NUMBER (6.8,10.0,.14,HTNY,0.,6) 
 1445   CALL SYMBOL (1.2,9.7,.14,MSGF,0.,6) 
        CALL NUMBER (2.0,9.7,.14,COL,0.,-1) 
        CALL SYMBOL (2.3,9.7,.14,MSGG,0.,16) 
        CALL NUMBER (4.6,9.7,.14,TIMES,0.,-1) 
        CALL SYMBOL (5.2,9.7,.14,MSGE,0.,12) 
        CALL NUMBER (7.0,9.7,.14,SF,0.,4) 
        CALL SYMBOL (1.2,9.4,.14,MSGH,0.,11) 
        CALL NUMBER (2.8,9.4,.14,AVGX,0.,6) 
        CALL SYMBOL (4.3,9.4,.14,MSGI,0.,11) 
        CALL NUMBER (5.9,9.4,.14,AVGY,0.,6) 
        CALL SYMBOL (1.2,9.1,.14,MSG4,0.,12) 
        CALL NUMBER (3.0,9.1,.14,SPOTP,0.,6) 
        CALL PLOT (1.2,2.5,-3) 
!             UFLAG=ONE IF UNITS ARE MILLIMETERS                        
!             UFLAG=TWO IF UNITS ARE CENTIMETERS                        
!             UFLAG=THREE IF UNITS ARE INCHES                           
!             UFLAG=FOUR IF UNITS ARE RADIANS                           
!             UFLAG=FIVE, IF PLOTTING ANGLES OF INCIDENCE               
        IF (UFLAG.EQ.ONE)   GO TO 1460 
        IF (UFLAG.EQ.TWO)   GO TO 1465 
        IF (UFLAG.EQ.THREE) GO TO 1470 
        IF (UFLAG.EQ.FOUR)  GO TO 1475 
        IF (UFLAG.EQ.FIVE)  GO TO 1480 
        GO TO 1485 
 1460   CALL AXIS (0.,0.,MSGT,21,6.,90.,YSMIN,YD) 
        GO TO 1490 
 1465   CALL AXIS (0.,0.,MSGD,21,6.,90.,YSMIN,YD) 
        GO TO 1490 
 1470   CALL AXIS (0.,0.,MSGB,16,6.,90.,YSMIN,YD) 
        GO TO 1490 
 1475   CALL AXIS (0.,0.,MSGZ,17,6.,90.,YSMIN,YD) 
        GO TO 1490 
 1480   CALL AXIS (0.,0.,MSG3,-18,6.,0.,XSMIN,XD) 
        GO TO 1530 
 1485   CALL AXIS (0.,0.,MSGW,15,6.,90.,YSMIN,YD) 
 1490   CALL AXIS (0.,6.,MSG2,1,6.,0.,XSMIN,XD) 
        CALL AXIS (6.,0.,MSG2,-1,6.,90.,YSMIN,YD) 
        IF (UFLAG.EQ.ONE)   GO TO 1495 
        IF (UFLAG.EQ.TWO)   GO TO 1500 
        IF (UFLAG.EQ.THREE) GO TO 1505 
        IF (UFLAG.EQ.FOUR)  GO TO 1510 
        IF (UFLAG.EQ.FIVE)  GO TO 1530 
        GO TO 1515 
 1495   CALL AXIS (0.,0.,MSGS,-21,6.,0.,XSMIN,XD) 
        GO TO 1520 
 1500   CALL AXIS (0.,0.,MSGC,-21,6.,0.,XSMIN,XD) 
        GO TO 1520 
 1505   CALL AXIS (0.,0.,MSGA,-16,6.,0.,XSMIN,XD) 
        GO TO 1520 
 1510   CALL AXIS (0.,0.,MSGY,-17,6.,0.,XSMIN,XD) 
        GO TO 1520 
 1515   CALL AXIS (0.,0.,MSGV,-15,6.,0.,XSMIN,XD) 
 1520   XP = .5 
        DO 1525 IND=1,25 
          INDD = IND+80 
          CALL SYMBOL (XP,-1.,.2,NAME(INDD),0.,1) 
          XP = XP+.2 
 1525   CONTINUE 
 1530   JSKIP=0 
 1540   NPTS=NTHRU 
        IF (IALLPL.EQ.1) NPTS=IJK 
        DO 1560 I=1,NPTS 
!                XT, YT ARE PLOT COORD                                  
          IF (IALLPL.EQ.1) GO TO 1544 
          XT=(XK(I)-XSMIN)/XD 
          YT=(YK(I)-YSMIN)/YD 
          GO TO 1546 
 1544     XT=(XKALL(I)-XSMIN)/XD 
          YT=(YKALL(I)-YSMIN)/YD 
 1546     IF (I.NE.1) GO TO 1550 
!                INITIALIZE MAX, MIN COORD                              
          XTMAX=XT 
          YTMAX=YT 
          XTMIN=XT 
          YTMIN=YT 
!                ACCUM MAX, MIN COORD CALCULATED                        
 1550     IF (XT.GT.XTMAX) XTMAX=XT 
          IF (YT.GT.YTMAX) YTMAX=YT 
          IF (YT.LT.YTMIN) YTMIN=YT 
          IF (XT.LT.XTMIN) XTMIN=XT 
!                IF POINT OUT OF GRID, SKIP PLOTTING IT                 
          IF (XT.LT.ZERO .OR. XT.GT.SIX .OR. YT.LT.ZERO .OR. YT.GT.SIX) &
     &        JSKIP = I                                                 
          IF ( JSKIP.EQ.I ) WRITE(6,207) 
          IF (JSKIP.EQ.I) GO TO 1560 
!                PLOT THE POINT                                         
          CALL SYMBOL (XT,YT,.05,1HX,0.,1) 
 1560   CONTINUE 
!                CLOSE UNIT PLOT, RESTORE ZEROS                         
        CALL PLOT (7.3,-2.5,-3) 
 1570   CONTINUE 
        GO TO 4000 
 1600 CALL PLOT(1.0,.4,-3) 
        ISW=0 
        RD=RSMAX/TEN 
!                WAS MAX RADIUS SPECIFIED                               
        IF (RD.EQ.ZERO) GO TO 1630 
        ISW=1 
        GO TO 1640 
!                IF NO MAX SPECIFIED, SCALE GRAPH                       
 1630   EPP(1)=ZERO 
        EPP(2)=XW(NTHRU) 
        RSMAX=EPP(4)*TEN 
!                EPP(4) IS ADJUSTED INCREMENT                           
        RD=EPP(4) 
 1640 CALL SCALE(EPP(1),10.,2,2) 
!          HEAD RED PLOT                                                
      IF (UFLAG.EQ.ONE) GO TO 1650 
      IF (UFLAG.EQ.TWO) GO TO 1660 
      IF (UFLAG.EQ.THREE) GO TO 1670 
      IF (UFLAG.EQ.FOUR) GO TO 1690 
      GO TO 1680 
 1650 CALL AXIS(0.,10.,MSGU,-21,10.,270.,RD,DIV) 
      GO TO 1700 
 1660 CALL AXIS(0.,10.,MSGK,-21,10.,270.,RD,DIV) 
      GO TO 1700 
 1670 CALL AXIS(0.,10.,MSGJ,-16,10.,270.,RD,DIV) 
      GO TO 1700 
 1680 CALL AXIS(0.,10.,MSGX,-15,10.,270.,RD,DIV) 
      GO TO 1700 
 1690 CALL AXIS(0.,10.,MSG1,-16,10.,270.,RD,DIV) 
 1700 CALL AXIS(0.,10.,MSGL,14,5.,0.,20.,DIV) 
      CALL PLOT(0.,10.,3) 
 1710 DO 1720 I=1,40 
       IPL=I*.025*TIMES 
       XT=I*.125 
       YT=TEN-(RP(IPL))/RD 
       IF (YT.LT.ZERO) GO TO 1720 
!          PLOT 40 POINTS                                               
       CALL PLOT(XT,YT,2) 
 1720 END DO 
      CALL SYMBOL(6.65,9.2,.2,MSGP,-90.,17) 
      YP=5.8 
      DO 1725 I=1,25 
         INDD=I+80 
         CALL SYMBOL(6.7,YP,.2,NAME(INDD),-90.,1) 
         YP=YP-.2 
 1725 END DO 
      CALL SYMBOL(6.25,8.4,.14,MSGQ,-90.,28) 
      CALL NUMBER(6.25,4.5,.14,HTNX,-90.,6) 
      CALL SYMBOL(6.25,3.3,.14,MSGR,-90.,3) 
      CALL NUMBER(6.25,2.8,.14,HTNY,-90.,6) 
!          CLOSE RED PLOT                                               
      CALL PLOT(7.5,-.4,-3) 
      IF (ISW.EQ.0) RSMAX=0. 
      GO TO 4000 
 1800 CALL PLOT(1.0,.4,-3) 
      CALL AXIS(0.,10.,MSGM,-41,10.,270.,0.,DELTPL,DIV) 
      CALL AXIS(0.,10.,MSGN,8,5.,0.,0.,.2,DIV) 
      CALL PLOT(0.,10.,3) 
      DO 1810 K=1,51 
      EN=(K-1)*DELTEN 
      XT=TX(K)*FIVE 
      YT=TEN-EN*DELTCC 
!            PLOT MTF VALUES                                            
      CALL PLOT(XT,YT,2) 
 1810 END DO 
      CALL SYMBOL(6.65,9.2,.2,MSGO,-90.,17) 
      YP=5.8 
      DO 1820 I=1,25 
        INDD=I+80 
        CALL SYMBOL(6.7,YP,.2,NAME(INDD),-90.,1) 
        YP=YP-.2 
 1820 END DO 
      CALL SYMBOL(6.25,8.4,.14,MSGQ,-90.,28) 
      CALL NUMBER(6.25,4.5,.14,HTNX,-90.,6) 
      CALL SYMBOL(6.25,3.3,.14,MSGR,-90.,3) 
      CALL NUMBER(6.25,2.8,.14,HTNY,-90.,6) 
!          CLOSE PLOT                                                   
      CALL PLOT(7.5,-.4,-3) 
      IF (ISW.EQ.1) SMAX=0 
      GO TO 4000 
 2000 CALL PLOTS(53,0,10) 
!           SET BLOCKSIZE = 1024 FOR ZETA ERROR CORRECTION MODE         
!           AT 1200 BAUD, M256 = 2 FOR 300 BAUD ERROR CORRECTION MODE   
      M256=IFOUR 
      CALL ZETOBS(M256) 
      GO TO 4000 
 3000 CALL PLOT(0.,0.,999) 
 4000 CONTINUE 
      RETURN 
      END                                           
!                                                                       
!******************************************************                 
      SUBROUTINE FNDSPT 
!******************************************************                 
!                ROUTINE SPOT PERFORMS SPOT PLOT AND PRINT              
!                CALCULATIONS                                           
      IMPLICIT REAL *8 (A-H,O-Z) 
      INTEGER *4 OPTNA,OPTNB,OPTNC,OPTND,OPTNE,OPTNF,OPTNG,ANULI,SECTRS 
!                                                                       
      COMMON /PMATX/  TRASH(10),S,D,RHO,UFLAG,RNOBJ,HYINIT,HYDEL,HXINIT,&
     &                HXDEL,APSTOP,SMAX,RSMAX,FOCL,OBJN(3),DELIMP,      &
     &                FPLANE,FAKEA,C(40),T(40),R(40),CONIC(40),FN(40,3),&
     &                FMASK(40),FAKEC(40),FAKEB(40),XDISP(40),YDISP(40),&
     &                TILTX(40),TILTY(40),TILTZ(40),ORDN(40,3),SIDE(40),&
     &                RDSPAC(40),Y0(40),SXY(40),SXNU(40),COEF(40,4),    &
     &                XMN(40),XMX(40),YMN(40),YMX(40),RX(40),CX(40),    &
     &                FREF(40),FREF0,WAVL(3)                            
      COMMON /COLLAT/ CLTRA(300),RADIMG,CVIMG,CONIMG,NPLANE,LATYPE,     &
     &                ICOL(3),NCOL,NSURF,IMODE,IPRINT,IPLTPR,IWVFLG(3), &
     &                IPR20,IREF,IJK,IALLPL                             
      COMMON /HEAD/   LINES,IPAGE,NSYS,LAMDA,NAME(160) 
      COMMON /PUPIL/  ENPUPR,ENPUPL,EXPUPR,EXPUPL 
      COMMON /SAGPAR/ YMAX,YMIN,DELY,REFCRV,CONST,WAVENM,ISRF 
      COMMON /CRED/ RPP(40),EPP(40),XW(812),HYSAVE,HXSAVE,              &
     &             TIMES,RD,OPTNE,ISW,IOPB                              
      COMMON /CMTF/   AN,RSC,DELTEN,DELTPL,DELTCC,                      &
     &                AR(40),DE(40),EN,SUM,AJ,TX(51),OPTNF,IOPC         
      COMMON /CSPOT/  XTMAX,XTMIN,YTMAX,YTMIN,AVGX,AVGY,RP(812),SPOTP,  &
     &                XK(812),YK(812),XKALL(7600),YKALL(7600),JSKIP,    &
     &                IOPA,NTHRU                                        
!                                                                       
!                                                                       
!                                                                       
      DIMENSION YW(812),TOR(40) 
      DIMENSION ZF(10),BLOB(10) 
      DIMENSION XSUM(10),YSUM(10),XSUMSQ(10),YSUMSQ(10) 
      DIMENSION QT(3),XX(3),XYZ(3) 
!                                                                       
      EQUIVALENCE (QT(1),QX),(QT(2),QY),(QT(3),QZ) 
      EQUIVALENCE (XX(1),XT),(XX(2),YT),(XX(3),ZT) 
      EQUIVALENCE (XYZ(1),X),(XYZ(2),Y),(XYZ(3),Z) 
      EQUIVALENCE (TRASH(1),XSMIN),(TRASH(2),XSMAX),(TRASH(3),YSMIN) 
      EQUIVALENCE (TRASH(4),YSMAX),(TRASH(5),FCODE),(TRASH(6),FXYJ) 
      EQUIVALENCE (TRASH(7),FXY),(TRASH(8),FXNU),(TRASH(9),FBY) 
      EQUIVALENCE (TRASH(10),FBNU),(Y0(1),TOR(1)) 
!                                                                       
!                                                                       
      DATA IZERO/0/,IONE/1/ 
      DATA ZERO/0.D0/,ONE/1.D0/,TWO/2.D0/,THREE/3.D0/,FOUR/4.D0/ 
      DATA FIVE/5.D0/,SIX/6.D0/,TEN/10.D0/,HUNDRD/100.D0/,THOUS/1000.D0/ 
      DATA EPS1/1.D-9/,EPS2/1.D-14/,PI/3.141592653589793D0/ 
!                                                                       
!                                                                       
        CALL GRAPH(2,1) 
        IF (ISW.EQ.1) GO TO 1580 
        XSMIN=ZERO 
        YSMIN=ZERO 
        XSMAX=ZERO 
        YSMAX=ZERO 
 1580  RETURN 
      END                                           
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!                                                                       
!     SUBROUTINE BESJN(X,N,BJN,ERRA, IER)                               
!                                                                       
!                                                                       
!                                                                       
!    PURPOSE                                                            
!        COMPUTE THE ORDINARY                                           
!        BESSEL FUNCTIONS OF ORDER N                                    
!              J (X)                                                    
!               N                                                       
!                                                                       
!    USAGE                                                              
!         CALL BESJN (X, N, BJN, ERRA, IER)                             
!                                                                       
!    PARAMETERS                                                         
!         BJN  - VALUE OF BESSEL FUNCTION                               
!         ERRA - ACCURACY                                               
!         IER  - ERROR FLAG                                             
!              = 1 => FATAL ERROR BESJN COULD NOT ACHIEVE DESIRED       
!                     PRECISION                                         
!                                                                       
!         N    - ORDER                                                  
!         X    - ARGUMENT                                               
!                                                                       
!    METHOD                                                             
!         USING THE RATIO X/DFLOAT(N)                                   
!         BESJN DIRECTS THE CALCULATION TO ONE OF 3 SUBROUTINES         
!     X/N < 0.50  => CALL BESJN0 FOR SUMMATION OF ASCENDING SERIES      
!     X/N > 100. => CALL BESJN1 FOR SUMMATION OF ASYMPTOTIC SERIES      
!      OTHERWISE  => CALL BESJN2 FOR DOWNWARDS RECURSION                
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    WRITTEN 4/12/84 BY JOHN PARKER OF SSAI UNDER CONTRACT WITH         
!    GSFC.                                                              
!                                                                       
!    REFERENCE: ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL         
!                FUNCTIONS (NATIONAL BUREAU OF STANDARDS)               
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BESJN(X,N,BJN,ERRA, IER) 
      IMPLICIT REAL*8 (A-H, O-Z) 
      DATA ERRMIN/1.D-13/, ZERO/0.0D0/, ONE/1.0D0/ 
      DATA TEN/10.0D0/ 
!                                                                       
      IER = 0 
      IER2 = 0 
      IER1 = 0 
      BJN = ZERO 
      ERR = DABS(ERRA) 
      IF(ERR .LE. ERRMIN) ERR= ERRMIN 
!                                                                       
!        DETERMINE PHASE                                                
!                                                                       
        PHASE = ONE 
      IF(((X .LT. ZERO) .OR. (N .LT. 0)) .AND.                          &
     &   (MOD(N,2) .NE. 0)) PHASE = -ONE                                
      IF ((X .LT. ZERO) .AND. (N.LT.0)) PHASE = ONE 
!                                                                       
!        CHANGE INPUT VBL NAME SO INPUT NOT DESTROYED                   
!                                                                       
      XA  = DABS(X) 
      NA  = IABS(N) 
      DNA1= DFLOAT(NA+1) 
!                                                                       
!        CHECK FOR SPECIAL VALUES                                       
!                                                                       
      IF( XA .GT. ZERO) GO TO 10 
!                                                                       
!         SPECIAL VALUES FOR X = 0                                      
!                                                                       
      BJN = ZERO 
      IF( NA .EQ. 0) BJN = ONE 
      GO TO 50 
!                                                                       
   10  CONTINUE 
       TEST = XA/DNA1 
       IF(TEST   .LE. 0.50D0) GO TO 20 
       IF(TEST   .GE. 1.D02 ) GO TO 30 
   12  CALL BESJN2 (XA, NA, BJN, ERR, IER) 
       IF(IER .EQ. 0) GO TO 50 
       IF( (IER2.NE.0) .OR. (IER3.NE.0)) GO TO 40 
!                                                                       
!    DOWNWARDS RECURSION FAILED TO CONVERGE.  TRY ASCENDING OR          
!    ASYMPTOTIC SERIES                                                  
!                                                                       
       IF(TEST.LE.ONE) GO TO 20 
       IF(TEST.GT.TEN) GO TO 30 
       GO TO 40 
   20  CALL BESJN0(XA, NA, BJN, ERR, IER2) 
       IF(IER2.EQ.0) GO TO 50 
       IF(IER.EQ.0) GO TO 12 
       IER = IER2 
       GO TO 40 
   30  CALL BESJN1(XA, NA, BJN, ERR, IER3) 
       IF(IER3.EQ.0) GO TO 50 
       IF(IER.EQ.0) GO TO 12 
       IER = IER3 
       GO TO 40 
   40  WRITE(6,10010) IER, X, N 
10010  FORMAT(/,' BESJN: ERROR. FAILURE TO CONVERGE',/                  &
     &          ' BESJN:                 ERROR FLAG IER = ',I5,         &
     &        /,' BESJN:',20X,'(X, N) =    ',1PD15.6,I5,/)              
   50  RETURN 
      END                                           
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            
      SUBROUTINE BESJN0 (X, N, BJN, ERR, IER) 
      IMPLICIT REAL*8 (A-H, O-Z) 
      DATA ZERO/0.0D0/, ONE/1.0D0/ 
!                                                                       
!       USE ASCENDING SERIES FOR SMALL ARGUMENTS                        
!                                                                       
      IF(IER .EQ. 100) WRITE(6, 1) 
    1 FORMAT(10X,' BESSEL FUNCTION ROUTINE: ASCENDING SERIES',          &
     & //, ' NUM. ', '   TERM VALUE  ', '          SUM  '/)             
           M= 0 
           XHALF = X/2.0D0 
           Y     = -XHALF*XHALF 
           DN    = DFLOAT(N) 
           DN1   = DFLOAT(N+1) 
           TERM  = ONE 
           BJN   = TERM 
   20 IF(M.GT.500) GO TO 50 
              M=M+1 
              DM = DFLOAT(M) 
              TERM  = TERM * Y /(DM*(DM + DN)) 
              BJN = BJN + TERM 
              IF(DABS(TERM/BJN).GT.ERR) GO TO 20 
      IF(N.EQ.0) GO TO 40 
      ATMP = ONE 
      DO 30 I=1,N 
          DI = DFLOAT(I) 
          ATMP = ATMP * (XHALF/DI) 
   30 END DO 
      BJN = BJN * ATMP 
   40 RETURN 
   50 IER = M 
      RETURN 
      END                                           
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                
      SUBROUTINE BESJN1 (X, N, BJN, ERR, IER) 
      IMPLICIT REAL*8 (A-H, O-Z) 
      DATA ZERO/0.0D0/, ONE/1.0D0/ 
      DATA PI/3.141592653589793238462643D0/ 
!                                                                       
!     ASYMPTOTIC SERIES FOR BJN                                         
!                                                                       
      PHASE = ONE 
       DN   = DFLOAT(N) 
       DMU = DN*DN *4.0D0 
       EIGHTX = 8.0D0*X 
       KTEST =  (N+N+1)/4 
       TERM = ONE 
       P= ONE 
       Q = ZERO 
       K = 0 
  230  K=K+1 
          OLDTRM = TERM 
          DK = DFLOAT(K) 
          TERM = TERM*(DMU-(DK+ONE)*(DK+ONE))/(DK*EIGHTX) 
          IF(MOD(K,2) .NE.0) GO TO 120 
            PHASE = -PHASE 
            P = P + PHASE * TERM 
          GO TO 130 
!                                                                       
!                                                                       
  120     Q = Q + PHASE*TERM 
          GO TO 130 
  130 CONTINUE 
!     WRITE(6,131) K,TERM,P,Q                                           
! 131 FORMAT('  K,TERM,P,Q=',I5,1P3D20.10,/)                            
      IF(K.LT.KTEST) GO TO 230 
      IF(P.NE.ZERO) ATMP = DABS(TERM/P) 
      IF(P.EQ.ZERO) ATMP = DABS(TERM) 
      IF ((ATMP .GT. ERR)                                               &
     &       .AND.(DABS(TERM) .LT. DABS(OLDTRM)))GO TO 230              
!                                                                       
!                                                                       
!                                                                       
      IF(DABS(TERM).GT.DABS(OLDTRM)) GO TO 140 
!                                                                       
!    ASYMPTOTIC SERIES SUCCESSFULLY CONVERGED                           
!                                                                       
       CHI = X - (DN +DN+ONE)*PI*0.25D0 
       BJN = P*DCOS(CHI)-Q*DSIN(CHI) 
       BJN = BJN*DSQRT(2.0D0/(X *PI)) 
       RETURN 
!                                                                       
!   ASYMPTOTIC SERIES FAILED TO CONVERGE                                
!                                                                       
  140  IF(K.LT.KTEST) GO TO 230 
       IER = 1 
       RETURN 
      END                                           
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
!                                                                       
!     SUBROUTINE BESJN2(X,N,BJN,ERRA, IER)                              
!                                                                       
!                                                                       
!                                                                       
!    PURPOSE                                                            
!        COMPUTE THE ORDINARY                                           
!        BESSEL FUNCTIONS OF ORDER N                                    
!              J (X)      BY DOWNWARDS RECURSION                        
!               N                                                       
!                                                                       
!    USAGE                                                              
!         CALL BESJN (X, N, BJN, ERRA, IER)                             
!                                                                       
!    PARAMETERS                                                         
!         BJN  - VALUE OF BESSEL FUNCTION                               
!         ERRA - ACCURACY                                               
!         IER  - ERROR FLAG                                             
!              = 1 => FATAL ERROR BESJN COULD NOT ACHIEVE DESIRED       
!                     PRECISION                                         
!              =100=> DETAILED PRINTS                                   
!         N    - ORDER                                                  
!         X    - ARGUMENT                                               
!                                                                       
!    METHOD                                                             
!        STARTING AT A LARGE ORDER, RECUR DOWNWARDS                     
!              USING                                                    
!      K (X) = 2*(N+1)*K   (X)/X - K(X)                                 
!       NN+1N+2                                                         
!                                                                       
!      THE K'S ARE EQUAL TO J'S WITHIN A CONSTANT FACTOR DETERMINED BY  
!       DNORM = K  + 2*(K  + K  + . . . )                               
!                0       2    4                                         
!                                                                       
!      THEN  J  = K  /DNORM                                             
!             N    N                                                    
!                                                                       
!    WRITTEN 4/12/84 BY JOHN PARKER OF SSAI UNDER CONTRACT WITH         
!    GSFC.                                                              
!                                                                       
!    REFERENCE: ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL         
!                FUNCTIONS (NATIONAL BUREAU OF STANDARDS)               
!                                                                       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BESJN2 (X, N, BJN, ERR, IER) 
      IMPLICIT REAL*8 (A-H, O-Z) 
      DATA ZERO/0.0D0/, ONE/1.0D0/ 
!                                                                       
      BJN    = ZERO 
      BSAVE  = ZERO 
      NSTART = MAX0(N + 10, IDINT(X + 10.0D0) ) 
      ICHECK = 0 
!                                                                       
! FOR A GIVEN NSTART, THIS LOOP IS PASSED THROUGH TWICE                 
!    TO CHECK ERRORS                                                    
!                                                                       
   40 M = NSTART 
      DNORM  = ZERO 
         BJMP1 = ZERO 
         BJMP2 = ZERO 
         BJM   = 1.D-28 
!               RECURSION                                               
   50    IF(M .LE. 0) GO TO 60 
              M=M-1 
              DM = DFLOAT(M) 
!  C                                                                    
              BJMP2 = BJMP1 
              BJMP1 = BJM 
!  C                                                                    
              BJM   =(2.0D0 * (DM+ONE ) * BJMP1/X)- BJMP2 
!                                                                       
              IF (M.EQ. N) BJN = BJM 
              IF (MOD(M,2).NE.0) GO TO 50 
              IF (M .NE. 0) DNORM = BJM + DNORM + BJM 
              IF (M .EQ. 0) DNORM = BJM + DNORM 
          GO TO 50 
!                NORMALIZE                                              
   60         BJN = BJN/DNORM 
              IF( ICHECK.EQ.1) GO TO 70 
                 ICHECK = 1 
                 BSAVE = BJN 
                 NSTART = NSTART + 2 
                 GO TO 40 
!                                                                       
!          ACCURATE???                                                  
!                                                                       
   70        IF(DABS(BSAVE-BJN) .LE. DABS(BJN*ERR)) GO TO 90 
!                                                                       
!           TRY AGAIN                                                   
!                                                                       
   80            BSAVE = BJN 
                 ICHECK = 1 
                 NSTART = NSTART  + 5 
                 IF(NSTART-N  .LT. 100 ) GO TO 40 
!                                                                       
!       FATAL  ERROR TRAP                                               
!                                                                       
                     IER = NSTART 
   90 RETURN 
      END                                           
