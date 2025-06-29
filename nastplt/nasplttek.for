      PROGRAM NasPltTek

      CHARACTER*35 FILNAM
      INTEGER PLTREC(750)
      LOGICAL LIBM, LSKIP
      DATA LSKIP/.FALSE./
C
      COMMON /NASTCOM/ ICMND,ICNTRL,IR,IS,IT,IU
      COMMON /SCOPE/ IBAUD
C
C     GET INPUT FILENAME & TYPE (VAX OR IBM)
10        CALL NAMTYP(FILNAM,LIBM)
C     OPEN PLOTFILE
!!!       OPEN(UNIT=13,NAME=FILNAM,TYPE='OLD',
!!!    *         ACCESS='SEQUENTIAL',RECORDTYPE='FIXED',
!!!    *         RECORDSIZE=3000,BLOCKSIZE=3000,
!!!    *         FORM='FORMATTED',ORGANIZATION='SEQUENTIAL',
!!!    *               READONLY,CARRIAGECONTROL='NONE', ERR=10)
!!!     *              READONLY,ERR=10)
      OPEN(UNIT=13,FILE=filnam,STATUS='OLD',ACTION='READ')

C     INITIALIZE PLOT SCREEN
      WRITE(*,*) 'ENTER TERMINAL BAUD RATE:'
      READ(*,*) IBAUD
          IBAUD = IBAUD/10  !! CHANGE TO CHAR/SEC
          CALL INITT(IBAUD)
C
      DO 200 IREC=1,1000000
C         GET NEXT PLOT RECORD
              READ(13,2000,END=250) PLTREC
C         CONVERT IBM RECORD
              IF (LIBM) CALL CNVIBM(PLTREC)
          DO 100 I=1,100
C             GET NEXT COMMAND (TOTAL OF 100 COMMANDS/RECORD)
                  CALL GETCMD(PLTREC)
C             EXECUTE COMMAND
                  IF (ICMND .EQ. 0)  CONTINUE     !! NULL COMMAND
                  IF (ICMND .EQ. 1)  CALL BGNPLT(LSKIP) !START NEW PLOT
                  IF (LSKIP)  GOTO 100
                  IF (ICMND .EQ. 2)  CONTINUE     !! SELECT CAMERA
                  IF (ICMND .EQ. 3)  CONTINUE     !! SKIP FRAME
                  IF (ICMND .EQ. 4)  CALL DRWCHR  !! TYPE A CHARACTER
                  IF (ICMND .EQ. 5)  CALL DRWLIN  !! DRAW A LINE
                  IF (ICMND .EQ. 6)  CALL DRWLIN  !! DRAW AN AXIS
                  IF (ICMND .EQ. 15) CALL DRWLIN  !! DRAW A LINE
                  IF (ICMND .EQ. 16) CALL DRWLIN  !! DRAW AN AXIS
100           CONTINUE
200       CONTINUE
C
C     END OF PLOT FILE
250       CALL RESET
          CALL SEETW(MINX,MAXX,MINY,MAXY)
          CALL FINITT(MINX,MAXY)
          STOP 'END OF PLOT FILE'
C
1000  FORMAT(1X,A25,1X,$)
2000  FORMAT(750A4)
      END Program NasPltTek


      SUBROUTINE BGNPLT(LSKIP)
C
      CHARACTER*1 YESNO
      LOGICAL LSKIP
      COMMON /NASTCOM/ ICMND,ICNTRL,IR,IS,IT,IU
      COMMON /SCOPE/ IBAUD
C
C     DUMP BUFFER FOR PREVIOUS PLOT
          CALL RESET
          CALL SEETW(IMINX,IMAXX,IMINY,IMAXY)
          CALL HOME
          CALL ANMODE
          PAUSE 'TYPE "CONTINUE" FOR NEXT PLOT:'
          LSKIP = .FALSE.
          WRITE(*,*) 'SKIP NEXT FRAME?  (Y OR N):'
          READ(*,*) YESNO
          IF (YESNO .EQ. 'Y')  LSKIP = .TRUE.
C
C     INITIALIZE FOR NEXT PLOT
          CALL INITT(IBAUD)
          CALL SEETW(IMINX,IMAXX,IMINY,IMAXY)
          IMAXY = NINT(.95 * FLOAT(IMAXY))
          XMAX = IS
          YMAX = IT
          RATIO = YMAX/XMAX
          LENY = IMAXY
          LENX = NINT( FLOAT(LENY) / RATIO )
          IF (LENX .GT. IMAXX)  THEN
              LENX = IMAXX
              LENY = NINT( FLOAT(LENX) * RATIO )
              END IF
C         SET SCREEN WINDOW
              MINX = (IMAXX-LENX)/2
              MINY = (IMAXY-LENY)/2
              CALL SWINDO(MINX,LENX,MINY,LENY)
C         SET VIRTUAL (USER) WINDOW
              CALL DWINDO(0.0,XMAX,0.0,YMAX)
      RETURN
C
1000  FORMAT(1X,A15,' (PLOT NO.',I4,')',A11,1X,$)
2000  FORMAT(A1)
      END

      SUBROUTINE CNVIBM(RECORD)
      INTEGER RECORD(750)
      BYTE BYTE(4), BTEMP(4)
      EQUIVALENCE (IWORD,BYTE(1)) , (ITEMP,BTEMP(1))
C
C     CONVERT FROM IBM FORMAT TO VAX FORMAT:
C     REVERSE THE ORDER OF THE BYTES IN EACH WORD OF THE RECORD
          DO 100 I=1,750
              ITEMP = RECORD(I)
              BYTE(1) = BTEMP(4)
              BYTE(2) = BTEMP(3)
              BYTE(3) = BTEMP(2)
              BYTE(4) = BTEMP(1)
              RECORD(I) = IWORD
100           CONTINUE
      RETURN
      END

      SUBROUTINE DRWCHR
C
C     THIS (CHARACTER-TYPING) CAPABILITY IS NOT YET IMPLEMENTED
C     FOR THIS PLOTTING PROGRAM.
C
      RETURN
      END

      SUBROUTINE DRWLIN
C
      COMMON /NASTCOM/ ICMND,ICNTRL,IR,IS,IT,IU
C
C     MOVE BEAM TO STARTING POINT OF LINE
          XSTART = IR
          YSTART = IS
          CALL MOVEA(XSTART,YSTART)
C     DRAW LINE FROM STARTING POINT TO END POINT
          XFINAL = IT
          YFINAL = IU
          CALL DRAWA(XFINAL,YFINAL)
      RETURN
      END

      SUBROUTINE GETCMD(PLTREC)
      IMPLICIT INTEGER (A-Z)
      DIMENSION PLTREC(750),Q(30),MASK(4)
      COMMON /NASTCOM/ CMND,CNTRL,R,S,T,U
      COMMON /RECPNT/ KWORD,KBYTE
C     EQUATE COMMAND INTEGERS TO MNEMONIC VARIABLES
          EQUIVALENCE
     *     (Q( 1),PC), (Q( 2),CI),
     *     (Q( 3),R4), (Q( 4),R3), (Q( 5),R2), (Q( 6),R1), (Q( 7),R0),
     *     (Q( 8),S4), (Q( 9),S3), (Q(10),S2), (Q(11),S1), (Q(12),S0),
     *     (Q(13),T4), (Q(14),T3), (Q(15),T2), (Q(16),T1), (Q(17),T0),
     *     (Q(18),U4), (Q(19),U3), (Q(20),U2), (Q(21),U1), (Q(22),U0)
C     INITIALIZE BYTE MASKS, AND PLTREC WORD & BYTE POINTERS
          DATA MASK(1) /'000000FF'X/ ,
     *         MASK(2) /'0000FF00'X/ ,
     *         MASK(3) /'00FF0000'X/ ,
     *         MASK(4) /'FF000000'X/
          DATA KWORD /1/, KBYTE /4/
C
C     BREAK OUT 30 INTEGERS OF NEXT COMMAND
          DO 100 I=1,30
              Q(I) = IAND( PLTREC(KWORD) , MASK(KBYTE) )
              Q(I) = ISHFT( Q(I) , -8*(KBYTE-1) )
C             SET POINTERS TO GET NEXT BYTE
                  KBYTE = KBYTE - 1
                  IF (KBYTE .EQ. 0) KWORD = KWORD + 1
                  IF (KBYTE .EQ. 0) KBYTE = 4
                  IF (KWORD .EQ. 751) KWORD = 1
100           CONTINUE
C     CALCULATE CMND,CNTRL,R,S,T,&U
          CMND = PC
          CNTRL = CI
          R = R0 + 10*R1 + 100*R2 + 1000*R3 + 10000*R4
          S = S0 + 10*S1 + 100*S2 + 1000*S3 + 10000*S4
          T = T0 + 10*T1 + 100*T2 + 1000*T3 + 10000*T4
          U = U0 + 10*U1 + 100*U2 + 1000*U3 + 10000*U4
      RETURN
      END

      SUBROUTINE NAMTYP(FILNAM,LIBM)
      CHARACTER(LEN=*),INTENT(OUT):: filnam
      LOGICAL,INTENT(OUT):: libm

      CHARACTER(LEN=1):: qflag

C     GET INPUT FILENAME
      WRITE(*,*) 'ENTER NAME OF INPUT PLT2 FILE: '
      READ(*,'(A)') FILNAM
C     IS THIS AN IBM FORMAT FILE ?
      WRITE(*,*) ' CONVERT FILE FROM IBM FORMAT ? (Y OR N): '
      READ(*,'(A)') QFLAG
      LIBM = (QFLAG .EQ. 'Y')
      RETURN
      END Subroutine NamTyp
