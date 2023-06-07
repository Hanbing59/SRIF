      SUBROUTINE TMPRNT (R,N,NAMES,TEXT,CODE,LABEL, 
     1                   PRTNAM,RNUM,IOUNIT)
C   
C        TO DISPLAY A VECTOR STORED UPPER TRIANGULAR MATRIX IN A
C        TWO-DIMENSIONAL TRIANGULAR FORMAT  
C   
C ** INPUTS:
C   
C     R(N*(N+1)/2)  VECTOR CONTAINING UPPER TRIANGULAR MATRIX   
C   
C     N             DIMENSION OF MATRIX 
C   
C **** NOTE: IF N IS 0 THEN ONLY THE TITLE LINE WILL
C            BE PRINTED. NO DIAGNOSTIC MESSAGE WILL BE PRODUCED.
C   
C   NAMES         ANY CHARACTER TYPE (CHARACTER*K)  
C                 BUT ONLY THE FIRST 16 ARE SIGNIFICANT.
C                 AESTHETICALLY  PLEASING RESULTS ARE OBTAINED  
C                 WHEN NAMES IS CHARACTER*16 AND (ALL 16 BYTES  
C                 ARE USED FOR EACH ELEMENT) OR (EACH ELEMENT   
C                 IS CENTERED WITHIN ITS 16 BYTE FIELD).
C   
C     TEXT( )     AN ARRAY OF ASCII CHARACTERS TO BE PRINTED
C                 A TITLE PRECEDING THE MATRIX. 
C   
C     CODE        (CHARACTER*(*)) IS A LIST OF PARAMETERS,  
C                 WHICH ARE INTERPRETED IN A POSITION INDEPENDENT   
C                 FASHION. ANY CHARACTERS NOT DEFINED IN THIS LIST  
C                 WILL HAVE NO EFFECT; BUT OTHER FEATURES MAY BE
C                 DEFINED AT A LATER DATE, SO EXTRANEOUS CHARACTERS 
C                 SHOULD BE AVOIDED.
C   
C      "P"        A PAGE EJECT WILL PRECEDE PRINTING
C   
C      "C"        CONDENSED PRINT   
C   
C                   (I.E. LESS PRECISION, MORE COLUMNS PER ROW) 
C                   THE DEFAULT (' ') CAUSES PRINTING TO BEGIN  
C                   AT THE NEXT LINE AND DISPLAYS 6 COLUMNS WITH
C                   8 DIGITS OF PRECISION. CONDENSED MODE CONSISTS  
C                   OF 10 COLUMNS OF 4 DIGIT
C                   PRECISION. IF THE NAMC CHARACTER NAMES ARE  
C                   LONGER THAN 8 CHARACTERS LONG OTHER CODE
C                   PARAMETERS CONTROLS HOW EACH NAME IS TO 
C                   BE PRINTED (SINCE ONLY 8 CHARACTERS CAN APPEAR  
C                   PER LINE).  
C   
C    "BK='C'"     WHERE 'C' IS ANY CHARACTER, INCLUDING SPACE.  
C   
C                   CONTROL FOR HANDLING NAMC NAMES WITH LENGTH 
C                   LONGER THAN 8. 'C' IS USED TO REQUEST   
C                   TO PUT THE REMAINING CHARACTERS ON THE  
C                   NEXT LINE. THE DEFAULT IS TO TREAT A
C                   BLANK AS THE REQUEST CHARACTER. NOTE:   
C                   THE REQUIRED SYNTAX IS BK=THING, BLANKS 
C                   ARE SIGNIFICANT. IT SHOULD ALSO BE NOTED
C                   THAT THIS OPTION IS IRRELEVANT UNLESS THE   
C                   CONDENSED PRINT OPTION IS IN EFFECT.
C   
C  EXAMPLES:   N = 3, CHARACTER*16 NAMES ='X-POSITION', 
C                     'X VELOCITY','X_ACCELERATION' 
C   
C    CODE =    'PC' WILL PRODUCE A NEW PAGE AND  CONDENSED PRINT
C    CODE =    'C'  WILL PRODUCE CONDENSED PRINT ONLY   
C   
C     THE COLUMN NAMES WILL BE: 
C   
C     X-POSITI   X            X_ACCELE  
C     ON         VELOCITY     RATION
C   
C  SINCE "BK=" WAS NOT SPECIFIED, SPACE TRIGGERED THE BREAK 
C   
C  IF CODE =   'C BK=_' THE RESULTING COLUMN NAMES WOULD BE 
C   
C     X-POSITI  X VELOCI      X 
C     ON        TY            ACCELERA  
C   
C   
C  SINCE THE "BREAK" CHARACTER IS NOT PRINTED (THE FIRST TIME IT
C  OCCURS) IT IS POSSIBLE (AND MEANINGFUL) TO USE CHARACTER*17  
C  NAMES, IF ONE OF THE CHARACTERS IS THE "BREAK" CHARACTER.
C   
C   
C     LABEL        TRUE TO PRINT PARAMETER NAMES
C   
C     PRTNAM   LABEL ASSOCIATED WITH REAL NUMBER, RNUM, TO BE   
C              PRINTED FOLLOWING THE TITLE  
C   
C                   IF .EQ. ' '  THEN   
C                      DO NOT PRINT ANYTHING FOLLOWING THE TITLE.   
C   
C     RNUM       NUMBER TO BE PRINTED WHEN PRTNAM   
C                .NE. ' '   
C   
C     IOUNIT       UNIT NUMBER OF OUTPUT DEVICE 
C                  IF IOUNIT .EQ. 0 NO OUTPUT IS PRODUCED.  
C   
C     REQUIRED SUBROUTINES: 
C   
C              MAX  <FORLIB>    INDEX  <FORLIB> 
C              MIN  <FORLIB>
C   
C     COGNIZANT PERSONS:
C   
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN 
C              FACTORIZED ESTIMATION APPLICATIONS   
C   
C              18 MARCH 1984     VERSION 0.00.00   FORTRAN  '77 VERSION 
C              17 JUNE  1984     VERSION 1.01.00
C              20 JUNE  1985     VERSION 2.00.00   16 CHARACTER VERSION 
C   
C ******************************************************************
C *                                                                *
C *   SUBROUTINE TMPRNT IS AN ADAPTATION OF AN EARLIER             *
C *   TRIANGULAR MATRIX PRINT SUBROUTINE                           *
C *   AUTHORED MY MELBA W. NEAD  OF THE JET  PROPULSION            *
C *   LABORATORY.                                                  *
C *                                                                *
C ******************************************************************
C *   This library is licensed under the Creative Commons          *
C *   Attribution 3.0 Unported License.                            *
C *   To view a copy of this license,                              *
C *   visit http://creativecommons.org/licenses/by/3.0/            *
C *   or send a letter to Creative Commons, 171 Second Street,     *
C *   Suite 300, San Francisco, California, 94105, USA.            *
C *                                                                *
C ******************************************************************
C   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   

      DIMENSION R(*),     B(10) 
      CHARACTER*(*) TEXT, PRTNAM, CODE, NAMES(N)
      CHARACTER*8 V(4), VFMT(6),    V6(2), V2I, 
     1                  VFMT10(10), V10(2),TMPNAM(2,10) 
      CHARACTER*9 V2A   
      CHARACTER*1 BREAK 
C   
      INTEGER L(10),LIST(10)
C   
      LOGICAL LABEL,CONDEN,LSEC,LBLNK   
C   
      DATA  V/'('' ''','A16,1PX,',' ','D17.7)'/,
     1 VFMT/'6','17X,5','34X,4','51X,3','68X,2','85X,1'/
     2 V2A/'A16,1P1X,'/,V2I/'I5,1P2X,'/ 
C   
      DATA V10/'('' '',','D11.3)'/,V6/'('' '',2X,','D17.7)'/
C   
      DATA VFMT10/'11','11X,10','22X,9','33X,8','44X,7','55X,6',
     1           '66X,5','77X,4','88X,3','99X,2'/   
C   
      IF (IOUNIT .EQ. 0) RETURN 
C   
      CONDEN = INDEX(CODE,'C') .NE. 0 .OR. INDEX(CODE,'c') .NE. 0   
      IF ((INDEX(CODE,'P') .NE. 0) .OR. ((INDEX(CODE,'p').NE. 0))) THEN 
         WRITE(IOUNIT,'(''1'')')
      ELSE  
         WRITE(IOUNIT,*) ' '
         WRITE(IOUNIT,*) ' '
      ENDIF 
      IF (PRTNAM.EQ.' ') THEN   
          WRITE (IOUNIT,*) TEXT 
      ELSE  
          WRITE (IOUNIT,9) TEXT,PRTNAM,RNUM 
 9        FORMAT(' ',A,'     ',A,' = ',1PD14.7) 
      END IF
C   
      IF (N .EQ. 0) RETURN  
C   
      I = INDEX(CODE,'BK=') 
      J = INDEX(CODE,'BK=') 
      I = MAX(I,J)  
      IF (I.EQ. 0) THEN 
         BREAK = ' '
      ELSE  
         BREAK = CODE(I+3:I+4)  
      ENDIF 
C   
C        M1,M2    ROW LIMITS FOR EACH PRINT SEQUENCE
C        N1,M2    COL LIMITS FOR EACH LINE OF PRINT 
C        L(I)     LOC OF EACH COLUMN IN A ROW   
C        KT       ROW COUNTER   
C   
C  * * * * *   INITIALIZE COUNTERS  
C   
C   
      IF (CONDEN) THEN  
C                                  SET FOR 10 COLUMN FORMAT 
        J0=10   
        V(1)=V10(1) 
        V(4)=V10(2) 
      ELSE  
C                                  SET FOR 6 COLUMN FORMAT  
        J0=6
        V(1)=V6(1)  
        V(4)=V6(2)  
      ENDIF 
C   
      M1=1  
      M2=J0 
      N1=1  
      KT=0  
      J0P1=J0+1 
      J0M1=J0-1 
C   
C     ASSUME THAT NAMES WERE REQUESTED  
C   
      V(2)=V2A  
C   
10    IF (M2.GT.N) M2=N 
C   
      WRITE(IOUNIT,*)' '
      IF (LABEL) THEN   
         IF (CONDEN) THEN   
            LSEC  = .FALSE. 
            LBLNK = .FALSE. 
            IPTR  = 0   
            DO 20 I = N1, M2
               J = INDEX(NAMES(I),BREAK)-1  
               J = MIN(J,8) 
               IF (J .LE. 0) THEN   
                  J = LEN(NAMES(I)) 
                  IF (J .GT. 8) J = J/2 
               ELSE 
                  LBLNK = .TRUE.
               ENDIF
               IPTR           = IPTR+1  
               TMPNAM(1,IPTR) = NAMES(I) (1:J)  
               TMPNAM(2,IPTR) = ' ' 
               JEND           = MIN(17,LEN(NAMES(I)))   
               IF (JEND .GE. J) THEN
                  IF (LBLNK .OR. (J .GE. JEND/2))        J = J+1
                  IF (INDEX(NAMES(I)(J:J),BREAK) .NE. 0) J = J+1
                  IF(JEND .GE. J)    THEN   
                     LSEC           = .TRUE.  
                     TMPNAM(2,IPTR) = NAMES(I)(J:JEND)  
                  ENDIF 
               ENDIF
  20         CONTINUE   
             WRITE(IOUNIT,*) ' '
             WRITE(IOUNIT,540)           (TMPNAM(1,I),I=1,IPTR) 
             IF (LSEC) WRITE(IOUNIT,540) (TMPNAM(2,I),I=1,IPTR) 
         ELSE   
             WRITE(IOUNIT,200) (NAMES(I),I=N1,M2)   
         ENDIF  
      ELSE  
C                      (NOT LABEL)  
         M         = N1 
         V(2)      = V2I
         L2        = M2-N1+1
         DO 30 I=1,L2   
           LIST(I) = M  
           M       = M+1
   30    CONTINUE   
         IF (CONDEN) THEN   
            WRITE(IOUNIT,210)  (LIST(I),I=1,L2) 
         ELSE   
            WRITE (IOUNIT,310) (LIST(I),I=1,L2) 
         ENDIF  
      ENDIF 
      WRITE(IOUNIT,*) ' '   
  200 FORMAT(' ',19X,6(A16,1X)) 
  210 FORMAT(' ',8X,10(I6,5X))  
  310 FORMAT(' ',14X,6(I6,11X)) 
  540 FORMAT(' ',16X,10(2X,A8,1X))  
C   
C   * * * * * * 
C   
      DO 190 IC=M1,M2   
        K=1 
        IF (IC.LE.(KT*J0)) GO TO 60 
        JJ=0
        DO 50 J=1,IC
   50     JJ=JJ+J   
        L(K)=JJ 
        I1=IC-KT*J0 
        IF (I1.EQ.J0) GO TO 90  
        GO TO 70
   60   CONTINUE
C   
        I1=1
        L(K)=L(K)+1 
   70   CONTINUE
        DO 80 I=I1,J0M1 
          K=K+1 
          II=I+KT*J0
C                                    @ OBTAIN COL INDEX FOR ROW 
   80     L(K)=L(K-1)+II
   90   CONTINUE
C   
        I2=MIN0(J0P1,(M2+1-KT*J0))-I1   
C   
C      ADD CORRECT NUMBER OF SPACES IN FRONT OF EACH ENTRY  
C          CORRECT NUMBER DEPENDS ON FORMAT (6 OR 10)   
      IF (CONDEN) THEN  
         V(3)=VFMT10(I1)
      ELSE  
         V(3)=VFMT(I1)  
      ENDIF 
C   
      DO 101 JB=1,I2
         LI=L(JB)   
         B(JB)=R(LI)
101   CONTINUE  
C   
C     WRITE PARAMETER NAMES OR NUMBERS DEPENDING ON NAMES   
C   
        IF( LABEL ) THEN
            WRITE (IOUNIT,FMT=V) NAMES(IC),(B(I),I=1,I2)
        ELSE
            WRITE (IOUNIT,FMT=V) IC,(B(I),I=1,I2)   
       ENDIF
C   
  190 CONTINUE  
      IF (M2.EQ.N) RETURN   
      N1=M2+1   
      M2=M2+J0  
      KT=KT+1   
      WRITE(IOUNIT,*) ' '   
C   
      IF (PRTNAM.EQ.' ') THEN   
          WRITE(IOUNIT,*) TEXT  
      ELSE  
          WRITE (IOUNIT,19) TEXT,PRTNAM,RNUM
      END IF
   19 FORMAT(' (CONTINUE) ',A,'     ',A,' = ',1PD14.7)  
C   
      GO TO 10  
      END   
