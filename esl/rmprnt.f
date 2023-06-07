      SUBROUTINE RMPRNT(A,IA,NR,NC,NAMR,NAMC,TEXT,CODE,LBLROW,  
     1                  LBLCOL,PRTNAM,RNUM,IOUNIT)  
C   
C         SUBROUTINE TO PRINT A RECTANGULAR VECTOR STORED MATRIX
C   
C *** INPUTS:   
C   
C     A(M*N)   VECTOR STORED RECTANGULAR MATRIX, PRINT AS A(M,N)
C   
C     IA       ROW DIMENSION OF A (AS DEFINED IN THE CALLING
C              PROGRAM) 
C   
C     NR       NUMBER OF ROWS OF A TO BE PRINTED
C   
C     NC       NUMBER OF COLS OF A TO BE PRINTED
C   
C **** NOTE: IF NR OR NC IS 0 THEN ONLY THE TITLE LINE WILL 
C            BE PRINTED. NO DIAGNOSTIC MESSAGE WILL BE PRODUCED 
C   
C     NAMR(NR) NAMES OF PARAMETERS FOR ROWS OF A
C   
C     NAMC(NC) NAMES OF PARAMETERS FOR COLS OF A
C   
C   *** NAMR AND NAMC MAY BE OF ANY CHARACTER TYPE (CHARACTER*K)
C       BUT ONLY THE FIRST 16 ARE SIGNIFICANT. NAMR AND NAMC
C       MAY BE OF DIFFERENT LENGTHS. AESTHETICALLY  PLEASING
C       RESULTS ARE OBTAINED WHEN NAMC IS CHARACTER*16 AND  
C       (ALL 16 BYTES ARE USED FOR EACH NAMC) OR (EACH NAMC 
C       IS CENTERED WITHIN ITS 16 BYTE FIELD).  
C   
C   
C     TEXT()   CHARACTER*(*) AN ARRAY OF ASCII  
C              CHARACTERS TO BE PRINTED AS A TITLE  
C              PRECEDING THE MATRIX.
C   
C     CODE     (CHARACTER*(*)) IS A LIST OF PARAMETERS, 
C              WHICH ARE INTERPRETED IN A POSITION INDEPENDENT  
C              FASHION. ANY CHARACTERS NOT DEFINED IN THIS LIST 
C              WILL HAVE NO EFFECT; BUT OTHER FEATURES MAY BE   
C              DEFINED AT A LATER DATE, SO EXTRANEOUS CHARACTERS
C              SHOULD BE AVOIDED.   
C   
C      "P"     A PAGE EJECT WILL PRECEDE PRINTING   
C   
C      "C"     CONDENSED PRINT  
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
C    "BK='C'"  WHERE 'C' IS ANY CHARACTER, INCLUDING SPACE. 
C              CONTROL FOR HANDLING NAMC NAMES WITH LENGTH LONGER   
C              THAN 8. 'C' IS USED TO REQUEST TO PUT THE
C              REMAINING CHARACTERS ON THE NEXT LINE.   
C              THE DEFAULT IS TO TREAT A BLANK AS THE   
C              REQUEST CHARACTER. NOTE: THE REQUIRED SYNTAX 
C              IS BK=THING, BLANKS ARE SIGNIFICANT. IT SHOULD   
C              ALSO BE NOTED THAT THIS OPTION IS IRRELEVANT 
C              UNLESS THE CONDENSED PRINT OPTION IS IN EFFECT.  
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
C     LBLROW   LOGICAL FLAG, IF TRUE PRINT ROW NAMES
C                            IF FALSE PRINT NUMBERED ROWS   
C                            IF FALSE AND THERE IS ONLY ONE 
C                               MATRIX ROW TO PRINT THE 
C                               NUMBER LABEL IS SUPPRESSED. 
C   
C     LBLCOL   LOGICAL FLAG, IF TRUE PRINT COL NAMES
C                            IF FALSE PRINT NUMBERED COLS.  
C   
C     PRTNAM   CHARACTER*(*) LABEL ASSOCIATED WITH REAL NUMBER, 
C              RNUM, TO BE PRINTED FOLLOWING THE TITLE. 
C              IF .EQ. ' ' THEN RNUM WILL NOT BE PRINTED.   
C   
C     RNUM     NUMBER TO BE PRINTED 
C   
C     IOUNIT   UNIT NUMBER OF OUTPUT DEVICE 
C              IF IOUNIT .EQ. 0 THEN NO OUTPUT IS PRODUCED  
C   
C *** OUTPUT:  
C   
C     NO VARIABLES ARE CHANGED FOR OUTPUT.  
C   
C  *** REQUIRED SUBROUTINE: 
C   
C     LEN (FORLIB)        INDEX (FORLIB)
C   
C     COGNIZANT PERSON: 
C   
C              KEITH H. BIERMAN 
C              FACTORIZED ESTIMATION APPLICATIONS   
C              JULY 5, 1981 
C   
C              REVISED FROM MNPRNT     18 MARCH 1984
C              VERSION 0.00.00         FORTRAN '77 VERSION  
C              REVISED 17 JUNE 1984  VERSION 0.01.00
C              REVISED 17 MAY  1985  VERSION 1.00.00
C              REVISED 11 JUNE 1985  VERSION 2.00.00
C   
C ******************************************************************
C *                                                                *
C *   SUBROUTINE RMPRNT IS AN ADAPTATION OF AN EARLIER             *
C *   RECTANGULAR MATRIX PRINT SUBROUTINE                          *
C *   AUTHORED BY MELBA W. NEAD  OF THE JET  PROPULSION            *
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

      DIMENSION A(*)
      LOGICAL LBLCOL,LBLROW,CONDEN,LSEC,LBLNK   
      CHARACTER*17 BLANKS   
      CHARACTER*8 TMPNAM(2,10)  
      CHARACTER*(*) NAMR(NR),NAMC(NC),PRTNAM,CODE   
      DOUBLE PRECISION V(10)
      CHARACTER*(*) TEXT
      CHARACTER*1 BREAK 
      INTEGER II(10)
C   
      DATA BLANKS /'                '/  
C   
      IF (IOUNIT .EQ. 0) RETURN 
C   
      CONDEN = INDEX(CODE,'C') .NE. 0 .OR. INDEX(CODE,'c') .NE. 0   
C   
      IF ((INDEX(CODE,'P').NE.0) .OR. ((INDEX(CODE,'p').NE.0))) THEN
         WRITE (IOUNIT,'(1H1)') 
      ELSE  
         WRITE (IOUNIT,*)   
         WRITE (IOUNIT,*)   
      END IF
      I = INDEX(CODE,'BK=') 
      J = INDEX(CODE,'BK=') 
      I = MAX(I,J)  
      IF (I.EQ. 0) THEN 
         BREAK = ' '
      ELSE  
         BREAK = CODE(I+3:I+4)  
      ENDIF 
C   
      IF(PRTNAM.EQ.' ') THEN
         WRITE(IOUNIT,*) TEXT   
      ELSE  
         WRITE(IOUNIT,9) TEXT,PRTNAM,RNUM   
9        FORMAT(' ',A,'     ',A,' = ',1PD14.7)  
      END IF
C   
      IF ((NR .EQ. 0) .OR. (NC .EQ. 0)) RETURN  
C   
      M   = NC  
      IC  =  1  
      IRR =  0  
C   
   10 L=M   
C   
      IF (CONDEN) THEN  
         IF (M.GT.10) L=10  
      ELSE  
          IF (M.GT.6) L=6   
      ENDIF 
      LL=IC+L-1 
      DO 30 I=1,L   
         II(I)=IC+I-1   
   30 CONTINUE  
C   
C     IF COLUMNS NAMES WERE REQUESTED, AND CONDENSED PRINT  
C     THEN WE NEED TO COPY NAMC INTO TEMP 8 CHAR ARRAY  
C   
      WRITE(IOUNIT,*) ' '   
      IF (LBLCOL) THEN  
         IF (CONDEN) THEN   
            LSEC  = .FALSE. 
            LBLNK = .FALSE. 
            IPTR  = 0   
            DO 35 I = IC, LL
               J = INDEX(NAMC(I),BREAK)-1   
               J = MIN(J,8) 
               IF (J .LE. 0) THEN   
                  J = LEN(NAMC(I))  
                  IF (J .GT. 8) J = J/2 
               ELSE 
                  LBLNK = .TRUE.
               ENDIF
               IPTR = IPTR + 1  
               TMPNAM(1,IPTR) = NAMC(I)(1:J)
               JEND = MIN(17,LEN(NAMC(I)))  
               TMPNAM(2,IPTR) = ' ' 
               IF (JEND .GE. J) THEN
                  IF (LBLNK .OR. (J .GE. JEND/2))       J = J+1 
                  IF (INDEX(NAMC(I)(J:J),BREAK) .NE. 0) J = J+1 
                  IF(JEND .GE. J)     THEN  
                     LSEC           = .TRUE.  
                     TMPNAM(2,IPTR) = NAMC(I)(J:JEND)   
                  ENDIF 
               ENDIF
 35         CONTINUE
            WRITE(IOUNIT,*) ' ' 
            WRITE(IOUNIT,540) (TMPNAM(1,I),I=1,IPTR)
            IF (LSEC) WRITE(IOUNIT,540) (TMPNAM(2,I),I=1,IPTR)  
         ELSE   
            WRITE(IOUNIT,200) (NAMC(I), I=IC,LL)
         ENDIF  
      ELSE  
C                          (NOT LBLCOL) 
         IF (CONDEN) THEN   
            WRITE (IOUNIT,210) (II(I),I=1,L)
         ELSE   
            WRITE (IOUNIT,310) (II(I),I=1,L)
         ENDIF  
      ENDIF 
C   
      WRITE(IOUNIT,*) ' '   
C   
  200 FORMAT (' ',17X,6(A16,3X))
  210 FORMAT (' ',17X,10(I6,5X))
  310 FORMAT (' ',18X,6(I6,13X))
  540 FORMAT (' ',16X,10(2X,A8,1X)) 
C   
      DO 50 I=1,NR  
         IR=IRR+I   
         DO 40 J=1,L
            V(J)=A(IR)  
            IR=IR+IA
   40    CONTINUE   
C   
C     IF ONLY ONE ROW IS TO BE PRINTED, AND NO NAMES ARE SELECTED   
C        PRINT BLANKS INSTEAD OF THE NUMBER 1.  
C   
      IF ( (.NOT. LBLROW) .AND. (NR .EQ. 1)) THEN   
         IF (CONDEN) THEN   
            WRITE (IOUNIT,220) BLANKS,(V(J),J=1,L)  
         ELSE   
            WRITE (IOUNIT,320) BLANKS,(V(J),J=1,L)  
         ENDIF  
      ELSE  
C   
C        PRINT CURRENT ROW WITH OR WITHOUT SIDE HDRS
C        DEPENDING ON LBLROW
C   
         IF (CONDEN) THEN   
            IF (LBLROW) THEN
               WRITE (IOUNIT,220) NAMR(I),(V(J),J=1,L)  
            ELSE
               WRITE (IOUNIT,230) I,(V(J),J=1,L)
            ENDIF   
         ELSE   
            IF (LBLROW) THEN
               WRITE (IOUNIT,320) NAMR(I),(V(J),J=1,L)  
            ELSE
               WRITE (IOUNIT,330) I,(V(J),J=1,L)
            ENDIF   
         ENDIF  
      ENDIF 
C   
  220 FORMAT (' ',A16,1P10D11.3)
  230 FORMAT (' ',I6,10X,1P10D11.3) 
  320 FORMAT (' ',A16,1P,6(2X,D14.7,3X))
  330 FORMAT (' ',I6,8X,1P6D19.7)   
C   
  50  CONTINUE  
      M   = M-L 
      IC  = IC+L
      IRR = IRR+L*IA
      IF (M.GT.0) GO TO 10  
      RETURN
      END   
