      SUBROUTINE PERMUT(A,IMAXA,IRA,JCIN,NAMIN,JCOUT,NAMOUT)
C
C        SUBROUTINE TO REARRANGE PARAMETERS OF A(IRA,JCIN), NAMIN
C        ORDER TO A(IRA,JCOUT), NAMOUT ORDER.  ZERO COLUMNS
C        ARE INSERTED CORRESPONDING TO THE NEWLY DEFINED NAMES.
C
C *** INPUTS:
C
C        A(IRA,L)      INPUT RECTANGULAR MATRIX, L=MAX(JCIN,JCOUT)
C
C        IMAXA         ROW DIMENSION OF A, IMAXA.GE.IRA
C               
C        IRA           NUMBER OF ROWS OF A THAT ARE TO BE REARRANGED
C
C        JCIN          NUMBER OF PARAMETER NAMES ASSOCIATED
C                      WITH THE INPUT  A MATRIX
C
C        NAMIN(L)     PARAMETER NAMES ASSOCIATED WITH A ON INPUT
C                      (ONLY THE FIRST JCIN ENTRIES APPLY TO THE
C                       INPUT A) NAMIN IS DESTROYED BY PERMUT.
C                      (UP TO CHARACTER*16).
C                      WHERE L = MAX(JCIN,JCOUT)   
C
C *** OUTPUTS:
C
C
C        A(IRA,L)      OUTPUT RECTANGULAR MATRIX.
C
C        JCOUT         NUMBER OF PARAMETER NAMES ASSOCIATED WITH
C                      THE OUTPUT A MATRIX
C
C        NAMOUT(JCOUT) PARAMETER NAMES ASSOCIATED WITH THE OUTPUT A
C                      (UP TO CHARACTER*16)
C
C     REQUIRED SUBROUTINE:
C
C      MAX0  <FORLIB>  
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23,1981
C
C              REVISED 20 AUGUST 1982         VERSION 0.00.01
C                                             VERSION 1.00.00
C              REVISED 20 MAY 1985            VERSION 2.00.00
C              (16 CHARACTER VERSION)              
C              REVISED 8  JUNE 1986           VERSION 2.01.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE PERMUT IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *          COPYRIGHT 1982, 1983, 1984, 1985, 1986                *
C *          FACTORIZED ESTIMATION APPLICATIONS INC.               *
C *   This library is licensed under the Creative Commons          *
C *   Attribution 3.0 Unported License.                            *
C *   To view a copy of this license,                              *
C *   visit http://creativecommons.org/licenses/by/3.0/            *
C *   or send a letter to Creative Commons, 171 Second Street,     *
C *   Suite 300, San Francisco, California, 94105, USA.            *
C *                                                                *
C ******************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)

      CHARACTER*(*) NAMIN(*), NAMOUT(JCOUT)
      CHARACTER*17 NAMTMP,BLANKS
      DIMENSION A(IMAXA,*)
      DATA BLANKS/'               '/
C
      ZERO=0.0
      L=MAX0(JCIN,JCOUT)
C
C        LOCATIONS BEYOND JCIN NOT REQUIRED, SO WE ARE FREE TO 
C        TO BLANK THEM
C
      IF (JCOUT .GT. JCIN) THEN
         DO 10 J = (JCIN+1),JCOUT
            NAMIN(J) = BLANKS
  10     CONTINUE
      ENDIF
      DO 30 J = 1, JCIN
         NAMTMP = NAMIN(J)
         DO 20 K = 1, JCOUT
            IF (NAMTMP .EQ. NAMOUT(K)) GO TO 30
  20     CONTINUE
         NAMIN(J) = BLANKS
  30  CONTINUE
      DO 120 J=1,JCOUT
         DO 60 I=1,L
            IF (NAMIN(I).EQ.NAMOUT(J))  GO TO 70
  60     CONTINUE
         GO TO 120
  70     CONTINUE
         IF (I .EQ. J) GO TO 120
         IF (I .LT. J) THEN
C                        ---------- FIND A BLANK
            DO 80 I2=1, L
               IF (NAMIN(I2) .EQ. BLANKS) GO TO 90
  80        CONTINUE
              PRINT*,' PERMUT: FATAL ERROR COULD NOT FIND A BLANK'
              PRINT*,'  COLUMN TO TEMPORARILY STORE COLUMN ',J,' IN'
C --------------         BY CONSTRUCTION THIS CANNOT HAPPEN
              STOP
  90        CONTINUE
            DO 100 K = 1, IRA
C                                       COPY COL J TO COL I2 (BLANK)
               A(K,I2) = A(K,J)
C                                       COPY COL I TO COL J
               A(K,J)  = A(K,I)
 100        CONTINUE
            NAMIN(I2)  = NAMIN(J)
            NAMIN(J)   = NAMIN(I) 
            NAMIN(I)   = ' '
         ELSE
C                                       INTERCHANGE COLS I AND J
            DO 110 K=1,IRA
               W=A(K,J)
               A(K,J)=A(K,I)
               A(K,I)=W
 110        CONTINUE
            NAMTMP   = NAMIN(I)
            NAMIN(I) = NAMIN(J)
            NAMIN(J) = NAMTMP
         ENDIF
 120  CONTINUE
C                           REPEAT TO FILL NEW COLS W/ZEROS
      DO 150 J=1,JCOUT
         DO 130 I = 1, L
            IF (NAMIN(I) .EQ. NAMOUT(J)) GO TO 150
 130     CONTINUE
         DO 140 K=1,IRA
            A(K,J)=ZERO
 140     CONTINUE
 150  CONTINUE
      RETURN
      END 
