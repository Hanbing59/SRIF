      SUBROUTINE P2P (P,LIN,NAMIN,LOUT,NAMOUT)
C
C        SUBROUTINE TO REARRANGE THE ROWS AND COLUMNS OF MATRIX
C        P(LIN,LIN) IN NAMIN ORDER AND PUT THE RESULT IN
C        P(LOUT,LOUT) IN NAMOUT ORDER. ZERO COLUMNS AND ROWS ARE
C        ASSOCIATED WITH OUTPUT DEFINED NAMES THAT ARE NOT
C        CONTAINED IN NAMIN.
C
C *** INPUTS:
C
C        P(LIN(LIN+1)/2)  INPUT MATRIX, VECTOR STORED.
C
C
C        LIN              NO. OF PARAMETER NAMES ASSOCIATED WITH THE
C                         INPUT P
C
C        NAMIN(L)         PARAMETER NAMES ASSOCIATED WITH C ON INPUT.
C                         L = MAX(LIN,LOUT).
C                         (ONLY THE FIRST LIN ENTRIES APPLY TO THE
C                         INPUT C) (NOTE: DIM(NAMIN) .GE. DIM(NAMOUT))
C                         UP TO 16 CHARACTER NAMES ARE ALLOWED.
C                         NAMIN IS DESTROYED DURING PROCESSING.
C
C        LOUT             NO. OF PARAMETER NAMES ASSOCIATED WITH THE
C                         OUTPUT P
C
C        NAMOUT(LOUT)     PARAMETER NAMES ASSOCIATED WITH THE OUTPUT P
C                         (ALSO UP TO 16 CHARACTER NAMES)
C
C *** OUTPUT:
C
C        P(LOUT*(LOUT+1/2)OUTPUT MATRIX (OVERWRITES INPUT MATRIX)
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23, 1981       (C2C)
C
C           REVISED 9 JULY 1982       VERSION 0.01.00
C                                     VERSION 1.00.00
C           REVISED 20 MAY 1985       VERSION 2.00.00
C           REVISED  2 FEB 1986       VERSION 3.00.00 END OF 2 DIM COV
C                                                     NAME CHANGED
C
C       REQUIRED SUBROUTINES:
C
C       MAX0   <FORLIB> INTEGER TO INTEGER
C       LISPOS <ESL>    LOCATE NAME IN A LIST
C       SWPCOV <ESL>    SWAP TWO COLS OF A VECTOR STORED
C                       COVARIANCE MATRIX
C
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE P2P    IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *                   COPYRIGHT 1982 - 1987                        *
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

      DIMENSION  P(*)
      CHARACTER*(*) NAMIN(*), NAMOUT(LOUT)
C
      N1 = LIN
C
C     STEP #1 REMOVE NAMIN COLUMNS NOT APPEARING IN NAMOUT
C
      DO 20 J = LIN,1,-1
         CALL LISPOS(NAMIN(J),NAMOUT,LOUT,IPOS)
         IF (IPOS .EQ. 0) THEN
C                                 NAMIN(J) IS NOT IN NAMOUT
            NTR = J*(J+1)/2
            DO 15 I = NTR-IPOS+1,NTR
               P(I) = 0.0
 15         CONTINUE
            NTR = J
            DO 16 I = J,LIN
               P(NTR) = 0.0
               NTR    = NTR+I
 16         CONTINUE
            NAMIN(J) = ' '
         ENDIF
 20   CONTINUE
      DO 22 J = LIN+1, LOUT
         NAMIN(J) = ' '
 22   CONTINUE
C
C     STEP #2 ADD COLUMNS OF ZEROS FOR EACH NAMOUT NOT IN NAMIN
C
      DO 30 J = 1, LOUT
         CALL LISPOS(NAMOUT(J),NAMIN,LIN,IPOS)
         IF (IPOS .EQ. 0) THEN
C                                            FIND A BLANK COLUMN
            CALL LISPOS(' ',NAMIN,LIN,IPOS)
            IF (IPOS .EQ. 0) THEN
                N1   = N1+1
                IPOS = N1
            ENDIF
            NAMIN(IPOS) = NAMOUT(J)
            NTR = IPOS*(IPOS+1)/2
            DO 25 I = NTR-IPOS+1,NTR
               P(I) = 0.0
 25         CONTINUE
         ENDIF
 30   CONTINUE
C
C     STEP #3 NOW P (AND NAMIN) HAS ALL THE NAMOUT STATES, BUT
C             IN THE WRONG ORDER
 
      DO 40 I = 1, LOUT
         CALL LISPOS(NAMOUT(I),NAMIN,LOUT,IPOS)
         CALL SWPCOV(P,LOUT,NAMIN,I,IPOS)
40    CONTINUE
      RETURN
      END
