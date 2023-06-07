C----------
      SUBROUTINE C2C (C,IMAXC,LIN,NAMIN,LOUT,NAMOUT)
C
C        SUBROUTINE TO REARRANGE THE ROWS AND COLUMNS OF MATRIX
C        C(LIN,LIN) IN NAMIN ORDER AND PUT THE RESULT IN
C        C(LOUT,LOUT) IN NAMOUT ORDER. ZERO COLUMNS AND ROWS ARE
C        ASSOCIATED WITH OUTPUT DEFINED NAMES THAT ARE NOT
C        CONTAINED IN NAMIN.
C
C *** INPUTS:
C
C        C(LIN,LIN)    INPUT MATRIX
C
C        IMAXC         ROW DIMENSION OF C, IMAXC.GE.L=MAX(LIN,LOUT)
C
C        LIN           NO. OF PARAMETER NAMES ASSOCIATED WITH THE
C                      INPUT C
C
C        NAMIN(L)      PARAMETER NAMES ASSOCIATED WITH C ON INPUT.
C                      (ONLY THE FIRST LIN ENTRIES APPLY TO THE
C                      INPUT C) (NOTE: DIM(NAMIN) .GE. DIM(NAMOUT))
C                      UP TO 16 CHARACTER NAMES ARE ALLOWED.
C
C        LOUT          NO. OF PARAMETER NAMES ASSOCIATED WITH THE
C                      OUTPUT C
C
C        NAMOUT(LOUT)  PARAMETER NAMES ASSOCIATED WITH THE OUTPUT C
C                      (ALSO UP TO 16 CHARACTER NAMES)
C
C *** OUTPUT:
C
C        C(LOUT,LOUT)  OUTPUT MATRIX (SHARES STORAGE WITH
C                                     INPUT MATRIX)
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23, 1981
C
C           REVISED 9 JULY 1982       VERSION 0.01.00
C                                     VERSION 1.00.00
C           REVISED 20 MAY 1985       VERSION 2.00.00
C
C       REQUIRED SUBROUTINES:
C
C         MAX0  <FORLIB>  INTEGER TO INTEGER
C
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE C2C    IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *             COPYRIGHT 1982, 1983, 1984, 1985                   *
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
      DIMENSION  C(IMAXC,LIN)
      CHARACTER*(*) NAMIN(*), NAMOUT(LOUT)
      CHARACTER*17 NAM,BLANKS
      DATA BLANKS/'                '/
C
      ZERO=0.0D0
      L=MAX0(LIN,LOUT)
      IF (L.LE.LIN) GO TO 5
      NM=LIN+1
      DO 1 K=NM,L
C                                   BLANK REMAINING NAMIN LOCNS
        NAMIN(K)=BLANKS
   1  CONTINUE
   5  DO 90 J=1,LOUT
         DO 10 I=1,L
             IF (NAMIN(I).EQ.NAMOUT(J)) GO TO 30
  10     CONTINUE
         GO TO 90
  30     IF (I.EQ.J) GO TO 90
         DO 40 K=1,L
C                                   INTERCHANGE COLUMNS I AND J
            H=C(K,J)
            C(K,J)=C(K,I)
            C(K,I)=H
  40     CONTINUE
         DO 80 K=1,L
C                                   INTERCHANGE ROWS I AND J
            H=C(J,K)
            C(J,K)=C(I,K)
            C(I,K)=H
  80     CONTINUE
C                                   INTERCHANGE LABELS I AND J
         NAM = NAMIN(I)
         NAMIN(I)=NAMIN(J)
         NAMIN(J) = NAM
  90    CONTINUE
C
C         FIND NAMOUT NAMES NOT IN NAMIN AND SET CORRESPONDING
C         ROWS AND COLUMNS TO ZERO
C
        DO 120 J=1,LOUT
           DO 100 I=1,L
               IF (NAMIN(I).EQ.NAMOUT(J)) GO TO 120
 100     CONTINUE
         DO 110 K=1,LOUT
            C(J,K)=ZERO
            C(K,J)=ZERO
 110     CONTINUE
 120    CONTINUE
      RETURN
      END
