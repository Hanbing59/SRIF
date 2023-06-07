      SUBROUTINE RINCON (RIN,ROUT,N,CNB)
C
C     TO COMPUTE THE INVERSE OF THE UPPER TRIANGULAR VECTOR STORED
C     INPUT MATRIX RIN AND STORE THE RESULT IN ROUT. (RIN=ROUT IS
C     PERMITTED) AND TO COMPUTE A CONDITION NUMBER ESTIMATE.
C     CNB=FROB.NORM(R)*FROB.NORM(R**-1).
C     THE FROBENIUS NORM IS THE SQUARE ROOT OF THE SUM OF SQUARES
C     OF THE ELEMENTS. THIS CONDITION NUMBER BOUND IS USED AS
C     AN UPPER BOUND AND IT ACTS AS A LOWER BOUND ON THE ACTUAL
C     CONDITION NUMBER OF THE PROBLEM. (SEE THE BOOK 'SOLVING LEAST
C     SQUARES', BY LAWSON AND HANSON)
C
C     CNB/N .LE. CONDITION NUMBER .LE. N*CNB
C
C     WHEN RIN IS SINGULAR, RINCON COMPUTES AN INVERSE
C     CORRESPONDING TO THE THE SYSTEM HAVING ZERO ROWS
C     AND COLUMNS ASSOCIATED WITH THE SINGULAR DIAGONAL
C     ENTRIES.
C
C *** INPUTS:
C
C     RIN(N*(N+1)/2)  INPUT VECTOR STORED UPPER TRIANGULAR MATRIX
C
C     N               DIMENSION OF R MATRICES, N .GE. 1
C
C     CNB             IF CNB IS .LT. 0, THE CONDITION NUMBER BOUND
C                     CALCULATION IS SKIPPED, AND CNB IS LEFT
C                     UNCHANGED.
C
C
C *** OUTPUT:
C
C     ROUT(N*(N+1)/2) OUTPUT VECTOR STORED UPPER TRIANGULAR MATRIX
C                     INVERSE IF THE COMPILER PERMITS ROUT MAY
C                     OVERWRITE RIN (COMPILER PERMITTING)
C
C     CNB             IF CNB WAS NON-NEGATIVE ON INPUT THEN
C                           CNB/N .LE. C .LE. CNB
C                        WHERE C IS THE CONDITION NUMBER OF RIN.
C                     ELSE
C                        UNCHANGED
C
C  REQUIRED SUBROUTINES:
C
C   DABS   <FORLIB>                  DSQRT    <FORLIB>
C
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              20 OCTOBER 1981
C
C              REVISED 13 APRIL 1982
C              REVISED 17 JUNE  1982
C              VERSION 0.01.00
C
C              REVISED 2 AUGUST 1982          VERSION 0.01.01
C                                             VERSION 1.00.00
C              REVISED 20      MAY 1985       VERSION 2.00.00
C              REVISED  3      DEC 1986       VERSION 2.01.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE RINCON IS A PART OF THE FEA. INC             *
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

      DOUBLE PRECISION RNM,DINV,SUM,EPS
      DIMENSION RIN(*), ROUT(*)
C
C     EPS MACHINE DEPENDENT, USED TO TEST FOR ZERO DIVIDES
C     THIS EPS IS ADEQUATE FOR VAX TYPE APPLICATIONS
C
C      DATA EPS/0.6D-37/ !suitable for single precision
      DATA EPS/1.0D-307/

C
      Z=0.D0
      ONE=1.D0
      NTOT=N*(N+1)/2
C
C     CHECK TO SEE IF THE RHS IS INCLUDED
C
C
      IF (CNB .GE. 0.0) THEN
         NBAR = NTOT
         IF (RIN(NTOT) .EQ. -1.0) NBAR = NTOT-N
         TEMP = ABS(RIN(1))
         DO 5 J = 2,  NBAR
             TEMP = MAX(TEMP,ABS(RIN(J)))
  5      CONTINUE
C                                 IF RIN IS THE ZERO MATRIX
         IF (TEMP .EQ. 0.0) RETURN
         RNM = 0.0
         DO 10 J=1,NBAR
            RNM = RNM+(RIN(J)/TEMP)**2
  10     CONTINUE
         CNB = TEMP*SQRT(RNM)
      ENDIF
      SUM     = RIN(1)
      ROUT(1) = Z
      IF (DABS(SUM) .GT. EPS) ROUT(1) = ONE/SUM
      IF ( N .EQ. 1 ) GO TO 50
      JJ=1
      DO 40 J=2,N
         JJOLD=JJ
         JJ=JJ+J
         SUM = RIN(JJ)
         DINV = Z
         IF (DABS(SUM) .GT. EPS) DINV = ONE/SUM
         ROUT(JJ)=DINV
         II=0
         IK=1
         JM1=J-1
         DO 40 I=1,JM1
            II  =II+I
            IK  =II
            SUM = Z
            DO 30 K=I,JM1
               SUM=SUM+ROUT(IK)*RIN(JJOLD+K)
               IK=IK+K
  30        CONTINUE
          ROUT(JJOLD+I)=-SUM*DINV
  40     CONTINUE
  50  CONTINUE
      IF (CNB .GE. 0.0) THEN
         TEMP = ABS(ROUT(1))
         DO 55 J = 2,  NBAR
             TEMP = MAX(TEMP,ABS(ROUT(J)))
  55     CONTINUE
         RNM = 0.0
         DO 60 J=1, NBAR
            RNM = RNM+(ROUT(J)/TEMP)**2
  60     CONTINUE
         CNB = CNB*TEMP*SQRT(RNM)
      ENDIF
      RETURN
      END
