      SUBROUTINE UTINV(R,N,IERR)
C
C        TO INVERT AN UPPER TRIANGULAR VECTOR STORED MATRIX AND
C        STORE THE RESULT IN VECTOR FORM. THE ALGORITHM IS SO
C        ARRANGED THAT THE RESULT OVERWRITES THE INPUT.
C        IN ADDITION TO  SOLVE RX=Z, SET RIN(N*(N+1)/2+1)=Z(1),
C        ETC., AND SET RIN((N+1)*(N+2)/2)=-1. CALL THE SUBROUTINE
C        USING N+1 INSTEAD OF N. ON RETURN THE FIRST N ENTRIES OF
C        COLUMN N+1 WILL CONTAIN X.
C        WHEN A DIAGONAL IS ZERO (OR MACHINE EPSILON CLOSE TO
C        ZERO) THE CORRESPONDING ROW AND COLUMN ARE SET TO ZERO
C
C *** INPUTS:
C
C        R  (N*(N+1)/2)  INPUT VECTOR STORED UPPER TRIANGULAR
C                        MATRIX
C
C        N               MATRIX DIMENSION
C
C
C *** OUTPUT:
C
C
C        R   (N*(N+1)/2) OUTPUT VECTOR STORED UPPER TRIANGULAR
C                        MATRIX INVERSE. OVERWRITES INPUT R
C
C     IERR            DIAGNOSTIC FLAG
C                     TO ALERT USER THAT HIS INPUT (AND OUTPUT)
C                     IS SINGULAR.
C
C                         .EQ. NUMBER OF ZERO COLUMNS IN
C                              THE OUTPUT R
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 28, 1981
C
C              REVISED MAY 26, 1982       VERSION 0.01.00
C                                         VERSION 1.00.00
C
C              REVISED 20   MAY 1985      VERSION 2.00.00
C
C    *** REQUIRED SUBROUTINE:
C
C        DABS  <FORLIB>
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE UTINV  IS A PART OF THE FEA. INC             *
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

      DOUBLE PRECISION ZERO, DINV, ONE, SUM, EPS
      DIMENSION R(*)
C
C      DATA  EPS/0.6D-37/! suitable for single precision
       DATA  ONE/1.0D0/, ZERO/0.0D0/
       DATA EPS/1.0D-307/

C
C     EPS -- MACHINE EPSILON, SUCH THAT RECIPROCAL EXISTS.
C     EPS =  0.6D-37 IS ACCEPTABLE ON VAX TYPE MACHINES
C
      IERR = 0
      SUM  = R(1)
      IF (DABS(SUM) .GT. EPS) GO TO 10
         IERR = IERR + 1
         DINV = ZERO
         GO TO 20
  10  DINV = ONE/SUM
  20  R(1) = DINV
      JJ = 1
      DO 60 J=2,N
         JJOLD = JJ
         JJ    = JJ+J
         SUM   = R(JJ)
         IF (DABS(SUM) .GT. EPS) GO TO 30
            IERR = IERR+1
            DINV = ZERO
            GO TO 40
  30     DINV = ONE/SUM
  40     R(JJ) = DINV
         II = 0
         IK = 1
         JM1 = J-1
         DO 60 I=1, JM1
            II = II+I
            IK = II
            SUM= ZERO
            DO 50 K=I, JM1
               SUM   = SUM + R(IK)*R(JJOLD+K)
               IK    = IK+K
  50        CONTINUE
            R(JJOLD+I) = -SUM*DINV
  60  CONTINUE
      RETURN
      END
