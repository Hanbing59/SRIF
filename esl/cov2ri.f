      SUBROUTINE COV2RI(CF,N)
C
C     TO CONSTRUCT THE UPPER TRIANGULAR CHOLESKY FACTOR
C     OF A POSITIVE SEMI-DEFINITE MATRIX. BOTH THE INPUT
C     COVARIANCE AND THE OUTPUT CHOLESKY FACTOR (SQUARE ROOT)
C     ARE VECTOR STORED. THE OUTPUT OVERWRITES THE INPUT.
C
C           COVARIANCE (INPUT) = CF*CF**T (CF IS OUTPUT)
C
C     IF THE INPUT COVARIANCE IS SINGULAR OR INDEFINITE THE
C     OUTPUT FACTOR HAS ZERO COLUMNS.
C
C *** INPUTS:
C
C     CF(N*(N+1)/2)        CONTAINS THE INPUT VECTOR STORED
C                          COVARIANCE MATRIX (ASSUMED POSITIVE
C                          DEFINITE)
C
C     N                    STATE DIMENSION
C
C *** OUTPUT:
C
C     CF(N*(N+1)/2)        THE UPPER TRIANGULAR SQUARE ROOT
C                          CHOLESKY FACTOR
C
C    REQUIRED SUBROUTINE:
C
C           DSQRT     <FORLIB>
C
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              FEBRURARY 20,  1978
C
C              REVISED 30 MARCH 1982
C
C                                               VERSION 1.00.00
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE COV2RI IS A PART OF THE FEA. INC             *
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
      DIMENSION CF(*)
      DOUBLE PRECISION ALPHA, EPS, S, ZERO
C
C     EPS IS, A MACHINE DEPENDENT LOWER BOUND FOR ALLOWABLE
C     DIVISION. 0.6E-37 IS SUITABLE FOR  VAX TYPE MACHINES
C
C      DATA EPS/0.6D-37/ ! suitable for single precision
      DATA EPS/1.0D-307/
      DATA ZERO/0.0D0/
C
      IF (N .EQ. 1) GO TO 60
         JJ = N*(N+1)/2
         JJN = JJ
         NP2 = N+2
         DO 50 L = 2, N
            J  = NP2-L
            S  = CF(JJ)
            IF (S .LT. ZERO) S = ZERO
            S      = DSQRT(S)
            CF(JJ) = S
            ALPHA  = ZERO
            IF (S .GT. EPS) ALPHA = 1.0/S
            KK  = 0
            JJN = JJ-J
            JM1 = J-1
            DO 40 K = 1, JM1
               JJNPK = JJN+K
               CF(JJNPK) = ALPHA*CF(JJNPK)
               S = CF(JJNPK)
               DO 30 I = 1, K
                  KKPI = KK + I
                  JJNPI= JJN + I
                  CF(KKPI) = CF(KKPI) - S*CF(JJNPI)
  30           CONTINUE
               KK = KK+K
  40        CONTINUE
            JJ = JJN
  50     CONTINUE
  60  CONTINUE
      S = CF(1)
      IF (S .LT. ZERO) S = ZERO
      CF(1) = DSQRT(S)
      RETURN
      END
