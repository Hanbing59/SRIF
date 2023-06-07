      SUBROUTINE TTHH(R,RA,N,IEST)
C
C     TO COMBINE TWO SINGLE SUBSCRIPTED SRIF ARRAYS
C     USING HOUSEHOLDER ORTHOGONAL TRANSFORMATIONS
C
C *** INPUTS:
C
C     R(N*(N+1)/2)       INPUT VECTOR STORED UPPER
C                        TRIANGULAR MATRIX
C
C     RA(N*(N+1)/2)      THE SECOND INPUT VECTOR STORED UPPER
C                        TRIANGULAR MATRIX, DESTROYED
C                        BY THE COMPUTATION
C
C     N                  DIMENSION OF THE ESTIMATED PARAMETER VECTOR
C
C     IEST               LOGICAL, IEST .TRUE. INDICATES THAT
C                        THE DATA ARRAYS INCLUDE A RIGHT HAND
C                        SIDE, NORMALIZED ESTIMATE, AND THE
C                        ARRAYS HAVE DIMENSION NBAR*(NBAR+1)/2
C                        WHERE NBAR = N+1
C *** OUTPUT:
C
C
C     R(NBAR*(NBAR+1)/2) THE OUTPUT OVERWRITES THE INPUT "R"
C
C
C  REQUIRED SUBROUTINES:
C
C    ABS        <FORLIB>
C    DSQRT      <FORLIB>
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23, 1981
C              REVISED  30 APRIL 1982
C
C              REVISED 9  AUGUST 1982      VERSION 0.00.01
C                                          VERSION 1.00.00
C              REVISED 5   MAY   1984      VERSION 2.00.00
C              REVISED 20  JUNE  1985      VERSION 2.01.00
C
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE TTHH   IS A PART OF THE FEA. INC             *
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

      DIMENSION RA(*), R(*)
      DOUBLE PRECISION SUM, DELTA, ZERO
      LOGICAL IEST
C
      DATA ZERO /0.0D0/
C
      IF (IEST) THEN
         NBAR = N+1
      ELSE
         NBAR = N
      ENDIF
      IJS = 1
      KK  = 0
      DO 100 J=1, N
C                       J-TH STEP OF HOUSEHOLDER REDUCTION
         KK = KK+J
         TST = 0.0
         DO 10 I = IJS, KK
            IF (TST .LT. ABS(RA(I))) TST = ABS(RA(I))
  10     CONTINUE
C
         IF (TST .LE. 0.0) GO TO 90
C
C        IF COLUMN J OF RA IS ZERO SKIP THIS STEP OF THE REDUCTION
C
         IF (ABS(R(KK)) .GT. TST) TST = ABS(R(KK))
C
         SUM = (R(KK)/TST)**2
         DO 20 I = IJS, KK
            SUM = SUM + (RA(I)/TST)**2
  20     CONTINUE
         SUM = TST * DSQRT(SUM)
         IF (R(KK) .GT. ZERO) SUM = -SUM
         DO 40  I = IJS, KK
            RA(I) = RA(I)/SUM
  40     CONTINUE
C
         DELTA = R(KK)/SUM - 1.0
         R(KK) = SUM
         JJ    = KK+J
         JP1   = J+1
         IKS   = KK+1
C
C ***      J - TH HOUSEHOLDER TRANS. DEFINED
C          "80" LOOP APPLIES TRANSFORM. TO COLS. J+1 TO NBAR
C
         DO 80 K = JP1, NBAR
            IK   = IKS
            SUM  = DELTA*R(JJ)
            DO 50 I = IJS, KK
               SUM  = SUM + RA(IK)*RA(I)
               IK   = IK+1
  50        CONTINUE
            IF (SUM .NE. ZERO) THEN
                R(JJ) = R(JJ) + SUM
                SUM   = SUM/DELTA
                IK    = IKS
                DO 60 I = IJS, KK
                   RA(IK) = RA(IK) + SUM*RA(I)
                   IK = IK+1
  60            CONTINUE
             ENDIF
             JJ  = JJ+K
             IKS = IKS+K
  80     CONTINUE
  90     CONTINUE
         IJS = KK+1
 100  CONTINUE
      RETURN
      END
