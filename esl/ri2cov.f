      SUBROUTINE RI2COV(RINV,N,SIG,COVOUT,KROW,KCOL)
C
C     TO COMPUTE THE STANDARD DEVIATIONS AND, IF DESIRED, THE
C     COVARIANCE MATRIX OF A VECTOR STORED UPPER TRIANGULAR SQUARE
C     ROOT COVARIANCE MATRIX. THE OUTPUT COVARIANCE MATRIX
C     IS ALSO VECTOR STORED.
C
C *** INPUTS:
C
C     RINV(N*(N+1)/2)   VECTOR STORED UPPER TRIANGULAR
C                       COVARIANCE SQUARE ROOT.
C                       (RINV = RINVERSE IS THE INVERSE OF
C                       THE SRIF MATRIX)
C
C     N                 DIMENSION OF THE RINV MATRIX, N .GE. 1
C
C
C     KROW .GT. 0       COMPUTES THE COVARIANCE AND SIGMAS
C
C          .LT. 0       COMPUTES ONLY THE SIGMAS OF THE FIRST
C                       KROW VARIABLES OF THE RINV MATRIX
C
C          .EQ. 0       NO COVARIANCE, BUT ALL SIGMAS
C                       (E.G. USE ALL N ROWS OF RINV)
C
C     KCOL              NUMBER OF COLUMNS OF COVOUT THAT ARE
C                       COMPUTED. IF KCOL .LE. 0 THEN KCOL=KROW.
C                       IF KROW .LE. 0 THIS INPUT IS IGNORED.
C
C *** OUTPUT
C
C     COVOUT(N*(N+1)/2) VECTOR STORED COVARIANCE
C                       (IF THE COMPILER PERMITS COVOUT MAY
C                        OVERWRITE RINV)
C
C     SIG(N)            OUTPUT VECTOR OF STANDARD DEVIATIONS
C
C     SUBROUTINES NEEDED:
C
C                 IABS       <FORLIB>
C                 DSQRT      <FORLIB>
C
C     COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              SEPTEMBER 18, 1976
C              VERSION 0.01.00
C
C              REVISED 2 AUGUST 1982          VERSION 0.02.00
C                                             VERSION 1.00.00
C              REVISED 20      MAY 1985       VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE RI2COV IS A PART OF THE FEA. INC             *
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

      DOUBLE PRECISION SUM, ZERO
      DIMENSION RINV(*), SIG(N), COVOUT(*)
C
      DATA ZERO/0.0D0/
C
      LIM  = N
      IF (KROW .NE. 0) LIM  = IABS(KROW)
C
C     COMPUTE SIGMAS
C
      IKS = 0
      DO 20 J = 1, LIM
         IKS = IKS+J
         SUM = ZERO
         IK  = IKS
         DO 10 K=J, N
            SUM = SUM+RINV(IK)**2
            IK = IK+K
  10     CONTINUE
         SIG(J) = DSQRT(SUM)
  20  CONTINUE
      IF (KROW .LE. 0) RETURN
      KKOL = KCOL
      IF (KKOL .LE. 0) KKOL = KROW
      JJ = 0
      IF (N .EQ. 1) GO TO 60
C
C *** COMPUTE COVARIANCE
C
         NM1 = LIM
         IF (KROW .EQ. N) NM1 = N-1
         DO 50 J=1, NM1
            JJ = JJ+J
            COVOUT(JJ) = SIG(J)**2
            IJS = JJ+J
            JP1 = J+1
            DO 40 I = JP1, KKOL
               IK  = IJS
               IMJ = I-J
               SUM = ZERO
               DO 30 K = I, N
                  IJK = IK + IMJ
                  SUM = SUM + RINV(IK)*RINV(IJK)
                  IK = IK+K
  30           CONTINUE
               COVOUT(IJS) = SUM
               IJS = IJS+I
  40        CONTINUE
  50     CONTINUE
  60  JJPN = JJ+N
      IF (KROW .EQ. N) COVOUT(JJPN) = SIG(N)**2
      RETURN
      END
