      SUBROUTINE INVTST(U,H,N,R,DZ,TST,ALPHA)
C
C     INVTST= " INNOVATIONS TEST "TO COMPUTE THE MEASUREMENT
C     INNOVATIONS VARIANCE ALPHA=H*P*(H**T)+R , USING THE U-D
C     COVARIANCE FACTORS OF P. THE RESULT IS USED TO FLAG THE
C     PREDICTED RESIDUAL VIA ( DZ ) ** 2 .LE. TST*ALPHA
C     IMPLIES THAT THE MEASUREMENT IS ACCEPTABLE.
C
C *** INPUTS:
C
C     U(N*(N+1)/2)         UPPER TRIANGULAR MATRIX, WITH "D"
C                          ELEMENTS STORED ON THE DIAGONAL.
C                          U IS VECTOR STORED AND CORRESPONDS
C                          TO THE A PRIORI COVARIANCE
C
C     H(N)                 OBSERVATION COEFFICIENTS
C
C     N                    STATE VECTOR DIMENSION (N .GE. 1)
C
C     DZ                   PREDICTED RESIDUAL
C
C     TST                  USED TO TEST DZ VALIDITY, I.E.
C                          ABS(DZ) .LE. SQRT(TST*ALPHA).
C                          E.G. TST= 9. CORRESPONDS TO A
C                          3 SIGMA TEST.
C
C *** OUTPUTS:
C
C     ALPHA                " SIGNED " INNOVATIONS VARIANCE.
C                          IF THE PREDICTED RESIDUAL FAILS
C                          THE TEST THEN NEGATIVE ALPHA
C                          IS OUTPUT.
C
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23, 1981
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE INVTST IS A PART OF THE FEA. INC             *
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
C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION U(*),H(N)
      DOUBLE PRECISION SUM,BETA
C
      SUM=R+U(1)*H(1)*H(1)
      IF (N .EQ. 1) GO TO 30
         JJ=1
         DO 20 J=2,N
            BETA=H(J)
            JM1=J-1
            DO 10 K=1,JM1
               JJPK = JJ+K
               BETA=BETA+H(K)*U(JJPK)
 10         CONTINUE
            JJ=JJ+J
            SUM=SUM+U(JJ)*BETA*BETA
 20      CONTINUE
 30   ALPHA = SUM
      IF (DZ*DZ .GT. TST*ALPHA) ALPHA = -SUM
      RETURN
      END
