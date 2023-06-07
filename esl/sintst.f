      SUBROUTINE SINTST(R,N,H,TST,DZ,ALPHA,V)
C
C
C     TO COMPUTE THE MEASUREMENT PREDICTION (INNOVATIONS)
C     VARIANCE:
C
C         ALPHA := H*[R**(-1)] * [R**(-T)] * (H**T) +1.0
C
C     AND IF TST .GE. ZERO TO COMPUTE:
C
C         DZ := Z(MEAS) - H**XEST
C             = H(N+1) - H * [R**(-1)]*Z
C
C     AND TO CARRY OUT THE CONFIDENCE TEST:
C
C         IF ( DZ**2 .LE. TST*ALPHA ) THEN
C            MEASUREMENT IS ACCEPTABLE, AND ALPHA POSITIVE IS
C            OUTPUT, OTHERWISE ALPHA NEGATIVE IS OUTPUT.
C
C
C *** INPUTS:
C
C     R(NBAR*(NBAR+1)/2)    SQUARE-ROOT INFORMATION ARRAY,
C                           IF TST .GE. 0.0 THEN
C                              COLUMN N+1 IS ASSUMED TO
C                              HOLD THE NORMALIZED ESTIMATE.
C                           R IS ASSUMED TO BE NON-SINGULAR.
C
C
C     N                     STATE VECTOR DIMENSION
C
C     TST                   CONFIDENCE TEST MULTIPLIER, VIZ.
C                           TST = 9.0 CORRESPONDS TO A 3 SIGMA
C                           CONFIDENCE TEST.
C
C     H(NBAR)               NORMALIZED MEASUREMENT COEFFICIENTS,
C                           AND IF TST .GE. 0.0 THEN
C                                  H(N+1) = Z(MEAS)   AND
C                                  NBAR   = N+1
C                                ELSE
C                                  ONLY THE FIRST N ENTRIES OF H
C                                  ARE USED.
C
C *** OUTPUTS:
C
C     DZ                    IF TST .GE. 0.0 THEN
C                              ERROR IS THE MEASUREMENT PREDICTION
C                              RESIDUAL:
C                                        Z(MEAS) - H*XEST
C
C     ALPHA                 PREDICTED MEASUREMENT RESIDUAL VARIANCE:
C                              H * P * (H**T) +1.0,
C                           WHERE P IS THE ESTIMATE ERROR
C                           COVARIANCE:
C                              [R**(-1)] * [R**(-T)].
C                           WHEN TST .LT. 0.0 OR WHEN DATA REJECTED
C                                THEN ALPHA  NEGATIVE IS RETURNED.
C
C     V(N)                  V**T = H*R**(-1).
C                           IF COMPILER ALLOWS, V CAN OVERWRITE H.
C
C
C
C  **** REQUIRED SUBROUTINES:
C
C       ABS   <FORLIB>       DABS <FORLIB>
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              MARCH 19, 1982
C
C              REVISED 5 AUGUST 1982      VERSION 0.01.00
C                                         VERSION 1.00.00
C              REVISED 20   MAY 1985      VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE SINTST IS A PART OF THE FEA. INC             *
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

      DOUBLE PRECISION TMP, SUM
      DOUBLE PRECISION EPS
      DIMENSION R(*), V(N), H(*)
C
C     EPS MACHINE DEPENDENT, USED TO TEST FOR ZERO DIVIDES
C     THIS EPS IS ADEQUATE FOR VAX TYPE APPLICATIONS
C
      DATA EPS/0.6D-37/
C
      IJ = 1
      V(1) = 0.0
      IF (ABS(H(1))*SQRT(EPS).LE.ABS(R(1)).AND.R(1).NE.0.0D0) 
     $   THEN
            V(1) = H(1)/R(1)                            
            TMP = 1.0 + V(1)*V(1)
      ELSE
            ALPHA=0.1/SQRT(EPS)    ! R NEAR SINGULAR -- ACCEPT DATA
            DZ=0.0
            RETURN
C                               ! MACHINE EQUIV. OF INFINITY
      END IF
      IF (N .EQ. 1) GO TO 30
         DO 20 J = 2, N
            JM1 = J-1
            SUM = H(J)
            DO 10 K = 1, JM1
               IJ = IJ+1
               SUM = SUM - V(K)*R(IJ)
  10        CONTINUE
            IJ = IJ + 1
            V(J) = 0.0
            IF (DABS(SUM)*SQRT(EPS).LE. ABS(R(IJ)).AND.R(IJ).NE.0.0D0) 
     $         THEN
                  V(J) = SUM/R(IJ)
                  IF(TMP.LE.(0.1/EPS)) THEN
                     TMP = TMP + V(J)*V(J)
                  ELSE
                     ALPHA=0.1/SQRT(EPS)    ! R NEAR SINGULAR -- ACCEPT DATA
                     DZ=0.0
                     RETURN  
                  END IF
            ELSE
               ALPHA=0.1/SQRT(EPS)    ! R NEAR SINGULAR -- ACCEPT DATA
               DZ=0.0
               RETURN
            END IF
  20     CONTINUE
  30  ALPHA = -TMP                
      IF (TST .LT. 0.0) RETURN
      NBAR = N+1
      NTR = N*(NBAR)/2
      SUM = H(NBAR)
      DO 40 J=1, N
         K = NTR+J
         SUM = SUM - V(J)*R(K)
  40  CONTINUE
      DZ = SUM     
C
C       OLD CODE: IF (SUM*SUM .LE. TST*TMP)   .... CHANGED 8/18/88 - S.LICHTEN
C
      IF (DABS(SUM).LE. SQRT(TST)*SQRT(TMP)) ALPHA = -ALPHA
      RETURN
      END
