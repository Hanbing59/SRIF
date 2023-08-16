!! Purpose    : COMPUTE THE MEASUREMENT PREDICTION (INNOVATIONS) VARIANCE,
!!              ALPHA := H*[R**(-1)] * [R**(-T)] * (H**T) +1.0
!!              IF TST.GE.ZERO, COMPUTE
!!                         DZ := Z(MEAS) - H**XEST
!!                             = H(N+1) - H * [R**(-1)]*Z
!!              AND CARRY OUT THE CONFIDENCE TEST:
!!                 IF(DZ**2 .LE. TST*ALPHA) THEN
!!                    MEASUREMENT IS ACCEPTABLE, AND ALPHA POSITIVE IS
!!                    OUTPUT
!!                 OTHERWISE
!!                    ALPHA NEGATIVE IS OUTPUT.
!!
!! Parameters :
!!               R --- SQUARE-ROOT INFORMATION ARRAY.
!!                     IF TST.GE.0, THEN COLUMN N+1 IS ASSUMED TO HOLD THE
!!                     NORMALIZED ESTIMATE. R IS ASSUMED TO BE NON-SINGULAR.
!!               N --- STATE VECTOR DIMENSION
!!               H --- NORMALIZED MEASUREMENT COEFFICIENTS.
!!                     IF TST.GE.0, THEN H(N+1)=Z(MEAS) AND NBAR = N+1
!!                     ELSE ONLY THE FIRST N ENTRIES OF H ARE USED.
!!             TST --- CONFIDENCE TEST MULTIPLIER, VIZ.
!!                     TST = 9.0 CORRESPONDS TO A 3 SIGMA CONFIDENCE TEST.
!!              DZ --- IF TST.GE.0, THEN ERROR IS THE MEASUREMENT PREDICTION
!!                     RESIDUAL, Z(MEAS) - H*XEST
!!           ALPHA --- PREDICTED MEASUREMENT RESIDUAL VARIANCE:
!!                     H * P * (H**T) +1.0, WHERE P IS THE ESTIMATE ERROR
!!                     COVARIANCE: [R**(-1)] * [R**(-T)].
!!                     WHEN TST.LT.0 OR WHEN DATA REJECTED, THEN ALPHA
!!                     NEGATIVE IS RETURNED.
!!               V --- V**T = H*R**(-1). IF COMPILER ALLOWS, V CAN OVERWRITE H.

SUBROUTINE SINTST(R,N,H,TST,DZ,ALPHA,V)

IMPLICIT none

real*8    R(*)
integer*4 N
real*8    H(*)
real*8    TST,DZ,ALPHA
real*8    V(N)
!local
integer*4 IJ,J,K,JM1,NBAR,NTR
real*8    TMP,SUM
real*8    EPS
DATA EPS/0.6D-37/

IJ = 1
V(1) = 0.0
IF(ABS(H(1))*SQRT(EPS).LE.ABS(R(1)) .AND. R(1).NE.0.0D0) THEN
    V(1) = H(1)/R(1)
    TMP = 1.0 + V(1)*V(1)
ELSE
    ! R NEAR SINGULAR -- ACCEPT DATA
    ALPHA=0.1/SQRT(EPS)
    DZ=0.0
    RETURN
ENDIF
IF(N.ne.1) then
    DO J = 2, N
        JM1 = J-1
        SUM = H(J)
        DO K = 1, JM1
            IJ = IJ+1
            SUM = SUM - V(K)*R(IJ)
        enddo
        IJ = IJ + 1
        V(J) = 0.0
        IF(DABS(SUM)*SQRT(EPS).LE.ABS(R(IJ)) .AND. R(IJ).NE.0.0D0) THEN
            V(J) = SUM/R(IJ)
            IF(TMP.LE.(0.1/EPS)) THEN
                TMP = TMP + V(J)*V(J)
            ELSE
                ALPHA=0.1/SQRT(EPS)
                DZ=0.0
                RETURN
            END IF
        ELSE
            ALPHA=0.1/SQRT(EPS)
            DZ=0.0
            RETURN
        ENDIF
    enddo
endif
ALPHA = -TMP
IF(TST.LT.0.0) RETURN
NBAR = N+1
NTR = N*(NBAR)/2
SUM = H(NBAR)
DO J=1, N
    K = NTR+J
    SUM = SUM - V(J)*R(K)
enddo
DZ = SUM
IF(DABS(SUM).LE.SQRT(TST)*SQRT(TMP)) ALPHA = -ALPHA

RETURN
END