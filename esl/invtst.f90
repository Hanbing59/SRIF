!! Purpose    : COMPUTE THE MEASUREMENT INNOVATIONS VARIANCE ALPHA=H*P*(H**T)+R, 
!!              USING THE U-D COVARIANCE FACTORS OF P. THE RESULT IS USED TO 
!!              FLAG THE PREDICTED RESIDUAL VIA ( DZ ) ** 2 .LE. TST*ALPHA
!!              IMPLIES THAT THE MEASUREMENT IS ACCEPTABLE.
!! Parameters : 
!!              U --- UPPER TRIANGULAR MATRIX, WITH "D" ELEMENTS STORED ON THE 
!!                    DIAGONAL. U IS VECTOR STORED AND CORRESPONDS TO THE 
!!                    A PRIORI COVARIANCE
!!              H --- OBSERVATION COEFFICIENTS
!!              N --- STATE VECTOR DIMENSION (N .GE. 1)
!!             DZ --- PREDICTED RESIDUAL
!!            TST --- USED TO TEST DZ VALIDITY, I.E. 
!!                    ABS(DZ) .LE. SQRT(TST*ALPHA). 
!!                    E.G. TST= 9 CORRESPONDS TO A 3 SIGMA TEST.
!!          ALPHA --- "SIGNED" INNOVATIONS VARIANCE.
!!                     IF THE PREDICTED RESIDUAL FAILS THE TEST, THEN NEGATIVE 
!!                     ALPHA IS OUTPUT.

SUBROUTINE INVTST(U,H,N,R,DZ,TST,ALPHA)

implicit none

real*8    U(*),H(N)
integer*4 N
real*8    R
real*8    DZ
real*8    TST
real*8    ALPHA
!local
real*8    SUM,BETA
integer*4 J,K,JJ,JM1,JJPK

SUM=R + U(1)*H(1)*H(1)
IF(N.ne.1) then
    JJ=1
    DO J=2,N
        BETA=H(J)
        JM1=J-1
        DO K=1,JM1
            JJPK = JJ+K
            BETA=BETA + H(K)*U(JJPK)
        enddo
        JJ=JJ+J
        SUM=SUM+U(JJ)*BETA*BETA
    enddo
endif
ALPHA = SUM
IF(DZ*DZ.GT.TST*ALPHA) ALPHA = -SUM

RETURN
END