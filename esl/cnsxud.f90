!! Purpose    : COLORED NOISE SMOOTHING ALGORITHM ESTIMATE UPDATE, SIMPLE
!!              PROCESS NOISE COMPONENT ADJUSTMENT, FROM TIME T BACK TO
!!              TIME T-DT
!! Parameters : 
!!              XSTAR --- SMOOTH ESTIMATE AT START OF SMOOTH CYCLE
!!                        For output, SMOOTH UPDATE ADJUSTED FOR THIS PROCESS NOISE
!!                  N --- STATE VECTOR DIMENSION
!!                VSM --- (  1:  N), BIERMAN MODIFIED R-T-S SMOOTHER GAIN VECTOR
!!                        (N+1:N+1), FILTER STATE VECTOR - SMOOTHER GAIN DOT PRODUCT
!!                                   THAT WAS COMPUTED AND STORED DURING THE FILTER
!!                                   OPERATION
!!
!!                BXJ --- PROCESS NOISE MODEL MULTIPLIER MATRIX
!!                 NX --- DIMENSION OF NON-TRIVIAL PORTION OF B
!!                        IF NX=0, SYSTEM DEGENERATE TO PREVIOUS SMOOTH COLORED
!!                        NOISE UPDATE CASE
!!                 EM --- EXPONENTIAL MULTIPLIERS
!!                  Q --- PROCESS NOISE VARIANCE
!!                 IP --- LOCATIONS OF PROCESS NOISE COMPONENT IN THE
!!                        STATE VECTOR (IP .GT. NX AND IP .GE. 1)
!!               IERR --- ERROR FLAG, IF INCORRECT INPUT (EM AND Q BOTH 0)
!!                        THEN IERR = IP ELSE IERR = 0. IF IERR .NE. 0
!!                        XSTAR RETURN UNCHANGED.
!!
!! Reference  : MODIFICATION OF THE BIERMAN COLORED NOISE SMOOTHER TO ACCOUNT 
!!              FOR A NON-TRIVIAL PROCESS NOISE COVARIANCE. BIERMAN ASSOCIATES 
!!              MEMO. 860429 (not publicly found)

SUBROUTINE CNSXUD(XSTAR,N,VSM,BXJ,NX,EM,Q,IP,IERR)

implicit none

real*8    XSTAR(N)
integer*4 N
real*8    VSM(*),BXJ(*)
integer*4 NX
real*8    EM,Q
integer*4 IP
integer*4 IERR
!local
integer*4 K
real*8    SUM,VXSTAR,VSMIP,THETA,D,C1,RLAMDA,DELTA
real*8    EPS
DATA EPS /1.0D-307/

IERR = 0
! SAVE VSM(IP) AND ZERO COMPONENT
VSMIP = VSM(IP)
VSM(IP) = 0.0

VXSTAR = 0.0
DO K=1,N
    VXSTAR = VXSTAR + VSM(K)*XSTAR(K)
enddo
SUM = 0.0
IF(NX.GT.0) THEN
    DO K=1,NX
        SUM = SUM + BXJ(K)*VSM(K)
    enddo
    SUM = Q*SUM
ENDIF
THETA = EM + SUM
D = EM*THETA + Q*VSMIP
IF(D.LT.EPS) THEN
    IF(Q.EQ.0.0 .AND. EM.NE.0.0) THEN
        C1 = 1.0/EM
        RLAMDA = 0.0
    ELSE
        IF(EM.EQ.0.0) THEN
            IERR = IP
            VSM(IP) = VSMIP
            RETURN
        ENDIF
    ENDIF
ELSE
    RLAMDA = Q/D
    C1 = THETA/D
ENDIF
DELTA = VSM(N+1) - VXSTAR
IF(NX.GT.0) THEN
    SUM = (VSMIP*XSTAR(IP)-EM*DELTA)*RLAMDA
    DO K=1,NX
        XSTAR(K) = XSTAR(K)-SUM*BXJ(K)
    enddo
ENDIF
XSTAR(IP) = C1*XSTAR(IP) + RLAMDA*DELTA
VSM(IP) = VSMIP

RETURN
END