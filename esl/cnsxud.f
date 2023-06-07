      SUBROUTINE CNSXUD(XSTAR,N,VSM,BXJ,NX,EM,Q,IP,IERR)
C
C    COLORED NOISE SMOOTHING ALGORITHM ESTIMATE UPDATE, SIMPLE
C    PROCESS NOISE COMPONENT ADJUSTMENT, FROM TIME T BACK TO
C    TIME T-DT
C
C *** INPUTS:
C
C   XSTAR(N)      SMOOTH ESTIMATE AT START OF SMOOTH CYCLE
C
C   N             STATE VECTOR DIMENSION
C
C   VSM(N)        BIERMAN MODIFIED R-T-S SMOOTHER GAIN VECTOR
C
C   VSM(N+1)      FILTER STATE VECTOR - SMOOTHER GAIN DOT PRODUCT
C                 THAT WAS COMPUTED AND STORED DURING THE FILTER
C                 OPERATION
C
C   BXJ(NX)       PROCESS NOISE MODEL MULTIPLIER MATRIX; I.E.
C                 X    = PHI  X  + B   W (1)  + ...+B    W (NP)
C                  J+1      X  J    X/1 J            X/NP J
C                 WHERE B   = BXJ AND W  ARE THE PROCESS NOISE
C                        X/J           J
C                 COMPONENTS
C
C   NX            DIMENSION OF NON-TRIVIAL PORTION OF B; IF NX=0
C                 SYSTEM DEGENERATE TO PREVIOUS SMOOTH COLORED
C                 NOISE UPDATE CASE
C
C   EM            EXPONENTIAL MULTIPLIERS
C
C   Q             PROCESS NOISE VARIANCE
C
C   IP            LOCATIONS OF PROCESS NOISE COMPONENT IN THE
C                 STATE VECTOR (IP .GT. NX AND IP .GE. 1)
C
C   F(N)          SCRATCH VECTOR
C
C
C *** OUTPUT:
C
C   XSTAR(N)      SMOOTH UPDATE ADJUSTED FOR THIS PROCESS NOISE
C
C   IERR          ERROR FLAG, IF INCORRECT INPUT (EM AND Q BOTH 0)
C                 THEN IERR = IP ELSE IERR = 0. IF IERR .NE. 0
C                 XSTAR RETURN UNCHANGED.
C
C
C   REFERENCE:
C                 MODIFICATION OF THE BIERMAN COLORED NOISE
C                 SMOOTHER TO ACCOUNT FOR A NON-TRIVIAL PROCESS
C                 NOISE COVARIANCE.
C                 BIERMAN ASSOCIATES MEMO. 860429
C
C  COGNIZANT PERSONS:
C
C         DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C         FACTORIZED ESTIMATION APPLICATIONS
C         JUNE 2, 1986
C
C *****************************************************************
C *                                                               *
C *         SUBROUTINE CNSXUD IS A PART OF THE FEA. INC.          *
C *              ESTIMATION SUBROUTINE LIBRARY                    *
C *                                                               *
C *                           (ESL)                               *
C *                                                               *
C *                     COPYRIGHT 1986                            *
C *          FACTORIZED ESTIMATION APLLICATIONS INC.              *
C *   This library is licensed under the Creative Commons         *
C *   Attribution 3.0 Unported License.                           *
C *   To view a copy of this license,                             *
C *   visit http://creativecommons.org/licenses/by/3.0/           *
C *   or send a letter to Creative Commons, 171 Second Street,    *
C *   Suite 300, San Francisco, California, 94105, USA.           *
C *                                                               *
C *****************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION SUM,VXSTAR
      DIMENSION XSTAR(N),VSM(*),BXJ(*)
C
C      DATA EPS/1.0E-37/ ! suitable for single precision
      DATA EPS/1.0D-307/
      IERR = 0
C
C   SAVE VSM(IP) AND ZERO COMPONENT
C
      VSMIP = VSM(IP)
      VSM(IP) = 0.0

      VXSTAR = 0.0
      DO 10 K=1,N
         VXSTAR = VXSTAR + VSM(K)*XSTAR(K)
   10 CONTINUE
      SUM = 0.0
      IF (NX .GT. 0) THEN
         DO 20 K=1,NX
            SUM = SUM+BXJ(K)*VSM(K)
   20    CONTINUE
         SUM = Q*SUM
      ENDIF
      THETA = EM+SUM
      D = EM*THETA+Q*VSMIP
      IF (D .LT. EPS) THEN
         IF (Q .EQ. 0.0 .AND. EM .NE. 0.0 ) THEN
            C1 = 1.0/EM
            RLAMDA = 0.0
         ELSE
            IF (EM .EQ. 0.0) THEN
               IERR = IP
               GO TO 500
            ENDIF
         ENDIF
      ELSE
         RLAMDA = Q/D
         C1 = THETA/D
      ENDIF
      DELTA = VSM(N+1)-VXSTAR
      IF (NX .GT. 0) THEN
         SUM = (VSMIP*XSTAR(IP)-EM*DELTA)*RLAMDA
         DO 30 K=1,NX
            XSTAR(K) = XSTAR(K)-SUM*BXJ(K)
   30    CONTINUE
      ENDIF
      XSTAR(IP) = C1*XSTAR(IP)+RLAMDA*DELTA
  500 CONTINUE
      VSM(IP) = VSMIP
      RETURN
      END
