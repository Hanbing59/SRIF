      SUBROUTINE UDOPT(U,N,R,V,G,ALPHA)
C
C
C     TO COMPUTE OPTIMAL U-D COVARIANCE FACTOR MEASUREMENT UDATE
C     ASSUMING THAT V = (H*U)**T AND G = D*V HAVE BEEN PRECOMPUTED
C
C
C *** INPUTS
C
C     U(N*(N+1)/2)     INPUT U-D COVARIANCE FACTORS, VECTOR STORED
C
C     N                STATE DIMENSION
C
C     R                MEASUREMENT VARIANCE
C
C     V(N)             (H*U)**T, PRECOMPUTED
C
C     G(N)             SCALED V, G = D*V PRECOMPUTED
C
C *** OUTPUTS
C
C     U(N*(N+1)/2)     OUTPUT U-D COVARIANCE FACTORS, OVERWRITES
C                      INPUT
C
C     G(N)             NORMALIZED KALMAN GAIN, U*D*(H*U)**T,
C                      OVERWRITES INPUT
C
C     ALPHA            INNOVATION VARIANCE
C
C     COGNIZANT PERSONS:
C
C                      DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C                      FACTORIZED ESTIMATION APPLICATIONS
C                      JUNE 27, 1981
C                      REVISED 25 SEPTEMBER 1981.
C                      REVISED 23 DECEMBER  1981.
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE UDOPT  IS A PART OF THE FEA. INC             *
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

      DOUBLE PRECISION SUM,GAMMA,BETA,RATIO,EPS
      DIMENSION U(*),V(N),G(N)
C
C      DATA  EPS/0.6D-37/ ! suitable for single precision
      DATA EPS/1.0D-307/

C
C     EPS -- MACHINE EPSILON, SUCH THAT RECIPROCAL EXISTS.
C     EPS =  0.6D-37 IS ACCEPTABLE ON VAX TYPE MACHINES
C
C     UDMEAS SUBROUTINE CODE
C
      SUM   = R+V(1)*G(1)
      GAMMA = 0.0
      IF (SUM .GT. EPS) GAMMA = 1.0/SUM
      IF (V(1).EQ. 0.0) GO TO 10
         RATIO = R*GAMMA
         U(1)  = U(1) * RATIO
  10  CONTINUE
C
      IF ( N .EQ. 1) GO TO 60
      KJ      = 2
      DO 50 J=2,N
         BETA = SUM
         TEMP = G(J)
         P    = V(J)
C                 P=V(J)*GAMMA   CORRESPONDS TO THE MORE
C                                EFFICIENT ALGORITHM MECHANIZATION
C
         JM1  = J-1
C
         DO 20 K=1,JM1
            S     = U(KJ)
            TEMPA = G(K)*GAMMA
            U(KJ) = S-P*TEMPA
C
C           U(KJ)=S-P*G(K)       THIS ARRANGEMENT SAVES
C                                N*(N-2) MULTIPLIES, BUT IN
C                                A CERTAIN SINGULAR SITUATION LEADS
C                                TO AN OVERFLOW CONDITION
            G(K)=G(K)+TEMP*S
            KJ=KJ+1
 20      CONTINUE
         SUM=SUM+TEMP*V(J)
         IF (SUM  .LE. EPS ) GO TO 30
         GAMMA=1.0/SUM
         RATIO = BETA * GAMMA
         U(KJ) = U(KJ)* RATIO
  30     IF ( U(KJ) .GT. 0.0 ) GO TO 50
C
C        40 LOOP ZEROS COLUMN IN SINGULAR CASES
C
         L = KJ-J+1
         DO 40 K=L,KJ
            U(K) = 0.0
  40     CONTINUE
C
  50    KJ=KJ+1
  60  CONTINUE
      ALPHA=SUM
      RETURN
      END
