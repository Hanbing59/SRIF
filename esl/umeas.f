      SUBROUTINE UMEAS(U,N,R,A,F,G,ALPHA,TST)
C
C       COMPUTES ESTIMATE AND U-D MEASUREMENT UPDATED
C       COVARIANCE FACTORS, U-D, WHERE  P=U*D*U**T
C
C       *** INPUTS ***
C
C       U(NU)  UPPER TRIANGULAR VECTOR STORED MATRIX, OF DIMENSION
C              NU=N*(N+1)/2,WITH D ELEMENTS STORED ON THE DIAGONAL.
C              U IS VECTOR STORED AND  CORRESPONDS  TO AN A PRIORI
C              COVARIANCE.  IF STATE ESTIMATES ARE  COMPUTED, THERE
C              IS AN (N+1)ST COLUMN OF U THAT CONTAINS X, AND IN
C              THAT CASE TST .GE. 0.0 AND  NU = (N+1)*(N+2)/2.
C
C       N      DIMENSION OF THE STATE ESTIMATE.
C
C       R      MEASUREMENT VARIANCE
C
C       A(N+1) VECTOR OF MEASUREMENT COEFFICIENTS, IF DATA
C              (VIZ. TST .GE.0.0) THEN A(N+1)=Z, THE RESIDUAL.
C
C
C       TST    IF TST .GE. 0.0 THEN ESTIMATES ARE
C              INCLUDED AND THEN COLUMN N+1  OF
C              THE U-MATRIX IS ASSUMED TO HOLD THE STATE ESTIMATE
C              XEST, AND  A(N+1)  IS  THE MEASUREMENT RESIDUAL.
C              TST IS USED TO CHECK THE VALIDITY OF THE PREDICTED
C              RESIDUAL, DLTA= Z-A**T*XEST.
C              EXAMPLE: TST=9 CORRESPONDS TO A 3-SIGMA
C              INNOVATIONS TEST; IF (DLTA**2 .LE. TST*ALPHA)
C              THEN ACCEPT MEASUREMENT
C              IF TST .LT. 0.0 THEN X AND Z NEED NOT BE INCLUDED.
C              ( U AND A THEN HAVE SMALLER DIMENSIONS,
C              N*(N+1)/2 AND  N, RESPECTIVELY.)
C
C       *** OUTPUTS ***
C
C       U(NU)  UPDATED, VECTOR STORED FACTORS AND IF TST .GE. 0.0
C              THEN  COLUMN (N+1)CONTAINS THE UPDATED ESTIMATE AND
C              THE LAST ELEMENT, U((N+1)(N+2)/2), CONTAINS
C              (Z-A**T*X).
C
C       ALPHA  INNOVATIONS VARIANCE OF THE MEASUREMENT RESIDUAL
C              IF TEST FAILS THEN OUTPUT ALPHA IS SET
C              NEGATIVE AND NO UPDATE OF X NOR U IS MADE.
C
C       G(N)   VECTOR OF UNWEIGHTED KALMAN GAINS. THE KALMAN
C              GAIN, K, IS EQUAL TO G/ALPHA
C              IF THE TEST FAILS, G CONTAINS D*(A*U)**T
C
C       F(N+1) FIRST N COMPONENTS CONTAIN U**T*A AND
C              IF TST .GE. 0.0 THE (N+1)ST ELEMENT CONTAINS
C              (Z-A**T*X)/ALPHA . IF THE COMPILER ALLOWS IT,
C              ONE CAN HAVE F OVERWRITE A TO SAVE STORAGE
C
C
C     COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              31 AUGUST  1981
C
C              REVISED 25 SEPTEMBER 1981
C              REVISED 10 OCTOBER   1981
C              REVISED 25 DECEMBER  1981
C              REVISED 13 AUGUST    1982  VERSION 0.00.01
C                                         VERSION 1.00.00
C              REVISED 13 JUNE      1985  VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE UMEAS   IS A PART OF THE FEA. INC            *
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

      DIMENSION U(*), A(*), F(*), G(N)
      DOUBLE PRECISION SUM,BETA,GAMMA,RATIO,EPS
C
C
C      DATA EPS/0.6D-37/ ! suitable for single precision
      DATA EPS/1.0D-307/

C
C     EPS -- MACHINE EPSILON, SUCH THAT RECIPROCAL EXISTS.
C     EPS =  0.6D-37 IS ACCEPTABLE ON VAX TYPE MACHINES
C
      ZERO = 0.D0
      ONE  = 1.D0
      NP1  = N+1
      NTOT = N*NP1/2
C
C     FIND FIRST NONZERO ELEMENT OF H
C
      IHNZ = 0
      DO 5 J = 1, N
         IF (A(J) .NE. 0.0) THEN
            IHNZ = J
            GO TO 6
         ENDIF
  5   CONTINUE
  6   CONTINUE
      IF (IHNZ .EQ. 0) THEN
         DO 7 K = 1, N
            G(K) = 0.0
            F(K) = 0.0
7        CONTINUE
         ALPHA = R
         IF ((TST .GE. 0.0) .AND. (R .GT. EPS)) U(NTOT+NP1)=A(NP1)/R
         RETURN
      ENDIF
      IF (N .EQ.1) GO TO 30
       JJN   =  NTOT
       MX2IH = MAX(2,IHNZ)
       DO 20 J=N, MX2IH, -1
          JJ=JJN-J
          SUM=A(J)
          JM1=J-1
          DO 10 K=1,JM1
            SUM=SUM+U(JJ+K)*A(K)
   10     CONTINUE
          F(J)=SUM
          G(J)=SUM*U(JJN)
          JJN=JJ
   20  CONTINUE
   30 CONTINUE
      IF (IHNZ .EQ. 1) THEN
         F(1)=A(1)
         G(1)=U(1)*F(1)
      ELSE
         DO 35 K = 1, (IHNZ-1)
            F(K) = 0.0
            G(K) = 0.0
  35     CONTINUE
      ENDIF
C
C         F=U**T*A AND G=D*(U**T*A)
C
      IF (TST .LT. 0.0) GO TO 50
      BETA = R
      SUM=0.0
      DO 40 J=IHNZ,N
          BETA = BETA+F(J)*G(J)
          SUM=SUM+A(J)*U(NTOT+J)
   40 CONTINUE
C                                       Z=Z-A**T*X
      SUM         = A(NP1)-SUM
      U(NTOT+NP1) = SUM
      IF (SUM*SUM .LE. TST*BETA) GO TO 50
         ALPHA    = -BETA
         RETURN
   50 CONTINUE
      GAMMA  = ZERO
      IF (IHNZ .EQ. 1) THEN
C                                        SUM(1)
         SUM   = R+G(1)*F(1)
C                                        FOR R=0 CASE
         IF (SUM.GT.EPS) GAMMA=ONE/SUM
C                                        D(1)
         RATIO = R*GAMMA
         U(1)  = U(1)*RATIO
         IF ( N .EQ. 1) GO TO 110
         KJ    = 2
      ELSE
         KJ    = IHNZ *(IHNZ-1)/2+1
         SUM   = R
         IF (SUM .GT. EPS) GAMMA = ONE/SUM
      ENDIF
      DO 100 J=MX2IH, N
C                                        BETA=SUM(J-1)
      BETA     = SUM
      TEMP     = G(J)
C                                        P=-F(J)*(1/SUM(J-1))
C                                         EQN(21)
      P=-F(J)*GAMMA
C                       P = -F(J)        FOR SINGULAR PROBLEMS
C                                         WHERE OVERFLOW IN EQN(23)
C                                         IS A PROBLEM
      DO 70 K=1,(J-1)
         S=U(KJ)
C                                        EQN(22)
         U(KJ)=S+P*G(K)
C                                        EQN(23)
C                 TEMPA = G(K)*GAMMA
C                 U(KJ) = S+P*TEMPA
C                                        NOTE THAT THIS
C                                         ARRANGEMENT OF EQN(23)
C                                         REQUIRES AN ADDITIONAL
C                                         N*(N-2) EXTRA MULTIPLIES
         G(K) = G(K)+TEMP*S
         KJ   = KJ+1
  70  CONTINUE
C                                        FOR R=0 CASE
C                                        SUM(J)
      SUM=SUM+TEMP*F(J)
      IF (SUM .GT. EPS) THEN
C                                        GAMMA=1/SUM(J)
          GAMMA = ONE/SUM
C                                        D(J)     EQN(19)
          RATIO = BETA*GAMMA
          U(KJ) = U(KJ)*RATIO
      ENDIF
      IF ( U(KJ) .EQ. ZERO) THEN
C
C     90 LOOP ZEROS U COLUMN IN SINGULARITY CASES
C
             L       = KJ
             DO 90 K = 1,J
                U(L) = ZERO
                L    = KJ-K
  90         CONTINUE
      ENDIF
 100  KJ    = KJ+1
 110  CONTINUE
      ALPHA = SUM
C
C         EQN. NOS. REFER TO BIERMAN'S 1975 CDC PAPER, PP. 337-346.
C
      IF (TST .LE. 0.0) RETURN
      SUM    = U(NTOT+NP1)*GAMMA
      F(NP1) = SUM
      DO 120 J=1,N
        U(NTOT+J)=U(NTOT+J)+G(J)*SUM
 120  CONTINUE
      RETURN
      END
