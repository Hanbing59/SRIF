      SUBROUTINE UDMVEC(U,N,H,MAXH,F,MAXF,UALPHA,M,V,MAXV,IEST,XEST,
     1 Z,RESTST)
C
C     U-D FACTOR MEASUREMENT UPDATE FOR VECTOR MEASUREMENT
C     PROCESSING. THE ALGORITHM IS AN APPLICATION OF MODIFIED
C     GIVENS (SQUARE-ROOT FREE) ROTATIONS TO THE MORF-KAILATH
C     PRE-ARRAY/POST-ARRAY MATRIX IDENTITY.
C
C  *** INPUT ***
C
C     U(NTR)            U-D COVARIANCE FACTORS; U HAS IMPLICIT UNIT
C                       DIAGONAL ENTRIES.
C
C     N                 FILTER DIMENSION; NTR=N*(N+1)/2.
C
C     H(M,N)            OBSERVATION MATRIX.
C
C     MAXH              ROW DIMENSION OF H MATRIX M .LE. MAXH
C
C     UALPHA(MTR)       U-D FACTORS OF MEASUREMENT NOISE COVARIANCE.
C
C     M                 NUMBER OF MEASUREMENTS TO BE SIMULTANEOUSLY
C                       PROCESSED; MTR = M*(M+1)/2.
C
C     IEST              (LOGICAL FLAG) WHEN
C
C       IEST = .TRUE.  ESTIMATES, XEST, AND DATA, Z, MUST BE INCLUDED.
C       IEST = .FALSE. ESTIMATES AND DATA ARE NOT USED.
C
C     XEST(N)           ESTIMATE STATE (USED ONLY IF IEST = .TRUE.)
C
C     Z(M)              DATA (Z=HX+NOISE); USED ONLY IF IEST = .TRUE.)
C                       WHEN IEST=.TRUE.Z IS MODIFED
C
C
C  *** OUTPUT ***
C
C     U(NTR)            MEASUREMENT UPDATED U-D COVARIANCE FACTORS.
C
C
C     UALPHA(MTR)       U-D FACTORS OF THE PREDICTION RESIDUAL
C                       (INNOVATIONS). THIS RESULT OVERWRITES THE
C                       INPUT MEASUREMENT NOISE COVARIANCE. THE
C                       ENTRIES CORRESPOND TO THE ORTHOGONALIZED
C                       RESIDUALS THAT ARE OUTPUT IN VECTOR Z.
C
C     V(N,M)            DALPHA NORMALIZED GAIN MATRIX, V=G*UALPHA
C                       WHERE G IS THE KALMAN FILTER GAIN AND UALPHA
C                       IS THE UNIT UPPER TRIANGULAR FACTOR OF THE
C                       INNOVATIONS COVARIANCE.
C
C  -------- WHEN IEST IS TRUE THEN:
C
C     Z :=              U(ALPHA)**-1(Z-H*XEST) WHERE
C                       HERE THE DIAGONALS OF U ARE IMPLICIT UNITY.
C                       THIS OUTPUT RESIDUAL IS THE GRAM-SCHMIDT
C                       ORTHOGONALIZATION OF THE PREDICTION RESIDUAL,
C                       AND OVERWRITES THE INPUT MEASUREMENT. THESE
C                       ARE THE RESIDUALS THAT WOULD RESULT IF THE
C                       MEASUREMENTS ARE SCALAR PROCESSED BUT IN
C                       REVERSE ORDER.
C
C     XEST :=           XEST(INPUT) + G*(Z(INPUT) - H*XEST(INPUT))
C           =           XEST(INPUT) + V*Z(OUTPUT)
C                       THIS RESULT OVERWRITES THE INPUT ESTIMATE.
C
C     RESTST            MAXIMUM LIKLIHOOD RESIDUAL TEST STATISTIC;
C                       USEFUL FOR DATA VALIDITY CHECKS.
C
C           =           (Z-H*XEST)**T * (P(ALPHA)**-1) * (Z-H*XEST)
C
C                       WHERE Z AND XEST HERE REPRESENT THE INPUTS.
C                       NOTE THAT THE U-D FACTORS AND ESTIMATES
C                       ARE UPDATED SO THAT IF THE USER WANTS TO
C                       USE RESTST TO AVOID BLUNDER POINTS (E.G.
C                       3-SIGMA TEST) HE WILL HAVE TO RESET U AND
C                       XEST IN THE PROGRAM THAT CALLS THIS UPDATE.
C  *** SCRATCH ***
C
C     F(M,N)            IF THE H MATRIX IS NOT NEEDED FURTHER, ONE
C                       CAN HAVE F AND H SHARE STORAGE.
C
C     MAXF              ROW DIMENSION OF F MATRIX, UNCHANGED.
C                       (IF REDUNDANT ARGUMENTS FOR H AND F ARE USED,
C                       F AND H SHARE STORAGE, THEN MAXF SHOULD BE
C                       EQUAL TO MAXH.)
C
C
C  **** SUBROUTINE REQUIRED:
C
C    PHIU               COMPUTES THE MATRIX PRODUCT PHI*U
C                       (IF THE OBSERVATION MATRIX, H, IS SPARSE
C                        ONE MAY PREFER TO USE THE SPARSE MULTIPLY
C                        SUBROUTINE,SFU).
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              AUGUST 22, 1983
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C              REVISED 19      MAY 1986         VERSION 2.01.00
C
C ******************************************************************
C *                                                                *
C *       SUBROUTINE UDMVEC IS A PART OF THE FEA. INC              *
C *              ESTIMATION  SUBROUTINE LIBRARY                    *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *           COPYRIGHT 1982, 1983, 1984, 1985, 1986               *
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

      LOGICAL IEST
C
C     MATRIX F IS  TREATED, INTERNAL TO THIS SUBROUTINE, AS
C     A VECTOR
C
      DIMENSION U(*), F(*), UALPHA(*), V(MAXV,M), H(MAXH,N)
      DIMENSION Z(M), XEST(N)
C
C     EPS = MACHINE APPROXIMATE ZERO
C
      DATA EPS/0.6D-37/
C
      MP1   = M+1
      MTR   = M*MP1/2
      JJOLD = 0
C
      IF ( .NOT. IEST) GO TO 30
C                                 COMPUTE Z = Z-H*XEST
         DO 20 J=1, N
            GAMMA = XEST(J)
            DO 10 I=1, M
               Z(I) = Z(I) - H(I,J)*GAMMA
  10        CONTINUE
  20     CONTINUE
  30  CONTINUE
C
C     CONSTRUCT F
C
      CALL PHIU(H,MAXH,M,N,U,N,F,MAXF)
      LJF = 0
      DO 100 J=1, N
         JJ = JJOLD+J
         LL = MTR
         DO 90  LB = 1, M
            L      = MP1-LB
            LJ     = LJF+L
            FLJ    = F(LJ)
C                                    ;(D(J)=U(JJ))
            Y      = FLJ*U(JJ)
C                                    ;(D_ALPHA(L)=U_ALPHA(LL))
            GAMMA  = Y*FLJ
            IF (GAMMA .GT. EPS) GO TO 50
C
C     D(J) * F(L,J)**2 IS EFFECTIVELY ZERO, SO IT NEED NOT BE NULLED
C
                LL     = LL-L
                V(J,L) = 0.0
                IF (U(JJ) .GT. 0) GO TO 90
                DO 40 K = 1, J
                   KJ    = JJOLD+K
                   U(KJ) = 0.0
   40           CONTINUE
                GO TO 90
   50       CONTINUE
            RHO2       = GAMMA + UALPHA(LL)
            GAMMA      = UALPHA(LL)/RHO2
            U(JJ)      = U(JJ)*GAMMA
            UALPHA(LL) = RHO2
            GAMMA      = Y/RHO2
            V(J,L)     = GAMMA
            IF (J .EQ. 1) GO TO 70
               JM1 = J-1
               DO 60  K  = 1, JM1
                  JK     = JJOLD+K
                  U(JK)  = U(JK)-V(K,L)*FLJ
                  V(K,L) = V(K,L)+GAMMA*U(JK)
   60         CONTINUE
   70       IF (L .EQ. 1) GO TO 90
C                                     ;LLN = LL(NEXT PRIOR DIAGONAL)
                LLN = LL-L
                LM1 = L-1
                DO 80  K = 1, LM1
                   LLNPK = LLN+K
                   KJ    = LJF + K
                   F(KJ) = F(KJ) - UALPHA(LLNPK)*FLJ
                   UALPHA(LLNPK) = UALPHA(LLNPK)+GAMMA*F(KJ)
  80            CONTINUE
                LL = LLN
  90     CONTINUE
         JJOLD = JJ
         LJF   = LJF + MAXF
100   CONTINUE
      IF (.NOT. IEST) RETURN
C                               COMPUTE Z=UALPHA**(-1)*Z
C      RESTST = LIKLIHOOD RESIDUAL (Z**T)*(PALPHA**-1)*Z
C
      RESTST = 0.0
      JJ    = MTR
      L     = M
      GAMMA = 0.0
      IF (M .EQ. 1) GO TO 180
         DO 170 I = 2, M
            Z(L)  = Z(L) - GAMMA
            IJ    = JJ - 1
            GAMMA = 0.0
            DO 160 J = L, M
               GAMMA = GAMMA+UALPHA(IJ) * Z(J)
               IJ    = IJ + J
 160        CONTINUE
            RESTST = RESTST + Z(L) * Z(L)/UALPHA(JJ)
                JJ = JJ-L
                L  = L-1
 170     CONTINUE
         Z(1) = Z(1) - GAMMA
 180  RESTST = RESTST + Z(1)*Z(1)/UALPHA(1)
C
C     UPDATE STATE ESTIMATE
C
      DO 200 J=1, M
         GAMMA = Z(J)
         DO 190 I = 1, N
            XEST(I) = XEST(I) + V(I,J)*GAMMA
 190     CONTINUE
 200  CONTINUE
      RETURN
      END
