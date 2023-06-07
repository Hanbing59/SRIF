      SUBROUTINE UCOLA(U,N,KS,V,EM,Q,IERR)
C
C     COLORED NOISE U-D COVARIANCE FACTOR UPDATING
C     U*D*(U**T)-OUTPUT=PHI*U*D*(U**T)*(PHI**T)+Q
C     PHI=DIAG(E(KS-1),EM,...,E(N-(KS-1+NCOLOR)))
C     E(K) IS A K-VECTOR WITH UNIT ELEMENTS.
C
C     THE ALGORITHM USED IS THE BIERMAN-THORNTON ONE COMPONENT
C     AT-A-TIME UPDATE. CF.BIERMAN "FACTORIZATION METHOD
C     FOR DISCRETE SEQUENTIAL ESTIMATION", ACADEMIC PRESS (1977)
C     PP.147-148
C
C *** INPUTS:
C
C     U(N*(N+1)/2)      INPUT U-D VECTOR STORED COVARIANCE FACTORS.
C                       THE COLORED NOISE UPDATE RESULT RESIDES
C                       IN U ON OUTPUT
C
C     N                 FILTER DIMENSION. WHEN THE LAST COLUMN OF U
C                       HOUSES THE FILTER ESTIMATES, THEN
C                       N=NUMBER FILTER VARIABLES + 1
C
C     KS                THE LOCATION OF THE FIRST COLORED NOISE TERM
C                       (KS.GE.1.AND.KS.LE.N)
C
C     NCOLOR            THE NUMBER OF COLORED NOISE TERMS
C                       (NCOLOR.GE. 1)
C
C     V(KS-1+NCOLOR)    PROCESS NOISE MULTIPLIER VECTOR
C                       (DESTROYED BY PROGRAM)
C                       IT IS ASSUMED THAT V(KS) IS IMPLICITLY
C                       UNITY, SO THAT Q WILL CORRESPOND TO
C                       THE COLORED NOISE VARIABLE PROCESS NOISE
C                       VARIANCE. WHEN V(KS) = -1.0 THE USER IS
C                       SPARED THE EFFORT OF DEFINING THE OTHER
C                       V-ELEMENTS. THEY ARE, IN THIS CASE
C                       ASSUMED TO BE ZERO.
C
C     EM                COLORED NOISE MAPPING TERM
C                       (UNALTERED BY PROGRAM)
C
C     Q                 PROCESS NOISE VARIANCE .GE. 0
C                       (UNALTERED BY PROGRAM)
C
C *** OUTPUTS:
C
C     U(N*(N+1)/2)      THE COLORED NOISE UPDATE RESULT RESIDES
C                       IN U ON OUTPUT
C
C     IERR            DIAGNOSTIC FLAG (FROM RNK1)
C
C                     IF  .EQ.  0       NORMAL RETURN
C
C                               J       ERROR RETURN DUE TO A
C                                       NEGATIVE COMPUTED J-TH
C                                       DIAGONAL ENTRY
C
C                              -1       UD RESULT IS SINGULAR
C
C ***  REQUIRED SUBROUTINE:
C
C               RNK1  <ESL>
C
C
C     COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JANUARY 27, 1985
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE UCOLA  IS A PART OF THE FEA. INC             *
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

      DOUBLE PRECISION S, EPS, TMP
      DIMENSION U(*),V(*)
C
C      DATA EPS/0.6D-37/ ! suitable for single precision
      DATA EPS/1.0D-307/

C
C     IN SINGLE PRECISION EPS(ILION) IS MACHINE ACCURACY
C     ADEQUATE FOR VAX TYPE MACHINES
C
C * * * * * * INITIALIZATION
C
      KSM1  = KS-1
      JJOLD = KS*KSM1/2
      JJ    = JJOLD+KS
      NM1   = N-1
C
C * * * * * *
C
      TMP = EM*U(JJ)
      S   = TMP*EM+Q
C
C     S = EM*EM*D(KS)+Q
C
      IF (KS .LT. N) THEN
C
C        MULTIPLY ROW KS OF U BY EM
C
         IJ  = JJ
         DO 10 J = KS,NM1
            IJ = IJ+J
            U(IJ) = U(IJ)*EM
  10     CONTINUE
      ENDIF
      IF (KS .GT. 1) THEN
         IF (S .LE. EPS) THEN
            DO 20 J = 1, KS
               U(J+JJOLD) = 0.0
  20        CONTINUE
         ELSE
            C     = Q*U(JJ)/S
            U(JJ) = S
            TMP   = TMP/S
C
C THIS ARRANGEMENT FOR C AND TMP IS WELL DEFINED FOR
C ALL Q .GE. 0 AND ABS(EM) .LE. 1)
C
C  SET V ELEMENTS IN THE CASE THAT V(KS) = -1.0
C
            IF (V(KS) .EQ. -1.0) THEN
               DO 30 J = 1, KSM1
                  L    = JJOLD+J
                  V(J) = U(L)
                  U(L) = TMP*V(J)
  30           CONTINUE
            ELSE
               DO 40 J = 1, KSM1
                  L    = JJOLD+J
                  S    = V(J)
                  V(J) = U(L)-EM*S
                  U(L) = S+TMP*V(J)
C
C    V(NEW) = COLUMN KS OF U(INPUT)-EM*V(OLD)
C             COLUMN KS OF U(OUT) = V(INPUT)+EM*RHO*V(NEW)
C             WHERE RHO = D(OLD)/D(NEW)
C
  40           CONTINUE
            ENDIF
            CALL RNK1(U,KSM1,C,V,IERR)
         ENDIF
      ELSE
        U(1) = S
      ENDIF
      RETURN
      END
