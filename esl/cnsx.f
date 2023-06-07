      SUBROUTINE CNSX(XSTAR,N,A,IMAXA,NP,XPRED,EM,Q,IP,F,DX)
C
C     COLORED NOISE SMOOTHING ALGORITHM ESTIMATE UPDATE, FROM TIME
C     T BACK TO TIME T-DT
C
C *** INPUTS:
C
C     XSTAR(N)      SMOOTH ESTIMATE AT TIME T
C
C     N             STATE VECTOR DIMENSION
C
C     A(N,NP)       SMOOTHER GAIN MULTPLIER MATRIX GENERATED DURING
C                   FILTER COLORED NOISE TIME PROPAGATION STEP
C
C     IMAXA         ROW DIMENSION OF A
C
C     NP            NUMBER OF COLORED NOISE PARAMETERS
C
C     XPRED(N)      FILTER DETERMINISTIC TIME PREDICTION AT TIME T
C
C     EM(NP)        COLORED NOISE EXPONENTIAL MODEL MULTIPLIERS
C
C     Q(NP)         PROCESS  NOISE VARIANCES ASSOCIATED WITH THE
C                   COLORED  NOISE UPDATE
C
C     IP(NP)        PROCESS NOISE COMPONENT LOCATIONS
C
C     F(NP),DX(N)   WORK VECTORS
C
C *** OUTPUT
C
C     XSTAR(N)      UPDATED SMOOTHED ESTIMATE AT TIME T-DT
C
C     XPRED(N)      FILTER TIME PREDICTION AT TIME T, I.E. THE
C                   MARKOV STATE FILTER ESTIMATES AT TIME T-DT
C                   ARE TIME PROPAGATED TO TIME T, AND THE OTHER
C                   FILTER PREDICTION ESTIMATES ARE LEFT UNCHANGED.
C
C
C     COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JULY  3, 1981
C
C              REVISED 6 SEPTEMBER 1982         VERSION 1.00.00
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE CNSX   IS A PART OF THE FEA. INC             *
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
      DIMENSION XSTAR(N),A(IMAXA,NP),XPRED(N),EM(NP),
     1 Q(NP),IP(NP), F(NP),DX(N)
C
      DOUBLE PRECISION SUM
C
      NPP1=NP+1
      DO 10 JP=1,NP
         K=IP(JP)
         F(JP)=XPRED(K)
         XPRED(K)=EM(JP)*XPRED(K)
10    CONTINUE
C
C     SAVE MARKOV ESTIMATES AT TIME T-DT AND ADJUST XPRED
C
      DO 20 I=1,N
        DX(I)=XSTAR(I)-XPRED(I)
20    CONTINUE
      DO 40 L=1,NP
         JP=NPP1-L
         C=Q(JP)
         EX=EM(JP)
C
C      FORM DX,                          IN   EQN (36)
C
         MPL=IP(JP)
         GAMMA=EX * EX + C * A(MPL,JP)
C                                       GAMMA DEFINED BY EQN (22P)
         DX(MPL) = -F(JP)
         SUM = 0.0
         DO 30 I=1,N
            SUM = SUM + A(I,JP) * DX(I)
 30      CONTINUE
C                                 SEE EQN (24P), SUM =(V**T)*DX
C  UPDATE MARKOV STATE
C
         XSTAR(MPL)=(EX*XSTAR(MPL)-C*SUM)/GAMMA
C
C        NOTE GAMMA IS POSITIVE UNLESS EM AND Q ARE BOTH ZERO
C
C        ADJUST DX ELEMENT FOR NEXT CYCLE OF 40 LOOP
C
         DX(MPL)=XSTAR(MPL)+DX(MPL)
C
 40   CONTINUE
      RETURN
      END
