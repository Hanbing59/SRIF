      SUBROUTINE UCON (U,N,H,R,G,NG,F,V,IERR)
C
C
C     TO COMPUTE THE GENERAL COVARIANCE FACTOR UPDATE
C     (I-GH) P (I-GH)**T + R * G*G ** T
C
C     WHERE THE VECTOR G HAS IMPLICIT ZEROS FOLLOWING ITS NG-TH
C     COMPONENT, AND P = U * D * U**T .
C
C *** INPUTS:
C
C     U(N*(N+1)/2)          VECTOR STORED U-D COVARIANCE FACTORS
C                           WITH D ELEMENTS STORED ON THE DIAGONAL
C
C     N                     FILTER DIMENSION, N .GE. 1
C
C     H(N)                  MEASUREMENT COEFFICIENTS
C
C     R                     MEASUREMENT VARIANCE, R .GE. 0.0
C
C     G(NG)                 NONTRIVIAL PORTION OF THE GAIN VECTOR,
C                           IMPLICITLY G(K) = 0.0 FOR K .GT. NG;
C                           NG .GE. 1
C
C *** OUTPUT:
C
C     U(N*(N+1)/2)          OUTPUT VECTOR STORED U-D COVARIANCE
C                           FACTORS WITH D ELEMENTS STORED
C                           ON THE DIAGONAL
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
C
C *** SCRATCH:
C
C     F(N)                  SCRATCH VECTOR
C
C     V(NG)                 SCRATCH VECTOR
C
C
C *** SUBROUTINES REQUIRED :
C
C              UDOPT    <ESL>
C              RNK1     <ESL>
C
C
C          COGNIZANT PERSONS:
C
C
C             DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C             FACTORIZED ESTIMATION APPLICATIONS INC.
C
C             AUGUST 9, 1981
C
C             REVISED  11 AUGUST 1982     VERSION  0.01.00
C                                         VERSION  1.00.00
C              REVISED 20    MAY 1985     VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE UCON   IS A PART OF THE FEA. INC             *
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

      DIMENSION U(*),H(N),G(NG),F(N),V(NG)
C
      DOUBLE PRECISION SUM
C
      IF (N .GT. 1) GO TO 10
         U(1) = (1.0 - G(1) * H(1))**2 * U(1) + R * G(1)**2
         GO TO 100
  10  CONTINUE
C
C     COMPUTE F**T = HU
C
      F(1)=H(1)
      JJ=1
      DO 30 K=2,N
         SUM = H(K)
         KM1 = K-1
         DO 20 J=1,KM1
            SUM = SUM + H(J) * U(JJ+J)
  20     CONTINUE
C
         F(K) = SUM
         JJ = JJ+K
  30  CONTINUE
C
C     COMPUTE V = D*F1
C
      JJ=0
      DO 40 J = 1,NG
         JJ=JJ+J
         V(J) = U(JJ) * F(J)
  40  CONTINUE
C
C     UPDATE U12
C
      IF (NG .GE. N) GO TO 70
         NGP1  = NG+1
         JJ    = NG * (NGP1)/2
         DO 60 J = NGP1,N
            SUM = F(J)
            DO 50 I = 1,NG
               U(JJ + I) = U(JJ + I) - G(I) * SUM
  50        CONTINUE
         JJ = JJ + J
  60     CONTINUE
C
C      UPDATE U1-D1 ELEMENTS
C
  70  CONTINUE
      CALL UDOPT(U,NG,R,F,V,ALPHA)
      IF (ALPHA .EQ. 0.0) GO TO 100
         DO 80 J=1,NG
            V(J) = G(J) - V(J)/ALPHA
  80     CONTINUE
         CALL RNK1 (U,NG,ALPHA,V,IERR)
 100  CONTINUE
      RETURN
      END
