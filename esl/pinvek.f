      SUBROUTINE PINVEK(U,N,K,V,IERR)
C
C     SUBROUTINE TO COMPUTE SOLUTION V TO P*V=E(SUB K) WHERE
C     P=U*D*U**T AND E(SUB K) HAS UNITY IN ELEMENT K
C     I.E. TO COMPUTE COLUMN K OF P INVERSE.
C     THE METHOD USED EMPLOYS BACK SUBSTITITION ON THE INPUT
C     VECTOR STORED U-D FACTORS.
C
C *** INPUTS:
C
C     U (N*N+1)/2)         VECTOR STORED U-D COVARIANCE FACTORS
C                          ( MUST BE NONSINGULAR )
C
C     N                    MATRIX DIMENSION
C
C     K                    LOCATION OF THE UNIT ELEMENT IN E(SUB K)
C
C
C
C *** OUTPUT:
C
C     V(N)                 COLUMN K OF P INVERSE
C
C     IERR                 ERROR RETURN FLAG
C
C                         .EQ. 0   NORMAL RETURN
C
C                            .EQ. J     ESTIMATE COMPUTATION
C                                       FAILED BECAUSE DIAGONAL "J"
C                                       OF U ARRAY IS TOO SMALL
C
C *** REQUIRED SUBROUTINE:
C
C      DABS    <FORLIB>
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JUNE 30, 1981
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C              REVISED  5      MAY 1986         VERSION 2.01.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE PINVEK IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
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
 
      DIMENSION     U(*),V(N)
      DOUBLE PRECISION S,EPS,TMP
C
C      DATA EPS /0.6D-37/ ! suitable for single precision
      DATA EPS/1.0D-307/

C
C      IN SINGLE PRECISION EPSILON IS MACHINE ACCURACY
C      ADEQUATE FOR VAX TYPE APPLICATIONS
C
C     IF N IS ONE SKIP THE BULK OF HE CODE
C
      IERR = 0
      IF (N .GT. 1) GO TO 1
         J = 1
         S = U(1)
         IF (S .LE. EPS) GO TO 90
         V(1)=1.0/S
         RETURN
C
C     STEP-1 SOLVE U*V = E(SUB K)
C
 1    V(K)=1.0
      IF ( K .EQ. 1) GO TO 35
      KM1=K-1
      JJ=K*KM1/2
C
      DO 10 I=1,KM1
         IT=JJ+I
         V(I)= -U(IT)
 10   CONTINUE
C
      IF (K .EQ. 2) GO TO 35
      DO 30 M=2,KM1
         JM1=K-M
         J=JM1+1
         JJ=JJ-J
         S=V(J)
         DO 20 I=1,JM1
            IT=JJ+I
            V(I)=V(I)-U(IT)*S
 20      CONTINUE
 30   CONTINUE
C
C     STEP-2 V=D**(-1)*V
C
 35   JJ=0
      DO 40 J=1,K
         JJ  = JJ+J
         S   = U(JJ)
         TMP = V(J)
         IF ((S .EQ. 0.0) .AND. (TMP .EQ. 0.0)) THEN
            V(J) = 0.0
         ELSE
           IF (S.LE. EPS*DABS(TMP)) GO TO 90
           V(J)= TMP/S
         END IF
 40   CONTINUE
C
C     STEP-3 V=U**(-T)*V
C
      IF (K .GE. N) GO TO 60
         KP1=K+1
         DO 50 I=KP1,N
            V(I)=0.0
 50      CONTINUE
 60   IJ=1
      DO 80 I=2,N
         S=V(I)
         IM1=I-1
         DO 70 J=1,IM1
            IJ=IJ+1
            S=S-U(IJ)*V(J)
 70      CONTINUE
         V(I)=S
         IJ=IJ+1
 80   CONTINUE
      RETURN
C
C ***  ERROR RETURN ***
C
 90   IERR = J
      RETURN
      END
