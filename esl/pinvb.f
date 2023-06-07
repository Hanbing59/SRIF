      SUBROUTINE PINVB(U,N,B,K,V,IERR)
C
C     SUBROUTINE TO COMPUTE SOLUTION V TO P*V=B  WHERE
C     P=U*D*U**T AND VECTOR B HAS ZEROS BELOW ELEMENT K
C     THE METHOD USED EMPLOYS BACK SUBSTITITION ON THE INPUT
C     VECTOR STORED U-D FACTORS.
C
C
C ** INPUTS:
C
C     U (N*N+1)/2)         VECTOR STORED U-D COVARIANCE FACTORS
C                          (SINGULAR D FACTOR COMPONENTS ARE
C                           ALLOWED PROVIDED CORRESPONDING B
C                           COMPONENTS AND ROW/COLUMN ENTRIES
C                           OF U ARE ALSO ZERO)
C
C     N                    MATRIX DIMENSION
C
C     B(K)                 B VECTOR INPUTS
C
C     K                    DIMENSION OF B
C
C
C ** OUTPUT:
C
C     V(N)                 V = P(INVERSE)*B
C
C
C     IERR                 ERROR RETURN FLAG
C
C                            .EQ. 0     NORMAL RETURN
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
C        DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C        FACTORIZED ESTIMATION APPLICATIONS INC.
C        NOVEMBER 23, 1981
C
C        REVISED 20      MAY        1985         VERSION 2.00.00
C        REVISED 10      SEPTEMBER  1985         VERSION 2.01.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE PINVB  IS A PART OF THE FEA. INC             *
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

      DIMENSION U(*),V(N),B(K)
      DOUBLE PRECISION S,EPS,TMP
C
C       DATA EPS /0.6D-37/ ! suitable for single precision
      DATA EPS/1.0D-307/

C
C
C     EPS -- MACHINE EPSILON, SUCH THAT RECIPROCAL EXISTS.
C     EPS =  0.6D-37 IS ACCEPTABLE ON VAX TYPE MACHINES
C
C     IF N IS ONE SKIP THE BULK OF THE CODE
C
      IERR = 0
      IF (N .EQ. 1) THEN
         J   = 1
         S   = U(1)
         TMP = B(1)
         IF (S .EQ. 0.0 .AND. TMP .EQ. 0.0) THEN
            V(1) = 0.0
         ELSE
            IF (S .LE. DABS(TMP)*EPS) GO TO 90
            V(1) = TMP/S
         ENDIF
         RETURN
      ENDIF
C
C     STEP-1 SOLVE U*V = B
C
      V(K)=B(K)
      IF ( K .EQ. 1) GO TO 35
      KM1=K-1
      JJ=K*KM1/2
      S = V(K)
      DO 10 I=1,KM1
         IT=JJ+I
         V(I)= B(I)-U(IT)*S
 10   CONTINUE
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
         JJ   = JJ+J
         S   = U(JJ)
         TMP = V(J)
         IF (S .EQ. 0.0 .AND. TMP .EQ. 0.0) THEN
            V(J) = 0.0
         ELSE
            IF (S .LE. EPS*DABS(TMP)) GO TO 90
            V(J) = TMP/S
         ENDIF
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
C *** ERROR RETURN ***
C
 90   IERR = J
C
      RETURN
      END
