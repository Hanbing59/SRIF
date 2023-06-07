      SUBROUTINE PINVBM (U,N,B,NX,K,EM,Q,IEST,XEST,V,IERR)
C
C   SUBROUTINE TO COMPUTE SOLUTION V TO P*V=BM WHERE
C   P=U*D*U**T, VECTOR B HAS ZEROS BELOW ELEMENT NX
C   BM(J)=EM*B(J),J=1 TO NX, BM HAS AN IMPLICIT UNIT VALUE
C   IN LOCATION K, AND THIS CORRESPONDENDS TO THE LOCATION
C   OF THE COLORED NOISE MULTIPLIER, EM.
C   THIS METHOD USED EMPLOYS BACK SUBSTITUTION ON THE INPUT
C   VECTOR STORED U-D FACTORS.
C
C
C ** INPUTS:
C
C    U (N*N+1)/2)         VECTOR STORED U-D COVARIANCE FACTORS
C                         SINGULAR D FACTORS ARE ALLOWED
C                         PROVIDED CORRESPONDING ROW & COLUMN
C                         ENTRIES OF U ARE ZERO
C
C    N                    MATRIX DIMENSION
C
C    B(NX)                B VECTOR INPUTS
C
C    NX                   DIMENSION OF NON TRIVIAL PORTION
C                         OF B
C
C    K                    LOCATION OF COLORED PROCESS NOISE STATE
C                         IN FILTER STATE VECTOR
C                         NX .LT. K .LE. N
C
C    EM                   EXPONENTIAL (COLORED NOISE) MULTIPLIER
C
C    Q                    COLORED NOISE VARIABLE PROCESS NOISE VARIANCE
C
C    IEST                 LOGICAL; IF(TRUE) COMPUTE BETA =
C                         (V**T)*XEST  ELSE INPUT XEST IS NOT USED
C
C    XEST(N)              FILTER ESTIMATE STATE VECTOR
C
C  ** OUTPUT:
C
C    V(N)                 V = P(INVERSE)*BM
C
C    V(N+1)               = (V**T)*XEST, WHEN IEST = TRUE
C                                        UNCHANGED OTHERWISE
C
C    IERR                 ERROR RETURN FLAG
C
C                            .EQ. 0    NORMAL RETURN
C
C                            .EQ. J    ESTIMATE COMPUTATION
C                                      FAILED BECAUSE DIAGONAL "J"
C                                      OF U ARRAY IS TOO SMALL
C                                      BETA NOT COMPUTED
C
C *** REQUIRED SUBROUTINE:
C
C     DABS    <FORLIB>
C
C    COGNIZANT PERSONS:
C
C       DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C       FACTORIZED ESTIMATION APPLICATIONS INC.
C       MAY 20, 1986
C
C                               VERSION 1.00.00       2 OCTOBER 1986
C
C ******************************************************************
C *                                                                *
C *           SUBROUTINE PINVBM IS PART OF THE FEA, INC            *
C *                 ESTIMATION SUBROUTINE LIBRARY                  *
C *                                                                *
C *                            (ESL)                               *
C *                                                                *
C *            COPYRIGHT 1982, 1983, 1984, 1985 , 1986             *
C *            FACTORIZED ESTIMATION APLLICATIONS INC.             *
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

      DIMENSION U(*),V(*),B(*),XEST(N)
      DOUBLE PRECISION S,EPS,TMP
      LOGICAL IEST
C
C      DATA EPS /0.6D-37/ !suitable for single precision
      DATA EPS/1.0D-307/
      DATA EPSCN/1.0D-7/
C
C    EPS -- MACHINE EPSILON, SUCH THAT RECIPROCAL EXISTS.
C    EPS =  0.6D-37 IS ACCEPTABLE ON VAX TYPE MACHINES
C
C    EPSCN IS USED TO TEST WHETHER THE COLORED NOISE PARAMETER
C           IS A CONSTANT.
C
C    WHEN THERE IS NO PROCESS NOISE, MODIFIED GAIN IS ESSENTIALLY ZERO
C
      XMONE=-1.0
      IF (EM .GT. (1.0D0-EPSCN)) THEN
         SIG = 1.0
         IF (Q .GT. 0.0) CALL U2SIG(U,N,' ',' ',XMONE,' ',.FALSE.,K,1,
     1                              SIG,' ',N,0)
         IF ((Q .EQ. 0.0) .OR. (Q .LT. (EPSCN*SIG)**2)) THEN
            DO 1 J = 1, N+1
               V(J) = 0.0
  1         CONTINUE
            RETURN
         ENDIF
      ENDIF
C
C    IF N IS ONE SKIP THE BULK OF THE CODE
C
      IERR = 0
      IF (N.EQ.1) THEN
         J   = 1
         S   = U(1)
         IF (S .LE. EPS) GO TO 90
             V(1) = 1.0/S
             V(N+1) = XEST(1)*V(1)
             RETURN
      ENDIF
C
C    STEP-1 SOLVE U*V = BM
C
      IF ((EM .EQ. 0.0) .OR. (NX .EQ. 0)) THEN
         DO 5 J=1,K
            V(J) = 0.0
    5    CONTINUE
      ELSE
         DO 10 J=1,NX
            V(J) = EM*B(J)
   10    CONTINUE
         DO 20 J = NX+1,K
            V(J) = 0.0
   20    CONTINUE
      ENDIF
      V(K) = 1.0
      IF (K .GT. 1) THEN
         JK = K*(K+1)/2
         DO 40 L=K,2,-1
            S = V(L)
            DO 30 J=L-1,1,-1
               JK = JK-1
               V(J) = V(J)-U(JK)*S
   30       CONTINUE
            JK = JK-1
   40    CONTINUE
      ENDIF
C
C    STEP-2 V=D**(-1)*V
C
      JJ = 0
      DO 50 J=1,K
         JJ = JJ+J
          S = U(JJ)
        TMP = V(J)
        IF (S .GT. 0.0) THEN
           IF (S .LE. EPS*DABS(TMP)) GO TO 90
              V(J) = TMP/S
        ELSE
           V(J) = 0.0
        ENDIF
   50 CONTINUE
C
C   STEP-3 V=U**(-T)*V
C
      IF (K .LT. N) THEN
         KP1 = K+1
         DO 60 I=KP1,N
            V(I) = 0.0
   60    CONTINUE
      ENDIF
      IJ = 1
      DO 80 I=2,N
         S = V(I)
         IM1 = I-1
         DO 70 J=1,IM1
            IJ = IJ+1
            S = S-U(IJ)*V(J)
   70    CONTINUE
         V(I) = S
         IJ = IJ+1
   80 CONTINUE
      IF (IEST) THEN
         S = 0.0
         DO 95 L=1,N
            S = S+V(L)*XEST(L)
   95    CONTINUE
         V(N+1) = S
      ENDIF
      RETURN
C
C   ERROR RETURN
C
   90 IERR = J
C
      RETURN
      END
