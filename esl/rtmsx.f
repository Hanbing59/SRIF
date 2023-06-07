      SUBROUTINE RTMSX(R,N,X,V)
C
C     SUBROUTINE TO MULTIPLY R (A VECTOR STORED MATRIX)
C     AND X (A VECTOR)
C
C
C *** INPUTS:
C
C     R(N*(N+1)/2)   INPUT VECTOR STORED MATRIX
C
C     N              DIMENSION OF INPUT VECTORS
C
C     X(N)           INPUT VECTOR
C
C *** OUTPUT:
C
C     V(N)           OUTPUT VECTOR (RESULT)
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              APRIL 10,1981
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE RTMSX  IS A PART OF THE FEA. INC             *
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

      DIMENSION R(*),X(N),V(N)
C
      II = 0
      DO 35 I=1,N
         SUM = 0.0
         II = II + I
C
C     IK =(I,K)
C
         IK = II
         DO 30 K = I,N
            SUM = SUM + R(IK)*X(K)
            IK=IK+K
 30      CONTINUE
         V(I) = SUM
 35   CONTINUE
      RETURN
      END
