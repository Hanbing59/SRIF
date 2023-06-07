      SUBROUTINE COVCOR(COV,N,RHO,V)
C
C     TO COMPUTE THE CORRELATION MATRIX RHO, FROM AN INPUT
C     COVARIANCE MATRIX COV. BOTH MATRICES ARE UPPER TRIANGULAR
C     VECTOR STORED. IF THE COMPILER ALLOWS THE CORRELATION MATRIX
C     RESULT CAN  OVERWRITE THE INPUT COVARIANCE
C
C *** INPUTS:
C
C     COV(N*(N+1)/2)  INPUT VECTOR STORED POSITIVE SEMI-DEFINITE
C                     COVARIANCE MATRIX
C
C     N               NUMBER OF PARAMETERS, N.GE.1
C
C *** OUTPUT:
C
C     RHO(N(N+1)/2)   OUTPUT VECTOR STORED CORRELATION MATRIX,
C                     RHO(IJ)=COV(IJ)/(SIGMA(I)*SIGMA(J))
C *** SCRATCH:
C
C     V(N)            WORK VECTOR
C
C
C *** REQUIRED SUBROUTINE:
C
C     DSQRT   <FORLIB>
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              18 OCTOBER 1981
C
C              REVISED 20 OCTOBER 1982        VERSION 1.00.00
C              REVISED 20     MAY 1985        VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE COVCOR IS A PART OF THE FEA. INC             *
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
      DIMENSION  COV(*), RHO(*), V(N)
      DOUBLE PRECISION EPS, S, Z, ONE
C
C       DATA EPS /0.6D-37/ ! suitable for single precision
	DATA Z/0.0D0/, ONE/1.0D0/
	DATA EPS/1.0D-307/

C
C      EPS  MACHINE EPSILON, SUCH THAT RECIPROCAL EXISITS
C      EPS = 0.6D-37 IS ACCEPTABLE ON VAX TYPE MACHINES
C
      JJ=0
      DO 10 J=1,N
         JJ=JJ+J
         V(J)=Z
         S = COV(JJ)
         IF (S .GT. EPS) V(J)=ONE/ DSQRT(S)
  10  CONTINUE
      IJ=0
      DO 30 J=1,N
         S=V(J)
         DO 20 I=1,J
            IJ=IJ+1
            RHO(IJ)=COV(IJ)*S*V(I)
  20     CONTINUE
      RHO(IJ) = 1.0
  30  CONTINUE
      RETURN
      END
