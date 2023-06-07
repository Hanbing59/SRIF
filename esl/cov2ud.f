      SUBROUTINE COV2UD (U,N,IERR,IOUNIT)
C
C        TO OBTAIN THE U-D FACTORS OF A POSITIVE
C        SEMI-DEFINITE MATRIX. THE INPUT VECTOR STORED MATRIX
C        IS OVERWRITTEN BY THE OUTPUT U-D FACTORS WHICH ARE
C        ALSO VECTOR STORED.
C
C *** INPUTS:
C
C        U(N*(N+1)/2)  VECTOR STORED COVARIANCE (POSITIVE
C                      SEMI-DEFINITE)MATRIX.
C
C        N             MATRIX DIMENSION, N .GE. 1
C
C        IOUNIT        UNIT NUMBER OF OUTPUT DEVICE
C
C *** OUTPUTS:
C
C        U(N*(N+1)/2)  VECTOR STORED U-D COVARIANCE FACTORS
C                      NOTE: OVERWRITING INPUT U
C
C        IERR          ERROR RETURN FLAG
C
C                        .EQ. 0       NORMAL RETURN
C
C                        .EQ. L       COMPUTED NEGATIVE DIAGONALS
C                                     THE ALGORITHM ASSUMES THAT
C                                     THE USER INTENDED THE INPUT
C                                     TO BE AT LEAST POSITIVE
C                                     SEMI-DEFINITE, AND THAT
C                                     NEGATIVE VALUES CORRESPOND
C                                     TO COMPUTER ROUNDOFF. THE
C                                     NEGATIVE DIAGONAL COMPUTED
C                                     VALUES ARE DISPLAYED AND
C                                     THE CORRESPONDING COLUMN
C                                     (INCLUDING THE DIAGONAL) IS
C                                     SET TO ZERO.
C
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              18 OCTOBER 1981
C
C                                      VERSION 1.00.00
C       REVISED 13 JUNE 1984           VERSION 1.01.00
C       REVISED 20  MAY 1985           VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE COV2UD IS A PART OF THE FEA. INC             *
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
      DIMENSION U(*)
      DOUBLE PRECISION S, EPS, Z, ONE, ALPHA
C
      DATA Z/0.0D0/, ONE /1.0D0/
C     DATA      EPS/0.6D-37/  ! single precision
C
C     EPS -- MACHINE EPSILON, SUCH THAT RECIPROCAL EXISTS.
C     EPS =  0.6D-37 IS ACCEPTABLE ON VAX TYPE MACHINES

      DATA EPS/1.0D-307/
C
      IERR = 0
      IF (N .EQ. 1) GO TO 60
      JJ=N*(N+1)/2
      NP2=N+2
      DO 50 L=2,N
        J=NP2-L
        ALPHA=Z
        S = U(JJ)
        IF (S.GE.Z) GO TO 10
           IF (IOUNIT .NE. 0) WRITE (IOUNIT,100) J,S
           IERR = IERR+1
           U(JJ)=Z
   10   IF (S.GT.EPS) ALPHA=ONE/S
        JJ=JJ-J
        KK=0
        KJ=JJ
        JM1=J-1
        DO 40 K=1,JM1
          KJ=KJ+1
          BETA=U(KJ)
          U(KJ)=ALPHA*U(KJ)
          IJ=JJ
          IK=KK
          DO 30 I=1,K
            IK=IK+1
            IJ=IJ+1
   30       U(IK)=U(IK)-BETA*U(IJ)
   40     KK=KK+K
   50   CONTINUE                                            
   60 CONTINUE
      S = U(1)
      IF (S.GE.Z) GO TO 70
         J = 1
         IF (IOUNIT .NE. 0) WRITE (IOUNIT,100) J, S
         IERR = IERR+1
         U(1)=Z
   70 RETURN
C
  100 FORMAT (1H0,20X,' AT STEP ',I4,' DIAGONAL ENTRY = ',1PD14.7)
      END
