      SUBROUTINE INF2R (R,N,IERR,IOUNIT)
C
C         TO CHOLESKY FACTOR AN INFORMATION MATRIX
C
C     COMPUTES A LOWER TRIANGULAR VECTOR STORED CHOLESKY
C     FACTORIZATION OF A POSITIVE SEMI-DEFINITE MATRIX.
C     R=R(**T)R, R UPPER TRIANGULAR. BOTH MATRICES ARE VECTOR
C     STORED  AND THE RESULT OVERWRITES THE INPUT
C
C *** INPUTS:
C
C     R(N*(N+1)/2)  ON INPUT THIS IS A POSITIVE SEMI-DEFINITE
C                   (INFORMATION) MATRIX.
C
C     N             DIMENSION OF MATRICES INVOLVED, N.GE.1
C
C     IOUNIT        UNIT NUMBER OF OUTPUT DEVICE
C
C *** OUTPUT:
C
C     R(N*(N+1)/2)  ON OUTPUT IT IS THE  TRANSPOSED LOWER
C                   TRIANGULAR CHOLESKY FACTOR.
C                   IF THE INPUT MATRIX IS SINGULAR THE OUTPUT
C                   MATRIX WILL HAVE ZERO DIAGONAL ENTRIES
C
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
C *** REQUIRED SUBROUTINE:
C
C       DSQRT   <FORLIB>
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23, 1981
C
C              REVISED 5 JULY 1982
C              REVISED 6 SEPTEMBER 1982       VERSION 1.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE INF2R  IS A PART OF THE FEA. INC             *
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

      DIMENSION R(*)
      DOUBLE PRECISION SUM, EPS, Z, ONE, ALPHA
C
      DATA Z/0.D0/, ONE /1.0D0/
C      EPS/0.6D-37/ ! suitable for single precision
C
C     EPS -- MACHINE EPSILON, SUCH THAT RECIPROCAL EXISTS.
C     EPS =  0.6D-37 IS ACCEPTABLE ON VAX TYPE MACHINES
      DATA EPS/1.0D-307/

C
      IERR = 0
      JJ=0
      NN=N*(N+1)/2
      IF (N .EQ. 1) GO TO 20
      NM1=N-1
      DO 10 J=1,NM1
C                                          @ JJ=(J,J)
        JJ=JJ+J
        SUM = R(JJ)
        IF (SUM.GE.Z)  GO TO 5
           WRITE (IOUNIT,100) J,SUM
           IERR = IERR + 1
           R(JJ)=Z
    5   SUM = DSQRT(SUM)
C
        ALPHA=Z
        IF (SUM .GT. EPS) ALPHA=ONE/SUM
C                                          @ JK=(J,K)
        JK    = NN+J
        JP1   = J+1
C                                          @ JIS=(J,I) START
        JIS   = JK
        NPJP1 = N+JP1
        DO 10 L=JP1,N
          K     = NPJP1-L
          JK    = JK-K
          R(JK) = ALPHA*R(JK)
          BETA  = R(JK)
          KI    = NN+K
          JI    = JIS
          NPK   = N+K
          DO 10 M=K,N
            I     = NPK-M
            KI    = KI-I
            JI    = JI-I
            R(KI) = R(KI)-R(JI)*BETA
   10     CONTINUE
C
   20 SUM = R(NN)
      IF (SUM .GE. Z) GO TO 30
         WRITE (IOUNIT,100) N,SUM
         R(NN) = Z
         IERR  = IERR+1
   30 R(NN)= DSQRT(SUM)
      RETURN
C
 100  FORMAT (1H0,20X,' AT STEP  ',I4,5X,' DIAGONAL ENTRY = ',
     1 1P1D14.7,2X,' IT IS RESET TO ZERO ')
      END
