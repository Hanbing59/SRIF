      SUBROUTINE RINZ(R,N,Z,X,IERR)
C
C     TO COMPUTE THE SOLUTION VECTOR X TO THE MATRIX EQUATION
C     R*X=Z, WHERE R IS AN N BY N UPPER-TRIANGULAR VECTOR
C     MATRIX. COMPUTATION IS MINIMIZED AND ACCURACY ENHANCED
C     BY NOT COMPUTING THE  MATRIX INVERSE. IF R IS SINGULAR,
C     WITH ITS JTH DIAGONAL ELEMENT AS ITS LAST DIAGONAL ZERO
C     VALUE, X(I) IS COMPUTED  ONLY FOR J .LT. I .LE. N,
C     THE OTHER COMPONENTS,   X(I), ARE LEFT UNCHANGED FOR
C     1 .LE. I .LE. J, AND THE VALUE OF J IS RETURNED AS IERR
C
C  *** INPUTS:
C
C     R(N*(N+1)/2)        VECTOR STORED UPPER TRIANGULAR
C                         MATRIX
C
C     N                   MATRIX VECTOR DIMENSION, N .GT. 0
C
C     Z(N)                RIGHT HAND SIDE VECTOR OF EQUATION
C                         R*X=Z
C
C
C  *** OUTPUT:
C
C     X(N)                SOLUTION VECTOR, WHICH MAY OVERWRITE Z,
C                         COMPILER PERMITTING.
C
C     IERR                ERROR RETURN FLAG
C
C                           .EQ. 0    NORMAL RETURN
C
C                           .EQ. J    SOLUTION, X, VALID ONLY FOR
C                                     COMPONENTS .GT. J
C
C
C *** REQUIRED SUBROUTINE:
C
C       DABS    <FORLIB>
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              01 MAY 1978
C
C              REVISED  3 AUGUST 1982         VERSION 0.01.00
C                                             VERSION 1.00.00
C              REVISED 20    MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE RINZ  IS A PART OF THE FEA. INC              *
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
C
      IMPLICIT REAL*8 (A-H,O-Z)

      DOUBLE PRECISION SUM, EPS, S, DX, ZERO
      DIMENSION R(*), Z(N), X(N)
C
C       DATA EPS /0.6D-37/ ! suitable for single precision
      DATA EPS/1.0D-307/
      DATA ZERO/0.0D0/
C
C      IN SINGLE PRECISION EPSILON IS MACHINE ACCURACY
C      ADEQUATE FOR VAX TYPE APPLICATIONS
C
      IERR = 0
      JJ   = N*(N+1)/2
      IX   = N
      SUM  = ZERO
      IF (N .EQ. 1) GO TO 30
         DO 20 I = 2,N
            S  = R(JJ)
            DX = Z(IX) - SUM
            IF (DABS(S) .LE. EPS * DABS(DX)) GO TO 90
            X(IX) = DX/S
            IJ    = JJ-1
            SUM   = ZERO
            DO 10 J = IX, N
               SUM  = SUM+R(IJ)*X(J)
               IJ   = IJ+J
  10        CONTINUE
            JJ = JJ-IX
            IX = IX-1
  20     CONTINUE
  30  CONTINUE
      S = R(1)
      DX = Z(1) - SUM
      IF (DABS(S) .LE. EPS * DABS(DX)) GO TO 90
      X(1) = DX/S
         RETURN
  90  CONTINUE
C
C *** ERROR RETURN ***
C
      IERR = IX
C
      RETURN
      END
