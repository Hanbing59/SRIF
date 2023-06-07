      SUBROUTINE GIVENS(A,B,SIG,S,C)
C
C     TO COMPUTE THE GIVENS ROTATION PARAMETERS SIG,S,C, SUCH THAT
C
C     |  C   S |  | A |   | SIG |
C     |        |  |   | = |     |
C     |  S  -C |  | B |   |  0  |
C
C     THE DEGENERATE CASE (WITH B=0) PRODUCES C = 1.0*SIGN(A),
C     S = 0.0, AND SIG = ABS(A)
C
C *** INPUT:
C
C     A,B             VECTOR COEFFICIENTS
C
C *** OUTPUT:
C
C     SIG             SQRT(A*A + B*B)
C     C               A/SIG
C     S               B/SIG
C
C       REQUIRED SUBROUTINES:
C
C         ABS  <FORLIB>
C         SIGN <FORLIB>
C         SQRT <FORLIB>
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              25 MARCH 1984
C
C              REVISED 2 JUNE 1985     (FORTRAN '77 LOGIC)
C                                       VERSION 2.00.00
C ******************************************************************
C *                                                                *
C *        SUBROUTINE GIVENS IS A PART OF THE FEA. INC             *
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
      ONE=1.0
      SIG = ABS(A)
      C   = ABS(B)
      IF (C .EQ. 0.0) THEN
         C = SIGN(ONE,A)
         S = 0.0
      ELSE
         IF (SIG .GE. C) THEN
            SIG = SIG*SQRT(ONE+(C/SIG)**2)
         ELSE
            SIG = C*SQRT(ONE+(SIG/C)**2)
         ENDIF
         C = A/SIG
         S = B/SIG
      ENDIF
      RETURN
      END
