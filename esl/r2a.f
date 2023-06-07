      SUBROUTINE R2A(R,N,NAMR,A,IMAXA,JA,NAMA)
C
C
C     TO PLACE THE TRIANGULAR VECTOR STORED MATRIX "R"
C     INTO THE MATRIX "A" AND TO ARRANGE THE COLUMNS
C     TO MATCH THE DESIRED NAMA PARAMETER LIST.
C     NAMES IN THE NAMA LIST THAT DO NOT CORRESPOND
C     TO ANY NAME IN THE NAMR LIST HAVE ZERO ENTRIES
C     IN THE CORRESPONDING "A" COLUMN
C
C *** INPUTS:
C
C     R(N*(N+1)/2)            INPUT UPPER TRIANGULAR VECTOR
C                             STORED ARRAY
C
C     N                       DIMENSION OF R
C
C     NAMR(N)                 PARAMETER NAMES ASSOCIATED WITH R
C
C *** OUTPUTS:
C
C     A(N,JA)                 MATRIX TO HOUSE THE REARRANGED
C                             R MATRIX
C
C     IMAXA                   ROW DIMENSION OF A
C
C     JA                      NUMBER OF PARAMETER NAMES ASSOCIATED
C                             WITH THE OUTPUT "A" MATRIX
C
C     NAMA(JA)                NAMES ASSOCIATED WITH "A"
C
C
C          COGNIZANT PERSONS:
C
C             DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C             FACTORIZED ESTIMATION APPLICATIONS INC.
C
C             01 SEPTEMBER  1976
C
C             REVISED 28 MARCH   1982
C             REVISED  2 AUGUST  1982
C             REVISED 20 AUGUST  1982   VERSION 0.00.02
C
C             REVISED 20      MAY 1985  VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE R2A    IS A PART OF THE FEA. INC             *
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

      DIMENSION R(*), A(IMAXA,JA)
      CHARACTER*(*) NAMR(N), NAMA(JA)
      CHARACTER*17 NAM
C
      DO 10 J = 1, JA
         DO 9 I =1, N
            A(I,J) = 0.0
   9     CONTINUE
  10  CONTINUE
      DO 40 J=1, JA
         NAM = NAMA(J)
         DO 15 I=1,N
            IF ( NAMR(I) .EQ. NAM)  GO TO 20
  15     CONTINUE
            GO TO 40
  20     JJ = I*(I-1)/2
         DO 30 K=1,I
            JJPK = JJ+K
            A(K,J) = R(JJPK)
  30     CONTINUE
  40  CONTINUE
      RETURN
      END
