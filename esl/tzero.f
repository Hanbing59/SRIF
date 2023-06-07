      SUBROUTINE  TZERO (R,N,ISTRT,ILST)
C
C        TO ZERO OUT ROWS IS (ISTART) TO IF (IFINAL) OF A VECTOR
C        STORED UPPER TRIANGULAR MATRIX
C
C        R(N*(N+1)/2)  INPUT VECTOR STORED UPPER TRIANGULAR MATRIX
C
C        N             DIMENSION OF R
C
C        ISTRT         FIRST ROW OF R THAT IS TO BE SET TO ZERO
C
C        ILST          LAST ROW OF R THAT IS TO BE SET TO ZERO
C
C     SUBROUTINES REQUIRED:
C
C     MAX0   <FORLIB>           MIN0   <FORLIB>
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23, 1981
C
C              REVISED  9 AUGUST 1982     VERSION 0.01.00
C                                         VERSION 1.00.00
C                      22 DECEMBER 1982   VERSION 1.01.00
C              REVISED 20 MAY      1985   VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE TZERO  IS A PART OF THE FEA. INC             *
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
C
      ZERO   = 0.D0
      IFINAL = MAX0(ISTRT,ILST)
      ISTART = MIN0(ISTRT,ILST)
      IJS    = ISTART*(ISTART-1)/2
      DO 10 I=ISTART, IFINAL
         IJS = IJS+I
         IJ  = IJS
         DO 10 J=I,N
            R(IJ) = ZERO
            IJ    = IJ+J
   10 CONTINUE
      RETURN
      END
