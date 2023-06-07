      SUBROUTINE TDHHT(S,IMAXS,IRS,JCS,JSTART,JSTOP,V)
C
C       TDHHT TRANSFORMS A RECTANGULAR DOUBLE SUBSCRIPTED MATRIX S
C       TO AN UPPER TRIANGULAR OR PARTIALLY UPPER TRIANGULAR FORM
C       BY THE APPLICATION OF HOUSEHOLDER ORTHOGONAL
C       TRANSFORMATIONS. IT IS ASSUMED THAT THE FIRST "JSTART"-1
C       COLUMNS OF S ARE ALREADY TRIANGULARIZED. THE ALGORITHM IS
C       DESCRIBED IN "FACTORIZATION METHODS FOR DISCRETE SEQUENTIAL
C        ESTIMATION" BY G.J.BIERMAN, ACADEMIC PRESS, 1977
C
C *** INPUTS:
C
C     S(IRS,JCS)  INPUT (POSSIBLY PARTIALLY) TRIANGULAR MATRIX. THE
C                 OUTPUT (POSSIBLY PARTIALLY) TRIANGULAR RESULT
C                 OVERWRITES THE INPUT.
C
C     IMAXS       ROW DIMENSION OF S
C
C     IRS         NUMBER OF ROWS IN S  (IRS.LE.IMAXS.AND.IRS.GE.1)
C
C     JCS         NUMBER OF COLUMNS IN S
C
C     JSTART      INDEX OF THE FIRST COLUMN TO BE TRIANGULARIZED. IF
C                 JSTART.LT.1 IT IS ASSUMED THAT JSTART=1, I.E.
C                 START TRIANGULARIZATION AT COLUMN 1.
C
C     JSTOP       INDEX OF LAST COLUMN TO BE TRIANGULARIZED.
C                 IF JSTOP.LT.JSTART.OR.JSTOP.GT.JCS THEN
C                   IF IRS.LE.JCS JSTOP IS SET EQUAL TO IRS-1
C                   IF IRS.GT.JCS JSTOP IS SET EQUAL TO JCS
C                 I.E. THE TRIANGULARIZATION IS COMPLETED AS FAR
C                 AS POSSIBLE
C
C     V(IRS)      WORK VECTOR
C
C  *** REQUIRED SUBROUTINES:
C
C      ABS     <FORLIB>
C      DABS    <FORLIB>
C      DSQRT   <FORLIB>
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              18 OCTOBER 1978
C              REVISED 29 MAY 1982        VERSION 0.01.00
C                                         VERSION 1.00.00
C              REVISED 05 MAY 1984        VERSION 2.00.00
C              REVISED 16 MAY 1986        VERSION 2.01.00
C              REVISED 17 OCT 1986        VERSION 2.02.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE TDHHT  IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *           COPYRIGHT 1982, 1983, 1984, 1985, 1986               *
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

      DIMENSION S(IMAXS,JCS), V(IRS)
      DOUBLE PRECISION SUM, DELTA,DEPS
C
      IF (IRS .LE. 1) RETURN
C
C     MACHINE EPSILON SUITABLE FOR VAX APPLICATIONS
C
C      EPS  = 0.6E-37 !suitable for single precision
      DATA EPS/1.0D-307/
      DEPS = EPS
C
      ONE  = 1.D0
      ZERO = 0.0D0
      JSTT = JSTART
      JSTP = JSTOP
      IF (JSTT .LT. 1) JSTT=1
      IF (JSTP .GE. JSTT .AND. (JSTP .LE. JCS)
     1                   .AND. (JSTP .LE. IRS)) GO TO 5
      IF (IRS .LE. JCS) JSTP=IRS-1
      IF (IRS .GT. JCS) JSTP=JCS
C
   5  DO 60 J=JSTT,JSTP
         TMPMAX = ZERO
         JP1    = J+1
         IF (JP1 .GT. IRS) RETURN
         DO 10 I=JP1,IRS
            TMP = ABS(S(I,J))
C
C       TMPMAX = MAX(S(I,J); I=(J+1),IRS)
C
            IF (TMPMAX .LT. TMP) TMPMAX = TMP
 10      CONTINUE
         IF (TMPMAX .LE. EPS) GO TO 60
C
C         IF TMPMAX = 0, COLUMN J ELEMENTS BELOW J-TH ARE
C         ZERO AND THIS STEP OF THE ALGORITHM IS OMITTED
C
         TMP = ABS(S(J,J))
         IF(TMP .GT. TMPMAX) TMPMAX = TMP
         SUM = ZERO
         DO 20 I = J, IRS
            V(I)   = S(I,J)
            S(I,J) = ZERO
            SUM    = SUM+(V(I)/TMPMAX)**2
  20     CONTINUE
         SUM = TMPMAX * DSQRT(SUM)
C
C       SUM = 2 NORM OF COLUMN J
C
         IF (V(J) .GT. ZERO) SUM = -SUM
         S(J,J)  = SUM
C
C     NORMALIZE HOUSEHOLDER VECTOR, V
C
         DO  30 I = J, IRS
            V(I) = V(I)/SUM
  30     CONTINUE
         V(J) = V(J) - ONE
         SUM  = V(J)
C
C         THE HOUSEHOLDER TRANSFORMATION IS T=I-(1/SUM)*V*V**T,
C         WHERE V IS NORM SCALED.
C
         IF (JP1.GT.JCS) GO TO 60
         DO 50  K = JP1,JCS
            DELTA = ZERO
            DO 40 I = J,IRS
               DELTA=DELTA+S(I,K)*V(I)
  40        CONTINUE
            IF (DABS(DELTA) .LE. DEPS) GO TO 50
            S(J,K)  = S(J,K) + DELTA
            DELTA   = DELTA/SUM
            DO 45 I = JP1,IRS
               S(I,K)  = S(I,K)+DELTA*V(I)
  45        CONTINUE
  50     CONTINUE
  60    CONTINUE
      RETURN
      END
