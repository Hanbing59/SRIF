      SUBROUTINE HHPOST(W,IMAXW,IRW,JCW,V,S)
C
C         TRIANGULARIZES RECTANGULAR W BY POST MULTIPLYING
C         IT BY AN ORTHOGONAL TRANSFORMATION T. THE RESULT IS IN S.
C
C *** INPUT:
C
C     W(IRW,JCW)       INPUT RECTANGULAR SQRT COVARIANCE MATRIX
C                      (W IS DESTROYED BY COMPUTATIONS).
C
C     IMAXW            ROW DIMENSION OF W.
C
C     IRW              NUMBER OF ROWS OF W TO BE TRIANGULARIZED
C                      AND THE DIMENSION OF S (IRW.GE.1).
C
C     JCW              NUMBER OF COLUMNS OF W (JCW.GE.IRW).
C
C *** OUTPUT:
C
C     S(IRW*(IRW+1)/2) OUTPUT UPPER TRIANGULAR VECTOR STORED SQRT
C                      COVARIANCE MATRIX.
C
C *** SCRATCH:
C
C     V(JCW)            WORK VECTOR.
C
C
C        REQUIRED SUBROUTINES:
C
C                 ABS      <FORLIB>
C                 DSQRT    <FORLIB>
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              10 OCTOBER 1981
C
C              REVISED 31 JANUARY 1982   VERSION 1.00.00
C              REVISED 2  MARCH   1983   VERSION 1.01.00
C              REVISED 5  MAY     1984   VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE HHPOST IS A PART OF THE FEA. INC             *
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

      DOUBLE PRECISION  SUM,BETA
      DIMENSION  S(*),W(IMAXW,JCW),V(JCW)
C
      ZERO = 0.0D0
      ONE  = 1.0D0
C
      JCOL = JCW
      IF (IRW .EQ. 1) GO TO 90
      NSYM = IRW*(IRW+1)/2
      JC   = IRW+2
C
      DO 80 L=2,IRW
        IROW   = JC-L
        TMPMAX = ZERO
C
        DO 10 K = 1,JCOL
           V(K)  = W(IROW,K)
           TMP   = ABS(V(K))
           IF (TMPMAX .LT. TMP) TMPMAX = TMP
  10    CONTINUE
        JCOLM1=JCOL-1
        IF (TMPMAX .EQ. 0.0) GO TO 70
        SUM = ZERO
        DO 20 K = 1, JCOL
           SUM = SUM + (V(K)/TMPMAX)**2
  20    CONTINUE
C
        SUM = TMPMAX*DSQRT(SUM)
C                                      DIAGONAL ENTRY (JCOL,JCOL)
C
        IF (V(JCOL) .GT. ZERO) SUM = -SUM
        S(NSYM)  = SUM
        NSYM     = NSYM-IROW
        DO 30 I  = 1, JCOL
           V(I)  = V(I)/SUM
  30    CONTINUE
        V(JCOL)  = V(JCOL)-1.0
        BETA     = ONE/V(JCOL)
C
C         T(ORTHOG. TRANS.)=I-BETA*V*V**T, V - SCALED
C
        IROWM1=IROW-1
C
        DO 60 I=1,IROWM1
           SUM=ZERO
           DO 40 K=1,JCOL
              SUM=SUM+V(K)*W(I,K)
  40       CONTINUE
C
           S(NSYM+I)=W(I,JCOL)+SUM
           SUM=BETA*SUM
           DO 50 K=1,JCOLM1
              W(I,K)=W(I,K)+SUM*V(K)
  50       CONTINUE
  60    CONTINUE
C
  70    CONTINUE
        JCOL=JCOLM1
  80  CONTINUE
C
  90  CONTINUE
      JC=JCW-IRW+1
      TMPMAX = ZERO
      DO 100 J=1,JC
         TMP = ABS(W(1,J))
         IF (TMPMAX .LT. TMP) TMPMAX = TMP
 100  CONTINUE
      IF (TMPMAX .NE. ZERO) GO TO 110
         S(1) = 0.0
         RETURN
 110  CONTINUE
      SUM = ZERO
      DO 120 J = 1, JC
         SUM = SUM + (W(1,J)/TMPMAX)**2
 120  CONTINUE
      S(1) = TMPMAX*DSQRT(SUM)
      RETURN
      END
