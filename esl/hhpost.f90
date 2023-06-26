!! Purpose    : TRIANGULARIZES RECTANGULAR W BY POST MULTIPLYING IT
!!              BY AN ORTHOGONAL TRANSFORMATION T. THE RESULT IS IN S.
!!
!! Parameters :
!!                W --- INPUT RECTANGULAR SQRT COVARIANCE MATRIX
!!                      (W IS DESTROYED BY COMPUTATIONS).
!!            IMAXW --- ROW DIMENSION OF W.
!!              IRW --- NUMBER OF ROWS OF W TO BE TRIANGULARIZED
!!                      AND THE DIMENSION OF S (IRW.GE.1).
!!              JCW --- NUMBER OF COLUMNS OF W (JCW.GE.IRW).
!!                S --- OUTPUT UPPER TRIANGULAR VECTOR STORED SQRT
!!                      COVARIANCE MATRIX.


SUBROUTINE HHPOST(W,IMAXW,IRW,JCW,V,S)

implicit none

real*8    W(IMAXW,JCW)
integer*4 IMAXW
integer*4 IRW,JCW
real*8    V(JCW)
real*8    S(*)
!local
real*8    SUM,BETA
real*8    ZERO,ONE
integer*4 IROW,JCOL,JC,I,J,L,K
integer*4 NSYM,JCOLM1,IROWM1
real*8    TMPMAX,TMP

ZERO = 0.0D0
ONE = 1.0D0
JCOL = JCW
IF(IRW.ne.1) then
    NSYM = IRW*(IRW+1)/2
    JC = IRW+2
    DO L=2,IRW
        IROW = JC-L
        TMPMAX = ZERO
        DO K = 1,JCOL
            V(K) = W(IROW,K)
            TMP = ABS(V(K))
            IF(TMPMAX.LT.TMP) TMPMAX = TMP
        enddo
        JCOLM1=JCOL-1
        IF(TMPMAX.ne.0.0) then
            SUM = ZERO
            DO K = 1, JCOL
                SUM = SUM + (V(K)/TMPMAX)**2
            enddo
            SUM = TMPMAX*DSQRT(SUM)
            ! DIAGONAL ENTRY (JCOL,JCOL)
            IF(V(JCOL).GT.ZERO) SUM = -SUM
            S(NSYM) = SUM
            NSYM = NSYM-IROW
            DO I  = 1, JCOL
                V(I)  = V(I)/SUM
            enddo
            V(JCOL) = V(JCOL)-1.0
            BETA = ONE/V(JCOL)
            ! T(ORTHOG. TRANS.)=I-BETA*V*V**T, V - SCALED
            IROWM1=IROW-1
            DO I=1,IROWM1
                SUM=ZERO
                DO K=1,JCOL
                    SUM=SUM+V(K)*W(I,K)
                enddo
                S(NSYM+I)=W(I,JCOL)+SUM
                SUM=BETA*SUM
                DO K=1,JCOLM1
                    W(I,K)=W(I,K)+SUM*V(K)
                enddo
            enddo
        endif
        JCOL=JCOLM1
    enddo
endif
JC=JCW-IRW+1
TMPMAX = ZERO
DO J=1,JC
    TMP = ABS(W(1,J))
    IF(TMPMAX.LT.TMP) TMPMAX = TMP
enddo
IF(TMPMAX.eq.ZERO) then
    S(1) = 0.0
    RETURN
endif
SUM = ZERO
DO J = 1, JC
    SUM = SUM + (W(1,J)/TMPMAX)**2
enddo
S(1) = TMPMAX*DSQRT(SUM)

RETURN
END