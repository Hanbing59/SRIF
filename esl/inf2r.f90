!! Purpose    : CHOLESKY FACTOR AN INFORMATION MATRIX.
!!              COMPUTES A LOWER TRIANGULAR VECTOR STORED CHOLESKY
!!              FACTORIZATION OF A POSITIVE SEMI-DEFINITE MATRIX.
!!              R=R(**T)R, R UPPER TRIANGULAR. BOTH MATRICES ARE VECTOR
!!              STORED  AND THE RESULT OVERWRITES THE INPUT
!!
!! Parameters :
!!                R --- ON INPUT THIS IS A POSITIVE SEMI-DEFINITE
!!                      (INFORMATION) MATRIX.
!!                      ON OUTPUT IT IS THE TRANSPOSED LOWER TRIANGULAR
!!                      CHOLESKY FACTOR. IF THE INPUT MATRIX IS SINGULAR,
!!                      THE OUTPUT MATRIX WILL HAVE ZERO DIAGONAL ENTRIES
!!                N --- DIMENSION OF MATRICES INVOLVED, N.GE.1
!!             IERR --- ERROR RETURN FLAG
!!                      0, NORMAL RETURN
!!                      L, COMPUTED NEGATIVE DIAGONALS
!!                         THE ALGORITHM ASSUMES THAT THE USER INTENDED THE
!!                         INPUT TO BE AT LEAST POSITIVE SEMI-DEFINITE, AND
!!                         THAT NEGATIVE VALUES CORRESPOND TO COMPUTER ROUNDOFF.
!!                         THE NEGATIVE DIAGONAL COMPUTED VALUES ARE DISPLAYED
!!                         AND THE CORRESPONDING COLUMN (INCLUDING THE DIAGONAL)
!!                         IS SET TO ZERO.
!!           IOUNIT --- UNIT NUMBER OF OUTPUT DEVICE

SUBROUTINE INF2R(R,N,IERR,IOUNIT)

implicit none

real*8    R(*)
integer*4 N
integer*4 IERR
integer*4 IOUNIT
!local
integer*4 I,J,L,K,M
integer*4 JJ,JK,NN,NM1,JIS,NPJP1,JP1,KI,JI,NPK
real*8    SUM,EPS,Z,ONE,ALPHA,BETA
DATA Z/0.D0/, ONE /1.0D0/
DATA EPS/1.0D-307/

100 FORMAT (1H0,20X,' AT STEP  ',I4,5X,' DIAGONAL ENTRY = ',1P1D14.7,2X,' IT IS RESET TO ZERO ')
    
IERR = 0
JJ=0
NN=N*(N+1)/2
IF(N.ne.1) then
    NM1=N-1
    DO J=1,NM1
        ! JJ=(J,J)
        JJ=JJ+J
        SUM = R(JJ)
        IF(SUM.lt.Z)  then
            WRITE(IOUNIT,100) J,SUM
            IERR = IERR + 1
            R(JJ)=Z
        endif
        SUM = DSQRT(SUM)
        ALPHA=Z
        IF(SUM.GT.EPS) ALPHA=ONE/SUM
        ! JK=(J,K)
        JK = NN+J
        JP1 = J+1
        ! JIS=(J,I) START
        JIS = JK
        NPJP1 = N+JP1
        DO L=JP1,N
            K = NPJP1-L
            JK = JK-K
            R(JK) = ALPHA*R(JK)
            BETA = R(JK)
            KI = NN+K
            JI = JIS
            NPK = N+K
            DO M=K,N
                I = NPK-M
                KI = KI-I
                JI = JI-I
                R(KI) = R(KI) - R(JI)*BETA
            enddo
        enddo
    enddo
endif
SUM = R(NN)
IF(SUM.lt.Z) then
    WRITE (IOUNIT,100) N,SUM
    R(NN) = Z
    IERR  = IERR+1
endif
R(NN)= DSQRT(SUM)

RETURN
END