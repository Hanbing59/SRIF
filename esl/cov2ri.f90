!! Purpose    : CONSTRUCT THE UPPER TRIANGULAR CHOLESKY FACTOR OF A POSITIVE 
!!              SEMI-DEFINITE MATRIX. BOTH THE INPUT COVARIANCE AND THE OUTPUT 
!!              CHOLESKY FACTOR (SQUARE ROOT) ARE VECTOR STORED. THE OUTPUT 
!!              OVERWRITES THE INPUT. IF THE INPUT COVARIANCE IS SINGULAR OR 
!!              INDEFINITE THE OUTPUT FACTOR HAS ZERO COLUMNS.
!!
!!              COVARIANCE (INPUT) = U*U**T (U IS OUTPUT)
!!
!! Parameters : 
!!               U --- CONTAINS THE INPUT VECTOR STORED COVARIANCE MATRIX 
!!                     (ASSUMED POSITIVEDEFINITE)
!!                     For output, THE UPPER TRIANGULAR SQUARE ROOT CHOLESKY FACTOR
!!               N --- STATE DIMENSION    


SUBROUTINE COV2RI(U,N)

implicit none

real*8    U(*)
integer*4 N
!local
real*8    ALPHA
real*8    S
integer*4 JJ,KK,JJN,JM1,NP2,L,J,K,JJNPK,I,KKPI,JJNPI
real*8    EPS
real*8    ZERO
DATA EPS /1.0D-307/
DATA ZERO /0.0D0/

IF(N.ne.1) then
    JJ = N*(N+1)/2
    JJN = JJ
    NP2 = N + 2
    DO L = 2, N
        J = NP2 - L
        S = U(JJ)
        IF(S.LT.ZERO) S = ZERO
        S = DSQRT(S)
        U(JJ) = S
        ALPHA = ZERO
        IF(S.GT.EPS) ALPHA = 1.0/S
        KK = 0
        JJN = JJ - J
        JM1 = J - 1
        DO K=1,JM1
            JJNPK = JJN + K
            U(JJNPK) = ALPHA*U(JJNPK)
            S = U(JJNPK)
            DO I = 1, K
                KKPI = KK + I
                JJNPI= JJN + I
                U(KKPI) = U(KKPI) - S*U(JJNPI)
            enddo
            KK = KK + K
        enddo
        JJ = JJN
    enddo
endif
S = U(1)
IF(S.LT.ZERO) S = ZERO
U(1) = DSQRT(S)
      
RETURN
END