!! Purpose    : TO OBTAIN THE U-D FACTORS OF A POSITIVE SEMI-DEFINITE MATRIX. 
!!              THE INPUT VECTOR STORED MATRIX IS OVERWRITTEN BY THE OUTPUT 
!!              U-D FACTORS WHICH ARE ALSO VECTOR STORED.    
!!
!! Parameters :
!!                U --- VECTOR STORED COVARIANCE (POSITIVE SEMI-DEFINITE) MATRIX
!!                      For output, VECTOR STORED U-D COVARIANCE FACTORS
!!                N --- MATRIX DIMENSION, N .GE. 1
!!             IERR --- ERROR RETURN FLAG
!!                      0, NORMAL RETURN
!!                      L, COMPUTED NEGATIVE DIAGONALS
!!                         THE ALGORITHM ASSUMES THAT THE USER INTENDED THE INPUT
!!                         TO BE AT LEAST POSITIVE SEMI-DEFINITE, AND THAT NEGATIVE 
!!                         VALUES CORRESPOND TO COMPUTER ROUNDOFF. THE NEGATIVE 
!!                         DIAGONAL COMPUTED VALUES ARE DISPLAYED AND THE 
!!                         CORRESPONDING COLUMN (INCLUDING THE DIAGONAL) IS SET TO 
!!                         ZERO.
!!           IOUNIT --- UNIT NUMBER OF OUTPUT DEVICE

SUBROUTINE COV2UD (U,N,IERR,IOUNIT)

implicit none

real*8    U(*)
integer*4 N
integer*4 IERR
integer*4 IOUNIT
!local
integer*4 I,JJ,NP2,L,J,K,KK,KJ,JM1,IJ,IK
real*8    S,EPS,Z,ONE,ALPHA,BETA
DATA Z /0.0D0/, ONE /1.0D0/
DATA EPS/1.0D-307/

100 FORMAT (1H0,20X,' AT STEP ',I4,' DIAGONAL ENTRY = ',1PD14.7)
    
IERR = 0
IF(N.ne.1) then
    JJ=N*(N+1)/2
    NP2=N+2
    DO L=2,N
        J=NP2-L
        ALPHA=Z
        S = U(JJ)
        IF(S.lt.Z) then
            IF(IOUNIT.NE.0) WRITE(IOUNIT,100) J,S
            IERR = IERR + 1
            U(JJ)=Z
        endif
        IF(S.GT.EPS) ALPHA=ONE/S
        JJ=JJ-J
        KK=0
        KJ=JJ
        JM1=J-1
        DO K=1,JM1
            KJ=KJ+1
            BETA=U(KJ)
            U(KJ)=ALPHA*U(KJ)
            IJ=JJ
            IK=KK
            DO I=1,K
                IK=IK+1
                IJ=IJ+1
                U(IK)=U(IK)-BETA*U(IJ)
            enddo
            KK=KK+K
        enddo
    enddo
endif
S = U(1)
IF(S.lt.Z) then
    J = 1
    IF(IOUNIT.NE.0) WRITE(IOUNIT,100) J,S
    IERR = IERR+1
    U(1)=Z
endif

RETURN
END