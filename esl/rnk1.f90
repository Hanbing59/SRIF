!! Purpose    : STABLE  U-D FACTOR RANK 1 UPDATE
!!
!! Parameters :
!!              U --- INPUT VECTOR STORED POSITIVE SEMI-DEFINITE U-D
!!                    ARRAY, WITH D ELEMENTS STORED ON THE DIAGONAL
!!                    OUTPUT(U) OVERWRITES INPUT.
!!              N --- MATRIX DIMENSION, N.GE.1
!!              C --- INPUT SCALAR. SHOULD BE NON-NEGATIVE
!!              V --- INPUT VECTOR FOR RANK ONE MODIFICATION. V IS 
!!                    DESTROYED DURING THE PROCESS
!!           IERR --- DIAGNOSTIC FLAG
!!                    0, NORMAL RETURN
!!                    J, ERROR RETURN DUE TO A NEGATIVE COMPUTED J-TH
!!                       DIAGONAL ENTRY
!!                   -1, UD RESULT IS SINGULAR

SUBROUTINE RNK1 (U,N,C,V,IERR)

implicit none

real*8    U(*)
integer*4 N
real*8    C
real*8    V(N)
integer*4 IERR
!local
integer*4 J,JJ,NP2,L,K,JM1,I
real*8    Z,CJ
real*8    ALPHA,BETA,S,D,EPS,TST
DATA EPS/1.0D-307/
! TST=1/16 IS USED FOR RANK1 ALGORITHM SWITCHING
DATA TST/.0625D0/

Z = 0.D0
IERR = 0
JJ = N*(N+1)/2
CJ = C
IF(CJ.LE.Z) RETURN
IF(N.ne.1) then
    NP2 =  N+2
    DO L=2,N
        J =  NP2-L
        S =  V(J)
        BETA =  CJ*S
        D =  U(JJ)+BETA*S
        IF(D.le.EPS) then
            IF(D.lt.Z) then
                IERR = J
                RETURN
            endif
            JJ = JJ-J
            IERR = -1
            DO K=1,J
                U(JJ+K)=Z
            enddo
            cycle
        endif
        BETA = BETA/D
        ALPHA = U(JJ)/D
        CJ = ALPHA*CJ
        U(JJ) = D
        JJ = JJ-J
        JM1 = J-1
        IF(ALPHA.ge.TST) then
            DO I=1,JM1
                V(I)=V(I)-S*U(JJ+I)
                U(JJ+I)=BETA*V(I)+U(JJ+I)
            enddo
            cycle
        endif
        DO I=1,JM1
            D = V(I)-S*U(JJ+I)
            U(JJ+I) = ALPHA*U(JJ+I)+BETA*V(I)
            V(I) =  D
        enddo
    enddo
endif
! CJ*V(1)**2
U(1) = U(1) + (CJ*V(1))*V(1)

RETURN
END