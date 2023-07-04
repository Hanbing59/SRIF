!! Purpose    : STABLE  U-D FACTOR RANK 1 UPDATE
!!              (UOUT)*DOUT*(UOUT)**T=(UIN)*DIN*(UIN)**T+C*V*V**T
!! Parameters :
!!               UIN --- INPUT VECTOR STORED POSITIVE SEMI-DEFINITE U-D
!!                       ARRAY, WITH D ELEMENTS STORED ON THE DIAGONAL
!!                 N --- MATRIX DIMENSION, N .GE. 1
!!                 C --- INPUT SCALAR. SHOULD BE NON-NEGATIVE
!!                 V --- INPUT VECTOR FOR RANK ONE MODIFICATION. V IS
!!                       DESTROYED DURING THE PROCESS
!!              UOUT --- OUTPUT VECTOR STORED POSITIVE (POSSIBLY) SEMI-DEFINITE
!!                       U-D RESULT. COMPILER PERMITTING UOUT=UIN IS ALLOWED
!!              IERR --- DIAGNOSTIC FLAG
!!                       0, NORMAL RETURN
!!                       J, ERROR RETURN DUE TO A NEGATIVE COMPUTED J-TH
!!                          DIAGONAL ENTRY
!!                      -1, UD RESULT IS SINGULAR


SUBROUTINE RANK1(UIN,UOUT,N,C,V,IERR)

implicit none

real*8    UIN(*)
real*8    UOUT(*)
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
IF(CJ.le.Z) then
    DO J=1,JJ
        UOUT(J)=UIN(J)
    enddo
    RETURN
endif
IF(N.gt.1) then
    NP2=N+2
    DO L=2,N
        J=NP2-L
        S=V(J)
        BETA=CJ*S
        D=UIN(JJ)+BETA*S
        IF(D.le.EPS) then
            IF(D.lt.Z) then
                IERR = J
                RETURN
            endif
            JJ = JJ-J
            IERR = -1
            DO K=1,J
                UOUT(JJ+K) = Z
            enddo
            cycle
        endif
        BETA = BETA/D
        ALPHA = UIN(JJ)/D
        CJ = ALPHA*CJ
        UOUT(JJ) = D
        JJ = JJ-J
        JM1 = J-1
        IF(ALPHA.ge.TST) then
            DO I=1,JM1
                V(I)  = V(I)-S*UIN(JJ+I)
                UOUT(JJ+I) = BETA*V(I) + UIN(JJ+I)
            enddo
            cycle
        endif
        DO I=1,JM1
            D = V(I) - S*UIN(JJ+I)
            UOUT(JJ+I) = ALPHA*UIN(JJ+I) + BETA*V(I)
            V(I) = D
        enddo
    enddo
endif
UOUT(1) = UIN(1) + (CJ*V(1))*V(1)

RETURN
END