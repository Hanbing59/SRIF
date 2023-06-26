!! Purpose    : COMPUTE SOLUTION V TO P*V=E(SUB K) WHERE P=U*D*U**T AND
!!              E(SUB K) HAS UNITY IN ELEMENT K, I.E., TO COMPUTE COLUMN
!!              K OF P INVERSE.
!!              THE METHOD USED EMPLOYS BACK SUBSTITITION ON THE INPUT
!!              VECTOR STORED U-D FACTORS.
!! Parameters :
!!              U --- VECTOR STORED U-D COVARIANCE FACTORS (MUST BE
!!                    NONSINGULAR)
!!              N --- MATRIX DIMENSION
!!              K --- LOCATION OF THE UNIT ELEMENT IN E(SUB K)
!!              V --- COLUMN K OF P INVERSE
!!           IERR --- ERROR RETURN FLAG
!!                    0, NORMAL RETURN
!!                    J, ESTIMATE COMPUTATION FAILED BECAUSE
!!                       DIAGONAL "J" OF U ARRAY IS TOO SMALL

SUBROUTINE PINVEK(U,N,K,V,IERR)

IMPLICIT none

real*8    U(*)
integer*4 N,K
real*8    V(N)
integer*4 IERR
!local
integer*4 I,J,M,JJ,IJ,IM1,JM1,KM1,IT,KP1
real*8    S,EPS,TMP
DATA EPS/1.0D-307/

IERR = 0
IF(N.le.1) then
    J = 1
    S = U(1)
    IF(S.LE.EPS) then
        IERR=J
        return
    endif
    V(1)=1.0/S
    RETURN
endif

! STEP-1 SOLVE U*V = E(SUB K)
V(K)=1.0
IF(K.ne.1) then
    KM1=K-1
    JJ=K*KM1/2
    DO I=1,KM1
        IT=JJ+I
        V(I)= -U(IT)
    enddo
    IF(K.ne.2) then
        DO M=2,KM1
            JM1=K-M
            J=JM1+1
            JJ=JJ-J
            S=V(J)
            DO I=1,JM1
                IT=JJ+I
                V(I)=V(I)-U(IT)*S
            enddo
        enddo
    endif
endif

! STEP-2 V=D**(-1)*V
JJ=0
DO J=1,K
    JJ = JJ+J
    S = U(JJ)
    TMP = V(J)
    IF(S.EQ.0.0 .AND. TMP.EQ.0.0) THEN
        V(J) = 0.0
    ELSE
        IF(S.LE.EPS*DABS(TMP)) then
            IERR=J
            return
        endif
        V(J)= TMP/S
    ENDIF
enddo

! STEP-3 V=U**(-T)*V
IF(K.lt.N) then
    KP1=K+1
    DO I=KP1,N
        V(I)=0.0
    enddo
endif

IJ=1
DO I=2,N
    S=V(I)
    IM1=I-1
    DO J=1,IM1
        IJ=IJ+1
        S=S-U(IJ)*V(J)
    enddo
    V(I)=S
    IJ=IJ+1
enddo

RETURN
END