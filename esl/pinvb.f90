!! Purpose    : COMPUTE SOLUTION V TO P*V=B WHERE P=U*D*U**T AND VECTOR B HAS
!!              ZEROS BELOW ELEMENT K.
!!              THE METHOD USED EMPLOYS BACK SUBSTITITION ON THE INPUT
!!              VECTOR STORED U-D FACTORS.
!! Parameters :
!!              U --- VECTOR STORED U-D COVARIANCE FACTORS (SINGULAR D FACTOR
!!                    COMPONENTS ARE ALLOWED PROVIDED CORRESPONDING B COMPONENTS
!!                    AND ROW/COLUMN ENTRIES OF U ARE ALSO ZERO)
!!              N --- MATRIX DIMENSION
!!              B --- B VECTOR INPUTS
!!              K --- DIMENSION OF B
!!              V --- V = P(INVERSE)*B
!!           IERR --- ERROR RETURN FLAG
!!                    0, NORMAL RETURN
!!                    J, ESTIMATE COMPUTATION FAILED BECAUSE DIAGONAL "J" OF
!!                       U ARRAY IS TOO SMALL

SUBROUTINE PINVB(U,N,B,K,V,IERR)

implicit none

real*8    U(*)
integer*4 N
real*8    B(K)
integer*4 K
real*8    V(N)
integer*4 IERR
!local
integer*4 I,J,JJ,IJ,KM1,IT,M,JM1,KP1,IM1
real*8    S,EPS,TMP
DATA EPS/1.0D-307/

IERR = 0
IF(N.EQ.1) THEN
    J=1
    S=U(1)
    TMP=B(1)
    IF(S.EQ.0.0 .AND. TMP.EQ.0.0) THEN
        V(1)=0.0
    ELSE
        IF(S.LE.DABS(TMP)*EPS) then
            IERR=J
            return
        endif
        V(1)=TMP/S
    ENDIF
    RETURN
ENDIF

! STEP-1 SOLVE U*V = B
V(K)=B(K)
IF(K.ne.1) then
    KM1=K-1
    JJ=K*KM1/2
    S = V(K)
    DO I=1,KM1
        IT=JJ+I
        V(I)= B(I) - U(IT)*S
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
        V(J) = TMP/S
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
        S=S - U(IJ)*V(J)
    enddo
    V(I)=S
    IJ=IJ+1
enddo

RETURN
END