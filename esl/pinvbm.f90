!! Purpose    : COMPUTE SOLUTION V TO P*V=BM WHERE P=U*D*U**T and VECTOR B HAS
!!              ZEROS BELOW ELEMENT NX.
!!              BM(J)=EM*B(J),J=1 TO NX, BM HAS AN IMPLICIT UNIT VALUE IN
!!              LOCATION K, AND THIS CORRESPONDENDS TO THE LOCATION OF THE
!!              COLORED NOISE MULTIPLIER, EM.
!!              THe METHOD USED EMPLOYS BACK SUBSTITUTION ON THE INPUT
!!              VECTOR STORED U-D FACTORS.
!! Parameters :
!!              U --- VECTOR STORED U-D COVARIANCE FACTORS. SINGULAR D FACTORS
!!                    ARE ALLOWED PROVIDED CORRESPONDING ROW & COLUMN ENTRIES
!!                    OF U ARE ZERO
!!              N --- MATRIX DIMENSION
!!              B --- B VECTOR INPUTS
!!             NX --- DIMENSION OF NON TRIVIAL PORTION OF B
!!              K --- LOCATION OF COLORED PROCESS NOISE STATE IN FILTER STATE
!!                    VECTOR, NX .LT. K .LE. N
!!             EM --- EXPONENTIAL (COLORED NOISE) MULTIPLIER
!!              Q --- COLORED NOISE VARIABLE PROCESS NOISE VARIANCE
!!           IEST --- TRUE, COMPUTE BETA = (V**T)*XEST
!!                    Otherwise, INPUT XEST IS NOT USED
!!           XEST --- FILTER ESTIMATE STATE VECTOR
!!              V --- V(1:N) = P(INVERSE)*BM
!!                    V(N+1) = (V**T)*XEST, WHEN IEST = TRUE UNCHANGED OTHERWISE
!!           IERR --- ERROR RETURN FLAG
!!                    0, NORMAL RETURN
!!                    J, ESTIMATE COMPUTATION FAILED BECAUSE DIAGONAL "J" OF
!!                       U ARRAY IS TOO SMALL. BETA NOT COMPUTED



SUBROUTINE PINVBM (U,N,B,NX,K,EM,Q,IEST,XEST,V,IERR)

implicit none

real*8    U(*)
integer*4 N
real*8    B(*)
integer*4 NX,K
real*8    EM,Q
LOGICAL   IEST
real*8    XEST(N)
real*8    V(*)
integer*4 IERR
!local
integer*4 I,J,L,JJ,JK,KP1,IJ,IM1
real*8    S,EPS,EPSCN,TMP,XMONE
real*8    SIG
DATA EPS/1.0D-307/
DATA EPSCN/1.0D-7/
EXTERNAL  U2SIG

! WHEN THERE IS NO PROCESS NOISE, MODIFIED GAIN IS ESSENTIALLY ZERO
XMONE=-1.0
! EPSCN IS USED TO TEST WHETHER THE COLORED NOISE PARAMETER IS A CONSTANT.
IF(EM.GT.(1.0D0-EPSCN)) THEN
    SIG = 1.0
    IF(Q.GT.0.0) CALL U2SIG(U,N,' ',' ',XMONE,' ',.FALSE.,K,1,SIG,' ',N,0)
    IF(Q.EQ.0.0 .OR. Q.LT.(EPSCN*SIG)**2) THEN
        DO J = 1, N+1
            V(J) = 0.0
        enddo
        RETURN
    ENDIF
ENDIF
! IF N IS ONE SKIP THE BULK OF THE CODE
IERR = 0
IF(N.EQ.1) THEN
    J = 1
    S = U(1)
    IF(S.LE.EPS) then
        IERR=J
        return
    endif
    V(1) = 1.0/S
    V(N+1) = XEST(1)*V(1)
    RETURN
ENDIF

! STEP-1 SOLVE U*V = BM
IF(EM.EQ.0.0 .OR. NX.EQ.0) THEN
    DO J=1,K
        V(J) = 0.0
    enddo
ELSE
    DO J=1,NX
        V(J) = EM*B(J)
    enddo
    DO J = NX+1,K
        V(J) = 0.0
    enddo
ENDIF
V(K) = 1.0
IF(K.GT.1) THEN
    JK = K*(K+1)/2
    DO L=K,2,-1
        S = V(L)
        DO J=L-1,1,-1
            JK = JK-1
            V(J) = V(J)-U(JK)*S
        enddo
        JK = JK-1
    enddo
ENDIF

! STEP-2 V=D**(-1)*V
JJ = 0
DO J=1,K
    JJ = JJ+J
    S = U(JJ)
    TMP = V(J)
    IF(S.GT.0.0) THEN
        IF(S.LE.EPS*DABS(TMP)) then
            IERR=J
            return
        endif
        V(J) = TMP/S
    ELSE
        V(J) = 0.0
    ENDIF
enddo

! STEP-3 V=U**(-T)*V
IF(K.LT.N) THEN
    KP1 = K+1
    DO I=KP1,N
        V(I) = 0.0
    enddo
ENDIF
IJ = 1
DO I=2,N
    S = V(I)
    IM1 = I-1
    DO J=1,IM1
        IJ = IJ+1
        S = S-U(IJ)*V(J)
    enddo
    V(I) = S
    IJ = IJ+1
enddo
IF(IEST) THEN
    S = 0.0
    DO L=1,N
        S = S+V(L)*XEST(L)
    enddo
    V(N+1) = S
ENDIF

RETURN
END