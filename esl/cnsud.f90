!! Purpose    : COLORED NOISE SMOOTH U-D COVARIANCE FACTOR UPDATE TO ACCOUNT
!!              FOR POSSIBLE OFF DIAGONAL PROCESS NOISE COMPONENTS.
!!              DYNAMIC MODEL:
!!                 !X   !     !I   0   0! !X !     !BX!
!!                 ! I+1!     !         ! ! I!     !  !
!!                 !P   !  =  !0  EM   0! !P !  +  ! 1! W
!!                 ! I+1!     !         ! ! I!     !  !  I
!!                 !Y   !     !0   0   I! !Y !     ! 0!
!!
!! Parameters :
!!              U --- VECTOR STORED U-D SMOOTH COVARIANCE FOR ESTIMATE AT
!!                    STAGE I+1.
!!                    For output, UPDATED U-D FACTOR SMOOTH COVARIANCE FACTORS
!!                    CORRESPONDING TO THE STATE AT TIME 'I'
!!              N --- FILTER MODEL DIMENSION
!!             BX --- STATE DEPENDENT PROCESS NOISE MULTIPLIER
!!             NX --- ONLY THE FIRST NX VARIABLES ARE ASSUMED TO BE DIRECTLY
!!                    AFFECTED BY THE NOISE COMPONENT. IF NX=0, THE MODEL IS
!!                    SIMPLY COLORED NOISE AND BX IS IGNORED. (I.E. IT IS NOT USED)
!!             EM --- COLORED NOISE MULTIPLIER
!!              Q --- PROCESS NOISE VARIANCE
!!             IP --- LOCATION OF THE COLORED NOISE COMPONENT,
!!                    IP.GE.1  (ELSE RETURN)
!!                    IP.GT.NX (ELSE PROGRAM TERMINATES, AS THIS IS A
!!                    VIOLATION OF THE MATH MODEL)
!!             VM --- SMOOTH GAIN INPUT (SMOOTH GAINS ARE MANUFACTURED FROM THIS)
!!         IOUNIT --- UNIT NUMBER FOR FAILURE MODE
!!
!! Reference  : G. BIERMAN - MODIFICATION OF THE BIERMAN COLORED NOISE
!!                           SMOOTHER TO ACCOUNT FOR A NON-DIAGONAL PROCESS
!!                           NOISE COVARIANCE B&A 860429 (not publicly found)

SUBROUTINE CNSUD(U,N,BX,NX,EM,Q,IP,VM,F,B,G,IOUNIT)

implicit none

real*8       U(*)
integer*4    N
real*8       BX(*)
integer*4    NX
real*8       EM,Q
integer*4    IP
real*8       VM(N)
real*8       F(N),B(N),G(N)
integer*4    IOUNIT
!local
real*8       TMP,S,SUM,C,D,THETA,ALPHA
integer*4    I,J,JJ,K,KL,KP1
integer*4    IERR
real*8       GAMMA,EPS
CHARACTER*80 ERRMSG
DATA EPS /1.0D-307/
DATA ERRMSG /'(10X,20(''*''),''CNSUD ERROR'',20(''*''))'/
EXTERNAL UDOPT,RNK1,UCON

IF(IP.LT.1) RETURN
IF(IP.LE.NX) THEN
    WRITE(IOUNIT,ERRMSG)
    WRITE(IOUNIT,*) ' IP .LE. NX, IP = ',IP,' NX = ',NX
    WRITE(IOUNIT,*) ' THIS IS A VIOLATION OF THE MATH MODEL'
    WRITE(IOUNIT,*) ' EXECUTION TERMINATED'
    STOP
ENDIF
TMP = Q*VM(IP)
GAMMA = TMP + EM**2
IF(GAMMA.LE.EPS) THEN
    WRITE(IOUNIT,ERRMSG)
    WRITE(IOUNIT,'(10X,'' ERROR INPUT TO SMOOTH SUBROUTINE: Q ='','//&
        '1PD14.7,'' EM ='',D14.7,''  SMOOTH GAIN('',I4,'') ='','//&
        'D14.7/40X,10(''*''),'' SMOOTH STEP SKIPPED '',10(''*''))') &
        Q,EM,IP,F(IP)
    RETURN
ENDIF
IF(VM(IP).EQ.0.0d0) THEN
    KL = IP*(IP-1)/2
    DO J = 1, IP
        U(KL+J) = 0.0d0
    enddo
    RETURN
ENDIF
IF(NX.EQ.0) THEN
    DO J=1,N
        F(J) = Q*VM(J)
    enddo
    F(IP) = GAMMA-EM
    S = 1.0d0/GAMMA
    ! SET B**T = (F**T)*U
    B(1) = F(1)
    IF(N.GT.1) THEN
        JJ = 1
        DO K=2,N
            SUM = F(K)
            DO I=1,K-1
                SUM = SUM + F(I)*U(JJ+I)
            enddo
            B(K) = SUM
            JJ = JJ+K
        enddo
        IF(IP.LT.N) THEN
            KP1 = IP+1
            KL = IP*KP1/2 + IP
            DO K  = KP1,N
                U(KL) = U(KL) - S*B(K)
                KL    = KL+K
            enddo
        ENDIF
    ENDIF
    ! UPDATE IP BY IP BLOCK OF U
    ! FIRST SET G = D*B
    JJ = 0
    DO I=1,IP
        JJ = JJ + I
        G(I) = U(JJ)*B(I)
    enddo
    TMP = Q*GAMMA
    CALL UDOPT(U,IP,TMP,B,G,ALPHA)
    IF(ALPHA.EQ.0.0D0) RETURN
    G(IP) = G(IP)-ALPHA*S
    C = 1.0d0/ALPHA
    CALL RNK1(U,IP,C,G,IERR)
    RETURN
ENDIF
! NX .NE. 0
TMP = Q*VM(IP)
GAMMA = TMP + EM**2
IF(GAMMA.LT.EPS) RETURN
S = TMP/GAMMA
SUM = 0.0d0
DO I=1,NX
    SUM = SUM + BX(I)*VM(I)
enddo
THETA = EM + Q*SUM
D = EM*THETA + TMP
IF(D.GT.0.0d0) THEN
    C = Q/D
ELSE
    WRITE(IOUNIT,ERRMSG)
    WRITE(IOUNIT,'(10X,'' INPUT ERROR: THETA ='',1PD14.7,'' D ='','//&
        'D14.7,'' EM ='',D14.7/ 40X,10(''*''),'' SKIPPED SMOOTH UPDATE'','//&
        '10(''*''))') THETA,D,EM
    RETURN
ENDIF
! FORM G,B (OBSERVATION)
DO K=1,N
    B(K) = 0.0d0
enddo
DO K=1,NX
    G(K) = S*BX(K)
enddo
IF((NX+2).LT.IP) THEN
    DO K=NX+1,IP-1
        G(K) = 0.0d0
    enddo
ENDIF
G(IP) = (TMP-EM*(1.0d0-EM))/GAMMA
B(IP) = 1.0d0
! PERFORM FIRST UPDATE
CALL UCON(U,N,B,0.0d0,G,IP,F,B,IERR)
IF (C.EQ.0.0d0) RETURN

! FORM SECOND G,B TERMS
DO K=1,IP-1
    G(K) = 0.0d0
enddo
S = EM*C
DO K=1,NX
    G(K) = S*BX(K)
enddo
G(IP) = C
TMP = VM(IP)
VM(IP) = 0.0d0
C = 1.0d0/C
CALL UCON(U,N,VM,C,G,IP,F,B,IERR)
VM(IP) = TMP

RETURN
END