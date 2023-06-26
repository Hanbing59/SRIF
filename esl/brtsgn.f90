!! Purpose    : COMPUTE THE BIERMAN MODIFIED R-T-S SMOOTH GAIN FACTORS FOR 
!!              COLORED NOISE MODEL X   = PHI  X   + BW
!!
!! Parameters :       
!!                   U --- VECTOR STORED U-D FILTER COVARIANCE FACTORS PRIOR 
!!                         TO THE COLORED NOISE COMPONENT UPDATE
!!                   N --- 
!!                  IP --- LOCATION OF COLORED NOISE PARAMETER 
!!                         IP.GE.1, AND IP.GT.(IBZ + NBX)
!!                  EM --- COLORED NOISE (EXPONENTIAL) MULTIPLIER
!!                  BX --- PROCESS NOISE SENSITIVITY MULTIPLIER TERMS.
!!                         WHEN 0 NOT USED
!!                 IBZ --- NUMBER OF LEADING ZEROS IN B VECTOR
!!                 NBX --- NUMBER OF ROWS IN BX (B = (0, BX ,0 ,1, 0))
!!                   Q --- PROCESS NOISE COMPONENT A PRIORI VARIANCE
!!                XEST --- FILTER STATE VECTOR PRIOR TO THE COLORED NOISE 
!!                         COMPONENT UPDATE, USED ONLY WHEN IEST = .TRUE.
!!               SGAIN --- AUGMENTED VECTOR OF SMOOTHER GAINS
!!                         WHEN IEST = .FALSE., THE (N+1)ST COMPONENT IS LEFT UNALTERED
!!              SGSTAR --- SCALAR MULTIPLIER USED IN THE SMOOTH U-D FACTOR 
!!                         COVARIANCE TIME UPDATE
!!                IEST --- IF .TRUE., THE FILTER ESTIMATE IS USED TO COMPUTE 
!!                         SGAIN(N+1) THAT IS USED AS PART OF THE SMOOTH 
!!                         ESTIMATE TIME UPDATE.
!!                IERR --- IF .EQ. -1, THE COMPUTED SMOOTH GAIN IS SINGULAR
!!                            .LE. -2, THE SMOOTH GAIN COMPUTATION FAILED.
      
SUBROUTINE BRTSGN (U,N,IP,EM,BX,IBZ,NBX,Q,XEST,SGAIN,SGSTAR,Y,C,S,IEST,IERR)

implicit none

!     TO COMPUTE THE BIERMAN MODIFIED R-T-S SMOOTH GAIN FACTORS FOR
!     COLORED NOISE MODEL X   = PHI  X   + BW
!                          T+      M  T-     T
!
!                               T
!     WHERE PHI  = I + (M-1)E  E
!              M             IP IP
!
!       (E  (J) = 0 J .NE. IP ; = 1 FOR J = IP)
!         IP
      
real*8    U(N*(N+1)/2)
integer*4 N
integer*4 IP
real*8    EM
real*8    BX(*)
integer*4 IBZ
integer*4 NBX
real*8    Q
real*8    XEST(N)
real*8    SGAIN(N+1)
real*8    SGSTAR
real*8    Y(N),C(N),S(N)
LOGICAL   IEST
integer*4 IERR
!local
integer*4 I,J,K
integer*4 IJ,JJ,JJS,IK
integer*4 IPIP,IPM1,IPP1
integer*4 IFV,NP1,JCOL
REAL*8    SIG,SN,CN
real*8    YNP1,SIGSTR
real*8    TMP,TMPA
real*8    EPSCN,EPS
DATA      EPS/1.0D-307/,EPSCN/1.0D-6/
EXTERNAL  GIVENS

IERR = 0
IF(EM.GT.(1d0-EPSCN) .AND. Q.EQ.0d0) THEN
    DO J = 1, N+1
        SGAIN(J) = 0d0
    enddo
    SGSTAR=0d0
    RETURN
ENDIF
IPM1 = IP - 1
IPP1 = IP + 1
IPIP = IP*IPP1/2
NP1  = N + 1

! SOLVE UY = BM; BM(IP) = 1, E(J) = M*B(J) J .LT. IP
DO I = 1, IPM1
    Y(I) = 0d0
enddo
DO I = 1, NBX
    Y(IBZ+I) = BX(I) * EM
enddo
Y(IP) = 1d0
IK = IPIP
DO J = IP,2,-1
    TMP = Y(J)
    DO I = J-1,1,-1
        IK = IK-1
        Y(I) = Y(I) - U(IK)*TMP
    enddo
    !DIAGONAL OF U MATRIX
    IK = IK-1
enddo
! NON-DIAGONAL Q
IF(NBX.GT.0) THEN
    ! FORM SGAIN**T = -((E**T)U*D;0)
    DO I = 1, IPM1
        SGAIN(I) = 0d0
    enddo
    JJ = IPIP
    SGAIN(IP) = -U(JJ)
    DO I = IPP1, N
        IJ = JJ+IP
        JJ = JJ+I
        SGAIN(I) = -U(JJ)*U(IJ)
        Y(I) = 0d0
    enddo
    SIGSTR = 0d0
    YNP1 = -EM
    ! SOLVE SYSTEM
    ! D*(U**T)*PHI*SGAIN + Y*SGSTAR = -D*(U**T)*E
    ! Q*(B**T) * SGAIN - EM*SGSTAR  = 0
    ! BY GIVENS TRANSFORMING IT TO TRIANGULAR FORM AND FORWARD
    ! SOLVING AS WE GO, BUT SAVING THE GIVENS DIRECTION SINES AND
    ! COSINES FOR LATER RECONSTRUCTION OF THE GAIN FROM THE GIVENS
    ! ROTATED SOLUTION.
    JJS = 0
    DO IFV = 1, N
        JJS = JJS+IFV
        JJ = JJS
        TMP = U(JJ)
        IF(IFV.EQ.IP) TMP = TMP*EM
        CALL GIVENS(TMP,Y(IFV),SIG,SN,CN)
        C(IFV) = CN
        S(IFV) = SN
        !TMP IS THE SOLVED VALUE FOR SGAIN(IFV)
        IF(SIG.GT.EPS) THEN
            TMP = SGAIN(IFV)/SIG
            ! APPLY GIVENS TO Y AND COLUMN B IFV OF THE COEFFICIENT MATRIX,
            ! ADJUST RHS CORRESPONDING TO THE FORWARD SOLVE.
            DO I = IFV+1, N
                IJ = JJ+IFV
                JJ = JJ+I
                TMPA = U(JJ) * U(IJ)
                IF(IFV.EQ.IP) TMPA = TMPA*EM
                SGAIN(I) = SGAIN(I) - (CN*TMPA + SN*Y(I))*TMP
                Y(I) = SN*TMPA  - CN*Y(I)
            enddo
            ! B HAS THE FORM 0000000 BBBBBBBBBBBB 0 IP 0000000
            !                |_IBZ_| |_IBZ+NBX_|
            ! CASE#             3         2       3  1  3
            IF(IFV.EQ.IP) THEN
                ! CASE #1, B IS "VIRTUALLY" ONE
                SIGSTR = SIGSTR - (CN*Q + SN*YNP1) * TMP
                YNP1 = SN*Q - CN*YNP1
            ELSE IF(IFV.GT.IBZ .AND. IFV.LE.(IBZ+NBX)) THEN
                ! CASE #2, B IS "REAL"
                TMPA = Q * BX(IFV-IBZ)
                SIGSTR = SIGSTR  - (CN*TMPA + SN*YNP1) * TMP
                YNP1 = SN*TMPA - CN*YNP1
            ELSE
                ! CASE #3, B IS "VIRTUALLY" ZERO
                SIGSTR = SIGSTR - SN*YNP1*TMP
                YNP1 = -CN*YNP1
            ENDIF
            SGAIN(IFV) = TMP
        ENDIF
    enddo
    ! NOW, CONSTRUCT SOLUTION BY APPLYING GIVENS TRANSFORMATION TO THE
    ! SGAIN VECTOR JUST COMPUTED
    IF(ABS(YNP1).LT.EPS) THEN
        IF(ABS(SIGSTR).LT.EPS) THEN
            SIGSTR = 0d0
            IERR = -1
        ELSE
            SIGSTR = 0d0
        ENDIF
    ELSE
        SIGSTR = SIGSTR/YNP1
    ENDIF
    DO J = N, 1, -1
        TMP = SGAIN(J)
        SGAIN(J) = C(J) * TMP + S(J) * SIGSTR
        SIGSTR   = S(J) * TMP - C(J) * SIGSTR
    enddo
    IF(SIGSTR.GT.0)  IERR = -3
    IF(IEST) THEN
        TMP = ((EM-1d0)*SGAIN(IP) + 1d0) * XEST(IP)
        DO J = 1, N
            TMP = TMP + SGAIN(J) *XEST(J)
        enddo
        SGAIN(NP1) = TMP
    ENDIF
    SGSTAR = -SIGSTR
    RETURN
ENDIF

! DIAGONAL Q
! Y(I) = 0 FOR I > IP, BUT THESE VALUES ARE NOT NEEDED
! (IMPLICITLY) FORM A MATRIX, A = (D(U**T)PHI+(Q/EM)(Y*E**T)PHI,
! TRIANGULARIZE A, SAVING Y AS WE GO AND THE DIRECTION COSINE TERMS OF THE 
! GIVENS ROTATIONS

! FORM COLUMN IP OF A, SCALED BY EM
DO I = 1, IPM1
    Y(I) = Q*Y(I)
enddo
TMP = EM*EM
JJ = IPIP
Y(IP) = U(JJ)*TMP+Q
DO I = IPP1, N
    IJ = JJ+IP
    JJ = JJ+I
    Y(I) = U(JJ)*U(IJ)*TMP
enddo
! IMPLICITLY APPLY GIVENS TO COLUMNS JCOL AND IP, SAVING GIVENS C AND S 
! AND ONLY EXPLICITLY MODIFYING Y
JJS = 0
DO JCOL = 1, IPM1
    SGAIN(JCOL) = 0d0
    JJS = JJS+JCOL
    JJ  = JJS
    CALL GIVENS(U(JJ),Y(JCOL),SIG,SN,CN)
    DO I = JCOL + 1 , N
        IJ   = JJ+JCOL
        JJ   = JJ+I
        Y(I) = SN*U(JJ)*U(IJ) - CN*Y(I)
    enddo
    C(JCOL) = CN
    S(JCOL) = SN
enddo

! SOLVE AV = -DU SUB(IP)
IF(ABS(Y(IP)).LE.EPS) THEN
    DO I = 1, NP1
        SGAIN(I) = 0d0
    enddo
    SGSTAR = 0d0
    RETURN
ENDIF
SGAIN(IP) = -U(IPIP)/Y(IP)
JJ = IPIP
DO I = IPP1, N
    IJ  = JJ+IP
    IK  = IJ
    JJ  = JJ+I
    TMP = -U(IJ)
    DO K = IPP1, I-1
        IK   = IK+1
        TMP  = TMP - U(IK)*SGAIN(K)
    enddo
    IF(U(JJ).GT.EPS) THEN
        SGAIN(I) = TMP - Y(I) * SGAIN(IP)/U(JJ)
    ELSE
        SGAIN(I) = 0d0
        IF(ABS(Y(I)).GE.EPS) THEN
            IERR = -1
        ENDIF
    ENDIF
enddo
! SOLVE FOR SGAIN BY APPLYING GIVENS TRANSFORMATION TO 
! SGAIN (OVERWRITING RESULT) AND THEN SCALING SGAIN(IP) BY EM
DO J = IPM1, 1, -1
    TMP       = SGAIN(J)
    SGAIN(J)  = C(J)*TMP + S(J)*SGAIN(IP)
    SGAIN(IP) = S(J)*TMP - C(J)*SGAIN(IP)
enddo
SIGSTR       = -Q*SGAIN(IP)
TMP          = EM*SGAIN(IP)
SGAIN(IP)    = TMP+1.0D0
! NOTE SGAIN(IP) IS REALLY (SGAIN(IP)-1)/EM
IF(IEST) THEN
    TMPA = (EM-1.0D0)*TMP*XEST(IP)
    DO I = 1, N
        TMPA = TMPA + SGAIN(I)*XEST(I)
    enddo
    SGAIN(NP1) = TMPA
ENDIF
SGSTAR = SIGSTR

RETURN
END