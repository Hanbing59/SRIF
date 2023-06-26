    !! Purpose    : COLORED NOISE UD SMOOTHER
    !!              U MATRIX UPDATE (WHEN UPDATU IS TRUE)
    !!              R-T-S SMOOTHER GAIN = I-E * (SGAIN)**T
    !!              ESTIMATE UPDATE (WHEN UPDATX IS TRUE)
    !!
    !! Parameters :
    !!                U --- VECTOR STORED U-D SMOOTH COVARIANCE FOR ESTIMATE AT STAGE I+1
    !!                      For output, UPDATED U-D FACTOR SMOOTH COVARIANCE FACTORS
    !!                      CORRESPONDING TO THE STATE AT TIME "K"
    !!                N --- FILTER MODEL DIMENSION
    !!               IP --- LOCATIONS OF COLORED NOISE STATES IN THE FILTER MODEL
    !!              IBZ --- IBZ(J) REPRESENTS THE NUMBER OF LEADING ZEROES IN COLUMN J OF BX
    !!               BX --- PROCESS NOISE SENSITIVITY MULTIPLIER MATRIX, NOT USED WHEN NBX() IS 0
    !!            NXMAX --- MAX ROW DIMENSION OF BX MATRIX
    !!              NBX --- NUMBER OF ROWS PER COLUMN OF BX
    !!              UBX --- WHEN UBX .EQ. 'SAME', IT IS ASSUMED THAT ALL THE COLUMNS
    !!                      HAVE THE SAME NUMBER OF ROWS (NBX(1)) AND THE SAME NUMBER
    !!                      OF LEADING ZEROS, VIZ IBZ(2..NP) = IBZ(1).
    !!                      WHEN UBX .EQ. 'NONE', NBX IS TREATED AS IDENTICALLY ZERO,
    !!                      I.E. DIAGONAL PROCESS NOISE.
    !!                      WHEN UBX .EQ. 'ALL', ALL NBX AND IBZ ELEMENTS ARE USED.
    !!                      'NONE' IS THE DEFAULT. WHEN 'NONE', OR 'SAME' ONE NEED NOT
    !!                      SPECIFY INTEGER VECTORS FOR NBX AND IBZ
    !!                      (E.G. CALL CNS(...,2,0,'SAME',...) MEANS THAT EACH COLUMN
    !!                      OF BX HAS 2 ROWS, AND THERE ARE NO (IMPLICIT) LEADING ZEROS.
    !!               NP --- NUMBER OF PROCESS NOISE STATES
    !!               EM --- EXPONENTIAL MULTIPLIERS
    !!           SGSTAR --- SMOOTH COEFFICIENT MULTIPLIERS
    !!            SGAIN --- R-T-S MODIFIED GAIN
    !!                        *
    !!                       G  = I-E(SGAIN)**T
    !!                      WHERE E IS A ZERO VECTOR, EXCEPT FOR UNITY IN ONE
    !!                      LOCATION, IP(J).
    !!                      ELEMENT N+1 CONTAINS: SGAIN**T * X(FILTER)
    !!            MAXSG --- MAX ROW DIMENSION OF SGAIN
    !!            XSTAR --- SMOOTH ESTIMATES, ON ENTRY CORRESPONDING TO TIME "K+1"
    !!                      For output, SMOOTH ESTIMATES AT TIME "K" WHEN UPDATX IS TRUE
    !!           UPDATU --- WHEN TRUE, COMPUTE UPDATED UD MATRIX
    !!           UPDATX --- WHEN TRUE, COMPUTE UPDATED ESTIMATES
    !!             IERR --- ERROR RETURN FROM RNK1 (CALLED BY UDGEN)
    !!                       0, NORMAL RETURN
    !!                       J, ERROR RETURN DUE TO A NEGATIVE COMPUTED J-TH
    !!                          DIAGONAL ENTRY
    !!                      -1, UD RESULT IS SINGULAR
    !!            ISING --- TOTAL NUMBER OF SINGULARITIES ENCOUNTERED

SUBROUTINE CNS(U,N,IP,IBZ,BX,NXMAX,NBX,UBX,NP,EM,SGSTAR,SGAIN,MAXSG,XSTAR,F,G,H,V,UPDATU,UPDATX,IERR,ISING)

implicit none

real*8       U(N*(N+1)/2)
integer*4    N
integer*4    IP(*),IBZ(*)
real*8       BX(NXMAX,*)
integer*4    NXMAX
integer*4    NBX(*)
CHARACTER*20 UBX
integer*4    NP
real*8       EM(*)
real*8       SGSTAR(*)
real*8       SGAIN(MAXSG,*)
integer*4    MAXSG
real*8       XSTAR(N)
real*8       F(N),G(N),H(N),V(N)
LOGICAL      UPDATU,UPDATX
integer*4    IERR,ISING
!local
INTEGER*4    I,J,K,L
integer*4    IBX,IZ
integer*4    NP1
real*8       SDOT1,SUM,TMP,TMPA
!USE ISW TO DETERMINE IF NBX, IBZ ARE IN EFFECT
LOGICAL      ISW
EXTERNAL     UDGEN

IERR = 0
NP1 = N+1
IF(INDEX(UBX,'ALL').NE.0 .OR. INDEX(UBX,'all').NE.0) THEN
    ISW = .TRUE.
ELSE IF(INDEX(UBX,'SAME').NE.0 .OR. INDEX(UBX,'same').NE.0) THEN
    IBX =  NBX(1)
    IZ  =  IBZ(1)
    ISW = .FALSE.
ELSE
    !NONE IS THE DEFAULT
    IBX =  0
    IZ  =  0
    ISW =  .FALSE.
ENDIF
!IF NO PROCESS NOISE RETURN
IF(NP.EQ.0) RETURN

IF(UPDATX) THEN
    !ESTIMATE UPDATE
    DO J = NP, 1, -1
        L = IP(J)
        IF(ISW) THEN
            IBX = NBX(J)
            IZ  = IBZ(J)
        ENDIF
        ! WHEN COLUMN J OF BX CONTAINS ONLY ZEROS, USE SIMPLER COMPUTATION
        SDOT1 = 0d0
        DO K = 1, N
            SDOT1 = SDOT1 + SGAIN(K,J)*XSTAR(K)
        enddo
        IF(IBX.EQ.0) THEN
            XSTAR(L) = XSTAR(L) - SDOT1 + SGAIN(NP1,J)
        ELSE
            ! GAMMA
            TMP = SDOT1
            DO I = 1, IBX
                XSTAR(IZ+I) = XSTAR(IZ+I)-BX(I,J) * (XSTAR(L) + (TMP-SGAIN(NP1,J)) * EM(J))
            enddo
            XSTAR(L) = SGAIN(NP1,J)-TMP
        ENDIF
    enddo
ENDIF

IF(UPDATU) THEN
    ! U UPDATE
    DO J = NP,1,-1
        IF(ISW) THEN
            IBX = NBX(J)
            IZ  = IBZ(J)
        ENDIF
        L = IP(J)
        ! WHEN COLUMN J OF BX CONTAINS ONLY ZEROS, USE SIMPLER COMPUTATION
        IF(IBX.EQ.0) THEN
            DO I = 1, N
                G(I) = 0d0
            enddo
            G(L) = 1.0d0
            CALL UDGEN(U,N,SGAIN(1,J),SGSTAR(J),G,L,F,H,IERR)
            IF(IERR.EQ.-1) THEN
                ISING = ISING + 1
            ELSE IF(IERR.NE.0) THEN
                WRITE(6,*)' ERROR RETURN FROM UDGEN = ',IERR,' STATE = ',J
            ENDIF
        ELSE
            ! NON DIAGONAL PROCESS NOISE
            ! COMPUTE P:=(I-G1*H1)P(I-G1H1)**T
            TMPA = 1.0-EM(J)
            DO I = 1, IZ
                G(I) = 0d0
            enddo
            DO I = 1, IBX
                G(I+IZ) = TMPA*BX(I,J)
            enddo
            DO I = 1, N
                H(I) = 0d0
            enddo
            H(L) = 1.0d0
            CALL UDGEN(U,N,H,0d0,G,IZ+IBX,F,V,IERR)
            ! P:= (I-G2*H2)P(I-G2H2)**T + SIGSTAR*G2*(G2)**T
            DO I = 1, L-1
                G(I) = 0d0
            enddo
            DO I = 1, IBX
                G(IZ+I) = EM(J)*BX(I,J)
            enddo
            G(L) = 1.0d0
            TMP = SGAIN(L,J)
            SUM = 0d0
            DO I = 1, IBX
                SUM = SUM + SGAIN(IZ+I,J)*BX(I,J)
            enddo
            SGAIN(L,J) = TMP + SUM*TMPA + 1.0d0
            CALL UDGEN(U,N,SGAIN(1,J),SGSTAR(J),G,L,F,V,IERR)
            SGAIN(L,J) = TMP
            IF(IERR.EQ.-1) THEN
                ISING = ISING+1
            ELSE IF(IERR.NE.0) THEN
                WRITE(6,*)' ERROR RETURN FROM UDGEN2 = ',IERR,' STATE = ',J
            ENDIF
        ENDIF
    enddo
ENDIF
RETURN
END