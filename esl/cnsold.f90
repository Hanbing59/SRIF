!! Purpose    : 
!!              SUBROUTINES NEEDED:
!!               UDGEN  <ESL>
!!               RNK1   <ESL>
!!               SDOT1  <ESL-BLAS>  SDOT FROM BLAS, BUT RETURNS DP
      
SUBROUTINE CNSOLD(U,N,BX,NXMAX,NBX,IBZ,UBX,NP,EM,SGSTAR,IP,SGAIN,MAXSG,XSTAR,F,G,H,V,UPDATU,UPDATX,IERR)

IMPLICIT NONE

real*8       U(N*(N+1)/2)
integer*4    N
real*8       BX(NXMAX,*)
integer*4    NXMAX
integer*4    NBX(*)
integer*4    IBZ(*)
CHARACTER*20 UBX
integer*4    NP
real*8       EM(*)
real*8       SGSTAR(*)
integer*4    IP(*)
real*8       SGAIN(MAXSG,*)
integer*4    MAXSG
real*8       XSTAR(N)
real*8       F(N),G(N),H(N),V(N)
LOGICAL      UPDATU,UPDATX
integer*4    IERR
!local
INTEGER*4    I,J,L
integer*4    IBX,IZ
integer*4    NP1
real*8       SDOT1,SUM,TMP,TMPA,TMPB
LOGICAL      ISW
EXTERNAL     UDGEN

IERR = 0
NP1 = N+1
IF(INDEX(UBX,'ALL').NE.0 .OR. INDEX(UBX,'all').NE.0) THEN
    ISW = .TRUE.
ELSE IF(INDEX(UBX,'SAME').NE.0 .OR. INDEX(UBX,'same').NE.0)THEN
    IBX =  NBX(1)
    IZ  =  IBZ(1)
    ISW = .FALSE.
ELSE
    IBX =  0
    IZ  =  0
    ISW =  .FALSE.
ENDIF
IF(NP.EQ.0) RETURN

IF(UPDATX) THEN
    DO J = NP, 1, -1
        L = IP(J)
        IF(ISW) THEN
            IBX = NBX(J)
            IZ  = IBZ(J)
        ENDIF
        ! WHEN COLUMN J OF BX CONTAINS ONLY ZEROS, USE SIMPLER COMPUTATION
        IF(IBX.EQ.0) THEN
            XSTAR(L) = XSTAR(L) - SDOT1(N,SGAIN(1,J),1,XSTAR,1) + SGAIN(N+1,J)
        ELSE
            ! GAMMA
            TMP        = SDOT1(N,SGAIN(1,J),1,XSTAR,1)
            ! GAMMAM
            TMPA       = TMP*EM(J)
            TMPB       = SGAIN(1,NP1)*EM(J)
            DO I = 1, IBX
                XSTAR(IZ+I) = XSTAR(IZ+I)-BX(I,J)*(TMPA+XSTAR(L)-TMPB)
            enddo
            XSTAR(L) = SGAIN(NP1,J)-TMP
        ENDIF
    enddo
ENDIF

IF (UPDATU) THEN
    DO J = NP,1,-1
        IF(ISW) THEN
            IBX = NBX(J)
            IZ  = IBZ(J)
        ENDIF
        L = IP(J)
        DO I=1,N
            G(I)=0.0D0
        ENDDO
        ! WHEN COLUMN J OF BX CONTAINS ONLY ZEROS, USE SIMPLER COMPUTATION
        IF(IBX.EQ.0) THEN
            G(L) = 1.0
            CALL UDGEN(U,N,SGAIN(1,J),SGSTAR(J),G,L,F,H,IERR)
            IF(IERR.NE.0) THEN
                WRITE(6,*)' ERROR RETURN FROM UDGEN = ',IERR,' STATE = ',J
            ENDIF
        ELSE
            TMPA = 1.0-EM(J)
            DO I = 1, IZ
                G(I) = 0.0
            enddo
            DO I = 1, IBX
                G(I+IZ) = TMPA*BX(I,J)
            enddo
            DO I = 1, N
                H(I) = 0.0
            enddo
            H(L) = 1.0
            CALL UDGEN(U,N,H,0d0,G,IZ+IBX,F,V,IERR)
            DO I = 1, L-1
                G(I) = 0.0
            enddo
            DO I = 1, IBX
                G(IZ+I) = EM(J)*BX(I,J)
            enddo
            G(L) = 1.0
            TMP = SGAIN(L,J)
            SUM = 0.0
            DO I = 1, IBX
                SUM = SUM + SGAIN(IZ+I,J)*BX(I,J)
            enddo
            SGAIN(L,J) = TMP + SUM*TMPA + 1.0
            CALL UDGEN(U,N,SGAIN(1,J),SGSTAR(J),G,L,F,V,IERR)
            SGAIN(L,J) = TMP
            IF (IERR.NE.0) THEN
                WRITE(6,*)' ERROR RETURN FROM UDGEN = ',IERR,' STATE = ',J
            ENDIF
        ENDIF
    enddo
ENDIF
RETURN
END
