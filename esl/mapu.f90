!! Purpose    : COMPUTE MAPPED COVARIANCE POUT = (PHI*U)*D*(PHI*U)**T
!!              IT IS ASSUMED THAT THE MAPPING MATRIX IS OF THE FORM:
!!                  PHI   0 WITH ONLY THE PHI ENTRIES EXPLICIT
!!
!! Parameters : 
!!              PHI --- MAPPING MATRIX
!!           MAXPHI --- ROW DIMENSION OF PHI, IRPHI.LE.MAXPHI
!!            IRPHI ---
!!            JCPHI ---
!!                U --- VECTOR STORED U-D COVARIANCE FACTORS
!!                N --- COLUMN DIMENSION OF U
!!                W --- SCRATCH MATRIX. IF THE COMPILER ALLOWS, W AND PHI 
!!                      CAN SHARE STORAGE.
!!             MAXW --- ROW DIMENSION OF W (IRPHI.LE.MAXW)
!!                D --- DIAGONAL ENTRIES OF INPUT U MATRIX
!!                      D(1), FLAG TO INDICATE WHETHER ONLY DIAGONAL 
!!                      ENTRIES (SIGMAS) OF THE COVARIANCE BEING REQUESTED.
!!                      IF D(1).GE.0, COMPUTE ENTIRE COVARIANCE MATRIX, 
!!                      ELSE COMPUTE ONLY DIAGONAL COVARIANCE SQUARE-ROOT VALUES.
!!             POUT --- IF D(1).GE.0, THEN VECTOR STORED OUTPUT COVARIANCE
!!                      MATRIX WHICH (IF THE COMPILER ALLOWS) MAY OVERWRITE U.
!!                      Otherwise, POUT HAS ONLY IRPHI STANDARD DEVIATIONS.          


SUBROUTINE MAPU(PHI,MAXPHI,IRPHI,JCPHI,U,N,W,MAXW,D,POUT)

implicit none

real*8    PHI(MAXPHI,IRPHI)
integer*4 MAXPHI,IRPHI,JCPHI
real*8    U(*)
integer*4 N
real*8    W(MAXW,N)
integer*4 MAXW
real*8    D(N)
real*8    POUT(*)
!local
LOGICAL   WNTVAR
real*8    SUM,TMP
integer*4 I,J,K,L,JJ,IJ
integer*4 NP2,IROW,IM1
EXTERNAL  PHIU

WNTVAR=D(1).LT.0.0
CALL PHIU(PHI,MAXPHI,IRPHI,JCPHI,U,N,W,MAXW)
! STORE DIAGONALS OF U IN VECTOR D
JJ=0
DO J=1,N
    JJ=JJ+J
    D(J)=U(JJ)
enddo
! CREATE POUT=W*D*(W**T)
IF(IRPHI.ne.1) then
    IJ = IRPHI*(IRPHI+1)/2+1
    NP2= IRPHI + 2
    DO IROW=2,IRPHI
        I=NP2-IROW
        SUM=0.0
        DO K=1,N
            TMP=D(K)*W(I,K)
            SUM=SUM+W(I,K)*TMP
            W(I,K)=TMP
        enddo
        IF(WNTVAR) THEN
            POUT(I) = DSQRT(SUM)
        ELSE
            IJ=IJ-1
            POUT(IJ)=SUM
            ! DIAGONAL ELEMENT OF OUTPUT COVARIANCE
            IM1=I-1
            DO L=1,IM1
                K=I-L
                SUM=0.0
                DO J=1,N
                    SUM=SUM+W(K,J)*W(I,J)
                enddo
                IJ=IJ-1
                POUT(IJ)=SUM
            enddo
        ENDIF
    enddo
endif
SUM=0.0
DO J=1,N
    SUM=SUM+D(J)*W(1,J)**2
enddo
IF(WNTVAR) THEN
    POUT(1) = DSQRT(SUM)
ELSE
    POUT(1)=SUM
ENDIF

RETURN
END