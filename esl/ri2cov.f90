!! Purpose    : COMPUTE THE STANDARD DEVIATIONS AND, IF DESIRED, THE COVARIANCE
!!              MATRIX OF A VECTOR STORED UPPER TRIANGULAR SQUARE ROOT
!!              COVARIANCE MATRIX. THE OUTPUT COVARIANCE MATRIX IS ALSO VECTOR STORED.
!! Parameters :
!!              RINV --- VECTOR STORED UPPER TRIANGULAR COVARIANCE SQUARE ROOT.
!!                       (RINV = RINVERSE IS THE INVERSE OF THE SRIF MATRIX)
!!                 N --- DIMENSION OF THE RINV MATRIX, N .GE. 1
!!               SIG --- OUTPUT VECTOR OF STANDARD DEVIATIONS
!!            COVOUT --- VECTOR STORED COVARIANCE (IF THE COMPILER PERMITS
!!                       COVOUT MAY OVERWRITE RINV)
!!              KROW --- .GT.0, COMPUTES THE COVARIANCE AND SIGMAS
!!                       .LT.0, COMPUTES ONLY THE SIGMAS OF THE FIRST
!!                              KROW VARIABLES OF THE RINV MATRIX
!!                       .EQ.0, NO COVARIANCE, BUT ALL SIGMAS
!!                              (E.G. USE ALL N ROWS OF RINV)
!!              KCOL --- NUMBER OF COLUMNS OF COVOUT THAT ARE COMPUTED.
!!                       IF KCOL .LE. 0 THEN KCOL=KROW.
!!                       IF KROW .LE. 0 THIS INPUT IS IGNORED.

SUBROUTINE RI2COV(RINV,N,SIG,COVOUT,KROW,KCOL)

implicit none

real*8    RINV(*)
integer*4 N
real*8    SIG(N)
real*8    COVOUT(*)
integer*4 KROW,KCOL
!local
integer*4 LIM,IKS,J,IK,K,KKOL,JJ,NM1,IJS,JP1,I,IMJ,IJK
real*8    SUM,ZERO
DATA ZERO/0.0D0/

LIM = N
IF(KROW.NE.0) LIM = IABS(KROW)
! COMPUTE SIGMAS
IKS = 0
DO J = 1, LIM
    IKS = IKS+J
    SUM = ZERO
    IK = IKS
    DO K=J, N
        SUM = SUM+RINV(IK)**2
        IK = IK+K
    enddo
    SIG(J) = DSQRT(SUM)
enddo
IF(KROW.LE.0) RETURN
KKOL = KCOL
IF(KKOL.LE.0) KKOL = KROW
JJ = 0
IF(N.ne.1) then
    ! COMPUTE COVARIANCE
    NM1 = LIM
    IF(KROW.EQ.N) NM1 = N-1
    DO J=1, NM1
        JJ = JJ+J
        COVOUT(JJ) = SIG(J)**2
        IJS = JJ+J
        JP1 = J+1
        DO I = JP1, KKOL
            IK  = IJS
            IMJ = I-J
            SUM = ZERO
            DO K = I, N
                IJK = IK + IMJ
                SUM = SUM + RINV(IK)*RINV(IJK)
                IK = IK+K
            enddo
            COVOUT(IJS) = SUM
            IJS = IJS+I
        enddo
    enddo
endif
IF(KROW.EQ.N) COVOUT(JJ+N) = SIG(N)**2

RETURN
END