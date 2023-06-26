!! Purpose    : COMPUTES W=PHI*U, WHERE PHI IS A RECTANGULAR MATRIX WITH
!!              IMPLICITLY DEFINED COLUMNS OF TRAILING ZEROS AND U IS A
!!              VECTOR STORED UPPER TRIANGULAR MATRIX.
!! Parameters :
!!              PHI --- INPUT RECTANGULAR MATRIX, IRPHI.LE.MAXPHI
!!           MAXPHI --- ROW DIMENSION OF PHI
!!            IRPHI --- NO. ROWS OF PHI
!!            JCPHI --- NO. COLS OF PHI
!!                U --- UPPER TRIANGULAR VECTOR STORED MATRIX
!!                N --- DIMENSION OF U MATRIX (JCPHI.LE.N)
!!                W --- OUTPUT, RESULT OF PHI*U, IF THE COMPILER PERMITS,
!!                      W CAN OVERWRITE PHI
!!            IMAXW --- ROW DIMENSION OF W

SUBROUTINE PHIU(PHI,MAXPHI,IRPHI,JCPHI,U,N,W,IMAXW)

implicit none

real*8    PHI(MAXPHI,JCPHI)
integer*4 MAXPHI
integer*4 IRPHI,JCPHI
real*8    U(*)
integer*4 N
real*8    W(IMAXW,N)
integer*4 IMAXW
!local
integer*4 I,J,K,L,NP2,KJS,JM1
real*8    SUM

DO I=1,IRPHI
    W(I,1)=PHI(I,1)
enddo
IF(N.EQ.1) RETURN

NP2=N+2
!CDIR$ IVDEP
DO L=2,N
    J=NP2-L
    KJS=(J-1)*J/2
    JM1=J-1
    DO I=1,IRPHI
        IF(J.LE.JCPHI) THEN
            SUM=PHI(I,J)
        ELSE
            SUM=0.D0
            JM1=JCPHI
        ENDIF
        DO K=1,JM1
            SUM=SUM+PHI(I,K)*U(KJS+K)
        enddo
        W(I,J)=SUM
    enddo
enddo

RETURN
END