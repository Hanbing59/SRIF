!! Purpose    : COMPUTE W(IRW,N)=F*U WHERE F IS SPARSE AND ONLY THE NON-ZERO
!!              ELEMENTS ARE DEFINED AND U IS VECTOR STORED, UPPER TRIANGULAR
!!              WITH IMPLICITLY DEFINED UNIT DIAGONAL ELEMENTS
!!
!! Parameters :
!!               FEL --- VALUES OF THE NON-ZERO ELEMENTS OF THE F MATRIX
!!              IROW --- ROW INDICES OF THE F ELEMENTS
!!              JCOL --- COLUMN INDICES OF THE F ELEMENTS, i.e.
!!                       F(IROW(K),JCOL(K))=FEL(K)
!!                NF --- NUMBER OF NON-ZERO ELEMENTS OF THE F MATRIX
!!                 U --- UPPER TRIANGULAR, VECTOR STORED MATRIX WITH
!!                       IMPLICITLY DEFINED UNIT DIAGONAL ELEMENTS
!!                       (U(J,J) ARE NOT, IN FACT, UNITY)
!!                 N --- DIMENSION OF U MATRIX
!!                 W --- OUTPUT RESULT
!!             IMAXW --- ROW DIMENSION OF W MATRIX
!!               IRW --- NUMBER OF ROWS IN W,
!!                       (IRW.LE.IMAXW .AND. IRW.GE.MAX(IROW(K)),
!!                       K=1,...,NF, I.E. W MUST HAVE AT LEAST AS
!!                       MANY ROWS AS DOES F.  ADDITIONAL ROWS OF W
!!                       COULD CORRESPOND TO ZERO ROWS OF F.
!!             JDIAG --- DIAGONAL ELEMENT INDICES OF A VECTOR STORED
!!                       UPPER TRIANGULAR MATRIX,
!!                       I.E. JDIAG(K)=K*(K+1)/2=JDIAG(K-1)+K

SUBROUTINE SFU(FEL,IROW,JCOL,NF,U,N,W,IMAXW,IRW,JDIAG)

IMPLICIT none

real*8    FEL(NF)
integer*4 IROW(NF)
integer*4 JCOL(NF)
integer*4 NF
real*8    U(*)
integer*4 N
real*8    W(IMAXW,N)
integer*4 IMAXW,IRW
integer*4 JDIAG(N)
!local
real*8    FIJ,ZERO
integer*4 I,J,K,NEL,JP1,IK

ZERO=0.D0

DO J=1,N
    DO I=1,IRW
        W(I,J)=ZERO
    enddo
enddo

DO NEL=1,NF
    I=IROW(NEL)
    J=JCOL(NEL)
    FIJ=FEL(NEL)
    W(I,J)=W(I,J)+FIJ
    ! THIS ACCOUNTS FOR THE IMPLICIT UNIT DIAGONAL U MATRIX
    ! ELEMENTS. WHEN NON-UNIT DIAGONALS ARE USED, DELETE
    ! THE ABOVE LINE AND USE J INSTEAD OF JP1 BELOW
    IF(J.EQ.N) cycle
    ! WHEN IT IS KNOWN THAT THE LAST COLUMN OF F IS ZERO
    ! THIS 'IF' TEST MAY BE OMITTED
    JP1=J+1
    IK=JDIAG(J)+J
    DO K=JP1,N
        W(I,K)=W(I,K)+FIJ*U(IK)
        IK=IK+K
    enddo
enddo

RETURN
END