!! Purpose    : COMPUTE THE SOLUTION VECTOR X TO THE MATRIX EQUATION
!!              R*X=Z, WHERE R IS AN N BY N UPPER-TRIANGULAR VECTOR
!!              MATRIX. COMPUTATION IS MINIMIZED AND ACCURACY ENHANCED
!!              BY NOT COMPUTING THE MATRIX INVERSE. IF R IS SINGULAR,
!!              WITH ITS JTH DIAGONAL ELEMENT AS ITS LAST DIAGONAL ZERO
!!              VALUE, X(I) IS COMPUTED  ONLY FOR J .LT. I .LE. N,
!!              THE OTHER COMPONENTS, X(I), ARE LEFT UNCHANGED FOR
!!              1 .LE. I .LE. J, AND THE VALUE OF J IS RETURNED AS IERR
!! Parameters :
!!              R --- VECTOR STORED UPPER TRIANGULAR MATRIX
!!              N --- MATRIX VECTOR DIMENSION, N .GT. 0
!!              Z --- RIGHT HAND SIDE VECTOR OF EQUATION R*X=Z
!!              X --- SOLUTION VECTOR, WHICH MAY OVERWRITE Z,
!!                    COMPILER PERMITTING.
!!           IERR --- ERROR RETURN FLAG
!!                    0, NORMAL RETURN
!!                    J, SOLUTION, X, VALID ONLY FOR COMPONENTS .GT. J

SUBROUTINE RINZ(R,N,Z,X,IERR)

implicit none

real*8    R(*)
integer*4 N
real*8    Z(N)
real*8    X(N)
integer*4 IERR
!local
integer*4 JJ,IX,I,IJ,J
real*8    SUM,EPS,S,DX,ZERO
DATA EPS/1.0D-307/
DATA ZERO/0.0D0/

IERR = 0
JJ = N*(N+1)/2
IX = N
SUM = ZERO
IF(N.ne.1) then
    DO I = 2,N
        S = R(JJ)
        DX = Z(IX) - SUM
        IF(DABS(S).LE.EPS*DABS(DX)) then
            IERR=IX
            return
        endif
        X(IX) = DX/S
        IJ = JJ-1
        SUM = ZERO
        DO J = IX, N
            SUM = SUM+R(IJ)*X(J)
            IJ = IJ+J
        enddo
        JJ = JJ-IX
        IX = IX-1
    enddo
endif
S = R(1)
DX = Z(1) - SUM
IF(DABS(S).gt.EPS*DABS(DX)) then
    X(1) = DX/S
else
    IERR = IX
endif

RETURN
END