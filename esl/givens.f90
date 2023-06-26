!! Purpose    : COMPUTE THE GIVENS ROTATION PARAMETERS SIG,S,C, SUCH THAT
!!                   |  C   S |  | A |   | SIG |
!!                   |        |  |   | = |     |
!!                   |  S  -C |  | B |   |  0  |
!!              THE DEGENERATE CASE (WITH B=0) PRODUCES C = 1.0*SIGN(A),
!!              S = 0d0, AND SIG = ABS(A)
!!
!! Parameters :
!!              A,B --- VECTOR COEFFICIENTS
!!              SIG --- SQRT(A*A + B*B)
!!                C --- A/SIG
!!                S --- B/SIG  

SUBROUTINE GIVENS(A,B,SIG,S,C)

implicit none
      
real*8 A,B
real*8 SIG
real*8 S,C
!local
real*8 ONE
         
ONE=1.0
SIG = ABS(A)
C   = ABS(B)
IF(C.EQ.0d0) THEN
    C = SIGN(ONE,A)
    S = 0d0
ELSE
    IF(SIG.GE.C) THEN
        SIG = SIG*SQRT(ONE+(C/SIG)**2)
    ELSE
        SIG = C*SQRT(ONE+(SIG/C)**2)
    ENDIF
    C = A/SIG
    S = B/SIG
ENDIF
RETURN
END
