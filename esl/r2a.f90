!! Purpose    : TO PLACE THE TRIANGULAR VECTOR STORED MATRIX "R" INTO THE
!!              MATRIX "A" AND TO ARRANGE THE COLUMNS TO MATCH THE DESIRED
!!              NAMA PARAMETER LIST. NAMES IN THE NAMA LIST THAT DO NOT
!!              CORRESPOND TO ANY NAME IN THE NAMR LIST HAVE ZERO ENTRIES
!!              IN THE CORRESPONDING "A" COLUMN
!! Parameters :
!!              R --- INPUT UPPER TRIANGULAR VECTOR STORED ARRAY
!!              N --- DIMENSION OF R
!!           NAMR --- PARAMETER NAMES ASSOCIATED WITH R
!!              A --- MATRIX TO HOUSE THE REARRANGED R MATRIX
!!          IMAXA --- ROW DIMENSION OF A
!!             JA --- NUMBER OF PARAMETER NAMES ASSOCIATED WITH THE OUTPUT
!!                    "A" MATRIX
!!           NAMA --- NAMES ASSOCIATED WITH "A"

SUBROUTINE R2A(R,N,NAMR,A,IMAXA,JA,NAMA)

IMPLICIT none

real*8       R(*)
integer*4    N
CHARACTER(*) NAMR(N)
real*8       A(IMAXA,JA)
integer*4    IMAXA,JA
CHARACTER(*) NAMA(JA)
!local
CHARACTER*17 NAM
integer*4    I,J,K,JJ,JJPK

DO J = 1, JA
    DO I =1, N
        A(I,J) = 0.0
    enddo
enddo
DO J=1, JA
    NAM = NAMA(J)
    DO I=1,N
        IF(NAMR(I).EQ.NAM) exit
    enddo
    if(I.gt.N) cycle
    JJ = I*(I-1)/2
    DO K=1,I
        JJPK = JJ+K
        A(K,J) = R(JJPK)
    enddo
enddo

RETURN
END