!! Purpose    : COPY VECTOR OF ELEMENTS (VIZ. A PRIORI VARIANCES) ONTO
!!              CONTIGUOUS DIAGONAL ENTRIES OF THE UD VECTOR STORED ARRAY
!!              STARTING WITH the JST-th diagonal entry
!!
!! Parameters :
!!                U --- VECTOR STORED U MATRIX
!!                      For output, VECTOR STORED U MATRIX WITH NEW DIAGONAL
!!                      ENTRIES
!!                N --- COLUMN DIMENSION OF U
!!                D --- INPUT VECTOR TO BE PLACED ON THE DIAGONAL OF U
!!              JST --- LOCATION OF FIRST ELEMENT TO BE COPIED, JSTRT.LE.N

SUBROUTINE D2U(U,N,D,JST)

implicit none

real*8    U(*)
integer*4 N
real*8    D(*)
integer*4 JST
!local
integer*4 JSTM1,JJ,J

JSTM1=JST-1
JJ=JST*JSTM1/2
DO J=JST,N
    JJ=JJ+J
    U(JJ)=D(J-JSTM1)
enddo

RETURN
END