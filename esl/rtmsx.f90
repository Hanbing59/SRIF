!! Purpose    : MULTIPLY R (A VECTOR STORED MATRIX) AND X (A VECTOR)
!!
!! Parameters :
!!              R --- INPUT VECTOR STORED MATRIX
!!              N --- DIMENSION OF INPUT VECTORS
!!              X --- INPUT VECTOR
!!              V --- OUTPUT VECTOR (RESULT)
      
SUBROUTINE RTMSX(R,N,X,V)

implicit none

real*8    R(*)
integer*4 N
real*8    X(N)
real*8    V(N)
!local
integer*4 II,I,IK,K
real*8    SUM

II = 0
DO I=1,N
    SUM = 0.0
    II = II + I
    ! IK =(I,K)
    IK = II
    DO K = I,N
    SUM = SUM + R(IK)*X(K)
    IK=IK+K
    enddo
    V(I) = SUM
enddo
      
RETURN
END