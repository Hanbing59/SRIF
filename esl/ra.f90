!! Purpose    : COMPUTE IRW ROWS OF THE PRODUCT W=R*A WHERE R IS UPPER 
!!              TRIANGULAR VECTOR SUBSCRIPTED AND OF DIMENSION N, A HAS 
!!              JCA COLUMNS AND IRA ROWS, THE TRAILING N-JCA COLUMNS ARE 
!!              IMPLICITY DEFINED AS ZERO. WHEN IRA.LT.N THE BOTTOM N-IRA 
!!              ROWS OF A ARE ASSUMED TO BE IMPLICITLY DEFINED AS THE 
!!              BOTTOM N-IRA ROWS OF THE N DIMENSION IDENTITY MATRIX.
!! Parameters :
!!               R --- UPPER TRIANGULAR VECTOR STORED INPUT MATRIX
!!               N --- DIMENSION OF R.
!!               A --- INPUT RECTANGULAR MATRIX
!!           IMAXA --- ROW DIMENSION OF A
!!             IRA --- NUMBER OF ROWS IN THE A MATRIX (IRA.LE.IMAXA AND 
!!                     IRA.LE.N)
!!             JCA --- NUMBER OF COLUMNS IN THE A MATRIX, JCA.LE.N
!!               W --- RESULTING RECTANGULAR MATRIX, = R*A, IF THE COMPILER 
!!                     ALLOWS IT W CAN OVERWRITE A
!!           IMAXW --- ROW DIMENSION OF W
!!             IRW --- NUMBER OF ROWS OF THE PRODUCT R*A THAT ARE COMPUTED 
!!                     (IRW.LE.IMAXW AND IRW.LE.N)

SUBROUTINE RA(R,N,A,IMAXA,IRA,JCA,W,IMAXW,IRW)

IMPLICIT none

real*8     R(*)
integer*4  N
real*8     A(IMAXA,JCA)
integer*4  IMAXA,IRA,JCA
real*8     W(IMAXW,N)
integer*4  IMAXW,IRW
!local
real*8     SUM
integer*4  I,J,K,II,IJ,IK

IJ=IRA*(IRA+1)/2
DO J=1,N
    II=0
    DO I=1,IRW
        II=II+I
        SUM=0.D0
        if(I.le.IRA .and. J.le.JCA) then
            IK=II
            DO K=I,IRA
                SUM=SUM+R(IK)*A(K,J)
                IK=IK+K
            enddo
        endif

        IF(J.gt.IRA .and. I.le.J) then
            SUM = SUM + R(IJ+I)
        endif
        W(I,J)=SUM
    enddo
    IF(J.GT.IRA) IJ=IJ+J
enddo

RETURN
END