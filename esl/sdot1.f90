!! Purpose    : RETURNS THE DOT PRODUCT OF SX AND SY.
!!               SDOT = SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
!!               WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!!               DEFINED IN A SIMILAR WAY USING INCY.
!!
!! Parameters :


DOUBLE PRECISION FUNCTION SDOT1(N,SX,INCX,SY,INCY)

IMPLICIT NONE

integer*4 N
REAL*8    SX(*)
integer*4 INCX
REAL*8    SY(*)
integer*4 INCY
!local
real*8    SUM
integer*4 IX,IY,I,M

SUM = 0.0
IF(N.LE.0) then
    SDOT1 = SUM
    return
endif

IF(INCX.ne.1 .or. INCY.ne.1) then
    !CODE FO UNEQUAL INCREMENTS OR EQUAL INCREMENTS NE 1
    IX = 1
    IY = 1
    IF(INCX.LT.0) IX = (-N+1)*INCX+1
    IF(INCY.LT.0) IY = (-N+1)*INCY+1
    DO I=1, N
        SUM = SUM + SX(IX)*SY(IY)
        IX = IX+INCX
        IY = IY+INCY
    enddo
    SDOT1 = SUM
    return
endif

!CODE FOR BOTH INCREMENTS = 1
M = MOD(N,5)
IF(M.ne.0) then
    DO I = 1,M
        SUM = SUM + SX(I)*SY(I)
    enddo
    IF(N.LT.5) then
        SDOT1 = SUM
        return
    endif
endif

DO I = M+1, N, 5
    SUM = SUM + SX(I)*SY(I) + SX(I+1)*SY(I+1) + SX(I+2)*SY(I+2) + &
        SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
enddo

SDOT1 = SUM
RETURN
END