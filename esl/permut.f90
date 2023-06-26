!! Purpose    : REARRANGE PARAMETERS OF A(IRA,JCIN), NAMIN ORDER TO
!!              A(IRA,JCOUT), NAMOUT ORDER. ZERO COLUMNS ARE INSERTED
!!              CORRESPONDING TO THE NEWLY DEFINED NAMES.
!! Parameters :
!!               A --- INPUT RECTANGULAR MATRIX, L=MAX(JCIN,JCOUT)
!!           IMAXA --- ROW DIMENSION OF A, IMAXA.GE.IRA
!!             IRA --- NUMBER OF ROWS OF A THAT ARE TO BE REARRANGED
!!            JCIN --- NUMBER OF PARAMETER NAMES ASSOCIATED WITH THE
!!                     INPUT A MATRIX
!!           NAMIN --- PARAMETER NAMES ASSOCIATED WITH A ON INPUT
!!                     (ONLY THE FIRST JCIN ENTRIES APPLY TO THE
!!                      INPUT A) NAMIN IS DESTROYED BY PERMUT.
!!                     (UP TO CHARACTER*16).
!!                     WHERE L = MAX(JCIN,JCOUT)
!!           JCOUT --- NUMBER OF PARAMETER NAMES ASSOCIATED WITH
!!                     THE OUTPUT A MATRIX
!!          NAMOUT --- PARAMETER NAMES ASSOCIATED WITH THE OUTPUT A
!!                     (UP TO CHARACTER*16)

SUBROUTINE PERMUT(A,IMAXA,IRA,JCIN,NAMIN,JCOUT,NAMOUT)

implicit none

real*8       A(IMAXA,*)
integer*4    IMAXA,IRA,JCIN
CHARACTER(*) NAMIN(*)
integer*4    JCOUT
CHARACTER(*) NAMOUT(JCOUT)
!local
integer*4    I,J,K,L,I2
real*8       ZERO,W
CHARACTER*17 NAMTMP,BLANKS
DATA BLANKS/'               '/

ZERO=0.0
L=MAX0(JCIN,JCOUT)
! LOCATIONS BEYOND JCIN NOT REQUIRED, SO WE ARE FREE TO BLANK THEM
IF(JCOUT.GT.JCIN) THEN
    DO J = (JCIN+1),JCOUT
        NAMIN(J) = BLANKS
    enddo
ENDIF
DO J=1,JCIN
    NAMTMP=NAMIN(J)
    DO K=1,JCOUT
        IF(NAMTMP.EQ.NAMOUT(K)) exit
    enddo
    if(K.gt.JCOUT) NAMIN(J) = BLANKS
enddo
DO J=1,JCOUT
    DO I=1,L
        IF(NAMIN(I).EQ.NAMOUT(J)) exit
    enddo
    if(I.gt.L) cycle
    IF(I.EQ.J) cycle
    IF(I.LT.J) THEN
        ! FIND A BLANK
        DO I2=1,L
            IF(NAMIN(I2).EQ.BLANKS) exit
        enddo
        if(I2.gt.L) then
            PRINT*,' PERMUT: FATAL ERROR COULD NOT FIND A BLANK'
            PRINT*,'  COLUMN TO TEMPORARILY STORE COLUMN ',J,' IN'
            ! BY CONSTRUCTION THIS CANNOT HAPPEN
            STOP
        endif
        DO K = 1, IRA
            ! COPY COL J TO COL I2 (BLANK)
            A(K,I2) = A(K,J)
            !COPY COL I TO COL J
            A(K,J)  = A(K,I)
        enddo
        NAMIN(I2)  = NAMIN(J)
        NAMIN(J)   = NAMIN(I)
        NAMIN(I)   = ' '
    ELSE
        !INTERCHANGE COLS I AND J
        DO K=1,IRA
            W=A(K,J)
            A(K,J)=A(K,I)
            A(K,I)=W
        enddo
        NAMTMP   = NAMIN(I)
        NAMIN(I) = NAMIN(J)
        NAMIN(J) = NAMTMP
    ENDIF
enddo
! REPEAT TO FILL NEW COLS W/ZEROS
DO J=1,JCOUT
    DO I=1, L
        IF(NAMIN(I).EQ.NAMOUT(J)) exit
    enddo
    if(I.gt.L) then
        DO K=1,IRA
            A(K,J)=ZERO
        enddo
    endif
enddo

RETURN
END