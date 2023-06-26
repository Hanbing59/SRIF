!! Purpose    : REARRANGE THE ROWS AND COLUMNS OF MATRIX C which is IN NAMIN 
!!              ORDER AND PUT THE RESULT IN C(LOUT,LOUT) IN NAMOUT ORDER. 
!!              ZERO COLUMNS AND ROWS ARE ASSOCIATED WITH parameters DEFINED 
!!              in NAMOUT but not CONTAINED IN NAMIN.
!! Parameters :
!!                  C --- INPUT MATRIX
!!              IMAXC --- ROW DIMENSION OF C, IMAXC.GE.L=MAX(LIN,LOUT)
!!                LIN --- NO. OF PARAMETER NAMES ASSOCIATED WITH THE INPUT C
!!              NAMIN --- PARAMETER NAMES ASSOCIATED WITH C ON INPUT.
!!                        (ONLY THE FIRST LIN ENTRIES APPLY TO THE
!!                        INPUT C) (NOTE: DIM(NAMIN) .GE. DIM(NAMOUT))
!!                        UP TO 16 CHARACTER NAMES ARE ALLOWED.
!!               LOUT --- NO. OF PARAMETER NAMES ASSOCIATED WITH THE OUTPUT C
!!             NAMOUT --- PARAMETER NAMES ASSOCIATED WITH THE OUTPUT C
!!                        (ALSO UP TO 16 CHARACTER NAMES)

SUBROUTINE C2C(C,IMAXC,LIN,NAMIN,LOUT,NAMOUT)

implicit none
      
real*8       C(IMAXC,LIN)
integer*4    IMAXC,LIN
CHARACTER(*) NAMIN(*)
integer*4    LOUT
CHARACTER(*) NAMOUT(LOUT)
!local
integer*4    I,J,K,L
real*8       xtmp
real*8       ZERO
CHARACTER*17 ctmp
CHARACTER*17 BLANKS
DATA BLANKS/'                '/

ZERO=0.0D0
L=MAX0(LIN,LOUT)

DO K=LIN+1,L
    !BLANK REMAINING NAMIN LOCNS
    NAMIN(K)=BLANKS
enddo

DO J=1,LOUT
    DO I=1,L
        IF(NAMIN(I).EQ.NAMOUT(J)) exit
    enddo
    if(i.gt.l .or. i.eq.j) cycle
    !INTERCHANGE COLUMNS I AND J
    DO K=1,L
        xtmp=C(K,J)
        C(K,J)=C(K,I)
        C(K,I)=xtmp
    enddo
    !INTERCHANGE ROWS I AND J
    DO K=1,L
        xtmp=C(J,K)
        C(J,K)=C(I,K)
        C(I,K)=xtmp
    enddo
    !INTERCHANGE LABELS I AND J
    ctmp=NAMIN(I)
    NAMIN(I)=NAMIN(J)
    NAMIN(J)=ctmp
enddo

!FIND NAMOUT NAMES NOT IN NAMIN AND SET CORRESPONDING
!ROWS AND COLUMNS TO ZERO
DO J=1,LOUT
    DO I=1,L
        IF(NAMIN(I).EQ.NAMOUT(J)) exit
    enddo
    if(i.le.l) cycle
    !Names in NAMOUT but not in NAMIN
    DO K=1,LOUT
        C(J,K)=ZERO
        C(K,J)=ZERO
    enddo
enddo
RETURN
END