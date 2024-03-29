!! Purpose    : COPY THE UPPER LEFT PORTION OF A VECTOR STORED UPPER TRIANGULAR
!!              MATRIX, RIN, INTO THE LOWER RIGHT PARTITION CORNER OF VECTOR STORED
!!              TRIANGULAR MATRIX ROUT, I.E. ROWS 1 THROUGH L AND COLUMNS 1 THROUGH
!!              L OF R BECOME ROWS N-L+1 THROUGH N AND COLUMNS N-L+1 THROUGH N OF
!!              ROUT. CORRESPONDING NAMOUT LABELS ARE TRANSFERRED TO NAMIN IF THE
!!              COMPILER ALLOWS IT ROUT MAY SHARE STORAGE WITH RIN
!!
!! Parameters :
!!              RIN --- VECTOR STORED UPPER TRIANGULAR MATRIX, WHERE L IS THE
!!                      INDEX THAT CORRESPONDS TO THE LABEL "LSTNAM".
!!
!!            NAMIN --- NAME LABELS ASSOCIATED WITH RIN.
!!                      NAMIN(L), NAMOUT(N), AND LSTNAM ARE CHARACTER*(*).
!!           LSTNAM --- LAST NAME OF THE GROUP OF CONTIGUOUS ENTRIES, STARTING
!!                      WITH THE FIRST, THAT ARE  TO BE MOVED. THIS NAME
!!                      CORRESPONDS TO THE LTH VARIABLE OF THE INPUT STATE VECTOR.
!!                      IF THE FIRST CHARACTER OF LSTNAM IS "#" THEN LSTNAM IS
!!                      INTERPRETED AS THE COLUMN NUMBER TO START COPYING FROM.
!!             ROUT --- VECTOR STORED UPPER TRIANGULAR MATRIX (ONLY ENTRIES IN
!!                      ROWS AND COLUMNS N-L+1 THROUGH N ARE CHANGED)
!!                N --- DIMENSION OF ROUT
!!           NAMOUT --- ONLY NAMES IN ENTRIES N-L+1 THROUGH N ARE CHANGED
!!             IERR --- ERROR FLAG
!!                      0, NORMAL RETURN
!!                      1, LABEL LSTNAM DOES NOT APPEAR IN THE NAMIN LIST

SUBROUTINE RU2RL (RIN,NAMIN,LSTNAM,ROUT,N,NAMOUT,IERR)

IMPLICIT none

real*8       RIN(*)
CHARACTER(*) NAMIN(*)
CHARACTER(*) LSTNAM
real*8       ROUT(*)
integer*4    N
CHARACTER(*) NAMOUT(N)
integer*4    IERR
!local
integer*4    LSTPOS,L,IJ,IJAS,KOLOUT,KOLIN,IJA
CHARACTER*20 CTEMP
LOGICAL      NAMED

IF(LSTNAM(1:1).EQ.'#') THEN
    DO LSTPOS = LEN(LSTNAM),1,-1
        IF(LSTNAM(LSTPOS:LSTPOS).NE.' ') exit
    enddo
    WRITE(CTEMP,'(''(1X,I'',I9.9,'')'')') LSTPOS-1
    READ(LSTNAM,CTEMP) LSTPOS
    NAMED = .FALSE.
    PRINT*,' CTEMP = ',CTEMP
    PRINT*,' LSTPOS = ',LSTPOS
    PRINT*,' LSTNAM = ',LSTNAM
else
    NAMED = .TRUE.
    DO LSTPOS=1, N
        IF(NAMIN(LSTPOS).EQ.LSTNAM) exit
    enddo
    if(LSTPOS.gt.N) then
        IERR = 1
        RETURN
    endif
ENDIF
L = LSTPOS+1
IJ = L*(LSTPOS)/2
IJAS = N*(N+1)/2
KOLOUT = N
DO KOLIN = LSTPOS,1,-1
    IF(NAMED) NAMOUT(KOLOUT) = NAMIN(KOLIN)
    DO IJA = IJAS,IJAS-KOLIN+1,-1
        ROUT(IJA) = RIN(IJ)
        IJ  = IJ-1
    enddo
    IJAS   = IJAS-KOLOUT
    KOLOUT = KOLOUT - 1
enddo
IERR = 0
RETURN
END
