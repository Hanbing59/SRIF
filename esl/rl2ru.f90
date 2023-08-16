!! Purpose    : COPY THE LOWER RIGHT PORTION OF A VECTOR STORED UPPER
!!              TRIANGULAR MATRIX, RIN, INTO THE UPPER LEFT CORNER BLOCK OF
!!              ROUT, A VECTOR STORED UPPER TRIANGULAR ARRAY, I.E.
!!              ROWS AND COLUMNS L THROUGH N OF RIN ARE COPIED INTO ROWS AND
!!              COLUMNS 1 THROUGH N-L+1 OF ROUT, THE NAMES ASSOCIATED WITH
!!              THE RIN COLUMNS ARE TRANSFERRED TO ROUT. IF THE COMPILER ALLOWS IT,
!!              ROUT CAN SHARE STORAGE WITH RIN.
!! Parameters :
!!              RIN --- VECTOR STORED UPPER TRIANGULAR MATRIX
!!                N --- DIMENSION OF RIN
!!            NAMIN --- NAME LABELS ASSOCIATED WITH RIN UP TO CHARACTER*16
!!           FSTNAM --- STARTING NAME OF THE GROUP OF CONTIGUOUS ENTRIES
!!                      THAT ARE TO BE TRANSFERRED. IF THE FIRST CHARACTER
!!                      OF FSTNAM IS "#" THEN FSTNAM IS INTERPRETED AS THE
!!                      COLUMN NUMBER TO START COPYING FROM.
!!             ROUT --- THE OUTPUT COPIED PORTION OF RIN;
!!                      ONLY THIS PORTION OF ROUT IS INVOLVED WHERE L IS THE
!!                      COLUMN INDEX ASSOCIATED WITH FSTNAM AND NOUT = N-L+1
!!           NAMOUT --- THE OUTPUT COPIED PORTION OF NAMIN (UP TO CHARACTER*16)
!!             IERR --- ERROR RETURN
!!                      0, NORMAL RETURN
!!                      1, FSTNAM IS NOT INCLUDED IN THE INPUT LIST.
!!                         IN THIS CASE, ROUT IS LEFT UNALTERED.
!!

SUBROUTINE RL2RU(RIN,N,NAMIN,FSTNAM,ROUT,NAMOUT,IERR)

IMPLICIT NONE

real*8       RIN(*)
integer*4    N
CHARACTER(*) NAMIN(N)
CHARACTER(*) FSTNAM
real*8       ROUT(*)
CHARACTER(*) NAMOUT(*)
integer*4    IERR
!local
CHARACTER*20 CTEMP
LOGICAL      NAMED
integer*4    J,L,IJS,LM1,IJA,KOLA,KOL,K,IR

IF(FSTNAM(1:1).EQ.'#') THEN
    DO J=LEN(FSTNAM),1,-1
        IF(FSTNAM(J:J).NE.' ') exit
    enddo
    WRITE(CTEMP,'(''(1X,I'',I9.9,'')'')') J-1
    READ(FSTNAM,CTEMP) J
    NAMED=.FALSE.
else
    NAMED = .TRUE.
    DO J=1,N
        IF(NAMIN(J).EQ.FSTNAM) exit
    enddo
    if(J.gt.N) then
        IERR = 1
        RETURN
    endif
ENDIF

L = J
IJS = L*(L+1)/2-1
LM1 = L-1
IJA = 0
KOLA = 0
DO KOL = L,N
    KOLA = KOLA + 1
    K = KOL-LM1
    IF(NAMED) NAMOUT(K) = NAMIN(KOL)
    DO IR = 1, KOLA
        IJA = IJA+1
        K = IJS+IR
        ROUT(IJA) = RIN(K)
    enddo
    IJS = IJS + KOL
enddo
IERR = 0

RETURN
END