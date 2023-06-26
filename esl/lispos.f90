!! Purpose    : PERFORMS A LINEAR SEARCH OF a LIST FOR a given TARGET. 
!!              IF IT IS FOUND, IPOS CONTAINS THE INDEX OF THE 
!!              ITEM IN LIST, OTHERWISE 0.
!! Parameters : 
!!              TARGET --- ITEM TO FIND IN LIST
!!                LIST --- LIST OF ITEMS
!!                LGTH --- LENGTH OF LIST
!!                IPOS --- LOCATION OF the TARGET


SUBROUTINE LISPOS(TARGET, LIST, LGTH,IPOS)

CHARACTER(*) TARGET
CHARACTER(*) LIST(*)
INTEGER*4    LGTH
INTEGER*4    IPOS
!local
integer*4    I

DO I=1, LGTH
    IF(TARGET.EQ.LIST(I)) THEN
        IPOS = I
        RETURN
    ENDIF
enddo
IPOS = 0

RETURN
END