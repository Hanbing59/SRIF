      SUBROUTINE LISPOS(TARGET, LIST, LGTH,IPOS)
      
C *   This library is licensed under the Creative Commons          *
C *   Attribution 3.0 Unported License.                            *
C *   To view a copy of this license,                              *
C *   visit http://creativecommons.org/licenses/by/3.0/            *
C *   or send a letter to Creative Commons, 171 Second Street,     *
C *   Suite 300, San Francisco, California, 94105, USA.            *
C *                                                                *
C ******************************************************************
C
C        LISPOS PERFORMS A LINEAR SEARCH OF LIST FOR THE ENTRY
C        WHICH IS EQUAL TO TARGET.  IF IT IS FOUND, IPOS CONTAINS
C        THE INDEX OF THE ITEM IN LIST, OTHERWISE 0.
C
C  ** INPUT:
C
      INTEGER LGTH
C                                  ! LENGTH OF LIST
      CHARACTER*(*) TARGET
C                                  ! ITEM TO FIND IN LIST
      CHARACTER*(*) LIST(*)
C                                  ! LIST OF ITEMS
C
C  ** OUTPUT:
C
      INTEGER IPOS
C                                  ! LOCATION OF TARGET
C
      DO 210 I = 1, LGTH
         IF (TARGET .EQ. LIST(I)) THEN
            IPOS = I
            RETURN
         ENDIF
210   CONTINUE
      IPOS = 0
      RETURN
      END
