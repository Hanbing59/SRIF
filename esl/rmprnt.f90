!! Purpose    : PRINT A RECTANGULAR VECTOR STORED MATRIX
!!
!! Parameters :
!!              A --- VECTOR STORED RECTANGULAR MATRIX, PRINT AS A(M,N)
!!             IA --- ROW DIMENSION OF A (AS DEFINED IN THE CALLING PROGRAM) 
!!             NR --- NUMBER OF ROWS OF A TO BE PRINTED
!!             NC --- NUMBER OF COLS OF A TO BE PRINTED
!!           NAMR --- NAMES OF PARAMETERS FOR ROWS OF A
!!           NAMC --- NAMES OF PARAMETERS FOR COLS OF A
!!                    NAMR AND NAMC MAY BE OF ANY CHARACTER TYPE (CHARACTER*K)
!!                    BUT ONLY THE FIRST 16 ARE SIGNIFICANT. NAMR AND NAMC
!!                    MAY BE OF DIFFERENT LENGTHS. AESTHETICALLY  PLEASING
!!                    RESULTS ARE OBTAINED WHEN NAMC IS CHARACTER*16 AND  
!!                    (ALL 16 BYTES ARE USED FOR EACH NAMC) OR (EACH NAMC 
!!                    IS CENTERED WITHIN ITS 16 BYTE FIELD).
!!           TEXT --- AN ARRAY OF ASCII CHARACTERS TO BE PRINTED AS A TITLE  
!!                    PRECEDING THE MATRIX.
!!           CODE --- A LIST OF PARAMETERS, WHICH ARE INTERPRETED IN A POSITION 
!!                    INDEPENDENT FASHION. ANY CHARACTERS NOT DEFINED IN THIS 
!!                    LIST WILL HAVE NO EFFECT; BUT OTHER FEATURES MAY BE   
!!                    DEFINED AT A LATER DATE, SO EXTRANEOUS CHARACTERS SHOULD 
!!                    BE AVOIDED.    
!!                    "P", A PAGE EJECT WILL PRECEDE PRINTING   
!!                    "C", CONDENSED PRINT (I.E. LESS PRECISION, MORE COLUMNS 
!!                         PER ROW). THE DEFAULT (' ') CAUSES PRINTING TO BEGIN  
!!                         AT THE NEXT LINE AND DISPLAYS 6 COLUMNS WITH 8 DIGITS 
!!                         OF PRECISION. CONDENSED MODE CONSISTS OF 10 COLUMNS 
!!                         OF 4 DIGIT PRECISION. IF THE NAMC CHARACTER NAMES ARE  
!!                         LONGER THAN 8 CHARACTERS LONG OTHER CODE PARAMETERS 
!!                         CONTROLS HOW EACH NAME IS TO BE PRINTED (SINCE ONLY 
!!                         8 CHARACTERS CAN APPEAR PER LINE).
!!                    "BK='C'", WHERE 'C' IS ANY CHARACTER, INCLUDING SPACE. 
!!                              CONTROL FOR HANDLING NAMC NAMES WITH LENGTH 
!!                              LONGER THAN 8. 'C' IS USED TO REQUEST TO PUT 
!!                              THE REMAINING CHARACTERS ON THE NEXT LINE.   
!!                              THE DEFAULT IS TO TREAT A BLANK AS THE REQUEST 
!!                              CHARACTER. NOTE: THE REQUIRED SYNTAX IS BK=THING, 
!!                              BLANKS ARE SIGNIFICANT. IT SHOULD ALSO BE NOTED 
!!                              THAT THIS OPTION IS IRRELEVANT UNLESS THE 
!!                              CONDENSED PRINT OPTION IS IN EFFECT.  
!! 
!!                    EXAMPLES: N = 3, CHARACTER*16 NAMES ='X-POSITION',
!!                              'X VELOCITY','X_ACCELERATION' 
!!                              CODE = 'PC' WILL PRODUCE A NEW PAGE AND CONDENSED PRINT
!!                              CODE = 'C'  WILL PRODUCE CONDENSED PRINT ONLY 
!!                              THE COLUMN NAMES WILL BE: 
!!                              X-POSITI   X            X_ACCELE  
!!                              ON         VELOCITY     RATION
!!                              SINCE "BK=" WAS NOT SPECIFIED, SPACE TRIGGERED 
!!                              THE BREAK   
!!                              IF CODE = 'C BK=_' THE RESULTING COLUMN NAMES WOULD BE  
!!                              X-POSITI  X VELOCI      X 
!!                              ON        TY            ACCELERA 
!!                              SINCE THE "BREAK" CHARACTER IS NOT PRINTED (THE FIRST 
!!                              TIME IT OCCURS) IT IS POSSIBLE (AND MEANINGFUL) TO USE 
!!                              CHARACTER*17 NAMES, IF ONE OF THE CHARACTERS IS THE 
!!                              "BREAK" CHARACTER.
!!         LBLROW --- LOGICAL FLAG, IF TRUE PRINT ROW NAMES
!!                                  IF FALSE PRINT NUMBERED ROWS   
!!                                  IF FALSE AND THERE IS ONLY ONE MATRIX ROW TO 
!!                                  PRINT THE NUMBER LABEL IS SUPPRESSED.  
!!         LBLCOL --- LOGICAL FLAG, IF TRUE PRINT COL NAMES
!!                                  IF FALSE PRINT NUMBERED COLS.  
!!         PRTNAM --- LABEL ASSOCIATED WITH REAL NUMBER, RNUM, TO BE PRINTED 
!!                    FOLLOWING THE TITLE. IF .EQ. ' ' THEN RNUM WILL NOT BE PRINTED.   
!!           RNUM --- NUMBER TO BE PRINTED  
!!         IOUNIT --- UNIT NUMBER OF OUTPUT DEVICE. 
!!                    IF IOUNIT .EQ. 0 THEN NO OUTPUT IS PRODUCED  
      
SUBROUTINE RMPRNT(A,IA,NR,NC,NAMR,NAMC,TEXT,CODE,LBLROW,LBLCOL,PRTNAM,RNUM,IOUNIT)  
 
implicit none

real*8       A(*)
integer*4    IA
integer*4    NR,NC
CHARACTER(*) NAMR(NR),NAMC(NC)
CHARACTER(*) TEXT
CHARACTER(*) CODE 
LOGICAL      LBLROW,LBLCOL
CHARACTER(*) PRTNAM
real*8       RNUM
integer*4    IOUNIT
!local
integer*4    i,j,m,l
integer*4    ic,ir,irr,ll,iptr,jEnd
LOGICAL      CONDEN,LSEC,LBLNK   
CHARACTER*8  TMPNAM(2,10)  
real*8       V(10)
CHARACTER*1  BREAK 
INTEGER      II(10)
CHARACTER*17 BLANKS   
DATA BLANKS /'                '/
      
IF(IOUNIT .EQ. 0) RETURN

CONDEN=INDEX(CODE,'C').NE.0 .OR. INDEX(CODE,'c').NE.0   
IF((INDEX(CODE,'P').NE.0) .OR. (INDEX(CODE,'p').NE.0)) THEN
    WRITE (IOUNIT,'(1H1)') 
ELSE  
    WRITE (IOUNIT,*)   
    WRITE (IOUNIT,*)   
ENDIF
I = INDEX(CODE,'BK=') 
J = INDEX(CODE,'BK=') 
I = MAX(I,J)  
IF(I.EQ.0) THEN 
    BREAK = ' '
ELSE  
    BREAK = CODE(I+3:I+4)  
ENDIF 

IF(PRTNAM.EQ.' ') THEN
    WRITE(IOUNIT,*) TEXT   
ELSE  
9   FORMAT(' ',A,'     ',A,' = ',1PD14.7) 
    WRITE(IOUNIT,9) TEXT,PRTNAM,RNUM   
ENDIF

IF (NR.EQ.0 .OR. NC.EQ.0) RETURN  

M   = NC  
IC  =  1  
IRR =  0  

10  L=M   
    IF(CONDEN) THEN
        IF (M.GT.10) L=10  
    ELSE
        IF (M.GT.6) L=6
    ENDIF 
    LL=IC+L-1 
    DO I=1,L
        II(I)=IC+I-1
    enddo
    ! IF COLUMNS NAMES WERE REQUESTED, AND CONDENSED PRINT THEN WE NEED TO COPY 
    ! NAMC INTO TEMP 8 CHAR ARRAY  
    WRITE(IOUNIT,*) ' '   
    IF(LBLCOL) THEN  
        IF(CONDEN) THEN   
            LSEC  = .FALSE. 
            LBLNK = .FALSE. 
            IPTR  = 0   
            DO I = IC, LL
                J = INDEX(NAMC(I),BREAK)-1   
                J = MIN(J,8) 
                IF (J .LE. 0) THEN   
                    J = LEN(NAMC(I))  
                    IF (J .GT. 8) J = J/2 
                ELSE 
                    LBLNK = .TRUE.
                ENDIF
                IPTR = IPTR + 1  
                TMPNAM(1,IPTR) = NAMC(I)(1:J)
                JEND = MIN(17,LEN(NAMC(I)))  
                TMPNAM(2,IPTR) = ' ' 
                IF(JEND .GE. J) THEN
                    IF(LBLNK .OR. (J .GE. JEND/2))       J = J+1 
                    IF(INDEX(NAMC(I)(J:J),BREAK) .NE. 0) J = J+1 
                    IF(JEND .GE. J) THEN  
                        LSEC = .TRUE.  
                        TMPNAM(2,IPTR) = NAMC(I)(J:JEND)   
                    ENDIF 
                ENDIF
            enddo
            WRITE(IOUNIT,*) ' ' 
            WRITE(IOUNIT,540) (TMPNAM(1,I),I=1,IPTR)
            IF (LSEC) WRITE(IOUNIT,540) (TMPNAM(2,I),I=1,IPTR)  
        ELSE   
            WRITE(IOUNIT,200) (NAMC(I), I=IC,LL)
        ENDIF  
    ELSE  
        IF(CONDEN) THEN   
            WRITE (IOUNIT,210) (II(I),I=1,L)
        ELSE   
            WRITE (IOUNIT,310) (II(I),I=1,L)
        ENDIF  
    ENDIF 
    WRITE(IOUNIT,*) ' '   
200 FORMAT (' ',17X,6(A16,3X))
210 FORMAT (' ',17X,10(I6,5X))
310 FORMAT (' ',18X,6(I6,13X))
540 FORMAT (' ',16X,10(2X,A8,1X)) 
    
220 FORMAT (' ',A16,1P10D11.3)
230 FORMAT (' ',I6,10X,1P10D11.3) 
320 FORMAT (' ',A16,1P,6(2X,D14.7,3X))
330 FORMAT (' ',I6,8X,1P6D19.7)
    DO I=1,NR  
        IR=IRR+I   
        DO J=1,L
            V(J)=A(IR)  
            IR=IR+IA
        enddo 
        ! IF ONLY ONE ROW IS TO BE PRINTED, AND NO NAMES ARE SELECTED PRINT 
        ! BLANKS INSTEAD OF THE NUMBER 1.    
        IF(.NOT.LBLROW .AND. NR.EQ.1) THEN   
            IF (CONDEN) THEN   
                WRITE (IOUNIT,220) BLANKS,(V(J),J=1,L)  
            ELSE   
                WRITE (IOUNIT,320) BLANKS,(V(J),J=1,L)  
            ENDIF  
        ELSE
            ! PRINT CURRENT ROW WITH OR WITHOUT SIDE HDRS DEPENDING ON LBLROW
            IF (CONDEN) THEN   
                IF (LBLROW) THEN
                    WRITE (IOUNIT,220) NAMR(I),(V(J),J=1,L)  
                ELSE
                    WRITE (IOUNIT,230) I,(V(J),J=1,L)
                ENDIF   
            ELSE   
                IF (LBLROW) THEN
                    WRITE (IOUNIT,320) NAMR(I),(V(J),J=1,L)  
                ELSE
                    WRITE (IOUNIT,330) I,(V(J),J=1,L)
                ENDIF   
            ENDIF  
        ENDIF
    enddo
    M = M-L 
    IC = IC+L
    IRR = IRR+L*IA
    IF(M.GT.0) GOTO 10 
    
RETURN
END   