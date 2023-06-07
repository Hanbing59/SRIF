      SUBROUTINE RU2RL (RIN,NAMIN,LSTNAM,ROUT,N,NAMOUT,IERR)
C
C     TO COPY THE UPPER LEFT PORTION OF A VECTOR STORED
C     UPPER TRIANGULAR MATRIX, RIN, INTO THE LOWER RIGHT
C     PARTITION CORNER OF VECTOR STORED TRIANGULAR MATRIX
C     ROUT, I.E. ROWS 1 THROUGH L AND COLUMNS 1 THROUGH
C     L OF R BECOME ROWS N-L+1 THROUGH N AND COLUMNS N-L+1
C     THROUGH N OF ROUT. CORRESPONDING NAMOUT LABELS ARE
C     TRANSFERRED TO NAMIN IF THE COMPILER ALLOWS IT ROUT MAY
C     SHARE STORAGE WITH RIN
C
C *** INPUTS:
C
C     RIN(L*(L+1)/2)     VECTOR STORED UPPER TRIANGULAR MATRIX,
C                        WHERE L IS THE INDEX THAT CORRESPONDS
C                        TO THE LABEL "LSTNAM".
C
C     NAMIN(L)           NAME LABELS ASSOCIATED WITH RIN.
C                        NAMIN(L), NAMOUT(N), AND LSTNAM
C                        ARE CHARACTER*(*).
C
C     LSTNAM             LAST NAME OF THE GROUP OF CONTIGUOUS
C                        ENTRIES, STARTING WITH THE FIRST,
C                        THAT ARE  TO BE MOVED. THIS NAME
C                        CORRESPONDS TO THE LTH VARIABLE OF
C                        THE INPUT STATE VECTOR. IF THE FIRST
C                        CHARACTER OF LSTNAM IS "#" THEN
C                        LSTNAM IS INTERPRETED AS THE COLUMN
C                        NUMBER TO START COPYING FROM.
C
C     N                  DIMENSION OF ROUT
C
C *** OUTPUTS:
C
C     ROUT(N*(N+1)/2     VECTOR STORED UPPER TRIANGULAR MATRIX
C                        (ONLY ENTRIES IN ROWS AND COLUMNS
C                         N-L+1 THROUGH N ARE CHANGED)
C
C     NAMOUT(N)          ONLY NAMES IN ENTRIES N-L+1 THROUGH
C                        N ARE CHANGED
C
C     IERR               ERROR FLAG
C                         .EQ. 0  NORMAL RETURN
C                         .EQ. 1  LABEL LSTNAM DOES NOT APPEAR
C                                 IN THE NAMIN  LIST
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23, 1981
C
C              REVISED  2 AUGUST 1982     VERSION 0.01.00
C              REVISED 13 AUGUST 1982     VERSION 0.01.01
C              REVISED 20 AUGUST 1982     VERSION 0.01.02
C                                         VERSION 1.00.00
C              REVISED  5 AUGUST    1984  VERSION 1.01.00
C              REVISED  5 SEPTEMBER 1984  VERSION 2.00.00
C              REVISED 20 MAY       1985  VERSION 2.01.00
C              REVISED 20 APRIL     1986  VERSION 2.02.00
C              REVISED 14 AUGUST    1986  VERSION 2.03.00
C
C *******************************************************************
C                                                                   *
C          SUBROUTINE RU2RL  IS A PART OF THE FEA. INC              *
C                ESTIMATION SUBROUTINE LIBRARY                      *
C                                                                   *
C             COPYRIGHT 1982, 1983, 1984, 1985, 1986                *
C            FACTORIZED ESTIMATION APPLICATIONS INC.                *
C *   This library is licensed under the Creative Commons          *
C *   Attribution 3.0 Unported License.                            *
C *   To view a copy of this license,                              *
C *   visit http://creativecommons.org/licenses/by/3.0/            *
C *   or send a letter to Creative Commons, 171 Second Street,     *
C *   Suite 300, San Francisco, California, 94105, USA.            *
C *                                                                *
C ******************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION  RIN(*), ROUT(*)
      CHARACTER*(*) NAMIN(*), NAMOUT(N), LSTNAM
      CHARACTER*20 CTEMP
      LOGICAL NAMED
C
      IF (LSTNAM(1:1) .EQ. '#') THEN
          DO 1 LSTPOS = LEN(LSTNAM),1,-1
             IF (LSTNAM(LSTPOS:LSTPOS).NE.' ') GO TO 2
1         CONTINUE
2         CONTINUE
          WRITE(CTEMP,'(''(1X,I'',I9.9,'')'')') LSTPOS-1
          READ(LSTNAM,CTEMP) LSTPOS
          NAMED = .FALSE.
          PRINT*,' CTEMP = ',CTEMP
          PRINT*,' LSTPOS = ',LSTPOS
          PRINT*,' LSTNAM = ',LSTNAM
          GO TO 20
      END IF
      NAMED = .TRUE.
      DO 10 LSTPOS=1, N
         IF (NAMIN(LSTPOS) .EQ. LSTNAM) GO TO 20
  10  CONTINUE
      IERR = 1
      RETURN
  20  L      = LSTPOS+1
      IJ     = L*(LSTPOS)/2
      IJAS   = N*(N+1)/2
      KOLOUT = N
      DO 40 KOLIN = LSTPOS,1,-1
         IF (NAMED) NAMOUT(KOLOUT) = NAMIN(KOLIN)
         DO 30 IJA    = IJAS,IJAS-KOLIN+1,-1
            ROUT(IJA) = RIN(IJ)
            IJ  = IJ-1
  30     CONTINUE
         IJAS   = IJAS-KOLOUT
         KOLOUT = KOLOUT - 1
  40  CONTINUE
      IERR = 0
      RETURN
      END
