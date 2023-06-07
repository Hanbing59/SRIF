      SUBROUTINE RL2RU(RIN,N,NAMIN,FSTNAM,ROUT,NAMOUT,IERR)
C
C     TO COPY THE LOWER RIGHT PORTION OF A VECTOR STORED
C     UPPER TRIANGULAR MATRIX, RIN, INTO THE UPPER LEFT
C     CORNER BLOCK OF ROUT, A VECTOR STORED UPPER TRIANGULAR
C     ARRAY, I.E. ROWS AND COLUMNS L THROUGH N OF RIN ARE
C     COPIED INTO ROWS AND COLUMNS 1 THROUGH N-L+1 OF ROUT,
C     THE NAMES ASSOCIATED WITH THE RIN COLUMNS ARE
C     TRANSFERRED TO ROUT. IF THE COMPILER ALLOWS IT,
C     ROUT CAN SHARE STORAGE WITH RIN.
C
C *** INPUTS:
C
C     RIN(N*(N+1)/2)     VECTOR STORED UPPER TRIANGULAR MATRIX
C
C     N                  DIMENSION OF RIN
C
C     NAMIN(N)           NAME LABELS ASSOCIATED WITH RIN
C                        UP TO CHARACTER*16
C
C     FSTNAM             STARTING NAME OF THE GROUP OF
C                        CONTIGUOUS ENTRIES THAT ARE TO
C                        BE TRANSFERRED. IF THE FIRST CHARACTER
C                        OF FSTNAM IS "#" THEN FSTNAM IS
C                        INTERPRETED AS THE COLUMN NUMBER TO
C                        START COPYING FROM.
C
C
C
C *** OUTPUT:
C
C     ROUT(NOUT*(NOUT+1)/2)   THE OUTPUT COPIED PORTION OF RIN;
C                        ONLY THIS PORTION OF ROUT IS INVOLVED
C                        WHERE L IS THE COLUMN INDEX ASSOCIATED
C                        WITH FSTNAM AND NOUT = N-L+1
C
C     NAMOUT(NOUT)       THE OUTPUT COPIED PORTION OF NAMIN
C                        (UP TO CHARACTER*16)
C
C
C     IERR               ERROR RETURN
C
C                          .EQ.  0      NORMAL RETURN
C
C                          .EQ.  1      FSTNAM IS NOT INCLUDED IN
C                                       THE INPUT LIST. IN THIS CAS
C                                       ROUT IS LEFT UNALTERED.
C
C
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              APRIL    1, 1982
C
C              REVISED 20 AUGUST 1982  VERSION 0.00.01
C                                      VERSION 1.00.00
C              REVISED 20    MAY 1985  VERSION 2.00.00
C              REVISED 14  AUGUST1986  VERSION 2.01.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE RL2RU  IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *         COPYRIGHT  1982, 1983,  1984, 1985, 1986               *
C *          FACTORIZED ESTIMATION APPLICATIONS INC.               *
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

      DIMENSION RIN(*), ROUT(*)
      CHARACTER*(*) NAMIN(N), NAMOUT(*), FSTNAM
      CHARACTER*20  CTEMP
      LOGICAL   NAMED
C
      IF (FSTNAM(1:1) .EQ. '#') THEN
          DO 1 J = LEN(FSTNAM),1,-1
             IF (FSTNAM(J:J).NE.' ') GO TO 2
1         CONTINUE
2         CONTINUE
          WRITE(CTEMP,'(''(1X,I'',I9.9,'')'')') J-1
          READ(FSTNAM,CTEMP) J
          NAMED = .FALSE.
          GO TO 20
      END IF
      NAMED = .TRUE.
      DO 10 J=1, N
         IF (NAMIN(J) .EQ. FSTNAM) GO TO 20
  10  CONTINUE
C
C   RETURN WITH ERROR
C
      IERR = 1
      RETURN
  20  L    = J
      IJS  = L*(L+1)/2-1
      LM1  = L-1
      IJA  = 0
      KOLA = 0
      DO 40 KOL = L,N
         KOLA = KOLA + 1
         K    = KOL-LM1
         IF (NAMED) NAMOUT(K) = NAMIN(KOL)
         DO 30 IR     = 1, KOLA
            IJA       = IJA+1
            K         = IJS+IR
            ROUT(IJA) = RIN(K)
  30     CONTINUE
         IJS = IJS + KOL
  40  CONTINUE
C
C   NORMAL RETURN
C
      IERR = 0
      RETURN
      END
