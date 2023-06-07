	SUBROUTINE STPHIU(PHIEL,INDX,NROW,JCOL,U,N)
C
C     TO COMPUTE U:=PHI*U WHERE PHI IS UNIT UPPER TRIANGULAR
C     AND SPARSE. ONLY THE NON-ZERO (AND NON-DIAGONAL) PHI
C     ELEMENTS ARE EXPLICITLY EMPLOYED.
C
C     U IS VECTOR STORED, UPPER TRIANGULAR WITH IMPLICITLY
C     DEFINED UNIT DIAGONAL ELEMENTS.
C
C     *** INPUTS:
C
C     PHIEL(NPHI)      VALUES OF THE NONZERO ELEMENTS OF
C                      PHI, NOT INCLUDING THE UNIT DIAGONAL
C                      ELEMENTS. ELEMENTS ARE STORED ROW-WISE
C                      VIZ. FIRST THE NONZERO ELEMENTS IN
C                      ROW #1 (NOT INCLUDING PHI(1,1)).
C                      THESE ARE FOLLOWED BY THE NONZERO ELEMENTS
C                      IN ROW 2 (IF THERE ARE ANY) AND NOT
C                      INCLUDING PHI(2,2). THE TOTAL NUMBER OF
C                      ELEMENTS, NPHI IS NOT EXPLICTLY USED.
C
C
C     INDX(2,NROW)     INDX(1,1)=ROW NUMBER OF FIRST ROW WITH
C                                NONZERO, NON-DIAGONAL ENTRIES.
C
C                      INDX(1,2)=ROW NUMBER OF THE SECOND ROW WITH
C                                NONZERO, NON-DIAGONAL ENTRIES.
C
C                      ETC.      ROW ORDERING MUST BE MONOTONIC
C
C                      INDX(2,1)=NUMBER OF NONZERO ELEMENTS
C                                (NOT INCLUDING THE DIAGONAL)
C                                 IN ROW INDX(1,1)
C
C                      INDX(2,2)=NUMBER OF NONZERO ELEMENTS
C                                (NOT INCLUDING THE DIAGONAL)
C                                 IN ROW INDX(1,2)
C
C     ETC.             NOTE THAT THE ROW INDICES MUST BE
C                      INCREASING IN THIS CODE FORMULATION.
C                      IT IS NOT REQUIRED THAT THE COLUMN INDICES
C                      BE MONOTONIC, SEE EXAMPLE.
C
C
C     NROW             THE NUMBER OF ROWS WITH NONZERO ENTRIES,
C                      NOT INCLUDING THE DIAGONALS.
C
C     JCOL(NPHI)       COLUMN LOCATIONS OF THE ELEMENTS IN THE
C                      VECTOR PHIEL.
C
C     U(N*(N+1)/2)     UPPER TRIANGULAR U MATRIX, VECTOR STORED
C                      WITH IMPLICIT UNIT DIAGONAL ENTRIES. THE
C                      OUTPUT U OVERWRITES THE INPUT U
C
C     N                COLUMN DIMENSION OF U
C
C
C     *** OUTPUT:
C
C     U(N*(N+1)/2)     UPPER TRIANGULAR U MATRIX, VECTOR STORED
C                      WITH IMPLICIT UNIT DIAGONAL ENTRIES. THE
C                      OUTPUT U OVERWRITES THE INPUT U.
C
C          EXAMPLE TO ILLUSTRATE PHI ELEMENT LABELING
C
C
C                   1.   0.  0.    5.0   0.     2.
C
C                   0.   1.  3.0   0.    4.     5.
C
C                   0.   0.  1.0   0.    0.     0.
C    PHI =
C                   0.   0.  0.    1.    0.     1.5
C
C                   0.   0.  0.    0.    1.    -2.
C
C                   0.   0.  0.    0.    0.     1.
C
C
C
C     PHIEL(7) = (5. , 2. , 5. , 3. , 4. , 1.5 , -2.)
C
C     INDX(1,1)= 1 , INDX(1,2)= 2 ,  INDX(1,3)= 4 ,  INDX(1,4)= 5
C     INDX(2,1)= 2 , INDX(2,2)= 3 ,  INDX(2,3)= 1 ,  INDX(2,4)= 1
C
C     JCOL(7)  = (4,6,6,3,5,6,6)
C
C     NROW     = 4
C
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              DECEMBER 28, 1981
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C              REVISED 02      JAN 1997 removed recursion on L for parallel
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE STPHIU IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *             COPYRIGHT 1982, 1983, 1984, 1985                   *
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

      PARAMETER( LCHUNK = 19 )
      DIMENSION PHIEL(*), INDX(2,NROW), JCOL(*), U(*)
C
      IF ( (N .EQ. 1) .OR. (NROW .EQ. 0)) RETURN
C
      NM1   = N-1
C Is this a case where NELEM is always 1, and all in same column?
      JFLG = 1
      IF (INDX(2,1).NE.1) JFLG = 0
      KOL = JCOL(1)
      DO 5 I=2,NROW
        IF(INDX(2,I).NE.1.OR.JCOL(I).NE.KOL) JFLG=0
 5    CONTINUE
      IF(JFLG .EQ. 1) THEN
C no need for recursive definition of jstart; may operate down columns
        KOL = JCOL(1)
        KOLDIAG = KOL*(KOL-1)/2
CDIR$ IVDEP
        DO 10 I = 1,NROW
C take care of first column of U, that uses diagonal unit U
          IK = INDX(1,I) + KOLDIAG
          U(IK) = U(IK) + PHIEL(I)
10      CONTINUE
CDIR$ IVDEP
        DO 20 K = KOL,NM1
          KDIAG = K*(K+1)/2
          DO 30 I=1,NROW
             IR = INDX(1,I)
             IK = IR + KDIAG
             JK = KOL + IK - IR
             U(IK) = U(IK) + PHIEL(I)*U(JK)
30        CONTINUE
20      CONTINUE

      ELSE

        ISTRT = 0
        DO 40 I=1, NROW
         IR = INDX(1,I)
C                                  IR = ROW INDEX
         NELEM = INDX(2,I)
C                                  NELEM = # ELEMENTS IN ROW IR
         DO 50 JC=1, NELEM
            L = ISTRT + JC
            A = PHIEL(L)
            KOL = JCOL(L)
            JK = KOL * (KOL-1)/2
            IK = JK+IR
            U(IK) = U(IK)+A
            IF (KOL .EQ. N) GO TO 50
            DO 60 K=KOL, NM1
               IK = IR + K*(K+1)/2
               JK = KOL + IK - IR
               U(IK)=U(IK)+A*U(JK)
  60        CONTINUE
  50     CONTINUE
         ISTRT = ISTRT+NELEM
  40    CONTINUE
      ENDIF
      RETURN
      END
