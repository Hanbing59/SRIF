      SUBROUTINE U2U(U,JCIN,NAMIN,JCOUT,NAMOUT,SCRATCH)
C
C        SUBROUTINE TO REARRANGE THE ROWS AND COLUMNS OF UD MATRIX
C        U(N*(N+1)/2) IN NAMIN ORDER AND PUT THE RESULT IN
C        U(N2*(N2+1)/2) IN NAMOUT ORDER. ZERO COLUMNS AND ROWS ARE
C        ASSOCIATED WITH OUTPUT DEFINED NAMES THAT ARE NOT
C        CONTAINED IN NAMIN.
C
C *** INPUTS:
C
C    U(L*(L+1)/2) U-D MATRIX, VECTOR STORED. (OUTPUT OVERWRITES INPUT)
C                 L = MAX(JCIN,JCOUT)
C
C    JCIN         NUMBER OF STATES (I.E. NAMES) IN INPUT U-D MATRIX.
C
C    NAMIN(L)     INPUT NAMES, WHERE L = MAX(JCIN,JCOUT).
C
C    JCOUT        NUMBER OF STATES IN OUTPUT U-D MATRIX.
C
C    NAMOUT       OUTPUT NAMES
C
C
C *** OUTPUT:
C
C    U        PERMUTED U-D MATRIX CORRESPONDING TO NAMOUT
C
C          NOTES:
C                   - JCOUT AND JCIN CAN BE DIFFERENT, AND BOTH
C                     NAMIN AND NAMOUT
C                     CAN HAVE PARAMETERS NOT IN THE OTHER. PARAMETERS IN
C                     NAMOUT NOT IN NAMIN GET SET TO ZERO COVARIANCE.
C                     PARAMETERS IN NAMIN NOT IN NAMOUT GET COMPRESSED OUT
C                     USING RANK 1 UPDATES.
C
C                   - REPEATED NAMES IN EITHER NAME LIST (NAMIN,NAMOUT)
C                     SHOULD BE AVOIDED, AS THE RESULTS ARE NOT WELL
C                     DEFINED.
C
C
C    THE CONTENTS OF NAMIN ARE DESTROYED
C
C *** SCRATCH:
C
C    SCRATCH(4,L)   INTERNALLY COLUMN 1 =G, 2 = H, ELSE SCRATCH
C
C**************************************************************************
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN / STEVE LICTHEN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JANUARY 23, 1987
C
C       REQUIRED SUBROUTINES:
C
C       LISPOS <ESL>  FIND A STATE NAME IN A LIST
C       SWPUCL <ESL>  SWAP TWO COLUMNS OF A VECTOR STORED UD ARRAY
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE U2U    IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *                  COPYRIGHT 1982 - 1987                         *
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

      DIMENSION     U(*),SCRATCH(*)
      CHARACTER*(*) NAMIN(*),NAMOUT(*)
      INTEGER       JCIN,JCOUT
C                                    LOCAL VARS
      INTEGER       I,J,K,L
C
      L = MAX(JCIN,JCOUT)
C                         ZERO FIRST TWO SCRATCH COLUMNS
      DO 10 I=1,L*2
         SCRATCH(I)=0.0
10    CONTINUE
C
C     STEP #1 REMOVE NAMIN COLUMNS NOT APPEARING IN NAMOUT
C
      K   = JCIN
      N1  = JCIN
      NTR = N1*(N1+1)/2
      DO 20 J=JCIN,1,-1
         CALL LISPOS(NAMIN(J),NAMOUT,JCOUT,IPOS)
C
C        IF NAMIN(J) IS NOT IN NAMOUT, MOVE IT TO END...
C
         IF (IPOS .EQ. 0) THEN
            IF (JCIN .NE. J) CALL SWPUCL(U,K,NAMIN,J,N1,SCRATCH)
C
C        THEN DELETE ITS EFFECTS
C
            CALL RNK1(U,K-1,U(NTR),U(NTR-K+1),IERR)
            NAMIN(N1) = ' '
C                            REDUCE THE EFFECTIVE MATRIX SIZE
            K   = K-1
            N1  = N1-1
            NTR = N1*(N1+1)/2
            IF (IERR .GT. 0) THEN
               WRITE(6,*)' NEGATIVE COMPUTED DIAGONAL',IERR
              STOP
           ENDIF
         ENDIF
 20   CONTINUE
C
C     STEP #2 ADDEND COLUMNS OF ZEROS FOR EACH NAMOUT NOT IN NAMIN
C
      DO 30 J = 1, JCOUT
         CALL LISPOS(NAMOUT(J),NAMIN,JCIN,IPOS)
         IF (IPOS .EQ. 0) THEN
            N1 = N1+1
            NAMIN(N1) = NAMOUT(J)
            NTR = N1*(N1+1)/2
            DO 25 I = NTR-N1+1,NTR
               U(I) = 0.0
 25         CONTINUE
         ENDIF
 30   CONTINUE
C
C     STEP #3 NOW U (AND NAMIN) HAS ALL THE NAMOUT STATES, BUT
C             IN THE WRONG ORDER
 
      DO 35 I = 1, JCOUT
         CALL LISPOS(NAMOUT(I),NAMIN,JCOUT,IPOS)
         WRITE(6,*)' SWAP COLS ',I,IPOS,' IN/OUT ',NAMIN(I),NAMOUT(I)
         CALL SWPUCL(U,JCOUT,NAMIN,I,IPOS,SCRATCH)
35    CONTINUE
      RETURN
      END
