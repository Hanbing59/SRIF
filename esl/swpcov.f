C----------------------------------------------
      SUBROUTINE SWPCOV(P,N,NAMES,ICOL,JCOL)
C
C     SWAP COVARIANCE  MATRIX COLUMNS
C
C   ** INPUT
C
C     P(N*(N+1)/2)     COVARIANCE MATRIX, VECTOR STORED.
C
C     N                STATE SIZE.
C
C     NAMES            STATE NAMES, IF NAMES(1) = '*' IT IS ASSUMED
C                      THAT THERE ARE NO STATE NAMES.
C
C     ICOL,JCOL        COLUMNS TO BE SWAPPED.
C
C ** OUTPUT **
C
C      P(N*(N+1)/2)   PERMUTED COVARIANCE ARRAY.
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE SWPCOV IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *                     COPYRIGHT 1987                             *
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

C
      DIMENSION     P(*)
      CHARACTER*(*) NAMES(*)
      CHARACTER*17  TMPNAM
      INTEGER       N,ICOL,JCOL
C
      IF (ICOL .EQ. JCOL) RETURN
C                                NOTHING TO DO
      I  = MIN(ICOL,JCOL)
      J  = MAX(ICOL,JCOL)
      II = I*(I-1)/2
      JJ = J*(J-1)/2
C
      IF (NAMES(1) .NE. '*') THEN
         TMPNAM      = NAMES(I)
         NAMES(I)    = NAMES(J)
         NAMES(J)    = TMPNAM
      ENDIF
C
C    INTERCHANGE P(K,I) AND P(K,J) K = 1, I (COL)
C
      DO 100 K = 1, I
         TMP     = P(II+K)
         P(II+K) = P(JJ+K)
         P(JJ+K) = TMP
100   CONTINUE
C
C    INTERCHANGE P(I,K) AND P(J,K) K = J, N (ROW)
C
      IK = JJ+I
      JK = JJ+J
      DO 200 K = J, N
         TMP   = P(IK)
         P(IK) = P(JK)
         P(JK) = TMP
         IK    = IK + K
         JK    = JK + K
200   CONTINUE
C
      IK = II+I+I
      II = II+I
      JJ = JJ+I
      JK = JJ+1
C
C    INTERCHANGE ROW ELEMENTS P(I,K) K = I+1 , J-1
C                COL ELEMENTS P(K,J) K = I+1 , J-1
C     I.E. P(I,J)<>P(J,I)
C
      DO 300 K = I+1, J-1
         TMP   = P(IK)
         P(IK) = P(JK)
         P(JK) = TMP
         IK    = IK+K
         JK    = JK+1
300   CONTINUE
C                      FIXUP, FOR LAST ELEMENT
      TMP     = P(JJ)
      P(JJ)   = P(II)
      P(II)   = TMP
      RETURN
      END
