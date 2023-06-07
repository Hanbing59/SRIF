      SUBROUTINE UCOLRDOLD(U,N,IP,BX,IBX,EM,Q,IERR)
C
C     COLORED NOISE U-D COVARIANCE FACTOR UPDATING
C     U*D*(U**T)-OUTPUT=PHI*U*D*(U**T)*(PHI**T)+Q*V*V**T
C     WHERE:
C       PHI = IDENTITY EXCEPT FOR IP-TH DIAGONAL, WHICH IS EM.
C       V**T = (BX**T,0**T,1,0**T)
C
C     THE ALGORITHM USED IS THE BIERMAN-THORNTON ONE COMPONENT
C     AT-A-TIME UPDATE. CF.BIERMAN "FACTORIZATION METHODs
C     FOR DISCRETE SEQUENTIAL ESTIMATION", ACADEMIC PRESS (1977)
C     PP.147-148
C
C *** INPUTS:
C
C     U(N*(N+1)/2)      INPUT U-D VECTOR STORED COVARIANCE FACTORS.
C                       THE COLORED NOISE UPDATE RESULT RESIDES
C                       IN U ON OUTPUT
C
C     N                 FILTER DIMENSION. WHEN THE LAST COLUMN OF U
C                       HOUSES THE FILTER ESTIMATES, THEN
C                       N=NUMBER FILTER VARIABLES + 1
C
C     IP                THE LOCATION OF THE COLORED NOISE TERM
C                       (IP.GE.1.AND.IP.LE.N)
C
C     BX(IBX)           PROCESS NOISE MULTIPLIER VECTOR
C                       (DESTROYED BY PROGRAM)
C                       IT IS ASSUMED THAT BX(IP) IS IMPLICITLY
C                       UNITY, SO THAT Q WILL CORRESPOND TO
C                       THE COLORED NOISE VARIABLE PROCESS NOISE
C                       VARIANCE.
C
C     IBX               LOCATION OF THE LAST NONZERO ENTRY IN THE
C                       PROCESS NOISE MULTIPLER VECTOR (BX). THE
C                       ASSUMPTION IS IP .GT. IBX. WHEN IBX IS ZERO
C                       THE ALGORITHM REUCES TO THE STANDARD CASE
C                       OF A DIAGONAL PROCESS NOISE COVARIANCE
C
C     EM                COLORED NOISE MAPPING TERM
C                       (UNALTERED BY PROGRAM)
C
C     Q                 PROCESS NOISE VARIANCE .GE. 0
C                       (UNALTERED BY PROGRAM)
C
C *** OUTPUTS:
C
C     U(N*(N+1)/2)      THE COLORED NOISE UPDATE RESULT RESIDES
C                       IN U ON OUTPUT
C
C     IERR            DIAGNOSTIC FLAG (FROM RNK1)
C
C                     IF  .EQ.  0       NORMAL RETURN
C
C                               J       ERROR RETURN DUE TO A
C                                       NEGATIVE COMPUTED J-TH
C                                       DIAGONAL ENTRY
C
C                              -1       UD RESULT IS SINGULAR
C
C ***  REQUIRED SUBROUTINE:
C
C               RNK1  <ESL>
C
C
C     COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JANUARY 27, 1985
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C              REVISED 10      NOV 1986         VERSION 3.00.00
C              REVISED  5      JAN 1986         VERSION 4.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE UCOLRD   IS A PART OF THE FEA. INC           *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *                   COPYRIGHT 1982 - 1987                        *
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

      DOUBLE PRECISION S, EPS, TMP
 
      DIMENSION U(*),BX(*)
C
C      DATA EPS/0.6D-37/! suitable for single precision
      DATA EPS/1.0D-307/

C
C     IN SINGLE PRECISION EPS(ILION) IS MACHINE ACCURACY
C     ADEQUATE FOR VAX TYPE MACHINES
C
C * * * * * * INITIALIZATION
C
      IPM1  = IP-1
      JJOLD = IP*IPM1/2
      JJ    = JJOLD+IP
      NM1   = N-1
      IERR  = 0
C
C * * * * * *
C
      TMP = EM*U(JJ)
      S   = TMP*EM+Q
C
C     S = EM*EM*D(IP)+Q
C
      IF (IP .LT. N) THEN
C
C        MULTIPLY ROW IP OF U BY EM
C
         IJ  = JJ
         DO 10 J = IP,NM1
            IJ = IJ+J
            U(IJ) = U(IJ)*EM
  10     CONTINUE
      ENDIF
      IF (IP .GT. 1) THEN
         IF (S .LE. EPS) THEN
            DO 20 J = 1, IP
               U(J+JJOLD) = 0.0
  20        CONTINUE
         ELSE
            C     = Q*U(JJ)/S
            U(JJ) = S
            TMP   = TMP/S
C
C THIS ARRANGEMENT FOR C AND TMP IS WELL DEFINED FOR
C ALL Q .GE. 0 AND ABS(EM) .LE. 1)
C
C  SET BX ELEMENTS IN THE CASE THAT BX(IP) = -1.0
C
            IF (IBX .LE. 0) THEN
               DO 30 J = 1, IPM1
                  L    = JJOLD+J
                  BX(J) = U(L)
                  U(L) = TMP*BX(J)
  30           CONTINUE
            ELSE
C                     SET UNSPECIFIED TERMS TO ZERO
C
               DO 35 J = IBX+1,IPM1
                  BX(J) = 0.0
  35           CONTINUE
               DO 40 J  = 1, IPM1
                  L     = JJOLD+J
                  S     = BX(J)
                  BX(J) = U(L)-EM*S
                  U(L)  = S+TMP*BX(J)
C
C    BX(NEW) = COLUMN IP OF U(INPUT)-EM*BX(OLD)
C             COLUMN IP OF U(OUT) = BX(INPUT)+EM*RHO*BX(NEW)
C             WHERE RHO = D(OLD)/D(NEW)
C
  40           CONTINUE
            ENDIF
            CALL RNK1(U,IPM1,C,BX,IERR)
         ENDIF
      ELSE
        U(1) = S
      ENDIF
      RETURN
      END
