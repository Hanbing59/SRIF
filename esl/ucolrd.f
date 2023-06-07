      SUBROUTINE UCOLRD(U,N,IP,IBZ,BX,IBX,EM,Q,V,IERR)
C
C     COLORED NOISE U-D COVARIANCE FACTOR UPDATING
C     U*D*(U**T)-OUTPUT=PHI*U*D*(U**T)*(PHI**T)+Q*V*V**T
C     WHERE:
C       PHI = IDENTITY EXCEPT FOR IP-TH DIAGONAL, WHICH IS EM.
C       V**T = (0**T,BX**T,0**T,1,0**T)
C
C     THE ALGORITHM USED IS THE BIERMAN-THORNTON ONE COMPONENT
C     AT-A-TIME UPDATE. CF.BIERMAN "FACTORIZATION METHOD
C     FOR DISCRETE SEQUENTIAL ESTIMATION", ACADEMIC PRESS (1977)
C     PP.147-148
C
C *** INPUTS:
C
C     U(N*(N+1)/2)      INPUT U-D VECTOR STORED COVARIANCE FACTORS.
C                       THE COLORED NOISE UPDATE RESULT RESIDES
C                       IN U ON OUTPUT
C
C     N                 U MATRIX DIMENSION. WHEN THE LAST COLUMN OF U
C                       HOUSES THE FILTER ESTIMATES, THEN
C                       N=NUMBER FILTER VARIABLES + 1
C
C     IP                THE LOCATION OF THE COLORED NOISE TERM
C                       (IP.GE.1.AND.IP.LE.N)
C
C     IBZ               NUMBER OF LEADING ZEROES PRECEDING BX, IBZ .GE. 0.
C
C     BX(IBX)           PROCESS NOISE MULTIPLIER VECTOR, NON-TRIVIAL
C                       PORTION. SEE EXAMPLE FOLLOWING IBX.
C
C     IBX               NUMBER OF NONZERO ENTRIES IN THE
C                       PROCESS NOISE MULTIPLER VECTOR (BX). THE
C                       ASSUMPTION IS IP .GE. IBX+IBZ=NBX. WHEN IBX IS ZERO
C                       THE ALGORITHM REDUCES TO THE STANDARD CASE WITH
C                       A DIAGONAL PROCESS NOISE COVARIANCE. WHEN NBX=IP
C                       THE VALUE OF BX(IBX) IS IGNORED; IT IS ASSUMED
C                       TO BE UNITY.
C
C      EXAMPLE: OUR THE TRUE PROCESS NOISE MULTIPLIER, V, BE
C
C                        V**T = ( 0.0,0.0,0.0,2.0,0.0,0.0,1.0)
C
C               INPUTS TO UCOLRD WOULD BE: IBZ = 3, BX = (2.0), IBX=1, IP = 7
C
C               NOTE THAT THE IP-TH TERM OF THE PROCESS NOISE MULTIPLER
C               IS IMPLICTLY UNITY. IT IS PERMISSIBLE TO WRITE
C
C                 IBZ=3, BX = (2.0,0.0,0.0,1.0), IBX=4
C
C               IT IS EVEN PERMISSIBLE (ALBEIT WASTEFUL) TO HAVE BX=V.
C
C     EM                COLORED NOISE MAPPING TERM
C                       (UNALTERED BY PROGRAM)
C
C     Q                 PROCESS NOISE VARIANCE .GE. 0
C                       (UNALTERED BY PROGRAM) IT IS ASSUMED THAT
C                       THE COLORED NOISE IP-TH COMPONENT OF THE THEORETICAL
C                       (IMPLICTLY DEFINED) NOISE MULTIPLER IS UNITY.
C
C *** OUTPUTS:
C
C     U(N*(N+1)/2)      THE COLORED NOISE UPDATE RESULT RESIDES
C                       IN U ON OUTPUT
C
C     IERR            DIAGNOSTIC FLAG
C
C                     IF  .EQ.  0       NORMAL RETURN
C
C                               J       ERROR RETURN DUE TO A
C                                       NEGATIVE COMPUTED J-TH
C                                       DIAGONAL ENTRY (FROM RNK1)
C
C                              -1       UD RESULT (FROM RNK1) IS SINGULAR
C
C                              -2       IBX+IBZ .GT. IP. NO PROCESSING
C                                       TOOK PLACE.
C
C                              -3       IBZ,IBX,IP, OR N LT 0
C *** SCRATCH:
C
C    V(IP-1)         THE INPUT BX IS EMBEDED IN V.
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
C              APRIL 20 , 1987
C
C              REVISED FROM VERSION 4.00 OF UCOL
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
 
      DIMENSION U(*),BX(*),V(*)
C
C      DATA EPS/0.6D-37/ ! suitable for single precision
      DATA EPS/1.0D-307/

C
C     IN SINGLE PRECISION EPS(ILON) IS MACHINE ACCURACY
C     ADEQUATE FOR VAX TYPE MACHINES
C
C * * * * * * INITIALIZATION
C
      NBX   = IBX+IBZ
      IF (NBX .GT. IP) THEN
         IERR = -2
         RETURN
      ENDIF
      IF (IBX .LT. 0 .OR. IBZ .LT. 0 .OR. IP .LT. 0 .OR. N .LT. 0) THEN
         IERR = -3
         RETURN
      ENDIF
      IERR  = 0
      IPM1  = IP-1
      JJOLD = IP*IPM1/2
      JJ    = JJOLD+IP
      NM1   = N-1
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
C  WHEN IBX .LE. 0 USE DIAGONAL LOGIC
C
            IF (IBX .LE. 0) THEN
               DO 30 J = 1, IPM1
                  L    = JJOLD+J
                  V(J) = U(L)
                  U(L) = TMP*V(J)
  30           CONTINUE
            ELSE
C                     SET UNSPECIFIED TERMS TO ZERO
C
               DO 35 J = 1,IBZ
                  V(J) = 0.0
  35           CONTINUE
               DO 36 J = 1,IBX
                  V(J+IBZ) = BX(J)
  36           CONTINUE
               DO 37 J = NBX+1,IPM1
                  V(J) = 0.0
  37           CONTINUE
               DO 40 J  = 1, IPM1
                  L     = JJOLD+J
                  S     = V(J)
                  V(J)  = U(L)-EM*S
                  U(L)  = S+TMP*V(J)
C
C    V(NEW) = COLUMN IP OF U(INPUT)-EM*V(OLD)
C             COLUMN IP OF U(OUT) = V(INPUT)+EM*RHO*BX(NEW)
C             WHERE RHO = D(OLD)/D(NEW)
C
  40           CONTINUE
            ENDIF
            CALL RNK1(U,IPM1,C,V,IERR)
         ENDIF
      ELSE
        U(1) = S
      ENDIF
      RETURN
      END
