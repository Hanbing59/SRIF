      SUBROUTINE THHC(R,N,A,IMAXA,M,RSOS,NSTRT,NSTOP)
C
C       THIS SUBROUTINE PERFORMS A TRIANGULARIZATION OF A
C       RECTANGULAR MATRIX INTO A SINGLY-SUBSCRIPTED
C       ARRAY BY  APPLICATION OF HOUSEHOLDER ORTHONORMAL
C       TRANSFORMATIONS.  A CONTROL FOR PARTIAL TRIANGULARIZATION
C       IS INCLUDED. SPECIAL EFFORT IS TAKEN TO AVOID POTENTIAL
C       UNDERFLOW/ZERO DIVIDE PROBLEMS ON MACHINES HAVING
C       A RESTRICTED EXPONENTIAL RANGE [VIZ. ON THE ORDER
C       OF 10**(+/-)38].
C
C *** INPUTS:
C
C       R(N*(N+3)/2) VECTOR STORED SQUARE ROOT INFORMATION MATRIX
C                    (LAST N LOCATIONS MAY CONTAIN A RIGHT HAND
C                     SIDE)
C
C       N            DIMENSION OF R MATRIX
C
C       A(M,NBAR)    MEASUREMENT MATRIX. WHEN A HAS LEADIN COLUMNS
C                    OF ZEROS, AND NSTRT .GT. 1, THESE COLUMNS ARE
C                    OMITTED FROM THE INPUT (I.E. COLUMN J OF THE
C                    INPUT A CORRESPONDS TO COLUMN NSTRT+J-1 OF
C                    THE ORIGINAL MEASUREMENT MATRIX.
C
C                    ALSO NBAR =  N+(NSTRT-1) COLUMNS WHEN THERE
C                    IS NOT DATA (I.E. RSOS .LT. 0) AND
C                    N+1-(NSTRT-1) COLUMNS OTHERWISE.
C
C       IMAXA        ROW DIMENSION OF A
C
C       M            NUMBER OF ROWS OF A THAT ARE TO BE COMBINED
C                    WITH R (M.LE.IMAXA)
C
C       RSOS         ACCUMULATED ROOT SUM OF SQUARES OF THE
C                    RESIDUALS SQRT(Z-A*X(EST)**2), INCLUDES A
C                    PRIORI RSOS MUST BE INPUT AS A VARIABLE, NOT
C                    AS A NUMERICAL VALUE. IF INPUT RSOS.LT.0, NO
C                    RSOS COMPUTATION OCCURS.
C
C       NSTRT        FIRST COL OF THE INPUT A MATRIX THAT HAS A
C                    NONZERO ENTRY. IF NSTRT.LE.1, IT IS SET TO 1.
C                    THIS OPTION  IS CONVENIENT WHEN PACKING A
C                    PRIORI BY BATCHES AND THE A MATRIX HAS
C                    LEADING COLUMNS OF ZEROS.
C
C       NSTOP        LAST COLUMN OF TRIANGULARIZATION.
C                    IF NSTOP .LT. NSTART, OR NSTOP .GT. N
C                    IT DEFAULTS TO N (THE NOMINAL CASE).
C
C     ON ENTRY R CONTAINS A PRIORI SQUARE ROOT INFORMATION FILTER
C        (SRIF)ARRAY, AND ON EXIT IT CONTAINS THE A POSTERIORI
C        (PACKED) ARRAY.
C     ON ENTRY A CONTAINS OBSERVATIONS WHICH ARE DESTROYED BY THE
C        INTERNAL COMPUTATIONS.
C     ON ENTRY IF RSOS IS .LT. ZERO ,PROGRAM WILL ASSUME THERE IS
C        NO RIGHT HAND SIDE DATA AND WILL NOT ALTER RSOS OR USE
C        LAST N LOCATIONS OF VECTOR R.
C
C    REQUIRED SUBROUTINES:
C
C       DSQRT   <FORLIB>
C       ABS     <FORLIB>
C       MAX     <FORLIB>
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23, 1981
C
C              REVISED 01 SEPTEMBER 1982      VERSION 1.00.00
C              REVISED 20  DECEMBER 1984      VERSION 1.01.00
C              REVISED  4      JULY 1984      VERSION 1.02.00
C              REVISED 20      MAY  1985      VERSION 2.00.00
C              REVISED  5      MAY  1986      VERSION 2.01.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE THHC   IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *               COPYRIGHT 1982, 1983, 1984                       *
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

      DIMENSION A(IMAXA,*),R(*)
      DOUBLE PRECISION SUM, BETA, DELTA, DEPS
C
C    MACHINE DEPENDENT ACCURACY TERM
C    ADEQUATE FOR VAX TYPE APPLICATIONS
C
C      EPS   =  0.6E-37 ! suitable for single precision
      DATA EPS/1.0D-307/

      DEPS  =  EPS
C
      NSTART = NSTRT
      NSTOP1 = NSTOP
C
      IF (NSTART.LE.0) NSTART=1
C                                 NO. COLUMNS OF R
      NP1=N+1
C                                 NO COLS. = N IF RSOS.LT.0
      IF(RSOS .LT. 0.0D0) NP1=N
C
      NSTM1 = NSTART - 1
C
      KK=NSTART*(NSTM1)/2
      IF ((NSTOP .LT. NSTART) .OR. (NSTOP .GT. N)) NSTOP1=N
C
C                          J-TH STEP OF HOUSEHOLDER REDUCTION
C
      DO 100 J=NSTART,NSTOP1
         KK    = KK+J
         JMNST = J-NSTM1
C
C     FIRST FIND LARGEST COLUMN ELEMENT.
C     TO BE USED FOR SCALING TO AVOID POTENTIAL NUMERICAL
C     PROBLEMS
C
         TEMP = 0.0
         DO 10 I = 1, M
            TEMP = MAX(TEMP,ABS(A(I,JMNST)))
  10     CONTINUE
C                         IF J-TH COL. OF A.EQ.0 GO TO STEP
C                                                      J+1
         IF (TEMP .LE. EPS) GO TO 100
         TEMP  =  MAX (TEMP,ABS(R(KK)))
         BETA  = 1.0/TEMP
         SUM   = (R(KK)*BETA)**2
         DO 20 I=1,M
            SUM = SUM + (A(I,JMNST)*BETA) **2
  20     CONTINUE
         SUM=DSQRT(SUM)*TEMP
         IF(R(KK) .GT. 0.0) SUM=-SUM
         DELTA = R(KK)-SUM
         R(KK) = SUM
         JP1   = J+1
         IF (JP1.GT.NP1) GO TO 105
         BETA  = 1.0/SUM
         JK    =  KK
         L     =  J
C
C     NORMALIZE HOUSEHOLDER VECTOR
C
         DELTA = DELTA*BETA
         DO 25 I=1, M
            A(I,JMNST) = A(I,JMNST)*BETA
  25     CONTINUE
C
C        ** READY TO APPLY J-TH HOUSEHOLDER TRANS.
C
         DO 40 K=JP1,NP1
            JK     = JK+L
            L      = L+1
            SUM    = DELTA*R(JK)
            KMNST  = K-NSTM1
            DO 30 I=1,M
               SUM = SUM+A(I,JMNST)*A(I,KMNST)
  30        CONTINUE
            IF(DABS(SUM) .LE. DEPS) GO TO 40
            R(JK) = R(JK) + SUM
            SUM   = SUM/DELTA
            DO 35 I=1,M
               A(I,KMNST)=A(I,KMNST)+SUM*A(I,JMNST)
  35        CONTINUE
  40     CONTINUE
 100  CONTINUE
 105  IF(RSOS .LT. 0.0) RETURN
C
C     CALCULATE RSOS
C
      SUM = RSOS*RSOS
      L   = NP1 - NSTM1
      DO 110 I=1,M
         SUM = SUM+A(I,L)**2
 110  CONTINUE
      RSOS = DSQRT(SUM)
      RETURN
      END
