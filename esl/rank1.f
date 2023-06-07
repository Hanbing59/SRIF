      SUBROUTINE RANK1 (UIN,UOUT,N,C,V,IERR)
C
C     STABLE  U-D FACTOR RANK 1 UPDATE
C
C              (UOUT)*DOUT*(UOUT)**T=(UIN)*DIN*(UIN)**T+C*V*V**T
C
C *** INPUTS:
C
C     UIN(N*(N+1)/2)  INPUT VECTOR STORED POSITIVE SEMI-DEFINITE U-D
C                     ARRAY, WITH D ELEMENTS STORED ON THE DIAGONAL
C
C     N               MATRIX DIMENSION, N .GE. 1
C
C     C               INPUT SCALAR. SHOULD BE NON-NEGATIVE
C
C     V(N)            INPUT VECTOR FOR RANK ONE MODIFICATION.
C                     V IS DESTROYED DURING THE PROCESS
C
C *** OUTPUTS:
C
C
C     UOUT(N*(N+1)/2) OUTPUT VECTOR STORED POSITIVE (POSSIBLY)
C                     SEMI-DEFINITE U-D RESULT.
C                     COMPILER PERMITTING UOUT=UIN IS ALLOWED
C
C     IERR            DIAGNOSTIC FLAG
C
C                     IF  .EQ.  0       NORMAL RETURN
C
C                               J       ERROR RETURN DUE TO A
C                                       NEGATIVE COMPUTED J-TH
C                                       DIAGONAL ENTRY
C
C                              -1       UD RESULT IS SINGULAR
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JUNE 23, 1981
C
C              REVISED JULY 24, 1982          VERSION 0.00.01
C                                             VERSION 1.00.00
C              REVISED 20    MAY   1985       VERSION 2.00.00
C              REVISED 20    APRIL 1986       VERSION 2.00.01
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE RANK1  IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *          COPYRIGHT 1982, 1983, 1984, 1985, 1986                *
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
 
      DIMENSION     UIN(*),   UOUT(*),  V(N)
      DOUBLE PRECISION ALPHA, BETA, S, D, EPS, TST
C
C      DATA EPS/0.6D-37/ !suitable for single precision
      DATA EPS/1.0D-307/
      DATA TST/.0625D0/
C
C      IN SINGLE PRECISION EPSILON IS MACHINE ACCURACY
C      ADEQUATE FOR VAX TYPE APPLICATIONS
C
C     TST=1/16 IS USED FOR RANK1 ALGORITHM SWITCHING
C
      Z    = 0.D0
      IERR = 0
C
      JJ = N*(N+1)/2
      CJ = C
      IF (CJ .GT. Z) GO TO 4
         DO 1 J=1,JJ
   1     UOUT(J)=UIN(J)
         RETURN
C
   4  IF (N .EQ. 1) GO TO 80
C
      NP2=N+2
      DO 70 L=2,N
         J=NP2-L
         S=V(J)
         BETA=CJ*S
         D=UIN(JJ)+BETA*S
         IF (D.GT.EPS) GO TO 30
         IF (D.GE.Z) GO TO 10
            IERR = J
            RETURN
   10   JJ = JJ-J
        IERR = -1
        DO 20 K=1,J
            UOUT(JJ+K) = Z
   20   CONTINUE
        GO TO 70
   30   BETA     = BETA/D
        ALPHA    = UIN(JJ)/D
        CJ       = ALPHA*CJ
        UOUT(JJ) = D
        JJ       = JJ-J
        JM1      = J-1
        IF (ALPHA .LT. TST) GO TO 50
           DO 40 I=1,JM1
              V(I)  = V(I)-S*UIN(JJ+I)
   40         UOUT(JJ+I) = BETA*V(I) + UIN(JJ+I)
           GO TO 70
   50      DO 60 I=1,JM1
              D = V(I) - S*UIN(JJ+I)
              UOUT(JJ+I) = ALPHA*UIN(JJ+I) + BETA*V(I)
   60      V(I) = D
   70   CONTINUE
C                           CJ*V(1)**2
   80   UOUT(1) = UIN(1) + (CJ*V(1))*V(1)
      RETURN
      END
