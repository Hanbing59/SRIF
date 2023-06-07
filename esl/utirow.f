      SUBROUTINE  UTIROW (RIN,N,ROUT,NRY,IERR)
C
C     TO COMPUTE THE INVERSE OF AN UPPER TRIANGULAR (VECTOR STORED)
C     MATRIX WHEN THE LOWER PORTION OF THE INVERSE IS GIVEN
C
C         ON INPUT:
C
C           RX   RXY            *     *                   RX   RXY
C     RIN=               ROUT=              WHERE     R=
C           *    *              0   RY**-1                0    RY
C
C         ON OUTPUT: RIN IS UNCHANGED AND ROUT=R**-1
C         COMPILER PERMITTING THE RESULT CAN OVER-WRITE
C         THE INPUT (I.E. RIN=ROUT)
C
C *** INPUTS:
C
C     RIN(N*(N+1)/2)    INPUT VECTOR STORED TRIANGULAR MATRIX
C                       THE BOTTOM NRY ROWS ARE IGNORED
C
C     N                 MATRIX DIMENSION
C
C
C *** OUTPUTS:
C
C     ROUT(N*(N+1)/2)   OUTPUT VECTOR STORED MATRIX. ON INPUT THE
C                       BOTTOM NRY ROWS CONTAIN THE LOWER PORTION
C                       OF R**-1. ON OUTPUT ROUT=R**-1
C
C     NRY               DIMENSION OF LOWER (ALREADY INVERTED)
C                       TRIANGULAR R. IF NRY=0, ORDINARY MATRIX
C                       INVERSION RESULTS.
C
C     IERR              DIAGNOSTIC FLAG TO ALERT USER THAT HIS
C                       INPUT RX MATRIX IS SINGULAR.
C                       IERR IS .EQ. TO THE NUMBER OF DIAGONAL
C                       ZEROS IN THE RX PORTION OF ROUT
C
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              MAY  26, 1982
C              VERSION 0.01.00
C              VERSION 0.02.01
C              REVISED 21 JUNE 1982       VERSION 0.03.00
C                                         VERSION 1.00.00
C              REVISED 20 MAY  1985       VERSION 2.00.00
C
C     REQUIRED SUBROUTINE:
C
C        DABS  <FORLIB>
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE UTIROW IS A PART OF THE FEA. INC             *
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

      DIMENSION RIN(*),ROUT(*)
      DOUBLE PRECISION SUM, ZERO, ONE, DINV, EPS
      DATA   ONE/1.D0/, ZERO/0.D0/
C      DATA EPS/0.6D-37/ ! suitable for single precision
      DATA EPS/1.0D-307/

C
C     EPS -- MACHINE EPSILON, SUCH THAT RECIPROCAL EXISTS.
C     EPS =  0.6D-37 IS ACCEPTABLE ON VAX TYPE MACHINES
C
C         INITIALIZATION
C
C    NR                 NO. ELEMENTS IN R
C    ISTRT              FIRST ROW TO BE INVERTED
C    IRLST              IRLST=PREVIOUS IROW INDEX
C    II                 II=DIAGONAL
C
      IERR  = 0
      NP1   = N+1
      NR=N*(NP1)/2
      ISTRT = N-NRY
      IRLST = ISTRT+1
      II    = ISTRT*IRLST/2
      IROW  = IRLST
      DO 40 IRW = 1, ISTRT
        IROW = IROW - 1
        SUM  = RIN(II)
        IF (DABS(SUM) .GT. EPS) GO TO 10
           IERR = IERR+1
           DINV  = ZERO
           GO TO 15
   10      DINV = ONE/SUM
   15   ROUT(II) = DINV
        KJS=NR+IROW
C                                  KJ (START) AND IK (START)
        IKS=II+IROW
C
        IF (IRLST.GT.N) GO TO 35
        J = NP1
        DO 30 J1=IRLST,N
          J = J-1
          KJS = KJS-J
          SUM = ZERO
          IK  = IKS
          KJ  = KJS
C
          DO 20 K=IRLST,J
            KJ=KJ+1
            SUM=SUM+RIN(IK)*ROUT(KJ)
   20       IK=IK+K
C
   30     ROUT(KJS)=-SUM*DINV
   35   IRLST = IROW
   40   II    = II - IROW
      RETURN
      END
