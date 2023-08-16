!! Purpose    : COMPUTE THE INVERSE OF THE UPPER TRIANGULAR VECTOR STORED
!!              INPUT MATRIX RIN AND STORE THE RESULT IN ROUT. (RIN=ROUT IS
!!              PERMITTED) AND TO COMPUTE A CONDITION NUMBER ESTIMATE,
!!              CNB=FROB.NORM(R)*FROB.NORM(R**-1).
!!              THE FROBENIUS NORM IS THE SQUARE ROOT OF THE SUM OF SQUARES
!!              OF THE ELEMENTS. THIS CONDITION NUMBER BOUND IS USED AS
!!              AN UPPER BOUND AND IT ACTS AS A LOWER BOUND ON THE ACTUAL
!!              CONDITION NUMBER OF THE PROBLEM. (SEE THE BOOK 'SOLVING LEAST
!!              SQUARES', BY LAWSON AND HANSON)
!!              CNB/N .LE. CONDITION NUMBER .LE. N*CNB
!!              WHEN RIN IS SINGULAR, RINCON COMPUTES AN INVERSE CORRESPONDING
!!              TO THE THE SYSTEM HAVING ZERO ROWS AND COLUMNS ASSOCIATED WITH
!!              THE SINGULAR DIAGONAL ENTRIES.
!! Parameters :
!!              RIN --- INPUT VECTOR STORED UPPER TRIANGULAR MATRIX
!!             ROUT --- OUTPUT VECTOR STORED UPPER TRIANGULAR MATRIX INVERSE
!!                      IF THE COMPILER PERMITS ROUT MAY OVERWRITE RIN
!!                N --- DIMENSION OF R MATRICES, N .GE. 1
!!              CNB --- IF CNB IS .LT. 0, THE CONDITION NUMBER BOUND
!!                      CALCULATION IS SKIPPED, AND CNB IS LEFT
!!                      UNCHANGED.
!!                      For output, IF CNB WAS NON-NEGATIVE ON INPUT THEN
!!                      CNB/N .LE. C .LE. CNB WHERE C IS THE CONDITION
!!                      NUMBER OF RIN. Otherwise, UNCHANGED.

SUBROUTINE RINCON(RIN,ROUT,N,CNB)

implicit none

real*8    RIN(*),ROUT(*)
integer*4 N
real*8    CNB
!local
integer*4 I,J,K,JJ,JJOLD,II,IK,JM1
integer*4 NTOT,NBAR
real*8    TEMP
real*8    RNM,DINV,SUM,EPS,Z,ONE
DATA EPS/1.0D-307/

Z=0.D0
ONE=1.D0
NTOT=N*(N+1)/2

! CHECK TO SEE IF THE RHS IS INCLUDED
IF(CNB.GE.0.0) THEN
    NBAR = NTOT
    IF(RIN(NTOT).EQ.-1.0) NBAR = NTOT-N
    TEMP = ABS(RIN(1))
    DO J = 2,NBAR
        TEMP = MAX(TEMP,ABS(RIN(J)))
    enddo
    ! IF RIN IS THE ZERO MATRIX
    IF(TEMP.EQ.0.0) RETURN
    RNM = 0.0
    DO J=1,NBAR
        RNM = RNM+(RIN(J)/TEMP)**2
    enddo
    CNB = TEMP*SQRT(RNM)
ENDIF
SUM = RIN(1)
ROUT(1) = Z
IF(DABS(SUM).GT.EPS) ROUT(1) = ONE/SUM
IF(N.ne.1) then
    JJ=1
    DO J=2,N
        JJOLD=JJ
        JJ=JJ+J
        SUM = RIN(JJ)
        DINV = Z
        IF(DABS(SUM).GT.EPS) DINV = ONE/SUM
        ROUT(JJ)=DINV
        II=0
        IK=1
        JM1=J-1
        DO I=1,JM1
            II =II+I
            IK =II
            SUM = Z
            DO K=I,JM1
                SUM=SUM+ROUT(IK)*RIN(JJOLD+K)
                IK=IK+K
            enddo
            ROUT(JJOLD+I)=-SUM*DINV
        enddo
    enddo
endif
IF(CNB.GE.0.0) THEN
    TEMP = ABS(ROUT(1))
    DO J = 2,  NBAR
        TEMP = MAX(TEMP,ABS(ROUT(J)))
    enddo
    RNM = 0.0
    DO J=1, NBAR
        RNM = RNM+(ROUT(J)/TEMP)**2
    enddo
    CNB = CNB*TEMP*SQRT(RNM)
ENDIF

RETURN
END