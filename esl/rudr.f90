!! Purpose    : FOR N.GT.0, TRANSFORMS AN UPPER TRIANGULAR VECTOR STORED SRIF
!!              MATRIX TO U-D FORM. AND WHEN N.LT.0, THE U-D VECTOR STORED
!!              ARRAY IS TRANSFORMED TO A VECTOR STORED SRIF ARRAY
!! Parameters :
!!              RIN --- INPUT VECTOR STORED SRIF OR U-D ARRAY SRIF INPUT
!!                      WITH ZERO ROWS AND COLUMNS IS PERMITTED
!!             ROUT --- OUTPUT IS THE CORRESPONDING U-D OR SRIF ARRAY
!!                      (RIN=ROUT IS PERMITTED)
!!                N --- ABS(N)= MATRIX DIMENSION .GE.1
!!                      N.GT.0, THE (INPUT) SRIF ARRAY IS (OUTPUT) IN U-D FORM
!!                      N.LT.0, THE (INPUT) U-D ARRAY IS (OUTPUT) IN SRIF FORM
!!           ESTFLG --- TRUE, THERE IS A RT. SIDE INPUT TO THE SRIF AND AN
!!                            ESTIMATE FOR THE U-D ARRAY. THESE RESIDE IN
!!                            COLUMN N+1.
!!                      FALSE, THERE IS NO RT. SIDE OR ESTIMATE STORED IN
!!                             COLUMN N+1, AND RIN NEED HAVE ONLY N COLUMNS,
!!                             I.E. RIN(N*(N+1)/2)
!!              CNB --- CONDITION NUMBER BOUND ASSOCIATED WITH THE TRIANGULAR
!!                      MATRIX INVERSION.

SUBROUTINE RUDR(RIN,ROUT,N,ESTFLG,CNB)

IMPLICIT none

real*8    RIN(*)
real*8    ROUT(*)
integer*4 N
LOGICAL   ESTFLG
real*8    CNB
!local
integer*4 NP1,JJ,IDIMR,J,JJN,JM1,I,NN
real*8    RNN,S
external  RINCON

IF(ESTFLG) THEN
    NP1 = IABS(N) +1
ELSE
    NP1 = IABS(N)
ENDIF

! INITIALIZE DIAGONAL INDEX
JJ = 1
IDIMR = NP1*(NP1 +1)/2
IF(ESTFLG) THEN
    RNN = RIN(IDIMR)
    RIN(IDIMR) = -1.0
    ROUT(IDIMR)= -1.0
ENDIF
IF(N.GE.0) THEN
    CALL RINCON(RIN,ROUT,NP1,CNB)
    ROUT(1) = ROUT(1)**2
    DO J=2, N
        JJN = JJ + J
        IF(ROUT(JJN).NE.0.0) THEN
            S = 1.0/ROUT(JJN)
        ELSE
            S = 0.0
        ENDIF
        ROUT(JJN) = ROUT(JJN)**2
        JM1 = J-1
        DO I=1,JM1
            ROUT(JJ+I)= ROUT(JJ+I)*S
        enddo
        JJ = JJN
    enddo
ELSE
    ! NN=NEGATIVE N
    NN = -N
    ROUT(1)= SQRT(RIN(1))
    DO J=2,NN
        JJN = JJ + J
        ROUT(JJN)= SQRT(RIN(JJN))
        S = ROUT(JJN)
        JM1 = J-1
        DO I=1,JM1
            ROUT(JJ+I)= RIN(JJ+I)*S
        enddo
        JJ=JJN
    enddo
    IF(ESTFLG) THEN
        DO I = 1, NN
            ROUT(IDIMR-I) = RIN(IDIMR-I)
        enddo
    ENDIF
    CALL RINCON(ROUT,ROUT,NP1,CNB)
ENDIF
IF(ESTFLG) THEN
    RIN(IDIMR) =  RNN
    ROUT(IDIMR) = -1.0
ENDIF

RETURN
END