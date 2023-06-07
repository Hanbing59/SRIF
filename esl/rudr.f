       SUBROUTINE RUDR(RIN,ROUT,N,ESTFLG,CNB)
C
C     FOR N.GT.0 THIS SUBROUTINE TRANSFORMS AN UPPER TRIANGULAR
C     VECTOR STORED SRIF MATRIX TO U-D FORM, AND WHEN N.LT.0 THE
C     U-D VECTOR STORED ARRAY IS TRANSFORMED TO A VECTOR STORED
C     SRIF ARRAY
C
C *** INPUTS:
C
C     RIN((N+1)*(N+2)/2)   INPUT VECTOR STORED SRIF OR U-D ARRAY
C                          SRIF INPUT WITH ZERO ROWS AND COLUMNS
C                          IS PERMITTED.
C
C     N                    ABS(N)= MATRIX DIMENSION .GE.1
C           N.GT.0         THE (INPUT) SRIF ARRAY IS (OUTPUT)
C                          IN U-D FORM
C           N.LT.0         THE (INPUT) U-D ARRAY IS (OUTPUT)
C                          IN SRIF FORM
C
C     IF ESTFLG  THEN      THERE IS A RT. SIDE INPUT TO THE SRIF AND
C                          AN ESTIMATE FOR THE U-D ARRAY. THESE
C                          RESIDE IN COLUMN N+1.
C                ELSE
C                          THERE IS NO RT. SIDE OR ESTIMATE STORED
C                          IN COLUMN N+1, AND RIN NEED HAVE ONLY
C                          N COLUMNS, I.E. RIN(N*(N+1)/2)
C
C *** OUTPUT:
C
C     ROUT((N+1)*(N+2)/2)  OUTPUT IS THE CORRESPONDING U-D OR SRIF
C                          ARRAY (RIN=ROUT IS PERMITTED)
C
C     CNB                  CONDITION NUMBER BOUND ASSOCIATED
C                          WITH THE TRIANGULAR MATRIX INVERSION.
C
C
C     REQUIRED SUBROUTINES:
C
C              RINCON <ESL>
C              IABS   <FORLIB>
C              SQRT   <FORLIB>
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              19 OCTOBER 1981
C
C              REVISED  AUGUST 2, 1982    VERSION 0.00.01
C                                         VERSION 1.00.00
C              REVISED 20      MAY 1985   VERSION 2.00.00
C              REVISED 30      JAN 1986   VERSION 2.01.00
C              REVISED 22     JULY 1986   VERSION 2.02.00
C              REVISED  2     NOV  1986   VERSION 2.03.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE RUDR   IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *             COPYRIGHT 1982, 1983, 1984, 1985, 1986             *
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

      DIMENSION RIN(*), ROUT(*)
      LOGICAL ESTFLG
C
      IF (ESTFLG) THEN
         NP1 = IABS(N) +1
      ELSE
         NP1 = IABS(N)
      ENDIF
C                                 INITIALIZE DIAGONAL INDEX
      JJ    = 1
      IDIMR = NP1*(NP1 +1)/2
      IF (ESTFLG) THEN
         RNN        = RIN(IDIMR)
         RIN(IDIMR) = -1.0
         ROUT(IDIMR)= -1.0
      END IF
      IF (N .GE. 0) THEN
         CALL RINCON(RIN,ROUT,NP1,CNB)
         ROUT(1) = ROUT(1)**2
         DO 20 J=2, N
            JJN = JJ + J
            IF (ROUT(JJN) .NE. 0.0) THEN
                S = 1.0/ROUT(JJN)
            ELSE
                S = 0.0
            END IF
            ROUT(JJN) = ROUT(JJN)**2
            JM1 = J-1
            DO 10 I=1,JM1
               ROUT(JJ+I)= ROUT(JJ+I)*S
   10       CONTINUE
            JJ = JJN
   20    CONTINUE
      ELSE
C                                 NN=NEGATIVE N
         NN = -N
         ROUT(1)= SQRT(RIN(1))
C
         DO 50 J=2,NN
            JJN = JJ + J
            ROUT(JJN)= SQRT(RIN(JJN))
            S   = ROUT(JJN)
            JM1 = J-1
            DO 40 I=1,JM1
               ROUT(JJ+I)= RIN(JJ+I)*S
   40       CONTINUE
            JJ=JJN
   50    CONTINUE
         IF (ESTFLG) THEN
            DO 60 I = 1, NN
               ROUT(IDIMR-I) = RIN(IDIMR-I)
   60       CONTINUE
         ENDIF
         CALL RINCON(ROUT,ROUT,NP1,CNB)
      END IF
      IF (ESTFLG) THEN
         RIN (IDIMR) =  RNN
         ROUT(IDIMR) = -1.0
      ENDIF
      RETURN
      END
