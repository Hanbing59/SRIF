      SUBROUTINE CNS (U,N,IP,IBZ,BX,NXMAX,NBX,UBX,NP,EM,SGSTAR,
     1                SGAIN,MAXSG,XSTAR,F,G,H,V,UPDATU,UPDATX,
     2                IERR,ISING)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C    COLORED NOISE UD SMOOTHER
C
C      U MATRIX UPDATE (WHEN UPDATU IS TRUE)
C
C        R-T-S SMOOTHER GAIN = I-E * (SGAIN)**T
C
C      ESTIMATE UPDATE (WHEN UPDATX IS TRUE)
C
C   INPUT:
C
C    U(N(N+1)/2)  VECTOR STORED U-D SMOOTH COVARIANCE
C                 FOR ESTIMATE AT STAGE I+1
C
C    N            FILTER MODEL DIMENSION
C
C    IP(NP)       LOCATIONS OF COLORED NOISE STATES IN THE
C                 FILTER MODEL.
C
C    IBZ(NP)      IBZ(J) REPRESENTS THE NUMBER OF LEADING ZEROES IN
C                 COLUMN J OF BX.
C
C    BX(NXMAX,NP) PROCESS NOISE SENSITIVITY MULTIPLIER MATRIX,
C                 NOT USED WHEN NBX() IS 0.
C
C    NXMAX        MAX ROW DIMENSION OF BX MATRIX.
C
C    NBX(NP)      NUMBER OF ROWS PER COLUMN OF BX
C
      CHARACTER*20 UBX
C
C                 WHEN UBX .EQ. 'SAME' IT IS ASSUMED THAT ALL
C                 THE COLUMNS HAVE THE SAME NUMBER OF ROWS (NBX(1))
C                 AND THE SAME NUMBER OF LEADING ZEROS, VIZ
C                 IBZ(2..NP) = IBZ(1).
C
C                 WHEN UBX .EQ. 'NONE' NBX IS TREATED AS
C                 IDENTICALLY ZERO, I.E. DIAGONAL PROCESS NOISE.
C
C                 WHEN UBX .EQ. 'ALL' ALL NBX AND IBZ ELEMENTS
C                 ARE USED.
C
C                 'NONE' IS THE DEFAULT. WHEN 'NONE', OR 'SAME'
C                 ONE NEED NOT SPECIFY INTEGER VECTORS FOR
C                 NBX AND IBZ (E.G. CALL CNS(...,2,0,'SAME',...)
C                 MEANS THAT EACH COLUMN OF BX HAS 2 ROWS, AND
C                 THERE ARE NO (IMPLICIT) LEADING ZEROS.
C
C    NP           NUMBER OF PROCESS NOISE STATES
C
C    EM(NP)         EXPONENTIAL MULTIPLIERS
C
C    SGSTAR(NP)     SMOOTH COEFFICIENT MULTIPLIERS
C
C    SGAIN(N+1,NP)  R-T-S MODIFIED GAIN
C                    *
C                   G  = I-E(SGAIN)**T
C
C                   WHERE E IS A ZERO VECTOR, EXCEPT FOR
C                   UNITY IN ONE LOCATION, IP(J). ELEMENT N+1
C                   CONTAINS: SGAIN**T * X(FILTER)
C
C    MAXSG        MAX ROW DIMENSION OF SGAIN
C
C    XSTAR(N)     SMOOTH ESTIMATES, ON ENTRY CORRESPONDING TO TIME "K+1"
C
C    UPDATU       (LOGICAL) WHEN TRUE COMPUTE UPDATED UD MATRIX
C
C    UPDATX       (LOGICAL) WHEN TRUE COMPUTE UPDATED ESTIMATES
C
C
C   OUTPUT:
C
C    U(N(N+1)/2)  UPDATED U-D FACTOR SMOOTH COVARIANCE FACTORS
C                 CORRESPONDING TO THE STATE AT TIME "K"
C
C    XSTAR(N)     SMOOTH ESTIMATES AT TIME "K", WHEN UPDATX IS TRUE.
C
C    IERR         ERROR RETURN FROM RNK1 (CALLED BY UDGEN)
C
C                     IF  .EQ.  0       NORMAL RETURN
C
C                               J       ERROR RETURN DUE TO A
C                                       NEGATIVE COMPUTED J-TH
C                                       DIAGONAL ENTRY
C
C                              -1       UD RESULT IS SINGULAR
C
C   ISING         TOTAL NUMBER OF SINGULARITIES ENCOUNTERED
C
C   SCRATCH:
C
C      F(N), G(N), H(N), V(N)
C
C    SUBROUTINES NEEDED:
C
C        UDGEN  <ESL>
C        RNK1   <ESL>
C
C    COGNIZANT PERSONS:
C
C           DR.GERALD J. BIERMAN/KEITH H. BIERMAN
C           FACTORIZED ESTIMATION APPLICATIONS
C           FEBRUARY 14, 1987
C                                VERSION 1.00.00
C
C ************************************************************
C *                                                          *
C *        SUBROUTINE CNS   IS A PART OF THE FEA. INC        *
C *              ESTIMATION SUBROUTINE LIBRARY               *
C *                                                          *
C *                          (ESL)                           *
C *                                                          *
C *                      COPYRIGHT 1987                      *
C *                  FACTORIZED ESTIMATION INC.              *
C *   This library is licensed under the Creative Commons    *
C *   Attribution 3.0 Unported License.                      *
C *   To view a copy of this license,                        *
C *   visit http://creativecommons.org/licenses/by/3.0/      *
C *   or send a letter to Creative Commons, 171 Second Street*
C *   Suite 300, San Francisco, California, 94105, USA.      *
C *                                                          *
C ************************************************************
C
      INTEGER IP(*),NP,L,NBX(*),IBZ(*)
      DIMENSION U(N*(N+1)/2),F(N),G(N),H(N),
     1          SGAIN(MAXSG,*),XSTAR(N),SGSTAR(*),EM(*),BX(NXMAX,*)
      LOGICAL UPDATU,UPDATX
      DOUBLE PRECISION  SDOT1,SUM,TMP,TMPA
C
C                       USE ISW TO DETERMINE IF NBX, IBZ ARE IN EFFECT
      LOGICAL ISW
C
C                              ASSUME NORMAL RETURN
c     CALL RMPRNT(BX,NXMAX,10,NP,N,N,' BX IN CNS','C',.FALSE.,.FALSE.,
c    1            ' ',N,6)
      IERR = 0
      NP1  = N+1
      IF (INDEX(UBX,'ALL').NE. 0 .OR. INDEX(UBX,'all') .NE. 0) THEN
          ISW = .TRUE.
      ELSEIF(INDEX(UBX,'SAME').NE. 0 .OR. INDEX(UBX,'same') .NE.0)THEN
          IBX =  NBX(1)
          IZ  =  IBZ(1)
          ISW = .FALSE.
      ELSE
C                        NONE IS THE DEFAULT
          IBX =  0
          IZ  =  0
          ISW =  .FALSE.
      ENDIF
C                              IF NO PROCESS NOISE RETURN
      IF (NP .EQ. 0) RETURN
C                              ESTIMATE UPDATE
      IF (UPDATX) THEN
        DO 200 J = NP, 1, -1
           L          = IP(J)
           IF (ISW) THEN
              IBX = NBX(J)
              IZ  = IBZ(J)
           ENDIF
C
C           WHEN COLUMN J OF BX CONTAINS ONLY ZEROS, USE
C           SIMPLER COMPUTATION
C
              SDOT1 = 0.0d0
              DO 146 K = 1, N
                 SDOT1 = SDOT1 + SGAIN(K,J)*XSTAR(K)
146           CONTINUE
           IF (IBX .EQ. 0) THEN
              XSTAR(L) = XSTAR(L) - SDOT1 + SGAIN(NP1,J)
           ELSE
              TMP        = SDOT1         ! GAMMA
              DO 150 I = 1, IBX
                 XSTAR(IZ+I) = XSTAR(IZ+I)-BX(I,J) * (XSTAR(L)+
     1                         (TMP-SGAIN(NP1,J)) * EM(J))
150           CONTINUE
              XSTAR(L) = SGAIN(NP1,J)-TMP
          ENDIF
200     CONTINUE
      ENDIF
C                               U UPDATE
      IF (UPDATU) THEN
          DO 300 J = NP,1,-1
             IF (ISW) THEN
                IBX = NBX(J)
                IZ  = IBZ(J)
             ENDIF
             L = IP(J)
C
C           WHEN COLUMN J OF BX CONTAINS ONLY ZEROS, USE
C           SIMPLER COMPUTATION
C
             IF (IBX .EQ. 0) THEN
                DO 205 I = 1, N
                   G(I) = 0.0d0
  205           CONTINUE
                G(L) = 1.0d0
                CALL UDGEN(U,N,SGAIN(1,J),SGSTAR(J),G,L,F,
     1                     H,IERR)
                IF (IERR .EQ. -1) THEN
                   ISING = ISING + 1
                ELSEIF(IERR .NE. 0) THEN
                   WRITE(6,*)' ERROR RETURN FROM UDGEN = ',IERR,
     1                    ' STATE = ',J
                ENDIF
             ELSE
C                    NON DIAGONAL PROCESS NOISE
C
C       COMPUTE P:=(I-G1*H1)P(I-G1H1)**T
C
                TMPA      = 1.0-EM(J)
                DO 210 I = 1, IZ
                   G(I) = 0.0d0
210             CONTINUE
                DO 220 I = 1, IBX
                   G(I+IZ) = TMPA*BX(I,J)
220             CONTINUE
                DO 230 I = 1, N
                   H(I) = 0.0d0
230             CONTINUE
                H(L) = 1.0d0
                CALL UDGEN(U,N,H,0.0d0,G,IZ+IBX,F,V,IERR)
C
C              P:= (I-G2*H2)P(I-G2H2)**T + SIGSTAR*G2*(G2)**T
C
                DO 250 I = 1, L-1
                   G(I) = 0.0d0
250             CONTINUE
                DO 260 I = 1, IBX
                   G(IZ+I) = EM(J)*BX(I,J)
260             CONTINUE
                G(L) = 1.0d0
                TMP = SGAIN(L,J)
                SUM = 0.0d0
                DO 270 I = 1, IBX
                   SUM = SUM + SGAIN(IZ+I,J)*BX(I,J)
270             CONTINUE
                SGAIN(L,J) = TMP + SUM*TMPA + 1.0d0
                CALL UDGEN(U,N,SGAIN(1,J),SGSTAR(J),G,L,F,V,IERR)
                SGAIN(L,J) = TMP
                IF (IERR .EQ. -1) THEN
                   ISING = ISING+1
                ELSEIF (IERR .NE. 0) THEN
                   WRITE(6,*)' ERROR RETURN FROM UDGEN2 = ',IERR,
     1                    ' STATE = ',J
                ENDIF
            ENDIF
300       CONTINUE
      ENDIF
      RETURN
      END
