       SUBROUTINE TRCOL(R,N,NPSTRT,NP,EM,B,MAXB,NX,RW,ZW,V,VB,
     1                  H,IMAXH,SGSTAR,ESTFLG,BFLAG)
C
C     TO INCLUDE MARKOV COLORED NOISE EFFECTS INTO A TRIANGULAR
C     VECTOR-STORED SQUARE ROOT INFORMATION ARRAY, AND TO
C     GENERATE WEIGHTING COEFFICIENTS FOR SMOOTHING. IT IS
C     ASSUMED THAT VARIABLES X(NPSTRT)...,X(NPSTRT+NP-1) ARE
C     MARKOV COLORED NOISE AND THAT EACH COMPONENT SATISFIES
C     A MODEL EQUATION OF THE FORM:
C
C          X(SUB-J+1)=EM*X(SUB-J)+W(SUB-J)
C
C     MOREOVER, THE NOISE DRIVING THE COLORED NOISE IS ASSUMED
C     TO ALSO DIRECTLY INFLUENCE THE DYNAMIC VARIABLES, ASSUMED
C     TO BE THE FIRST NX VARIABLES, VIZ
C
C     X(I)SUM-J+1) = X(I)(SUB-J)+SUM[B(I,K)W(K);
C     K = 1 TO NP] ; I = 1 TO NX
C
C     THE ALGORITHM USED IS A GIVENS TRANSFORMATION
C     IMPLEMENTATION OF THE BIERMAN "ONE-AT-A-TIME" MARKOV UPDATE;
C     SEE  "FACTORIZATION METHODS FOR DISCRETE SEQUENTIAL
C     ESTIMATION", G.J. BIERMAN, ACADEMIC PRESS(1977)
C     FOR ALGORITHM DESCRIPTION.
C
C     USE OF THE SMOOTH COEFFICIENTS IS DESCRIBED IN THE
C     MONOGRAPH, AND THE MECHANICS OF THE SMOOTHER IMPLEMENTATION
C     ARE DESCRIBED IN THE FEA ESTIMATION LIBRARY USER REPORT,
C     FUNCTIONAL DESCRIPTION
C
C *** INPUTS:
C
C     R(NTOT)             UPPER TRIANGULAR VECTOR-STORED SRIF
C                         ARRAY. NBAR = NBAR*(NBAR+1)/2
C                         IF ESTFLG THEN NBAR = NBAR+1
C                            OTHERWISE   NBAR = N
C
C     N                   STATE VECTOR DIMENSION
C
C     NPSTRT              LOCATION OF FIRST COLORED NOISE
C                         VARIABLE.
C
C     NP                  NUMBER OF CONTIGUOUS MARKOV PROCESS
C                         NOISE VARIABLES.
C
C     EM(NP)              COLORED NOISE MAPPING MULTIPLIERS
C                         (USUALLY OF THE FORM EXP(-DT/TAU))
C
C     B(NX,NP)            COLORED NOISE INFLUENCE COEFFICIENTS
C                         THAT DIRECTLY COUPLE THE DRIVING NOISE
C                         TO THE NX-DYNAMIC STATE VARAIBLES.
C
C     MAXB                ROW DIMENSION OF MATRIX B, NX .LE. MAXB
C
C     BFLAG               LOGICAL, IF FALSE PROCESS NOISE IS ASSUMED
C                         NOT TO DIRECTLY EFFECT THE DYNAMIC
C                         VARIABLES.
C
C     VB(NX)              SCRATCH VECTOR, TO HOLD R(SUBX)*K-TH
C                         COLUMN OF B-MATRIX.
C
C
C     RW(NP)              RECIPROCAL PROCESS NOISE STANDARD
C                         DEVIATIONS; THESE SHOULD ALL BE
C                         POSITIVE
C
C     ZW(NP)              ZW = RW*(W-ESTIMATE); PROCESS NOISE
C                         ARE GENERALLY ZERO MEAN.
C
C     V(NBAR)             SCRATCH VECTOR.
C
C     ESTFLG (LOGICAL)    IF ESTFLG IS .TRUE. THEN
C                            ESTIMATES ARE INCLUDED (AND NBAR=N+1)
C                         IF ESTFLG IS .FALSE. THEN
C                                ZW IS IGNORED, V(NBAR)
C                                IS NOT USED, AND ONLY
C                                N-COLUMNS OF H ARE USED.
C
C     SGSTAR(1)           IF SGSTAR(1) .GE. 0.0
C                            SMOOTHING COEFFICIENTS ARE STORED
C                            IN SIGSTAR AND H. NOTE THAT THIS
C                            FLAG SETTING HAS NO EFFECT ON
C                            THE COMPUTATIONS.
C
C *** OUTPUTS:
C
C
C     R(NTOT)             RESULTING MARKOV PARAMETER TIME UPDATED
C                         SRIF ARRAY. THE OUTPUT OVERWRITES THE
C                         INPUT.
C
C     H(NP,NBAR)          MATRIX OF SMOOTHING COEFFICIENTS
C                         IF SGSTAR IS NEGATIVE THEN
C                            H IS NOT USED. (NP .LE. IMAXH
C
C     IMAXH               ROW DIMENSION OF H
C
C     SGSTAR(NP)          SMOOTHING COEFFICIENT MULTIPLIERS.
C                         IF SGSTAR(1) IS NEGATIVE THEN
C                            SGSTAR IS NOT FURTHER USED
C
C                         NOTE: IF SMOOTHNG COEFFICIENTS ARE
C                               STORED, THE SGSTAR COEFFICENTS
C                               ARE POSITIVE. THERE IS THEREFORE
C                               NO POSSIBILITY THAT ON A LATER
C                               SUBROUTINE CALL THAT THE SMOOTH
C                               FLAG WOULD BE TURNED OFF BECAUSE
C                               OF A SGSTAR(1) RESETING.
C
C       REQUIRED SUBROUTINES:
C
C         ABS    <FORLIB>   USED IN SUBROUTINE GIVENS
C         SIGN   <FORLIB>   USED IN SUBROUTINE GIVENS
C         SQRT   <FORLIB>   USED IN SUBROUTINE GIVENS
C         GIVENS <ESL>
C         RTMSX  <ESL>
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              23 OCTOBER 1985
C
C                                REVISED  15 JANUARY 1986
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE TRCOL  IS A PART OF THE FEA. INC             *
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

      DIMENSION R(*), EM(NP), RW(NP),ZW(NP), V(*), H(IMAXH,*),
     1          SGSTAR(NP),B(MAXB,NP),VB(NX)
C
      LOGICAL ESTFLG, BFLAG
C
      IF (ESTFLG) THEN
         NBAR = N+1
      ELSE
         NBAR = N
      ENDIF
C
C     JMP  = INDEX OF MARKOV PROCESS NOISE TERM
C
      JMP = NPSTRT
      DO 100 JCLRD=1, NP
         IF (BFLAG) CALL RTMSX(R,NX,B(1,JCLRD),VB)
         JMPP1= JMP+1
         JJP  = JMP*(JMPP1)/2
C
C        JJP = INDEX OF R MATRIX DIAGONAL CORRESPONDING TO JMP
C
         RWN   =  RW(JCLRD)
         CALL GIVENS((-RWN*EM(JCLRD)),R(JJP),SIG,S,C)
C
C                     C,S ARE THE GIVENS REFLECTON DIRECTION
C                     COSINES
C
         V(JMP)= C*RWN
         R(JJP)= S*RWN
         IJ    = JJP+JMP
         IF (JMP .LT. N) THEN
            DO 10 J=JMPP1, N
               V(J) = S*R(IJ)
               R(IJ)=-C*R(IJ)
               IJ = IJ+J
  10        CONTINUE
         ENDIF
C
         IF (ESTFLG) THEN
            V(NBAR) = C*ZW(JCLRD)+S*R(IJ)
            R(IJ)   = S*ZW(JCLRD)-C*R(IJ)
         ENDIF
C
C     IF ZW IS IDENTICALLY ZERO, ONE CAN SIMPLIFY V(NBAR)
C     AND R(IJ) INSTRUCTIONS
C
         IF (JJP .EQ. 1) GO TO 70
            IJXP  = JJP
            JJ    = JJP-JMP
            JMPM1 = JMP-1
            DO 60 L=1,JMPM1
               J     = JMP-L
               V(J)  = 0.0
               IJXP  = IJXP-1
               RWN   = SIG
               IF (BFLAG .AND. (J .LE. NX)) THEN
                  TMP = R(IJXP) + VB(J)*EM(JCLRD)
               ELSE
                  TMP = R(IJXP)
               ENDIF
               CALL GIVENS(RWN,TMP,SIG,S,C)
               IK    = JJ
               DO 40 K=J, JMPM1
                  SUM  = C*V(K) + S*R(IK)
                  R(IK)= S*V(K) - C*R(IK)
                  V(K) = SUM
                  IK   = IK+K
  40           CONTINUE
C
C     40 LOOP UPDATES THE ROW J ELEMENTS OF R THAT
C     ARE TO THE LEFT OF THE MARKOV STATE AT JMP
C
               R(IK)  = S * V(JMP)
               V(JMP) = C * V(JMP)
               IF (BFLAG .AND. (J .LE. NX)) THEN
                   V(JMP) = V(JMP) - S*VB(J)
                   R(IK)  = R(IK)  + C*VB(J)
               ENDIF
               IF (JMP .EQ. NBAR) GO TO 60
                  IK = IK+JMP
                  DO 50 K = JMPP1,NBAR
                     SUM   = C*V(K) + S*R(IK)
                     R(IK) = S*V(K) - C*R(IK)
                     V(K)  = SUM
                     IK    = IK+K
  50              CONTINUE
                  JJ = JJ-J
  60        CONTINUE
C
C     60 LOOP UPDATES ALL THE ROWS OF R ABOVE THE
C     MARKOV STATE AT JMP, AND COMPLETES THE COMPUTATION
C     OF THE SMOOTH COEFFICIENTS
C
  70     CONTINUE
         JMP = JMPP1
         IF (SGSTAR(1) .LT. 0.0) GO TO 100
            SGSTAR(JCLRD) = SIG
            DO 80 K=1, NBAR
               H(JCLRD,K) = V(K)
  80        CONTINUE
 100  CONTINUE
      RETURN
      END
