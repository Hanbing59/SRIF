C         Introduction
C
C The Estimation Subroutine Library (ESL) was originally authored by Gerald J.
C Bierman and Keith H. Bierman. The algorithms largely correspond closely
C to those found the Research Monograph "Factorization Methods for Discrete
C Sequential Estimation" originally published by Academic Press. It has
C been reprinted by Dover with ISBN-10: 0486449815 and ISBN-13: 978-0486449814
C
C Originally published as commerical software by Factorized Estimation
C Applications Inc. it has largely been unavailable to the general public
C for over 20 years. This library is now licensed under the Creative Commons 
C Attribution 3.0 Unported License. 
C To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/ 
C or send a letter to Creative Commons, 171 Second Street, Suite 300, 
C San Francisco, California, 94105, USA.
C
C As most are aware, Gerald J. Bierman (Ph.D., IEEE Fellow) died suddenly in 
C his prime in 1987. This library is re-released in his memory.
C
C Keith H. Bierman 2008
C khbkhb@gmail.com
C Cherry Hills Village, CO 
C
C------------------------------------------------------------------------------- 
C
C                      AUG 17 1987    (SL)
C                                RENAMED: UCOL => UCOLRD
C                                         UCOLOLD => UCOL
C
C                     JUNE 1987      (SL)
C                                       BRTSGN, LISPOS, SDOT1, UCOL,
C                                       GIVENS,SWPCOV,SWPUCL
C                                       CNS
C                     6 NOV  1986    (SL)
C                                    NEW VERSIONS OF:
C                                       RL2RU
C                                       RU2RL
C                                       XY2U
C                                       PHIU
C                                       TDHHT
C                                       TRCOL
C                                       RUDR
C
C                     9 OCT  1986    (SL)
C
C                                    NEW SUBROUTINES:
C                                       RUDR (NEW VERSION)
C                                       CNSUD
C                                       CNSXUD
C                                       PINVBM
C-------------------------------------------------------------------------
C                     5 MAY  1986
C                                     --- IMPLICIT
C                                     --- MODEST RU2RL MOD
C                                     --- MODEST RANK1,RNK1 MODIFICATION
C  IMPLEMENTATION NOTES:
C
C  THIS COPY OF THE LIBRARY IS INTENDED FOR USE WITH FORTRAN'77
C  (ANSI X3.9-1978, ISO 1539-1980(E) "FULL" COMPILERS.
C
C  IF FOR SOME REASON IT IS NECESSARY TO TRANSPORT IT TO ANSI '66
C  COMPILERS SPECIAL EFFORT WILL HAVE TO BE MADE. CONTACT FEA FOR
C  DETAILS.
C
C  IF A FORTRAN '77 SUBSET COMPILER IS USED (LIKE MSFORTRAN, OR
C  THE P-SYSTEM FORTRAN) ALL CHARACTER*(*) WILL HAVE TO BE
C  MADE FIXED LENGTH (E.G. CHARACTER*16)
C
C  THE FOLOWING ROUTINES HAVE A MACHINE DEPENDENT EPSILION TERM
C  (USUALLY CALLED EPS):
C
C  CNS, COV2RI, COV2UD, COVCOR, INF2R, RNK1, RANK1, PINVB, PINVEK,
C  RINCON, RINZ, SINTST, THHC, UDOPT, UDMVEC, UMEAS, UTINV,
C  UTIROW, WGSG
C
C  IN GENERAL IT IS THE SMALLEST NUMBER SUCH THAT 1/EPS IS MACHINE
C  REPRESENTABLE (IN THE PRECISION OF THE LIBRARY. THE ESL IS
C  DISTRIBUTED AS A SYSTEM OF SINGLE PRECISION MACHINES, SO EVEN
C  THOUGH EPS IS DOUBLE PRECISION 1/EPS MUST FIT IN A SINGLE
C  PRECISION VARIABLE). THE DISTRIBUTED VALUE 0.6E-37 IS SUITABLE FOR
C  MOST MACHINES WITH 32-BIT SINGLE PRECISION.
C  ---- 
C  Since this version of the code is "pure" double precision, 
C  a value of 1.0d-307 has been chosen which should be suitable
C  for most modern (2008+) systems. 
C
C  THE ENTIRE LIBRARY HAS BEEN COMPILED AND TESTED ON A WIDE VARIETY
C  OF MACHINES.
C
C SINCE PUBLICATION OF THE USER GUIDE THE CALLING SEQUENCE OF THE
C OUTPUT ROUTINES RMPRNT AND TMPRNT HAVE CHANGED. RMPRNT WAS ENHANCED
C TO INCLUDE AN ABBREVIATED PRINT OPTION, TMPRNT WAS CHANGED IN ORDER
C TO BE CONSISTENT WITH THE NEW RMPRNT. THE NEW RMPRNT ARGUMENT
C LIST IS AS FOLLOWS:
C
C
C      SUBROUTINE RMPRNT(A,IA,NR,NC,NAMR,NAMC,TEXT,CODE,LBLROW,
C     1                  LBLCOL,PRTNAM,RNUM,IOUNIT)
C
C
C  ALL OF THE ARGUMENTS ARE AS DEFINED IN THE GUIDE EXCEPT THAT
C  NCHAR HAS BEEN REPLACED BY CODE:
C
C     CODE        (CHARACTER*(*)) IS A LIST OF PARAMETERS,
C                 WHICH ARE INTERPRETED IN A POSITION INDEPENDENT
C                 FASHION. ANY CHARACTERS NOT DEFINED IN THIS LIST
C                 WILL HAVE NO EFFECT; BUT OTHER FEATURES MAY BE
C                 DEFINED AT A LATER DATE, SO EXTRANEOUS CHARACTERS
C                 SHOULD BE AVOIDED.
C
C      "P"        A PAGE EJECT WILL PRECEDE PRINTING
C
C      "C"        CONDENSED PRINT
C
C                   (I.E. LESS PRECISION, MORE COLUMNS PER ROW)
C                   THE DEFAULT (' ') CAUSES PRINTING TO BEGIN
C                   AT THE NEXT LINE AND DISPLAYS 6 COLUMNS WITH
C                   8 DIGITS OF PRECISION. CONDENSED MODE CONSISTS
C                   OF 10 COLUMNS OF 4 DIGIT
C                   PRECISION. IF THE NAMC CHARACTER NAMES ARE
C                   LONGER THAN 8 CHARACTERS LONG OTHER CODE
C                   PARAMETERS CONTROLS HOW EACH NAME IS TO
C                   BE PRINTED (SINCE ONLY 8 CHARACTERS CAN APPEAR
C                   PER VARIABLE PER LINE).
C
C    "BK='C'"     WHERE 'C' IS ANY CHARACTER, INCLUDING SPACE.
C
C                   CONTROL FOR HANDLING NAMC NAMES WITH LENGTH
C                   LONGER THAN 8. 'C' IS USED TO REQUEST
C                   TO PUT THE REMAINING CHARACTERS ON THE
C                   NEXT LINE. THE DEFAULT IS TO TREAT A
C                   BLANK AS THE REQUEST CHARACTER. NOTE:
C                   THE REQUIRED SYNTAX IS BK=THING, BLANKS
C                   ARE SIGNIFICANT. IT SHOULD ALSO BE NOTED
C                   THAT THIS OPTION IS IRRELEVANT UNLESS THE
C                   CONDENSED PRINT OPTION IS IN EFFECT.
C
C  EXAMPLES:   N = 3, CHARACTER*16 NAMES ='X-POSITION',
C                     'X VELOCITY','X_ACCELERATION'
C
C    CODE =    'PC' WILL PRODUCE A NEW PAGE AND  CONDENSED PRINT
C    CODE =    'C'  WILL PRODUCE CONDENSED PRINT ONLY
C
C     THE COLUMN NAMES WILL BE:
C
C     X-POSITI   X            X_ACCELE
C     ON         VELOCITY     RATION
C
C  SINCE "BK=" WAS NOT SPECIFIED, SPACE TRIGGERED THE BREAK
C
C  IF CODE =   'C BK=_' THE RESULTING COLUMN NAMES WOULD BE
C
C     X-POSITI  X VELOCI      X
C     ON        TY            ACCELERA
C
C
C  SINCE THE "BREAK" CHARACTER IS NOT PRINTED (THE FIRST TIME IT
C  OCCURS) IT IS POSSIBLE (AND MEANINGFUL) TO USE CHARACTER*17
C  NAMES, IF ONE OF THE CHARACTERS IS THE "BREAK" CHARACTER.
C
C--------------------------------------------------
C
C    TMPRNT HAS BEEN CHANGED TO
C
C         SUBROUTINE TMPRNT (R,N,NAMES,TEXT,CODE,LABEL,
C     1                   PRTNAM,RNUM,IOUNIT)
C
C   THOSE ESL ROUTINES (U2SIG,SGDCMP) WHICH CALL THESE ROUTINES
C      HAVE HAD THEIR ARGUMENT LISTS CHANGED ACCORDINGLY.
C
C   FEA IS NO LONGER SUPPORTING ANSI '66 FORTRAN. THOSE REMARKS
C   IN THE USER GUIDE REFERING TO ANSI '66 DEFAULTS, ETC. SHOULD
C   BE IGNORED. (SUCH REMARKS APPEAR WITH RESPECT TO NAME PARAMETERS
C   E.G. NAMIN, NAMOUT, NAMR, ETC.)
C
C        SUBROUTINE PINVB HAS BEEN ENHANCED. SOMETIMES
C        ONE WISHES TO DELETE A VARIABLE (OR SET OF
C        VARIABLES) FROM THE PROBLEM, WITHOUT CODE
C        MODIFICATION. THIS IS COMMONLY DONE BY SIMPLY
C        ZEROING OUT THE DIAGONAL, AND ITS ENTIRE ROW
C        AND COLUMN. PINVB NOW ALLOWS THIS, PINVEK MAY
C        BE ENHANCED SIMILARILY.
C
C JAN 86
C
C        MAPU HAS BEEN ENHANCED, IT IS NOW POSSIBLE TO COMPUTE ONLY
C        STANDARD DEVIATIONS; D(1) (PREVIOUSLY ONLY AN OUTPUT)
C         CONTROLS THIS BEHAVIOR. IF D(1) .GE. 0 THEN MAPU OPERATES
C         NORMALLY.
C
C        RUDR HAS HAD "IS" CHANGED TO ESTFLG, "IS" WAS AN INTEGER WHICH
C        SIGNALED WHETHER THERE IS A RIGHT HAND SIDE, ESTFLG
C        IS A LOGICAL VARIABLE, SO IN THE ESL ALL RHS'S ARE
C        INPUT WITH A LOGICAL VARIABLE.
C
C ---
      SUBROUTINE BRTSGN (U,N,IP,EM,BX,IBZ,NBX,Q,XEST,SGAIN,SGSTAR,
     1                   Y,C,S,IEST,IERR)
      IMPLICIT REAL*8 (A-H,O-Z)

C
C     TO COMPUTE THE BIERMAN MODIFIED R-T-S SMOOTH GAIN FACTORS FOR
C     COLORED NOISE MODEL X   = PHI  X   + BW
C                          T+      M  T-     T
C
C                               T
C     WHERE PHI  = I + (M-1)E  E
C              M             IP IP
C
C       (E  (J) = 0 J .NE. IP ; = 1 FOR J = IP)
C         IP
C
C
      INTEGER N
C                   FILTER DIMENSION
C
      DIMENSION U(N*(N+1)/2)
C
C                   VECTOR STORED U-D FILTER COVARIANCE FACTORS PRIOR
C                   TO THE COLORED NOISE COMPONENT UPDATE.
C
      INTEGER IP
C                   LOCATION OF COLORED NOISE PARAMETER; IP .GE. 1,
C                   AND IP .GT. (IBZ + NBX)
C
C     EM            COLORED NOISE (EXPONENTIAL) MULTIPLIER
C
      INTEGER NBX
C                                          T     T   T  T    T
C                   NUMBER OF ROWS IN BX (B  = (0, BX ,0 ,1,0 ))
C
      DIMENSION BX(*)
C
C                   (NBX) PROCESS NOISE SENSITIVITY MULTIPLIER
C                   TERMS WHEN 0 NOT USED.
C
      INTEGER IBZ
C                   NUMBER OF LEADING ZEROS IN B VECTOR
C
C     Q             PROCESS NOISE COMPONENT A PRIORI VARIANCE
C
      DIMENSION XEST(N)
C
C                   FILTER STATE VECTOR PRIOR TO THE COLORED NOISE
C                   COMPONENT UPDATE; USED ONLY WHEN IEST = .TRUE.
C
      LOGICAL IEST
C                   IF TRUE THE FILTER ESTIMATE IS USED
C                   TO COMPUTE SGAIN(N+1) THAT IS USED AS
C                   PART OF THE SMOOTH ESTIMATE TIME UPDATE.
C OUTPUT:
C
      DIMENSION SGAIN(N+1)
C
C                   AUGMENTED VECTOR OF SMOOTHER GAINS, WHEN
C                   IEST = .FALSE. THE (N+1)ST COMPONENT IS LEFT
C                   UNALTERED.
C
C      SGSTAR       SCALAR MULTIPLIER USED IN THE SMOOTH U-D FACTOR
C                   COVARIANCE TIME UPDATE
C
       INTEGER IERR
C                   IF .EQ. -1 THE COMPUTED SMOOTH GAIN IS SINGULAR
C                      .LE. -2 THE SMOOTH GAIN COMPUTATION FAILED.
C SCRATCH:
C
      DIMENSION Y(N),C(N),S(N)
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              17 MARCH 1987
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE BRTSGN IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *                      COPYRIGHT 1987                            *
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
C         LOCAL VARIABLES
C
      INTEGER I,J,K,IJ,JJ,JJS,IK,IPIP,IPM1,IPP1,IFV,NP1,JCOL
      REAL*8 SIG,SN,CN
      DOUBLE PRECISION YNP1,SIGSTR
C
      DOUBLE PRECISION TMP,TMPA,EPSCN,EPS
c      DATA EPS/0.6D-37/,EPSCN/1.0D-6/  ! single precision typical values
      DATA EPS /1.0D-307/,EPSCN/1.0D-6/
C
      IERR = 0
      IF ((EM .GT. 1.0 - EPSCN) .AND. Q .EQ. 0.0) THEN
         DO 5 J = 1, N+1
            SGAIN(J) = 0.0D0
 5       CONTINUE
         SGSTAR=0.0D0
         RETURN
C                      COMMON EXIT
      ENDIF
      IPM1 = IP - 1
      IPP1 = IP + 1
      IPIP = IP*IPP1/2
      NP1  = N + 1
C**************************************************************
C***  SOLVE UY = BM; BM(IP) = 1, E(J) = M*B(J) J .LT. IP    ***
C**************************************************************
C
      DO 10 I = 1, IPM1
         Y(I) = 0.0
10    CONTINUE
      DO 20 I = 1, NBX
         Y(IBZ+I) = BX(I) * EM
20    CONTINUE
      Y(IP)  = 1.0
      IK     = IPIP
      DO 40 J = IP,2,-1
         TMP = Y(J)
         DO 30 I = J-1,1,-1
            IK   = IK-1
            Y(I) = Y(I) - U(IK)*TMP
30       CONTINUE
C                              DIAGONAL OF U MATRIX
          IK = IK-1
40     CONTINUE
C                              IF NON-DIAGONAL Q
       IF (NBX .GT. 0) THEN
C                              FORM SGAIN**T = -((E**T)U*D;0)
           DO 50 I = 1, IPM1
              SGAIN(I) = 0.0
50         CONTINUE
           JJ = IPIP
           SGAIN(IP) = -U(JJ)
           DO 60 I = IPP1, N
              IJ       = JJ+IP
              JJ       = JJ+I
              SGAIN(I) = -U(JJ)*U(IJ)
              Y(I)     = 0.0
 60        CONTINUE
           SIGSTR = 0.0
           YNP1   = -EM
C********************************************************************
C**            SOLVE SYSTEM                                        **
C**----------------------------------------------------------------**
C** D*(U**T)*PHI*SGAIN + Y*SGSTAR = -D*(U**T)*E                    **
C** Q*(B**T) * SGAIN - EM*SGSTAR  = 0                              **
C**                                                                **
C** BY GIVENS TRANSFORMING IT TO TRIANGULAR FORM AND FORWARD       **
C** SOLVING AS WE GO, BUT SAVING THE GIVENS DIRECTION SINES AND    **
C** COSINES FOR LATER RECONSTRUCTION OF THE GAIN FROM THE GIVENS   **
C** ROTATED SOLUTION.                                              **
C********************************************************************
C
           JJS = 0
           DO 100 IFV = 1, N
              JJS    = JJS+IFV
              JJ     = JJS
              TMP    = U(JJ)
              IF (IFV .EQ. IP) TMP = TMP*EM
              CALL GIVENS(TMP,Y(IFV),SIG,SN,CN)
              C(IFV) = CN
              S(IFV) = SN
C                              TMP IS THE SOLVED VALUE FOR SGAIN(IFV)
C
              IF (SIG .GT. EPS) THEN
                 TMP = SGAIN(IFV)/SIG
C**********************************************************************
C**   APPLY GIVENS TO Y AND COLUMN B IFV OF THE COEFFICIENT MATRIX,  **
C**   ADJUST RHS CORRESPONDING TO THE FORWARD SOLVE.                 **
C**********************************************************************
C
                 DO 80 I = IFV+1, N
                    IJ       = JJ+IFV
                    JJ       = JJ+I
                    TMPA     = U(JJ) * U(IJ)
                    IF (IFV .EQ. IP) TMPA = TMPA*EM
                    SGAIN(I) = SGAIN(I) - (CN*TMPA + SN*Y(I)) * TMP
                    Y(I)     = SN*TMPA  - CN*Y(I)
 80              CONTINUE
C
C        B HAS THE FORM 0000000 BBBBBBBBBBBB 0 IP 0000000
C                       |_IBZ_| |_IBZ+NBX_|
C        CASE#             3         2       3  1  3
C
                 IF (IFV .EQ. IP) THEN
C
C                    B IS "VIRTUALLY" ONE  CASE #1
C
                    SIGSTR = SIGSTR - (CN*Q + SN*YNP1) * TMP
                    YNP1   = SN*Q - CN*YNP1
                 ELSEIF ((IFV .GT. IBZ) .AND.
     1                   (IFV .LE. (IBZ + NBX))) THEN
C
C                   B IS "REAL"    CASE #2
C
                    TMPA   = Q * BX(IFV-IBZ)
                    SIGSTR = SIGSTR  - (CN*TMPA + SN*YNP1) * TMP
                    YNP1   = SN*TMPA - CN*YNP1
                 ELSE
C
C                   B IS "VIRTUALLY" ZERO  CASE#3
C
                    SIGSTR = SIGSTR - SN*YNP1*TMP
                    YNP1   = -CN*YNP1
                 ENDIF
                 SGAIN(IFV) = TMP
             ENDIF
100        CONTINUE
C
C**********************************************************************
C*  NOW, CONSTRUCT SOLUTION BY APPLYING GIVENS TRANSFORMATION TO THE  *
C*  SGAIN VECTOR JUST COMPUTED                                        *
C**********************************************************************
C
           IF (ABS(YNP1) .LT. EPS) THEN
              IF (ABS(SIGSTR) .LT. EPS) THEN
                 SIGSTR = 0.0
                 IERR   = -1
              ELSE
                 SIGSTR = 0.0
              ENDIF
           ELSE
              SIGSTR = SIGSTR/YNP1
           ENDIF
           DO 110 J = N, 1, -1
                TMP = SGAIN(J)
                SGAIN(J) = C(J) * TMP + S(J) * SIGSTR
                SIGSTR   = S(J) * TMP - C(J) * SIGSTR
110        CONTINUE
           IF (SIGSTR .GT. 0)  IERR = -3
           IF (IEST) THEN
              TMP = ((EM-1.0)*SGAIN(IP) + 1.0) * XEST(IP)
              DO 120 J = 1, N
                 TMP = TMP + SGAIN(J) *XEST(J)
120           CONTINUE
              SGAIN(NP1) = TMP
           ENDIF
C--------------------------------------
           SGSTAR = -SIGSTR
           RETURN
C--------------------------------------
       ENDIF
C-------------- END OF NEW CODE (NON -DIAGONAL Q CASE)
C               BACK TO DIAGONAL Q CASE
C
C                  Y(I) = 0 FOR I > IP, BUT THESE VALUES ARE NOT NEEDED
C
C  II  (IMPLICITLY) FORM A MATRIX, A = (D(U**T)PHI+(Q/EM)(Y*E**T)PHI,
C
C  III TRIANGULARIZE A, SAVING Y AS WE GO AND THE DIRECTION
C      COSINE TERMS OF THE GIVENS ROTATIONS
C
C   FORM COLUMN IP OF A, SCALED BY EM
C
      DO 300 I = 1, IPM1
         Y(I) = Q*Y(I)
300   CONTINUE
      TMP   = EM*EM
      JJ    = IPIP
      Y(IP) = U(JJ)*TMP+Q
      DO 400 I = IPP1, N
         IJ   = JJ+IP
         JJ   = JJ+I
         Y(I) = U(JJ)*U(IJ)*TMP
400   CONTINUE
C
C   IMPLICITLY APPLY GIVENS TO COLUMNS JCOL AND IP, SAVING GIVENS
C   C AND S AND ONLY EXPLICITLY MODIFYING Y
C
      JJS = 0
C                             ! JJS(TART)
      DO 500 JCOL = 1, IPM1
         SGAIN(JCOL) = 0.0
         JJS = JJS+JCOL
         JJ  = JJS
         CALL GIVENS(U(JJ),Y(JCOL),SIG,SN,CN)
         DO 450 I = JCOL + 1 , N
            IJ   = JJ+JCOL
            JJ   = JJ+I
            Y(I) = SN*U(JJ)*U(IJ) - CN*Y(I)
450      CONTINUE
         C(JCOL) = CN
         S(JCOL) = SN
500   CONTINUE
C                     T
C   IV  SOLVE AV = -DU SUB(IP)
C
      IF (ABS(Y(IP)) .LE. EPS) THEN
         DO 510 I = 1, NP1
            SGAIN(I) = 0.0
510      CONTINUE
         SGSTAR = 0.0
         RETURN
      ENDIF
      SGAIN(IP) = -U(IPIP)/Y(IP)
      JJ        = IPIP
      DO 600 I = IPP1, N
         IJ  = JJ+IP
         IK  = IJ
         JJ  = JJ+I
         TMP = -U(IJ)
         DO 550 K = IPP1, I-1
            IK   = IK+1
            TMP  = TMP - U(IK)*SGAIN(K)
550      CONTINUE
         IF (U(JJ) .GT. EPS) THEN
            SGAIN(I) = TMP - Y(I) * SGAIN(IP)/U(JJ)
         ELSE
            SGAIN(I) = 0.0
            IF ( ABS(Y(I)) .GE. EPS) THEN
                IERR = -1
            ENDIF
         ENDIF
600   CONTINUE
C
C    SOLVE FOR SGAIN BY APPLYING GIVENS TRANSFORMATION TO SGAIN
C    (OVERWRITING RESULT) AND THEN SCALING SGAIN(IP) BY EM
C
      DO 700 J = IPM1, 1, -1
         TMP       = SGAIN(J)
         SGAIN(J)  = C(J)*TMP + S(J)*SGAIN(IP)
         SGAIN(IP) = S(J)*TMP - C(J)*SGAIN(IP)
700   CONTINUE
      SIGSTR       = -Q*SGAIN(IP)
      TMP          = EM*SGAIN(IP)
      SGAIN(IP)    = TMP+1.0D0
C                      NOTE SGAIN(IP) IS REALLY (SGAIN(IP)-1)/EM
      IF (IEST) THEN
         TMPA       = (EM-1.0D0)*TMP*XEST(IP)
         DO 950 I = 1, N
            TMPA = TMPA + SGAIN(I)*XEST(I)
950      CONTINUE
         SGAIN(NP1) = TMPA
      ENDIF
      SGSTAR        = SIGSTR
      RETURN
      END
