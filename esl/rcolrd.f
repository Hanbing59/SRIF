      SUBROUTINE RCOLRD(S,MAXS,IRS,JCS,NPSTRT,NP,EM,RW,ZW,V,SGSTAR)
C
C     TO ADD IN MARKOV COLORED NOISE MODEL EFFECTS INTO THE SQUARE ROOT
C     INFORMATION FILTER, AND TO GENERATE WEIGHTING COEFFICIENTS
C     FOR SMOOTHING. IT IS ASSUMED THAT VARIABLES X(NPSTRT),
C     X(NPSTRT+1),...,X(NPSTRT+NP-1) ARE COLORED NOISE AND THAT
C     EACH COMPONENT SATISFIES A MODEL EQUATION OF THE FORM
C     X(SUB)(J+1)=EM*X(SUB)(J)+W(SUB)(J). FOR DETAILS, SEE
C     'FACTORIZATION METHODS FOR DISCRETE SEQUENTIAL
C     ESTIMATION', G.J.BIERMAN, ACADEMIC PRESS (1977)
C     FOR SMOOTHING, COEFFICIENTS TO BE GENERATED
C     SGSTAR(1) .GT. 0.0 .  THE SIGNIFICANCE OF THE
C     SMOOTHING COEFFICIENTS IS EXPLAINED IN THE FUNCTIONAL
C     DESCRIPTION.
C
C *** INPUTS:
C
C     S(IRS,JCS)  INPUT SQUARE ROOT INFORMATION ARRAY. OUTPUT
C                 COLORED NOISE ARRAY HOUSED HERE TOO. WHEN
C                 THERE IS SMOOTHING, NR ADDITIONAL ROWS MUST BE
C                 INCLUDED IN S
C
C     MAXS        ROW DIMENSION OF S. WHEN THERE ARE SMOOTHING
C                 COMPUTATIONS IT IS NECESSARY THAT MAXS.GE.IRS+NP
C                 BECAUSE THE BOTTOM NP ROWS OF S HOUSE THE
C                 SMOOTHING INFORMATION
C
C     IRS         NUMBER OF ROWS OF S (.LE. NUMBER OF FILTER
C                 VARIABLES) (IRS.GE.1)
C
C     JCS         NUMBER OF COLUMNS OF S (EQUALS NUMBER OF FILTER
C                 VARIABLES + POSSIBLY A RIGHT SIDE), WHICH CONTAINS
C                 THE DATA EQUATION NORMALIZED ESTIMATE (JCS.GE.1)
C
C     NPSTRT      LOCATION OF THE FIRST COLORED NOISE VARIABLE
C                 (1.LE.NPSTRT.LE.JCS)
C
C     NP          NUMBER OF CONTIGUOUS COLORED NOISE VARIABLES
C                 (NP.GE.1)
C
C     EM(NP)      COLORED NOISE MAPPING COEFFICIENTS (GENERALLY
C                 OF EXPONENTIAL FORM, EM=EXP(-DT/TAU))
C
C     RW(NP)      RECIPROCAL PROCESS NOISE STANDARD DEVIATIONS
C                 (MUST BE POSITIVE)
C
C     ZW(NP)      ZW=RW*W-ESTIMATE (PROCESS NOISE ESTIMATES ARE
C                 GENERALLY ZERO MEAN).  WHEN THE RIGHT HAND SIDE
C                 COLUMN IS OMITTED, ZW IS NOT USED.
C
C     SGSTAR(1)   IF SGSTAR(1) .GT. 0.0
C                    SMOOTHING COEFFICIENTS ARE GENERATED.
C
C                 IF SGSTAR(1) .LE. 0.0
C                    NO USE IS MADE OF SGSTAR OR OF ROWS
C                    IRS+1,...IRS+NP OF S.
C
C *** OUTPUTS:
C
C     V(IRS)      WORK VECTOR
C
C     S(IRS,JCS)  RESULTING SQUARE-ROOT INFORMATION ARRAY
C
C     SGSTAR(NP)  VECTOR OF SMOOTHING COEFFICIENTS. WHEN
C                 SGSTAR(1) .LE. 0.0 SGSTAR IS NOT USED.
C                 INCIDENTALLY, WHEN SMOOTHING THE SGSTAR
C                 MULTIPLIERS ARE ALWAYS POSITIVE.
C
C *** REQUIRED SUBROUTINES:
C
C     DSQRT       < FORLIB >
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              APRIL 10,1981
C
C              REVISED 30 MARCH  1982
C              REVISED 28 MAY    1982     VERSION 0.01.00
C              REVISED 10 AUGUST 1982     VERSION 0.01.01
C                                         VERSION 1.00.00
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE RCOLRD IS A PART OF THE FEA. INC             *
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

      DIMENSION S(MAXS,JCS),EM(NP),RW(NP),ZW(NP), V(IRS),SGSTAR(NP)
      DOUBLE PRECISION ALPHA,SIGMA,BETA,GAMMA, TMP, ZERO, ONE
      LOGICAL ISMTH
C
      DATA ZERO/0.0D0/, ONE/1.0D0/
C
      ISMTH = .FALSE.
      IF (SGSTAR(1) .LE. 0.0) ISMTH = .TRUE.
C
C                                  COL NO OF COLORED
C                                  NOISE TERM TO BE OPERATED ON
      NPCOL=NPSTRT
C
      DO 70 JCOLRD=1,NP
         TMP   = -RW(JCOLRD)*EM(JCOLRD)
         ALPHA = TMP
         SIGMA = ZERO
         DO 10 K=1,IRS
C
C                                 FIRST IRS ELEMENTS OF HOUSEHOLDER
C                                 TRANSFORMATION VECTOR
C
            V(K)=S(K,NPCOL)
   10     SIGMA=SIGMA+V(K)**2
C
C     IPJC  = I+JC
C
      IPJC = IRS + JCOLRD
      IF ((SIGMA .GT. ZERO) .OR. (.NOT. ISMTH)) GO TO 19
         DO 15 KOL = 1, JCS
            S(IPJC,KOL) = ZERO
  15     CONTINUE
         S(IPJC,NPCOL) = RW(JCOLRD)
         GO TO 65
  19  CONTINUE
          SIGMA=DSQRT(SIGMA+ALPHA*ALPHA)
C
C                                  LAST ELEMENT OF HOUSEHOLDER
C                                  TRANSFORMATION VECTOR
          ALPHA=ALPHA-SIGMA
C
C                                  HOUSEHOLDER=I+BETA*V*V**T
          BETA=ONE/(SIGMA*ALPHA)
C
C  HOUSEHOLDER TRANSFORMATION DEFINED, NOW APPLY IT TO S,
C  I.E. "60"  LOOP
C
        DO 60 KOL=1,JCS
           IF (KOL.NE.NPCOL) GO TO 30
              GAMMA= RW(JCOLRD)*ALPHA*BETA
              IF (ISMTH)   S(IPJC,NPCOL)=TMP*GAMMA
              DO 20 K=1,IRS
   20            S(K,NPCOL)=GAMMA*V(K)
              GO TO 60
   30      GAMMA=ZERO
           IF (KOL.EQ.JCS) GAMMA=ZW(JCOLRD)*ALPHA
C
C         IF ZW ALWAYS ZERO OR THERE ARE NO ESTIMATES
C         ONE CAN COMMENT OUT THE ABOVE IF STATEMENT
C
          DO 40 K=1,IRS
   40       GAMMA=GAMMA+S(K,KOL)*V(K)
            GAMMA= GAMMA*BETA
          DO 50 K=1,IRS
   50       S(K,KOL)=S(K,KOL)+GAMMA*V(K)
            IF (ISMTH) S(IPJC,KOL)=GAMMA*ALPHA
   60   CONTINUE
        IF ( .NOT. ISMTH ) GO TO 70
   65      SGSTAR(JCOLRD) = SIGMA
           S(IPJC,JCS) = S(IPJC,JCS) + ZW(JCOLRD)
C
C     IF ZW IS ALWAYS ZERO OR THERE ARE NO ESTIMATES INVOLVED,
C     ONE CAN COMMENT OUT THE ABOVE STATEMENT
C
   70 NPCOL=NPCOL+1
C
      RETURN
      END
