      SUBROUTINE UCOL(U,N,KS,NCOLOR,V,EM,Q,IERR)
C
C
C     COLORED NOISE UPDATING OF THE U-D COVARIANCE FACTORS, I.E.
C       U*D*(U**T)-OUTPUT=PHI*U*D*(U**T)*(PHI**T)+Q
C       PHI=DIAG(E(KS-1),EM(1),...,EM(NCOLOR),E(N-(KS-1+NCOLOR)))
C     Q=DIAG(0(KS-1),Q(1),...,Q(NCOLOR),0(N-(KS-1+NCOLOR)))
C     E(K) IS A K  VECTOR OF ZEROS, AND E(K) IS UNITY IN ELEMENT
C     K, ZEROS ELSEWHERE
C
C     THE ALGORITHM USED IS THE BIERMAN-THORNTON ONE COMPONENT
C     AT-A-TIME UPDATE. CF.BIERMAN "FACTORIZATION METHOD
C     FOR DISCRETE SEQUENTIAL ESTIMATION", ACADEMIC PRESS (1977)
C     PP.147-148
C
C *** INPUTS:
C
C     U(N*(N+1)/2)      INPUT U-D VECTOR STORED COVARIANCE FACTORS.
C                       THE COLORED NOISE UPDATE RESULT RESIDES
C                       IN U ON OUTPUT
C
C     N                 FILTER DIMENSION. IF THE LAST COLUMN OF U
C                       HOUSES THE FILTER ESTIMATES, THEN
C                       N=NUMBER FILTER VARIABLES + 1
C
C     KS                THE LOCATION OF THE FIRST COLORED NOISE TERM
C                       (KS.GE.1.AND.KS.LE.N)
C
C     NCOLOR            THE NUMBER OF COLORED NOISE TERMS
C                       (NCOLOR.GE. 1)
C
C     V(KS-1+NCOLOR)    WORK VECTOR
C
C     EM(NCOLOR)        INPUT VECTOR OF COLORED NOISE MAPPING TERMS
C                       (UNALTERED BY PROGRAM)
C
C     Q(NCOLOR)         INPUT VECTOR OF PROCESS NOISE VARIANCES
C                       (UNALTERED BY PROGRAM)
C
C *** OUTPUTS:
C
C     U(N*(N+1)/2)      THE COLORED NOISE UPDATE RESULT RESIDES
C                       IN U ON OUTPUT
C
C     IERR            DIAGNOSTIC FLAG (FROM RNK1)
C
C                     IF  .EQ.  0       NORMAL RETURN
C
C                               J       ERROR RETURN DUE TO A
C                                       NEGATIVE COMPUTED J-TH
C                                       DIAGONAL ENTRY
C
C                              -1       UD RESULT IS SINGULAR
C
C ***  REQUIRED SUBROUTINE:
C
C               RNK1  <ESL>
C
C
C     COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              SEPTEMBER 8, 1981
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE UCOL   IS A PART OF THE FEA. INC             *
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

      DOUBLE PRECISION TMP,S
      DIMENSION U(*),V(*),EM(NCOLOR),Q(NCOLOR)
C
C * * * * * * INITIALIZATION
C
      NM1   = N-1
      KSM1  = KS-1
      JJOLD = KS*KSM1/2
      KOL   = KSM1
C
C * * * * * *
C
      DO 50 K=1,NCOLOR
        KOLM1 = KOL
        KOL   = KOL+1
        JJ    = JJOLD+KOL
        C     = U(JJ)
        A     = EM(K)
        TMP   = C*A
        S=TMP*A+Q(K)
C                           D(J) UPDATE
        U(JJ)=S
C
        IF (KOL.GE.N) GO TO 20
           IJ=JJ
           DO 10 J=KOL,NM1
              IJ   = IJ+J
              U(IJ)= U(IJ)*A
   10      CONTINUE
C                     UPDATING ROW KOL ENTRIES
C
   20   IF (JJ.EQ.1) GO TO 50
C                     (WHEN KS=1)
           IF (S.LE.0.D0) GO TO 30
              C   = C/S
              TMP = C*A
C                        TMP=EM(K)*D(KOL)-OLD/D(KOL)-NEW
   30      C   = C*Q(K)
C                        C=Q(K)*D(KOL)-OLD/D(KOL)-NEW
           DO 40 I=1,KOLM1
              IFIX=JJOLD+I
              V(I)=U(IFIX)
   40      U(IFIX)=TMP*V(I)
   45      CALL RNK1(U,KOLM1,C,V,IERR)
   50     JJOLD=JJ
      RETURN
      END
