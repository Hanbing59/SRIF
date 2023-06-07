      SUBROUTINE MAPU(PHI,MAXPHI,IRPHI,JCPHI,U,N,W,MAXW,D,POUT)
C
C     TO COMPUTE MAPPED COVARIANCE POUT = (PHI*U)*D*(PHI*U)**T
C
C        IT IS ASSUMED THAT THE MAPPING MATRIX IS OF THE FORM:
C
C
C
C                   PHI   0         WITH ONLY THE PHI
C                                   ENTRIES EXPLICIT
C    *** INPUTS:
C
C     PHI(IRPHI,JCPHI)        MAPPING MATRIX
C
C     MAXPHI                  ROW DIMENSION OF PHI,
C                             IRPHI .LE. MAXPHI
C
C     U(N*(N+1)/2)            VECTOR STORED U-D COVARIANCE FACTORS
C
C     N                       COLUMN DIMENSION OF U
C
C     W(IRPHI,N)              SCRATCH MATRIX. IF THE COMPILER ALLOWS
C                             W AND PHI CAN SHARE STORAGE.
C
C     MAXW                    ROW DIMENSION OF W (IRPHI .LE. MAXW)
C
C     D(1)                    FLAG TO INDICATE WHETHER ONLY DIAGONAL
C                             ENTRIES (SIGMAS) OF THE COVARIANCE
C                             BEING REQUESTED.
C
C                             IF (D(1) .GE. 0) COMPUTE ENTIRE
C                             COVARIANCE MATRIX, ELSE COMPUTE ONLY
C                             DIAGONAL COVARIANCE SQUARE-ROOT VALUES.
C
C    *** OUTPUTS:
C
C     POUT(IRPHI*(IRPHI+1)/2) IF (D(1) .GE. 0) THEN
C                                VECTOR STORED OUTPUT COVARIANCE
C                                MATRIX WHICH (IF THE COMPILER ALLOWS)
C                                MAY OVERWRITE U.
C                             ELSE
C                                POUT HAS ONLY IRPHI STANDARD
C                                DEVIATIONS.
C                             ENDIF
C
C  *** SCRATCH:
C
C     D(N)               DIAGONAL ENTRIES OF INPUT U MATRIX
C
C          COGNIZANT PERSONS:
C
C             DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C             FACTORIZED ESTIMATION APPLICATIONS INC.
C
C             AUGUST  10, 1981
C
C             REVISED 18 OCTOBER 1981
C             REVISED 31 JANUARY 1982
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C              REVISED 13      JAN 1986         VERSION 2.01.00
C
C     SUBROUTINE REQUIRED:
C
C        PHIU
C
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE MAPU   IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *           COPYRIGHT 1982, 1983, 1984, 1985, 1986               *
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

      DIMENSION PHI(MAXPHI,IRPHI),U(*),W(MAXW,N),POUT(*),D(N)
      LOGICAL WNTVAR
C
      DOUBLE PRECISION SUM
C
      WNTVAR = D(1) .LT. 0.0
C
      CALL PHIU(PHI,MAXPHI,IRPHI,JCPHI,U,N,W,MAXW)
C
C     STORE DIAGONALS OF U IN VECTOR D
C
      JJ=0
      DO 10 J=1,N
         JJ=JJ+J
         D(J)=U(JJ)
  10  CONTINUE
C
C     CREATE POUT=W*D*(W**T)
C
      IF (IRPHI .EQ. 1) GO TO 60
      IJ = IRPHI*(IRPHI+1)/2+1
      NP2= IRPHI + 2
      DO 50 IROW=2,IRPHI
         I=NP2-IROW
         SUM=0.0
         DO 20 K=1,N
            TMP=D(K)*W(I,K)
            SUM=SUM+W(I,K)*TMP
            W(I,K)=TMP
  20     CONTINUE
         IF (WNTVAR) THEN
            POUT(I) = DSQRT(SUM)
         ELSE
            IJ=IJ-1
            POUT(IJ)=SUM
C
C          DIAGONAL ELEMENT OF OUTPUT COVARIANCE
C
            IM1=I-1
            DO 40 L=1,IM1
               K=I-L
               SUM=0.0
               DO 30 J=1,N
                  SUM=SUM+W(K,J)*W(I,J)
  30           CONTINUE
               IJ=IJ-1
               POUT(IJ)=SUM
  40        CONTINUE
         ENDIF
  50  CONTINUE
  60  CONTINUE
      SUM=0.0
      DO 70 J=1,N
         SUM=SUM+D(J)*W(1,J)**2
  70  CONTINUE
      IF (WNTVAR) THEN
         POUT(1) = DSQRT(SUM)
      ELSE
         POUT(1)=SUM
      ENDIF
      RETURN
      END
