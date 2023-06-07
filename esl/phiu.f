      SUBROUTINE PHIU(PHI,MAXPHI,IRPHI,JCPHI,U,N,W,IMAXW)
C
C     THIS SUBROUTINE  COMPUTES  W=PHI*U  WHERE PHI IS A
C     RECTANGULAR MATRIX WITH IMPLICITLY DEFINED COLUMNS
C     OF  TRAILING  ZEROS AND U IS A VECTOR STORED UPPER
C     TRIANGULAR MATRIX.
C
C *** INPUTS:
C
C     PHI(IRPHI,JCPHI) INPUT RECTANGULAR MATRIX, IRPHI.LE.MAXPHI
C
C     MAXPHI           ROW DIMENSION OF PHI
C
C     IRPHI            NO. ROWS OF PHI
C
C     JCPHI            NO. COLS OF PHI
C
C     U(N*(N+1)/2)     UPPER TRIANGULAR VECTOR STORED MATRIX
C
C     N                DIMENSION OF U MATRIX (JCPHI.LE.N)
C
C *** OUTPUTS:
C
C     W(IRPHI,N)       OUTPUT, RESULT OF PHI*U, IF THE
C                      COMPILER PERMITS, W CAN OVERWRITE
C                      PHI
C
C     IMAXW            ROW DIMENSION OF W
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JUNE 11, 1981
C              REVISED 18 OCTOBER 1981
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C              REVISED 20      OCT 1986         VERSION 2.01.00
C              REVISED 02      JAN 1997: removed recursion on KJS for parallel
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE PHIU   IS A PART OF THE FEA. INC             *
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

      DIMENSION PHI(MAXPHI,JCPHI),U(*),W(IMAXW,N)
      DOUBLE PRECISION SUM
C
      DO 10 I=1,IRPHI
         W(I,1)=PHI(I,1)
  10  CONTINUE
      IF (N .EQ. 1) RETURN
      NP2=N+2
CDIR$ IVDEP
      DO 40 L=2,N
        J=NP2-L
        KJS=(J-1)*J/2
        JM1=J-1
        DO 30 I=1,IRPHI
           IF (J.LE.JCPHI) THEN
               SUM=PHI(I,J)
           ELSE
               SUM=0.D0
               JM1=JCPHI
           ENDIF
           DO 20 K=1,JM1
              SUM=SUM+PHI(I,K)*U(KJS+K)
  20       CONTINUE
           W(I,J)=SUM
  30    CONTINUE
  40  CONTINUE
      RETURN
      END
