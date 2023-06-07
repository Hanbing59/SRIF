      SUBROUTINE TPHIU(PHI,MRPHI,IRPHI,JCPHI,U,N)
C
C
C       TO COMPUTE PHI*U AND OVERWRITE U WITH THE RESULTS.  PHI
C       IS ASSUMED TO BE UNIT UPPER TRIANGULAR WITH IMPLICIT
C       TRAILING COLUMNS OF ZEROS.  THE LOWER PORTION OF PHI
C       (IE. THE ELEMENTS BELOW ROW JCPHI) ARE ONLY IMPLICITLY
C            DEFINED, AS (0 I).
C
C ** INPUTS:
C
C       PHI(IRPHI,JCPHI)   INPUT MATRIX OF TRANSITON MULTIPLIER
C                   INPUTS. ONLY THE ELEMENTS PHI(I,J) FOR
C                          J .GT. I ARE EXPLICITLY USED.
C                          PHI IS LEFT UNCHANGED.
C
C       MRPHI              MAX DIMENSION OF PHI MATRIX
C
C       IRPHI              NUMBER OF ROWS IN PHI (IRPHI .LE. N)
C
C       JCPHI              NUMBER OF COLUMNS IN PHI (JCPHI .LE. N)
C
C       U(N*(N+1)/2)       VECTOR STORED U-D MATRIX WITH
C                          IMPLICITLY DEFINED UNIT DIAGONAL
C                          ELEMENTS
C
C       N                  DIMENSION OF U WHEN VIEWED AS A DOUBLE
C                          SUBSCRIPTED SQUARE ARRAY.
C
C ** OUTPUT:
C
C       U(N*(N+1)/2)       VECTOR STORED U-D MATRIX WITH
C                          IMPLICITLY DEFINED UNIT DIAGONAL
C                          ELEMENTS
C
C
C     COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              AUGUST 31, 1981
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE TPHIU  IS A PART OF THE FEA. INC             *
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

      DIMENSION PHI(MRPHI,JCPHI),U(*)
C
        II=0
        ILAST=IRPHI
        IF (ILAST .EQ. N) ILAST=N-1
C
        DO 40 I=1,ILAST
C
C     II=(I,I) ELEMENT OF VECTOR STORED U MATRIX
C
           II=II+I
C
C     KJU=(I,I+1) ELEMENT OF VECTOR STORED U MATRIX
C
           KJU=II+I
           IP1=I+1
           DO 30 J=IP1,N
C
C  S=U(I,J)+SUM(PHI(I,K)*U(K,J);K=I+1,J-1)+PHI(I,J)
C
              S=U(KJU)
              IF (J .EQ. IP1) GO TO 20
C
C  SUM(PHI(I,K)*U(K,J);K=I+1,J-1)
C
              KJ=KJU
              JM1=J-1
              IF (J .GT. JCPHI) JM1=JCPHI
C
C   JM1 RESET TO ACCOUNT FOR TRAILING ZERO
C   COLUMNS IN PHI MATRIX
C
         IF (IP1 .GT. JM1) GO TO 20
C
           DO 10 K=IP1,JM1
                 KJ=KJ+1
                 S=S+PHI(I,K)*U(KJ)
  10          CONTINUE
  20          IF (J .LE. JCPHI) S=S+PHI(I,J)
C
C     IF J .GT. JCPHI THEN PHI*U PRODUCT DOES
C            NOT CONTAIN UNIT U(J,J) TERM
C
              U(KJU)=S
              KJU=KJU+J
  30    CONTINUE
  40   CONTINUE
       RETURN
       END
