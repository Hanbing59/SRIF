      SUBROUTINE RA (R,N,A,IMAXA,IRA,JCA,W,IMAXW,IRW)
C
C         TO COMPUTE  IRW ROWS OF THE PRODUCT W=R*A
C
C     WHERE R IS UPPER TRIANGULAR VECTOR SUBSCRIPTED AND OF
C     DIMENSION N, A HAS JCA COLUMNS AND IRA ROWS, THE TRAILING
C     N-JCA COLUMNS ARE IMPLICITY DEFINED AS ZERO. WHEN
C     IRA .LT. N THE BOTTOM N-IRA ROWS OF A ARE ASSUMED TO BE
C     IMPLICITLY DEFINED AS THE BOTTOM N-IRA ROWS OF THE
C     N DIMENSION IDENTITY MATRIX.
C
C *** INPUTS:
C
C     R(N*(N+1)/2)  UPPER TRIANGULAR VECTOR STORED INPUT MATRIX
C
C     N             DIMENSION OF R.
C
C     A(IRA,JCA)    INPUT RECTANGULAR MATRIX
C
C     IMAXA         ROW DIMENSION OF A
C
C     IRA           NUMBER OF ROWS IN THE A MATRIX
C                   (IRA.LE.IMAXA AND IRA .LE. N)
C
C     JCA           NUMBER OF COLUMNS IN THE A MATRIX, JCA .LE. N
C
C     IMAXW         ROW DIMENSION OF W
C
C     IRW           NUMBER OF ROWS OF THE PRODUCT R*A THAT ARE
C                   COMPUTED (IRW.LE.IMAXW AND IRW .LE. N)
C
C *** OUTPUT:
C
C     W(IRW,N)      RESULTING RECTANGULAR MATRIX, = R*A
C                   IF THE COMPILER ALLOWS IT W CAN OVERWRITE A
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23, 1981
C
C              REVISED 20 SEPTEMBER 1982    VERSION 1.00.00
C              REVISED 20      MAY 1985     VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE RA     IS A PART OF THE FEA. INC             *
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

      DIMENSION  R(*),A(IMAXA,JCA),W(IMAXW,N)
      DOUBLE PRECISION SUM
C
      IJ=IRA*(IRA+1)/2
C
      DO 40 J=1,N
        II=0
        DO 30 I=1,IRW
          II=II+I
          SUM=0.D0
          IF (I.GT.IRA .OR. J .GT. JCA) GO TO 15
              IK=II
              DO 10 K=I,IRA
                 SUM=SUM+R(IK)*A(K,J)
                 IK=IK+K
  10          CONTINUE
  15      CONTINUE
          IF (J .LE. IRA .OR. I .GT. J) GO TO 20
             SUM = SUM + R(IJ+I)
  20      CONTINUE
          W(I,J)=SUM
  30    CONTINUE
        IF (J.GT.IRA) IJ=IJ+J
  40  CONTINUE
      RETURN
      END
