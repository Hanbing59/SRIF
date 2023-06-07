      SUBROUTINE UD2COV (UIN,POUT,N)
C
C     TO OBTAIN A COVARIANCE FROM ITS U-D FACTORIZATION. BOTH
C     MATRICES ARE VECTOR STORED AND THE OUTPUT COVARIANCE CAN
C     OVERWRITE THE INPUT U-D ARRAY. UIN=U-D IS RELATED TO POUT
C     VIA POUT=UDU(**T)
C
C     UIN(N*(N+1)/2)  INPUT U-D FACTORS, VECTOR STORED WITH
C                     THE D ENTRIES STORED ON THE DIAGONAL OF
C                     UIN
C
C     POUT(N*(N+1)/2) OUTPUT COVARIANCE, VECTOR STORED.
C                     IF THE COMPILER ALLOWS POUT MAY
C                     OVERWRITE UIN.
C
C     N               DIMENSION OF THE MATRICES INVOLVED,
C                     N .GE. 1
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JUNE 7, 1981
C
C              REVISED 13 AUGUST 1982
C                                         VERSION 1.00.00
C              REVISED 20    MAY 1985     VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE UD2COV IS A PART OF THE FEA. INC             *
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

      DIMENSION UIN(*),   POUT(*)
C
      POUT(1)=UIN(1)
      IF (N .EQ. 1) RETURN
      JJ=1
      DO 20 J=2,N
        JJL=JJ
C                              (J-1,J-1)
        JJ=JJ+J
        POUT(JJ)=UIN(JJ)
        S=POUT(JJ)
        II=0
        JM1=J-1
        DO 20 I=1,JM1
          II=II+I
C                              JJL+I=(I,J)
          ALPHA=S*UIN(JJL+I)
          IK=II
          DO 10 K=I,JM1
C                                               JJL+K=(K,J)
            POUT(IK)=POUT(IK)+ALPHA*UIN(JJL+K)
   10       IK=IK+K
   20   POUT(JJL+I)=ALPHA
C
      RETURN
      END
