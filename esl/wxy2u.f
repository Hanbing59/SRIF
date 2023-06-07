      SUBROUTINE WXY2U(WXY,IWMAX,IROW,U,N)
C     TO COPY WXY INTO UXY WHERE UXY IS A SYMBOLIC PARTITION OF
C     VECTOR STORED U AND WXY, IS A SYMBOLIC PARTITION OF  MATRIX W
C     IT IS ASSUMED THAT WXY IS BEING COPIED INTO 'COLUMNS'
C     IROW+1,..N OF U
C
C *** INPUTS
C
C     WXY(IROW,N-IROW) MATRIX TO BE COPIED. IN PRACTICE, THE CALLING
C                    SEQUENCE WOULD USE W(1,IROWP1),IROWP1=IROW+1; AND
C                    THIS  IS THE REASON THAT WE CALL WXY A
C                    "SYMBOLIC" PARTITION.
C
C     IWMAX          ROW DIMENSION OF WXY
C
C     IROW           NUMBER OF ROWS TO BE COPIED IN U.
C                    IROW .LT. N
C
C     U(N*(N+1)/2)   VECTOR STORED U-D COVARIANCE FACTORS WITH D
C                    ELEMENTS STORED ON THE DIAGONAL
C
C     N              COLUMN DIMENSION OF U
C
C
C *** OUTPUT
C
C     U(N*(N+1)/2)   U, WITH THE CORRESPONDING COLUMNS OF WXY SET
C                    INTO PLACE
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JUNE 12, 1981
C
C              REVISED  17 AUGUST 1982    VERSION 0.00.01
C                                         VERSION 1.00.00
C              REVISED 20  MAY    1985    VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE WXY2U  IS A PART OF THE FEA. INC             *
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
      DIMENSION U(*),WXY(IWMAX,*)
C
      IF (IROW .GE. N) RETURN
      IRP1 = IROW+1
      JJ   = IROW*IRP1/2
      DO 20 J=IRP1,N
         DO 10 I=1,IROW
            JJI    = JJ+I
            JMIR   =J-IROW
            U(JJI) = WXY(I,JMIR)
10    CONTINUE
      JJ = JJ+J
20    CONTINUE
      RETURN
      END
