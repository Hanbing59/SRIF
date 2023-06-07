      SUBROUTINE U2D(U,N,D,JST)
C     TO COPY DIAGONAL ELEMENTS OF VECTOR STORED  MATRIX U
C     INTO VECTOR D, STARTING WITH ELEMENT (JSTRT,JSTRT)
C
C *** INPUTS
C
C     U(N*(N+1)/2)   VECTOR STORED U MATRIX
C
C     N              COLUMN DIMENSION OF U
C
C     JST            LOCATION OF FIRST ELEMENT TO BE COPIED.
C                    JSTRT .LE.N
C
C *** OUTPUT
C
C     D(N-JST+1)     VECTOR OF OUTPUT DIAGONAL ENTRIES
C                    U(JST,JST)...U(N,N)
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JUNE 11, 1981
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE D2U    IS A PART OF THE FEA. INC             *
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
 
      DIMENSION U(*),D(*)
C
      JSTM1=JST-1
      JJ=JST*(JSTM1)/2
C
      DO 10 J=JST,N
         JJ   = JJ+J
         JTMP = J-JSTM1
         D(JTMP)=U(JJ)
 10   CONTINUE
      RETURN
      END
