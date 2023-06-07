      SUBROUTINE SFU(FEL,IROW,JCOL,NF,U,N,W,IMAXW,IRW,JDIAG)
C
C         TO COMPUTE W(IRW,N)=F*U WHERE F IS SPARSE AND ONLY THE
C         NON-ZERO ELEMENTS ARE DEFINED AND U IS VECTOR STORED,
C         UPPER TRIANGULAR WITH IMPLICITLY DEFINED UNIT DIAGONAL
C         ELEMENTS
C
C *** INPUTS:
C
C     FEL(NF)       VALUES OF THE NON-ZERO ELEMENTS OF THE F MATRIX
C
C     IROW(NF)      ROW INDICES OF THE F ELEMENTS
C
C     JCOL(NF)      COLUMN INDICES OF THE F ELEMENTS
C                   F(IROW(K),JCOL(K))=FEL(K)
C
C     NF            NUMBER OF NON-ZERO ELEMENTS OF THE F MATRIX
C
C     U(N*(N+1)/2)  UPPER TRIANGULAR, VECTOR STORED MATRIX WITH
C                   IMPLICITLY DEFINED UNIT DIAGONAL ELEMENTS
C                   (U(J,J) ARE NOT, IN FACT, UNITY)
C
C     N             DIMENSION OF U MATRIX
C
C
C
C     JDIAG(N)      DIAGONAL ELEMENT INDICES OF A VECTOR STORED
C                   UPPER TRIANGULAR MATRIX,
C                   I.E. JDIAG(K)=K*(K+1)/2=JDIAG(K-1)+K
C
C *** OUTPUT:
C
C     W(IRW,N)      OUTPUT RESULT
C
C     IMAXW         ROW DIMENSION OF W MATRIX
C
C     IRW           NUMBER OF ROWS IN W,
C                   (IRW.LE.IMAXW.AND.IRW.GE.MAX(IROW(K)),
C                   K=1,...,NF, I.E. W MUST HAVE AT LEAST AS
C                   MANY ROWS AS DOES F.  ADDITIONAL ROWS OF W
C                   COULD CORRESPOND TO ZERO  ROWS OF F.
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              OCTOBER 23,1981
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE SFU    IS A PART OF THE FEA. INC             *
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

      DIMENSION  FEL(NF),U(*),W(IMAXW,N),IROW(NF)
      DIMENSION  JCOL(NF),JDIAG(N)
C
      ZERO=0.D0
C
C * * * * INITIALIZE W
C
      DO 10 J=1,N
        DO 10 I=1,IRW
   10     W(I,J)=ZERO
C
C         IF IMAXW=IRW, IT IS MORE EFFICIENT TO REPLACE THIS LOOP BY
C
C             DO 10 IJ=1,IFUN
C          10   W(IJ,1)=ZERO
C                                 WHERE  IFUN=IRW*N
C
      DO 30 NEL=1,NF
C
C         NEL REPRESENTS THE ELEMENT NUMBER OF THE F MATRIX
C
         I=IROW(NEL)
         J=JCOL(NEL)
         FIJ=FEL(NEL)
         W(I,J)=W(I,J)+FIJ
C
C           THIS ACCOUNTS FOR THE IMPLICIT UNIT DIAGONAL U MATRIX
C           ELEMENTS. WHEN NON-UNIT DIAGONALS ARE USED, DELETE
C           THE ABOVE LINE AND USE J INSTEAD OF JP1 BELOW
C
         IF (J.EQ.N) GO TO 30
C
C           WHEN IT IS KNOWN THAT THE LAST COLUMN OF F IS ZERO
C           THIS 'IF' TEST MAY BE OMITTED
C
        JP1=J+1
        IK=JDIAG(J)+J
        DO 20 K=JP1,N
           W(I,K)=W(I,K)+FIJ*U(IK)
   20      IK=IK+K
   30 CONTINUE
      RETURN
      END
