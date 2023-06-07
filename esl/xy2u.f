      SUBROUTINE XY2U(XY2UF,U,N,ISTART,ISTOP,JSTART,W2D,IWMAX)
C
C   WHEN XY2UF(LAG) IS TRUE W2D (DOUBLE SUBSCRIPT) IS COPIED INTO
C   THE SYMBOLIC PARTITION UXBY OF VECTOR STORED U (SEE FIGURE).
C   WHEN XY2UF(LAG) IS FALSE THE OPPOSITE COPY TAKES PLACE
C   I.E. THE SYMBOLIC U ARRAY PARTITION IS COPIED INTO W2D.
C
C                           ISTOP        JSTART            N
C    *************************$*************$***************
C     *                       $             $              *
C       *                     $             $              *
C         *                   $             $              * XA
C           *                 $             $              *
C             *               $             $              *
C               *             $             $              *
C                 %%%%%%%%%%%%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX ISTART
C                 % *         X             X   U          X
C                 %   *       X             X    XBY       X
C                 %     *     X             X              X XB
C                 %       *   X             X   W2D        X
C                 %         * X             X              X
C                 %%%%%%%%%%%%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX ISTOP
C                               *                          *
C                                 *                        *
C                                   *                      *
C                                     *                    *
C                                       *                  *
C                                         *                *
C                                           *              * Y
C                                             *            *
C                                               *          *
C                                                 *        *
C                                                   *      *
C                                                     *    *
C                                                       *  *
C                                                          *
C
C *** INPUTS
C   XY2UF            LOGICAL
C
C                    "TRUE"  IF DOUBLE SUBSCRIPT ARRAY W2D IS TO
C                            BE COPIED INTO APPROPRIATE BLOCK OF U
C                    "FALSE" IF VECTOR STORED BLOCK OF U IS TO BE
C                            COPIED INTO DOUBLE SUBSCRIPT ARRAY W2D
C IF (XY2UF IS TRUE)
C
C   W2D(IROW,JCOL)   MATRIX TO BE COPIED INTO U. IROWS=ISTOP-ISTART,
C                    JCOL=N-JSTART+1. IN PRACTICE, THE CALLING
C                    SEQUENCE WOULD USE W(ISTART,JSTART), AND THIS IS
C                    THE REASON THAT WXBY IS CALLED A "SYMBOLIC"
C                    PARTITION.
C
C   IWMAX            ROW DIMENSION OF W2D
C
C IF (XY2UF IS FALSE)
C
C   N                COLUMN DIMENSION OF U
C
C   U(N*(N+1)/2)     VECTOR STORED UPPER TRIANGULAR INPUT
C                    SYMBOLIC PARTITION TO GO INTO W2D
C
C   ISTART           START ROW OF SYMBOLIC PARTITION
C
C   ISTOP            LOCATION OF THE BOTTOMMOST ROW OF
C                    SYMBOLIC PARTITION
C
C   JSTART           STARTING COLUMN OF SYMBOLIC PARTITION
C
C
C *** OUTPUT
C
C IF XY2UF = .TRUE.
C
C   U(N*(N+1)/2)     U, WITH THE CORRESPONDING COLUMNS OF W2D SET
C                    INTO PLACE
C
C ELSE
C
C   W2D              CORRESPONDING XBY BLOCK OF U IS COPIED INTO W2D
C
C
C     COGNIZANT PERSONS:
C
C         DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C         FACTORIZED ESTIMATION APPLICATIONS
C         JULY 3, 1986
C
C         REVISED 1 DEC 86                      VERSION 1.00.00
C
C ***************************************************************
C *                                                             *
C *           SUBROUTINE XY2U IS PART OF THE FEA. INC           *
C *                ESTIMATION SUBROUTINE LIBRARY                *
C *                                                             *
C *                           (ESL)                             *
C *                                                             *
C *            COPYRIGHT 1982, 1983, 1984, 1985, 1986           *
C *            FACTORIZED ESTIMATION APPLICATIONS INC.          *
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
      DIMENSION U(*), W2D(IWMAX,*)
      LOGICAL XY2UF
C
      IROW  = ISTOP  - ISTART+1
      JSTM1 = JSTART -   1
      JJ    = JSTART*JSTM1/2
      IF (XY2UF) THEN
         DO 20 J = JSTART,N
            IJ   = ISTART+JJ
            DO 10 I=1,IROW
               U(IJ) = W2D(I,J-JSTM1)
               IJ    = IJ+1
   10       CONTINUE
            JJ = JJ+J
   20    CONTINUE
      ELSE
         DO 40 J = JSTART,N
            IJ   = ISTART+JJ
            DO 30 I=1,IROW
               W2D(I,J-JSTM1) = U(IJ)
               IJ = IJ+1
   30       CONTINUE
            JJ = JJ+J
   40    CONTINUE
      ENDIF
      RETURN
      END
