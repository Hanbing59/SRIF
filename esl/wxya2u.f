      SUBROUTINE WXYA2U(WXYA,IWMAX,IROW,JYA,U,N)
C     TO COPY WXYA INTO UXYA WHERE UXYA IS A SYMBOLIC PARTITION OF
C     VECTOR STORED U AND WXYA, IS A SYMBOLIC PARTITION OF  MATRIX W
C     IT IS ASSUMED THAT WXYA IS BEING COPIED INTO 'ROWS'
C     (JYA-IROW+1, ..., JYA) AND 'COLUMNS' JYA+1,..N OF U
C
C                           JYA                     N
C     ***********************$***********************
C      *                     $                      *
C        *                   $                      *
C          *                 $                      *
C            *               $                      *
C              *             $                      *
C                *           $                      *
C                  %%%%%%%%%%XXXXXXXXXXXXXXXXXXXXXXXX JYA-IROW+1
C                  % *       X                      X
C                  %   *     X                      X
C                  %     *   X       WXYA           X
C                  %       * X                      X
C                  %         X                      X
C                  %%%%%%%%%%XXXXXXXXXXXXXXXXXXXXXXXX JYA
C                              *                    *
C                                *                  *
C                                  *                *
C                                    *              *
C                                      *            *
C                                        *          *
C                                          *        *
C                                            *      *
C                                              *    *
C                                                *  *
C                                                   *
C
C *** INPUTS
C
C     WXYA(IROW,N-JYA) MATRIX TO BE COPIED. IN PRACTICE, THE CALLING
C                      SEQUENCE WOULD USE W(1,JYA+1),AND THIS  IS
C                      THE REASON THAT WE CALL WXYA A "SYMBOLIC"
C                      PARTITION.
C
C     IWMAX          ROW DIMENSION OF WXYA
C
C     IROW           NUMBER OF ROWS TO BE COPIED INTO U.
C                    IROW .LT. N
C
C     JYA            LOCATION OF THE BOTTOMMOST ROW OF WXYA;
C                    JYA .GE. IROW
C
C     U(N*(N+1)/2)   VECTOR STORED U-D COVARIANCE FACTORS WITH D
C                    ELEMENTS STORED ON THE DIAGONAL
C
C     N              COLUMN DIMENSION OF U
C
C
C *** OUTPUT
C
C     U(N*(N+1)/2)   U, WITH THE CORRESPONDING COLUMNS OF WXYA SET
C                    INTO PLACE
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              NOVEMBER 18, 1984
C
C              REVISED 20      MAY 1985         VERSION 2.00.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE WXYA2U  IS A PART OF THE FEA. INC            *
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
      DIMENSION U(*), WXYA(IWMAX,*)
C
      IF ((JYA .GE. N) .OR. (JYA .LT. IROW)) RETURN
      JYAP1   = JYA + 1
      JJ      = JYA *(JYAP1)/2
      ISTRT   = JYA - IROW
      DO 20 J = JYAP1,N
         IJ   = ISTRT + JJ
         DO 10 I = 1, IROW
            IJ   = IJ+1
            JMIR = J-JYA
            U(IJ)= WXYA(I,JMIR)
10       CONTINUE
         JJ=JJ+J
20    CONTINUE
      RETURN
      END
