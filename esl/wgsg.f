      SUBROUTINE WGSG(W,IMAXW,IW,JW,DW,F,V,U)
C
C
C     MODIFIED WEIGHTED GRAM-SCHMIDT GENERALIZED TRIANGULARIZATION.
C     MODIFIED WEIGHTED GRAM-SCHMIDT ALGORITHM FOR TRANSFORMING
C     W-DW TO U-D WHERE UDU**T=W*DW*W**T AND EITHER:
C
C
C          **********                   *****
C      W = *        * IW --------> U =   *  *  IW .LE. JW
C          **********                     * *
C             JW                            *
C
C
C            OR
C
C
C                                       *********
C                                       *       * IW-JW   WHEN
C          ********                     *       *         IW .GT. JW
C          *      *                     *       *
C      W = *      * IW -------> U =     *--------
C          *      *                      *      *
C          *      *                       *     *
C          ********                         *   *
C             JW                              * *
C                                               *
C
C
C     REFERENCES:
C
C                (1)   EFFICIENT TIME PROPAGATION OF U-D COVARIANCE
C                      FACTORS, G.J. BIERMAN, IEEE TRANS. A/C
C                      VOL. AC-26, NO. 4,8-81, PP 890-894
C
C                (2)   FACTORIZATION METHODS FOR DESCRETE SEQUENTIAL
C                      ESTIMATION ACADEMIC PRESS, 1977, G.J. BIERMAN
C
C
C *** INPUTS:
C
C     W(IW,JW)         INPUT MATRIX TO BE REDUCED TO TRIANGULAR
C                      FORM. THIS MATRIX IS DESTROYED BY THE
C                      CALCULATION.
C
C     IMAXW            ROW DIMENSION OF W MATRIX.
C
C     IW               NUMBER OF ROWS OF W MATRIX = COLUMN DIMENSION
C                      OF OUTPUT U MATRIX.
C                      ( (IW .GE. 1) .AND. (IW .LE. IMAXW) )
C
C     JW                NUMBER OF COLUMNS OF W MATRIX. ( JW .GE. 1)
C
C     DW(JW)           VECTOR OF NON-NEGATIVE WEIGHTS FOR THE
C                      OTHOGONALIZATION PROCESS. THESE (DW) ELEMENTS
C                      ARE UNCHANGED BY THE CALCULATION.
C
C     F(JW)            SCRATCH VECTOR
C
C     V(JW)            SCRATCH VECTOR
C
C *** OUTPUT:
C
C     U(IW*(IW+1)/2)   OUTPUT VECTOR STORED U-D COVARIANCE FACTORS
C                      IF ( IW .GT. JW ) THEN ( SEE REF. #1) IT IS
C                       ASSUMED THAT COLUMNS 1 UP TO IW-JW (OF U)
C                      ARE LEFT UNALTERED AND THE RESULT OF THE
C                      WGS PROCESS IS (VECTOR) STORED IN COLUMNS
C                      IW-JW+1,...,IW.
C
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              19 SEPTEMBER  1981
C
C              REVISED   10 OCTOBER  1981  VERSION 1.00.00
C              REVISED   05 DECEMBER 1982  VERSION 1.01.00
C              REVISED   20 MAY      1985  VERSION 2.00.00
C              REVISED   01 JANUARY  1992  VERSION 2.01.00 (RJM)
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE WGSG   IS A PART OF THE FEA. INC             *
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
  
      DOUBLE PRECISION SUM,DINV,Z, EPS
      DIMENSION W(IMAXW,JW),DW(JW),F(JW),V(JW),U(*)
C
      DATA Z /0.0D0/
C      DATA EPS/0.6D-37/ !suitable for single precision
      DATA EPS/1.0D-307/

C
C     EPS --- IS BOUND ON MACHINE ZERO
C     EPS =   0.6D-37 FOR VAX TYPE MACHINES
C
C     CALCULATE SOME CONSTANTS
C
      IWP1 = IW + 1
      JJ   = IW * IWP1/2
C
C     IF (IW = 1) SKIP MOST OF THE COMPUTATION
C
      IF (IW .EQ. 1) GO TO 200
C                                  MORE CONSTANTS
      L = IW - JW
      LP1 = L + 1
      IF (L .LE. 0) LP1=2
      LBACK = IW + LP1
      DO 100 LROW=LP1,IW
         J     = LBACK-LROW
         JM1   = J-1 
         JJNXT = JJ-J
         SUM=Z
         DO 10 K=1,JW
            V(K)=W(J,K)
            F(K)=V(K)*DW(K)
            SUM =SUM+V(K)*F(K)
  10     CONTINUE
C
C     SUM IS NEW D ELEMENT(J)
C
         DINV = SUM
         IF (DINV .GT. EPS) GO TO 30
C
C     WHEN D-ELEMENT IS EFFECTIVELY ZERO THE CORRESPONDING
C     COLUMN OF U ALSO SET TO ZERO
C
C           Also need to zero out the diagonal term of the UD array (RJM).
C
C           DO 20 K=1,JM1  
            DO 20 K=1,J
               KJ = JJNXT + K
               U(KJ) = Z
  20        CONTINUE
            GO TO 70
C
  30     CONTINUE
         U(JJ)= SUM
         DO 60 K=1, JM1
            SUM=Z
            DO 40 I=1,JW
               SUM = SUM + W(K,I) * F(I)
  40        CONTINUE
            SUM = SUM/DINV
C
C FOR NUMERICAL REASONS DIVIDE HERE IS BETTER THAN USING A
C RECIPROCAL
C
            DO 50 I=1,JW
               W(K,I) = W(K,I) - SUM * V(I)
  50        CONTINUE
C
            ITMP = JJNXT + K
            U(ITMP) = SUM
  60     CONTINUE
C
  70     CONTINUE
         JJ = JJNXT
 100  CONTINUE
C
C     NOW HANDLE FIRST ROW ( IF APPLICABLE )
C
      IF ( L .GT. 0) RETURN
 200  CONTINUE
C
      SUM = Z
      DO 210 K=1,JW
         SUM = SUM + DW(K) * W(1,K)*W(1,K)
 210  CONTINUE
C
      U(JJ) = SUM
C
      RETURN
      END
