      SUBROUTINE SGDCMP (UIN,UOUT,N,WNTSIG,SIGMA,NAMES,NAMFLG,
     1                   TITLE,CODE,PRTNAM,RNUM,IOUNIT)
C
C     TO COMPUTE AND DISPLAY THE U-D SIGMA DECOMPOSITION ASSOCIATED
C     WITH THE COVARIANCE MATRIX UDU**T
C
C *** INPUT:
C
C     UIN(N*(N+1)/2)  INPUT U-D VECTOR STORED UPPER TRIANGULAR
C                     ARRAY WITH  D  ENTRIES STORED ALONG THE
C                     DIAGONAL  OF UIN.
C
C
C     N               DIMENSION OF THE ASSOCIATED ESTIMATE, N > 1.
C
C     WNTSIG          LOGICAL FLAG TO INCLUDE STANDARD DEVIATION
C                     PRINTED OUTPUT.
C
C     NAMES(N)        CHARACTER*(*) VECTOR OF VARIABLE NAMES.
C
C     NAMFLG          LOGICAL FLAG.
C                     IF .TRUE.  NAMES WILL BE USED ON OUTPUT,
C                       NAMES(N) WILL BE ASSUMED TO HAVE
C                       ALPHANUMERIC DATA.
C                     IF .FALSE. NUMBERS WILL BE USED IN PLACE OF
C                                LABELS.
C
C     TITLE           CHARACTER*(*) OUTPUT TITLE.
C
C
C     CODE            (CHARACTER*(*)) IF THE LETTER "P" APPEARS
C                     THEN A PAGE EJECT WILL PRECEDE PRINTING,
C                     IF THE LETTER "C" APPEARS CONDENSED PRINT
C                     WILL BE PRODUCED (I.E. LESS PRECISION, MORE
C                     COLUMNS PER ROW) THE DEFAULT (' ') CAUSES
C                     PRINTING TO BEGIN AT THE NEXT LINE AND
C                     DISPLAYS 6 COLUMNS WITH  8 DIGITS OF
C                     PRECISION. CONDENSED MODE CONSISTS OF
C                     11 COLUMNS OF 4 DIGIT PRECISION.
C
C     PRTNAM          LABEL FOR RNUM (SIMPLY PASSED TO
C                     RMPRNT AND TMPRNT)
C
C     RNUM            REAL NUMBER   (SIMPLY PASSED
C                     TO RMPRNT AND TMPRNT)
C
C     IOUNIT          UNIT NUMBER OF OUTPUT DEVICE.
C                     WHEN IOUNIT .EQ. 0 THEN OUTPUT
C                     IS SUPPRESSED.
C
C
C *** OUTPUT:
C
C     UOUT(N*(N+1)/2) OUTPUT U-D VECTOR STORED UPPER TRIANGULAR
C                     SIGMA DECOMPOSITION. UOUT MAY OVERWRITE UIN.
C                     NOTE: DIAGONAL ZEROS CAUSE THAT COLUMN TO
C                           BECOME ZERO FILLED.
C
C     SIGMA(N)        OUTPUT VECTOR OF ESTIMATE ERROR STANDARD
C                     DEVIATIONS
C
C
C     SIGMA(I)=SQRT(D(I)+D(I+1)*U(I,I+1)**2+...
C              +D(N)*U(I,N)**2)
C
C     UOUT(I,J)=U(I,J)*SQRT(D(J))/SIGMA(I)
C
C     UOUT(I,I)=SQRT(D(I))/SIGMA(I)
C
C     (IF SIGMA(I)=0., THE ASSOCIATED ROW OF UOUT IS SET TO ZERO
C
C     THE ROW SUM OF SQUARES OF UOUT SUMS TO UNITY, AND UOUT(1,I)
C     IS THE FRACTIONAL PART OF OUTPUT SIGMA(1) THAT IS ASSOCIATED
C     WITH THE MODELING OF VARIABLES 1 UP TO AND INCLUDING I.
C     UOUT(I,J), IS THE FRACTIONAL PART OF UOUT SIGMA(I)
C     ATTRIBUTABLE TO ESTIMATE STATE(J). THE UOUT ARRAY, DISPLAYED
C     IN TRIANGULAR FORM, IS A KIND OF CORRELATION MATRIX.
C
C  ***     REQUIRED SUBROUTINES:
C
C                 RMPRNT   <ESL>
C                 TMPRNT   <ESL>
C
C                 SQRT   <FORLIB>   IABS <FORLIB>
C                 DSQRT  <FORLIB>   MOD  <FORLIB>
C
C
C         COGNIZANT PERSONS:
C
C              DR. GERALD J. BIERMAN / KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              JUNE 2, 1981
C
C              REVISED 19 OCTOBER  1981
C              REVISED 02 FEBRUARY 1982
C              REVISED 02 AUGUST   1982   VERSION 0.01.00
C              REVISED 20 AUGUST   1982   VERSION 0.01.01
C                                         VERSION 1.00.00
C              REVISED 30  JUNE    1984   VERSION 1.01.00
C              REVISED 12  JULY    1984   VERSION 1.02.00
C              REVISED 16   MAY    1985   VERSION 2.00.00
C              REVISED 15  JANUARY 1986   VERSION 2.01.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE SGDCMP IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *            COPYRIGHT 1982, 1983, 1984, 1985                    *
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
  
      LOGICAL WNTSIG, NAMFLG
      DOUBLE PRECISION SUM
      DIMENSION UIN(*),UOUT(*),SIGMA(N)
C
      CHARACTER*(*) TITLE,NAMES(N)
      CHARACTER*(*) PRTNAM,CODE
C
C                        PRINT TITLE
C
      IF (IOUNIT .NE. 0) THEN
         IF     ((INDEX(CODE,'P').NE.0)
     1     .OR. ((INDEX(CODE,'P').NE.0))) THEN
                WRITE (IOUNIT,'(1H1)')
         ELSE
             WRITE (IOUNIT,*)
             WRITE (IOUNIT,*)
          END IF
         WRITE(IOUNIT,*) TITLE
      ENDIF
      NM1=N-1
      JJ=0
      DO 10 J=1,N
         JJ=JJ+J
C
C     NOTE:  HERE SIGMAS ARE USED FOR TEMPORARY STORAGE
C
         SIGMA(J)=UIN(JJ)
         UOUT(JJ)=SQRT(SIGMA(J))
  10  CONTINUE
      JJ=0
      IF (NM1 .GT. 0) GO TO 20
         KK = 1
         GO TO 110
  20  CONTINUE
      DO 100 J=1,NM1
         JJ=JJ+J
         SUM=SIGMA(J)
C
C        SUM= SCALED DIAGONAL OF UD
C
         JKS=JJ+J
C
C        JKS = STARTING INDEX (J,K)=(J,J+1)
C
         KK=JJ
         JP1=J+1
         DO 30 K=JP1,N
            SUM=SUM+SIGMA(K)*UIN(JKS)**2
            KK=KK+K
            UOUT(JKS)=UOUT(KK)*UIN(JKS)
            JKS=JKS+K
  30     CONTINUE
C
         SUM=DSQRT(SUM)
         SIGMA(J)=SUM
         JKS=JJ
C
         IF (SUM .LE. 0.0) THEN
            DO 40 K=J,N
               UOUT(JKS)=0.0
               JKS=JKS+K
  40        CONTINUE
            GO TO 100
         ENDIF
         DO 50 K=J,N
            UOUT(JKS)=UOUT(JKS)/SUM
            JKS=JKS+K
  50     CONTINUE
 100  CONTINUE
 110  CONTINUE
      SIGMA(N)= UOUT(KK)
      IF (UOUT(KK) .NE. 0.0) UOUT(KK) = 1.0
C
C   PRINT SIGMA'S DEPENDING ON FLAG
C
      IF (WNTSIG) CALL RMPRNT(SIGMA,1,1,N,NAMES,NAMES,
     1 ' STANDARD DEVIATIONS  ',CODE,.FALSE.,NAMFLG,PRTNAM,RNUM,IOUNIT)
C
      CALL TMPRNT(UOUT,N,NAMES,' FRACTIONAL SIGMA DECOMPOSITION '
     1     ,CODE,NAMFLG,PRTNAM,RNUM,IOUNIT)
      RETURN
      END
