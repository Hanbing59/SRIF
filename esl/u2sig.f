      SUBROUTINE U2SIG(U,N,TEXT,CODE,SF,NAMES,PRNAM,
     1                 JSTRT,NSIG,SIG,PRTNAM,RNUM,IOUNIT)
C
C     COMPUTE STANDARD DEVIATIONS (SIGMAS) FROM
C     U-D COVARIANCE FACTORS
C
C *** INPUTS:
C
C     U(N*(N+1)/2)  INPUT VECTOR STORED ARRAY CONTAINING THE U-D
C                   FACTORS. THE D (DIAGONAL) ELEMENTS ARE STORED
C                   ON THE DIAGONAL
C
C     N             STATE DIMENSION  (N .GE. 1)
C
C     TEXT( )       ARRAY OF ASCII CHARACTERS TO BE PRINTED
C                   PRECEDING THE VECTOR OF SIGMAS
C
C     CODE          (CHARACTER*(*)) IF THE LETTER "P" APPEARS
C                   THEN A PAGE EJECT WILL PRECEDE PRINTING,
C                   IF THE LETTER "C" APPEARS CONDENSED PRINT
C                   WILL BE PRODUCED (I.E. LESS PRECISION, MORE
C                   COLUMNS PER ROW) THE DEFAULT (' ') CAUSES
C                   PRINTING TO BEGIN AT THE NEXT LINE AND DISPLAYS
C                   6 COLUMNS WITH  8 DIGITS OF PRECISION.
C                   CONDENSED MODE CONSISTS OF 11 COLUMNS OF 4 DIGIT
C                   PRECISION.
C
C     NAMES(N)      CHARACTER*(*) VECTOR OF NAME LABELS ASSOCIATED
C                   WITH THE STATE VECTOR COMPONENTS THAT ARE
C                   TO BE PRINTED.
C
C
C     PRNAM         LOGICAL FLAG. IF .TRUE. NAMES LABELS
C                   ARE USED; IF .FALSE., DEFAULT OF
C                   INTEGER INDICES ARE USED.
C
C     SF(N)         VECTOR OF SCALE FACTORS
C                      IF SF(1) .LT. 0.0, SCALE FACTORS ARE NOT USED
C                      I.E. THEY ARE IMPLICITLY TREATED AS UNITY
C
C     JSTRT         START LOCATION OF SIGMAS THAT ARE
C                   TO BE COMPUTED
C
C     NSIG          NUMBER OF SIGMAS
C
C     PRTNAM        CHARACTER*(*) AN ALPHANUMERIC VARIABLE
C                   OR TITLE, PRINTED AFTER TEXT(); TO BE USED
C                   AS A LABEL FOR THE SCALAR RNUM.
C                   IF .EQ. ' ' THEN PRTNAM AND RNUM ARE
C                   NOT PRINTED. SEE SUBROUTINE RMPRNT FOR DETAILS
C
C     RNUM          A (REAL) NUMBER TO BE PRINTED FOLLOWING
C                   TEXT() AND PRTNAM, DEPENDING ON THE VALUE
C                   OF PRTNAM.
C
C     IOUNIT        UNIT NUMBER OF OUTPUT DEVICE.
C                   IF .EQ. 0 THEN NO OUTPUT IS PRODUCED
C
C *** OUTPUT:
C
C     SIG(NSIG)     VECTOR OF OUTPUT STANDARD DEVIATIONS.
C
C    REQUIRED SUBROUTINES:
C
C             RMPRNT   <ESL>       SQRT <FORLIB>
C
C              DR. GERALD J. BIERMAN/KEITH H. BIERMAN
C              FACTORIZED ESTIMATION APPLICATIONS
C              DECEMBER 7, 1981
C
C              REVISED   9 AUGUST 1982    VERSION  0.01.00
C              REVISED  13 AUGUST 1982    VERSION  0.01.01
C              REVISED  20 AUGUST 1982    VERSION  0.01.02
C                                         VERSION  1.00.00
C              REVISED  23 APRIL  1984    VERSION  2.00.00
C              REVISED  16   MAY  1985    VERSION  2.01.00
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE U2SIG  IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *              COPYRIGHT 1982, 1983, 1984, 1985                  *
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

      DIMENSION U(*), SIG(NSIG), SF(N)
      LOGICAL PRNAM
      CHARACTER*(*) NAMES(N)
      CHARACTER*(*) PRTNAM
      CHARACTER*(*) TEXT,CODE
C
      JSM1  = JSTRT-1
      JSP1  = JSTRT+1
      JJ    = JSTRT*JSP1/2
      K     = 1
      SIG(1)=U(JJ)
      IF (JSP1  .GT. N) GO TO 20
      DO 15 J=JSP1,N
         K   = K+1
         JJL = JJ+JSM1
         JJ  = JJ+J
         S   = U(JJ)
         KM1 = NSIG
         IF (K .GT. NSIG) GO TO 5
            SIG(K) = S
            KM1    = K-1
   5     CONTINUE
         DO 10 I=1,KM1
            SIG(I)=SIG(I)+S*U(JJL+I)**2
  10     CONTINUE
  15  CONTINUE
  20  CONTINUE
C
C           WE NOW HAVE VARIANCES
C
      DO 30 J=1,NSIG
         SIG(J) = SQRT(SIG(J))
  30  CONTINUE
      IF (SF(1) .LT. 0.0) GO TO 50
         DO 40 J=1, NSIG
            SIG(J) = SIG(J) * SF(J)
  40     CONTINUE
  50  CONTINUE
      IF (IOUNIT .NE. 0) THEN
        CALL RMPRNT(SIG,1,1,NSIG,NAMES,NAMES,TEXT,CODE,.FALSE.,
     2               PRNAM,PRTNAM,RNUM,IOUNIT)
      ENDIF
      RETURN
      END
