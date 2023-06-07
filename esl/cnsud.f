
      SUBROUTINE CNSUD (U,N,BX,NX,EM,Q,IP,VM,F,B,G,IOUNIT)
C
C   COLORED NOISE SMOOTH U-D COVARIANCE FACTOR UPDATE TO ACCOUNT
C   FOR POSSIBLE OFF DIAGONAL PROCESS NOISE COMPONENTS.
C   DYNAMIC MODEL:
C
C      !X   !     !I   0   0! !X !     !BX!
C      ! I+1!     !         ! ! I!     !  !
C      !P   !  =  !0  EM   0! !P !  +  ! 1! W
C      ! I+1!     !         ! ! I!     !  !  I
C      !Y   !     !0   0   I! !Y !     ! 0!
C
C   REFERENCE:
C     G. BIERMAN - MODIFICATION OF THE BIERMAN COLORED NOISE
C                  SMOOTHER TO ACCOUNT FOR A NON-DIAGONAL
C                  PROCESS NOISE COVARIANCE  B&A 860429
C
C   INPUT:
C
C      U(N(N+1)/2)  VECTOR STORED U-D SMOOTH COVARIANCE
C                   FOR ESTIMATE AT STAGE I+1
C
C      N            FILTER MODEL DIMENSION
C
C      BX(NX)       STATE DEPENDENT PROCESS NOISE MULTIPLIER
C
C      NX           ONLY THE FIRST NX VARIABLES ARE ASSUMED
C                   TO BE DIRECTLY AFFECTED BY THE NOISE COMPONENT.
C                   IF NX=0 THE MODEL IS SIMPLY COLORED NOISE
C                   AND BX IS IGNORED. (I.E. IT IS NOT USED)
C
C      EM           COLORED NOISE MULTIPLIER
C
C      Q            PROCESS NOISE VARIANCE
C
C      IP           LOCATION OF THE COLORED NOISE COMPONENT,
C                   IP.GE.1  (ELSE RETURN)
C                   IP.GT.NX (ELSE PROGRAM TERMINATES, AS THIS IS A
C                   VIOLATION OF THE MATH MODEL)
C
C      VM(N)        SMOOTH GAIN INPUT (SMOOTH GAINS ARE
C                   MANUFACTURED FROM THIS)
C
C      F(N),B(N),G(N) SCRATCH VECTORS
C
C      IOUNIT       UNIT NUMBER FOR FAILURE MODE
C
C   OUTPUT:
C
C      U(N(N+1)/2)  UPDATED U-D FACTOR SMOOTH COVARIANCE FACTORS
C                   CORRESPONDING TO THE STATE AT TIME 'I'
C
C    SUBROUTINES NEEDED:
C
C        UDOPT
C        UCON
C        RNK1
C
C    COGNIZANT PERSONS:
C
C           DR.GERALD J. BIERMAN/KEITH H. BIERMAN
C           FACTORIZED ESTIMATION APPLICATIONS
C           MAY 14, 1986
C
C ************************************************************
C *                                                          *
C *        SUBROUTINE CNSUD IS A PART OF THE FEA. INC        *
C *              ESTIMATION SUBROUTINE LIBRARY               *
C *                                                          *
C *                          (ESL)                           *
C *                                                          *
C *                     COPYRIGHT 1986                       *
C *           FACTORIZED ESTIMATION APPLICATIONS INC.        *
C *   This library is licensed under the Creative Commons    *
C *   Attribution 3.0 Unported License.                      *
C *   To view a copy of this license,                        *
C *   visit http://creativecommons.org/licenses/by/3.0/      *
C *   or send a letter to Creative Commons, 171 Second Street*
C *   Suite 300, San Francisco, California, 94105, USA.      *
C *                                                          *
C ************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(*),BX(*),VM(N),F(N),B(N),G(N)
      DOUBLE PRECISION GAMMA,EPS
      CHARACTER*80 ERRMSG
C
C      DATA EPS/1.D-37/ ! suitable for single precision
      DATA EPS /1.0D-307/
      DATA ERRMSG/'(10X,20(''*''),''CNSUD ERROR'',20(''*''))'/
C      
      IF (IP .LT. 1 ) RETURN
      IF (IP .LE. NX) THEN
          WRITE(IOUNIT,ERRMSG) 
          WRITE(IOUNIT,*) ' IP .LE. NX, IP = ',IP,' NX = ',NX
          WRITE(IOUNIT,*) ' THIS IS A VIOLATION OF THE MATH MODEL'
          WRITE(IOUNIT,*) ' EXECUTION TERMINATED'
          STOP
      ENDIF
      TMP   = Q*VM(IP)
      GAMMA = TMP+EM**2
      IF (GAMMA .LE. EPS) THEN
         WRITE(IOUNIT,ERRMSG) 
         WRITE (IOUNIT,'(10X,'' ERROR INPUT TO SMOOTH SUBROUTINE: Q ='',
     1   1PD14.7,'' EM ='',D14.7,''  SMOOTH GAIN('',I4,'') ='',
     2   D14.7/40X,10(''*''),'' SMOOTH STEP SKIPPED '',10(''*''))')
     3        Q,EM,IP,F(IP)
         RETURN
      ENDIF
      IF (VM(IP) .EQ. 0.0d0) THEN
         KL = IP*(IP-1)/2
         DO 5 J = 1, IP
            U(KL+J) = 0.0d0
    5    CONTINUE
         RETURN
      END IF
      IF (NX .EQ. 0) THEN
         DO 10 J=1,N
             F(J) = Q*VM(J)
   10    CONTINUE
         F(IP) = GAMMA-EM
         S     = 1.0d0/GAMMA
C
C   SET B**T = (F**T)*U
C
C
         B(1) = F(1)
         IF (N .GT. 1) THEN
            JJ = 1
            DO 30 K=2,N
               SUM = F(K)
               DO 20 I=1,K-1
                  SUM = SUM+F(I)*U(JJ+I)
   20          CONTINUE
               B(K) = SUM
               JJ   = JJ+K
   30       CONTINUE
            IF (IP .LT. N) THEN
               KP1 = IP+1
               KL  = IP*KP1/2+IP
               DO 40 K  = KP1,N
                  U(KL) = U(KL)-S*B(K)
                  KL    = KL+K
   40          CONTINUE
            ENDIF
         ENDIF
C
C   UPDATE IP BY IP BLOCK OF U
C   FIRST SET G = D*B
C
         JJ = 0
         DO 50 I=1,IP
            JJ   = JJ+I
            G(I) = U(JJ)*B(I)
   50    CONTINUE
         TMP   = Q*GAMMA
         CALL UDOPT (U,IP,TMP,B,G,ALPHA)
         IF(ALPHA.EQ.0.0D0) RETURN
         G(IP) = G(IP)-ALPHA*S
         C     = 1.0d0/ALPHA
         CALL RNK1 (U,IP,C,G,IERR)
         RETURN
      END IF
C                         NX .NE. 0
      TMP    = Q*VM(IP)
      GAMMA  = TMP+EM**2
      IF (GAMMA .LT. EPS) RETURN
      S      = TMP/GAMMA
      SUM    = 0.0d0
      DO 60 I=1,NX
         SUM = SUM+BX(I)*VM(I)
   60 CONTINUE
      THETA = EM + Q*SUM
      D     = EM*THETA + TMP
      IF (D .GT. 0.0d0) THEN
         C  = Q/D
      ELSE
        WRITE(IOUNIT,ERRMSG) 
        WRITE (IOUNIT,'(10X,'' INPUT ERROR: THETA ='',1PD14.7,
     1  '' D ='',D14.7,'' EM ='',D14.7/ 40X,10(''*''),
     2  '' SKIPPED SMOOTH UPDATE'',10(''*''))') THETA,D,EM
         RETURN
      ENDIF
C
C   FORM G,B (OBSERVATION)
C
      DO 70 K=1,N
         B(K) = 0.0d0
   70 CONTINUE
      DO 75 K=1,NX
         G(K) = S*BX(K)
   75 CONTINUE
      IF ((NX+2) .LT. IP) THEN
         DO 80 K = NX+1, IP-1
            G(K) = 0.0d0
   80    CONTINUE
      ENDIF
      G(IP) = (TMP-EM*(1.0d0-EM))/GAMMA
      B(IP) = 1.0d0
C
C   PERFORM FIRST UPDATE
C
      CALL UCON (U,N,B,0.0d0,G,IP,F,B,IERR)
      IF (C .EQ. 0.0d0) RETURN
C
C   FORM SECOND G,B TERMS
C
      DO 85 K=1,IP-1
         G(K) = 0.0d0
   85 CONTINUE
      S = EM*C
      DO 90 K=1,NX
         G(K) = S*BX(K)
   90 CONTINUE
      G(IP)  = C
      TMP    = VM(IP)
      VM(IP) = 0.0d0
      C      = 1.0d0/C
      CALL UCON (U,N,VM,C,G,IP,F,B,IERR)
C
C   RESTORE VM TERM
C
      VM(IP) = TMP
      RETURN
      END
