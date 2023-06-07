C-------------------------------------------------
      SUBROUTINE SWPUCL(U,N,NAMES,ICOL,JCOL,SCRATCH)
C
C     SWAP U MATRIX COLUMNS
C
C   ** INPUT
C
C     U(N*(N+1)/2)     UD FACTORS, VECTOR STORED.
C
C     N                STATE SIZE.
C
C     NAMES(N)         STATE NAMES, IF NAMES(1) = '*', IT IS ASSUMED THAT
C                      THERE ARE NO STATE NAMES.
C
C     ICOL,JCOL        COLUMNS TO BE SWAPPED.
C
C     SCRATCH(N,4)     COLUMNS 1 AND 2 OF SCRATCH MUST BE SET TO ZERO
C                      BY THE CALLER. ALTHOUGH IT WOULD BE A BIT MORE
C                      CONVIENIENT IF SCRATCH WERE TO BE ZEROED
C                      INTERNALLY; BUT SWPUCL IS TYPICALLY CALLED
C                      REPEATEDLY, IN A TIGHT LOOP. HAVING SCRATCH
C                      ZEROED EXTERNALLY SAVES CPU CYCLES.
C ** OUTPUT **
C
C      U(N*(N+1)/2)   PERMUTED U ARRAY.
C
C ******************************************************************
C *                                                                *
C *        SUBROUTINE SWPUCL IS A PART OF THE FEA. INC             *
C *              ESTIMATION SUBROUTINE LIBRARY                     *
C *                                                                *
C *                         (ESL)                                  *
C *                                                                *
C *                     COPYRIGHT 1987                             *
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

C
      DIMENSION     U(*)
      INTEGER       N,ICOL,JCOL
      DIMENSION     SCRATCH(N,4)
      CHARACTER*(*) NAMES(*)
C                              LOCAL VARS
      CHARACTER*17  NAMTMP
      INTEGER       IERR
C
      IF (ICOL .EQ. JCOL) RETURN
C                                NOTHING TO DO
      ZERO=0.0
      SCRATCH(ICOL,1) = 1.0
      SCRATCH(ICOL,2) = 1.0
      SCRATCH(JCOL,1) =-1.0
      SCRATCH(JCOL,2) =-1.0
      CALL UDGEN(U,N,SCRATCH(1,2),ZERO,SCRATCH(1,1),N,
     1           SCRATCH(1,3),SCRATCH(1,4),IERR)
      SCRATCH(ICOL,1) = 0.0
      SCRATCH(ICOL,2) = 0.0
      SCRATCH(JCOL,1) = 0.0
      SCRATCH(JCOL,2) = 0.0
      IF (IERR .GT. 0) THEN
         WRITE(6,*)' ERROR IN SWPUCL, UDGEN  ',IERR,
     1             '-TH DIAGONAL IS NEGATIVE'
         STOP
      ENDIF
      IF (NAMES(1) .NE. '*') THEN
         NAMTMP      = NAMES(ICOL)
         NAMES(ICOL) = NAMES(JCOL)
         NAMES(JCOL) = NAMTMP
      ENDIF
      RETURN
      END
