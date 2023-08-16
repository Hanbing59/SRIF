    !! Purpose    : COMPUTE U:=PHI*U WHERE PHI IS UNIT UPPER TRIANGULAR AND SPARSE.
    !!              ONLY THE NON-ZERO (AND NON-DIAGONAL) PHI ELEMENTS ARE
    !!              EXPLICITLY EMPLOYED. U IS VECTOR STORED, UPPER TRIANGULAR WITH
    !!              IMPLICITLY DEFINED UNIT DIAGONAL ELEMENTS.
    !!
    !! Parameters :
    !!              PHIEL --- VALUES OF THE NONZERO ELEMENTS OF PHI, NOT INCLUDING
    !!                        THE UNIT DIAGONAL ELEMENTS. ELEMENTS ARE STORED ROW-WISE
    !!                        VIZ. FIRST THE NONZERO ELEMENTS IN ROW #1 (NOT
    !!                        INCLUDING PHI(1,1)). THESE ARE FOLLOWED BY THE NONZERO
    !!                        ELEMENTS IN ROW 2 (IF THERE ARE ANY) AND NOT
    !!                        INCLUDING PHI(2,2). THE TOTAL NUMBER OF ELEMENTS, NPHI,
    !!                        IS NOT EXPLICTLY USED.
    !!               INDX --- INDX(1,1)=ROW NUMBER OF the FIRST ROW WITH NONZERO,
    !!                        NON-DIAGONAL ENTRIES.
    !!                        INDX(1,2)=ROW NUMBER OF THE SECOND ROW WITH NONZERO,
    !!                        NON-DIAGONAL ENTRIES, ETC.
    !!                        ROW ORDERING MUST BE MONOTONIC
    !!                        INDX(2,1)=NUMBER OF NONZERO ELEMENTS (NOT INCLUDING
    !!                        THE DIAGONAL) IN ROW INDX(1,1)
    !!                        INDX(2,2)=NUMBER OF NONZERO ELEMENTS (NOT INCLUDING
    !!                        THE DIAGONAL) IN ROW INDX(1,2), ETC.
    !!                        NOTE THAT THE ROW INDICES MUST BE INCREASING IN THIS
    !!                        CODE FORMULATION. IT IS NOT REQUIRED THAT THE COLUMN
    !!                        INDICES BE MONOTONIC, SEE EXAMPLE.
    !!
    !!                        EXAMPLE TO ILLUSTRATE PHI ELEMENT LABELING
    !!
    !!                                       1.   0.  0.    5.0   0.     2.
    !!
    !!                                       0.   1.  3.0   0.    4.     5.
    !!
    !!                                       0.   0.  1.0   0.    0.     0.
    !!                        PHI =
    !!                                       0.   0.  0.    1.    0.     1.5
    !!
    !!                                       0.   0.  0.    0.    1.    -2.
    !!
    !!                                       0.   0.  0.    0.    0.     1.
    !!                        PHIEL(7) = (5. , 2. , 5. , 3. , 4. , 1.5 , -2.)
    !!                        INDX(1,1)= 1 , INDX(1,2)= 2 ,  INDX(1,3)= 4 ,  INDX(1,4)= 5
    !!                        INDX(2,1)= 2 , INDX(2,2)= 3 ,  INDX(2,3)= 1 ,  INDX(2,4)= 1
    !!                        JCOL(7)  = (4,6,6,3,5,6,6)
    !!                        NROW     = 4

    !!               NROW --- THE NUMBER OF ROWS WITH NONZERO ENTRIES, NOT INCLUDING
    !!                        THE DIAGONALS.
    !!               JCOL --- COLUMN LOCATIONS OF THE ELEMENTS IN THE VECTOR PHIEL.
    !!                  U --- UPPER TRIANGULAR U MATRIX, VECTOR STORED WITH IMPLICIT
    !!                        UNIT DIAGONAL ENTRIES. THE OUTPUT U OVERWRITES THE INPUT U
    !!                  N --- COLUMN DIMENSION OF U
    !!
    SUBROUTINE STPHIU(PHIEL,INDX,NROW,JCOL,U,N)

    IMPLICIT none

    real*8    PHIEL(*)
    integer*4 INDX(2,NROW)
    integer*4 NROW
    integer*4 JCOL(*)
    real*8    U(*)
    integer*4 N
    !local
    integer*4 NM1,JFLG,KOL,I,KOLDIAG,IK,K,KDIAG,IR,JK,ISTRT,NELEM,JC,L
    real*8    A

    IF(N.EQ.1 .OR. NROW.EQ.0) RETURN

    NM1 = N-1
    ! Is this a case where NELEM is always 1, and all in same column?
    JFLG = 1
    IF(INDX(2,1).NE.1) JFLG = 0
    KOL = JCOL(1)
    DO I=2,NROW
        IF(INDX(2,I).NE.1 .OR. JCOL(I).NE.KOL) JFLG=0
    enddo
    IF(JFLG .EQ. 1) THEN
        ! no need for recursive definition of jstart; may operate down columns
        KOL = JCOL(1)
        KOLDIAG = KOL*(KOL-1)/2
        !CDIR$ IVDEP
        DO I = 1,NROW
            ! take care of first column of U, that uses diagonal unit U
            IK = INDX(1,I) + KOLDIAG
            U(IK) = U(IK) + PHIEL(I)
        enddo
        !CDIR$ IVDEP
        DO K = KOL,NM1
            KDIAG = K*(K+1)/2
            DO I=1,NROW
                IR = INDX(1,I)
                IK = IR + KDIAG
                JK = KOL + IK - IR
                U(IK) = U(IK) + PHIEL(I)*U(JK)
            enddo
        enddo
    ELSE
        ISTRT = 0
        DO I=1, NROW
            IR = INDX(1,I)
            NELEM = INDX(2,I)
            DO JC=1, NELEM
                L = ISTRT + JC
                A = PHIEL(L)
                KOL = JCOL(L)
                JK = KOL * (KOL-1)/2
                IK = JK+IR
                U(IK) = U(IK)+A
                IF(KOL.EQ.N) cycle
                DO K=KOL, NM1
                    IK = IR + K*(K+1)/2
                    JK = KOL + IK - IR
                    U(IK)=U(IK)+A*U(JK)
                enddo
            enddo
            ISTRT = ISTRT+NELEM
        enddo
    ENDIF

    RETURN
    END