!! Purpose    : ADD IN MARKOV COLORED NOISE MODEL EFFECTS INTO THE SQUARE ROOT
!!              INFORMATION FILTER, AND TO GENERATE WEIGHTING COEFFICIENTS
!!              FOR SMOOTHING.
!!              IT IS ASSUMED THAT VARIABLES X(NPSTRT), X(NPSTRT+1),...,
!!              X(NPSTRT+NP-1) ARE COLORED NOISE AND THAT EACH COMPONENT SATISFIES
!!              A MODEL EQUATION OF THE FORM X(SUB)(J+1)=EM*X(SUB)(J)+W(SUB)(J).
!!              FOR SMOOTHING, COEFFICIENTS TO BE GENERATED SGSTAR(1).GT.0.0 .
!!              THE SIGNIFICANCE OF THE SMOOTHING COEFFICIENTS IS EXPLAINED IN
!!              THE FUNCTIONAL DESCRIPTION.
!! Parameters :
!!                 S --- INPUT SQUARE ROOT INFORMATION ARRAY. OUTPUT COLORED NOISE
!!                       ARRAY HOUSED HERE TOO. WHEN THERE IS SMOOTHING, NR ADDITIONAL
!!                       ROWS MUST BE INCLUDED IN S
!!              MAXS --- ROW DIMENSION OF S. WHEN THERE ARE SMOOTHING COMPUTATIONS
!!                       IT IS NECESSARY THAT MAXS.GE.IRS+NP BECAUSE THE BOTTOM NP
!!                       ROWS OF S HOUSE THE SMOOTHING INFORMATION
!!               IRS --- NUMBER OF ROWS OF S (.LE. NUMBER OF FILTER VARIABLES) (IRS.GE.1)
!!               JCS --- NUMBER OF COLUMNS OF S (EQUALS NUMBER OF FILTER VARIABLES +
!!                       POSSIBLY A RIGHT SIDE), WHICH CONTAINS THE DATA EQUATION
!!                       NORMALIZED ESTIMATE (JCS.GE.1)
!!            NPSTRT --- LOCATION OF THE FIRST COLORED NOISE VARIABLE (1.LE.NPSTRT.LE.JCS)
!!                NP --- NUMBER OF CONTIGUOUS COLORED NOISE VARIABLES (NP.GE.1)
!!                EM --- COLORED NOISE MAPPING COEFFICIENTS (GENERALLY OF EXPONENTIAL
!!                       FORM, EM=EXP(-DT/TAU))
!!                RW --- RECIPROCAL PROCESS NOISE STANDARD DEVIATIONS (MUST BE POSITIVE)
!!                ZW --- ZW=RW*W-ESTIMATE (PROCESS NOISE ESTIMATES ARE GENERALLY ZERO MEAN).
!!                       WHEN THE RIGHT HAND SIDE COLUMN IS OMITTED, ZW IS NOT USED.
!!                 V --- WORK VECTOR
!!            SGSTAR --- SGSTAR(1).GT.0.0, SMOOTHING COEFFICIENTS ARE GENERATED.
!!                       SGSTAR(1).LE.0.0, NO USE IS MADE OF SGSTAR OR OF ROWS IRS+1,...IRS+NP OF S.
!!                       SGSTAR(NP), VECTOR OF SMOOTHING COEFFICIENTS. WHEN SGSTAR(1).LE.0.0,
!!                                   SGSTAR IS NOT USED. INCIDENTALLY, WHEN SMOOTHING THE SGSTAR
!!                                   MULTIPLIERS ARE ALWAYS POSITIVE.

SUBROUTINE RCOLRD(S,MAXS,IRS,JCS,NPSTRT,NP,EM,RW,ZW,V,SGSTAR)

implicit none

real*8    S(MAXS,JCS)
integer*4 MAXS
integer*4 IRS,JCS
integer*4 NPSTRT,NP
real*8    EM(NP)
real*8    RW(NP)
real*8    ZW(NP)
real*8    V(IRS)
real*8    SGSTAR(NP)
!local
integer*4 NPCOL,K
integer*4 JCOLRD,IPJC,KOL
LOGICAL   ISMTH
real*8    ALPHA,SIGMA,BETA,GAMMA,TMP,ZERO,ONE
DATA ZERO/0.0D0/, ONE/1.0D0/

ISMTH = .FALSE.
IF(SGSTAR(1).LE.0.0) ISMTH=.TRUE.

! COL NO OF COLORED NOISE TERM TO BE OPERATED ON
NPCOL=NPSTRT
DO JCOLRD=1,NP
    TMP = -RW(JCOLRD)*EM(JCOLRD)
    ALPHA = TMP
    SIGMA = ZERO
    DO K=1,IRS
        ! FIRST IRS ELEMENTS OF HOUSEHOLDER TRANSFORMATION VECTOR
        V(K)=S(K,NPCOL)
        SIGMA=SIGMA+V(K)**2
    enddo
    ! IPJC  = I+JC
    IPJC = IRS + JCOLRD
    IF(SIGMA.le.ZERO .and. ISMTH) then
        DO KOL = 1, JCS
            S(IPJC,KOL) = ZERO
        enddo
        S(IPJC,NPCOL) = RW(JCOLRD)
        SGSTAR(JCOLRD) = SIGMA
        S(IPJC,JCS) = S(IPJC,JCS) + ZW(JCOLRD)
        NPCOL=NPCOL+1
        cycle
    endif
    SIGMA=DSQRT(SIGMA+ALPHA*ALPHA)
    ! LAST ELEMENT OF HOUSEHOLDER TRANSFORMATION VECTOR
    ALPHA=ALPHA-SIGMA
    ! HOUSEHOLDER=I+BETA*V*V**T
    BETA=ONE/(SIGMA*ALPHA)
    ! HOUSEHOLDER TRANSFORMATION DEFINED, NOW APPLY IT TO S,

    DO KOL=1,JCS
        IF(KOL.eq.NPCOL) then
            GAMMA=RW(JCOLRD)*ALPHA*BETA
            IF(ISMTH) S(IPJC,NPCOL)=TMP*GAMMA
            DO K=1,IRS
                S(K,NPCOL)=GAMMA*V(K)
            enddo
            cycle
        endif
        GAMMA=ZERO
        IF(KOL.EQ.JCS) GAMMA=ZW(JCOLRD)*ALPHA
        ! IF ZW ALWAYS ZERO OR THERE ARE NO ESTIMATES ONE
        ! CAN COMMENT OUT THE ABOVE IF STATEMENT
        DO K=1,IRS
            GAMMA=GAMMA+S(K,KOL)*V(K)
        enddo
        GAMMA= GAMMA*BETA
        DO K=1,IRS
            S(K,KOL)=S(K,KOL)+GAMMA*V(K)
        enddo
        IF(ISMTH) S(IPJC,KOL)=GAMMA*ALPHA
    ENDDO
    IF(ISMTH) then
        SGSTAR(JCOLRD) = SIGMA
        S(IPJC,JCS) = S(IPJC,JCS) + ZW(JCOLRD)
        ! IF ZW IS ALWAYS ZERO OR THERE ARE NO ESTIMATES INVOLVED,
        ! ONE CAN COMMENT OUT THE ABOVE STATEMENT
    endif
    NPCOL=NPCOL+1
enddo

RETURN
END