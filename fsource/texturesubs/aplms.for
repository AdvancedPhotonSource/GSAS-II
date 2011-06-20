      SUBROUTINE APLMS(L,M,S,AP)

!PURPOSE: Compute A'(l,m,s)

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:             

      INTEGER*4     L,M,S               !Order & subindices
      REAL*4        AP                  !Output value

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

      REAL*4        LNORM,ALM0S,A1,A2   

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

!CODE:                     

      LNORM = SQRT((2.0*L+1.0)/2.0)
      CALL QLMN(L,M,S,A1)
      CALL QLMN(L,0,S,A2)
      ALM0S = A1*A2
      IF ( MOD(ABS(M),2).EQ.0 ) THEN
        IF ( S.EQ.0 ) THEN
          AP = ALM0S*LNORM
        ELSE
          AP = 2.0*ALM0S*LNORM
        END IF
        IF ( MOD(ABS(M),4).EQ.2 ) AP = -AP
      ELSE
        IF ( S.EQ.0 ) THEN
          AP = 0.0
        ELSE
          AP = 2.0*ALM0S*LNORM
        END IF      
        IF ( MOD(ABS(M),4).EQ.1 ) AP = -AP
      END IF
      IF ( ABS(AP).LT.1.0E-6 ) AP = 0.0
      RETURN                                
      END                   

