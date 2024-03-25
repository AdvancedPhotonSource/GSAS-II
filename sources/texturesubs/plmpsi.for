      SUBROUTINE PLMPSI(L,M,PSI,P,DPDPS)

!PURPOSE: Compute P(l,m,psi)

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:

      INTEGER*4     L,M                 !Order & index
      REAL*4        PSI                 !Angle (in deg)
      REAL*4        P                   !Value returned

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:
                                                 
      INTEGER*4     S                   
      REAL*4        APR,RS              
                                                 
!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

!CODE:                 
         
      P = 0.0
      DPDPS = 0.0
      IF ( MOD(ABS(M),2).EQ.0 ) THEN
        DO S=0,L,2
          CALL APLMS(L,M,S,APR)
          RS = S
          P = P+APR*COSD(RS*PSI)
          DPDPS = DPDPS-RS*APR*SIND(RS*PSI)
        END DO
      ELSE                                            
        DO S=2,L,2        
          CALL APLMS(L,M,S,APR)      
          RS = S
          P = P+APR*SIND(RS*PSI)        
          DPDPS = DPDPS+RS*APR*COSD(RS*PSI)
        END DO      
      END IF
      RETURN
      END

