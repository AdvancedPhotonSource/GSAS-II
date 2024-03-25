      FUNCTION FACTLN(N)

!PURPOSE: Compute ln(factorial) - adapted from Numerical Recipes

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:

      INTEGER*4     N                   !Order of factorial
      REAL*8        FACTLN              !ln(N!) returned

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:
                                                                
      REAL*8        A(100)              
      INTEGER*4     IA(100)             
                                                                
!FUNCTION DEFINITIONS:                          

      REAL*8        DGAMMLN             

!DATA STATEMENTS:                     
         
      DATA          A/100*-1.0D0/,IA/100*0/

!CODE:              

      IF ( N.LE.99 ) THEN
        IF ( IA(N+1).EQ.0 ) THEN
          IA(N+1) = 1
          A(N+1) = DGAMMLN(N+1.0)
        END IF
        FACTLN = A(N+1)
      ELSE
        FACTLN = DGAMMLN(N+1.0)
      END IF
      RETURN
      END      

