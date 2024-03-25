      FUNCTION DGAMMLN(XX)

!PURPOSE: ln(gamma(xx)), taken from Numerical Recipes, p156-157

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!PSEUDOCODE:

!CALLING ARGUMENTS:

      REAL*8        DGAMMLN              !ln(gamma(xx)) with xx>0
      REAL*4        XX                  !argument must be >0

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

      REAL*8        COF(6),STP,HALF,ONE,FPF,X,TMP,SER 

!SUBROUTINES CALLED:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA      COF,STP/76.18009172947146D0,-86.50532032941677D0,
     1  24.01409824083091D0,-1.231739572450155D0,0.1208650973866D-2,
     1  -0.5395239384953D-5,2.5066282746310005D0/
      DATA      HALF,ONE,FPF/0.5D0,1.000000000190015D0,5.5D0/

!CODE:
                                            
      X = XX-ONE
      TMP = X+FPF
      TMP = (X+HALF)*DLOG(TMP)-TMP
      SER = ONE
      DO J=1,6
        X = X+ONE
        SER = SER+COF(J)/X
      END DO
      DGAMMLN = TMP+DLOG(STP*SER)
      RETURN
      END
