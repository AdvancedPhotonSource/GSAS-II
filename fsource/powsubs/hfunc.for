      FUNCTION HFUNC(X,Y,IFLAG,DHDY)

!PURPOSE: Compute exp(x)exp(y**2)erfc(y) and partial derivatives

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!PSEUDOCODE:

!CALLING ARGUMENTS:

      REAL*4        X,Y                 !Arguments
      INTEGER*4     IFLAG               !Control for calculation method
      REAL*4        DHDY                !Partial derivative wrt y
      REAL*4        HFUNC               !Value of function

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

!SUBROUTINES CALLED:

!FUNCTION DEFINITIONS:

      REAL*4        GERFC               !Complementary error function

!DATA STATEMENTS:

      DATA CONST/1.128379167/,CONS2/1.273239545/                  !2/SQRT(PI),4/PI

!CODE:

      T = Y**2+X
      HFUNC = 0.0
      DHDY = 0.0
      IF ( Y.GE.6.000 .OR. IFLAG.EQ.1 ) THEN                        !Compute function - protect against underflows
        T1 = SQRT(Y**2+CONS2)
        IF ( X.GT.-75.0 ) HFUNC = EXP(X)                        !Upper inequality
        HFUNC = HFUNC*CONST/(Y+T1)
        DHDY = -HFUNC/T1
      ELSE IF ( Y.LE.-6.375 ) THEN
        IF ( T.GT.-75.0 ) HFUNC = 2.0*EXP(T)                        !Upper inequality
        DHDY = 2.0*Y*HFUNC
      ELSE
        IF ( T.GT.-75.0 ) HFUNC = EXP(T)
        HFUNC = HFUNC*GERFC(Y)
        IF ( X.GT.-75.0 ) DHDY = EXP(X)
        DHDY = 2.0*Y*HFUNC-CONST*DHDY
      END IF
      RETURN
      END
