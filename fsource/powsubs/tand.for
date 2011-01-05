      REAL*4 FUNCTION TAND(ARG)

!PURPOSE: Calculate tangent from angle in deg.

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:

      REAL*4        ARG                 !Tangent argument in degrees

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA  RPD/0.017453292519943/

!CODE:

      TAND = TAN(ARG*RPD)
      RETURN
      END
