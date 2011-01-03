      REAL*4 FUNCTION SIND(ARG)

!PURPOSE: Calculate sine from angle in deg.

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:

      REAL*4        ARG                 !Sine argument in degrees

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA  RPD/0.017453292519943/

!CODE:

      SIND = SIN(ARG*RPD)
      RETURN
      END
