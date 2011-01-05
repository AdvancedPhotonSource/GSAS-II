      REAL*4 FUNCTION COSD(ARG)

!PURPOSE: Calculate cosine from angle in deg.

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:

      REAL*4        ARG                 !Cosine argument in degrees

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA  RPD/0.017453292519943/

!CODE:

      COSD = COS(ARG*RPD)
      RETURN
      END
