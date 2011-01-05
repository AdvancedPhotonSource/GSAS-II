      SUBROUTINE LORENTZ(DT,GAM,FUNC,DLDT,DLDG)

!PURPOSE:Calculate Lorentzian function & derivatives

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!PSEUDOCODE:

!CALLING ARGUMENTS:

      REAL*4        DT                  !Delta
      REAL*4        GAM                 !Coefficient
      REAL*4        FUNC                !Function
      REAL*4        DLDT                !df/dt
      REAL*4        DLDG                !df/dg

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

!SUBROUTINES CALLED:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA            PI/3.1415292654/
!CODE:

      FNORM = 2.0*GAM/PI
      DENOM = GAM**2+4.0*DT**2
      IF ( ABS(GAM).GT.5.0E-6 ) THEN
        FUNC = FNORM/DENOM
        DLDT = -8.0*DT*FUNC/DENOM
        DLDG = (8.0*DT**2-2.0*GAM**2)/(PI*DENOM**2)
      ELSE
        FUNC = 0.0
        DLDT = 0.0
        IF ( ABS(DT).LT.1.0E-5 ) THEN
          DLDG = 0.0
        ELSE
          DLDG = 2.0/(PI*DT**2)
        END IF
      END IF
      RETURN
      END
