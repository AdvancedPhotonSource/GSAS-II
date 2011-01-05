      SUBROUTINE PSVOIGT(DX,SIG,GAM,FUNC,DFDX,DFDS,DFDG)

!PURPOSE: Compute function & derivatives pseudovoigt
!pseudo Voigt P.Tompson, D.E. Cox & J.B. Hastings (1987). J. Appl. Cryst.,20,79-83

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!PSEUDOCODE:

!CALLING ARGUMENTS:

      REAL*4        DX                  !Delta-x from center
      REAL*4        SIG                 !Gaussian variance
      REAL*4        GAM                 !Lorentzian FWHM
      REAL*4        FUNC                !Value of pseudo-Voigt at DX
      REAL*4        DFDX                !dF/dx
      REAL*4        DFDS                !dF/ds
      REAL*4        DFDG                !dF/dg

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

      REAL*4        COFG(6),COFL(6)     !Linear combination coeffs
      REAL*4        ACOFG(7),ACOFL(7)   
      REAL*4        GNORM               !Gaussian Normalization constant
      REAL*4        COFT(6),COFN(3)     

!SUBROUTINES CALLED:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA STOFW/2.35482005/            !2*SQRT(2LN2)
      DATA SQ2PI/2.506628275/            !SQRT(2PI)
      DATA COFT/1.0,2.69269,2.42843,4.47163,0.07842,1.0/
      DATA COFN/1.36603,-0.47719,0.11116/
      DATA ITYPE/1/

!CODE:

      SQSG = MAX(SQRT(SIG),0.001)
      FWHG = STOFW*SQSG
      PGL = FWHG**5
      SUMHM = PGL
      DSDL = 0.0
      DSDG = 0.0
      DO ITRM=1,5
        PGL = PGL/FWHG
        DSDL = DSDL+FLOAT(ITRM)*COFT(ITRM+1)*PGL
        DSDG = DSDG+FLOAT(6-ITRM)*COFT(ITRM)*PGL
        PGL = PGL*GAM
        SUMHM = SUMHM+COFT(ITRM+1)*PGL
      END DO
      FWHM = EXP(0.2*LOG(SUMHM))
      FRAC = GAM/FWHM
      DEDF = 0.0
      PF = 1.0
      ETA = 0.0
      DO ITRM=1,3
        DEDF = DEDF+FLOAT(ITRM)*COFN(ITRM)*PF
        PF = PF*FRAC
        ETA = ETA+COFN(ITRM)*PF
      END DO
      CALL LORENTZ(DX,FWHM,TL,DTLDT,DTLDFW)
      SIGP = (FWHM/STOFW)**2
      EX = MAX(-20.0,-0.5*DX**2/SIGP)
      TG = STOFW*EXP(EX)/(SQ2PI*FWHM)
      FUNC = ETA*TL+(1.0-ETA)*TG

      TS = -2.0*(1.0-ETA)*TG*(EX+0.5)/FWHM
      DFDX = ETA*DTLDT-(1.0-ETA)*TG*DX/SIGP

      DFWDG = 0.2*DSDG*FWHM/SUMHM
      DFRDG = -FRAC*DFWDG/FWHM
      DFDS = DEDF*DFRDG*(TL-TG)+(ETA*DTLDFW+TS)*DFWDG
      DFDS = 0.5*DFDS*STOFW/SQSG

      DFWDL = 0.2*DSDL*FWHM/SUMHM
      DFRDL = (1.0-FRAC*DFWDL)/FWHM
      DFDG = DEDF*DFRDL*(TL-TG)+(ETA*DTLDFW+TS)*DFWDL

      RETURN
      END
