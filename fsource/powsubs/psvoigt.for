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
      REAL*4        EPS                 !Are values different
! Local variables saved between calls
      REAL*4        prev_sig,prev_gam
      REAL*4        eta,fwhm,frac,dsdg,dsdl,sumhm,dedf,sqsg
      save          eta,fwhm,frac,prev_sig,prev_gam,dsdg,dsdl,sumhm,
     1                   dedf,sqsg
      
!SUBROUTINES CALLED:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA STOFW/2.35482005/            !2*SQRT(2LN2)
      DATA SQ2PI/2.506628275/            !SQRT(2PI)
      DATA COFT/1.0,2.69269,2.42843,4.47163,0.07842,1.0/
      DATA COFN/1.36603,-0.47719,0.11116/
      DATA prev_sig/-1.0/
      DATA prev_gam/-1.0/
      DATA eps/1.0e-10/           !Threshold for recalculation

!CODE:

! Check for repeat call
      if (abs(prev_sig-sig) .gt. eps .or.   
     1  (abs(prev_gam-gam).gt.eps)) then !need to recalculate
         prev_sig = sig
         prev_gam = gam
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
      end if   !end of recalculation step 
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
      
      SUBROUTINE PSVOIGT2(DX,SIG,GAM,FUNC,DFDX,DFDS,DFDG)

!PURPOSE: Compute function & derivatives pseudovoigt - unfinished; 
!   no derivatives
!pseudo Voigt W.I.F. David, J. Appl. Cryst. 19, 63-64 (1986)

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

      REAL*4        GNORM               !Gaussian Normalization constant
      REAL*4        COFEG(7),COFEL(7),COFGG(6),COFGL(6)    

!SUBROUTINES CALLED:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA STOFW/2.35482005/            !2*SQRT(2LN2)
      DATA SQ2PI/2.506628275/            !SQRT(2PI)
      DATA ITYPE/1/
      DATA COFEG/0.00268,0.75458,2.88898,-3.85144,-.55765,3.03824,
     1  -1.27539/
      DATA COFEL/1.35248,0.41168,-2.18731,6.42452,-10.29036,6.88093,
     1  -1.59194/
      DATA COFGG/-.50734,-.22744,1.63804,-2.28532,1.31943,0.0/
      DATA COFGL/-.99725,1.14594,2.56150,-6.52088,5.82647,-1.91086/

!CODE:

      SQSG = MAX(SQRT(SIG),0.001)
      FWHG = STOFW*SQSG
      FW = FWHG+GAM
      R1 = FWHG/FW
      R17 = R1
      R2 = GAM/FWHG
      R27 = R2
      DSDL = 0.0
      DSDG = 0.0
      ETAG = 0.0
      ETAL = 0.0
      DO ITRM=1,7
        ETAG = ETAG+R17*COFEG(ITRM)
        ETAL = ETAL+R27*COFEL(ITRM)
        R17 = R17*R1
        R27 = R27*R2
      END DO
      WG = 1.0
      WL = 1.0
      R16 = R1
      R26 = R2
      DO ITRM=1,6
        WG = WG+R26*COFGG(ITRM)
        WL = WL+R16*COFGL(ITRM)
        R16 = R16*R1
        R26 = R26*R2
      END DO
      CALL LORENTZ(DX,WL,TL,DTLDT,DTLDFW)
      SIGP = (WG/STOFW)**2
      EX = MAX(-20.0,-0.5*DX**2/SIGP)
      TG = STOFW*EXP(EX)/(SQ2PI*FWHM)
      FUNC = ETAL*TL+ETAG*TG
        
! Unfinished - no derivatives

      RETURN
      END
