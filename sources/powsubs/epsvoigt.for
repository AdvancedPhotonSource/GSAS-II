      SUBROUTINE EPSVOIGT(DT,ALP,BET,SIG,GAM,FUNC,DFDX,DFDA,DFDB,
     1  DFDS,DFDG)

!PURPOSE: Compute function & derivatives exponential X pseudovoigt
!P.Tompson, D.E. Cox & J.B. Hastings (1987). J. Appl. Cryst.,20,79-83

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!PSEUDOCODE:

!CALLING ARGUMENTS:

      REAL*4        DT                  !Delta-TOF from center
      REAL*4        ALP                 !Exponential rise
      REAL*4        BET                 !Exponential decay
      REAL*4        SIG                 !Gaussian variance
      REAL*4        GAM                 !Lorentzian FWHM
      REAL*4        FUNC                !Value of pseudo-Voigt at DX
      REAL*4        DFDX                !dF/dta
      REAL*4        DFDA                !dF/da
      REAL*4        DFDB                !dF/db
      REAL*4        DFDS                !dF/ds
      REAL*4        DFDG                !dF/dg

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

      REAL*4        COFG(6),COFL(6)     !Linear combination coeffs
      REAL*4        ACOFG(7),ACOFL(7)   
      REAL*4        GNORM               !Gaussian Normalization constant
      REAL*4        COFT(6),COFN(3)     
      REAL*4        NORM                !Normalization constant
      REAL*4        RXA,IXA,RXB,IXB     !Exp-integral arguements
      REAL*4        RFA,IFA,RFB,IFB     !Exp-integral results

!SUBROUTINES CALLED:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA STOFW/2.35482005/            !2*SQRT(2LN2)
      DATA SQ2PI/2.506628275/            !SQRT(2PI)
      DATA PI/3.141592654/            !PI
      DATA COFT/1.0,2.69269,2.42843,4.47163,0.07842,1.0/
      DATA COFN/1.36603,-0.47719,0.11116/

!CODE:

      SQSGT = MAX(SQRT(SIG),0.00001)
      GAM = MAX(GAM,0.00001)
      FWHG = STOFW*SQSGT
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

      NORM = 0.5*ALP*BET/(ALP+BET)
      SIGP = (FWHM/STOFW)**2
      TJ1 = ALP*SIGP
      TJ2 = BET*SIGP
      SQSG = SQRT(2.0*SIGP)
      S2SG = 2.0*SQSG*SIGP
      Y1 = (TJ1+DT)/SQSG
      DY1DS = (TJ1-DT)/S2SG
      Y2 = (TJ2-DT)/SQSG
      DY2DS = (TJ2+DT)/S2SG
      EX1 = -0.5*(DT**2)/SIGP
      DEXDS = -EX1/SIGP
      H1 = HFUNC(EX1,Y1,0,DH1DY)
      H2 = HFUNC(EX1,Y2,0,DH2DY)

      TG = NORM*(H1+H2)

      DTGDA = 2.0*NORM*TG/(ALP*ALP)
      DTGDA = DTGDA+0.5*NORM*SQSG*DH1DY                  !OK
      DTGDB = 2.0*NORM*TG/(BET*BET)
      DTGDB = DTGDB+0.5*NORM*SQSG*DH2DY                  !OK
      DTGDFW = TG*DEXDS+NORM*(DH1DY*DY1DS+DH2DY*DY2DS)
      DTGDFW = 2.0*DTGDFW*FWHM/STOFW**2                  !OK
      DTGDT = -NORM*(ALP*H1-BET*H2)                        !OK

      RXA = -ALP*DT
      RXB = -BET*DT
      IXA = ALP*FWHM/2.0
      IXB = BET*FWHM/2.0
      CALL EXPINT(RXA,IXA,RFA,IFA)
      CALL EXPINT(RXB,IXB,RFB,IFB)
      TL = -2.0*NORM*(IFA+IFB)/PI
      DIVSOR = DT**2+FWHM**2/4.0
      DTLDT = -2.0*NORM*(ALP*IFA+BET*IFB+FWHM/DIVSOR)/PI      !OK

      DTLDFW = NORM*(-ALP*RFA-BET*RFB-2.0*DT/DIVSOR)/PI

      DTLDA = 2.0*NORM*TL/ALP**2
      DTLDA = DTLDA+2.0*NORM*(DT*IFA-FWHM*RFA/2.0)/PI            !OK
      DTLDB = 2.0*NORM*TL/BET**2
      DTLDB = DTLDB+2.0*NORM*(DT*IFB-FWHM*RFB/2.0)/PI            !OK

      FUNC = ETA*TL+(1.0-ETA)*TG                        !OK

      DFDX = ETA*DTLDT+(1.0-ETA)*DTGDT

      DFWDG = 0.2*DSDG*FWHM/SUMHM
      DFRDG = -FRAC*DFWDG/FWHM

      DFDS = DEDF*DFRDG*(TL-TG)
      DFDS = DFDS+(ETA*DTLDFW+(1.0-ETA)*DTGDFW)*DFWDG
      DFDS = 0.5*DFDS*STOFW/SQSGT                        !OK

      DFWDL = 0.2*DSDL*FWHM/SUMHM
      DFRDL = (1.0-FRAC*DFWDL)/FWHM

      DFDG = DEDF*DFRDL*(TL-TG)
      DFDG = DFDG+(ETA*DTLDFW+(1.0-ETA)*DTGDFW)*DFWDL            !OK

      DFDA = ETA*DTLDA+(1.0-ETA)*DTGDA                  !OK

      DFDB = ETA*DTLDB+(1.0-ETA)*DTGDB                  !OK

      RETURN
      END
