      SUBROUTINE PYPSVFCJ(NPTS,DTT,TTHETA,SIG,GAM,SPH,PRFUNC)
C DTT in degrees
C TTHETA in degrees
C SPH is S/L + H/L
Cf2py intent(in) NPTS
Cf2py intent(in) DTT
cf2py depend(NPTS) DTT
Cf2py intent(in) TTHETA
Cf2py intent(in) SIG
Cf2py intent(in) GAM
Cf2py intent(in) SPH
Cf2py intent(out) PRFUNC
Cf2py depend(NPTS) PRFUNC

      REAL*4 DTT(0:NPTS-1),PRFUNC(0:NPTS-1)
      SL = SPH/2.0
      FW = (2.355*SQRT(SIG)+GAM)/100.0
      FMIN = 10.0*(-FW-SPH*COSD(TTHETA))
      FMAX = 15.0*FW
      DO I=0,NPTS-1
        CALL PSVFCJ(DTT(I)*100.,TTHETA*100.,SL,SL,SIG,GAM,
     1    PRFUNC(I),DPRDT,SLPART,HLPART,SIGPART,GAMPART)
      END DO
      RETURN
      END

      SUBROUTINE PYPSVFCJ2(NPTS,DTT,TTHETA,SIG,GAM,SPH,PRFUNC)
C DTT in degrees
C TTHETA in degrees
C SPH is S/L + H/L
Cf2py intent(in) NPTS
Cf2py intent(in) DTT
cf2py depend(NPTS) DTT
Cf2py intent(in) TTHETA
Cf2py intent(in) SIG
Cf2py intent(in) GAM
Cf2py intent(in) SPH
Cf2py intent(out) PRFUNC
Cf2py depend(NPTS) PRFUNC

      REAL*4 DTT(0:NPTS-1),PRFUNC(0:NPTS-1)
      SL = SPH/2.0
      FW = (2.355*SQRT(SIG)+GAM)/100.0
      FMIN = 10.0*(-FW-SPH*COSD(TTHETA))
      FMAX = 15.0*FW
      DO I=0,NPTS-1
        CALL PSVFCJ2(DTT(I)*100.,TTHETA*100.,SL,SL,SIG,GAM,
     1    PRFUNC(I))
      END DO
      RETURN
      END
