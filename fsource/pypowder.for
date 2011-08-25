      SUBROUTINE PYPSVFCJ(NPTS,DTT,TTHETA,SIG,GAM,SPH,PRFUNC)
C DTT in degrees
C TTHETA in degrees
C SPH is S/L + H/L
C RETURNS FUNCTION ONLY
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
      FW = (2.355*SQRT(SIG)+GAM)/100.0
      FMIN = 10.0*(-FW-SPH*COSD(TTHETA))
      FMAX = 15.0*FW
      DO I=0,NPTS-1
        CALL PSVFCJ(DTT(I)*100.,TTHETA*100.,SIG,GAM,SPH,
     1    PRFUNC(I),DPRDT,SIGPART,GAMPART,SLPART)
      END DO
      RETURN
      END

      SUBROUTINE PYDPSVFCJ(NPTS,DTT,TTHETA,SIG,GAM,SPH,PRFUNC,
     1  DPRDT,SIGPART,GAMPART,SLPART)
C DTT in degrees
C TTHETA in degrees
C SPH is S/L + H/L
C RETURNS FUNCTION & DERIVATIVES
Cf2py intent(in) NPTS
Cf2py intent(in) DTT
cf2py depend(NPTS) DTT
Cf2py intent(in) TTHETA
Cf2py intent(in) SIG
Cf2py intent(in) GAM
Cf2py intent(in) SPH
Cf2py intent(out) PRFUNC
Cf2py depend(NPTS) PRFUNC
Cf2py intent(out) DPRDT
Cf2py depend(NPTS) DPRDT
Cf2py intent(out) SIGPART
Cf2py depend(NPTS) SIGPART
Cf2py intent(out) GAMPART
Cf2py depend(NPTS) GAMPART
Cf2py intent(out) SLPART
Cf2py depend(NPTS) SLPART

      REAL*4 DTT(0:NPTS-1),DPRDT(0:NPTS-1),SIGPART(0:NPTS-1),
     1  GAMPART(0:NPTS-1),SLPART(0:NPTS-1),PRFUNC(0:NPTS-1)
      FW = (2.355*SQRT(SIG)+GAM)/100.0
      FMIN = 10.0*(-FW-SPH*COSD(TTHETA))
      FMAX = 15.0*FW
      DO I=0,NPTS-1
        CALL PSVFCJ(DTT(I)*100.,TTHETA*100.,SIG,GAM,SPH,
     1    PRFUNC(I),DPRDT(I),SIGPART(I),GAMPART(I),SLPART(I))
        DPRDT(I) = DPRDT(I)*100.
      END DO
      RETURN
      END

      SUBROUTINE PYPSVFCJO(NPTS,DTT,TTHETA,SIG,GAM,SPH,PRFUNC)
C DTT in degrees
C TTHETA in degrees
C SPH is S/L + H/L
C RETURNS FUNCTION ONLY
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
      FW = (2.355*SQRT(SIG)+GAM)/100.0
      FMIN = 10.0*(-FW-SPH*COSD(TTHETA))
      FMAX = 15.0*FW
      DO I=0,NPTS-1
        CALL PSVFCJO(DTT(I)*100.,TTHETA*100.,SIG,GAM,SPH/2.0,SPH/2.0,
     1    PRFUNC(I),DPRDT,SIGPART,GAMPART,SLPART,HLPART)
      END DO
      RETURN
      END

      SUBROUTINE PYDPSVFCJO(NPTS,DTT,TTHETA,SIG,GAM,SHL,PRFUNC,
     1  DPRDT,SIGPART,GAMPART,SLPART)
C DTT in degrees
C TTHETA in degrees
C SPH is S/L + H/L
C RETURNS FUNCTION & DERIVATIVES
Cf2py intent(in) NPTS
Cf2py intent(in) DTT
cf2py depend(NPTS) DTT
Cf2py intent(in) TTHETA
Cf2py intent(in) SIG
Cf2py intent(in) GAM
Cf2py intent(in) SHL
Cf2py intent(out) PRFUNC
Cf2py depend(NPTS) PRFUNC
Cf2py intent(out) DPRDT
Cf2py depend(NPTS) DPRDT
Cf2py intent(out) SIGPART
Cf2py depend(NPTS) SIGPART
Cf2py intent(out) GAMPART
Cf2py depend(NPTS) GAMPART
Cf2py intent(out) SLPART
Cf2py depend(NPTS) SLPART

      REAL*4 DTT(0:NPTS-1),DPRDT(0:NPTS-1),SIGPART(0:NPTS-1),
     1  GAMPART(0:NPTS-1),SLPART(0:NPTS-1),PRFUNC(0:NPTS-1)
      FW = (2.355*SQRT(SIG)+GAM)/100.0
      FMIN = 10.0*(-FW-SPH*COSD(TTHETA))
      FMAX = 15.0*FW
      DO I=0,NPTS-1
        CALL PSVFCJO(DTT(I)*100.,TTHETA*100.,SIG,GAM,SHL/2.,SHL/2.,
     1    PRFUNC(I),DPRDT(I),SIGPART(I),GAMPART(I),SPART,HPART)
          SLPART(I) = SPART
        DPRDT(I) = DPRDT(I)*100.
      END DO
      RETURN
      END

