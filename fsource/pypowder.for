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
      REAL*4 TTHETA,SIG,GAM,SPH
      INTEGER*4 NPTS,I
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

      INTEGER*4 NPTS
      REAL*4 TTHETA,SIG,GAM,SPH
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

      INTEGER*4 NPTS
      REAL*4 TTHETA,SIG,GAM,SPH
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

      INTEGER*4 NPTS
      REAL*4 TTHETA,SIG,GAM,SHL
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

C Fortran (fast) linear interpolation -- B.H. Toby 9/2011
      SUBROUTINE PYFINTERP(NIN,XIN,YIN,NOUT,XOUT,YOUT)
C XIN(1:NIN) and YIN(1:NIN) are arrays of (x,y) points to be interpolated
C     Values must be sorted increasing in XIN
C XOUT(1:NOUT) is an array of x values, must also be sorted increasing in x
C    XOUT may contain values smaller than XIN(1) or larger than XIN(NIN)
C RETURNS interpolated y values corresponding to XOUT. Values outside the 
C   range of XIN are set to zero.
C Needs a way to signal an error if XIN or XOUT is not sorted -- for now stops
Cf2py intent(in) NIN
Cf2py intent(in)  XIN
cf2py depend(NIN) XIN
Cf2py intent(in)  YIN
cf2py depend(NIN) YIN
Cf2py intent(in) NOUT
Cf2py intent(in)   XOUT
cf2py depend(NOUT) XOUT
Cf2py intent(out)  YOUT
cf2py depend(NOUT) YOUT

      INTEGER NIN,NOUT
      REAL XIN(NIN),YIN(NIN)
      REAL XOUT(NOUT),YOUT(NOUT)
      INTEGER IERROR
      REAL X,F
      INTEGER IIN,I

      IERROR = 1
      IIN = 1
      X = XOUT(1)
      DO I=1,NOUT
         IF (X .GT. XOUT(I)) STOP ! test if Xout not sorted
         X = XOUT(I)
         IF (X .LT. XIN(1) .OR. X .GT. XIN(NIN) ) THEN
            YOUT(I) = 0.0
         ELSE
            DO WHILE (X .GT.  XIN(IIN+1))
               IF (XIN(IIN) .GT. XIN(IIN+1)) STOP ! test if Xin not sorted
               IIN = IIN + 1
             ENDDO
             F = (X - XIN(IIN)) / (XIN(IIN+1) - XIN(IIN))
             YOUT(I) = (1.-F)*YIN(IIN) + F*YIN(IIN+1)
          ENDIF
          !write (*,*) xout(i),iin,f,yout(i)
      END DO
      IERROR = 0
      RETURN
      END
