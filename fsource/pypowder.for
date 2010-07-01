      SUBROUTINE PYPSVFCJ(YCALC,DTT,TTHETA,SIG,GAM,SPH,SMH,
     1  PRFUNC,DPRDT,DPRDY,SIGPART,GAMPART,SPHPART,SMHPART)
C DTT in degrees
C TTHETA in degrees
C SPH is S/L + H/L
C SML is S/L - H/L (frequently zero)
Cf2py intent(in) YCALC
Cf2py intent(in) DTT
Cf2py intent(in) TTHETA
Cf2py intent(in) SIG
Cf2py intent(in) GAM
Cf2py intent(in) SPH
Cf2py intent(in) SMH
Cf2py intent(out) PRFUNC
Cf2py intent(out) DPRDT
Cf2py intent(out) DPRDY
Cf2py intent(out) SIGPART
Cf2py intent(out) GAMPART
Cf2py intent(out) SPHPART
Cf2py intent(out) SMHPART
      SL = (SPH+SMH)/2.0
      HL = (SPH-SMH)/2.0
      FW = (2.355*SQRT(SIG)+GAM)/100.0
      FMIN = 10.0*(-FW-SPH*COSD(TTHETA))
      FMAX = 15.0*FW
      IF ( DTT .GE. FMIN .AND. DTT .LE. FMAX ) THEN
        CALL PSVFCJ(DTT*100.,TTHETA*100.,SL,HL,SIG,GAM,
     1    PRFUNC,DPRDT,SLPART,HLPART,SIGPART,GAMPART)
        DPRDT = DPRDT*YCALC*100.0
        DPRDY = PRFUNC
        SIGPART = SIGPART*YCALC
        GAMPART = GAMPART*YCALC
        SPHPART = 0.5*(SLPART+HLPART)*YCALC
        SMHPART = 0.5*(SLPART-HLPART)*YCALC
      ELSE
        PRFUNC = 0.0
        DPRDT = 0.0
        DPRDY = 0.0
        SIGPART = 0.0
        GAMPART = 0.0
        SPHPART = 0.0
        SMHPART = 0.0
      END IF
      RETURN
      END

      SUBROUTINE BUILDMV(WDELT,W,M,DP,A,V)
Cf2py intent(in) WDELT
Cf2py intent(in) W
Cf2py intent(in) M
Cf2py intent(in) DP
Cf2py depend(M) DP
Cf2py intent(in,out) A
Cf2py depend(M) a
Cf2py intent(in,out) V
Cf2py depend(M) V
      REAL*4 DP(M),A(M,M),V(M)
      DO I=1,M
        V(I) = V(I)+WDELT*DP(I)
        DO J=1,M
          A(I,J) = A(I,J)+W*DP(I)*DP(J)
        END DO
      END DO
      RETURN
      END
      
