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
      
      SUBROUTINE SGFORPY(SPG,LAUE,SGINV,SGLATT,SGUNIQ,SGPOL,
     1  SGNOPS,SGMTRX,SGTRNS,IERR)
Cf2py intent(in)  SPG
Cf2py intent(out) LAUE
Cf2py intent(out) SGINV
Cf2py intent(out) SGLATT
Cf2py intent(out) SGUNIQ
Cf2py intent(out) SGPOL
Cf2py intent(out) SGNOPS
Cf2py intent(out) SGMTRX
Cf2py intent(out) SGTRNS
Cf2py intent(out) IERR

      CHARACTER*(20) SPG
      INTEGER*4     LAUE,SGINV,SGLATT,SGUNIQ,SGNOPS,IERR
      REAL*4        SGMTRX(24,3,3),SGTRNS(24,3)
      REAL*4        RT(5,4,25),CEN(3,4)
      INTEGER*4     JRT(3,5,24)


      CALL SGROUPNP(SPG,LAUE,SGUNIQ,SGINV,SGLATT,SGNOPS,SGPOL,JRT,
     1  CEN,SGNCEN,RT,IERR)

      DO K=1,SGNOPS
        DO I=1,3
          DO J=1,3
            SGMTRX(K,I,J) = JRT(I,J,K)
            SGTRNS(K,I) = JRT(I,4,K)/12.
          END DO
        END DO
      END DO
      RETURN
      END
