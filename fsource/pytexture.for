      SUBROUTINE PYPLMPSI(L,I,NPHI,PHI,PCRS,DPDPS)
Cf2py intent(in) L
Cf2Py intent(in) I
Cf2py intent(in) NPHI
Cf2py intent(in) PHI
Cf2py intent(out) PCRS
Cf2py intent(out) DPDPS
Cf2py depend(NPHI) PHI,PCRS,DPDPS

      INTEGER*4     L
      INTEGER*4     I
      REAL*4        PHI(0:NPHI-1)
      REAL*4        PCRS(0:NPHI-1)
      REAL*4        DPDPS(0:NPHI-1)

      DO K = 0,NPHI-1
         CALL PLMPSI(L,I,PHI(K),PCRS(K),DPDPS(K))
      END DO
      RETURN
      END

      SUBROUTINE PYQLMNINIT()
      CALL QLMNINIT
      RETURN
      END


