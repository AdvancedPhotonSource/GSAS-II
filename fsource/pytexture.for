      SUBROUTINE PYQLMNINIT()
      CALL QLMNINIT
      RETURN
      END

      SUBROUTINE PYPLMPSI(L,I,NPHI,PHI,PCRS)
Cf2py intent(in) L
Cf2Py intent(in) I
Cf2py intent(in) NPHI
Cf2py intent(in) PHI
Cf2py intent(out) PCRS
Cf2py depend(NPHI) PHI,PCRS

      INTEGER*4     L
      INTEGER*4     I
      REAL*4        PHI(0:NPHI-1)
      REAL*4        PCRS(0:NPHI-1)


      DO K = 0,NPHI-1
         CALL PLMPSI(L,I,PHI(K),PCRS(K))
      END DO
      RETURN
      END

