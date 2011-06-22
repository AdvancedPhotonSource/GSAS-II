      SUBROUTINE PYQLMNINIT()
      CALL QLMNINIT
      RETURN
      END

      SUBROUTINE PYPLMPSI(L,I,PHI,PCRS)
Cf2py intent(in) L
Cf2Py intent(in) I
Cf2py intent(in) PHI
Cf2py intent(out) PCRS

      INTEGER*4     L
      INTEGER*4     I
      REAL*4        PHI
      REAL*4        PCRS

      CALL PLMPSI(L,I,PHI,PCRS)
      RETURN
      END

