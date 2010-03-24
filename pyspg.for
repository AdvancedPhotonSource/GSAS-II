C Space group access routines for python

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
