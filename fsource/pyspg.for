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
      INTEGER*4     LAUE,SGINV,SGLATT,SGUNIQ,SGNOPS,IERR,SGNCEN
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

      SUBROUTINE GENHKLPY(XH,NSYM,SGMTRX,SGTRNS,ICEN,NCV,SGCEN,JHK,
     1  HKL,IABSNT,MULP)
Cf2py intent(in)  XH
Cf2py intent(in)  NSYM
Cf2py intent(in)  SGMTRX
Cf2py intent(in)  SGTRNS
Cf2py depend(NSYM) SGMTRX,SGTRNS
Cf2py intent(in)  ICEN
Cf2py intent(in)  NCV
Cf2py intent(in)  SGCEN
Cf2py depend(NCV) SGCEN
Cf2py intent(out) JHK
Cf2py intent(out) HKL
Cf2py intent(out) IABSNT
Cf2py intent(out) MULP

      INTEGER*4     ICEN,NSYM
      REAL*4        SGMTRX(NSYM,3,3),SGTRNS(NSYM,3),SGCEN(NCV,3)
      REAL*4        CEN(3,4),HKL(4,24),XH(4)
      INTEGER*4     JRT(3,5,24),JHK,NCV

      DO J=1,NCV
        DO I=1,3
          CEN(I,J) = SGCEN(J,I)
        END DO
      END DO 
      DO K=1,NSYM
        DO I=1,3
          DO J=1,3
            JRT(I,J,K) = SGMTRX(K,I,J)*1.
            JRT(I,4,K) = SGTRNS(K,I)*12.
          END DO
        END DO
      END DO
      CALL GENHKL(XH,NSYM,JRT,ICEN,NCV,CEN,JHK,HKL,IABSNT,MULP)
      RETURN
      END