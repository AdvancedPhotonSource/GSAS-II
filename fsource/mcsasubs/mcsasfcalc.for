      SUBROUTINE MCSASFCALC(INV,NTD,TDATA,MDATA,XDATA,MUL,NFFS,FFS,
     1  NUNIQ,UNIQ,PHI,ICALC)
     
      LOGICAL*4 INV
      INTEGER*4 NTD,MUL,NFFS,NUNIQ,I,J,K,TDATA(0:NTD-1)
      REAL*8 XDATA(0:3*NTD-1),UNIQ(0:3*NUNIQ-1)
      REAL*8 MDATA(0:NTD-1),FFS(0:NFFS-1)
      REAL*8 ICALC,PHI(0:NUNIQ-1)
      REAL*8 PHASE,FF,FAS,FBS,TWOPI

      TWOPI = 8.0*ATAN(1.0)
      FAS = 0.
      FBS = 0.
      DO I=0,NTD-1
        FF = FFS(TDATA(I))*MDATA(I)/NUNIQ
        DO J=0,NUNIQ-1
          PHASE = 0.
          DO K=0,2
            PHASE = PHASE+UNIQ(3*J+K)*XDATA(3*I+K)
          END DO
          PHASE = PHASE+PHI(J)
          FAS = FAS+FF*COS(TWOPI*PHASE)
          IF ( .NOT. INV ) FBS = FBS+FF*SIN(TWOPI*PHASE)
        END DO
      END DO
      ICALC = (FAS**2+FBS**2)*MUL
      RETURN
      END
