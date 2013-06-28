      SUBROUTINE SETUPSINTAB
       
      REAL*4 SINTAB
      COMMON /SINTAB/SINTAB(0:24999)
      
      TWOPI = 8.0*ATAN(1.0)
      DO I=0,24999
        SINTAB(I) = SIN(TWOPI*I/10000.)
      END DO
        
      RETURN
      END
