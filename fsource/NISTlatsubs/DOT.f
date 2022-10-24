      SUBROUTINE DOT (UX,VX,WX,UY,VY,WY,AI,BI,CI,AR,BR,GR,DOTXY)
C
C
C     SUBROUTINE  'DOT' ...
C        CALCULATE DOT PRODUCT OF VECTORS X AND Y
C
C
C     --- VECTOR X = (UX A + VX B + WX C)
C         VECTOR Y = (UY A + VY B + WY C)
C
      DOTXY = UX*UY*(AI**2) + VX*VY*(BI**2) + WX*WY*(CI**2) +
     $       (VX*WY + WX*VY)*BI*CI*COS(AR) +
     $       (UX*WY + WX*UY)*AI*CI*COS(BR) +
     $       (UX*VY + VX*UY)*AI*BI*COS(GR)
      RETURN
      END
