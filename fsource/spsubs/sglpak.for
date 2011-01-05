      SUBROUTINE SGLPAK(L,IER)

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

C       THIS PROGRAM WAS DEVELOPED FOR
C                    THE DIVISION OF CHEMISTRY
C                               OF
C               THE NATIONAL RESEARCH COUNCIL OF CANADA
C                               BY
C       ALLEN C. LARSON, P.O.BOX 5898, SANTA FE, NM 87502,USA

      DIMENSION     L(4)                

      IF ( L(2).LT.12 ) IER=4
      IF ( L(2).GT.17 ) IER=4
      L(1) = L(2)
      L(2) = 18
      RETURN
      END
