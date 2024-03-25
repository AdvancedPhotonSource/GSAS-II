      SUBROUTINE MULTIP
      COMMON /MATR1/ U1,V1,W1,U2,V2,W2,U3,V3,W3
      COMMON /MATR2/ U(9),T(9)
C
C
C     SUBROUTINE 'MULTIP' ...
C        MULTIPLY TWO 3X3 MATRICES
C        WHEN THIS SUBROUTINE IS CALLED (RSS FUNCTION),
C        THE TOTAL TRANSFORMATION MATRIX IS UPDATED
C
C
C
C     --- (PRODUCT MATRIX) = (U-MATRIX) (T-MATRIX)
      U1 = T(1)*U(1) + T(4)*U(2) + T(7)*U(3)
      V1 = T(2)*U(1) + T(5)*U(2) + T(8)*U(3)
      W1 = T(3)*U(1) + T(6)*U(2) + T(9)*U(3)
      U2 = T(1)*U(4) + T(4)*U(5) + T(7)*U(6)
      V2 = T(2)*U(4) + T(5)*U(5) + T(8)*U(6)
      W2 = T(3)*U(4) + T(6)*U(5) + T(9)*U(6)
      U3 = T(1)*U(7) + T(4)*U(8) + T(7)*U(9)
      V3 = T(2)*U(7) + T(5)*U(8) + T(8)*U(9)
      W3 = T(3)*U(7) + T(6)*U(8) + T(9)*U(9)
      RETURN
      END
