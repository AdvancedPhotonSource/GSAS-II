      FUNCTION PRODVMV(X,A,Y)

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 
!
! ROUTINE TO COMPUTE THE DOT PRODUCT OF TWO 3 VECTORS

      REAL          X(3),Y(3),A(3,3),PRODVMV 

      T1 = 0.0
      DO J=1,3
        T1 = T1+X(J)*(Y(1)*A(J,1)+Y(2)*A(J,2)+Y(3)*A(J,3))
      END DO
      PRODVMV = T1
      RETURN
      END
