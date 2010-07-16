      SUBROUTINE HISTOGRAM2D(N,X,Y,Z,NXBINS,NYBINS,XLIM,YLIM,DX,DY,
     1  NST,HST)

Cf2py intent(in) N
Cf2py intent(in) X
Cf2py depend(N) X
Cf2py intent(in) Y
Cf2py depend(N) Y
Cf2py intent(in) Z
Cf2py depend(N) Z
Cf2py intent(in) NXBINS
Cf2py intent(in) NYBINS
Cf2py intent(in) XLIM
Cf2py intent(in) YLIM
Cf2py intent(in) DX
Cf2py intent(in) DY
Cf2py intent(in,out) NST
Cf2py depend(NXBINS,NYBINS) NST
Cf2py intent(in,out) HST
Cf2py depend(NXBINS,NYBINS) HST

      IMPLICIT NONE
      INTEGER*4   N
      REAL*4      X(0:N-1),Y(0:N-1),Z(0:N-1)
      INTEGER*4   NXBINS,NYBINS
      REAL*8      XLIM(0:1),YLIM(0:1)
      REAL*8      DX,DY
      INTEGER*8   NST(0:NXBINS-1,0:NYBINS-1)
      REAL*4      HST(0:NXBINS-1,0:NYBINS-1)

      INTEGER*4   I,J,K

      DO K=0,N-1
        IF ( ( X(K) .GE. XLIM(0) .AND. X(K) .LE. XLIM(1)) .AND.
     1    (Y(K) .GE. YLIM(0) .AND. Y(K). LE. YLIM(1)) ) THEN
          I = NINT((X(K)-XLIM(0))/DX+0.5)-1
          J = NINT((Y(K)-YLIM(0))/DY+0.5)-1
          NST(I,J) = NST(I,J)+1
          HST(I,J) = HST(I,J)+Z(K)
        END IF
      END DO
      RETURN
      END
