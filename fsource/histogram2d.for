      SUBROUTINE HISTOGRAM2D(N,X,Y,Z,NXBINS,NYBINS,XLIM,YLIM,
     1  NST,HST,HSTX,HSTY)

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
Cf2py intent(in,out) NST
Cf2py depend(NXBINS,NYBINS) 
Cf2py intent(in,out) HST
Cf2py depend(NXBINS,NYBINS) 
Cf2py intent(in,out) HSTX
Cf2py depend(NXBINS) 
Cf2py intent(in,out) HSTY
Cf2py depend(NYBINS) 

      IMPLICIT NONE
      INTEGER*4   N
      REAL*4      X(0:N-1),Y(0:N-1),Z(0:N-1)
      INTEGER*4   NXBINS,NYBINS
      REAL*8      XLIM(0:1),YLIM(0:1)
      INTEGER*8   NST(0:NXBINS-1,0:NYBINS-1)
      REAL*8      HST(0:NXBINS-1,0:NYBINS-1)
      REAL*8      HSTX(0:NXBINS),HSTY(0:NYBINS)

      INTEGER*4   I,J,K
      REAL*8      DX,DY

      DX = (XLIM(1)-XLIM(0))/FLOAT(NXBINS)
      DY = (YLIM(1)-YLIM(0))/FLOAT(NYBINS)

      DO I=0,NXBINS
        HSTX(I) = XLIM(0)+FLOAT(I)*DX
      END DO
      DO J=0,NYBINS
        HSTY(J) = YLIM(0)+FLOAT(J)*DY
      END DO

      DO K=0,N
        IF ( ( X(K) .GE. XLIM(0) .AND. X(K) .LE. XLIM(1)) .AND.
     1    (Y(K) .GE. YLIM(0) .AND. Y(K). LE. YLIM(1)) ) THEN
          I = INT((X(K)-XLIM(0))/DX)
          J = INT((Y(K)-YLIM(0))/DY)
          NST(I,J) = NST(I,J)+1
          HST(I,J) = HST(I,J)+Z(K)
        END IF
      END DO
      RETURN
      END
