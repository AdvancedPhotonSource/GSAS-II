      SUBROUTINE HISTOGRAM2D(N,X,Y,Z,NXBINS,NYBINS,XLIM,YLIM,
     1  NST,HST,HSTX,HSTY)

Cf2py intent(in) n
Cf2py intent(in) x
Cf2py depend(n) x
Cf2py intent(in) y
Cf2py depend(n) y
Cf2py intent(in) z
Cf2py depend(n) z
Cf2py intent(in) nxbins
Cf2py intent(in) nybins
Cf2py intent(in) xlim       
Cf2py intent(in) ylim       
Cf2py intent(inout) nst
Cf2py depend(nxbins,nybins) nst
Cf2py intent(inout) hst
Cf2py depend(nxbins,nybins) hst
Cf2py intent(inout) hstx
Cf2py depend(nxbins) hstx
Cf2py intent(inout) hsty
Cf2py depend(nybins) hsty

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
