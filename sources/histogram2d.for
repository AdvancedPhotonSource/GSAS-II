      SUBROUTINE HISTOGRAM2D(N,X,Y,Z,NXBINS,NYBINS,XLIM,YLIM,DX,DY,
     1  NST,HST)

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
Cf2py intent(in) dx
Cf2py intent(in) dy
Cf2py intent(in,out) nst
Cf2py depend(nxbins,nybins) nst
Cf2py intent(in,out) hst
Cf2py depend(nxbins,nybins) hst

      IMPLICIT NONE
      INTEGER*8   N
      REAL*4      X(0:N-1),Y(0:N-1),Z(0:N-1)
      INTEGER*8   NXBINS,NYBINS
      REAL*8      XLIM(0:1),YLIM(0:1)
      REAL*4      NST(0:NXBINS-1,0:NYBINS-1)
      REAL*4      HST(0:NXBINS-1,0:NYBINS-1)

      INTEGER*4   I,J,K
      REAL*8      DX,DY
      REAL*4      DDX,DDY
      DO K=0,N-1
C        if ( mod(k,8000) .eq. 0 )
C     1    print *,k,x(k),xlim,y(k),ylim
        IF ( ( X(K).GE.XLIM(0) .AND. X(K).LT.XLIM(1) ) .AND.
     1    ( Y(K).GE.YLIM(0) .AND. Y(K).LT.YLIM(1) )) THEN
          DDX = (X(K)-XLIM(0))/DX
          I = MIN(MAX(0,INT(DDX)),NXBINS-1)
          DDY = (Y(K)-YLIM(0))/DY
          J = MIN(MAX(0,INT(DDY)),NYBINS-1)
          NST(I,J) = NST(I,J)+1.0
          HST(I,J) = HST(I,J)+Z(K)
C          if ( mod(k,8000) .eq. 0 )
C     1      print *,i,j,nst(i,j),hst(i,j)
        END IF
      END DO
      RETURN
      END
