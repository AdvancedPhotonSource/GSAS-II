      SUBROUTINE SPOTMASK(N,X,Y,M,SPOTS,MASK)

Cf2py intent(in) N
Cf2py intent(in) X
Cf2py depend(N) X
Cf2py intent(in) Y
Cf2py depend(N) Y
Cf2py intent(in) M
Cf2py intent(in) SPOTS
Cf2py depend(M) SPOTS
Cf2py intent(in,out) MASK

      IMPLICIT NONE
      INTEGER*4    N,M
      REAL*4       X(0:N-1),Y(0:N-1)
      REAL*8       SPOTS(0:M-1,0:2)
      LOGICAL*1    MASK(0:1024*1024-1)

      INTEGER*4    I,K
      REAL*4       XYRAD2,XINTERS
      
      DO K=0,N-1
        MASK(K) = .FALSE.
        DO I=0,M-1
          XYRAD2 = (X(K)-SPOTS(I,0))**2+(Y(K)-SPOTS(I,1))**2
          IF ( XYRAD2 .LE. SPOTS(I,2) ) THEN
                  MASK(K) = .NOT.MASK(K)
          END IF
        END DO
      END DO

      RETURN
      END
