      SUBROUTINE POLYMASK(N,X,Y,M,POLY,MASK)

Cf2py intent(in) N
Cf2py intent(in) X
Cf2py depend(N) X
Cf2py intent(in) Y
Cf2py depend(N) Y
Cf2py intent(in) M
Cf2py intent(in) POLY
Cf2py depend(M) POLY
Cf2py intent(in,out) MASK

      IMPLICIT NONE
      INTEGER*4    N,M
      REAL*4       X(0:N-1),Y(0:N-1)
      REAL*8       POLY(0:M-1,0:1)
      LOGICAL*1    MASK(0:1024*1024-1)

      INTEGER*4    I,K
      REAL*4       P1X,P1Y,P2X,P2Y,XINTERS
      
      DO K=0,N-1
        MASK(K) = .FALSE.
        DO I=0,M-1
          P2X = POLY(I,0)
          P2Y = POLY(I,1)
          IF (Y(K) .GT. MIN(P1Y,P2Y)) THEN
            IF (Y(K) .LE. MAX(P1Y,P2Y)) THEN
              IF (X(K) .LE. MAX(P1X,P2X)) THEN
                IF (P1Y .NE.P2Y) THEN
                  XINTERS = (Y(K)-P1Y)*(P2X-P1X)/(P2Y-P1Y)+P1X
                END IF
                IF ( (P1X .EQ. P2X) .OR. (X(K) .LE. XINTERS) ) THEN
                  MASK(K) = .NOT.MASK(K)
                END IF
              END IF
            END IF
          END IF
          P1X = P2X
          P1Y = P2Y
        END DO
      END DO

      RETURN
      END

