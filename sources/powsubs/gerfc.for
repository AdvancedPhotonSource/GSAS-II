      FUNCTION GERFC(Y)

!PURPOSE: This routine returns the Complementary ERROR Function of a single precision variable (Y)

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:

      REAL*4        Y                   !
      REAL*4        GERFC                !Result

!INCLUDE STATEMENTS:

      DIMENSION     P(5),Q(4),P1(9),Q1(8),P2(6),Q2(5) 
      REAL*4        P,Q,P1,Q1,P2,Q2,XMIN,XLARGE,SQRPI,X, 
     1                RES,XSQ,XNUM,XDEN,XI
      INTEGER       ISW,I               

!DATA STATEMENTS:

      DATA P(1)/113.86415/,
     1  P(2)/377.48524/,
     1  P(3)/3209.3776/,
     1  P(4)/.18577771/,
     1  P(5)/3.1611237/
      DATA Q(1)/244.02464/,
     1  Q(2)/1282.6165/,
     1  Q(3)/2844.2368/,
     1  Q(4)/23.601291/
      DATA P1(1)/8.8831498/,
     1  P1(2)/66.119190/,
     1  P1(3)/298.63514/,
     1  P1(4)/881.95222/,
     1  P1(5)/1712.0476/,
     1  P1(6)/2051.0784/,
     1  P1(7)/1230.3394/,
     1  P1(8)/2.1531154E-8/,
     1  P1(9)/.56418850/
      DATA Q1(1)/117.69395/,
     1  Q1(2)/537.18110/,
     1  Q1(3)/1621.3896/,
     1  Q1(4)/3290.7992/,
     1  Q1(5)/4362.6191/,
     1  Q1(6)/3439.3677/,
     1  Q1(7)/1230.3394/,
     1  Q1(8)/15.744926/
      DATA P2(1)/-3.6034490E-01/,
     1  P2(2)/-1.2578173E-01/,
     1  P2(3)/-1.6083785E-02/,
     1  P2(4)/-6.5874916E-04/,
     1  P2(5)/-1.6315387E-02/,
     1  P2(6)/-3.0532663E-01/
      DATA Q2(1)/1.8729528/,
     1  Q2(2)/5.2790510E-01/,
     1  Q2(3)/6.0518341E-02/,
     1  Q2(4)/2.3352050E-03/,
     1  Q2(5)/2.5685202/
      DATA XMIN/1.0E-10/,XLARGE/6.375/
      DATA SQRPI/.56418958/

!CODE:

      IF ( Y.GT.0.0 ) THEN
        X = Y
        ISW = 1
      ELSE
        X = -Y
        ISW = -1
      END IF
      XSQ = X*X
      IF ( X.LT.0.477 ) THEN
        IF ( X.LT.XMIN ) THEN
          RES = X*P(3)/Q(3)
        ELSE
          XNUM = P(4)*XSQ+P(5)
          XDEN = XSQ+Q(4)
          DO I = 1,3
            XNUM = XNUM*XSQ+P(I)
            XDEN = XDEN*XSQ+Q(I)
          END DO
          RES = X*XNUM/XDEN
        END IF
        IF ( ISW.EQ.-1 ) RES = -RES
        RES = 1.0-RES
      ELSE
        IF ( X.LT.XLARGE ) THEN
          IF ( X.LE.4.0 ) THEN
            XNUM = P1(8)*X+P1(9)
            XDEN = X+Q1(8)
            DO I=1,7
              XNUM = XNUM*X+P1(I)
              XDEN = XDEN*X+Q1(I)
            END DO
            RES = XNUM/XDEN
          ELSE
            XI = 1.0/XSQ
            XNUM = P2(5)*XI+P2(6)
            XDEN = XI+Q2(5)
            DO I = 1,4
              XNUM = XNUM*XI+P2(I)
              XDEN = XDEN*XI+Q2(I)
            END DO
            RES = (SQRPI+XI*XNUM/XDEN)/X
          END IF
          RES = RES*EXP(-XSQ)
        ELSE
          RES = 0.0
        END IF
        IF ( ISW.EQ.-1 ) RES = 2.0-RES
      END IF
      GERFC = RES
      RETURN
      END
