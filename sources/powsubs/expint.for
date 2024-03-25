      SUBROUTINE EXPINT(AX,AY,ANSX,ANSY)

!PURPOSE: Compute the real & imaginary components and derivatives of
!the function exp(z)E(z) where E(z) is the exponential integral function.
!Originally written by       W.I.F.David      1-SEP-84
!Minor changes by R.B. Von Dreele 1-Aug-1986

!PSEUDOCODE:

!CALLING ARGUMENTS:

      REAL*4        AX,AY               !Real & imaginary parts of argument z
      REAL*4        ANSX,ANSY           !Real & imaginary parts of result

!INCLUDE STATEMENTS:

!LOCAL VARIABLES:

      REAL*4        ZX,ZY               
      REAL*4        CLNX,CLNY           
      REAL*4        TEMX,TEMY           
      REAL*4        SUBX,SUBY           
      REAL*4        EULERX              
      REAL*4        ZANSX,ZANSY         
      REAL*4        DEDZX,DEDZY         
      REAL*4        RATX,RATY           
      REAL*4        ADDX,ADDY           

!SUBROUTINES CALLED:

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA EULERX/0.5772156649/

!CODE:

      ZX = AX                              !Z=CMPLX(AX,AY)
      ZY = AY
      AYCH= -1.
      IF (AX .LT. 0.) AYCH= 0.125*AX+4.
      A= SQRT(AX*AX+AY*AY)
      IF (A .LT. 4.0 .OR. ABS(AY) .LT. AYCH) THEN
        ITER= 10 + 4*NINT(A)
        CLNX = 0.5*LOG(ZX**2+ZY**2)            !CLN=CLOG(Z)
        CLNY = ATAN2(ZY,ZX)
        TEMX = FLOAT(ITER)/FLOAT((ITER+1)*(ITER+1))
        TEMY = ZY*TEMX
        TEMX = ZX*TEMX
        SUMX = TEMX                        !SUM=TEM
        SUMY = TEMY
        DO J= ITER,2,-1
          T= 1./FLOAT(J)
          SUBX = 1.- SUMX                  !SUB=1.0-SUM
          SUBY = -SUMY
          T = (1.+T)*(1-T*T)
          TEMX = TEMX*T
          TEMY = TEMY*T
          SUMX = TEMX*SUBX-TEMY*SUBY            !SUM=TEM*SUB
          SUMY = TEMX*SUBY+TEMY*SUBX
        END DO
        TEMX = 1.0-SUMX
        TEMY = -SUMY
        SUMX = TEMX*ZX-TEMY*ZY
        SUMY = TEMX*ZY+TEMY*ZX
        SUMX = -EULERX-CLNX+SUMX
        SUMY = -CLNY+SUMY
        EX = 0.0
        IF ( ZX.GT.-75.0 ) EX = EXP(ZX)
        EY = EX*SIN(ZY)
        EX = EX*COS(ZY)
        ANSX = SUMX*EX-SUMY*EY            !ZANS=SUM*EXP(Z)
        ANSY = SUMX*EY+SUMY*EX
      ELSE
        IF (AX .LT. 0.) A = SQRT((AX+29.)*(AX+29.)/9.0+AY*AY)
        ITER = 4 + NINT(128./A)
        TEMX = FLOAT(ITER)
        ADDX = TEMX
        ADDY = 0.0
        DO I=1,ITER-1
          TEMX = TEMX - 1.
          SUMX = ZX + ADDX                  !SUM = Z+ADD
          SUMY = ZY + ADDY
          X2PY2 = SUMX**2+SUMY**2
          RATX = TEMX*SUMX/X2PY2      !RAT = TEM/SUM
          RATY = -TEMX*SUMY/X2PY2
          RATX = 1.0+RATX
          X2PY2 = RATX**2+RATY**2
          ADDX = TEMX*RATX/X2PY2      !ADD = TEM/(1.+RAT))
          ADDY = -TEMX*RATY/X2PY2
        END DO
        ZX = ZX+ADDX                        !Z=A+ADD
        ZY = ZY+ADDY
        X2PY2 = ZX**2+ZY**2
        ANSX      = ZX/X2PY2                  !ZANS= 1./Z
        ANSY = -ZY/X2PY2
      END IF

      RETURN
      END
