      SUBROUTINE QLMN(L,MM,NN,Q)

!PURPOSE: Compute Ql,m,n for spherical harmonics from lookup table

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:         

      INTEGER*4     L,MM,NN             !order & subindices (m may be <0)
      REAL*4        Q                   !Output value

!INCLUDE STATEMENTS:

      REAL*4        QT                  
      COMMON /QLMNVAL/QT(2109)


!LOCAL VARIABLES:

      REAL*8        SUM,TEMP,TEMP1      
      INTEGER*4     LMN,I,J,M,N         

!FUNCTION DEFINITIONS:                     

      REAL*8        FACTLN              !Compute ln-factorial & binominal coeffs.

!DATA STATEMENTS:

!CODE:  
                          
      M = ABS(MM)            
      N = ABS(NN)
      IF ( N.GT.M ) THEN
        I = M
        M = N
        N = I
      END IF
      IF ( MOD(N,2).EQ.0 .AND. MOD(L,2).EQ.0 ) THEN         !Even L,N - do lookup
        J = 0
        DO I=2,L,2
          J = J+(I/2)**2      !points to last term in L-2 block
        END DO
        J = J+1
        DO I=0,M-1
          J = J+(I+2)/2       !points to 1st term in M block
        END DO
        J = J+N/2             !offset to N term
!       PRINT '(A,I4,F12.8)',' J,Q ',J,QT(J)
        Q = QT(J)                          
      ELSE                       !Odd L or N - calculate Q
        LMN = L-M-N
        TEMP = 0.5D0*(FACTLN(L+N)+FACTLN(L+M)+
     1    FACTLN(L-M)+FACTLN(L-N))
        SUM = 0.0D0
        DO I=0,L-N
          IF ( (L-M-I).GE.0 .AND. (M+N+I).GE.0 ) THEN
            TEMP1 = TEMP-FACTLN(I)-FACTLN(L-M-I)-
     1        FACTLN(L-N-I)-FACTLN(M+N+I)
            TEMP1 = DEXP(TEMP1)
            IF ( MOD(I,2).NE.0 ) TEMP1 = -TEMP1
            SUM = SUM+TEMP1
          END IF
        END DO
        Q = SUM/2.**(1.*L)
        IF ( MOD(LMN,2).NE.0 ) Q = -Q
!       PRINT '(A,3I4,F12.8)',' l,m,n,Q(lmn)',L,M,N,Q      
      END IF
      IF ( MM.LT.0 .AND. MOD(L+NN,2).NE.0 ) Q = -Q
      IF ( NN.LT.0 .AND. MOD(L+MM,2).NE.0 ) Q = -Q
      RETURN
      END

