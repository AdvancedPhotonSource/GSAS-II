      SUBROUTINE QLMNINIT

!PURPOSE: Compute Ql,m,n for spherical harmonics up to L=34 - only does even orders 
!     and only even N terms - by R.I. Sheldon & modified by R. Von Dreele for GSAS

      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:         

!INCLUDE STATEMENTS:
                             
      REAL*4        QT                  
      COMMON /QLMNVAL/QT(2109)

!LOCAL VARIABLES:

      REAL*8        SUM,TEMP,TEMP1      
      INTEGER*4     I,LMN,M,MM          

!FUNCTION DEFINITIONS:                     

      REAL*8        FACTLN              !Compute ln-factorial & binominal coeffs.

!DATA STATEMENTS:

!CODE:  
      
      J = 1
      QT(1) = 1.0      
      DO L=2,34,2
        DO M=0,L
          DO N=0,M,2
            J = J+1
            LMN = L-M-N
            TEMP = 0.5D0*(FACTLN(L+N)+FACTLN(L+M)+
     1        FACTLN(L-M)+FACTLN(L-N))
            SUM = 0.0D0
            DO I=0,L-N
              IF ( (L-M-I).GE.0 .AND. (M+N+I).GE.0 ) THEN
                TEMP1 = TEMP-FACTLN(I)-FACTLN(L-M-I)-
     1            FACTLN(L-N-I)-FACTLN(M+N+I)
                TEMP1 = DEXP(TEMP1)
                IF ( MOD(I,2).NE.0 ) TEMP1 = -TEMP1
                SUM = SUM+TEMP1
              END IF
            END DO
            QT(J) = SUM/2.**(1.*L)
            IF ( MOD(LMN,2).NE.0 ) QT(J) = -QT(J)
!            PRINT '(A,3I4,F12.8)',' l,m,n,Q(lmn)',L,M,N,QT(J)
          END DO
        END DO
      END DO
      RETURN
      END

