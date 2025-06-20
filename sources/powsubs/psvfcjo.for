      SUBROUTINE PSVFCJO(DTT,TTHETA,SIG,GAM,SL,HL,PRFUNC,DPRDT,
     1  SIGPART,GAMPART,SLPART,HLPART)

!PURPOSE: Compute function & derivatives for Pseudovoigt profile
!   [W.I.F.David (1986), J. Appl. Cryst. 19, 63-64 &
!    P. Thompson, D.E. Cox & J.B. Hastings (1987) J. Appl. Cryst.,20,79-83.]
! Finger-Cox-Jephcoat (FCJ94) asymmetry correction 
!   [L.W. Finger, D.E. Cox & A.P. Jephcoat (1994) J. Appl. Cryst.,27,892-900.]
! coded 11/95 by B. H. Toby (NIST). revised version
! parameterized as asym1=S/L asym2=H/L


      INCLUDE       '../INCLDS/COPYRIGT.FOR' 

!CALLING ARGUMENTS:

      REAL*4        DTT                 !delta 2-theta in centidegrees                 
      REAL*4        TTHETA              !2-theta in centidegrees              
      REAL*4        SIG,GAM             
      REAL*4        SL,HL               !S/L & H/L               
      REAL*4        PRFUNC              
      REAL*4        DPRDT               
      REAL*4        SIGPART,GAMPART     
      REAL*4        SLPART,HLPART       

!INCLUDE STATEMENTS:
      REAL*4 SIND,COSD,TAND,ACOSD

!LOCAL VARIABLES:

      REAL*4        R                   ! pseudo-voight intensity
      REAL*4        DRDT                ! deriv R w/r theta
      REAL*4        DRDS                ! deriv R w/r sig
      REAL*4        DRDG                ! deriv R w/r gam
      REAL*4        F                   
      REAL*4        G                   
      REAL*4        DFDA                
      REAL*4        DGDA                
      REAL*4        DGDB                
      REAL*4        DYDA                
      REAL*4        DYDB                
      REAL*4        SIN2THETA2          ! sin(2theta)**2
      REAL*4        COS2THETA           ! cos(2theta)
      REAL*4        COS2THETA2          ! cos(2theta)**2
      REAL*4        SIN2THETA           ! sin(2THETA)
      REAL*4        SINDELTA            ! sin(Delta)
      REAL*4        COSDELTA            ! cos(Delta)
      REAL*4        RCOSDELTA           ! 1/cos(Delta)
      REAL*4        TANDELTA            ! tan(Delta)
      REAL*4        COSDELTA2           ! cos(Delta)**2
      REAL*4        A                   ! asym1 [coff(7)]
      REAL*4        B                   ! asym2 [coff(8)]
      REAL*4        APB                 ! (A+B)
      REAL*4        AMB                 ! (A-B)
      REAL*4        APB2                ! (A+B)**2
      REAL*4        TTHETAD             ! Two Theta in degrees

! Intermediate variables

      REAL*4        RSUMWG2             !      1.0/(sum of w G)**2
      REAL*4        SUMWG               !      sum of w G
      REAL*4        WG                  !      w G
      REAL*4        RSUMWG              !      1.0/sum of w G
      REAL*4        SUMWRG              !      sum of w G
      REAL*4        SUMWDGDA            !      sum of w dGdA
      REAL*4        SUMWRDGDA           !      sum of w R dGdA
      REAL*4        SUMWDGDB            !      sum of w dGdB
      REAL*4        SUMWRDGDB           !      sum of w R dGdB
      REAL*4        SUMWGDRD2T          !      sum of w G dRd(2theta)
      REAL*4        SUMWGDRDSIG         !      sum of w G dRdp(n)
      REAL*4        SUMWGDRDGAM         !      sum of w G dRdp(n)
      REAL*4        SUMWGDRDA           
      REAL*4        SUMWGDRDB           
      REAL*4        EMIN                ! 2phi minimum
      REAL*4        EINFL               ! 2phi of inflection point
      REAL*4        DEMINDA             ! Derivative of Emin wrt A
      REAL*4        DELTA               ! Angle of integration for convolution
      REAL*4        DDELTADA            
      REAL*4        TMP,TMP1,TMP2       ! intermediates
      INTEGER*4     I,K,IT              ! Miscellaneous loop variables
c       
c       Local Variables for Gaussian Integration
c        
      INTEGER*4     NGT                 !NUMBER OF TERMS IN GAUSSIAN QUADRATURE
      INTEGER*4     NUMTAB              ! NUMBER OF PRE-COMPUTED GAUSSIAN TABLES
      PARAMETER      (NUMTAB=34)
      INTEGER*4     NTERMS(NUMTAB)      ! NUMBER OF TERMS IN EACH TABLE - MUST BE EVEN
      INTEGER*4     FSTTERM(NUMTAB)     ! LOCATION OF 1ST TERM: N.B. N/2 TERMS
      LOGICAL*4     CALCLFG(NUMTAB)     ! TRUE IF TABLE HAS PREVIOUSLY BEEN CALCULATED
      INTEGER*4     ARRAYNUM            ! NUMBER OF SELECTED ARRAY
      INTEGER*4     ARRAYSIZE           ! SIZE OF COMPLETE ARRAY
      PARAMETER      (ARRAYSIZE=1670)
      REAL*4        XP(ARRAYSIZE)       !GAUSSIAN ABSCISSAS
      REAL*4        WP(ARRAYSIZE)       !GAUSSIAN WEIGHTS
      REAL*4        XPT(400)            !TEMPORARY GAUSSIAN ABSCISSAS
      REAL*4        WPT(400)            !TEMPORARY GAUSSIAN WEIGHTS
      REAL*4        STOFW               
      PARAMETER (STOFW=2.35482005)
      REAL*4        TODEG               
      PARAMETER (TODEG=57.2957795)
      SAVE      CALCLFG,XP,WP          !VALUES TO BE SAVED ACROSS CALLS

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA      NTERMS                ! NUMBER OF TERMS IN EACH TABLE - MUST BE EVEN
     1  /2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,40,50,60,70,80,90	! NUMBERS OF TERMS CHANGED MAY 2004
     2   ,100,110,120,140,160,180,200,220,240,260,280,300,400/
C       NOTE THAT NTERMS DETERMINES BOTH ARRAYSIZE AND FSTTERM
      DATA      FSTTERM                ! LOC. OF 1ST TERM: N.B. N/2 TERMS ARE SAVED
     1  /0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,140,165,195,230, !FSTERM(1) SHOULD BE 0 - INDEXING STARTS AT 1+THIS NUMBER!
     2  270,315,365,420,480,550,630,720,820,930,1050,1180,1320,1470/
      DATA      CALCLFG/NUMTAB*.FALSE./ ! TRUE IF TABLE ENTRY HAS BEEN CALCULATED 

!CODE:

c
C f(2theta) intermediates
c
      TTHETAD = 0.01 * TTHETA
      SIN2THETA = SIND(TTHETAD)
      COS2THETA = COSD(TTHETAD)
      SIN2THETA2 = SIN2THETA * SIN2THETA
      COS2THETA2 = COS2THETA * COS2THETA
c
C Asymmetry terms
c
      A = SL            ! A = S/L IN FCJ
      B = HL            ! B = H/L IN FCJ
      APB = A+B
      AMB = A-B
      APB2 = APB*APB
c
C handle the case where there is asymmetry
c
      IF (A .NE. 0.0 .OR. B .NE. 0.0) THEN
        EINFL = ACOSD(SQRT(1.0 + AMB**2)*COS2THETA) ! 2PHI(INFL) FCJ EQ 5 (DEGREES)
        TMP2 = 1.0 + APB2
        TMP = SQRT(TMP2)*COS2THETA
c
C Treat case where A or B is zero - Set Einfl = 2theta
c
        IF (A.EQ.0.0 .OR. B .EQ. 0.0)EINFL = ACOSD(COS2THETA)
        IF (ABS(TMP) .LE. 1.0) THEN
          EMIN = ACOSD(TMP)      ! 2PHI(MIN) FCJ EQ 4 (DEGREES)
          TMP1 = TMP2*(1.0 - TMP2*(1.0-SIN2THETA2))
        ELSE
          TMP1 = 0.0
          IF (TMP .GT. 0.0) THEN
            EMIN = 0.0
          ELSE
            EMIN = 180.0
          ENDIF
        ENDIF
        IF (TMP1 .GT. 0 .AND. ABS(TMP) .LE. 1.0) THEN
          DEMINDA = -APB*COS2THETA/SQRT(TMP1) ! N. B. SYMM W/R A,B
        ELSE
          DEMINDA = 0.0
        ENDIF

c       
c Compute Gaussian Quadrature interval 
C Note that the number of points must be large enough so that the
c interval between 2phi(min) and 2theta must be divided into steps no
c larger than 0.005 degrees.  LWF  8/10/95
c       
c Determine which Gauss-Legendre Table to use
c       
        ARRAYNUM = 1
c       
c Calculate desired number of intervals
c
        TMP = ABS(TTHETAD - EMIN)               ! TMP IN DEGREES
        IF ( GAM.LE.0.0 ) THEN
          K = INT(TMP*200.0)                    !RVD FORMULATION
        ELSE
          K = INT(300.0*TMP/GAM)                !NEW FORMULATION OF MAY 2004
        END IF
c
C Find the next largest set
c
        DO WHILE (ARRAYNUM .LT. NUMTAB .AND. K.GT.NTERMS(ARRAYNUM))
          ARRAYNUM = ARRAYNUM + 1
        ENDDO
        NGT = NTERMS(ARRAYNUM)
c
C calculate the terms, if they have not been used before
c
        IF (.NOT. CALCLFG(ARRAYNUM)) THEN
          CALCLFG(ARRAYNUM) = .TRUE.
          CALL GAULEG(-1.,1.,XPT,WPT,NGT)
          IT = FSTTERM(ARRAYNUM)-NGT/2
c
C copy the ngt/2 terms from our working array to the stored array
c
          DO K=NGT/2+1,NGT
            XP(K+IT) = XPT(K)
            WP(K+IT) = WPT(K)
          ENDDO
        ENDIF
        SUMWG = 0.
        SUMWRG = 0.
        SUMWDGDA = 0.
        SUMWRDGDA = 0.
        SUMWDGDB = 0.
        SUMWRDGDB = 0.
        SUMWGDRD2T = 0.
        SUMWGDRDSIG = 0.
        SUMWGDRDGAM = 0.
        SUMWGDRDA = 0.
        SUMWGDRDB = 0.
c       
c Compute Convolution integral for 2phi(min) <= delta <= 2theta
c       
        IT = FSTTERM(ARRAYNUM)-NGT/2
        DO K=NGT/2+1,NGT
          DELTA = EMIN + (TTHETAD - EMIN) * XP(K+IT) ! DELTA IN DEGREES
          CALL PSVOIGT(DTT+TTHETA-DELTA*100.,SIG,GAM,R,DRDT,DRDS,DRDG)

          DDELTADA = (1. - XP(K+IT))*DEMINDA ! N. B. SYMM W/R A,B
          SINDELTA = SIND(DELTA)
          COSDELTA = COSD(DELTA)
          RCOSDELTA = 1. / COSDELTA
          TANDELTA = TAND(DELTA)
          COSDELTA2 = COSDELTA*COSDELTA      
          TMP1 = COSDELTA2 - COS2THETA2
          TMP2 = SIN2THETA2 - SINDELTA * SINDELTA
          TMP = TMP2
          IF ( TTHETA.GT.4500.0 ) TMP = TMP1 
          IF (TMP .GT. 0) THEN
            TMP1 = 1.0/SQRT(TMP)
            F = ABS(COS2THETA) * TMP1
            DFDA = COSDELTA*COS2THETA*SINDELTA*DDELTADA * 
     1           (TMP1*TMP1*TMP1)
          ELSE
            F = 1.0
            DFDA = 0.0
          ENDIF
c       
c Calculate G(Delta,2theta) [G = W /(h cos(delta) ] [ FCJ eq. 7(a) and 7(b) ]
c       
          IF(ABS(DELTA-EMIN) .GT. ABS(EINFL-EMIN))THEN
            IF ( A.GE.B) THEN
c
C N.B. this is the only place where d()/dA <> d()/dB
c
              G = 2.0*B*F*RCOSDELTA
              DGDA = 2.0*B*RCOSDELTA*(DFDA + F*TANDELTA*DDELTADA)
              DGDB = DGDA + 2.0*F*RCOSDELTA
            ELSE
              G = 2.0*A*F*RCOSDELTA
              DGDB = 2.0*A*RCOSDELTA*(DFDA + F*TANDELTA*DDELTADA)
              DGDA = DGDB + 2.0*F*RCOSDELTA
            ENDIF
          ELSE                                            ! DELTA .LE. EINFL .OR. MIN(A,B) .EQ. 0
            G = (-1.0 + APB*F) * RCOSDELTA
            DGDA = RCOSDELTA*(F - TANDELTA*DDELTADA
     1             + APB*F*TANDELTA*DDELTADA + APB*DFDA)
            DGDB = DGDA
          ENDIF

          WG = WP(K+IT) * G
          SUMWG = SUMWG + WG
          SUMWRG = SUMWRG + WG * R
          SUMWDGDA = SUMWDGDA + WP(K+IT) * DGDA
          SUMWRDGDA = SUMWRDGDA + WP(K+IT) * R * DGDA
          SUMWDGDB = SUMWDGDB + WP(K+IT) * DGDB
          SUMWRDGDB = SUMWRDGDB + WP(K+IT) * R * DGDB
          SUMWGDRD2T = SUMWGDRD2T + WG * DRDT ! N.B. 1/CENTIDEGREES
          SUMWGDRDSIG = SUMWGDRDSIG + WG * DRDS
          SUMWGDRDGAM = SUMWGDRDGAM + WG * DRDG
          SUMWGDRDA = SUMWGDRDA + WG * DRDT * DDELTADA ! N. B. SYMM W/R A,B
        ENDDO
        RSUMWG = 1.0/SUMWG
        RSUMWG2 = RSUMWG * RSUMWG
        PRFUNC = SUMWRG * RSUMWG

        DYDA = (-(SUMWRG*SUMWDGDA) +
     1    SUMWG*(SUMWRDGDA- 100.0 * TODEG * SUMWGDRDA)) * RSUMWG2
        DYDB = (-(SUMWRG*SUMWDGDB) +
     1    SUMWG*(SUMWRDGDB- 100.0 * TODEG * SUMWGDRDA)) * RSUMWG2
        SIGPART = SUMWGDRDSIG * RSUMWG
        GAMPART = SUMWGDRDGAM * RSUMWG
        DPRDT = -SUMWGDRD2T * RSUMWG
      ELSE
c
C no asymmetry -- nice and simple!
c
        CALL PSVOIGT(DTT,SIG,GAM,R,DRDT,DRDS,DRDG)
        PRFUNC = R
        DYDA = 0.002 * SIGN(1.0,TTHETA - DTT)
        DYDB = DYDA
        SIGPART = DRDS
        GAMPART = DRDG
        DPRDT = -DRDT
      END IF
      SLPART = DYDA
      HLPART = DYDB
      
      RETURN
      END
      
