      SUBROUTINE PSVFCJEXPO(DTT,TTHETA,ALP,BET,SIG,GAM,SL,HL,PRFUNC,
     1 DPRDT,ALPART,BEPART,SIGPART,GAMPART,SLPART,HLPART)

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
      REAL*4        ALP,BET,SIG,GAM             
      REAL*4        SL,HL               !S/L & H/L               
      REAL*4        PRFUNC              
      REAL*4        DPRDT               
      REAL*4        ALPART,BEPART,SIGPART,GAMPART     
      REAL*4        SLPART,HLPART       

!INCLUDE STATEMENTS:
      REAL*4 SIND,COSD,TAND,ACOSD

!LOCAL VARIABLES:

      REAL*4        R                   ! pseudo-voight intensity
      REAL*4        DRDT                ! deriv R w/r theta
      REAL*4        DRDA                ! deriv R w/r alpha
      REAL*4        DRDB                ! deriv R w/r beta
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
      REAL*4        SIN2THETA           ! sin(2THETA)
      REAL*4        SINDELTA            ! sin(Delta)
      REAL*4        COSDELTA            ! cos(Delta)
      REAL*4        RCOSDELTA           ! 1/cos(Delta)
      REAL*4        TANDELTA            ! tan(Delta)
      REAL*4        COSDELTA2           ! cos(Delta)**2
      REAL*4        A                   ! asym1 [coff(9)]
      REAL*4        B                   ! asym2 [coff(10)]
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
      REAL*4        SUMWGDRDALP         !      sum of w G dRdp(n)
      REAL*4        SUMWGDRDBET         !      sum of w G dRdp(n)
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
      INTEGER*4     NGT                 !number of terms in Gaussian quadrature
      INTEGER*4     NUMTAB              ! number of pre-computed Gaussian tables
      PARAMETER      (NUMTAB=34)
      INTEGER*4     NTERMS(NUMTAB)      ! number of terms in each table - must be even
      INTEGER*4     FSTTERM(NUMTAB)     ! location of 1st term: N.B. N/2 terms
      LOGICAL*4     CALCLFG(NUMTAB)     ! true if table has previously been calculated
      INTEGER*4     ARRAYNUM            ! number of selected array
      INTEGER*4     ARRAYSIZE           ! size of complete array
      PARAMETER      (ARRAYSIZE=1670)
      REAL*4        XP(ARRAYSIZE)       !Gaussian abscissas
      REAL*4        WP(ARRAYSIZE)       !Gaussian weights
      REAL*4        XPT(400)            !temporary Gaussian abscissas
      REAL*4        WPT(400)            !temporary Gaussian weights
      REAL*4        STOFW               
      PARAMETER (STOFW=2.35482005)
      REAL*4        TODEG               
      PARAMETER (todeg=57.2957795)
      SAVE      CALCLFG,XP,WP          !VALUES TO BE SAVED ACROSS CALLS

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA      NTERMS                ! number of terms in each table - must be even
     1  /2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,40,50,60,70,80,90	! Numbers of terms changed May 2004
     2   ,100,110,120,140,160,180,200,220,240,260,280,300,400/
C       note that nterms determines both arraysize and fstterm
      DATA      FSTTERM                ! loc. of 1st term: N.B. N/2 terms are saved
     1  /0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,140,165,195,230, !FSTERM(1) should be 0 - indexing starts at 1+this number!
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
C
C ASYMMETRY TERMS
C
      A = SL            ! A = S/L IN FCJ
      B = HL            ! B = H/L IN FCJ
      APB = A+B
      AMB = A-B
      APB2 = APB*APB
C
C HANDLE THE CASE WHERE THERE IS ASYMMETRY
C
      IF (A .NE. 0.0 .OR. B .NE. 0.0) THEN
        EINFL = ACOSD(SQRT(1.0 + AMB**2)*COS2THETA) ! 2PHI(INFL) FCJ EQ 5 (DEGREES)
        TMP2 = 1.0 + APB2
        TMP = SQRT(TMP2)*COS2THETA
C
C TREAT CASE WHERE A OR B IS ZERO - SET EINFL = 2THETA
C
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

C       
C COMPUTE GAUSSIAN QUADRATURE INTERVAL 
C NOTE THAT THE NUMBER OF POINTS MUST BE LARGE ENOUGH SO THAT THE
C INTERVAL BETWEEN 2PHI(MIN) AND 2THETA MUST BE DIVIDED INTO STEPS NO
C LARGER THAN 0.005 DEGREES.  LWF  8/10/95
C       
C DETERMINE WHICH GAUSS-LEGENDRE TABLE TO USE
C       
        ARRAYNUM = 1
C       
C CALCULATE DESIRED NUMBER OF INTERVALS
C
        TMP = ABS(TTHETAD - EMIN)               ! TMP IN DEGREES
        IF ( GAM.LE.0.0 ) THEN
          K = INT(TMP*200.0)                    !RVD FORMULATION
        ELSE
          K = INT(300.0*TMP/GAM)                !NEW FORMULATION OF MAY 2004
        END IF
C
C FIND THE NEXT LARGEST SET
C
        DO WHILE (ARRAYNUM .LT. NUMTAB .AND. K.GT.NTERMS(ARRAYNUM))
          ARRAYNUM = ARRAYNUM + 1
        ENDDO
        NGT = NTERMS(ARRAYNUM)
C
C CALCULATE THE TERMS, IF THEY HAVE NOT BEEN USED BEFORE
C
        IF (.NOT. CALCLFG(ARRAYNUM)) THEN
          CALCLFG(ARRAYNUM) = .TRUE.
          CALL GAULEG(-1.,1.,XPT,WPT,NGT)
          IT = FSTTERM(ARRAYNUM)-NGT/2
C
C COPY THE NGT/2 TERMS FROM OUR WORKING ARRAY TO THE STORED ARRAY
C
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
        SUMWGDRDALP = 0.
        SUMWGDRDBET = 0.
        SUMWGDRDSIG = 0.
        SUMWGDRDGAM = 0.
        SUMWGDRDA = 0.
        SUMWGDRDB = 0.
C       
C COMPUTE CONVOLUTION INTEGRAL FOR 2PHI(MIN) <= DELTA <= 2THETA
C       
        IT = FSTTERM(ARRAYNUM)-NGT/2
        DO K=NGT/2+1,NGT
          DELTA = EMIN + (TTHETAD - EMIN) * XP(K+IT) ! DELTA IN DEGREES
          CALL EPSVOIGT(DTT+TTHETA-DELTA*100.,ALP,BET,SIG,GAM,R,DRDT,
     1      DRDA,DRDB,DRDS,DRDG)

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
C       
C CALCULATE G(DELTA,2THETA) [G = W /(H COS(DELTA) ] [ FCJ EQ. 7(A) AND 7(B) ]
C       
          IF(ABS(DELTA-EMIN) .GT. ABS(EINFL-EMIN))THEN
            IF ( A.GE.B) THEN
C
C N.B. THIS IS THE ONLY PLACE WHERE D()/DA <> D()/DB
C
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
          SUMWGDRDALP = SUMWGDRDALP + WG * DRDA
          SUMWGDRDBET = SUMWGDRDBET + WG * DRDB
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
        ALPART = SUMWGDRDALP * RSUMWG
        BEPART = SUMWGDRDBET * RSUMWG
        SIGPART = SUMWGDRDSIG * RSUMWG
        GAMPART = SUMWGDRDGAM * RSUMWG
        DPRDT = SUMWGDRD2T * RSUMWG
      ELSE
C
C NO ASYMMETRY -- NICE AND SIMPLE! (never invoked!)
C
        CALL EPSVOIGT(DTT,ALP,BET,SIG,GAM,R,DRDT,DRDA,DRDB,DRDS,DRDG)
        
        PRFUNC = R
        DYDA = 0.002 * SIGN(1.0,TTHETA - DTT)
        DYDB = DYDA
        ALPART = DRDA
        BEPART = DRDB
        SIGPART = DRDS
        GAMPART = DRDG
        DPRDT = DRDT
      END IF
      SLPART = DYDA
      HLPART = DYDB
      RETURN
      END
      
