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
      real*4 sind,cosd,tand,acosd

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
      INTEGER*4     NGT                 !number of terms in Gaussian quadrature
      INTEGER*4     NUMTAB              ! number of pre-computed Gaussian tables
      parameter      (numtab=34)
      INTEGER*4     NTERMS(NUMTAB)      ! number of terms in each table - must be even
      INTEGER*4     FSTTERM(NUMTAB)     ! location of 1st term: N.B. N/2 terms
      LOGICAL*4     CALCLFG(NUMTAB)     ! true if table has previously been calculated
      INTEGER*4     ARRAYNUM            ! number of selected array
      INTEGER*4     ARRAYSIZE           ! size of complete array
      parameter      (arraysize=1670)
      REAL*4        XP(ARRAYSIZE)       !Gaussian abscissas
      REAL*4        WP(ARRAYSIZE)       !Gaussian weights
      REAL*4        XPT(400)            !temporary Gaussian abscissas
      REAL*4        WPT(400)            !temporary Gaussian weights
      REAL*4        STOFW               
      PARAMETER (STOFW=2.35482005)
      REAL*4        TODEG               
      PARAMETER (todeg=57.2957795)
      save      calclfg,xp,wp          !Values to be saved across calls

!FUNCTION DEFINITIONS:

!DATA STATEMENTS:

      DATA      NTERMS                ! number of terms in each table - must be even
     1  /2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,40,50,60,70,80,90	! Numbers of terms changed May 2004
     2   ,100,110,120,140,160,180,200,220,240,260,280,300,400/
C       note that nterms determines both arraysize and fstterm
      DATA      FSTTERM                ! loc. of 1st term: N.B. N/2 terms are saved
     1  /0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,140,165,195,230, !FSTERM(1) should be 0 - indexing starts at 1+this number!
     2  270,315,365,420,480,550,630,720,820,930,1050,1180,1320,1470/
      DATA      calclfg/numtab*.false./ ! true if table entry has been calculated 

!CODE:

c
C f(2theta) intermediates
c
      TTHETAD = 0.01 * TTHETA
      SIN2THETA = SIND(TTHETAD)
      cos2THETA = COSD(TTHETAD)
      sin2theta2 = sin2THETA * sin2THETA
      cos2THETA2 = cos2THETA * cos2THETA
c
C Asymmetry terms
c
      A = SL            ! A = S/L in FCJ
      B = HL            ! B = H/L in FCJ
      ApB = A+B
      AmB = A-B
      ApB2 = ApB*ApB
c
C handle the case where there is asymmetry
c
      IF (A .ne. 0.0 .or. B .ne. 0.0) then
        Einfl = Acosd(SQRT(1.0 + AmB**2)*cos2THETA) ! 2phi(infl) FCJ eq 5 (degrees)
        tmp2 = 1.0 + ApB2
        tmp = SQRT(tmp2)*cos2THETA
c
C Treat case where A or B is zero - Set Einfl = 2theta
c
        if (A.eq.0.0 .or. B .eq. 0.0)Einfl = Acosd(cos2THETA)
        if (abs(tmp) .le. 1.0) then
          Emin = Acosd(tmp)      ! 2phi(min) FCJ eq 4 (degrees)
          tmp1 = tmp2*(1.0 - tmp2*(1.0-sin2THETA2))
        else
          tmp1 = 0.0
          if (tmp .gt. 0.0) then
            Emin = 0.0
          else
            Emin = 180.0
          endif
        endif
        if (tmp1 .gt. 0 .and. abs(tmp) .le. 1.0) then
          dEmindA = -ApB*cos2THETA/SQRT(tmp1) ! N. B. symm w/r A,B
        ELSE
          dEmindA = 0.0
        ENDIF

c       
c Compute Gaussian Quadrature interval 
C Note that the number of points must be large enough so that the
c interval between 2phi(min) and 2theta must be divided into steps no
c larger than 0.005 degrees.  LWF  8/10/95
c       
c Determine which Gauss-Legendre Table to use
c       
        arraynum = 1
c       
c Calculate desired number of intervals
c
        tmp = abs(TTHETAD - emin)               ! tmp in degrees
        IF ( GAM.LE.0.0 ) THEN
          K = INT(TMP*200.0)                    !RVD formulation
        ELSE
          k = int(300.0*tmp/GAM)                !New formulation of May 2004
        END IF
c
C Find the next largest set
c
        do while (arraynum .lt. numtab .and. k.gt.nterms(arraynum))
          arraynum = arraynum + 1
        enddo
        ngt = nterms(arraynum)
c
C calculate the terms, if they have not been used before
c
        if (.not. calclfg(arraynum)) then
          calclfg(arraynum) = .true.
          call gauleg(-1.,1.,xpt,wpt,ngt)
          it = fstterm(arraynum)-ngt/2
c
C copy the ngt/2 terms from our working array to the stored array
c
          do k=ngt/2+1,ngt
            xp(k+it) = xpt(k)
            wp(k+it) = wpt(k)
          enddo
        endif
        sumWG = 0.
        sumWRG = 0.
        sumWdGdA = 0.
        sumWRdGdA = 0.
        sumWdGdB = 0.
        sumWRdGdB = 0.
        sumWGdRd2t = 0.
        sumWGdRdsig = 0.
        sumWGdRdgam = 0.
        sumWGdRdA = 0.
        sumWGdRdB = 0.
c       
c Compute Convolution integral for 2phi(min) <= delta <= 2theta
c       
        it = fstterm(arraynum)-ngt/2
        do k=ngt/2+1,ngt
          delta = emin + (TTHETAD - emin) * xp(k+it) ! Delta in degrees
          CALL PSVOIGT(DTT+TTHETA-delta*100.,SIG,GAM,R,dRdT,dRdS,dRdG)

          dDELTAdA = (1. - xp(k+it))*dEmindA ! N. B. symm w/r A,B
          sinDELTA = sind(Delta)
          cosDELTA = cosd(Delta)
          RcosDELTA = 1. / cosDELTA
          tanDELTA = tand(Delta)
          cosDELTA2 = cosDELTA*cosDELTA      
          tmp1 = cosDELTA2 - cos2THETA2
          tmp2 = sin2THETA2 - sinDelta * sinDELTA
          tmp = tmp2
          if ( ttheta.gt.4500.0 ) tmp = tmp1 
          if (tmp .gt. 0) then
            tmp1 = 1.0/SQRT(tmp)
            F = abs(cos2THETA) * tmp1
            dFdA = cosDELTA*cos2THETA*sinDELTA*dDELTAdA * 
     1           (tmp1*tmp1*tmp1)
          else
            F = 1.0
            dFdA = 0.0
          endif
c       
c Calculate G(Delta,2theta) [G = W /(h cos(delta) ] [ FCJ eq. 7(a) and 7(b) ]
c       
          if(abs(delta-emin) .gt. abs(einfl-emin))then
            if ( A.ge.B) then
c
C N.B. this is the only place where d()/dA <> d()/dB
c
              G = 2.0*B*F*RcosDELTA
              dGdA = 2.0*B*RcosDELTA*(dFdA + F*tanDELTA*dDELTAdA)
              dGdB = dGdA + 2.0*F*RcosDELTA
            else
              G = 2.0*A*F*RcosDELTA
              dGdB = 2.0*A*RcosDELTA*(dFdA + F*tanDELTA*dDELTAdA)
              dGdA = dGdB + 2.0*F*RcosDELTA
            endif
          else                                            ! delta .le. einfl .or. min(A,B) .eq. 0
            G = (-1.0 + ApB*F) * RcosDELTA
            dGdA = RcosDELTA*(F - tanDELTA*dDELTAdA
     1             + ApB*F*tanDELTA*dDELTAdA + ApB*dFdA)
            dGdB = dGdA
          endif

          WG = wp(k+it) * G
          sumWG = sumWG + WG
          sumWRG = sumWRG + WG * R
          sumWdGdA = sumWdGdA + wp(k+it) * dGdA
          sumWRdGdA = sumWRdGdA + wp(k+it) * R * dGdA
          sumWdGdB = sumWdGdB + wp(k+it) * dGdB
          sumWRdGdB = sumWRdGdB + wp(k+it) * R * dGdB
          sumWGdRd2t = sumWGdRd2t + WG * dRdT ! N.B. 1/centidegrees
          sumWGdRdsig = sumWGdRdsig + WG * dRdS
          sumWGdRdgam = sumWGdRdgam + WG * dRdG
          sumWGdRdA = sumWGdRdA + WG * dRdT * dDELTAdA ! N. B. symm w/r A,B
        enddo
        RsumWG = 1.0/sumWG
        RsumWG2 = RsumWG * RsumWG
        PRFUNC = sumWRG * RsumWG

        dydA = (-(sumWRG*sumWdGdA) +
     1    sumWG*(sumWRdGdA- 100.0 * todeg * sumWGdRdA)) * RsumWG2
        dydB = (-(sumWRG*sumWdGdB) +
     1    sumWG*(sumWRdGdB- 100.0 * todeg * sumWGdRdA)) * RsumWG2
        sigpart = sumWGdRdsig * RsumWG
        gampart = sumWGdRdgam * RsumWG
        DPRDT = -SumWGdRd2T * RsumWG
      else
c
C no asymmetry -- nice and simple!
c
        CALL PSVOIGT(DTT,SIG,GAM,R,dRdT,dRdS,dRdG)
        PRFUNC = R
        dydA = 0.002 * sign(1.0,TTHETA - DTT)
        dydB = dydA
        sigpart = dRdS
        gampart = dRdG
        DPRDT = -dRdT
      END IF
      SLPART = DYDA
      HLPART = DYDB
      
      RETURN
      END
      
