 Module Profile_Finger
    !    Here the functions provided in the original code by L. Finger have been
    !    implemented as subroutines as is the common practice in Fortran 90 when
    !    several dummy arguments have output intent attribute.
    !    The subroutine Prof_Val returns the value of the profile function at twoth of
    !    a peak of centre twoth0 as well as the derivatives wrt profile parameters.
    !    Asymmetry due to axial divergence using the method of Finger, Cox and Jephcoat,
    !    J. Appl. Cryst. 27, 892, 1992.
    !    This version based on code provided by Finger, Cox and Jephcoat as modified
    !    by J Hester (J. Appl. Cryst. (2013),46,1219-1220)
    !    with a new derivative calculation method and then optimised,
    !    adapted to Fortran 90 and further improved by J. Rodriguez-Carvajal.
    !    Further changed by J. Hester to include optimisations found in GSAS-II version
    !    of original FCJ code, and streamlined for incorporation into GSAS-II
   
    implicit none

    !---- List of Public Subroutines ----!
    public :: prof_val

    !---- List of Private Functions ----!
    private ::  dfunc_int, extra_int

    integer,       parameter :: sp = selected_real_kind(6,30)
    integer,       parameter :: dp = selected_real_kind(14,150)
    integer,       parameter :: cp = sp    !switch to double precision by putting cp=dp
    real(kind=dp), parameter :: pi = 3.141592653589793238463_dp
    real(kind=dp), parameter :: to_DEG  = 180.0_dp/pi
    real(kind=dp), parameter :: to_RAD  = pi/180.0_dp
    integer,       parameter :: numtab = 34
    integer, private,dimension(numtab) :: nterms =(/2,4,6,8,10,12,14,16,18,20,22,24,26,28, &
         30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,240,260,280,300,400/)
    integer, private,dimension(numtab) :: fstterm=(/0,1,3,6,10,15,21,28,36,45,55,66,78, &
         91,105,120,140,165,195,230,270,315,365,420,480,550,630, &
         720,820,930,1050,1180,1320,1470/) !fstterm(1) should be zero, indexing starts at 1+ this number
    real(kind=cp),private,dimension(0:1883) :: wp  !extra space for expansion
    real(kind=cp),private,dimension(0:1883) :: xp  !extra space for expansion

    !
    ! Variables to switch to new calculations of some variables that depend
    ! only on (twoth0,asym1,asym2) and not on the particular point of the profile.
    ! When the subroutine is invoked for the same reflection some variables
    ! are not calculated.
    !
    !
    real(kind=cp), private, save :: twoth0_prev = 0.0_cp
    real(kind=cp), private, save ::  asym1_prev = 0.0_cp
    real(kind=cp), private, save ::  asym2_prev = 0.0_cp

    ! Fixed constants
    real(kind=cp), private, parameter :: pi_over_two=0.5_cp*pi
    real(kind=cp), private, parameter :: eps=1.0e-6_cp
    ! Following gives the notional number of steps per degree for quadrature
    integer,       private, parameter :: ctrl_nsteps= 300
    logical,       private, dimension(numtab) :: calclfg
    !calclfg true if table entry calculated
    data calclfg/numtab*.false./
 Contains

    Subroutine Prof_Val( SIG, gamma, asym1, asym2, twoth, twoth0, dprdt, dprdg,  &
         dprdlz , dprds , dprdd , profval)
      real(kind=cp), Intent(In)    :: sig        ! Gaussian variance (i.e. sig squared)
      real(kind=cp), Intent(In)    :: gamma      ! Lorenzian FWHM (i.e. 2 HWHM)
      real(kind=cp), Intent(In)    :: asym1      ! D_L+S_L
      real(kind=cp), Intent(In)    :: asym2      ! D_L-S_L
      real(kind=cp), Intent(In)    :: twoth      ! point at which to evaluate the profile
      real(kind=cp), Intent(In)    :: twoth0     ! two_theta value for peak
      real(kind=cp), Intent(Out)   :: dprdt      ! derivative of profile wrt TwoTH0
      real(kind=cp), Intent(Out)   :: dprdg      ! derivative of profile wrt Gaussian sig
      real(kind=cp), Intent(Out)   :: dprdlz     ! derivative of profile wrt Lorenzian
      real(kind=cp), Intent(Out)   :: dprds      ! derivative of profile wrt asym1
      real(kind=cp), Intent(Out)   :: dprdd      ! derivative of profile wrt asym2
      real(kind=cp), Intent(Out)   :: profval    ! Value of the profile at point twoth
      !The variables below have the "save" attribute in order to save calculation
      !time when the subroutine is invoked for different points of the same peak
      real(kind=cp),save :: s_l , d_l, half_over_dl
      real(kind=cp),save :: df_dh_factor, df_ds_factor
      real(kind=cp),save :: dfi_emin, dfi_einfl
      real(kind=cp),save :: normv_analytic
      real(kind=cp),save :: einfl              ! 2phi value for inflection point
      real(kind=cp),save :: emin               ! 2phi value for minimum
      real(kind=cp),save :: cstwoth            ! cos(twoth)
      real(kind=cp),save :: coseinfl           ! cos(Einfl)
      real(kind=cp),save :: apb                ! (S + H)/L
      real(kind=cp),save :: amb                ! (S - H)/L
      real(kind=cp),save :: apb2               ! (ApB) **2
      Integer,save       :: arraynum, ngt, ngt2, it
      Logical,save       :: s_eq_d
      

      ! Variables not conserving their value between calls
      Integer       :: side, k
      real(kind=cp) :: tmp , tmp1 , tmp2  ! intermediate values
      real(kind=cp) :: delta              ! Angle of integration for comvolution
      real(kind=cp) :: sindelta           ! sine of DELTA
      real(kind=cp) :: cosdelta           ! cosine of DELTA
      real(kind=cp) :: rcosdelta          ! 1/cos(DELTA)
      real(kind=cp) :: f,g, einflr,eminr,twoth0r
      real(kind=cp) :: sumwg, sumwrg, sumwrdgda ,sumwdgdb , sumwrdgdb
      real(kind=cp) :: sumwgdrdg, sumwgdrdlz, sumwgdrd2t
      real(kind=cp) :: sumwx
      real(kind=cp) :: xpt(1000)          !temporary storage
      real(kind=cp) :: wpt(1000)          !temporary storage
      logical       :: re_calculate

      ! First simple calculation of Pseudo-Voigt if asymmetry is not used
      if(asym1 == 0.0) then
        call Psvoigt(twoth-twoth0,sig,gamma,tmp,tmp1,dprdg,dprdlz)
        profval = tmp
        dprds = 0.4*sign(1.0,2.0*twoth0-twoth)
        dprdd = 0.0
        dprdt = -tmp1    !derivative relative to centre position
        return
      end if

      !From here to the end of the procedure asymmetry is used.
      !Make the calculations of some variables only if twoth0,asym1,asym2
      !are different from previous values. This saves calculation time if the
      !different points of a peak are calculated sequentially for the same values
      !of twotheta and asymmetry parameters.

      re_calculate= abs(twoth0_prev-twoth0) > eps .or.  &
                    abs(asym1_prev-asym1)   > eps .or.  &
                    abs(asym2_prev-asym2)   > eps

      if(re_calculate) then
        twoth0_prev=twoth0
         asym1_prev=asym1
         asym2_prev=asym2

        twoth0r=twoth0*to_rad
        cstwoth = Cos(twoth0r)
        s_l = 0.5*(asym1 - asym2)  ! 1/2(s_l+d_l - (d_l-s_l))
        d_l = 0.5*(asym1 + asym2)  ! 1/2(s_l+d_l + (d_l-s_l))
        apb = asym1
        amb = asym2
        ! Catch special case of S_L = D_L
        If (Abs(amb) < 0.00001) Then
          s_eq_d = .TRUE.
        Else
          s_eq_d = .FALSE.
        End If
        apb2 = apb*apb

        tmp = Sqrt(1.0 + amb*amb)*cstwoth
        If ((Abs(tmp) > 1.0) .or. (Abs(tmp) <= Abs(cstwoth))) Then
          einfl = twoth0
          einflr=einfl*to_rad
          dfi_einfl = pi_over_two
        Else
          einflr = Acos(tmp)
          einfl=einflr*to_deg
          dfi_einfl = dfunc_int(einflr,twoth0r)
        End If
        coseinfl = Cos(einflr)
        tmp2 = 1.0 + apb2
        tmp = Sqrt(tmp2) * cstwoth

        ! If S_L or D_L are zero, set Einfl = 2theta
        ! If S_L equals D_L, set Einfl = 2theta

        If ((s_l == 0.0) .OR. (d_l == 0.0) .OR. s_eq_d) then
          einfl = twoth0
          einflr=einfl*to_rad
        End if

        If (Abs(tmp) <= 1.0) Then
          eminr = Acos(tmp)
          emin = eminr * to_deg
          tmp1 = tmp2 * (1.0 - tmp2 * cstwoth*cstwoth)
        Else
          tmp1 = 0.0
          If (tmp > 0.0) Then
            emin = 0.0
            eminr= 0.0
          Else
            emin = 180.0
            eminr= pi
          End If
        End If

        dfi_emin = dfunc_int(eminr,twoth0r)
        !
        ! Simplifications if S_L equals D_L
        !
        half_over_dl=0.5_cp/d_l
        If (s_eq_d) Then
          dfi_einfl = pi_over_two
          normv_analytic = (pi_over_two - dfi_emin)  &
              - 2.0_cp*half_over_dl*(extra_int(einflr)-extra_int(eminr))
          df_dh_factor =  half_over_dl * (pi_over_two - dfi_emin)
          df_ds_factor =  half_over_dl * (pi_over_two - dfi_emin)
          df_dh_factor = df_dh_factor - 2.0_cp*half_over_dl * normv_analytic
        Else
          dfi_einfl = dfunc_int(einflr,twoth0r)
          normv_analytic = Min(s_l,d_l)/d_l*(pi_over_two - dfi_einfl)
          normv_analytic = normv_analytic + apb*half_over_dl*(dfi_einfl-dfi_emin)   &
                       - 2.0_cp*half_over_dl*(extra_int(einflr)-extra_int(eminr))
          tmp= half_over_dl*(pi - dfi_einfl - dfi_emin)
          tmp1=half_over_dl*(dfi_einfl - dfi_emin)
          If(d_l < s_l) Then
            df_dh_factor = tmp
            df_ds_factor = tmp1
          Else
            df_dh_factor = tmp1
            df_ds_factor = tmp
          End If
          df_dh_factor = df_dh_factor - 2.0_cp*half_over_dl * normv_analytic
        End If
        arraynum = 1
        ! Number of terms needed, GSAS-II formulation
        tmp = abs(twoth0 - emin)
        if (gamma <= 0.0) then
           k = ctrl_nsteps*tmp/2
        else
           k = ctrl_nsteps*tmp/(100*gamma)
        endif
        Do
           if ( .not. ( arraynum < numtab  .And.  k > nterms(arraynum) ) ) exit
           arraynum = arraynum + 1
        End Do
        ngt = nterms(arraynum)              ! Save the number of terms
        ngt2 = ngt / 2
        ! Calculate Gauss-Legendre quadrature terms the first time they
        ! are required
        if (.not. calclfg(arraynum)) then
           calclfg(arraynum) = .true.
           call gauleg(-1.,1.,xpt,wpt,ngt)
           it = fstterm(arraynum)-ngt2
           !
           ! copy the ngt/2 terms from our working array to the stored array
           !
           do k=ngt2+1,ngt
              xp(k+it) = xpt(k)
              wp(k+it) = wpt(k)
           enddo
      end if
      it = fstterm(arraynum)-ngt2  !in case skipped initialisation

   End if   !re_calculate
        ! Clear terms needed for summations
      sumwg = 0.0
      sumwrg = 0.0
      sumwrdgda = 0.0
      sumwdgdb = 0.0
      sumwrdgdb = 0.0
      sumwgdrd2t = 0.0
      sumwgdrdg = 0.0
      sumwgdrdlz = 0.0
      sumwx = 0.0
      ! Compute the convolution integral with the pseudovoight.
      ! using Gauss-Legendre quadrature.
      ! In theory, we should use the weighted w_i values of the
      ! product of the asymmetry
      ! profile with the pseudovoight at a set of points x_i in the
      ! interval [-1,1]. However, 
      ! the following adopts the GSAS-II approach of instead integrating
      ! between 0 and 1, which can be imagined as extending the range
      ! of integration be from Emin - twoth0 up to twoth0.
      ! This seems to be preferable because the most important area
      ! for integration is the rapidly increasing peak, and so the
      ! higher density of points at the end of the interval works in our
      ! favour.
      Do k = ngt2+1 , ngt
          delta = emin + (twoth0 - emin) * xp(k + it)
          sindelta = Sin(delta*to_rad)
          cosdelta = Cos(delta*to_rad)
          If (Abs(cosdelta) < 1.0E-15) cosdelta = 1.0E-15
          rcosdelta = Abs(1.0 / cosdelta)
          tmp = cosdelta*cosdelta - cstwoth*cstwoth
          If (tmp > 0.0) Then
            tmp1 = Sqrt(tmp)
            f = Abs(cstwoth) / tmp1           !h-function in FCJ
          Else
            f = 0.0
          End If
          !  calculate G(Delta,2theta) , FCJ eq. 7a and 7b
          If ( Abs(delta - emin) > Abs(einfl - emin)) Then
            If (s_l > d_l) Then
              g = 2.0 * d_l * f * rcosdelta
            Else
              g = 2.0 * s_l * f * rcosdelta
            End If
          Else
            g = (-1.0 + apb * f) * rcosdelta
          End If
          call Psvoigt(twoth-delta,sig,gamma,tmp,dprdt,dprdg,dprdlz)
          sumwg = sumwg + wp(k+it) * g
          sumwrg = sumwrg + wp(k+it) * g * tmp
          If ( Abs(cosdelta) > Abs(coseinfl)) Then
            sumwrdgda = sumwrdgda + wp(k+it) * f * rcosdelta * tmp
            sumwrdgdb = sumwrdgdb + wp(k+it) * f * rcosdelta * tmp
          Else
            If (s_l < d_l) Then
              sumwrdgdb = sumwrdgdb + 2.0*wp(k+it)*f* rcosdelta*tmp
            Else
              sumwrdgda = sumwrdgda + 2.0*wp(k+it)*f* rcosdelta*tmp
            End If
          End If
          sumwgdrd2t = sumwgdrd2t + wp(k+it) * g * dprdt
          sumwgdrdg = sumwgdrdg + wp(k+it) * g * dprdg
          sumwgdrdlz = sumwgdrdlz + wp(k+it) * g * dprdlz
      End Do

      If (sumwg == 0.0) sumwg = 1.0_cp
      profval = sumwrg / sumwg
      ! Minus sign in following as Psvoight returns derivs against x, not
      ! against the centre position.
      dprdt = -sumwgdrd2t/ sumwg
      dprdg = sumwgdrdg / sumwg
      dprdlz = sumwgdrdlz / sumwg
      !
      If(normv_analytic <= 0.0) normv_analytic=1.0_cp
      dprdd = sumwrdgda / sumwg - df_dh_factor*profval/normv_analytic - profval/d_l
      dprds = sumwrdgdb / sumwg - df_ds_factor*profval/normv_analytic

      dprds = 0.5_cp*(dprdd + dprds)  !S is really D+S
      dprdd = 0.5_cp*(dprdd - dprds)  !D is really D-S
      Return
    End Subroutine Prof_Val

!  Function to give the analytical value of the normalisation constant

    Function dfunc_int(twopsi, twoth0) result(dfunc)
      Real(kind=cp), Intent(In)  :: twopsi
      Real(kind=cp), Intent(In)  :: twoth0
      Real(kind=cp)              :: dfunc
      !--- Local variables
      Real(kind=cp) :: sintp        !Sin twopsi
      Real(kind=cp) :: sin2t,sin2t2,csp,csm,ssp,ssm,a,b ! sin2Theta, (sin2Theta)^2

      If(Abs(twopsi-twoth0) < 1.0E-5) Then
        dfunc=pi_over_two
      Else
        sin2t=Sin(twoth0)
        sin2t2=sin2t*sin2t
        sintp = Sin(twopsi)
        csp=sintp+sin2t2
        csm=sintp-sin2t2
        ssp=Abs((sintp+1.0_cp)*sin2t)
        ssm=Abs((sintp-1.0_cp)*sin2t)
        a=csm/ssm; b=-csp/ssp
        If(a > 1.0_cp) a=1.0_cp
        If(b <-1.0_cp) b=-1.0_cp
        dfunc=0.5_cp*(Asin(a)-Asin(b))
      End If
    End Function dfunc_int

    !  Function to calculate 1/4(log(|sin(x)+1|)-log(|sin(x)-1|))
    Function extra_int(x) result(extra)
      Real(kind=cp), Intent(In) :: x
      Real(kind=cp)             :: extra
      !--- Local variables
      Real(kind=cp)             :: sinx

      sinx = Sin(x)
      extra = 0.25_cp*(Log(Abs(sinx+1.0_cp))-Log(Abs(sinx-1.0_cp)))
    End Function extra_int

End Module Profile_Finger

! Linkage for F2Py as direct calling of Fortran module contents
! from outside file containing module fails when separate object files combined into
! single Python loadable module
Subroutine get_prof_val( SIG, gamma, asym1, asym2, twoth, twoth0, dprdt, dprdg,  &
     dprdlz , dprds , dprdd , profval)
  use Profile_Finger, only:Prof_Val
      integer,       parameter :: sp = selected_real_kind(6,30)
      integer,       parameter :: dp = selected_real_kind(14,150)
      integer,       parameter :: cp = sp    !switch to double precision by putting cp=dp
      real(kind=cp), Intent(In)    :: sig        ! Gaussian variance (i.e. sig squared)
      real(kind=cp), Intent(In)    :: gamma      ! Lorenzian FWHM (i.e. 2 HWHM)
      real(kind=cp), Intent(In)    :: asym1      ! D_L+S_L
      real(kind=cp), Intent(In)    :: asym2      ! D_L-S_L
      real(kind=cp), Intent(In)    :: twoth      ! point at which to evaluate the profile
      real(kind=cp), Intent(In)    :: twoth0     ! two_theta value for peak
      real(kind=cp), Intent(Out)   :: dprdt      ! derivative of profile wrt TwoTH0
      real(kind=cp), Intent(Out)   :: dprdg      ! derivative of profile wrt Gaussian sig
      real(kind=cp), Intent(Out)   :: dprdlz     ! derivative of profile wrt Lorenzian
      real(kind=cp), Intent(Out)   :: dprds      ! derivative of profile wrt asym1
      real(kind=cp), Intent(Out)   :: dprdd      ! derivative of profile wrt asym2
      real(kind=cp), Intent(Out)   :: profval    ! Value of the profile at point twoth
      Call Prof_Val(SIG,GAMMA,asym1,asym2,twoth,twoth0,DPRDT,DPRDG,DPRDLZ,DPRDS, &
           DPRDD,profval)
      RETURN
    END subroutine get_prof_val


