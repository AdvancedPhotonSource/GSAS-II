************************************************************************
*                                                                      *
*     Copyright 1987-1995 Michael M. J. Treacy and Michael W. Deem     *
*                                                                      *
************************************************************************
*******************      Source file DIFFaX.inc       ******************
************************************************************************
*                                                                      *
* This file contains the declaration of COMMON variables used          *
* throughout DIFFaX.                                                   *
*                                                                      *
* Some DIFFaX routines require static storage of variables. Not all    *
* compilers will automatically declare variables to be static (in fact *
* the fortran77 standard implies that non static allocation should be  *
* the default). To be safe, use a 'static' compiler option, if         *
* available. Alternatively, uncomment the 'save' statement at the      *
* beginning of the declaration section of this file, and also at the   *
* beginning of those subroutines and functions which do not include    *
* this file, such as RAN3. If variables are not 'saved', results will  *
* be unpredictable.                                                    *
*                                                                      *
************************************************************************
*************      Description of variables in COMMON       ************
************************************************************************
*                                                                      *
*            d-> indicates that the variable is read in                *
*                from the user's data file.                            *
*                                                                      *
*            i-> indicates that the variable is acquired               *
*                interactively at run_time, either by                  *
*                prompting the user for keyboard input, or             *
*                (if present) by reading the control file.             *
*                                                                      *
*            s-> indicates that the variable is read in                *
*                from the atomic scattering factor data file.          *
*                Normally, this data should never change.              *
*                                                                      *
************************************************************************
***********************     logical variables     **********************
************************************************************************
*    any_sharp    -  TRUE if DIFFaX suspects there are any sharp       *
*                    peaks.                                            *
*    all_Bs_zero  -  TRUE when all Bs_zero(i,j) = TRUE                 *
*    Bs_zero(MAX_L, MAX_L)                                             *
*                 -  TRUE if all layer stacking uncertainty factors in *
*                    one transition are all zero.                      *
*    CFile        -  TRUE if there is a file named 'control.dif'       *
*                    in the current directory. This is the so-called   *
*                    'control file' that automates the running of      *
*                    DIFFaX.                                           *
*    check_sym    -  TRUE if user specified a point group symmetry     *
*                    to test against. If the user's choice is          *
*                    incompatible with the input data, then this       *
*                    will be reset to FALSE.                           *
*    debug        -  TRUE for printouts to appear on screen (RVD)      *
*    DoDatdump    -  TRUE if the user wants a dump of the data file    *
*    DoSymDump    -  TRUE if the user wants to dump the output of      *
*                    the symmetry tests                                *
*d-> finite_width -  TRUE if layer widths are finite.                  *
*    has_l_mirror    TRUE if the diffraction data has a mirror         *
*                    perpendicular to c* (along fault axis)            *
*    h_mirror     -  TRUE if diffraction data has a mirror across      *
*                    the a*-c* plane                                   *
*    hk_mirror    -  TRUE if diffraction data has a mirror across      *
*                    the a*=b*, c* plane                               *
*    k_mirror     -  TRUE if diffraction data has a mirror across      *
*                    the b*-c* plane                                   *
*    intp_F       -  TRUE if the form factors are to be interpolated   *
*                    for the purposes of the Gauss-Legendre adaptive   *
*                    quadrature integration.                           *
*    inf_thick    -  TRUE if crystal is to be treated as if it has     *
*                    an infinite number of layers                      *
*    one_B(MAX_L) -  TRUE if all Debye-Waller factors are the          *
*                    same.                                             *
*    one_occup(MAX_L)                                                  *
*                 -  TRUE if all occupancy factors are the             *
*                    same.                                             *
*    only_real       TRUE if all layers are centrosymmetric            *
*d-> recrsv       -  TRUE if the recursive model is to be used to      *
*                    calculate diffraction from a statistical          *
*                    ensemble of crystallites.                         *
*d-> rndm         -  TRUE if 'xplcit' is TRUE, and if the user wishes  *
*                    the computer to generate a random sequence of     *
*                    layers biassed by the transition probabilities    *
*                    'l_alpha'.                                        *
*    rot_only     -  TRUE if diffraction point group symmetry          *
*                    has no vertical mirrors.                          *
*    same_Bs      -  TRUE if all layer stacking uncertainty            *
*                    factors are identical.                            *
*    same_layer   -  TRUE if all of the explicitly defined             *
*                    layers are identical.                             *
*    same_rz      -  TRUE if all stacking vectors have the same        *
*                    z-component.                                      *
*    there(MAX_L,MAX_L)                                                *
*                 -  TRUE if the transition j to i is non-zero         *
*d-> xplcit       -  TRUE if the user specifies a layer sequence       *
*                    explicitly.                                       *
*                    'recrsv' and 'xplcit' cannot both be TRUE         *
*    trim_origin  -  If TRUE, intensity near the origin will be        *
*                    ignored for the purposes of applying              *
*                    instrumental broadening.                          *
************************************************************************
**********************     integer*4 variables     *********************
************************************************************************
*d-> a_number(MAX_A,MAX_L)                                             *
*                 -  numeric label of atom in the layer.               *
*d-> a_type(MAX_A,MAX_L)                                               *
*                 -  type of each atom in each layer.                  *
*                    (a_type <= MAX_TA).                               *
*i-> bitdepth     -  The bit-depth of the selected area diffraction.   *
*                    pattern (sadp) data. Can equal 8 or 16.           *
*d-> blurring     -  Type of instrumental broadening to simulate.      *
*                    Choices are; PS_VGT, GAUSS, LORENZ                *
*    CENTRO       -  numeric constant (= 1) indicating layer has       *
*                    a center of symmetry.                             *
*    cntrl        -  The device number to read input from.             *
*                    If the file 'control.dif' is present, this file   *
*                    becomes the default input. Otherwise, it is the   *
*                    keyboard.                                         *
*s-> e_sf(MAX_TA) -  Electron scattering factors.                      *
*    ELECTN       -  numeric constant (= 2) indicating radiation       *
*                    type is electrons.                                *
*i-> full_brd     -  1 if full adaptive integration is to be carried   *
*                    out on the sharp spots.                           *
*i-> full_shrp    -  1 if full adaptive integration is to be carried   *
*                    out on the streaks.                               *
*    GAUSS        -  numeric constant (= 1) indicating user wishes     *
*                    to simulate Gaussian instrumental                 *
*                    broadening, with a constant half width            *
*    h_bnd        -  Maximum value of h to explore.                    *
*    k_bnd        -  Maximum value of k to explore.                    *
*    l_actual(MAX_L)                                                   *
*                 -  Contains the layer number that layer i is         *
*                    structurally identical to. If all layers are      *
*                    unique, l_actual(i) = i;                          *
*                    else, l_actual(i) <= i.                           *
*    l_cnt        -  Number of layers in explicit sequence. This is    *
*                    tallied by DIFFaX, by counting the layers.        *
*    l_n_atoms(MAX_L)                                                  *
*                 -  number of atoms in each layer                     *
*d-> l_symmetry(MAX_L)                                                 *
*                 -  symmetry of layer (NONE or CENTRO)                *
*d-> l_seq(XP_MAX)                                                     *
*                 -  array containing the explicitly defined           *
*                    sequence of layers. Used only if                  *
*                    'xplcit' = TRUE                                   *
*i-> loglin       -  1 if logarithmic scaling of SADP data is          *
*                     required                                         *
*    LORENZ       -  numeric constant (= 2) indicating user wishes     *
*                    to simulate Lorentzian instrumental               *
*                    broadening, with a constant half width            *
*    maxsad       -  Maximum intensity of sadp patterns.               *
*    n_actual     -  Number of unique layers ( <= n_layers).           *
*    n_atoms      -  Temporary variable holding the number of unique   *
*                    atoms in a given layer.                           *
*d-> n_layers     -  Number of user-defined layers. n_layers <= MAX_L  *
*    NEUTRN       -  Numeric constant (= 1) indicating radiation       *
*                    type is neutrons.                                 *
*    no_trials    -  Number of reciprocal space points to sample.      *
*                    whilst evaluating diffraction symmetry.           *
*    NONE         -  numeric constant (= 0) indicating layer has       *
*                    no center of symmetry.                            *
*    PS_VGT       -  numeric constant (= 3) indicating user wishes     *
*                    to simulate pseudo-Voigt instrumental             *
*                    broadening.                                       *
*    PV_LRN       -  numeric constant (= 5) indicating user wishes     *
*                    to simulate Lorentzian instrumental               *
*                    broadening, with a variable half width            *
*    PV_GSS       -  numeric constant (= 4) indicating user wishes     *
*                    to simulate Gaussian instrumental                 *
*                    broadening, with a variable half width            *
*d-> rad_type     -  Type of radiation for which to calculate          *
*                    diffraction intensities. Choices are;             *
*                    X_RAY, NEUTRN, ELECTN                             *
*    sadblock     -  Length of a row of SADP data                      *
*d-> SymGrpNo     -  Numeric descriptor of diffraction symmetry.       *
*    X_RAY        -  numeric constant (= 0) indicating radiation       *
*                    type is X-rays.                                   *
*                                                                      *
************************************************************************
**********************     character variables     *********************
************************************************************************
*d-> a_name(MAX_A,MAX_L)*4                                             *
*                 -  Name of each atom. DIFFaX expects 4               *
*                    characters. See file 'data.sfc' for               *
*                    allowed names.                                    *
*    atom_l(MAX_TA)                                                    *
*                 -  The name of each type of atom found in the        *
*                    structure data file                               *
*    cfname*16    -  The name of the control file containing the       *
*                    automated responses to DIFFaX's prompts           *
*                    (usually set to 'control.dif').                   *
*d-> pnt_grp*12   -  Symbolic name of the point group symmetry of      *
*                    the diffraction data.                             *
*    sfname*16    -  The name of the data file containing the          *
*                    atomic scattering factor data (usually set to     *
*                    'data.sfc')                                       *
*                                                                      *
************************************************************************
***********************     real*8 variables     ***********************
************************************************************************
*    a0           -  One of seven reciprocal lattice constants         *
*d-> a_B(MAX_A,MAX_L)                                                  *
*                 -  isotropic Debye-Waller factor for each            *
*                    atom in each layer.                               *
*d-> a_occup(MAX_A,MAX_L)                                              *
*                 -  Occupancy of each atom in each layer.             *
*                    (Normally this will lie between 0 and 1).         *
*d-> a_pos(3,MAX_A,MAX_L)                                              *
*                 -  x,y,z relative coordinates of each atom           *
*                    in each layer.                                    *
*    a_B11,a_B22,a_B33,a_B12,a_B23,a_B31                               *
*                 -  The average values of the r_Bij arrays            *
*    ab0          -  One of seven reciprocal lattice constants         *
*    b0           -  One of seven reciprocal lattice constants         *
*    bc0          -  One of seven reciprocal lattice constants         *
*    bnds_wt      -  Equals 1.0 if rot_only is TRUE, otherwise         *
*                    equals 0.5 if rot_only is FALSE (ie there is      *
*                    a vertical mirror                                 *
*    brd_spc(MAX_SP)                                                   *
*                 -  Array holding the powder diffraction data         *
*                    after instrumental broadening has been            *
*                    applied.                                          *
*    c0           -  One of seven reciprocal lattice constants         *
*    ca0          -  One of seven reciprocal lattice constants         *
*d-> cell_a       -  Unit cell a axis dimension.                       *
*d-> cell_b       -  Unit cell b axis dimension.                       *
*d-> cell_c       -  Unit cell c axis dimension.                       *
*d-> cell_gamma   -  Angle between a and b axes. Angle between b and   *
*                    c, and a and c axes is 90 degrees by default.     *
*    d0           -  One of seven reciprocal lattice constants         *
*d-> d_theta      -  Angular increment in PXD spectrum.                *
*    DEG2RAD      -  Conversion factor for degrees to radians          *
*    detune(MAX_L,MAX_L)                                               *
*                 -  Array of small positive numbers whose             *
*                    purpose is to prevent the determinant of the      *
*                    recursion array 'mat' from becoming zero at the   *
*                    sharp peaks. This produces a singularity which    *
*                    is hard to integrate accurately over. In essence, *
*                    the 'detune' parameters are small stacking        *
*                    uncertainty factors. The result is to reduce      *
*                    the value of l_alpha(j,i) by an amount            *
*                    detune(j,i), such that the sum over the alphas    *
*                    for stacking from a given layer do not quite      *
*                    add to unity.                                     *
*    fatsWalla_hk -  temporary storage for Fats-Waller factor          *
*    ffact_scale  -  Angular scale (radians) of array 'formfactor'.    *
*    ffhkcnst     -  Constant associated with the form-factor half-    *
*                    width. Depends on reflection indices h and k, as  *
*                    well as Wx and Wy.                                *
*    ffwdth       -  Form factor half width in reciprocal Angstroms.   *
*    formfactor(FFACT_SIZE)                                            *
*                 -  Array containing a normalized Lorentzian profile, *
*                    the form factor due to in-plane size broadening.  *
*                    The profile is Lorentzian out to N_SIGMAS half-   *
*                    widths, and linear to zero from there. The linear *
*                    portion has the same gradient as the last point   *
*                    of the Lorentzian portion, thus the sampling step *
*                    is governed by N_SIGMAS as well as FFACT_SIZE.    *
*    h_end        -  (h_end,k_end) is a vector in the                  *
*                    h-k plane of reciprocal space defining the        *
*                    upper boundary of the wedge in reciprocal space   *
*                    to integrate within. This is defined by the       *
*                    symmetry of the diffraction data.                 *
*    high_atom(MAX_L)                                                  *
*                 -  The highest atomic z-rel position in each layer   *
*    h_start      -  (h_start,k_start) is a vector in the              *
*                    h-k plane of reciprocal space defining the        *
*                    lower boundary of the wedge in reciprocal space   *
*                    to integrate within. This is defined by the       *
*                    symmetry of the diffraction data.                 *
*d-> FWHM         -  Full width half maximum of instrumental           *
*                    broadening.                                       *
*    hx_ky(MAX_A,MAX_L)                                                *
*                 -  Temporary storage of h*Rx + k*Ry, whilst          *
*                    l*Rz is being computed along the streaks.         *
*    k_end        -  (h_end,k_end) is a vector in the                  *
*                    h-k plane of reciprocal space defining the        *
*                    upper boundary of the wedge in reciprocal space   *
*                    to integrate within. This is defined by the       *
*                    symmetry of the diffraction data.                 *
*    k_start      -  (h_start,k_start) is a vector in the              *
*                    h-k plane of reciprocal space defining the        *
*                    lower boundary of the wedge in reciprocal space   *
*                    to integrate within. This is defined by the       *
*                    symmetry of the diffraction data.                 *
*d-> l_alpha(MAX_L,MAX_L)                                              *
*                 -  Array of layer transition probabilities.          *
*                    The order is (column, row).                       *
*    l_bnd        -  Maximum value of l to explore.                    *
*    l_g(MAX_L)   -  Array of layer existence probabilities.           *
*                    These are determined by the transition            *
*                    probabilities 'l_alpha' entered by the user.      *
*    low_atom(MAX_L)                                                   *
*                 -  The lowest atomic z-rel position in each layer    *
*d-> l_r(3,MAX_L,MAX_L)                                                *
*                 -  Array of layer stacking vectors.                  *
*                    The order is (column, row).                       *
*    l_rz         -  Value of Rz if same_rz = TRUE.                    *
*d-> lambda       -  Radiation wavelength.                             *
*    max_angle    -  Maximum angle that intensity information is to    *
*                    be taken from for the purposes of evaluating      *
*                    diffraction symmetry.                             *
*    max_var      -  Maximum mean variation of intensities when a      *
*                    given symmetry operator was applied.              *
*    mltplcty     -  1/mltplcty is the fraction of reciprocal space    *
*                    necessary to integrate over, as determined by     *
*                    the diffraction point group symmetry.             *
*s-> n_sf(MAX_TA) -  Neutron scattering factors.                       *
*    PI           -  The value of pi, 3.141592653589793.....           *
*    PI2          -  The value of 2*pi                                 *
*d-> pv_gamma     -  Pseudo-Voigt gamma parameter.                     *
*d-> pv_u         -  Pseudo-Voigt u parameter.                         *
*d-> pv_v         -  Pseudo-Voigt v parameter.                         *
*d-> pv_w         -  Pseudo-Voigt w parameter.                         *
*d-> r_B11(MAX_L,MAX_L)                                                *
*d-> r_B22(MAX_L,MAX_L)                                                *
*d-> r_B33(MAX_L,MAX_L)                                                *
*d-> r_B12(MAX_L,MAX_L)                                                *
*d-> r_B23(MAX_L,MAX_L)                                                *
*d-> r_B31(MAX_L,MAX_L)                                                *
*                 -  The 6 components of the anisotropic layer         *
*                    stacking uncertainties. These are                 *
*                    equivalent to the atomic Debye-Waller             *
*                    factors, except they apply to the stacking        *
*                    vectors. These parameters allow for               *
*                    'turbostratic' disorder such as is found          *
*                    in liquid crystals. These parameters are          *
*                    optional, and can be entered by the user          *
*                    in the PARAMETERS section enclosed in             *
*                    parentheses.                                      *
*    RAD2DEG      -  Conversion factor for radians to degrees          *
*    scaleint     -  Intensity scaling factor used in calculating      *
*                    the selected area diffraction patterns.           *
*    spec(MAX_SP) -  Array holding the unbroadened powder              *
*                    diffraction data. This array also holds           *
*                    the SADP image data.                              *
*d-> th2_max      -  Upper bound of angle in PXD spectrum.             *
*d-> th2_min      -  Lower bound of angle in PXD spectrum.             *
*    theta1       -  angle relative to (1,0,0) of lower wedge bound    *
*    theta2       -  angle relative to (1,0,0) of upper wedge bound    *
*d-> tolerance    -  Maximum deviation that intensities can have       *
*                    from symmetry related points if intensities are   *
*                    to be considered equal.                           *
*    tiny_inty    -  a small intensity value used in the diffraction   *
*                    symmetry checking routines. Intensities lower     *
*                    than tiny_inty are treated as being close to zero.*
*                 -  X-ray scattering factors.                         *
*s-> x_sf(9,MAX_TA)                                                    *
*                 -  X-ray scattering factors.                         *
*i-> Wa           -  In-plane width of crystal along a-direction.      *
*i-> Wb           -  In-plane width of crystal perpendicular to        *
*                    a-direction. Wx and Wy in Angstroms.              *
*                                                                      *
************************************************************************
*********************     complex*16 variables     *********************
************************************************************************
*    l_phi(MAX_L,MAX_L)                                                *
*                 -  Phases of components of 'mat'                     *
*    mat(MAX_L,MAX_L)                                                  *
*                 -  Recursion matrix relating the scattering          *
*                    from crystals centered on different layers        *
*    mat1(MAX_L,MAX_L)                                                 *
*                 -  Storage for intermediate 'mat' results.           *
*    wavefn       -  Coherent wavefunction calculated for an           *
*                    explicitly defined sequence of layers (if         *
*                    requested)                                        *
*                                                                      *
************************************************************************
***************      Declaration of COMMON variables      **************
************************************************************************
*
*     save
*
      character*4 a_name(MAX_A,MAX_L), atom_l(MAX_TA)
      character*12 pnt_grp
*      character*16 sfname, cfname
*
      logical one_B(MAX_L), one_occup(MAX_L), Bs_zero(MAX_L, MAX_L),
     |        there(MAX_L,MAX_L)
      logical only_real, same_Bs, all_Bs_zero, rot_only, CFile,
     |        DoDatdump, DoSymDump, intp_F, trim_origin, recrsv,
     |        xplcit, rndm, inf_thick, has_l_mirror, h_mirror,
     |        k_mirror, hk_mirror, check_sym, same_rz, any_sharp,
     |        same_layer, finite_width,debug
*
      integer*4 l_seq(XP_MAX), pow(MAX_BIN), a_type(MAX_A,MAX_L),
     |          l_n_atoms(MAX_L), l_symmetry(MAX_L), l_actual(MAX_L),
     |          a_number(MAX_A,MAX_L), e_sf(MAX_TA)
      integer*4 SymGrpNo, no_trials, h_bnd, k_bnd, cntrl, max_pow,
     |          l_cnt,full_brd, full_shrp, sadblock, loglin, bitdepth,
     |          rad_type, n_layers, n_actual, blurring, n_atoms, maxsad
      integer*4 NONE, CENTRO, GAUSS, LORENZ, PS_VGT, PV_GSS, PV_LRN,
     |          X_RAY, NEUTRN, ELECTN
*
      real*8 l_alpha(MAX_L,MAX_L), l_r(3,MAX_L,MAX_L), l_g(MAX_L),
     |       a_pos(3,MAX_A,MAX_L), a_B(MAX_A,MAX_L),
     |       a_occup(MAX_A,MAX_L), high_atom(MAX_L), low_atom(MAX_L),
     |       r_B11(MAX_L,MAX_L),r_B22(MAX_L,MAX_L),r_B33(MAX_L,MAX_L),
     |       r_B12(MAX_L,MAX_L),r_B23(MAX_L,MAX_L),r_B31(MAX_L,MAX_L),
     |       hx_ky(MAX_A,MAX_L), spec(MAX_SP), brd_spc(MAX_SP),
     |       detune(MAX_L,MAX_L), x_sf(9,MAX_TA)
      real*8 a_B11,a_B22,a_B33,a_B12,a_B23,a_B31
      real*8 tolerance, max_var, max_angle, l_bnd, l_rz,
     |       PI, PI2, RAD2DEG, DEG2RAD, scaleint, brightness, lambda,
     |       th2_min, th2_max, d_theta, h_start, k_start,h_end, k_end,
     |       cell_a, cell_b, cell_c, cell_gamma, pv_u, pv_v, pv_w,
     |       pv_gamma, FWHM, mltplcty, bnds_wt, theta1, theta2, a0, b0,
     |       c0, d0, ab0, bc0, ca0, tiny_inty, fatsWalla_hk
      real*8 formfactor(FFACT_SIZE),ffact_scale,Wa,Wb,ffhkcnst,ffwdth
      real*8 n_sf(MAX_TA)
*
      complex*16 mat(MAX_L,MAX_L), mat1(MAX_L,MAX_L),
     |           l_phi(MAX_L,MAX_L), wavefn
*
      common /chars1/ a_name, atom_l
      common /chars2/ pnt_grp
*      common /chars3/ sfname, cfname
*
      common /logic1/ one_B, one_occup, Bs_zero, there
      common /logic2/ only_real, same_Bs, all_Bs_zero, rot_only, CFile,
     |                DoDatdump, DoSymDump, intp_F, trim_origin,
     |                recrsv, xplcit, rndm, inf_thick, has_l_mirror,
     |                h_mirror, k_mirror, hk_mirror, check_sym,
     |                same_rz, any_sharp, same_layer, finite_width,
     |                debug
*
      common /integ1/ l_seq, pow, a_type, l_n_atoms, l_symmetry,
     |                l_actual, a_number, e_sf
      common /integ2/ SymGrpNo, no_trials, h_bnd, k_bnd, cntrl,
     |                max_pow, l_cnt, full_brd, full_shrp, sadblock,
     |                loglin, bitdepth, rad_type, n_layers, n_actual,
     |                blurring, n_atoms, maxsad
*
*The common block /consts/ is initialized in the block data section
      common /consts/ NONE, CENTRO, GAUSS, LORENZ, PS_VGT, PV_GSS,
     |                PV_LRN, X_RAY, NEUTRN, ELECTN
*
      common /reals1/ l_alpha, l_r, l_g, a_pos, a_B,a_occup, high_atom,
     |                low_atom, r_B11,r_B22,r_B33, r_B12,r_B23,r_B31,
     |                hx_ky,spec, brd_spc, detune, x_sf
      common /reals2/ a_B11,a_B22,a_B33,a_B12,a_B23,a_B31
      common /reals3/ tolerance, max_var, max_angle, l_bnd, l_rz,
     |                PI, PI2, RAD2DEG, DEG2RAD, scaleint,
     |                brightness, lambda, th2_min, th2_max, d_theta,
     |                h_start, k_start,h_end, k_end, cell_a, cell_b,
     |                cell_c, cell_gamma, pv_u, pv_v, pv_w, pv_gamma,
     |                FWHM, mltplcty, bnds_wt, theta1, theta2, a0, b0,
     |                c0, d0, ab0, bc0, ca0, tiny_inty, fatsWalla_hk
      common /reals4/ formfactor, ffact_scale, Wa, Wb, ffhkcnst, ffwdth
*
      common /cmplx1/ mat, mat1, l_phi, wavefn
*
      equivalence (n_sf, x_sf)
*
