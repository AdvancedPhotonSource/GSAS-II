.. 
    This file is created using the makeVarTbl.py. Edit that not this file.
         
.. list-table:: Naming for GSAS-II parameter names, ``p:h:<var>:n``
   :widths: 35 65
   :header-rows: 1

   * - ``<var>``
     - usage
   * - \ :math:`\scriptstyle K`\  (example: ``a``)
     - Lattice parameter, \ :math:`\scriptstyle K`\ , from Ai and Djk; where \ :math:`\scriptstyle K`\  is one of the characters a, b or c.
   * - α
     - Lattice parameter, α, computed from both Ai and Djk.
   * - β
     - Lattice parameter, β, computed from both Ai and Djk.
   * - γ
     - Lattice parameter, γ, computed from both Ai and Djk.
   * - Scale
     - Phase fraction (as p:h:Scale) or Histogram scale factor (as :h:Scale).
   * - A\ :math:`\scriptstyle I`\  (example: ``A0``)
     - Reciprocal metric tensor component \ :math:`\scriptstyle I`\ ; where \ :math:`\scriptstyle I`\  is a digit between 0 and 5.
   * - \ :math:`\scriptstyle L`\ ol (example: ``vol``)
     - Unit cell volume; where \ :math:`\scriptstyle L`\  is one of the characters v or V.
   * - dA\ :math:`\scriptstyle M`\  (example: ``dAx``)
     - Refined change to atomic coordinate, \ :math:`\scriptstyle M`\ ; where \ :math:`\scriptstyle M`\  is one of the characters x, y or z.
   * - A\ :math:`\scriptstyle M`\  (example: ``Ax``)
     - Fractional atomic coordinate, \ :math:`\scriptstyle M`\ ; where \ :math:`\scriptstyle M`\  is one of the characters x, y or z.
   * - AUiso
     - Atomic isotropic displacement parameter.
   * - AU\ :math:`\scriptstyle N_0`\ \ :math:`\scriptstyle N_1`\  (example: ``AU11``)
     - Atomic anisotropic displacement parameter U\ :math:`\scriptstyle N_0`\ \ :math:`\scriptstyle N_1`\ ; where \ :math:`\scriptstyle N_0`\  is one of the characters 1, 2 or 3 and \ :math:`\scriptstyle N_1`\  is one of the characters 1, 2 or 3.
   * - Afrac
     - Atomic site fraction parameter.
   * - Amul
     - Atomic site multiplicity value.
   * - AM\ :math:`\scriptstyle M`\  (example: ``AMx``)
     - Atomic magnetic moment parameter, \ :math:`\scriptstyle M`\ ; where \ :math:`\scriptstyle M`\  is one of the characters x, y or z.
   * - Back\ :math:`\scriptstyle J`\  (example: ``Back11``)
     - Background term #\ :math:`\scriptstyle J`\ ; where \ :math:`\scriptstyle J`\  is the background term number.
   * - BkPkint;\ :math:`\scriptstyle J`\  (example: ``BkPkint;11``)
     - Background peak #\ :math:`\scriptstyle J`\  intensity; where \ :math:`\scriptstyle J`\  is the background peak number.
   * - BkPkpos;\ :math:`\scriptstyle J`\  (example: ``BkPkpos;11``)
     - Background peak #\ :math:`\scriptstyle J`\  position; where \ :math:`\scriptstyle J`\  is the background peak number.
   * - BkPksig;\ :math:`\scriptstyle J`\  (example: ``BkPksig;11``)
     - Background peak #\ :math:`\scriptstyle J`\  Gaussian width; where \ :math:`\scriptstyle J`\  is the background peak number.
   * - BkPkgam;\ :math:`\scriptstyle J`\  (example: ``BkPkgam;11``)
     - Background peak #\ :math:`\scriptstyle J`\  Cauchy width; where \ :math:`\scriptstyle J`\  is the background peak number.
   * - BF mult
     - Background file multiplier.
   * - Bab\ :math:`\scriptstyle O`\  (example: ``BabA``)
     - Babinet solvent scattering coef. \ :math:`\scriptstyle O`\ ; where \ :math:`\scriptstyle O`\  is one of the characters A or U.
   * - D\ :math:`\scriptstyle N_0`\ \ :math:`\scriptstyle N_1`\  (example: ``D11``)
     - Anisotropic strain coef. \ :math:`\scriptstyle N_0`\ \ :math:`\scriptstyle N_1`\ ; where \ :math:`\scriptstyle N_0`\  is one of the characters 1, 2 or 3 and \ :math:`\scriptstyle N_1`\  is one of the characters 1, 2 or 3.
   * - Extinction
     - Extinction coef.
   * - MD
     - March-Dollase coef.
   * - Mustrain;\ :math:`\scriptstyle J`\  (example: ``Mustrain;11``)
     - Microstrain coefficient (delta Q/Q x 10**6); where \ :math:`\scriptstyle J`\  can be i for isotropic or equatorial and a is axial (uniaxial broadening), a number for generalized (Stephens) broadening or mx for the Gaussian/Lorentzian mixing term (LGmix).
   * - Size;\ :math:`\scriptstyle J`\  (example: ``Size;11``)
     - Crystallite size value (in microns); where \ :math:`\scriptstyle J`\  can be i for isotropic or equatorial, and a is axial (uniaxial broadening), a number between 0 and 5 for ellipsoidal broadening or mx for the Gaussian/Lorentzian mixing term (LGmix).
   * - eA
     - Cubic mustrain value.
   * - Ep
     - Primary extinction.
   * - Es
     - Secondary type II extinction.
   * - Eg
     - Secondary type I extinction.
   * - Flack
     - Flack parameter.
   * - TwinFr
     - Twin fraction.
   * - Layer Disp
     - Layer displacement along beam.
   * - Absorption
     - Absorption coef.
   * - LayerDisp
     - Bragg-Brentano Layer displacement.
   * - Displace\ :math:`\scriptstyle P`\  (example: ``DisplaceX``)
     - Debye-Scherrer sample displacement \ :math:`\scriptstyle P`\ ; where \ :math:`\scriptstyle P`\  is one of the characters X or Y.
   * - Lam
     - Wavelength.
   * - I(L2)\\/I(L1)
     - Ka2/Ka1 intensity ratio.
   * - Polariz.
     - Polarization correction.
   * - SH/L
     - FCJ peak asymmetry correction.
   * - \ :math:`\scriptstyle Q`\  (example: ``U``)
     - Gaussian instrument broadening \ :math:`\scriptstyle Q`\ ; where \ :math:`\scriptstyle Q`\  is one of the characters U, V or W.
   * - \ :math:`\scriptstyle R`\  (example: ``X``)
     - Cauchy instrument broadening \ :math:`\scriptstyle R`\ ; where \ :math:`\scriptstyle R`\  is one of the characters X, Y or Z.
   * - Zero
     - Debye-Scherrer zero correction.
   * - Shift
     - Bragg-Brentano sample displ.
   * - SurfRoughA
     - Bragg-Brenano surface roughness A.
   * - SurfRoughB
     - Bragg-Brenano surface roughness B.
   * - Transparency
     - Bragg-Brentano sample tranparency.
   * - DebyeA
     - Debye model amplitude.
   * - DebyeR
     - Debye model radius.
   * - DebyeU
     - Debye model Uiso.
   * - RBV\ :math:`\scriptstyle J`\  (example: ``RBV11``)
     - Vector rigid body parameter.
   * - RBVO\ :math:`\scriptstyle S`\  (example: ``RBVOa``)
     - Vector rigid body orientation parameter \ :math:`\scriptstyle S`\ ; where \ :math:`\scriptstyle S`\  is one of the characters a, i, j or k.
   * - RBVP\ :math:`\scriptstyle M`\  (example: ``RBVPx``)
     - Vector rigid body \ :math:`\scriptstyle M`\  position parameter; where \ :math:`\scriptstyle M`\  is one of the characters x, y or z.
   * - RBVf
     - Vector rigid body site fraction.
   * - RBV\ :math:`\scriptstyle T_0`\ \ :math:`\scriptstyle U_0`\ \ :math:`\scriptstyle U_1`\  (example: ``RBVT11``)
     - Residue rigid body group disp. param.; where \ :math:`\scriptstyle T_0`\  is one of the characters T, L or S and \ :math:`\scriptstyle U_0`\  is one of the characters 1, 2, 3, A or B and \ :math:`\scriptstyle U_1`\  is one of the characters 1, 2, 3, A or B.
   * - RBVU
     - Residue rigid body group Uiso param.
   * - RBRO\ :math:`\scriptstyle S`\  (example: ``RBROa``)
     - Residue rigid body orientation parameter \ :math:`\scriptstyle S`\ ; where \ :math:`\scriptstyle S`\  is one of the characters a, i, j or k.
   * - RBRP\ :math:`\scriptstyle M`\  (example: ``RBRPx``)
     - Residue rigid body \ :math:`\scriptstyle M`\  position parameter; where \ :math:`\scriptstyle M`\  is one of the characters x, y or z.
   * - RBRTr;\ :math:`\scriptstyle J`\  (example: ``RBRTr;11``)
     - Residue rigid body torsion parameter.
   * - RBRf
     - Residue rigid body site fraction.
   * - RBR\ :math:`\scriptstyle T_0`\ \ :math:`\scriptstyle U_0`\ \ :math:`\scriptstyle U_1`\  (example: ``RBRT11``)
     - Residue rigid body group disp. param.; where \ :math:`\scriptstyle T_0`\  is one of the characters T, L or S and \ :math:`\scriptstyle U_0`\  is one of the characters 1, 2, 3, A or B and \ :math:`\scriptstyle U_1`\  is one of the characters 1, 2, 3, A or B.
   * - RBRU
     - Residue rigid body group Uiso param.
   * - RBSAtNo
     - Atom number for spinning rigid body.
   * - RBSO\ :math:`\scriptstyle S`\  (example: ``RBSOa``)
     - Spinning rigid body orientation parameter \ :math:`\scriptstyle S`\ ; where \ :math:`\scriptstyle S`\  is one of the characters a, i, j or k.
   * - RBSP\ :math:`\scriptstyle M`\  (example: ``RBSPx``)
     - Spinning rigid body \ :math:`\scriptstyle M`\  position parameter; where \ :math:`\scriptstyle M`\  is one of the characters x, y or z.
   * - RBSShC\ :math:`\scriptstyle V`\  (example: ``RBSShC1``)
     - Spinning rigid body sph. harmonics term; where \ :math:`\scriptstyle V`\  is one of the characters 1, -, 2 or 0 ,, 1, -, 2 or 0.
   * - constr\ :math:`\scriptstyle G`\  (example: ``constr10``)
     - Generated degree of freedom from constraint; where \ :math:`\scriptstyle G`\  is one or more digits (0, 1,... 9).
   * - nv-(.+)
     - New variable assignment with name \1.
   * - mV\ :math:`\scriptstyle H`\  (example: ``mV0``)
     - Modulation vector component \ :math:`\scriptstyle H`\ ; where \ :math:`\scriptstyle H`\  is the digits 0, 1, or 2.
   * - Fsin
     - Sin site fraction modulation.
   * - Fcos
     - Cos site fraction modulation.
   * - Fzero
     - Crenel function offset.
   * - Fwid
     - Crenel function width.
   * - Tmin
     - ZigZag/Block min location.
   * - Tmax
     - ZigZag/Block max location.
   * - \ :math:`\scriptstyle R`\ max (example: ``Xmax``)
     - ZigZag/Block max value for \ :math:`\scriptstyle R`\ ; where \ :math:`\scriptstyle R`\  is one of the characters X, Y or Z.
   * - \ :math:`\scriptstyle R`\ sin (example: ``Xsin``)
     - Sin position wave for \ :math:`\scriptstyle R`\ ; where \ :math:`\scriptstyle R`\  is one of the characters X, Y or Z.
   * - \ :math:`\scriptstyle R`\ cos (example: ``Xcos``)
     - Cos position wave for \ :math:`\scriptstyle R`\ ; where \ :math:`\scriptstyle R`\  is one of the characters X, Y or Z.
   * - U\ :math:`\scriptstyle N_0`\ \ :math:`\scriptstyle N_1`\ sin (example: ``U11sin``)
     - Sin thermal wave for U\ :math:`\scriptstyle N_0`\ \ :math:`\scriptstyle N_1`\ ; where \ :math:`\scriptstyle N_0`\  is one of the characters 1, 2 or 3 and \ :math:`\scriptstyle N_1`\  is one of the characters 1, 2 or 3.
   * - U\ :math:`\scriptstyle N_0`\ \ :math:`\scriptstyle N_1`\ cos (example: ``U11cos``)
     - Cos thermal wave for U\ :math:`\scriptstyle N_0`\ \ :math:`\scriptstyle N_1`\ ; where \ :math:`\scriptstyle N_0`\  is one of the characters 1, 2 or 3 and \ :math:`\scriptstyle N_1`\  is one of the characters 1, 2 or 3.
   * - M\ :math:`\scriptstyle R`\ sin (example: ``MXsin``)
     - Sin mag. moment wave for \ :math:`\scriptstyle R`\ ; where \ :math:`\scriptstyle R`\  is one of the characters X, Y or Z.
   * - M\ :math:`\scriptstyle R`\ cos (example: ``MXcos``)
     - Cos mag. moment wave for \ :math:`\scriptstyle R`\ ; where \ :math:`\scriptstyle R`\  is one of the characters X, Y or Z.
   * - PDFpos
     - PDF peak position.
   * - PDFmag
     - PDF peak magnitude.
   * - PDFsig
     - PDF peak std. dev.
   * - Aspect ratio
     - Particle aspect ratio.
   * - Length
     - Cylinder length.
   * - Diameter
     - Cylinder/disk diameter.
   * - Thickness
     - Disk thickness.
   * - Shell thickness
     - Multiplier to get inner(<1) or outer(>1) sphere radius.
   * - Dist
     - Interparticle distance.
   * - VolFr
     - Dense scatterer volume fraction.
   * - epis
     - Sticky sphere epsilon.
   * - Sticky
     - Stickyness.
   * - Depth
     - Well depth.
   * - Width
     - Well width.
   * - Volume
     - Particle volume.
   * - Radius
     - Sphere/cylinder/disk radius.
   * - Mean
     - Particle mean radius.
   * - StdDev
     - Standard deviation in Mean.
   * - G
     - Guinier prefactor.
   * - Rg
     - Guinier radius of gyration.
   * - B
     - Porod prefactor.
   * - P
     - Porod power.
   * - Cutoff
     - Porod cutoff.
   * - PkInt
     - Bragg peak intensity.
   * - PkPos
     - Bragg peak position.
   * - PkSig
     - Bragg peak sigma.
   * - PkGam
     - Bragg peak gamma.
   * - e\ :math:`\scriptstyle W_0`\ \ :math:`\scriptstyle W_1`\  (example: ``e11``)
     - strain tensor e\ :math:`\scriptstyle W_0`\ \ :math:`\scriptstyle W_1`\ ; where \ :math:`\scriptstyle W_0`\  is one of the characters 1 or 2 and \ :math:`\scriptstyle W_1`\  is one of the characters 1 or 2.
   * - Dcalc
     - Calc. d-spacing.
   * - Back
     - background parameter.
   * - pos
     - peak position.
   * - int
     - peak intensity.
   * - WgtFrac
     - phase weight fraction.
   * - alpha
     - TOF profile term.
   * - alpha-\ :math:`\scriptstyle X`\  (example: ``alpha-0``)
     - Pink profile term; where \ :math:`\scriptstyle X`\  is one of the characters 0 or 1.
   * - beta-\ :math:`\scriptstyle Y`\  (example: ``beta-0``)
     - TOF/Pink profile term; where \ :math:`\scriptstyle Y`\  is one of the characters 0, 1 or q.
   * - sig-\ :math:`\scriptstyle Z`\  (example: ``sig-0``)
     - TOF profile term; where \ :math:`\scriptstyle Z`\  is one of the characters 0, 1, 2 or q.
   * - dif\ :math:`\scriptstyle [`\  (example: ``difA``)
     - TOF to d-space calibration; where \ :math:`\scriptstyle [`\  is one of the characters A, B or C.
   * - C\ :math:`\scriptstyle G_0`\ ,\ :math:`\scriptstyle G_1`\  (example: ``C10,10``)
     - spherical harmonics preferred orientation coef.; where \ :math:`\scriptstyle G_0`\  is one or more digits (0, 1,... 9) and \ :math:`\scriptstyle G_1`\  is one or more digits (0, 1,... 9).
   * - Pressure
     - Pressure level for measurement in MPa.
   * - Temperature
     - T value for measurement, K.
   * - FreePrm\ :math:`\scriptstyle N`\  (example: ``FreePrm1``)
     - User defined measurement parameter \ :math:`\scriptstyle N`\ ; where \ :math:`\scriptstyle N`\  is one of the characters 1, 2 or 3.
   * - Gonio. radius
     - Distance from sample to detector, mm.
