"""
Complete list of GSAS-II documentation sources.

All Home page, Help & Tutorial URLs resolve under:
    https://advancedphotonsource.github.io/GSAS-II-tutorials/
    URL lists in HOME_SOURCES, HELP_SOURCES and TUTORIAL_SOURCES

Additional HTML sources:
  - GSAS-II Programmer's Guide — 24 pages from readthedocs (https://gsas-ii.readthedocs.io/en/latest)
    URL list in READTHEDOCS_SOURCES
    
  - Powder Diffraction Crystallography book — 185 HTML pages (https://briantoby.github.io/PowderCrystallography)
    URL list in BOOK_HTML_SOURCES

Sources are selected based on the arguments supplied when ingest.py is run (e.g. --book includes BOOK_HTML_SOURCES)

"""

BASE_URL = "https://advancedphotonsource.github.io/GSAS-II-tutorials"

# ── Home / installation pages (22) ────────────────────────────────────────────

HOME_SOURCES = [
    "GSAS-II Home pages",
    {"title": "GSAS-II Home",                    "url": f"{BASE_URL}/index.html",             "category": "Home"},
    {"title": "About GSAS-II",                   "url": f"{BASE_URL}/AboutGSASII.html",        "category": "Home"},
    {"title": "Documentation Overview",          "url": f"{BASE_URL}/documentation.html",      "category": "Home"},
    {"title": "Tutorials Index",                 "url": f"{BASE_URL}/tutorials.html",           "category": "Home"},
    {"title": "Help Overview",                   "url": f"{BASE_URL}/help.html",               "category": "Home"},
    {"title": "Miscellaneous Notes",             "url": f"{BASE_URL}/misc.html",               "category": "Home"},
    {"title": "Installation Overview",           "url": f"{BASE_URL}/install.html",            "category": "Installation"},
    {"title": "Install with pip",                "url": f"{BASE_URL}/install-pip.html",        "category": "Installation"},
    {"title": "Install with pixi",               "url": f"{BASE_URL}/install_pixi.html",       "category": "Installation"},
    {"title": "Install GSAS2Full (Mac)",         "url": f"{BASE_URL}/install-g2f-mac.html",    "category": "Installation"},
    {"title": "Install GSAS2Full (Linux)",       "url": f"{BASE_URL}/install-g2f-linux.html",  "category": "Installation"},
    {"title": "Install GSAS2Full (Windows)",     "url": f"{BASE_URL}/install-g2f-win.html",    "category": "Installation"},
    {"title": "Install for Developers",          "url": f"{BASE_URL}/install_dev.html",        "category": "Installation"},
    {"title": "Install External Programs",       "url": f"{BASE_URL}/install-external.html",   "category": "Installation"},
    {"title": "Install Python Manually",         "url": f"{BASE_URL}/install-python.html",     "category": "Installation"},
    {"title": "GSAS-II Options",                 "url": f"{BASE_URL}/options.html",            "category": "Home"},
    {"title": "Proxy Configuration",             "url": f"{BASE_URL}/proxy.html",              "category": "Installation"},
    {"title": "Developer Notes",                 "url": f"{BASE_URL}/developers.html",         "category": "Development"},
    {"title": "Compiling Extensions",            "url": f"{BASE_URL}/compile.html",            "category": "Development"},
    {"title": "Mailing List",                    "url": f"{BASE_URL}/mailinglist.html",        "category": "Home"},
    {"title": "Bug Reporting",                   "url": f"{BASE_URL}/bug.html",               "category": "Home"},
    {"title": "General Index",                   "url": f"{BASE_URL}/genindex.html",           "category": "Home"},
]
    
# ── Help pages (42, skipping 404.html) ────────────────────────────────────────

HELP_SOURCES = [
    "GSAS-II Help pages",
    {"title": "Help: Application Window",            "url": f"{BASE_URL}/help/applicationwindow.html",   "category": "Help"},
    {"title": "Help: Main Menu",                     "url": f"{BASE_URL}/help/mainmenu.html",            "category": "Help"},
    {"title": "Help: Data Tree",                     "url": f"{BASE_URL}/help/datatree.html",            "category": "Help"},
    {"title": "Help: Common Tree Items",             "url": f"{BASE_URL}/help/commontreeitems.html",     "category": "Help"},
    {"title": "Help: Preface",                       "url": f"{BASE_URL}/help/preface.html",             "category": "Help"},
    {"title": "Help: Others",                        "url": f"{BASE_URL}/help/others.html",              "category": "Help"},
    {"title": "Help: Powder Diffraction Overview",   "url": f"{BASE_URL}/help/powder.html",              "category": "Help: Powder"},
    {"title": "Help: Powder Parent",                 "url": f"{BASE_URL}/help/powderparent.html",        "category": "Help: Powder"},
    {"title": "Help: Powder Comments",               "url": f"{BASE_URL}/help/powdercomments.html",      "category": "Help: Powder"},
    {"title": "Help: Powder Instrument Parameters",  "url": f"{BASE_URL}/help/powderinst.html",          "category": "Help: Powder"},
    {"title": "Help: Powder Sample Parameters",      "url": f"{BASE_URL}/help/powdersample.html",        "category": "Help: Powder"},
    {"title": "Help: Powder Limits",                 "url": f"{BASE_URL}/help/powderlimits.html",        "category": "Help: Powder"},
    {"title": "Help: Powder Background",             "url": f"{BASE_URL}/help/powderbkg.html",           "category": "Help: Powder"},
    {"title": "Help: Powder Peaks",                  "url": f"{BASE_URL}/help/powderpeaks.html",         "category": "Help: Powder"},
    {"title": "Help: Powder Peak Indexing",          "url": f"{BASE_URL}/help/powderindexppeaks.html",   "category": "Help: Powder"},
    {"title": "Help: Powder Cells",                  "url": f"{BASE_URL}/help/powdercells.html",         "category": "Help: Powder"},
    {"title": "Help: Powder Reflections",            "url": f"{BASE_URL}/help/powderrefs.html",          "category": "Help: Powder"},
    {"title": "Help: Peak List",                     "url": f"{BASE_URL}/help/peaks.html",               "category": "Help: Powder"},
    {"title": "Help: Sequential Refinement",         "url": f"{BASE_URL}/help/sequential.html",          "category": "Help: Powder"},
    {"title": "Help: Phase Overview",                "url": f"{BASE_URL}/help/phaseoverview.html",       "category": "Help: Phase"},
    {"title": "Help: Phase General",                 "url": f"{BASE_URL}/help/phasegeneral.html",        "category": "Help: Phase"},
    {"title": "Help: Phase Atoms",                   "url": f"{BASE_URL}/help/phaseatoms.html",          "category": "Help: Phase"},
    {"title": "Help: Phase Data",                    "url": f"{BASE_URL}/help/phasedata.html",           "category": "Help: Phase"},
    {"title": "Help: Phase Draw Options",            "url": f"{BASE_URL}/help/phasedrawopts.html",       "category": "Help: Phase"},
    {"title": "Help: Phase Draw Atoms",              "url": f"{BASE_URL}/help/phasedrawatoms.html",      "category": "Help: Phase"},
    {"title": "Help: Phase Map Peaks",               "url": f"{BASE_URL}/help/phasemappeaks.html",       "category": "Help: Phase"},
    {"title": "Help: Phase Pawley",                  "url": f"{BASE_URL}/help/phasepawley.html",         "category": "Help: Phase"},
    {"title": "Help: Phase Texture",                 "url": f"{BASE_URL}/help/phasetexture.html",        "category": "Help: Phase"},
    {"title": "Help: Phase Layers (Stacking)",       "url": f"{BASE_URL}/help/phaselayers.html",         "category": "Help: Phase"},
    {"title": "Help: Phase Waves (Modulated)",       "url": f"{BASE_URL}/help/phasewave.html",           "category": "Help: Phase"},
    {"title": "Help: Phase Rigid Bodies",            "url": f"{BASE_URL}/help/phaseRB.html",             "category": "Help: Phase"},
    {"title": "Help: Phase RMC",                     "url": f"{BASE_URL}/help/phaseRMC.html",            "category": "Help: Phase"},
    {"title": "Help: Phase MCSA",                    "url": f"{BASE_URL}/help/phasemcsa.html",           "category": "Help: Phase"},
    {"title": "Help: Phase ISODISTORT",              "url": f"{BASE_URL}/help/phaseisodistort.html",     "category": "Help: Phase"},
    {"title": "Help: Phase DYSNOMIA",                "url": f"{BASE_URL}/help/phasedysnomia.html",       "category": "Help: Phase"},
    {"title": "Help: Single Crystal",                "url": f"{BASE_URL}/help/singlecrystal.html",       "category": "Help: Single Crystal"},
    {"title": "Help: Image Processing",              "url": f"{BASE_URL}/help/image.html",               "category": "Help: Image"},
    {"title": "Help: Small Angle Scattering",        "url": f"{BASE_URL}/help/smallanglescattering.html","category": "Help: SAXS"},
    {"title": "Help: Pair Distribution Function",    "url": f"{BASE_URL}/help/pairdistribution.html",    "category": "Help: PDF"},
    {"title": "Help: Reflectometry",                 "url": f"{BASE_URL}/help/reflectometry.html",       "category": "Help: Reflectometry"},
    {"title": "Help: Cluster Analysis",              "url": f"{BASE_URL}/help/cluster.html",             "category": "Help: Analysis"},
    {"title": "Help: Index",                         "url": f"{BASE_URL}/help/index.html",               "category": "Help"},
]

# ── Tutorials (62, skipping tutorial_template) ────────────────────────────────

TUTORIAL_SOURCES = [
    "GSAS-II Tutorials",
    # Getting Started
    {"title": "Starting GSAS-II",
     "url": f"{BASE_URL}/StartingGSASII/Starting%20GSAS.htm",
     "category": "Getting Started"},

    # Rietveld Refinement
    {"title": "Fitting CW Neutron Powder Data (YIG)",
     "url": f"{BASE_URL}/CWNeutron/Neutron%20CW%20Powder%20Data.htm",
     "category": "Rietveld Refinement"},
    {"title": "Fitting Laboratory X-ray Powder Data (Fluoroapatite)",
     "url": f"{BASE_URL}/LabData/Laboratory%20X.htm",
     "category": "Rietveld Refinement"},
    {"title": "Combined X-ray and CW-Neutron Refinement (PbSO4)",
     "url": f"{BASE_URL}/CWCombined/Combined%20refinement.htm",
     "category": "Rietveld Refinement"},
    {"title": "Combined X-ray and TOF-Neutron Rietveld Refinement",
     "url": f"{BASE_URL}/TOF-CW%20Joint%20Refinement/TOF%20combined%20XN%20Rietveld%20refinement%20in%20GSAS.htm",
     "category": "Rietveld Refinement"},
    {"title": "Simulating Powder Diffraction with GSAS-II",
     "url": f"{BASE_URL}/Simulation/SimTutorial.htm",
     "category": "Rietveld Refinement"},

    # Background and Profile
    {"title": "Fitting the Background using Fixed Points",
     "url": f"{BASE_URL}/BkgFit/FitBkgTut.htm",
     "category": "Background & Profile"},
    {"title": "Using the Auto Background Feature",
     "url": f"{BASE_URL}/AutoBkg/AutoBkg.html",
     "category": "Background & Profile"},
    {"title": "Le Bail Intensity Extraction (Sucrose)",
     "url": f"{BASE_URL}/LeBail/LeBailSucrose.htm",
     "category": "Background & Profile"},
    {"title": "Determining Profile Parameters with Fundamental Parameters",
     "url": f"{BASE_URL}/FPAfit/FPAfit.htm",
     "category": "Background & Profile"},
    {"title": "Create Instrument Parameter File (CW Profile from Standard)",
     "url": f"{BASE_URL}/CWInstDemo/FindProfParamCW.html",
     "category": "Background & Profile"},
    {"title": "Use of Parameter Limits",
     "url": f"{BASE_URL}/ParameterLimits/ParameterLimitsUse.html",
     "category": "Background & Profile"},
    {"title": "Rietveld Fitting with Rigid Bodies",
     "url": f"{BASE_URL}/RigidBody/RigidBodyRef.html",
     "category": "Background & Profile"},

    # Sequential Refinement
    {"title": "Sequential Refinement of Multiple Datasets",
     "url": f"{BASE_URL}/SeqRefine/SequentialTutorial.htm",
     "category": "Sequential Refinement"},
    {"title": "Parametric Fitting and Pseudo Variables for Sequential Fits",
     "url": f"{BASE_URL}/SeqParametric/ParametricFitting.htm",
     "category": "Sequential Refinement"},
    {"title": "Sequential Fitting of Single Peaks and Strain Analysis",
     "url": f"{BASE_URL}/TOF%20Sequential%20Single%20Peak%20Fit/TOF%20Sequential%20Single%20Peak%20Fit.htm",
     "category": "Sequential Refinement"},
    {"title": "Sequential Refinement with Small Angle Scattering Data",
     "url": f"{BASE_URL}/SAseqref/Sequential%20Refinement%20of%20Small%20Angle%20Scattering%20Data.htm",
     "category": "Sequential Refinement"},

    # Magnetic Structures
    {"title": "Simple Magnetic Structure Analysis",
     "url": f"{BASE_URL}/SimpleMagnetic/SimpleMagnetic.htm",
     "category": "Magnetic Structures"},
    {"title": "Register for Bilbao Crystallographic Server",
     "url": f"{BASE_URL}/RegisterBilbao/RegisterBilbao.html",
     "category": "Magnetic Structures"},
    {"title": "Magnetic Structure Analysis I",
     "url": f"{BASE_URL}/Magnetic-I/Magnetic%20Structures-I.htm",
     "category": "Magnetic Structures"},
    {"title": "Magnetic Structure Analysis II",
     "url": f"{BASE_URL}/Magnetic-II/Magnetic-II.htm",
     "category": "Magnetic Structures"},
    {"title": "Magnetic Structure Analysis III",
     "url": f"{BASE_URL}/Magnetic-III/Magnetic-III.htm",
     "category": "Magnetic Structures"},
    {"title": "Magnetic Structure Analysis IV",
     "url": f"{BASE_URL}/Magnetic-IV/Magnetic-IV.htm",
     "category": "Magnetic Structures"},
    {"title": "Magnetic Structure Analysis V",
     "url": f"{BASE_URL}/Magnetic-V/Magnetic-V.htm",
     "category": "Magnetic Structures"},
    {"title": "k-vector Searching in GSAS-II (zero vector)",
     "url": f"{BASE_URL}/k_vec_tutorial/k_vec_tutorial.html",
     "category": "Magnetic Structures"},
    {"title": "k-vector Searching in GSAS-II (non-zero vector)",
     "url": f"{BASE_URL}/k_vec_tutorial_non_zero/k_vec_tutorial_non_zero.html",
     "category": "Magnetic Structures"},
    {"title": "Use of ISODISTORT with k-vector from GSAS-II",
     "url": f"{BASE_URL}/k_vec_isodistort/k_vec_isodistort.html",
     "category": "Magnetic Structures"},

    # Structure Solution
    {"title": "Fitting Individual Peaks and Autoindexing",
     "url": f"{BASE_URL}/FitPeaks/Fit%20Peaks.htm",
     "category": "Structure Solution"},
    {"title": "Charge Flipping Structure Solution (Jadarite)",
     "url": f"{BASE_URL}/CFjadarite/Charge%20Flipping%20in%20GSAS.htm",
     "category": "Structure Solution"},
    {"title": "Charge Flipping Structure Solution (Sucrose)",
     "url": f"{BASE_URL}/CFsucrose/Charge%20Flipping%20-%20sucrose.htm",
     "category": "Structure Solution"},
    {"title": "Charge Flipping from X-ray Single Crystal Data",
     "url": f"{BASE_URL}/CFXraySingleCrystal/CFSingleCrystal.htm",
     "category": "Structure Solution"},
    {"title": "Charge Flipping from Neutron TOF Single Crystal Data",
     "url": f"{BASE_URL}/TOF%20Charge%20Flipping/Charge%20Flipping%20with%20TOF%20single%20crystal%20data%20in%20GSASII.htm",
     "category": "Structure Solution"},
    {"title": "Monte-Carlo Simulated Annealing Structure Determination",
     "url": f"{BASE_URL}/MCsimanneal/MCSA%20in%20GSAS.htm",
     "category": "Structure Solution"},
    {"title": "Merohedral Twin Refinements",
     "url": f"{BASE_URL}/MerohedralTwins/Merohedral%20twin%20refinement%20in%20GSAS.htm",
     "category": "Structure Solution"},
    {"title": "Single Crystal Refinement from TOF Data",
     "url": f"{BASE_URL}/TOF%20Single%20Crystal%20Refinement/TOF%20single%20crystal%20refinement%20in%20GSAS.htm",
     "category": "Structure Solution"},

    # PDF: RMCProfile
    {"title": "RMC Modeling with RMCProfile I",
     "url": f"{BASE_URL}/RMCProfile-I/RMCProfile-I.htm",
     "category": "PDF: RMCProfile"},
    {"title": "RMC Modeling with RMCProfile II",
     "url": f"{BASE_URL}/RMCProfile-II/RMCProfile-II.htm",
     "category": "PDF: RMCProfile"},
    {"title": "RMC Modeling with RMCProfile III",
     "url": f"{BASE_URL}/RMCProfile-III/RMCProfile-III.htm",
     "category": "PDF: RMCProfile"},
    {"title": "RMC Modeling with RMCProfile IV",
     "url": f"{BASE_URL}/RMCProfile-IV/RMCProfile-IV.htm",
     "category": "PDF: RMCProfile"},

    # PDF: PDFfit
    {"title": "Small Box PDF Modeling with PDFfit I",
     "url": f"{BASE_URL}/PDFfit-I/PDFfit-I.htm",
     "category": "PDF: PDFfit"},
    {"title": "Small Box PDF Modeling with PDFfit II",
     "url": f"{BASE_URL}/PDFfit-II/PDFfit-II.htm",
     "category": "PDF: PDFfit"},
    {"title": "Sequential PDF Fitting with PDFfit III",
     "url": f"{BASE_URL}/PDFfit-III/PDFfit-III.htm",
     "category": "PDF: PDFfit"},
    {"title": "Nanoparticle PDF Fitting with PDFfit IV",
     "url": f"{BASE_URL}/PDFfit-IV/PDFfit-IV.htm",
     "category": "PDF: PDFfit"},

    # PDF: fullrmc
    {"title": "RMC and Rigid Body Modeling with fullrmc (Ni)",
     "url": f"{BASE_URL}/fullrmc-Ni/fullrmc-Ni.html",
     "category": "PDF: fullrmc"},
    {"title": "RMC and Rigid Body Modeling with fullrmc (SF6)",
     "url": f"{BASE_URL}/fullrmc-SF6/fullrmc-SF6.html",
     "category": "PDF: fullrmc"},

    # Stacking Faults
    {"title": "Stacking Fault Simulations (Diamond)",
     "url": f"{BASE_URL}/StackingFaults-I/Stacking%20Faults-I.htm",
     "category": "Stacking Faults"},
    {"title": "Stacking Fault Simulations (Keokuk Kaolinite)",
     "url": f"{BASE_URL}/StackingFaults-II/Stacking%20Faults%20II.htm",
     "category": "Stacking Faults"},
    {"title": "Stacking Fault Simulations (Georgia Kaolinite)",
     "url": f"{BASE_URL}/StackingFaults-III/Stacking%20Faults-III.htm",
     "category": "Stacking Faults"},

    # TOF Calibration
    {"title": "Calibration of a Neutron TOF Diffractometer",
     "url": f"{BASE_URL}/TOF%20Calibration/Calibration%20of%20a%20TOF%20powder%20diffractometer.htm",
     "category": "Calibration"},

    # 2D Image Processing
    {"title": "Calibration of an Area Detector",
     "url": f"{BASE_URL}/2DCalibration/Calibration%20of%20an%20area%20detector%20in%20GSAS.htm",
     "category": "2D Image Processing"},
    {"title": "Integration of Area Detector Data",
     "url": f"{BASE_URL}/2DIntegration/Integration%20of%20area%20detector%20data%20in%20GSAS.htm",
     "category": "2D Image Processing"},
    {"title": "Strain Fitting of 2D Data",
     "url": f"{BASE_URL}/2DStrain/Strain%20fitting%20of%202D%20data%20in%20GSAS-II.htm",
     "category": "2D Image Processing"},
    {"title": "Texture Analysis of 2D Data",
     "url": f"{BASE_URL}/2DTexture/Texture%20analysis%20of%202D%20data%20in%20GSAS-II.htm",
     "category": "2D Image Processing"},
    {"title": "Area Detector Calibration: Determine Wavelength",
     "url": f"{BASE_URL}/DeterminingWavelength/DeterminingWavelength.html",
     "category": "2D Image Processing"},
    {"title": "Area Detector Calibration: Detector Distances",
     "url": f"{BASE_URL}/CalibrationTutorial/CalibrationTutorial.html",
     "category": "2D Image Processing"},

    # Small Angle Scattering
    {"title": "Small Angle X-ray Data Size Distribution",
     "url": f"{BASE_URL}/SAsize/Small%20Angle%20Size%20Distribution.htm",
     "category": "Small Angle Scattering"},
    {"title": "Fitting Small Angle X-ray Data",
     "url": f"{BASE_URL}/SAfit/Fitting%20Small%20Angle%20Scattering%20Data.htm",
     "category": "Small Angle Scattering"},
    {"title": "Image Processing of Small Angle X-ray Data",
     "url": f"{BASE_URL}/SAimages/Small%20Angle%20Image%20Processing.htm",
     "category": "Small Angle Scattering"},

    # Python Scripting
    {"title": "Scripting a GSAS-II Refinement from Python",
     "url": f"{BASE_URL}/PythonScript/Scripting.htm",
     "category": "Python Scripting"},
    {"title": "Running a GSAS-II Refinement from the Command Line",
     "url": f"{BASE_URL}/PythonScript/CommandLine.htm",
     "category": "Python Scripting"},

    # Publication
    {"title": "Create a CIF for Publication",
     "url": f"{BASE_URL}/CIFtutorial/CIFtutorial.html",
     "category": "Publication"},
    {"title": "Create a Publication-Ready Rietveld Plot",
     "url": f"{BASE_URL}/RietPlot/PublicationPlot.htm",
     "category": "Publication"},

    # Misc
    {"title": "Cluster and Outlier Analysis",
     "url": f"{BASE_URL}/ClusterAnalysis/Cluster and Outlier Analysis.htm",
     "category": "Analysis"},
    {"title": "Changing the GSAS-II Font Size",
     "url": f"{BASE_URL}/FontSize/FontSize.html",
     "category": "Miscellaneous"},
]

# ── GSAS-II Programmer's Guide — readthedocs HTML (preferred over PDF) ────────

_RTD = "https://gsas-ii.readthedocs.io/en/latest"

READTHEDOCS_SOURCES = [
    "GSAS-II programmers' manual",
    {"title": "GSAS-II Packages Overview",          "url": f"{_RTD}/packages.html",          "category": "Programmer's Guide"},
    {"title": "GSAS-II Versioning",                 "url": f"{_RTD}/versioning.html",         "category": "Programmer's Guide"},
    {"title": "Object/Variable Organization",       "url": f"{_RTD}/objvarorg.html",          "category": "Programmer's Guide"},
    {"title": "GSASII module",                      "url": f"{_RTD}/GSASII.html",             "category": "Programmer's Guide"},
    {"title": "GSASIIobj module",                   "url": f"{_RTD}/GSASIIobj.html",          "category": "Programmer's Guide"},
    {"title": "GSASIIutil module",                  "url": f"{_RTD}/GSASIIutil.html",         "category": "Programmer's Guide"},
    {"title": "GSASIIGUIr module",                  "url": f"{_RTD}/GSASIIGUIr.html",         "category": "Programmer's Guide"},
    {"title": "GSASIIGUI module",                   "url": f"{_RTD}/GSASIIGUI.html",          "category": "Programmer's Guide"},
    {"title": "GSASIIdata module",                  "url": f"{_RTD}/GSASIIdata.html",         "category": "Programmer's Guide"},
    {"title": "GSASIIstruc module",                 "url": f"{_RTD}/GSASIIstruc.html",        "category": "Programmer's Guide"},
    {"title": "GSASIImapvars module",               "url": f"{_RTD}/GSASIImapvars.html",      "category": "Programmer's Guide"},
    {"title": "GSASIIimage module",                 "url": f"{_RTD}/GSASIIimage.html",        "category": "Programmer's Guide"},
    {"title": "GSASIImath module",                  "url": f"{_RTD}/GSASIImath.html",         "category": "Programmer's Guide"},
    {"title": "GSAS-II Index",                      "url": f"{_RTD}/GSASIIindex.html",        "category": "Programmer's Guide"},
    {"title": "Graphics modules",                   "url": f"{_RTD}/graphics.html",           "category": "Programmer's Guide"},
    {"title": "GSASIIpwd module",                   "url": f"{_RTD}/GSASIIpwd.html",          "category": "Programmer's Guide"},
    {"title": "Small Angle Scattering module",      "url": f"{_RTD}/SAS.html",                "category": "Programmer's Guide"},
    {"title": "GSASIIscriptable module",            "url": f"{_RTD}/GSASIIscriptable.html",   "category": "Programmer's Guide"},
    {"title": "GSASIIscripts module",               "url": f"{_RTD}/GSASIIscripts.html",      "category": "Programmer's Guide"},
    {"title": "GSASIIweb module",                   "url": f"{_RTD}/GSASIIweb.html",          "category": "Programmer's Guide"},
    {"title": "Import modules",                     "url": f"{_RTD}/imports.html",            "category": "Programmer's Guide"},
    {"title": "Export modules",                     "url": f"{_RTD}/exports.html",            "category": "Programmer's Guide"},
    {"title": "G2tools module",                     "url": f"{_RTD}/G2tools.html",            "category": "Programmer's Guide"},
    #{"title": "GSAS-II General Index",              "url": f"{_RTD}/indices.html",            "category": "Programmer's Guide"},
]

# ── PDF sources ────────────────────────────────────────────────────────────────

PDF_SOURCES = []  # Book and Programmer's Guide now available as HTML

# ── Powder Diffraction Crystallography book — 185 HTML pages (Brian Toby) ─────
# https://briantoby.github.io/PowderCrystallography/
# Pages are accessible by direct URL; include via `gsas2-query --setup --book`.

_BOOK_BASE = "https://briantoby.github.io/PowderCrystallography"

# Placeholder — populated lazily by get_book_sources() at ingest time, not import time.
BOOK_HTML_SOURCES = [
    "Powder Diff. Cryst. Book book",
    {"title": "Powder Diff. Cryst. Book (Contents)", "url": f"{_BOOK_BASE}/HTML-template.html", "category": "Powder Crystallography Book"},
]


def _url_exists(url):
    import requests
    try:
        response = requests.head(url, allow_redirects=True, timeout=5)
        return response.status_code == 200
    except requests.RequestException:
        return False


def get_book_sources():
    """Probe the Powder Crystallography book URLs and return the full source list.

    Called at ingest time (not import time) so that importing sources.py never
    makes network requests.  Returns a list with the string label as the first
    element, followed by one dict per discovered chapter/section page.
    """
    sources = [
        "Powder Diff. Cryst. Book book",
        {"title": "Powder Diff. Cryst. Book (Contents)", "url": f"{_BOOK_BASE}/HTML-template.html", "category": "Powder Crystallography Book"},
    ]

    _i = 0
    lastURL = ''
    while True:
        _i += 1
        url = f"{_BOOK_BASE}/HTML-templatech{_i}.html"
        if _url_exists(url):
            lastURL = url
            sources.append({
                "title": f"Powder Diff. Cryst. Book Chapter {_i}",
                "url": url,
                "category": "Powder Crystallography Book",
            })
        else:
            print(f'Last Powder Diff. Cryst. Book Chapter is {_i-1} ({lastURL})')
            break

    _i = 0
    lastURL = ''
    while True:
        _i += 1
        url = f"{_BOOK_BASE}/HTML-templatese{_i}.html"
        if _url_exists(url):
            lastURL = url
            sources.append({
                "title": f"Powder Diff. Cryst. Book Section {_i}",
                "url": url,
                "category": "Powder Crystallography Book",
            })
        else:
            print(f'Last Powder Diff. Cryst. Book Section is {_i-1}  ({lastURL})')
            break

    return sources
        
# ── Dynamic tutorial list from GSAS-II's tutorialIndex.py ────────────────────
# Fetched at ingest time so new tutorials are picked up automatically.
# Falls back to the hardcoded TUTORIAL_SOURCES on any network or parse error.

_TUTORIAL_INDEX_URL = (
    "https://raw.githubusercontent.com/AdvancedPhotonSource/GSAS-II"
    "/main/GSASII/tutorialIndex.py"
)


def get_tutorial_sources():
    """Return tutorial sources from GSAS-II's canonical tutorialIndex.py.

    Uses ast.literal_eval (not eval) so untrusted file content is never executed.
    Falls back to the hardcoded TUTORIAL_SOURCES list on any error.
    """
    import ast
    import urllib.request

    try:
        req = urllib.request.Request(
            _TUTORIAL_INDEX_URL, headers={"User-Agent": "gsas2-query"}
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            content = resp.read().decode()

        tree = ast.parse(content)
        index_value = None
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id == "tutorialIndex":
                        index_value = ast.literal_eval(node.value)
                        break
                if index_value is not None:
                    break

        if not index_value:
            return TUTORIAL_SOURCES

        sources = ['GSAS-II tutorials']
        for entry in index_value:
            if len(entry) == 4:
                directory, filename, title, _ = entry
                sources.append({
                    "title": title.strip(),
                    "url": f"{BASE_URL}/{directory}/{filename}",
                    "category": "Tutorial",
                })
        return sources if sources else TUTORIAL_SOURCES

    except Exception:
        return TUTORIAL_SOURCES
