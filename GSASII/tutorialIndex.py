'''A catalog of GSAS-II tutorials with headings. This is the master list 
of GSAS-II tutorials and must be updated when tutorials are added. 
Each item has either one or three items. Titles are single item in a 
list or tuple. Tutorials have four items: 

     (a) the name of the directory,
     (b) the name of the web page,
     (c) a title for the tutorial and
     (d) a short text description (optional). 

Tutorials that depend on a previous tutorial being completed should 
have the title for the tutorial indented by five spaces.
    
Note that :data:`GSASIIpath.tutorialCatalog` is generated from this table using the
`makeGitTutorial.py` script (see
https://gsas-ii.readthedocs.io/en/latest/GSASIIscripts.html#other-scripts) in the 
GSAS-II tutorials repo (https://github.com/AdvancedPhotonSource/GSAS-II-tutorials) and 
which creates the tutorials.html file in that repo.
'''

tutorialIndex = (
    # tutorial dir,      web page file name,      title for page,  description
['Getting started'], #######################################################
    ['StartingGSASII', 'Starting GSAS.htm', 'Starting GSAS-II',
     '''An introduction to GSAS-II with starting instructions and a brief description of the displays.'''],

['Rietveld refinement'], #######################################################
    ['CWNeutron', 'Neutron CW Powder Data.htm', 'CW Neutron Powder fit for Yttrium-Iron Garnet',
     '''This shows a simple Rietveld refinement with constraints from CW neutron powder diffraction data.'''],
     
    ['LabData', 'Laboratory X.htm', 'Fitting laboratory X-ray powder data for fluoroapatite',
     '''This shows a simple Rietveld refinement with CuKa lab Bragg-Brentano powder data.'''],
     
    ['CWCombined', 'Combined refinement.htm', 'Combined X-ray/CW-neutron refinement of PbSO4',
     '''This shows Rietveld refinement of a structure with room temperature lab CuKa data and low temperature CW neutron data; 
     use is made of the lattice parameter offsets to account for thermal expansion.'''],
     
    ['TOF-CW Joint Refinement', 'TOF combined XN Rietveld refinement in GSAS.htm', 'Combined X-ray/TOF-neutron Rietveld refinement',
     '''This shows Rietveld refinement with high resolution synchrotron powder data and neutron TOF data'''],
     
    ['Simulation', 'SimTutorial.htm',  'Simulating Powder Diffraction with GSAS-II',
     '''This show how to create a simulated powder pattern from a lab diffractometer.'''],
     
['Advanced Rietveld'], #######################################################
    ['BkgFit', 'FitBkgTut.htm',  'Fitting the Starting Background using Fixed Points',
     '''This shows how to get an initial estimate of background parameters from a suite of fixed points 
     before beginning Rietveld refinement.'''],
     
    ['AutoBkg', 'AutoBkg.html',  'Using the "Auto Background" feature',
     '''This shows how to use the "Auto Background" feature in GSAS-II to get an estimate of background parameters for a 
     series of histograms with quite significant background levels. This estimate can be used to define a set of fixed points 
     or to define a "Fixed background histogram."'''],
     
    ['LeBail', 'LeBailSucrose.htm', 'Le Bail Intensity Extraction in GSAS-II - Sucrose',
     '''Shows the process of setting up a Le Bail fit, where reflection 
     intensities are treated as arbitrary, and how to converge the Le Bail
     intensities before a combined Le Bail/Least-Squares fit that 
     optimizes lattice, peak shape and background parameters.'''],
     
    ['CIFtutorial', 'CIFtutorial.html', 'Create a CIF for Publication',
     '''Shows how to create a full CIF for a project that includes the structure(s),
     bond distances/angles/ the observed and computed data, etc as well as documentary 
     information about the sample(s), instrument(s) and the project in a way that allows
     for updating the CIF without having to reenter any of that information. The tutorial 
     also explains how creation of template file can allow for reuse of that information.'''],

    ['RietPlot', 'PublicationPlot.htm', 'Create a Publication-Ready Rietveld Plot',
     '''Shows how to create a customized version of a plot from a fit, 
     with enlarged letters, different colors or symbols which can be written 
     as a bitmap file, a pdf file or be exported to the Grace or Igor Pro 
     plotting programs.'''],

    ['ParameterLimits', 'ParameterLimitsUse.html', 'Use of Parameter Limits',
     '''Shows how to set lower and upper limits on selected parameters to keep refinements from refining unreasonably.'''],
     
    ['RigidBody', 'RigidBodyRef.html', 'Rietveld Fitting with Rigid Bodies',
     '''Shows how to set up and refine with rigid bodies to simplify and improve the crystal structure model.'''],
     
['Magnetic Structure Analysis'], #######################################################
    ['SimpleMagnetic', 'SimpleMagnetic.htm',"Simple Magnetic Structure Analysis",
     '''Analysis of a simple antiferromagnet and a simple ferromagnet from CW neutron powder data'''],
     
    ['Magnetic-I', 'Magnetic Structures-I.htm',"Magnetic Structure Analysis-I",
     '''Analysis of a simple antiferromagnet using Bilbao k-SUBGROUPSMAG from CW neutron powder data'''],
     
    ['Magnetic-II', 'Magnetic-II.htm',"Magnetic Structure Analysis-II",
     '''Analysis of a antiferromagnet with change of space group using Bilbao k-SUBGROUPSMAG from CW neutron powder data'''],
     
    ['Magnetic-III', 'Magnetic-III.htm',"Magnetic Structure Analysis-III",
     '''Analysis of a Type IV antiferromagnet with a cell axis doubling using Bilbao k-SUBGROUPSMAG from CW neutron powder data'''],
         
    ['Magnetic-IV', 'Magnetic-IV.htm',"Magnetic Structure Analysis-IV",
     '''Analysis of a Type IV antiferromagnet with a lattice centering change using Bilbao k-SUBGROUPSMAG from CW neutron powder data'''],
         
    ['Magnetic-V', 'Magnetic-V.htm',"Magnetic Structure Analysis-V",
     '''Analysis of a complex Type IV antiferromagnet with two propagation vectorse using Bilbao k-SUBGROUPSMAG from TOF neutron powder data'''],

    ['k_vec_tutorial', 'k_vec_tutorial.html', 'k-vector searching in GSAS-II #1',
    '''Search neutron diffraction data of Er2Ge2O7 for a all-zero magnetic k-vector'''],

    ['k_vec_tutorial_non_zero', 'k_vec_tutorial_non_zero.html', 'k-vector searching in GSAS-II #2',
    '''Search neutron diffraction data used in Magnetic Structure Analysis-III for a non-zero magnetic k-vector'''],

    #['ExampleDir', 'ExamplePage.html', 'Example Tutorial Title', '''Example descriptive text'''],
         
['Parametric sequential fitting'], #######################################################
    ['SeqRefine', 'SequentialTutorial.htm', 'Sequential refinement of multiple datasets',
     '''This shows the fitting of a structural model to multiple data sets collected as a function of temperature (7-300K). 
     This tutorial is the prerequisite for the next one.'''],
     
    ['SeqParametric', 'ParametricFitting.htm', '     Parametric Fitting and Pseudo Variables for Sequential Fits',
     '''This explores the results of the sequential refinement obtained in the previous tutorial; includes 
     plotting of variables and fitting the changes with simple equations.'''],
     
     ['TOF Sequential Single Peak Fit','TOF Sequential Single Peak Fit.htm','Sequential fitting of single peaks and strain analysis of result',
      '''This shows the fitting of single peaks in a sequence of TOF powder patterns from a sample under load; includes
      fitting of the result to get Hookes Law coefficients for elastic deformations.'''],

['Structure solution'], #######################################################
    ['FitPeaks', 'Fit Peaks.htm', 'Fitting individual peaks & autoindexing',
     '''This covers two examples of selecting individual powder diffraction peaks, fitting them and then 
     indexing to determine the crystal lattice and possible space group. This is the prerequisite for the next two tutorials.'''],
     
    ['CFjadarite', 'Charge Flipping in GSAS.htm', '     Charge Flipping structure solution for jadarite',
     '''Solving the structure of jadarite (HLiNaSiB3O8) by charge flipping from Pawley extracted intensities
     from a high resolution synchrotron powder pattern.'''],
     
    ['CFsucrose', 'Charge Flipping - sucrose.htm','     Charge Flipping structure solution for sucrose',
          '''Solving the structure of sucrose (C12H22O11) by charge flipping from Pawley extracted intensities
     from a high resolution synchrotron powder pattern.'''],
     
    ['CFXraySingleCrystal', 'CFSingleCrystal.htm', 'Charge Flipping structure solution with Xray single crystal data',
     '''Solving the structure of dipyridyl disulfate by charge flipping and then refine the structure by least-squares.'''],
       
    ['TOF Charge Flipping', 'Charge Flipping with TOF single crystal data in GSASII.htm', 
     'Charge flipping with neutron TOF single crystal data',
     '''Solving the crystal structure or rubrene (C42H28) from single crystal neutron data 
     via charge flipping and then refine the structure by least squares.'''],
     
    ['MCsimanneal', 'MCSA in GSAS.htm', 'Monte-Carlo simulated annealing structure determination',
     '''Solving the structures of 3-aminoquinoline and α-d-lactose monohydrate from powder diffraction data 
     via Monte Carlo/Simulated Annealing (MC/SA).'''],
     
['PDF 1: RMCProfile Reverse Monte-Carlo PDF & S(Q) Modeling'], #######################################################
    ['RMCProfile-I', 'RMCProfile-I.htm','RMC Modeling with RMCProfile-I',
     '''Big box modeling for real and reciprocal space diffraction data for SF6'''],
    ['RMCProfile-II', 'RMCProfile-II.htm','RMC Modeling with RMCProfile-II',
     '''Big box modeling for real and reciprocal space diffraction data for SrTiO3'''],
    ['RMCProfile-III', 'RMCProfile-III.htm','RMC Modeling with RMCProfile-III',
     '''Combined x-ray/neutron big box modeling for real and reciprocal space diffraction data for GaPO4'''],
    ['RMCProfile-IV', 'RMCProfile-IV.htm','RMC Modeling with RMCProfile-IV',
     '''x-ray big box modeling with potential energy restraints for real and reciprocal space diffraction data for GaPO4'''],
    
['PDF 2: PDFfit Pair Distribution Function Modeling'], #######################################################
    ['PDFfit-I','PDFfit-I.htm','Small Box PDF modeling with PDFfit-I',
     '''Small box modeling of G(r); introduction to PDFfit'''],
    ['PDFfit-II','PDFfit-II.htm','Small Box PDF modeling with PDFfit-II',
     '''Small box modeling of G(r); using ISODISTORT mode analysis'''],
    ['PDFfit-III','PDFfit-III.htm','Sequential PDF fitting with PDFfit-III',
     '''Small box modeling of G(r); sequential fitting of a temperature series of G(r)'''],
    ['PDFfit-IV','PDFfit-IV.htm','Nanoparticle PDF fitting with PDFfit-IV',
     '''Small box modeling of G(r); fitting G(r) from nanoparticles'''],

['PDF 3: fullrmc Stochastic PDF & S(Q) Modeling'], #######################################################
    ['fullrmc-Ni', 'fullrmc-Ni.html','RMC & Rigid Body Modeling with fullrmc-I',
     '''Big box modeling with real space diffraction data for Ni'''],
    ['fullrmc-SF6', 'fullrmc-SF6.html','RMC & Rigid Body Modeling with fullrmc-II',
     '''Multiple approaches to big box modeling for real and reciprocal space diffraction data for SF6'''],

['Stacking Fault Modeling'], #######################################################
    ['StackingFaults-I', 'Stacking Faults-I.htm', 'Stacking fault simulations for diamond',
     '''This shows how to simulate the diffraction patterns from faulted diamond.'''],
     
    ['StackingFaults-II', 'Stacking Faults II.htm', 'Stacking fault simulations for Keokuk kaolinite',
     '''This shows how to simulate some diffraction patterns from well ordered Keokuk kaolinite (Al2Si2O5(OH)4) clay.'''],
     
    ['StackingFaults-III', 'Stacking Faults-III.htm', 'Stacking fault simulations for Georgia kaolinite',
     '''This shows how to simulate some diffraction patterns from poorly ordered Georgia kaolinite (Al2Si2O5(OH)4) clay.'''],

['Powder diffractometer calibration'], #######################################################
    ['CWInstDemo', 'FindProfParamCW.htm',  'Create Instrument Parameter File: Determine Starting Profile from a Standard',
     '''This shows how to determine profile parameters by fitting individual peaks
        with data collected on a standard using a lab diffractometer and then same them for reuse.'''],
    ['FPAfit', 'FPAfit.htm',  'Determining Profile Parameters with Fundamental Parameters',
     '''This shows how to determine profile parameters by fitting 
     peaks that are computed using the NIST Fundamental Parameters Python
     code. 
     Input is formulated to use FPA values similar to those in Topas.'''],      
    ['TOF Calibration', 'Calibration of a TOF powder diffractometer.htm', 'Calibration of a Neutron TOF diffractometer',
     '''This uses the fitted positions of all visible peaks in a pattern of NIST SRM 660b La11B6 
     (a=4.15689Å) obtained in a multiple single peak fit. The positions are compared to those expected from the 
     known lattice parameters to establish the diffractometer constants (difC, difA, difB and Zero) used for 
     calculating TOF peak positions from d-spacings. In addition, the peak fitting includes the various profile 
     coefficients thus fully describing the instrument contribution to the peak profiles.''' ],

['2D Image Processing'], #######################################################
    ['2DCalibration', 'Calibration of an area detector in GSAS.htm', 'Calibration of an area detector',
     '''A demonstration of calibrating a Perkin-Elmer area detector,  where the detector was intentionally tilted at 45 degrees.
     This exercise is the prerequisite for the next one.'''],
     
    ['2DIntegration', 'Integration of area detector data in GSAS.htm', '     Integration of area detector data',
     '''Integration of the image from a Perkin-Elmer area detector, where the detector was intentionally tilted at 45 degrees.'''],
     
    ['2DStrain', 'Strain fitting of 2D data in GSAS-II.htm', 'Strain fitting of 2D data',
     '''This show how to determine 3 strain tensor values using the method of He & Smith (Adv. in X-ray Anal. 41, 501, 1997) 
     directly froom a sequence of 2D imges from a loaded sample.'''],
    
    ['2DTexture', 'Texture analysis of 2D data in GSAS-II.htm', 'Texture analysis of 2D data',
     '''This shows 3 different methods for determining texture via spherical harmonics from 2D x-ray diffraction images. '''],
     
    ['DeterminingWavelength', 'DeterminingWavelength.html', 'Area Detector Calibration with Multiple Distances: Determine Wavelength',
     '''To get an accurate wavelength, without knowing the sample-to-detector distance accurately, images recorded with
     several different distances can be used. This exercise shows how to determine the wavelength from such a series. 
     This exercise is the prerequisite for the next one.'''],
     
    ['CalibrationTutorial', 'CalibrationTutorial.html', '    Area Detector Calibration with Multiple Distances: Calibrate Detector Distances',
     '''To get an accurate wavelength, without knowing the sample-to-detector distance accurately, images recorded with
     several different distances can be used. After using the previous exercise to determine the wavelength,
     this exercise calibrates the detector distances and shows examples of how to mask, integrate, and save those parameters
     for future reuse.'''],
                   
['Small-Angle Scattering'], #######################################################      
    ['SAsize', 'Small Angle Size Distribution.htm', 'Small angle x-ray data size distribution (alumina powder)',
     '''This shows how to determine the size distribution of particles using data from a constant 
     wavelength synchrotron X-ray USAXS instrument. This is the prerequisite for the next tutorial'''],
     
    ['SAfit', 'Fitting Small Angle Scattering Data.htm', '     Fitting small angle x-ray data (alumina powder)',
     '''This shows how to fit small angle scattering data using data from a constant wavelength synchrotron X-ray USAXS instrument. '''],
     
    ['SAimages', 'Small Angle Image Processing.htm', 'Image Processing of small angle x-ray data',
     '''This shows how to  reduce 2D SAXS data to create 1D absolute scaled data. '''],
     
    ['SAseqref', 'Sequential Refinement of Small Angle Scattering Data.htm', 'Sequential refinement with small angle scattering data',
     '''This shows how to fit USAXS small angle scattering data for a suite of samples to demonstrate the 
     sequential refinement technique in GSAS-II for SASD and demonstrates fitting with a hard sphere structure 
     factor for non-dilute systems. '''],

['Other'],    #######################################################
    ['MerohedralTwins', 'Merohedral twin refinement in GSAS.htm', 'Merohedral twin refinements',
     '''This shows how to use GSAS-II to refine the structure of a few single crystal structures where there is merohedral twinning. '''],
     
    ['TOF Single Crystal Refinement', 'TOF single crystal refinement in GSAS.htm', 'Single crystal refinement from TOF data',
     '''This shows how to refine the structure of sapphire (really corundum, Al2O3) from single crystal diffraction data 
     collected at the SNS on the TOPAZ instrument at room temperature.  '''],
     
    ['PythonScript','Scripting.htm','Scripting a GSAS-II Refinement from Python',
     '''This demonstrates the use of the GSASIIscriptable module. This uses a Python script to perform a refinement or 
     computation, but without use of the GSAS-II graphical user interface. This is a prerequisite for the next tutorial.'''],
     
    ['PythonScript','CommandLine.htm','     Running a GSAS-II Refinement from the Command Line',
     '''This shows a unix script that duplicates the previous Python Scripting GSAS-II tutorial. '''],
    
    ['ClusterAnalysis','Cluster and Outlier Analysis.htm', 'Cluster and Outlier Analysis',
     '''This gives an example of using Cluster and Outlier Analysis with PWDR data.'''],
    
    ['FontSize', 'FontSize.html', 'Changing the GSAS-II Font Size',
    '''Configuration variables allow GSAS-II to be customized. In this example the configuration variable that changes the font size 
    used in the GUI is changed.'''],

    #['ExampleDir', 'ExamplePage.html', 'Example Tutorial Title', '''Example descriptive text'''],
    )
    
