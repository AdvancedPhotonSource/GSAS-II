.. _GSASIIscriptable:

=========================================
 *GSASIIscriptable: Scripting Interface*
=========================================


*Summary/Contents*
==================

Routines to use an increasing amount of GSAS-II's capabilities from scripts, 
without use of the graphical user interface (GUI). GSASIIscriptable can create and access
GSAS-II project (.gpx) files and can directly perform image handling and refinements.  
The module defines wrapper classes (inheriting from :class:`~GSASIIscriptable.G2ObjectWrapper`) for a growing number 
of data tree items.

GSASIIscriptable can be used in two ways. It offers a command-line mode, but 
the more widely used and more powerful mode of GSASIIscriptable is 
used is via Python scripts that 
call the module's application interface (API), these are summarized immediately below and are documented in the :ref:`complete API documentation <API>` section.

While the command-line mode 
provides access a number of features without writing Python scripts 
via shell/batch commands (see :ref:`CommandlineInterface`), use in practice
seems somewhat clumsy. Command-line mode
is no longer being developed and its use is discouraged.

.. contents:: Scripting Documentation Contents
   :depth: 2
   :backlinks: none
	       
Installation of GSASIIscriptable
================================

GSASIIscriptable is included as part of a standard GSAS-II installation that includes the GSAS-II GUI (as described in the `installation instructions <https://advancedphotonsource.github.io/GSAS-II-tutorials/install.html>`_). People who will will use scripting extensively will still need access to the GUI
for some activities, since the scripting API has not yet been extended to all
features of GSAS-II and even if that is ever completed, there will still be some things that GSAS-II does with the GUI would be almost impossible to implement without a interactive graphical view of the data.

Nonetheless, there may be times where it does make sense to install GSAS-II without all of the GUI components, for example on a compute server.
The minimal requirements for use of GSASIIscriptable are only Python, numpy and scipy, but additional optional packages that can be utilized are described in 
the :ref:`ScriptingRequirements` section of the requirements chapter, which also provides some installation instructions.

In a standard GSAS-II installation, no changes are made to Python. When the GUI is invoked, a small script or Windows batch file is used to start GSAS-II inside Python. When
GSASIIscriptable is used, Python must be provided with the location of the GSAS-II files. There are two ways this can be done:

 #. define the GSAS-II installation location in the Python ``sys.path``, or
 #. install a reference to GSAS-II inside Python. 

The latter method requires an extra installation step, but has the advantage that
it allows writing portable GSAS-II scripts. This is discussed further in the
:ref:`ScriptingShortcut` section of this chapter.

Application Interface (API) Summary
===================================
This section of the documentation provides an overview to API, with full documentation 
in the :ref:`API` section. The typical API use will be with a Python script, such as 
what is found in :ref:`CodeExamples`. Most functionality is provided via the objects and methods
summarized below.

Overview of Classes
-------------------

.. tabularcolumns:: |l|p{4in}|

========================================   ===============================================================================================================
Scripting class name                              Description 
========================================   ===============================================================================================================
:ref:`G2Project <Class_G2Project>`         :class:`~GSASIIscriptable.G2Project`:
                                           A GSAS-II project file; provides references to objects below,
                                           each corresponding to a tree item 
                                           (exception is :class:`~GSASIIscriptable.G2AtomRecord`)

:ref:`G2Phase <Class_G2Phase>`             :class:`~GSASIIscriptable.G2Phase`:
                                           Provides phase information access 
                                           (also provides access to atom info via :class:`~GSASIIscriptable.G2AtomRecord`)

:ref:`G2AtomRecord <Class_G2AtomRecord>`   :class:`~GSASIIscriptable.G2AtomRecord`:
                                           Access to an atom within a phase

:ref:`G2PwdrData <Class_G2PwdrData>`       :class:`~GSASIIscriptable.G2PwdrData`:
                                           Access to powder histogram info

:ref:`G2Single <Class_G2Single>`           :class:`~GSASIIscriptable.G2Single`: Access to single crystal histogram info

:ref:`G2Image <Class_G2Image>`             :class:`~GSASIIscriptable.G2Image`: Access to image info

:ref:`G2PDF <Class_G2PDF>`                 :class:`~GSASIIscriptable.G2PDF`: PDF histogram info

:ref:`G2SmallAngle <Class_G2SmallAngle>`          :class:`~GSASIIscriptable.G2SmallAngle`: Small Angle scattering histogram info

:ref:`G2SeqRefRes <Class_G2SeqRefRes>`     :class:`~GSASIIscriptable.G2SeqRefRes`:
                                           The sequential results table
========================================   ===============================================================================================================

Independent Functions
----------------------

A small number of Scriptable routines do not require existence of a G2Project object. 

.. tabularcolumns:: |l|p{4in}|

===================================================   ===============================================================================================================
method                                                Use
===================================================   ===============================================================================================================
:func:`~GSASIIscriptable.ShowVersions`                Shows Python and GSAS-II version information
:func:`~GSASIIscriptable.GenerateReflections`         Generates a list of unique powder reflections 
:func:`~GSASIIscriptable.SetPrintLevel`               Sets the amount of output generated when running a script 
:func:`~GSASIIscriptable.installScriptingShortcut`    Installs GSASIIscriptable within Python as G2script
===================================================   ===============================================================================================================

.. _Class_G2Project:

Class :class:`~GSASIIscriptable.G2Project`
---------------------------------------------

  All GSASIIscriptable scripts will need to create a :class:`~GSASIIscriptable.G2Project` object 
  either for a new GSAS-II project or to read in an existing project (.gpx) file. 
  The most commonly used routines in this object are:

.. tabularcolumns:: |l|p{3.in}|

====================================================================    ===============================================================================================================
method                                                                  Use
====================================================================    ===============================================================================================================
:meth:`~GSASIIscriptable.G2Project.save`                                Writes the current project to disk.

:meth:`~GSASIIscriptable.G2Project.add_powder_histogram`                Used to read in powder diffraction data into a project file.

:meth:`~GSASIIscriptable.G2Project.add_simulated_powder_histogram`      Defines a "dummy" powder diffraction data that will be simulated after a refinement step.

:meth:`~GSASIIscriptable.G2Project.add_image`                           Reads in an image into a project.

:meth:`~GSASIIscriptable.G2Project.add_phase`                           Adds a phase to a project

:meth:`~GSASIIscriptable.G2Project.add_PDF`                             Adds a PDF entry to a project (does not compute it)

:meth:`~GSASIIscriptable.G2Project.add_single_histogram`                Used to read in a single crystal diffraction dataset into a project file.

:meth:`~GSASIIscriptable.G2Project.histogram`                           Finds a histogram from an object, name or random id reference, returning a
                                                                        a :class:`~GSASIIscriptable.G2PwdrData` or :class:`~GSASIIscriptable.G2Single` object.

:meth:`~GSASIIscriptable.G2Project.histograms`                          Provides a list of histograms in the current project, as :class:`~GSASIIscriptable.G2PwdrData` or
                                                                        as :class:`~GSASIIscriptable.G2Single` objects.

:meth:`~GSASIIscriptable.G2Project.histType`                            Determines the histogram type from an object, name or random id reference.
                                                                        
:meth:`~GSASIIscriptable.G2Project.phases`                              Provides a list of phases defined in the current project, as :class:`~GSASIIscriptable.G2Phase` objects

:meth:`~GSASIIscriptable.G2Project.images`                              Provides a list of images in the current project, as :class:`~GSASIIscriptable.G2Image` objects

:meth:`~GSASIIscriptable.G2Project.pdfs`                                Provides a list of PDFs in the current project, as :class:`~GSASIIscriptable.G2PDF` objects

:meth:`~GSASIIscriptable.G2Project.seqref`                              Returns a :class:`~GSASIIscriptable.G2SeqRefRes` object if there are Sequential Refinement results

:meth:`~GSASIIscriptable.G2Project.do_refinements`                      This is passed a list of dictionaries, where each dict defines a refinement step. 
                                                                        Passing a list with a single empty dict initiates a refinement with the current
                                                                        parameters and flags. A refinement dict sets up a single refinement step 
                                                                        (as described in :ref:`Project_dicts`). Also see :ref:`Refinement_recipe`.

:meth:`~GSASIIscriptable.G2Project.set_refinement`                      This is passed a single dict which is used to set parameters and flags.
                                                                        These actions can be performed also in :meth:`~GSASIIscriptable.G2Project.do_refinements`. 
:meth:`~GSASIIscriptable.G2Project.get_Variable`                        Retrieves the value and esd for a parameter
:meth:`~GSASIIscriptable.G2Project.get_Covariance`                      Retrieves values and covariance for a set of refined parameters
:meth:`~GSASIIscriptable.G2Project.set_Controls`                        Set overall GSAS-II control settings such as number of cycles and to set up a sequential
                                                                        fit. (Also see :meth:`~GSASIIscriptable.G2Project.get_Controls` to read values.)
:meth:`~GSASIIscriptable.G2Project.imageMultiDistCalib`                 Performs a global calibration fit with images at multiple distance settings.
:meth:`~GSASIIscriptable.G2Project.get_Constraints`                     Retrieves :ref:`constraint definition <Constraint_definitions_table>` entries.
:meth:`~GSASIIscriptable.G2Project.add_HoldConstr`                      Adds a hold constraint on one or more variables
:meth:`~GSASIIscriptable.G2Project.add_EquivConstr`                     Adds an equivalence constraint on two or more variables
:meth:`~GSASIIscriptable.G2Project.add_EqnConstr`                       Adds an equation-type constraint on two or more variables
:meth:`~GSASIIscriptable.G2Project.add_NewVarConstr`                    Adds an new variable as a constraint on two or more variables
:meth:`~GSASIIscriptable.G2Project.ComputeWorstFit`                     Determines the parameters that will have the greatest impact on the fit if refined
====================================================================    ===============================================================================================================

.. _Class_G2Phase:

Class :class:`~GSASIIscriptable.G2Phase`
-------------------------------------------


  Another common object in GSASIIscriptable scripts is :class:`~GSASIIscriptable.G2Phase`, used to encapsulate each phase in a project, with commonly used methods:

.. tabularcolumns:: |l|p{3.5in}|

========================================================    ===============================================================================================================
method                                                      Use
========================================================    ===============================================================================================================
:meth:`~GSASIIscriptable.G2Phase.set_refinements`           Provides a mechanism to set values and refinement flags for the phase. See :ref:`Phase_parameters_table` 
                                                            for more details. This information also can be supplied within a call
                                                            to :meth:`~GSASIIscriptable.G2Project.do_refinements` 
                                                            or :meth:`~GSASIIscriptable.G2Project.set_refinement`.
:meth:`~GSASIIscriptable.G2Phase.clear_refinements`         Unsets refinement flags for the phase. 
:meth:`~GSASIIscriptable.G2Phase.set_HAP_refinements`       Provides a mechanism to set values and refinement flags for parameters specific to both this phase and 
                                                            one of its histograms. See :ref:`HAP_parameters_table`. This information also can be supplied within 
                                                            a call to :meth:`~GSASIIscriptable.G2Project.do_refinements` or
                                                            :meth:`~GSASIIscriptable.G2Project.set_refinement`.
:meth:`~GSASIIscriptable.G2Phase.clear_HAP_refinements`     Clears refinement flags specific to both this phase and one of its histograms.
:meth:`~GSASIIscriptable.G2Phase.getHAPvalues`              Returns values of parameters specific to both this phase and one of its histograms.
:meth:`~GSASIIscriptable.G2Phase.copyHAPvalues`             Copies HAP settings between from one phase/histogram and to other histograms in same phase.
:meth:`~GSASIIscriptable.G2Phase.HAPvalue`                  Sets or retrieves values for some of the parameters specific to both this phase and 
                                                            one or more of its histograms. 
:meth:`~GSASIIscriptable.G2Phase.atoms`                     Returns a list of atoms in the phase
:meth:`~GSASIIscriptable.G2Phase.atom`                      Returns an atom from its label 
:meth:`~GSASIIscriptable.G2Phase.add_atom`                  Adds an atom to a phase
:meth:`~GSASIIscriptable.G2Phase.histograms`                Returns a list of histograms linked to the phase
:meth:`~GSASIIscriptable.G2Phase.get_cell`                  Returns unit cell parameters (also see :meth:`~GSASIIscriptable.G2Phase.get_cell_and_esd`)
:meth:`~GSASIIscriptable.G2Phase.export_CIF`                Writes a CIF for the phase
:meth:`~GSASIIscriptable.G2Phase.setSampleProfile`          Sets sample broadening parameters
:meth:`~GSASIIscriptable.G2Phase.clearDistRestraint`        Clears any previously defined bond distance restraint(s) for the selected phase
:meth:`~GSASIIscriptable.G2Phase.addDistRestraint`          Finds and defines new bond distance restraint(s) for the selected phase
:meth:`~GSASIIscriptable.G2Phase.setDistRestraintWeight`    Sets the weighting factor for the bond distance restraints
========================================================    ===============================================================================================================

.. _Class_G2PwdrData:

Class :class:`~GSASIIscriptable.G2PwdrData`
---------------------------------------------

  Another common object in GSASIIscriptable scripts is :class:`~GSASIIscriptable.G2PwdrData`, which encapsulate each powder diffraction histogram in a project, with commonly used methods:

.. tabularcolumns:: |l|p{3.5in}|

=======================================================  ===============================================================================================================
method                                                     Use
=======================================================  ===============================================================================================================
:meth:`~GSASIIscriptable.G2PwdrData.set_refinements`     Provides a mechanism to set values and refinement flags for the powder histogram. See 
                                                         :ref:`Histogram_parameters_table` for details.  
:meth:`~GSASIIscriptable.G2PwdrData.clear_refinements`   Unsets refinement flags for the powder histogram.
:meth:`~GSASIIscriptable.G2PwdrData.residuals`           Reports R-factors etc. for the powder histogram (also see :meth:`~GSASIIscriptable.G2PwdrData.get_wR`) 
:meth:`~GSASIIscriptable.G2PwdrData.add_back_peak`       Adds a background peak to the histogram. Also see :meth:`~GSASIIscriptable.G2PwdrData.del_back_peak`
                                                         and :meth:`~GSASIIscriptable.G2PwdrData.ref_back_peak`.
:meth:`~GSASIIscriptable.G2PwdrData.fit_fixed_points`    Fits background to the specified fixed points.
:meth:`~GSASIIscriptable.G2PwdrData.set_background`      Sets a background histogram that will be subtracted (point by point) from the current histogram.
:meth:`~GSASIIscriptable.G2PwdrData.calc_autobkg`        Estimates the background and sets the fixed background points from that.
:meth:`~GSASIIscriptable.G2PwdrData.getdata`             Provides access to the diffraction data associated with the histogram.
:meth:`~GSASIIscriptable.G2PwdrData.reflections`         Provides access to the reflection lists for the histogram.
:meth:`~GSASIIscriptable.G2PwdrData.Export`              Writes the diffraction data or reflection list into a file
:meth:`~GSASIIscriptable.G2PwdrData.add_peak`            Adds a peak to the peak list. Also see :ref:`PeakRefine`.
:meth:`~GSASIIscriptable.G2PwdrData.set_peakFlags`       Sets refinement flags for peaks
:meth:`~GSASIIscriptable.G2PwdrData.refine_peaks`        Starts a peak/background fitting cycle, returns refinement results
:attr:`~GSASIIscriptable.G2PwdrData.Peaks`               Provides access to the peak list data structure
:attr:`~GSASIIscriptable.G2PwdrData.PeakList`            Provides the peak list parameter values 
:meth:`~GSASIIscriptable.G2PwdrData.Export_peaks`        Writes the peak parameters to a text file 
:meth:`~GSASIIscriptable.G2PwdrData.Limits`              Reads or sets the region of data used in fitting (histogram limits)
:meth:`~GSASIIscriptable.G2PwdrData.Excluded`            Reads or sets regions of powder data that will be ignored
=======================================================  ===============================================================================================================

.. _Class_G2Single:

Class :class:`~GSASIIscriptable.G2Single`
---------------------------------------------

  A less commonly-used object in GSASIIscriptable scripts is :class:`~GSASIIscriptable.G2Single`, which will encapsulate each single crystal diffraction
  histogram in a project. At present, very few methods are provided:

.. tabularcolumns:: |l|p{3.5in}|

=======================================================  ===============================================================================================================
method                                                     Use
=======================================================  ===============================================================================================================
:meth:`~GSASIIscriptable.G2Single.set_refinements`        Provides a mechanism to set refinement flags for the single crystal histogram. See 
                                                          :ref:`Histogram_parameters_table` for details.  
:meth:`~GSASIIscriptable.G2Single.clear_refinements`      Unsets refinement flags for the single crystal powder histogram.
:meth:`~GSASIIscriptable.G2Single.Export`                 Writes the reflection list into a file
=======================================================  ===============================================================================================================


.. _Class_G2Image:

Class :class:`~GSASIIscriptable.G2Image`
-----------------------------------------

  When working with images, there will be a :class:`~GSASIIscriptable.G2Image` object for each image (also see :meth:`~GSASIIscriptable.G2Project.add_image`  and :meth:`~GSASIIscriptable.G2Project.images`).

.. tabularcolumns:: |l|p{3.5in}|

====================================================  ===============================================================================================================
method                                                Use
====================================================  ===============================================================================================================
:meth:`~GSASIIscriptable.G2Image.Recalibrate`         Invokes a recalibration fit starting from the current Image Controls calibration coefficients.
:meth:`~GSASIIscriptable.G2Image.Integrate`           Invokes an image integration All parameters Image Controls will have previously been set.
:meth:`~GSASIIscriptable.G2Image.GeneratePixelMask`   Searches for "bad" pixels creating a pixel mask. 
:meth:`~GSASIIscriptable.G2Image.setControl`          Set an Image Controls parameter in the current image. 
:meth:`~GSASIIscriptable.G2Image.getControl`          Return an Image Controls parameter in the current image.
:meth:`~GSASIIscriptable.G2Image.findControl`         Get the names of Image Controls parameters.
:meth:`~GSASIIscriptable.G2Image.loadControls`        Load controls from a .imctrl file (also see :meth:`~GSASIIscriptable.G2Image.saveControls`).
:meth:`~GSASIIscriptable.G2Image.loadMasks`           Load masks from a .immask file.
:meth:`~GSASIIscriptable.G2Image.setVary`             Set a refinement flag for Image Controls parameter in the current image.
                                                      (Also see :meth:`~GSASIIscriptable.G2Image.getVary`)
:meth:`~GSASIIscriptable.G2Image.setCalibrant`        Set a calibrant type (or show choices) for the current image.
:meth:`~GSASIIscriptable.G2Image.setControlFile`      Set a image to be used as a background/dark/gain map image.
:meth:`~GSASIIscriptable.G2Image.getControls`         Returns the Image Controls dict for the current image. 
:meth:`~GSASIIscriptable.G2Image.setControls`         Updates the Image Controls dict for the current image with specified key/value pairs.
:meth:`~GSASIIscriptable.G2Image.getMasks`            Returns the Masks dict for the current image. 
:meth:`~GSASIIscriptable.G2Image.setMasks`            Updates the Masks dict for the current image with specified key/value pairs.
:meth:`~GSASIIscriptable.G2Image.IntThetaAzMap`       Computes the set of 2theta-azimuth mapping matrices to integrate the current image. 
:meth:`~GSASIIscriptable.G2Image.IntMaskMap`          Computes the masking map for the current image for integration. 
:meth:`~GSASIIscriptable.G2Image.MaskThetaMap`        Computes the 2theta mapping matrix to determine a pixel mask. 
:meth:`~GSASIIscriptable.G2Image.MaskFrameMask`       Computes the Frame mask needed to determine a pixel mask. 
:meth:`~GSASIIscriptable.G2Image.TestFastPixelMask`   Returns True if fast pixel masking is available.
:meth:`~GSASIIscriptable.G2Image.clearImageCache`     Clears a saved image from memory, if one is present. 
:meth:`~GSASIIscriptable.G2Image.clearPixelMask`      Clears a saved Pixel map from the project, if one is present. 
====================================================  ===============================================================================================================

.. _Class_G2PDF:

Class :class:`~GSASIIscriptable.G2PDF`
-----------------------------------------

  To work with PDF entries, object :class:`~GSASIIscriptable.G2PDF`, encapsulates a PDF entry with methods:

.. tabularcolumns:: |l|p{3.5in}|

==================================================    ===============================================================================================================
method                                                Use
==================================================    ===============================================================================================================
:meth:`~GSASIIscriptable.G2PDF.export`                                   Used to write G(r), etc. as a file
:meth:`~GSASIIscriptable.G2PDF.calculate`                                Computes the PDF using parameters in the object
:meth:`~GSASIIscriptable.G2PDF.optimize`                                 Optimizes selected PDF parameters
:meth:`~GSASIIscriptable.G2PDF.set_background`                           Sets the histograms used for sample background, container, etc. 
:meth:`~GSASIIscriptable.G2PDF.set_formula`                              Sets the chemical formula for the sample
==================================================    ===============================================================================================================

.. _Class_G2SmallAngle:

Class :class:`~GSASIIscriptable.G2SmallAngle`
-----------------------------------------

  To work with Small Angle (currently only SASD entries), object :class:`~GSASIIscriptable.G2SmallAngle`, encapsulates a SASD entry.
  At present no methods are provided.

.. _Class_G2SeqRefRes:

Class :class:`~GSASIIscriptable.G2SeqRefRes`
-----------------------------------------------

  To work with Sequential Refinement results, object :class:`~GSASIIscriptable.G2SeqRefRes`, encapsulates the sequential refinement table with methods:

.. tabularcolumns:: |l|p{3.5in}|

======================================================    ===============================================================================================================
method                                                    Use
======================================================    ===============================================================================================================
:meth:`~GSASIIscriptable.G2SeqRefRes.histograms`           Provides a list of histograms used in the Sequential Refinement 
:meth:`~GSASIIscriptable.G2SeqRefRes.get_cell_and_esd`     Returns cell dimensions and standard uncertainties for a phase and histogram from the Sequential Refinement 
:meth:`~GSASIIscriptable.G2SeqRefRes.get_Variable`         Retrieves the value and esd for a parameter from a particular histogram in the Sequential Refinement 
:meth:`~GSASIIscriptable.G2SeqRefRes.get_Covariance`       Retrieves values and covariance for a set of refined parameters for a particular histogram 
======================================================    ===============================================================================================================

.. _Class_G2AtomRecord:

Class :class:`~GSASIIscriptable.G2AtomRecord`
-----------------------------------------------

  When working with phases, :class:`~GSASIIscriptable.G2AtomRecord` methods provide access to the contents of each atom in a phase. This provides access to atom
  values via class "properties" that can be used to get values of much of the atoms associated settings, as below. Most can also be used to set values via
  "setter" methods.
  See the :class:`~GSASIIscriptable.G2AtomRecord` docs and source code.

========================================================    ===============================================================================================================
method/property                                                    Use
========================================================    ===============================================================================================================
:data:`~GSASIIscriptable.G2AtomRecord.label`                 Reference as ``<atom>.label``` to get or set label value for atom

:data:`~GSASIIscriptable.G2AtomRecord.type`                  Reference as ``<atom>.G2AtomRecord.type`` to get or set the atom type 

:data:`~GSASIIscriptable.G2AtomRecord.element`               Reference as ``<atom>.G2AtomRecord.element`` to get the element symbol
                                                             associated with an atom (change with ``<atom>.G2AtomRecord.type``,
                                                             see :data:`~GSASIIscriptable.G2AtomRecord.type`)

:data:`~GSASIIscriptable.G2AtomRecord.refinement_flags`      Reference class property ``<atom>.G2AtomRecord.refinement_flags`` to get or set
                                                             the refinement flags associated with an atom
                                                             
:data:`~GSASIIscriptable.G2AtomRecord.coordinates`           Reference as ``<atom>.G2AtomRecord.coordinates`` to get or set the three coordinates
                                                             associated with an atom

:data:`~GSASIIscriptable.G2AtomRecord.occupancy`             Reference class property ``<atom>.G2AtomRecord.occupancy`` to get or set the 
                                                             site occupancy associated with an atom
                                                             
:data:`~GSASIIscriptable.G2AtomRecord.mult`                  Reference as ``<atom>.G2AtomRecord.mult`` to get an atom site multiplicity
                                                             (value cannot be changed in script)

:data:`~GSASIIscriptable.G2AtomRecord.ranId`                 Reference as ``<atom>.G2AtomRecord.ranId`` to get an atom random Id number
                                                             (value cannot be changed in script)
                                                             
:data:`~GSASIIscriptable.G2AtomRecord.adp_flag`              Reference as ``<atom>.G2AtomRecord.adp_flag`` to get either 'U' or 'I'
                                                             specifying that an atom is set as anisotropic or isotropic
                                                             (value cannot be changed in script)

:data:`~GSASIIscriptable.G2AtomRecord.uiso`                  Reference pseudo class variable ``<atom>.G2AtomRecord.uiso`` to get 
                                                             or set the Uiso value associated with an atom
                                                             
========================================================    ===============================================================================================================
  
.. _Refinement_dicts:

Refinement parameters
=====================
While scripts can be written that setup refinements by changing individual parameters 
through calls to the methods associated with objects that wrap each data tree item, 
many of these actions can be combined into fairly complex dict structures to conduct refinement
steps. Use of these dicts is required with the :ref:`CommandlineInterface`. This section of the 
documentation describes these dicts. 

.. _Project_dicts:

Project-level Parameter Dict
----------------------------

As noted below (:ref:`Refinement_parameters_kinds`), there are three types of refinement parameters,
which can be accessed individually by the objects that encapsulate individual phases and histograms
but it will often be simplest to create a composite dictionary
that is used at the project-level. A dict is created with keys
"set" and "clear" that can be supplied to :meth:`~GSASIIscriptable.G2Project.set_refinement`
(or :meth:`~GSASIIscriptable.G2Project.do_refinements`, see :ref:`Refinement_recipe` below) that will
determine parameter values and will determine which parameters will be refined. 

The specific keys and subkeys that can be used are defined in tables 
:ref:`Histogram_parameters_table`, :ref:`Phase_parameters_table` and :ref:`HAP_parameters_table`.

Note that optionally a list of histograms and/or phases can be supplied in the call to 
:meth:`~GSASIIscriptable.G2Project.set_refinement`, but if not specified, the default is to use all defined
phases and histograms. 

As an example: 

.. code-block::  python

    pardict = {'set': { 'Limits': [0.8, 12.0],
                       'Sample Parameters': ['Absorption', 'Contrast', 'DisplaceX'],
                       'Background': {'type': 'chebyschev-1', 'refine': True,
                                      'peaks':[[0,True],[1,1,1]] }},
              'clear': {'Instrument Parameters': ['U', 'V', 'W']}}
    my_project.set_refinement(pardict)

.. _Refinement_recipe:

Refinement recipe
-----------------

Building on the :ref:`Project_dicts`,
it is possible to specify a sequence of refinement actions as a list of
these dicts and supplying this list 
as an argument to :meth:`~GSASIIscriptable.G2Project.do_refinements`.

As an example, this code performs the same actions as in the example in the section above: 

.. code-block::  python

    pardict = {'set': { 'Limits': [0.8, 12.0],
                       'Sample Parameters': ['Absorption', 'Contrast', 'DisplaceX'],
                       'Background': {'type': 'chebyschev-1', 'refine': True}},
              'clear': {'Instrument Parameters': ['U', 'V', 'W']}}
    my_project.do_refinements([pardict])

However, in addition to setting a number of parameters, this example will perform a refinement as well,
after setting the parameters. More than one refinement can be performed by including more 
than one dict in the list. 

In this example, two refinement steps will be performed:

.. code-block::  python

    my_project.do_refinements([pardict,pardict1])


The keys defined in the following table
may be used in a dict supplied to :meth:`~GSASIIscriptable.G2Project.do_refinements`. Note that keys ``histograms``
and ``phases`` are used to limit actions to specific sets of parameters within the project. 

.. tabularcolumns:: |l|p{4in}|

========== ============================================================================
key         explanation
========== ============================================================================
set                    Specifies a dict with keys and subkeys as described in the
                       :ref:`Refinement_parameters_fmt` section. Items listed here
                       will be set to be refined.
clear                  Specifies a dict, as above for set, except that parameters are
                       cleared and thus will not be refined.
once                   Specifies a dict as above for set, except that parameters are
                       set for the next cycle of refinement and are cleared once the
                       refinement step is completed.
skip                   Normally, once parameters are processed with a set/clear/once
                       action(s), a refinement is started. If skip is defined as True
                       (or any other value) the refinement step is not performed.
output                 If a file name is specified for output is will be used to save
                       the current refinement. 
histograms             Should contain a list of histogram(s) to be used for the
                       set/clear/once action(s) on :ref:`Histogram_parameters_table` or
                       :ref:`HAP_parameters_table`. Note that this will be
                       ignored for :ref:`Phase_parameters_table`. Histograms may be
                       specified as a list of strings [('PWDR ...'),...], indices
                       [0,1,2] or as list of objects [hist1, hist2]. 
phases                 Should contain a list of phase(s) to be used for the
                       set/clear/once action(s) on :ref:`Phase_parameters_table` or 
                       :ref:`HAP_parameters_table`. Note that this will be
                       ignored for :ref:`Histogram_parameters_table`.
                       Phases may be specified as a list of strings
                       [('Phase name'),...], indices [0,1,2] or as list of objects
                       [phase0, phase2]. 
call                   Specifies a function to call after a refinement is completed.
                       The value supplied can be the object (typically a function)
                       that will be called or a string that will evaluate (in the
                       namespace inside
		       :meth:`~GSASIIscriptable.G2Project.iter_refinements` where
                       ``self`` references the project.)
                       Nothing is called if this is not specified.
callargs               Provides a list of arguments that will be passed to the function
                       in call (if any). If call is defined and callargs is not, the
                       current <tt>G2Project</tt> is passed as a single argument. 
========== ============================================================================

An example that performs a series of refinement steps follows:

.. code-block::  python

    reflist = [
            {"set": { "Limits": { "low": 0.7 },
                      "Background": { "no. coeffs": 3,
                                      "refine": True }}},
            {"set": { "LeBail": True,
                      "Cell": True }},
            {"set": { "Sample Parameters": ["DisplaceX"]}},
            {"set": { "Instrument Parameters": ["U", "V", "W", "X", "Y"]}},
            {"set": { "Mustrain": { "type": "uniaxial",
                                    "refine": "equatorial",
                                    "direction": [0, 0, 1]}}},
            {"set": { "Mustrain": { "type": "uniaxial",
                                    "refine": "axial"}}},
            {"clear": { "LeBail": True},
             "set": { "Atoms": { "Mn": "X" }}},
            {"set": { "Atoms": { "O1": "X", "O2": "X" }}},]
    my_project.do_refinements(reflist)


In this example, a separate refinement step will be performed for each dict in the list. The keyword 
"skip" can be used to specify a dict that should not include a refinement. 
Note that in the second from last refinement step, parameters are both set and cleared. 

.. _Refinement_parameters_kinds:

Refinement parameter types
--------------------------

Note that parameters and refinement flags used in GSAS-II fall into three classes:

    * **Histogram**: There will be a set of these for each dataset loaded into a
      project file. The parameters available depend on the type of histogram
      (Bragg-Brentano, Single-Crystal, TOF,...). Typical Histogram parameters
      include the overall scale factor, background, instrument and sample parameters;
      see the :ref:`Histogram_parameters_table` table for a list of the histogram
      parameters where access has been provided.

    * **Phase**: There will be a set of these for each phase loaded into a
      project file. While some parameters are found in all types of phases,
      others are only found in certain types (modulated, magnetic, protein...).
      Typical phase parameters include unit cell lengths and atomic positions; see the
      :ref:`Phase_parameters_table` table for a list of the phase      
      parameters where access has been provided.

    * **Histogram-and-phase** (HAP): There is a set of these for every histogram
      that is associated with each phase, so that if there are ``N`` phases and ``M``
      histograms, there can be ``N*M`` total sets of "HAP" parameters sets (fewer if all
      histograms are not linked to all phases.) Typical HAP parameters include the
      phase fractions, sample microstrain and crystallite size broadening terms,
      hydrostatic strain perturbations of the unit cell and preferred orientation
      values.
      See the :ref:`HAP_parameters_table` table for the HAP parameters where access has
      been provided. 

.. _Refinement_parameters_fmt:


Specifying Refinement Parameters
================================

Refinement parameter values and flags to turn refinement on and off are specified within dictionaries,
where the details of these dicts are organized depends on the
type of parameter (see :ref:`Refinement_parameters_kinds`), with a different set
of keys (as described below) for each of the three types of parameters.

.. _Histogram_parameters_table:

Histogram parameters
--------------------

This table describes the dictionaries supplied to :func:`~GSASIIscriptable.G2PwdrData.set_refinements`
and :func:`~GSASIIscriptable.G2PwdrData.clear_refinements`. As an example, 

.. code-block::  python

   hist.set_refinements({"Background": {"no. coeffs": 3, "refine": True},
                         "Sample Parameters": ["Scale"],
                         "Limits": [10000, 40000]})

With :meth:`~GSASIIscriptable.G2Project.do_refinements`, these parameters should be placed inside a dict with a key
``set``, ``clear``, or ``once``. Values will be set for all histograms, unless the ``histograms``
key is used to define specific histograms. As an example: 

.. code-block::  python

  gsas_proj.do_refinements([
      {'set': {
          'Background': {'no. coeffs': 3, 'refine': True},
          'Sample Parameters': ['Scale'],
          'Limits': [10000, 40000]},
      'histograms': [1,2]}
                            ])

Note that below in the Instrument Parameters section, 
related profile parameters (such as U and V) are grouped together but
separated by commas to save space in the table.

.. tabularcolumns:: |l|l|p{3.5in}|

===================== ====================  ==============================================================
key                   subkey                explanation
===================== ====================  ==============================================================
Limits                                      The range of 2-theta (degrees) or TOF (in 
                                            microsec) range of values to use. Can
                                            be either a dictionary of 'low' and/or 'high',
                                            or a list of 2 items [low, high]
                                            Available for powder histograms only.
\                     low                   Sets the low limit
\                     high                  Sets the high limit

Sample Parameters                           Should be provided as a **list** of subkeys
                                            to set or clear refinement flags for,
                                            e.g. ['DisplaceX', 'Scale']
                                            Available for powder histograms only.
\                     Absorption
\                     Contrast
\                     DisplaceX             Sample displacement along the X direction (Debye-Scherrer)
\                     DisplaceY             Sample displacement along the Y direction (Debye-Scherrer)
\                     Shift                 Bragg-Brentano sample displacement 
\                     Scale                 Histogram Scale factor

Background                                  Sample background. Value will be a dict or 
                                            a boolean. If True or False, the refine 
                                            parameter for background is set to that.
                                            Available for powder histograms only.
                                            Note that background peaks are not handled
                                            via this; see 
                                            :meth:`~GSASIIscriptable.G2PwdrData.ref_back_peak` 
                                            instead. When value is a dict,
                                            supply any of the following keys:
\                     type                  The background model, e.g. 'chebyschev-1'
\                     refine                The value of the refine flag, boolean
\                     'no. coeffs'          Number of coefficients to use, integer
\                     coeffs                List of floats, literal values for background
\                     FixedPoints           List of (2-theta, intensity) values for fixed points
\                     'fit fixed points'    If True, triggers a fit to the fixed points to
                                            be calculated. It is calculated when this key is
                                            detected, regardless of calls to refine.
\                     peaks                 Specifies a set of flags for refining 
                                            background peaks as a nested list. There may
                                            be an item for each defined background peak
                                            (or fewer) and each item is a list with the flag 
                                            values for pos,int,sig & gam (fewer than 4 values 
                                            are allowed). 

Instrument Parameters                       As in Sample Parameters, provide as a **list** of
                                            subkeys to set or clear refinement flags,
                                            e.g. ['X', 'Y', 'Zero', 'SH/L']
                                            Available for powder histograms only.
\                     U, V, W               Gaussian peak profile terms
\                     X, Y, Z               Lorentzian peak profile terms
\                     alpha, beta-0,        TOF profile terms 
                      beta-1, beta-q,
\                     sig-0, sig-1,         TOF profile terms
                      sig-2, sig-q
\                     difA, difB, difC      TOF Calibration constants
\                     Zero                  Zero shift
\                     SH/L                  Finger-Cox-Jephcoat low-angle peak asymmetry
\                     Polariz.              Polarization parameter
\                     Lam                   Lambda, the incident wavelength

Single xtal                                 As in Sample Parameters, provide as a **list** of
                                            subkeys to set or clear refinement flags,
                                            e.g. [...].
                                            Available for single crystal histograms only.
\                     Scale                 Single crystal scale factor
\                     BabA, BabU            Babinet A & U parameters
\                     Eg, Es, Ep            Extinction parameters
\                     Flack                 Flack absolute configuration parameter
===================== ====================  ==============================================================

.. _Phase_parameters_table:

Phase parameters
----------------

This table describes the dictionaries supplied to :func:`~GSASIIscriptable.G2Phase.set_refinements`
and :func:`~GSASIIscriptable.G2Phase.clear_refinements`. With :meth:`~GSASIIscriptable.G2Project.do_refinements`,
these parameters should be placed inside a dict with a key
``set``, ``clear``, or ``once``. Values will be set for all phases, unless the ``phases``
key is used to define specific phase(s). 


.. tabularcolumns:: |l|p{4.5in}|

======= ==========================================================
key                   explanation
======= ==========================================================
Cell                  Whether or not to refine the unit cell.
Atoms                 Dictionary of atoms and refinement flags.
                      Each key should be an atom label, e.g.
                      'O3', 'Mn5', and each value should be
                      a string defining what values to refine.
                      Values can be any combination of 'F'
                      for site fraction, 'X' for position,
                      and 'U' for Debye-Waller factor
LeBail                Enables LeBail intensity extraction.
======= ==========================================================


.. _HAP_parameters_table:


Histogram-and-phase parameters
------------------------------

This table describes the dictionaries supplied to :func:`~GSASIIscriptable.G2Phase.set_HAP_refinements`
and :func:`~GSASIIscriptable.G2Phase.clear_HAP_refinements`. When supplied to
:meth:`~GSASIIscriptable.G2Project.do_refinements`, these parameters should be placed inside a dict with a key
``set``, ``clear``, or ``once``. Values will be set for all histograms used in each phase,
unless the ``histograms`` and ``phases`` keys are used to define specific phases and histograms.

.. tabularcolumns:: |l|l|p{3.5in}|

=============  ==========  =========================================================================
key             subkey                 explanation
=============  ==========  =========================================================================
Babinet                                Should be a **list** of the following
                                       subkeys. If not, assumes both
                                       BabA and BabU
\              BabA
\              BabU
Extinction                             Boolean, True to refine.
HStrain                                Boolean or list/tuple, True to refine all 
                                       appropriate D\ :sub:`ij` terms or False
                                       to not refine any. If a list/tuple, will
                                       be a set of True & False values for each 
                                       D\ :sub:`ij` term; number of items must 
                                       match number of terms.
Mustrain
\              type                   Mustrain model. One of 'isotropic',
                                      'uniaxial', or 'generalized'. This should
                                      be specified to change the model.
\              direction              For uniaxial only. A list of three
                                      integers,
                                      the [hkl] direction of the axis.
\              refine                 Usually boolean, set to True to refine.
                                      or False to clear. 
                                      For uniaxial model, can specify a value
                                      of 'axial' or 'equatorial' to set that flag
                                      to True or a single
                                      boolean sets both axial and equatorial.
Size                                   
\              type                   Size broadening model. One of 'isotropic',
                                      'uniaxial', or 'ellipsoid'. This should 
                                      be specified to change from the current.
\              direction              For uniaxial only. A list of three
                                      integers,
                                      the [hkl] direction of the axis.
\              refine                 Boolean, True to refine.
\              value                  float, size value in microns 
Pref.Ori.                             Boolean, True to refine
Show                                  Boolean, True to refine
Use                                   Boolean, True to refine
Scale                                 Phase fraction; Boolean, True to refine
PhaseFraction                         PhaseFraction can also be used in place of
                                      Scale for the routines that access HAP
                                      parameters:
                                      :func:`~GSASIIscriptable.G2Phase.HAPvalue`,
                                      :func:`~GSASIIscriptable.G2Phase.setHAPvalues`,
                                      :func:`~GSASIIscriptable.G2Phase.copyHAPvalues`,
                                      :meth:`~GSASIIscriptable.G2Project.set_refinement`,
                                      :meth:`~GSASIIscriptable.G2Project.do_refinements`,
                                      :func:`~GSASIIscriptable.G2Phase.clear_HAP_refinements`
                                      and :func:`~GSASIIscriptable.G2Phase.set_HAP_refinements`.
=============  ==========  =========================================================================

Histogram/Phase objects
-----------------------
Each phase and powder histogram in a :class:`~GSASIIscriptable.G2Project` object has an associated
object. Parameters within each individual object can be turned on and off by calling
:meth:`~GSASIIscriptable.G2PwdrData.set_refinements` or :meth:`~GSASIIscriptable.G2PwdrData.clear_refinements`
for histogram parameters;
:meth:`~GSASIIscriptable.G2Phase.set_refinements` or :meth:`~GSASIIscriptable.G2Phase.clear_refinements`
for phase parameters; and :meth:`~GSASIIscriptable.G2Phase.set_HAP_refinements` or
:meth:`~GSASIIscriptable.G2Phase.clear_HAP_refinements`. As an example, if some_histogram is a histogram object (of type :class:`~GSASIIscriptable.G2PwdrData`), use this to set parameters in that histogram:

.. code-block::  python

    params = { 'Limits': [0.8, 12.0],
               'Sample Parameters': ['Absorption', 'Contrast', 'DisplaceX'],
               'Background': {'type': 'chebyschev-1', 'refine': True}}
    some_histogram.set_refinements(params)

Likewise to turn refinement flags on, use code such as this:

.. code-block::  python

    params = { 'Instrument Parameters': ['U', 'V', 'W']} 
    some_histogram.set_refinements(params)

and to turn these refinement flags, off use this (Note that the
``.clear_refinements()`` methods will usually will turn off refinement even
if a refinement parameter is set in the dict to True.):

.. code-block::  python

    params = { 'Instrument Parameters': ['U', 'V', 'W']} 
    some_histogram.clear_refinements(params)

For phase parameters, use code such as this:

.. code-block::  python

    params = { 'LeBail': True, 'Cell': True,
               'Atoms': { 'Mn1': 'X',
                          'O3': 'XU',
                          'V4': 'FXU'}}
    some_histogram.set_refinements(params)

and here is an example for HAP parameters:

.. code-block::  python

    params = { 'Babinet': 'BabA',
               'Extinction': True,
               'Mustrain': { 'type': 'uniaxial',
                             'direction': [0, 0, 1],
                             'refine': True}}
    some_phase.set_HAP_refinements(params)

Note that the parameters must match the object type and method (phase vs. histogram vs. HAP).

.. _AccessingOtherItems:

Access to other parameter settings
==================================

There are several hundred different types of values that can be stored in a 
GSAS-II project (.gpx) file. All can be changed from the GUI but only a 
subset have direct mechanism implemented for change from the GSASIIscriptable 
API. In practice all parameters in a .gpx file can be edited via scripting, 
but sometimes determining what should be set to implement a parameter 
change can be complex. 
Several routines, :meth:`~GSASIIscriptable.G2Phase.getHAPentryList`, 
:meth:`~GSASIIscriptable.G2Phase.getPhaseEntryList` and :meth:`~GSASIIscriptable.G2PwdrData.getHistEntryList` 
(and their related get...Value and set.Value entries), 
provide a mechanism to discover what the GUI is changing inside a .gpx file. 

As an example, a user in changing the data type for a histogram from Debye-Scherrer 
mode to Bragg-Brentano. This capability is not directly exposed in the API. To 
find out what changes when the histogram type is changed we can create a short script 
that displays the contents of all the histogram settings:

.. code-block::  python

    from __future__ import division, print_function
    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    from . import GSASIIscriptable as G2sc
    gpx = G2sc.G2Project('/tmp/test.gpx')
    h = gpx.histograms()[0]
    for h in h.getHistEntryList():
        print(h)

This can be run with a command like this::

       python test.py > before.txt

(This will create file ``before.txt``, which will contain hundreds of lines.) 

At this point open the project file, ``test.gpx`` in the GSAS-II GUI and 
change in Histogram/Sample Parameters the diffractometer type from Debye-Scherrer 
mode to Bragg-Brentano and then save the file. 

Rerun the previous script creating a new file::

       python test.py > after.txt

Finally look for the differences between files ``before.txt`` and ``after.txt`` using a tool 
such as diff (on Linux/OS X) or fc (in Windows).

in Windows:: 

    Z:\>fc before.txt after.txt
    Comparing files before.txt and after.txt
    ***** before.txt
           fill_value = 1e+20)
    , 'PWDR Co_PCP_Act_d900-00030.fxye Bank 1', 'PWDR Co_PCP_Act_d900-00030.fxye Ban
    k 1'])
    (['Comments'], <class 'list'>, ['Co_PCP_Act_d900-00030.tif #0001 Azm= 180.00'])
    ***** AFTER.TXT
           fill_value = 1e+20)
    , 'PWDR Co_PCP_Act_d900-00030.fxye Bank 1', 'PWDR Co_PCP_Act_d900-00030.fxye Ban
    k 1', 'PWDR Co_PCP_Act_d900-00030.fxye Bank 1']

    (['Comments'], <class 'list'>, ['Co_PCP_Act_d900-00030.tif #0001 Azm= 180.00'])
    *****

    ***** before.txt
    (['Sample Parameters', 'Scale'], <class 'list'>, [1.276313196832068, True])
    (['Sample Parameters', 'Type'], <class 'str'>, 'Debye-Scherrer')
    (['Sample Parameters', 'Absorption'], <class 'list'>, [0.0, False])
    ***** AFTER.TXT
    (['Sample Parameters', 'Scale'], <class 'list'>, [1.276313196832068, True])
    (['Sample Parameters', 'Type'], <class 'str'>, 'Bragg-Brentano')
    (['Sample Parameters', 'Absorption'], <class 'list'>, [0.0, False])
    *****

in Linux/Mac:: 

    bht14: toby$ diff before.txt after.txt 
    103c103
    < , 'PWDR Co_PCP_Act_d900-00030.fxye Bank 1', 'PWDR Co_PCP_Act_d900-00030.fxye Bank 1'])
    ---
    > , 'PWDR Co_PCP_Act_d900-00030.fxye Bank 1', 'PWDR Co_PCP_Act_d900-00030.fxye Bank 1', 'PWDR Co_PCP_Act_d900-00030.fxye Bank 1'])
    111c111
    < (['Sample Parameters', 'Type'], <class 'str'>, 'Debye-Scherrer')
    ---
    > (['Sample Parameters', 'Type'], <class 'str'>, 'Bragg-Brentano')

From this we can see there are two changes that took place. One is fairly obscure, 
where the histogram name is added to a list, which can be ignored, but the second change 
occurs in a straight-forward way and we discover that a simple call::

    h.setHistEntryValue(['Sample Parameters', 'Type'], 'Bragg-Brentano')

can be used to change the histogram type. 

.. _CodeExamples:

Code Examples
=============

.. contents:: Contents for Scripting Examples
   :local: 

.. _ScriptingShortcut:

Shortcut for Scripting Access
-----------------------------

As is seen in a number of the code examples, the location where GSAS-II is
specified in the GSAS-II script using commands such as 

.. code-block::  python

    import sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII') # needed to "find" GSAS-II modules
    from . import GSASIIscriptable as G2sc

An alternative to this is to "install" the current GSAS-II installation into the current
Python interpreter. Once this has been done a single time, this single command can be used to replace
the three commands listed above for all future uses of GSASIIscripting:

.. code-block::  python

    import G2script as G2sc

There are two ways this installation can be done. The most easy way is to invoke the
"Install GSASIIscriptable shortcut" command in the GSAS-II GUI
File menu. Alternatively it can be accomplished from within GSASIIscriptable
using these commands:

.. code-block::  python

    import sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII') # update this for your installation
    from . import GSASIIscriptable as G2sc
    G2sc.installScriptingShortcut()

An even simpler way to do this is from the command-line, from the GSAS-II directory.
A full path for Python is only needed if if the Python to be used with GSAS-II is not in the
path. 

.. code-block::  bash

		 terrier:toby> cd /home/beams1/TOBY/gsas2full/GSASII/
		 terrier:toby> /mypath/bin/python -c "import GSASIIscriptable as G2sc; G2sc.installScriptingShortcut()"
		 GSAS-II binary directory: /home/beams1/TOBY/gsas2full/GSASII/bindist
		 Created file /home/beams1/TOBY/gsas2full/lib/python3.10/site-packages/G2script.py
		 setting up GSASIIscriptable from /home/beams1/TOBY/gsas2full/GSASII
		 success creating /home/beams1/TOBY/gsas2full/lib/python3.10/site-packages/G2script.py

Note the shortcut only installs use of GSAS-II with the current Python
installation. If more than one Python installation will be used with GSAS-II
(for example because different conda environments are used), a shortcut
should be created from within each Python environment.

If more than one GSAS-II installation will be used with a Python installation, 
a shortcut can only be used with one of them.

Status Information
-----------------------------

To find information on Python, Python packages and the GSAS-II version, one can call the
:func:`~GSASIIscriptable.ShowVersions` function. This will show versions and
install locations. 

.. code-block::  python

    import G2script as G2sc
    print(f'Version information:\n{G2sc.ShowVersions()}')

which produces output like this::

  setting up GSASIIscriptable from /Users/toby/G2/git/g2full/GSAS-II/GSASII
  Version information:
    Python      3.11.9:  from /Users/toby/py/mf3/envs/py311/bin/python
    numpy       1.26.4:  
    scipy       1.13.0:  
    IPython     8.22.2:  
    GSAS-II:    641a65, 24-May-2024 10:16 (0.5 days old). Last tag: #5789

  GSAS-II location: /Users/toby/G2/git/g2full/GSAS-II/GSASII
  Binary location:  /Users/toby/G2/git/g2full/GSAS-II/GSASII-bin/mac_arm_p3.11_n1.26

    
.. _PeakRefine:  
 
Peak Fitting
------------

Peak refinement is performed with routines 
:meth:`~GSASIIscriptable.G2PwdrData.add_peak`, :meth:`~GSASIIscriptable.G2PwdrData.set_peakFlags` and
:meth:`~GSASIIscriptable.G2PwdrData.refine_peaks`. Method :meth:`~GSASIIscriptable.G2PwdrData.Export_peaks` and
properties :attr:`~GSASIIscriptable.G2PwdrData.Peaks` and :attr:`~GSASIIscriptable.G2PwdrData.PeakList` 
provide ways to access the results. Note that when peak parameters are 
refined with :meth:`~GSASIIscriptable.~G2PwdrData.refine_peaks`, the background may also
be refined. Use :meth:`~GSASIIscriptable.G2PwdrData.set_refinements` to change background 
settings and the range of data used in the fit. See below for an example
peak refinement script, where the data files are taken from the 
"Rietveld refinement with CuKa lab Bragg-Brentano powder data" tutorial 
(in https://advancedphotonsource.github.io/GSAS-II-tutorials/LabData/data/).

.. code-block::  python

    from __future__ import division, print_function
    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII') # needed to "find" GSAS-II modules
    from . import GSASIIscriptable as G2sc
    datadir = os.path.expanduser("~/Scratch/peakfit")
    PathWrap = lambda fil: os.path.join(datadir,fil)
    gpx = G2sc.G2Project(newgpx=PathWrap('pkfit.gpx'))
    hist = gpx.add_powder_histogram(PathWrap('FAP.XRA'), PathWrap('INST_XRY.PRM'),
                                    fmthint='GSAS powder')
    hist.set_refinements({'Limits': [16.,24.],
          'Background': {"no. coeffs": 2,'type': 'chebyschev-1', 'refine': True}
                         })
    peak1 = hist.add_peak(1, ttheta=16.8)
    peak2 = hist.add_peak(1, ttheta=18.9)
    peak3 = hist.add_peak(1, ttheta=21.8)
    peak4 = hist.add_peak(1, ttheta=22.9)
    hist.set_peakFlags(area=True)
    hist.refine_peaks()
    hist.set_peakFlags(area=True,pos=True)
    hist.refine_peaks()
    hist.set_peakFlags(area=True, pos=True, sig=True, gam=True)
    res = hist.refine_peaks()
    print('peak positions: ',[i[0] for i in hist.PeakList])
    for i in range(len(hist.Peaks['peaks'])):
        print('peak',i,'pos=',hist.Peaks['peaks'][i][0],'sig=',hist.Peaks['sigDict']['pos'+str(i)])
    hist.Export_peaks('pkfit.txt')
    #gpx.save()  # gpx file is not written without this

Pattern Simulation
------------------

This shows two examples where a structure is read from a CIF, a 
pattern is computed using a instrument parameter file to specify the 
probe type (neutrons here) and wavelength. 

The first example uses a CW neutron instrument parameter file. 
The pattern is computed over a 2θ range of 5 to 120 degrees 
with 1000 points. 
The pattern and reflection list are written into files. 
Data files are found in the 
`Scripting Tutorial <https://advancedphotonsource.github.io/GSAS-II-tutorials/PythonScript/data/>`_.

.. code-block::  python

    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    from . import GSASIIscriptable as G2sc
    datadir = "/Users/toby/software/G2/Tutorials/PythonScript/data"
    PathWrap = lambda fil: os.path.join(datadir,fil)
    gpx = G2sc.G2Project(newgpx='PbSO4sim.gpx') # create a project    
    phase0 = gpx.add_phase(PathWrap("PbSO4-Wyckoff.cif"),
             phasename="PbSO4",fmthint='CIF') # add a phase to the project
    # add a simulated histogram and link it to the previous phase(s)
    hist1 = gpx.add_simulated_powder_histogram("PbSO4 simulation",
                PathWrap("inst_d1a.prm"),5.,120.,Npoints=1000,
                phases=gpx.phases(),scale=500000.)
    gpx.do_refinements()   # calculate pattern
    gpx.save()
    # save results
    gpx.histogram(0).Export('PbSO4data','.csv','hist') # data
    gpx.histogram(0).Export('PbSO4refl','.csv','refl') # reflections

This example uses bank#2 from a TOF neutron instrument parameter file. 
The pattern is computed over a TOF range of 14 to 35 milliseconds with 
the default of 2500 points. 
This uses the same CIF as in the example before, but the instrument is found in the  
`TOF-CW Joint Refinement Tutorial <https://advancedphotonsource.github.io/GSAS-II-tutorials/TOF-CW Joint Refinement/data>`_
tutorial. 

.. code-block::  python

    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    from . import GSASIIscriptable as G2sc
    cifdir = "/Users/toby/software/G2/Tutorials/PythonScript/data"
    datadir = "/Users/toby/software/G2/Tutorials/TOF-CW Joint Refinement/data"
    gpx = G2sc.G2Project(newgpx='/tmp/PbSO4simT.gpx') # create a project
    phase0 = gpx.add_phase(os.path.join(cifdir,"PbSO4-Wyckoff.cif"),
             phasename="PbSO4",fmthint='CIF') # add a phase to the project
    hist1 = gpx.add_simulated_powder_histogram("PbSO4 simulation",
                os.path.join(datadir,"POWGEN_1066.instprm"),14.,35.,
                phases=gpx.phases(),ibank=2)
    gpx.do_refinements([{}])
    gpx.save()

Simple Refinement
-----------------

GSASIIscriptable can be used to setup and perform simple refinements. 
This example reads in an existing project (.gpx) file, adds a background
peak, changes some refinement flags and performs a refinement.

.. code-block::  python

    from __future__ import division, print_function
    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII') # needed to "find" GSAS-II modules
    from . import GSASIIscriptable as G2sc
    datadir = "/Users/Scratch/"
    gpx = G2sc.G2Project(os.path.join(datadir,'test2.gpx'))
    gpx.histogram(0).add_back_peak(4.5,30000,5000,0)
    pardict = {'set': {'Sample Parameters': ['Absorption', 'Contrast', 'DisplaceX'],
                       'Background': {'type': 'chebyschev-1', 'refine': True,
                                      'peaks':[[0,True]]}}}
    gpx.set_refinement(pardict)

Sequential Refinement
---------------------

GSASIIscriptable can be used to setup and perform sequential refinements. This example script 
is used to take the single-dataset fit at the end of Step 1 of the 
`Sequential Refinement tutorial <https://advancedphotonsource.github.io/GSAS-II-tutorials/SeqRefine/SequentialTutorial.htm>`_
and turn on and off refinement flags, add histograms and setup the sequential fit, which is then run:

.. code-block::  python

    import os,sys,glob
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    from . import GSASIIscriptable as G2sc
    datadir = os.path.expanduser("~/Scratch/SeqTut2019Mar")
    PathWrap = lambda fil: os.path.join(datadir,fil)
    # load and rename project
    gpx = G2sc.G2Project(PathWrap('7Konly.gpx'))
    gpx.save(PathWrap('SeqRef.gpx'))
    # turn off some variables; turn on Dijs
    for p in gpx.phases():
        p.set_refinements({"Cell": False})
    gpx.phase(0).set_HAP_refinements(
        {'Scale': False,
         "Size": {'type':'isotropic', 'refine': False},
         "Mustrain": {'type':'uniaxial', 'refine': False},
         "HStrain":True,})
    gpx.phase(1).set_HAP_refinements({'Scale': False})
    gpx.histogram(0).clear_refinements({'Background':False,
                     'Sample Parameters':['DisplaceX'],})
    gpx.histogram(0).ref_back_peak(0,[])
    gpx.phase(1).set_HAP_refinements({"HStrain":(1,1,1,0)})
    for fil in sorted(glob.glob(PathWrap('*.fxye'))): # load in remaining fxye files
        if '00' in fil: continue
        gpx.add_powder_histogram(fil, PathWrap('OH_00.prm'), fmthint="GSAS powder",phases='all')
    # copy HAP values, background, instrument params. & limits, not sample params. 
    gpx.copyHistParms(0,'all',['b','i','l']) 
    for p in gpx.phases(): p.copyHAPvalues(0,'all')
    # setup and launch sequential fit
    gpx.set_Controls('sequential',gpx.histograms())
    gpx.set_Controls('cycles',10)
    gpx.set_Controls('seqCopy',True)
    gpx.refine()  

.. _ImageProc:

Image Processing
----------------

A sample script where an image is read, assigned calibration values from a file 
and then integrated follows. 
The data files are found in the 
`Scripting Tutorial <https://advancedphotonsource.github.io/GSAS-II-tutorials/PythonScript/data/>`_.

.. code-block::  python

    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    from . import GSASIIscriptable as G2sc
    datadir = "/tmp"
    PathWrap = lambda fil: os.path.join(datadir,fil)

    gpx = G2sc.G2Project(newgpx=PathWrap('inttest.gpx'))
    imlst = gpx.add_image(PathWrap('Si_free_dc800_1-00000.tif'),fmthint="TIF")
    imlst[0].loadControls(PathWrap('Si_free_dc800_1-00000.imctrl'))
    pwdrList = imlst[0].Integrate()
    gpx.save()

This example shows a computation similar to what is done in tutorial 
`Area Detector Calibration with Multiple Distances <https://advancedphotonsource.github.io/GSAS-II-tutorials/DeterminingWavelength/DeterminingWavelength.html>`_

.. code-block::  python

    import os,sys,glob
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    from . import GSASIIscriptable as G2sc
    PathWrap = lambda fil: os.path.join(
        "/Users/toby/wp/Active/MultidistanceCalibration/multimg",
        fil)

    gpx = G2sc.G2Project(newgpx='/tmp/img.gpx')
    for f in glob.glob(PathWrap('*.tif')): 
        im = gpx.add_image(f,fmthint="TIF")
    # image parameter settings
    defImgVals = {'wavelength': 0.24152, 'center': [206., 205.],
      'pixLimit': 2,  'cutoff': 5.0, 'DetDepth': 0.055,'calibdmin': 1.,}
    # set controls and vary options, then fit
    for img in gpx.images():
        img.setCalibrant('Si    SRM640c')
        img.setVary('*',False)
        img.setVary(['det-X', 'det-Y', 'phi', 'tilt', 'wave'], True)
        img.setControls(defImgVals)
        img.Recalibrate()
        img.Recalibrate() # 2nd run better insures convergence
    gpx.save()
    # make dict of images for sorting
    images = {img.getControl('setdist'):img for img in gpx.images()}
    # show values
    for key in sorted(images.keys()):
        img = images[key]
        c = img.getControls()
        print(c['distance'],c['wavelength'])

.. _MultiDist_Example:

Image Calibration
-----------------

This example performs a number of cycles of constrained fitting. 
A project is created with the images found in a directory, setting initial
parameters as the images are read. The initial values 
for the calibration are not very good, so a :meth:`~GSASIIscriptable.G2Image.Recalibrate` is done
to quickly improve the fit. Once that is done, a fit of all images is performed
where the wavelength, an offset and detector orientation are constrained to 
be the same for all images. The detector penetration correction is then added. 
Note that as the calibration values improve, the algorithm is able to find more 
points on diffraction rings to use for calibration and the number of "ring picks" 
increase. The calibration is repeated until that stops increasing significantly (<10%). 
Detector control files are then created. 
The files used for this exercise are found in the
`Area Detector Calibration Tutorial <https://advancedphotonsource.github.io/GSAS-II-tutorials/DeterminingWavelength/data/>`_
(see 
`Area Detector Calibration with Multiple Distances <https://advancedphotonsource.github.io/GSAS-II-tutorials/DeterminingWavelength/DeterminingWavelength.html>`_ ). 

.. code-block::  python

    import os,sys,glob
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    from . import GSASIIscriptable as G2sc
    PathWrap = lambda fil: os.path.join(
        "/Users/toby/wp/Active/MultidistanceCalibration/multimg",
        fil)

    gpx = G2sc.G2Project(newgpx='/tmp/calib.gpx')
    for f in glob.glob(PathWrap('*.tif')): 
        im = gpx.add_image(f,fmthint="TIF")
    # starting image parameter settings
    defImgVals = {'wavelength': 0.240, 'center': [206., 205.],
      'pixLimit': 2,  'cutoff': 5.0, 'DetDepth': 0.03,'calibdmin': 0.5,}
    # set controls and vary options, then initial fit
    for img in gpx.images():
        img.setCalibrant('Si    SRM640c')
        img.setVary('*',False)
        img.setVary(['det-X', 'det-Y', 'phi', 'tilt', 'wave'], True)
        img.setControls(defImgVals)
        if img.getControl('setdist') > 900:
            img.setControls({'calibdmin': 1.,})
        img.Recalibrate()
    G2sc.SetPrintLevel('warn') # cut down on output
    result,covData = gpx.imageMultiDistCalib()
    print('1st global fit: initial ring picks',covData['obs'])
    print({i:result[i] for i in result if '-' not in i})
    # add parameter to all images & refit multiple times
    for img in gpx.images(): img.setVary('dep',True)
    ringpicks = covData['obs']
    delta = ringpicks
    while delta > ringpicks/10:
        result,covData = gpx.imageMultiDistCalib(verbose=False)
        delta = covData['obs'] - ringpicks
        print('ring picks went from',ringpicks,'to',covData['obs'])
        print({i:result[i] for i in result if '-' not in i})
        ringpicks = covData['obs']
    # once more for good measure & printout
    result,covData = gpx.imageMultiDistCalib(verbose=True)
    # create image control files
    for img in gpx.images():
        img.saveControls(os.path.splitext(img.name)[0]+'.imctrl')
    gpx.save()

.. _OptInteg_Example:

Optimized Image Integration
---------------------------

This example shows how image integration, including pixel masking of outliers,
can be accomplished for a series of images where the calibration and 
other masking (Frame, Spots, etc) are the same for all images. This code has been optimized
significantly so that computations are cached and are not repeated where possible. For one
set of test data, processing of the first image takes ~5 seconds, but processing of subsequent
takes on the order of 0.7 sec. 

This code uses an ``import G2script as G2sc`` statement to access GSASIIscriptable
without referencing the GSAS-II installation directory. This requires installing a reference to
the GSAS-II location into the current a Python installation, which can be done from the GUI
or with scripting commands, as is discussed in :ref:`ScriptingShortcut`. Here 
function :func:`~GSASIIscriptable.installScriptingShortcut` was used to create
the :mod:`G2script` module. That code has been retained here as comments to show what was done. 

To simplify use of this script, it is assumed that the script will be placed in the same
directory as where the data files will be collected. Other customization is done
in variables at the beginning of the code. Note that the beamline where these data are collected
opens the output .tif files before the data collection for that image is complete. Once the .metadata
file has been created, the image may be read.

Processing progresses as follows:
 * Once a set of images are found, a project is created. This is never written and will be deleted after the images are processed.
 * For each image file, routine :meth:`~GSASIIscriptable.G2Project.add_image` is used to add image(s) from that file to the project. The .tif format can only hold one image, but others can have more than one. 
 * When the first image is processed, calibration and mask info is read; a number of computations are performed and cached.
 * For subsequent images cached information is used.
 * Pixel masking is performed in :meth:`~GSASIIscriptable.G2Image.GeneratePixelMask` and the mask is saved into the image.
 * Image integration is performed in :meth:`~GSASIIscriptable.G2Image.Integrate`.
 * Note that multiple powder patterns could be created from one image, so creation of data files is done in a loop with :meth:`~GSASIIscriptable.G2PwdrData.Export`.
 * To reduce memory demands, cached versions of the Pixel map and the Image are deleted and the image file is moved to a separate directory so note that it has been processed. 
 * The project (.gpx file) is deleted and recreated periodically so that the memory footprint for this script does not grow.

The speed of this code will depend on many things, but the number of pixels in the
image is primary, as well as CPU speed. With ~9 Mb images, I have seen average times in the range of 0.7 to 0.9 sec/image, after the first image is processed and the cached arrays are computed. With the Apple M1 chip the time is closer to 0.6 sec/image.
There is also a possible tuning parameter that may change speed based on the speed of the CPU vs. memory
constraints in variable :data:`GSASIIscriptable.blkSize`. This value should be a power of two and defaults to
128. You might find that a larger or smaller value will improve performance for you. 
   
.. code-block::  python

    import os,glob,time,shutil

    #### Create G2script: do this once ################################################
    #import sys
    #sys.path.insert(0,'/Users/toby/software/G2/GSASII') # update with your install loc
    #import GSASIIscriptable as G2sc
    #G2sc.installScriptingShortcut()
    ###################################################################################

    import G2script as G2sc
    G2sc.blkSize = 2**8  # computer-dependent tuning parameter
    G2sc.SetPrintLevel('warn')   # reduces output

    cache = {}  # place to save intermediate computations
    # define location & names of files
    dataLoc = os.path.abspath(os.path.split(__file__)[0]) # data in location of this file
    PathWrap = lambda fil: os.path.join(dataLoc,fil) # convenience function for file paths 
    imgctrl = PathWrap('Si_ch3_d700-00000.imctrl')
    imgmask = PathWrap('Si_ch3_d700-00000.immask')
    globPattern = PathWrap("*_d700-*.tif")

    def wait_for_metadata(tifname):
        '''A .tif file is created before it can be read. Wait for the
        metadata file to be created before trying to read both.
        '''
        while not os.path.exists(tifname + '.metadata'):
            time.sleep(0.05)

    # make a subfolder to store integrated images & integrated patterns
    pathImg = os.path.join(dataLoc,'img') 
    if not os.path.exists(pathImg): os.mkdir(pathImg)
    pathxye = os.path.join(dataLoc,'xye')
    if not os.path.exists(pathxye): os.mkdir(pathxye)

    while True:	 # Loop will never end, stop with ctrl+C
        tiflist = sorted(glob.glob(globPattern),key=lambda x: os.path.getctime(x)) # get images sorted by creation time, oldest 1st
        if not tiflist:
            time.sleep(0.1)
            continue
        gpx = G2sc.G2Project(newgpx=PathWrap('integration.gpx')) # temporary use
        for tifname in tiflist:
            starttime = time.time()
            wait_for_metadata(tifname)
            for img in gpx.add_image(tifname,fmthint="TIF",cacheImage=True):  # loop unneeded for TIF (1 image/file)
                if not cache: # load & compute controls & 2theta values once
                    img.loadControls(imgctrl)	# set controls/calibrations/masks
                    img.loadMasks(imgmask)
                    cache['Image Controls'] = img.getControls() # save controls & masks contents for quick reload
                    cache['Masks'] = img.getMasks()
                    cache['intMaskMap'] = img.IntMaskMap() # calc mask & TA arrays to save for integrations
                    cache['intTAmap'] = img.IntThetaAzMap()
                    cache['FrameMask'] = img.MaskFrameMask() # calc Frame mask & T array to save for Pixel masking
                    cache['maskTmap'] = img.MaskThetaMap() 
                else:
                    img.setControls(cache['Image Controls'])
                    img.setMasks(cache['Masks'],True)  # True: reset threshold masks
                img.GeneratePixelMask(esdMul=3,ThetaMap=cache['maskTmap'],FrameMask=cache['FrameMask'])
                for pwdr in img.Integrate(MaskMap=cache['intMaskMap'],ThetaAzimMap=cache['intTAmap']):
                    pwdr.Export(os.path.join(pathxye,os.path.split(tifname)[1]),'.xye')  # '.tif in name ignored         
                img.clearImageCache()  # save some space
                img.clearPixelMask()
            shutil.move(tifname, pathImg)	# move file after integration so that it is not searchable
            shutil.move(tifname + '.metadata', pathImg)
            print('*=== processing complete, time=',time.time()-starttime,'sec\n')
        del gpx
	
.. _MultiCoreInteg_Example:

Multicore Image Integration
---------------------------


The previous example (:ref:`OptInteg_Example`) can be accelerated even further
on a multicore computer using the following script. In this example, 
the image integration is moved to a function, `integrate_tif`, that accepts
a filename to integrate. Note that with the multiprocessing module is used, 
the script will be read on each core that will be used, but only on the primary 
(controller) process will this ``__name__ == '__main__'`` be True. 
Thus the code following the if statement runs on the primary process. 
The primary process uses the mp.Pool() statement to create a set of 
secondary (worker) processes that are intended to run on other cores. 
The primary process locates .tif files, if the corresponding 
.tif.metadata is also found, both are moved to a separate directory where they 
will be processed in a secondary process. When the secondary process starts, 
the script is imported and then `integrate_tif` is called with the name of the 
image file from the primary process. The `integrate_tif` routine 
will initially have an empty cache and thus the code preceeded by 
"load & compute controls & 2theta values" will be computed once for every 
secondary process, which should be on an independent core. The size of the pool
determines how many images will be processed simultaneously. 

The script as given below uses the first argument on the command 
line to specify the number of cores to be used, where 0 is used to 
mean run `integrate_tif` directly rather than through a pool. This 
facilitates timing comparisons. 
This code seems to have a maximum speed using slightly less than the 
total number of available cores and does benefit partially from 
hyperthreading. A two- to three-fold speedup is seen with four cores and a 
six-fold speedup has been seen with 16 cores. 

.. code-block::  python

    import os,sys,glob,time,shutil
    scriptstart = time.time()

    if len(sys.argv) >= 2:
        nodes = int(sys.argv[1])
    else:
        nodes = 4

    if nodes == 0:
        print('no multiprocessing')
    else:
        print(f'multiprocessing with {nodes} cores')

    import G2script as G2sc
    G2sc.blkSize = 2**8  # computer-dependent tuning parameter
    #G2sc.SetPrintLevel('warn')

    cache = {}  # place to save intermediate computations

    # define location & names of files
    dataLoc = '/dataserv/inttest'  # images found here
    globPattern = os.path.join(dataLoc,"*_d700-*.tif")
    calibLoc = os.path.abspath(os.path.split(__file__)[0]) # calib in location of this file
    imgctrl = os.path.join(calibLoc,'Si_ch3_d700-00000.imctrl')
    imgmask = os.path.join(calibLoc,'Si_ch3_d700-00000.immask')
    # locations to put processed files
    pathImg = os.path.join(dataLoc,'img') 
    pathxye = os.path.join(dataLoc,'xye')

    def integrate_tif(tifname):
        starttime = time.time()
        gpx = G2sc.G2Project(newgpx='integration.gpx') # temporary use, not written
        for img in gpx.add_image(tifname,fmthint="TIF",cacheImage=True):  # loop unneeded for TIF (1 image/file)
            img.setControl('pixelSize',[150,150])
            if not cache: # load & compute controls & 2theta values once
                print('Initializing cache for',tifname)
                img.loadControls(imgctrl)	# set controls/calibrations/masks
                img.loadMasks(imgmask)
                cache['Image Controls'] = img.getControls() # save file contents for quick reload
                cache['Masks'] = img.getMasks()
                cache['intMaskMap'] = img.IntMaskMap() # calc mask & TA arrays to save for integrations
                cache['intTAmap'] = img.IntThetaAzMap()
                cache['FrameMask'] = img.MaskFrameMask() # calc Frame mask & T array to save for Pixel masking
                cache['maskTmap'] = img.MaskThetaMap() 
            else:
                img.setControls(cache['Image Controls'])
                img.setMasks(cache['Masks'],True)  # not using threshold masks
            img.GeneratePixelMask(esdMul=3,ThetaMap=cache['maskTmap'],FrameMask=cache['FrameMask'])
            for pwdr in img.Integrate(MaskMap=cache['intMaskMap'],ThetaAzimMap=cache['intTAmap']):
                pwdr.Export(os.path.join(pathxye,os.path.split(tifname)[1]),'.xye')  # '.tif in name ignored
            img.clearImageCache()  # save some space
            img.clearPixelMask()

        print(f'*=== image processed, time={time.time()-starttime:.3f} sec\n')
        del gpx

    if __name__ == '__main__':
        if nodes > 0: import multiprocessing as mp

        # make folder to store integrated images & integrated patterns if needed
        if not os.path.exists(pathImg): os.mkdir(pathImg)
        if not os.path.exists(pathxye):	os.mkdir(pathxye)

        if nodes > 0: pool = mp.Pool(nodes)

        while True:	 # Loop will never end, stop with ctrl+C
            tiflist = sorted(glob.glob(globPattern),key=lambda x: os.path.getctime(x)) # get images sorted by creation time, oldest 1st
            if not tiflist:
                time.sleep(0.1)
                continue
            intlist = []  # list of images read to process
            for tifname in tiflist:
                if not os.path.exists(tifname + '.metadata'): continue
                shutil.move(tifname, pathImg)	# move file before integration so that it is not found in another search
                shutil.move(tifname + '.metadata', pathImg)
                intlist.append(os.path.join(pathImg,os.path.split(tifname)[1]))
            if nodes == 0:
                for newtifname in intlist: integrate_tif(newtifname)
            else:
                pool.map(integrate_tif,intlist)

        if nodes > 0: pool.close()
        print(f'Total elapsed time={time.time()-scriptstart:.3f} sec')

.. _HistExport:

Histogram Export
----------------

This example shows how to export a series of histograms from a collection of 
.gpx (project) files. The Python ``glob()`` function is used to find all files 
matching a wildcard in the specified directory (``dataloc``). For each file 
there is a loop over histograms in that project and for each histogram 
:meth:`~GSASIIscriptable.G2PwdrData.Export` is called to write out the contents of that histogram
as CSV (comma-separated variable) file that contains data positions, 
observed, computed and background intensities as well as weighting for each 
point and Q. Note that for the Export call, there is more than one choice of
exporter that can write ``.csv`` extension files, so the export hint must 
be specified. 

.. code-block::  python

    import os,sys,glob
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')  # change this
    from . import GSASIIscriptable as G2sc

    dataloc = "/Users/toby/Scratch/"                 # where to find data 
    PathWrap = lambda fil: os.path.join(dataloc,fil) # EZ way 2 add dir to filename

    for f in glob.glob(PathWrap('bkg*.gpx')):  # put filename prefix here
        print(f)
        gpx = G2sc.G2Project(f)
        for i,h in enumerate(gpx.histograms()):
            hfil = os.path.splitext(f)[0]+'_'+str(i) # file to write
            print('\t',h.name,hfil+'.csv')
            h.Export(hfil,'.csv','histogram CSV')

.. _AutoBackground:

Automatic Background
---------------------

This example shows how to use the automatic background feature in GSAS-II to
compute an approximate background and set fixed background points from that
background. This approximately example follows that of the
`Autobackground Tutorial <https://advancedphotonsource.github.io/GSAS-II-tutorials/AutoBkg/AutoBkg.html>`_. In this example, a new project is created and
the data files from the tutorial are read. Note that scripting is not able
to read files from inside a zip archive or use defaulted instrument parameters.
The histograms are then processed in turn.
The first step is to use `calc_autobkg` to compute the fixed background points.
The refinement flag is then set for the Chebyschev polynomial terms and three
background peaks are added with the width flag set for refinement. The first
call to `fit_fixed_points()` will refine the three Chebyschev terms and
the intensities of the three background peaks to fit the
fixed background points. The refinement flags for 
the widths of the three background peaks are then set as well and the
refinement is repeated. The location of the third background peaks is added
and the refinement is repeated.
Finally, the number of Chebyschev polynomial terms is increased to six
and the refinement is repeated.
            
.. code-block::  python

    import os,glob
    import G2script as G2sc
    PathWrap = lambda fil: os.path.join('/tmp',fil)
    gpx = G2sc.G2Project(newgpx=PathWrap('autobkg.gpx'))
    for i in glob.glob(PathWrap('test_RampDown-*.xye')):
        hist = gpx.add_powder_histogram(i,PathWrap('testData.instprm'))
    for hist in gpx.histograms('PWDR'):
        hist.calc_autobkg(logLam=3.5)
        hist.set_refinements({"Background": {"no. coeffs": 3, "refine": True}})
        for pk in [2.4,3.1,4.75]:
            hist.add_back_peak(pk,1000,1000,0,[False,True,False,False])
        hist.fit_fixed_points()
        for i in [0,1,2]: hist.ref_back_peak(i,[False,True,True,False])
        hist.fit_fixed_points()
        hist.ref_back_peak(2,[True,True,True,False])
        hist.fit_fixed_points()
        hist.set_refinements({"Background": {"no. coeffs": 6, "refine": True}})
        hist.fit_fixed_points()
        gpx.save()

            
.. _CommandlineInterface:

GSASIIscriptable Command-line Interface
=======================================

The routines described above are intended to be called from a Python script, but an
alternate way to access some of the same functionality is to 
invoke the ``GSASIIscriptable.py`` script from 
the command line usually from within a shell script or batch file. 
This mode of accessing GSAS-II scripting does not appear to get much use and 
is no longer being developed. Please do communicate to the developers if 
keeping this mode of access would be of value in your work.

To use the command-line mode is done with a command like this::

       python <path/>GSASIIscriptable.py <subcommand> <file.gpx> <options>

The following subcommands are defined:

        * create, see :func:`~GSASIIscriptable.create`
        * add, see :func:`~GSASIIscriptable.add`
        * dump, see :func:`~GSASIIscriptable.dump`
        * refine, see :func:`~GSASIIscriptable.refine`
        * export, :func:`~GSASIIscriptable.export`
        * browse, see :func:`~GSASIIscriptable.IPyBrowse`

Run::

   python GSASIIscriptable.py --help

to show the available subcommands, and inspect each subcommand with
`python GSASIIscriptable.py <subcommand> --help` or see the documentation for each of the above routines.

.. _JsonFormat:

Parameters in JSON files
------------------------

The refine command requires two inputs: an existing GSAS-II project (.gpx) file and
a JSON format file
(see `Introducing JSON <http://json.org/>`_) that contains a single dict.
This dict may have two keys:

refinements:
  This defines the a set of refinement steps in a JSON representation of a
  :ref:`Refinement_recipe` list. 

code:
  This optionally defines Python code that will be executed after the project is loaded,
  but before the refinement is started. This can be used to execute Python code to change
  parameters that are not accessible via a :ref:`Refinement_recipe` dict (note that the
  project object is accessed with variable ``proj``) or to define code that will be called
  later (see key ``call`` in the :ref:`Refinement_recipe` section.)

JSON website: `Introducing JSON <http://json.org/>`_.

.. _API:

API: Complete Documentation
===========================

.. automodule:: GSASIIscriptable
    :members: 

