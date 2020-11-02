#!/usr/bin/env python
# -*- coding: utf-8 -*-
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
#
"""
*GSASIIscriptable: Scripting Interface*
=======================================

Routines to use an increasing amount of GSAS-II's capabilities from scripts, 
without use of the graphical user interface (GUI). GSASIIscriptable can create and access
GSAS-II project (.gpx) files and can directly perform image handling and refinements.  
The module defines wrapper classes (inheriting from :class:`G2ObjectWrapper`) for a growing number 
of data tree items.

GSASIIscriptable can be used in two ways. It offers a command-line mode 
(see :ref:`CommandlineInterface`) that 
provides access a number of features without writing Python scripts 
via shell/batch commands. The more widely used and more powerful mode of GSASIIscriptable is 
use is through Python scripts that 
call the module's application interface (API), see API summary that follows or the :ref:`API` 
section.

==================================================
Application Interface (API) Summary
==================================================
This section of the documentation provides an overview to API, with full documentation 
in the :ref:`API` section. The typical API use will be with a Python script, such as 
what is found in :ref:`CodeExamples`. Most functionality is provided via the objects and methods
summarized below.

---------------------
Overview of Classes 
---------------------

===============================    ===============================================================================================================
class                              Encapsulates 
===============================    ===============================================================================================================
:class:`G2Project`                  a GSAS-II project file, provides references to objects below,
                                    each corresponding to a tree item 
                                    (excepting :class:`G2AtomRecord`)
:class:`G2Phase`                    a phase
:class:`G2PwdrData`                 a powder histogram
:class:`G2Image`                    an image
:class:`G2PDF`                      a PDF histogram
:class:`G2SeqRefRes`                the sequential results table
:class:`G2AtomRecord`               an atom within a phase
===============================    ===============================================================================================================

---------------------
Functions
---------------------

A small amount of the Scriptable code does not require use of objects. 

==================================================    ===============================================================================================================
method                                                Use
==================================================    ===============================================================================================================
:func:`GenerateReflections`                            Generates a list of unique powder reflections 
:func:`SetPrintLevel`                                  Sets the amout of output generated when running a script
==================================================    ===============================================================================================================

---------------------------------
Class :class:`G2Project`
---------------------------------

  All GSASIIscriptable scripts will need to create a :class:`G2Project` object 
  either for a new GSAS-II project or to read in an existing project (.gpx) file. 
  The most commonly used routines in this object are:

.. tabularcolumns:: |l|p{3.5in}|

==================================================    ===============================================================================================================
method                                                Use
==================================================    ===============================================================================================================
:meth:`G2Project.save`                                Writes the current project to disk.

:meth:`G2Project.add_powder_histogram`                Used to read in powder diffraction data into a project file.

:meth:`G2Project.add_simulated_powder_histogram`      Defines a "dummy" powder diffraction data that will be simulated after a refinement step.

:meth:`G2Project.add_image`                           Reads in an image into a project.

:meth:`G2Project.add_phase`                           Adds a phase to a project

:meth:`G2Project.add_PDF`                             Adds a PDF entry to a project (does not compute it)

:meth:`G2Project.histograms`                          Provides a list of histograms in the current project, as :class:`G2PwdrData` objects

:meth:`G2Project.phases`                              Provides a list of phases defined in the current project, as :class:`G2Phase` objects

:meth:`G2Project.images`                              Provides a list of images in the current project, as :class:`G2Image` objects

:meth:`G2Project.pdfs`                                Provides a list of PDFs in the current project, as :class:`G2PDF` objects

:meth:`G2Project.seqref`                              Returns a :class:`G2SeqRefRes` object if there are Sequential Refinement results

:meth:`G2Project.do_refinements`                      This is passed a list of dictionaries, where each dict defines a refinement step. 
                                                      Passing a list with a single empty dict initiates a refinement with the current
                                                      parameters and flags. A refinement dict sets up a single refinement step 
                                                      (as described in :ref:`Project_dicts`). Also see :ref:`Refinement_recipe`.

:meth:`G2Project.set_refinement`                      This is passed a single dict which is used to set parameters and flags.
                                                      These actions can be performed also in :meth:`G2Project.do_refinements`. 
:meth:`G2Project.get_Variable`                        Retrieves the value and esd for a parameter
:meth:`G2Project.get_Covariance`                      Retrieves values and covariance for a set of refined parameters
:meth:`G2Project.set_Controls`                        Set overall GSAS-II control settings such as number of cycles and to set up a sequential
                                                      fit. (Also see :meth:`G2Project.get_Controls` to read values.)
:meth:`G2Project.imageMultiDistCalib`                 Performs a global calibration fit with images at multiple distance settings.
==================================================    ===============================================================================================================

---------------------------------
Class :class:`G2Phase`
---------------------------------


  Another common object in GSASIIscriptable scripts is :class:`G2Phase`, used to encapsulate each phase in a project, with commonly used methods:

.. tabularcolumns:: |l|p{3.5in}|

==================================================    ===============================================================================================================
method                                                Use
==================================================    ===============================================================================================================
:meth:`G2Phase.set_refinements`                       Provides a mechanism to set values and refinement flags for the phase. See :ref:`Phase_parameters_table` 
                                                      for more details. This information also can be supplied within a call to :meth:`G2Project.do_refinements` 
                                                      or :meth:`G2Project.set_refinement`.
:meth:`G2Phase.clear_refinements`                     Unsets refinement flags for the phase. 
:meth:`G2Phase.set_HAP_refinements`                   Provides a mechanism to set values and refinement flags for parameters specific to both this phase and 
                                                      one of its histograms. See :ref:`HAP_parameters_table`. This information also can be supplied within 
                                                      a call to :meth:`G2Project.do_refinements` or :meth:`G2Project.set_refinement`.
:meth:`G2Phase.clear_HAP_refinements`                 Clears refinement flags specific to both this phase and one of its histograms.
:meth:`G2Phase.getHAPvalues`                          Returns values of parameters specific to both this phase and one of its histograms.
:meth:`G2Phase.copyHAPvalues`                         Copies HAP settings between from one phase/histogram and to other histograms in same phase.
:meth:`G2Phase.atoms`                                 Returns a list of atoms in the phase
:meth:`G2Phase.atom`                                  Returns an atom from its label 
:meth:`G2Phase.histograms`                            Returns a list of histograms linked to the phase
:meth:`G2Phase.get_cell`                              Returns unit cell parameters (also see :meth:`G2Phase.get_cell_and_esd`)
:meth:`G2Phase.export_CIF`                            Writes a CIF for the phase
:meth:`G2Phase.setSampleProfile`                      Sets sample broadening parameters
==================================================    ===============================================================================================================

---------------------------------
Class :class:`G2PwdrData`
---------------------------------

  Another common object in GSASIIscriptable scripts is :class:`G2PwdrData`, which encapsulate each powder diffraction histogram in a project, with commonly used methods:

.. tabularcolumns:: |l|p{3.5in}|

==================================================    ===============================================================================================================
method                                                Use
==================================================    ===============================================================================================================
:meth:`G2PwdrData.set_refinements`                    Provides a mechanism to set values and refinement flags for the powder histogram. See 
                                                      :ref:`Histogram_parameters_table` for details.  
:meth:`G2PwdrData.clear_refinements`                  Unsets refinement flags for the the powder histogram.
:meth:`G2PwdrData.residuals`                          Reports R-factors etc. for the the powder histogram (also see :meth:`G2PwdrData.get_wR`) 
:meth:`G2PwdrData.add_back_peak`                      Adds a background peak to the histogram. Also see :meth:`G2PwdrData.del_back_peak` and 
                                                      :meth:`G2PwdrData.ref_back_peak`.
:meth:`G2PwdrData.fit_fixed_points`                   Fits background to the specified fixed points.
:meth:`G2PwdrData.getdata`                            Provides access to the diffraction data associated with the histogram.
:meth:`G2PwdrData.reflections`                        Provides access to the reflection lists for the histogram.
:meth:`G2PwdrData.Export`                             Writes the diffraction data or reflection list into a file
:meth:`G2PwdrData.add_peak`                           Adds a peak to the peak list. Also see :ref:`PeakRefine`.
:meth:`G2PwdrData.set_peakFlags`                      Sets refinement flags for peaks
:meth:`G2PwdrData.refine_peaks`                       Starts a peak/background fitting cycle
:attr:`G2PwdrData.Peaks`                              Provides access to the peak list data structure
:attr:`G2PwdrData.PeakList`                           Provides the peak list parameter values 
:meth:`G2PwdrData.Export_peaks`                       Writes the peak parameters to a text file 
==================================================    ===============================================================================================================

---------------------------------
Class :class:`G2Image`
---------------------------------

  When working with images, there will be a :class:`G2Image` object for each image (also see :meth:`G2Project.add_image`  and :meth:`G2Project.images`).

.. tabularcolumns:: |l|p{3.5in}|

==================================================    ===============================================================================================================
method                                                Use
==================================================    ===============================================================================================================
:meth:`G2Image.Recalibrate`                           Invokes a recalibration fit starting from the current Image Controls calibration coefficients.
:meth:`G2Image.Integrate`                             Invokes an image integration All parameters Image Controls will have previously been set.
:meth:`G2Image.setControl`                            Set an Image Controls parameter in the current image.
:meth:`G2Image.getControl`                            Return an Image Controls parameter in the current image.
:meth:`G2Image.findControl`                           Get the names of Image Controls parameters.
:meth:`G2Image.loadControls`                          Load controls from a .imctrl file (also see :meth:`G2Image.saveControls`).
:meth:`G2Image.loadMasks`                             Load masks from a .immask file.
:meth:`G2Image.setVary`                               Set a refinement flag for Image Controls parameter in the current image. (Also see :meth:`G2Image.getVary`)
:meth:`G2Image.setCalibrant`                          Set a calibrant type (or show choices) for the current image.
:meth:`G2Image.setControlFile`                        Set a image to be used as a background/dark/gain map image.
==================================================    ===============================================================================================================


---------------------------------
Class :class:`G2PDF`
---------------------------------

  To work with PDF entries, object :class:`G2PDF`, encapsulates a PDF entry with methods:

.. tabularcolumns:: |l|p{3.5in}|

==================================================    ===============================================================================================================
method                                                Use
==================================================    ===============================================================================================================
:meth:`G2PDF.export`                                   Used to write G(r), etc. as a file
:meth:`G2PDF.calculate`                                Computes the PDF using parameters in the object
:meth:`G2PDF.optimize`                                 Optimizes selected PDF parameters
:meth:`G2PDF.set_background`                           Sets the histograms used for sample background, container, etc. 
:meth:`G2PDF.set_formula`                              Sets the chemical formula for the sample
==================================================    ===============================================================================================================

---------------------------------
Class :class:`G2SeqRefRes`
---------------------------------

  To work with Sequential Refinement results, object :class:`G2SeqRefRes`, encapsulates the sequential refinement table with methods:

.. tabularcolumns:: |l|p{3.5in}|

==================================================    ===============================================================================================================
method                                                Use
==================================================    ===============================================================================================================
:meth:`G2SeqRefRes.histograms`                         Provides a list of histograms used in the Sequential Refinement 
:meth:`G2SeqRefRes.get_cell_and_esd`                   Returns cell dimensions and standard uncertainies for a phase and histogram from the Sequential Refinement 
:meth:`G2SeqRefRes.get_Variable`                       Retrieves the value and esd for a parameter from a particular histogram in the Sequential Refinement 
:meth:`G2SeqRefRes.get_Covariance`                     Retrieves values and covariance for a set of refined parameters for a particular histogram 
==================================================    ===============================================================================================================

---------------------------------
Class :class:`G2AtomRecord`
---------------------------------

  When working with phases, :class:`G2AtomRecord` objects provide access to the contents of each atom in a phase. This provides access to "properties" that can be 
  used to get values of much of the atoms associated settings: label, type, refinement_flags, coordinates, occupancy, ranId, adp_flag, and uiso. In addition, 
  refinement_flags, occupancy and uiso can be used to set values. See the :class:`G2AtomRecord` docs and source code.

.. _Refinement_dicts:

=====================
Refinement parameters
=====================
While scripts can be written that setup refinements by changing individual parameters 
through calls to the methods associated with objects that wrap each data tree item, 
many of these actions can be combined into fairly complex dict structures to conduct refinement
steps. Use of these dicts is required with the :ref:`CommandlineInterface`. This section of the 
documentation describes these dicts. 

.. _Project_dicts:

-----------------------------
Project-level Parameter Dict
-----------------------------

As noted below (:ref:`Refinement_parameters_kinds`), there are three types of refinement parameters,
which can be accessed individually by the objects that encapsulate individual phases and histograms
but it will often be simplest to create a composite dictionary
that is used at the project-level. A dict is created with keys
"set" and "clear" that can be supplied to :meth:`G2Project.set_refinement`
(or :meth:`G2Project.do_refinements`, see :ref:`Refinement_recipe` below) that will
determine parameter values and will determine which parameters will be refined. 

The specific keys and subkeys that can be used are defined in tables 
:ref:`Histogram_parameters_table`, :ref:`Phase_parameters_table` and :ref:`HAP_parameters_table`.

Note that optionally a list of histograms and/or phases can be supplied in the call to 
:meth:`G2Project.set_refinement`, but if not specified, the default is to use all defined
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

------------------------
Refinement recipe
------------------------

Building on the :ref:`Project_dicts`,
it is possible to specify a sequence of refinement actions as a list of
these dicts and supplying this list 
as an argument to :meth:`G2Project.do_refinements`.

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
may be used in a dict supplied to :meth:`G2Project.do_refinements`. Note that keys ``histograms``
and ``phases`` are used to limit actions to specific sets of parameters within the project. 

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
                       namespace inside :meth:`G2Project.iter_refinements` where
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

----------------------------
Refinement parameter types
----------------------------

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

=================================
Specifying Refinement Parameters
=================================

Refinement parameter values and flags to turn refinement on and off are specified within dictionaries,
where the details of these dicts are organized depends on the
type of parameter (see :ref:`Refinement_parameters_kinds`), with a different set
of keys (as described below) for each of the three types of parameters.

.. _Histogram_parameters_table:

--------------------
Histogram parameters
--------------------

This table describes the dictionaries supplied to :func:`G2PwdrData.set_refinements`
and :func:`G2PwdrData.clear_refinements`. As an example, 

.. code-block::  python

   hist.set_refinements({"Background": {"no.coeffs": 3, "refine": True},
                         "Sample Parameters": ["Scale"],
                         "Limits": [10000, 40000]})

With :meth:`G2Project.do_refinements`, these parameters should be placed inside a dict with a key
``set``, ``clear``, or ``once``. Values will be set for all histograms, unless the ``histograms``
key is used to define specific histograms. As an example: 

.. code-block::  python

  gsas_proj.do_refinements([
      {'set': {
          'Background': {'no.coeffs': 3, 'refine': True},
          'Sample Parameters': ['Scale'],
          'Limits': [10000, 40000]},
      'histograms': [1,2]}
                            ])

Note that below in the Instrument Parameters section, 
related profile parameters (such as U and V) are grouped together but
separated by commas to save space in the table.

.. tabularcolumns:: |l|l|p{3.5in}|

===================== ====================  =================================================
key                   subkey                explanation
===================== ====================  =================================================
Limits                                      The range of 2-theta (degrees) or TOF (in 
                                            microsec) range of values to use. Can
                                            be either a dictionary of 'low' and/or 'high',
                                            or a list of 2 items [low, high]
\\                     low                   Sets the low limit
\\                     high                  Sets the high limit

Sample Parameters                           Should be provided as a **list** of subkeys
                                            to set or clear,\\e.g. ['DisplaceX', 'Scale']
\\                     Absorption
\\                     Contrast
\\                     DisplaceX             Sample displacement along the X direction
\\                     DisplaceY             Sample displacement along the Y direction
\\                     Scale                 Histogram Scale factor

Background                                  Sample background. Value will be a dict or 
                                            a boolean. If True or False, the refine 
                                            parameter for background is set to that.
                                            Note that background peaks are not handled
                                            via this; see 
                                            :meth:`G2PwdrData.ref_back_peak` instead.
                                            When value is a dict,
                                            supply any of the following keys:
\\                     type                  The background model, e.g. 'chebyschev-1'
\\                     refine                The value of the refine flag, boolean
\\                     'no. coeffs'          Number of coefficients to use, integer
\\                     coeffs                List of floats, literal values for background
\\                     FixedPoints           List of (2-theta, intensity) values for fixed points
\\                     'fit fixed points'    If True, triggers a fit to the fixed points to
                                            be calculated. It is calculated when this key is
                                            detected, regardless of calls to refine.
\\                     peaks                 Specifies a set of flags for refining 
                                            background peaks as a nested list. There may
                                            be an item for each defined background peak
                                            (or fewer) and each item is a list with the flag 
                                            values for pos,int,sig & gam (fewer than 4 values 
                                            are allowed). 

Instrument Parameters                       As in Sample Paramters, provide as a **list** of
                                            subkeys to
                                            set or clear, e.g. ['X', 'Y', 'Zero', 'SH/L']
\\                     U, V, W               Gaussian peak profile terms
\\                     X, Y, Z               Lorentzian peak profile terms
\\                     alpha, beta-0,        TOF profile terms 
                      beta-1, beta-q,
\\                     sig-0, sig-1,         TOF profile terms
                      sig-2, sig-q
\\                     difA, difB, difC      TOF Calibration constants
\\                     Zero                  Zero shift
\\                     SH/L                  Finger-Cox-Jephcoat low-angle peak asymmetry
\\                     Polariz.              Polarization parameter
\\                     Lam                   Lambda, the incident wavelength
===================== ====================  =================================================

.. _Phase_parameters_table:

----------------
Phase parameters
----------------

This table describes the dictionaries supplied to :func:`G2Phase.set_refinements`
and :func:`G2Phase.clear_refinements`. With :meth:`G2Project.do_refinements`,
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
                      for fractional occupancy, 'X' for position,
                      and 'U' for Debye-Waller factor
LeBail                Enables LeBail intensity extraction.
======= ==========================================================


.. _HAP_parameters_table:


Histogram-and-phase parameters
------------------------------

This table describes the dictionaries supplied to :func:`G2Phase.set_HAP_refinements`
and :func:`G2Phase.clear_HAP_refinements`. When supplied to
:meth:`G2Project.do_refinements`, these parameters should be placed inside a dict with a key
``set``, ``clear``, or ``once``. Values will be set for all histograms used in each phase,
unless the ``histograms`` and ``phases`` keys are used to define specific phases and histograms.

.. tabularcolumns:: |l|l|p{3.5in}|

=============  ==========  ============================================================
key             subkey                 explanation
=============  ==========  ============================================================
Babinet                                Should be a **list** of the following
                                       subkeys. If not, assumes both
                                       BabA and BabU
\\              BabA
\\              BabU
Extinction                             Boolean, True to refine.
HStrain                                Boolean or list/tuple, True to refine all 
                                       appropriate D\\ :sub:`ij` terms or False
                                       to not refine any. If a list/tuple, will
                                       be a set of True & False values for each 
                                       D\\ :sub:`ij` term; number of items must 
                                       match number of terms.
Mustrain
\\              type                   Mustrain model. One of 'isotropic',
                                       'uniaxial', or 'generalized'. **Should always
                                       be included when Mustrain is used.**
\\              direction              For uniaxial only. A list of three
                                       integers,
                                       the [hkl] direction of the axis.
\\              refine                 Usually boolean, set to True to refine.
                                       or False to clear. 
                                       For uniaxial model, can specify a value
                                       of 'axial' or 'equatorial' to set that flag
                                       to True or a single
                                       boolean sets both axial and equatorial.
Size                                   
\\              type                   Size broadening model. One of 'isotropic',
                                       'uniaxial', or 'ellipsoid'. **Should always
                                       be specified when Size is used.**
\\              direction              For uniaxial only. A list of three
                                       integers,
                                       the [hkl] direction of the axis.
\\              refine                 Boolean, True to refine.
\\              value                  float, size value in microns 
Pref.Ori.                              Boolean, True to refine
Show                                   Boolean, True to refine
Use                                    Boolean, True to refine
Scale                                  Phase fraction; Boolean, True to refine
=============  ==========  ============================================================

------------------------
Histogram/Phase objects
------------------------
Each phase and powder histogram in a :class:`G2Project` object has an associated
object. Parameters within each individual object can be turned on and off by calling
:meth:`G2PwdrData.set_refinements` or :meth:`G2PwdrData.clear_refinements`
for histogram parameters;
:meth:`G2Phase.set_refinements` or :meth:`G2Phase.clear_refinements`
for phase parameters; and :meth:`G2Phase.set_HAP_refinements` or
:meth:`G2Phase.clear_HAP_refinements`. As an example, if some_histogram is a histogram object (of type :class:`G2PwdrData`), use this to set parameters in that histogram:

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

===================================
Access to other parameter settings
===================================

There are several hundred different types of values that can be stored in a 
GSAS-II project (.gpx) file. All can be changed from the GUI but only a 
subset have direct mechanism implemented for change from the GSASIIscriptable 
API. In practice all parameters in a .gpx file can be edited via scripting, 
but sometimes determining what should be set to implement a parameter 
change can be complex. 
Several routines, :meth:`G2Phase.getHAPentryList`, 
:meth:`G2Phase.getPhaseEntryList` and :meth:`G2PwdrData.getHistEntryList` 
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
    import GSASIIscriptable as G2sc
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

    Z:\\>fc before.txt after.txt
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

=================================
Code Examples
=================================

.. _PeakRefine:

--------------------
Peak Fitting
--------------------

Peak refinement is performed with routines 
:meth:`G2PwdrData.add_peak`, :meth:`G2PwdrData.set_peakFlags` and
:meth:`G2PwdrData.refine_peaks`. Method :meth:`G2PwdrData.Export_peaks` and
properties :attr:`G2PwdrData.Peaks` and :attr:`G2PwdrData.PeakList` 
provide ways to access the results. Note that when peak parameters are 
refined with :meth:`~G2PwdrData.refine_peaks`, the background may also
be refined. Use :meth:`G2PwdrData.set_refinements` to change background 
settings and the range of data used in the fit. See below for an example
peak refinement script, where the data files are taken from the 
"Rietveld refinement with CuKa lab Bragg-Brentano powder data" tutorial 
(in https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/LabData/data/).

.. code-block::  python

    from __future__ import division, print_function
    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII') # needed to "find" GSAS-II modules
    import GSASIIscriptable as G2sc
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
    hist.refine_peaks()
    print('peak positions: ',[i[0] for i in hist.PeakList])
    for i in range(len(hist.Peaks['peaks'])):
        print('peak',i,'pos=',hist.Peaks['peaks'][i][0],'sig=',hist.Peaks['sigDict']['pos'+str(i)])
    hist.Export_peaks('pkfit.txt')
    #gpx.save()  # gpx file is not written without this

-----------------------
Pattern Simulation
-----------------------

This shows two examples where a structure is read from a CIF, a 
pattern is computed using a instrument parameter file to specify the 
probe type (neutrons here) and wavelength. 

The first example uses a CW neutron instrument parameter file. 
The pattern is computed over a 2θ range of 5 to 120 degrees 
with 1000 points. 
The pattern and reflection list are written into files. 
Data files are found in the 
`Scripting Tutorial <https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/PythonScript/data/>`_.

.. code-block::  python

    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    import GSASIIscriptable as G2sc
    datadir = "/Users/toby/software/G2/Tutorials/PythonScript/data"
    PathWrap = lambda fil: os.path.join(datadir,fil)
    gpx = G2sc.G2Project(filename='PbSO4sim.gpx') # create a project    
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
`TOF-CW Joint Refinement Tutorial <https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/TOF-CW Joint Refinement/data>`_
tutorial. 

.. code-block::  python

    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    import GSASIIscriptable as G2sc
    cifdir = "/Users/toby/software/G2/Tutorials/PythonScript/data"
    datadir = "/Users/toby/software/G2/Tutorials/TOF-CW Joint Refinement/data"
    gpx = G2sc.G2Project(filename='/tmp/PbSO4simT.gpx') # create a project
    phase0 = gpx.add_phase(os.path.join(cifdir,"PbSO4-Wyckoff.cif"),
             phasename="PbSO4",fmthint='CIF') # add a phase to the project
    hist1 = gpx.add_simulated_powder_histogram("PbSO4 simulation",
                os.path.join(datadir,"POWGEN_1066.instprm"),14.,35.,
                phases=gpx.phases(),ibank=2)
    gpx.do_refinements([{}])
    gpx.save()

----------------------
Simple Refinement
----------------------

GSASIIscriptable can be used to setup and perform simple refinements. 
This example reads in an existing project (.gpx) file, adds a background
peak, changes some refinement flags and performs a refinement.

.. code-block::  python

    from __future__ import division, print_function
    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII') # needed to "find" GSAS-II modules
    import GSASIIscriptable as G2sc
    datadir = "/Users/Scratch/"
    gpx = G2sc.G2Project(os.path.join(datadir,'test2.gpx'))
    gpx.histogram(0).add_back_peak(4.5,30000,5000,0)
    pardict = {'set': {'Sample Parameters': ['Absorption', 'Contrast', 'DisplaceX'],
                       'Background': {'type': 'chebyschev-1', 'refine': True,
                                      'peaks':[[0,True]]}}}
    gpx.set_refinement(pardict)

----------------------
Sequential Refinement
----------------------

GSASIIscriptable can be used to setup and perform sequential refinements. This example script 
is used to take the single-dataset fit at the end of Step 1 of the 
`Sequential Refinement tutorial <https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/SeqRefine/SequentialTutorial.htm>`_
and turn on and off refinement flags, add histograms and setup the sequential fit, which is then run:

.. code-block::  python

    import os,sys,glob
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    import GSASIIscriptable as G2sc
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

----------------------
Image Processing
----------------------

A sample script where an image is read, assigned calibration values from a file 
and then integrated follows. 
The data files are found in the 
`Scripting Tutorial <https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/PythonScript/data/>`_.

.. code-block::  python

    import os,sys
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    import GSASIIscriptable as G2sc
    datadir = "/tmp"
    PathWrap = lambda fil: os.path.join(datadir,fil)

    gpx = G2sc.G2Project(filename=PathWrap('inttest.gpx'))
    imlst = gpx.add_image(PathWrap('Si_free_dc800_1-00000.tif'),fmthint="TIF")
    imlst[0].loadControls(PathWrap('Si_free_dc800_1-00000.imctrl'))
    pwdrList = imlst[0].Integrate()
    gpx.save()

This example shows a computation similar to what is done in tutorial 
`Area Detector Calibration with Multiple Distances <https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/DeterminingWavelength/DeterminingWavelength.html>`_

.. code-block::  python

    import os,sys,glob
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    import GSASIIscriptable as G2sc
    PathWrap = lambda fil: os.path.join(
        "/Users/toby/wp/Active/MultidistanceCalibration/multimg",
        fil)

    gpx = G2sc.G2Project(filename='/tmp/img.gpx')
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

----------------------
Image Calibration
----------------------

This example performs a number of cycles of constrained fitting. 
A project is created with the images found in a directory, setting initial
parameters as the images are read. The initial values 
for the calibration are not very good, so a :meth:`G2Image.Recalibrate` is done
to quickly improve the fit. Once that is done, a fit of all images is performed
where the wavelength, an offset and detector orientation are constrained to 
be the same for all images. The detector penetration correction is then added. 
Note that as the calibration values improve, the algorithm is able to find more 
points on diffraction rings to use for calibration and the number of "ring picks" 
increase. The calibration is repeated until that stops increasing significantly (<10%). 
Detector control files are then created. 
The files used for this exercise are found in the
`Area Detector Calibration Tutorial <https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/DeterminingWavelength/data/>`_
(see 
`Area Detector Calibration with Multiple Distances <https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/DeterminingWavelength/DeterminingWavelength.html>`_ ). 

.. code-block::  python

    import os,sys,glob
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    import GSASIIscriptable as G2sc
    PathWrap = lambda fil: os.path.join(
        "/Users/toby/wp/Active/MultidistanceCalibration/multimg",
        fil)

    gpx = G2sc.G2Project(filename='/tmp/calib.gpx')
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

.. _HistExport:

--------------------
Histogram Export
--------------------

This example shows how to export a series of histograms from a collection of 
.gpx (project) files. The Python ``glob()`` function is used to find all files 
matching a wildcard in the specified directory (``dataloc``). For each file 
there is a loop over histograms in that project and for each histogram 
:meth:`G2PwdrData.Export` is called to write out the contents of that histogram
as CSV (comma-separated variable) file that contains data positions, 
observed, computed and backgroun intensities as well as weighting for each 
point and Q. Note that for the Export call, there is more than one choice of
exporter that can write ``.csv`` extension files, so the export hint must 
be specified. 

.. code-block::  python

    import os,sys,glob
    sys.path.insert(0,'/Users/toby/software/G2/GSASII')  # change this
    import GSASIIscriptable as G2sc

    dataloc = "/Users/toby/Scratch/"                 # where to find data 
    PathWrap = lambda fil: os.path.join(dataloc,fil) # EZ way 2 add dir to filename

    for f in glob.glob(PathWrap('bkg*.gpx')):  # put filename prefix here
        print(f)
        gpx = G2sc.G2Project(f)
        for i,h in enumerate(gpx.histograms()):
            hfil = os.path.splitext(f)[0]+'_'+str(i) # file to write
            print('\t',h.name,hfil+'.csv')
            h.Export(hfil,'.csv','histogram CSV')

.. _CommandlineInterface:

=======================================
GSASIIscriptable Command-line Interface
=======================================

The routines described above are intended to be called from a Python script, but an
alternate way to access some of the same functionality is to 
invoke the ``GSASIIscriptable.py`` script from 
the command line usually from within a shell script or batch file. This
will usually be done with a command such as::

       python <path/>GSASIIscriptable.py <subcommand> <file.gpx> <options>

The following subcommands are defined:

        * create, see :func:`create`
        * add, see :func:`add`
        * dump, see :func:`dump`
        * refine, see :func:`refine`
        * export, :func:`export`
        * browse, see :func:`IPyBrowse`

Run::

   python GSASIIscriptable.py --help

to show the available subcommands, and inspect each subcommand with
`python GSASIIscriptable.py <subcommand> --help` or see the documentation for each of the above routines.

.. _JsonFormat:

-------------------------
Parameters in JSON files
-------------------------

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

============================================================
API: Complete Documentation
============================================================

The classes and modules in this module are described below.
A script will create one or more :class:`G2Project` objects by reading 
a GSAS-II project (.gpx) file or creating a new one and will then
perform actions such as adding a histogram (method :meth:`G2Project.add_powder_histogram`),
adding a phase (method :meth:`G2Project.add_phase`),
or setting parameters and performing a refinement
(method :meth:`G2Project.do_refinements`).

To change settings within histograms, images and phases, one usually needs to use
methods inside :class:`G2PwdrData`, :class:`G2Image` or :class:`G2Phase`. 
"""

#============================================================================
# Notes for adding a new object type
# 1) add a new object class (e.g. G2PDF)
# 2) add the wrapper into G2Project (e.g. _pdfs, pdf, pdfs)
# 3) add a new method to add the object into a project (G2Project.add_PDF)
# 4) add to documentation in section :class:`G2Project`
# 5) add a new documentation section for the new class
#============================================================================

from __future__ import division, print_function
import argparse
import os.path as ospath
#import datetime as dt
import sys
import platform
if '2' in platform.python_version_tuple()[0]:
    import cPickle
    strtypes = (str,unicode)
else:
    import pickle as cPickle
    strtypes = (str,bytes)
#import imp
import copy
import os
import random as ran

import numpy.ma as ma
import scipy.interpolate as si
import numpy as np
import scipy as sp

import GSASIIpath
GSASIIpath.SetBinaryPath(True)  # for now, this is needed before some of these modules can be imported
import GSASIIobj as G2obj
import GSASIIpwd as G2pwd
import GSASIIstrMain as G2strMain
import GSASIIstrIO as G2strIO
import GSASIIspc as G2spc
import GSASIIElem as G2elem
import GSASIIfiles as G2fil
import GSASIIimage as G2img

# Delay imports to not slow down small scripts that don't need them
Readers = {'Pwdr':[], 'Phase':[], 'Image':[]}
'''Readers by reader type'''
exportersByExtension = {}
'''Specifies the list of extensions that are supported for Powder data export'''
npsind = lambda x: np.sin(x*np.pi/180.)

def SetPrintLevel(level):
    '''Set the level of output from calls to :func:`GSASIIfiles.G2Print`, 
    which should be used in place of print() where possible. This is a 
    wrapper for :func:`GSASIIfiles.G2SetPrintLevel` so that this routine is
    documented here. 
    
    :param str level: a string used to set the print level, which may be 
      'all', 'warn', 'error' or 'none'.
      Note that capitalization and extra letters in level are ignored, so 
      'Warn', 'warnings', etc. will all set the mode to 'warn'
    '''
    G2fil.G2SetPrintLevel(level)
    global printLevel
    for mode in  'all', 'warn', 'error', 'none':
        if mode in level.lower():
            printLevel = mode
            return
        
def LoadG2fil():
    '''Setup GSAS-II importers. 
    Delay importing this module when possible, it is slow.
    Multiple calls are not. Only the first does anything.
    '''
    if len(Readers['Pwdr']) > 0: return
    # initialize imports
    Readers['Pwdr'] = G2fil.LoadImportRoutines("pwd", "Powder_Data")
    Readers['Phase'] = G2fil.LoadImportRoutines("phase", "Phase")
    Readers['Image'] = G2fil.LoadImportRoutines("img", "Image")

    # initialize exports
    for obj in G2fil.LoadExportRoutines(None):
        try:
            obj.Writer
        except AttributeError:
            continue
        for typ in obj.exporttype:
            if typ not in exportersByExtension:
                exportersByExtension[typ] = {obj.extension:obj}
            elif obj.extension in exportersByExtension[typ]:
                if type(exportersByExtension[typ][obj.extension]) is list:
                    exportersByExtension[typ][obj.extension].append(obj)
                else:
                    exportersByExtension[typ][obj.extension] = [
                        exportersByExtension[typ][obj.extension],
                        obj]
            else:
                exportersByExtension[typ][obj.extension] = obj

def LoadDictFromProjFile(ProjFile):
    '''Read a GSAS-II project file and load items to dictionary
    
    :param str ProjFile: GSAS-II project (name.gpx) full file name
    :returns: Project,nameList, where

      * Project (dict) is a representation of gpx file following the GSAS-II tree structure
        for each item: key = tree name (e.g. 'Controls','Restraints',etc.), data is dict
        data dict = {'data':item data whch may be list, dict or None,'subitems':subdata (if any)}
      * nameList (list) has names of main tree entries & subentries used to reconstruct project file

    Example for fap.gpx::

      Project = {                 #NB:dict order is not tree order
        'Phases':{'data':None,'fap':{phase dict}},
        'PWDR FAP.XRA Bank 1':{'data':[histogram data list],'Comments':comments,'Limits':limits, etc},
        'Rigid bodies':{'data': {rigid body dict}},
        'Covariance':{'data':{covariance data dict}},
        'Controls':{'data':{controls data dict}},
        'Notebook':{'data':[notebook list]},
        'Restraints':{'data':{restraint data dict}},
        'Constraints':{'data':{constraint data dict}}]
        }
      nameList = [                #NB: reproduces tree order
        ['Notebook',],
        ['Controls',],
        ['Covariance',],
        ['Constraints',],
        ['Restraints',],
        ['Rigid bodies',],
        ['PWDR FAP.XRA Bank 1',
             'Comments',
             'Limits',
             'Background',
             'Instrument Parameters',
             'Sample Parameters',
             'Peak List',
             'Index Peak List',
             'Unit Cells List',
             'Reflection Lists'],
        ['Phases', 'fap']
        ]
    '''
    # Let IOError be thrown if file does not exist
    if not ospath.exists(ProjFile):
        G2fil.G2Print ('\n*** Error attempt to open project file that does not exist: \n    {}'.
                   format(ProjFile))
        raise IOError('GPX file {} does not exist'.format(ProjFile))
    try:
        Project, nameList = G2strIO.GetFullGPX(ProjFile)
    except Exception as msg:
        raise IOError(msg)
    return Project,nameList

def SaveDictToProjFile(Project,nameList,ProjFile):
    '''Save a GSAS-II project file from dictionary/nameList created by LoadDictFromProjFile

    :param dict Project: representation of gpx file following the GSAS-II
        tree structure as described for LoadDictFromProjFile
    :param list nameList: names of main tree entries & subentries used to reconstruct project file
    :param str ProjFile: full file name for output project.gpx file (including extension)
    '''
    file = open(ProjFile,'wb')
    try:
        for name in nameList:
            data = []
            item = Project[name[0]]
            data.append([name[0],item['data']])
            for item2 in name[1:]:
                data.append([item2,item[item2]])
            cPickle.dump(data,file,1)
    finally:
        file.close()
    G2fil.G2Print('gpx file saved as %s'%ProjFile)

# def ImportPowder(reader,filename):
#     '''Use a reader to import a powder diffraction data file

#     :param str reader: a scriptable reader
#     :param str filename: full name of powder data file; can be "multi-Bank" data

#     :returns: list rdlist: list of reader objects containing powder data, one for each
#         "Bank" of data encountered in file. Items in reader object of interest are:

#           * rd.comments: list of str: comments found on powder file
#           * rd.dnames: list of str: data nammes suitable for use in GSASII data tree NB: duplicated in all rd entries in rdlist
#           * rd.powderdata: list of numpy arrays: pos,int,wt,zeros,zeros,zeros as needed for a PWDR entry in  GSASII data tree.
#     '''
#     rdfile,rdpath,descr = imp.find_module(reader)
#     rdclass = imp.load_module(reader,rdfile,rdpath,descr)
#     rd = rdclass.GSAS_ReaderClass()
#     if not rd.scriptable:
#         G2fil.G2Print(u'**** ERROR: '+reader+u' is not a scriptable reader')
#         return None
#     rdlist = []
#     if rd.ContentsValidator(filename):
#         repeat = True
#         rdbuffer = {} # create temporary storage for file reader
#         block = 0
#         while repeat: # loop if the reader asks for another pass on the file
#             block += 1
#             repeat = False
#             rd.objname = ospath.basename(filename)
#             flag = rd.Reader(filename,None,buffer=rdbuffer,blocknum=block,)
#             if flag:
#                 rdlist.append(copy.deepcopy(rd)) # save the result before it is written over
#                 if rd.repeat:
#                     repeat = True
#         return rdlist
#     G2fil.G2Print(rd.errors)
#     return None

def SetDefaultDData(dType,histoName,NShkl=0,NDij=0):
    '''Create an initial Histogram dictionary

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    '''
    if dType in ['SXC','SNC']:
        return {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
            'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
            'Extinction':['Lorentzian','None', {'Tbar':0.1,'Cos2TM':0.955,
            'Eg':[1.e-10,False],'Es':[1.e-10,False],'Ep':[1.e-10,False]}],
            'Flack':[0.0,False]}
    elif dType == 'SNT':
        return {'Histogram':histoName,'Show':False,'Scale':[1.0,True],
            'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]},
            'Extinction':['Lorentzian','None', {
            'Eg':[1.e-10,False],'Es':[1.e-10,False],'Ep':[1.e-10,False]}]}
    elif 'P' in dType:
        return {'Histogram':histoName,'Show':False,'Scale':[1.0,False],
            'Pref.Ori.':['MD',1.0,False,[0,0,1],0,{},[],0.1],
            'Size':['isotropic',[1.,1.,1.],[False,False,False],[0,0,1],
                [1.,1.,1.,0.,0.,0.],6*[False,]],
            'Mustrain':['isotropic',[1000.0,1000.0,1.0],[False,False,False],[0,0,1],
                NShkl*[0.01,],NShkl*[False,]],
            'HStrain':[NDij*[0.0,],NDij*[False,]],
            'Extinction':[0.0,False],'Babinet':{'BabA':[0.0,False],'BabU':[0.0,False]}}


def PreSetup(data):
    '''Create part of an initial (empty) phase dictionary

    from GSASIIphsGUI.py, near end of UpdatePhaseData

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    '''
    if 'RBModels' not in data:
        data['RBModels'] = {}
    if 'MCSA' not in data:
        data['MCSA'] = {'Models':[{'Type':'MD','Coef':[1.0,False,[.8,1.2],],'axis':[0,0,1]}],'Results':[],'AtInfo':{}}
    if 'dict' in str(type(data['MCSA']['Results'])):
        data['MCSA']['Results'] = []
    if 'Modulated' not in data['General']:
        data['General']['Modulated'] = False
#    if 'modulated' in data['General']['Type']:
#        data['General']['Modulated'] = True
#        data['General']['Type'] = 'nuclear'


def SetupGeneral(data, dirname):
    """Helps initialize phase data.

    From GSASIIphsGui.py, function of the same name. Minor changes for imports etc.

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    mapDefault = {'MapType':'','RefList':'','Resolution':0.5,'Show bonds':True,
                'rho':[],'rhoMax':0.,'mapSize':10.0,'cutOff':50.,'Flip':False}
    generalData = data['General']
    atomData = data['Atoms']
    generalData['AtomTypes'] = []
    generalData['Isotopes'] = {}

    if 'Isotope' not in generalData:
        generalData['Isotope'] = {}
    if 'Data plot type' not in generalData:
        generalData['Data plot type'] = 'Mustrain'
    if 'POhkl' not in generalData:
        generalData['POhkl'] = [0,0,1]
    if 'Map' not in generalData:
        generalData['Map'] = mapDefault.copy()
    if 'Flip' not in generalData:
        generalData['Flip'] = {'RefList':'','Resolution':0.5,'Norm element':'None',
            'k-factor':0.1,'k-Max':20.,}
    if 'testHKL' not in generalData['Flip']:
        generalData['Flip']['testHKL'] = [[0,0,2],[2,0,0],[1,1,1],[0,2,0],[1,2,3]]
    if 'doPawley' not in generalData:
        generalData['doPawley'] = False     #ToDo: change to ''
    if 'Pawley dmin' not in generalData:
        generalData['Pawley dmin'] = 1.0
    if 'Pawley neg wt' not in generalData:
        generalData['Pawley neg wt'] = 0.0
    if 'Algolrithm' in generalData.get('MCSA controls',{}) or \
        'MCSA controls' not in generalData:
        generalData['MCSA controls'] = {'Data source':'','Annealing':[50.,0.001,50],
        'dmin':2.0,'Algorithm':'log','Jump coeff':[0.95,0.5],'boltzmann':1.0,
        'fast parms':[1.0,1.0,1.0],'log slope':0.9,'Cycles':1,'Results':[],'newDmin':True}
    if 'AtomPtrs' not in generalData:
        generalData['AtomPtrs'] = [3,1,7,9]
        if generalData['Type'] == 'macromolecular':
            generalData['AtomPtrs'] = [6,4,10,12]
        elif generalData['Type'] == 'magnetic':
            generalData['AtomPtrs'] = [3,1,10,12]
    if generalData['Modulated']:
        generalData['Type'] = 'nuclear'
        if 'Super' not in generalData:
            generalData['Super'] = 1
            generalData['SuperVec'] = [[0,0,.1],False,4]
            generalData['SSGData'] = {}
        if '4DmapData' not in generalData:
            generalData['4DmapData'] = mapDefault.copy()
            generalData['4DmapData'].update({'MapType':'Fobs'})
    if 'Modulated' not in generalData:
        generalData['Modulated'] = False
    if 'HydIds' not in generalData:
        generalData['HydIds'] = {}
    cx,ct,cs,cia = generalData['AtomPtrs']
    generalData['NoAtoms'] = {}
    generalData['BondRadii'] = []
    generalData['AngleRadii'] = []
    generalData['vdWRadii'] = []
    generalData['AtomMass'] = []
    generalData['Color'] = []
    if generalData['Type'] == 'magnetic':
        generalData['MagDmin'] = generalData.get('MagDmin',1.0)
        landeg = generalData.get('Lande g',[])
    generalData['Mydir'] = dirname
    badList = {}
    for iat,atom in enumerate(atomData):
        atom[ct] = atom[ct].lower().capitalize()              #force to standard form
        if generalData['AtomTypes'].count(atom[ct]):
            generalData['NoAtoms'][atom[ct]] += atom[cx+3]*float(atom[cs+1])
        elif atom[ct] != 'UNK':
            Info = G2elem.GetAtomInfo(atom[ct])
            if not Info:
                if atom[ct] not in badList:
                    badList[atom[ct]] = 0
                badList[atom[ct]] += 1
                atom[ct] = 'UNK'
                continue
            atom[ct] = Info['Symbol'] # N.B. symbol might be changed by GetAtomInfo
            generalData['AtomTypes'].append(atom[ct])
            generalData['Z'] = Info['Z']
            generalData['Isotopes'][atom[ct]] = Info['Isotopes']
            generalData['BondRadii'].append(Info['Drad'])
            generalData['AngleRadii'].append(Info['Arad'])
            generalData['vdWRadii'].append(Info['Vdrad'])
            if atom[ct] in generalData['Isotope']:
                if generalData['Isotope'][atom[ct]] not in generalData['Isotopes'][atom[ct]]:
                    isotope = list(generalData['Isotopes'][atom[ct]].keys())[-1]
                    generalData['Isotope'][atom[ct]] = isotope
                generalData['AtomMass'].append(Info['Isotopes'][generalData['Isotope'][atom[ct]]]['Mass'])
            else:
                generalData['Isotope'][atom[ct]] = 'Nat. Abund.'
                if 'Nat. Abund.' not in generalData['Isotopes'][atom[ct]]:
                    isotope = list(generalData['Isotopes'][atom[ct]].keys())[-1]
                    generalData['Isotope'][atom[ct]] = isotope
                generalData['AtomMass'].append(Info['Mass'])
            generalData['NoAtoms'][atom[ct]] = atom[cx+3]*float(atom[cs+1])
            generalData['Color'].append(Info['Color'])
            if generalData['Type'] == 'magnetic':
                if len(landeg) < len(generalData['AtomTypes']):
                    landeg.append(2.0)
    if generalData['Type'] == 'magnetic':
        generalData['Lande g'] = landeg[:len(generalData['AtomTypes'])]

    if badList:
        msg = 'Warning: element symbol(s) not found:'
        for key in badList:
            msg += '\n\t' + key
            if badList[key] > 1:
                msg += ' (' + str(badList[key]) + ' times)'
        raise G2ScriptException("Phase error:\n" + msg)
        # wx.MessageBox(msg,caption='Element symbol error')
    F000X = 0.
    F000N = 0.
    for i,elem in enumerate(generalData['AtomTypes']):
        F000X += generalData['NoAtoms'][elem]*generalData['Z']
        isotope = generalData['Isotope'][elem]
        F000N += generalData['NoAtoms'][elem]*generalData['Isotopes'][elem][isotope]['SL'][0]
    generalData['F000X'] = F000X
    generalData['F000N'] = F000N
    import GSASIImath as G2mth
    generalData['Mass'] = G2mth.getMass(generalData)


def make_empty_project(author=None, filename=None):
    """Creates an dictionary in the style of GSASIIscriptable, for an empty
    project.

    If no author name or filename is supplied, 'no name' and
    <current dir>/test_output.gpx are used , respectively.

    Returns: project dictionary, name list

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    if not filename:
        filename = 'test_output.gpx'
    filename = os.path.abspath(filename)
    gsasii_version = str(GSASIIpath.GetVersionNumber())
    LoadG2fil()
    try:
        import matplotlib as mpl
        python_library_versions = G2fil.get_python_versions([mpl, np, sp])
    except ImportError:
        python_library_versions = G2fil.get_python_versions([np, sp])

    controls_data = dict(G2obj.DefaultControls)
    controls_data['LastSavedAs'] = filename
    controls_data['LastSavedUsing'] = gsasii_version
    controls_data['PythonVersions'] = python_library_versions
    if author:
        controls_data['Author'] = author

    output = {'Constraints': {'data': {'HAP': [], 'Hist': [], 'Phase': [],
                                       'Global': []}},
              'Controls': {'data': controls_data},
              u'Covariance': {'data': {}},
              u'Notebook': {'data': ['']},
              u'Restraints': {'data': {}},
              u'Rigid bodies': {'data': {'RBIds': {'Residue': [], 'Vector': []},
                                'Residue': {'AtInfo': {}},
                                'Vector':  {'AtInfo': {}}}}}

    names = [[u'Notebook'], [u'Controls'], [u'Covariance'],
             [u'Constraints'], [u'Restraints'], [u'Rigid bodies']]

    return output, names

def GenerateReflections(spcGrp,cell,Qmax=None,dmin=None,TTmax=None,wave=None):
    """Generates the crystallographically unique powder diffraction reflections
    for a lattice and space group (see :func:`GSASIIlattice.GenHLaue`). 

    :param str spcGrp: A GSAS-II formatted space group (with spaces between 
       axial fields, e.g. 'P 21 21 21' or 'P 42/m m c'). Note that non-standard
       space groups, such as 'P 21/n' or 'F -1' are allowed (see 
       :func:`GSASIIspc.SpcGroup`). 
    :param list cell: A list/tuple with six unit cell constants, 
      (a, b, c, alpha, beta, gamma) with values in Angstroms/degrees.
      Note that the cell constants are not checked for consistency 
      with the space group.
    :param float Qmax: Reflections up to this Q value are computed 
       (do not use with dmin or TTmax)
    :param float dmin: Reflections with d-space above this value are computed 
       (do not use with Qmax or TTmax)
    :param float TTmax: Reflections up to this 2-theta value are computed 
       (do not use with dmin or Qmax, use of wave is required.)
    :param float wave: wavelength in Angstroms for use with TTmax (ignored
       otherwise.)
    :returns: a list of reflections, where each reflection contains four items:
       h, k, l, d, where d is the d-space (Angstroms)

    Example:

    >>> import os,sys
    >>> sys.path.insert(0,'/Users/toby/software/G2/GSASII')
    >>> import GSASIIscriptable as G2sc
    GSAS-II binary directory: /Users/toby/software/G2/GSASII/bin
    17 values read from config file /Users/toby/software/G2/GSASII/config.py
    >>> refs = G2sc.GenerateReflections('P 1',
    ...                     (5.,6.,7.,90.,90.,90),
    ...                     TTmax=20,wave=1)
    >>> for r in refs: print(r)
    ... 
    [0, 0, 1, 7.0]
    [0, 1, 0, 6.0]
    [1, 0, 0, 5.0]
    [0, 1, 1, 4.55553961419178]
    [0, 1, -1, 4.55553961419178]
    [1, 0, 1, 4.068667356033675]
    [1, 0, -1, 4.068667356033674]
    [1, 1, 0, 3.8411063979868794]
    [1, -1, 0, 3.8411063979868794]
    """

    import GSASIIlattice as G2lat
    if len(cell) != 6:
        raise G2ScriptException("GenerateReflections: Invalid unit cell:" + str(cell))
    opts = (Qmax is not None) + (dmin is not None) + (TTmax is not None)
    if Qmax:
        dmin = 2 * np.pi / Qmax
        #print('Q,d',Qmax,dmin)
    elif TTmax and wave is None:
        raise G2ScriptException("GenerateReflections: specify a wavelength with TTmax")
    elif TTmax:
        dmin = wave / (2.0 * np.sin(np.pi*TTmax/360.))
        #print('2theta,d',TTmax,dmin)
    if opts != 1:
        raise G2ScriptException("GenerateReflections: specify one Qmax, dmin or TTmax")
    err,SGData = G2spc.SpcGroup(spcGrp)
    if err != 0:
        print('GenerateReflections space group error:',G2spc.SGErrors(err))
        raise G2ScriptException("GenerateReflections: Invalid space group: " + str(spcGrp))
    A = G2lat.cell2A(cell)
    return G2lat.GenHLaue(dmin,SGData,A)

class G2ImportException(Exception):
    pass

class G2ScriptException(Exception):
    pass

def import_generic(filename, readerlist, fmthint=None, bank=None):
    """Attempt to import a filename, using a list of reader objects.

    Returns the first reader object which worked."""
    # Translated from OnImportGeneric method in GSASII.py
    primaryReaders, secondaryReaders = [], []
    for reader in readerlist:
        if fmthint is not None and fmthint not in reader.formatName: continue
        flag = reader.ExtensionValidator(filename)
        if flag is None:
            secondaryReaders.append(reader)
        elif flag:
            primaryReaders.append(reader)
    if not secondaryReaders and not primaryReaders:
        raise G2ImportException("Could not read file: ", filename)

    with open(filename, 'r'):
        rd_list = []

        for rd in primaryReaders + secondaryReaders:
            # Initialize reader
            rd.selections = []
            if bank is None:
                rd.selections = []
            else:
                rd.selections = [bank-1]
            rd.dnames = []
            rd.ReInitialize()
            # Rewind file
            rd.errors = ""
            if not rd.ContentsValidator(filename):
                # Report error
                G2fil.G2Print("Warning: File {} has a validation error, continuing".format(filename))
            if len(rd.selections) > 1:
                raise G2ImportException("File {} has {} banks. Specify which bank to read with databank param."
                                .format(filename,len(rd.selections)))

            block = 0
            rdbuffer = {}
            repeat = True
            while repeat:
                repeat = False
                block += 1
                rd.objname = os.path.basename(filename)
                try:
                    flag = rd.Reader(filename,buffer=rdbuffer, blocknum=block)
                except:
                    flag = False
                if flag:
                    # Omitting image loading special cases
                    rd.readfilename = filename
                    rd_list.append(copy.deepcopy(rd))
                    repeat = rd.repeat
                else:
                    G2fil.G2Print("Warning: {} Reader failed to read {}".format(rd.formatName,filename))
            if rd_list:
                if rd.warnings:
                    G2fil.G2Print("Read warning by", rd.formatName, "reader:",
                          rd.warnings)
                elif bank is None:
                    G2fil.G2Print("{} read by Reader {}"
                              .format(filename,rd.formatName))
                else:
                    G2fil.G2Print("{} block # {} read by Reader {}"
                              .format(filename,bank,rd.formatName))
                return rd_list
    raise G2ImportException("No reader could read file: " + filename)


def load_iprms(instfile, reader, bank=None):
    """Loads instrument parameters from a file, and edits the
    given reader.

    Returns a 2-tuple of (Iparm1, Iparm2) parameters
    """
    LoadG2fil()
    ext = os.path.splitext(instfile)[1]

    if ext.lower() == '.instprm':
        # New GSAS File, load appropriate bank
        with open(instfile) as f:
            lines = f.readlines()
        if bank is None: 
            bank = reader.powderentry[2] 
        numbanks = reader.numbanks
        iparms = G2fil.ReadPowderInstprm(lines, bank, numbanks, reader)
        reader.instfile = instfile
        reader.instmsg = '{} (G2 fmt) bank {}'.format(instfile,bank)
        return iparms
    elif ext.lower() not in ('.prm', '.inst', '.ins'):
        raise ValueError('Expected .prm file, found: ', instfile)

    # It's an old GSAS file, load appropriately
    Iparm = {}
    with open(instfile, 'r') as fp:
        for line in fp:
            if '#' in line:
                continue
            Iparm[line[:12]] = line[12:-1]
    ibanks = int(Iparm.get('INS   BANK  ', '1').strip())
    if bank is not None:
        # pull out requested bank # bank from the data, and change the bank to 1
        Iparm,IparmC = {},Iparm
        for key in IparmC:
            if 'INS' not in key[:3]: continue   #skip around rubbish lines in some old iparm
            if key[4:6] == "  ":
                Iparm[key] = IparmC[key]
            elif int(key[4:6].strip()) == bank:
                Iparm[key[:4]+' 1'+key[6:]] = IparmC[key]            
        reader.instbank = bank
    elif ibanks == 1:
        reader.instbank = 1
    else: 
        raise G2ImportException("Instrument parameter file has {} banks, select one with instbank param."
                                    .format(ibanks))
    reader.powderentry[2] = 1
    reader.instfile = instfile
    reader.instmsg = '{} bank {}'.format(instfile,reader.instbank)
    return G2fil.SetPowderInstParms(Iparm, reader)

def load_pwd_from_reader(reader, instprm, existingnames=[],bank=None):
    """Loads powder data from a reader object, and assembles it into a GSASII data tree.

    :returns: (name, tree) - 2-tuple of the histogram name (str), and data

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    HistName = 'PWDR ' + G2obj.StripUnicode(reader.idstring, '_')
    HistName = G2obj.MakeUniqueLabel(HistName, existingnames)

    try:
        Iparm1, Iparm2 = instprm
    except ValueError:
        Iparm1, Iparm2 = load_iprms(instprm, reader, bank=bank)
        G2fil.G2Print('Instrument parameters read:',reader.instmsg)
    Ymin = np.min(reader.powderdata[1])
    Ymax = np.max(reader.powderdata[1])
    valuesdict = {'wtFactor': 1.0,
                  'Dummy': False,
                  'ranId': ran.randint(0, sys.maxsize),
                  'Offset': [0.0, 0.0], 'delOffset': 0.02*Ymax,
                  'refOffset': -0.1*Ymax, 'refDelt': 0.1*Ymax,
                  'Yminmax': [Ymin, Ymax]}
    reader.Sample['ranId'] = valuesdict['ranId']
    if 'T' in Iparm1['Type'][0]:
        if not reader.clockWd and reader.GSAS:
            reader.powderdata[0] *= 100.        #put back the CW centideg correction

    # Ending keys:
    # [u'Reflection Lists',
    #  u'Limits',
    #  'data',
    #  u'Index Peak List',
    #  u'Comments',
    #  u'Unit Cells List',
    #  u'Sample Parameters',
    #  u'Peak List',
    #  u'Background',
    #  u'Instrument Parameters']
    Tmin = np.min(reader.powderdata[0])
    Tmax = np.max(reader.powderdata[0])

    default_background = [['chebyschev-1', False, 3, 1.0, 0.0, 0.0],
                          {'nDebye': 0, 'debyeTerms': [], 'nPeaks': 0, 'peaksList': []}]

    output_dict = {u'Reflection Lists': {},
                   u'Limits': reader.pwdparms.get('Limits', [(Tmin, Tmax), [Tmin, Tmax]]),
                   u'data': [valuesdict, reader.powderdata, HistName],
                   u'Index Peak List': [[], []],
                   u'Comments': reader.comments,
                   u'Unit Cells List': [],
                   u'Sample Parameters': reader.Sample,
                   u'Peak List': {'peaks': [], 'sigDict': {}},
                   u'Background': reader.pwdparms.get('Background', default_background),
                   u'Instrument Parameters': [Iparm1, Iparm2],
                   }

    names = [u'Comments',
             u'Limits',
             u'Background',
             u'Instrument Parameters',
             u'Sample Parameters',
             u'Peak List',
             u'Index Peak List',
             u'Unit Cells List',
             u'Reflection Lists']

    # TODO controls?? GSASII.py:1664-7

    return HistName, [HistName] + names, output_dict

def _deep_copy_into(from_, into):
    """Helper function for reloading .gpx file. See G2Project.reload()

    :author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    if isinstance(from_, dict) and isinstance(into, dict):
        combined_keys = set(from_.keys()).union(into.keys())
        for key in combined_keys:
            if key in from_ and key in into:
                both_dicts = (isinstance(from_[key], dict)
                              and isinstance(into[key], dict))
                both_lists = (isinstance(from_[key], list)
                              and isinstance(into[key], list))
                if both_dicts or both_lists:
                    _deep_copy_into(from_[key], into[key])
                else:
                    into[key] = from_[key]
            elif key in from_:
                into[key] = from_[key]
            else:  # key in into
                del into[key]
    elif isinstance(from_, list) and isinstance(into, list):
        if len(from_) == len(into):
            for i in range(len(from_)):
                both_dicts = (isinstance(from_[i], dict)
                              and isinstance(into[i], dict))
                both_lists = (isinstance(from_[i], list)
                              and isinstance(into[i], list))
                if both_dicts or both_lists:
                    _deep_copy_into(from_[i], into[i])
                else:
                    into[i] = from_[i]
        else:
            into[:] = from_

def _getCorrImage(ImageReaderlist,proj,imageRef):
    '''Gets image & applies dark, background & flat background corrections.
    based on :func:`GSASIIimgGUI.GetImageZ`. Expected to be for internal
    use only.

    :param list ImageReaderlist: list of Reader objects for images
    :param object proj: references a :class:`G2Project` project
    :param imageRef: A reference to the desired image in the project. 
      Either the Image tree name (str), the image's index (int) or
      a image object (:class:`G2Image`)

    :return: array sumImg: corrected image for background/dark/flat back
    '''
    ImgObj = proj.image(imageRef)
    Controls = ImgObj.data['Image Controls']
    formatName = Controls.get('formatName','')
    imagefile = ImgObj.data['data'][1]
    if isinstance(imagefile, tuple) or isinstance(imagefile, list):
        imagefile, ImageTag =  imagefile # fix for multiimage files
    else:
        ImageTag = None # single-image file
    sumImg = G2fil.RereadImageData(ImageReaderlist,imagefile,ImageTag=ImageTag,FormatName=formatName)
    if sumImg is None:
        return []
    sumImg = np.array(sumImg,dtype='int32')
    darkImg = False
    if 'dark image' in Controls:
        darkImg,darkScale = Controls['dark image']
        if darkImg:
            dImgObj = proj.image(darkImg)
            formatName = dImgObj.data['Image Controls'].get('formatName','')
            imagefile = dImgObj.data['data'][1]
            if type(imagefile) is tuple:
                imagefile,ImageTag  = imagefile
            darkImage = G2fil.RereadImageData(ImageReaderlist,imagefile,ImageTag=ImageTag,FormatName=formatName)
            if darkImg is None:
                raise Exception('Error reading dark image {}'.format(imagefile))
            sumImg += np.array(darkImage*darkScale,dtype='int32')
    if 'background image' in Controls:
        backImg,backScale = Controls['background image']            
        if backImg:     #ignores any transmission effect in the background image
            bImgObj = proj.image(backImg)
            formatName = bImgObj.data['Image Controls'].get('formatName','')
            imagefile = bImgObj.data['data'][1]
            ImageTag = None # fix this for multiimage files
            backImage = G2fil.RereadImageData(ImageReaderlist,imagefile,ImageTag=ImageTag,FormatName=formatName)
            if backImage is None:
                raise Exception('Error reading background image {}'.format(imagefile))
            if darkImg:
                backImage += np.array(darkImage*darkScale/backScale,dtype='int32')
            else:
                sumImg += np.array(backImage*backScale,dtype='int32')
    if 'Gain map' in Controls:
        gainMap = Controls['Gain map']
        if gainMap:
            gImgObj = proj.image(gainMap)
            formatName = gImgObj.data['Image Controls'].get('formatName','')
            imagefile = gImgObj.data['data'][1]
            ImageTag = None # fix this for multiimage files
            GMimage = G2fil.RereadImageData(ImageReaderlist,imagefile,ImageTag=ImageTag,FormatName=formatName)
            if GMimage is None:
                raise Exception('Error reading Gain map image {}'.format(imagefile))
            sumImg = sumImg*GMimage/1000
    sumImg -= int(Controls.get('Flat Bkg',0))
    Imax = np.max(sumImg)
    Controls['range'] = [(0,Imax),[0,Imax]]
    return np.asarray(sumImg,dtype='int32')

def patchControls(Controls):
    '''patch routine to convert variable names used in parameter limits 
    to G2VarObj objects 
    (See :ref:`Parameter Limits<ParameterLimits>` description.)
    '''
    #patch (added Oct 2020) convert variable names for parm limits to G2VarObj
    for d in 'parmMaxDict','parmMinDict':
        if d not in Controls: Controls[d] = {}
        for k in Controls[d]:  
            if type(k) is str:
                print("Applying patch to Controls['{}']".format(d))
                Controls[d] = {G2obj.G2VarObj(k):v for k,v in Controls[d].items()}
                break
    conv = False
    if 'parmFrozen' not in Controls: Controls['parmFrozen'] = {}
    for k in Controls['parmFrozen']:
        for item in Controls['parmFrozen'][k]:
            if type(item) is str:
                conv = True
                Controls['parmFrozen'][k] = [G2obj.G2VarObj(i) for i in Controls['parmFrozen'][k]]
                break
    if conv: print("Applying patch to Controls['parmFrozen']")
    # end patch

class G2ObjectWrapper(object):
    """Base class for all GSAS-II object wrappers.

    The underlying GSAS-II format can be accessed as `wrapper.data`. A number
    of overrides are implemented so that the wrapper behaves like a dictionary.

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    def __init__(self, datadict):
        self.data = datadict

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __contains__(self, key):
        return key in self.data

    def get(self, k, d=None):
        return self.data.get(k, d)

    def keys(self):
        return self.data.keys()

    def values(self):
        return self.data.values()

    def items(self):
        return self.data.items()


class G2Project(G2ObjectWrapper):    
    """Represents an entire GSAS-II project. The object contains these 
    class variables:

     * G2Project.filename: contains the .gpx filename
     * G2Project.names: contains the contents of the project "tree" as a list 
       of lists. Each top-level entry in the tree is an item in the list. The 
       name of the top-level item is the first item in the inner list. Children
       of that item, if any, are subsequent entries in that list.
     * G2Project.data: contains the entire project as a dict. The keys 
       for the dict are the top-level names in the project tree (initial items
       in the G2Project.names inner lists) and each top-level item is stored
       as a dict. 

         * The contents of Top-level entries will be found in the item 
           named 'data', as an example, ``G2Project.data['Notebook']['data']`` 

         * The contents of child entries will be found in the item
           using the names of the parent and child, for example 
           ``G2Project.data['Phases']['NaCl']``

    :param str gpxfile: Existing .gpx file to be loaded. If nonexistent,
            creates an empty project.
    :param str author: Author's name (not yet implemented)
    :param str newgpx: The filename the project should be saved to in
            the future. If both newgpx and gpxfile are present, the project is
            loaded from the gpxfile, then when saved will be written to newgpx.
    :param str filename: Name to be used to save the project. Has same function as
            parameter newgpx (do not use both gpxfile and filename). Use of newgpx
            is preferred over filename.

    There are two ways to initialize this object:

    >>> # Load an existing project file
    >>> proj = G2Project('filename.gpx')
    
    >>> # Create a new project
    >>> proj = G2Project(newgpx='new_file.gpx')
    
    Histograms can be accessed easily.

    >>> # By name
    >>> hist = proj.histogram('PWDR my-histogram-name')
    
    >>> # Or by index
    >>> hist = proj.histogram(0)
    >>> assert hist.id == 0
    
    >>> # Or by random id
    >>> assert hist == proj.histogram(hist.ranId)

    Phases can be accessed the same way.

    >>> phase = proj.phase('name of phase')

    New data can also be loaded via :meth:`~G2Project.add_phase` and
    :meth:`~G2Project.add_powder_histogram`.

    >>> hist = proj.add_powder_histogram('some_data_file.chi',
                                         'instrument_parameters.prm')
    >>> phase = proj.add_phase('my_phase.cif', histograms=[hist])

    Parameters for Rietveld refinement can be turned on and off as well.
    See :meth:`~G2Project.set_refinement`, :meth:`~G2Project.clear_refinements`,
    :meth:`~G2Project.iter_refinements`, :meth:`~G2Project.do_refinements`.
    """
    def __init__(self, gpxfile=None, author=None, filename=None, newgpx=None):
        if filename is not None and newgpx is not None:
            raise G2ScriptException('Do not use filename and newgpx together')
        elif newgpx:
            filename = newgpx
        if not gpxfile:
            filename = os.path.abspath(os.path.expanduser(filename))
            self.filename = filename
            self.data, self.names = make_empty_project(author=author, filename=filename)
        elif os.path.exists(os.path.expanduser(gpxfile)):
            # TODO set author
            self.data, self.names = LoadDictFromProjFile(gpxfile)
            self.update_ids()
            if filename:
                filename = os.path.abspath(os.path.expanduser(filename))
                dr = os.path.split(filename)[0]
                if not os.path.exists(dr):
                    raise Exception("Directory {} for filename/newgpx does not exist".format(dr))
                self.filename = filename
            else: 
                self.filename = os.path.abspath(os.path.expanduser(gpxfile))
        else:
            raise ValueError("Not sure what to do with gpxfile {}. Does not exist?".format(gpxfile))

    @classmethod
    def from_dict_and_names(cls, gpxdict, names, filename=None):
        """Creates a :class:`G2Project` directly from
        a dictionary and a list of names. If in doubt, do not use this.

        :returns: a :class:`G2Project`
        """
        out = cls()
        if filename:
            filename = os.path.abspath(os.path.expanduser(filename))
            out.filename = filename
            gpxdict['Controls']['data']['LastSavedAs'] = filename
        else:
            try:
                out.filename = gpxdict['Controls']['data']['LastSavedAs']
            except KeyError:
                out.filename = None
        out.data = gpxdict
        out.names = names

    def save(self, filename=None):
        """Saves the project, either to the current filename, or to a new file.

        Updates self.filename if a new filename provided"""
        # TODO update LastSavedUsing ?
        if filename:
            filename = os.path.abspath(os.path.expanduser(filename))
            self.data['Controls']['data']['LastSavedAs'] = filename
            self.filename = filename
        elif not self.filename:
            raise AttributeError("No file name to save to")
        SaveDictToProjFile(self.data, self.names, self.filename)

    def add_powder_histogram(self, datafile, iparams, phases=[], fmthint=None,
                                 databank=None, instbank=None):
        """Loads a powder data histogram into the project.

        Automatically checks for an instrument parameter file, or one can be
        provided. Note that in unix fashion, "~" can be used to indicate the
        home directory (e.g. ~/G2data/data.fxye).

        Note that the data type (x-ray/CW neutron/TOF) for the histogram 
        will be set from the instrument parameter file. The instrument
        geometry is assumed to be Debye-Scherrer except for 
        dual-wavelength x-ray, where Bragg-Brentano is assumed. 

        :param str datafile: The powder data file to read, a filename.
        :param str iparams: The instrument parameters file, a filename.
        :param list phases: A list of phases to link to the new histogram,
           phases can be references by object, name, rId or number.
           Alternately, use 'all' to link to all phases in the project. 
        :param str fmthint: If specified, only importers where the format name
          (reader.formatName, as shown in Import menu) contains the
          supplied string will be tried as importers. If not specified, all
          importers consistent with the file extension will be tried
          (equivalent to "guess format" in menu).
        :param int databank: Specifies a dataset number to read, if file contains 
          more than set of data. This should be 1 to read the first bank in 
          the file (etc.) regardless of the number on the Bank line, etc.
          Default is None which means there should only be one dataset in the 
          file. 
        :param int instbank: Specifies an instrument parameter set to read, if 
          the instrument parameter file contains more than set of parameters. 
          This will match the INS # in an GSAS type file so it will typically 
          be 1 to read the first parameter set in the file (etc.) 
          Default is None which means there should only be one parameter set 
          in the file.

        :returns: A :class:`G2PwdrData` object representing
            the histogram
        """
        LoadG2fil()
        datafile = os.path.abspath(os.path.expanduser(datafile))
        iparams = os.path.abspath(os.path.expanduser(iparams))
        pwdrreaders = import_generic(datafile, Readers['Pwdr'],fmthint=fmthint,bank=databank)
        histname, new_names, pwdrdata = load_pwd_from_reader(
                                          pwdrreaders[0], iparams,
                                          [h.name for h in self.histograms()],bank=instbank)
        if histname in self.data:
            G2fil.G2Print("Warning - redefining histogram", histname)
        elif self.names[-1][0] == 'Phases':
            self.names.insert(-1, new_names)
        else:
            self.names.append(new_names)
        self.data[histname] = pwdrdata
        self.update_ids()

        if phases == 'all':
            phases = self.phases()
        for phase in phases:
            phase = self.phase(phase)
            self.link_histogram_phase(histname, phase)

        return self.histogram(histname)

    def clone_powder_histogram(self, histref, newname, Y, Yerr=None):
        '''Creates a copy of a powder diffraction histogram with new Y values.
        The X values are not changed. The number of Y values must match the
        number of X values.         

        :param histref: The histogram object, the name of the histogram (str), or ranId 
           or histogram index.
        :param str newname: The name to be assigned to the new histogram
        :param list Y: A set of intensity values
        :param list Yerr: A set of uncertainties for the intensity values (may be None,
           sets all weights to unity)
        :returns: the new histogram object (type G2PwdrData)
        '''
        hist = self.histogram(histref)
        for i in self.names:
            if i[0] == hist.name:
                subkeys = i[1:]
                break
        else:
            raise Exception("error in self.names, hist not found")
        orighist = hist.name
        newhist = 'PWDR '+newname
        if len(Y) != len(self[orighist]['data'][1][0]):
            raise Exception("clone error: length of Y does not match number of X values ({})"
                                .format(len(self[orighist]['data'][1][0])))            
        if Yerr is not None and len(Yerr) != len(self[orighist]['data'][1][0]):
            raise Exception("clone error: length of Yerr does not match number of X values ({})"
                                .format(len(self[orighist]['data'][1][0])))
        
        self[newhist] = copy.deepcopy(self[orighist])
        # intensities
        yo = self[newhist]['data'][1][1] = ma.MaskedArray(Y,mask=self[orighist]['data'][1][1].mask)
        
        Ymin,Ymax = yo.min(),yo.max()
        # set to zero: weights, calc, bkg, obs-calc
        for i in [2,3,4,5]:
            self[newhist]['data'][1][i] *= 0
        # weights
        if Yerr is not None:
            self[newhist]['data'][1][2] += 1./np.array(Yerr)**2
        else:
            self[newhist]['data'][1][2] += 1            # set all weights to 1
        self[newhist]['data'][0] = {'wtFactor': 1.0, 'Dummy': False,
                  'ranId': ran.randint(0, sys.maxsize),
                  'Offset': [0.0, 0.0], 'delOffset': 0.02*Ymax,
                  'refOffset': -0.1*Ymax, 'refDelt': 0.1*Ymax,
                  'Yminmax': [Ymin, Ymax]}
        self[newhist]['Comments'].insert(0,'Cloned from '+orighist)
        self[newhist]['Reflection Lists'] = {}
        self[newhist]['Index Peak List'] = [[], []]
        self[newhist]['Unit Cells List'] = []
        self[newhist]['Peak List'] = {'peaks': [], 'sigDict': {}}
        self.names.append([newhist]+subkeys)
        self.update_ids()
        return self.histogram(newhist)
        
    def add_simulated_powder_histogram(self, histname, iparams, Tmin, Tmax, Tstep=None,
                                       wavelength=None, scale=None, phases=[], ibank=None,
                                           Npoints=None):
        """Create a simulated powder data histogram for the project.

        Requires an instrument parameter file. 
        Note that in unix fashion, "~" can be used to indicate the
        home directory (e.g. ~/G2data/data.prm). The instrument parameter file
        will determine if the histogram is x-ray, CW neutron, TOF, etc. as well
        as the instrument type. 

        :param str histname: A name for the histogram to be created.
        :param str iparams: The instrument parameters file, a filename.
        :param float Tmin: Minimum 2theta or TOF (millisec) for dataset to be simulated
        :param float Tmax: Maximum 2theta or TOF (millisec) for dataset to be simulated
        :param float Tstep: Step size in 2theta or deltaT/T (TOF) for simulated dataset. 
           Default is to compute this from Npoints.
        :param float wavelength: Wavelength for CW instruments, overriding the value
           in the instrument parameters file if specified.
        :param float scale: Histogram scale factor which multiplies the pattern. Note that
           simulated noise is added to the pattern, so that if the maximum intensity is
           small, the noise will mask the computed pattern. The scale 
           needs to be a large number for CW neutrons.
           The default, None, provides a scale of 1 for x-rays and TOF; 10,000 for CW neutrons
           and 100,000 for TOF.
        :param list phases: Phases to link to the new histogram. Use proj.phases() to link to
           all defined phases.
        :param int ibank: provides a bank number for the instrument parameter file. The 
           default is None, corresponding to load the first bank.
        :param int Νpoints: the number of data points to be used for computing the 
            diffraction pattern. Defaults as None, which sets this to 2500. Do not specify
            both Npoints and Tstep. Due to roundoff the actual nuber of points used may differ
            by +-1 from Npoints. Must be below 25,000. 

        :returns: A :class:`G2PwdrData` object representing the histogram
        """
        LoadG2fil()
        iparams = os.path.abspath(os.path.expanduser(iparams))
        if not os.path.exists(iparams):
            raise G2ScriptException("File does not exist:"+iparams)
        rd = G2obj.ImportPowderData( # Initialize a base class reader
            extensionlist=tuple(),
            strictExtension=False,
            formatName = 'Simulate dataset',
            longFormatName = 'Compute a simulated pattern')
        rd.powderentry[0] = '' # no filename
        rd.powderentry[2] = 1 # only one bank
        rd.comments.append('This is a dummy dataset for powder pattern simulation')
        rd.idstring = histname
        #Iparm1, Iparm2 = load_iprms(iparams, rd)
        if Tmax < Tmin:
            Tmin,Tmax = Tmax,Tmin
        if Tstep is not None and Npoints is not None:
            raise G2ScriptException("Error: Tstep and Npoints both specified")
        elif Tstep is not None:
            Tstep = abs(Tstep)
        elif Npoints is None:
            Npoints = 2500
        Iparm1, Iparm2 = load_iprms(iparams, rd, bank=ibank)
        #G2fil.G2Print('Instrument parameters read:',reader.instmsg)
        if 'T' in Iparm1['Type'][0]:
            # patch -- anticipate TOF values in microsec from buggy version
            if Tmax > 200.:
                print('Error: Tmax is too large. Note that input for TOF Tmin & Tmax has changed.')
                print('       Tmin & Tmax are now in milliseconds not microsec. Step is now deltaT/T.')
                raise G2ScriptException("Error: Tmax is too large")
            if Npoints:
                N = Npoints
                Tstep = (np.log(Tmax)-np.log(Tmin))/N
            else:
                N = (np.log(Tmax)-np.log(Tmin))/Tstep
            if N > 25000:
                raise G2ScriptException("Error: Tstep is too small. Would need "+str(N)+" points.")
            x = np.exp((np.arange(0,N))*Tstep+np.log(Tmin*1000.))
            N = len(x)
            unit = 'millisec'
        else:            
            if Npoints:
                N = Npoints
            else:
                N = int((Tmax-Tmin)/Tstep)+1
            if N > 25000:
                raise G2ScriptException("Error: Tstep is too small. Would need "+str(N)+" points.")
            x = np.linspace(Tmin,Tmax,N,True)
            N = len(x)
            unit = 'degrees 2theta'
        if N < 3:
            raise G2ScriptException("Error: Range is too small or step is too large, <3 points")
        G2fil.G2Print('Simulating {} points from {} to {} {}'.format(N,Tmin,Tmax,unit))
        rd.powderdata = [
            np.array(x), # x-axis values
            np.zeros_like(x), # powder pattern intensities
            np.ones_like(x), # 1/sig(intensity)^2 values (weights)
            np.zeros_like(x), # calc. intensities (zero)
            np.zeros_like(x), # calc. background (zero)
            np.zeros_like(x), # obs-calc profiles
            ]
        Tmin = rd.powderdata[0][0]
        Tmax = rd.powderdata[0][-1]
        histname, new_names, pwdrdata = load_pwd_from_reader(rd, iparams,
                                            [h.name for h in self.histograms()],ibank)
        if histname in self.data:
            G2fil.G2Print("Warning - redefining histogram", histname)
        elif self.names[-1][0] == 'Phases':
            self.names.insert(-1, new_names)
        else:
            self.names.append(new_names)
        if scale is not None:
            pwdrdata['Sample Parameters']['Scale'][0] = scale
        elif pwdrdata['Instrument Parameters'][0]['Type'][0].startswith('PNC'):
            pwdrdata['Sample Parameters']['Scale'][0] = 10000.
        elif pwdrdata['Instrument Parameters'][0]['Type'][0].startswith('PNT'):
            pwdrdata['Sample Parameters']['Scale'][0] = 100000.
        self.data[histname] = pwdrdata
        self.update_ids()

        for phase in phases:
            phase = self.phase(phase)
            self.link_histogram_phase(histname, phase)
            
        self.set_Controls('cycles', 0)
        self.data[histname]['data'][0]['Dummy'] = True
        return self.histogram(histname)
    
    def add_phase(self, phasefile, phasename=None, histograms=[], fmthint=None):
        """Loads a phase into the project from a .cif file

        :param str phasefile: The CIF file from which to import the phase.
        :param str phasename: The name of the new phase, or None for the default
        :param list histograms: The names of the histograms to associate with
            this phase. Use proj.histograms() to add to all histograms.
        :param str fmthint: If specified, only importers where the format name
          (reader.formatName, as shown in Import menu) contains the
          supplied string will be tried as importers. If not specified, all
          importers consistent with the file extension will be tried
          (equivalent to "guess format" in menu).

        :returns: A :class:`G2Phase` object representing the
            new phase.
        """
        LoadG2fil()
        histograms = [self.histogram(h).name for h in histograms]
        phasefile = os.path.abspath(os.path.expanduser(phasefile))

        # TODO handle multiple phases in a file
        phasereaders = import_generic(phasefile, Readers['Phase'], fmthint=fmthint)
        phasereader = phasereaders[0]
        
        phasename = phasename or phasereader.Phase['General']['Name']
        phaseNameList = [p.name for p in self.phases()]
        phasename = G2obj.MakeUniqueLabel(phasename, phaseNameList)
        phasereader.Phase['General']['Name'] = phasename

        if 'Phases' not in self.data:
            self.data[u'Phases'] = { 'data': None }
        assert phasename not in self.data['Phases'], "phase names should be unique"
        self.data['Phases'][phasename] = phasereader.Phase

        # process constraints, currently generated only from ISODISTORT CIFs
        if phasereader.Constraints:
            Constraints = self.data['Constraints']
            for i in phasereader.Constraints:
                if isinstance(i, dict):
                    if '_Explain' not in Constraints:
                        Constraints['_Explain'] = {}
                    Constraints['_Explain'].update(i)
                else:
                    Constraints['Phase'].append(i)

        data = self.data['Phases'][phasename]
        generalData = data['General']
        SGData = generalData['SGData']
        NShkl = len(G2spc.MustrainNames(SGData))
        NDij = len(G2spc.HStrainNames(SGData))
        Super = generalData.get('Super', 0)
        if Super:
            SuperVec = np.array(generalData['SuperVec'][0])
        else:
            SuperVec = []
        UseList = data['Histograms']

        for hist in histograms:
            self.link_histogram_phase(hist, phasename)

        for obj in self.names:
            if obj[0] == 'Phases':
                phasenames = obj
                break
        else:
            phasenames = [u'Phases']
            self.names.append(phasenames)
        phasenames.append(phasename)

        # TODO should it be self.filename, not phasefile?
        SetupGeneral(data, os.path.dirname(phasefile))
        self.index_ids()

        self.update_ids()
        return self.phase(phasename)

    def link_histogram_phase(self, histogram, phase):
        """Associates a given histogram and phase.

        .. seealso::

            :meth:`~G2Project.histogram`
            :meth:`~G2Project.phase`"""
        hist = self.histogram(histogram)
        phase = self.phase(phase)

        generalData = phase['General']

        if hist.name.startswith('HKLF '):
            raise NotImplementedError("HKLF not yet supported")
        elif hist.name.startswith('PWDR '):
            hist['Reflection Lists'][generalData['Name']] = {}
            UseList = phase['Histograms']
            SGData = generalData['SGData']
            NShkl = len(G2spc.MustrainNames(SGData))
            NDij = len(G2spc.HStrainNames(SGData))
            UseList[hist.name] = SetDefaultDData('PWDR', hist.name, NShkl=NShkl, NDij=NDij)
            UseList[hist.name]['hId'] = hist.id
            for key, val in [('Use', True), ('LeBail', False),
                             ('newLeBail', True),
                             ('Babinet', {'BabA': [0.0, False],
                                          'BabU': [0.0, False]})]:
                if key not in UseList[hist.name]:
                    UseList[hist.name][key] = val
        else:
            raise RuntimeError("Unexpected histogram" + hist.name)

    def reload(self):
        """Reload self from self.filename"""
        data, names = LoadDictFromProjFile(self.filename)
        self.names = names
        # Need to deep copy the new data file data into the current tree,
        # so that any existing G2Phase, or G2PwdrData objects will still be
        # valid
        _deep_copy_into(from_=data, into=self.data)

    def refine(self, newfile=None, printFile=None, makeBack=False):
        '''Invoke a refinement for the project. The project is written to
        the currently selected gpx file and then either a single or sequential refinement 
        is performed depending on the setting of 'Seq Data' in Controls
        (set in :meth:`get_Controls`). 
        '''
        seqSetting = self.data['Controls']['data'].get('Seq Data',[])
        if not seqSetting:
            self.index_ids()    # index_ids will automatically save the project
            # TODO: migrate to RefineCore G2strMain does not properly use printFile
            # G2strMain.RefineCore(Controls,Histograms,Phases,restraintDict,rigidbodyDict,parmDict,varyList,
            #      calcControls,pawleyLookup,ifPrint,printFile,dlg)
            G2strMain.Refine(self.filename, makeBack=makeBack)
        else:
            self._seqrefine()
        self.reload() # get file from GPX

    def _seqrefine(self):
        '''Perform a sequential refinement.
        '''

        self.data['Controls']['data']['ShowCell'] = True
        # add to tree item to project, if not present
        if 'Sequential results' not in self.data:
            self.data['Sequential results'] = {'data':{}}
            self.names.append(['Sequential results'])
        self.index_ids()    # index_ids will automatically save the project
        #GSASIIpath.IPyBreak_base()

        # check that constraints are OK
        errmsg, warnmsg = G2strIO.ReadCheckConstraints(self.filename)
        if errmsg:
            G2fil.G2Print('Refinement error',errmsg)
            raise Exception('Constraint error')
        if warnmsg:
            G2fil.G2Print(u'Warning: Conflict between refinment flag settings and constraints:\n'+
                  warnmsg+u'\nRefinement not possible')
            raise Exception('Constraint error')
        OK,Msg = G2strMain.SeqRefine(self.filename,None)

    def histogram(self, histname):
        """Returns the histogram named histname, or None if it does not exist.

        :param histname: The name of the histogram (str), or ranId or index.
        :returns: A :class:`G2PwdrData` object, or None if
            the histogram does not exist

        .. seealso::
            :meth:`~G2Project.histograms`
            :meth:`~G2Project.phase`
            :meth:`~G2Project.phases`
        """
        if isinstance(histname, G2PwdrData):
            if histname.proj == self:
                return histname
            else:
                raise Exception('Histogram object (G2PwdrData) is not in current project')
        if histname in self.data:
            return G2PwdrData(self.data[histname], self, histname)
        try:
            # see if histname is an id or ranId
            histname = int(histname)
        except ValueError:
            return

        for histogram in self.histograms():
            if histogram.id == histname or histogram.ranId == histname:
                return histogram

    def histograms(self, typ=None):
        """Return a list of all histograms, as :class:`G2PwdrData` objects

        For now this only finds Powder/Single Xtal histograms, since that is all that is
        currently implemented in this module.

        :param ste typ: The prefix (type) the histogram such as 'PWDR '. If None
          (the default) all known histograms types are found. 
        :returns: a list of objects

        .. seealso::
            :meth:`~G2Project.histogram`
            :meth:`~G2Project.phase`
            :meth:`~G2Project.phases`
        """
        output = []
        # loop through each tree entry. If it is more than one level (more than one item in the
        # list of names). then it must be a histogram, unless labeled Phases or Restraints
        if typ is None:
            for obj in self.names:
                if obj[0].startswith('PWDR ') or obj[0].startswith('HKLF '): 
                    output.append(self.histogram(obj[0]))
        else:
            for obj in self.names:
                if len(obj) > 1 and obj[0].startswith(typ): 
                    output.append(self.histogram(obj[0]))
        return output

    def phase(self, phasename):
        """
        Gives an object representing the specified phase in this project.

        :param str phasename: A reference to the desired phase. Either the phase 
            name (str), the phase's ranId, the phase's index (both int) or
            a phase object (:class:`G2Phase`)
        :returns: A :class:`G2Phase` object
        :raises: KeyError

        .. seealso::
            :meth:`~G2Project.histograms`
            :meth:`~G2Project.phase`
            :meth:`~G2Project.phases`
            """
        if isinstance(phasename, G2Phase):
            if phasename.proj == self:
                return phasename
        phases = self.data['Phases']
        if phasename in phases:
            return G2Phase(phases[phasename], phasename, self)

        try:
            # phasename should be phase index or ranId
            phasename = int(phasename)
        except ValueError:
            return

        for phase in self.phases():
            if phase.id == phasename or phase.ranId == phasename:
                return phase

    def phases(self):
        """
        Returns a list of all the phases in the project.

        :returns: A list of :class:`G2Phase` objects

        .. seealso::
            :meth:`~G2Project.histogram`
            :meth:`~G2Project.histograms`
            :meth:`~G2Project.phase`
            """
        for obj in self.names:
            if obj[0] == 'Phases':
                return [self.phase(p) for p in obj[1:]]
        return []

    def _images(self):
        """Returns a list of all the phases in the project.
        """
        return [i[0] for i in self.names if i[0].startswith('IMG ')]
    
    def image(self, imageRef):
        """
        Gives an object representing the specified image in this project.

        :param str imageRef: A reference to the desired image. Either the Image 
          tree name (str), the image's index (int) or
          a image object (:class:`G2Image`)
        :returns: A :class:`G2Image` object
        :raises: KeyError

        .. seealso::
            :meth:`~G2Project.images`
        """
        if isinstance(imageRef, G2Image):
            if imageRef.proj == self:
                return imageRef
            else:
                raise Exception("Image {} not in current selected project".format(imageRef.name))
        if imageRef in self._images():
            return G2Image(self.data[imageRef], imageRef, self)

        try:
            # imageRef should be an index
            num = int(imageRef)
            imageRef = self._images()[num] 
            return G2Image(self.data[imageRef], imageRef, self)
        except ValueError:
            raise Exception("imageRef {} not an object, name or image index in current selected project"
                                .format(imageRef))
        except IndexError:
            raise Exception("imageRef {} out of range (max={}) in current selected project"
                                .format(imageRef,len(self._images())-1))
            
    def images(self):
        """
        Returns a list of all the images in the project.

        :returns: A list of :class:`G2Image` objects
        """
        return [G2Image(self.data[i],i,self) for i in self._images()]
    
    def _pdfs(self):
        """Returns a list of all the PDF entries in the project.
        """
        return [i[0] for i in self.names if i[0].startswith('PDF ')]
    
    def pdf(self, pdfRef):
        """
        Gives an object representing the specified PDF entry in this project.

        :param pdfRef: A reference to the desired image. Either the PDF 
          tree name (str), the pdf's index (int) or
          a PDF object (:class:`G2PDF`)
        :returns: A :class:`G2PDF` object
        :raises: KeyError

        .. seealso::
            :meth:`~G2Project.pdfs`
            :class:`~G2PDF`
        """
        if isinstance(pdfRef, G2PDF):
            if pdfRef.proj == self:
                return pdfRef
            else:
                raise Exception("PDF {} not in current selected project".format(pdfRef.name))
        if pdfRef in self._pdfs():
            return G2PDF(self.data[pdfRef], pdfRef, self)

        try:
            # pdfRef should be an index
            num = int(pdfRef)
            pdfRef = self._pdfs()[num] 
            return G2PDF(self.data[pdfRef], pdfRef, self)
        except ValueError:
            raise Exception("pdfRef {} not an object, name or PDF index in current selected project"
                                .format(pdfRef))
        except IndexError:
            raise Exception("pdfRef {} out of range (max={}) in current selected project"
                                .format(pdfRef,len(self._images())-1))
    def pdfs(self):
        """
        Returns a list of all the PDFs in the project.

        :returns: A list of :class:`G2PDF` objects
        """
        return [G2PDF(self.data[i],i,self) for i in self._pdfs()]
        
    def copy_PDF(self, PDFobj, histogram):
        '''Creates a PDF entry that can be used to compute a PDF 
        as a copy of settings in an existing PDF (:class:`G2PDF`)
        object. 
        This places an entry in the project but :meth:`G2PDF.calculate` 
        must be used to actually perform the PDF computation. 

        :param PDFobj: A :class:`G2PDF` object which may be 
          in a separate project or the dict associated with the
          PDF object (G2PDF.data).
        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number.
        :returns: A :class:`G2PDF` object for the PDF entry
        '''
        LoadG2fil()
        PDFname = 'PDF ' + self.histogram(histogram).name[4:]
        PDFdict = {'data':None}
        for i in 'PDF Controls', 'PDF Peaks':
            PDFdict[i] = copy.deepcopy(PDFobj[i])
        self.names.append([PDFname]+['PDF Controls', 'PDF Peaks'])
        self.data[PDFname] = PDFdict
        for i in 'I(Q)','S(Q)','F(Q)','G(R)','g(r)':
            self.data[PDFname]['PDF Controls'][i] = []
        G2fil.G2Print('Adding "{}" to project'.format(PDFname))
        return G2PDF(self.data[PDFname], PDFname, self)
        
    def add_PDF(self, prmfile, histogram):
        '''Creates a PDF entry that can be used to compute a PDF. 
        Note that this command places an entry in the project,
        but :meth:`G2PDF.calculate` must be used to actually perform
        the computation. 

        :param str datafile: The powder data file to read, a filename.
        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number.
        :returns: A :class:`G2PDF` object for the PDF entry
        '''
        
        LoadG2fil()
        PDFname = 'PDF ' + self.histogram(histogram).name[4:]
        peaks = {'Limits':[1.,5.],'Background':[2,[0.,-0.2*np.pi],False],'Peaks':[]}
        Controls = {
            'Sample':{'Name':self.histogram(histogram).name,'Mult':1.0},
            'Sample Bkg.':{'Name':'','Mult':-1.0,'Refine':False},
            'Container':{'Name':'','Mult':-1.0,'Refine':False},
            'Container Bkg.':{'Name':'','Mult':-1.0},
            'ElList':{},
            'Geometry':'Cylinder','Diam':1.0,'Pack':0.50,'Form Vol':0.0,'Flat Bkg':0,
            'DetType':'Area detector','ObliqCoeff':0.2,'Ruland':0.025,'QScaleLim':[20,25],
            'Lorch':False,'BackRatio':0.0,'Rmax':100.,'noRing':False,'IofQmin':1.0,'Rmin':1.0,
            'I(Q)':[],'S(Q)':[],'F(Q)':[],'G(R)':[],'g(r)':[]}

        fo = open(prmfile,'r')
        S = fo.readline()
        while S:
            if '#' in S:
                S = fo.readline()
                continue
            key,val = S.split(':',1)
            try:
                Controls[key] = eval(val)
            except:
                Controls[key] = val.strip()
            S = fo.readline()
        fo.close()
        Controls['Sample']['Name'] = self.histogram(histogram).name
        for i in 'Sample Bkg.','Container','Container Bkg.':
            Controls[i]['Name'] = ''
        PDFdict = {'data':None,'PDF Controls':Controls, 'PDF Peaks':peaks}
        self.names.append([PDFname]+['PDF Controls', 'PDF Peaks'])
        self.data[PDFname] = PDFdict
        G2fil.G2Print('Adding "{}" to project'.format(PDFname))
        return G2PDF(self.data[PDFname], PDFname, self)

    def seqref(self):
        """
        Returns a sequential refinement results object, if present

        :returns: A :class:`G2SeqRefRes` object or None if not present
        """
        if 'Sequential results' not in self.data: return
        return G2SeqRefRes(self.data['Sequential results']['data'], self)
    
    def update_ids(self):
        """Makes sure all phases and histograms have proper hId and pId"""
        # Translated from GetUsedHistogramsAndPhasesfromTree,
        #   GSASIIdataGUI.py:4107
        for i, h in enumerate(self.histograms()):
            h.id = i
        for i, p in enumerate(self.phases()):
            p.id = i

    def do_refinements(self, refinements=[{}], histogram='all', phase='all',
                       outputnames=None, makeBack=False):
        """Conducts one or a series of refinements according to the
           input provided in parameter refinements. This is a wrapper
           around :meth:`iter_refinements`

        :param list refinements: A list of dictionaries specifiying changes to be made to
            parameters before refinements are conducted. 
            See the :ref:`Refinement_recipe` section for how this is defined. 
            If not specified, the default value is ``[{}]``, which performs a single 
            refinement step is performed with the current refinement settings. 
        :param str histogram: Name of histogram for refinements to be applied
            to, or 'all'; note that this can be overridden for each refinement
            step via a "histograms" entry in the dict.
        :param str phase: Name of phase for refinements to be applied to, or
            'all'; note that this can be overridden for each refinement
            step via a "phases" entry in the dict.
        :param list outputnames: Provides a list of project (.gpx) file names
            to use for each refinement step (specifying None skips the save step).
            See :meth:`save`. 
            Note that this can be overridden using an "output" entry in the dict.
        :param bool makeBack: determines if a backup ).bckX.gpx) file is made
            before a refinement is performed. The default is False.
            
        To perform a single refinement without changing any parameters, use this
        call:

        .. code-block::  python
        
            my_project.do_refinements([])
        """
        
        for proj in self.iter_refinements(refinements, histogram, phase,
                                          outputnames, makeBack):
            pass
        return self

    def iter_refinements(self, refinements, histogram='all', phase='all',
                         outputnames=None, makeBack=False):
        """Conducts a series of refinements, iteratively. Stops after every
        refinement and yields this project, to allow error checking or
        logging of intermediate results. Parameter use is the same as for
        :meth:`do_refinements` (which calls this method). 

        >>> def checked_refinements(proj):
        ...     for p in proj.iter_refinements(refs):
        ...         # Track intermediate results
        ...         log(p.histogram('0').residuals)
        ...         log(p.phase('0').get_cell())
        ...         # Check if parameter diverged, nonsense answer, or whatever
        ...         if is_something_wrong(p):
        ...             raise Exception("I need a human!")

            
        """
        if outputnames:
            if len(refinements) != len(outputnames):
                raise ValueError("Should have same number of outputs to"
                                 "refinements")
        else:
            outputnames = [None for r in refinements]

        for output, refinedict in zip(outputnames, refinements):
            if 'histograms' in refinedict:
                hist = refinedict['histograms']
            else:
                hist = histogram
            if 'phases' in refinedict:
                ph = refinedict['phases']
            else:
                ph = phase
            if 'output' in refinedict:
                output = refinedict['output']
            self.set_refinement(refinedict, hist, ph)
            # Handle 'once' args - refinements that are disabled after this
            # refinement
            if 'once' in refinedict:
                temp = {'set': refinedict['once']}
                self.set_refinement(temp, hist, ph)

            if output:
                self.save(output)

            if 'skip' not in refinedict:
                self.refine(makeBack=makeBack)
            yield self

            # Handle 'once' args - refinements that are disabled after this
            # refinement
            if 'once' in refinedict:
                temp = {'clear': refinedict['once']}
                self.set_refinement(temp, hist, ph)
            if 'call' in refinedict:
                fxn = refinedict['call']
                if callable(fxn):
                    fxn(*refinedict.get('callargs',[self]))
                elif callable(eval(fxn)):
                    eval(fxn)(*refinedict.get('callargs',[self]))
                else:
                    raise G2ScriptException("Error: call value {} is not callable".format(fxn))

    def set_refinement(self, refinement, histogram='all', phase='all'):
        """Apply specified refinements to a given histogram(s) or phase(s).

        :param dict refinement: The refinements to be conducted
        :param histogram: Specifies either 'all' (default), a single histogram or
          a list of histograms. Histograms may be specified as histogram objects
          (see :class:`G2PwdrData`), the histogram name (str) or the index number (int)
          of the histogram in the project, numbered starting from 0.
          Omitting the parameter or the string 'all' indicates that parameters in
          all histograms should be set. 
        :param phase: Specifies either 'all' (default), a single phase or
          a list of phases. Phases may be specified as phase objects
          (see :class:`G2Phase`), the phase name (str) or the index number (int)
          of the phase in the project, numbered starting from 0.
          Omitting the parameter or the string 'all' indicates that parameters in
          all phases should be set. 

        Note that refinement parameters are categorized as one of three types:

        1. Histogram parameters
        2. Phase parameters
        3. Histogram-and-Phase (HAP) parameters
        
        .. seealso::
            :meth:`G2PwdrData.set_refinements`
            :meth:`G2PwdrData.clear_refinements`
            :meth:`G2Phase.set_refinements`
            :meth:`G2Phase.clear_refinements`
            :meth:`G2Phase.set_HAP_refinements`
            :meth:`G2Phase.clear_HAP_refinements`"""

        if histogram == 'all':
            hists = self.histograms()
        elif isinstance(histogram, list) or isinstance(histogram, tuple):
            hists = []
            for h in histogram:
                if isinstance(h, str) or isinstance(h, int):
                    hists.append(self.histogram(h))
                else:
                    hists.append(h)
        elif isinstance(histogram, str) or isinstance(histogram, int):
            hists = [self.histogram(histogram)]
        else:
            hists = [histogram]

        if phase == 'all':
            phases = self.phases()
        elif isinstance(phase, list) or isinstance(phase, tuple):
            phases = []
            for ph in phase:
                if isinstance(ph, str) or isinstance(ph, int):
                    phases.append(self.phase(ph))
                else:
                    phases.append(ph)
        elif isinstance(phase, str) or isinstance(phase, int):
            phases = [self.phase(phase)]
        else:
            phases = [phase]

        pwdr_set = {}
        phase_set = {}
        hap_set = {}
        for key, val in refinement.get('set', {}).items():
            # Apply refinement options
            if G2PwdrData.is_valid_refinement_key(key):
                pwdr_set[key] = val
            elif G2Phase.is_valid_refinement_key(key):
                phase_set[key] = val
            elif G2Phase.is_valid_HAP_refinement_key(key):
                hap_set[key] = val
            else:
                raise ValueError("Unknown refinement key", key)

        for hist in hists:
            hist.set_refinements(pwdr_set)
        for phase in phases:
            phase.set_refinements(phase_set)
        for phase in phases:
            phase.set_HAP_refinements(hap_set, hists)

        pwdr_clear = {}
        phase_clear = {}
        hap_clear = {}
        for key, val in refinement.get('clear', {}).items():
            # Clear refinement options
            if G2PwdrData.is_valid_refinement_key(key):
                pwdr_clear[key] = val
            elif G2Phase.is_valid_refinement_key(key):
                phase_clear[key] = val
            elif G2Phase.is_valid_HAP_refinement_key(key):
                hap_set[key] = val
            else:
                raise ValueError("Unknown refinement key", key)

        for hist in hists:
            hist.clear_refinements(pwdr_clear)
        for phase in phases:
            phase.clear_refinements(phase_clear)
        for phase in phases:
            phase.clear_HAP_refinements(hap_clear, hists)

    def index_ids(self):
        self.save()
        return G2strIO.GetUsedHistogramsAndPhases(self.filename)

    def add_constraint_raw(self, cons_scope, constr):
        """Adds a constraint of type consType to the project.
        cons_scope should be one of "Hist", "Phase", "HAP", or "Global".

        WARNING it does not check the constraint is well-constructed"""
        constrs = self.data['Constraints']['data']
        if 'Global' not in constrs:
            constrs['Global'] = []
        constrs[cons_scope].append(constr)

    def hold_many(self, vars, type):
        """Apply holds for all the variables in vars, for constraint of a given type.

        type is passed directly to add_constraint_raw as consType

        :param list vars: A list of variables to hold. Either :class:`GSASIIobj.G2VarObj` objects,
            string variable specifiers, or arguments for :meth:`make_var_obj`
        :param str type: A string constraint type specifier. See
            :class:`G2Project.add_constraint_raw`

        """
        for var in vars:
            if isinstance(var, str):
                var = self.make_var_obj(var)
            elif not isinstance(var, G2obj.G2VarObj):
                var = self.make_var_obj(*var)
            self.add_constraint_raw(type, [[1.0, var], None, None, 'h'])

    def make_var_obj(self, phase=None, hist=None, varname=None, atomId=None,
                     reloadIdx=True):
        """Wrapper to create a G2VarObj. Takes either a string representation ("p:h:name:a")
        or individual names of phase, histogram, varname, and atomId.

        Automatically converts string phase, hist, or atom names into the ID required
        by G2VarObj.

        Note that this will cause the project to be saved if not 
        already done so.
        """

        if reloadIdx:
            self.index_ids()
        elif G2obj.TestIndexAll():
            self.index_ids()

        # If string representation, short circuit
        if hist is None and varname is None and atomId is None:
            if isinstance(phase, str) and ':' in phase:
                return G2obj.G2VarObj(phase)

        # Get phase index
        phaseObj = None
        if isinstance(phase, G2Phase):
            phaseObj = phase
            phase = G2obj.PhaseRanIdLookup[phase.ranId]
        elif phase in self.data['Phases']:
            phaseObj = self.phase(phase)
            phaseRanId = phaseObj.ranId
            phase = G2obj.PhaseRanIdLookup[phaseRanId]

        # Get histogram index
        if isinstance(hist, G2PwdrData):
            hist = G2obj.HistRanIdLookup[hist.ranId]
        elif hist in self.data:
            histRanId = self.histogram(hist).ranId
            hist = G2obj.HistRanIdLookup[histRanId]

        # Get atom index (if any)
        if isinstance(atomId, G2AtomRecord):
            atomId = G2obj.AtomRanIdLookup[phase][atomId.ranId]
        elif phaseObj:
            atomObj = phaseObj.atom(atomId)
            if atomObj:
                atomRanId = atomObj.ranId
                atomId = G2obj.AtomRanIdLookup[phase][atomRanId]

        return G2obj.G2VarObj(phase, hist, varname, atomId)

    def add_image(self, imagefile, fmthint=None, defaultImage=None, indexList=None):
        """Load an image into a project

        :param str imagefile: The image file to read, a filename.
        :param str fmthint: If specified, only importers where the format name
          (reader.formatName, as shown in Import menu) contains the
          supplied string will be tried as importers. If not specified, all
          importers consistent with the file extension will be tried
          (equivalent to "guess format" in menu).
        :param str defaultImage: The name of an image to use as a default for 
          setting parameters for the image file to read. 
        :param list indexList: specifies the image numbers (counting from zero) 
          to be used from the file when a file has multiple images. A value of
          ``[0,2,3]`` will cause the only first, third and fourth images in the file
          to be included in the project. 

        :returns: a list of :class:`G2Image` object(s) for the added image(s) 
        """
        LoadG2fil()
        imagefile = os.path.abspath(os.path.expanduser(imagefile))
        readers = import_generic(imagefile, Readers['Image'], fmthint=fmthint)
        objlist = []
        for i,rd in enumerate(readers):
            if indexList is not None and i not in indexList:
                G2fil.G2Print("Image {} skipped".format(i))
                continue
            if rd.SciPy:        #was default read by scipy; needs 1 time fixes
                G2fil.G2Print('Warning: Image {} read by generic SciPy import. Image parameters likely wrong'.format(imagefile))
                #see this: G2IO.EditImageParms(self,rd.Data,rd.Comments,rd.Image,imagefile)
                rd.SciPy = False
            rd.readfilename = imagefile
            if rd.repeatcount == 1 and not rd.repeat: # skip image number if only one in set
                rd.Data['ImageTag'] = None
            else:
                rd.Data['ImageTag'] = rd.repeatcount
            rd.Data['formatName'] = rd.formatName
            if rd.sumfile:
                rd.readfilename = rd.sumfile
            # Load generic metadata, as configured
            G2fil.GetColumnMetadata(rd)
            # Code from G2IO.LoadImage2Tree(rd.readfilename,self,rd.Comments,rd.Data,rd.Npix,rd.Image)
            Imax = np.amax(rd.Image)
            ImgNames = [i[0] for i in self.names if i[0].startswith('IMG ')]
            TreeLbl = 'IMG '+os.path.basename(imagefile)
            ImageTag = rd.Data.get('ImageTag')
            if ImageTag:
                TreeLbl += ' #'+'%04d'%(ImageTag)
                imageInfo = (imagefile,ImageTag)
            else:
                imageInfo = imagefile
            TreeName = G2obj.MakeUniqueLabel(TreeLbl,ImgNames)
            # MT dict to contain image info
            ImgDict = {}
            ImgDict['data'] = [rd.Npix,imageInfo]
            ImgDict['Comments'] = rd.Comments
            if defaultImage:
                if isinstance(defaultImage, G2Image):
                    if defaultImage.proj == self:
                        defaultImage = self.data[defaultImage.name]['data']
                    else:
                        raise Exception("Image {} not in current selected project".format(defaultImage.name))
                elif defaultImage in self._images():
                    defaultImage = self.data[defaultImage]['data']
                else:
                    defaultImage = None
            Data = rd.Data
            if defaultImage:
                Data = copy.copy(defaultImage)
                Data['showLines'] = True
                Data['ring'] = []
                Data['rings'] = []
                Data['cutoff'] = 10.
                Data['pixLimit'] = 20
                Data['edgemin'] = 100000000
                Data['calibdmin'] = 0.5
                Data['calibskip'] = 0
                Data['ellipses'] = []
                Data['calibrant'] = ''
                Data['GonioAngles'] = [0.,0.,0.]
                Data['DetDepthRef'] = False
            else:
                Data['type'] = 'PWDR'
                Data['color'] = GSASIIpath.GetConfigValue('Contour_color','Paired')
                if 'tilt' not in Data:          #defaults if not preset in e.g. Bruker importer
                    Data['tilt'] = 0.0
                    Data['rotation'] = 0.0
                    Data['pixLimit'] = 20
                    Data['calibdmin'] = 0.5
                    Data['cutoff'] = 10.
                Data['showLines'] = False
                Data['calibskip'] = 0
                Data['ring'] = []
                Data['rings'] = []
                Data['edgemin'] = 100000000
                Data['ellipses'] = []
                Data['GonioAngles'] = [0.,0.,0.]
                Data['DetDepth'] = 0.
                Data['DetDepthRef'] = False
                Data['calibrant'] = ''
                Data['IOtth'] = [5.0,50.0]
                Data['LRazimuth'] = [0.,180.]
                Data['azmthOff'] = 0.0
                Data['outChannels'] = 2500
                Data['outAzimuths'] = 1
                Data['centerAzm'] = False
                Data['fullIntegrate'] = GSASIIpath.GetConfigValue('fullIntegrate',True)
                Data['setRings'] = False
                Data['background image'] = ['',-1.0]                            
                Data['dark image'] = ['',-1.0]
                Data['Flat Bkg'] = 0.0
                Data['Oblique'] = [0.5,False]
            Data['varyList'] = {'dist':True,'det-X':True,'det-Y':True,'tilt':True,'phi':True,'dep':False,'wave':False}
            Data['setDefault'] = False
            Data['range'] = [(0,Imax),[0,Imax]]
            ImgDict['Image Controls'] = Data
            ImgDict['Masks'] = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],
                'Frames':[],'Thresholds':[(0,Imax),[0,Imax]],'SpotMask':{'esdMul':3.,'spotMask':None},}
            ImgDict['Stress/Strain']  = {'Type':'True','d-zero':[],'Sample phi':0.0,
                'Sample z':0.0,'Sample load':0.0}
            self.names.append([TreeName]+['Comments','Image Controls','Masks','Stress/Strain'])
            self.data[TreeName] = ImgDict
            del rd.Image
            objlist.append(G2Image(self.data[TreeName], TreeName, self))
        return objlist
    
    def imageMultiDistCalib(self,imageList=None,verbose=False):
        '''Invokes a global calibration fit (same as Image Controls/Calibration/Multi-distance Recalibrate 
        menu command) with images as multiple distance settings. 
        Note that for this to work properly, the initial calibration parameters 
        (center, wavelength, distance & tilts) must be close enough to converge.
        This may produce a better result if run more than once.

        See :ref:`MultiDist_Example` for example code.

        :param str imageList: the images to include in the fit, if not specified 
          all images in the project will be included.

        :returns: parmDict,covData where parmDict has the refined parameters 
          and their values and covData is a dict containing the covariance matrix ('covMatrix'),
          the number of ring picks ('obs') the reduced Chi-squared ('chisq'),
          the names of the variables ('varyList') and their values ('variables')
        '''
        if imageList is None:
            imageList = self.images()

        # code based on GSASIIimgGUI..OnDistRecalib
        obsArr = np.array([]).reshape(0,4)
        parmDict = {}
        varList = []
        HKL = {}

        for img in imageList:
            name = img.name
            G2fil.G2Print ('getting rings for',name)
            Data = img.data['Image Controls']
            key = str(int(Data['setdist']))
            # create a parameter dict for combined fit
            if 'wavelength' not in parmDict:
                parmDict['wavelength'] = Data['wavelength']
                if Data['varyList']['wave']:
                    varList += ['wavelength']
                if Data['varyList']['dist']:
                    raise Exception(
                                'You cannot vary individual detector positions and the global wavelength.\n\nChange flags for 1st image.',
                                'Conflicting vars')
                parmDict['dep'] = Data['DetDepth']
                if Data['varyList']['dep']:
                    varList += ['dep']
                # distance flag determines if individual values are refined
                if not Data['varyList']['dist']:
                    # starts as zero, single variable, always refined
                    parmDict['deltaDist'] = 0.
                    varList += ['deltaDist']
                parmDict['phi'] = Data['rotation']
                if Data['varyList']['phi']:
                    varList += ['phi']
                parmDict['tilt'] = Data['tilt']
                if Data['varyList']['tilt']:
                    varList += ['tilt']

            ImageZ = _getCorrImage(Readers['Image'],self,img)
            Data['setRings'] = True
            Masks = img.data['Masks']
            result = G2img.ImageRecalibrate(None,ImageZ,Data,Masks,getRingsOnly=True)
            if not len(result):
                raise Exception('calibrant missing from local image calibrants files')
            rings,HKL[key] = result
            # add detector set dist into data array, create a single really large array
            distarr = np.zeros_like(rings[:,2:3])
            if 'setdist' not in Data:
                raise Exception('Distance (setdist) not in image metadata')
            distarr += Data['setdist']
            obsArr = np.concatenate((
                        obsArr,
                        np.concatenate((rings[:,0:2],distarr,rings[:,2:3]),axis=1)),axis=0)
            if 'deltaDist' not in parmDict:
                # starts as zero, variable refined for each image
                parmDict['delta'+key] = 0
                varList += ['delta'+key]
            for i,z in enumerate(['X','Y']):
                v = 'det-'+z
                if v+key in parmDict:
                    raise Exception('Error: two images with setdist ~=',key)
                parmDict[v+key] = Data['center'][i]
                if Data['varyList'][v]:
                    varList += [v+key]
        #GSASIIpath.IPyBreak()
        G2fil.G2Print('\nFitting',obsArr.shape[0],'ring picks and',len(varList),'variables...')
        result = G2img.FitMultiDist(obsArr,varList,parmDict,covar=True,Print=verbose)
        
        for img in imageList: # update GPX info with fit results
            name = img.name
            #print ('updating',name)
            Data = img.data['Image Controls']
            Data['wavelength'] = parmDict['wavelength']
            key = str(int(Data['setdist']))
            Data['center'] = [parmDict['det-X'+key],parmDict['det-Y'+key]]
            if 'deltaDist' in parmDict:
                Data['distance'] = Data['setdist'] - parmDict['deltaDist']
            else:
                Data['distance'] = Data['setdist'] - parmDict['delta'+key]
            Data['rotation'] = np.mod(parmDict['phi'],360.0)
            Data['tilt'] = parmDict['tilt']
            Data['DetDepth'] = parmDict['dep']
            N = len(Data['ellipses'])
            Data['ellipses'] = []           #clear away individual ellipse fits
            for H in HKL[key][:N]:
                ellipse = G2img.GetEllipse(H[3],Data)
                Data['ellipses'].append(copy.deepcopy(ellipse+('b',)))
                
        covData = {'title':'Multi-distance recalibrate','covMatrix':result[3],
                       'obs':obsArr.shape[0],'chisq':result[0],
                       'varyList':varList,'variables':result[1]}
        return parmDict,covData
                
    def set_Frozen(self, variable=None, histogram=None, mode='remove'):
        '''Removes one or more Frozen variables (or adds one)
        (See :ref:`Parameter Limits<ParameterLimits>` description.)
        Note that use of this 
        will cause the project to be saved if not already done so.

        :param str variable: a variable name as a str or 
          (as a :class:`GSASIIobj.G2VarObj` object). Should
          not contain wildcards. 
          If None (default), all frozen variables are deleted
          from the project, unless a sequential fit and 
          a histogram is specified. 
        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number.
          Used for sequential fits only. 
        :param str mode: The default mode is to remove variables
          from the appropriate Frozen list, but if the mode
          is specified as 'add', the variable is added to the 
          list. 
        :returns: True if the variable was added or removed, False 
          otherwise. Exceptions are generated with invalid requests.
        '''
        Controls = self.data['Controls']['data']
        for key in ('parmMinDict','parmMaxDict','parmFrozen'):
            if key not in Controls: Controls[key] = {}
        if G2obj.TestIndexAll(): self.index_ids()
        patchControls(Controls)

        if mode == 'remove':
            if variable is None and histogram is None:
                Controls['parmFrozen'] = {}
                return True
        elif mode == 'add':
            if variable is None:
                raise Exception('set_Frozen error: variable must be specified')
        else:
            raise Exception('Undefined mode ({}) in set_Frozen'.format(mode))
        
        if histogram is None:
            h = 'FrozenList'
        else:
            hist = self.histogram(histogram)
            if hist:
                h = hist.name
                if h not in Controls['parmFrozen']: # no such list, not found
                    return False
            else:
                raise Exception('set_Frozen error: histogram {} not found'.format(histogram))
        if mode == 'remove':
            if variable is None:
                Controls['parmFrozen'][h] = []
                return True
            if h not in Controls['parmFrozen']:
                return True
            delList = []
            for i,v in enumerate(Controls['parmFrozen'][h]):
                if v == variable: delList.append(i)
            if delList:
                for i in reversed(delList):
                    del Controls['parmFrozen'][h][i]
                return True
            return False
        elif mode == 'add':
            if type(variable) is str:
                variable = G2obj.G2VarObj(variable)
            elif type(v) is not G2obj.G2VarObj:
                raise Exception(
                    'set_Frozen error: variable {} wrong type ({})'
                    .format(variable,type(variable)))
            if h not in Controls['parmFrozen']: Controls['parmFrozen'][h] = []
            Controls['parmFrozen'][h].append(variable)
            return True

    def get_Frozen(self, histogram=None):
        '''Gets a list of Frozen variables. 
        (See :ref:`Parameter Limits<ParameterLimits>` description.)
        Note that use of this 
        will cause the project to be saved if not already done so.

        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number. Used
          for sequential fits only. If left as the default (None)
          for a sequential fit, all Frozen variables in all
          histograms are returned.
        :returns: a list containing variable names, as str values
        '''
        Controls = self.data['Controls']['data']
        for key in ('parmMinDict','parmMaxDict','parmFrozen'):
            if key not in Controls: Controls[key] = {}
        if G2obj.TestIndexAll(): self.index_ids()
        patchControls(Controls)
        if histogram:
            hist = self.histogram(histogram)
            if hist:
                h = hist.name
                return [str(i) for i in Controls['parmFrozen'][h]]
            elif histogram.lower() == 'all':
                names = set()
                for h in self.histograms():
                    if h.name not in Controls['parmFrozen']: continue
                    names = names.union([str(i) for i in Controls['parmFrozen'][h.name]])
                return list(names)
            else:
                raise Exception('histogram {} not recognized'.format(histogram))
        elif 'FrozenList' in Controls['parmFrozen']:
            return [str(i) for i in Controls['parmFrozen']['FrozenList']]
        else:
            return []
        
    def get_Controls(self, control, variable=None):
        '''Return project controls settings

        :param str control: the item to be returned. See below for allowed values.
        :param str variable: a variable name as a str or 
          (as a :class:`GSASIIobj.G2VarObj` object). 
          Used only with control set to "parmMin" or "parmMax".
        :returns: The value for the control.

        Allowed values for parameter control:

            * cycles: the maximum number of cycles (returns int)
            * sequential: the histograms used for a sequential refinement as a list
              of histogram names or an empty list when in non-sequential mode.
            * Reverse Seq: returns True or False. True indicates that fitting of the
              sequence of histograms proceeds in reversed order.
            * seqCopy: returns True or False. True indicates that results from 
              each sequential fit are used as the starting point for the next
              histogram.
            * parmMin & parmMax: retrieves a maximum or minimum value for 
              a refined parameter. Note that variable will be a GSAS-II 
              variable name, optionally with * specified for a histogram 
              or atom number. Return value will be a float.
              (See :ref:`Parameter Limits<ParameterLimits>` description.)
            * Anything else returns the value in the Controls dict, if present. An 
              exception is raised if the control value is not present. 

        .. seealso::
            :meth:`set_Controls`
        '''
        if control.startswith('parmM'):
            if not variable:
                raise Exception('set_Controls requires a value for variable for control=parmMin or parmMax')
            for key in ('parmMinDict','parmMaxDict','parmFrozen'):
                if key not in self.data['Controls']['data']: self.data['Controls']['data'][key] = {}
            if G2obj.TestIndexAll(): self.index_ids()
            patchControls(self.data['Controls']['data'])
        if control == 'cycles':
            return self.data['Controls']['data']['max cyc']
        elif control == 'sequential':
            return self.data['Controls']['data']['Seq Data']
        elif control == 'Reverse Seq':
            return self.data['Controls']['data']['Reverse Seq']
        elif control == 'seqCopy':
            return self.data['Controls']['data']['Copy2Next']
        elif control == 'parmMin' or control == 'parmMax':
            key = G2obj.G2VarObj(variable)
            return G2obj.prmLookup(
                variable,
                self.data['Controls']['data'][control+'Dict'])[1]
        elif control in self.data['Controls']['data']:
            return self.data['Controls']['data'][control]
        else:
            G2fil.G2Print('Defined Controls:',self.data['Controls']['data'].keys())
            raise Exception('{} is an invalid control value'.format(control))
        
    def set_Controls(self, control, value, variable=None):
        '''Set project controls.

        Note that use of this with control set to parmMin or parmMax 
        will cause the project to be saved if not already done so.

        :param str control: the item to be set. See below for allowed values. 
        :param value: the value to be set.
        :param str variable: used only with control set to "parmMin" or "parmMax"

        Allowed values for *control* parameter:

        * ``'cycles'``: sets the maximum number of cycles (value must be int)
        * ``'sequential'``: sets the histograms to be used for a sequential refinement. 
          Use an empty list to turn off sequential fitting. 
          The values in the list may be the name of the histogram (a str), or 
          a ranId or index (int values), see :meth:`histogram`.
        * ``'seqCopy'``: when True, the results from each sequential fit are used as
          the starting point for the next. After each fit is is set to False. 
          Ignored for non-sequential fits. 
        * ``'Reverse Seq'``: when True, sequential refinement is performed on the
          reversed list of histograms.
        * ``'parmMin'`` & ``'parmMax'``: set a maximum or minimum value for a refined 
          parameter. Note that variable will be a GSAS-II variable name, 
          optionally with * specified for a histogram or atom number and
          value must be a float. 
          (See :ref:`Parameter Limits<ParameterLimits>` description.)

        .. seealso::
            :meth:`get_Controls`
        '''
        if control.startswith('parmM'):
            if not variable:
                raise Exception('set_Controls requires a value for variable for control=parmMin or parmMax')
            for key in ('parmMinDict','parmMaxDict','parmFrozen'):
                if key not in self.data['Controls']['data']: self.data['Controls']['data'][key] = {}
            if G2obj.TestIndexAll(): self.index_ids()
            patchControls(self.data['Controls']['data'])
        if control == 'cycles':
            self.data['Controls']['data']['max cyc'] = int(value)
        elif control == 'seqCopy':
            self.data['Controls']['data']['Copy2Next'] = bool(value)
        elif control == 'Reverse Seq':
            self.data['Controls']['data']['Reverse Seq'] = bool(value)
        elif control == 'sequential':
            histlist = []
            for i,j in enumerate(value):
                h = self.histogram(j)
                if h:
                    histlist.append(h.name)
                else:
                    raise Exception('item #{} ({}) is an invalid histogram value'
                                        .format(i,j))
            self.data['Controls']['data']['Seq Data'] = histlist
        elif control == 'parmMin' or control == 'parmMax':
            key = G2obj.G2VarObj(variable)
            self.data['Controls']['data'][control+'Dict'][key] = float(value)
        else:
            raise Exception('{} is an invalid control value'.format(control))
        
    def copyHistParms(self,sourcehist,targethistlist='all',modelist='all'):
        '''Copy histogram information from one histogram to others

        :param sourcehist: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram
        :param list targethistlist: a list of histograms where each item in the
            list can be a histogram object (:class:`G2PwdrData`), 
            a histogram name or the index number of the histogram.
            if the string 'all' (default value), then all histograms in 
            the project are used. 
        :param list modelist: May be a list of sections to copy, which may 
           include 'Background', 'Instrument Parameters', 'Limits' and 
           'Sample Parameters' (items may be shortened to uniqueness and
           capitalization is ignored, so ['b','i','L','s'] will work.)
           The default value, 'all' causes the listed sections to 

        '''
        sections = ('Background','Instrument Parameters','Limits',
                        'Sample Parameters')
        hist_in = self.histogram(sourcehist)
        if not hist_in:
            raise Exception('{} is not a valid histogram'.format(sourcehist))
        if targethistlist == "all":
            targethistlist = self.histograms()
        if 'all' in modelist:
            copysections = sections
        else:
            copysections = set()
            for s in sections:
                for m in modelist:
                    if s.lower().startswith(m.lower()):
                        copysections.add(s)
        for h in targethistlist:
            hist_out = self.histogram(h)
            if not hist_out:
                raise Exception('{} is not a valid histogram'.format(h))
            for key in copysections: 
                hist_out[key] = copy.deepcopy(hist_in[key])

    def get_VaryList(self):
        '''Returns a list of the refined variables in the 
        last refinement cycle

        :returns: a list of variables or None if no refinement has been
          performed. 
        '''
        try:
            return self['Covariance']['data']['varyList']
        except:
            return

    def get_ParmList(self):
        '''Returns a list of all the parameters defined in the 
        last refinement cycle

        :returns: a list of parameters or None if no refinement has been
          performed. 
        '''
        try:
            return list(self['Covariance']['data']['parmDict'].keys())
        except:
            return
        
    def get_Variable(self,var):
        '''Returns the value and standard uncertainty (esd) for a variable 
        parameters, as defined in the last refinement cycle

        :param str var: a variable name of form '<p>:<h>:<name>', such as 
          ':0:Scale'
        :returns: (value,esd) if the parameter is refined or 
          (value, None) if the variable is in a constraint or is not 
          refined or None if the parameter is not found. 
        '''
        if var not in self['Covariance']['data']['parmDict']:
            return None
        val = self['Covariance']['data']['parmDict'][var]
        try:
            pos = self['Covariance']['data']['varyList'].index(var)
            esd = np.sqrt(self['Covariance']['data']['covMatrix'][pos,pos])
            return (val,esd)
        except ValueError:
            return (val,None)

    def get_Covariance(self,varList):
        '''Returns the values and covariance matrix for a series of variable 
        parameters. as defined in the last refinement cycle

        :param tuple varList: a list of variable names of form '<p>:<h>:<name>'
        :returns: (valueList,CovMatrix) where valueList contains the (n) values
          in the same order as varList (also length n) and CovMatrix is a 
          (n x n) matrix. If any variable name is not found in the varyList
          then None is returned. 

        Use this code, where sig provides standard uncertainties for 
        parameters and where covArray provides the correlation between 
        off-diagonal terms::

            sig = np.sqrt(np.diag(covMatrix))
            xvar = np.outer(sig,np.ones_like(sig))
            covArray = np.divide(np.divide(covMatrix,xvar),xvar.T)

        '''
        missing = [i for i in varList if i not in self['Covariance']['data']['varyList']]
        if missing:
            G2fil.G2Print('Warning: Variable(s) {} were not found in the varyList'.format(missing))
            return None
        vals = [self['Covariance']['data']['parmDict'][i] for i in varList]
        import GSASIImath as G2mth
        cov = G2mth.getVCov(varList,
                            self['Covariance']['data']['varyList'],
                            self['Covariance']['data']['covMatrix'])
        return (vals,cov)

class G2AtomRecord(G2ObjectWrapper):
    """Wrapper for an atom record. Has convenient accessors via @property:
    label, type, refinement_flags, coordinates, occupancy, ranId, id, adp_flag, uiso

    Example:

    >>> atom = some_phase.atom("O3")
    >>> # We can access the underlying data format
    >>> atom.data
    ['O3', 'O-2', '', ... ]
    >>> # We can also use wrapper accessors
    >>> atom.coordinates
    (0.33, 0.15, 0.5)
    >>> atom.refinement_flags
    u'FX'
    >>> atom.ranId
    4615973324315876477
    >>> atom.occupancy
    1.0
    """
    def __init__(self, data, indices, proj):
        self.data = data
        self.cx, self.ct, self.cs, self.cia = indices
        self.proj = proj

    @property
    def label(self):
        '''Get the associated atom's label
        '''
        return self.data[self.ct-1]

    @property
    def type(self):
        '''Get the associated atom's type
        '''
        return self.data[self.ct]

    @property
    def element(self):
        '''Get the associated atom's element symbol
        '''
        import re
        try:
            return re.match('^([A-Z][a-z]?)',self.data[self.ct]).group(1)
        except:
            raise Exception("element parse error with type {}".
                                format(self.data[self.ct]))

    @property
    def refinement_flags(self):
        '''Get or set refinement flags for the associated atom
        '''
        return self.data[self.ct+1]

    @refinement_flags.setter
    def refinement_flags(self, other):
        # Automatically check it is a valid refinement
        for c in other:
            if c not in ' FXU':
                raise ValueError("Invalid atom refinement: ", other)
        self.data[self.ct+1] = other

    @property
    def coordinates(self):
        '''Get the associated atom's coordinates
        '''
        return tuple(self.data[self.cx:self.cx+3])

    @property
    def occupancy(self):
        '''Get or set the associated atom's occupancy fraction
        '''
        return self.data[self.cx+3]
    
    @occupancy.setter
    def occupancy(self, val):
        self.data[self.cx+3] = float(val)

    @property
    def mult(self):
        '''Get the associated atom's multiplicity value
        '''
        return self.data[self.cs+1]
    
    @property
    def ranId(self):
        '''Get the associated atom's Random Id number
        '''
        return self.data[self.cia+8]

    @property
    def adp_flag(self):
        '''Get the associated atom's iso/aniso setting, 'I' or 'A'
        '''
        # Either 'I' or 'A'
        return self.data[self.cia]

    @property
    def uiso(self):
        '''Get or set the associated atom's Uiso or Uaniso value(s)
        '''
        if self.adp_flag == 'I':
            return self.data[self.cia+1]
        else:
            return self.data[self.cia+2:self.cia+8]

    @uiso.setter
    def uiso(self, value):
        if self.adp_flag == 'I':
            self.data[self.cia+1] = float(value)
        else:
            assert len(value) == 6
            self.data[self.cia+2:self.cia+8] = [float(v) for v in value]

class G2PwdrData(G2ObjectWrapper):
    """Wraps a Powder Data Histogram. 
    The object contains these class variables:

        * G2PwdrData.proj: contains a reference to the :class:`G2Project`
          object that contains this histogram 
        * G2PwdrData.name: contains the name of the histogram
        * G2PwdrData.data: contains the histogram's associated data in a dict,
          as documented for the :ref:`Powder Diffraction Tree<Powder_table>`.
          The actual histogram values are contained in the 'data' dict item,
          as documented for Data. 

    """
    def __init__(self, data, proj, name):
        self.data = data
        self.proj = proj
        self.name = name

    @staticmethod
    def is_valid_refinement_key(key):
        valid_keys = ['Limits', 'Sample Parameters', 'Background',
                      'Instrument Parameters']
        return key in valid_keys

    #@property
    #def name(self):
    #    return self.data['data'][-1]

    @property
    def ranId(self):
        return self.data['data'][0]['ranId']

    @property
    def residuals(self):
        '''Provides a dictionary with with the R-factors for this histogram.
        Includes the weighted and unweighted profile terms (R, Rb, wR, wRb, wRmin)
        as well as the Bragg R-values for each phase (ph:H:Rf and ph:H:Rf^2).
        '''
        data = self.data['data'][0]
        return {key: data[key] for key in data
                if key in ['R', 'Rb', 'wR', 'wRb', 'wRmin']
                   or ':' in key}

    @property
    def InstrumentParameters(self):
        '''Provides a dictionary with with the Instrument Parameters
        for this histogram.
        '''
        return self.data['Instrument Parameters'][0]

    @property
    def SampleParameters(self):
        '''Provides a dictionary with with the Sample Parameters
        for this histogram.
        '''
        return self.data['Sample Parameters']

    @property
    def Background(self):
        '''Provides a list with with the Background parameters
        for this histogram.

        :returns: list containing a list and dict with background values
        '''
        return self.data['Background']

    def add_back_peak(self,pos,int,sig,gam,refflags=[]):
        '''Adds a background peak to the Background parameters
        
        :param float pos: position of peak, a 2theta or TOF value
        :param float int: integrated intensity of background peak, usually large
        :param float sig: Gaussian width of background peak, usually large
        :param float gam: Lorentzian width of background peak, usually unused (small)
        :param list refflags: a list of 1 to 4 boolean refinement flags for 
            pos,int,sig & gam, respectively (use [0,1] to refine int only). 
            Defaults to [] which means nothing is refined. 
        '''
        if 'peaksList' not in self.Background[1]:
            self.Background[1]['peaksList'] = []
        flags = 4*[False]
        for i,f in enumerate(refflags):
            if i>3: break
            flags[i] = bool(f)
        bpk = []
        for i,j in zip((pos,int,sig,gam),flags):
            bpk += [float(i),j]
        self.Background[1]['peaksList'].append(bpk)
        self.Background[1]['nPeaks'] = len(self.Background[1]['peaksList'])

    def del_back_peak(self,peaknum):
        '''Removes a background peak from the Background parameters
        
        :param int peaknum: the number of the peak (starting from 0)
        '''
        npks = self.Background[1].get('nPeaks',0)
        if peaknum >= npks:
            raise Exception('peak {} not found in histogram {}'.format(peaknum,self.name))
        del self.Background[1]['peaksList'][peaknum]
        self.Background[1]['nPeaks'] = len(self.Background[1]['peaksList'])
        
    def ref_back_peak(self,peaknum,refflags=[]):
        '''Sets refinement flag for a background peak
        
        :param int peaknum: the number of the peak (starting from 0)
        :param list refflags: a list of 1 to 4 boolean refinement flags for 
            pos,int,sig & gam, respectively. If a flag is not specified
            it defaults to False (use [0,1] to refine int only). 
            Defaults to [] which means nothing is refined. 
        '''
        npks = self.Background[1].get('nPeaks',0)
        if peaknum >= npks:
            raise Exception('peak {} not found in histogram {}'.format(peaknum,self.name))
        flags = 4*[False]
        for i,f in enumerate(refflags):
            if i>3: break
            flags[i] = bool(f)
        for i,f in enumerate(flags):
            self.Background[1]['peaksList'][peaknum][2*i+1] = f
                    
    @property
    def id(self):
        self.proj.update_ids()
        return self.data['data'][0]['hId']

    @id.setter
    def id(self, val):
        self.data['data'][0]['hId'] = val

    def fit_fixed_points(self):
        """Attempts to apply a background fit to the fixed points currently specified."""
        def SetInstParms(Inst):
            dataType = Inst['Type'][0]
            insVary = []
            insNames = []
            insVals = []
            for parm in Inst:
                insNames.append(parm)
                insVals.append(Inst[parm][1])
                if parm in ['U','V','W','X','Y','Z','SH/L','I(L2)/I(L1)','alpha',
                    'beta-0','beta-1','beta-q','sig-0','sig-1','sig-2','sig-q',] and Inst[parm][2]:
                        Inst[parm][2] = False
            instDict = dict(zip(insNames, insVals))
            instDict['X'] = max(instDict['X'], 0.01)
            instDict['Y'] = max(instDict['Y'], 0.01)
            if 'SH/L' in instDict:
                instDict['SH/L'] = max(instDict['SH/L'], 0.002)
            return dataType, instDict, insVary

        bgrnd = self.data['Background']

        # Need our fixed points in order
        bgrnd[1]['FixedPoints'].sort(key=lambda pair: pair[0])
        X = [x for x, y in bgrnd[1]['FixedPoints']]
        Y = [y for x, y in bgrnd[1]['FixedPoints']]

        limits = self.data['Limits'][1]
        if X[0] > limits[0]:
            X = [limits[0]] + X
            Y = [Y[0]] + Y
        if X[-1] < limits[1]:
            X += [limits[1]]
            Y += [Y[-1]]

        # Some simple lookups
        controls = self.proj['Controls']['data']
        inst, inst2 = self.data['Instrument Parameters']
        pwddata = self.data['data'][1]

        # Construct the data for background fitting
        xBeg = np.searchsorted(pwddata[0], limits[0])
        xFin = np.searchsorted(pwddata[0], limits[1])
        xdata = pwddata[0][xBeg:xFin]
        ydata = si.interp1d(X,Y)(ma.getdata(xdata))

        W = [1]*len(xdata)
        Z = [0]*len(xdata)

        dataType, insDict, insVary = SetInstParms(inst)
        bakType, bakDict, bakVary = G2pwd.SetBackgroundParms(bgrnd)

        # Do the fit
        data = np.array([xdata, ydata, W, Z, Z, Z])
        G2pwd.DoPeakFit('LSQ', [], bgrnd, limits, inst, inst2, data,
                        prevVaryList=bakVary, controls=controls)

        # Post-fit
        parmDict = {}
        bakType, bakDict, bakVary = G2pwd.SetBackgroundParms(bgrnd)
        parmDict.update(bakDict)
        parmDict.update(insDict)
        pwddata[3][xBeg:xFin] *= 0
        pwddata[5][xBeg:xFin] *= 0
        pwddata[4][xBeg:xFin] = G2pwd.getBackground('', parmDict, bakType, dataType, xdata)[0]

        # TODO adjust pwddata? GSASIIpwdGUI.py:1041
        # TODO update background
        self.proj.save()

    def getdata(self,datatype):
        '''Provides access to the histogram data of the selected data type

        :param str datatype: must be one of the following values (case is ignored)
        
           * 'X': the 2theta or TOF values for the pattern
           * 'Yobs': the observed intensity values 
           * 'Yweight': the weights for each data point (1/sigma**2)
           * 'Ycalc': the computed intensity values 
           * 'Background': the computed background values
           * 'Residual': the difference between Yobs and Ycalc (obs-calc)

        :returns: an numpy MaskedArray with data values of the requested type
        
        '''
        enums = ['x', 'yobs', 'yweight', 'ycalc', 'background', 'residual']
        if datatype.lower() not in enums:
            raise G2ScriptException("Invalid datatype = "+datatype+" must be one of "+str(enums))
        return self.data['data'][1][enums.index(datatype.lower())]
        
    def y_calc(self):
        '''Returns the calculated intensity values; better to 
        use :meth:`getdata`
        '''
        return self.data['data'][1][3]

    def reflections(self):
        '''Returns a dict with an entry for every phase in the 
        current histogram. Within each entry is a dict with keys
        'RefList' (reflection list, see 
        :ref:`Powder Reflections <PowderRefl_table>`), 
        'Type' (histogram type), 'FF' 
        (form factor information), 'Super' (True if this is superspace 
        group). 
        '''
        return self.data['Reflection Lists']
    
    def Export(self,fileroot,extension,fmthint=None):
        '''Write the histogram into a file. The path is specified by fileroot and
        extension.
        
        :param str fileroot: name of the file, optionally with a path (extension is
           ignored)
        :param str extension: includes '.', must match an extension in global
           exportersByExtension['powder'] or a Exception is raised.
        :param str fmthint: If specified, the first exporter where the format 
           name (obj.formatName, as shown in Export menu) contains the
           supplied string will be used. If not specified, an error 
           will be generated showing the possible choices.
        :returns: name of file that was written
        '''
        LoadG2fil()
        if extension not in exportersByExtension.get('powder',[]):
            print('Defined exporters are:')
            print('  ',list(exportersByExtension.get('powder',[])))
            raise G2ScriptException('No Writer for file type = "'+extension+'"')
        fil = os.path.abspath(os.path.splitext(fileroot)[0]+extension)
        obj = exportersByExtension['powder'][extension]
        if type(obj) is list:
            if fmthint is None:
                print('Defined ',extension,'exporters are:')
                for o in obj:
                    print('\t',o.formatName)
                raise G2ScriptException('No format hint for file type = "'+extension+'"')
            for o in obj:
              if fmthint.lower() in o.formatName.lower():
                  obj = o
                  break
            else:
                print('Hint ',fmthint,'not found. Defined ',extension,'exporters are:')
                for o in obj:
                    print('\t',o.formatName)
                raise G2ScriptException('Bad format hint for file type = "'+extension+'"')
        self._SetFromArray(obj)
        obj.Writer(self.name,fil)
        return fil
    
    def _SetFromArray(self,expObj):
        '''Load a histogram into the exporter in preparation for use of 
        the .Writer method in the object. 

        :param Exporter expObj: Exporter object
        '''
        expObj.Histograms[self.name] =  {}
        expObj.Histograms[self.name]['Data'] = self.data['data'][1]
        for key in 'Instrument Parameters','Sample Parameters','Reflection Lists':
            expObj.Histograms[self.name][key] = self.data[key]
            
    def plot(self, Yobs=True, Ycalc=True, Background=True, Residual=True):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            G2fil.G2Print('Warning: Unable to import matplotlib, skipping plot')
            return
        data = self.data['data'][1]
        if Yobs:
            plt.plot(data[0], data[1], label='Yobs')
        if Ycalc:
            plt.plot(data[0], data[3], label='Ycalc')
        if Background:
            plt.plot(data[0], data[4], label='Background')
        if Residual:
            plt.plot(data[0], data[5], label="Residual")

    def get_wR(self):
        """returns the overall weighted profile R factor for a histogram
        
        :returns: a wR value as a percentage or None if not defined
        """
        return self['data'][0].get('wR')

    def _decodeHist(self,hist):
        '''Convert a histogram reference to a histogram name string
        '''
        if isinstance(hist, G2PwdrData):
            return hist.name
        elif hist in [h.name for h in self.proj.histograms()]:
            return hist
        elif type(hist) is int:
            return self.proj.histograms()[hist].name
        else:
            raise G2ScriptException("Invalid histogram reference: "+str(hist))
        
    def set_background(self, key, value):
        '''Set background parameters (this serves a similar function as in 
        :meth:`set_refinements`, but with a simplified interface). 

        :param str key: a string that defines the background parameter that will 
           be changed. Must appear in the table below.

           =================   ==============   ===========================================
           key name            type of value     meaning of value
           =================   ==============   ===========================================
           fixedHist           int, str,         reference to a histogram in the current
                               None or           project or None to remove the reference.
                               G2PwdrData        
           fixedFileMult       float             multiplier applied to intensities in 
                                                 the background histogram where a value 
                                                 of -1.0 means full subtraction of 
                                                 the background histogram.
           =================   ==============   ===========================================

        :param value: a value to set the selected background parameter. The meaning 
           and type for this parameter is listed in the table above.

        '''
        bkgPrms, bkgDict = self.data['Background']
        if key == 'fixedHist':
            if value is None:
                bkgDict['background PWDR'][0] = ''
                return
            bkgDict['background PWDR'][0] = self._decodeHist(value)
        elif key == 'fixedFileMult':
            bkgDict['background PWDR'][1] = float(value)
        else:
            raise ValueError("Invalid key in set_background:", key)
        
    def set_refinements(self, refs):
        """Sets the histogram refinement parameter 'key' to the specification 'value'.

        :param dict refs: A dictionary of the parameters to be set. See the 
                          :ref:`Histogram_parameters_table` table for a description of
                          what these dictionaries should be.

        :returns: None

        """
        do_fit_fixed_points = False
        for key, value in refs.items():
            if key == 'Limits':
                old_limits = self.data['Limits'][1]
                new_limits = value
                if isinstance(new_limits, dict):
                    if 'low' in new_limits:
                        old_limits[0] = new_limits['low']
                    if 'high' in new_limits:
                        old_limits[1] = new_limits['high']
                else:
                    old_limits[0], old_limits[1] = new_limits
            elif key == 'Sample Parameters':
                sample = self.data['Sample Parameters']
                for sparam in value:
                    if sparam not in sample:
                        raise ValueError("Unknown refinement parameter, "
                                         + str(sparam))
                    sample[sparam][1] = True
            elif key == 'Background':
                bkg, peaks = self.data['Background']

                # If True or False, just set the refine parameter
                if value in (True, False):
                    bkg[1] = value
                    return

                if 'type' in value:
                    bkg[0] = value['type']
                if 'refine' in value:
                    bkg[1] = value['refine']
                if 'no. coeffs' in value:
                    cur_coeffs = bkg[2]
                    n_coeffs = value['no. coeffs']
                    if n_coeffs > cur_coeffs:
                        for x in range(n_coeffs - cur_coeffs):
                            bkg.append(0.0)
                    else:
                        for _ in range(cur_coeffs - n_coeffs):
                            bkg.pop()
                    bkg[2] = n_coeffs
                if 'coeffs' in value:
                    bkg[3:] = value['coeffs']
                if 'FixedPoints' in value:
                    peaks['FixedPoints'] = [(float(a), float(b))
                                            for a, b in value['FixedPoints']]
                if value.get('fit fixed points', False):
                    do_fit_fixed_points = True
                if 'peaks' in value:
                    for i,flags in enumerate(value['peaks']):
                        self.ref_back_peak(i,flags)

            elif key == 'Instrument Parameters':
                instrument, secondary = self.data['Instrument Parameters']
                for iparam in value:
                    try:
                        instrument[iparam][2] = True
                    except IndexError:
                        raise ValueError("Invalid key:", iparam)
            else:
                raise ValueError("Unknown key:", key)
        # Fit fixed points after the fact - ensure they are after fixed points
        # are added
        if do_fit_fixed_points:
            # Background won't be fit if refinement flag not set
            orig = self.data['Background'][0][1]
            self.data['Background'][0][1] = True
            self.fit_fixed_points()
            # Restore the previous value
            self.data['Background'][0][1] = orig

    def clear_refinements(self, refs):
        """Clears the refinement parameter 'key' and its associated value.

        :param dict refs: A dictionary of parameters to clear.
          See the :ref:`Histogram_parameters_table` table for what can be specified. 
        """
        for key, value in refs.items():
            if key == 'Limits':
                old_limits, cur_limits = self.data['Limits']
                cur_limits[0], cur_limits[1] = old_limits
            elif key == 'Sample Parameters':
                sample = self.data['Sample Parameters']
                for sparam in value:
                    sample[sparam][1] = False
            elif key == 'Background':
                bkg, peaks = self.data['Background']

                # If True or False, just set the refine parameter
                if value in (True, False):
                    bkg[1] = False
                    return

                bkg[1] = False
                if 'FixedPoints' in value:
                    if 'FixedPoints' in peaks:
                        del peaks['FixedPoints']
                if 'peaks' in value:
                    for i in range(len(self.Background[1]['peaksList'])):
                        self.ref_back_peak(i,[])
            elif key == 'Instrument Parameters':
                instrument, secondary = self.data['Instrument Parameters']
                for iparam in value:
                    instrument[iparam][2] = False
            else:
                raise ValueError("Unknown key:", key)

    def add_peak(self,area,dspace=None,Q=None,ttheta=None):
        '''Adds a single peak to the peak list
        :param float area: peak area
        :param float dspace: peak position as d-space (A)
        :param float Q: peak position as Q (A-1)
        :param float ttheta: peak position as 2Theta (deg)

        Note: only one of the parameters: dspace, Q or ttheta may be specified.
        See :ref:`PeakRefine` for an example.
        '''
        import GSASIIlattice as G2lat
        import GSASIImath as G2mth
        if (not dspace) + (not Q) + (not ttheta) != 2:
            G2fil.G2Print('add_peak error: too many or no peak position(s) specified')
            return
        pos = ttheta
        Parms,Parms2 = self.data['Instrument Parameters']
        if Q:
            pos = G2lat.Dsp2pos(Parms,2.*np.pi/Q)
        elif dspace:
            pos = G2lat.Dsp2pos(Parms,dspace)
        peaks = self.data['Peak List']
        peaks['sigDict'] = {}        #no longer valid
        peaks['peaks'].append(G2mth.setPeakparms(Parms,Parms2,pos,area))

    def set_peakFlags(self,peaklist=None,area=None,pos=None,sig=None,gam=None):
        '''Set refinement flags for peaks
        
        :param list peaklist: a list of peaks to change flags. If None (default), changes
          are made to all peaks.
        :param bool area: Sets or clears the refinement flag for the peak area value.
          If None (the default), no change is made.
        :param bool pos: Sets or clears the refinement flag for the peak position value.
          If None (the default), no change is made.
        :param bool sig: Sets or clears the refinement flag for the peak sig (Gaussian width) value.
          If None (the default), no change is made.
        :param bool gam: Sets or clears the refinement flag for the peak sig (Lorentzian width) value.
          If None (the default), no change is made.
          
        Note that when peaks are first created the area flag is on and the other flags are
        initially off.

        Example::
        
           set_peakFlags(sig=False,gam=True)

        causes the sig refinement flag to be cleared and the gam flag to be set, in both cases for
        all peaks. The position and area flags are not changed from their previous values.
        '''
        peaks = self.data['Peak List']
        if peaklist is None:
            peaklist = range(len(peaks['peaks']))
        for i in peaklist:
            for var,j in [(area,3),(pos,1),(sig,5),(gam,7)]:
                if var is not None:
                    peaks['peaks'][i][j] = var
            
    def refine_peaks(self):
        '''Causes a refinement of peak position, background and instrument parameters
        '''
        import GSASIIpwd as G2pwd
        controls = self.proj.data.get('Controls',{})
        controls = controls.get('data',
                            {'deriv type':'analytic','min dM/M':0.001,}     #fill in defaults if needed
                            )
        peaks = self.data['Peak List']
        Parms,Parms2 = self.data['Instrument Parameters']
        background = self.data['Background']
        limits = self.data['Limits'][1]
        bxye = np.zeros(len(self.data['data'][1][1]))
        peaks['sigDict'] = G2pwd.DoPeakFit('LSQ',peaks['peaks'],background,limits,
                                           Parms,Parms2,self.data['data'][1],bxye,[],
                                           False,controls,None)[0]

    @property
    def Peaks(self):
        '''Provides a dict with the Peak List parameters
        for this histogram.

        :returns: dict with two elements where item
          'peaks' is a list of peaks where each element is 
          [pos,pos-ref,area,area-ref,sig,sig-ref,gam,gam-ref], 
          where the -ref items are refinement flags and item
          'sigDict' is a dict with possible items 'Back;#', 
          'pos#', 'int#', 'sig#', 'gam#'
        '''
        return self.data['Peak List']

    @property
    def PeakList(self):
        '''Provides a list of peaks parameters
        for this histogram.

        :returns: a list of peaks, where each peak is a list containing
          [pos,area,sig,gam] 
          (position, peak area, Gaussian width, Lorentzian width)
           
        '''
        return [i[::2] for i in self.data['Peak List']['peaks']]
    
    def Export_peaks(self,filename):
        '''Write the peaks file. The path is specified by filename
        extension.
        
        :param str filename: name of the file, optionally with a path, 
            includes an extension
        :returns: name of file that was written
        '''
        import GSASIIlattice as G2lat
        import math
        nptand = lambda x: np.tan(x*math.pi/180.)
        fil = os.path.abspath(filename)
        fp = open(filename,'w')
        Inst,Inst2 = self.data['Instrument Parameters']
        Type = Inst['Type'][0]
        if 'T' not in Type:
            import GSASIImath as G2mth
            wave = G2mth.getWave(Inst)
        else:
            wave = None
        pkdata = self.data['Peak List']
        peaks = pkdata['peaks']
        sigDict = pkdata['sigDict']
        # code taken from GSASIIdataGUI OnExportPeakList
        fp.write("#%s \n" % (self.name+' Peak List'))
        if wave:
            fp.write('#wavelength = %10.6f\n'%(wave))
        if 'T' in Type:
            fp.write('#%9s %10s %10s %12s %10s %10s %10s %10s %10s\n'%('pos','dsp','esd','int','alp','bet','sig','gam','FWHM'))                                    
        else:
            fp.write('#%9s %10s %10s %12s %10s %10s %10s\n'%('pos','dsp','esd','int','sig','gam','FWHM'))
        for ip,peak in enumerate(peaks):
            dsp = G2lat.Pos2dsp(Inst,peak[0])
            if 'T' in Type:  #TOF - more cols
                esds = {'pos':0.,'int':0.,'alp':0.,'bet':0.,'sig':0.,'gam':0.}
                for name in list(esds.keys()):
                    esds[name] = sigDict.get('%s%d'%(name,ip),0.)
                sig = np.sqrt(peak[8])
                gam = peak[10]
                esddsp = G2lat.Pos2dsp(Inst,esds['pos'])
                FWHM = G2pwd.getgamFW(gam,sig) +(peak[4]+peak[6])*np.log(2.)/(peak[4]*peak[6])     #to get delta-TOF from Gam(peak)
                fp.write("%10.2f %10.5f %10.5f %12.2f %10.3f %10.3f %10.3f %10.3f %10.3f\n" % \
                    (peak[0],dsp,esddsp,peak[2],peak[4],peak[6],peak[8],peak[10],FWHM))
            else:               #CW
                #get esds from sigDict for each peak & put in output - esds for sig & gam from UVWXY?
                esds = {'pos':0.,'int':0.,'sig':0.,'gam':0.}
                for name in list(esds.keys()):
                    esds[name] = sigDict.get('%s%d'%(name,ip),0.)
                sig = np.sqrt(peak[4]) #var -> sig
                gam = peak[6]
                esddsp = 0.5*esds['pos']*dsp/nptand(peak[0]/2.)
                FWHM = G2pwd.getgamFW(gam,sig)      #to get delta-2-theta in deg. from Gam(peak)
                fp.write("%10.4f %10.5f %10.5f %12.2f %10.5f %10.5f %10.5f \n" % \
                    (peak[0],dsp,esddsp,peak[2],np.sqrt(max(0.0001,peak[4]))/100.,peak[6]/100.,FWHM/100.)) #convert to deg
        fp.close()
        return fil

    def SaveProfile(self,filename):
        '''Writes a GSAS-II (new style) .instprm file
        '''
        data,Parms2 = self.data['Instrument Parameters']
        filename = os.path.splitext(filename)[0]+'.instprm'         # make sure extension is .instprm
        File = open(filename,'w')
        File.write("#GSAS-II instrument parameter file; do not add/delete items!\n")
        for item in data:
            File.write(item+':'+str(data[item][1])+'\n')
        File.close()
        G2fil.G2Print ('Instrument parameters saved to: '+filename)

    def LoadProfile(self,filename,bank=0):
        '''Reads a GSAS-II (new style) .instprm file and overwrites the current 
        parameters

        :param str filename: instrument parameter file name, extension ignored if not
          .instprm
        :param int bank: bank number to read, defaults to zero
        '''
        filename = os.path.splitext(filename)[0]+'.instprm'         # make sure extension is .instprm
        File = open(filename,'r')
        S = File.readline()
        newItems = []
        newVals = []
        Found = False
        while S:
            if S[0] == '#':
                if Found:
                    break
                if 'Bank' in S:
                    if bank == int(S.split(':')[0].split()[1]):
                        S = File.readline()
                        continue
                    else:
                        S = File.readline()
                        while S and '#Bank' not in S:
                            S = File.readline()
                        continue
                else:   #a non #Bank file
                    S = File.readline()
                    continue
            Found = True
            [item,val] = S[:-1].split(':')
            newItems.append(item)
            try:
                newVals.append(float(val))
            except ValueError:
                newVals.append(val)                        
            S = File.readline()                
        File.close()
        LoadG2fil()
        self.data['Instrument Parameters'][0] = G2fil.makeInstDict(newItems,newVals,len(newVals)*[False,])

    def EditSimulated(self,Tmin, Tmax, Tstep=None, Npoints=None):
        '''Change the parameters for an existing simulated powder histogram. 
        This will reset the previously computed "observed" pattern.

        :param float Tmin: Minimum 2theta or TOF (microsec) for dataset to be simulated
        :param float Tmax: Maximum 2theta or TOF (usec) for dataset to be simulated
        :param float Tstep: Step size in 2theta or TOF (usec) for dataset to be simulated       
           Default is to compute this from Npoints.
        :param int Νpoints: the number of data points to be used for computing the 
            diffraction pattern. Defaults as None, which sets this to 2500. Do not specify
            both Npoints and Tstep. Due to roundoff the actual nuber of points used may differ
            by +-1 from Npoints. Must be below 25,000.
         '''
        if not self.data['data'][0]['Dummy']:
            raise G2ScriptException("Error: histogram for G2PwdrData.EditSimulated is not simulated")            
        if Tmax < Tmin:
            Tmin,Tmax = Tmax,Tmin
        if Tstep is not None and Npoints is not None:
            raise G2ScriptException("Error: Tstep and Npoints both specified")
        elif Tstep is not None:
            Tstep = abs(Tstep)
        elif Npoints is None:
            Npoints = 2500
            
        if 'T' in self.data['Instrument Parameters'][0]['Type'][0]:
            if Tmax > 200.:
                raise G2ScriptException("Error: Tmax is too large")
            if Npoints:
                N = Npoints
                Tstep = (np.log(Tmax)-np.log(Tmin))/N
            else:
                N = (np.log(Tmax)-np.log(Tmin))/Tstep
            if N > 25000:
                raise G2ScriptException("Error: Tstep is too small. Would need "+str(N)+" points.")
            x = np.exp((np.arange(0,N))*Tstep+np.log(Tmin*1000.))
            N = len(x)
            unit = 'millisec'
            limits = [(1000*Tmin, 1000*Tmax), [1000*Tmin, 1000*Tmax]]
        else:
            if Npoints:
                N = Npoints
            else:
                N = int((Tmax-Tmin)/Tstep)+1
            if N > 25000:
                raise G2ScriptException("Error: Tstep is too small. Would need "+str(N)+" points.")
            x = np.linspace(Tmin,Tmax,N,True)
            N = len(x)
            unit = 'degrees 2theta'
            limits = [(Tmin, Tmax), [Tmin, Tmax]]
        if N < 3:
            raise G2ScriptException("Error: Range is too small or step is too large, <3 points")
        G2fil.G2Print('Simulating {} points from {} to {} {}'.format(N,Tmin,Tmax,unit))
        self.data['data'][1] = [
            np.array(x), # x-axis values
            np.zeros_like(x), # powder pattern intensities
            np.ones_like(x), # 1/sig(intensity)^2 values (weights)
            np.zeros_like(x), # calc. intensities (zero)
            np.zeros_like(x), # calc. background (zero)
            np.zeros_like(x), # obs-calc profiles
            ]
        self.data['Limits'] = limits

    def getHistEntryList(self, keyname=''):
        """Returns a dict with histogram setting values. 

        :param str keyname: an optional string. When supplied only entries
          where at least one key contains the specified string are reported.
          Case is ignored, so 'sg' will find entries where one of the keys 
          is 'SGdata', etc. 
        :returns: a set of histogram dict keys. 

        See :meth:`G2Phase.getHAPentryList` for a related example.

        .. seealso::
            :meth:`getHistEntryValue`
            :meth:`setHistEntryValue`

        """
        return [i for i in dictDive(self.data,keyname) if i[0] != ['Histograms']]

    def getHistEntryValue(self, keylist):
        """Returns the histogram control value associated with a list of keys. 
        Where the value returned is a list, it may be used as the target of 
        an assignment (as in 
        ``getHistEntryValue(...)[...] = val``) 
        to set a value inside a list.        

        :param list keylist: a list of dict keys, typically as returned by 
          :meth:`getHistEntryList`. 

        :returns: a histogram setting; may be a int, float, bool, list,...

        See :meth:`G2Phase.getHAPentryValue` for a related example.

        """
        d = self.data
        for key in keylist:
            d = d[key]
        return d

    def setHistEntryValue(self, keylist, newvalue):
        """Sets a histogram control value associated with a list of keys. 

        See :meth:`G2Phase.setHAPentryValue` for a related example.

       :param list keylist: a list of dict keys, typically as returned by 
          :meth:`getHistEntryList`. 

        :param newvalue: a new value for the hist setting. The type must be
          the same as the initial value, but if the value is a container 
          (list, tuple, np.array,...) the elements inside are not checked.

        """
        oldvalue = self.getHistEntryValue(keylist)
        if type(oldvalue) is not type(newvalue):
            raise G2ScriptException("getHistEntryValue error: types do not agree for keys {}".format(keylist))
        d = self.data
        for key in keylist:
            dlast = d
            d = d[key]
        dlast[key] = newvalue
        
class G2Phase(G2ObjectWrapper):
    """A wrapper object around a given phase.
    The object contains these class variables:

        * G2Phase.proj: contains a reference to the :class:`G2Project`
          object that contains this phase 
        * G2Phase.name: contains the name of the phase
        * G2Phase.data: contains the phases's associated data in a dict,
          as documented for the :ref:`Phase Tree items<Phase_table>`.

    Author: Jackson O'Donnell (jacksonhodonnell .at. gmail.com)
    """
    def __init__(self, data, name, proj):
        self.data = data
        self.name = name
        self.proj = proj

    @staticmethod
    def is_valid_refinement_key(key):
        valid_keys = ["Cell", "Atoms", "LeBail"]
        return key in valid_keys

    @staticmethod
    def is_valid_HAP_refinement_key(key):
        valid_keys = ["Babinet", "Extinction", "HStrain", "Mustrain",
                      "Pref.Ori.", "Show", "Size", "Use", "Scale"]
        return key in valid_keys

    def atom(self, atomlabel):
        """Returns the atom specified by atomlabel, or None if it does not
        exist.

        :param str atomlabel: The name of the atom (e.g. "O2")
        :returns: A :class:`G2AtomRecord` object
            representing the atom.
        """
        # Consult GSASIIobj.py for the meaning of this
        cx, ct, cs, cia = self.data['General']['AtomPtrs']
        ptrs = [cx, ct, cs, cia]
        atoms = self.data['Atoms']
        for atom in atoms:
            if atom[ct-1] == atomlabel:
                return G2AtomRecord(atom, ptrs, self.proj)

    def atoms(self):
        """Returns a list of atoms present in the current phase.

        :returns: A list of :class:`G2AtomRecord` objects.

        .. seealso::
            :meth:`~G2Phase.atom`
            :class:`G2AtomRecord`
        """
        ptrs = self.data['General']['AtomPtrs']
        return [G2AtomRecord(atom, ptrs, self.proj) for atom in self.data['Atoms']]

    def histograms(self):
        '''Returns a list of histogram names associated with the current phase ordered 
        as they appear in the tree (see :meth:`G2Project.histograms`).
        '''
        return [i.name for i in self.proj.histograms() if i.name in self.data.get('Histograms', {})]
    
    @property
    def composition(self):
        '''Provides a dict where keys are atom types and values are the number of 
        atoms of that type in cell (such as {'H': 2.0, 'O': 1.0})
        '''
        out = {}
        for a in self.atoms():
            typ = a.element
            if typ in out:
                out[typ] += a.mult*a.occupancy
            else:
                out[typ] = a.mult*a.occupancy
        return out

    def mu(self,wave):
        '''Provides mu values for a phase at the supplied wavelength in A.
        Uses GSASIImath.XScattDen which seems to be off by an order of 
        magnitude, which has been corrected here.
        '''
        import GSASIImath as G2mth
        vol = self.data['General']['Cell'][7]
        out = {}
        for typ in self.data['General']['NoAtoms']:
            if typ in out:
                out[typ]['Num'] += self.data['General']['NoAtoms'][typ]
            else:
                out[typ] = {}
                out[typ]['Num'] = self.data['General']['NoAtoms'][typ]
                out[typ]['Z'] = 0 # wrong but not needed
        return 10*G2mth.XScattDen(out,vol,wave)[1]
    
    @property
    def density(self):
        '''Provides a scalar with the density of the phase. In case of a 
        powder this assumes a 100% packing fraction.
        '''
        import GSASIImath as G2mth
        density,mattCoeff = G2mth.getDensity(self.data['General'])
        return density
    
    @property
    def ranId(self):
        return self.data['ranId']

    @property
    def id(self):
        return self.data['pId']

    @id.setter
    def id(self, val):
        self.data['pId'] = val

    def get_cell(self):
        """Returns a dictionary of the cell parameters, with keys:
            'length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma', 'volume'

        :returns: a dict

        .. seealso::
           :meth:`~G2Phase.get_cell_and_esd`

        """
        cell = self.data['General']['Cell']
        return {'length_a': cell[1], 'length_b': cell[2], 'length_c': cell[3],
                'angle_alpha': cell[4], 'angle_beta': cell[5], 'angle_gamma': cell[6],
                'volume': cell[7]}

    def get_cell_and_esd(self):
        """
        Returns a pair of dictionaries, the first representing the unit cell, the second
        representing the estimated standard deviations of the unit cell.

        :returns: a tuple of two dictionaries

        .. seealso::
           :meth:`~G2Phase.get_cell`

        """
        # translated from GSASIIstrIO.ExportBaseclass.GetCell
        import GSASIIlattice as G2lat
        import GSASIImapvars as G2mv
        try:
            pfx = str(self.id) + '::'
            sgdata = self['General']['SGData']
            covDict = self.proj['Covariance']['data']

            parmDict = dict(zip(covDict.get('varyList',[]),
                                covDict.get('variables',[])))
            sigDict = dict(zip(covDict.get('varyList',[]),
                               covDict.get('sig',[])))

            if covDict.get('covMatrix') is not None:
                sigDict.update(G2mv.ComputeDepESD(covDict['covMatrix'],
                                                  covDict['varyList'],
                                                  parmDict))

            A, sigA = G2strIO.cellFill(pfx, sgdata, parmDict, sigDict)
            cellSig = G2strIO.getCellEsd(pfx, sgdata, A, self.proj['Covariance']['data'])
            cellList = G2lat.A2cell(A) + (G2lat.calc_V(A),)
            cellDict, cellSigDict = {}, {}
            for i, key in enumerate(['length_a', 'length_b', 'length_c',
                                     'angle_alpha', 'angle_beta', 'angle_gamma',
                                     'volume']):
                cellDict[key] = cellList[i]
                cellSigDict[key] = cellSig[i]
            return cellDict, cellSigDict
        except KeyError:
            cell = self.get_cell()
            return cell, {key: 0.0 for key in cell}

    def export_CIF(self, outputname, quickmode=True):
        """Write this phase to a .cif file named outputname

        :param str outputname: The name of the .cif file to write to
        :param bool quickmode: Currently ignored. Carryover from exports.G2export_CIF"""
        # This code is all taken from exports/G2export_CIF.py
        # Functions copied have the same names
        import GSASIImath as G2mth
        import GSASIImapvars as G2mv
        from exports import G2export_CIF as cif

#        CIFdate = dt.datetime.strftime(dt.datetime.now(),"%Y-%m-%dT%H:%M")
        CIFname = os.path.splitext(self.proj.filename)[0]
        CIFname = os.path.split(CIFname)[1]
        CIFname = ''.join([c if ord(c) < 128 else ''
                           for c in CIFname.replace(' ', '_')])
        # try:
        #     author = self.proj['Controls']['data'].get('Author','').strip()
        # except KeyError:
        #     pass
        # oneblock = True

        covDict = self.proj['Covariance']['data']
        parmDict = dict(zip(covDict.get('varyList',[]),
                            covDict.get('variables',[])))
        sigDict = dict(zip(covDict.get('varyList',[]),
                           covDict.get('sig',[])))

        if covDict.get('covMatrix') is not None:
            sigDict.update(G2mv.ComputeDepESD(covDict['covMatrix'],
                                              covDict['varyList'],
                                              parmDict))

        with open(outputname, 'w') as fp:
            fp.write(' \n' + 70*'#' + '\n')
            cif.WriteCIFitem(fp, 'data_' + CIFname)
            # from exports.G2export_CIF.WritePhaseInfo
            cif.WriteCIFitem(fp, '\n# phase info for '+str(self.name) + ' follows')
            cif.WriteCIFitem(fp, '_pd_phase_name', self.name)
            # TODO get esds
            cellDict = self.get_cell()
            defsigL = 3*[-0.00001] + 3*[-0.001] + [-0.01] # significance to use when no sigma
            names = ['length_a','length_b','length_c',
                     'angle_alpha','angle_beta ','angle_gamma',
                     'volume']
            for key, val in cellDict.items():
                cif.WriteCIFitem(fp, '_cell_' + key, G2mth.ValEsd(val))

            cif.WriteCIFitem(fp, '_symmetry_cell_setting',
                         self.data['General']['SGData']['SGSys'])

            spacegroup = self.data['General']['SGData']['SpGrp'].strip()
            # regularize capitalization and remove trailing H/R
            spacegroup = spacegroup[0].upper() + spacegroup[1:].lower().rstrip('rh ')
            cif.WriteCIFitem(fp, '_symmetry_space_group_name_H-M', spacegroup)

            # generate symmetry operations including centering and center of symmetry
            SymOpList, offsetList, symOpList, G2oprList, G2opcodes = G2spc.AllOps(
                self.data['General']['SGData'])
            cif.WriteCIFitem(fp, 'loop_\n    _space_group_symop_id\n    _space_group_symop_operation_xyz')
            for i, op in enumerate(SymOpList,start=1):
                cif.WriteCIFitem(fp, '   {:3d}  {:}'.format(i,op.lower()))

            # TODO skipped histograms, exports/G2export_CIF.py:880

            # report atom params
            if self.data['General']['Type'] in ['nuclear','macromolecular']:        #this needs macromolecular variant, etc!
                cif.WriteAtomsNuclear(fp, self.data, self.name, parmDict, sigDict, [])
                # self._WriteAtomsNuclear(fp, parmDict, sigDict)
            else:
                raise G2ScriptException("no export for "+str(self.data['General']['Type'])+" coordinates implemented")
            # report cell contents
            cif.WriteComposition(fp, self.data, self.name, parmDict)
            if not quickmode and self.data['General']['Type'] == 'nuclear':      # report distances and angles
                # WriteDistances(fp,self.name,SymOpList,offsetList,symOpList,G2oprList)
                raise NotImplementedError("only quickmode currently supported")
            if 'Map' in self.data['General'] and 'minmax' in self.data['General']['Map']:
                cif.WriteCIFitem(fp,'\n# Difference density results')
                MinMax = self.data['General']['Map']['minmax']
                cif.WriteCIFitem(fp,'_refine_diff_density_max',G2mth.ValEsd(MinMax[0],-0.009))
                cif.WriteCIFitem(fp,'_refine_diff_density_min',G2mth.ValEsd(MinMax[1],-0.009))


    def set_refinements(self, refs):
        """Sets the phase refinement parameter 'key' to the specification 'value'

        :param dict refs: A dictionary of the parameters to be set. See the
                          :ref:`Phase_parameters_table` table for a description of
                          this dictionary.

        :returns: None"""
        for key, value in refs.items():
            if key == "Cell":
                self.data['General']['Cell'][0] = value

            elif key == "Atoms":
                for atomlabel, atomrefinement in value.items():
                    if atomlabel == 'all':
                        for atom in self.atoms():
                            atom.refinement_flags = atomrefinement
                    else:
                        atom = self.atom(atomlabel)
                        if atom is None:
                            raise ValueError("No such atom: " + atomlabel)
                        atom.refinement_flags = atomrefinement

            elif key == "LeBail":
                hists = self.data['Histograms']
                for hname, hoptions in hists.items():
                    if 'LeBail' not in hoptions:
                        hoptions['newLeBail'] = bool(True)
                    hoptions['LeBail'] = bool(value)
            else:
                raise ValueError("Unknown key:", key)

    def clear_refinements(self, refs):
        """Clears a given set of parameters.

        :param dict refs: The parameters to clear.
          See the :ref:`Phase_parameters_table` table for what can be specified. 
        """
        for key, value in refs.items():
            if key == "Cell":
                self.data['General']['Cell'][0] = False
            elif key == "Atoms":
                cx, ct, cs, cia = self.data['General']['AtomPtrs']

                for atomlabel in value:
                    atom = self.atom(atomlabel)
                    # Set refinement to none
                    atom.refinement_flags = ' '
            elif key == "LeBail":
                hists = self.data['Histograms']
                for hname, hoptions in hists.items():
                    if 'LeBail' not in hoptions:
                        hoptions['newLeBail'] = True
                    hoptions['LeBail'] = False
            else:
                raise ValueError("Unknown key:", key)

    def set_HAP_refinements(self, refs, histograms='all'):
        """Sets the given HAP refinement parameters between the current phase and
        the specified histograms.

        :param dict refs: A dictionary of the parameters to be set. See
                          the :ref:`HAP_parameters_table` table for a description of this
                          dictionary.
        :param histograms: Either 'all' (default) or a list of the histograms by index, name
            or object. The index number is relative to all histograms in the tree, not to 
            those in the phase. 
            Histograms not associated with the current phase will be ignored. 
            whose HAP parameters will be set with this phase. Histogram and phase
            must already be associated. 
        :returns: None
        """
        if not self.data.get('Histograms',[]):
            G2fil.G2Print("Error likely: Phase {} has no linked histograms".format(self.name))
            return
            
        if histograms == 'all':
            histograms = self.data['Histograms'].keys()
        else:
            histograms = [self._decodeHist(h) for h in histograms 
                          if self._decodeHist(h) in self.data['Histograms']]    
        if not histograms:
            G2fil.G2Print("Warning: Skipping HAP set for phase {}, no selected histograms".format(self.name))
            return
        for key, val in refs.items():
            if key == 'Babinet':
                try:
                    sets = list(val)
                except ValueError:
                    sets = ['BabA', 'BabU']
                for param in sets:
                    if param not in ['BabA', 'BabU']:
                        raise ValueError("Not sure what to do with" + param)
                    for h in histograms:
                        self.data['Histograms'][h]['Babinet'][param][1] = True
            elif key == 'Extinction':
                for h in histograms:
                    self.data['Histograms'][h]['Extinction'][1] = bool(val)
            elif key == 'HStrain':
                if isinstance(val,list) or isinstance(val,tuple):
                    for h in histograms:
                        if len(self.data['Histograms'][h]['HStrain'][1]) != len(val):
                            raise Exception('Need {} HStrain terms for phase {} hist {}'
                                .format(len(self.data['Histograms'][h]['HStrain'][1]),self.name,h))
                        for i,v in enumerate(val):
                            self.data['Histograms'][h]['HStrain'][1][i] = bool(v)
                else:
                    for h in histograms:
                        self.data['Histograms'][h]['HStrain'][1] = [bool(val) for p in self.data['Histograms'][h]['HStrain'][1]]
            elif key == 'Mustrain':
                for h in histograms:
                    mustrain = self.data['Histograms'][h]['Mustrain']
                    newType = None
                    direction = None
                    if isinstance(val, strtypes):
                        if val in ['isotropic', 'uniaxial', 'generalized']:
                            newType = val
                        else:
                            raise ValueError("Not a Mustrain type: " + val)
                    elif isinstance(val, dict):
                        newType = val.get('type', None)
                        direction = val.get('direction', None)

                    if newType:
                        mustrain[0] = newType
                        if newType == 'isotropic':
                            mustrain[2][0] = True == val.get('refine',False)
                            mustrain[5] = [False for p in mustrain[4]]
                        elif newType == 'uniaxial':
                            if 'refine' in val:
                                mustrain[2][0] = False
                                types = val['refine']
                                if isinstance(types, strtypes):
                                    types = [types]
                                elif isinstance(types, bool):
                                    mustrain[2][1] = types
                                    mustrain[2][2] = types
                                    types = []
                                else:
                                    raise ValueError("Not sure what to do with: "
                                                     + str(types))
                            else:
                                types = []

                            for unitype in types:
                                if unitype == 'equatorial':
                                    mustrain[2][0] = True
                                elif unitype == 'axial':
                                    mustrain[2][1] = True
                                else:
                                    msg = 'Invalid uniaxial mustrain type'
                                    raise ValueError(msg + ': ' + unitype)
                        else:  # newtype == 'generalized'
                            mustrain[2] = [False for p in mustrain[1]]
                            if 'refine' in val:
                                mustrain[5] = [True == val['refine']]*len(mustrain[5])

                    if direction:
                        if len(direction) != 3:
                            raise ValueError("Expected hkl, found", direction)
                        direction = [int(n) for n in direction]
                        mustrain[3] = direction
            elif key == 'Size':
                newSize = None
                if 'value' in val:
                    newSize = float(val['value'])
                for h in histograms:
                    size = self.data['Histograms'][h]['Size']
                    newType = None
                    direction = None
                    if isinstance(val, strtypes):
                        if val in ['isotropic', 'uniaxial', 'ellipsoidal']:
                            newType = val
                        else:
                            raise ValueError("Not a valid Size type: " + val)
                    elif isinstance(val, dict):
                        newType = val.get('type', None)
                        direction = val.get('direction', None)

                    if newType:
                        size[0] = newType
                        refine = bool(val.get('refine'))
                        if newType == 'isotropic' and refine is not None:
                            size[2][0] = bool(refine)
                            if newSize: size[1][0] = newSize
                        elif newType == 'uniaxial' and refine is not None:
                            size[2][1] = bool(refine)
                            size[2][2] = bool(refine)
                            if newSize: size[1][1] = size[1][2] =newSize
                        elif newType == 'ellipsoidal' and refine is not None:
                            size[5] = [bool(refine) for p in size[5]]
                            if newSize: size[4] = [newSize for p in size[4]]

                    if direction:
                        if len(direction) != 3:
                            raise ValueError("Expected hkl, found", direction)
                        direction = [int(n) for n in direction]
                        size[3] = direction
            elif key == 'Pref.Ori.':
                for h in histograms:
                    self.data['Histograms'][h]['Pref.Ori.'][2] = bool(val)
            elif key == 'Show':
                for h in histograms:
                    self.data['Histograms'][h]['Show'] = bool(val)
            elif key == 'Use':
                for h in histograms:
                    self.data['Histograms'][h]['Use'] = bool(val)
            elif key == 'Scale':
                for h in histograms:
                    self.data['Histograms'][h]['Scale'][1] = bool(val)
            else:
                G2fil.G2Print(u'Warning: Unknown HAP key: '+key)

    def clear_HAP_refinements(self, refs, histograms='all'):
        """Clears the given HAP refinement parameters between this phase and
        the given histograms.

        :param dict refs: A dictionary of the parameters to be cleared.
            See the the :ref:`HAP_parameters_table` table for what can be specified. 
        :param histograms: Either 'all' (default) or a list of the histograms by index, name
            or object. 
            The index number is relative to all histograms in the tree, not to 
            those in the phase.
            Histograms not associated with the current phase will be ignored. 
            whose HAP parameters will be set with this phase. Histogram and phase
            must already be associated
        :returns: None
        """
        if histograms == 'all':
            histograms = self.data['Histograms'].keys()
        else:
            histograms = [self._decodeHist(h) for h in histograms 
                          if self._decodeHist(h) in self.data['Histograms']]    

        for key, val in refs.items():
            for h in histograms:
                if key == 'Babinet':
                    try:
                        sets = list(val)
                    except ValueError:
                        sets = ['BabA', 'BabU']
                    for param in sets:
                        if param not in ['BabA', 'BabU']:
                            raise ValueError("Not sure what to do with" + param)
                        for h in histograms:
                            self.data['Histograms'][h]['Babinet'][param][1] = False
                elif key == 'Extinction':
                    for h in histograms:
                        self.data['Histograms'][h]['Extinction'][1] = False
                elif key == 'HStrain':
                    for h in histograms:
                        self.data['Histograms'][h]['HStrain'][1] = [False for p in self.data['Histograms'][h]['HStrain'][1]]
                elif key == 'Mustrain':
                    for h in histograms:
                        mustrain = self.data['Histograms'][h]['Mustrain']
                        mustrain[2] = [False for p in mustrain[2]]
                        mustrain[5] = [False for p in mustrain[4]]
                elif key == 'Pref.Ori.':
                    for h in histograms:
                        self.data['Histograms'][h]['Pref.Ori.'][2] = False
                elif key == 'Show':
                    for h in histograms:
                        self.data['Histograms'][h]['Show'] = False
                elif key == 'Size':
                    for h in histograms:
                        size = self.data['Histograms'][h]['Size']
                        size[2] = [False for p in size[2]]
                        size[5] = [False for p in size[5]]
                elif key == 'Use':
                    for h in histograms:
                        self.data['Histograms'][h]['Use'] = False
                elif key == 'Scale':
                    for h in histograms:
                        self.data['Histograms'][h]['Scale'][1] = False
                else:
                    G2fil.G2Print(u'Warning: Unknown HAP key: '+key)

    def _decodeHist(self,hist):
        '''Convert a histogram reference to a histogram name string
        '''
        if isinstance(hist, G2PwdrData):
            return hist.name
        elif hist in self.data['Histograms']:
            return hist
        elif type(hist) is int:
            return self.proj.histograms()[hist].name
        else:
            raise G2ScriptException("Invalid histogram reference: "+str(hist))
        
    def getHAPvalues(self, histname):
        """Returns a dict with HAP values for the selected histogram

        :param histogram: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram.
            The index number is relative to all histograms in the tree, not to 
            those in the phase.
        :returns: HAP value dict
        """
        return self.data['Histograms'][self._decodeHist(histname)]

    def copyHAPvalues(self, sourcehist, targethistlist='all', skip=[], use=None):
        """Copies HAP parameters for one histogram to a list of other histograms.
        Use skip or use to select specific entries to be copied or not used.

        :param sourcehist: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram to copy 
            parameters from.
            The index number is relative to all histograms in the tree, not to 
            those in the phase.
        :param list targethistlist: a list of histograms where each item in the
            list can be a histogram object (:class:`G2PwdrData`), 
            a histogram name or the index number of the histogram.
            If the string 'all' (default), then all histograms in the phase 
            are used.
        :param list skip: items in the HAP dict that should not be 
            copied. The default is an empty list, which causes all items
            to be copied. To see a list of items in the dict, use
            :meth:`getHAPvalues` or use an invalid item, such as '?'.
        :param list use: specifies the items in the HAP dict should be 
            copied. The default is None, which causes all items
            to be copied. 

        examples::

            ph0.copyHAPvalues(0,[1,2,3])
            ph0.copyHAPvalues(0,use=['HStrain','Size'])

        The first example copies all HAP parameters from the first histogram to 
        the second, third and fourth histograms (as listed in the project tree). 
        The second example copies only the 'HStrain' (Dij parameters and 
        refinement flags) and the 'Size' (crystallite size settings, parameters
        and refinement flags) from the first histogram to all histograms.
        """
        sourcehist = self._decodeHist(sourcehist)
        if targethistlist == 'all':
            targethistlist = self.histograms()
        
        copydict = copy.deepcopy(self.data['Histograms'][sourcehist])
        for item in skip:
            if item in list(copydict.keys()):
                del copydict[item]
            else:
                G2fil.G2Print('items in HAP dict are: {}'.format(
                    list(self.data['Histograms'][sourcehist])))
                raise Exception('HAP skip list entry {} invalid'.format(item))
        if use:
            for item in list(copydict.keys()):
                if item not in use:
                    del copydict[item]

        G2fil.G2Print('Copying item(s) {} from histogram {}'.format(list(copydict.keys()),sourcehist))
        G2fil.G2Print(' to histogram(s) {}'.format([self._decodeHist(h) for h in targethistlist]))
        for h in targethistlist:
            h = self._decodeHist(h)
            if h not in self.data['Histograms']:
                G2fil.G2Print('Unexpected Warning: histogram {} not in phase {}'.format(h,self.name))
                continue
            self.data['Histograms'][h].update(copy.deepcopy(copydict))
            
    def setSampleProfile(self, histname, parmType, mode, val1, val2=None, axis=None, LGmix=None):
        """Sets sample broadening parameters for a histogram associated with the 
        current phase. This currently supports isotropic and uniaxial broadening 
        modes only. 

        :param histogram: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram.
            The index number is relative to all histograms in the tree, not to 
            those in the phase.
        :param str parmType: should be 'size' or 'microstrain' (can be abbreviated to 's' or 'm')
        :param str mode: should be 'isotropic' or 'uniaxial' (can be abbreviated to 'i' or 'u')
        :param float val1: value for isotropic size (in :math:`\\mu m`) or  
           microstrain (unitless, :math:`\\Delta Q/Q \\times 10^6`) or the equatorial value in the uniaxial case
        :param float val2: value for axial size (in :math:`\\mu m`) or  
           axial microstrain (unitless, :math:`\\Delta Q/Q \\times 10^6`) 
           in uniaxial case; not used for isotropic 
        :param list axis: tuple or list with three values indicating the preferred direction
          for uniaxial broadening; not used for isotropic 
        :param float LGmix: value for broadening type (1=Lorentzian, 0=Gaussian or a value
          between 0 and 1. Default value (None) is ignored.

        Examples::

            phase0.setSampleProfile(0,'size','iso',1.2)
            phase0.setSampleProfile(0,'micro','isotropic',1234)
            phase0.setSampleProfile(0,'m','u',1234,4567,[1,1,1],.5) 
            phase0.setSampleProfile(0,'s','uni',1.2,2.3,[0,0,1])
        """
        if parmType.lower().startswith('s'):
            key = 'Size'
        elif parmType.lower().startswith('m'):
            key = 'Mustrain'
        else:
            G2fil.G2Print('setSampleProfile Error: value for parmType of {} is not size or microstrain'.
                              format(parmType))
            raise Exception('Invalid parameter in setSampleProfile')
        if mode.lower().startswith('i'):
            iso = True
        elif mode.lower().startswith('u'):
            iso = False
            if val2 is None:
                G2fil.G2Print('setSampleProfile Error: value for val2 is required with mode of uniaxial')
                raise Exception('Invalid val2 parameter in setSampleProfile')
            if axis is None:
                G2fil.G2Print('setSampleProfile Error: value for axis is required with mode of uniaxial')
                raise Exception('Invalid axis parameter in setSampleProfile')
        else:
            G2fil.G2Print('setSampleProfile Error: value for mode of {} is not isotropic or uniaxial'.
                              format(mode))
            raise Exception('Invalid parameter in setSampleProfile')
        
        d = self.data['Histograms'][self._decodeHist(histname)][key]
        if iso:
            d[0] = 'isotropic'
            d[1][0] = float(val1)
            if LGmix is not None: d[1][2] = float(LGmix)
        else:
            d[3] = [int(axis[0]),int(axis[1]),int(axis[2])]            
            d[0] = 'uniaxial'
            d[1][0] = float(val1)
            d[1][1] = float(val2)
            if LGmix is not None: d[1][2] = float(LGmix)            

    def setHAPvalues(self, HAPdict, targethistlist='all', skip=[], use=None):
        """Copies HAP parameters for one histogram to a list of other histograms.
        Use skip or use to select specific entries to be copied or not used.
        Note that ``HStrain`` and sometimes ``Mustrain`` values can be specific to 
        a Laue class and should be copied with care between phases of different 
        symmetry. A "sanity check" on the number of Dij terms is made if ``HStrain``
        values are copied.

        :param dict HAPdict: is a dict returned by :meth:`getHAPvalues` containing 
            HAP parameters.
        :param list targethistlist: a list of histograms where each item in the
            list can be a histogram object (:class:`G2PwdrData`), 
            a histogram name or the index number of the histogram.
            The index number is relative to all histograms in the tree, not to 
            those in the phase.
            If the string 'all' (default), then all histograms in the phase 
            are used.
        :param list skip: items in the HAP dict that should not be 
            copied. The default is an empty list, which causes all items
            to be copied. To see a list of items in the dict, use
            :meth:`getHAPvalues` or use an invalid item, such as '?'.
        :param list use: specifies the items in the HAP dict should be 
            copied. The default is None, which causes all items
            to be copied. 

        example::

            HAPdict = ph0.getHAPvalues(0)
            ph1.setHAPvalues(HAPdict,use=['HStrain','Size'])

        This copies the Dij (hydrostatic strain) HAP parameters and the 
        crystallite size broadening terms from the first histogram in 
        phase ``ph0`` to all histograms in phase ``ph1``. 
        """
        if targethistlist == 'all':
            targethistlist = self.histograms()
        copydict = copy.deepcopy(HAPdict)
        for item in skip:
            if item in list(copydict.keys()):
                del copydict[item]
            # else:
            #     G2fil.G2Print('items in HAP dict are: {}'.format(
            #         list(self.data['Histograms'][sourcehist])))
            #     raise Exception('HAP skip list entry {} invalid'.format(item))
        if use:
            for item in list(copydict.keys()):
                if item not in use:
                    del copydict[item]

        first = True
        for h in targethistlist:
            h = self._decodeHist(h)
            if h not in self.data['Histograms']:
                G2fil.G2Print('Warning: histogram {} not in phase {}'.format(h,self.name))
                continue
            if first:
                first = False
                if 'HStrain' in self.data['Histograms'][h] and 'HStrain' in copydict:
                    if len(copydict['HStrain'][0]) != len(self.data['Histograms'][h]['HStrain'][0]):
                        G2fil.G2Print('Error: HStrain has differing numbers of terms. Input: {}, phase {}: {}'.
                                  format(len(copydict['HStrain'][0]),
                                        self.name,len(self.data['Histograms'][h]['HStrain'][0])))
                        raise Exception('HStrain has differing numbers of terms.')
            self.data['Histograms'][h].update(copy.deepcopy(copydict))            
        G2fil.G2Print('Copied item(s) {} from dict'.format(list(copydict.keys())))
        G2fil.G2Print(' to histogram(s) {}'.format([self._decodeHist(h) for h in targethistlist]))

    def getPhaseEntryList(self, keyname=''):
        """Returns a dict with control values. 

        :param str keyname: an optional string. When supplied only entries
          where at least one key contains the specified string are reported.
          Case is ignored, so 'sg' will find entries where one of the keys 
          is 'SGdata', etc. 
        :returns: a set of phase dict keys. Note that HAP items, while 
          technically part of the phase entries, are not included. 

        See :meth:`getHAPentryList` for a related example.

        .. seealso::
            :meth:`getPhaseEntryValue`
            :meth:`setPhaseEntryValue`

        """
        return [i for i in dictDive(self.data,keyname) if i[0] != ['Histograms']]

    def getPhaseEntryValue(self, keylist):
        """Returns the value associated with a list of keys. 
        Where the value returned is a list, it may be used as the target of 
        an assignment (as in 
        ``getPhaseEntryValue(...)[...] = val``) 
        to set a value inside a list.        

        :param list keylist: a list of dict keys, typically as returned by 
          :meth:`getPhaseEntryList`. 

        :returns: a phase setting; may be a int, float, bool, list,...

        See :meth:`getHAPentryValue` for a related example.

        """
        d = self.data
        for key in keylist:
            d = d[key]
        return d

    def setPhaseEntryValue(self, keylist, newvalue):
        """Sets a phase control value associated with a list of keys. 

        :param list keylist: a list of dict keys, typically as returned by 
          :meth:`getPhaseEntryList`. 

        :param newvalue: a new value for the phase setting. The type must be
          the same as the initial value, but if the value is a container 
          (list, tuple, np.array,...) the elements inside are not checked.

        See :meth:`setHAPentryValue` for a related example.

        """
        oldvalue = self.getPhaseEntryValue(keylist)
        if type(oldvalue) is not type(newvalue):
            raise G2ScriptException("getPhaseEntryValue error: types do not agree for keys {}".format(keylist))
        d = self.data
        for key in keylist:
            dlast = d
            d = d[key]
        dlast[key] = newvalue
        
    def getHAPentryList(self, histname=None, keyname=''):
        """Returns a dict with HAP values. Optionally a histogram 
        may be selected.

        :param histname: is a histogram object (:class:`G2PwdrData`) or
            a histogram name or the index number of the histogram.
            The index number is relative to all histograms in the tree, not to 
            those in the phase. If no histogram is specified, all histograms 
            are selected. 
        :param str keyname: an optional string. When supplied only entries
          where at least one key contains the specified string are reported.
          Case is ignored, so 'sg' will find entries where one of the keys 
          is 'SGdata', etc. 
        :returns: a set of HAP dict keys. 

        Example: 

        >>> p.getHAPentryList(0,'Scale')
        [(['PWDR test Bank 1', 'Scale'], list, [1.0, False])]

        .. seealso::
            :meth:`getHAPentryValue`
            :meth:`setHAPentryValue`

        """
        if histname:
            h = [self._decodeHist(histname)]
        else:
            h = []
        return dictDive(self.data['Histograms'],keyname,h)

    def getHAPentryValue(self, keylist):
        """Returns the HAP value associated with a list of keys. Where the
        value returned is a list, it may be used as the target of 
        an assignment (as in 
        ``getHAPentryValue(...)[...] = val``) 
        to set a value inside a list.

        :param list keylist: a list of dict keys, typically as returned by 
          :meth:`getHAPentryList`. Note the first entry is a histogram name. 
          Example: ``['PWDR hist1.fxye Bank 1', 'Scale']``

        :returns: HAP value

        Example: 

        >>> sclEnt = p.getHAPentryList(0,'Scale')[0]                                    >>> sclEnt
        [(['PWDR test Bank 1', 'Scale'], list, [1.0, False])]
        >>> p.getHAPentryValue(sclEnt[0])
        [1.0, False]
        >>> p.getHAPentryValue(sclEnt[0])[1] = True
        >>> p.getHAPentryValue(sclEnt[0])
        [1.0, True]

        """
        d = self.data['Histograms']
        for key in keylist:
            d = d[key]
        return d

    def setHAPentryValue(self, keylist, newvalue):
        """Sets an HAP value associated with a list of keys. 

        :param list keylist: a list of dict keys, typically as returned by 
          :meth:`getHAPentryList`. Note the first entry is a histogram name. 
          Example: ``['PWDR hist1.fxye Bank 1', 'Scale']``

        :param newvalue: a new value for the HAP setting. The type must be
          the same as the initial value, but if the value is a container 
          (list, tuple, np.array,...) the elements inside are not checked.

        Example: 

        >>> sclEnt = p.getHAPentryList(0,'Scale')[0]
        >>> p.getHAPentryValue(sclEnt[0])
        [1.0, False]
        >>> p.setHAPentryValue(sclEnt[0], (1, True))
        GSASIIscriptable.G2ScriptException: setHAPentryValue error: types do not agree for keys ['PWDR test.fxye Bank 1', 'Scale']
        >>> p.setHAPentryValue(sclEnt[0], [1, True])
        >>> p.getHAPentryValue(sclEnt[0])
        [1, True]

        """
        oldvalue = self.getHAPentryValue(keylist)
        if type(oldvalue) is not type(newvalue):
            raise G2ScriptException("setHAPentryValue error: types do not agree for keys {}".format(keylist))
        d = self.data['Histograms']
        for key in keylist:
            dlast = d
            d = d[key]
        dlast[key] = newvalue
    
class G2SeqRefRes(G2ObjectWrapper):
    '''Wrapper for a Sequential Refinement Results tree entry, containing the 
    results for a refinement

    As an example:: 

        from __future__ import division, print_function
        import os,sys
        sys.path.insert(0,'/Users/toby/software/G2/GSASII')
        PathWrap = lambda fil: os.path.join('/Users/toby/Scratch/SeqTut2019Mar',fil)
        import GSASIIscriptable as G2sc
        gpx = G2sc.G2Project(PathWrap('scr4.gpx'))
        seq = gpx.seqref()
        lbl = ('a','b','c','alpha','beta','gamma','Volume')
        for j,h in enumerate(seq.histograms()):
            cell,cellU,uniq = seq.get_cell_and_esd(1,h)
            print(h)
            print([cell[i] for i in list(uniq)+[6]])
            print([cellU[i] for i in list(uniq)+[6]])
            print('')
        print('printed',[lbl[i] for i in list(uniq)+[6]])

    .. seealso::
        :meth:`G2Project.seqref`
    '''
    def __init__(self, data, proj):
        self.data = data
        self.proj = proj
        self.newCellDict = {} # dict with recp. cell tensor & Dij name
        # newAtomDict = {} # dict with atom positions; relative & absolute
        for name in self.data['histNames']:
            self.newCellDict.update(self.data[name].get('newCellDict',{}))
            #newAtomDict.update(self.data[name].get('newAtomDict',{}) # dict with atom positions; relative & absolute
        #ESDlookup = {self.newCellDict[item][0]:item for item in self.newCellDict}
        #Dlookup = {item:self.newCellDict[item][0] for item in self.newCellDict}
        # Possible error: the next might need to be data[histNames[0]]['varyList']
        #atomLookup = {newAtomDict[item][0]:item for item in newAtomDict if item in self.data['varyList']}
        #Dlookup.update({atomLookup[parm]:parm for parm in atomLookup}
        #ESDlookup.update({parm:atomLookup[parm] for parm in atomLookup})

#    @property
    def histograms(self):
        '''returns a list of histograms in the squential fit
        '''
        return self.data['histNames']

    def get_cell_and_esd(self,phase,hist):
        '''Returns a vector of cell lengths and esd values

        :param phase: A phase, which may be specified as a phase object
          (see :class:`G2Phase`), the phase name (str) or the index number (int)
          of the phase in the project, numbered starting from 0.
        :param hist: Specify a histogram or using the histogram name (str) 
          or the index number (int) of the histogram in the sequential 
          refinement (not the project), numbered as in in the project tree
          starting from 0.
        :returns: cell,cellESD,uniqCellIndx where cell (list) 
          with the unit cell parameters (a,b,c,alpha,beta,gamma,Volume);
          cellESD are the standard uncertainties on the 7 unit cell 
          parameters; and uniqCellIndx is a tuple with indicies for the 
          unique (non-symmetry determined) unit parameters (e.g. 
          [0,2] for a,c in a tetragonal cell)
        '''
        def striphist(var,insChar=''):
            'strip a histogram number from a var name'
            sv = var.split(':')
            if len(sv) <= 1: return var
            if sv[1]:
                sv[1] = insChar
            return ':'.join(sv)
        import GSASIIlattice as G2lat
        import GSASIIstrIO as G2stIO
        
        uniqCellLookup = [
        [['m3','m3m'],(0,)],
        [['3R','3mR'],(0,3)],
        [['3','3m1','31m','6/m','6/mmm','4/m','4/mmm'],(0,2)],
        [['mmm'],(0,1,2)],
        [['2/m'+'a'],(0,1,2,3)],
        [['2/m'+'b'],(0,1,2,4)],
        [['2/m'+'c'],(0,1,2,5)],
        [['-1'],(0,1,2,3,4,5)],
        ]

        seqData,histData = self.RefData(hist)
        hId = histData['data'][0]['hId']
        phasedict = self.proj.phase(phase).data
        pId = phasedict['pId']
        pfx = str(pId)+'::' # prefix for A values from phase
        phfx = '%d:%d:'%(pId,hId) # Dij prefix
        # get unit cell & symmetry for phase
        RecpCellTerms = G2lat.cell2A(phasedict['General']['Cell'][1:7])
        zeroDict = {pfx+'A'+str(i):0.0 for i in range(6)}
        SGdata = phasedict['General']['SGData']
        # determine the cell items not defined by symmetry
        laue = SGdata['SGLaue'][:]
        if laue == '2/m':
            laue += SGdata['SGUniq']
        for symlist,celllist in uniqCellLookup:
            if laue in symlist:
                uniqCellIndx = celllist
                break
        else: # should not happen
            uniqCellIndx = list(range(6))
        initialCell = {}
        for i in uniqCellIndx:
            initialCell[str(pId)+'::A'+str(i)] =  RecpCellTerms[i]

        esdLookUp = {}
        dLookup = {}
        # note that varyList keys are p:h:Dij while newCellDict keys are p::Dij
        for nKey in seqData['newCellDict']:
            p = nKey.find('::')+1
            vKey = nKey[:p] + str(hId) + nKey[p:]
            if vKey in seqData['varyList']:
                esdLookUp[self.newCellDict[nKey][0]] = nKey
                dLookup[nKey] = self.newCellDict[nKey][0]
        covData = {'varyList': [dLookup.get(striphist(v),v) for v in seqData['varyList']],
                'covMatrix': seqData['covMatrix']}
        A = RecpCellTerms[:] # make copy of starting A values

        for i,j in enumerate(('D11','D22','D33','D12','D13','D23')):
            var = pfx+'A'+str(i)
            Dvar = phfx+j
            # apply Dij value if non-zero
            if Dvar in seqData['parmDict']:
                A[i] += seqData['parmDict'][Dvar]
            # override with fit result if is Dij varied
            try:
                A[i] = seqData['newCellDict'][esdLookUp[var]][1] # get refined value 
            except KeyError:
                pass
        Albls = [pfx+'A'+str(i) for i in range(6)]
        cellDict = dict(zip(Albls,A))

        
        A,zeros = G2stIO.cellFill(pfx,SGdata,cellDict,zeroDict)
        # convert to direct cell
        c = G2lat.A2cell(A)
        vol = G2lat.calc_V(A)
        cE = G2stIO.getCellEsd(pfx,SGdata,A,covData)
        return list(c)+[vol],cE,uniqCellIndx
    
    def get_VaryList(self,hist):
        '''Returns a list of the refined variables in the 
        last refinement cycle for the selected histogram

        :param hist: Specify a histogram or using the histogram name (str) 
          or the index number (int) of the histogram in the sequential 
          refinement (not the project), numbered starting from 0.
        :returns: a list of variables or None if no refinement has been
          performed. 
        '''
        try:
            seqData,histData = self.RefData(hist)
            return seqData['varyList']
        except:
            return

    def get_ParmList(self,hist):
        '''Returns a list of all the parameters defined in the 
        last refinement cycle for the selected histogram

        :param hist: Specify a histogram or using the histogram name (str) 
          or the index number (int) of the histogram in the sequential 
          refinement (not the project), numbered as in the project tree 
          starting from 0.
        :returns: a list of parameters or None if no refinement has been
          performed. 
        '''
        try:
            seqData,histData = self.RefData(hist)
            return list(seqData['parmDict'].keys())
        except:
            return

    def get_Variable(self,hist,var):
        '''Returns the value and standard uncertainty (esd) for a variable 
        parameters, as defined for the selected histogram 
        in the last sequential refinement cycle

        :param hist: Specify a histogram or using the histogram name (str) 
          or the index number (int) of the histogram in the sequential 
          refinement (not the project), numbered as in the project tree 
          starting from 0.
        :param str var: a variable name of form '<p>:<h>:<name>', such as 
          ':0:Scale'
        :returns: (value,esd) if the parameter is refined or 
          (value, None) if the variable is in a constraint or is not 
          refined or None if the parameter is not found. 
        '''
        try:
            seqData,histData = self.RefData(hist)
            val = seqData['parmDict'][var]
        except:
            return
        try:
            pos = seqData['varyList'].index(var)
            esd = np.sqrt(seqData['covMatrix'][pos,pos])
            return (val,esd)
        except ValueError:
            return (val,None)

    def get_Covariance(self,hist,varList):
        '''Returns the values and covariance matrix for a series of variable 
        parameters, as defined for the selected histogram 
        in the last sequential refinement cycle

        :param hist: Specify a histogram or using the histogram name (str) 
          or the index number (int) of the histogram in the sequential 
          refinement (not the project), numbered as in the project tree 
          starting from 0.
        :param tuple varList: a list of variable names of form '<p>:<h>:<name>'
        :returns: (valueList,CovMatrix) where valueList contains the (n) values
          in the same order as varList (also length n) and CovMatrix is a 
          (n x n) matrix. If any variable name is not found in the varyList
          then None is returned. 

        Use this code, where sig provides standard uncertainties for 
        parameters and where covArray provides the correlation between 
        off-diagonal terms::

            sig = np.sqrt(np.diag(covMatrix))
            xvar = np.outer(sig,np.ones_like(sig))
            covArray = np.divide(np.divide(covMatrix,xvar),xvar.T)

        '''
        try:
            seqData,histData = self.RefData(hist)
        except:
            G2fil.G2Print('Warning: Histogram {} not found in the sequential fit'.format(hist))
            return
        missing = [i for i in varList if i not in seqData['varyList']]
        if missing:
            G2fil.G2Print('Warning: Variable(s) {} were not found in the varyList'.format(missing))
            return None
        vals = [seqData['parmDict'][i] for i in varList]
        import GSASIImath as G2mth
        cov = G2mth.getVCov(varList,seqData['varyList'],seqData['covMatrix'])
        return (vals,cov)
    
    def RefData(self,hist):
        '''Provides access to the output from a particular histogram

        :param hist: Specify a histogram or using the histogram name (str) 
          or the index number (int) of the histogram in the sequential 
          refinement (not the project), numbered as in the project tree 
          starting from 0.
        :returns: a list of dicts where the first element has sequential 
          refinement results and the second element has the contents of 
          the histogram tree items.
        '''
        try:
            hist = self.data['histNames'][hist]
        except IndexError:
            raise Exception('Histogram #{} is out of range from the Sequential Refinement'
                                .format(hist))
        except TypeError:
            pass
        if hist not in self.data['histNames']:
            raise Exception('Histogram {} is not included in the Sequential Refinement'
                                .format(hist))
        return self.data[hist],self.proj.histogram(hist).data

class G2PDF(G2ObjectWrapper):
    """Wrapper for a PDF tree entry, containing the information needed to 
    compute a PDF and the S(Q), G(r) etc. after the computation is done. 
    Note that in a GSASIIscriptable script, instances of G2PDF will be created by 
    calls to :meth:`G2Project.add_PDF` or :meth:`G2Project.pdf`, not via calls 
    to :meth:`G2PDF.__init__`.

    Example use of :class:`G2PDF`::

       gpx.add_PDF('250umSiO2.pdfprm',0)
       pdf.set_formula(['Si',1],['O',2])
       pdf.set_background('Container',1,-0.21)
       for i in range(5):
           if pdf.optimize(): break
       pdf.calculate()
       pdf.export(gpx.filename,'S(Q), pdfGUI')
       gpx.save('pdfcalc.gpx')

    .. seealso::
        :meth:`G2Project.pdf`
        :meth:`G2Project.pdfs`
    """
    def __init__(self, data, name, proj):
        self.data = data
        self.name = name
        self.proj = proj
    def set_background(self,btype,histogram,mult=-1.,refine=False):
        '''Sets a histogram to be used as the 'Sample Background',
        the 'Container' or the 'Container Background.'

        :param str btype: Type of background to set, must contain 
          the string 'samp' for Sample Background', 'cont' and 'back' 
          for the 'Container Background' or only 'cont' for the 
          'Container'. Note that capitalization and extra characters 
          are ignored, so the full strings (such as 'Sample 
          Background' & 'Container Background') can be used.
        
        :param histogram: A reference to a histogram,
          which can be reference by object, name, or number.
        :param float mult: a multiplier for the histogram; defaults 
          to -1.0
        :param bool refine: a flag to enable refinement (only 
          implemented for 'Sample Background'); defaults to False
        '''
        if 'samp' in btype.lower():
            key = 'Sample Bkg.'
        elif 'cont' in btype.lower() and 'back' in btype.lower():
            key = 'Container Bkg.'
        elif 'cont' in btype.lower():
            key = 'Container'
        else:
            raise Exception('btype = {} is invalid'.format(btype))
        self.data['PDF Controls'][key]['Name'] = self.proj.histogram(histogram).name
        self.data['PDF Controls'][key]['Mult'] = mult
        self.data['PDF Controls'][key]['Refine'] = refine

    def set_formula(self,*args):
        '''Set the chemical formula for the PDF computation. 
        Use pdf.set_formula(['Si',1],['O',2]) for SiO2.

        :param list item1: The element symbol and number of atoms in formula for first element
        :param list item2: The element symbol and number of atoms in formula for second element,...

        repeat parameters as needed for all elements in the formula.
        '''
        powderHist = self.proj.histogram(self.data['PDF Controls']['Sample']['Name'])
        inst = powderHist.data['Instrument Parameters'][0]
        ElList = self.data['PDF Controls']['ElList']
        ElList.clear()
        sumVol = 0.
        for elem,mult in args:
            ElList[elem] = G2elem.GetElInfo(elem,inst)
            ElList[elem]['FormulaNo'] = mult
            Avol = (4.*np.pi/3.)*ElList[elem]['Drad']**3
            sumVol += Avol*ElList[elem]['FormulaNo']
        self.data['PDF Controls']['Form Vol'] = max(10.0,sumVol)
        
    def calculate(self,xydata=None,limits=None,inst=None):
        '''Compute the PDF using the current parameters. Results are set
        in the PDF object arrays (self.data['PDF Controls']['G(R)'] etc.). 
        Note that if ``xydata``, is specified, the background histograms(s)
        will not be accessed from the project file associated with the current 
        PDF entry. If ``limits`` and ``inst`` are both specified, no histograms
        need be in the current project. However, the self.data['PDF Controls']
        sections ('Sample', 'Sample Bkg.','Container Bkg.') must be  
        non-blank for the corresponding items to be used from``xydata``.

        :param dict xydata: an array containing the Sample's I vs Q, and 
          any or none of the Sample Background, the Container scattering and 
          the Container Background. If xydata is None (default), the values are
          taken from histograms, as named in the PDF's self.data['PDF Controls']
          entries with keys 'Sample', 'Sample Bkg.','Container Bkg.' & 
          'Container'.
        :param list limits: upper and lower Q values to be used for PDF 
          computation. If None (default), the values are
          taken from the Sample histogram's .data['Limits'][1] values. 
        :param dict inst: The Sample histogram's instrument parameters
          to be used for PDF computation. If None (default), the values are
          taken from the Sample histogram's .data['Instrument Parameters'][0]
          values. 
        '''
        data = self.data['PDF Controls']
        if xydata is None:
            xydata = {}
            for key in 'Sample Bkg.','Container Bkg.','Container','Sample':
                name = data[key]['Name'].strip()
                if name:
                    xydata[key] = self.proj.histogram(name).data['data']
        if limits is None:
            name = data['Sample']['Name'].strip()
            limits = self.proj.histogram(name).data['Limits'][1]
        if inst is None:
            name = data['Sample']['Name'].strip()
            inst = self.proj.histogram(name).data['Instrument Parameters'][0]
        G2pwd.CalcPDF(data,inst,limits,xydata)
        data['I(Q)'] = xydata['IofQ']
        data['S(Q)'] = xydata['SofQ']
        data['F(Q)'] = xydata['FofQ']
        data['G(R)'] = xydata['GofR']
        data['g(r)'] = xydata['gofr']

    def optimize(self,showFit=True,maxCycles=5,
                     xydata=None,limits=None,inst=None):
        '''Optimize the low R portion of G(R) to minimize selected 
        parameters. Note that this updates the parameters in the settings 
        (self.data['PDF Controls']) but does not update the PDF object 
        arrays (self.data['PDF Controls']['G(R)'] etc.) with the computed 
        values, use :meth:`calculate` after a fit to do that. 

        :param bool showFit: if True (default) the optimized parameters
          are shown before and after the fit, as well as the RMS value 
          in the minimized region.
        :param int maxCycles: the maximum number of least-squares cycles;
          defaults to 5. 
        :returns: the result from the optimizer as True or False, depending 
          on if the refinement converged.
        :param dict xydata: an array containing the Sample's I vs Q, and 
          any or none of the Sample Background, the Container scattering and 
          the Container Background. If xydata is None (default), the values are
          taken from histograms, as named in the PDF's self.data['PDF Controls']
          entries with keys 'Sample', 'Sample Bkg.','Container Bkg.' & 
          'Container'.
        :param list limits: upper and lower Q values to be used for PDF 
          computation. If None (default), the values are
          taken from the Sample histogram's .data['Limits'][1] values. 
        :param dict inst: The Sample histogram's instrument parameters
          to be used for PDF computation. If None (default), the values are
          taken from the Sample histogram's .data['Instrument Parameters'][0]
          values. 
        '''
        data = self.data['PDF Controls']
        if xydata is None:
            xydata = {}
            for key in 'Sample Bkg.','Container Bkg.','Container','Sample':
                name = data[key]['Name'].strip()
                if name:
                    xydata[key] = self.proj.histogram(name).data['data']
        if limits is None:
            name = data['Sample']['Name'].strip()
            limits = self.proj.histogram(name).data['Limits'][1]
        if inst is None:
            name = data['Sample']['Name'].strip()
            inst = self.proj.histogram(name).data['Instrument Parameters'][0]

        res = G2pwd.OptimizePDF(data,xydata,limits,inst,showFit,maxCycles)
        return res['success']
        
    def export(self,fileroot,formats):
        '''Write out the PDF-related data (G(r), S(Q),...) into files

        :param str fileroot: name of file(s) to be written. The extension 
          will be ignored and set to .iq, .sq, .fq or .gr depending 
          on the formats selected.
        :param str formats: string specifying the file format(s) to be written,
          should contain at least one of the following keywords:
          I(Q), S(Q), F(Q), G(r) and/or PDFgui (capitalization and
          punctuation is ignored). Note that G(r) and PDFgui should not 
          be specifed together.
        '''
        PDFsaves = 5*[False]
        PDFentry = self.name
        name = self.data['PDF Controls']['Sample']['Name'].strip()
        limits = self.proj.histogram(name).data['Limits']
        inst = self.proj.histogram(name).data['Instrument Parameters'][0]
        for i,lbl in enumerate(['I(Q)', 'S(Q)', 'F(Q)', 'G(r)', 'PDFgui']):
            PDFsaves[i] = lbl.lower() in formats.lower()
        G2fil.PDFWrite(PDFentry,fileroot,PDFsaves,self.data['PDF Controls'],inst,limits)

blkSize = 256   #256 seems to be optimal; will break in polymask if >1024
'Integration block size; 256 seems to be optimal, must be <=1024 (for polymask)'

def calcMaskMap(imgprms,mskprms):
    '''Computes the mask array for a set of image controls and mask parameters
    '''
    return G2img.MakeUseMask(imgprms,mskprms,blkSize)

def calcThetaAzimMap(imgprms):
    '''Computes the array for theta-azimuth mapping for a set of image controls
    '''
    return G2img.MakeUseTA(imgprms,blkSize)

class G2Image(G2ObjectWrapper):
    '''Wrapper for an IMG tree entry, containing an image and associated metadata. 

    Note that in a GSASIIscriptable script, instances of G2Image will be created by 
    calls to :meth:`G2Project.add_image` or :meth:`G2Project.images`. 
    Scripts will not use ``G2Image()`` to call :meth:`G2Image.__init__` directly. 
    The object contains these class variables:

        * G2Image.proj: contains a reference to the :class:`G2Project`
          object that contains this image 
        * G2Image.name: contains the name of the image
        * G2Image.data: contains the image's associated data in a dict,
          as documented for the :ref:`Image Data Structure<Image_table>`.

    Example use of G2Image:

    >>> gpx = G2sc.G2Project(filename='itest.gpx')
    >>> imlst = gpx.add_image(idata,fmthint="TIF")
    >>> imlst[0].loadControls('stdSettings.imctrl')
    >>> imlst[0].setCalibrant('Si    SRM640c')
    >>> imlst[0].loadMasks('stdMasks.immask')
    >>> imlst[0].Recalibrate()
    >>> imlst[0].setControl('outAzimuths',3)
    >>> pwdrList = imlst[0].Integrate()
 
    More detailed image processing examples are shown at :ref:`ImageProc`.

    '''
    # parameters in that can be accessed via setControl. This may need future attention
    ControlList = {
        'int': ['calibskip', 'pixLimit', 'edgemin', 'outChannels',
                    'outAzimuths'],
        'float': ['cutoff', 'setdist', 'wavelength', 'Flat Bkg',
                      'azmthOff', 'tilt', 'calibdmin', 'rotation',
                      'distance', 'DetDepth'],
        'bool': ['setRings', 'setDefault', 'centerAzm', 'fullIntegrate',
                     'DetDepthRef', 'showLines'],
        'str': ['SampleShape', 'binType', 'formatName', 'color',
                    'type', ],
        'list': ['GonioAngles', 'IOtth', 'LRazimuth', 'Oblique', 'PolaVal',
                   'SampleAbs', 'center', 'ellipses', 'linescan',
                    'pixelSize', 'range', 'ring', 'rings', 'size', ],
        'dict': ['varyList'],
        }
    '''Defines the items known to exist in the Image Controls tree section 
    and the item's data types. A few are not included here
    ('background image', 'dark image', 'Gain map', and 'calibrant') because
    these items have special set routines,
    where references to entries are checked to make sure their values are
    correct.
    ''' 
        
    def __init__(self, data, name, proj):
        self.data = data
        self.name = name
        self.proj = proj

    def setControl(self,arg,value):
        '''Set an Image Controls parameter in the current image.
        If the parameter is not found an exception is raised.

        :param str arg: the name of a parameter (dict entry) in the 
          image. The parameter must be found in :data:`ControlList`
          or an exception is raised.
        :param value: the value to set the parameter. The value is 
          cast as the appropriate type from :data:`ControlList`.
        '''
        for typ in self.ControlList:
            if arg in self.ControlList[typ]: break
        else:
            G2fil.G2Print('Allowed args:\n',[nam for nam,typ in self.findControl('')])
            raise Exception('arg {} not defined in G2Image.setControl'
                                .format(arg))
        try:
            if typ == 'int':
                self.data['Image Controls'][arg] = int(value)
            elif typ == 'float':
                self.data['Image Controls'][arg] = float(value)
            elif typ == 'bool':
                self.data['Image Controls'][arg] = bool(value)
            elif typ == 'str':
                self.data['Image Controls'][arg] = str(value)
            elif typ == 'list':
                self.data['Image Controls'][arg] = list(value)
            elif typ == 'dict':
                self.data['Image Controls'][arg] = dict(value)
            else:
                raise Exception('Unknown type {} for arg {} in  G2Image.setControl'
                                    .format(typ,arg))
        except:
            raise Exception('Error formatting value {} as type {} for arg {} in  G2Image.setControl'
                                    .format(value,typ,arg))

    def getControl(self,arg):
        '''Return an Image Controls parameter in the current image.
        If the parameter is not found an exception is raised.

        :param str arg: the name of a parameter (dict entry) in the 
          image. 
        :returns: the value as a int, float, list,...
        '''
        if arg in self.data['Image Controls']:
            return self.data['Image Controls'][arg]
        G2fil.G2Print(self.findControl(''))
        raise Exception('arg {} not defined in G2Image.getControl'.format(arg))

    def findControl(self,arg=''):
        '''Finds the Image Controls parameter(s) in the current image
        that match the string in arg. Default is '' which returns all 
        parameters.

            Example: 

            >>> findControl('calib')
            [['calibskip', 'int'], ['calibdmin', 'float'], ['calibrant', 'str']]

        :param str arg: a string containing part of the name of a 
          parameter (dict entry) in the image's Image Controls. 
        :returns: a list of matching entries in form 
          [['item','type'], ['item','type'],...] where each 'item' string 
          contains the sting in arg.
        '''
        matchList = []
        for typ in self.ControlList:
            for item in self.ControlList[typ]:
                if arg in item:
                    matchList.append([item,typ])
        return matchList

    def setCalibrant(self,calib):
        '''Set a calibrant for the current image

        :param str calib: specifies a calibrant name which must be one of
          the entries in file ImageCalibrants.py. This is validated and
          an error provides a list of valid choices.
        '''
        import ImageCalibrants as calFile
        if calib in calFile.Calibrants.keys():
            self.data['Image Controls']['calibrant'] = calib
            return
        G2fil.G2Print('Calibrant {} is not valid. Valid calibrants'.format(calib))
        for i in calFile.Calibrants.keys():
            if i: G2fil.G2Print('\t"{}"'.format(i))
        
    def setControlFile(self,typ,imageRef,mult=None):
        '''Set a image to be used as a background/dark/gain map image

        :param str typ: specifies image type, which must be one of:
           'background image', 'dark image', 'gain map'; N.B. only the first
           four characters must be specified and case is ignored.
        :param imageRef: A reference to the desired image. Either the Image 
          tree name (str), the image's index (int) or
          a image object (:class:`G2Image`)
        :param float mult: a multiplier to be applied to the image (not used 
          for 'Gain map'; required for 'background image', 'dark image'
        '''
        if 'back' in typ.lower():
            key = 'background image'
            mult = float(mult)
        elif 'dark' in typ.lower():
            key = 'dark image'
            mult = float(mult)
        elif 'gain' in typ.lower():
            #key = 'Gain map'
            if mult is not None:
                G2fil.G2Print('Warning: Ignoring multiplier for Gain map')
            mult = None
        else:
            raise Exception("Invalid typ {} for setControlFile".format(typ))
        imgNam = self.proj.image(imageRef).name
        if mult is None:
            self.data['Image Controls']['Gain map'] = imgNam
        else:
            self.data['Image Controls'][key] = [imgNam,mult]

    def loadControls(self,filename=None,imgDict=None):
        '''load controls from a .imctrl file

        :param str filename: specifies a file to be read, which should end 
          with .imctrl (defaults to None, meaning parameters are input 
          with imgDict.)
        :param dict imgDict: contains a set of image parameters (defaults to 
          None, meaning parameters are input with filename.)
        '''
        if filename:
            File = open(filename,'r')
            Slines = File.readlines()
            File.close()
            G2fil.LoadControls(Slines,self.data['Image Controls'])
            G2fil.G2Print('file {} read into {}'.format(filename,self.name))
        elif imgDict:
            self.data['Image Controls'].update(imgDict)
            G2fil.G2Print('Image controls set')
        else:
            raise Exception("loadControls called without imgDict or filename specified")

    def saveControls(self,filename):
        '''write current controls values to a .imctrl file

        :param str filename: specifies a file to write, which should end 
          with .imctrl
        '''
        G2fil.WriteControls(filename,self.data['Image Controls'])
        G2fil.G2Print('file {} written from {}'.format(filename,self.name))

    def getControls(self,clean=False):
        '''returns current Image Controls as a dict

        :param bool clean: causes the calbration information to be deleted
        '''
        ImageControls = copy.deepcopy(self.data['Image Controls'])
        if clean:
            ImageControls['showLines'] = True
            ImageControls['ring'] = []
            ImageControls['rings'] = []
            ImageControls['ellipses'] = []
            ImageControls['setDefault'] = False
            for i in 'range','size','GonioAngles':
                if i in ImageControls: del ImageControls[i]
        return ImageControls
    
    def setControls(self,controlsDict):
        '''uses dict from :meth:`getControls` to set Image Controls for current image
        '''
        self.data['Image Controls'].update(copy.deepcopy(controlsDict))
    
    
    def loadMasks(self,filename,ignoreThreshold=False):
        '''load masks from a .immask file

        :param str filename: specifies a file to be read, which should end 
          with .immask
        :param bool ignoreThreshold: If True, masks are loaded with
          threshold masks. Default is False which means any Thresholds 
          in the file are ignored.
        '''
        G2fil.readMasks(filename,self.data['Masks'],ignoreThreshold)
        G2fil.G2Print('file {} read into {}'.format(filename,self.name))

    def initMasks(self):
        '''Initialize Masks, including resetting the Thresholds values
        '''
        self.data['Masks'] = {'Points':[],'Rings':[],'Arcs':[],'Polygons':[],'Frames':[]}
        ImageZ = _getCorrImage(Readers['Image'],self.proj,self)
        Imin = max(0.,np.min(ImageZ))
        Imax = np.max(ImageZ)
        self.data['Masks']['Thresholds'] = [(0,Imax),[Imin,Imax]]
        
    def getMasks(self):
        '''load masks from an IMG tree entry
        '''
        return self.data['Masks']
        
    def setMasks(self,maskDict,resetThresholds=False):
        '''load masks dict (from :meth:`getMasks`) into current IMG record

        :param dict maskDict: specifies a dict with image parameters, 
          from :meth:`getMasks`
        :param bool resetThresholds: If True, Threshold Masks in the 
          dict are ignored. The default is False which means Threshold 
          Masks are retained.
        '''
        self.data['Masks'] = copy.deepcopy(maskDict)
        if resetThresholds:
            ImageZ = _getCorrImage(Readers['Image'],self.proj,self)
            Imin = max(0.,np.min(ImageZ))
            Imax = np.max(ImageZ)
            self.data['Masks']['Thresholds'] [(0,Imax),[Imin,Imax]]        

    def getVary(self,*args):
        '''Return the refinement flag(s) for Image Controls parameter(s)
        in the current image.
        If the parameter is not found, an exception is raised.

        :param str arg: the name of a refinement parameter in the 
          varyList for the image. The name should be one of 
          'dep', 'det-X', 'det-Y', 'dist', 'phi', 'tilt', or 'wave'
        :param str arg1: the name of a parameter (dict entry) as before,
          optional


        :returns: a list of bool value(s)
        '''
        res = []
        for arg in args:
            if arg in self.data['Image Controls']['varyList']:
                res.append(self.data['Image Controls']['varyList'][arg])
            else:
                raise Exception('arg {} not defined in G2Image.getVary'.format(arg))
        return res
    
    def setVary(self,arg,value):
        '''Set a refinement flag for Image Controls parameter in the 
        current image.
        If the parameter is not '*' or found, an exception is raised.

        :param str arg: the name of a refinement parameter in the 
          varyList for the image. The name should be one of 
          'dep', 'det-X', 'det-Y', 'dist', 'phi', 'tilt', or 'wave',
          or it may be a list or tuple of names, 
          or it may be '*' in which all parameters are set accordingly.
        :param value: the value to set the parameter. The value is 
          cast as the appropriate type from :data:`ControlList`.
        '''
        if arg == '*':
            for a in self.data['Image Controls']['varyList']:
                self.data['Image Controls']['varyList'][a] = bool(value)
            return
        if not isinstance(arg,tuple) and not isinstance(arg,list):
            arg = [arg]
        for a in arg:
            if a in self.data['Image Controls']['varyList']:
                self.data['Image Controls']['varyList'][a] = bool(value)
            else:
                raise Exception('arg {} not defined in G2Image.setVary'.format(a))

    def Recalibrate(self):
        '''Invokes a recalibration fit (same as Image Controls/Calibration/Recalibrate 
        menu command). Note that for this to work properly, the calibration 
        coefficients (center, wavelength, distance & tilts) must be fairly close.
        This may produce a better result if run more than once.
        '''
        LoadG2fil()
        ImageZ = _getCorrImage(Readers['Image'],self.proj,self)
        G2img.ImageRecalibrate(None,ImageZ,self.data['Image Controls'],self.data['Masks'])

    def Integrate(self,name=None,MaskMap=None,ThetaAzimMap=None):
        '''Invokes an image integration (same as Image Controls/Integration/Integrate
        menu command). All parameters will have previously been set with Image Controls
        so no input is needed here. However, the optional parameters MaskMap 
        and ThetaAzimMap may be supplied to save computing these items more than
        once, speeding integration of multiple images with the same 
        image/mask parameters.

        Note that if integration is performed on an 
        image more than once, histogram entries may be overwritten. Use the name
        parameter to prevent this if desired. 

        :param str name: base name for created histogram(s). If None (default), 
          the histogram name is taken from the image name. 
        :param list MaskMap: from :func:`calcMaskMap` 
        :param list ThetaAzimMap: from :func:`calcThetaAzimMap`
        :returns: a list of created histogram (:class:`G2PwdrData`) objects.
        '''
        ImageZ = _getCorrImage(Readers['Image'],self.proj,self)
        # do integration
        ints,azms,Xvals,cancel = G2img.ImageIntegrate(ImageZ,
                self.data['Image Controls'],self.data['Masks'],blkSize=blkSize,
                useMask=MaskMap,useTA=ThetaAzimMap)
        # code from here on based on G2IO.SaveIntegration, but places results in the current
        # project rather than tree
        X = Xvals[:-1]
        N = len(X)

        data = self.data['Image Controls']
        Comments = self.data['Comments']
        # make name from image, unless overridden
        if name:
            if not name.startswith(data['type']+' '):
                name = data['type']+' '+name
        else:
            name = self.name.replace('IMG ',data['type']+' ')
        if 'PWDR' in name:
            if 'target' in data:
                names = ['Type','Lam1','Lam2','I(L2)/I(L1)','Zero','Polariz.','U','V','W','X','Y','Z','SH/L','Azimuth'] 
                codes = [0 for i in range(14)]
            else:
                names = ['Type','Lam','Zero','Polariz.','U','V','W','X','Y','Z','SH/L','Azimuth'] 
                codes = [0 for i in range(12)]
        elif 'SASD' in name:
            names = ['Type','Lam','Zero','Azimuth'] 
            codes = [0 for i in range(4)]
            X = 4.*np.pi*npsind(X/2.)/data['wavelength']    #convert to q
        Xminmax = [X[0],X[-1]]
        Azms = []
        dazm = 0.
        if data['fullIntegrate'] and data['outAzimuths'] == 1:
            Azms = [45.0,]                              #a poor man's average?
        else:
            for i,azm in enumerate(azms[:-1]):
                if azm > 360. and azms[i+1] > 360.:
                    Azms.append(G2img.meanAzm(azm%360.,azms[i+1]%360.))
                else:    
                    Azms.append(G2img.meanAzm(azm,azms[i+1]))
            dazm = np.min(np.abs(np.diff(azms)))/2.
        # pull out integration results and make histograms for each
        IntgOutList = []
        for i,azm in enumerate(azms[:-1]):
            Aname = name+" Azm= %.2f"%((azm+dazm)%360.)
            # MT dict to contain histogram
            HistDict = {}
            histItems = [Aname]
            Sample = G2obj.SetDefaultSample()       #set as Debye-Scherrer
            Sample['Gonio. radius'] = data['distance']
            Sample['Omega'] = data['GonioAngles'][0]
            Sample['Chi'] = data['GonioAngles'][1]
            Sample['Phi'] = data['GonioAngles'][2]
            Sample['Azimuth'] = (azm+dazm)%360.    #put here as bin center 
            polariz = 0.99    #set default polarization for synchrotron radiation!
            for item in Comments:
                if 'polariz' in item:
                    try:
                        polariz = float(item.split('=')[1])
                    except:
                        polariz = 0.99
                for key in ('Temperature','Pressure','Time','FreePrm1','FreePrm2','FreePrm3','Omega',
                    'Chi','Phi'):
                    if key.lower() in item.lower():
                        try:
                            Sample[key] = float(item.split('=')[1])
                        except:
                            pass
                # if 'label_prm' in item.lower():
                #     for num in ('1','2','3'):
                #         if 'label_prm'+num in item.lower():
                #             Controls['FreePrm'+num] = item.split('=')[1].strip()
            if 'PWDR' in Aname:
                if 'target' in data:    #from lab x-ray 2D imaging data
                    waves = {'CuKa':[1.54051,1.54433],'TiKa':[2.74841,2.75207],'CrKa':[2.28962,2.29351],
                                 'FeKa':[1.93597,1.93991],'CoKa':[1.78892,1.79278],'MoKa':[0.70926,0.713543],
                                 'AgKa':[0.559363,0.563775]}
                    wave1,wave2 = waves[data['target']]
                    parms = ['PXC',wave1,wave2,0.5,0.0,polariz,290.,-40.,30.,6.,-14.,0.0,0.0001,Azms[i]]
                else:
                    parms = ['PXC',data['wavelength'],0.0,polariz,1.0,-0.10,0.4,0.30,1.0,0.0,0.0001,Azms[i]]
            elif 'SASD' in Aname:
                Sample['Trans'] = data['SampleAbs'][0]
                parms = ['LXC',data['wavelength'],0.0,Azms[i]]
            Y = ints[i]
            Ymin = np.min(Y)
            Ymax = np.max(Y)
            W = np.where(Y>0.,1./Y,1.e-6)                    #probably not true
            section = 'Comments'
            histItems += [section]
            HistDict[section] = Comments
            section = 'Limits'
            histItems += [section]
            HistDict[section] = copy.deepcopy([tuple(Xminmax),Xminmax])
            if 'PWDR' in Aname:
                section = 'Background'
                histItems += [section]
                HistDict[section] = [['chebyschev-1',1,3,1.0,0.0,0.0],
                    {'nDebye':0,'debyeTerms':[],'nPeaks':0,'peaksList':[]}]
            inst = [dict(zip(names,zip(parms,parms,codes))),{}]
            for item in inst[0]:
                inst[0][item] = list(inst[0][item])
            section = 'Instrument Parameters'
            histItems += [section]
            HistDict[section] = inst
            if 'PWDR' in Aname:
                section = 'Sample Parameters'
                histItems += [section]
                HistDict[section] = Sample
                section = 'Peak List'
                histItems += [section]
                HistDict[section] = {'sigDict':{},'peaks':[]}
                section = 'Index Peak List'
                histItems += [section]
                HistDict[section] = [[],[]]
                section = 'Unit Cells List'
                histItems += [section]
                HistDict[section] = []
                section = 'Reflection Lists'
                histItems += [section]
                HistDict[section] = {}
            # elif 'SASD' in Aname:             
            #     section = 'Substances'
            #     histItems += [section]
            #     HistDict[section] = G2pdG.SetDefaultSubstances()  # this needs to be moved
            #     section = 'Sample Parameters'
            #     histItems += [section]
            #     HistDict[section] = Sample
            #     section = 'Models'
            #     histItems += [section]
            #     HistDict[section] = G2pdG.SetDefaultSASDModel() # this needs to be moved
            valuesdict = {
                'wtFactor':1.0,'Dummy':False,'ranId':ran.randint(0,sys.maxsize),'Offset':[0.0,0.0],'delOffset':0.02*Ymax,
                'refOffset':-0.1*Ymax,'refDelt':0.1*Ymax,'Yminmax':[Ymin,Ymax]}
            # if Aname is already in the project replace it
            for j in self.proj.names:
                if j[0] == Aname: 
                    G2fil.G2Print('Replacing "{}" in project'.format(Aname))
                    break
            else:
                G2fil.G2Print('Adding "{}" to project'.format(Aname))
                self.proj.names.append([Aname]+
                        [u'Comments',u'Limits',u'Background',u'Instrument Parameters',
                         u'Sample Parameters', u'Peak List', u'Index Peak List',
                         u'Unit Cells List', u'Reflection Lists'])
            HistDict['data'] = [valuesdict,
                    [np.array(X),np.array(Y),np.array(W),np.zeros(N),np.zeros(N),np.zeros(N)]]
            self.proj.data[Aname] = HistDict
            IntgOutList.append(self.proj.histogram(Aname))
        return IntgOutList

##########################
# Command Line Interface #
##########################
# Each of these takes an argparse.Namespace object as their argument,
# representing the parsed command-line arguments for the relevant subcommand.
# The argument specification for each is in the subcommands dictionary (see
# below)

commandhelp={}
commandhelp["create"] = "creates a GSAS-II project, optionally adding histograms and/or phases"
def create(args):
    """Implements the create command-line subcommand. This creates a GSAS-II project, optionally adding histograms and/or phases::

  usage: GSASIIscriptable.py create [-h] [-d HISTOGRAMS [HISTOGRAMS ...]]
                                  [-i IPARAMS [IPARAMS ...]]
                                  [-p PHASES [PHASES ...]]
                                  filename
                                  
positional arguments::

  filename              the project file to create. should end in .gpx

optional arguments::

  -h, --help            show this help message and exit
  -d HISTOGRAMS [HISTOGRAMS ...], --histograms HISTOGRAMS [HISTOGRAMS ...]
                        list of datafiles to add as histograms
  -i IPARAMS [IPARAMS ...], --iparams IPARAMS [IPARAMS ...]
                        instrument parameter file, must be one for every
                        histogram
  -p PHASES [PHASES ...], --phases PHASES [PHASES ...]
                        list of phases to add. phases are automatically
                        associated with all histograms given.

    """
    proj = G2Project(gpxfile=args.filename)

    hist_objs = []
    if args.histograms:
        for h,i in zip(args.histograms,args.iparams):
            G2fil.G2Print("Adding histogram from",h,"with instparm ",i)
            hist_objs.append(proj.add_powder_histogram(h, i))

    if args.phases: 
        for p in args.phases:
            G2fil.G2Print("Adding phase from",p)
            proj.add_phase(p, histograms=hist_objs)
        G2fil.G2Print('Linking phase(s) to histogram(s):')
        for h in hist_objs:
            G2fil.G2Print ('   '+h.name)

    proj.save()

commandhelp["add"] = "adds histograms and/or phases to GSAS-II project"
def add(args):
    """Implements the add command-line subcommand. This adds histograms and/or phases to GSAS-II project::

  usage: GSASIIscriptable.py add [-h] [-d HISTOGRAMS [HISTOGRAMS ...]]
                               [-i IPARAMS [IPARAMS ...]]
                               [-hf HISTOGRAMFORMAT] [-p PHASES [PHASES ...]]
                               [-pf PHASEFORMAT] [-l HISTLIST [HISTLIST ...]]
                               filename


positional arguments::

  filename              the project file to open. Should end in .gpx

optional arguments::

  -h, --help            show this help message and exit
  -d HISTOGRAMS [HISTOGRAMS ...], --histograms HISTOGRAMS [HISTOGRAMS ...]
                        list of datafiles to add as histograms
  -i IPARAMS [IPARAMS ...], --iparams IPARAMS [IPARAMS ...]
                        instrument parameter file, must be one for every
                        histogram
  -hf HISTOGRAMFORMAT, --histogramformat HISTOGRAMFORMAT
                        format hint for histogram import. Applies to all
                        histograms
  -p PHASES [PHASES ...], --phases PHASES [PHASES ...]
                        list of phases to add. phases are automatically
                        associated with all histograms given.
  -pf PHASEFORMAT, --phaseformat PHASEFORMAT
                        format hint for phase import. Applies to all phases.
                        Example: -pf CIF
  -l HISTLIST [HISTLIST ...], --histlist HISTLIST [HISTLIST ...]
                        list of histgram indices to associate with added
                        phases. If not specified, phases are associated with
                        all previously loaded histograms. Example: -l 2 3 4
    
    """
    proj = G2Project(args.filename)

    if args.histograms:
        for h,i in zip(args.histograms,args.iparams):
            G2fil.G2Print("Adding histogram from",h,"with instparm ",i)
            proj.add_powder_histogram(h, i, fmthint=args.histogramformat)

    if args.phases: 
        if not args.histlist:
            histlist = proj.histograms()
        else:
            histlist = [proj.histogram(i) for i in args.histlist]

        for p in args.phases:
            G2fil.G2Print("Adding phase from",p)
            proj.add_phase(p, histograms=histlist, fmthint=args.phaseformat)
            
        if not args.histlist:
            G2fil.G2Print('Linking phase(s) to all histogram(s)')
        else:
            G2fil.G2Print('Linking phase(s) to histogram(s):')
            for h in histlist:
                G2fil.G2Print('   '+h.name)
    proj.save()


commandhelp["dump"] = "Shows the contents of a GSAS-II project"
def dump(args):
    """Implements the dump command-line subcommand, which shows the contents of a GSAS-II project::

       usage: GSASIIscriptable.py dump [-h] [-d] [-p] [-r] files [files ...]

positional arguments::

  files

optional arguments::

  -h, --help        show this help message and exit
  -d, --histograms  list histograms in files, overrides --raw
  -p, --phases      list phases in files, overrides --raw
  -r, --raw         dump raw file contents, default
  
    """
    if not args.histograms and not args.phases:
        args.raw = True
    if args.raw:
        import IPython.lib.pretty as pretty

    for fname in args.files:
        if args.raw:
            proj, nameList = LoadDictFromProjFile(fname)
            print("file:", fname)
            print("names:", nameList)
            for key, val in proj.items():
                print(key, ":")
                pretty.pprint(val)
        else:
            proj = G2Project(fname)
            if args.histograms:
                hists = proj.histograms()
                for h in hists:
                    print(fname, "hist", h.id, h.name)
            if args.phases:
                phase_list = proj.phases()
                for p in phase_list:
                    print(fname, "phase", p.id, p.name)


commandhelp["browse"] = "Load a GSAS-II project and then open a IPython shell to browse it"
def IPyBrowse(args):
    """Load a .gpx file and then open a IPython shell to browse it::

  usage: GSASIIscriptable.py browse [-h] files [files ...]

positional arguments::

  files       list of files to browse

optional arguments::

  -h, --help  show this help message and exit

    """
    for fname in args.files:
        proj, nameList = LoadDictFromProjFile(fname)
        msg = "\nfname {} loaded into proj (dict) with names in nameList".format(fname)
        GSASIIpath.IPyBreak_base(msg)
        break


commandhelp["refine"] = '''
Conducts refinements on GSAS-II projects according to a list of refinement
steps in a JSON dict
'''
def refine(args):
    """Implements the refine command-line subcommand:
    conducts refinements on GSAS-II projects according to a JSON refinement dict::

        usage: GSASIIscriptable.py refine [-h] gpxfile [refinements]

positional arguments::

  gpxfile      the project file to refine
  refinements  json file of refinements to apply. if not present refines file
               as-is

optional arguments::

  -h, --help   show this help message and exit
  
    """
    proj = G2Project(args.gpxfile)
    if args.refinements is None:
        proj.refine()
    else:
        import json
        with open(args.refinements) as refs:
            refs = json.load(refs)
        if type(refs) is not dict:
            raise G2ScriptException("Error: JSON object must be a dict.")
        if "code" in refs:
            print("executing code:\n|  ",'\n|  '.join(refs['code']))
            exec('\n'.join(refs['code']))
        proj.do_refinements(refs['refinements'])


commandhelp["export"] = "Export phase as CIF"
def export(args):
    """Implements the export command-line subcommand: Exports phase as CIF::

      usage: GSASIIscriptable.py export [-h] gpxfile phase exportfile

positional arguments::

  gpxfile     the project file from which to export
  phase       identifier of phase to export
  exportfile  the .cif file to export to

optional arguments::

  -h, --help  show this help message and exit

    """
    proj = G2Project(args.gpxfile)
    phase = proj.phase(args.phase)
    phase.export_CIF(args.exportfile)


def _args_kwargs(*args, **kwargs):
    return args, kwargs

# A dictionary of the name of each subcommand, and a tuple
# of its associated function and a list of its arguments
# The arguments are passed directly to the add_argument() method
# of an argparse.ArgumentParser

subcommands = {"create":
               (create, [_args_kwargs('filename',
                                      help='the project file to create. should end in .gpx'),

                         _args_kwargs('-d', '--histograms',
                                      nargs='+',
                                      help='list of datafiles to add as histograms'),
                                      
                         _args_kwargs('-i', '--iparams',
                                      nargs='+',
                                      help='instrument parameter file, must be one'
                                           ' for every histogram'
                                      ),

                         _args_kwargs('-p', '--phases',
                                      nargs='+',
                                      help='list of phases to add. phases are '
                                           'automatically associated with all '
                                           'histograms given.')]),
               "add": (add, [_args_kwargs('filename',
                                      help='the project file to open. Should end in .gpx'),

                         _args_kwargs('-d', '--histograms',
                                      nargs='+',
                                      help='list of datafiles to add as histograms'),
                                      
                         _args_kwargs('-i', '--iparams',
                                      nargs='+',
                                      help='instrument parameter file, must be one'
                                           ' for every histogram'
                                      ),
                                      
                         _args_kwargs('-hf', '--histogramformat',
                                      help='format hint for histogram import. Applies to all'
                                           ' histograms'
                                      ),

                         _args_kwargs('-p', '--phases',
                                      nargs='+',
                                      help='list of phases to add. phases are '
                                           'automatically associated with all '
                                           'histograms given.'),

                         _args_kwargs('-pf', '--phaseformat',
                                      help='format hint for phase import. Applies to all'
                                           ' phases. Example: -pf CIF'
                                      ),
                                      
                         _args_kwargs('-l', '--histlist',
                                      nargs='+',
                                      help='list of histgram indices to associate with added'
                                           ' phases. If not specified, phases are'
                                           ' associated with all previously loaded'
                                           ' histograms. Example: -l 2 3 4')]),
                                           
               "dump": (dump, [_args_kwargs('-d', '--histograms',
                                     action='store_true',
                                     help='list histograms in files, overrides --raw'),

                               _args_kwargs('-p', '--phases',
                                            action='store_true',
                                            help='list phases in files, overrides --raw'),

                               _args_kwargs('-r', '--raw',
                                      action='store_true', help='dump raw file contents, default'),

                               _args_kwargs('files', nargs='+')]),

               "refine":
               (refine, [_args_kwargs('gpxfile', help='the project file to refine'),
                         _args_kwargs('refinements',
                                      help='JSON file of refinements to apply. if not present'
                                           ' refines file as-is',
                                      default=None,
                                      nargs='?')]),

               "export": (export, [_args_kwargs('gpxfile',
                                                help='the project file from which to export'),
                                   _args_kwargs('phase', help='identifier of phase to export'),
                                   _args_kwargs('exportfile', help='the .cif file to export to')]),
               "browse": (IPyBrowse, [_args_kwargs('files', nargs='+',
                                                   help='list of files to browse')])}


def main():
    '''The command-line interface for calling GSASIIscriptable as a shell command,
    where it is expected to be called as::

       python GSASIIscriptable.py <subcommand> <file.gpx> <options>

    The following subcommands are defined:

        * create, see :func:`create`
        * add, see :func:`add`
        * dump, see :func:`dump`
        * refine, see :func:`refine`
        * export, :func:`export`
        * browse, see :func:`IPyBrowse`

    .. seealso::
        :func:`create`
        :func:`add`
        :func:`dump`
        :func:`refine`
        :func:`export`
        :func:`IPyBrowse`
    '''
    parser = argparse.ArgumentParser(description=
        "Use of "+os.path.split(__file__)[1]+" Allows GSAS-II actions from command line."
        )
    subs = parser.add_subparsers()

    # Create all of the specified subparsers
    for name, (func, args) in subcommands.items():
        new_parser = subs.add_parser(name,help=commandhelp.get(name),
                                     description='Command "'+name+'" '+commandhelp.get(name))
        for listargs, kwds in args:
            new_parser.add_argument(*listargs, **kwds)
        new_parser.set_defaults(func=func)

    # Parse and trigger subcommand
    result = parser.parse_args()
    result.func(result)

def dictDive(d,search="",keylist=[],firstcall=True,l=None):
    '''Recursive routine to scan a nested dict. Reports a list of keys 
    and the associated type and value for that key. 

    :param dict d: a dict that will be scanned
    :param str search: an optional search string. If non-blank,
       only entries where one of the keys constains search (case ignored)
    :param list keylist: a list of keys to apply to the dict.
    :param bool firstcall: do not specify
    :param list l: do not specify

    :returns: a list of keys located by this routine
       in form [([keylist], type, value),...] where if keylist is ['a','b','c']
       then d[['a']['b']['c'] will have the value.

    This routine can be called in a number of ways, as are shown in a few 
    examples: 

    >>> for i in G2sc.dictDive(p.data['General'],'paw'): print(i)
    ... 
    (['Pawley dmin'], <class 'float'>, 1.0)
    (['doPawley'], <class 'bool'>, False)
    (['Pawley dmax'], <class 'float'>, 100.0)
    (['Pawley neg wt'], <class 'float'>, 0.0)
    >>>
    >>> for i in G2sc.dictDive(p.data,'paw',['General']): print(i)
    ... 
    (['General', 'Pawley dmin'], <class 'float'>, 1.0)
    (['General', 'doPawley'], <class 'bool'>, False)
    (['General', 'Pawley dmax'], <class 'float'>, 100.0)
    (['General', 'Pawley neg wt'], <class 'float'>, 0.0)
    >>>
    >>> for i in G2sc.dictDive(p.data,'',['General','doPawley']): print(i)
    ... 
    (['General', 'doPawley'], <class 'bool'>, False)

    '''
    if firstcall:
        for k in keylist:
            d = d[k]
        l = []
#        first = True
    if type(d) is dict:
        [dictDive(d[key],search,keylist+[key],False,l) for key in d]
    elif search:
        for key in keylist:
            if search.lower() in key.lower():
                l.append((keylist,type(d),d))
                break
        return
    else:
        l.append((keylist,type(d),d))
    if firstcall:
        return l
        
if __name__ == '__main__':
    #fname='/tmp/corundum-template.gpx'
    #prj = G2Project(fname)
    main()
