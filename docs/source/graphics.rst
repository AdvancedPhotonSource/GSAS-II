*Visualization routines*
===============================

*Summary/Contents*
----------------------------

Routines for visualization, using matplotlib and OpenGL graphics. 
Note that the plot toolbar is customized with :class:`GSASIItoolbar` 

.. contents:: Section Contents 

*List of Graphics routines*
-----------------------------

The following plotting routines are defined: 

============================================  ===========================================================================
plotting routine                               action        
============================================  ===========================================================================
:func:`~GSASIIpwdplot.PlotPatterns`            Powder pattern plotting
:func:`~GSASIIpwdplot.PublishRietveldPlot`     Create publication-quality Rietveld plots from :func:`PlotPatterns` plot
:func:`~GSASIIplot.PlotImage`                  Plots of 2D detector images
:func:`~GSASIIplot.PlotPeakWidths`             Plot instrument broadening terms as function of 2-theta/TOF
:func:`~GSASIIplot.PlotCovariance`             Show covariance terms in 2D 
:func:`~GSASIIplot.PlotStructure`              Crystal structure plotting with balls, sticks, lines,
                                               ellipsoids, polyhedra and magnetic moments
:func:`~GSASIIplot.PlotBeadModel`              Plots representation of protein shape from small angle scattering
:func:`~GSASIIplot.Plot1DSngl`                 1D stick plots of structure factors                              
:func:`~GSASIIplot.PlotSngl`                   Structure factor plotting
:func:`~GSASIIplot.Plot3DSngl`                 3D Structure factor plotting
:func:`~GSASIIplot.PlotDeltSig`                Normal probability plot (powder or single crystal)
:func:`~GSASIIplot.PlotISFG`                   PDF analysis: displays I(Q), S(Q), F(Q) and G(r)
:func:`~GSASIIplot.PlotCalib`                  CW or TOF peak calibration
:func:`~GSASIIplot.PlotXY`                     Simple plot of xy data
:func:`~GSASIIplot.PlotXYZ`                    Simple contour plot of xyz data
:func:`~GSASIIplot.PlotXYZvect`                Quiver Plot for 3D cartesian vectors
:func:`~GSASIIplot.Plot3dXYZ`                  Surface Plot for 3D vectors
:func:`~GSASIIplot.PlotAAProb`                 Protein "quality" plot 
:func:`~GSASIIplot.PlotStrain`                 Plot of strain data, used for diagnostic purposes
:func:`~GSASIIplot.PlotSASDSizeDist`           Small angle scattering size distribution plot
:func:`~GSASIIplot.PlotPowderLines`            Plot powder pattern as a stick plot (vertical lines)
:func:`~GSASIIplot.PlotSizeStrainPO`           Plot 3D mustrain/size/preferred orientation figure
:func:`~GSASIIplot.PlotTexture`                Pole figure, inverse pole figure plotting
:func:`~GSASIIplot.ModulationPlot`             Plots modulation information
:func:`~GSASIIplot.PlotTorsion`                Plots MC torsion angles
:func:`~GSASIIplot.PlotRama`                   Ramachandran of energetically allowed regions for dihedral
                                               angles in protein
:func:`~GSASIIplot.PlotSelectedSequence`       Plot one or more sets of values selected from the sequential
                                               refinement table
:func:`~GSASIIplot.PlotIntegration`            Rectified plot of 2D image after image integration with 2-theta and
                                               azimuth as coordinates
:func:`~GSASIIplot.PlotTRImage`                test plot routine
:func:`~GSASIIplot.PlotRigidBody`              show rigid body structures as balls & sticks
:func:`~GSASIIplot.PlotLayers`                 show layer structures as balls & sticks
:func:`~GSASIIplot.PlotFPAconvolutors`         plots the convolutors from Fundamental Parameters
:func:`~GSASIIplot.PlotClusterXYZ`             plots the result of cluster analysis
============================================  ===========================================================================

*Window management routines*
--------------------------------------

The above plotting routines place their graphics in the GSAS-II Plot Window, which contains a
:class:`GSASIIplot.G2PlotNoteBook` tabbed panel allowing multiple plots to be viewed. Methods 
:meth:`GSASIIplot.G2PlotNoteBook.addMpl` (2-D matplotlib), 
:meth:`GSASIIplot.G2PlotNoteBook.add3D` (3-D matplotlib), and 
:meth:`GSASIIplot.G2PlotNoteBook.addOgl` (OpenGL) are used to
create tabbed plot objects to hold plots of the following classes:
:class:`GSASIIplot.G2PlotMpl` (2-D matplotlib), 
:class:`GSASIIplot.G2Plot3D` (3-D matplotlib), and 
:class:`GSASIIplot.G2PlotOgl` (OpenGL). Note that two :class:`GSASIIplot.G2PlotNoteBook` methods are
potentially used to determine how plot updates after a refinement are handled: 

=========================================================     ========================================================
class method                                                  description
=========================================================     ========================================================
:meth:`GSASIIplot.G2PlotNoteBook.RegisterRedrawRoutine`       This specifies a function 
                                                              to redraw the plot after the data tree has been
                                                              reloaded. Be sure this updates data
                                                              objects with new values from the tree, when needed.
                                                  
:meth:`GSASIIplot.G2PlotNoteBook.SetNoDelete`                 Use this to indicate that a plot does not need to be
                                                              updated after a refinement and should not be closed.
=========================================================     ========================================================

These two methods define the following attributes (variables) in the plot tab classes: 

======================    ===============     ============================================================
variable                   default             use
======================    ===============     ============================================================
replotFunction              None               Defines a routine to be called to update the plot 
                                               after a refinement (unless None). Use
                                               :meth:`GSASIIplot.G2PlotNoteBook.RegisterRedrawRoutine`
                                               to define this (and replotArgs & replotKwArgs). 
                                               Plotting functions that take significant time
                                               to complete should probably not use this.)
replotArgs                  []                 Defines the positional arguments to be supplied to
                                               the replotFunction function or method.
replotKwArgs                {}                 Defines the keyword arguments to be supplied to
                                               the replotFunction function or method. 
plotRequiresRedraw         True                If set to True, after a refinement, the plot will be
                                               closed (in :func:`GSASIIdataGUI.GSASII.ResetPlots`)
                                               if it was not updated after the refinement. Set this to
                                               False using
                                               :meth:`GSASIIplot.G2PlotNoteBook.SetNoDelete`
                                               for plots that should not be deleted or do
                                               not change based on refinement results.
plotInvalid                 False              Used to track if a plot has been updated. Set to False
                                               in :meth:`GSASIIplot.G2PlotNoteBook.FindPlotTab` when a
                                               plot is drawn. After a refinement is completed, method
                                               :func:`GSASIIdataGUI.GSASII.ResetPlots` sets
                                               plotInvalid to False for all plots before any routines
                                               are called. 
======================    ===============     ============================================================

*GSASIIplot Classes and Routines*
------------------------------------

.. automodule:: GSASII.GSASIIplot
    :members: 
    :private-members:
    :special-members:

*GSASIIpwdplot Routines*
---------------------------------------

.. automodule:: GSASII.GSASIIpwdplot
    :members: 
    :private-members:
    :special-members:
