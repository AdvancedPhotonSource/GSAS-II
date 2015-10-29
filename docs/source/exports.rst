*GSAS-II Export Modules*
====================================

Exports are implemented by deriving a class from 
:class:`GSASIIIO.ExportBaseClass`. Initialization of 
``self.exporttype`` determines the type of export that will be performed
('project', 'phase', 'single', 'powder', 'image', 'map' or (someday)
'pdf') and of ``self.multiple``
determines if only a single phase, data set, etc. can be exported at a
time (when False) or more than one can be selected.

Powder export routines may optionally define a ``Writer()``
method that accepts the histogram tree name as well as a file name to
be written. This allows :func:`ExportPowder` to use the exporter
independent of the GUI. 

.. automodule:: G2export_examples
    :members: 
    :synopsis: Demonstrates sample code that exports a phase or dataset to
      a text file. 

.. automodule:: G2export_csv
    :members: 
    :synopsis: Exports a phase or dataset to a spreadsheet via a 
       comma-separated-variable (csv) format file.

.. automodule:: G2export_PDB
    :members: 
    :synopsis: Cartesian coordinate export, including PDB format

.. automodule:: G2export_image
    :members: 
    :synopsis: Exports images
 
.. automodule:: G2export_map
    :members: 
    :synopsis: Export Fourier and charge-flip atomic density maps

.. automodule:: G2export_shelx 
    :members: 
    :synopsis: Export a phase in Shelx format

.. automodule:: G2export_CIF
    :members: 
    :synopsis: Export a project in CIF format

.. automodule:: G2export_pwdr
    :members: 
    :synopsis: Export powder data in GSAS and Topas formats

.. automodule:: G2export_FIT2D 
    :members: 
    :synopsis: Export powder data in Fit2D (.chi) format
       
