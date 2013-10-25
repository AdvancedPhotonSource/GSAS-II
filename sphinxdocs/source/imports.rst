*GSAS-II Import Modules*
====================================

Imports are implemented by deriving a class from 
:class:`GSASIIIO.ImportPhase`, :class:`GSASIIIO.ImportStructFactor`
or :class:`GSASIIIO.ImportPowderData` to implement import of 
a phase, a single crystal or a powder dataset, respectively. 
Module names are used to determine which menu an import routine should
be placed in. 

.. automodule:: G2phase
    :members: 
    :synopsis: Uses previously implemented code: PDB and GSAS .EXP


.. automodule:: G2phase_GPX
    :members: 
    :synopsis: Reads phase information from a GSAS-II project (.gpx) file
      a text file. 

.. automodule:: G2phase_CIF
    :members: 
    :synopsis: Reads phase information from a CIF

.. automodule:: G2pwd_GPX
    :members: 
    :synopsis: Reads powder data from from a GSAS-II project (.gpx) file

.. automodule:: G2pwd_fxye
    :members: 
    :synopsis: Reads powder data in all of the GSAS formats

.. automodule:: G2pwd_xye
    :members: 
    :synopsis: Reads powder data from a Topas format file

.. automodule:: G2pwd_CIF
    :members: 
    :synopsis: Reads powder data from a CIF

.. automodule:: G2sfact
    :members: 
    :synopsis: Reads single crystal data from simple hkl files

.. automodule:: G2sfact_CIF
    :members: 
    :synopsis: Reads single crystal data from CIF files
