*GSAS-II Web Modules*
=========================

These modules are used to access external web sites.

---------------------------------------------------------------------
*SUBGROUPS: Interface Bilbao SUBGROUPS & k-SUBGROUPSMAG web pages*
---------------------------------------------------------------------

Extraction of  space subgroups for a given space group and a propagation vector
from the GSAS version of SUBGROUPS & k-SUBGROUPSMAG web page on the Bilbao Crystallographic server. Note that the web pages are special to GSAS-II. 
This uses the :func:`GSASIIpath.postURL` function for web access. 

.. automodule:: GSASII.SUBGROUPS
    :members: 


------------------------------------------------------
*ISODISTORT: Interface to BYU ISODISTORT web pages*
------------------------------------------------------

Uses the BYU ISODISTORT web site to search over all k-points for a structure or to relate a parent and child structure by irreps.
This uses the mod:`requests` package.

.. automodule:: GSASII.ISODISTORT
    :members: 
