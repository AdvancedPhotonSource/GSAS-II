*GSASIIimage: Image calc module*
================================

*Summary/Contents*
----------------------------

Image calibration, masking & image integration routines.

Note that the GSAS-II image coordinate system is defined as follows:
standing facing the x-ray source (e.g. behind the beamstop),
the synchrotron ring will be to the left (for a left handed
synchrotron – almost all are left handed). That left-right direction
defines X. Y is up and thus Z is toward the source/sample.
The resulting 2D image is then viewed from the sample position
(e.g. between the x-ray source and the detector).
The detector is addressed in units of pixels (or distances using the
pixel size) with the origin as the lower left corner.
The beam center is measured from this point; usually somewhere near
the center of the image, thus both Xc & Yc will be greater than zero
unless the beam center is not on the image. Note that when images are
displayed in image viewers, most software puts the origin in the upper
left corner.

.. contents:: Section Contents 

*GSASIIimage Routines*
------------------------------------


.. automodule:: GSASII.GSASIIimage
    :members: 
    :private-members:
    :special-members:
