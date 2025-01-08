from warnings import warn

from GSASII.GSASIIimage import *

warn(
    "Importing GSASIIimage as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
