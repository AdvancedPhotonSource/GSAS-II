from warnings import warn

from GSASII.GSASIIdata import *

warn(
    "Importing GSASIIdata as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
