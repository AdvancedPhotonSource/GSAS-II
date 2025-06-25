from warnings import warn

from GSASII.GSASIIstrMath import *

warn(
    "Importing GSASIIstrMath as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
