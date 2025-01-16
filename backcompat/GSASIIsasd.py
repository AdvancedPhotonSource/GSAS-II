from warnings import warn

from GSASII.GSASIIsasd import *

warn(
    "Importing GSASIIsasd as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
