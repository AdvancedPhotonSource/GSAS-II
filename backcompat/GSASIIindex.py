from warnings import warn

from GSASII.GSASIIindex import *

warn(
    "Importing GSASIIindex as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
