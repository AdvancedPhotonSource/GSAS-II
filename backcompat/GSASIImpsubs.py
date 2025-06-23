from warnings import warn

from GSASII.GSASIImpsubs import *

warn(
    "Importing GSASIImpsubs as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
