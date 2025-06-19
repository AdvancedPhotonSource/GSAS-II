from warnings import warn

from GSASII.GSASIIpwd import *

warn(
    "Importing GSASIIpwd as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
