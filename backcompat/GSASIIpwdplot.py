from warnings import warn

from GSASII.GSASIIpwdplot import *

warn(
    "Importing GSASIIpwdplot as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
