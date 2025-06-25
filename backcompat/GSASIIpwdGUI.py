from warnings import warn

from GSASII.GSASIIpwdGUI import *

warn(
    "Importing GSASIIpwdGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
