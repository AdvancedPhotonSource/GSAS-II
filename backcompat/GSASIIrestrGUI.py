from warnings import warn

from GSASII.GSASIIrestrGUI import *

warn(
    "Importing GSASIIrestrGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
