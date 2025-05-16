from warnings import warn

from GSASII.GSASIIimgGUI import *

warn(
    "Importing GSASIIimgGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
