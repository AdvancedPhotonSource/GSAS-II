from warnings import warn

from GSASII.GSASIIElemGUI import *

warn(
    "Importing GSASIIElemGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
