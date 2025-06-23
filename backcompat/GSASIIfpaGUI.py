from warnings import warn

from GSASII.GSASIIfpaGUI import *

warn(
    "Importing GSASIIfpaGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
