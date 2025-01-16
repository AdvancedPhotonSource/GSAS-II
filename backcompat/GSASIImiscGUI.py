from warnings import warn

from GSASII.GSASIImiscGUI import *

warn(
    "Importing GSASIImiscGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
