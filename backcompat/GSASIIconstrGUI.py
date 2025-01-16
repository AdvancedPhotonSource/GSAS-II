from warnings import warn

from GSASII.GSASIIconstrGUI import *

warn(
    "Importing GSASIIconstrGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
