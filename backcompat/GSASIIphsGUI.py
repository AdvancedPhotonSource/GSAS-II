from warnings import warn

from GSASII.GSASIIphsGUI import *

warn(
    "Importing GSASIIphsGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
