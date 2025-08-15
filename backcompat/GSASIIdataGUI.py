from warnings import warn

from GSASII.GSASIIdataGUI import *  # noqa: F403

warn(
    "Importing GSASIIdataGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
