from warnings import warn

from GSASII.GSASIIctrlGUI import *  # noqa: F403

warn(
    "Importing GSASIIctrlGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
