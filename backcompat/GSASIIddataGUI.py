from warnings import warn

from GSASII.GSASIIddataGUI import *  # noqa: F403

warn(
    "Importing GSASIIddataGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
