from warnings import warn

from GSASII.GSASIIrestrGUI import *  # noqa: F403

warn(
    "Importing GSASIIrestrGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
