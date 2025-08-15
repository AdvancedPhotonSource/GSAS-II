from warnings import warn

from GSASII.GSASIIseqGUI import *  # noqa: F403

warn(
    "Importing GSASIIseqGUI as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
