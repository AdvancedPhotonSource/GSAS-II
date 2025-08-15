from warnings import warn

from GSASII.GSASIIdata import *  # noqa: F403

warn(
    "Importing GSASIIdata as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
