from warnings import warn

from GSASII.GSASIIobj import *  # noqa: F403

warn(
    "Importing GSASIIobj as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
