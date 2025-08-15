from warnings import warn

from GSASII.GSASIIindex import *  # noqa: F403

warn(
    "Importing GSASIIindex as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
