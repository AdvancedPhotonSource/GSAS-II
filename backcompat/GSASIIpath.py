from warnings import warn

from GSASII.GSASIIpath import *  # noqa: F403

warn(
    "Importing GSASIIpath as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
