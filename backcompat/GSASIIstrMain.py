from warnings import warn

from GSASII.GSASIIstrMain import *  # noqa: F403

warn(
    "Importing GSASIIstrMain as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
