from warnings import warn

from GSASII.GSASIIfiles import *  # noqa: F403

warn(
    "Importing GSASIIfiles as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
