from warnings import warn

from GSASII.GSASIImath import *  # noqa: F403

warn(
    "Importing GSASIImath as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
