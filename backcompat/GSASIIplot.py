from warnings import warn

from GSASII.GSASIIplot import *  # noqa: F403

warn(
    "Importing GSASIIplot as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
