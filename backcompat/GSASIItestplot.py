from warnings import warn

from GSASII.GSASIItestplot import *  # noqa: F403

warn(
    "Importing GSASIItestplot as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
