from warnings import warn

from GSASII.GSASIIspc import *  # noqa: F403

warn(
    "Importing GSASIIspc as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
