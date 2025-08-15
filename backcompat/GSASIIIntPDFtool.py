from warnings import warn

from GSASII.GSASIIIntPDFtool import *  # noqa: F403

warn(
    "Importing GSASIIIntPDFtool as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
