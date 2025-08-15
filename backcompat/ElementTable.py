from warnings import warn

from GSASII.data.ElementTable import *  # noqa: F403

warn(
    "Importing GSASIIElem as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
