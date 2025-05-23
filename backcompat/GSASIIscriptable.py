from warnings import warn

from GSASII.GSASIIscriptable import *

warn(
    "Importing GSASIIscriptable as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
