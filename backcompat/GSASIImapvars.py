from warnings import warn

from GSASII.GSASIImapvars import *

warn(
    "Importing GSASIImapvars as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
