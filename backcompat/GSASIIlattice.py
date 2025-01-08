from warnings import warn

from GSASII.GSASIIlattice import *

warn(
    "Importing GSASIIlattice as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
