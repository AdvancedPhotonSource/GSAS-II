from warnings import warn

from GSASII.GSASIIlattice import *  # noqa: F403

warn(
    "Importing GSASIIlattice as a top level module is deprecated, please import "
    + "it as a sub-module of GSASII",
    stacklevel=2,
)
