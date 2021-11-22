#!/usr/bin/env python
##############################################################################
#
# pdffit2           by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""
Definition of __version__, __date__, __timestamp__, __git_commit__.

Notes
-----
Variable `__gitsha__` is deprecated as of version 1.2.
Use `__git_commit__` instead.
"""

__all__ = ['__date__', '__git_commit__', '__timestamp__', '__version__']

import os.path

from pkg_resources import resource_filename


# obtain version information from the version.cfg file
cp = dict(version='', date='', commit='', timestamp='0')
fcfg = resource_filename(__name__, 'version.cfg')
if not os.path.isfile(fcfg):    # pragma: no cover
    from warnings import warn
    warn('Package metadata not found, execute "./setup.py egg_info".')
    fcfg = os.devnull
with open(fcfg) as fp:
    kwords = [[w.strip() for w in line.split(' = ', 1)]
              for line in fp if line[:1].isalpha() and ' = ' in line]
assert all(w[0] in cp for w in kwords), "received unrecognized keyword"
cp.update(kwords)

__version__ = cp['version']
__date__ = cp['date']
__git_commit__ = cp['commit']
__timestamp__ = int(cp['timestamp'])

# TODO remove deprecated __gitsha__ in version 1.3.
__gitsha__ = __git_commit__

del cp, fcfg, fp, kwords

# End of file
