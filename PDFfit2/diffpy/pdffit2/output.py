#!/usr/bin/env python
##############################################################################
#
# pdffit2           by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Take care of sending engine output to given file-like object.
The output file is stored in local module variable stdout.
"""


# create module variable stdout

from sys import stdout as stdout
# silence pyflakes checker
assert stdout


def redirect_stdout(dst):
    """Redirect PDFfit2 standard output to a file-like object dst.
    The dst value is stored in module variable stdout.
    """
    from diffpy.pdffit2.pdffit2 import redirect_stdout
    redirect_stdout(dst)
    global stdout
    stdout = dst
    return

#  End of file
