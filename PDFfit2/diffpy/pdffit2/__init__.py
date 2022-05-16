#!/usr/bin/env python
##############################################################################
#
# pdffit2           by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""PDFfit2 - real space structure refinement program.
Classes:
    PdfFit
Routines:
    redirect_stdout
"""


from diffpy.pdffit2.version import __version__, __date__
from diffpy.pdffit2.pdffit import PdfFit
from diffpy.pdffit2.output import redirect_stdout
try: 
    from diffpy.pdffit2.pdffit2 import is_element
except ImportError:
    from pdffit2 import is_element  # added for GSAS-II to relocate binary file

# silence pyflakes checker
assert __version__ or True
assert __date__ or True
assert all((PdfFit, redirect_stdout, is_element))

# End of file
