'''
*GSASIIpy3: Python 3.x Routines*
================================

Module to hold python 3-compatible code, to keep it separate from
code that will break with __future__ options.

'''
from __future__ import division
import numpy as np
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
# de
sind = sin = s = lambda x: np.sin(x*np.pi/180.)
cosd = cos = c = lambda x: np.cos(x*np.pi/180.)
tand = tan = t = lambda x: np.tan(x*np.pi/180.)
sqrt = sq = lambda x: np.sqrt(x)
pi = np.pi

def FormulaEval(string):
    '''Evaluates a algebraic formula into a float, if possible. Works
    properly on fractions e.g. 2/3 only with python 3.0+ division.

    Expressions such as 2/3, 3*pi, sin(45)/2, 2*sqrt(2), 2**10 can all
    be evaluated.

    :param str string: Character string containing a Python expression
      to be evaluated.

    :returns: the value for the expression as a float or None if the expression does not
      evaluate to a valid number. 
    
    '''
    try:
        val = float(eval(string))
        if np.isnan(val) or np.isinf(val): return None
    except:
        return None
    return val

def FormatValue(val,maxdigits=[10,2]):
    '''Format a float to fit in ``maxdigits[0]`` spaces with maxdigits[1] after decimal.

    :param float val: number to be formatted.

    :param list maxdigits: the number of digits & places after decimal to be used for display of the
      number (defaults to [10,2]).

    :returns: a string with <= maxdigits characters (I hope).  
    '''
    # does the standard str() conversion fit?
    string = str(val)
    if len(string) <= maxdigits[0]: return string.strip()
    # negative numbers, leave room for a sign
    if val < 0: maxdigits[0] -= 1
    decimals = maxdigits[0] - maxdigits[1]
    if abs(val) < 1e-99 or abs(val) > 1e99:
        decimals = min(maxdigits[0]-6,maxdigits[1])
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits[0],decimals))+"}" # create format string
    elif abs(val) < 1e-9 or abs(val) > 1e9:
        decimals = min(maxdigits[0]-5,maxdigits[1])
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits[0],decimals))+"}"
    elif abs(val) < 10**(4-decimals): # make sure at least 4 decimals show
        decimals = min(maxdigits[0]-5,maxdigits[1])
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits[0],decimals))+"}"
    elif abs(val) >= 10**decimals: # deal with large numbers in smaller spaces
        decimals = min(maxdigits[0]-5,maxdigits[1])
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits[0],decimals))+"}"
    elif abs(val) < 1: # use f format for small numbers
        decimals = min(maxdigits[0]-3,maxdigits[1])
        fmt = "{" + (":{:d}.{:d}f".format(maxdigits[0],decimals))+"}"
    else: # in range where g formatting should do what I want
        decimals = maxdigits[0] - 1
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits[0],decimals))+"}"
    try:
        return fmt.format(val).strip()
    except ValueError as err:
        print 'FormatValue Error with val,maxdigits,fmt=',val,maxdigits,fmt
        return str(val)
