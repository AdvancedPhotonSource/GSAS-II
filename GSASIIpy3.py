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
# declare symbol (pi) and functions allowed in expressions
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

def FormatPadValue(val,maxdigits=None):
    '''Format a float to fit in ``maxdigits[0]`` spaces with maxdigits[1] after decimal.

    :param float val: number to be formatted.

    :param list maxdigits: the number of digits & places after decimal to be used for display of the
      number (defaults to [10,2]).

    :returns: a string with exactly maxdigits[0] characters (except under error conditions),
      but last character will always be a space
    '''
    if maxdigits is None:
        digits = [10,2]
    else:
        digits = list(maxdigits)
    fmt = '{:'+str(digits[0])+'}'
    s = fmt.format(FormatValue(val,digits))
    if s[-1] == ' ':
        return s
    else:
        return s+' '
    

def FormatValue(val,maxdigits=None):
    '''Format a float to fit in at most ``maxdigits[0]`` spaces with maxdigits[1] after decimal.
    Note that this code has been hacked from FormatSigFigs and may have unused sections. 

    :param float val: number to be formatted.

    :param list maxdigits: the number of digits & places after decimal to be used for display of the
      number (defaults to [10,2]).

    :returns: a string with <= maxdigits characters (usually).  
    '''
    if 'str' in str(type(val)) and (val == '?' or val == '.'):
        return val        
    if maxdigits is None:
        digits = [10,2]
    else:
        digits = list(maxdigits)
    fmt="{:."+str(digits[1])+"f}"
    string = fmt.format(float(val)).strip() # will standard .f formatting work?
    if len(string) <= digits[0]:
        if ':' in string: # deal with weird bug where a colon pops up in a number when formatting (EPD 7.3.2!)
            string = str(val)
        if digits[1] > 0: # strip off extra zeros on right side
            string = string.rstrip('0')
            if string[-1] == '.': string += "0"
        return string
    if val < 0: digits[0] -= 1 # negative numbers, reserve space for the sign
    decimals = digits[0] - digits[1]
    if abs(val) > 1e99: # for very large numbers, use scientific notation and use all digits
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],digits[0]-6))+"}"
    elif abs(val) > 1e9:
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],digits[0]-5))+"}"
    elif abs(val) < 10**(4-decimals): # make sure at least 4 decimals show
        # this clause is probably no longer needed since the number probably shows as 0.0
        decimals = min(digits[0]-5,digits[1])
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],decimals))+"}"
    elif abs(val) >= 10**(decimals-1): # deal with large numbers in smaller spaces
        decimals = max(0,digits[0]-5)
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],decimals))+"}"
    elif abs(val) < 1: # use f format for small numbers
        # this clause is probably no longer needed since the number probably shows as 0.0
        decimals = min(digits[0]-3,digits[1])
        fmt = "{" + (":{:d}.{:d}f".format(digits[0],decimals))+"}"
    else: # in range where g formatting should do what I want
        # used? 
        decimals = digits[0] - 1
        fmt = "{" + (":{:d}.{:d}g".format(digits[0],decimals))+"}"
    try:
        return fmt.format(float(val)).strip()
    except ValueError as err:
        print 'FormatValue Error with val,maxdigits,fmt=',val,maxdigits,fmt
        return str(val)

def FormatSigFigs(val, maxdigits=10, sigfigs=5, treatAsZero=1e-20):
    '''Format a float to use ``maxdigits`` or fewer digits with ``sigfigs``
    significant digits showing (if room allows).

    :param float val: number to be formatted.

    :param int maxdigits: the number of digits to be used for display of the
       number (defaults to 10).

    :param int sigfigs: the number of significant figures to use, if room allows

    :param float treatAsZero: numbers that are less than this in magnitude
      are treated as zero. Defaults to 1.0e-20, but this can be disabled
      if set to None. 

    :returns: a string with <= maxdigits characters (I hope).  
    '''
    if 'str' in str(type(val)) and (val == '?' or val == '.'):
        return val        
    if treatAsZero is not None:
        if abs(val) < treatAsZero:
            return '0.0'
    # negative numbers, leave room for a sign
    if val < 0: maxdigits -= 1
    if abs(val) < 1e-99 or abs(val) > 9.999e99:
        decimals = min(maxdigits-6,sigfigs)
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits,decimals))+"}" # create format string
    elif abs(val) < 1e-9 or abs(val) > 9.999e9:
        decimals = min(maxdigits-5,sigfigs)
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits,decimals))+"}"
    elif abs(val) < 9.9999999*10**(sigfigs-maxdigits):
        decimals = min(maxdigits-5,sigfigs)
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits,decimals))+"}"
    elif abs(val) >= 10**sigfigs: # deal with large numbers in smaller spaces
        decimals = min(maxdigits-5,sigfigs)
        fmt = "{" + (":{:d}.{:d}g".format(maxdigits,decimals))+"}"
    elif abs(val) < 1: # small numbers, add to decimal places
        decimals = sigfigs - int(np.log10(abs(val)))
        fmt = "{" + (":{:d}.{:d}f".format(maxdigits,decimals))+"}"
    else: # larger numbers, remove decimal places
        decimals = sigfigs - 1 - int(np.log10(abs(val)))
        if decimals <= 0: 
            fmt = "{" + (":{:d}.0f".format(maxdigits))+"}."
        else:
            fmt = "{" + (":{:d}.{:d}f".format(maxdigits,decimals))+"}"
    try:
        return fmt.format(float(val)).strip()
    except ValueError as err:
        print 'FormatValue Error with val,maxdigits, sigfigs, fmt=',val, maxdigits,sigfigs, fmt
        return str(val)

if __name__ == '__main__':
    for i in (1.23456789e-129,1.23456789e129,1.23456789e-99,1.23456789e99,-1.23456789e-99,-1.23456789e99):
        print FormatSigFigs(i),i
    for i in (1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000):
        print FormatSigFigs(1.23456789e-9*i),1.23456789e-9*i
    for i in (1,10,100,1000,10000,100000,1000000,10000000,100000000):
        print FormatSigFigs(1.23456789e9/i),1.23456789e9/i

    print FormatSigFigs(200,10,3)
