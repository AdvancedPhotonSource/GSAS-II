# -*- coding: utf-8 -*-
"""
Routines to move information between Python interpreters for code 
development/debugging purposes. This is helpful to move numpy 
data objects to Jupyter notebooks. 

This will need more work to be used on Windows [where it will need 
to use C:/TEMP or the equilvalent, see tempfile.gettempdir()]

Not currently used anywhere in GSAS-II but is easy to insert if needed.
This module is kept separate and small so that it can be copied 
to locations outside of the GSASII project directory where it 
may be imported easily.

Created on Fri Jul 21 10:40:08 2023 by Brian Toby
"""

import sys
import pickle

def dump2tmp(uselocals=True, useglobals=True, usefunctions=False):
    '''
    Places variables from current interpreter into a scratch pickle file 
    that can be read into another Python interpreter. 
    
    Parameters
    ----------
    uselocals : bool, optional
        If True, include objects defined at the local level. 
        The default is True.
    useglobals : bool, optional
        If True, include objects defined at the global level. 
        The default is True.
    usefunctions : bool, optional
        Normally callable functions will not included in the pickled file, 
        but if set to True these functions will be included. Reading the 
        pickle file may then require that the sys.path include references
        to the location of the modules containing these functions.
        The default is False. 

    Returns
    -------
    None.

    '''
    if sys.platform.startswith('win'):
        raise ImportError("Module dmp needs work for Windows")
    fp = open('/tmp/pickledump.pickle', 'wb')
    if uselocals:     # get locals in caller
        import inspect
        frame = inspect.currentframe()
        callerLocals = frame.f_back.f_locals
        for o in callerLocals:
            if o.startswith('_'): continue
            if (not usefunctions) and callable(callerLocals[o]): continue
            try:
                pickle.dump((o,callerLocals[o]), fp)
                print('dumpped',o)
            except:
                print('no dump for ', o, type(callerLocals[o]))
        del frame
    if useglobals:
        for o in globals():
            if o.startswith('_'): continue
            if (not usefunctions) and callable(globals()[o]): continue
            try:
                pickle.dump((o,globals()[o]), fp)
                print('dumpped',o)
            except:
                print('no dump for ', o, type(globals()[o]))
    fp.close()
#dumpStuff()

def undumptmp(setglobals=True):
    '''
    Reads variables saved from another Python interpreter via a 
    scratch pickle file into the current Python interpreter.     

    Parameters
    ----------
    setglobals : bool, optional
        When True variables read will be declared as global. 
        The default is True.

    Returns
    -------
    A dict with all the variables read from the pickle file.

    '''
    import inspect
    frame = inspect.currentframe().f_back
    vars = {}
    fp = open('/tmp/pickledump.pickle', 'rb')
    while True:
        try:
            nam,obj = pickle.load(fp)
            vars[nam] = obj
            if setglobals:
                frame.f_globals[nam] = obj
                #exec(f'global {nam}; {nam} = obj')
#                print('global loaded',nam)
        except EOFError:
            break
        except ModuleNotFoundError:
            print('Ending read due to ModuleNotFoundError')
            break
        except Exception as msg:
            print(nam,'error',msg)
            pass
    fp.close()
    return vars
