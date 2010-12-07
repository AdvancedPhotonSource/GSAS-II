#GSASIIsolve - structure solving routines
import sys
import os.path as ospath
import numpy as np
import cPickle
import time
import math
import GSASIIpath
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIstruct as G2str

def ShowBanner():
    print 80*'*'
    print '    General Structure Analysis System-II Crystal Structure Solution'
    print '     by Robert B. Von Dreele, Argonne National Laboratory(C), 2010'
    print ' This product includes software developed by the UChicago Argonne, LLC,' 
    print '            as Operator of Argonne National Laboratory.'
    print 80*'*','\n'
    
def ShowControls(Controls):
    print ' Controls:'
    
def Solve(GPXfile):
    ShowBanner()
    Controls = G2str.GetControls(GPXfile)
    ShowControls(Controls)
        
def main():
    arg = sys.argv
    if len(arg) > 1:
        GPXfile = arg[1]
        if not ospath.exists(GPXfile):
            print 'ERROR - ',GPXfile," doesn't exist!"
            exit()
        GPXpath = ospath.dirname(arg[1])
        Solve(GPXfile)
    else:
        print 'ERROR - missing filename'
        exit()
    print "Done"
         
if __name__ == '__main__':
    main()

