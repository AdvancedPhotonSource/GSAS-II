# -*- coding: utf-8 -*-
'''
*GSASIIsolve - structure solving routines*
==========================================

'''
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import sys
import os.path as ospath
import numpy as np
import cPickle
import time
import math
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIstrIO as G2stIO

def ShowBanner():
    'Print authorship, copyright and citation notice'
    print 80*'*'
    print '    General Structure Analysis System-II Crystal Structure Solution'
    print '              by Robert B. Von Dreele & Brian H. Toby'
    print '                Argonne National Laboratory(C), 2010'
    print ' This product includes software developed by the UChicago Argonne, LLC,' 
    print '            as Operator of Argonne National Laboratory.'
    print 80*'*','\n'
    
def ShowControls(Controls):
    'Print controls information'
    print ' Controls:'
    
def Solve(GPXfile):
    'perform the computation'
    ShowBanner()
    Controls = G2stIO.GetControls(GPXfile)
    ShowControls(Controls)
        
def main():
    'needs doc string'
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

