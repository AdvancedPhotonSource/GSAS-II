# -*- coding: utf-8 -*-
"""
*MAXMAGN: Interface to Bilbao MAXMAGN web page*
-------------------------------

Extraction of Maximal magnetic space groups for a given space group and a propagation vector
from the MAXMAGN web page on the Bilbao Crystallographic server

"""
########### SVN repository information ###################
# $Date: 2018-07-10 11:41:00 -0500 (Tue, 10 Jul 2018) $
# $Author: vondreele $
# $Revision: 3465 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/MAXMAGN.py $
# $Id: MAXMAGN.py 3465 2018-07-10 16:41:00Z vondreele $
########### SVN repository information ###################
from __future__ import division, print_function
import requests
try:
    import HTMLParser as HTML
except ImportError: 
    import html.parser as HTML # Python 3
    
site='http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-msglist2'

class TableParser(HTML.HTMLParser):
    def __init__(self):
        HTML.HTMLParser.__init__(self)
        self.spgp = ''
        self.in_sp = False
        self.in_pre = False
        self.in_sub = False
        self.MV = ''
        self.BNS = ''
        self.beg = False
        self.SPGPs = []
        self.MVs = []
        self.BNSs = []
    
    def handle_starttag(self, tag, attrs):
#        print('tag:',tag)
        if tag == 'i' and self.beg:
            self.in_sp = True
            self.spgp = ''
            self.BNS = ''
        elif tag == 'pre':
            self.in_pre = True
            self.MV = ''
        elif tag == 'sub':
            self.in_sub = True
    
    def handle_data(self, data):
#        print('*',data)
        if self.in_sp:
            if self.in_sub and len(self.spgp) == 1:
                self.spgp += '_'
            self.spgp += data
            if len(self.spgp) == 3 and '_' in self.spgp:
                self.spgp += ' '
                self.BNS = data
        if self.in_pre:
            self.MV = data
    
    def handle_endtag(self, tag):
#        print ('end tag:',tag)
        if tag == 'h2':
            self.beg = True
        elif tag == 'i' and self.beg:
            self.in_sp = False
            if self.spgp[:3] not in ['Sys','Ten']:
                self.SPGPs.append(self.spgp)
                self.BNSs.append(self.BNS)
#                print('space group:',self.spgp,' BNS: ',self.BNS)
            self.spgp = ''
        elif tag == 'pre':
            self.in_pre = False
            self.MVs.append(self.MV.replace('\n',' '))
#            print('MV:')
#            print(self.MV)
        elif tag == 'sub':
            self.in_sub = False
    
def MAXMAGN(sg,kvec):
    print('''
    For use of MAXMAGN, please cite:
      Symmetry-Based Computational Tools for Magnetic Crystallography,
      J.M. Perez-Mato, S.V. Gallego, E.S. Tasci, L. Elcoro, G. de la Flor, and M.I. Aroyo
      Annu. Rev. Mater. Res. 2015. 45,217–48.
      doi: 10.1146/annurev-matsci-070214-021008
    ''')
    postdict = {'gua':'','gor':'','grha':'','kb':'conventional+dual+%28ITA%29',
                'bnsog':'BNS','list':'Submit'}
    postdict['gnum'] = str(sg)
    for i,k in zip(('x','y','z'),kvec):
        postdict['k'+i] = str(k)
    r = requests.post(site,postdict)
    if r.status_code == 200:
        print('request OK')
        page = r.text
    else:
        page = ''
        print('request failed. Reason=',r.reason)
        return []
    r.close()
    
    p = TableParser()
    p.feed(page)
    return zip(p.SPGPs,p.BNSs,p.MVs)

def test():
    results = MAXMAGN(141,(0.5,0,0))
    if results:
        for spgp,bns,mv in results:
            print('Space group:',spgp, 'BNS:',bns)
            print('MV')
            print(mv)
        

if __name__ == '__main__':
    # run self-tests
    selftestquiet = False
    test()
    print ("OK")
