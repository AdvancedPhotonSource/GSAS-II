# -*- coding: utf-8 -*-
"""
*kSUBGROUPSMAG: Interface to Bilbao k-SUBGROUPSMAG web page*
-------------------------------

Extraction of magnetic space groups for a given space group and a propagation vector
from the k-SUBGROUPSMAG web page on the Bilbao Crystallographic server

"""
########### SVN repository information ###################
# $Date: 2018-07-10 11:41:00 -0500 (Tue, 10 Jul 2018) $
# $Author: vondreele $
# $Revision: 3465 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/kSUBGROUPSMAG.py $
# $Id: kSUBGROUPSMAG.py 3465 2018-07-10 16:41:00Z vondreele $
########### SVN repository information ###################
from __future__ import division, print_function
import requests
try:
    import HTMLParser as HTML
except ImportError: 
    import html.parser as HTML # Python 3
import GSASIIspc as G2spc
import GSASIIpath
GSASIIpath.SetBinaryPath()
submagSite = 'http://www.cryst.ehu.es/cgi-bin/cryst/programs/subgrmag1_general.pl'

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
#        print('tag:',tag,self.beg,self.in_sp)
        if tag == 'i' and self.beg:
            if not self.in_sp:
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
            if 'No.' in data:
                self.spgp += data.split('(')[0]        #pick up trailing number!
                self.in_sp = False
                self.SPGPs.append(self.spgp)
                self.BNSs.append(self.BNS)
#                print('space group:',self.spgp,' BNS: ',self.BNS)
            else:
                if self.in_sub and len(self.spgp) == 1:
                    self.spgp += '_'
                self.spgp += data
                if len(self.spgp) == 3 and '_' in self.spgp:
                    self.spgp += ' '
                    self.BNS = data
                    if 'P_S' in self.spgp:      #fix ambiguous std. symbol for triclinics
                        self.spgp.replace('_S','_c')
                        self.BNS = 'c'
#                print(self.spgp,self.BNS)
        if self.in_pre:
            self.MV = data
    
    def handle_endtag(self, tag):
#        print ('    end tag: %s'%(tag))
        if tag == 'script':
            self.beg = True
        elif tag == 'pre':
            self.in_pre = False
            MV = self.MV.split()
            MV = ['[[%s,%s,%s],[%s,%s,%s],[%s,%s,%s]]'%(MV[0],MV[4],MV[8],MV[1],MV[5],MV[9],MV[2],MV[6],MV[10]),'[%s.,%s.,%s.]'%(MV[3],MV[7],MV[11])]
            self.MVs.append(MV)
#            print('MV:',self.MV)
        elif tag == 'sub':
            self.in_sub = False
            
def GetNonStdSubgroupsmag(SGData, kvec,star=False,landau=False,maximal=False):
    '''Run Bilboa's k-Subgroupsmag for a non-standard space group. 
    This requires doing a post to the Bilboa site, which returns all
    magnetic subgroups of the entered subgroup as the text of a web page 
    with a table containing the BNS magnetic space group symbol, the 
    transformation matrix and index for each subgroup.

    :params list kvec: propogation vector as a list of three numbers
    :params SGData: space group object (see :ref:`Space Group object<SGData_table>`) 

    :returns: (error,text) error: if True no error or False; where 
      text containts a possible web page text
    '''
    print('''
    For use of k-SUBGROUPSMAG, please cite:
      Symmetry-Based Computational Tools for Magnetic Crystallography,
      J.M. Perez-Mato, S.V. Gallego, E.S. Tasci, L. Elcoro, G. de la Flor, and M.I. Aroyo
      Annu. Rev. Mater. Res. 2015. 45,217-48.
      doi: 10.1146/annurev-matsci-070214-021008
    ''')
    starmag = 'no'
    if star:
        starmag = 'yes'
    land = 'no'
    if landau:
        land = 'yes'
    celtodas = 'no'
    limite = 'spgroup'
    if maximal:
        limite = 'maximal'
    postdict = {'centrosymmetry':'0','crystalsystem':'0','landau':land,
               'eleccion':'subgrmag1_k','inicio':'nostandard','celtodas':celtodas,
               'limite':limite,'list':'Submit','listado':'lista','starmagnetica':starmag,
               'pointgroup':'0','polarity':'0','sub':'1.1',
               'super':'','tipog':'gmag','wyckoffstrain':''}
    text,table = G2spc.SGPrint(SGData)
    OpList = G2spc.TextOps(text,table,reverse=True)
#    GenList = G2spc.TextGen(SGData,reverse=True)
    for item in OpList:
        item += '\n'
    sym = ""
    for i in OpList:
        if sym: sym += '\n'
        #if sym: sym += ' ' # use this for testing to generate an error in place of previous
        sym += i.lower()
    postdict['generators'] = sym
    for j in [1,2,3]:
        if kvec[3*j-3] == ' ':
            break
        for i,k in zip(('x','y','z'),kvec[3*j-3:3*j]):
            postdict['km%d%s'%(j,i)] = k
    try:
        r = requests.post(submagSite,postdict)
    except:     #ConnectionError?
        page = ''
        print('connection error - not on internet')
        return None,None
    if r.status_code == 200:
        print('request OK')
        page = r.text
        page = page.replace('<font style= "text-decoration: overline;">','<font>-')
    else:
        page = ''
        print('request failed. Reason=',r.reason)
        return None,None
    r.close()

    p = TableParser()
    p.feed(page)
    nItms = len(p.MVs)
    result = list(zip(p.SPGPs,p.BNSs,p.MVs,range(nItms),nItms*[[],],nItms*['[]',]))
    return result,range(nItms)

def GetNonStdSubgroups(SGData, kvec,star=False,landau=False,maximal=False):
    '''Run Bilboa's SUBGROUPS for a non-standard space group. 
    This requires doing a post to the Bilboa site, which returns all
    subgroups of the entered space group as the text of a web page 
    with a table containing the space group symbol, the 
    transformation matrix and index for each subgroup.

    :params list kvec: propogation vector as a list of nine string fractions or blank
    :params SGData: space group object (see :ref:`Space Group object<SGData_table>`) 

    :returns: (error,text) error: if True no error or False; where 
      text containts a possible web page text
    '''
    print('''
    For use of SUBGROUPS, please cite:
      Symmetry-Based Computational Tools for Magnetic Crystallography,
      J.M. Perez-Mato, S.V. Gallego, E.S. Tasci, L. Elcoro, G. de la Flor, and M.I. Aroyo
      Annu. Rev. Mater. Res. 2015. 45,217-48.
      doi: 10.1146/annurev-matsci-070214-021008
    ''')
    starmag = 'no'
    if star:
        starmag = 'yes'
    land = 'no'
    if landau:
        land = 'yes'
    celtodas = 'no'
    limite = 'spgroup'
    if maximal:
        limite = 'maximal'
    postdict = {'centrosymmetry':'0','crystalsystem':'0','landau':land,
               'eleccion':'subgrmag1_k','inicio':'nostandard','celtodas':celtodas,
               'limite':limite,'list':'Submit','listado':'lista','starmagnetica':starmag,
               'pointgroup':'0','polarity':'0','sub':'1',
               'super':'','tipog':'gesp','wyckoffstrain':''}
    text,table = G2spc.SGPrint(SGData)
    OpList = G2spc.TextOps(text,table,reverse=True)
#    GenList = G2spc.TextGen(SGData,reverse=True)
    for item in OpList:
        item += '\n'
    sym = ""
    for i in OpList:
        if sym: sym += '\n'
        #if sym: sym += ' ' # use this for testing to generate an error in place of previous
        sym += i.lower()
    postdict['generators'] = sym
    for j in [1,2,3]:
        if kvec[3*j-3] == ' ':
            break
        for i,k in zip(('x','y','z'),kvec[3*j-3:3*j]):
            postdict['knm%d%s'%(j,i)] = k
    try:
        r = requests.post(submagSite,postdict)
    except:     #ConnectionError?
        page = ''
        print('connection error - not on internet')
        return None
    if r.status_code == 200:
        print('request OK')
        page = r.text
        page = page.replace('<font style= "text-decoration: overline;">','<font>-')
    else:
        page = ''
        print('request failed. Reason=',r.reason)
        return None
    r.close()

    p = TableParser()
    p.feed(page)
    nItms = len(p.MVs)
    result = list(zip(p.SPGPs,p.MVs,range(nItms),range(nItms),nItms*[0,]))
    return result,range(nItms)

def test():
    SGData = G2spc.SpcGroup('p -3 m 1')[1]
    results,baseList = GetNonStdSubgroupsmag(SGData,('1/3','1/3','1/2',' ',' ',' ',' ',' ',' ',' '))
    if results:
        for [spgp,mv,bns,gid,altList,supList] in results:
            if gid in baseList:
                print('Space group:',spgp, 'BNS:',bns)
                print('MV')
                print(mv)
    results,baseList = GetNonStdSubgroups(SGData,('1/3','1/3','1/2',' ',' ',' ',' ',' ',' ',' '))
    if results:
        for [spgp,mv,gid,altList,supList] in results:
            if gid in baseList:
                print('Space group:',spgp)
                print('MV')
                print(mv)
        

if __name__ == '__main__':
    # run self-tests
    selftestquiet = False
    test()
    print ("OK")
