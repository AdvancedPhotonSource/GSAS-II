# -*- coding: utf-8 -*-
"""
*SUBGROUPS: Interface to special GSAS Bilbao SUBGROUPS & k-SUBGROUPSMAG web pages*
-------------------------------

Extraction of  space subgroups for a given space group and a propagation vector
from the GSAS version of SUBGROUPS & k-SUBGROUPSMAG web page on the Bilbao Crystallographic server

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
import numpy as np
import numpy.linalg as nl
import GSASIIspc as G2spc
import GSASIIpath
GSASIIpath.SetBinaryPath()
submagSite = 'http://www.cryst.ehu.es/cgi-bin/cryst/programs/subgrmag1_general_GSAS.pl'

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
        
        
    def getSpGrp(item):
        return item.replace('<i>','').replace('</i>','').replace('<sub>','').replace('</sub>','')
    
    def getMatVec(item):
        return item.replace('{','[').replace('}',']')
               
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
    
    result = page.replace('&','\n')
    result = result.split('\n')
    SPGPs = []
    MVs = []
    baseList = []
    itemList = []
    superList = []
    altList = []
    start = 0
    for line in result:    #work around bug report from Bilbao
        start += 1
        if 'yesz' in line:
            break
    for line in result[start:]:
        if 'GGG' in line:
            lines = line.split('GGG')
            line = lines[0]
            alts = []
            beg = True
            for sline in lines:
                items = sline.split('z')
                gid = int(items[0])
                if beg:
                    baseList.append(gid)
                    beg = False
                alts.append(gid)
                itemList.append(gid)
                superList.append(getMatVec(items[7]))
                SPGPs.append(getSpGrp(items[4]))
                MVs.append([getMatVec(items[5]),getMatVec(items[6])])
            altList.append(alts)
            for sline in lines[1:]:
                altList.append([])
        else:
            items = line.split('z')
            gid = int(items[0])
            altList.append([gid,])
            baseList.append(gid)
            itemList.append(gid)
            superList.append(getMatVec(items[7]))
            SPGPs.append(getSpGrp(items[4]))
            MVs.append([getMatVec(items[5]),getMatVec(items[6])])
    result = list(zip(SPGPs,MVs,itemList,altList,superList))
    return result,baseList

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

    def getSpGrp(item):
        return item.replace('<i>','').replace('</i>','').replace('<sub>','').replace('</sub>','')
    
    def getBNS(item):
        spgrp = getSpGrp(item)
        bns = ''
        sid = item.find('<sub>')
        if sid == 8:
            bns = spgrp[1]
            spgrp = '%s_%s %s'%(spgrp[0],bns,spgrp[2:])
        return spgrp,bns
    
    def getMatVec(item):
        return item.replace('{','[').replace('}',']')
               
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
#    OpList = G2spc.TextGen(SGData,reverse=True)
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

    result = page.replace('&','\n')
    result = result.split('\n')
    start = 0
    for line in result:    #work around bug report from Bilbao
        start += 1
        if 'yesz' in line:
            break
    SPGPs = []
    BNSs = []
    MVs = []
    baseList = []
    itemList = []
    superList = []
    altList = []
    for line in result[start:]:
        if 'GGG' in line:
            lines = line.split('GGG')
            alts = []
            beg = True
            for sline in lines:
                items = sline.split('z')
                gid = int(items[0])
                if beg:
                    baseList.append(gid)
                    beg = False
                alts.append(gid)
                itemList.append(gid)
                superList.append(getMatVec(items[7]))
                spgrp,bns = getBNS(items[4])
                SPGPs.append(spgrp)
                BNSs.append(bns)
                MVs.append([getMatVec(items[5]),getMatVec(items[6])])
            altList.append(alts)
            for sline in lines[1:]:
                altList.append([])
        else:
            items = line.split('z')
            gid = int(items[0])
            altList.append([gid,])
            baseList.append(gid)
            itemList.append(gid)
            superList.append(getMatVec(items[7]))
            spgrp,bns = getBNS(items[4])
            SPGPs.append(spgrp)
            BNSs.append(bns)
            MVs.append([getMatVec(items[5]),getMatVec(items[6])])
    result = list(zip(SPGPs,BNSs,MVs,itemList,altList,superList))
    return result,baseList

def subBilbaoCheckLattice(spgNum,cell,tol=5):
    '''submit a unit cell to  Bilbao PseudoLattice
    '''
    psSite = "http://www.cryst.ehu.es/cgi-bin/cryst/programs/pseudosym/nph-pseudolattice"
    cellstr = '+'.join(['{:.5f}'.format(i) for i in cell])
    datastr = "sgr={:}&cell={:}&tol={:}&submit=Show".format(
        str(int(spgNum)),cellstr,str(int(tol)))
    try:
        r = requests.post(psSite,data=datastr)
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
    return page

def parseBilbaoCheckLattice(page):
    '''find the cell options from the web page returned by Bilbao PseudoLattice
    '''
    cellopts = [i for i in page.split('<tr>') if '<td><pre>' in i]
    found = []
    for c in cellopts:
        cells = c.split("pre")[1].split('<')[0].replace('>','').split('\n') # list of cells, 1st is approx
        try:
            acell = [float(i) for i in cells[0].split()]
            xmatA = [c.split('[')[i].split(']')[0].split() for i in (1,2,3)]
            xmat =  np.array([[eval(i) for i in j] for j in xmatA])
            cellmat = nl.inv(xmat).T
        except:
            print('Error processing cell in',c)
            continue
        found.append((acell,cellmat))
    return found


def test():
    SGData = G2spc.SpcGroup('f d -3 m')[1]
    
    print('test SUBGROUPSMAG')            
    results,baseList = GetNonStdSubgroupsmag(SGData,('0','0','0',' ',' ',' ',' ',' ',' ',' '))
    if results:
        for [spgp,bns,mv,gid,altList,supList] in results:
            if gid in baseList:
                print('Space group: %d %s BNS: %s'%(gid,spgp,bns))
                print('MV',mv)
                print('altList:',altList)
                print('superList: ',supList)
                
    print('test SUBGROUPS')
    results,baseList = GetNonStdSubgroups(SGData,('1/3','1/3','1/2',' ',' ',' ',' ',' ',' ',' '))
    if results:
        for [spgp,mv,gid,altList,supList] in results:
            if gid in baseList:
                print('Space group: %d %s'%(gid,spgp))
                print('MV',mv)
                print('altList:',altList)
                print('superList: ',supList)
                
        

if __name__ == '__main__':
    # run self-tests
    selftestquiet = False
    test()
    print ("OK")
