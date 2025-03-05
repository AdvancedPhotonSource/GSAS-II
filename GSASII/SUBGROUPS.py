# -*- coding: utf-8 -*-
from __future__ import division, print_function
import re
import copy
import random as ran
import sys
import os
import numpy as np
import numpy.linalg as nl
from . import GSASIIpath
GSASIIpath.SetBinaryPath()
from . import GSASIIspc as G2spc
from . import GSASIIlattice as G2lat
from . import GSASIIElem as G2elem
from . import GSASIIctrlGUI as G2G
bilbaoSite = 'https://www.cryst.ehu.es/cgi-bin/cryst/programs/'
submagSite = bilbaoSite + 'subgrmag1_general_GSAS.pl?'
pseudosym = 'pseudosym/nph-pseudosym'
timeout=150  # time to wait for www.cryst.ehu.es to respond; 2.5 minutes

def GetNonStdSubgroups(SGData, kvec,star=False,landau=False,maximal=False):
    '''Run Bilbao's SUBGROUPS for a non-standard space group. 
    This requires doing a post to the Bilbao site, which returns all
    subgroups of the entered space group as the text of a web page 
    with a table containing the space group symbol, the 
    transformation matrix and index for each subgroup.

    :params list kvec: propogation vector as a list of nine string fractions or blank
    :params SGData: space group object (see :ref:`Space Group object<SGData_table>`) 

    :returns: (error,text) error: if True no error or False; where 
      text containts a possible web page text
    '''
    print(f'''
    For use of SUBGROUPS, please cite:
{G2G.GetCite('Bilbao: k-SUBGROUPSMAG',wrap=70,indent=5)}''')

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
    page = GSASIIpath.postURL(submagSite,postdict)
    if not page:
        print('connection error - not on internet?')
        return None,None
    page = page.replace('<font style= "text-decoration: overline;">','<font>-')
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
    '''Run Bilbao's k-Subgroupsmag for a non-standard space group. 
    This requires doing a post to the Bilbao site, which returns all
    magnetic subgroups of the entered subgroup as the text of a web page 
    with a table containing the BNS magnetic space group symbol, the 
    transformation matrix and index for each subgroup.

    :params list kvec: propogation vector as a list of three numbers
    :params SGData: space group object (see :ref:`Space Group object<SGData_table>`) 

    :returns: (error,text) error: if True no error or False; where 
      text containts a possible web page text
    '''
    print(f'''
    For use of k-SUBGROUPSMAG, please cite:
{G2G.GetCite('Bilbao: k-SUBGROUPSMAG',wrap=70,indent=5)}''')

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
    page = GSASIIpath.postURL(submagSite,postdict)
    if not page:
        print('connection error - not on internet?')
        return None,None
    page = page.replace('<font style= "text-decoration: overline;">','<font>-')
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
    psSite = bilbaoSite + "pseudosym/nph-pseudolattice"
    cellstr = '+'.join(['{:.5f}'.format(i) for i in cell])
    datastr = "sgr={:}&cell={:}&tol={:}&submit=Show".format(
        str(int(spgNum)),cellstr,str(int(tol)))
    page = GSASIIpath.postURL(psSite,datastr,timeout=timeout)
    if not page:
        print('connection error - not on internet?')
        return None
    page = page.replace('<font style= "text-decoration: overline;">','<font>-')
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

def GetStdSGset(SGData=None, oprList=[]):
    '''Determine the standard setting for a space group from either
    a list of symmetry operators or a space group object using the 
    Bilbao Crystallographic Server utility IDENTIFY GROUP 

    :param list oprList: a list of symmetry operations (example: 
       ['x,y,z', '-x,-y,-z']). Supply either this or SGData, not both. 
    :param SGData: from :func:`GSASIIspc.SpcGroup` 
       Supply either this or oprList, not both.

    :returns: (sgnum, sgnam, xformM, offsetV) where:

      * sgnum is the space group number, 
      * sgnam is the space group short H-M name (as used in GSAS-II)
      * xformM is 3x3 numpy matrix relating the old cell & coordinates to the new
      * offsetV is the vector of offset to be applied to the coordinates

      Note that the new cell is given by G2lat.TransformCell([a,b,...],xformM)
    '''
    import re
    Site = bilbaoSite + 'checkgr.pl'

    if not bool(oprList) ^ bool(SGData):
        raise ValueError('GetStdSGset: Must specify oprList or SGData and not both')
    elif SGData:
        SymOpList,offsetList,symOpList,G2oprList,G2opcodes = G2spc.AllOps(SGData)
        oprList = [x.lower() for x in SymOpList]
    print('Using Bilbao Crystallographic Server utility IDENTIFY GROUP. '+
              'Please cite:\n'+
              G2G.GetCite('Bilbao: k-SUBGROUPSMAG',wrap=70,indent=5))
    postdict = {'tipog':'gesp','generators':'\n'.join(oprList)}
    page = GSASIIpath.postURL(Site,postdict,timeout=timeout)
    if not page:
        print('error:','No response')
        return [None,None,None,None]

    # scrape the HTML output for the new space group # and the xform info
    try:
        sgnum = int(re.search(r'\(No. (\d+)\)',page).group(1))
    except Exception as msg:
        print('error:',msg)
        return [None,None,None,None]
    sgnam = G2spc.spgbyNum[sgnum]
    
    xform = re.split(r'parclose\.png',re.split(r'paropen\.png',page)[1])[0] # pull out section w/Xform matrix
    mat = re.split(r'</pre>',re.split('<pre>',xform)[1])[0].split('\n')
    offsetV = [eval(m.split()[3]) for m in mat]
    xformM = np.array([[float(i) for i in m.split()[:3]] for m in mat])
    return sgnum, sgnam, xformM.T, offsetV

def GetSupergroup(SGnum,dlg=None):
    '''Get supergroups for a space group in a standard setting 
    using the Bilbao Crystallographic Server utility "Minimal 
    Supergroups of Space Groups" (minsup)

    This routine is not fully tested and is not currently implemented.

    :param int SGnum: a space group number (1-230)

    :returns: a list of supergroups where each entry in the list contains 
      [sgnum, sgnam, index, xtype, url, xformlist) where: 

      * sgnum is the space group number, 
      * sgnam is the space group short H-M name (no spacing, not GSAS-II usage)
      * index (int) is the change in asymmetric unit size
      * xtype (char) is the transformation type (t,k) 
      * url is a link to compute the subgroup transformations
      * xformlist is a list containing all the subgroup transformations. 
        xformlist contains [(M,V), (M,V),... ]   where:

        * M is 3x3 numpy matrix relating the old cell & coordinates to the new
        * V is the vector of offset to be applied to the coordinates

        Note that the new cell is given by G2lat.TransformCell([a,b,...],M)
    '''
    import re
    Site = bilbaoSite + 'nph-minsup'
    if dlg: dlg.Update(0,newmsg='Waiting for initial web response')
    out = GSASIIpath.postURL(Site,{'gnum':f'{SGnum:}'},timeout=timeout)
    if not out: return None
        
    if dlg: dlg.Update(1,newmsg='Initial table of supergroups returned')
    lines = out.split('HM symbol')[1].split('/table')[0].split('<tr')
    xforms = []
    for l in lines[1:]:
        ls = l.split('<td>')
        if len(ls) < 8: continue
        spg = re.sub(r'</?[a-zA-Z]+>','',ls[3])
        try:
            spgnum = int(ls[4].split('>')[1].split('<')[0])
        except:
            spgnum = 0
        try:
            index = int(re.sub(r'</?[a-zA-Z]+>','',ls[5]))
        except:
            index = 0
        xformtype = re.sub(r'</?[a-zA-Z]+>','',ls[6])
        click = l.split('click')[1].split('value="')[1].split('"')[0]
        xforms.append([spgnum,spg,index,xformtype,click])

    if dlg: dlg.SetRange(2+len(xforms))
    for i,line in enumerate(xforms):
        click = line[-1]
        print(SGnum,click)
        out1 = GSASIIpath.postURL(Site,{'gnum':SGnum,'show':'show','click':click}
                                ,timeout=timeout)
        if not out1: return None
        #with open(f'/tmp/{click}.html','w') as fp:
        #    fp.write(out1)
        #with open(f'/tmp/last.html','w') as fp:
        #    fp.write(out1)
        if dlg: dlg.Update(2+i,newmsg=f'processing {line[1]}, #{i+1} of {len(xforms)}')
        mvlist = []
        for s in [i.split('</pre>')[0] for i in out1.split('<pre>')[1::2]]:
            mvc = np.array(s.split()).reshape(3,4) 
            mc,vc = mvc[:,:-1],mvc[:,-1]  # separate matrix & vector
            # cast as float if possible
            try:
                m = np.array([float(eval(i)) for i in mc.flatten()]).reshape(3,3)
            except:
                m = mc
            try:
                v = np.array([float(eval(i)) for i in vc.flatten()])
            except:
                v = vc
            mvlist.append((m,v))
        line.append(mvlist)
    return xforms

def applySym(xform,cell):
    '''Determine a unit cell according to a supergroup transformation
    computed with the Bilbao Crystallographic Server utility "Minimal 
    Supergroups of Space Groups" (minsup) utility. 

    The required symmetry is applied to the cell and the cell is 
    scaled so that the unit cell volume is unchanged beyond the 
    scaling in the transformation matrix. 

    This is not used in GSAS-II at present. Some additional
    thought is likely needed to to drop unit cells that are too 
    far from the required lattice symmetry.
    
    :param list xform: a list of transformations from :func:`GetSupergroup`
    :param list cell: the unit cell to be transformmed. 

    :returns: a list of two cells for each transformation matrix in xform. 
      The first cell in each pair is the scaled cell where lattice 
      symmetry has been applied, while the second cell is the direct 
      transform of the input cell. 
    '''
    
    SGData = G2spc.SpcGroup(G2spc.spgbyNum[xform[0]])[1]
    v = G2lat.calc_V(G2lat.cell2A(cell[:6])) # current cell volume
    gmat = G2lat.cell2Gmat(cell[:6])[0] # reciprocal cell tensor
    
    cellsList = []
    for x in reversed(xform[-1]):
        m = x[0]
        vrat = abs(np.linalg.det(m)) # size ratio for new cell
        newg = G2lat.prodMGMT(gmat,m)  # xform the inverse cell & get the A vector
        origA = G2lat.Gmat2A(newg)
        A = copy.copy(origA)
        print('raw',[f"{i:.4f}" for i in G2lat.A2cell(A)])
        # apply lattice symmetry to A
        if SGData['SGLaue'] in ['2/m',]:
            if SGData['SGUniq'] == 'a':
                newA = [A[0],A[1],A[2],0,0,A[5]]
            elif SGData['SGUniq'] == 'b':
                newA = [A[0],A[1],A[2],0,A[4],0]
            else:
                newA = [A[0],A[1],A[2],A[3],0,0]
        elif SGData['SGLaue'] in ['mmm',]:
            newA = [A[0],A[1],A[2],0,0,0]
        elif SGData['SGLaue'] in ['4/m','4/mmm']:
            A[0] = (A[0]+A[1])/2.
            newA = [A[0],A[0],A[2],0,0,0]
        elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
            A[0] = (A[0]+A[1]+A[3])/3.
            newA = [A[0],A[0],A[2],A[0],0,0]
        elif SGData['SGLaue'] in ['3R', '3mR']:
            A[0] = (A[0]+A[1]+A[2])/3.
            A[3] = (A[3]+A[4]+A[5])/3.
            newA = [A[0],A[0],A[0],A[3],A[3],A[3]]
        elif SGData['SGLaue'] in ['m3m','m3']:
            newA = [A[0],A[0],A[0],0,0,0]
        else:
            newA = copy.copy(A)
        # scale the symmetry-enforced A unit cell to match the expected size
        vadj = (G2lat.calc_V(newA)*vrat/v)**(2/3)
        scalA = [i*vadj for i in newA]
        cellsList.append((
            list(G2lat.A2cell(scalA)) + [G2lat.calc_V(scalA)],
            list(G2lat.A2cell(origA)) + [G2lat.calc_V(origA)],
            ))
    return cellsList

def BilbaoSymSearch1(sgnum, phase, maxdelta=2, angtol=None,
                         pagelist=None, keepCell=False):
    '''Perform a search for a supergroup consistent with a phase
    using the Bilbao Pseudosymmetry search (PSEUDO) program, see
    C. Capillas, E.S. Tasci, G. de la Flor, D. Orobengoa, J.M. Perez-Mato 
    and M.I. Aroyo. "A new computer tool at the Bilbao Crystallographic 
    Server to detect and characterize pseudosymmetry". Z. Krist. (2011), 
    226(2), 186-196 DOI:10.1524/zkri.2011.1321.

    The phase must be in a standard setting. 

    :param int sgnum: A space group number (1-230)
    :param dict phase: a GSAS-II phase object (see 
        :ref:`Phase Information<Phase_Information>`). Note that the 
        phase must be in a standard setting (see :func:`GetStdSGset`).
    :param float maxdelta: Allowed distance tolerance in pseudosym search (default 2) 
    :param float angtol: Allowed tolerance for cell angles, used for finding 
      possible unit cells in from triclinic or monoclinic cells, ignored 
      otherwise. Defaults to None, which will cause 5 degrees to be used. 
    :param list pagelist: a list to contain references to the text of web
        pages created by the Bilbao web site. If None (default) the web 
        pages are not saved.
    :param bool keepCell: if False (default) and the cell is monoclinic 
        or triclinic, a search is made for higher symmetry cells. If True,
        the search is made with the current cell. 
    :returns: valsdict,csdict,rowdict,savedcookies where the contents
      will change depending on the space group, but valsdict
      will contain values to be used in the next call to Bilbao and 
      savedcookies will contain a cookie needed for this as well. 
        
      * For monoclinic and triclinic unit cells: csdict will be None and 
        rowdict (rowlist) will be a list containing unit cells of higher 
        symmetry matching the input unit cell to be used for searching
        for supergroups. 
        
      * For higher symmetry unit cells, csdict will be used to select 
        which entries will be used in the next search and rowdict 
        contain possible supergroup settings.
    '''
    
    print(f'''\nUsing the Bilbao Crystallographic Server Pseudosymmetry search (PSEUDO) 
program; Please cite:
{G2G.GetCite('Bilbao: PSEUDO',wrap=70,indent=5)}\n''')

    postdict = {
        "formulae":'',
        "cifile":'',
        "filename":"",
        "what":'minsup',
        "maxik":'1',
        "onlythisik":'1',
        "mysuper2":'211',
        "x1":'1',
        "x2":'0',
        "x3":'0',
        "y1":'0',
        "y2":'1',
        "y3":'0',
        "z1":'0',
        "z2":'0',
        "z3":'1',
        "x4":'0',
        "y4":'0',
        "z4":'0',
        "angtol":'5',
        "submit":'Show',}

    postdict["maxdelta"] = f'{maxdelta:.1f}'
    if sgnum <= 14 and not keepCell: 
        postdict["what"] = 'pseudocell'
        if angtol: # and monoclinic/triclinic
            postdict["angtol"] = f'{angtol:.1f}'
    postdict["stru"] = f'# Space Group ITA number\n{sgnum}\n'
    postdict["stru"] += '# Lattice parameters\n'
    postdict["stru"] += ' '.join([f"{i}" for i in phase['General']['Cell'][1:7]])
    postdict["stru"] += f"\n# number of atoms & atoms\n{len(phase['Atoms'])}\n"
    cx,ct,cs,cia = phase['General']['AtomPtrs']
    for i,atom in enumerate(phase['Atoms'],1):
        el = ''.join([i for i in atom[ct] if i.isalpha()])   # strip to element symbol
        postdict["stru"] += f"{el:4s} {i} - {atom[cx]:.5f} {atom[cx+1]:.5f} {atom[cx+2]:.5f}\n"
#    if GSASIIpath.GetConfigValue('debug'): print(postdict["stru"])
    savedcookies = {}
    page0 = GSASIIpath.postURL(bilbaoSite+pseudosym,postdict,
                             getcookie=savedcookies,timeout=timeout)
    if not page0: return None
    if pagelist is not None:
        pagelist[0] = page0
    if page0 is None: return None
    return scanBilbaoSymSearch1(page0,postdict)+[savedcookies]
        
def scanBilbaoSymSearch1(page0,postdict):
    #open(f'/tmp/pseudosym0.html','w') as fp:
    #    fp.write(page0)
    valsdict = {} # base for all supergroups 
    csdict = {}   # supergroups w/default selection value
    rowdict = {}  # information in table row for each supergroup
    if postdict["what"] == 'pseudocell':
        form = page0.split('Lattice Pseudosymmetry')[1].split('<form')[1].split('</form')[0]
        valsdict["id"] = form.split('name="id"')[1].split('value="')[1].split('"')[0]
        valsdict["tol"] = postdict["angtol"]
        valsdict["idcell"] = form.split('name="idcell"')[1].split('value="')[1].split('"')[0]
        valsdict["what"] = 'minsup-lattice'
        valsdict["formulae"] = postdict["formulae"]
        valsdict["maxdelta"] = postdict["maxdelta"]
        valsdict["submit"] = 'Show'
        rowList = []
        for line in form.split('<tr')[2:]:
            items = line.split('<td')
            # parse table row
            row = [items[2].split('value=')[1].split('"')[1]]  # "lattice" #
            row.append(items[3].split('>')[1].split('<')[0]) # lattice type
            row += items[4].split('<pre>')[1].split('<')[0].split('\n') # ideal & as xformed lattice
            row.append(np.array([
                float(eval(i)) for i in 
                    items[5].split('<pre>')[1].split('<')[0]
                    .replace('[','').replace(']','').split()
                    ]).reshape(3,3))                
            row.append(items[6].split('>')[1].split('<')[0]) # strain
            row.append(items[7].split('>')[1].split('<')[0]) # tol
            rowList.append(row)
        return [valsdict,None,rowList]
            
    # scan output for values to be used next
    for row in page0.split('<input'):
        field = row.split('>')[0]
        if 'name=' not in field: continue
        if 'value=' not in field: continue
        name = value = None
        for item in field.split():
            if '=' not in item: continue
            key,val = item.split('=')
            if key.lower() == 'name':
                name = val.replace('"','')
            if key.lower() == 'value':
                value = val.replace('"','')
            if name == 'cs' and value is not None:
                csdict[value] = "checked" in field
                tr = [i.split('>',1)[1].split('</td')[0] for i in row.split('<td')]
                tr[6] = tr[6].replace('<br>','\n')
                tr = [re.sub(r'\<.*?\>','',i).strip() for i in tr]
                rowdict[value] = tr[1:7] + [not 'invalid' in tr[7]]
                break
            elif name is not None and value is not None:
                if name in valsdict and type(valsdict[name]) is str:
                    valsdict[name] = [valsdict[name],value]
                elif name in valsdict:
                    valsdict[name] += [value]
                else:
                    valsdict[name] = value
                break
    return [valsdict,csdict,rowdict]

def BilbaoLowSymSea1(valsdict,row,savedcookies,pagelist=None):
    '''Using a candidate higher symmetry unit cell from 
    :func:`BilbaoSymSearch1` for monoclinic and triclinic cells, 
    create a list of possible supergroups. 
    Those that match the possible lattice types
    are marked for potential follow-up to see if coordinates can be
    are consistent with that symmetry. 
    
    :returns: latticeList,valsdict,tbl where 
    
     * latticeList: a list of the possible Bravais lattice types
     * valsdict: a dict with values needed for the next web form
     * tbl a list of supergroups with four values per entry, 
         True/False if the lattice type matches, 
         a label with the space group number and the index (sg@ind),
         the space group number and a lattice type (cell & centering)
    '''
    postdict = {}
    postdict.update(valsdict)
    num = row[0]
    if GSASIIpath.GetConfigValue('debug'): print(f"processing cell #{num}")
    postdict['lattice'] = num
    page1 = GSASIIpath.postURL(bilbaoSite+pseudosym,postdict,
                                     usecookie=savedcookies,timeout=timeout)
    if not page1: return None,None,None,None

    lbl = f'cell{num}'
    if pagelist is not None:
        pagelist[lbl] = page1
    if 'Possible lattices:' not in page1: return 3*[None]
    if '<form' not in page1: return 3*[None]
    latticeList = page1.split('Possible lattices:')[1].split('<')[0].strip().split(',')
    form = page1.split('<form')[1].split('</form')[0]
    valsdict = {}
    for key in ('id','polares','idcell','what','formulae','maxdelta','lattice','subcosets'):
        valsdict[key] = form.split(f'name={key}')[1].split('"')[1]
    tbl = []
    for l in form.split('<tr'):
        line = l.split('</tr')[0]
        if 'super_numind' not in line: continue
        items = line.split('<td')
        super_numind = items[2].split('value=')[1].split('>')[0]
        super_numind = super_numind.replace('checked','').strip()
        sg = items[3].split('center">')[1].split('</td')[0]
        sg = sg.replace('<i>','').replace('</i>','')
        sg = sg.replace('<sub>','').replace('</sub>','')
        lattice = items[5].split('>')[1].split('<')[0]
        tbl.append([(lattice in latticeList),super_numind,sg,lattice])
    return lbl,latticeList,valsdict,tbl

def BilbaoLowSymSea2(num,valsdict,row,savedcookies,pagelist=None):
    '''For a selected cell & supergroup from :func:`BilbaoLowSymSea1`, 
    see if the coordinates are consistent with the supergroup
    '''
    postdict = {}
    postdict.update(valsdict)
    postdict['super_numind'] = row[1]
    if GSASIIpath.GetConfigValue('debug'): print(f"processing cell #{num} supergroup {row[1]}")
    page1 = GSASIIpath.postURL(bilbaoSite+pseudosym,postdict,
                                     usecookie=savedcookies,timeout=timeout)
    if page1 is None: return '',None
    lbl = f'cell{num}_{row[1]}'
    if pagelist is not None: 
        pagelist[lbl] = page1
    superdat = page1.split('Idealized structure (supergroup setting')
    if len(superdat) >= 2:
        return lbl,superdat[1].split('<pre>')[1].split('</pre')[0].replace('\t',' ').split('\n')
    return lbl,None

def BilbaoSymSearch2(valsdict,csdict,rowdict,savedcookies,
                         pagelist=None,dlg=None,ophsnam='?'):
    '''Perform a supergroup search from the list of identified supergroups 
    found in BilbaoSymSearch1 (typically where the cell is higher symmetry 
    than monoclinic or triclinic; see :func:`BilbaoLowSymSea1` and
    :func:`BilbaoLowSymSea2` for the low symmetry case.) 
    '''
    structures = {}
    for i,num in enumerate(sorted(csdict)):
        if dlg:
            GoOn = dlg.Update(i+1,newmsg=
                        f'Searching for supergroup(s) consistent with phase {ophsnam}'+
                            f'\nProcessing supergroup #{num}')
            if not GoOn: return {}
        if csdict[num]:
            if GSASIIpath.GetConfigValue('debug'): print(f"processing {num}")
            postdict = {}
            postdict.update(valsdict)
            postdict['cs'] = num
            page1 = GSASIIpath.postURL(bilbaoSite+pseudosym,postdict,
                                     usecookie=savedcookies,timeout=timeout)
            if pagelist is not None:
                pagelist[num] = page1
            if page1 is None:
                structures[num] = "No response, likely web timeout"
            else:
                searchRes = page1.split('Case #')[1].split('/table')[0]         
                superdat = page1.split('Idealized structure (supergroup setting')
                if len(superdat) >= 2:
                    structures[num] = superdat[1].split('<pre>')[1].split('</pre'
                                    )[0].replace('\t',' ').split('\n')
                elif '>tol' in searchRes: 
                    structures[num] = f"coordinates inconsistent with symmetry {rowdict[num][0]}"
    return structures

def find2SearchAgain(pagelist,req='@'):
    '''Scan through web pages from supergroup tests and pull out 
    where coordinates pass the tests to be potentially used to search 
    entries to be used to search for a higher symmetry setting.
    '''
    dicts = {}
    for key,page in pagelist.items():
        if key == 0: continue
        if req and req not in key: continue
        if page is None: continue
        if '<form' not in page: continue
        if 'Continue to search for pseudosymmetry' not in page: continue
        for f in page.split('<form')[1:]:
            if 'Continue to search for pseudosymmetry' not in f: continue
            d = {}
            for field in f.split('<input')[1:]:
                k = field.split('name=')[1].split()[0]
                if 'value="' in field:
                    v = field.split('value="')[1].split('"')[0]
                else:
                    v = field.split('value=')[1].split('>')[0]
                d[k] = v
            k = key
            i = 0
            while k in dicts:
                i += 1
                k = key + '_' + str(i)
            dicts[k] = d
    return dicts

def BilbaoReSymSearch(key,postdict,pagelist=None):
    '''Perform a supergroup search on a result from previously 
    identified supergroup that was found in :func:`find2SearchAgain`
    from the returned web pages. Provides results about the same as from 
    :func:`BilbaoSymSearch1`

    :returns: valsdict,csdict,rowdict,savedcookies where valsdict
      will contain values to be used in the next call to Bilbao and 
      savedcookies will contain a cookie needed for this as well. 

        csdict will be used to select 
        which entries will be used in the next search and rowdict 
        contains possible supergroup settings.

    '''
    savedcookies = {}
    page1 = GSASIIpath.postURL(bilbaoSite+pseudosym,postdict
                             ,getcookie=savedcookies,timeout=timeout)
    if pagelist is not None:
        pagelist[key] = page1
    if page1 is None: return {},{},{},savedcookies
    valsdict,csdict,rowdict = scanBilbaoSymSearch1(page1,postdict)
    return valsdict,csdict,rowdict,savedcookies


def createStdSetting(cifFile,rd):
    '''Use the Bilbao "CIF to Standard Setting" web service to obtain a 
    structure in a standard setting. Then update the reader object with
    the space group, cell and atom positions from this.
    '''
    try:
        import requests # delay this until now, since rarely needed
    except:
        # this import seems to fail with the Anaconda pythonw on
        # macs; it should not!
        print('Warning: failed to import requests. Python config error')
        return None
    cif2std = 'nph-cif2std'
    if not os.path.exists(cifFile):
        print(f'createStdSetting error: file {cifFile} not found')
        return False
    files = {'cifile': open(cifFile,'rb')}
    values = {'strtidy':''}
    print(f'''Submitting structure to Bilbao "CIF to Standard Setting" (strtidy)
web service. Please cite:
{G2G.GetCite('Bilbao: PSEUDOLATTICE',wrap=70,indent=5)}''')
    r0 = requests.post(bilbaoSite+cif2std, files=files, data=values)
    structure = r0.text[r0.text.lower().find('<pre>')+5:r0.text.lower().find('</pre>')].strip()
    spnum,celllist,natom = structure.split('\n')[:3]
    #spgNam = G2spc.spgbyNum[int(spnum)]
    cell = [float(i) for i in celllist.split()]
    # replace cell, space group and atom info with info from Bilbao
    # could try to xfer Uiso (Uij needs xform), but that would be too involved
    rd.Phase['General']['SGData'] = SGData = G2spc.SpcGroup(G2spc.spgbyNum[int(spnum)])[1]
    rd.Phase['General']['Cell'] = [False] + list(cell) + [G2lat.calc_V(G2lat.cell2A(cell))]
    rd.Phase['Atoms'] = []
    for i,line in enumerate(structure.split('\n')[3:]):
        atomlist = ['','Xe','',0.,0.,0.,1.0,'',0.,'I',0.01, 0.,0.,0.,0.,0.,0.,]
        elem,lbl,wyc,x,y,z = line.split()
        elem = elem.rstrip('0123456789-+')
        atomlist[0] = elem + lbl
        if G2elem.CheckElement(elem):
            atomlist[1] = elem
        atomlist[3:6] = [float(i) for i in (x,y,z)]
        atomlist[7],atomlist[8] = G2spc.SytSym(atomlist[3:6],SGData)[:2]
        atomlist[1] = G2elem.FixValence(atomlist[1])

        atomlist.append(ran.randint(0,sys.maxsize)) # add a random Id
        rd.Phase['Atoms'].append(atomlist)
        if i == int(natom)-1: break
    del rd.SymOps['xyz'] # as-read sym ops now obsolete

def test():
    '''This tests that routines in Bilbao Crystallographic Server 
    are accessible and produce output that we can parse. The output 
    is displayed but not checked to see that it agrees with what 
    has been provided previously.
    '''    
    SGData = G2spc.SpcGroup('f d -3 m')[1]
    
    print('test SUBGROUPSMAG')
    results,baseList = GetNonStdSubgroupsmag(SGData,('0','0','0',' ',' ',' ',' ',' ',' ',' '))
    print(results)
    if results:
        for [spgp,bns,mv,gid,altList,supList] in results:
            if gid in baseList:
                print('Space group: %d %s BNS: %s'%(gid,spgp,bns))
                print('MV',mv)
                print('altList:',altList)
                print('superList: ',supList)
                
    print('\n\ntest SUBGROUPS')
    results,baseList = GetNonStdSubgroups(SGData,('1/3','1/3','1/2',' ',' ',' ',' ',' ',' ',' '))
    print(results)
    if results:
        for [spgp,mv,gid,altList,supList] in results:
            if gid in baseList:
                print('Space group: %d %s'%(gid,spgp))
                print('MV',mv)
                print('altList:',altList)
                print('superList: ',supList)
                
    print('\n\ntest Bilbao IDENTIFY GROUP')
    sgnum,sgsym,xmat,xoff = GetStdSGset(G2spc.SpcGroup('R 3 C r')[1])
    if sgnum:
        print(f'Space group R3c (rhomb) transforms to std setting: {sgsym} (#{sgnum})')
        print('  xform matrix',xmat)
        print('  coord offset:',xoff)

if __name__ == '__main__':
    # run self-tests
    selftestquiet = False
    test()
    print ("OK")
