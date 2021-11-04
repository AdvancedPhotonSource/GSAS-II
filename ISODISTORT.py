# -*- coding: utf-8 -*-
"""
*ISODISTORT: Interface to BYU ISODISTORT web pages*
-------------------------------


"""
########### SVN repository information ###################
# $Date: 2018-07-10 11:41:00 -0500 (Tue, 10 Jul 2018) $
# $Author: vondreele $
# $Revision: 3465 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/kSUBGROUPSMAG.py $
# $Id: kSUBGROUPSMAG.py 3465 2018-07-10 16:41:00Z vondreele $
########### SVN repository information ###################
from __future__ import division, print_function
import subprocess as subp
import os.path
import requests
import copy
isouploadsite = 'https://stokes.byu.edu/iso/isodistortuploadfile.php'
isoformsite = 'https://iso.byu.edu/iso/isodistortform.php'

def HandleError(out):
    open('out.html','wb').write(out.encode("utf-8"))
    url = os.path.realpath('out.html')
    try:
        os.startfile(url)
    except AttributeError:
        try: # should work on MacOS and most linux versions
            subp.call(['open', url])
        except:
            print('Could not open URL')

def GetISODISTORT(Phase,parentcif):
    '''Run Stokes & Campbell ISODISTORT. 
    This requires doing a post to the BYU upload site with a cif file, which returns a BYU local
    copy. This is then sent to the BYU form site with various options, which returns all
    subgroups of the entered space group as the text of a web page with a table containing the space 
    group symbol, the transformation matrix and index for each subgroup. Selection of one of these is 
    returned to the BYU form site which returns the text of a cif file to be used to create the new phase
    which can apply the distortion mode constraints

    :params dict Phase: GSAS-II phase data
    :params str parentcif: parent cif file name - should be local to working directory

    :returns: radio: dict of possible distortion structures
    :returns: data2: list of str input for next run of isositortform for extracting cif file
    '''

    print('''
    For use of ISODISTORT, please cite:
      H. T. Stokes, D. M. Hatch, and B. J. Campbell, ISODISTORT, ISOTROPY Software Suite, iso.byu.edu.
      B. J. Campbell, H. T. Stokes, D. E. Tanner, and D. M. Hatch, "ISODISPLACE: An Internet Tool for Exploring Structural Distortions." 
      J. Appl. Cryst. 39, 607-614 (2006).
    ''')
                   
    
    #upload cif file to BYU web site 
      
    up1 = {'toProcess':(parentcif,open(parentcif,'rb')),}
    out1 = requests.post(isouploadsite,files=up1).text
    
    #retrieve BYU temp file name for cif file    
    
    pos = out1.index('<INPUT')+7
    pos1 = out1.index('VALUE=',pos)+7
    filename = out1[pos1:out1.index('"',pos1)]
    
    print('filename=',filename)
    
    #submit cif for processing by ISODISTORT
    
    up2 = {'filename':filename,'input':'uploadparentcif'}
    out2 = requests.post(isoformsite,up2).text
    
    #recover required info for the distortion search; includes info from cif file (no longer needed)
    
    try:
        pos = out2.index('<p><FORM')
    except ValueError:
        HandleError(out2)
        return [],[]
    data = {}
    while True:
        try:
            posB = out2[pos:].index('INPUT TYPE')+pos
            posF = out2[posB:].index('>')+posB
            items = out2[posB:posF].split('=',3)
            name = items[2].split()[0].replace('"','')
            if 'isosystem' in name:
                break
            vals = items[3].replace('"','')
            data[name] = vals
            pos = posF
        except ValueError:
            break
    #save copy for future use
    data2 = copy.deepcopy(data)
            
    #no limits on space group or lattice
    
    data['isosubgroup'] = 'no choice'
    data['isolattice'] = 'no choice'
    data['isoplattice'] = 'no choice'
    
    #do the distortion search - result is radio button list of choices
    
    out3 = requests.post(isoformsite,data=data).text
    
    #extract the radio button collection
    
    radio = {}
    num = 0
    try:
        pos = out3.index('RADIO')
    except ValueError:
        HandleError(out3)
        return [],[]
        
    while True:
        try:
            posF = out3[pos:].index('<BR>')+pos
            num += 1
            items = out3[pos:posF].split('=',2)
            radio['orderparam%d'%num] = items[2].replace('"','')
            pos = out3[posF:].index('RADIO')+posF
        except ValueError:
            break
    
    return radio,data2

def GetISODISTORTcif(Phase):
    '''Run Stokes & Campbell ISODISTORT. 
    Selection of one of the order parameter disrections is returned to the BYU 
    form site which returns the text of a cif file to be used to create the new phase
    which can apply the distortion mode constraints
     
    :params dict Phase: GSAS-II phase data; contains result of GetISODISTORT above & selection
    
    :returns: CIFfile str: name of cif file created by this in local directory
    '''
    
    ISOdata = Phase['ISODISTORT']
    data2 = ISOdata['rundata']
    #choose one & resubmit
    data2['origintype'] = 'method1'
    data2['orderparam'] = ISOdata['selection'][1]
    data2['input'] = 'distort'
    # for item in data2:
    #     print(item,data2[item])
    out4 = requests.post(isoformsite,data=data2).text
    #print(out4)
    #open('pyout4.html','wb').write(out4.encode("utf-8"))
     
    #retrieve data needed for next(last) step

    try:
        pos = out4.index('<FORM ACTION')
    except ValueError:
        HandleError(out4)
    data3 = {}
    while True:
        try:
            posB = out4[pos:].index('INPUT TYPE')+pos
            posF = out4[posB:].index('>')+posB
            items = out4[posB:posF].split('=',3)
            name = items[2].split()[0].replace('"','')
            if 'subsetting' in name:
                data3[name] = ''
                pos = posF
                continue
            elif 'atomsfile' in name:
                data3[name] = ' '
                pos = posF
                continue
            vals = items[3].replace('"','')
            data3[name] = vals
            pos = posF
            if 'lattparamsub' in name:
                break
        except ValueError:
            break
       
   #request a cif file    
       
    data3['origintype'] = 'structurefile'
    data3['inputvalues'] = 'false'
    data3['atomicradius'] = '0.4'
    data3['bondlength'] = '2.50'
    data3['modeamplitude'] = '1.0'
    data3['strainamplitude'] = '0.1'
    # for item in data3:
    #     print(item,data3[item])
    k = requests.post(isoformsite,data=data3)
    out5 = k.text   #this is output cif!
    #print(out5)
    names = ISOdata['selection'][1].split()
    cifFile = '%s_%s%s%s.cif'%(Phase['General']['Name'],names[1],names[2].replace('*','_'),names[3])
    fl = open(cifFile,'wb')
    fl.write(out5.encode("utf-8"))
    fl.close()
    return cifFile