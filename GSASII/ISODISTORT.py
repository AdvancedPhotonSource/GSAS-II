# -*- coding: utf-8 -*-
from __future__ import division, print_function
import subprocess as subp
import os
import os.path
try:
    import requests
except:
    print('Module requests not installed, access to ISODISTORT not possible')
import copy
from . import GSASIIscriptable as G2sc
from . import GSASIIctrlGUI as G2G
#import tempfile
isouploadsite = 'https://stokes.byu.edu/iso/isodistortuploadfile.php'
isoformsite = 'https://iso.byu.edu/iso/isodistortform.php'

def HandleError(out):
    with open('out.html','wb') as fp:
        fp.write(out.encode("utf-8"))
    url = os.path.realpath('out.html')
    try:
        os.startfile(url)
    except AttributeError:
        try: # should work on MacOS and most linux versions
            subp.call(['open', url])
        except:
            print('Could not open URL')

def UploadCIF(cifname):
       #upload cif file to BYU web site
    ciffile = open(cifname,'rb')
    up1 = {'toProcess':(cifname,ciffile),}
    out1 = requests.post(isouploadsite,files=up1).text
    ciffile.close()

    #retrieve BYU temp file name for cif file

    pos = out1.index('<INPUT')+7
    pos1 = out1.index('VALUE=',pos)+7
    filename = out1[pos1:out1.index('"',pos1)]

    print('ciffile %s uploaded to ISODISTORT to make %s'%(cifname,filename))
    return filename

def MakePhaseCif(data):
    data['pId'] = data.get('pId',0) # needs a pId
    proj = G2sc.G2Project(newgpx='tmp4cif.gpx')
    ph = G2sc.G2Phase(data,data['General']['Name'],proj)
    tempcif = 'ISOin.cif'
    ph.export_CIF(tempcif)
    return tempcif

def GetISODISTORT(Phase):
    '''Run Stokes & Campbell ISODISTORT.
    This requires doing a post to the BYU upload site with a cif file, which returns a BYU local
    copy. This is then sent to the BYU form site with various options, which returns all
    subgroups of the entered space group as the text of a web page with a table containing the space
    group symbol, the transformation matrix and index for each subgroup. Selection of one of these is
    returned to the BYU form site which returns the text of a cif file to be used to create the new phase
    which can apply the distortion mode constraints

    :params dict Phase: GSAS-II phase data

    :returns: radio: dict of possible distortion structures
    :returns: data2: list of str input for next run of isositortform for extracting cif file
    '''

    print(f'''
    For use of ISODISTORT, please cite:

{G2G.GetCite('ISOTROPY, ISODISTORT, ISOCIF...',wrap=60,indent=5)}

{G2G.GetCite('ISODISPLACE',wrap=60,indent=5)}
    ''')

    ISOdata = Phase['ISODISTORT']
    parentcif = ISOdata['ParentCIF']
    childcif = None
    if 'Use this phase' in parentcif:
        parentcif = MakePhaseCif(Phase)
    print(' Run ISODISTORT with %s as parent cif'%parentcif)
    ISOparentcif = UploadCIF(parentcif)

    #submit cif for processing by ISODISTORT

    up2 = {'filename':ISOparentcif,'input':'uploadparentcif'}
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

    if ISOdata['ISOmethod'] == 1:
        data['isosubgroup'] = 'no choice'
        data['isolattice'] = 'no choice'
        data['isoplattice'] = 'no choice'

    elif ISOdata['ISOmethod'] == 3:
        print('method  3 TBD')
        return [],[]
    elif ISOdata['ISOmethod'] == 4:
        childcif = ISOdata['ChildCIF']
        if 'Use this phase' in childcif:
            childcif = MakePhaseCif(Phase)
        print(' Run ISODISTORT with %s as child cif'%childcif)
        ISOchildcif = UploadCIF(childcif)
        data['input'] = 'uploadsubgroupcif'
        data['filename'] = ISOchildcif
        out24 = requests.post(isoformsite,data=data).text
        posB = out24.index('OPTION VALUE=')
        posF = out24[posB:].index('>')+posB
        value = out24[posB+13:posF]
        data['input'] = 'distort'
        data['origintype'] = 'method4'
        data['inputbasis'] = 'list'
        data['basisselect'] = value[1:-1]
        data['chooseorigin'] = False
        data['trynearest'] = True
        data['dmax'] = '1'
        out25 = requests.post(isoformsite,data=data).text
        cifout = GetISOcif(out25,4)
        if cifout is None:
            return None,None
        cifFile = '%s_%s.cif'%(Phase['General']['Name'],'child')
        fl = open(cifFile,'wb')
        fl.write(cifout.encode("utf-8"))
        fl.close()
        return [],cifFile

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
            items = out3[pos:posF].split('=',2)[2].split('>')[0]
            radio['orderparam%d'%num] = items.replace('"','').replace('CHECKED','')
            pos = out3[posF:].index('RADIO')+posF
        except ValueError:
            break

    if parentcif == 'ISOin.cif' or childcif == 'ISOin.cif':
        os.remove('ISOin.cif')

    return radio,data2

def GetISOcif(out4,method):

    try:
        pos = out4.index('<FORM ACTION')
    except ValueError:
        HandleError(out4)
        return None
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
            if 'origintype' in name:
                break
        except ValueError:
            break
    if method == 4:
        try:
            pos = out4[posF:].index('Enter mode')+posF
            pos = out4[pos:].index('<p>')+pos
        except ValueError:
            HandleError(out4)
        while True:
            try:
                posB = out4[pos:].index('name')+pos
                posF = out4[posB:].index('>')+posB
                items = out4[posB:posF].split('=')
                name = items[1].split()[0].replace('"','')
                if name == 'atomicradius':
                    break
                vals = items[2].split()[0].replace('"','')
                data3[name] = vals
                pos = posF
            except ValueError:
                break
        data3['zeromodes'] = False

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
    return out5


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
    #with open('pyout4.html','wb') as fp:
        #fp.write(out4.encode("utf-8"))

    #retrieve data needed for next(last) step

    out5 = GetISOcif(out4,1)
    names = ISOdata['selection'][1].split()
    cifFile = '%s_%s%s%s.cif'%(Phase['General']['Name'],names[1],names[2].replace('*','_'),names[3])
    fl = open(cifFile,'wb')
    fl.write(out5.encode("utf-8"))
    fl.close()

    return cifFile
