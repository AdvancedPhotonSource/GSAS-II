#GSASIIstructure - structure computation routines
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import sys
import numpy as np
import numpy.linalg as nl
import time
import math
import GSASIIpath
import GSASIIIO as G2IO
import GSASIIElem as G2el
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIIpwd as G2pwd
import GSASIImapvars as G2mv
import GSASIImath as G2mth
import scipy.optimize as so

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi


def ShowBanner():
    print 80*'*'
    print '   General Structure Analysis System-II Crystal Structure Refinement'
    print '     by Robert B. Von Dreele, Argonne National Laboratory(C), 2010'
    print ' This product includes software developed by the UChicago Argonne, LLC,' 
    print '            as Operator of Argonne National Laboratory.'
    print 80*'*','\n'

def ShowControls(Controls):
    print ' Least squares controls:'
    print ' Derivative type: ',Controls['deriv type']
    print ' Minimum delta-M/M for convergence: ','%.2g'%(Controls['min dM/M'])
    print ' Initial shift factor: ','%.3f'%(Controls['shift factor'])
    
def GetFFtable(General):
    ''' returns a dictionary of form factor data for atom types found in General
    input:
        General = dictionary of phase info.; includes AtomTypes
    return:
        FFtable = dictionary of form factor data; key is atom type
    '''
    atomTypes = General['AtomTypes']
    FFtable = {}
    for El in atomTypes:
        FFs = G2el.GetFormFactorCoeff(El.split('+')[0].split('-')[0])
        for item in FFs:
            if item['Symbol'] == El.upper():
                FFtable[El] = item
    return FFtable
    
def GetBLtable(General):
    ''' returns a dictionary of neutron scattering length data for atom types & isotopes found in General
    input:
        General = dictionary of phase info.; includes AtomTypes & Isotopes
    return:
        BLtable = dictionary of scattering length data; key is atom type
    '''
    atomTypes = General['AtomTypes']
    BLtable = {}
    isotopes = General['Isotopes']
    isotope = General['Isotope']
    for El in atomTypes:
        BLtable[El] = [isotope[El],isotopes[El][isotope[El]]]
    return BLtable
        
def GetPawleyConstr(SGLaue,PawleyRef,pawleyVary):
    if SGLaue in ['-1','2/m','mmm']:
        return                      #no Pawley symmetry required constraints
    for i,varyI in enumerate(pawleyVary):
        refI = int(varyI.split(':')[-1])
        ih,ik,il = PawleyRef[refI][:3]
        for varyJ in pawleyVary[0:i]:
            refJ = int(varyJ.split(':')[-1])
            jh,jk,jl = PawleyRef[refJ][:3]
            if SGLaue in ['4/m','4/mmm']:
                isum = ih**2+ik**2
                jsum = jh**2+jk**2
                if abs(il) == abs(jl) and isum == jsum:
                    G2mv.StoreEquivalence(varyJ,(varyI,))
            elif SGLaue in ['3R','3mR']:
                isum = ih**2+ik**2+il**2
                jsum = jh**2+jk**2*jl**2
                isum2 = ih*ik+ih*il+ik*il
                jsum2 = jh*jk+jh*jl+jk*jl
                if isum == jsum and isum2 == jsum2:
                    G2mv.StoreEquivalence(varyJ,(varyI,))
            elif SGLaue in ['3','3m1','31m','6/m','6/mmm']:
                isum = ih**2+ik**2+ih*ik
                jsum = jh**2+jk**2+jh*jk
                if abs(il) == abs(jl) and isum == jsum:
                    G2mv.StoreEquivalence(varyJ,(varyI,))
            elif SGLaue in ['m3','m3m']:
                isum = ih**2+ik**2+il**2
                jsum = jh**2+jk**2+jl**2
                if isum == jsum:
                    G2mv.StoreEquivalence(varyJ,(varyI,))
                    
def cellVary(pfx,SGData): 
    if SGData['SGLaue'] in ['-1',]:
        return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3',pfx+'A4',pfx+'A5']
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3']
        elif SGData['SGUniq'] == 'b':
            return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A4']
        else:
            return [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A5']
    elif SGData['SGLaue'] in ['mmm',]:
        return [pfx+'A0',pfx+'A1',pfx+'A2']
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        return [pfx+'A0',pfx+'A2']
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        return [pfx+'A0',pfx+'A2']
    elif SGData['SGLaue'] in ['3R', '3mR']:
        return [pfx+'A0',pfx+'A3']                       
    elif SGData['SGLaue'] in ['m3m','m3']:
        return [pfx+'A0',]
                    
def GetPhaseData(PhaseData,Print=True):
            
    def PrintFFtable(FFtable):
        print '\n X-ray scattering factors:'
        print '   Symbol     fa                                      fb                                      fc'
        print 99*'-'
        for Ename in FFtable:
            ffdata = FFtable[Ename]
            fa = ffdata['fa']
            fb = ffdata['fb']
            print ' %8s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f' %  \
                (Ename.ljust(8),fa[0],fa[1],fa[2],fa[3],fb[0],fb[1],fb[2],fb[3],ffdata['fc'])
                
    def PrintBLtable(BLtable):
        print '\n Neutron scattering factors:'
        print '   Symbol   isotope       mass       b       resonant terms'
        print 99*'-'
        for Ename in BLtable:
            bldata = BLtable[Ename]
            isotope = bldata[0]
            mass = bldata[1][0]
            blen = bldata[1][1]
            bres = []
            if len(bldata[1]) > 2:
                bres = bldata[1][2:]
            line = ' %8s%11s %10.3f %8.3f'%(Ename.ljust(8),isotope.center(11),mass,blen)
            for item in bres:
                line += '%10.5g'%(item)
            print line
                
    def PrintAtoms(General,Atoms):
        print '\n Atoms:'
        line = '   name    type  refine?   x         y         z    '+ \
            '  frac site sym  mult I/A   Uiso     U11     U22     U33     U12     U13     U23'
        if General['Type'] == 'magnetic':
            line += '   Mx     My     Mz'
        elif General['Type'] == 'macromolecular':
            line = ' res no  residue  chain '+line
        print line
        if General['Type'] == 'nuclear':
            print 135*'-'
            for i,at in enumerate(Atoms):
                line = '%7s'%(at[0])+'%7s'%(at[1])+'%7s'%(at[2])+'%10.5f'%(at[3])+'%10.5f'%(at[4])+ \
                    '%10.5f'%(at[5])+'%8.3f'%(at[6])+'%7s'%(at[7])+'%5d'%(at[8])+'%5s'%(at[9])
                if at[9] == 'I':
                    line += '%8.4f'%(at[10])+48*' '
                else:
                    line += 8*' '
                    for j in range(6):
                        line += '%8.4f'%(at[11+j])
                print line
        
    def PrintTexture(textureData):
        topstr = '\n Spherical harmonics texture: Order:' + \
            str(textureData['Order'])
        if textureData['Order']:
            print topstr+' Refine? '+str(textureData['SH Coeff'][0])
        else:
            print topstr
            return
        names = ['omega','chi','phi']
        line = '\n'
        for name in names:
            line += ' SH '+name+':'+'%12.4f'%(textureData['Sample '+name][1])+' Refine? '+str(textureData['Sample '+name][0])
        print line
        print '\n Texture coefficients:'
        ptlbls = ' names :'
        ptstr =  ' values:'
        SHcoeff = textureData['SH Coeff'][1]
        for item in SHcoeff:
            ptlbls += '%12s'%(item)
            ptstr += '%12.4f'%(SHcoeff[item]) 
        print ptlbls
        print ptstr    
        
    if Print: print ' Phases:'
    phaseVary = []
    phaseDict = {}
    phaseConstr = {}
    pawleyLookup = {}
    FFtables = {}                   #scattering factors - xrays
    BLtables = {}                   # neutrons
    Natoms = {}
    AtMults = {}
    AtIA = {}
    shModels = ['cylindrical','none','shear - 2/m','rolling - mmm']
    SamSym = dict(zip(shModels,['0','-1','2/m','mmm']))
    for name in PhaseData:
        General = PhaseData[name]['General']
        pId = PhaseData[name]['pId']
        pfx = str(pId)+'::'
        FFtable = GetFFtable(General)
        BLtable = GetBLtable(General)
        FFtables.update(FFtable)
        BLtables.update(BLtable)
        Atoms = PhaseData[name]['Atoms']
        try:
            PawleyRef = PhaseData[name]['Pawley ref']
        except KeyError:
            PawleyRef = []
        SGData = General['SGData']
        SGtext = G2spc.SGPrint(SGData)
        cell = General['Cell']
        A = G2lat.cell2A(cell[1:7])
        phaseDict.update({pfx+'A0':A[0],pfx+'A1':A[1],pfx+'A2':A[2],pfx+'A3':A[3],pfx+'A4':A[4],pfx+'A5':A[5]})
        if cell[0]:
            phaseVary += cellVary(pfx,SGData)
        Natoms[pfx] = 0
        if Atoms:
            if General['Type'] == 'nuclear':
                Natoms[pfx] = len(Atoms)
                for i,at in enumerate(Atoms):
                    phaseDict.update({pfx+'Atype:'+str(i):at[1],pfx+'Afrac:'+str(i):at[6],pfx+'Amul:'+str(i):at[8],
                        pfx+'Ax:'+str(i):at[3],pfx+'Ay:'+str(i):at[4],pfx+'Az:'+str(i):at[5],
                        pfx+'dAx:'+str(i):0.,pfx+'dAy:'+str(i):0.,pfx+'dAz:'+str(i):0.,         #refined shifts for x,y,z
                        pfx+'AI/A:'+str(i):at[9],})
                    if at[9] == 'I':
                        phaseDict[pfx+'AUiso:'+str(i)] = at[10]
                    else:
                        phaseDict.update({pfx+'AU11:'+str(i):at[11],pfx+'AU22:'+str(i):at[12],pfx+'AU33:'+str(i):at[13],
                            pfx+'AU12:'+str(i):at[14],pfx+'AU13:'+str(i):at[15],pfx+'AU23:'+str(i):at[16]})
                    if 'F' in at[2]:
                        phaseVary.append(pfx+'Afrac:'+str(i))
                    if 'X' in at[2]:
                        xId,xCoef = G2spc.GetCSxinel(at[7])
                        delnames = [pfx+'dAx:'+str(i),pfx+'dAy:'+str(i),pfx+'dAz:'+str(i)]
                        for j in range(3):
                            if xId[j] > 0:                               
                                phaseVary.append(delnames[j])
                                for k in range(j):
                                    if xId[j] == xId[k]:
                                        G2mv.StoreEquivalence(delnames[k],((delnames[j],xCoef[j]),)) 
                    if 'U' in at[2]:
                        if at[9] == 'I':
                            phaseVary.append(pfx+'AUiso:'+str(i))
                        else:
                            uId,uCoef = G2spc.GetCSuinel(at[7])[:2]
                            names = [pfx+'AU11:'+str(i),pfx+'AU22:'+str(i),pfx+'AU33:'+str(i),
                                pfx+'AU12:'+str(i),pfx+'AU13:'+str(i),pfx+'AU23:'+str(i)]
                            for j in range(6):
                                if uId[j] > 0:                               
                                    phaseVary.append(names[j])
                                    for k in range(j):
                                        if uId[j] == uId[k]:
                                            G2mv.StoreEquivalence(names[k],((names[j],uCoef[j]),))
#            elif General['Type'] == 'magnetic':
#            elif General['Type'] == 'macromolecular':

                
            if 'SH Texture' in General:
                textureData = General['SH Texture']
                phaseDict[pfx+'SHmodel'] = SamSym[textureData['Model']]
                phaseDict[pfx+'SHorder'] = textureData['Order']
                for name in ['omega','chi','phi']:
                    phaseDict[pfx+'SH '+name] = textureData['Sample '+name][1]
                    if textureData['Sample '+name][0]:
                        phaseVary.append(pfx+'SH '+name)
                for name in textureData['SH Coeff'][1]:
                    phaseDict[pfx+name] = textureData['SH Coeff'][1][name]
                    if textureData['SH Coeff'][0]:
                        phaseVary.append(pfx+name)
                
            if Print:
                print '\n Phase name: ',General['Name']
                print 135*'-'
                PrintFFtable(FFtable)
                PrintBLtable(BLtable)
                print ''
                for line in SGtext: print line
                PrintAtoms(General,Atoms)
                print '\n Unit cell: a =','%.5f'%(cell[1]),' b =','%.5f'%(cell[2]),' c =','%.5f'%(cell[3]), \
                    ' alpha =','%.3f'%(cell[4]),' beta =','%.3f'%(cell[5]),' gamma =', \
                    '%.3f'%(cell[6]),' volume =','%.3f'%(cell[7]),' Refine?',cell[0]
                if 'SH Texture' in General:
                    PrintTexture(textureData)
                    
        elif PawleyRef:
            pawleyVary = []
            for i,refl in enumerate(PawleyRef):
                phaseDict[pfx+'PWLref:'+str(i)] = refl[6]
                pawleyLookup[pfx+'%d,%d,%d'%(refl[0],refl[1],refl[2])] = i
                if refl[5]:
                    pawleyVary.append(pfx+'PWLref:'+str(i))
            GetPawleyConstr(SGData['SGLaue'],PawleyRef,pawleyVary)      #does G2mv.StoreEquivalence
            phaseVary += pawleyVary
                
    return Natoms,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables
    
def cellFill(pfx,SGData,parmDict,sigDict): 
    if SGData['SGLaue'] in ['-1',]:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
            parmDict[pfx+'A3'],parmDict[pfx+'A4'],parmDict[pfx+'A5']]
        sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
            sigDict[pfx+'A3'],sigDict[pfx+'A4'],sigDict[pfx+'A5']]
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
                parmDict[pfx+'A3'],0,0]
            sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                sigDict[pfx+'A3'],0,0]
        elif SGData['SGUniq'] == 'b':
            A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
                0,parmDict[pfx+'A4'],0]
            sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                0,sigDict[pfx+'A4'],0]
        else:
            A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],
                0,0,parmDict[pfx+'A5']]
            sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],
                0,0,sigDict[pfx+'A5']]
    elif SGData['SGLaue'] in ['mmm',]:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A1'],parmDict[pfx+'A2'],0,0,0]
        sigA = [sigDict[pfx+'A0'],sigDict[pfx+'A1'],sigDict[pfx+'A2'],0,0,0]
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A2'],0,0,0]
        sigA = [sigDict[pfx+'A0'],0,sigDict[pfx+'A2'],0,0,0]
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A2'],
            parmDict[pfx+'A0'],0,0]
        sigA = [sigDict[pfx+'A0'],0,sigDict[pfx+'A2'],0,0,0]
    elif SGData['SGLaue'] in ['3R', '3mR']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A0'],
            parmDict[pfx+'A3'],parmDict[pfx+'A3'],parmDict[pfx+'A3']]
        sigA = [sigDict[pfx+'A0'],0,0,sigDict[pfx+'A3'],0,0]
    elif SGData['SGLaue'] in ['m3m','m3']:
        A = [parmDict[pfx+'A0'],parmDict[pfx+'A0'],parmDict[pfx+'A0'],0,0,0]
        sigA = [sigDict[pfx+'A0'],0,0,0,0,0]
    return A,sigA
        
def getCellEsd(pfx,SGData,A,covData):
    dpr = 180./np.pi
    rVsq = G2lat.calc_rVsq(A)
    G,g = G2lat.A2Gmat(A)       #get recip. & real metric tensors
    cell = np.array(G2lat.Gmat2cell(g))   #real cell
    cellst = np.array(G2lat.Gmat2cell(G)) #recip. cell
    scos = cosd(cellst[3:6])
    ssin = sind(cellst[3:6])
    scot = scos/ssin
    rcos = cosd(cell[3:6])
    rsin = sind(cell[3:6])
    rcot = rcos/rsin
    RMnames = [pfx+'A0',pfx+'A1',pfx+'A2',pfx+'A3',pfx+'A4',pfx+'A5']
    varyList = covData['varyList']
    covMatrix = covData['covMatrix']
    vcov = G2mth.getVCov(RMnames,varyList,covMatrix)
    Ax = np.array(A)
    Ax[3:] /= 2.
    drVdA = np.array([Ax[1]*Ax[2]-Ax[5]**2,Ax[0]*Ax[2]-Ax[4]**2,Ax[0]*Ax[1]-Ax[3]**2,
        Ax[4]*Ax[5]-Ax[2]*Ax[3],Ax[3]*Ax[5]-Ax[1]*Ax[4],Ax[3]*Ax[4]-Ax[0]*Ax[5]])
    srcvlsq = np.inner(drVdA,np.inner(vcov,drVdA.T))
    Vol = 1/np.sqrt(rVsq)
    sigVol = Vol**3*np.sqrt(srcvlsq)/2.
    R123 = Ax[0]*Ax[1]*Ax[2]
    dsasdg = np.zeros((3,6))
    dadg = np.zeros((6,6))
    for i0 in range(3):         #0  1   2
        i1 = (i0+1)%3           #1  2   0
        i2 = (i1+1)%3           #2  0   1
        i3 = 5-i2               #3  5   4
        i4 = 5-i1               #4  3   5
        i5 = 5-i0               #5  4   3
        dsasdg[i0][i1] = 0.5*scot[i0]*scos[i0]/Ax[i1]
        dsasdg[i0][i2] = 0.5*scot[i0]*scos[i0]/Ax[i2]
        dsasdg[i0][i5] = -scot[i0]/np.sqrt(Ax[i1]*Ax[i2])
        denmsq = Ax[i0]*(R123-Ax[i1]*Ax[i4]**2-Ax[i2]*Ax[i3]**2+(Ax[i4]*Ax[i3])**2)
        denom = np.sqrt(denmsq)
        dadg[i5][i0] = -Ax[i5]/denom-rcos[i0]/denmsq*(R123-0.5*Ax[i1]*Ax[i4]**2-0.5*Ax[i2]*Ax[i3]**2)
        dadg[i5][i1] = -0.5*rcos[i0]/denmsq*(Ax[i0]**2*Ax[i2]-Ax[i0]*Ax[i4]**2)
        dadg[i5][i2] = -0.5*rcos[i0]/denmsq*(Ax[i0]**2*Ax[i1]-Ax[i0]*Ax[i3]**2)
        dadg[i5][i3] = Ax[i4]/denom+rcos[i0]/denmsq*(Ax[i0]*Ax[i2]*Ax[i3]-Ax[i3]*Ax[i4]**2)
        dadg[i5][i4] = Ax[i3]/denom+rcos[i0]/denmsq*(Ax[i0]*Ax[i1]*Ax[i4]-Ax[i3]**2*Ax[i4])
        dadg[i5][i5] = -Ax[i0]/denom
    for i0 in range(3):
        i1 = (i0+1)%3
        i2 = (i1+1)%3
        i3 = 5-i2
        for ij in range(6):
            dadg[i0][ij] = cell[i0]*(rcot[i2]*dadg[i3][ij]/rsin[i2]-dsasdg[i1][ij]/ssin[i1])
            if ij == i0:
                dadg[i0][ij] = dadg[i0][ij]-0.5*cell[i0]/Ax[i0]
            dadg[i3][ij] = -dadg[i3][ij]*rsin[2-i0]*dpr
    sigMat = np.inner(dadg,np.inner(vcov,dadg.T))
    var = np.diag(sigMat)
    CS = np.where(var>0.,np.sqrt(var),0.)
    cellSig = [CS[0],CS[1],CS[2],CS[5],CS[4],CS[3],sigVol]  #exchange sig(alp) & sig(gam) to get in right order
    return cellSig            
    
def SetPhaseData(parmDict,sigDict,Phases,covData):
    
    def PrintAtomsAndSig(General,Atoms,atomsSig):
        print '\n Atoms:'
        line = '   name      x         y         z      frac   Uiso     U11     U22     U33     U12     U13     U23'
        if General['Type'] == 'magnetic':
            line += '   Mx     My     Mz'
        elif General['Type'] == 'macromolecular':
            line = ' res no  residue  chain '+line
        print line
        if General['Type'] == 'nuclear':
            print 135*'-'
            fmt = {0:'%7s',1:'%7s',3:'%10.5f',4:'%10.5f',5:'%10.5f',6:'%8.3f',10:'%8.5f',
                11:'%8.5f',12:'%8.5f',13:'%8.5f',14:'%8.5f',15:'%8.5f',16:'%8.5f'}
            noFXsig = {3:[10*' ','%10s'],4:[10*' ','%10s'],5:[10*' ','%10s'],6:[8*' ','%8s']}
            for i,at in enumerate(Atoms):
                name = fmt[0]%(at[0])+fmt[1]%(at[1])+':'
                valstr = ' values:'
                sigstr = ' sig   :'
                for ind in [3,4,5,6]:
                    sigind = str(i)+':'+str(ind)
                    valstr += fmt[ind]%(at[ind])                    
                    if sigind in atomsSig:
                        sigstr += fmt[ind]%(atomsSig[sigind])
                    else:
                        sigstr += noFXsig[ind][1]%(noFXsig[ind][0])
                if at[9] == 'I':
                    valstr += fmt[10]%(at[10])
                    if str(i)+':10' in atomsSig:
                        sigstr += fmt[10]%(atomsSig[str(i)+':10'])
                    else:
                        sigstr += 8*' '
                else:
                    valstr += 8*' '
                    sigstr += 8*' '
                    for ind in [11,12,13,14,15,16]:
                        sigind = str(i)+':'+str(ind)
                        valstr += fmt[ind]%(at[ind])
                        if sigind in atomsSig:                        
                            sigstr += fmt[ind]%(atomsSig[sigind])
                        else:
                            sigstr += 8*' '
                print name
                print valstr
                print sigstr
                
    def PrintSHtextureAndSig(textureData,SHtextureSig):
        print '\n Spherical harmonics texture: Order:' + str(textureData['Order'])
        names = ['omega','chi','phi']
        namstr = '  names :'
        ptstr =  '  values:'
        sigstr = '  esds  :'
        for name in names:
            namstr += '%12s'%(name)
            ptstr += '%12.3f'%(textureData['Sample '+name][1])
            if 'Sample '+name in SHtextureSig:
                sigstr += '%12.3f'%(SHtextureSig['Sample '+name])
            else:
                sigstr += 12*' '
        print namstr
        print ptstr
        print sigstr
        print '\n Texture coefficients:'
        namstr = '  names :'
        ptstr =  '  values:'
        sigstr = '  esds  :'
        SHcoeff = textureData['SH Coeff'][1]
        for name in SHcoeff:
            namstr += '%12s'%(name)
            ptstr += '%12.3f'%(SHcoeff[name])
            if name in SHtextureSig:
                sigstr += '%12.3f'%(SHtextureSig[name])
            else:
                sigstr += 12*' '
        print namstr
        print ptstr
        print sigstr
        
            
    print '\n Phases:'
    for phase in Phases:
        print ' Result for phase: ',phase
        Phase = Phases[phase]
        General = Phase['General']
        SGData = General['SGData']
        Atoms = Phase['Atoms']
        cell = General['Cell']
        pId = Phase['pId']
        pfx = str(pId)+'::'
        if cell[0]:
            A,sigA = cellFill(pfx,SGData,parmDict,sigDict)
            cellSig = getCellEsd(pfx,SGData,A,covData)  #includes sigVol
            print ' Reciprocal metric tensor: '
            ptfmt = "%15.9f"
            names = ['A11','A22','A33','A12','A13','A23']
            namstr = '  names :'
            ptstr =  '  values:'
            sigstr = '  esds  :'
            for name,a,siga in zip(names,A,sigA):
                namstr += '%15s'%(name)
                ptstr += ptfmt%(a)
                if siga:
                    sigstr += ptfmt%(siga)
                else:
                    sigstr += 15*' '
            print namstr
            print ptstr
            print sigstr
            cell[1:7] = G2lat.A2cell(A)
            cell[7] = G2lat.calc_V(A)
            print ' New unit cell:'
            ptfmt = ["%12.6f","%12.6f","%12.6f","%12.4f","%12.4f","%12.4f","%12.3f"]
            names = ['a','b','c','alpha','beta','gamma','Volume']
            namstr = '  names :'
            ptstr =  '  values:'
            sigstr = '  esds  :'
            for name,fmt,a,siga in zip(names,ptfmt,cell[1:8],cellSig):
                namstr += '%12s'%(name)
                ptstr += fmt%(a)
                if siga:
                    sigstr += fmt%(siga)
                else:
                    sigstr += 12*' '
            print namstr
            print ptstr
            print sigstr
            
        if 'Pawley' in Phase['General']['Type']:
            pawleyRef = Phase['Pawley ref']
            for i,refl in enumerate(pawleyRef):
                key = pfx+'PWLref:'+str(i)
                refl[6] = abs(parmDict[key])        #suppress negative Fsq
                if key in sigDict:
                    refl[7] = sigDict[key]
                else:
                    refl[7] = 0
        else:
            atomsSig = {}
            if General['Type'] == 'nuclear':
                for i,at in enumerate(Atoms):
                    names = {3:pfx+'Ax:'+str(i),4:pfx+'Ay:'+str(i),5:pfx+'Az:'+str(i),6:pfx+'Afrac:'+str(i),
                        10:pfx+'AUiso:'+str(i),11:pfx+'AU11:'+str(i),12:pfx+'AU22:'+str(i),13:pfx+'AU33:'+str(i),
                        14:pfx+'AU12:'+str(i),15:pfx+'AU13:'+str(i),16:pfx+'AU23:'+str(i)}
                    for ind in [3,4,5,6]:
                        at[ind] = parmDict[names[ind]]
                        if ind in [3,4,5]:
                            name = names[ind].replace('A','dA')
                        else:
                            name = names[ind]
                        if name in sigDict:
                            atomsSig[str(i)+':'+str(ind)] = sigDict[name]
                    if at[9] == 'I':
                        at[10] = parmDict[names[10]]
                        if names[10] in sigDict:
                            atomsSig[str(i)+':10'] = sigDict[names[10]]
                    else:
                        for ind in [11,12,13,14,15,16]:
                            at[ind] = parmDict[names[ind]]
                            if names[ind] in sigDict:
                                atomsSig[str(i)+':'+str(ind)] = sigDict[names[ind]]
            PrintAtomsAndSig(General,Atoms,atomsSig)
        
        if 'SH Texture' in General:
            textureData = General['SH Texture']    
            if textureData['Order']:
                SHtextureSig = {}
                for name in ['omega','chi','phi']:
                    aname = pfx+'SH '+name
                    textureData['Sample '+name][1] = parmDict[aname]
                    if aname in sigDict:
                        SHtextureSig['Sample '+name] = sigDict[aname]
                for name in textureData['SH Coeff'][1]:
                    aname = pfx+name
                    textureData['SH Coeff'][1][name] = parmDict[aname]
                    if aname in sigDict:
                        SHtextureSig[name] = sigDict[aname]
                PrintSHtextureAndSig(textureData,SHtextureSig)

def GetHistogramPhaseData(Phases,Histograms,Print=True):
    
    def PrintSize(hapData):
        if hapData[0] in ['isotropic','uniaxial']:
            line = '\n Size model    : %9s'%(hapData[0])
            line += ' equatorial:'+'%12.3f'%(hapData[1][0])+' Refine? '+str(hapData[2][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:'+'%12.3f'%(hapData[1][1])+' Refine? '+str(hapData[2][1])
            print line
        else:
            print '\n Size model    : %s'%(hapData[0])
            Snames = ['S11','S22','S33','S12','S13','S23']
            ptlbls = ' names :'
            ptstr =  ' values:'
            varstr = ' refine:'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                varstr += '%12s' % (str(hapData[5][i]))
            print ptlbls
            print ptstr
            print varstr
        
    def PrintMuStrain(hapData,SGData):
        if hapData[0] in ['isotropic','uniaxial']:
            line = '\n Mustrain model: %9s'%(hapData[0])
            line += ' equatorial:'+'%12.1f'%(hapData[1][0])+' Refine? '+str(hapData[2][0])
            if hapData[0] == 'uniaxial':
                line += ' axial:'+'%12.1f'%(hapData[1][1])+' Refine? '+str(hapData[2][1])
            print line
        else:
            print '\n Mustrain model: %s'%(hapData[0])
            Snames = G2spc.MustrainNames(SGData)
            ptlbls = ' names :'
            ptstr =  ' values:'
            varstr = ' refine:'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                varstr += '%12s' % (str(hapData[5][i]))
            print ptlbls
            print ptstr
            print varstr

    def PrintHStrain(hapData,SGData):
        print '\n Hydrostatic/elastic strain: '
        Hsnames = G2spc.HStrainNames(SGData)
        ptlbls = ' names :'
        ptstr =  ' values:'
        varstr = ' refine:'
        for i,name in enumerate(Hsnames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.6f' % (hapData[0][i])
            varstr += '%12s' % (str(hapData[1][i]))
        print ptlbls
        print ptstr
        print varstr

    def PrintSHPO(hapData):
        print '\n Spherical harmonics preferred orientation: Order:' + \
            str(hapData[4])+' Refine? '+str(hapData[2])
        ptlbls = ' names :'
        ptstr =  ' values:'
        for item in hapData[5]:
            ptlbls += '%12s'%(item)
            ptstr += '%12.3f'%(hapData[5][item]) 
        print ptlbls
        print ptstr
    
    hapDict = {}
    hapVary = []
    controlDict = {}
    poType = {}
    poAxes = {}
    spAxes = {}
    spType = {}
    
    for phase in Phases:
        HistoPhase = Phases[phase]['Histograms']
        SGData = Phases[phase]['General']['SGData']
        cell = Phases[phase]['General']['Cell'][1:7]
        A = G2lat.cell2A(cell)
        pId = Phases[phase]['pId']
        histoList = HistoPhase.keys()
        histoList.sort()
        for histogram in histoList:
            try:
                Histogram = Histograms[histogram]
            except KeyError:                        
                #skip if histogram not included e.g. in a sequential refinement
                continue
            hapData = HistoPhase[histogram]
            hId = Histogram['hId']
            limits = Histogram['Limits'][1]
            inst = Histogram['Instrument Parameters']
            inst = dict(zip(inst[3],inst[1]))
            Zero = inst['Zero']
            if 'C' in inst['Type']:
                try:
                    wave = inst['Lam']
                except KeyError:
                    wave = inst['Lam1']
                dmin = wave/(2.0*sind(limits[1]/2.0))
            pfx = str(pId)+':'+str(hId)+':'
            for item in ['Scale','Extinction']:
                hapDict[pfx+item] = hapData[item][0]
                if hapData[item][1]:
                    hapVary.append(pfx+item)
            names = G2spc.HStrainNames(SGData)
            for i,name in enumerate(names):
                hapDict[pfx+name] = hapData['HStrain'][0][i]
                if hapData['HStrain'][1][i]:
                    hapVary.append(pfx+name)
            controlDict[pfx+'poType'] = hapData['Pref.Ori.'][0]
            if hapData['Pref.Ori.'][0] == 'MD':
                hapDict[pfx+'MD'] = hapData['Pref.Ori.'][1]
                controlDict[pfx+'MDAxis'] = hapData['Pref.Ori.'][3]
                if hapData['Pref.Ori.'][2]:
                    hapVary.append(pfx+'MD')
            else:                           #'SH' spherical harmonics
                controlDict[pfx+'SHord'] = hapData['Pref.Ori.'][4]
                controlDict[pfx+'SHncof'] = len(hapData['Pref.Ori.'][5])
                for item in hapData['Pref.Ori.'][5]:
                    hapDict[pfx+item] = hapData['Pref.Ori.'][5][item]
                    if hapData['Pref.Ori.'][2]:
                        hapVary.append(pfx+item)
            for item in ['Mustrain','Size']:
                controlDict[pfx+item+'Type'] = hapData[item][0]
                if hapData[item][0] in ['isotropic','uniaxial']:
                    hapDict[pfx+item+':i'] = hapData[item][1][0]
                    if hapData[item][2][0]:
                        hapVary.append(pfx+item+':i')
                    if hapData[item][0] == 'uniaxial':
                        controlDict[pfx+item+'Axis'] = hapData[item][3]
                        hapDict[pfx+item+':a'] = hapData[item][1][1]
                        if hapData[item][2][1]:
                            hapVary.append(pfx+item+':a')
                else:       #generalized for mustrain or ellipsoidal for size
                    if item == 'Mustrain':
                        names = G2spc.MustrainNames(SGData)
                        pwrs = []
                        for name in names:
                            h,k,l = name[1:]
                            pwrs.append([int(h),int(k),int(l)])
                        controlDict[pfx+'MuPwrs'] = pwrs
                    for i in range(len(hapData[item][4])):
                        sfx = ':'+str(i)
                        hapDict[pfx+item+sfx] = hapData[item][4][i]
                        if hapData[item][5][i]:
                            hapVary.append(pfx+item+sfx)
                            
            if Print: 
                print '\n Phase: ',phase,' in histogram: ',histogram
                print 135*'-'
                print ' Phase fraction  : %10.4f'%(hapData['Scale'][0]),' Refine?',hapData['Scale'][1]
                print ' Extinction coeff: %10.4f'%(hapData['Extinction'][0]),' Refine?',hapData['Extinction'][1]
                if hapData['Pref.Ori.'][0] == 'MD':
                    Ax = hapData['Pref.Ori.'][3]
                    print ' March-Dollase PO: %10.4f'%(hapData['Pref.Ori.'][1]),' Refine?',hapData['Pref.Ori.'][2], \
                        ' Axis: %d %d %d'%(Ax[0],Ax[1],Ax[2])
                else: #'SH' for spherical harmonics
                    PrintSHPO(hapData['Pref.Ori.'])
                PrintSize(hapData['Size'])
                PrintMuStrain(hapData['Mustrain'],SGData)
                PrintHStrain(hapData['HStrain'],SGData)
            HKLd = np.array(G2lat.GenHLaue(dmin,SGData,A))
            refList = []
            for h,k,l,d in HKLd:
                ext,mul,Uniq,phi = G2spc.GenHKLf([h,k,l],SGData)
                if ext:
                    continue
                if 'C' in inst['Type']:
                    pos = 2.0*asind(wave/(2.0*d))+Zero
                    if limits[0] < pos < limits[1]:
                        refList.append([h,k,l,mul,d,pos,0.0,0.0,0.0,0.0,0.0,Uniq,phi,0.0])
                else:
                    raise ValueError 
            Histogram['Reflection Lists'][phase] = refList
    return hapVary,hapDict,controlDict
    
def SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms,Print=True):
    
    def PrintSizeAndSig(hapData,sizeSig):
        line = '\n Size model:     %9s'%(hapData[0])
        refine = False
        if hapData[0] in ['isotropic','uniaxial']:
            line += ' equatorial:%12.3f'%(hapData[1][0])
            if sizeSig[0][0]:
                line += ', sig: %8.3f'%(sizeSig[0][0])
                refine = True
            if hapData[0] == 'uniaxial':
                line += ' axial:%12.3f'%(hapData[1][1])
                if sizeSig[0][1]:
                    refine = True
                    line += ', sig: %8.3f'%(sizeSig[0][1])
            if refine:
                print line
        else:
            Snames = ['S11','S22','S33','S12','S13','S23']
            ptlbls = ' name  :'
            ptstr =  ' value :'
            sigstr = ' sig   :'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                if sizeSig[1][i]:
                    refine = True
                    sigstr += '%12.6f' % (sizeSig[1][i])
                else:
                    sigstr += 12*' '
            if refine:
                print line
                print ptlbls
                print ptstr
                print sigstr
        
    def PrintMuStrainAndSig(hapData,mustrainSig,SGData):
        line = '\n Mustrain model: %9s'%(hapData[0])
        refine = False
        if hapData[0] in ['isotropic','uniaxial']:
            line += ' equatorial:%12.1f'%(hapData[1][0])
            if mustrainSig[0][0]:
                line += ', sig: %8.1f'%(mustrainSig[0][0])
                refine = True
            if hapData[0] == 'uniaxial':
                line += ' axial:%12.1f'%(hapData[1][1])
                if mustrainSig[0][1]:
                     line += ', sig: %8.1f'%(mustrainSig[0][1])
            if refine:
                print line
        else:
            Snames = G2spc.MustrainNames(SGData)
            ptlbls = ' name  :'
            ptstr =  ' value :'
            sigstr = ' sig   :'
            for i,name in enumerate(Snames):
                ptlbls += '%12s' % (name)
                ptstr += '%12.6f' % (hapData[4][i])
                if mustrainSig[1][i]:
                    refine = True
                    sigstr += '%12.6f' % (mustrainSig[1][i])
                else:
                    sigstr += 12*' '
            if refine:
                print line
                print ptlbls
                print ptstr
                print sigstr
            
    def PrintHStrainAndSig(hapData,strainSig,SGData):
        Hsnames = G2spc.HStrainNames(SGData)
        ptlbls = ' name  :'
        ptstr =  ' value :'
        sigstr = ' sig   :'
        refine = False
        for i,name in enumerate(Hsnames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.6g' % (hapData[0][i])
            if name in strainSig:
                refine = True
                sigstr += '%12.6g' % (strainSig[name])
            else:
                sigstr += 12*' '
        if refine:
            print '\n Hydrostatic/elastic strain: '
            print ptlbls
            print ptstr
            print sigstr
        
    def PrintSHPOAndSig(hapData,POsig):
        print '\n Spherical harmonics preferred orientation: Order:'+str(hapData[4])
        ptlbls = ' names :'
        ptstr =  ' values:'
        sigstr = ' sig   :'
        for item in hapData[5]:
            ptlbls += '%12s'%(item)
            ptstr += '%12.3f'%(hapData[5][item])
            if item in POsig:
                sigstr += '%12.3f'%(POsig[item])
            else:
                sigstr += 12*' ' 
        print ptlbls
        print ptstr
        print sigstr
    
    for phase in Phases:
        HistoPhase = Phases[phase]['Histograms']
        SGData = Phases[phase]['General']['SGData']
        pId = Phases[phase]['pId']
        histoList = HistoPhase.keys()
        histoList.sort()
        for histogram in histoList:
            try:
                Histogram = Histograms[histogram]
            except KeyError:                        
                #skip if histogram not included e.g. in a sequential refinement
                continue
            print '\n Phase: ',phase,' in histogram: ',histogram
            print 130*'-'
            hapData = HistoPhase[histogram]
            hId = Histogram['hId']
            pfx = str(pId)+':'+str(hId)+':'
            print ' Final refinement RF, RF^2 = %.2f%%, %.2f%% on %d reflections'   \
                %(Histogram[pfx+'Rf'],Histogram[pfx+'Rf^2'],Histogram[pfx+'Nref'])
            
            PhFrExtPOSig = {}
            for item in ['Scale','Extinction']:
                hapData[item][0] = parmDict[pfx+item]
                if pfx+item in sigDict:
                    PhFrExtPOSig[item] = sigDict[pfx+item]
            if hapData['Pref.Ori.'][0] == 'MD':
                hapData['Pref.Ori.'][1] = parmDict[pfx+'MD']
                if pfx+'MD' in sigDict:
                    PhFrExtPOSig['MD'] = sigDict[pfx+'MD']
            else:                           #'SH' spherical harmonics
                for item in hapData['Pref.Ori.'][5]:
                    hapData['Pref.Ori.'][5][item] = parmDict[pfx+item]
                    if pfx+item in sigDict:
                        PhFrExtPOSig[item] = sigDict[pfx+item]
            if Print:
                if 'Scale' in PhFrExtPOSig:
                    print ' Phase fraction  : %10.4f, sig %10.4f'%(hapData['Scale'][0],PhFrExtPOSig['Scale'])
                if 'Extinction' in PhFrExtPOSig:
                    print ' Extinction coeff: %10.4f, sig %10.4f'%(hapData['Extinction'][0],PhFrExtPOSig['Extinction'])
                if hapData['Pref.Ori.'][0] == 'MD':
                    if 'MD' in PhFrExtPOSig:
                        print ' March-Dollase PO: %10.4f, sig %10.4f'%(hapData['Pref.Ori.'][1],PhFrExtPOSig['MD'])
                else:
                    PrintSHPOAndSig(hapData['Pref.Ori.'],PhFrExtPOSig)
            SizeMuStrSig = {'Mustrain':[[0,0],[0 for i in range(len(hapData['Mustrain'][4]))]],
                'Size':[[0,0],[0 for i in range(len(hapData['Size'][4]))]],
                'HStrain':{}}                  
            for item in ['Mustrain','Size']:
                if hapData[item][0] in ['isotropic','uniaxial']:                    
                    hapData[item][1][0] = parmDict[pfx+item+':i']
                    if item == 'Size':
                        hapData[item][1][0] = min(10.,max(0.001,hapData[item][1][0]))
                    if pfx+item+':i' in sigDict: 
                        SizeMuStrSig[item][0][0] = sigDict[pfx+item+':i']
                    if hapData[item][0] == 'uniaxial':
                        hapData[item][1][1] = parmDict[pfx+item+':a']
                        if item == 'Size':
                            hapData[item][1][1] = min(10.,max(0.001,hapData[item][1][1]))                        
                        if pfx+item+':a' in sigDict:
                            SizeMuStrSig[item][0][1] = sigDict[pfx+item+':a']
                else:       #generalized for mustrain or ellipsoidal for size
                    for i in range(len(hapData[item][4])):
                        sfx = ':'+str(i)
                        hapData[item][4][i] = parmDict[pfx+item+sfx]
                        if pfx+item+sfx in sigDict:
                            SizeMuStrSig[item][1][i] = sigDict[pfx+item+sfx]
            names = G2spc.HStrainNames(SGData)
            for i,name in enumerate(names):
                hapData['HStrain'][0][i] = parmDict[pfx+name]
                if pfx+name in sigDict:
                    SizeMuStrSig['HStrain'][name] = sigDict[pfx+name]
            if Print:
                PrintSizeAndSig(hapData['Size'],SizeMuStrSig['Size'])
                PrintMuStrainAndSig(hapData['Mustrain'],SizeMuStrSig['Mustrain'],SGData)
                PrintHStrainAndSig(hapData['HStrain'],SizeMuStrSig['HStrain'],SGData)
    
def GetHistogramData(Histograms,Print=True):
    
    def GetBackgroundParms(hId,Background):
        Back = Background[0]
        Debye = Background[1]
        bakType,bakFlag = Back[:2]
        backVals = Back[3:]
        backNames = [':'+str(hId)+':Back:'+str(i) for i in range(len(backVals))]
        backDict = dict(zip(backNames,backVals))
        backVary = []
        if bakFlag:
            backVary = backNames
        backDict[':'+str(hId)+':nDebye'] = Debye['nDebye']
        debyeDict = {}
        debyeList = []
        for i in range(Debye['nDebye']):
            debyeNames = [':'+str(hId)+':DebyeA:'+str(i),':'+str(hId)+':DebyeR:'+str(i),':'+str(hId)+':DebyeU:'+str(i)]
            debyeDict.update(dict(zip(debyeNames,Debye['debyeTerms'][i][::2])))
            debyeList += zip(debyeNames,Debye['debyeTerms'][i][1::2])
        debyeVary = []
        for item in debyeList:
            if item[1]:
                debyeVary.append(item[0])
        backDict.update(debyeDict)
        backVary += debyeVary    
        return bakType,backDict,backVary            
        
    def GetInstParms(hId,Inst):
        insVals,insFlags,insNames = Inst[1:4]
        dataType = insVals[0]
        instDict = {}
        insVary = []
        pfx = ':'+str(hId)+':'
        for i,flag in enumerate(insFlags):
            insName = pfx+insNames[i]
            instDict[insName] = insVals[i]
            if flag:
                insVary.append(insName)
        instDict[pfx+'X'] = max(instDict[pfx+'X'],0.001)
        instDict[pfx+'Y'] = max(instDict[pfx+'Y'],0.001)
        instDict[pfx+'SH/L'] = max(instDict[pfx+'SH/L'],0.0005)
        return dataType,instDict,insVary
        
    def GetSampleParms(hId,Sample):
        sampVary = []
        hfx = ':'+str(hId)+':'        
        sampDict = {hfx+'Gonio. radius':Sample['Gonio. radius'],hfx+'Omega':Sample['Omega'],
            hfx+'Chi':Sample['Chi'],hfx+'Phi':Sample['Phi']}
        Type = Sample['Type']
        if 'Bragg' in Type:             #Bragg-Brentano
            for item in ['Scale','Shift','Transparency']:       #surface roughness?, diffuse scattering?
                sampDict[hfx+item] = Sample[item][0]
                if Sample[item][1]:
                    sampVary.append(hfx+item)
        elif 'Debye' in Type:        #Debye-Scherrer
            for item in ['Scale','Absorption','DisplaceX','DisplaceY']:
                sampDict[hfx+item] = Sample[item][0]
                if Sample[item][1]:
                    sampVary.append(hfx+item)
        return Type,sampDict,sampVary
        
    def PrintBackground(Background):
        Back = Background[0]
        Debye = Background[1]
        print '\n Background function: ',Back[0],' Refine?',bool(Back[1])
        line = ' Coefficients: '
        for i,back in enumerate(Back[3:]):
            line += '%10.3f'%(back)
            if i and not i%10:
                line += '\n'+15*' '
        print line
        if Debye['nDebye']:
            print '\n Debye diffuse scattering coefficients'
            parms = ['DebyeA','DebyeR','DebyeU']
            line = ' names :'
            for parm in parms:
                line += '%16s'%(parm)
            print line
            for j,term in enumerate(Debye['debyeTerms']):
                line = ' term'+'%2d'%(j)+':'
                for i in range(3):
                    line += '%10.4g %5s'%(term[2*i],bool(term[2*i+1]))                    
                print line
        
    def PrintInstParms(Inst):
        print '\n Instrument Parameters:'
        ptlbls = ' name  :'
        ptstr =  ' value :'
        varstr = ' refine:'
        instNames = Inst[3][1:]
        for i,name in enumerate(instNames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.6f' % (Inst[1][i+1])
            if name in ['Lam1','Lam2','Azimuth']:
                varstr += 12*' '
            else:
                varstr += '%12s' % (str(bool(Inst[2][i+1])))
        print ptlbls
        print ptstr
        print varstr
        
    def PrintSampleParms(Sample):
        print '\n Sample Parameters:'
        print ' Goniometer omega = %.2f, chi = %.2f, phi = %.2f'% \
            (Sample['Omega'],Sample['Chi'],Sample['Phi'])
        ptlbls = ' name  :'
        ptstr =  ' value :'
        varstr = ' refine:'
        if 'Bragg' in Sample['Type']:
            for item in ['Scale','Shift','Transparency']:
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                varstr += '%14s'%(str(bool(Sample[item][1])))
            
        elif 'Debye' in Type:        #Debye-Scherrer
            for item in ['Scale','Absorption','DisplaceX','DisplaceY']:
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                varstr += '%14s'%(str(bool(Sample[item][1])))

        print ptlbls
        print ptstr
        print varstr
        

    histDict = {}
    histVary = []
    controlDict = {}
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        Histogram = Histograms[histogram]
        hId = Histogram['hId']
        pfx = ':'+str(hId)+':'
        controlDict[pfx+'Limits'] = Histogram['Limits'][1]
        
        Background = Histogram['Background']
        Type,bakDict,bakVary = GetBackgroundParms(hId,Background)
        controlDict[pfx+'bakType'] = Type
        histDict.update(bakDict)
        histVary += bakVary
        
        Inst = Histogram['Instrument Parameters']
        Type,instDict,insVary = GetInstParms(hId,Inst)
        controlDict[pfx+'histType'] = Type
        if pfx+'Lam1' in instDict:
            controlDict[pfx+'keV'] = 12.397639/instDict[pfx+'Lam1']
        else:
            controlDict[pfx+'keV'] = 12.397639/instDict[pfx+'Lam']            
        histDict.update(instDict)
        histVary += insVary
        
        Sample = Histogram['Sample Parameters']
        Type,sampDict,sampVary = GetSampleParms(hId,Sample)
        controlDict[pfx+'instType'] = Type
        histDict.update(sampDict)
        histVary += sampVary

        if Print: 
            print '\n Histogram: ',histogram,' histogram Id: ',hId
            print 135*'-'
            Units = {'C':' deg','T':' msec'}
            units = Units[controlDict[pfx+'histType'][2]]
            Limits = controlDict[pfx+'Limits']
            print ' Instrument type: ',Sample['Type']
            print ' Histogram limits: %8.2f%s to %8.2f%s'%(Limits[0],units,Limits[1],units)     
            PrintSampleParms(Sample)
            PrintInstParms(Inst)
            PrintBackground(Background)
        
    return histVary,histDict,controlDict
    
def SetHistogramData(parmDict,sigDict,Histograms,Print=True):
    
    def SetBackgroundParms(pfx,Background,parmDict,sigDict):
        Back = Background[0]
        Debye = Background[1]
        lenBack = len(Back[3:])
        backSig = [0 for i in range(lenBack+3*Debye['nDebye'])]
        for i in range(lenBack):
            Back[3+i] = parmDict[pfx+'Back:'+str(i)]
            if pfx+'Back:'+str(i) in sigDict:
                backSig[i] = sigDict[pfx+'Back:'+str(i)]
        if Debye['nDebye']:
            for i in range(Debye['nDebye']):
                names = [pfx+'DebyeA:'+str(i),pfx+'DebyeR:'+str(i),pfx+'DebyeU:'+str(i)]
                for j,name in enumerate(names):
                    Debye['debyeTerms'][i][2*j] = parmDict[name]
                    if name in sigDict:
                        backSig[lenBack+3*i+j] = sigDict[name]            
        return backSig
        
    def SetInstParms(pfx,Inst,parmDict,sigDict):
        insVals,insFlags,insNames = Inst[1:4]
        instSig = [0 for i in range(len(insVals))]
        for i,flag in enumerate(insFlags):
            insName = pfx+insNames[i]
            insVals[i] = parmDict[insName]
            if insName in sigDict:
                instSig[i] = sigDict[insName]
        return instSig
        
    def SetSampleParms(pfx,Sample,parmDict,sigDict):
        if 'Bragg' in Sample['Type']:             #Bragg-Brentano
            sampSig = [0 for i in range(3)]
            for i,item in enumerate(['Scale','Shift','Transparency']):       #surface roughness?, diffuse scattering?
                Sample[item][0] = parmDict[pfx+item]
                if pfx+item in sigDict:
                    sampSig[i] = sigDict[pfx+item]
        elif 'Debye' in Sample['Type']:        #Debye-Scherrer
            sampSig = [0 for i in range(4)]
            for i,item in enumerate(['Scale','Absorption','DisplaceX','DisplaceY']):
                Sample[item][0] = parmDict[pfx+item]
                if pfx+item in sigDict:
                    sampSig[i] = sigDict[pfx+item]
        return sampSig
        
    def PrintBackgroundSig(Background,backSig):
        Back = Background[0]
        Debye = Background[1]
        lenBack = len(Back[3:])
        valstr = ' value : '
        sigstr = ' sig   : '
        refine = False
        for i,back in enumerate(Back[3:]):
            valstr += '%10.4g'%(back)
            if Back[1]:
                refine = True
                sigstr += '%10.4g'%(backSig[i])
            else:
                sigstr += 10*' '
        if refine:
            print '\n Background function: ',Back[0]
            print valstr
            print sigstr 
        if Debye['nDebye']:
            ifAny = False
            ptfmt = "%12.5f"
            names =  ' names :'
            ptstr =  ' values:'
            sigstr = ' esds  :'
            for item in sigDict:
                if 'Debye' in item:
                    ifAny = True
                    names += '%12s'%(item)
                    ptstr += ptfmt%(parmDict[item])
                    sigstr += ptfmt%(sigDict[item])
            if ifAny:
                print '\n Debye diffuse scattering coefficients'
                print names
                print ptstr
                print sigstr
        
    def PrintInstParmsSig(Inst,instSig):
        ptlbls = ' names :'
        ptstr =  ' value :'
        sigstr = ' sig   :'
        instNames = Inst[3][1:]
        refine = False
        for i,name in enumerate(instNames):
            ptlbls += '%12s' % (name)
            ptstr += '%12.6f' % (Inst[1][i+1])
            if instSig[i+1]:
                refine = True
                sigstr += '%12.6f' % (instSig[i+1])
            else:
                sigstr += 12*' '
        if refine:
            print '\n Instrument Parameters:'
            print ptlbls
            print ptstr
            print sigstr
        
    def PrintSampleParmsSig(Sample,sampleSig):
        ptlbls = ' names :'
        ptstr =  ' values:'
        sigstr = ' sig   :'
        refine = False
        if 'Bragg' in Sample['Type']:
            for i,item in enumerate(['Scale','Shift','Transparency']):
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                if sampleSig[i]:
                    refine = True
                    sigstr += '%14.4f'%(sampleSig[i])
                else:
                    sigstr += 14*' '
            
        elif 'Debye' in Sample['Type']:        #Debye-Scherrer
            for i,item in enumerate(['Scale','Absorption','DisplaceX','DisplaceY']):
                ptlbls += '%14s'%(item)
                ptstr += '%14.4f'%(Sample[item][0])
                if sampleSig[i]:
                    refine = True
                    sigstr += '%14.4f'%(sampleSig[i])
                else:
                    sigstr += 14*' '

        if refine:
            print '\n Sample Parameters:'
            print ptlbls
            print ptstr
            print sigstr
        
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            pfx = ':'+str(hId)+':'
            Background = Histogram['Background']
            backSig = SetBackgroundParms(pfx,Background,parmDict,sigDict)
            
            Inst = Histogram['Instrument Parameters']
            instSig = SetInstParms(pfx,Inst,parmDict,sigDict)
        
            Sample = Histogram['Sample Parameters']
            sampSig = SetSampleParms(pfx,Sample,parmDict,sigDict)

            print '\n Histogram: ',histogram,' histogram Id: ',hId
            print 135*'-'
            print ' Final refinement wRp = %.2f%% on %d observations in this histogram'%(Histogram['wRp'],Histogram['Nobs'])
            if Print:
                print ' Instrument type: ',Sample['Type']
                PrintSampleParmsSig(Sample,sampSig)
                PrintInstParmsSig(Inst,instSig)
                PrintBackgroundSig(Background,backSig)

def GetAtomFXU(pfx,FFtables,BLtables,calcControls,parmDict):
    Natoms = calcControls['Natoms'][pfx]
    Tdata = Natoms*[' ',]
    Mdata = np.zeros(Natoms)
    IAdata = Natoms*[' ',]
    Fdata = np.zeros(Natoms)
    FFdata = []
    BLdata = []
    Xdata = np.zeros((3,Natoms))
    dXdata = np.zeros((3,Natoms))
    Uisodata = np.zeros(Natoms)
    Uijdata = np.zeros((6,Natoms))
    keys = {'Atype:':Tdata,'Amul:':Mdata,'Afrac:':Fdata,'AI/A:':IAdata,
        'dAx:':dXdata[0],'dAy:':dXdata[1],'dAz:':dXdata[2],
        'Ax:':Xdata[0],'Ay:':Xdata[1],'Az:':Xdata[2],'AUiso:':Uisodata,
        'AU11:':Uijdata[0],'AU22:':Uijdata[1],'AU33:':Uijdata[2],
        'AU12:':Uijdata[3],'AU13:':Uijdata[4],'AU23:':Uijdata[5]}
    for iatm in range(Natoms):
        for key in keys:
            parm = pfx+key+str(iatm)
            if parm in parmDict:
                keys[key][iatm] = parmDict[parm]
        FFdata.append(FFtables[Tdata[iatm]])
        BLdata.append(BLtables[Tdata[iatm]][1])
    return FFdata,BLdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata
    
def StructureFactor(refList,G,hfx,pfx,SGData,calcControls,parmDict):
    ''' Compute structure factors for all h,k,l for phase
    input:
        refList: [ref] where each ref = h,k,l,m,d,...,[equiv h,k,l],phase[equiv] 
        G:      reciprocal metric tensor
        pfx:    phase id string
        SGData: space group info. dictionary output from SpcGroup
        calcControls:
        ParmDict:
    puts result F^2 in each ref[8] in refList
    '''        
    twopi = 2.0*np.pi
    twopisq = 2.0*np.pi**2
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    FFdata,BLdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata = GetAtomFXU(pfx,FFtables,BLtables,calcControls,parmDict)
    if 'N' in parmDict[hfx+'Type']:
        FP,FPP = G2el.BlenRes(BLdata,parmDict[hfx+'Lam'])
    else:
        FP = np.array([El[hfx+'FP'] for El in FFdata])
        FPP = np.array([El[hfx+'FPP'] for El in FFdata])
    maxPos = len(SGData['SGOps'])
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    for refl in refList:
        fbs = np.array([0,0])
        H = refl[:3]
        SQ = 1./(2.*refl[4])**2
        if 'N' in parmDict[hfx+'Type']:
            FF = np.array([El[1] for El in BLdata])
        else:       #'X'
            FF = np.array([G2el.ScatFac(El,SQ)[0] for El in FFdata])
        SQfactor = 4.0*SQ*twopisq
        Uniq = refl[11]
        phi = refl[12]
        phase = twopi*(np.inner(Uniq,(dXdata.T+Xdata.T))+phi[:,np.newaxis])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        occ = Mdata*Fdata/len(Uniq)
        biso = -SQfactor*Uisodata
        Tiso = np.where(biso<1.,np.exp(biso),1.0)
        HbH = np.array([-np.inner(h,np.inner(bij,h)) for h in Uniq])
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)
        Tcorr = Tiso*Tuij
        fa = np.array([(FF+FP)*occ*cosp*Tcorr,-FPP*occ*sinp*Tcorr])
        fas = np.sum(np.sum(fa,axis=1),axis=1)        #real
        if not SGData['SGInv']:
            fb = np.array([(FF+FP)*occ*sinp*Tcorr,FPP*occ*cosp*Tcorr])
            fbs = np.sum(np.sum(fb,axis=1),axis=1)
        fasq = fas**2
        fbsq = fbs**2        #imaginary
        refl[9] = np.sum(fasq)+np.sum(fbsq)
        refl[10] = atan2d(fbs[0],fas[0])
    return refList
    
def StructureFactorDerv(refList,G,hfx,pfx,SGData,calcControls,parmDict):
    twopi = 2.0*np.pi
    twopisq = 2.0*np.pi**2
    ast = np.sqrt(np.diag(G))
    Mast = twopisq*np.multiply.outer(ast,ast)
    FFtables = calcControls['FFtables']
    BLtables = calcControls['BLtables']
    FFdata,BLdata,Mdata,Fdata,Xdata,dXdata,IAdata,Uisodata,Uijdata = GetAtomFXU(pfx,FFtables,BLtables,calcControls,parmDict)
    if 'N' in parmDict[hfx+'Type']:
        FP = 0.
        FPP = 0.
    else:
        FP = np.array([El[hfx+'FP'] for El in FFdata])
        FPP = np.array([El[hfx+'FPP'] for El in FFdata])
    maxPos = len(SGData['SGOps'])       
    Uij = np.array(G2lat.U6toUij(Uijdata))
    bij = Mast*Uij.T
    dFdvDict = {}
    dFdfr = np.zeros((len(refList),len(Mdata)))
    dFdx = np.zeros((len(refList),len(Mdata),3))
    dFdui = np.zeros((len(refList),len(Mdata)))
    dFdua = np.zeros((len(refList),len(Mdata),6))
    for iref,refl in enumerate(refList):
        H = np.array(refl[:3])
        SQ = 1./(2.*refl[4])**2          # or (sin(theta)/lambda)**2
        if 'N' in parmDict[hfx+'Type']:
            FF = np.array([El[1] for El in BLdata])
        else:       #'X'
            FF = np.array([G2el.ScatFac(El,SQ)[0] for El in FFdata])
        SQfactor = 8.0*SQ*np.pi**2
        Uniq = refl[11]
        phi = refl[12]
        phase = twopi*(np.inner((dXdata.T+Xdata.T),Uniq)+phi[np.newaxis,:])
        sinp = np.sin(phase)
        cosp = np.cos(phase)
        occ = Mdata*Fdata/len(Uniq)
        biso = -SQfactor*Uisodata
        Tiso = np.where(biso<1.,np.exp(biso),1.0)
#        HbH = np.array([-np.inner(h,np.inner(bij,h)) for h in Uniq])
        HbH = -np.inner(H,np.inner(bij,H))
        Hij = np.array([Mast*np.multiply.outer(U,U) for U in Uniq])
        Hij = np.array([G2lat.UijtoU6(Uij) for Uij in Hij])
        Tuij = np.where(HbH<1.,np.exp(HbH),1.0)
        Tcorr = Tiso*Tuij
        fot = (FF+FP)*occ*Tcorr
        fotp = FPP*occ*Tcorr
        fa = np.array([fot[:,np.newaxis]*cosp,fotp[:,np.newaxis]*cosp])       #non positions
        fb = np.array([fot[:,np.newaxis]*sinp,-fotp[:,np.newaxis]*sinp])
        
        fas = np.sum(np.sum(fa,axis=1),axis=1)
        fbs = np.sum(np.sum(fb,axis=1),axis=1)
        fax = np.array([-fot[:,np.newaxis]*sinp,-fotp[:,np.newaxis]*sinp])   #positions
        fbx = np.array([fot[:,np.newaxis]*cosp,-fot[:,np.newaxis]*cosp])
        #sum below is over Uniq
        dfadfr = np.sum(fa/occ[:,np.newaxis],axis=2)
        dfadx = np.sum(twopi*Uniq*fax[:,:,:,np.newaxis],axis=2)
        dfadui = np.sum(-SQfactor*fa,axis=2)
        dfadua = np.sum(-Hij*fa[:,:,:,np.newaxis],axis=2)
        #NB: the above have been checked against PA(1:10,1:2) in strfctr.for      
        dFdfr[iref] = 2.*(fas[0]*dfadfr[0]+fas[1]*dfadfr[1])*Mdata/len(Uniq)
        dFdx[iref] = 2.*(fas[0]*dfadx[0]+fas[1]*dfadx[1])
        dFdui[iref] = 2.*(fas[0]*dfadui[0]+fas[1]*dfadui[1])
        dFdua[iref] = 2.*(fas[0]*dfadua[0]+fas[1]*dfadua[1])
        if not SGData['SGInv']:
            dfbdfr = np.sum(fb/occ[:,np.newaxis],axis=2)        #problem here if occ=0 for some atom
            dfbdx = np.sum(twopi*Uniq*fbx[:,:,:,np.newaxis],axis=2)          
            dfbdui = np.sum(-SQfactor*fb,axis=2)
            dfbdua = np.sum(-Hij*fb[:,:,:,np.newaxis],axis=2)
            dFdfr[iref] += 2.*(fbs[0]*dfbdfr[0]-fbs[1]*dfbdfr[1])*Mdata/len(Uniq)
            dFdx[iref] += 2.*(fbs[0]*dfbdx[0]+fbs[1]*dfbdx[1])
            dFdui[iref] += 2.*(fbs[0]*dfbdui[0]-fbs[1]*dfbdui[1])
            dFdua[iref] += 2.*(fbs[0]*dfbdua[0]+fbs[1]*dfbdua[1])
        #loop over atoms - each dict entry is list of derivatives for all the reflections
    for i in range(len(Mdata)):     
        dFdvDict[pfx+'Afrac:'+str(i)] = dFdfr.T[i]
        dFdvDict[pfx+'dAx:'+str(i)] = dFdx.T[0][i]
        dFdvDict[pfx+'dAy:'+str(i)] = dFdx.T[1][i]
        dFdvDict[pfx+'dAz:'+str(i)] = dFdx.T[2][i]
        dFdvDict[pfx+'AUiso:'+str(i)] = dFdui.T[i]
        dFdvDict[pfx+'AU11:'+str(i)] = dFdua.T[0][i]
        dFdvDict[pfx+'AU22:'+str(i)] = dFdua.T[1][i]
        dFdvDict[pfx+'AU33:'+str(i)] = dFdua.T[2][i]
        dFdvDict[pfx+'AU12:'+str(i)] = 2.*dFdua.T[3][i]
        dFdvDict[pfx+'AU13:'+str(i)] = 2.*dFdua.T[4][i]
        dFdvDict[pfx+'AU23:'+str(i)] = 2.*dFdua.T[5][i]
    return dFdvDict
        
def Dict2Values(parmdict, varylist):
    '''Use before call to leastsq to setup list of values for the parameters 
    in parmdict, as selected by key in varylist'''
    return [parmdict[key] for key in varylist] 
    
def Values2Dict(parmdict, varylist, values):
    ''' Use after call to leastsq to update the parameter dictionary with 
    values corresponding to keys in varylist'''
    parmdict.update(zip(varylist,values))
    
def GetNewCellParms(parmDict,varyList):
    newCellDict = {}
    Ddict = dict(zip(['D11','D22','D33','D12','D13','D23'],['A'+str(i) for i in range(6)]))
    for item in varyList:
        keys = item.split(':')
        if keys[2] in Ddict:
            key = keys[0]+'::'+Ddict[keys[2]]
            parm = keys[0]+'::'+keys[2]
            newCellDict[parm] = [key,parmDict[key]+parmDict[item]]
    return newCellDict
    
def ApplyXYZshifts(parmDict,varyList):
    ''' takes atom x,y,z shift and applies it to corresponding atom x,y,z value
        input:
            parmDict - parameter dictionary
            varyList - list of variables
        returns:
            newAtomDict - dictitemionary of new atomic coordinate names & values; 
                key is parameter shift name
    '''
    newAtomDict = {}
    for item in parmDict:
        if 'dA' in item:
            parm = ''.join(item.split('d'))
            parmDict[parm] += parmDict[item]
            newAtomDict[item] = [parm,parmDict[parm]]
    return newAtomDict
    
def SHTXcal(refl,g,pfx,hfx,SGData,calcControls,parmDict):
    IFCoup = 'Bragg' in calcControls[hfx+'instType']
    odfCor = 1.0
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangls = [parmDict[pfx+'SH omega'],parmDict[pfx+'SH chi'],parmDict[pfx+'SH phi']]
    Gangls = [parmDict[hfx+'Omega'],parmDict[hfx+'Chi'],parmDict[hfx+'Phi'],parmDict[hfx+'Azimuth']]
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,x,x = G2lat.SamAng(refl[5]/2.,Gangls,Sangls,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],parmDict[pfx+'SHmodel'],parmDict[pfx+'SHorder'])
    for item in SHnames:
        L,M,N = eval(item.strip('C'))
        Kcl = G2lat.GetKcl(L,N,SGData['SGLaue'],phi,beta)
        Ksl,x,x = G2lat.GetKsl(L,M,parmDict[pfx+'SHmodel'],psi,gam)
        Lnorm = G2lat.Lnorm(L)
        odfCor += parmDict[pfx+item]*Lnorm*Kcl*Ksl
    return odfCor
    
def SHTXcalDerv(refl,g,pfx,hfx,SGData,calcControls,parmDict):
    FORPI = 12.5663706143592
    IFCoup = 'Bragg' in calcControls[hfx+'instType']
    odfCor = 1.0
    dFdODF = {}
    dFdSA = [0,0,0]
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangls = [parmDict[pfx+'SH omega'],parmDict[pfx+'SH chi'],parmDict[pfx+'SH phi']]
    Gangls = [parmDict[hfx+'Omega'],parmDict[hfx+'Chi'],parmDict[hfx+'Phi'],parmDict[hfx+'Azimuth']]
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,dPSdA,dGMdA = G2lat.SamAng(refl[5]/2.,Gangls,Sangls,IFCoup)
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],parmDict[pfx+'SHmodel'],parmDict[pfx+'SHorder'])
    for item in SHnames:
        L,M,N = eval(item.strip('C'))
        Kcl = G2lat.GetKcl(L,N,SGData['SGLaue'],phi,beta)
        Ksl,dKsdp,dKsdg = G2lat.GetKsl(L,M,parmDict[pfx+'SHmodel'],psi,gam)
        Lnorm = G2lat.Lnorm(L)
        odfCor += parmDict[pfx+item]*Lnorm*Kcl*Ksl
        dFdODF[pfx+item] = Lnorm*Kcl*Ksl
        for i in range(3):
            dFdSA[i] += parmDict[pfx+item]*Lnorm*Kcl*(dKsdp*dPSdA[i]+dKsdg*dGMdA[i])
    return odfCor,dFdODF,dFdSA
    
def SHPOcal(refl,g,phfx,hfx,SGData,calcControls,parmDict):
    odfCor = 1.0
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangl = [0.,0.,0.]
    if 'Bragg' in calcControls[hfx+'instType']:
        Gangls = [0.,90.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = True
    else:
        Gangls = [0.,0.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = False
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,x,x = G2lat.SamAng(refl[5]/2.,Gangls,Sangl,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',calcControls[phfx+'SHord'],False)
    for item in SHnames:
        L,N = eval(item.strip('C'))
        Kcsl,Lnorm = G2lat.GetKclKsl(L,N,SGData['SGLaue'],psi,phi,beta) 
        odfCor += parmDict[phfx+item]*Lnorm*Kcsl
    return odfCor
    
def SHPOcalDerv(refl,g,phfx,hfx,SGData,calcControls,parmDict):
    FORPI = 12.5663706143592
    odfCor = 1.0
    dFdODF = {}
    H = refl[:3]
    cell = G2lat.Gmat2cell(g)
    Sangl = [0.,0.,0.]
    if 'Bragg' in calcControls[hfx+'instType']:
        Gangls = [0.,90.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = True
    else:
        Gangls = [0.,0.,0.,parmDict[hfx+'Azimuth']]
        IFCoup = False
    phi,beta = G2lat.CrsAng(H,cell,SGData)
    psi,gam,x,x = G2lat.SamAng(refl[5]/2.,Gangls,Sangl,IFCoup) #ignore 2 sets of angle derivs.
    SHnames = G2lat.GenSHCoeff(SGData['SGLaue'],'0',calcControls[phfx+'SHord'],False)
    for item in SHnames:
        L,N = eval(item.strip('C'))
        Kcsl,Lnorm = G2lat.GetKclKsl(L,N,SGData['SGLaue'],psi,phi,beta) 
        odfCor += parmDict[phfx+item]*Lnorm*Kcsl
        dFdODF[phfx+item] = Kcsl*Lnorm
    return odfCor,dFdODF
    
def GetPrefOri(refl,G,g,phfx,hfx,SGData,calcControls,parmDict):
    if calcControls[phfx+'poType'] == 'MD':
        MD = parmDict[phfx+'MD']
        MDAxis = calcControls[phfx+'MDAxis']
        sumMD = 0
        for H in refl[11]:            
            cosP,sinP = G2lat.CosSinAngle(H,MDAxis,G)
            A = 1.0/np.sqrt((MD*cosP)**2+sinP**2/MD)
            sumMD += A**3
        POcorr = sumMD/len(refl[11])
    else:   #spherical harmonics
        POcorr = SHPOcal(refl,g,phfx,hfx,SGData,calcControls,parmDict)
    return POcorr
    
def GetPrefOriDerv(refl,G,g,phfx,hfx,SGData,calcControls,parmDict):
    POderv = {}
    if calcControls[phfx+'poType'] == 'MD':
        MD = parmDict[phfx+'MD']
        MDAxis = calcControls[phfx+'MDAxis']
        sumMD = 0
        sumdMD = 0
        for H in refl[11]:            
            cosP,sinP = G2lat.CosSinAngle(H,MDAxis,G)
            A = 1.0/np.sqrt((MD*cosP)**2+sinP**2/MD)
            sumMD += A**3
            sumdMD -= (1.5*A**5)*(2.0*MD*cosP**2-(sinP/MD)**2)
        POcorr = sumMD/len(refl[11])
        POderv[phfx+'MD'] = sumdMD/len(refl[11])
    else:   #spherical harmonics
        POcorr,POderv = SHPOcalDerv(refl,g,phfx,hfx,SGData,calcControls,parmDict)
    return POcorr,POderv
    
def GetIntensityCorr(refl,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict):
    Icorr = parmDict[phfx+'Scale']*parmDict[hfx+'Scale']*refl[3]               #scale*multiplicity
    if 'X' in parmDict[hfx+'Type']:
        Icorr *= G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5],parmDict[hfx+'Azimuth'])[0]
    Icorr *= GetPrefOri(refl,G,g,phfx,hfx,SGData,calcControls,parmDict)
    if pfx+'SHorder' in parmDict:
        Icorr *= SHTXcal(refl,g,pfx,hfx,SGData,calcControls,parmDict)
    refl[13] = Icorr        
    
def GetIntensityDerv(refl,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict):
    dIdsh = 1./parmDict[hfx+'Scale']
    dIdsp = 1./parmDict[phfx+'Scale']
    if 'X' in parmDict[hfx+'Type']:
        pola,dIdPola = G2pwd.Polarization(parmDict[hfx+'Polariz.'],refl[5],parmDict[hfx+'Azimuth'])
        dIdPola /= pola
    else:       #'N'
        dIdPola = 0.0
    POcorr,dIdPO = GetPrefOriDerv(refl,G,g,phfx,hfx,SGData,calcControls,parmDict)
    for iPO in dIdPO:
        dIdPO[iPO] /= POcorr
    dFdODF = {}
    dFdSA = [0,0,0]
    if pfx+'SHorder' in parmDict:
        odfCor,dFdODF,dFdSA = SHTXcalDerv(refl,g,pfx,hfx,SGData,calcControls,parmDict)
        for iSH in dFdODF:
            dFdODF[iSH] /= odfCor
        for i in range(3):
            dFdSA[i] /= odfCor
    return dIdsh,dIdsp,dIdPola,dIdPO,dFdODF,dFdSA
        
def GetSampleGam(refl,wave,G,GB,phfx,calcControls,parmDict):
    costh = cosd(refl[5]/2.)
    #crystallite size
    if calcControls[phfx+'SizeType'] == 'isotropic':
        gam = 1.8*wave/(np.pi*parmDict[phfx+'Size:i']*costh)
    elif calcControls[phfx+'SizeType'] == 'uniaxial':
        H = np.array(refl[:3])
        P = np.array(calcControls[phfx+'SizeAxis'])
        cosP,sinP = G2lat.CosSinAngle(H,P,G)
        gam = (1.8*wave/np.pi)/(parmDict[phfx+'Size:i']*parmDict[phfx+'Size:a']*costh)
        gam *= np.sqrt((sinP*parmDict[phfx+'Size:a'])**2+(cosP*parmDict[phfx+'Size:i'])**2)
    else:           #ellipsoidal crystallites
        Sij =[parmDict[phfx+'Size:%d'%(i)] for i in range(6)]
        H = np.array(refl[:3])
        lenR = G2pwd.ellipseSize(H,Sij,GB)
        gam = 1.8*wave/(np.pi*costh*lenR)
    #microstrain                
    if calcControls[phfx+'MustrainType'] == 'isotropic':
        gam += 0.018*parmDict[phfx+'Mustrain:i']*tand(refl[5]/2.)/np.pi
    elif calcControls[phfx+'MustrainType'] == 'uniaxial':
        H = np.array(refl[:3])
        P = np.array(calcControls[phfx+'MustrainAxis'])
        cosP,sinP = G2lat.CosSinAngle(H,P,G)
        Si = parmDict[phfx+'Mustrain:i']
        Sa = parmDict[phfx+'Mustrain:a']
        gam += 0.018*Si*Sa*tand(refl[5]/2.)/(np.pi*np.sqrt((Si*cosP)**2+(Sa*sinP)**2))
    else:       #generalized - P.W. Stephens model
        pwrs = calcControls[phfx+'MuPwrs']
        sum = 0
        for i,pwr in enumerate(pwrs):
            sum += parmDict[phfx+'Mustrain:'+str(i)]*refl[0]**pwr[0]*refl[1]**pwr[1]*refl[2]**pwr[2]
        gam += 0.018*refl[4]**2*tand(refl[5]/2.)*sum            
    return gam
        
def GetSampleGamDerv(refl,wave,G,GB,phfx,calcControls,parmDict):
    gamDict = {}
    costh = cosd(refl[5]/2.)
    tanth = tand(refl[5]/2.)
    #crystallite size derivatives
    if calcControls[phfx+'SizeType'] == 'isotropic':
        gamDict[phfx+'Size:i'] = -1.80*wave/(np.pi*costh)
    elif calcControls[phfx+'SizeType'] == 'uniaxial':
        H = np.array(refl[:3])
        P = np.array(calcControls[phfx+'SizeAxis'])
        cosP,sinP = G2lat.CosSinAngle(H,P,G)
        Si = parmDict[phfx+'Size:i']
        Sa = parmDict[phfx+'Size:a']
        gami = (1.8*wave/np.pi)/(Si*Sa)
        sqtrm = np.sqrt((sinP*Sa)**2+(cosP*Si)**2)
        gam = gami*sqtrm/costh            
        gamDict[phfx+'Size:i'] = gami*Si*cosP**2/(sqtrm*costh)-gam/Si
        gamDict[phfx+'Size:a'] = gami*Sa*sinP**2/(sqtrm*costh)-gam/Sa         
    else:           #ellipsoidal crystallites
        const = 1.8*wave/(np.pi*costh)
        Sij =[parmDict[phfx+'Size:%d'%(i)] for i in range(6)]
        H = np.array(refl[:3])
        R,dRdS = G2pwd.ellipseSizeDerv(H,Sij,GB)
        for i,item in enumerate([phfx+'Size:%d'%(j) for j in range(6)]):
            gamDict[item] = -(const/R**2)*dRdS[i]
    #microstrain derivatives                
    if calcControls[phfx+'MustrainType'] == 'isotropic':
        gamDict[phfx+'Mustrain:i'] =  0.018*tanth/np.pi            
    elif calcControls[phfx+'MustrainType'] == 'uniaxial':
        H = np.array(refl[:3])
        P = np.array(calcControls[phfx+'MustrainAxis'])
        cosP,sinP = G2lat.CosSinAngle(H,P,G)
        Si = parmDict[phfx+'Mustrain:i']
        Sa = parmDict[phfx+'Mustrain:a']
        gami = 0.018*Si*Sa*tanth/np.pi
        sqtrm = np.sqrt((Si*cosP)**2+(Sa*sinP)**2)
        gam = gami/sqtrm
        gamDict[phfx+'Mustrain:i'] = gam/Si-gami*Si*cosP**2/sqtrm**3
        gamDict[phfx+'Mustrain:a'] = gam/Sa-gami*Sa*sinP**2/sqtrm**3
    else:       #generalized - P.W. Stephens model
        pwrs = calcControls[phfx+'MuPwrs']
        const = 0.018*refl[4]**2*tanth
        for i,pwr in enumerate(pwrs):
            gamDict[phfx+'Mustrain:'+str(i)] = const*refl[0]**pwr[0]*refl[1]**pwr[1]*refl[2]**pwr[2]
    return gamDict
        
def GetReflPos(refl,wave,G,hfx,calcControls,parmDict):
    h,k,l = refl[:3]
    dsq = 1./G2lat.calc_rDsq2(np.array([h,k,l]),G)
    d = np.sqrt(dsq)

    refl[4] = d
    pos = 2.0*asind(wave/(2.0*d))+parmDict[hfx+'Zero']
    const = 9.e-2/(np.pi*parmDict[hfx+'Gonio. radius'])                  #shifts in microns
    if 'Bragg' in calcControls[hfx+'instType']:
        pos -= const*(4.*parmDict[hfx+'Shift']*cosd(pos/2.0)+ \
            parmDict[hfx+'Transparency']*sind(pos)*100.0)            #trans(=1/mueff) in cm
    else:               #Debye-Scherrer - simple but maybe not right
        pos -= const*(parmDict[hfx+'DisplaceX']*cosd(pos)+parmDict[hfx+'DisplaceY']*sind(pos))
    return pos

def GetReflPosDerv(refl,wave,A,hfx,calcControls,parmDict):
    dpr = 180./np.pi
    h,k,l = refl[:3]
    dstsq = G2lat.calc_rDsq(np.array([h,k,l]),A)
    dst = np.sqrt(dstsq)
    pos = refl[5]-parmDict[hfx+'Zero']
    const = dpr/np.sqrt(1.0-wave**2*dstsq/4.0)
    dpdw = const*dst
    dpdA = np.array([h**2,k**2,l**2,h*k,h*l,k*l])
    dpdA *= const*wave/(2.0*dst)
    dpdZ = 1.0
    const = 9.e-2/(np.pi*parmDict[hfx+'Gonio. radius'])                  #shifts in microns
    if 'Bragg' in calcControls[hfx+'instType']:
        dpdSh = -4.*const*cosd(pos/2.0)
        dpdTr = -const*sind(pos)*100.0
        return dpdA,dpdw,dpdZ,dpdSh,dpdTr,0.,0.
    else:               #Debye-Scherrer - simple but maybe not right
        dpdXd = -const*cosd(pos)
        dpdYd = -const*sind(pos)
        return dpdA,dpdw,dpdZ,0.,0.,dpdXd,dpdYd
            
def GetHStrainShift(refl,SGData,phfx,parmDict):
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    h,k,l = refl[:3]
    if laue in ['m3','m3m']:
        Dij = parmDict[phfx+'D11']*(h**2+k**2+l**2)+ \
            refl[4]**2*parmDict[phfx+'eA']*((h*k)**2+(h*l)**2+(k*l)**2)/(h**2+k**2+l**2)**2
    elif laue in ['6/m','6/mmm','3m1','31m','3']:
        Dij = parmDict[phfx+'D11']*(h**2+k**2+h*k)+parmDict[phfx+'D33']*l**2
    elif laue in ['3R','3mR']:
        Dij = parmDict[phfx+'D11']*(h**2+k**2+l**2)+parmDict[phfx+'D12']*(h*k+h*l+k*l)
    elif laue in ['4/m','4/mmm']:
        Dij = parmDict[phfx+'D11']*(h**2+k**2)+parmDict[phfx+'D33']*l**2
    elif laue in ['mmm']:
        Dij = parmDict[phfx+'D11']*h**2+parmDict[phfx+'D22']*k**2+parmDict[phfx+'D33']*l**2
    elif laue in ['2/m']:
        Dij = parmDict[phfx+'D11']*h**2+parmDict[phfx+'D22']*k**2+parmDict[phfx+'D33']*l**2
        if uniq == 'a':
            Dij += parmDict[phfx+'D23']*k*l
        elif uniq == 'b':
            Dij += parmDict[phfx+'D13']*h*l
        elif uniq == 'c':
            Dij += parmDict[phfx+'D12']*h*k
    else:
        Dij = parmDict[phfx+'D11']*h**2+parmDict[phfx+'D22']*k**2+parmDict[phfx+'D33']*l**2+ \
            parmDict[phfx+'D12']*h*k+parmDict[phfx+'D13']*h*l+parmDict[phfx+'D23']*k*l
    return Dij*refl[4]**2*tand(refl[5]/2.0)
            
def GetHStrainShiftDerv(refl,SGData,phfx):
    laue = SGData['SGLaue']
    uniq = SGData['SGUniq']
    h,k,l = refl[:3]
    if laue in ['m3','m3m']:
        dDijDict = {phfx+'D11':h**2+k**2+l**2,
            phfx+'eA':((h*k)**2+(h*l)**2+(k*l)**2)/(h**2+k**2+l**2)**2}
    elif laue in ['6/m','6/mmm','3m1','31m','3']:
        dDijDict = {phfx+'D11':h**2+k**2+h*k,phfx+'D33':l**2}
    elif laue in ['3R','3mR']:
        dDijDict = {phfx+'D11':h**2+k**2+l**2,phfx+'D12':h*k+h*l+k*l}
    elif laue in ['4/m','4/mmm']:
        dDijDict = {phfx+'D11':h**2+k**2,phfx+'D33':l**2}
    elif laue in ['mmm']:
        dDijDict = {phfx+'D11':h**2,phfx+'D22':k**2,phfx+'D33':l**2}
    elif laue in ['2/m']:
        dDijDict = {phfx+'D11':h**2,phfx+'D22':k**2,phfx+'D33':l**2}
        if uniq == 'a':
            dDijDict[phfx+'D23'] = k*l
        elif uniq == 'b':
            dDijDict[phfx+'D13'] = h*l
        elif uniq == 'c':
            dDijDict[phfx+'D12'] = h*k
            names.append()
    else:
        dDijDict = {phfx+'D11':h**2,phfx+'D22':k**2,phfx+'D33':l**2,
            phfx+'D12':h*k,phfx+'D13':h*l,phfx+'D23':k*l}
    for item in dDijDict:
        dDijDict[item] *= refl[4]**2*tand(refl[5]/2.0)
    return dDijDict
    
def GetFprime(controlDict,Histograms):
    FFtables = controlDict['FFtables']
    if not FFtables:
        return
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            keV = controlDict[hfx+'keV']
            for El in FFtables:
                Orbs = G2el.GetXsectionCoeff(El.split('+')[0].split('-')[0])
                FP,FPP,Mu = G2el.FPcalc(Orbs, keV)
                FFtables[El][hfx+'FP'] = FP
                FFtables[El][hfx+'FPP'] = FPP                
            
def getPowderProfile(parmDict,x,varylist,Histogram,Phases,calcControls,pawleyLookup):
    
    def GetReflSIgGam(refl,wave,G,GB,hfx,phfx,calcControls,parmDict):
        U = parmDict[hfx+'U']
        V = parmDict[hfx+'V']
        W = parmDict[hfx+'W']
        X = parmDict[hfx+'X']
        Y = parmDict[hfx+'Y']
        tanPos = tand(refl[5]/2.0)
        sig = U*tanPos**2+V*tanPos+W        #save peak sigma
        sig = max(0.001,sig)
        gam = X/cosd(refl[5]/2.0)+Y*tanPos+GetSampleGam(refl,wave,G,GB,phfx,calcControls,parmDict) #save peak gamma
        gam = max(0.001,gam)
        return sig,gam
                
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    bakType = calcControls[hfx+'bakType']
    yb = G2pwd.getBackground(hfx,parmDict,bakType,x)
    yc = np.zeros_like(yb)
        
    if 'C' in calcControls[hfx+'histType']:    
        shl = max(parmDict[hfx+'SH/L'],0.002)
        Ka2 = False
        if hfx+'Lam1' in parmDict.keys():
            wave = parmDict[hfx+'Lam1']
            Ka2 = True
            lamRatio = 360*(parmDict[hfx+'Lam2']-parmDict[hfx+'Lam1'])/(np.pi*parmDict[hfx+'Lam1'])
            kRatio = parmDict[hfx+'I(L2)/I(L1)']
        else:
            wave = parmDict[hfx+'Lam']
    else:
        print 'TOF Undefined at present'
        raise ValueError
    for phase in Histogram['Reflection Lists']:
        refList = Histogram['Reflection Lists'][phase]
        Phase = Phases[phase]
        pId = Phase['pId']
        pfx = '%d::'%(pId)
        phfx = '%d:%d:'%(pId,hId)
        hfx = ':%d:'%(hId)
        SGData = Phase['General']['SGData']
        A = [parmDict[pfx+'A%d'%(i)] for i in range(6)]
        G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
        GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
        Vst = np.sqrt(nl.det(G))    #V*
        if 'Pawley' not in Phase['General']['Type']:
            refList = StructureFactor(refList,G,hfx,pfx,SGData,calcControls,parmDict)
        for refl in refList:
            if 'C' in calcControls[hfx+'histType']:
                h,k,l = refl[:3]
                refl[5] = GetReflPos(refl,wave,G,hfx,calcControls,parmDict)         #corrected reflection position
                Lorenz = 1./(2.*sind(refl[5]/2.)**2*cosd(refl[5]/2.))           #Lorentz correction
                refl[5] += GetHStrainShift(refl,SGData,phfx,parmDict)               #apply hydrostatic strain shift
                refl[6:8] = GetReflSIgGam(refl,wave,G,GB,hfx,phfx,calcControls,parmDict)    #peak sig & gam
                GetIntensityCorr(refl,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)    #puts corrections in refl[13]
                refl[13] *= Vst*Lorenz
                if 'Pawley' in Phase['General']['Type']:
                    try:
                        refl[9] = abs(parmDict[pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])])
                    except KeyError:
#                        print ' ***Error %d,%d,%d missing from Pawley reflection list ***'%(h,k,l)
                        continue
                Wd,fmin,fmax = G2pwd.getWidths(refl[5],refl[6],refl[7],shl)
                iBeg = np.searchsorted(x,refl[5]-fmin)
                iFin = np.searchsorted(x,refl[5]+fmax)
                if not iBeg+iFin:       #peak below low limit - skip peak
                    continue
                elif not iBeg-iFin:     #peak above high limit - done
                    break
                yc[iBeg:iFin] += refl[13]*refl[9]*G2pwd.getFCJVoigt3(refl[5],refl[6],refl[7],shl,x[iBeg:iFin])    #>90% of time spent here
                if Ka2:
                    pos2 = refl[5]+lamRatio*tand(refl[5]/2.0)       # + 360/pi * Dlam/lam * tan(th)
                    Wd,fmin,fmax = G2pwd.getWidths(pos2,refl[6],refl[7],shl)
                    iBeg = np.searchsorted(x,pos2-fmin)
                    iFin = np.searchsorted(x,pos2+fmax)
                    if not iBeg+iFin:       #peak below low limit - skip peak
                        continue
                    elif not iBeg-iFin:     #peak above high limit - done
                        return yc,yb
                    yc[iBeg:iFin] += refl[13]*refl[9]*kRatio*G2pwd.getFCJVoigt3(pos2,refl[6],refl[7],shl,x[iBeg:iFin])        #and here
            elif 'T' in calcControls[hfx+'histType']:
                print 'TOF Undefined at present'
                raise Exception    #no TOF yet
    return yc,yb
    
def GetFobsSq(Histograms,Phases,parmDict,calcControls):
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            Limits = calcControls[hfx+'Limits']
            shl = max(parmDict[hfx+'SH/L'],0.002)
            Ka2 = False
            kRatio = 0.0
            if hfx+'Lam1' in parmDict.keys():
                Ka2 = True
                lamRatio = 360*(parmDict[hfx+'Lam2']-parmDict[hfx+'Lam1'])/(np.pi*parmDict[hfx+'Lam1'])
                kRatio = parmDict[hfx+'I(L2)/I(L1)']
            x,y,w,yc,yb,yd = Histogram['Data']
            ymb = np.array(y-yb)
            ycmb = np.array(yc-yb)
            ratio = np.where(ycmb!=0.,ymb/ycmb,0.0)            
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])
            refLists = Histogram['Reflection Lists']
            for phase in refLists:
                Phase = Phases[phase]
                pId = Phase['pId']
                phfx = '%d:%d:'%(pId,hId)
                refList = refLists[phase]
                sumFo = 0.0
                sumdF = 0.0
                sumFosq = 0.0
                sumdFsq = 0.0
                for refl in refList:
                    if 'C' in calcControls[hfx+'histType']:
                        yp = np.zeros_like(yb)
                        Wd,fmin,fmax = G2pwd.getWidths(refl[5],refl[6],refl[7],shl)
                        iBeg = np.searchsorted(x[xB:xF],refl[5]-fmin)
                        iFin = np.searchsorted(x[xB:xF],refl[5]+fmax)
                        iFin2 = iFin
                        yp[iBeg:iFin] = refl[13]*refl[9]*G2pwd.getFCJVoigt3(refl[5],refl[6],refl[7],shl,x[iBeg:iFin])    #>90% of time spent here                            
                        if Ka2:
                            pos2 = refl[5]+lamRatio*tand(refl[5]/2.0)       # + 360/pi * Dlam/lam * tan(th)
                            Wd,fmin,fmax = G2pwd.getWidths(pos2,refl[6],refl[7],shl)
                            iBeg2 = np.searchsorted(x,pos2-fmin)
                            iFin2 = np.searchsorted(x,pos2+fmax)
                            yp[iBeg2:iFin2] += refl[13]*refl[9]*kRatio*G2pwd.getFCJVoigt3(pos2,refl[6],refl[7],shl,x[iBeg2:iFin2])        #and here
                        refl[8] = np.sum(np.where(ratio[iBeg:iFin2]>0.,yp[iBeg:iFin2]*ratio[iBeg:iFin2]/(refl[13]*(1.+kRatio)),0.0))
                    elif 'T' in calcControls[hfx+'histType']:
                        print 'TOF Undefined at present'
                        raise Exception    #no TOF yet
                    Fo = np.sqrt(np.abs(refl[8]))
                    Fc = np.sqrt(np.abs(refl[9]))
                    sumFo += Fo
                    sumFosq += refl[8]**2
                    sumdF += np.abs(Fo-Fc)
                    sumdFsq += (refl[8]-refl[9])**2
                Histogram[phfx+'Rf'] = min(100.,(sumdF/sumFo)*100.)
                Histogram[phfx+'Rf^2'] = min(100.,np.sqrt(sumdFsq/sumFosq)*100.)
                Histogram[phfx+'Nref'] = len(refList)
                
def getPowderProfileDerv(parmDict,x,varylist,Histogram,Phases,calcControls,pawleyLookup):
    
    def cellVaryDerv(pfx,SGData,dpdA): 
        if SGData['SGLaue'] in ['-1',]:
            return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],
                [pfx+'A3',dpdA[3]],[pfx+'A4',dpdA[4]],[pfx+'A5',dpdA[5]]]
        elif SGData['SGLaue'] in ['2/m',]:
            if SGData['SGUniq'] == 'a':
                return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A3',dpdA[3]]]
            elif SGData['SGUniq'] == 'b':
                return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A4',dpdA[4]]]
            else:
                return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A5',dpdA[5]]]
        elif SGData['SGLaue'] in ['mmm',]:
            return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]]]
        elif SGData['SGLaue'] in ['4/m','4/mmm']:
#            return [[pfx+'A0',dpdA[0]+dpdA[1]],[pfx+'A2',dpdA[2]]]
            return [[pfx+'A0',dpdA[0]],[pfx+'A2',dpdA[2]]]
        elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
#            return [[pfx+'A0',dpdA[0]+dpdA[1]+dpdA[3]],[pfx+'A2',dpdA[2]]]
            return [[pfx+'A0',dpdA[0]],[pfx+'A2',dpdA[2]]]
        elif SGData['SGLaue'] in ['3R', '3mR']:
            return [[pfx+'A0',dpdA[0]+dpdA[1]+dpdA[2]],[pfx+'A3',dpdA[3]+dpdA[4]+dpdA[5]]]                       
        elif SGData['SGLaue'] in ['m3m','m3']:
#            return [[pfx+'A0',dpdA[0]+dpdA[1]+dpdA[2]]]
            return [[pfx+'A0',dpdA[0]]]
    # create a list of dependent variables and set up a dictionary to hold their derivatives
    dependentVars = G2mv.GetDependentVars()
    depDerivDict = {}
    for j in dependentVars:
        depDerivDict[j] = np.zeros(shape=(len(x)))
    #print 'dependent vars',dependentVars
    lenX = len(x)                
    hId = Histogram['hId']
    hfx = ':%d:'%(hId)
    bakType = calcControls[hfx+'bakType']
    dMdv = np.zeros(shape=(len(varylist),len(x)))
    dMdb,dMddb = G2pwd.getBackgroundDerv(hfx,parmDict,bakType,x)
    if hfx+'Back:0' in varylist: # for now assume that Back:x vars to not appear in constraints
        bBpos =varylist.index(hfx+'Back:0')
        dMdv[bBpos:bBpos+len(dMdb)] = dMdb
    names = [hfx+'DebyeA',hfx+'DebyeR',hfx+'DebyeU']
    for name in varylist:
        if 'Debye' in name:
            id = int(name.split(':')[-1])
            parm = name[:int(name.rindex(':'))]
            ip = names.index(parm)
            dMdv[varylist.index(name)] = dMddb[3*id+ip]
    if 'C' in calcControls[hfx+'histType']:    
        dx = x[1]-x[0]
        shl = max(parmDict[hfx+'SH/L'],0.002)
        Ka2 = False
        if hfx+'Lam1' in parmDict.keys():
            wave = parmDict[hfx+'Lam1']
            Ka2 = True
            lamRatio = 360*(parmDict[hfx+'Lam2']-parmDict[hfx+'Lam1'])/(np.pi*parmDict[hfx+'Lam1'])
            kRatio = parmDict[hfx+'I(L2)/I(L1)']
        else:
            wave = parmDict[hfx+'Lam']
    else:
        print 'TOF Undefined at present'
        raise ValueError
    for phase in Histogram['Reflection Lists']:
        refList = Histogram['Reflection Lists'][phase]
        Phase = Phases[phase]
        SGData = Phase['General']['SGData']
        pId = Phase['pId']
        pfx = '%d::'%(pId)
        phfx = '%d:%d:'%(pId,hId)
        A = [parmDict[pfx+'A%d'%(i)] for i in range(6)]
        G,g = G2lat.A2Gmat(A)       #recip & real metric tensors
        GA,GB = G2lat.Gmat2AB(G)    #Orthogonalization matricies
        if 'Pawley' not in Phase['General']['Type']:
            dFdvDict = StructureFactorDerv(refList,G,hfx,pfx,SGData,calcControls,parmDict)
        for iref,refl in enumerate(refList):
            if 'C' in calcControls[hfx+'histType']:        #CW powder
                h,k,l = refl[:3]
                dIdsh,dIdsp,dIdpola,dIdPO,dFdODF,dFdSA = GetIntensityDerv(refl,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
                if 'Pawley' in Phase['General']['Type']:
                    try:
                        refl[9] = abs(parmDict[pfx+'PWLref:%d'%(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])])
                    except KeyError:
#                        print ' ***Error %d,%d,%d missing from Pawley reflection list ***'%(h,k,l)
                        continue
                Wd,fmin,fmax = G2pwd.getWidths(refl[5],refl[6],refl[7],shl)
                iBeg = np.searchsorted(x,refl[5]-fmin)
                iFin = np.searchsorted(x,refl[5]+fmax)
                if not iBeg+iFin:       #peak below low limit - skip peak
                    continue
                elif not iBeg-iFin:     #peak above high limit - done
                    break
                pos = refl[5]
                tanth = tand(pos/2.0)
                costh = cosd(pos/2.0)
                dMdpk = np.zeros(shape=(6,len(x)))
                dMdipk = G2pwd.getdFCJVoigt3(refl[5],refl[6],refl[7],shl,x[iBeg:iFin])
                for i in range(1,5):
                    dMdpk[i][iBeg:iFin] += 100.*dx*refl[13]*refl[9]*dMdipk[i]
                dMdpk[0][iBeg:iFin] += 100.*dx*refl[13]*refl[9]*dMdipk[0]
                dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4]}
                if Ka2:
                    pos2 = refl[5]+lamRatio*tanth       # + 360/pi * Dlam/lam * tan(th)
                    kdelt = int((pos2-refl[5])/dx)               
                    iBeg2 = min(lenX,iBeg+kdelt)
                    iFin2 = min(lenX,iFin+kdelt)
                    if iBeg2-iFin2:
                        dMdipk2 = G2pwd.getdFCJVoigt3(pos2,refl[6],refl[7],shl,x[iBeg2:iFin2])
                        for i in range(1,5):
                            dMdpk[i][iBeg2:iFin2] += 100.*dx*refl[13]*refl[9]*kRatio*dMdipk2[i]
                        dMdpk[0][iBeg2:iFin2] += 100.*dx*refl[13]*refl[9]*kRatio*dMdipk2[0]
                        dMdpk[5][iBeg2:iFin2] += 100.*dx*refl[13]*dMdipk2[0]
                        dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4],'L1/L2':dMdpk[5]*refl[9]}
                if 'Pawley' in Phase['General']['Type']:
                    try:
                        idx = varylist.index(pfx+'PWLref:'+str(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)]))
                        dMdv[idx] = dervDict['int']/refl[9]
                        # Assuming Pawley variables not in constraints
                    except ValueError:
                        pass
                dpdA,dpdw,dpdZ,dpdSh,dpdTr,dpdX,dpdY = GetReflPosDerv(refl,wave,A,hfx,calcControls,parmDict)
                names = {hfx+'Scale':[dIdsh,'int'],hfx+'Polariz.':[dIdpola,'int'],phfx+'Scale':[dIdsp,'int'],
                    hfx+'U':[tanth**2,'sig'],hfx+'V':[tanth,'sig'],hfx+'W':[1.0,'sig'],
                    hfx+'X':[1.0/costh,'gam'],hfx+'Y':[tanth,'gam'],hfx+'SH/L':[1.0,'shl'],
                    hfx+'I(L2)/I(L1)':[1.0,'L1/L2'],hfx+'Zero':[dpdZ,'pos'],hfx+'Lam':[dpdw,'pos'],
                    hfx+'Shift':[dpdSh,'pos'],hfx+'Transparency':[dpdTr,'pos'],hfx+'DisplaceX':[dpdX,'pos'],
                    hfx+'DisplaceY':[dpdY,'pos'],}
                for name in names:
                    item = names[name]
                    if name in varylist:
                        dMdv[varylist.index(name)] += item[0]*dervDict[item[1]]
                    elif name in dependentVars:
                        depDerivDict[name] += item[0]*dervDict[item[1]]

                for iPO in dIdPO:
                    if iPO in varylist:
                        dMdv[varylist.index(iPO)] += dIdPO[iPO]*dervDict['int']
                    elif iPO in dependentVars:
                        depDerivDict[iPO] = dIdPO[iPO]*dervDict['int']

                for i,name in enumerate(['omega','chi','phi']):
                    aname = pfx+'SH '+name
                    if aname in varylist:
                        dMdv[varylist.index(aname)] += dFdSA[i]*dervDict['int']
                    elif aname in dependentVars:
                        depDerivDict[aname] += dFdSA[i]*dervDict['int']
                for iSH in dFdODF:
                    if iSH in varylist:
                        dMdv[varylist.index(iSH)] += dFdODF[iSH]*dervDict['int']
                    elif iSH in dependentVars:
                        depDerivDict[iSH] += dFdODF[iSH]*dervDict['int']
                cellDervNames = cellVaryDerv(pfx,SGData,dpdA)
                for name,dpdA in cellDervNames:
                    if name in varylist:
                        dMdv[varylist.index(name)] += dpdA*dervDict['pos']
                    elif name in dependentVars:
                        depDerivDict[name] += dpdA*dervDict['pos']
                dDijDict = GetHStrainShiftDerv(refl,SGData,phfx)
                for name in dDijDict:
                    if name in varylist:
                        dMdv[varylist.index(name)] += dDijDict[name]*dervDict['pos']
                    elif name in dependentVars:
                        depDerivDict[name] += dDijDict[name]*dervDict['pos']
                gamDict = GetSampleGamDerv(refl,wave,G,GB,phfx,calcControls,parmDict)
                for name in gamDict:
                    if name in varylist:
                        dMdv[varylist.index(name)] += gamDict[name]*dervDict['gam']
                    elif name in dependentVars:
                        depDerivDict[name] += gamDict[name]*dervDict['gam']
                                               
            elif 'T' in calcControls[hfx+'histType']:
                print 'TOF Undefined at present'
                raise Exception    #no TOF yet
            #do atom derivatives -  for F,X & U so far              
            corr = dervDict['int']/refl[9]
            for name in varylist+dependentVars:
                try:
                    aname = name.split(pfx)[1][:2]
                    if aname not in ['Af','dA','AU']: continue # skip anything not an atom param
                except IndexError:
                    continue
                if name in varylist:
                    dMdv[varylist.index(name)] += dFdvDict[name][iref]*corr
                elif name in dependentVars:
                    depDerivDict[name] += dFdvDict[name][iref]*corr
    # now process derivatives in constraints
    G2mv.Dict2Deriv(varylist,depDerivDict,dMdv)
    return dMdv

def dervRefine(values,HistoPhases,parmdict,varylist,calcControls,pawleyLookup,dlg):
    parmdict.update(zip(varylist,values))
    G2mv.Dict2Map(parmdict,varylist)
    Histograms,Phases = HistoPhases
    nvar = len(varylist)
    dMdv = np.empty(0)
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            Limits = calcControls[hfx+'Limits']
            x,y,w,yc,yb,yd = Histogram['Data']
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])
            dMdvh = np.sqrt(w[xB:xF])*getPowderProfileDerv(parmdict,x[xB:xF],
                varylist,Histogram,Phases,calcControls,pawleyLookup)
            if len(dMdv):
                dMdv = np.concatenate((dMdv.T,dMdvh.T)).T
            else:
                dMdv = dMdvh
    return dMdv

def HessRefine(values,HistoPhases,parmdict,varylist,calcControls,pawleyLookup,dlg):
    parmdict.update(zip(varylist,values))
    G2mv.Dict2Map(parmdict,varylist)
    Histograms,Phases = HistoPhases
    nvar = len(varylist)
    Hess = np.empty(0)
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            Limits = calcControls[hfx+'Limits']
            x,y,w,yc,yb,yd = Histogram['Data']
            dy = y-yc
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])
            dMdvh = np.sqrt(w[xB:xF])*getPowderProfileDerv(parmdict,x[xB:xF],
                varylist,Histogram,Phases,calcControls,pawleyLookup)
            if dlg:
                dlg.Update(Histogram['wRp'],newmsg='Hessian for histogram %d Rwp=%8.3f%s'%(hId,Histogram['wRp'],'%'))[0]
            if len(Hess):
                Vec += np.sum(dMdvh*np.sqrt(w[xB:xF])*dy[xB:xF],axis=1)
                Hess += np.inner(dMdvh,dMdvh)
            else:
                Vec = np.sum(dMdvh*np.sqrt(w[xB:xF])*dy[xB:xF],axis=1)
                Hess = np.inner(dMdvh,dMdvh)
    return Vec,Hess

def errRefine(values,HistoPhases,parmdict,varylist,calcControls,pawleyLookup,dlg):        
    parmdict.update(zip(varylist,values))
    Values2Dict(parmdict, varylist, values)
    G2mv.Dict2Map(parmdict,varylist)
    Histograms,Phases = HistoPhases
    M = np.empty(0)
    sumwYo = 0
    Nobs = 0
    histoList = Histograms.keys()
    histoList.sort()
    for histogram in histoList:
        if 'PWDR' in histogram[:4]:
            Histogram = Histograms[histogram]
            hId = Histogram['hId']
            hfx = ':%d:'%(hId)
            Limits = calcControls[hfx+'Limits']
            x,y,w,yc,yb,yd = Histogram['Data']
            yc *= 0.0                           #zero full calcd profiles
            yb *= 0.0
            yd *= 0.0
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])
            Histogram['Nobs'] = xF-xB
            Nobs += Histogram['Nobs']
            Histogram['sumwYo'] = np.sum(w[xB:xF]*y[xB:xF]**2)
            sumwYo += Histogram['sumwYo']
            yc[xB:xF],yb[xB:xF] = getPowderProfile(parmdict,x[xB:xF],
                varylist,Histogram,Phases,calcControls,pawleyLookup)
            yc[xB:xF] += yb[xB:xF]
            yd[xB:xF] = y[xB:xF]-yc[xB:xF]
            Histogram['sumwYd'] = np.sum(np.sqrt(w[xB:xF])*(yd[xB:xF]))
            wdy = -np.sqrt(w[xB:xF])*(yd[xB:xF])
            Histogram['wRp'] = min(100.,np.sqrt(np.sum(wdy**2)/Histogram['sumwYo'])*100.)
            if dlg:
                dlg.Update(Histogram['wRp'],newmsg='For histogram %d Rwp=%8.3f%s'%(hId,Histogram['wRp'],'%'))[0]
            M = np.concatenate((M,wdy))
    Histograms['sumwYo'] = sumwYo
    Histograms['Nobs'] = Nobs
    Rwp = min(100.,np.sqrt(np.sum(M**2)/sumwYo)*100.)
    if dlg:
        GoOn = dlg.Update(Rwp,newmsg='%s%8.3f%s'%('Powder profile Rwp =',Rwp,'%'))[0]
        if not GoOn:
            parmDict['saved values'] = values
            raise Exception         #Abort!!
    return M
    
                    
def Refine(GPXfile,dlg):
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics
    
    ShowBanner()
    varyList = []
    parmDict = {}
    calcControls = {}
    G2mv.InitVars()    
    Controls = G2IO.GetControls(GPXfile)
    ShowControls(Controls)            
    constrDict,constrFlag,fixedList = G2IO.GetConstraints(GPXfile)
    Histograms,Phases = G2IO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print ' *** ERROR - you have no histograms to refine! ***'
        print ' *** Refine aborted ***'
        raise Exception
    if not Histograms:
        print ' *** ERROR - you have no data to refine with! ***'
        print ' *** Refine aborted ***'
        raise Exception        
    Natoms,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables = GetPhaseData(Phases)
    calcControls['Natoms'] = Natoms
    calcControls['FFtables'] = FFtables
    calcControls['BLtables'] = BLtables
    hapVary,hapDict,controlDict = GetHistogramPhaseData(Phases,Histograms)
    calcControls.update(controlDict)
    histVary,histDict,controlDict = GetHistogramData(Histograms)
    calcControls.update(controlDict)
    varyList = phaseVary+hapVary+histVary
    parmDict.update(phaseDict)
    parmDict.update(hapDict)
    parmDict.update(histDict)
    GetFprime(calcControls,Histograms)
    # do constraint processing
    try:
        groups,parmlist = G2mv.GroupConstraints(constrDict)
        G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,constrFlag,fixedList)
    except:
        print ' *** ERROR - your constraints are internally inconsistent ***'
        raise Exception(' *** Refine aborted ***')
    # check to see which generated parameters are fully varied
    msg = G2mv.SetVaryFlags(varyList)
    if msg:
        print ' *** ERROR - you have not set the refine flags for constraints consistently! ***'
        print msg
        raise Exception(' *** Refine aborted ***')
    G2mv.Map2Dict(parmDict,varyList)
#    print G2mv.VarRemapShow(varyList)

    while True:
        begin = time.time()
        values =  np.array(Dict2Values(parmDict, varyList))
        Ftol = Controls['min dM/M']
        Factor = Controls['shift factor']
        maxCyc = Controls['max cyc']
        if 'Jacobian' in Controls['deriv type']:            
            result = so.leastsq(errRefine,values,Dfun=dervRefine,full_output=True,
                ftol=Ftol,col_deriv=True,factor=Factor,
                args=([Histograms,Phases],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/2)
        elif 'Hessian' in Controls['deriv type']:
            result = G2mth.HessianLSQ(errRefine,values,Hess=HessRefine,ftol=Ftol,maxcyc=maxCyc,
                args=([Histograms,Phases],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = result[2]['num cyc']+1                           
        else:           #'numeric'
            result = so.leastsq(errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                args=([Histograms,Phases],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/len(varyList))
#        table = dict(zip(varyList,zip(values,result[0],(result[0]-values))))
#        for item in table: print item,table[item]               #useful debug - are things shifting?
        runtime = time.time()-begin
        chisq = np.sum(result[2]['fvec']**2)
        Values2Dict(parmDict, varyList, result[0])
        G2mv.Dict2Map(parmDict,varyList)
        
        Rwp = np.sqrt(chisq/Histograms['sumwYo'])*100.      #to %
        GOF = chisq/(Histograms['Nobs']-len(varyList))
        print '\n Refinement results:'
        print 135*'-'
        print ' Number of function calls:',result[2]['nfev'],' Number of observations: ',Histograms['Nobs'],' Number of parameters: ',len(varyList)
        print ' Refinement time = %8.3fs, %8.3fs/cycle, for %d cycles'%(runtime,runtime/ncyc,ncyc)
        print ' wRp = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rwp,chisq,GOF)
        try:
            covMatrix = result[1]*GOF
            sig = np.sqrt(np.diag(covMatrix))
            if np.any(np.isnan(sig)):
                print '*** Least squares aborted - some invalid esds possible ***'
#            table = dict(zip(varyList,zip(values,result[0],(result[0]-values)/sig)))
#            for item in table: print item,table[item]               #useful debug - are things shifting?
            break                   #refinement succeeded - finish up!
        except TypeError:          #result[1] is None on singular matrix
            print '**** Refinement failed - singular matrix ****'
            if 'Hessian' in Controls['deriv type']:
                for i in result[2]['psing'].reverse():
                        print 'Removing parameter: ',varyList[i]
                        del(varyList[i])                    
            else:
                Ipvt = result[2]['ipvt']
                for i,ipvt in enumerate(Ipvt):
                    if not np.sum(result[2]['fjac'],axis=1)[i]:
                        print 'Removing parameter: ',varyList[ipvt-1]
                        del(varyList[ipvt-1])
                        break

#    print 'dependentParmList: ',G2mv.dependentParmList
#    print 'arrayList: ',G2mv.arrayList
#    print 'invarrayList: ',G2mv.invarrayList
#    print 'indParmList: ',G2mv.indParmList
#    print 'fixedDict: ',G2mv.fixedDict
#    print 'test1'
    GetFobsSq(Histograms,Phases,parmDict,calcControls)
#    print 'test2'
    sigDict = dict(zip(varyList,sig))
    newCellDict = GetNewCellParms(parmDict,varyList)
    newAtomDict = ApplyXYZshifts(parmDict,varyList)
    covData = {'variables':result[0],'varyList':varyList,'sig':sig,
        'covMatrix':covMatrix,'title':GPXfile,'newAtomDict':newAtomDict,'newCellDict':newCellDict}
    # add the uncertainties into the esd dictionary (sigDict)
    sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList,parmDict))
    SetPhaseData(parmDict,sigDict,Phases,covData)
    SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms)
    SetHistogramData(parmDict,sigDict,Histograms)
    G2mv.PrintIndependentVars(parmDict,varyList,sigDict)
    G2IO.SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases,covData)
    
#for testing purposes!!!
#    import cPickle
#    file = open('structTestdata.dat','wb')
#    cPickle.dump(parmDict,file,1)
#    cPickle.dump(varyList,file,1)
#    for histogram in Histograms:
#        if 'PWDR' in histogram[:4]:
#            Histogram = Histograms[histogram]
#    cPickle.dump(Histogram,file,1)
#    cPickle.dump(Phases,file,1)
#    cPickle.dump(calcControls,file,1)
#    cPickle.dump(pawleyLookup,file,1)
#    file.close()

    if dlg:
        return Rwp

def SeqRefine(GPXfile,dlg):
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics
    
    ShowBanner()
    print ' Sequential Refinement'
    G2mv.InitVars()    
    Controls = G2IO.GetControls(GPXfile)
    ShowControls(Controls)            
    constrDict,constrFlag,fixedList = G2IO.GetConstraints(GPXfile)
    Histograms,Phases = G2IO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print ' *** ERROR - you have no histograms to refine! ***'
        print ' *** Refine aborted ***'
        raise Exception
    if not Histograms:
        print ' *** ERROR - you have no data to refine with! ***'
        print ' *** Refine aborted ***'
        raise Exception
    Natoms,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables = GetPhaseData(Phases,False)
    if 'Seq Data' in Controls:
        histNames = Controls['Seq Data']
    else:
        histNames = G2IO.GetHistogramNames(GPXfile,['PWDR',])
    if 'Reverse Seq' in Controls:
        if Controls['Reverse Seq']:
            histNames.reverse()
    SeqResult = {'histNames':histNames}
    makeBack = True
    for ihst,histogram in enumerate(histNames):
        ifPrint = False
        if dlg:
            dlg.SetTitle('Residual for histogram '+str(ihst))
        calcControls = {}
        calcControls['Natoms'] = Natoms
        calcControls['FFtables'] = FFtables
        calcControls['BLtables'] = BLtables
        varyList = []
        parmDict = {}
        Histo = {histogram:Histograms[histogram],}
        hapVary,hapDict,controlDict = GetHistogramPhaseData(Phases,Histo,False)
        calcControls.update(controlDict)
        histVary,histDict,controlDict = GetHistogramData(Histo,False)
        calcControls.update(controlDict)
        varyList = phaseVary+hapVary+histVary
        if not ihst:
            saveVaryList = varyList[:]
            for i,item in enumerate(saveVaryList):
                items = item.split(':')
                if items[1]:
                    items[1] = ''
                item = ':'.join(items)
                saveVaryList[i] = item
            SeqResult['varyList'] = saveVaryList
        else:
            newVaryList = varyList[:]
            for i,item in enumerate(newVaryList):
                items = item.split(':')
                if items[1]:
                    items[1] = ''
                item = ':'.join(items)
                newVaryList[i] = item
            if newVaryList != SeqResult['varyList']:
                print newVaryList
                print SeqResult['varyList']
                print '**** ERROR - variable list for this histogram does not match previous'
                raise Exception
        parmDict.update(phaseDict)
        parmDict.update(hapDict)
        parmDict.update(histDict)
        GetFprime(calcControls,Histo)
        constrDict,constrFlag,fixedList = G2mv.InputParse([])        #constraints go here?
        groups,parmlist = G2mv.GroupConstraints(constrDict)
        G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,constrFlag,fixedList)
        G2mv.Map2Dict(parmDict,varyList)
    
        while True:
            begin = time.time()
            values =  np.array(Dict2Values(parmDict, varyList))
            Ftol = Controls['min dM/M']
            Factor = Controls['shift factor']
            maxCyc = Controls['max cyc']

            if 'Jacobian' in Controls['deriv type']:            
                result = so.leastsq(errRefine,values,Dfun=dervRefine,full_output=True,
                    ftol=Ftol,col_deriv=True,factor=Factor,
                    args=([Histograms,Phases],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = int(result[2]['nfev']/2)
            elif 'Hessian' in Controls['deriv type']:
                result = G2mth.HessianLSQ(errRefine,values,Hess=HessRefine,ftol=Ftol,maxcyc=maxCyc,
                    args=([Histograms,Phases],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = result[2]['num cyc']+1                           
            else:           #'numeric'
                result = so.leastsq(errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                    args=([Histograms,Phases],parmDict,varyList,calcControls,pawleyLookup,dlg))
                ncyc = int(result[2]['nfev']/len(varyList))



            runtime = time.time()-begin
            chisq = np.sum(result[2]['fvec']**2)
            Values2Dict(parmDict, varyList, result[0])
            G2mv.Dict2Map(parmDict,varyList)
            
            Rwp = np.sqrt(chisq/Histo['sumwYo'])*100.      #to %
            GOF = chisq/(Histo['Nobs']-len(varyList))
            print '\n Refinement results for histogram: v'+histogram
            print 135*'-'
            print ' Number of function calls:',result[2]['nfev'],' Number of observations: ',Histo['Nobs'],' Number of parameters: ',len(varyList)
            print ' Refinement time = %8.3fs, %8.3fs/cycle, for %d cycles'%(runtime,runtime/ncyc,ncyc)
            print ' wRp = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f'%(Rwp,chisq,GOF)
            try:
                covMatrix = result[1]*GOF
                sig = np.sqrt(np.diag(covMatrix))
                if np.any(np.isnan(sig)):
                    print '*** Least squares aborted - some invalid esds possible ***'
                    ifPrint = True
                break                   #refinement succeeded - finish up!
            except TypeError:          #result[1] is None on singular matrix
                print '**** Refinement failed - singular matrix ****'
                if 'Hessian' in Controls['deriv type']:
                    for i in result[2]['psing'].reverse():
                            print 'Removing parameter: ',varyList[i]
                            del(varyList[i])                    
                else:
                    Ipvt = result[2]['ipvt']
                    for i,ipvt in enumerate(Ipvt):
                        if not np.sum(result[2]['fjac'],axis=1)[i]:
                            print 'Removing parameter: ',varyList[ipvt-1]
                            del(varyList[ipvt-1])
                            break
    
        GetFobsSq(Histo,Phases,parmDict,calcControls)
        sigDict = dict(zip(varyList,sig))
        newCellDict = GetNewCellParms(parmDict,varyList)
        newAtomDict = ApplyXYZshifts(parmDict,varyList)
        covData = {'variables':result[0],'varyList':varyList,'sig':sig,
            'covMatrix':covMatrix,'title':histogram,'newAtomDict':newAtomDict,'newCellDict':newCellDict}
        SetHistogramPhaseData(parmDict,sigDict,Phases,Histo,ifPrint)
        SetHistogramData(parmDict,sigDict,Histo,ifPrint)
        SeqResult[histogram] = covData
        G2IO.SetUsedHistogramsAndPhases(GPXfile,Histo,Phases,covData,makeBack)
        makeBack = False
    G2IO.SetSeqResult(GPXfile,Histograms,SeqResult)

def DistAngle(DisAglCtls,DisAglData):
    import numpy.ma as ma
    
    def ShowBanner(name):
        print 80*'*'
        print '   Interatomic Distances and Angles for phase '+name
        print 80*'*','\n'

    ShowBanner(DisAglCtls['Name'])
    SGData = DisAglData['SGData']
    SGtext = G2spc.SGPrint(SGData)
    for line in SGtext: print line
    Cell = DisAglData['Cell']
    
    Amat,Bmat = G2lat.cell2AB(Cell[:6])
    covData = {}
    if 'covData' in DisAglData:   
        covData = DisAglData['covData']
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        pfx = str(DisAglData['pId'])+'::'
        A = G2lat.cell2A(Cell[:6])
        cellSig = getCellEsd(pfx,SGData,A,covData)
        names = [' a = ',' b = ',' c = ',' alpha = ',' beta = ',' gamma = ',' Volume = ']
        valEsd = [G2mth.ValEsd(Cell[i],cellSig[i],True) for i in range(7)]
        line = '\n Unit cell:'
        for name,vals in zip(names,valEsd):
            line += name+vals  
        print line
    else: 
        print '\n Unit cell: a = ','%.5f'%(Cell[0]),' b = ','%.5f'%(Cell[1]),' c = ','%.5f'%(Cell[2]), \
            ' alpha = ','%.3f'%(Cell[3]),' beta = ','%.3f'%(Cell[4]),' gamma = ', \
            '%.3f'%(Cell[5]),' volume = ','%.3f'%(Cell[6])
    Factor = DisAglCtls['Factors']
    Radii = dict(zip(DisAglCtls['AtomTypes'],zip(DisAglCtls['BondRadii'],DisAglCtls['AngleRadii'])))
    Units = np.array([                   #is there a nicer way to make this?
        [-1,-1,-1],[-1,-1,0],[-1,-1,1],[-1,0,-1],[-1,0,0],[-1,0,1],[-1,1,-1],[-1,1,0],[-1,1,1],
        [0,-1,-1],[0,-1,0],[0,-1,1],[0,0,-1],[0,0,0],[0,0,1],[0,1,-1],[0,1,0],[0,1,1],
        [1,-1,-1],[1,-1,0],[1,-1,1],[1,0,-1],[1,0,0],[1,0,1],[1,1,-1],[1,1,0],[1,1,1]])
    origAtoms = DisAglData['OrigAtoms']
    targAtoms = DisAglData['TargAtoms']
    for Oatom in origAtoms:
        OxyzNames = ''
        IndBlist = []
        Dist = []
        Vect = []
        VectA = []
        angles = []
        for Tatom in targAtoms:
            Xvcov = []
            TxyzNames = ''
            if 'covData' in DisAglData:
                OxyzNames = [pfx+'dAx:%d'%(Oatom[0]),pfx+'dAy:%d'%(Oatom[0]),pfx+'dAz:%d'%(Oatom[0])]
                TxyzNames = [pfx+'dAx:%d'%(Tatom[0]),pfx+'dAy:%d'%(Tatom[0]),pfx+'dAz:%d'%(Tatom[0])]
                Xvcov = G2mth.getVCov(OxyzNames+TxyzNames,varyList,covMatrix)
            result = G2spc.GenAtom(Tatom[3:6],SGData,False)
            BsumR = (Radii[Oatom[2]][0]+Radii[Tatom[2]][0])*Factor[0]
            AsumR = (Radii[Oatom[2]][1]+Radii[Tatom[2]][1])*Factor[1]
            for Txyz,Top,Tunit in result:
                Dx = (Txyz-np.array(Oatom[3:6]))+Units
                dx = np.inner(Amat,Dx)
                dist = ma.masked_less(np.sqrt(np.sum(dx**2,axis=0)),0.5)
                IndB = ma.nonzero(ma.masked_greater(dist-BsumR,0.))
                if np.any(IndB):
                    for indb in IndB:
                        for i in range(len(indb)):
                            if str(dx.T[indb][i]) not in IndBlist:
                                IndBlist.append(str(dx.T[indb][i]))
                                unit = Units[indb][i]
                                tunit = '[%2d%2d%2d]'%(unit[0]+Tunit[0],unit[1]+Tunit[1],unit[2]+Tunit[2])
                                pdpx = G2mth.getDistDerv(Oatom[3:6],Tatom[3:6],Amat,unit,Top,SGData)
                                sig = 0.0
                                if len(Xvcov):
                                    sig = np.sqrt(np.inner(pdpx,np.inner(Xvcov,pdpx)))
                                Dist.append([Oatom[1],Tatom[1],tunit,Top,ma.getdata(dist[indb])[i],sig])
                                if (Dist[-1][-1]-AsumR) <= 0.:
                                    Vect.append(dx.T[indb][i]/Dist[-1][-2])
                                    VectA.append([OxyzNames,np.array(Oatom[3:6]),TxyzNames,np.array(Tatom[3:6]),unit,Top])
                                else:
                                    Vect.append([0.,0.,0.])
                                    VectA.append([])
        Vect = np.array(Vect)
        angles = np.zeros((len(Vect),len(Vect)))
        angsig = np.zeros((len(Vect),len(Vect)))
        for i,veca in enumerate(Vect):
            if np.any(veca):
                for j,vecb in enumerate(Vect):
                    if np.any(vecb):
                        angles[i][j],angsig[i][j] = G2mth.getAngSig(VectA[i],VectA[j],Amat,SGData,covData)
        line = ''
        for i,x in enumerate(Oatom[3:6]):
            if len(Xvcov):
                line += '%12s'%(G2mth.ValEsd(x,np.sqrt(Xvcov[i][i]),True))
            else:
                line += '%12.5f'%(x)
        print '\n Distances & angles for ',Oatom[1],' at ',line
        print 80*'*'
        line = ''
        for dist in Dist[:-1]:
            line += '%12s'%(dist[1].center(12))
        print '  To       cell +(sym. op.)      dist.  ',line
        for i,dist in enumerate(Dist):
            line = ''
            for j,angle in enumerate(angles[i][0:i]):
                sig = angsig[i][j]
                if angle:
                    if sig:
                        line += '%12s'%(G2mth.ValEsd(angle,sig,True).center(12))
                    else:
                        val = '%.3f'%(angle)
                        line += '%12s'%(val.center(12))
                else:
                    line += 12*' '
            if dist[5]:            #sig exists!
                val = G2mth.ValEsd(dist[4],dist[5])
            else:
                val = '%8.4f'%(dist[4])
            print '  %8s%10s+(%4d) %12s'%(dist[1].ljust(8),dist[2].ljust(10),dist[3],val.center(12)),line

def Torsion(TorsionData):

    def ShowBanner(name):
        print 80*'*'
        print '   Torsion angle for phase '+name
        print 80*'*','\n'

    ShowBanner(TorsionData['Name'])
    SGData = TorsionData['SGData']
    SGtext = G2spc.SGPrint(SGData)
    for line in SGtext: print line
    Cell = TorsionData['Cell']
    
    Amat,Bmat = G2lat.cell2AB(Cell[:6])
    covData = {}
    if 'covData' in TorsionData:   
        covData = TorsionData['covData']
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        pfx = str(TorsionData['pId'])+'::'
        A = G2lat.cell2A(Cell[:6])
        cellSig = getCellEsd(pfx,SGData,A,covData)
        names = [' a = ',' b = ',' c = ',' alpha = ',' beta = ',' gamma = ',' Volume = ']
        valEsd = [G2mth.ValEsd(Cell[i],cellSig[i],True) for i in range(7)]
        line = '\n Unit cell:'
        for name,vals in zip(names,valEsd):
            line += name+vals  
        print line
    else: 
        print '\n Unit cell: a = ','%.5f'%(Cell[0]),' b = ','%.5f'%(Cell[1]),' c = ','%.5f'%(Cell[2]), \
            ' alpha = ','%.3f'%(Cell[3]),' beta = ','%.3f'%(Cell[4]),' gamma = ', \
            '%.3f'%(Cell[5]),' volume = ','%.3f'%(Cell[6])
    #find one end of 4 atom string - involved in longest distance
    dist = {}
    for i,X1 in enumerate(TorsionData['Datoms']):
        for j,X2 in enumerate(TorsionData['Datoms'][:i]):
            dist[np.sqrt(np.sum(np.inner(Amat,np.array(X2[3:6])-np.array(X1[3:6]))**2))] = [i,j]
    sortdist = dist.keys()
    sortdist.sort()
    end = dist[sortdist[-1]][0]
    #order atoms in distance from end - defines sequence of atoms for the torsion angle
    dist = {}
    X1 = TorsionData['Datoms'][end]
    for i,X2 in enumerate(TorsionData['Datoms']):                
        dist[np.sqrt(np.sum(np.inner(Amat,np.array(X2[3:6])-np.array(X1[3:6]))**2))] = i
    sortdist = dist.keys()
    sortdist.sort()
    Datoms = []
    for d in sortdist:
        atom = TorsionData['Datoms'][dist[d]]
        symop = atom[-1].split('+')
        if len(symop) == 1:
            symop.append('0,0,0')        
        symop[0] = int(symop[0])
        symop[1] = eval(symop[1])
        atom[-1] = symop
        print atom
        Datoms.append(atom)
    Tors,sig = G2mth.GetTorsionSig(Datoms,Amat,SGData,covData={})
    print ' Torsion: ',G2mth.ValEsd(Tors,sig)
        
        
def BestPlane(PlaneData):

    def ShowBanner(name):
        print 80*'*'
        print '   Best plane result for phase '+name
        print 80*'*','\n'

    ShowBanner(PlaneData['Name'])

    Cell = PlaneData['Cell']    
    Amat,Bmat = G2lat.cell2AB(Cell[:6])        
    Atoms = PlaneData['Atoms']
    sumXYZ = np.zeros(3)
    XYZ = []
    Natoms = len(Atoms)
    for atom in Atoms:
        xyz = np.array(atom[3:6])
        XYZ.append(xyz)
        sumXYZ += xyz
    sumXYZ /= Natoms
    XYZ = np.array(XYZ)-sumXYZ
    XYZ = np.inner(Amat,XYZ).T
    Zmat = np.zeros((3,3))
    for i,xyz in enumerate(XYZ):
        Zmat += np.outer(xyz.T,xyz)
    print ' Selected atoms centered at %10.5f %10.5f %10.5f'%(sumXYZ[0],sumXYZ[1],sumXYZ[2])
    Evec,Emat = nl.eig(Zmat)
    Evec = np.sqrt(Evec)/(Natoms-3)
    Order = np.argsort(Evec)
    XYZ = np.inner(XYZ,Emat.T).T
    XYZ = np.array([XYZ[Order[2]],XYZ[Order[1]],XYZ[Order[0]]]).T
    print ' Atoms in Cartesian best plane coordinates:'
    print ' Name         X         Y         Z'
    for i,xyz in enumerate(XYZ):
        print ' %6s%10.3f%10.3f%10.3f'%(Atoms[i][1].ljust(6),xyz[0],xyz[1],xyz[2])
    print '\n Best plane RMS X =%8.3f, Y =%8.3f, Z =%8.3f'%(Evec[Order[2]],Evec[Order[1]],Evec[Order[0]])   

            
def main():
    arg = sys.argv
    if len(arg) > 1:
        GPXfile = arg[1]
        if not ospath.exists(GPXfile):
            print 'ERROR - ',GPXfile," doesn't exist!"
            exit()
        GPXpath = ospath.dirname(arg[1])
        Refine(GPXfile,None)
    else:
        print 'ERROR - missing filename'
        exit()
    print "Done"
         
if __name__ == '__main__':
    main()
