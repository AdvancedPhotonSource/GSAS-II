# -*- coding: utf-8 -*-
'''
*GSASIImpsubs - routines used in multiprocessing*
-------------------------------------------------

The routines here are called either directly when GSAS-II is used without multiprocessing
or in separate cores when multiprocessing is used.

These routines are designed to be used in one of two ways:

 * when multiprocessing is
   enabled (see global variable useMP) the computational routines are called in
   separate Python interpreter that is created and then deleted after use.

 * when useMP is False, these routines are called directly from the main "thread".

Note that :func:`GSASIImpsubs.InitMP` should be called before any of the other routines
in this module are used. 
'''
########### SVN repository information ###################
# $Date: $
# $Author: $
# $Revision: $
# $URL: $
# $Id: $
########### SVN repository information ###################
import multiprocessing as mp
import numpy as np
import numpy.ma as ma
import numpy.linalg as nl
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 2895 $")
import GSASIIpwd as G2pwd
import GSASIIstrMath as G2stMth

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
#asind = lambda x: 180.*np.arcsin(x)/np.pi
#acosd = lambda x: 180.*np.arccos(x)/np.pi
#atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
    
ncores = None

def InitMP(allowMP=True):
    '''Called in routines to initialize use of Multiprocessing
    '''
    global useMP,ncores
    if ncores is not None: return
    useMP = False
    if not allowMP:
        print('Multiprocessing disabled')
        ncores = 0
        return
    ncores = GSASIIpath.GetConfigValue('Multiprocessing_cores',-1)
    if ncores < 0: ncores = mp.cpu_count()
    if ncores > 1:
        useMP = True
    #if GSASIIpath.GetConfigValue('debug'):
    if True:
        print('Multiprocessing with {} cores enabled'.format(ncores))

################################################################################
# derivative computation
################################################################################        
def InitDerivGlobals(im1,calcControls1,SGMT1,hfx1,phfx1,pfx1,G1,GB1,g1,SGData1,
                     parmDict1,wave1,shl1,x1,cw1,Ka21,A1,varylist1,dependentVars1,
                     dFdvDict1,lamRatio1,kRatio1,doPawley1,pawleyLookup1):
    '''Initialize for the computation of derivatives. Puts lots of junk into the global
    namespace in this module, including the arrays for derivatives (when needed.)
    '''
    global im,calcControls,SGMT,hfx,phfx,pfx,G,GB,g,SGData,parmDict,wave,shl,x,cw,Ka2,A
    global varylist,dependentVars,dFdvDict,lamRatio,kRatio,doPawley,pawleyLookup
    im = im1
    calcControls = calcControls1
    SGMT = SGMT1
    hfx = hfx1
    phfx = phfx1
    pfx = pfx1
    G = G1
    GB = GB1
    g = g1
    SGData = SGData1
    parmDict = parmDict1
    wave = wave1
    shl = shl1
    x = ma.getdata(x1)
    cw = cw1
    Ka2 = Ka21
    A = A1
    varylist = varylist1
    dependentVars = dependentVars1
    dFdvDict = dFdvDict1
    lamRatio = lamRatio1
    kRatio = kRatio1
    doPawley = doPawley1
    pawleyLookup = pawleyLookup1
    # determine the parameters that will have derivatives computed only at end
    global nonatomvarylist
    nonatomvarylist = []
    for name in varylist:
        if '::RBV;' not in name:
            try:
                aname = name.split(pfx)[1][:2]
                if aname not in ['Af','dA','AU','RB','AM','Xs','Xc','Ys','Yc','Zs','Zc',    \
                    'Tm','Xm','Ym','Zm','U1','U2','U3']: continue # skip anything not an atom or rigid body param
            except IndexError:
                continue
        nonatomvarylist.append(name)
    global nonatomdependentVars
    nonatomdependentVars = []
    for name in dependentVars:
        if '::RBV;' not in name:
            try:
                aname = name.split(pfx)[1][:2]
                if aname not in ['Af','dA','AU','RB','AM','Xs','Xc','Ys','Yc','Zs','Zc',    \
                    'Tm','Xm','Ym','Zm','U1','U2','U3']: continue # skip anything not an atom or rigid body param
            except IndexError:
                continue
        nonatomdependentVars.append(name)
    # create local copies of derivative arrays, if multiprocessing will be used
    if useMP:
        global dMdv
        dMdv = np.zeros(shape=(len(varylist),len(x)))
        global depDerivDict
        depDerivDict = {j:np.zeros(shape=len(x)) for j in dependentVars}
            

def cellVaryDerv(pfx,SGData,dpdA): 
    if SGData['SGLaue'] in ['-1',]:
        return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],
            [pfx+'A3',dpdA[3]],[pfx+'A4',dpdA[4]],[pfx+'A5',dpdA[5]]]
    elif SGData['SGLaue'] in ['2/m',]:
        if SGData['SGUniq'] == 'a':
            return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A5',dpdA[5]]]
        elif SGData['SGUniq'] == 'b':
            return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A4',dpdA[4]]]
        else:
            return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]],[pfx+'A3',dpdA[3]]]
    elif SGData['SGLaue'] in ['mmm',]:
        return [[pfx+'A0',dpdA[0]],[pfx+'A1',dpdA[1]],[pfx+'A2',dpdA[2]]]
    elif SGData['SGLaue'] in ['4/m','4/mmm']:
        return [[pfx+'A0',dpdA[0]],[pfx+'A2',dpdA[2]]]
    elif SGData['SGLaue'] in ['6/m','6/mmm','3m1', '31m', '3']:
        return [[pfx+'A0',dpdA[0]],[pfx+'A2',dpdA[2]]]
    elif SGData['SGLaue'] in ['3R', '3mR']:
        return [[pfx+'A0',dpdA[0]+dpdA[1]+dpdA[2]],[pfx+'A3',dpdA[3]+dpdA[4]+dpdA[5]]]                       
    elif SGData['SGLaue'] in ['m3m','m3']:
        return [[pfx+'A0',dpdA[0]]]

def ComputeDerivMPbatch(reflsList):
    '''Computes a the derivatives for a batch of reflections and sums the them into
    global arrays dMdv & depDerivDict. These arrays are returned once the computation
    is completed.
    '''
    for refl,iref,fmin,fmax,iBeg,iFin in reflsList:
        if ComputeDeriv(refl,iref,fmin,fmax,iBeg,iFin,dMdv,depDerivDict): break
    return dMdv,depDerivDict

def ComputeDeriv(refl,iref,fmin,fmax,iBeg,iFin,dMdv,depDerivDict):
    '''Compute the parameter derivatives for a single reflection and add the results
    into either array dMdv or depDerivDict
    '''
    global wave
    if im:
        h,k,l,m = refl[:4]
    else:
        h,k,l = refl[:3]
    Uniq = np.inner(refl[:3],SGMT)
    if 'T' in calcControls[hfx+'histType']:
        wave = refl[14+im]
    dIdsh,dIdsp,dIdpola,dIdPO,dFdODF,dFdSA,dFdAb,dFdEx = G2stMth.GetIntensityDerv(refl,im,wave,Uniq,G,g,pfx,phfx,hfx,SGData,calcControls,parmDict)
    pos = refl[5+im]
    calcKa2 = False
    if 'C' in calcControls[hfx+'histType']:
        tanth = tand(pos/2.0)
        costh = cosd(pos/2.0)
        if Ka2:
            pos2 = refl[5+im]+lamRatio*tanth       # + 360/pi * Dlam/lam * tan(th)
            iBeg2 = np.searchsorted(x,pos2-fmin)
            iFin2 = np.searchsorted(x,pos2+fmax)
            if iBeg2-iFin2:
                calcKa2 = True
                iFin = iFin2        
        lenBF = iFin-iBeg
        dMdpk = np.zeros(shape=(6,lenBF))
        dMdipk = G2pwd.getdFCJVoigt3(refl[5+im],refl[6+im],refl[7+im],shl,x[iBeg:iFin])
        for i in range(5):
            dMdpk[i] += 100.*cw[iBeg:iFin]*refl[11+im]*refl[9+im]*dMdipk[i]
        dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'sig':dMdpk[2],'gam':dMdpk[3],'shl':dMdpk[4],'L1/L2':np.zeros_like(dMdpk[0])}
        if calcKa2: 
            dMdpk2 = np.zeros(shape=(6,lenBF))
            dMdipk2 = G2pwd.getdFCJVoigt3(pos2,refl[6+im],refl[7+im],shl,x[iBeg:iFin])
            for i in range(5):
                dMdpk2[i] = 100.*cw[iBeg:iFin]*refl[11+im]*refl[9+im]*kRatio*dMdipk2[i]
            dMdpk2[5] = 100.*cw[iBeg:iFin]*refl[11+im]*dMdipk2[0]
            dervDict2 = {'int':dMdpk2[0],'pos':dMdpk2[1],'sig':dMdpk2[2],'gam':dMdpk2[3],'shl':dMdpk2[4],'L1/L2':dMdpk2[5]*refl[9]}
    else:   #'T'OF
        lenBF = iFin-iBeg
        if lenBF < 0: return True  #bad peak coeff
        dMdpk = np.zeros(shape=(6,lenBF))
        dMdipk = G2pwd.getdEpsVoigt(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im],x[iBeg:iFin])
        for i in range(6):
            dMdpk[i] += refl[11+im]*refl[9+im]*dMdipk[i]      #cw[iBeg:iFin]*
        dervDict = {'int':dMdpk[0],'pos':dMdpk[1],'alp':dMdpk[2],'bet':dMdpk[3],'sig':dMdpk[4],'gam':dMdpk[5]}            
    if doPawley:
        try:
            if im:
                pIdx = pfx+'PWLref:'+str(pawleyLookup[pfx+'%d,%d,%d,%d'%(h,k,l,m)])
            else:
                pIdx = pfx+'PWLref:'+str(pawleyLookup[pfx+'%d,%d,%d'%(h,k,l)])
            idx = varylist.index(pIdx)
            dMdv[idx][iBeg:iFin] = dervDict['int']/refl[9+im]
            if Ka2: #not for TOF either
                dMdv[idx][iBeg:iFin] += dervDict2['int']/refl[9+im]
        except: # ValueError:
            pass
    if 'C' in calcControls[hfx+'histType']:
        dpdA,dpdw,dpdZ,dpdSh,dpdTr,dpdX,dpdY,dpdV = G2stMth.GetReflPosDerv(refl,im,wave,A,pfx,hfx,calcControls,parmDict)
        names = {hfx+'Scale':[dIdsh,'int'],hfx+'Polariz.':[dIdpola,'int'],phfx+'Scale':[dIdsp,'int'],
            hfx+'U':[tanth**2,'sig'],hfx+'V':[tanth,'sig'],hfx+'W':[1.0,'sig'],
            hfx+'X':[1.0/costh,'gam'],hfx+'Y':[tanth,'gam'],hfx+'SH/L':[1.0,'shl'],
            hfx+'I(L2)/I(L1)':[1.0,'L1/L2'],hfx+'Zero':[dpdZ,'pos'],hfx+'Lam':[dpdw,'pos'],
            hfx+'Shift':[dpdSh,'pos'],hfx+'Transparency':[dpdTr,'pos'],hfx+'DisplaceX':[dpdX,'pos'],
            hfx+'DisplaceY':[dpdY,'pos'],}
        if 'Bragg' in calcControls[hfx+'instType']:
            names.update({hfx+'SurfRoughA':[dFdAb[0],'int'],
                hfx+'SurfRoughB':[dFdAb[1],'int'],})
        else:
            names.update({hfx+'Absorption':[dFdAb,'int'],})
    else:   #'T'OF
        dpdA,dpdZ,dpdDC,dpdDA,dpdDB,dpdV = G2stMth.GetReflPosDerv(refl,im,0.0,A,pfx,hfx,calcControls,parmDict)
        names = {hfx+'Scale':[dIdsh,'int'],phfx+'Scale':[dIdsp,'int'],
            hfx+'difC':[dpdDC,'pos'],hfx+'difA':[dpdDA,'pos'],hfx+'difB':[dpdDB,'pos'],
            hfx+'Zero':[dpdZ,'pos'],hfx+'X':[refl[4+im],'gam'],hfx+'Y':[refl[4+im]**2,'gam'],
            hfx+'alpha':[1./refl[4+im],'alp'],hfx+'beta-0':[1.0,'bet'],hfx+'beta-1':[1./refl[4+im]**4,'bet'],
            hfx+'beta-q':[1./refl[4+im]**2,'bet'],hfx+'sig-0':[1.0,'sig'],hfx+'sig-1':[refl[4+im]**2,'sig'],
            hfx+'sig-2':[refl[4+im]**4,'sig'],hfx+'sig-q':[1./refl[4+im]**2,'sig'],
            hfx+'Absorption':[dFdAb,'int'],phfx+'Extinction':[dFdEx,'int'],}
    for name in names:
        item = names[name]
        if name in varylist:
            dMdv[varylist.index(name)][iBeg:iFin] += item[0]*dervDict[item[1]]
            if calcKa2:
                dMdv[varylist.index(name)][iBeg:iFin] += item[0]*dervDict2[item[1]]
        elif name in dependentVars:
            depDerivDict[name][iBeg:iFin] += item[0]*dervDict[item[1]]
            if calcKa2:
                depDerivDict[name][iBeg:iFin] += item[0]*dervDict2[item[1]]
    for iPO in dIdPO:
        if iPO in varylist:
            dMdv[varylist.index(iPO)][iBeg:iFin] += dIdPO[iPO]*dervDict['int']
            if calcKa2:
                dMdv[varylist.index(iPO)][iBeg:iFin] += dIdPO[iPO]*dervDict2['int']
        elif iPO in dependentVars:
            depDerivDict[iPO][iBeg:iFin] += dIdPO[iPO]*dervDict['int']
            if calcKa2:
                depDerivDict[iPO][iBeg:iFin] += dIdPO[iPO]*dervDict2['int']
    for i,name in enumerate(['omega','chi','phi']):
        aname = pfx+'SH '+name
        if aname in varylist:
            dMdv[varylist.index(aname)][iBeg:iFin] += dFdSA[i]*dervDict['int']
            if calcKa2:
                dMdv[varylist.index(aname)][iBeg:iFin] += dFdSA[i]*dervDict2['int']
        elif aname in dependentVars:
            depDerivDict[aname][iBeg:iFin] += dFdSA[i]*dervDict['int']
            if calcKa2:
                depDerivDict[aname][iBeg:iFin] += dFdSA[i]*dervDict2['int']
    for iSH in dFdODF:
        if iSH in varylist:
            dMdv[varylist.index(iSH)][iBeg:iFin] += dFdODF[iSH]*dervDict['int']
            if calcKa2:
                dMdv[varylist.index(iSH)][iBeg:iFin] += dFdODF[iSH]*dervDict2['int']
        elif iSH in dependentVars:
            depDerivDict[iSH][iBeg:iFin] += dFdODF[iSH]*dervDict['int']
            if calcKa2:
                depDerivDict[iSH][iBeg:iFin] += dFdODF[iSH]*dervDict2['int']
    cellDervNames = cellVaryDerv(pfx,SGData,dpdA)
    for name,dpdA in cellDervNames:
        if name in varylist:
            dMdv[varylist.index(name)][iBeg:iFin] += dpdA*dervDict['pos']
            if calcKa2:
                dMdv[varylist.index(name)][iBeg:iFin] += dpdA*dervDict2['pos']
        elif name in dependentVars: #need to scale for mixed phase constraints?
            depDerivDict[name][iBeg:iFin] += dpdA*dervDict['pos']
            if calcKa2:
                depDerivDict[name][iBeg:iFin] += dpdA*dervDict2['pos']
    dDijDict = G2stMth.GetHStrainShiftDerv(refl,im,SGData,phfx,hfx,calcControls,parmDict)
    for name in dDijDict:
        if name in varylist:
            dMdv[varylist.index(name)][iBeg:iFin] += dDijDict[name]*dervDict['pos']
            if calcKa2:
                dMdv[varylist.index(name)][iBeg:iFin] += dDijDict[name]*dervDict2['pos']
        elif name in dependentVars:
            depDerivDict[name][iBeg:iFin] += dDijDict[name]*dervDict['pos']
            if calcKa2:
                depDerivDict[name][iBeg:iFin] += dDijDict[name]*dervDict2['pos']
    for i,name in enumerate([pfx+'mV0',pfx+'mV1',pfx+'mV2']):
        if name in varylist:
            dMdv[varylist.index(name)][iBeg:iFin] += dpdV[i]*dervDict['pos']
            if calcKa2:
                dMdv[varylist.index(name)][iBeg:iFin] += dpdV[i]*dervDict2['pos']
        elif name in dependentVars:
            depDerivDict[name][iBeg:iFin] += dpdV[i]*dervDict['pos']
            if calcKa2:
                depDerivDict[name][iBeg:iFin] += dpdV[i]*dervDict2['pos']
    if 'C' in calcControls[hfx+'histType']:
        sigDict,gamDict = G2stMth.GetSampleSigGamDerv(refl,im,wave,G,GB,SGData,hfx,phfx,calcControls,parmDict)
    else:   #'T'OF
        sigDict,gamDict = G2stMth.GetSampleSigGamDerv(refl,im,0.0,G,GB,SGData,hfx,phfx,calcControls,parmDict)
    for name in gamDict:
        if name in varylist:
            dMdv[varylist.index(name)][iBeg:iFin] += gamDict[name]*dervDict['gam']
            if calcKa2:
                dMdv[varylist.index(name)][iBeg:iFin] += gamDict[name]*dervDict2['gam']
        elif name in dependentVars:
            depDerivDict[name][iBeg:iFin] += gamDict[name]*dervDict['gam']
            if calcKa2:
                depDerivDict[name][iBeg:iFin] += gamDict[name]*dervDict2['gam']
    for name in sigDict:
        if name in varylist:
            dMdv[varylist.index(name)][iBeg:iFin] += sigDict[name]*dervDict['sig']
            if calcKa2:
                dMdv[varylist.index(name)][iBeg:iFin] += sigDict[name]*dervDict2['sig']
        elif name in dependentVars:
            depDerivDict[name][iBeg:iFin] += sigDict[name]*dervDict['sig']
            if calcKa2:
                depDerivDict[name][iBeg:iFin] += sigDict[name]*dervDict2['sig']
    for name in ['BabA','BabU']:
        if refl[9+im]:
            if phfx+name in varylist:
                dMdv[varylist.index(phfx+name)][iBeg:iFin] += parmDict[phfx+'Scale']*dFdvDict[phfx+name][iref]*dervDict['int']/refl[9+im]
                if calcKa2:
                    dMdv[varylist.index(phfx+name)][iBeg:iFin] += parmDict[phfx+'Scale']*dFdvDict[phfx+name][iref]*dervDict2['int']/refl[9+im]
            elif phfx+name in dependentVars:                    
                depDerivDict[phfx+name][iBeg:iFin] += parmDict[phfx+'Scale']*dFdvDict[phfx+name][iref]*dervDict['int']/refl[9+im]
                if calcKa2:
                    depDerivDict[phfx+name][iBeg:iFin] += parmDict[phfx+'Scale']*dFdvDict[phfx+name][iref]*dervDict2['int']/refl[9+im]                  
    if not doPawley and not parmDict[phfx+'LeBail']:
        #do atom derivatives -  for RB,F,X & U so far - how do I scale mixed phase constraints?
        corr = 0.
        #corr2 = 0.
        if refl[9+im]:             
            corr = dervDict['int']/refl[9+im]
            #if calcKa2:  # commented out in Bob's code. Why?
            #    corr2 = dervDict2['int']/refl[9+im]
        for name in nonatomvarylist:
            dMdv[varylist.index(name)][iBeg:iFin] += dFdvDict[name][iref]*corr
            #if calcKa2:
            #   dMdv[varylist.index(name)][iBeg:iFin] += dFdvDict[name][iref]*corr2 # unneeded w/o above
        for name in nonatomdependentVars:
           depDerivDict[name][iBeg:iFin] += dFdvDict[name][iref]*corr
           #if calcKa2:
           #    depDerivDict[name][iBeg:iFin] += dFdvDict[name][iref]*corr2

################################################################################
# Fobs Squared computation
################################################################################        
#x,ratio,shl,xB,xF,im,lamRatio,kRatio,xMask,Ka2
def InitFobsSqGlobals(x1,ratio1,shl1,xB1,xF1,im1,lamRatio1,kRatio1,xMask1,Ka21):
    '''Initialize for the computation of Fobs Squared for powder histograms.
    Puts lots of junk into the global namespace in this module.
    '''
    global x,ratio,shl,xB,xF,im,lamRatio,kRatio,xMask,Ka2
    x = ma.getdata(x1)
    ratio = ratio1
    shl = shl1
    xB = xB1
    xF = xF1
    im = im1
    lamRatio = lamRatio1
    kRatio = kRatio1
    xMask = xMask1
    Ka2 = Ka21

def ComputeFobsSqCWbatch(profList):
    sInt = 0
    resList = []
    for refl,iref in profList:
        icod = ComputeFobsSqCW(refl,iref)
        if type(icod) is tuple:
            resList.append((icod[0],iref))
            sInt += icod[1]
        elif icod == -1:
            res.append((None,iref))
        elif icod == -2:
            break
    return sInt,resList

def ComputeFobsSqTOFbatch(profList):
    sInt = 0
    resList = []
    for refl,iref in profList:
        icod = ComputeFobsSqTOF(refl,iref)
        if type(icod) is tuple:
            resList.append((icod[0],iref))
            sInt += icod[1]
        elif icod == -1:
            res.append((None,iref))
        elif icod == -2:
            break
    return sInt,resList
        
def ComputeFobsSqCW(refl,iref):
    yp = np.zeros(len(x)) # not masked
    sInt = 0
    refl8im = 0
    Wd,fmin,fmax = G2pwd.getWidthsCW(refl[5+im],refl[6+im],refl[7+im],shl)
    iBeg = max(xB,np.searchsorted(x,refl[5+im]-fmin))
    iFin = max(xB,min(np.searchsorted(x,refl[5+im]+fmax),xF))
    iFin2 = iFin
    if not iBeg+iFin:       #peak below low limit - skip peak
        return 0
    if ma.all(xMask[iBeg:iFin]):    #peak entirely masked - skip peak
        return -1
    elif not iBeg-iFin:     #peak above high limit - done
        return -2
    elif iBeg < iFin:
        yp[iBeg:iFin] = refl[11+im]*refl[9+im]*G2pwd.getFCJVoigt3(
            refl[5+im],refl[6+im],refl[7+im],shl,x[iBeg:iFin])
        sInt = refl[11+im]*refl[9+im]
        if Ka2:
            pos2 = refl[5+im]+lamRatio*tand(refl[5+im]/2.0)       # + 360/pi * Dlam/lam * tan(th)
            Wd,fmin,fmax = G2pwd.getWidthsCW(pos2,refl[6+im],refl[7+im],shl)
            iBeg2 = max(xB,np.searchsorted(x,pos2-fmin))
            iFin2 = min(np.searchsorted(x,pos2+fmax),xF)
            if iFin2 > iBeg2: 
                yp[iBeg2:iFin2] += refl[11+im]*refl[9+im]*kRatio*G2pwd.getFCJVoigt3(
                    pos2,refl[6+im],refl[7+im],shl,x[iBeg2:iFin2])
                sInt *= 1.+kRatio
    refl8im = np.sum(np.where(ratio[iBeg:iFin2]>0.,yp[iBeg:iFin2]*ratio[iBeg:iFin2]/(refl[11+im]*(1.+kRatio)),0.0))
    return refl8im,sInt

def ComputeFobsSqTOF(refl,iref):
    yp = np.zeros(len(x)) # not masked
    refl8im = 0
    Wd,fmin,fmax = G2pwd.getWidthsTOF(refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im])
    iBeg = max(xB,np.searchsorted(x,refl[5+im]-fmin))
    iFin = max(xB,min(np.searchsorted(x,refl[5+im]+fmax),xF))
    if not iBeg+iFin:       #peak below low limit - skip peak
        return 0
    if ma.all(xMask[iBeg:iFin]):    #peak entirely masked - skip peak
        return -1
    elif not iBeg-iFin:     #peak above high limit - done
        return -2
    if iBeg < iFin:
        yp[iBeg:iFin] = refl[11+im]*refl[9+im]*G2pwd.getEpsVoigt(
            refl[5+im],refl[12+im],refl[13+im],refl[6+im],refl[7+im],x[iBeg:iFin])
    refl8im = np.sum(np.where(ratio[iBeg:iFin]>0.,yp[iBeg:iFin]*ratio[iBeg:iFin]/refl[11+im],0.0))
    return refl8im,refl[11+im]*refl[9+im]
################################################################################
# Powder Profile computation
################################################################################        
def InitPwdrProfGlobals(im1,shl1,x1):
    '''Initialize for the computation of Fobs Squared for powder histograms.
    Puts lots of junk into the global namespace in this module.
    '''
    global im,shl,x
    im = im1
    shl = shl1
    x = ma.getdata(x1)
    global cw
    cw = np.diff(x)
    cw = np.append(cw,cw[-1])
    # create local copies of ycalc array
    if useMP:
        global yc
        yc = np.zeros_like(x1)


def ComputePwdrProfCW(profList):
    'Compute the peaks profile for a set of CW peaks and add into the yc array'
    for pos,refl,iBeg,iFin,kRatio in profList:
        yc[iBeg:iFin] += refl[11+im]*refl[9+im]*kRatio*G2pwd.getFCJVoigt3(
            pos,refl[6+im],refl[7+im],shl,x[iBeg:iFin])
    return yc

def ComputePwdrProfTOF(profList):
    'Compute the peaks profile for a set of TOF peaks and add into the yc array'
    for pos,refl,iBeg,iFin in profList:
        yc[iBeg:iFin] += refl[11+im]*refl[9+im]*G2pwd.getEpsVoigt(
            pos,refl[12+im],refl[13+im],refl[6+im],refl[7+im],x[iBeg:iFin])/cw[iBeg:iFin]
    return yc
    
