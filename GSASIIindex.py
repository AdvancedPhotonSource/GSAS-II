# -*- coding: utf-8 -*-
#GSASII cell indexing program: variation on that of A. Coehlo
#   includes cell refinement from peak positions (not zero as yet)
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import math
import wx
import time
import numpy as np
import numpy.linalg as nl
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision$")
import GSASIIlattice as G2lat
import scipy.optimize as so

# trig functions in degrees
sind = lambda x: math.sin(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
tand = lambda x: math.tan(x*math.pi/180.)
atand = lambda x: 180.*math.atan(x)/math.pi
atan2d = lambda y,x: 180.*math.atan2(y,x)/math.pi
cosd = lambda x: math.cos(x*math.pi/180.)
acosd = lambda x: 180.*math.acos(x)/math.pi
rdsq2d = lambda x,p: round(1.0/math.sqrt(x),p)
#numpy versions
npsind = lambda x: np.sin(x*np.pi/180.)
npasind = lambda x: 180.*np.arcsin(x)/math.pi
npcosd = lambda x: np.cos(x*math.pi/180.)
nptand = lambda x: np.tan(x*math.pi/180.)
npatand = lambda x: 180.*np.arctan(x)/np.pi
npatan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi
    
def scaleAbyV(A,V):
    v = G2lat.calc_V(A)
    scale = math.exp(math.log(v/V)/3.)**2
    for i in range(6):
        A[i] *= scale
    
def ranaxis(dmin,dmax):
    import random as rand
    return rand.random()*(dmax-dmin)+dmin
    
def ran2axis(k,N):
    import random as rand
    T = 1.5+0.49*k/N
#    B = 0.99-0.49*k/N
#    B = 0.99-0.049*k/N
    B = 0.99-0.149*k/N
    R = (T-B)*rand.random()+B
    return R
    
#def ranNaxis(k,N):
#    import random as rand
#    T = 1.0+1.0*k/N
#    B = 1.0-1.0*k/N
#    R = (T-B)*rand.random()+B
#    return R
    
def ranAbyV(Bravais,dmin,dmax,V):
    cell = [0,0,0,0,0,0]
    bad = True
    while bad:
        bad = False
        cell = rancell(Bravais,dmin,dmax)
        G,g = G2lat.cell2Gmat(cell)
        A = G2lat.Gmat2A(G)
        if G2lat.calc_rVsq(A) < 1:
            scaleAbyV(A,V)
            cell = G2lat.A2cell(A)
            for i in range(3):
                bad |= cell[i] < dmin
    return A
    
def ranAbyR(Bravais,A,k,N,ranFunc):
    R = ranFunc(k,N)
    if Bravais in [0,1,2]:          #cubic - not used
        A[0] = A[1] = A[2] = A[0]*R
        A[3] = A[4] = A[5] = 0.
    elif Bravais in [3,4]:          #hexagonal/trigonal
        A[0] = A[1] = A[3] = A[0]*R
        A[2] *= R
        A[4] = A[5] = 0.        
    elif Bravais in [5,6]:          #tetragonal
        A[0] = A[1] = A[0]*R
        A[2] *= R
        A[3] = A[4] = A[5] = 0.        
    elif Bravais in [7,8,9,10]:     #orthorhombic
        A[0] *= R
        A[1] *= R
        A[2] *= R
        A[3] = A[4] = A[5] = 0.        
    elif Bravais in [11,12]:        #monoclinic
        A[0] *= R
        A[1] *= R
        A[2] *= R
        A[4] *= R
        A[3] = A[5] = 0.        
    else:                           #triclinic
        A[0] *= R
        A[1] *= R
        A[2] *= R
        A[3] *= R
        A[4] *= R
        A[5] *= R
    return A
    
def rancell(Bravais,dmin,dmax):
    if Bravais in [0,1,2]:          #cubic
        a = b = c = ranaxis(dmin,dmax)
        alp = bet = gam = 90
    elif Bravais in [3,4]:          #hexagonal/trigonal
        a = b = ranaxis(dmin,dmax)
        c = ranaxis(dmin,dmax)
        alp = bet =  90
        gam = 120
    elif Bravais in [5,6]:          #tetragonal
        a = b = ranaxis(dmin,dmax)
        c = ranaxis(dmin,dmax)
        alp = bet = gam = 90
    elif Bravais in [7,8,9,10]:       #orthorhombic - F,I,C,P - a<b<c convention
        abc = [ranaxis(dmin,dmax),ranaxis(dmin,dmax),ranaxis(dmin,dmax)]
        abc.sort()
        a = abc[0]
        b = abc[1]
        c = abc[2]
        alp = bet = gam = 90
    elif Bravais in [11,12]:        #monoclinic - C,P - a<c convention
        ac = [ranaxis(dmin,dmax),ranaxis(dmin,dmax)]
        ac.sort()
        a = ac[0]
        b = ranaxis(dmin,dmax)
        c = ac[1]
        alp = gam = 90
        bet = ranaxis(90.,130.)
    else:                           #triclinic - a<b<c convention
        abc = [ranaxis(dmin,dmax),ranaxis(dmin,dmax),ranaxis(dmin,dmax)]
        abc.sort()
        a = abc[0]
        b = abc[1]
        c = abc[2]
        r = 0.5*b/c
        alp = ranaxis(acosd(r),acosd(-r))
        r = 0.5*a/c
        bet = ranaxis(acosd(r),acosd(-r))
        r = 0.5*a/b
        gam = ranaxis(acosd(r),acosd(-r))  
    return [a,b,c,alp,bet,gam]
    
def calc_M20(peaks,HKL):
    diff = 0
    X20 = 0
    for Nobs20,peak in enumerate(peaks):
        if peak[3]:
            Qobs = 1.0/peak[7]**2
            Qcalc = 1.0/peak[8]**2
            diff += abs(Qobs-Qcalc)
        elif peak[2]:
            X20 += 1
        if Nobs20 == 19: 
            d20 = peak[7]
            break
    else:
        d20 = peak[7]
        Nobs20 = len(peaks)
    for N20,hkl in enumerate(HKL):
        if hkl[3] < d20:
            break                
    eta = diff/Nobs20
    Q20 = 1.0/d20**2
    if diff:
        M20 = Q20/(2.0*diff)
    else:
        M20 = 0
    M20 /= (1.+X20)
    return M20,X20
    
def sortM20(cells):
    #cells is M20,X20,Bravais,a,b,c,alp,bet,gam
    #sort highest M20 1st
    T = []
    for i,M in enumerate(cells):
        T.append((M[0],i))
    D = dict(zip(T,cells))
    T.sort()
    T.reverse()
    X = []
    for key in T:
        X.append(D[key])
    return X
                
def sortVolume(cells):
    #cells is M20,X20,Bravais,a,b,c,alp,bet,gam,volume
    #sort smallest volume 1st
    T = []
    for i,M in enumerate(cells):
        T.append((M[9],i))
    D = dict(zip(T,cells))
    T.sort()
    X = []
    for key in T:
        X.append(D[key])
    return X
                
def IndexPeaks(peaks,HKL):
    import bisect
    N = len(HKL)
    if N == 0: return False
    hklds = list(np.array(HKL).T[3])+[1000.0,0.0,]
    hklds.sort()                                        # ascending sort - upper bound at end
    hklmax = [0,0,0]
    for ipk,peak in enumerate(peaks):
        if peak[2]:
            i = bisect.bisect_right(hklds,peak[7])          # find peak position in hkl list
            dm = peak[7]-hklds[i-1]                         # peak to neighbor hkls in list
            dp = hklds[i]-peak[7]
            pos = N-i                                       # reverse the order
            if dp > dm: pos += 1                            # closer to upper than lower
            if pos >= N:
                print pos,N
                break
            hkl = HKL[pos]                                 # put in hkl
            if hkl[4] >= 0:                                 # peak already assigned - test if this one better
                opeak = peaks[hkl[4]]
                dold = abs(opeak[7]-hkl[3])
                dnew = min(dm,dp)
                if dold > dnew:                             # new better - zero out old
                    opeak[4:7] = [0,0,0]
                    opeak[8] = 0.
                else:                                       # old better - do nothing
                    continue                
            hkl[4] = ipk
            peak[4:7] = hkl[:3]
            peak[8] = hkl[3]                                # fill in d-calc
    for peak in peaks:
        peak[3] = False
        if peak[2]:
            if peak[8] > 0.:
                for j in range(3):
                    if abs(peak[j+4]) > hklmax[j]: hklmax[j] = abs(peak[j+4])
                peak[3] = True
    if hklmax[0]*hklmax[1]*hklmax[2] > 0:
        return True
    else:
        return False
        
def Values2A(ibrav,values):
    if ibrav in [0,1,2]:
        return [values[0],values[0],values[0],0,0,0]
    elif ibrav in [3,4]:
        return [values[0],values[0],values[1],values[0],0,0]
    elif ibrav in [5,6]:
        return [values[0],values[0],values[1],0,0,0]
    elif ibrav in [7,8,9,10]:
        return [values[0],values[1],values[2],0,0,0]
    elif ibrav in [11,12]:
        return [values[0],values[1],values[2],0,values[3],0]
    else:
        return list(values)
        
def A2values(ibrav,A):
    if ibrav in [0,1,2]:
        return [A[0],]
    elif ibrav in [3,4,5,6]:
        return [A[0],A[2]]
    elif ibrav in [7,8,9,10]:
        return [A[0],A[1],A[2]]
    elif ibrav in [11,12]:
        return [A[0],A[1],A[2],A[4]]
    else:
        return A

def FitHKL(ibrav,peaks,A,Pwr):
                
    def errFit(values,ibrav,d,H,Pwr):
        A = Values2A(ibrav,values)
        Qo = 1./d**2
        Qc = G2lat.calc_rDsq(H,A)
        return (Qo-Qc)*d**Pwr
    
    Peaks = np.array(peaks).T
    
    values = A2values(ibrav,A)
    result = so.leastsq(errFit,values,full_output=True,ftol=0.0001,
        args=(ibrav,Peaks[7],Peaks[4:7],Pwr))
    A = Values2A(ibrav,result[0])
    return True,np.sum(errFit(result[0],ibrav,Peaks[7],Peaks[4:7],Pwr)**2),A,result
           
def FitHKLZ(wave,ibrav,peaks,A,Z,Zref,Pwr):
    
    def errFit(values,ibrav,d,H,tth,wave,Zref,Pwr):
        Zero = Z
        if Zref:    
            Zero = values[-1]
        A = Values2A(ibrav,values[:6])
        Qo = 1./d**2
        Qc = G2lat.calc_rDsqZ(H,A,Zero,tth,wave)
        return (Qo-Qc)*d**Pwr
    
    Peaks = np.array(peaks).T
    
    values = A2values(ibrav,A)
    if Zref:
        values.append(Z)
    result = so.leastsq(errFit,values,full_output=True,ftol=0.0001,factor=0.001,
        args=(ibrav,Peaks[7],Peaks[4:7],Peaks[0],wave,Zref,Pwr))
    A = Values2A(ibrav,result[0][:6])
    if Zref:
        Z = result[0][-1]
    
    return True,np.sum(errFit(result[0],ibrav,Peaks[7],Peaks[4:7],Peaks[0],wave,Zref,Pwr)**2),A,Z,result
               
def rotOrthoA(A):
    return [A[1],A[2],A[0],0,0,0]
    
def swapMonoA(A):
    return [A[2],A[1],A[0],0,A[4],0]
    
def oddPeak(indx,peaks):
    noOdd = True
    for peak in peaks:
        H = peak[4:7]
        if H[indx] % 2:
            noOdd = False
    return noOdd
    
def halfCell(ibrav,A,peaks):
    if ibrav in [0,1,2]:
        if oddPeak(0,peaks):
            A[0] *= 2
            A[1] = A[2] = A[0]
    elif ibrav in [3,4,5,6]:
        if oddPeak(0,peaks):
            A[0] *= 2
            A[1] = A[0]
        if oddPeak(2,peaks):
            A[2] *=2
    else:
        if oddPeak(0,peaks):
            A[0] *=2
        if oddPeak(1,peaks):
            A[1] *=2
        if oddPeak(2,peaks):
            A[2] *=2
    return A
    
def getDmin(peaks):
    return peaks[-1][7]
    
def getDmax(peaks):
    return peaks[0][7]
    
def refinePeaksZ(peaks,wave,ibrav,A,Zero,ZeroRef):
    dmin = getDmin(peaks)
    OK,smin,Aref,Z,result = FitHKLZ(wave,ibrav,peaks,A,Zero,ZeroRef,0)
    Peaks = np.array(peaks).T
    H = Peaks[4:7]
    Peaks[8] = 1./np.sqrt(G2lat.calc_rDsqZ(H,Aref,Z,Peaks[0],wave))
    peaks = Peaks.T    
    HKL = G2lat.GenHBravais(dmin,ibrav,Aref)
    M20,X20 = calc_M20(peaks,HKL)
    return len(HKL),M20,X20,Aref,Z
    
def refinePeaks(peaks,ibrav,A):
    dmin = getDmin(peaks)
    smin = 1.0e10
    pwr = 8
    maxTries = 10
    OK = False
    tries = 0
    HKL = G2lat.GenHBravais(dmin,ibrav,A)
    while len(HKL) > 2 and IndexPeaks(peaks,HKL):
        Pwr = pwr - (tries % 2)
        HKL = []
        tries += 1
        osmin = smin
        oldA = A[:]
        Vold = G2lat.calc_V(oldA)
        OK,smin,A,result = FitHKL(ibrav,peaks,A,Pwr)
        Vnew = G2lat.calc_V(A)
        if Vnew > 2.0*Vold or Vnew < 2.:
            A = ranAbyR(ibrav,oldA,tries+1,maxTries,ran2axis)
            OK = False
            continue
        try:
            HKL = G2lat.GenHBravais(dmin,ibrav,A)
        except FloatingPointError:
            A = oldA
            OK = False
            break
        if len(HKL) == 0: break                         #absurd cell obtained!
        rat = (osmin-smin)/smin
        if abs(rat) < 1.0e-5 or not OK: break
        if tries > maxTries: break
    if OK:
        OK,smin,A,result = FitHKL(ibrav,peaks,A,2)
        Peaks = np.array(peaks).T
        H = Peaks[4:7]
        try:
            Peaks[8] = 1./np.sqrt(G2lat.calc_rDsq(H,A))
            peaks = Peaks.T
        except FloatingPointError:
            A = oldA
        
    M20,X20 = calc_M20(peaks,HKL)
    return len(HKL),M20,X20,A
        
def findBestCell(dlg,ncMax,A,Ntries,ibrav,peaks,V1):
# dlg & ncMax are used for wx progress bar 
# A != 0 find the best A near input A,
# A = 0 for random cell, volume normalized to V1;
# returns number of generated hkls, M20, X20 & A for best found
    mHKL = [3,3,3, 5,5, 5,5, 7,7,7,7, 9,9, 10]
    dmin = getDmin(peaks)-0.05
    amin = 2.5
    amax = 5.*getDmax(peaks)
    Asave = []
    GoOn = True
    if A:
        HKL = G2lat.GenHBravais(dmin,ibrav,A[:])
        if len(HKL) > mHKL[ibrav]:
            IndexPeaks(peaks,HKL)
            Asave.append([calc_M20(peaks,HKL),A[:]])
    tries = 0
    while tries < Ntries:
        if A:
            Abeg = ranAbyR(ibrav,A,tries+1,Ntries,ran2axis)
            if ibrav in [11,12,13]:         #monoclinic & triclinic
                Abeg = ranAbyR(ibrav,A,tries/10+1,Ntries,ran2axis)
        else:
            Abeg = ranAbyV(ibrav,amin,amax,V1)
        HKL = G2lat.GenHBravais(dmin,ibrav,Abeg)
        
        if IndexPeaks(peaks,HKL) and len(HKL) > mHKL[ibrav]:
            Lhkl,M20,X20,Aref = refinePeaks(peaks,ibrav,Abeg)
            Asave.append([calc_M20(peaks,HKL),Aref[:]])
            if ibrav == 9:                          #C-centered orthorhombic
                for i in range(2):
                    Abeg = rotOrthoA(Abeg[:])
                    Lhkl,M20,X20,Aref = refinePeaks(peaks,ibrav,Abeg)
                    HKL = G2lat.GenHBravais(dmin,ibrav,Aref)
                    IndexPeaks(peaks,HKL)
                    Asave.append([calc_M20(peaks,HKL),Aref[:]])
            elif ibrav == 11:                      #C-centered monoclinic
                Abeg = swapMonoA(Abeg[:])
                Lhkl,M20,X20,Aref = refinePeaks(peaks,ibrav,Abeg)
                HKL = G2lat.GenHBravais(dmin,ibrav,Aref)
                IndexPeaks(peaks,HKL)
                Asave.append([calc_M20(peaks,HKL),Aref[:]])
        else:
            break
        Nc = len(HKL)
        if Nc >= ncMax:
            GoOn = False
        elif dlg:
            GoOn = dlg.Update(Nc)[0]
            if not GoOn:
                break
        tries += 1    
    X = sortM20(Asave)
    if X:
        Lhkl,M20,X20,A = refinePeaks(peaks,ibrav,X[0][1])
        return GoOn,Lhkl,M20,X20,A
        
    else:
        return GoOn,0,0,0,0
        
def monoCellReduce(ibrav,A):
    a,b,c,alp,bet,gam = G2lat.A2cell(A)
    G,g = G2lat.A2Gmat(A)
    if ibrav in [11]:
        u = [0,0,-1]
        v = [1,0,2]
        anew = math.sqrt(np.dot(np.dot(v,g),v))
        if anew < a:
            cang = np.dot(np.dot(u,g),v)/(anew*c)
            beta = acosd(-abs(cang))
            A = G2lat.cell2A([anew,b,c,90,beta,90])
    else:
        u = [-1,0,0]
        v = [1,0,1]
        cnew = math.sqrt(np.dot(np.dot(v,g),v))
        if cnew < c:
            cang = np.dot(np.dot(u,g),v)/(a*cnew)
            beta = acosd(-abs(cang))
            A = G2lat.cell2A([a,b,cnew,90,beta,90])
    return A

def DoIndexPeaks(peaks,inst,controls,bravais):
    
    delt = 0.005                                     #lowest d-spacing cushion - can be fixed?
    amin = 2.5
    amax = 5.0*getDmax(peaks)
    dmin = getDmin(peaks)-delt
    bravaisNames = ['Cubic-F','Cubic-I','Cubic-P','Trigonal-R','Trigonal/Hexagonal-P',
        'Tetragonal-I','Tetragonal-P','Orthorhombic-F','Orthorhombic-I','Orthorhombic-C',
        'Orthorhombic-P','Monoclinic-C','Monoclinic-P','Triclinic']
    tries = ['1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th']
    N1s = [1,1,1,   5,5,  5,5, 50,50,50,50,  50,50, 200]
    N2s = [1,1,1,   2,2,  2,2,     2,2,2,2,   2,2,   4]
    Nm  = [1,1,1,   1,1,  1,1,     1,1,1,1,   2,2,   4]
    Nobs = len(peaks)
    wave = inst[1]
    zero,ncno = controls[1:3]
    ncMax = Nobs*ncno
    print "%s %8.3f %8.3f" % ('lattice parameter range = ',amin,amax)
    print "%s %8.5f %s %.4f %s %d %s %d" % ('Wavelength =',wave,'Zero =',zero,'Nc/No max =',ncno,' Max Nc =',ncno*Nobs)
    cells = []
    for ibrav in range(14):
        begin = time.time()
        if bravais[ibrav]:
            print 'cell search for ',bravaisNames[ibrav]
            print '      M20  X20  Nc       a          b          c        alpha       beta      gamma     volume      V-test'
            V1 = controls[3]
            bestM20 = 0
            topM20 = 0
            cycle = 0
            while cycle < 5:
                dlg = wx.ProgressDialog("Generated reflections",tries[cycle]+" cell search for "+bravaisNames[ibrav],ncMax, 
                    style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
                screenSize = wx.ClientDisplayRect()
                Size = dlg.GetSize()
                dlg.SetPosition(wx.Point(screenSize[2]-Size[0]-305,screenSize[1]+5))
                try:
                    GoOn = True
                    while GoOn:                                                 #Loop over increment of volume
                        N2 = 0
                        while N2 < N2s[ibrav]:                                  #Table 2 step (iii)               
                            if ibrav > 2:
                                if not N2:
                                    A = []
                                    GoOn,Nc,M20,X20,A = findBestCell(dlg,ncMax,A,Nm[ibrav]*N1s[ibrav],ibrav,peaks,V1)
                                if A:
                                    GoOn,Nc,M20,X20,A = findBestCell(dlg,ncMax,A[:],N1s[ibrav],ibrav,peaks,0)
                            else:
                                GoOn,Nc,M20,X20,A = findBestCell(dlg,ncMax,0,Nm[ibrav]*N1s[ibrav],ibrav,peaks,V1)
                            if Nc >= ncMax:
                                GoOn = False
                                break
                            elif 3*Nc < Nobs:
                                N2 = 10
                                break
                            else:
                                if not GoOn:
                                    break
                                if M20 > 1.0:
                                    bestM20 = max(bestM20,M20)
                                    A = halfCell(ibrav,A[:],peaks)
                                    if ibrav in [12]:
                                        A = monoCellReduce(ibrav,A[:])
                                    HKL = G2lat.GenHBravais(dmin,ibrav,A)
                                    IndexPeaks(peaks,HKL)
                                    a,b,c,alp,bet,gam = G2lat.A2cell(A)
                                    V = G2lat.calc_V(A)
                                    print "%10.3f %3d %3d %10.5f %10.5f %10.5f %10.3f %10.3f %10.3f %10.2f %10.2f" % (M20,X20,Nc,a,b,c,alp,bet,gam,V,V1)
                                    if M20 >= 2.0:
                                        cells.append([M20,X20,ibrav,a,b,c,alp,bet,gam,V,False,False])
                            if not GoOn:
                                break
                            N2 += 1
                        if ibrav < 11:
                            V1 *= 1.1
                        elif ibrav in range(11,14):
                            V1 *= 1.05
                        if not GoOn:
                            if bestM20 > topM20:
                                topM20 = bestM20
                                if cells:
                                    V1 = cells[0][9]
                                else:
                                    V1 = 25
                                ncMax += Nobs
                                cycle += 1
                                print 'Restart search, new Max Nc = ',ncMax
                            else:
                                cycle = 10
                finally:
                    dlg.Destroy()
            print '%s%s%s%s'%('finished cell search for ',bravaisNames[ibrav], \
                ', elapsed time = ',G2lat.sec2HMS(time.time()-begin))
            
    if cells:
        cells = sortM20(cells)
        cells[0][-1] = True
        return True,dmin,cells
    else:
        return False,0,[]
        
        
NeedTestData = True
def TestData():
    array = np.array
    global NeedTestData
    NeedTestData = False
    global TestData
    TestData = [12, [7.,8.70,10.86,90.,102.95,90.], [7.76006,8.706215,10.865679,90.,102.947,90.],3,
        [[2.176562137832974, 761.60902227696033, True, True, 0, 0, 1, 10.591300714328161, 10.589436], 
        [3.0477561489789498, 4087.2956049071572, True, True, 1, 0, 0, 7.564238997554908, 7.562777], 
        [3.3254921120068524, 1707.0253890991009, True, True, 1, 0, -1, 6.932650301411212, 6.932718], 
        [3.428121546163426, 2777.5082170150563, True, True, 0, 1, 1, 6.725163158013632, 6.725106], 
        [4.0379791325512118, 1598.4321673135987, True, True, 1, 1, 0, 5.709789097440156, 5.70946], 
        [4.2511182350743937, 473.10955149057577, True, True, 1, 1, -1, 5.423637972781876, 5.42333], 
        [4.354684330373451, 569.88528280256071, True, True, 0, 0, 2, 5.2947091882172534, 5.294718],
        [4.723324574319177, 342.73882372499997, True, True, 1, 0, -2, 4.881681587039431, 4.881592], 
        [4.9014773581253994, 5886.3516356615492, True, True, 1, 1, 1, 4.704350709093183, 4.70413], 
        [5.0970774474587275, 3459.7541692903033, True, True, 0, 1, 2, 4.523933797797693, 4.523829], 
        [5.2971997607389518, 1290.0229964239879, True, True, 0, 2, 0, 4.353139557169142, 4.353108], 
        [5.4161306205553847, 1472.5726977257755, True, True, 1, 1, -2, 4.257619398422479, 4.257944], 
        [5.7277364698554161, 1577.8791668322888, True, True, 0, 2, 1, 4.026169751907777, 4.026193], 
        [5.8500213058834163, 420.74210142657131, True, True, 1, 0, 2, 3.9420803081518443, 3.942219],
        [6.0986764166731708, 163.02160537058708, True, True, 2, 0, 0, 3.7814965150452537, 3.781389], 
        [6.1126665157702753, 943.25461245706833, True, True, 1, 2, 0, 3.772849962062199, 3.772764], 
        [6.2559260555056957, 250.55355015505376, True, True, 1, 2, -1, 3.6865353266375283, 3.686602], 
        [6.4226243128279892, 5388.5560141098349, True, True, 1, 1, 2, 3.5909481979190283, 3.591214], 
        [6.5346132446561134, 1951.6070344509026, True, True, 0, 0, 3, 3.5294722429440584, 3.529812], 
        [6.5586952135236443, 259.65938178131034, True, True, 2, 1, -1, 3.516526936765838, 3.516784], 
        [6.6509216222783722, 93.265376597376573, True, True, 2, 1, 0, 3.4678179073694952, 3.468369], 
        [6.7152737044107722, 289.39386813803162, True, True, 1, 2, 1, 3.4346235125812807, 3.434648], 
        [6.8594130457361899, 603.54959764648322, True, True, 0, 2, 2, 3.362534044860622, 3.362553], 
        [7.0511627728884454, 126.43246447656593, True, True, 0, 1, 3, 3.2712038721790675, 3.271181], 
        [7.077700845503319, 125.49742760019636, True, True, 1, 1, -3, 3.2589538988480626, 3.259037], 
        [7.099393757363675, 416.55444885434633, True, True, 1, 2, -2, 3.2490085228959193, 3.248951], 
        [7.1623933278642742, 369.27397110921817, True, True, 2, 1, -2, 3.2204673608202383, 3.220487], 
        [7.4121734953058924, 482.84120827021826, True, True, 2, 1, 1, 3.1120858221599876, 3.112308]]
        ]
    global TestData2
    TestData2 = [12, [0.15336547830008007, 0.017345499139401827, 0.008122368657493792, 0, 0.02893538955687591, 0], 3,
        [[2.176562137832974, 761.6090222769603, True, True, 0, 0, 1, 10.591300714328161, 11.095801], 
        [3.0477561489789498, 4087.295604907157, True, True, 0, 1, 0, 7.564238997554908, 7.592881], 
        [3.3254921120068524, 1707.025389099101, True, False, 0, 0, 0, 6.932650301411212, 0.0], 
        [3.428121546163426, 2777.5082170150563, True, True, 0, 1, 1, 6.725163158013632, 6.266192], 
        [4.037979132551212, 1598.4321673135987, True, False, 0, 0, 0, 5.709789097440156, 0.0], 
        [4.251118235074394, 473.10955149057577, True, True, 0, 0, 2, 5.423637972781876, 5.5479], 
        [4.354684330373451, 569.8852828025607, True, True, 0, 0, 2, 5.2947091882172534, 5.199754], 
        [4.723324574319177, 342.738823725, True, False, 0, 0, 0, 4.881681587039431, 0.0], 
        [4.901477358125399, 5886.351635661549, True, False, 0, 0, 0, 4.704350709093183, 0.0], 
        [5.0970774474587275, 3459.7541692903033, True, True, 0, 1, 2, 4.523933797797693, 4.479534], 
        [5.297199760738952, 1290.022996423988, True, True, 0, 1, 0, 4.353139557169142, 4.345087],
        [5.416130620555385, 1472.5726977257755, True, False, 0, 0, 0, 4.257619398422479, 0.0], 
        [5.727736469855416, 1577.8791668322888, True, False, 0, 0, 0, 4.026169751907777, 0.0], 
        [5.850021305883416, 420.7421014265713, True, False, 0, 0, 0, 3.9420803081518443, 0.0], 
        [6.098676416673171, 163.02160537058708, True, True, 0, 2, 0, 3.7814965150452537, 3.796441], 
        [6.112666515770275, 943.2546124570683, True, False, 0, 0, 0, 3.772849962062199, 0.0], 
        [6.255926055505696, 250.55355015505376, True, True, 0, 0, 3, 3.6865353266375283, 3.6986], 
        [6.422624312827989, 5388.556014109835, True, True, 0, 2, 1, 3.5909481979190283, 3.592005], 
        [6.534613244656113, 191.6070344509026, True, True, 1, 0, -1, 3.5294722429440584, 3.546166], 
        [6.558695213523644, 259.65938178131034, True, True, 0, 0, 3, 3.516526936765838, 3.495428], 
        [6.650921622278372, 93.26537659737657, True, True, 0, 0, 3, 3.4678179073694952, 3.466503], 
        [6.715273704410772, 289.3938681380316, True, False, 0, 0, 0, 3.4346235125812807, 0.0], 
        [6.85941304573619, 603.5495976464832, True, True, 0, 1, 3, 3.362534044860622, 3.32509], 
        [7.051162772888445, 126.43246447656593, True, True, 0, 1, 2, 3.2712038721790675, 3.352121], 
        [7.077700845503319, 125.49742760019636, True, False, 0, 0, 0, 3.2589538988480626, 0.0], 
        [7.099393757363675, 416.55444885434633, True, False, 0, 0, 0, 3.2490085228959193, 0.0], 
        [7.162393327864274, 369.27397110921817, True, False, 0, 0, 0, 3.2204673608202383, 0.0], 
        [7.412173495305892, 482.84120827021826, True, True, 0, 2, 2, 3.112085822159976, 3.133096]]
        ]
#
def test0():
    if NeedTestData: TestData()
    msg = 'test FitHKL'
    ibrav,cell,bestcell,Pwr,peaks = TestData
    print 'best cell:',bestcell
    print 'old cell:',cell
    Peaks = np.array(peaks)
    HKL = Peaks[4:7]
    print calc_M20(peaks,HKL)
    A = G2lat.cell2A(cell)
    OK,smin,A,result = FitHKL(ibrav,peaks,A,Pwr)
    print 'new cell:',G2lat.A2cell(A)    
    print 'x:',result[0]
    print 'cov_x:',result[1]
    print 'infodict:'
    for item in result[2]:
        print item,result[2][item]
    print 'msg:',result[3]
    print 'ier:',result[4]
    result = refinePeaks(peaks,ibrav,A)
    N,M20,X20,A = result
    print 'refinePeaks:',N,M20,X20,G2lat.A2cell(A)
    print 'compare bestcell:',bestcell
#
def test1():
    if NeedTestData: TestData()
    msg = 'test FitHKL'
    ibrav,A,Pwr,peaks = TestData2
    print 'bad cell:',G2lat.A2cell(A)
    print 'FitHKL'
    OK,smin,A,result = FitHKL(ibrav,peaks,A,Pwr)
    result = refinePeaks(peaks,ibrav,A)
    N,M20,X20,A = result
    print 'refinePeaks:',N,M20,X20,A
#    Peaks = np.array(peaks)
#    HKL = Peaks[4:7]
#    print calc_M20(peaks,HKL)
#    OK,smin,A,result = FitHKL(ibrav,peaks,A,Pwr)
#    print 'new cell:',G2lat.A2cell(A)    
#    print 'x:',result[0]
#    print 'cov_x:',result[1]
#    print 'infodict:'
#    for item in result[2]:
#        print item,result[2][item]
#    print 'msg:',result[3]
#    print 'ier:',result[4]
    
#
if __name__ == '__main__':
    test0()
    test1()
#    test2()
#    test3()
#    test4()
#    test5()
#    test6()
#    test7()
#    test8()
    print "OK"
