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
import pypowder as pyp              #assumes path has been amended to include correctr bin directory
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

def FitHKL(ibrav,peaks,A,Pwr):
    
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
            return values
            
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
    
    def errFit(values,ibrav,d,H,Pwr):
        A = Values2A(ibrav,values)
        Qo = 1./d**2
        Qc = G2lat.calc_rDsq(H,A)
        return (Qo-Qc)*d**Pwr
    
    Peaks = np.array(peaks).T
    values = A2values(ibrav,A)    
    result = so.leastsq(errFit,values,args=(ibrav,Peaks[7],Peaks[4:7],Pwr),full_output=True)
    A = Values2A(ibrav,result[0])
    return True,np.sum(errFit(result[0],ibrav,Peaks[7],Peaks[4:7],Pwr)**2),A
           
def FitHKLZ(ibrav,peaks,Z,A):
    return A,Z
    
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
    
def refinePeaks(peaks,ibrav,A):
    dmin = getDmin(peaks)
    smin = 1.0e10
    pwr = 3
    maxTries = 10
    if ibrav == 13:
        pwr = 4
        maxTries = 10
    OK = False
    tries = 0
    HKL = G2lat.GenHBravais(dmin,ibrav,A)
    while IndexPeaks(peaks,HKL):
        Pwr = pwr - (tries % 2)
        HKL = []
        tries += 1
        osmin = smin
        oldA = A
        OK,smin,A = FitHKL(ibrav,peaks,A,Pwr)
        if min(A[:3]) <= 0:
            A = oldA
            OK = False
            break
        if OK:
            HKL = G2lat.GenHBravais(dmin,ibrav,A)
        if len(HKL) == 0: break                         #absurd cell obtained!
        rat = (osmin-smin)/smin
        if abs(rat) < 1.0e-5 or not OK: break
        if tries > maxTries: break
    if OK:
        OK,smin,A = FitHKL(ibrav,peaks,A,2)
        Peaks = np.array(peaks).T
        H = Peaks[4:7]
        Peaks[8] = 1./np.sqrt(G2lat.calc_rDsq(H,A))
        peaks = Peaks.T
        
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
            Anew = ranAbyR(ibrav,A[:],tries+1,Ntries,ran2axis)
            if ibrav in [11,12,13]:
                Anew = ranAbyR(ibrav,A[:],tries/10+1,Ntries,ran2axis)
        else:
            Anew = ranAbyV(ibrav,amin,amax,V1)
        HKL = G2lat.GenHBravais(dmin,ibrav,Anew)
        
        if IndexPeaks(peaks,HKL) and len(HKL) > mHKL[ibrav]:
            Lhkl,M20,X20,Anew = refinePeaks(peaks,ibrav,Anew)
            Asave.append([calc_M20(peaks,HKL),Anew[:]])
            if ibrav == 9:                          #C-centered orthorhombic
                for i in range(2):
                    Anew = rotOrthoA(Anew[:])
                    Lhkl,M20,X20,Anew = refinePeaks(peaks,ibrav,Anew)
                    HKL = G2lat.GenHBravais(dmin,ibrav,Anew)
                    IndexPeaks(peaks,HKL)
                    Asave.append([calc_M20(peaks,HKL),Anew[:]])
            elif ibrav == 11:                      #C-centered monoclinic
                Anew = swapMonoA(Anew[:])
                Lhkl,M20,X20,Anew = refinePeaks(peaks,ibrav,Anew)
                HKL = G2lat.GenHBravais(dmin,ibrav,Anew)
                IndexPeaks(peaks,HKL)
                Asave.append([calc_M20(peaks,HKL),Anew[:]])
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
        return GoOn,Lhkl,M20,X20,X[0][1]
    else:
        return GoOn,0,0,0,Anew
        
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
    if len(inst) > 10:
        zero = inst[3]
    else:
        zero = inst[2]
    print "%s %8.5f %6.3f" % ('wavelength, zero =',wave,zero)
    print "%s %8.3f %8.3f" % ('lattice parameter range = ',amin,amax)
    ifzero,maxzero,ncno = controls[:3]
    ncMax = Nobs*ncno
    print "%s %d %s %d %s %d" % ('change zero =',ifzero,'Nc/No max =',ncno,' Max Nc =',ncno*Nobs)
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
                                        cells.append([M20,X20,ibrav,a,b,c,alp,bet,gam,V,False])
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
        return False,0,0
        