#GSASII computational module
import sys
import math
import wx
import time
import numpy as np
import numpy.linalg as nl
import os.path as ospath
import GSASIIpath
import pypowder as pyp              #assumes path has been amended to include correctr bin directory
import GSASIIplot as G2plt

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

def sec2HMS(sec):
    H = int(sec/3600)
    M = int(sec/60-H*60)
    S = sec-3600*H-60*M
    return '%d:%2d:%.2f'%(H,M,S)

def ValEsd(value,esd=0,nTZ=False):
    # returns value(esd) string; nTZ=True for no trailing zeros
    # use esd < 0 for level of precision shown e.g. esd=-0.01 gives 2 places beyond decimal
    #get the 2 significant digits in the esd 
    edig = lambda esd: int(round(10**(math.log10(esd) % 1+1)))
    #get the number of digits to represent them 
    epl = lambda esd: 2+int(1.545-math.log10(10*edig(esd)))
    
    mdec = lambda esd: -int(math.log10(abs(esd)))
    ndec = lambda esd: int(1.545-math.log10(abs(esd)))
    if esd > 0:
        fmt = '"%.'+str(ndec(esd))+'f(%d)"'
        print fmt,ndec(esd),esd*10**(mdec(esd)+1)
        return fmt%(value,int(esd*10**(mdec(esd)+1)))
    elif esd < 0:
         return str(round(value,mdec(esd)))
    else:
        text = "%F"%(value)
        if nTZ:
            return text.rstrip('0')
        else:
            return text

        
        

def DoPeakFit(peaks,background,limits,inst,data):
    
    def backgroundPrint(background,sigback):
        if background[1]:
            print 'Background coefficients for',background[0],'function'
            ptfmt = "%12.5f"
            ptstr =  'values:'
            sigstr = 'esds  :'
            for i,back in enumerate(background[3:]):
                ptstr += ptfmt % (back)
                sigstr += ptfmt % (sigback[i+3])
            print ptstr
            print sigstr
        else:
            print 'Background not refined'
    
    def instPrint(instVal,siginst,insLabels):
        print 'Instrument Parameters:'
        ptfmt = "%12.6f"
        ptlbls = 'names :'
        ptstr =  'values:'
        sigstr = 'esds  :'
        for i,value in enumerate(instVal):
            ptlbls += "%s" % (insLabels[i].center(12))
            ptstr += ptfmt % (value)
            if siginst[i]:
                sigstr += ptfmt % (siginst[i])
            else:
                sigstr += 12*' '
        print ptlbls
        print ptstr
        print sigstr
    
    def peaksPrint(peaks,sigpeaks):
        print 'Peak list='

    begin = time.time()
    Np = len(peaks[0])
    DataType = inst[1][0]
    instVal = inst[1][1:]
    Insref = inst[2][1:]
    insLabels = inst[3][1:]
    Ka2 = False
    Ioff = 3
    if len(instVal) == 12:
        lamratio = instVal[1]/instVal[0]
        Ka2 = True
        Ioff = 5
    insref = Insref[len(Insref)-7:-1]               #just U,V,W,X,Y,SH/L
    for peak in peaks:
        dip = []
        dip.append(tand(peak[0]/2.0)**2)
        dip.append(tand(peak[0]/2.0))
        dip.append(1.0/cosd(peak[0]/2.0))
        dip.append(tand(peak[0]/2.0))
        peak.append(dip)
    B = background[2]
    bcof = background[3:3+B]
    Bv = 0
    if background[1]:
        Bv = B
    x,y,w,yc,yb,yd = data               #these are numpy arrays!
    V = []
    A = []
    swobs = 0.0
    smin = 0.0
    first = True
    GoOn = True
    Go = True
    dlg = wx.ProgressDialog("Elapsed time","Fitting peaks to pattern",len(x), \
        style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE|wx.PD_REMAINING_TIME|wx.PD_CAN_ABORT)
    screenSize = wx.DisplaySize()
    Size = dlg.GetSize()
    dlg.SetPosition(wx.Point(screenSize[0]-Size[0]-300,0))
    try:
        i = 0
        for xi in x:
            Go = dlg.Update(i)[0]
            if GoOn:
                GoOn = Go
            if limits[0] <= xi <= limits[1]:
                yb[i] = 0.0
                dp = []
                for j in range(B):
                    t = (xi-limits[0])**j
                    yb[i] += t*bcof[j]
                    if background[1]:
                        dp.append(t)
                yc[i] = yb[i]
                Iv = 0
                for j in range(6):
                    if insref[j]:
                        dp.append(0.0)
                        Iv += 1
                for peak in peaks:
                    dip = peak[-1]
                    f = pyp.pypsvfcj(peak[2],xi-peak[0],peak[0],peak[4],peak[6],instVal[-2],0.0)
                    yc[i] += f[0]*peak[2]
                    if f[0] > 0.0:
                        j = 0
                        if insref[0]:              #U
                            dp[Bv+j] += f[3]*dip[0]
                            j += 1
                        if insref[1]:              #V
                            dp[Bv+j] += f[3]*dip[1]
                            j += 1
                        if insref[2]:              #W
                            dp[Bv+j] += f[3]
                            j += 1
                        if insref[3]:              #X
                            dp[Bv+j] += f[4]*dip[2]
                            j += 1
                        if insref[4]:              #Y
                            dp[Bv+j] += f[4]*dip[3]
                            j += 1
                        if insref[5]:              #SH/L
                            dp[Bv+j] += f[5]
                    if Ka2:
                       pos2 = 2.0*asind(lamratio*sind(peak[0]/2.0))
                       f2 = pyp.pypsvfcj(peak[2],xi-pos2,peak[0],peak[4],peak[6],instVal[-2],0.0)
                       yc[i] += f2[0]*peak[2]*instVal[3]
                       if f[0] > 0.0:
                           j = 0
                           if insref[0]:              #U
                               dp[Bv+j] += f2[3]*dip[0]*instVal[3]
                               j += 1
                           if insref[1]:              #V
                               dp[Bv+j] += f2[3]*dip[1]*instVal[3]
                               j += 1
                           if insref[2]:              #W
                               dp[Bv+j] += f2[3]*instVal[3]
                               j += 1
                           if insref[3]:              #X
                               dp[Bv+j] += f2[4]*dip[2]*instVal[3]
                               j += 1
                           if insref[4]:              #Y
                               dp[Bv+j] += f2[4]*dip[3]*instVal[3]
                               j += 1
                           if insref[5]:              #SH/L
                               dp[Bv+j] += f2[5]*instVal[3]                       
                    for j in range(0,Np,2):
                        if peak[j+1]: dp.append(f[j/2+1])
                yd[i] = y[i]-yc[i]
                swobs += w[i]*y[i]**2
                t2 = w[i]*yd[i]
                smin += t2*yd[i]                 
                if first:
                    first = False
                    M = len(dp)
                    A = np.zeros(shape=(M,M))
                    V = np.zeros(shape=(M))
                A,V = pyp.buildmv(t2,w[i],M,dp,A,V)
            i += 1
    finally:
        dlg.Destroy()
    print GoOn
    Rwp = smin/swobs
    Rwp = math.sqrt(Rwp)*100.0
    norm = np.diag(A)
    for i,elm in enumerate(norm):
        if elm <= 0.0:
            print norm
            return False,0,0,0,False
    for i in xrange(len(V)):
        norm[i] = 1.0/math.sqrt(norm[i])
        V[i] *= norm[i]
        a = A[i]
        for j in xrange(len(V)):
            a[j] *= norm[i]
    A = np.transpose(A)
    for i in xrange(len(V)):
        a = A[i]
        for j in xrange(len(V)):
            a[j] *= norm[i]
    b = nl.solve(A,V)
    A = nl.inv(A)
    sig = np.diag(A)
    for i in xrange(len(V)):
        b[i] *= norm[i]
        sig[i] *= norm[i]*norm[i]
        sig[i] = math.sqrt(abs(sig[i]))
    sigback = [0,0,0]
    for j in range(Bv):
        background[j+3] += b[j]
        sigback.append(sig[j])
    backgroundPrint(background,sigback)
    k = 0
    delt = []
    if Ka2:
        siginst = [0,0,0,0,0]
    else:
        siginst = [0,0,0]
    for j in range(6):
        if insref[j]:
            instVal[j+Ioff] += b[Bv+k]*0.5
            siginst.append(sig[Bv+k])
            delt.append(b[Bv+k]*0.5)
            k += 1
        else:
            delt.append(0.0)
            siginst.append(0.0)
    delt.append(0.0)                    #dummies for azm
    siginst.append(0.0)
    instPrint(instVal,siginst,insLabels)
    inst[1] = [DataType,]
    for val in instVal:
        inst[1].append(val)
    B = Bv+Iv
    for peak in peaks:
        del peak[-1]                        # remove dip from end
        delsig = delt[0]*tand(peak[0]/2.0)**2+delt[1]*tand(peak[0]/2.0)+delt[2]
        delgam = delt[3]/cosd(peak[0]/2.0)+delt[4]*tand(peak[0]/2.0)
        for j in range(0,len(peak[:-1]),2):
            if peak[j+1]: 
                peak[j] += b[B]*0.5
                B += 1
        peak[4] += delsig
        if peak[4] < 0.0:
            print 'ERROR - negative sigma'
            return False,0,0,0,False            
        peak[6] += delgam
        if peak[6] < 0.0:
            print 'ERROR - negative gamma'
            return False,0,0,0,False
    runtime = time.time()-begin    
    data = [x,y,w,yc,yb,yd]
    return True,smin,Rwp,runtime,GoOn
    
#some cell utilities
#for these: H = [h,k,l]; A is as used in calc_rDsq; G - inv metric tensor, g - metric tensor; 
#           cell - a,b,c,alp,bet,gam in A & deg
   
def calc_rDsq(H,A):
    rdsq = H[0]*H[0]*A[0]+H[1]*H[1]*A[1]+H[2]*H[2]*A[2]+H[0]*H[1]*A[3]+H[0]*H[2]*A[4]+H[1]*H[2]*A[5]
    return rdsq
    
def calc_rDsqZ(H,A,Z,tth,lam):
    rpd = math.pi/180.
    rdsq = calc_rDsq(H,A)+Z*math.sin(tth*rpd)*2.0*rpd/(lam*lam)
    return rdsq
    
def calc_rVsq(A):
    rVsq = A[0]*A[1]*A[2]+0.25*(A[3]*A[4]*A[5]-A[0]*A[5]**2-A[1]*A[4]**2-A[2]*A[3]**2)
    if rVsq < 0:
        return 1
    return rVsq
    
def calc_rV(A):
    return math.sqrt(calc_rVsq(A))
    
def calc_V(A):
    return 1./calc_rV(A)
   
def scaleAbyV(A,V):
    v = calc_V(A)
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
    
def ranNaxis(k,N):
    import random as rand
    T = 1.0+1.0*k/N
    B = 1.0-1.0*k/N
    R = (T-B)*rand.random()+B
    return R
    
def ranAbyV(Bravais,dmin,dmax,V):
    cell = [0,0,0,0,0,0]
    bad = True
    while bad:
        bad = False
        cell = rancell(Bravais,dmin,dmax)
        G,g = cell2Gmat(cell)
        A = Gmat2A(G)
        if calc_rVsq(A) < 1:
            scaleAbyV(A,V)
            cell = A2cell(A)
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
    
def A2Gmat(A):
    G = np.zeros(shape=(3,3))
    G = [[A[0],A[3]/2.,A[4]/2.], [A[3]/2.,A[1],A[5]/2.], [A[4]/2.,A[5]/2.,A[2]]]
    g = nl.inv(G)
    return G,g
    
def fillgmat(cell):
    a,b,c,alp,bet,gam = cell
    g = np.array([[a*a,a*b*cosd(gam),a*c*cosd(bet)],[a*b*cosd(gam),b*b,b*c*cosd(alp)],
        [a*c*cosd(bet),b*c*cosd(alp),c*c]])
    return g
           
def cell2Gmat(cell):
    #returns reciprocal (G) & real (g) metric tensors
    g = fillgmat(cell)
    G = nl.inv(g)        
    return G,g
    
def invcell2Gmat(invcell):
    G = fillgmat(invcell)
    g = nl.inv(G)
    return G,g
    
def Gmat2cell(g):
    #returns lattice parameters from real metric tensor (g)
    a = math.sqrt(max(0,g[0][0]))
    b = math.sqrt(max(0,g[1][1]))
    c = math.sqrt(max(0,g[2][2]))
    alp = acosd(g[2][1]/(b*c))
    bet = acosd(g[2][0]/(a*c))
    gam = acosd(g[0][1]/(a*b))
    return a,b,c,alp,bet,gam
    
def Gmat2A(G):
    return [G[0][0],G[1][1],G[2][2],2.*G[0][1],2.*G[0][2],2.*G[1][2]]
    
def cell2A(cell):
    G,g = cell2Gmat(cell)
    return Gmat2A(G)
    
def A2cell(A):
    G,g = A2Gmat(A)
    return Gmat2cell(g)
    
def A2invcell(A):
    ainv = math.sqrt(max(0.,A[0]))
    binv = math.sqrt(max(0.,A[1]))
    cinv = math.sqrt(max(0.,A[2]))
    gaminv = acosd(max(-0.5,min(0.5,0.5*A[3]/(ainv*binv))))
    betinv = acosd(max(-0.5,min(0.5,0.5*A[4]/(ainv*cinv))))
    alpinv = acosd(max(-0.5,min(0.5,0.5*A[5]/(binv*cinv))))
    return ainv,binv,cinv,alpinv,betinv,gaminv
    
def cell2AB(cell):
    #from real lattice parameters - cell
    # returns A for Cartesian to crystal transformations A*X = x 
    # and inverse B for crystal to Cartesian transformation B*x = X
    G,g = cell2Gmat(cell)       #reciprocal & real metric tensors
    cosAlpStar = G[2][1]/math.sqrt(G[1][1]*G[2][2])
    sinAlpStar = math.sqrt(1.0-cosAlpStar**2)
    B = np.eye(3)
    B *= cell[:3]
    A = np.zeros(shape=(3,3))
    A[0][0] = 1.0
    A[0][1] = cosd(cell[5])
    A[1][1] = sinAlpStar*sind(cell[5])
    A[1][2] = -cosAlpStar*sind(cell[5])
    A[0][2] = cosd(cell[4])
    A[2][2] = sind(cell[4])
    B = np.dot(A,B)
    A = nl.inv(B)
    return A,B
    
def makeMat(Angle,Axis):
    #Make rotation matrix from Angle in degrees,Axis =0 for rotation about x, =1 for about y, etc.
    cs = cosd(Angle)
    ss = sind(Angle)
    M = np.array(([1.,0.,0.],[0.,cs,-ss],[0.,ss,cs]),dtype=np.float32)
    return np.roll(np.roll(M,Axis,axis=0),Axis,axis=1)
                    
def MaxIndex(dmin,A):
    #finds maximum allowed hkl for given A within dmin
    Hmax = [0,0,0]
    try:
        cell = A2cell(A)
    except:
        cell = [1,1,1,90,90,90]
    for i in range(3):
        Hmax[i] = int(round(cell[i]/dmin))
    return Hmax
    
def sortHKLd(HKLd,ifreverse,ifdup):
    #HKLd is a list of [h,k,l,d,...]; ifreverse=True for largest d first
    #ifdup = True if duplicate d-spacings allowed 
    T = []
    for i,H in enumerate(HKLd):
        if ifdup:
            T.append((H[3],i))
        else:
            T.append(H[3])            
    D = dict(zip(T,HKLd))
    T.sort()
    if ifreverse:
        T.reverse()
    X = []
    okey = ''
    for key in T: 
        if key != okey: X.append(D[key])    #remove duplicate d-spacings
        okey = key
    return X
    
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
                
 
def GenHBravais(dmin,Bravais,A):
# dmin - minimum d-spacing
# Bravais in range(14) to indicate Bravais lattice; 0-2 cubic, 3,4 - hexagonal/trigonal,
# 5,6 - tetragonal, 7-10 - orthorhombic, 11,12 - monoclinic, 13 - triclinic
# A - as defined in calc_rDsq
# returns HKL = [h,k,l,d,0] sorted so d largest first 
    import math
    if Bravais in [9,11]:
        Cent = 'C'
    elif Bravais in [1,5,8]:
        Cent = 'I'
    elif Bravais in [0,7]:
        Cent = 'F'
    elif Bravais in [3]:
        Cent = 'R'
    else:
        Cent = 'P'
    Hmax = MaxIndex(dmin,A)
    dminsq = 1./(dmin**2)
    HKL = []
    if Bravais == 13:                       #triclinic
        for l in range(-Hmax[2],Hmax[2]+1):
            for k in range(-Hmax[1],Hmax[1]+1):
                hmin = 0
                if (k < 0): hmin = 1
                if (k ==0 and l < 0): hmin = 1
                for h in range(hmin,Hmax[0]+1):
                    H=[h,k,l]
                    rdsq = calc_rDsq(H,A)
                    if 0 < rdsq <= dminsq:
                        HKL.append([h,k,l,rdsq2d(rdsq,6),-1])
    elif Bravais in [11,12]:                #monoclinic - b unique
        Hmax = SwapIndx(2,Hmax)
        for h in range(Hmax[0]+1):
            for k in range(-Hmax[1],Hmax[1]+1):
                lmin = 0
                if k < 0:lmin = 1
                for l in range(lmin,Hmax[2]+1):
                    [h,k,l] = SwapIndx(-2,[h,k,l])
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,rdsq2d(rdsq,6),-1])
                    [h,k,l] = SwapIndx(2,[h,k,l])
    elif Bravais in [7,8,9,10]:            #orthorhombic
        for h in range(Hmax[0]+1):
            for k in range(Hmax[1]+1):
                for l in range(Hmax[2]+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,rdsq2d(rdsq,6),-1])
    elif Bravais in [5,6]:                  #tetragonal
        for l in range(Hmax[2]+1):
            for k in range(Hmax[1]+1):
                for h in range(k,Hmax[0]+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,rdsq2d(rdsq,6),-1])
    elif Bravais in [3,4]:
        lmin = 0
        if Bravais == 3: lmin = -Hmax[2]                  #hexagonal/trigonal
        for l in range(lmin,Hmax[2]+1):
            for k in range(Hmax[1]+1):
                hmin = k
                if l < 0: hmin += 1
                for h in range(hmin,Hmax[0]+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,rdsq2d(rdsq,6),-1])

    else:                                   #cubic
        for l in range(Hmax[2]+1):
            for k in range(l,Hmax[1]+1):
                for h in range(k,Hmax[0]+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,rdsq2d(rdsq,6),-1])
    return sortHKLd(HKL,True,False)
    
def SwapXY(x,y):
    return [y,x]
    
def SwapIndx(Axis,H):
    if Axis in [1,-1]:
        return H
    elif Axis in [2,-3]:
        return [H[1],H[2],H[0]]
    else:
        return [H[2],H[0],H[1]]
        
def CentCheck(Cent,H):
    h,k,l = H
    if Cent == 'A' and (k+l)%2:
        return False
    elif Cent == 'B' and (h+l)%2:
        return False
    elif Cent == 'C' and (h+k)%2:
        return False
    elif Cent == 'I' and (h+k+l)%2:
        return False
    elif Cent == 'F' and ((h+k)%2 or (h+l)%2 or (k+l)%2):
        return False
    elif Cent == 'R' and (-h+k+l)%3:
        return False
    else:
        return True
                                    
def GenHLaue(dmin,Laue,Cent,Axis,A):
# dmin - minimum d-spacing
# Laue - Laue group symbol = '-1','2/m','mmmm','4/m','6/m','4/mmm','6/mmm',
#                            '3m1', '31m', '3', '3R', '3mR', 'm3', 'm3m'
# Cent - lattice centering = 'P','A','B','C','I','F'
# Axis - code for unique monoclinic axis = 'a','b','c'
# A - 6 terms as defined in calc_rDsq
# returns - HKL = list of [h,k,l,d] sorted with largest d first and is unique 
# part of reciprocal space ignoring anomalous dispersion
    import math
    Hmax = MaxIndex(dmin,A)
    dminsq = 1./(dmin**2)
    HKL = []
    if Laue == '-1':                       #triclinic
        for l in range(-Hmax[2],Hmax[2]+1):
            for k in range(-Hmax[1],Hmax[1]+1):
                hmin = 0
                if (k < 0) or (k ==0 and l < 0): hmin = 1
                for h in range(hmin,Hmax[0]+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    rdsq = calc_rDsq(H,A)
                    if 0 < rdsq <= dminsq:
                        HKL.append([h,k,l,1/math.sqrt(rdsq)])
    elif Laue == '2/m':                #monoclinic
        Hmax = SwapIndx(Axis,Hmax)
        for h in range(Hmax[0]+1):
            for k in range(-Hmax[1],Hmax[1]+1):
                lmin = 0
                if k < 0:lmin = 1
                for l in range(lmin,Hmax[2]+1):
                    [h,k,l] = SwapIndx(-Axis,[h,k,l])
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1/math.sqrt(rdsq)])
                    [h,k,l] = SwapIndx(Axis,[h,k,l])
    elif Laue in ['mmm','4/m','6/m']:            #orthorhombic
        for l in range(Hmax[2]+1):
            for h in range(Hmax[0]+1):
                kmin = 1
                if Laue == 'mmm' or h ==0: kmin = 0
                for k in range(kmin,Hmax[2]+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1/math.sqrt(rdsq)])
    elif Laue in ['4/mmm','6/mmm']:                  #tetragonal & hexagonal
        for l in range(Hmax[2]+1):
            for h in range(Hmax[0]+1):
                for k in range(h+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1/math.sqrt(rdsq)])
    elif Laue in ['3m1','31m','3','3R','3mR']:                  #trigonals
        for l in range(-Hmax[2],Hmax[2]+1):
            hmin = 0
            if l < 0: hmin = 1
            for h in range(hmin,Hmax[0]+1):
                if Laue in ['3R','3']:
                    kmax = h
                    kmin = -int((h-1)/2.)
                else:
                    kmin = 0
                    kmax = h
                    if Laue in ['3m1','3mR'] and l < 0: kmax = h-1
                    if Laue == '31m' and l < 0: kmin = 1
                for k in range(kmin,kmax+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1/math.sqrt(rdsq)])
    else:                                   #cubic
        for h in range(Hmax[0]+1):
            for k in range(h+1):
                lmin = 0
                lmax = k
                if Laue =='m3':
                    lmax = h-1
                    if h == k: lmax += 1
                for l in range(lmin,lmax+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1/math.sqrt(rdsq)])
    return sortHKLd(HKL,True,True)
    
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
    
def IndexPeaks(peaks,HKL):
    import bisect
    hklds = [1000.0]                                    #make bounded list of available d-spacings
    N = len(HKL)
    if N == 0: return False
    for hkl in HKL:
        hklds.append(hkl[3])
    hklds.append(0.0)
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
                    for j in range(3):
                        opeak[j+4] = 0
                    opeak[8] = 0.
                else:                                       # old better - do nothing
                    continue                
            hkl[4] = ipk
            for j in range(3):
                peak[j+4] = hkl[j]
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
    
def FitHKL(ibrav,peaks,A,wtP):
    def ShiftTest(a,b):
        if b < -0.1*a: 
            b = -0.0001*a
        return b
    smin = 0.
    first = True
    for peak in peaks:
        if peak[2] and peak[3]:
            h,k,l = H = peak[4:7]
            Qo = 1./peak[7]**2
            Qc = calc_rDsq(H,A)
            try:
                peak[8] = 1./math.sqrt(Qc)
            except:
                print A2invcell(A)
            delt = Qo-Qc
            smin += delt**2
            dp = []
            if ibrav in [0,1,2]:            #m3m
                dp.append(h*h+k*k+l*l)
            elif ibrav in [3,4]:            #R3H, P3/m & P6/mmm
                dp.append(h*h+k*k+h*k)
                dp.append(l*l)
            elif ibrav in [5,6]:            #4/mmm
                dp.append(h*h+k*k)
                dp.append(l*l)
            elif ibrav in [7,8,9,10]:       #mmm
                dp.append(h*h)
                dp.append(k*k)
                dp.append(l*l)
            elif ibrav in [11,12]:          #2/m
                dp.append(h*h)
                dp.append(k*k)
                dp.append(l*l)
                dp.append(h*l)
            else:                           #1
#    derivatives for a0*h^2+a1*k^2+a2*l^2+a3*h*k+a4*h*l+a5*k*l
                dp.append(h*h)
                dp.append(k*k)
                dp.append(l*l)
                dp.append(h*k)
                dp.append(h*l)
                dp.append(k*l)
            if first:
                first = False
                M = len(dp)
                B = np.zeros(shape=(M,M))
                V = np.zeros(shape=(M))
            dwt = peak[7]**wtP
            B,V = pyp.buildmv(delt*dwt,dwt,M,dp,B,V)
    if nl.det(B) > 0.0:
        try:
            b = nl.solve(B,V)
            B = nl.inv(B)
            sig = np.diag(B)
        except SingularMatrix:
            return False,0
        if ibrav in [0,1,2]:            #m3m
            A[0] += ShiftTest(A[0],b[0])
            A[1] = A[2] = A[0]
        elif ibrav in [3,4]:            #R3H, P3/m & P6/mmm
            A[0] += ShiftTest(A[0],b[0])
            A[1] = A[3] = A[0]
            A[2] += ShiftTest(A[2],b[1])
        elif ibrav in [5,6]:            #4/mmm
            A[0] += ShiftTest(A[0],b[0])
            A[1] = A[0]
            A[2] += ShiftTest(A[2],b[1])
        elif ibrav in [7,8,9,10]:       #mmm
            for i in range(3):
                A[i] += ShiftTest(A[i],b[i])
        elif ibrav in [11,12]:          #2/m
            for i in range(3):
                A[i] += ShiftTest(A[i],b[i])
            A[4] += ShiftTest(A[4],b[3])
            A[4] = min(1.4*math.sqrt(A[0]*A[2]),A[4])   #min beta star = 45
        else:                           #1
            oldV = math.sqrt(1./calc_rVsq(A))
            oldA = A[:]
            for i in range(6):
                A[i] += b[i]*0.2
            A[3] = min(1.1*math.sqrt(max(0,A[1]*A[2])),A[3])
            A[4] = min(1.1*math.sqrt(max(0,A[0]*A[2])),A[4])
            A[5] = min(1.1*math.sqrt(max(0,A[0]*A[1])),A[5])
            ratio = math.sqrt(1./calc_rVsq(A))/oldV
            if 0.9 > ratio or ratio > 1.1:
                A = oldA
#    else:
#        print B
#        print V
#        for peak in peaks:
#            print peak
    return True,smin
       
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
    pwr = 6
    maxTries = 10
    if ibrav == 13:
        pwr = 4
        maxTries = 10
    OK = False
    tries = 0
    HKL = GenHBravais(dmin,ibrav,A)
    while IndexPeaks(peaks,HKL):
        Pwr = pwr - 2*(tries % 2)
        HKL = []
        tries += 1
        osmin = smin
        oldA = A
        OK,smin = FitHKL(ibrav,peaks,A,Pwr)
        for a in A[:3]:
            if a < 0:
                A = oldA
                OK = False
                break
        if OK:
            HKL = GenHBravais(dmin,ibrav,A)
        if len(HKL) == 0: break                         #absurd cell obtained!
        rat = (osmin-smin)/smin
        if abs(rat) < 1.0e-5 or not OK: break
        if tries > maxTries: break
    if OK:
        OK,smin = FitHKL(ibrav,peaks,A,4)
    M20,X20 = calc_M20(peaks,HKL)
    return len(HKL),M20,X20
        
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
        HKL = GenHBravais(dmin,ibrav,A[:])
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
        HKL = GenHBravais(dmin,ibrav,Anew)
        if len(HKL) > mHKL[ibrav] and IndexPeaks(peaks,HKL):
            Lhkl,M20,X20 = refinePeaks(peaks,ibrav,Anew)
            Asave.append([calc_M20(peaks,HKL),Anew[:]])
            if ibrav == 9:                          #C-centered orthorhombic
                for i in range(2):
                    Anew = rotOrthoA(Anew[:])
                    Lhkl,M20,X20 = refinePeaks(peaks,ibrav,Anew)
                    HKL = GenHBravais(dmin,ibrav,Anew)
                    IndexPeaks(peaks,HKL)
                    Asave.append([calc_M20(peaks,HKL),Anew[:]])
            elif ibrav == 11:                      #C-centered monoclinic
                Anew = swapMonoA(Anew[:])
                Lhkl,M20,X20 = refinePeaks(peaks,ibrav,Anew)
                HKL = GenHBravais(dmin,ibrav,Anew)
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
        Lhkl,M20,X20 = refinePeaks(peaks,ibrav,X[0][1])
        return GoOn,Lhkl,M20,X20,X[0][1]
    else:
        return GoOn,0,0,0,0
        
def monoCellReduce(ibrav,A):
    a,b,c,alp,bet,gam = A2cell(A)
    G,g = A2Gmat(A)
    if ibrav in [11]:
        u = [0,0,-1]
        v = [1,0,2]
        anew = math.sqrt(np.dot(np.dot(v,g),v))
        if anew < a:
            cang = np.dot(np.dot(u,g),v)/(anew*c)
            beta = acosd(-abs(cang))
            A = cell2A([anew,b,c,90,beta,90])
    else:
        u = [-1,0,0]
        v = [1,0,1]
        cnew = math.sqrt(np.dot(np.dot(v,g),v))
        if cnew < c:
            cang = np.dot(np.dot(u,g),v)/(a*cnew)
            beta = acosd(-abs(cang))
            A = cell2A([a,b,cnew,90,beta,90])
    return A

def DoIndexPeaks(peaks,inst,controls,bravais):
    
    def peakDspace(peaks,A):
        for peak in peaks:
            if peak[3]:
                dsq = calc_rDsq(peak[4:7],A)
                if dsq > 0:
                    peak[8] = 1./math.sqrt(dsq)
        return
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
                screenSize = wx.DisplaySize()
                Size = dlg.GetSize()
                dlg.SetPosition(wx.Point(screenSize[0]-Size[0]-300,0))
                try:
                    GoOn = True
                    while GoOn:                                                 #Loop over increment of volume
                        N2 = 0
                        while N2 < N2s[ibrav]:                                  #Table 2 step (iii)               
                            if ibrav > 2:
                                if not N2:
                                    GoOn,Nc,M20,X20,A = findBestCell(dlg,ncMax,0,Nm[ibrav]*N1s[ibrav],ibrav,peaks,V1)
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
                                    HKL = GenHBravais(dmin,ibrav,A)
                                    IndexPeaks(peaks,HKL)
                                    a,b,c,alp,bet,gam = A2cell(A)
                                    V = calc_V(A)
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
                ', elapsed time = ',sec2HMS(time.time()-begin))
            
    if cells:
        cells = sortM20(cells)
        cells[0][-1] = True
        return True,dmin,cells
    else:
        return False,0,0
        
def FitRing(ring):
    Err,parms = FitCircle(ring)
    Err /= len(ring)
#    print 'circle error:','%8f'%(Err)
    if Err > 20000.:
        eparms = FitEllipse(ring)
        if eparms:
            parms = eparms
    return parms
        
def FitCircle(ring):
    import numpy.linalg as nl
    
    def makeParmsCircle(B):
        cent = [-B[0]/2,-B[1]/2]
        phi = 0.
        sr1 = sr2 = math.sqrt(cent[0]**2+cent[1]**2-B[2])
        return cent,phi,[sr1,sr2]
        
    ring = np.array(ring)
    x = np.asarray(ring.T[0])
    y = np.asarray(ring.T[1])
    
    M = np.array((x,y,np.ones_like(x)))
    B = np.array(-(x**2+y**2))
    result = nl.lstsq(M.T,B)
    return result[1],makeParmsCircle(result[0])
        
def FitEllipse(ring):
    import numpy.linalg as nl
            
    def makeParmsEllipse(B):
        det = 4.*(1.-B[0]**2)-B[1]**2
        if det < 0.:
            print 'hyperbola!'
            return 0
        elif det == 0.:
            print 'parabola!'
            return 0
        cent = [(B[1]*B[3]-2.*(1.-B[0])*B[2])/det, \
            (B[1]*B[2]-2.*(1.+B[0])*B[3])/det]
        phi = 0.5*atand(0.5*B[1]/B[0])
        
        a = (1.+B[0])/cosd(2*phi)
        b = 2.-a
        f = (1.+B[0])*cent[0]**2+(1.-B[0])*cent[1]**2+B[1]*cent[0]*cent[1]-B[4]
        if f/a < 0 or f/b < 0:
            return 0
        sr1 = math.sqrt(f/a)
        sr2 = math.sqrt(f/b)
        if sr1 > sr2:
            sr1,sr2 = SwapXY(sr1,sr2)
            phi -= 90.
            if phi < -90.:
                phi += 180.
        return cent,phi,[sr1,sr2]
                
    ring = np.array(ring)
    x = np.asarray(ring.T[0])
    y = np.asarray(ring.T[1])
    M = np.array((x**2-y**2,x*y,x,y,np.ones_like(x)))
    B = np.array(-(x**2+y**2))
    result = nl.lstsq(M.T,B)
    return makeParmsEllipse(result[0])
    
def FitDetector(rings,p0,wave):
    from scipy.optimize import leastsq
    def ellipseCalc(B,xyd,wave):
        x = xyd[0]
        y = xyd[1]
        dsp = xyd[2]
        dist,x0,y0,phi,tilt = B
        tth = 2.0*npasind(wave/(2.*dsp))
        ttth = nptand(tth)
        radius = dist*ttth
        stth = npsind(tth)
        cosb = npcosd(tilt)
        R1 = dist*stth*npcosd(tth)*cosb/(cosb**2-stth**2)
        R0 = np.sqrt(R1*radius*cosb)
        zdis = R1*ttth*nptand(tilt)
        X = x-x0+zdis*npsind(phi)
        Y = y-y0-zdis*npcosd(phi)
        XR = X*npcosd(phi)-Y*npsind(phi)
        YR = X*npsind(phi)+Y*npcosd(phi)
        return (XR/R0)**2+(YR/R1)**2-1
    result = leastsq(ellipseCalc,p0,args=(rings.T,wave))
    return result[0]
            
def ImageLocalMax(image,w,Xpix,Ypix):
    w2 = w*2
    size = len(image)
    xpix = int(Xpix)            #get reference corner of pixel chosen
    ypix = int(Ypix)
    if (w < xpix < size-w) and (w < ypix < size-w) and image[ypix,xpix]:
        Z = image[ypix-w:ypix+w,xpix-w:xpix+w]
        Zmax = np.argmax(Z)
        Zmin = np.argmin(Z)
        xpix += Zmax%w2-w
        ypix += Zmax/w2-w
        return xpix,ypix,np.ravel(Z)[Zmax],np.ravel(Z)[Zmin]
    else:
        return 0,0,0,0
    
def makeRing(dsp,ellipse,pix,reject,scalex,scaley,image):
    cent,phi,radii = ellipse
    cphi = cosd(phi)
    sphi = sind(phi)
    ring = []
    for a in range(-180,180,2):
        x = radii[0]*cosd(a)
        y = radii[1]*sind(a)
        X = (cphi*x-sphi*y+cent[0])*scalex      #convert mm to pixels
        Y = (sphi*x+cphi*y+cent[1])*scaley
        X,Y,I,J = ImageLocalMax(image,pix,X,Y)      
        if I and J and I/J > reject:
            X += .5                             #set to center of pixel
            Y += .5
            X /= scalex                         #convert to mm
            Y /= scaley
            ring.append([X,Y,dsp])
    if len(ring) < 45:             #want more than 1/4 of a circle
        return []
    return ring
    
def makeIdealRing(ellipse):
    cent,phi,radii = ellipse
    cphi = cosd(phi)
    sphi = sind(phi)
    ring = []
    for a in range(0,360,2):
        x = radii[0]*cosd(a)
        y = radii[1]*sind(a)
        X = (cphi*x-sphi*y+cent[0])
        Y = (sphi*x+cphi*y+cent[1])
        ring.append([X,Y])
    return ring
    
def calcDist(radii,tth):
    stth = sind(tth)
    ctth = cosd(tth)
    ttth = tand(tth)
    return math.sqrt(radii[0]**4/(ttth**2*((radii[0]*ctth)**2+(radii[1]*stth)**2)))
    
def calcZdisCosB(radius,tth,radii):
    cosB = sinb = radii[0]**2/(radius*radii[1])
    if cosB > 1.:
        return 0.,1.
    else:
        cosb = math.sqrt(1.-sinb**2)
        ttth = tand(tth)
        zdis = radii[1]*ttth*cosb/sinb
        return zdis,cosB
    
def GetEllipse(dsp,data):
    dist = data['distance']
    cent = data['center']
    tilt = data['tilt']
    phi = data['rotation']
    radii = [0,0]
    tth = 2.0*asind(data['wavelength']/(2.*dsp))
    ttth = tand(tth)
    stth = sind(tth)
    ctth = cosd(tth)
    cosb = cosd(tilt)
    radius = dist*ttth
    radii[1] = dist*stth*ctth*cosb/(cosb**2-stth**2)
    if radii[1] > 0:
        radii[0] = math.sqrt(radii[1]*radius*cosb)
        zdis = radii[1]*ttth*tand(tilt)
        elcent = [cent[0]-zdis*sind(phi),cent[1]+zdis*cosd(phi)]
        return elcent,phi,radii
    else:
        return False
        
def GetDetectorXY(dsp,azm,data):
    from scipy.optimize import fsolve
    def func(xy,*args):
       azm,phi,R0,R1,A,B = args
       cp = cosd(phi)
       sp = sind(phi)
       x,y = xy
       out = []
       out.append(y-x*tand(azm))
       out.append(R0**2*((x+A)*sp-(y+B)*cp)**2+R1**2*((x+A)*cp+(y+B)*sp)**2-(R0*R1)**2)
       return out
    elcent,phi,radii = GetEllipse(dsp,data)
    cent = data['center']
    tilt = data['tilt']
    phi = data['rotation']
    wave = data['wavelength']
    dist = data['distance']
    tth = 2.0*asind(wave/(2.*dsp))
    ttth = tand(tth)
    radius = dist*ttth
    stth = sind(tth)
    cosb = cosd(tilt)
    R1 = dist*stth*cosd(tth)*cosb/(cosb**2-stth**2)
    R0 = math.sqrt(R1*radius*cosb)
    zdis = R1*ttth*tand(tilt)
    A = zdis*sind(phi)
    B = -zdis*cosd(phi)
    xy0 = [radius*cosd(azm),radius*sind(azm)]
    xy = fsolve(func,xy0,args=(azm,phi,R0,R1,A,B))+cent
    return xy
                    
def GetTthAzmDsp(x,y,data):
    wave = data['wavelength']
    dist = data['distance']
    cent = data['center']
    tilt = data['tilt']
    phi = data['rotation']
    dx = np.array(x-cent[0],dtype=np.float32)
    dy = np.array(y-cent[1],dtype=np.float32)
    X = np.array(([dx,dy,np.zeros_like(dx)]),dtype=np.float32).T
    X = np.dot(X,makeMat(phi,2))
    Z = np.dot(X,makeMat(tilt,0)).T[2]
    tth = npatand(np.sqrt(dx**2+dy**2-Z**2)/(dist-Z))
    dsp = wave/(2.*npsind(tth/2.))
    azm = npatan2d(dy,dx)
    return tth,azm,dsp
    
def GetTth(x,y,data):
    return GetTthAzmDsp(x,y,data)[0]
    
def GetTthAzm(x,y,data):
    return GetTthAzmDsp(x,y,data)[0:2]
    
def GetDsp(x,y,data):
    return GetTthAzmDsp(x,y,data)[2]
       
def ImageCompress(image,scale):
    if scale == 1:
        return image
    else:
        return image[::scale,::scale]
        
def ImageCalibrate(self,data):
    import copy
    import ImageCalibrants as calFile
    print 'image calibrate'
    ring = data['ring']
    pixelSize = data['pixelSize']
    scalex = 1000./pixelSize[0]
    scaley = 1000./pixelSize[1]
    cutoff = data['cutoff']
    if len(ring) < 5:
        print 'not enough inner ring points for ellipse'
        return False
        
    #fit start points on inner ring
    data['ellipses'] = []
    outE = FitRing(ring)
    if outE:
        print 'start ellipse:',outE
        ellipse = outE
    else:
        return False
        
    #setup 180 points on that ring for "good" fit
    Ring = makeRing(1.0,ellipse,20,cutoff,scalex,scaley,self.ImageZ)
    if Ring:
        ellipse = FitRing(Ring)
        Ring = makeRing(1.0,ellipse,20,cutoff,scalex,scaley,self.ImageZ)    #do again
        ellipse = FitRing(Ring)
    else:
        print '1st ring not sufficiently complete to proceed'
        return False
    print 'inner ring:',ellipse
    data['center'] = copy.copy(ellipse[0])           #not right!! (but useful for now)
    data['ellipses'].append(ellipse[:]+('r',))
    G2plt.PlotImage(self)
    
    #setup for calibration
    data['rings'] = []
    data['ellipses'] = []
    if not data['calibrant']:
        print 'no calibration material selected'
        return True
        
    Bravais,cell = calFile.Calibrants[data['calibrant']]
    A = cell2A(cell)
    wave = data['wavelength']
    cent = data['center']
    pixLimit = data['pixLimit']
    elcent,phi,radii = ellipse
    HKL = GenHBravais(0.5,Bravais,A)
    dsp = HKL[0][3]
    tth = 2.0*asind(wave/(2.*dsp))
    ttth = tand(tth)
    data['distance'] = dist = calcDist(radii,tth)
    radius = dist*tand(tth)
    zdis,cosB = calcZdisCosB(radius,tth,radii)
    cent1 = []
    cent2 = []
    xSum = 0
    ySum = 0
    zxSum = 0
    zySum = 0
    phiSum = 0
    tiltSum = 0
    distSum = 0
    Zsum = 0
    for i,H in enumerate(HKL):
        dsp = H[3]
        tth = 2.0*asind(0.5*wave/dsp)
        stth = sind(tth)
        ctth = cosd(tth)
        ttth = tand(tth)
        radius = dist*ttth
        elcent,phi,radii = ellipse
        radii[1] = dist*stth*ctth*cosB/(cosB**2-stth**2)
        radii[0] = math.sqrt(radii[1]*radius*cosB)
        zdis,cosB = calcZdisCosB(radius,tth,radii)
        zsinp = zdis*sind(phi)
        zcosp = zdis*cosd(phi)
        cent = data['center']
        elcent = [cent[0]+zsinp,cent[1]-zcosp]
        ratio = radii[1]/radii[0]
        Ring = makeRing(dsp,ellipse,pixLimit,cutoff,scalex,scaley,self.ImageZ)
        if Ring:
            numZ = len(Ring)
            data['rings'].append(np.array(Ring))
            ellipse = FitRing(Ring)
            elcent,phi,radii = ellipse                
            if abs(phi) > 45. and phi < 0.:
                phi += 180.
            dist = calcDist(radii,tth)
            distR = 1.-dist/data['distance']
            if distR > 0.001:
                print 'Wavelength too large?'
            elif distR < -0.001:
                print 'Wavelength too small?'
            else:
                if abs((radii[1]/radii[0]-ratio)/ratio) > 0.01:
                    print 'Bad fit for ring # %i. Try reducing Pixel search range'%(i)
                    return False
            zdis,cosB = calcZdisCosB(radius,tth,radii)
            Tilt = acosd(cosB)          # 0 <= tilt <= 90
            zsinp = zdis*sind(ellipse[1])
            zcosp = zdis*cosd(ellipse[1])
            cent1.append(np.array([elcent[0]+zsinp,elcent[1]-zcosp]))
            cent2.append(np.array([elcent[0]-zsinp,elcent[1]+zcosp]))
            if i:
                d1 = cent1[-1]-cent1[-2]        #get shift of 2 possible center solutions
                d2 = cent2[-1]-cent2[-2]
                if np.dot(d2,d2) > np.dot(d1,d1):  #right solution is the larger shift
                    data['center'] = cent1[-1]
                else:
                    data['center'] = cent2[-1]
                Zsum += numZ
                phiSum += numZ*phi
                distSum += numZ*dist
                xSum += numZ*data['center'][0]
                ySum += numZ*data['center'][1]
                tiltSum += numZ*abs(Tilt)
            cent = data['center']
            print ('for ring # %2i dist %.3f rotate %6.2f tilt %6.2f Xcent %.3f Ycent %.3f Npts %d' 
                %(i,dist,phi,Tilt,cent[0],cent[1],numZ))
            data['ellipses'].append(copy.deepcopy(ellipse+('r',)))
            G2plt.PlotImage(self)
        else:
            break
    fullSize = len(self.ImageZ)/scalex
    if 2*radii[1] < .9*fullSize:
        print 'Are all usable rings (>25% visible) used? Try reducing Min ring I/Ib'
    if not Zsum:
        print 'Only one ring fitted. Check your wavelength.'
        return False
    cent = data['center'] = [xSum/Zsum,ySum/Zsum]
    wave = data['wavelength']
    dist = data['distance'] = distSum/Zsum
    
    #possible error if no. of rings < 3! Might need trap here
    d1 = cent1[-1]-cent1[1]             #compare last ring to 2nd ring
    d2 = cent2[-1]-cent2[1]
    Zsign = 1
    len1 = math.sqrt(np.dot(d1,d1))
    len2 = math.sqrt(np.dot(d2,d2))
    t1 = d1/len1
    t2 = d2/len2
    if len2 > len1:
        if -135. < atan2d(t2[1],t2[0]) < 45.:
            Zsign = -1
    else:
        if -135. < atan2d(t1[1],t1[0]) < 45.:
            Zsign = -1
    
    tilt = data['tilt'] = Zsign*tiltSum/Zsum
    phi = data['rotation'] = phiSum/Zsum
    rings = np.concatenate((data['rings']),axis=0)
    p0 = [dist,cent[0],cent[1],phi,tilt]
    result = FitDetector(rings,p0,wave)
    data['distance'] = result[0]
    data['center'] = result[1:3]
    data['rotation'] = np.mod(result[3],360.0)
    data['tilt'] = result[4]
    N = len(data['ellipses'])
    data['ellipses'] = []           #clear away individual ellipse fits
    for H in HKL[:N]:
        ellipse = GetEllipse(H[3],data)
        data['ellipses'].append(copy.deepcopy(ellipse+('b',)))
    G2plt.PlotImage(self)        
    return True
    
def Make2ThetaAzimuthMap(data,imageN):
    #transforms 2D image from x,y space to 2-theta,azimuth space based on detector orientation
    pixelSize = data['pixelSize']
    scalex = pixelSize[0]/1000.
    scaley = pixelSize[1]/1000.
    tax,tay = np.mgrid[0.5:imageN+.5,0.5:imageN+.5]         #bin centers not corners
    tax = np.asfarray(tax*scalex,dtype=np.float32)
    tay = np.asfarray(tay*scaley,dtype=np.float32)
    return GetTthAzm(tay,tax,data)           #2-theta & azimuth arrays

def Fill2ThetaAzimuthMap(data,TA,image):
    import numpy.ma as ma
    LUtth = data['IOtth']
    if data['fullIntegrate']:
        LRazm = [-180,180]
    else:
        LRazm = data['LRazimuth']
    imageN = len(image)
    TA = np.reshape(TA,(2,imageN,imageN))
    TA = np.dstack((ma.getdata(TA[1]),ma.getdata(TA[0])))    #azimuth, 2-theta
    tax,tay = np.dsplit(TA,2)    #azimuth, 2-theta
    tax = ma.masked_outside(tax.flatten(),LRazm[0],LRazm[1])
    tay = ma.masked_outside(tay.flatten(),LUtth[0],LUtth[1])
    tam = ma.getmask(tax)+ma.getmask(tay)
    taz = ma.masked_where(tam,image.flatten())
    return tax,tay,taz,tam
    
def Bin2ThetaAzimuthMap(data,tax,tay,taz):
    import numpy.ma as ma
    LUtth = data['IOtth']
    if data['fullIntegrate']:
        LRazm = [-180,180]
    else:
        LRazm = data['LRazimuth']
    numAzms = data['outAzimuths']
    numChans = data['outChannels']
    NST = np.histogram2d(tax,tay,normed=False,bins=(numAzms,numChans),range=[LRazm,LUtth])
    HST = np.histogram2d(tax,tay,normed=False,bins=(numAzms,numChans),range=[LRazm,LUtth],weights=taz)
    return NST,HST

def ImageIntegrate(self,data):
    dlg = wx.ProgressDialog("Elapsed time","2D image integration",5,
        style = wx.PD_ELAPSED_TIME|wx.PD_AUTO_HIDE)
    try:
        print 'Begin image integration'
        print 'Create 2-theta,azimuth map'
        t0 = time.time()
        dlg.Update(0)
        imageN = len(self.ImageZ)
        TA = Make2ThetaAzimuthMap(data,imageN)           #2-theta & azimuth arrays
        dlg.Update(1)
        print 'Fill map with 2-theta/azimuth values'
        tax,tay,taz,tam = Fill2ThetaAzimuthMap(data,TA,self.ImageZ)
        del TA
        dlg.Update(2)
        print 'Bin image by 2-theta/azimuth intervals'
        NST,HST = Bin2ThetaAzimuthMap(data,tax,tay,taz)
        del tax,tay,taz
        dlg.Update(3)
        print 'Form normalized 1D pattern(s)'
        self.Integrate = [HST[0]/NST[0],HST[1],HST[2]]
        del NST,HST
        dlg.Update(4)
        t1 = time.time()
        print 'Integration complete'
        print "Elapsed time:","%8.3f"%(t1-t0), "s"
    finally:
        dlg.Destroy()
       
def test():
    cell = [5,5,5,90,90,90]
    A = cell2A(cell)
    assert ( calc_V(A) == 125. )
    
    print 'test passed'

print __name__
if __name__ == '__main__':
    test()
    
        