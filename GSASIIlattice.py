'''Perform lattice-related computations'''

import numpy as np
import numpy.linalg as nl

# trig functions in degrees
sind = lambda x: np.sin(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
tand = lambda x: np.tan(x*np.pi/180.)
atand = lambda x: 180.*np.arctan(x)/np.pi
atan2d = lambda y,x: 180.*np.atan2(y,x)/np.pi
cosd = lambda x: np.cos(x*np.pi/180.)
acosd = lambda x: 180.*np.arccos(x)/np.pi
rdsq2d = lambda x,p: round(1.0/np.sqrt(x),p)

def sec2HMS(sec):
    H = int(sec/3600)
    M = int(sec/60-H*60)
    S = sec-3600*H-60*M
    return '%d:%2d:%.2f'%(H,M,S)

def fillgmat(cell):
    '''Compute lattice metric tensor from unit cell constants
    cell is tuple with a,b,c,alpha, beta, gamma (degrees)
    returns 3x3 numpy array
    '''
    a,b,c,alp,bet,gam = cell
    g = np.array([
        [a*a,  a*b*cosd(gam),  a*c*cosd(bet)],
        [a*b*cosd(gam),  b*b,  b*c*cosd(alp)],
        [a*c*cosd(bet) ,b*c*cosd(alp),   c*c]])
    return g
           
def cell2Gmat(cell):
    '''Compute real and reciprocal lattice metric tensor from unit cell constants
    cell is tuple with a,b,c,alpha, beta, gamma (degrees)
    returns reciprocal (G) & real (g) metric tensors (list of two 3x3 arrays)
    '''
    g = fillgmat(cell)
    G = nl.inv(g)        
    return G,g

def A2Gmat(A):
    '''Fill reciprocal metric tensor (G) from A
    returns reciprocal (G) & real (g) metric tensors (list of two 3x3 arrays)
    '''
    G = np.zeros(shape=(3,3))
    G = [
        [A[0],  A[3]/2.,  A[4]/2.], 
        [A[3]/2.,A[1],    A[5]/2.], 
        [A[4]/2.,A[5]/2.,    A[2]]]
    g = nl.inv(G)
    return G,g

def Gmat2A(G):
    'Extract A from reciprocal metric tensor (G)'
    return [G[0][0],G[1][1],G[2][2],2.*G[0][1],2.*G[0][2],2.*G[1][2]]
    
def cell2A(cell):
    G,g = cell2Gmat(cell)
    return Gmat2A(G)

def A2cell(A):
    '''Compute unit cell constants from A tensor
    returns tuple with a,b,c,alpha, beta, gamma (degrees)
    '''
    G,g = A2Gmat(A)
    return Gmat2cell(g)

def Gmat2cell(g):
    '''Compute lattice parameters from real metric tensor (g)
    returns tuple with a,b,c,alpha, beta, gamma (degrees)
    Alternatively,compute reciprocal lattice parameters from inverse metric tensor (G)
    returns tuple with a*,b*,c*,alpha*, beta*, gamma* (degrees)
    '''
    a = np.sqrt(max(0,g[0][0]))
    b = np.sqrt(max(0,g[1][1]))
    c = np.sqrt(max(0,g[2][2]))
    alp = acosd(g[2][1]/(b*c))
    bet = acosd(g[2][0]/(a*c))
    gam = acosd(g[0][1]/(a*b))
    return a,b,c,alp,bet,gam

def invcell2Gmat(invcell):
    '''Compute real and reciprocal lattice metric tensor from reciprocal 
       unit cell constants
    invcell is tuple with a*,b*,c*,alpha*, beta*, gamma* (degrees)
    returns reciprocal (G) & real (g) metric tensors (list of two 3x3 arrays)
    '''
    G = fillgmat(invcell)
    g = nl.inv(G)
    return G,g
        
def calc_rVsq(A):
    'Compute the square of the reciprocal lattice volume (V* **2) from A'
    G,g = A2Gmat(A)
    rVsq = nl.det(G)
    if rVsq < 0:
        return 1
    return rVsq
    
def calc_rV(A):
    'Compute the reciprocal lattice volume (V*) from A'
    return np.sqrt(calc_rVsq(A))
    
def calc_V(A):
    'Compute the real lattice volume (V) from A'
    return 1./calc_rV(A)

def A2invcell(A):
    '''Compute reciprocal unit cell constants from A
    returns tuple with a*,b*,c*,alpha*, beta*, gamma* (degrees)
    '''
    G,g = A2Gmat(A)
    return Gmat2cell(G)

def cell2AB(cell):
    '''Computes orthogonalization matrix from unit cell constants
    cell is tuple with a,b,c,alpha, beta, gamma (degrees)
    returns list of two 3x3 numpy arrays
       A for crystal to Cartesian transformations A*x = X 
       B (inverse) for Cartesian to crystal transformation B*X = x
    '''
    G,g = cell2Gmat(cell) 
    cellstar = Gmat2cell(G)
    A = np.zeros(shape=(3,3))
    # from Giacovazzo (Fundamentals 2nd Ed.) p.75
    A[0][0] = cell[0]                # a
    A[0][1] = cell[1]*cosd(cell[5])  # b cos(gamma)
    A[0][2] = cell[2]*cosd(cell[4])  # c cos(beta)
    A[1][1] = cell[1]*sind(cell[5])  # b sin(gamma)
    A[1][2] = -cell[2]*cosd(cellstar[3])*sind(cell[4]) # - c cos(alpha*) sin(beta)
    A[2][2] = 1/cellstar[2]         # 1/c*
    B = nl.inv(A)
    return A,B

#reflection generation routines
#for these: H = [h,k,l]; A is as used in calc_rDsq; G - inv metric tensor, g - metric tensor; 
#           cell - a,b,c,alp,bet,gam in A & deg
   
def calc_rDsq(H,A):
    rdsq = H[0]*H[0]*A[0]+H[1]*H[1]*A[1]+H[2]*H[2]*A[2]+H[0]*H[1]*A[3]+H[0]*H[2]*A[4]+H[1]*H[2]*A[5]
    return rdsq
    
def calc_rDsqZ(H,A,Z,tth,lam):
    rpd = math.pi/180.
    rdsq = calc_rDsq(H,A)+Z*math.sin(tth*rpd)*2.0*rpd/(lam*lam)
    return rdsq
       
def MaxIndex(dmin,A):
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
    
def SwapIndx(Axis,H):
    if Axis in [1,-1]:
        return H
    elif Axis in [2,-3]:
        return [H[1],H[2],H[0]]
    else:
        return [H[2],H[0],H[1]]
        
def Rh2Hx(Rh):
    Hx = [0,0,0]
    Hx[0] = Rh[0]-Rh[1]
    Hx[1] = Rh[1]-Rh[2]
    Hx[2] = np.sum(Rh)
    return Hx
    
def Hx2Rh(Hx):
        Rh = [0,0,0]
        itk = -Hx[0]+Hx[1]+Hx[2]
        if itk%3 != 0:
            return 0        #error - not rhombohedral reflection
        else:
            Rh[1] = itk/3
            Rh[0] = Rh[1]+Hx[0]
            Rh[2] = Rh[1]-Hx[1]
            if Rh[0] < 0:
                for i in range(3):
                    Rh[i] = -Rh[i]
            return Rh
        
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
                                    
def GetBraviasNum(center,system):
    '''Determine the Bravais lattice number, as used in GenHBravais
         center = one of: P, C, I, F, R (see SGLatt from GSASIIspc.SpcGroup)
         lattice = is cubic, hexagonal, tetragonal, orthorhombic, trigonal (R)
             monoclinic, triclinic (see SGSys from GSASIIspc.SpcGroup)
       Returns a number between 0 and 13 
          or throws an exception if the setting is non-standard
       '''
    if center.upper() == 'F' and system.lower() == 'cubic':
        return 0
    elif center.upper() == 'I' and system.lower() == 'cubic':
        return 1
    elif center.upper() == 'P' and system.lower() == 'cubic':
        return 2
    elif center.upper() == 'R' and system.lower() == 'trigonal':
        return 3
    elif center.upper() == 'P' and system.lower() == 'hexagonal':
        return 4
    elif center.upper() == 'I' and system.lower() == 'tetragonal':
        return 5
    elif center.upper() == 'P' and system.lower() == 'tetragonal':
        return 6
    elif center.upper() == 'F' and system.lower() == 'orthorhombic':
        return 7
    elif center.upper() == 'I' and system.lower() == 'orthorhombic':
        return 8
    elif center.upper() == 'C' and system.lower() == 'orthorhombic':
        return 9
    elif center.upper() == 'P' and system.lower() == 'orthorhombic':
        return 10
    elif center.upper() == 'C' and system.lower() == 'monoclinic':
        return 11
    elif center.upper() == 'P' and system.lower() == 'monoclinic':
        return 12
    elif center.upper() == 'P' and system.lower() == 'triclinic':
        return 13
    raise ValueError,'non-standard Bravais lattice center=%s, cell=%s' % (center,system)

def GenHBravais(dmin,Bravais,A):
    '''Generate the positionally unique powder diffraction reflections 
    input:
       dmin is minimum d-space
       Bravais is 0-13 to indicate lattice type (see GetBraviasNum)
       A is reciprocal cell tensor (see Gmat2A or cell2A)
    returns:
       a list of tuples containing: h,k,l,d-space,-1   
    '''
# Bravais in range(14) to indicate Bravais lattice:
#   0 F cubic
#   1 I cubic
#   2 P cubic
#   3 R hexagonal (trigonal not rhombohedral)
#   4 P hexagonal
#   5 I tetragonal
#   6 P tetragonal
#   7 F orthorhombic
#   8 I orthorhombic
#   9 C orthorhombic
#  10 P orthorhombic
#  11 C monoclinic
#  12 P monoclinic
#  13 P triclinic
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
    
def GenHLaue(dmin,Laue,Cent,Axis,A):
    '''Generate the crystallographically unique powder diffraction reflections
    for a lattice and Bravais type
    '''
# dmin - minimum d-spacing
# Laue - Laue group symbol = '-1','2/m','mmm','4/m','6/m','4/mmm','6/mmm',
#                            '3m1', '31m', '3', '3R', '3mR', 'm3', 'm3m'
# Cent - lattice centering = 'P','A','B','C','I','F'
# Axis - code for unique monoclinic axis = 'a','b','c'
# A - 6 terms as defined in calc_rDsq
# returns - HKL = list of [h,k,l,d] sorted with largest d first and is unique 
# part of reciprocal space ignoring anomalous dispersion
    import math
    #finds maximum allowed hkl for given A within dmin
    if Laue in ['3R','3mR']:        #Rhombohedral axes
        Hmax = [0,0,0]
        cell = A2cell(A)
        aHx = cell[0]*math.sqrt(2.0*(1.0-cosd(cell[3])))
        cHx = cell[0]*math.sqrt(3.0*(1.0+2.0*cosd(cell[3])))
        Hmax[0] = Hmax[1] = int(round(aHx/dmin))
        Hmax[2] = int(round(cHx/dmin))
        print Hmax,aHx,cHx
    else:                           # all others
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
        axisnum = 1 + ['a','b','c'].index(Axis)
        Hmax = SwapIndx(axisnum,Hmax)
        for h in range(Hmax[0]+1):
            for k in range(-Hmax[1],Hmax[1]+1):
                lmin = 0
                if k < 0:lmin = 1
                for l in range(lmin,Hmax[2]+1):
                    [h,k,l] = SwapIndx(-axisnum,[h,k,l])
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1/math.sqrt(rdsq)])
                    [h,k,l] = SwapIndx(axisnum,[h,k,l])
    elif Laue in ['mmm','4/m','6/m']:            #orthorhombic
        for l in range(Hmax[2]+1):
            for h in range(Hmax[0]+1):
                kmin = 1
                if Laue == 'mmm' or h ==0: kmin = 0
                for k in range(kmin,Hmax[1]+1):
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
#                    kmin = -int((h-1.)/2.)
                    kmin = -(h-1)/2
                else:
                    kmin = 0
                    kmax = h
                    if Laue in ['3m1','3mR'] and l < 0: kmax = h-1
                    if Laue == '31m' and l < 0: kmin = 1
                for k in range(kmin,kmax+1):
                    H = []
                    if CentCheck(Cent,[h,k,l]): H=[h,k,l]
                    if Laue in ['3R','3mR']:
                        H = Hx2Rh(H)
                    if H:
                        rdsq = calc_rDsq(H,A)
                        if 0 < rdsq <= dminsq:
                            HKL.append([h,k,l,1/math.sqrt(rdsq)])
                            print H,1/math.sqrt(rdsq)
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
    
# output from uctbx computed on platform darwin on 2010-05-28
NeedTestData = True
def TestData():
    array = np.array
    global NeedTestData
    NeedTestData = False
    global CellTestData
    CellTestData = [
# cell, g, G, cell*, V, V*
  [(4, 4, 4, 90, 90, 90), 
   array([[  1.60000000e+01,   9.79717439e-16,   9.79717439e-16],
       [  9.79717439e-16,   1.60000000e+01,   9.79717439e-16],
       [  9.79717439e-16,   9.79717439e-16,   1.60000000e+01]]), array([[  6.25000000e-02,   3.82702125e-18,   3.82702125e-18],
       [  3.82702125e-18,   6.25000000e-02,   3.82702125e-18],
       [  3.82702125e-18,   3.82702125e-18,   6.25000000e-02]]), (0.25, 0.25, 0.25, 90.0, 90.0, 90.0), 64.0, 0.015625],
# cell, g, G, cell*, V, V*
  [(4.0999999999999996, 5.2000000000000002, 6.2999999999999998, 100, 80, 130), 
   array([[ 16.81      , -13.70423184,   4.48533243],
       [-13.70423184,  27.04      ,  -5.6887143 ],
       [  4.48533243,  -5.6887143 ,  39.69      ]]), array([[ 0.10206349,  0.05083339, -0.00424823],
       [ 0.05083339,  0.06344997,  0.00334956],
       [-0.00424823,  0.00334956,  0.02615544]]), (0.31947376387537696, 0.25189277536327803, 0.16172643497798223, 85.283666420376008, 94.716333579624006, 50.825714168082683), 100.98576357983838, 0.0099023858863968445],
# cell, g, G, cell*, V, V*
  [(3.5, 3.5, 6, 90, 90, 120), 
   array([[  1.22500000e+01,  -6.12500000e+00,   1.28587914e-15],
       [ -6.12500000e+00,   1.22500000e+01,   1.28587914e-15],
       [  1.28587914e-15,   1.28587914e-15,   3.60000000e+01]]), array([[  1.08843537e-01,   5.44217687e-02,   3.36690552e-18],
       [  5.44217687e-02,   1.08843537e-01,   3.36690552e-18],
       [  3.36690552e-18,   3.36690552e-18,   2.77777778e-02]]), (0.32991443953692895, 0.32991443953692895, 0.16666666666666669, 90.0, 90.0, 60.000000000000021), 63.652867178156257, 0.015710211406520427],
  ]
    global CoordTestData
    CoordTestData = [
# cell, ((frac, ortho),...)
  ((4,4,4,90,90,90,), [
 ((0.10000000000000001, 0.0, 0.0),(0.40000000000000002, 0.0, 0.0)),
 ((0.0, 0.10000000000000001, 0.0),(2.4492935982947065e-17, 0.40000000000000002, 0.0)),
 ((0.0, 0.0, 0.10000000000000001),(2.4492935982947065e-17, -2.4492935982947065e-17, 0.40000000000000002)),
 ((0.10000000000000001, 0.20000000000000001, 0.29999999999999999),(0.40000000000000013, 0.79999999999999993, 1.2)),
 ((0.20000000000000001, 0.29999999999999999, 0.10000000000000001),(0.80000000000000016, 1.2, 0.40000000000000002)),
 ((0.29999999999999999, 0.20000000000000001, 0.10000000000000001),(1.2, 0.80000000000000004, 0.40000000000000002)),
 ((0.5, 0.5, 0.5),(2.0, 1.9999999999999998, 2.0)),
]),
# cell, ((frac, ortho),...)
  ((4.1,5.2,6.3,100,80,130,), [
 ((0.10000000000000001, 0.0, 0.0),(0.40999999999999998, 0.0, 0.0)),
 ((0.0, 0.10000000000000001, 0.0),(-0.33424955703700043, 0.39834311042186865, 0.0)),
 ((0.0, 0.0, 0.10000000000000001),(0.10939835193016617, -0.051013289294572106, 0.6183281045774256)),
 ((0.10000000000000001, 0.20000000000000001, 0.29999999999999999),(0.069695941716497567, 0.64364635296002093, 1.8549843137322766)),
 ((0.20000000000000001, 0.29999999999999999, 0.10000000000000001),(-0.073350319180835066, 1.1440160419710339, 0.6183281045774256)),
 ((0.29999999999999999, 0.20000000000000001, 0.10000000000000001),(0.67089923785616512, 0.74567293154916525, 0.6183281045774256)),
 ((0.5, 0.5, 0.5),(0.92574397446582857, 1.7366491056364828, 3.0916405228871278)),
]),
# cell, ((frac, ortho),...)
  ((3.5,3.5,6,90,90,120,), [
 ((0.10000000000000001, 0.0, 0.0),(0.35000000000000003, 0.0, 0.0)),
 ((0.0, 0.10000000000000001, 0.0),(-0.17499999999999993, 0.3031088913245536, 0.0)),
 ((0.0, 0.0, 0.10000000000000001),(3.6739403974420595e-17, -3.6739403974420595e-17, 0.60000000000000009)),
 ((0.10000000000000001, 0.20000000000000001, 0.29999999999999999),(2.7675166561703527e-16, 0.60621778264910708, 1.7999999999999998)),
 ((0.20000000000000001, 0.29999999999999999, 0.10000000000000001),(0.17500000000000041, 0.90932667397366063, 0.60000000000000009)),
 ((0.29999999999999999, 0.20000000000000001, 0.10000000000000001),(0.70000000000000018, 0.6062177826491072, 0.60000000000000009)),
 ((0.5, 0.5, 0.5),(0.87500000000000067, 1.5155444566227676, 3.0)),
]),
]

def test0():
    if NeedTestData: TestData()
    msg = 'test cell2Gmat, fillgmat, Gmat2cell'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = cell2Gmat(cell)
        assert np.allclose(G,tG),msg
        assert np.allclose(g,tg),msg
        tcell = Gmat2cell(g)
        assert np.allclose(cell,tcell),msg
        tcell = Gmat2cell(G)
        assert np.allclose(tcell,trcell),msg

def test1():
    if NeedTestData: TestData()
    msg = 'test cell2A and A2Gmat'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = A2Gmat(cell2A(cell))
        assert np.allclose(G,tG),msg
        assert np.allclose(g,tg),msg

def test2():
    if NeedTestData: TestData()
    msg = 'test Gmat2A, A2cell, A2Gmat, Gmat2cell'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = cell2Gmat(cell)
        tcell = A2cell(Gmat2A(G))
        assert np.allclose(cell,tcell),msg

def test3():
    if NeedTestData: TestData()
    msg = 'test invcell2Gmat'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = invcell2Gmat(trcell)
        assert np.allclose(G,tG),msg
        assert np.allclose(g,tg),msg

def test4():
    if NeedTestData: TestData()
    msg = 'test calc_rVsq, calc_rV, calc_V'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        assert np.allclose(calc_rV(cell2A(cell)),trV), msg
        assert np.allclose(calc_V(cell2A(cell)),tV), msg

def test5():
    if NeedTestData: TestData()
    msg = 'test A2invcell'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        rcell = A2invcell(cell2A(cell))
        assert np.allclose(rcell,trcell),msg

def test6():
    if NeedTestData: TestData()
    msg = 'test cell2AB'
    for (cell,coordlist) in CoordTestData:
        A,B = cell2AB(cell)
        for (frac,ortho) in coordlist:
            to = np.inner(A,frac)
            tf = np.inner(B,to)
            assert np.allclose(ortho,to), msg
            assert np.allclose(frac,tf), msg
            to = np.sum(A*frac,axis=1)
            tf = np.sum(B*to,axis=1)
            assert np.allclose(ortho,to), msg
            assert np.allclose(frac,tf), msg

# test GetBraviasNum(...) and GenHBravais(...)
def test7():
    import os.path
    import sys
    import GSASIIspc as spc
    testdir = os.path.join(os.path.split(os.path.abspath( __file__ ))[0],'testinp')
    if os.path.exists(testdir):
        if testdir not in sys.path: sys.path.insert(0,testdir)
    import sgtbxlattinp
    derror = 1e-4
    def indexmatch(hklin, hkllist, system):
        for hklref in hkllist:
            hklref = list(hklref)
            # these permutations are far from complete, but are sufficient to 
            # allow the test to complete
            if system == 'cubic':
                permlist = [(1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1),]
            elif system == 'monoclinic':
                permlist = [(1,2,3),(-1,2,-3)]
            else:
                permlist = [(1,2,3)]

            for perm in permlist:
                hkl = [abs(i) * hklin[abs(i)-1] / i for i in perm]
                if hkl == hklref: return True
                if [-i for i in hkl] == hklref: return True
        else:
            return False

    for key in sgtbxlattinp.sgtbx7:
        spdict = spc.SpcGroup(key)
        cell = sgtbxlattinp.sgtbx7[key][0]
        system = spdict[1]['SGSys']
        center = spdict[1]['SGLatt']

        bravcode = GetBraviasNum(center, system)

        g2list = GenHBravais(sgtbxlattinp.dmin, bravcode, cell2A(cell))

        assert len(sgtbxlattinp.sgtbx7[key][1]) == len(g2list), 'Reflection lists differ for %s' % key
        for h,k,l,d,num in g2list:
            for hkllist,dref in sgtbxlattinp.sgtbx7[key][1]: 
                if abs(d-dref) < derror:
                    if indexmatch((h,k,l,), hkllist, system):
                        break
            else:
                assert 0,'No match for %s at %s (%s)' % ((h,k,l),d,key)

def test8():
    import GSASIIspc as spc
    import sgtbxlattinp
    derror = 1e-4
    dmin = sgtbxlattinp.dmin

    def indexmatch(hklin, hklref, system, axis):
        # these permutations are far from complete, but are sufficient to 
        # allow the test to complete
        if system == 'cubic':
            permlist = [(1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1),]
        elif system == 'monoclinic' and axis=='b':
            permlist = [(1,2,3),(-1,2,-3)]
        elif system == 'monoclinic' and axis=='a':
            permlist = [(1,2,3),(1,-2,-3)]
        elif system == 'monoclinic' and axis=='c':
            permlist = [(1,2,3),(-1,-2,3)]
        elif system == 'trigonal':
            permlist = [(1,2,3),(2,1,3),(-1,-2,3),(-2,-1,3)]
        elif system == 'rhombohedral':
            permlist = [(1,2,3),(2,3,1),(3,1,2),(-1,-2,-3),(-2,-3,-1),(-3,-1,-2)]
        else:
            permlist = [(1,2,3)]

        hklref = list(hklref)
        for perm in permlist:
            hkl = [abs(i) * hklin[abs(i)-1] / i for i in perm]
            if hkl == hklref: return True
            if [-i for i in hkl] == hklref: return True
        return False

    for key in sgtbxlattinp.sgtbx8:
        spdict = spc.SpcGroup(key)[1]
        cell = sgtbxlattinp.sgtbx8[key][0]
        center = spdict['SGLatt']
        Laue = spdict['SGLaue']
        Axis = spdict['SGUniq']
        system = spdict['SGSys']

        g2list = GenHLaue(dmin,Laue,center,Axis,cell2A(cell))
        #if len(g2list) != len(sgtbxlattinp.sgtbx8[key][1]):
        #    print 'failed',key,':' ,len(g2list),'vs',len(sgtbxlattinp.sgtbx8[key][1])
        #    print 'GSAS-II:'
        #    for h,k,l,d in g2list: print '  ',(h,k,l),d
        #    print 'SGTBX:'
        #    for hkllist,dref in sgtbxlattinp.sgtbx8[key][1]: print '  ',hkllist,dref
        assert len(g2list) == len(sgtbxlattinp.sgtbx8[key][1]), (
            'Reflection lists differ for %s' % key
            )
        #match = True
        for h,k,l,d in g2list:
            for hkllist,dref in sgtbxlattinp.sgtbx8[key][1]: 
                if abs(d-dref) < derror:
                    if indexmatch((h,k,l,), hkllist, system, Axis): break
            else:
                assert 0,'No match for %s at %s (%s)' % ((h,k,l),d,key)
                #match = False
        #if not match: 
            #for hkllist,dref in sgtbxlattinp.sgtbx8[key][1]: print '  ',hkllist,dref
            #print center, Laue, Axis, system

if __name__ == '__main__':
    test0()
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
    print "OK"
