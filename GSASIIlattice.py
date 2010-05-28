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

def fillgmat(cell):
    '''Compute lattice metric tensor from unit cell constants
    cell is tuple with a,b,c,alpha, beta, gamma (degrees)
    returns 3x3 numpy array
    '''
    a,b,c,alp,bet,gam = cell
    g = np.array([[a*a,a*b*cosd(gam),a*c*cosd(bet)],[a*b*cosd(gam),b*b,b*c*cosd(alp)],
        [a*c*cosd(bet),b*c*cosd(alp),c*c]])
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
    '''Compute reciprocal metric tensor (G) from A tensor
    returns reciprocal (G) & real (g) metric tensors (list of two 3x3 arrays)
    '''
    G = np.zeros(shape=(3,3))
    G = [[A[0],A[3]/2.,A[4]/2.], [A[3]/2.,A[1],A[5]/2.], [A[4]/2.,A[5]/2.,A[2]]]
    g = nl.inv(G)
    return G,g

def Gmat2A(G):
    'Compute A tensor from reciprocal metric tensor (G)'
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
    'Compute the square of the reciprocal lattice volume (V* **2) from A tensor'
    rVsq = A[0]*A[1]*A[2]+0.25*(A[3]*A[4]*A[5]-A[0]*A[5]**2-A[1]*A[4]**2-A[2]*A[3]**2)
    if rVsq < 0:
        return 1
    return rVsq
    
def calc_rV(A):
    'Compute the reciprocal lattice volume (V*) from A tensor'
    return np.sqrt(calc_rVsq(A))
    
def calc_V(A):
    'Compute the real lattice volume (V) from A tensor'
    return 1./calc_rV(A)

def A2invcell(A):
    '''Compute reciprocal unit cell constants from A tensor
    returns tuple with a*,b*,c*,alpha*, beta*, gamma* (degrees)
    '''
    G,g = A2Gmat(A)
    return Gmat2cell(G)
    # Code below is broken
    ainv = np.sqrt(max(0.,A[0]))
    binv = np.sqrt(max(0.,A[1]))
    cinv = np.sqrt(max(0.,A[2]))
    gaminv = acosd(max(-0.5,min(0.5,0.5*A[3]/(ainv*binv))))
    betinv = acosd(max(-0.5,min(0.5,0.5*A[4]/(ainv*cinv))))
    alpinv = acosd(max(-0.5,min(0.5,0.5*A[5]/(binv*cinv))))
    return ainv,binv,cinv,alpinv,betinv,gaminv

# not tested yet (broken?)
def cell2AB(cell):
    '''Computes orthogonalization matrix from unit cell constants
    cell is tuple with a,b,c,alpha, beta, gamma (degrees)
    returns list of two 3x3 numpy arrays
       A for Cartesian to crystal transformations A*X = x 
       B (inverse) for crystal to Cartesian transformation B*x = X
    '''
    G,g = cell2Gmat(cell)       #reciprocal & real metric tensors
    cosAlpStar = G[2][1]/np.sqrt(G[1][1]*G[2][2])
    sinAlpStar = np.sqrt(1.0-cosAlpStar**2)
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


# output from uctbx computed on platform darwin on 2010-05-28
array = np.array
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

def test0():
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
    msg = 'test cell2A and A2Gmat'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = A2Gmat(cell2A(cell))
        assert np.allclose(G,tG),msg
        assert np.allclose(g,tg),msg

def test2():
    msg = 'test Gmat2A, A2cell, A2Gmat, Gmat2cell'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = cell2Gmat(cell)
        tcell = A2cell(Gmat2A(G))
        assert np.allclose(cell,tcell),msg

def test3():
    msg = 'test invcell2Gmat'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        G, g = invcell2Gmat(trcell)
        assert np.allclose(G,tG),msg
        assert np.allclose(g,tg),msg

def test4():
    msg = 'test calc_rVsq, calc_rV, calc_V'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        assert np.allclose(calc_rV(cell2A(cell)),trV), msg
        assert np.allclose(calc_V(cell2A(cell)),tV), msg

def test5():
    msg = 'test A2invcell'
    for (cell, tg, tG, trcell, tV, trV) in CellTestData:
        rcell = A2invcell(cell2A(cell))
        assert np.allclose(rcell,trcell),msg

if __name__ == '__main__':
    test0()
    test1()
    test2()
    test3()
    test4()
    test5()
    print "OK"
