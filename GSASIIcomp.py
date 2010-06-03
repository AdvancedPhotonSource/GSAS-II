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
import GSASIIlattice as G2lat

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

def ValEsd(value,esd=0,nTZ=False):                  #NOT complete - don't use
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

        
#GSASII peak fitting routine: Thompson, Cox & Hastings; Finger, Cox & Jephcoat model        

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
    
           
        
