#testGSASIIstruct.py

import GSASIIpath
import GSASIIstruct as G2st
import numpy as np
import pytexture as ptx
ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics

#testing data
NeedTestData = True
def TestData():
    import cPickle
    file = open('structTestdata.dat','rb')
    global parmDict
    parmDict = cPickle.load(file)
    global varylist
    varylist = cPickle.load(file)
    global Histogram
    Histogram = cPickle.load(file)
    global Phases
    Phases = cPickle.load(file)
    global RBData
    RBData = cPickle.load(file)
    global calcControls
    calcControls = cPickle.load(file)
    global pawleyLookup
    pawleyLookup = cPickle.load(file)
    file.close()

    global NeedTestData
    NeedTestData = False
    
def test1():
    if NeedTestData: TestData()
    limits = Histogram['Limits'][1]
    data = Histogram['Data']
    xdata = data[0]
    xB = np.searchsorted(xdata,limits[0])
    xF = np.searchsorted(xdata,limits[1])
    fplot = plotter.add('function test').gca()
    yc,yb = G2st.getPowderProfile(parmDict,xdata[xB:xF],varylist,Histogram,Phases,calcControls,pawleyLookup)
    fplot.plot(xdata[xB:xF],yc+yb,'r',label='calc+bkg')
    fplot.legend()

def test2(name,delt):
    if NeedTestData: TestData()
    varyList = [name,]
    limits = Histogram['Limits'][1]
    data = Histogram['Data']
    xdata = data[0]
    xB = np.searchsorted(xdata,limits[0])
    xF = np.searchsorted(xdata,limits[1])
    hplot = plotter.add('derivatives test for '+name).gca()
    ya = G2st.getPowderProfileDerv(parmDict,xdata[xB:xF],varyList,Histogram,Phases,RBData,calcControls,pawleyLookup)[0]
    hplot.plot(xdata[xB:xF],ya,'b',label='analytic deriv')
    if 'dA' in name:
        name = ''.join(name.split('d'))
        varyList = [name,]
    parmDict[name] -= delt
    y0,yb = G2st.getPowderProfile(parmDict,xdata[xB:xF],varyList,Histogram,Phases,calcControls,pawleyLookup)
    y0 += yb
    parmDict[name] += 2.*delt
    y1,yb = G2st.getPowderProfile(parmDict,xdata[xB:xF],varyList,Histogram,Phases,calcControls,pawleyLookup)
    y1 += yb 
    yn = (y1-y0)/(2.*delt)
    hplot.plot(xdata[xB:xF],yn,'r+',label='numeric deriv')
    hplot.plot(xdata[xB:xF],ya-yn,'g',label='diff')
    hplot.legend()
    
if __name__ == '__main__':
    if NeedTestData: TestData()
    import GSASIItestplot as plot
#    parmDict['0:0:Scale'] = 1.0
    
    global plotter
    plotter = plot.PlotNotebook()
    test1()
    for name in parmDict:
        print name,parmDict[name]
    names = [
        ['0::AUiso:0',0.001],
        ['::RBV;0:0',0.001],
        ['0::RBVT11:0:0',0.1],
        ['0::RBVL11:0:0',0.1],
        ['0::RBVPz:0:0',0.0001],
        ['0::RBVOa:0:0',0.001],
        ]
    for [name,delt] in names:
        test2(name,delt)
    print "OK"
    plotter.StartEventLoop()
