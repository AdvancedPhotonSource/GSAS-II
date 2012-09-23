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
    fplot = plotter.add('function test').gca()
    yc,yb = G2st.getPowderProfile(parmDict,xdata,varylist,Histogram,Phases,calcControls,pawleyLookup)
    fplot.plot(xdata,yc+yb,'r')

def test2(name,delt):
    if NeedTestData: TestData()
    varyList = [name,]
    limits = Histogram['Limits'][1]
    data = Histogram['Data']
    xdata = data[0]
    hplot = plotter.add('derivatives test for '+name).gca()
    ya = G2st.getPowderProfileDerv(parmDict,xdata,varyList,Histogram,Phases,calcControls,pawleyLookup)[0]
    hplot.plot(xdata,ya)
    if 'dA' in name:
        name = ''.join(name.split('d'))
        varyList = [name,]
    parmDict[name] -= delt
    y0,yb = G2st.getPowderProfile(parmDict,xdata,varyList,Histogram,Phases,calcControls,pawleyLookup)
    y0 += yb
    parmDict[name] += 2.*delt
    y1,yb = G2st.getPowderProfile(parmDict,xdata,varyList,Histogram,Phases,calcControls,pawleyLookup)
    y1 += yb 
    yn = (y1-y0)/(2.*delt)
    hplot.plot(xdata,yn,'r+')
    hplot.plot(xdata,ya-yn)
    
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
        ['0:0:Size:mx',0.001],
        ['0:0:Mustrain:mx',0.001],
        ['0:0:Size:i',0.001],
        ['0:0:Mustrain:i',0.1],
        ]
    for [name,delt] in names:
        test2(name,delt)
    print "OK"
    plotter.StartEventLoop()
