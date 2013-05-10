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
    global values
    values = cPickle.load(file)
    global HistoPhases
    HistoPhases = cPickle.load(file)
    global parmDict
    parmDict = cPickle.load(file)
    global varylist
    varylist = cPickle.load(file)
    global calcControls
    calcControls = cPickle.load(file)
    global pawleyLookup
    pawleyLookup = cPickle.load(file)
    file.close()

    global NeedTestData
    NeedTestData = False
    
def test1():
    if NeedTestData: TestData()
    fplot = plotter.add('function test').gca()
    M = G2st.errRefine(values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup,None)
    fplot.plot(M,'r',label='M')
    fplot.legend(loc='best')

def test2(name,delt):
    if NeedTestData: TestData()
    hplot = plotter.add('derivatives test for '+name).gca()
    dMdV = G2st.dervRefine(values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup,None)
    hplot.plot(dMdV[varylist.index(name)],'b',label='analytic deriv')
    if name in varylist:
        values[varylist.index(name)] -= delt
        M0 = G2st.errRefine(values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup,None)
        values[varylist.index(name)] += 2.*delt
        M1 = G2st.errRefine(values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup,None)
        values[varylist.index(name)] -= delt    
        Mn = (M1-M0)/(2.*delt)
        hplot.plot(Mn,'r+',label='numeric deriv')
        hplot.plot(dMdV[varylist.index(name)]-Mn,'g',label='diff')
    hplot.legend(loc='best')
    
if __name__ == '__main__':
    if NeedTestData: TestData()
    import GSASIItestplot as plot
#    parmDict['0:0:Scale'] = 1.0
    
    global plotter
    plotter = plot.PlotNotebook()
    test1()
    for name in varylist:
        print ' %15s    %f'%(name,parmDict[name])
    names = [
    
#FAP derivatives
        ['0::RBVOa:0:0',0.0001],       #almost perfect
        ['0::RBVPx:0:0',0.000001],     #perfect
        ['0::RBVPy:0:0',0.000001],     #perfect
        ['::RBV;0:0',0.00001],     #perfect
        ['0::RBVT11:0:0',0.0001],   #not good
        ['0::RBVT22:0:0',0.0001],
        ['0::RBVT33:0:0',0.0001],
        ['0::RBVT12:0:0',0.0001],
        ['0::RBVL11:0:0',0.01],     #not good
        ['0::RBVL22:0:0',0.01],
        ['0::RBVL33:0:0',0.01],
        ['0::RBVL12:0:0',0.01],
    
# quinuclidine derivatives
#        ['0::RBVOi:0:0',0.0001],  #bad
#        ['0::RBVOj:0:0',0.0001],  #bad
#        ['0::RBVOk:0:0',0.0001],  #bad
#        ['0::RBVPx:0:0',0.000001],  #perfect
#        ['0::RBVPz:0:0',0.000001],  #perfect
#        ['0::RBVT11:0:0',0.0001],   #not good
#        ['0::RBVT22:0:0',0.0001],
#        ['0::RBVT33:0:0',0.0001],
#        ['0::RBVT13:0:0',0.0001],
#        ['0::RBVL11:0:0',0.01],     #not good
#        ['0::RBVL22:0:0',0.01],
#        ['0::RBVL33:0:0',0.01],
#        ['0::RBVL13:0:0',0.01],
#        ['0::RBVS12:0:0',0.001],        #very close
#        ['0::RBVS21:0:0',0.001],
#        ['0::RBVS23:0:0',0.001],
#        ['0::RBVS32:0:0',0.001],
#        ['::RBV;0:0',0.00001],  #wrong
#        ['::RBV;0:2',0.00001],
#        ['::RBV;1:1',0.00001],
#        ['::RBV;1:2',0.00001],
#        ['::RBV;2:1',0.00001],
#        ['::RBV;2:2',0.00001],
#        ['::RBV;3:3',0.00001],
        
    
# Na BenzRB derivatives
#        ['0::RBROi:0:0',0.0001],       #almost perfect
#        ['0::RBROj:0:0',0.0001],       #almost perfect
#        ['0::RBROk:0:0',0.0001],       #almost perfect
#        ['0::RBRPx:0:0',0.000001],     #perfect
#        ['0::RBRPy:0:0',0.000001],     #perfect
#        ['0::RBRPz:0:0',0.000001],     #perfect
#        ['0::AUiso:0',0.0001],         #perfect
#        [':0:BkPkint:2',0.1],          #perfect
#        ['0:0:Mustrain:6',0.00001],    #perfect
#        ['0:0:Size;a',0.0001],         #perfect
#        ['0::dAz:0',0.0001],           #perfect
#        [':0:DisplaceX',0.001],        #OK but flaky
#        ['0::dAx:0',1.e-7],            #perfect
#        ['0::RBRTr;0:0:0',0.001],      #perfect
#        ['0::RBRU:0:0',0.00001],       #perfect
# Na Benz - restraints derivatives
#        ['0::dAz:0',0.0001],           #perfect
#        ['0::dAx:10',0.0001],           #perfect
#        ['0::dAy:10',0.0001],           #perfect
#        ['0::dAz:10',0.0001],           #perfect
#        ['0::dAx:12',0.0001],           #perfect
#        ['0::dAy:12',0.0001],           #perfect
#        ['0::dAz:12',0.0001],           #perfect
        
        ]
    for [name,delt] in names:
        test2(name,delt)
    print "OK"
    plotter.StartEventLoop()
