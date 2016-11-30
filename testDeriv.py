# -*- coding: utf-8 -*-
#testDeriv.py
'''
*testDeriv: Check derivative computation*
=========================================

Use this to check derivatives used in structure least squares
refinement against numerical values computed in this script.

To use set ``DEBUG=True`` in GSASIIstrMain.py (line 40, as of version
2546); run the least squares - zero cycles is sufficient.  Do the "Save
Results"; this will write the file testDeriv.dat in the local
directory.

Then run this program to see plots of derivatives for all
parameters refined in the last least squares.  Shown will be numerical
derivatives generated over all observations (including penalty terms)
and the corresponding analytical ones produced in the least
squares. They should match. Profiling is also done for function 
calculation & for the 1st selected derivative (rest should be the same).
'''

import sys
import cPickle
import cProfile,pstats,StringIO
import wx
import numpy as np
import GSASIIpath
import GSASIIstrMath as G2stMth
import GSASIItestplot as plot
import GSASIImapvars as G2mv
import pytexture as ptx
ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics

def create(parent):
    return testDeriv(parent)
    
[wxID_FILEEXIT, wxID_FILEOPEN, wxID_MAKEPLOTS,
] = [wx.NewId() for _init_coll_File_Items in range(3)]

def FileDlgFixExt(dlg,file):            #this is needed to fix a problem in linux wx.FileDialog
    ext = dlg.GetWildcard().split('|')[2*dlg.GetFilterIndex()+1].strip('*')
    if ext not in file:
        file += ext
    return file
    
class testDeriv(wx.Frame):

    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='testDeriv', parent=parent,
            size=wx.Size(800, 250),style=wx.DEFAULT_FRAME_STYLE, title='Test Jacobian Derivatives')
        self.testDerivMenu = wx.MenuBar()
        self.File = wx.Menu(title='')
        self.File.Append(help='Open testDeriv.dat', id=wxID_FILEOPEN,
             kind=wx.ITEM_NORMAL,text='Open testDeriv.dat file')
        self.File.Append(help='Make derivative plots',id=wxID_MAKEPLOTS,
            kind=wx.ITEM_NORMAL,text='Make plots')
        self.File.Append(help='Exit from testDeriv', id=wxID_FILEEXIT, kind=wx.ITEM_NORMAL,
            text='Exit')
        self.Bind(wx.EVT_MENU, self.OnTestRead, id=wxID_FILEOPEN)
        self.Bind(wx.EVT_MENU,self.OnMakePlots,id=wxID_MAKEPLOTS)
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=wxID_FILEEXIT)
        self.testDerivMenu.Append(menu=self.File, title='File')
        self.SetMenuBar(self.testDerivMenu)
        self.testDerivPanel = wx.ScrolledWindow(self)        
        self.plotNB = plot.PlotNotebook()
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)    
        self.dirname = ''
        self.testfile = []
        self.dataFrame = None

    def ExitMain(self, event):
        sys.exit()
        
    def OnFileExit(self,event):
        if self.dataFrame:
            self.dataFrame.Clear() 
            self.dataFrame.Destroy()
        self.Close()

    def OnTestRead(self,event):
        dlg = wx.FileDialog(self, 'Open testDeriv.dat file',defaultFile='testDeriv.dat',
            wildcard='testDeriv.dat')
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.dirname = dlg.GetDirectory()
                testFile = dlg.GetPath()
                file = open(testFile,'rb')
                self.values = cPickle.load(file)
                self.HistoPhases = cPickle.load(file)
                (self.constrDict,self.fixedList,self.depVarList) = cPickle.load(file)
                self.parmDict = cPickle.load(file)
                self.varylist = cPickle.load(file)
                self.calcControls = cPickle.load(file)
                self.pawleyLookup = cPickle.load(file)
                self.use = [False for i in range(len(self.varylist+self.depVarList))]
                self.delt = [max(abs(self.parmDict[name])*0.0001,1e-6) for name in self.varylist+self.depVarList]
                file.close()
                groups,parmlist = G2mv.GroupConstraints(self.constrDict)
                G2mv.GenerateConstraints(groups,parmlist,self.varylist,self.constrDict,self.fixedList,self.parmDict)
                self.UpdateControls(event)
                print G2mv.VarRemapShow(self.varylist)
                print 'Dependent Vary List:',self.depVarList
        finally:
            dlg.Destroy()
            
    def UpdateControls(self,event):
        
        def OnItemCk(event):
            Obj = event.GetEventObject()
            item = ObjInd[Obj.GetId()]
            self.use[item] = Obj.GetValue()
            
        def OnDelValue(event):
            event.Skip()
            Obj = event.GetEventObject()
            item = ObjInd[Obj.GetId()]
            try:
                value = float(Obj.GetValue())
            except ValueError:
                value = self.delt[item]
            self.delt[item] = value
            Obj.SetValue('%g'%(value))
        
        self.testDerivPanel.DestroyChildren()
        ObjInd = {}
        varylist = self.varylist
        depVarList = self.depVarList
        use = self.use
        delt = self.delt
        mainSizer = wx.FlexGridSizer(0,8,5,5)
        for id,[ck,name,d] in enumerate(zip(use,varylist+depVarList,delt)):
            useVal = wx.CheckBox(self.testDerivPanel,label=name)
            useVal.SetValue(ck)
            ObjInd[useVal.GetId()] = id
            useVal.Bind(wx.EVT_CHECKBOX, OnItemCk)
            mainSizer.Add(useVal,0)
            delVal = wx.TextCtrl(self.testDerivPanel,wx.ID_ANY,'%g'%(d),style=wx.TE_PROCESS_ENTER)
            ObjInd[delVal.GetId()] = id
            delVal.Bind(wx.EVT_TEXT_ENTER,OnDelValue)
            delVal.Bind(wx.EVT_KILL_FOCUS,OnDelValue)
            mainSizer.Add(delVal,0)
        self.testDerivPanel.SetSizer(mainSizer)    
        Size = mainSizer.GetMinSize()
        self.testDerivPanel.SetScrollbars(10,10,Size[0]/10-4,Size[1]/10-1)
        Size[1] = min(200,Size[1])
        self.testDerivPanel.SetSize(Size)

    def OnMakePlots(self,event):
        
        def test1():
            fplot = self.plotNB.add('function test').gca()
            pr = cProfile.Profile()
            pr.enable()
            M = G2stMth.errRefine(self.values,self.HistoPhases,
                self.parmDict,self.varylist,self.calcControls,
                self.pawleyLookup,None)
            pr.disable()
            s = StringIO.StringIO()
            sortby = 'tottime'
            ps = pstats.Stats(pr, stream=s).strip_dirs().sort_stats(sortby)
            print 'Profiler of function calculation; top 50% of routines:'
            ps.print_stats("GSASII",.5)
            print s.getvalue()
            fplot.plot(M,'r',label='M')
            fplot.legend(loc='best')
            
        def test2(name,delt,doProfile):
            Title = 'derivatives test for '+name
            varyList = self.varylist+self.depVarList
            hplot = self.plotNB.add(Title).gca()
            if doProfile:
                pr = cProfile.Profile()
                pr.enable()
            dMdV = G2stMth.dervRefine(self.values,self.HistoPhases,self.parmDict,
                varyList,self.calcControls,self.pawleyLookup,None)
            if doProfile:
                pr.disable()
                s = StringIO.StringIO()
                sortby = 'tottime'
                ps = pstats.Stats(pr, stream=s).strip_dirs().sort_stats(sortby)
                ps.print_stats("GSASII",.5)
                print 'Profiler of '+name+' derivative calculation; top 50% of routines:'
                print s.getvalue()
            M2 = dMdV[varyList.index(name)]
            hplot.plot(M2,'b',label='analytic deriv')
            mmin = np.min(dMdV[varyList.index(name)])
            mmax = np.max(dMdV[varyList.index(name)])
            print 'parameter:',name,self.parmDict[name],delt,mmin,mmax
            if name in self.varylist:
                self.values[self.varylist.index(name)] -= delt
                M0 = G2stMth.errRefine(self.values,self.HistoPhases,self.parmDict,
                    varyList,self.calcControls,self.pawleyLookup,None)
                self.values[self.varylist.index(name)] += 2.*delt
                M1 = G2stMth.errRefine(self.values,self.HistoPhases,self.parmDict,
                    varyList,self.calcControls,self.pawleyLookup,None)
                self.values[self.varylist.index(name)] -= delt
            elif name in self.depVarList:   #in depVarList
                if 'dA' in name:
                    name = name.replace('dA','A')
                    delt *= -1
                self.parmDict[name] -= delt
                M0 = G2stMth.errRefine(self.values,self.HistoPhases,self.parmDict,
                    varyList,self.calcControls,self.pawleyLookup,None)
                self.parmDict[name] += 2.*delt
                M1 = G2stMth.errRefine(self.values,self.HistoPhases,self.parmDict,
                    varyList,self.calcControls,self.pawleyLookup,None)
                self.parmDict[name] -= delt    
            Mn = (M1-M0)/(2.*abs(delt))
            hplot.plot(Mn,'r',label='numeric deriv')
            hplot.plot(M2-Mn,'g',label='diff')
#            GSASIIpath.IPyBreak()
            hplot.legend(loc='best')            
            
        while self.plotNB.nb.GetPageCount():
            self.plotNB.nb.DeletePage(0)
            
        test1()

        doProfile = True
        for use,name,delt in zip(self.use,self.varylist+self.depVarList,self.delt):
            if use:
                test2(name,delt,doProfile)
                doProfile = False
        
        self.plotNB.Show()
        
class testDerivmain(wx.App):
    def OnInit(self):
        self.main = testDeriv(None)
        self.main.Show()
        self.SetTopWindow(self.main)
        return True

def main():
    'Starts main application to compute and plot derivatives'
    application = testDerivmain(0)
    application.MainLoop()
    
if __name__ == '__main__':
    GSASIIpath.InvokeDebugOpts()
    main()
