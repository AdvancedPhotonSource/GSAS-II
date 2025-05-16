# -*- coding: utf-8 -*-
#testDeriv.py
'''
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
import os
import copy
import pickle
import io as StringIO
import cProfile,pstats
import wx
import numpy as np
# path hack for restart, when needed
import importlib.util
try:
    importlib.util.find_spec('GSASII.GSASIIGUI')
except ModuleNotFoundError:
    print('GSAS-II not installed in Python; Hacking sys.path')
    sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))
from GSASII import GSASIIpath
GSASIIpath.SetBinaryPath()
from GSASII import GSASIIstrMath as G2stMth
from GSASII import GSASIItestplot as plot
from GSASII import GSASIImapvars as G2mv
try:  # fails on doc build
    from . import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics
except ImportError:
    pass

try:
    NewId = wx.NewIdRef
except AttributeError:
    NewId = wx.NewId
[wxID_FILEEXIT, wxID_FILEOPEN, wxID_MAKEPLOTS, wxID_CLEARSEL,wxID_SELECTALL,
] = [NewId() for _init_coll_File_Items in range(5)]

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
        self.File.Append(wxID_FILEOPEN,'Open testDeriv file\tCtrl+O','Open testDeriv')
        self.File.Append(wxID_MAKEPLOTS,'Make plots\tCtrl+P','Make derivative plots')
        self.File.Append(wxID_SELECTALL,'Select all\tCtrl+S')
        self.File.Append(wxID_CLEARSEL,'Clear selections\tCtrl+C')
        self.File.Append(wxID_FILEEXIT,'Exit\tALT+F4','Exit from testDeriv')
        self.Bind(wx.EVT_MENU,self.OnTestRead, id=wxID_FILEOPEN)
        self.Bind(wx.EVT_MENU,self.OnMakePlots,id=wxID_MAKEPLOTS)
        self.Bind(wx.EVT_MENU,self.ClearSelect,id=wxID_CLEARSEL)
        self.Bind(wx.EVT_MENU,self.SelectAll,id=wxID_SELECTALL)
        self.Bind(wx.EVT_MENU,self.OnFileExit, id=wxID_FILEEXIT)
        self.testDerivMenu.Append(menu=self.File, title='File')
        self.SetMenuBar(self.testDerivMenu)
        self.testDerivPanel = wx.ScrolledWindow(self)
        self.plotNB = plot.PlotNotebook()
        self.testFile = ''
        arg = sys.argv
        if len(arg) > 1 and arg[1]:
            try:
                self.testFile = os.path.splitext(arg[1])[0]+u'.testDeriv'
            except:
                self.testFile = os.path.splitext(arg[1])[0]+'.testDeriv'
            self.TestRead()
            self.UpdateControls(None)
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)    
        self.dirname = ''
        self.testfile = []
        self.dataFrame = None
        self.timingOn = False

    def ExitMain(self, event):
        sys.exit()
        
    def OnFileExit(self,event):
        if self.dataFrame:
            self.dataFrame.Clear() 
            self.dataFrame.Destroy()
        self.Close()
        
    def SelectAll(self,event):
        self.use = [True for name in self.names]
        for i,name in enumerate(self.names):
            if 'Back' in name:
                self.use[i] = False
        self.UpdateControls(event)
        
    def ClearSelect(self,event):
        self.use = [False for i in range(len(self.names))]
        self.UpdateControls(event)

    def OnTestRead(self,event):
        dlg = wx.FileDialog(self, 'Open *.testDeriv file',defaultFile='*.testDeriv',
            wildcard='*.testDeriv')
        if self.dirname:
            dlg.SetDirectory(self.dirname)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.dirname = dlg.GetDirectory()
                self.testFile = dlg.GetPath()
                self.TestRead()
                self.UpdateControls(event)
        finally:
            dlg.Destroy()
            
    def TestRead(self):
        file = open(self.testFile,'rb')
        self.values = pickle.load(file,encoding='Latin-1')
        self.HistoPhases = pickle.load(file,encoding='Latin-1')
        (self.constrDict,self.fixedList,self.depVarList) = pickle.load(file,encoding='Latin-1')
        self.parmDict = pickle.load(file,encoding='Latin-1')
        self.varylist = pickle.load(file,encoding='Latin-1')
        self.calcControls = pickle.load(file,encoding='Latin-1')
        self.pawleyLookup = pickle.load(file,encoding='Latin-1')
        self.names = self.varylist+self.depVarList
        self.use = [False for i in range(len(self.names))]
        self.delt = [max(abs(self.parmDict[name])*0.0001,1e-6) for name in self.names]
        for iname,name in enumerate(self.names):
            if name.split(':')[-1] in ['Shift','DisplaceX','DisplaceY',]:
                self.delt[iname] = 0.1
        file.close()
        G2mv.InitVars()
        msg = G2mv.EvaluateMultipliers(self.constrDict,self.parmDict)
        if msg:
            print('Unable to interpret multiplier(s): '+msg)
            raise Exception
        G2mv.GenerateConstraints(self.varylist,self.constrDict,self.fixedList,self.parmDict)
        print(G2mv.VarRemapShow(self.varylist))
        print('Dependent Vary List:',self.depVarList)
        G2mv.Map2Dict(self.parmDict,copy.copy(self.varylist))   # compute independent params, N.B. changes varylist
        G2mv.Dict2Map(self.parmDict) # imposes constraints on dependent values

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
        
        if self.testDerivPanel.GetSizer():
            self.testDerivPanel.GetSizer().Clear(True)
        ObjInd = {}
        use = self.use
        delt = self.delt
        topSizer = wx.BoxSizer(wx.VERTICAL)
        self.timingVal = wx.CheckBox(self.testDerivPanel,label='Show Execution Profiling')
        topSizer.Add(self.timingVal,0)
        topSizer.Add((-1,10))
        mainSizer = wx.FlexGridSizer(0,8,5,5)
        for id,[ck,name,d] in enumerate(zip(use,self.names,delt)):
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
        topSizer.Add(mainSizer,0)
        self.testDerivPanel.SetSizer(topSizer)    
        Size = topSizer.GetMinSize()
        self.testDerivPanel.SetScrollbars(10,10,int(Size[0]/10-4),int(Size[1]/10-1))
        Size[1] = max(200,Size[1])
        Size[0] += 20
        self.SetSize(Size)

    def OnMakePlots(self,event):
        
        def test1():
            fplot = self.plotNB.add('function test').gca()
            if self.timingOn:
                pr = cProfile.Profile()
                pr.enable()
            M = G2stMth.errRefine(self.values,self.HistoPhases,
                self.parmDict,self.varylist,self.calcControls,
                self.pawleyLookup,None)
            if self.timingOn:
                pr.disable()
                s = StringIO.StringIO()
                sortby = 'tottime'
                ps = pstats.Stats(pr, stream=s).strip_dirs().sort_stats(sortby)
                print('Profiler of function calculation; top 50% of routines:')
                ps.print_stats("GSASII",.5)
                print(s.getvalue())
            fplot.plot(M,'r',label='M')
            fplot.legend(loc='best')
            
        def test2(name,delt,doProfile):
            Title = 'derivatives test for '+name
            ind = self.names.index(name)
            hplot = self.plotNB.add(Title).gca()
            if doProfile and self.timingOn:
                pr = cProfile.Profile()
                pr.enable()
            #regenerate minimization fxn
            G2stMth.errRefine(self.values,self.HistoPhases,
                self.parmDict,self.varylist,self.calcControls,
                self.pawleyLookup,None)
            dMdV = G2stMth.dervRefine(self.values,self.HistoPhases,self.parmDict,
                self.names,self.calcControls,self.pawleyLookup,None)
            if doProfile and self.timingOn:
                pr.disable()
                s = StringIO.StringIO()
                sortby = 'tottime'
                ps = pstats.Stats(pr, stream=s).strip_dirs().sort_stats(sortby)
                ps.print_stats("GSASII",.5)
                print('Profiler of '+name+' derivative calculation; top 50% of routines:')
                print(s.getvalue())
            M2 = dMdV[ind]
            hplot.plot(M2,'b',label='analytic deriv')
            mmin = np.min(dMdV[ind])
            mmax = np.max(dMdV[ind])
            if name in self.varylist:
                ind = self.varylist.index(name)
                orig = copy.copy(self.parmDict)  # save parmDict before changes
                self.parmDict[name] = self.values[ind] =  self.values[ind] - delt
                G2mv.Dict2Map(self.parmDict)
                first = True
                for i in self.parmDict:
                    if 'UVmat' in i:
                        continue
                    if orig[i] != self.parmDict[i] and i != name:
                        if first:
                            print('Propagated changes from this shift')
                            print(name,orig[name],self.parmDict[name],orig[name]-self.parmDict[name])
                            print('are:')
                            first = False
                        print(i,orig[i],self.parmDict[i],orig[i]-self.parmDict[i])
                M0 = G2stMth.errRefine(self.values,self.HistoPhases,self.parmDict,
                    self.names,self.calcControls,self.pawleyLookup,None)
                self.parmDict[name] = self.values[ind] =  self.values[ind] + 2.*delt
                G2mv.Dict2Map(self.parmDict)
                M1 = G2stMth.errRefine(self.values,self.HistoPhases,self.parmDict,
                    self.names,self.calcControls,self.pawleyLookup,None)
                self.parmDict[name] = self.values[ind] =  self.values[ind] - delt
                G2mv.Dict2Map(self.parmDict)
            elif name in self.depVarList:   #in depVarList
                if 'dA' in name:
                    name = name.replace('dA','A')
                    #delt *= -1  # why???
                self.parmDict[name] -= delt
                G2mv.Dict2Map(self.parmDict)
                M0 = G2stMth.errRefine(self.values,self.HistoPhases,self.parmDict,
                        self.names,self.calcControls,self.pawleyLookup,None)
                self.parmDict[name] += 2.*delt
                G2mv.Dict2Map(self.parmDict)
                M1 = G2stMth.errRefine(self.values,self.HistoPhases,self.parmDict,
                        self.names,self.calcControls,self.pawleyLookup,None)
                self.parmDict[name] -= delt    
                G2mv.Dict2Map(self.parmDict)
            Mn = (M1-M0)/(2.*abs(delt))
            print('parameter:',name,self.parmDict[name],delt,mmin,mmax,np.sum(M0),np.sum(M1),np.sum(Mn))
            hplot.plot(Mn,'r',label='numeric deriv')
            hplot.legend(loc='best')            
            
        while self.plotNB.nb.GetPageCount():
            self.plotNB.nb.DeletePage(0)
            
        test1()
        self.timingOn = self.timingVal.GetValue()

        doProfile = True
        for use,name,delt in zip(self.use,self.names,self.delt):
            if use:
                test2(name,delt,doProfile)
                doProfile = False
        
        self.plotNB.Show()
        
def main():
    'Starts main application to compute and plot derivatives'
    application = wx.App(0)
    application.main = testDeriv(None)
    application.main.Show()
    application.SetTopWindow(application.main)
    application.MainLoop()
    
if __name__ == '__main__':
    GSASIIpath.InvokeDebugOpts()
    main()
