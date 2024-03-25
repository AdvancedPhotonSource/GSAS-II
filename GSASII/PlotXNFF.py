# -*- coding: utf-8 -*-
#PlotXNFF.py
'''
*PlotXNFF: Check x-ray & neutron form factor computation*
=========================================

Use this to check form factors used in x-ray & neutron scattering

'''

import sys
import wx
# the next line removes the need for pythonw. Thanks to Matt Newville!
# appears unneaded from wx 4.2.1 on
if sys.platform.lower() == 'darwin': wx.PyApp.IsDisplayAvailable = lambda _: True
import numpy as np
import GSASIIpath
GSASIIpath.SetBinaryPath()
import GSASIItestplot as plot
import GSASIIElem as G2el
import GSASIIElemGUI as G2elG
import atmdata
WACV = wx.ALIGN_CENTER_VERTICAL

try:
    wx.NewIdRef
    wx.NewId = wx.NewIdRef
except AttributeError:
    pass

[wxID_FILEEXIT, wxID_PICKELEM,
] = [wx.NewId() for _init_coll_File_Items in range(2)]

class PlotXNFF(wx.Frame):
    
    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='PlotXNFF', parent=parent,
            size=wx.Size(460, 250),style=wx.DEFAULT_FRAME_STYLE, title='Plot Xray, Neutron & Magnetic Form Factors')
        self.PlotFFMenu = wx.MenuBar()
        self.File = wx.Menu(title='')
        self.File.Append(wxID_PICKELEM,'Pick element','Pick element')
        self.File.Append(wxID_FILEEXIT,'Exit','Exit from PlotXNFF')
        self.Bind(wx.EVT_MENU,self.OnPickElement,id=wxID_PICKELEM)
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=wxID_FILEEXIT)
        self.PlotFFMenu.Append(menu=self.File, title='File')
        self.SetMenuBar(self.PlotFFMenu)
        self.PlotFFPanel = wx.ScrolledWindow(self)        
        self.plotNB = plot.PlotNotebook()
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)

    def ExitMain(self, event):
        sys.exit()
        
    def OnFileExit(self,event):
        self.Close()

    def OnPickElement(self,Parms):

        PE = G2elG.PickElement(self.PlotFFPanel,oneOnly=True)
        if PE.ShowModal() == wx.ID_OK:
            self.xrayFFs = G2el.GetFormFactorCoeff(PE.Elem)
            self.elecFFs = G2el.GetEFormFactorCoeff(PE.Elem)
            self.Elem = PE.Elem
            self.atmInfo = G2el.GetAtomInfo(PE.Elem)
            self.magFFs = G2el.GetMagFormFacCoeff(PE.Elem)
        PE.Destroy()
        self.MakePlot()

    def MakePlot(self):
        
        def test1():
            fplot = self.plotNB.add('Neutron scattering lengths').gca()
            lams = np.linspace(0.3,10.0,1000,True)
            
            BLtable = {}
            El = self.Elem.capitalize()
            for isotope in self.atmInfo['Isotopes']:
                if 'Nat' in isotope:
                    BLtable[isotope] = ['',atmdata.AtmBlens[El+'_']]
                else:
                    BLtable[isotope] = ['',atmdata.AtmBlens[El+'_'+isotope]]
            for isotope in BLtable:
                if 'BW-LS' in BLtable[isotope][1]:
                    b0 = np.ones_like(lams)*BLtable[isotope][1]['BW-LS'][0]                    
                else:
                    b0 = np.ones_like(lams)*BLtable[isotope][1]['SL'][0]
                bp,bpp = G2el.BlenResTOF([isotope,],BLtable,lams)
                fplot.plot(lams,b0,label=isotope+El+' b0')
                fplot.plot(lams,bp[0]-b0,label=isotope+El+" b'")
                fplot.plot(lams,bpp[0],label=isotope+El+' b"')
                
            fplot.legend(loc='best')
            fplot.set_xlabel('wavelength, A')
            fplot.set_ylabel('b')

        def test2():
            fplot = self.plotNB.add('X-ray form factors').gca()
            sq = np.linspace(0,2.,1000,True)
            for El in self.xrayFFs:
                F = G2el.ScatFac(El, sq**2)
                fplot.plot(sq,F,label=El['Symbol'])
            fplot.legend(loc='best')
            fplot.set_xlabel('sin-theta/lambda')
            fplot.set_ylabel('form factor')
            
        def test3():
            fplot = self.plotNB.add('Electron form factor').gca()
            sq = np.linspace(0,2.,1000,True)
            for El in self.elecFFs:
                F = G2el.ScatFac(El, sq**2)
                fplot.plot(sq,F,label=El['Symbol'])
            fplot.legend(loc='best')
            fplot.set_xlabel('sin-theta/lambda')
            fplot.set_ylabel('form factor')
            
        def test4():
            fplot = self.plotNB.add('magnetic form factors').gca()
            sq = np.linspace(0,1.,1000,True)
            for El in self.magFFs:
                El['gfac'] = 2.0
                F = G2el.MagScatFac(El, sq**2)
                fplot.plot(sq,F,label=El['Symbol'])
            fplot.legend(loc='best')
            fplot.set_xlabel('sin-theta/lambda')
            fplot.set_ylabel('form factor')
            
        while self.plotNB.nb.GetPageCount():
            self.plotNB.nb.DeletePage(0)
        test1()
        test2()
        test3()
        if len(self.magFFs):
            test4()
        
        self.plotNB.Show()
        
if __name__ == "__main__":
    import GSASIIplot as G2plt
    app = wx.App()
    GSASIIpath.InvokeDebugOpts()
    frm = wx.Frame(None) # create a frame
    frm.Show(False)
    win = PlotXNFF(frm)
    win.Show()
    win.Bind(wx.EVT_WINDOW_DESTROY,lambda event: sys.exit())
    app.MainLoop()
