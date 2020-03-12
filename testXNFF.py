# -*- coding: utf-8 -*-
#testXNFF.py
'''
*testXNFF: Check x-ray & neutron form factor computation*
=========================================

Use this to check form factors used in neutron scattering

'''

import sys
import wx
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

def create(parent):
    return testXNFF(parent)
    
[wxID_FILEEXIT, wxID_PICKELEM,
] = [wx.NewId() for _init_coll_File_Items in range(2)]

class testXNFF(wx.Frame):
    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, name='testXNFF', parent=parent,
            size=wx.Size(460, 250),style=wx.DEFAULT_FRAME_STYLE, title='Test Xray & Neutron Form Factors')
        self.testFFMenu = wx.MenuBar()
        self.File = wx.Menu(title='')
        self.File.Append(wxID_PICKELEM,'Pick element','Pick element')
        self.File.Append(wxID_FILEEXIT,'Exit','Exit from testXNFF')
        self.Bind(wx.EVT_MENU,self.OnPickElement,id=wxID_PICKELEM)
        self.Bind(wx.EVT_MENU, self.OnFileExit, id=wxID_FILEEXIT)
        self.testFFMenu.Append(menu=self.File, title='File')
        self.SetMenuBar(self.testFFMenu)
        self.testFFPanel = wx.ScrolledWindow(self)        
        self.plotNB = plot.PlotNotebook()
        
    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Bind(wx.EVT_CLOSE, self.ExitMain)

    def ExitMain(self, event):
        sys.exit()
        
    def OnFileExit(self,event):
        self.Close()

    def OnPickElement(self,Parms):

        PE = G2elG.PickElement(self.testFFPanel,oneOnly=True)
        if PE.ShowModal() == wx.ID_OK:
            self.xrayFFs = G2el.GetFormFactorCoeff(PE.Elem)
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
                fplot.plot(lams,bp[0],label=isotope+El+" b'")
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
        
        self.plotNB.Show()
        
class testFFmain(wx.App):
    def OnInit(self):
        self.main = testXNFF(None)
        self.main.Show()
        self.SetTopWindow(self.main)
        return True

def main():
    'Starts main application to compute and plot form factors'
    application = testFFmain(0)
    application.MainLoop()
    
if __name__ == '__main__':
    GSASIIpath.InvokeDebugOpts()
    main()
