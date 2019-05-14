#!/usr/bin/env python

"""main Absorb routines
   Copyright: 2009, Robert B. Von Dreele (Argonne National Laboratory)
"""
from __future__ import division, print_function
import platform
import math
import wx
import numpy as np
import sys
import matplotlib as mpl
import GSASIIpath
GSASIIpath.SetVersionNumber("$Revision: 3765 $")
import GSASIIElem as G2elem
import GSASIIElemGUI as G2elemGUI

try:
    wx.NewIdRef
    wx.NewId = wx.NewIdRef
except AttributeError:
    pass

if '2' in platform.python_version_tuple()[0]:
    Gktheta = unichr(0x3b8)
    Gklambda = unichr(0x3bb)
    GkDelta = unichr(0x0394)
    Pwr10 = unichr(0x0b9)+unichr(0x2070)
    Pwr20 = unichr(0x0b2)+unichr(0x2070)
    Pwrm1 = unichr(0x207b)+unichr(0x0b9)
    Pwrm2 = unichr(0x207b)+unichr(0x0b2)
    Pwrm6 = unichr(0x207b)+unichr(0x2076)
    Pwrm4 = unichr(0x207b)+unichr(0x2074)
    Angstr = unichr(0x00c5)
    Gkmu = unichr(0x3bc)
    Pwr3 = unichr(0x0b3)
    Pwr4 = unichr(0x2074)
    Pwr20 = unichr(0x0b2)+unichr(0x0b0)
    Pwrm1 = unichr(0x207b)+unichr(0x0b9)

else:
    Gktheta = chr(0x3b8)
    Gklambda = chr(0x3bb)
    GkDelta = chr(0x0394)
    Pwr10 = chr(0x0b9)+chr(0x2070)
    Pwr20 = chr(0x0b2)+chr(0x2070)
    Pwrm1 = chr(0x207b)+chr(0x0b9)
    Pwrm2 = chr(0x207b)+chr(0x0b2)
    Pwrm6 = chr(0x207b)+chr(0x2076)
    Pwrm4 = chr(0x207b)+chr(0x2074)
    Angstr = chr(0x00c5)   
    Gkmu = chr(0x3bc)
    Pwr3 = chr(0x0b3)
    Pwr4 = chr(0x2074)
    Pwr20 = chr(0x0b2)+chr(0x0b0)
    Pwrm1 = chr(0x207b)+chr(0x0b9)

[wxID_CHOICE1, wxID_SPINTEXT1, wxID_SPINTEXT2, wxID_SPINTEXT3, wxID_SPINTEXT4,
 wxID_RESULTS,wxID_SLIDER1, wxID_SPINBUTTON, wxID_NUMELEM, wxID_SPINTEXT5,wxID_SPINTEXT6,
] = [wx.NewId() for _init_ctrls in range(11)]

[wxID_EXIT, wxID_DELETE, wxID_NEW, 
] = [wx.NewId() for _init_coll_ABSORB_Items in range(3)]
    
[wxID_KALPHAAGKA, wxID_KALPHACOKA, wxID_KALPHACRKA, 
 wxID_KALPHACUKA, wxID_KALPHAFEKA, wxID_KALPHAMNKA, 
 wxID_KALPHAMOKA, wxID_KALPHANIKA, wxID_KALPHAZNKA, 
] = [wx.NewId() for _init_coll_KALPHA_Items in range(9)]

[wxID_ABSORBABOUT] = [wx.NewId() for _init_coll_ABOUT_Items in range(1)]

class Absorb(wx.Frame):
    ''' '''
    Elems = []
    Wave = 1.5405      #CuKa default
    Kev = 12.397639    #keV for 1A x-rays
    for arg in sys.argv:
        if '-w' in arg:
            Wave = float(arg.split('-w')[1])
        elif '-e' in arg:
            E = float(arg.split('-e')[1])
            Wave = Kev/E
        elif '-h' in arg:
            print ('''
Absorb.py can take the following arguments:
-h   -  this help listing
-wv  -  set default wavelength to v, e.g. -w1.54 sets wavelength to 1.54A
-ev  -  set default energy to v, e.g. -e27 sets energy to 27keV
without arguments Absorb uses CuKa as default (Wave=1.54052A, E=8.0478keV)
''')
            sys.exit()
    Wmin = 0.05        #wavelength range
    Wmax = 3.0
    Wres = 0.004094    #plot resolution step size as const delta-lam/lam - gives 1000 steps for Wmin to Wmax
    Eres = 1.5e-4      #typical energy resolution for synchrotron x-ray sources
    Energy = Kev/Wave
    ifWave = True
    Volume = 0
    ifVol = False
    Zcell = 1
    Pack = 0.50
    Radius = 0.4
    def _init_coll_ABOUT_Items(self, parent):

        parent.Append(wxID_ABSORBABOUT,'About')
        self.Bind(wx.EVT_MENU, self.OnABOUTItems0Menu, id=wxID_ABSORBABOUT)

    def _init_coll_menuBar1_Menus(self, parent):

        parent.Append(menu=self.ABSORB, title='Absorb')
        parent.Append(menu=self.KALPHA, title='Kalpha')
        parent.Append(menu=self.ABOUT, title='About')

    def _init_coll_KALPHA_Items(self, parent):
        "Set of characteristic radiation from sealed tube sources"
        def OnCrkaMenu(event):
            self.SetWaveEnergy(2.28962)
    
        def OnMnkaMenu(event):
            self.SetWaveEnergy(2.10174)
    
        def OnFekaMenu(event):
            self.SetWaveEnergy(1.93597)
    
        def OnCokaMenu(event):
            self.SetWaveEnergy(1.78896)
    
        def OnNikaMenu(event):
            self.SetWaveEnergy(1.65784)
    
        def OnCukaMenu(event):
            self.SetWaveEnergy(1.54052)
    
        def OnZnkaMenu(event):
            self.SetWaveEnergy(1.43510)
    
        def OnMokaMenu(event):
            self.SetWaveEnergy(0.70926)
    
        def OnAgkaMenu(event):
            self.SetWaveEnergy(0.55936)
            
        parent.Append(wxID_KALPHACRKA, 'CrKa')
        parent.Append(wxID_KALPHAMNKA, 'MnKa')
        parent.Append(wxID_KALPHAFEKA, 'FeKa')
        parent.Append(wxID_KALPHACOKA, 'CoKa')
        parent.Append(wxID_KALPHANIKA, 'NiKa')
        parent.Append(wxID_KALPHACUKA, 'CuKa')
        parent.Append(wxID_KALPHAZNKA, 'ZnKa')
        parent.Append(wxID_KALPHAMOKA, 'MoKa')
        parent.Append(wxID_KALPHAAGKA, 'AgKa')
        self.Bind(wx.EVT_MENU, OnCrkaMenu, id=wxID_KALPHACRKA)
        self.Bind(wx.EVT_MENU, OnMnkaMenu, id=wxID_KALPHAMNKA)
        self.Bind(wx.EVT_MENU, OnFekaMenu, id=wxID_KALPHAFEKA)
        self.Bind(wx.EVT_MENU, OnCokaMenu, id=wxID_KALPHACOKA)
        self.Bind(wx.EVT_MENU, OnNikaMenu, id=wxID_KALPHANIKA)
        self.Bind(wx.EVT_MENU, OnCukaMenu, id=wxID_KALPHACUKA)
        self.Bind(wx.EVT_MENU, OnZnkaMenu, id=wxID_KALPHAZNKA)
        self.Bind(wx.EVT_MENU, OnMokaMenu, id=wxID_KALPHAMOKA)
        self.Bind(wx.EVT_MENU, OnAgkaMenu, id=wxID_KALPHAAGKA)

    def _init_coll_ABSORB_Items(self, parent):
        parent.Append(wxID_NEW,'&New Element','Add new element')
        self.Delete = parent.Append(wxID_DELETE,'&Delete Element','Delete an element')
        self.Delete.Enable(False)
        parent.Append(wxID_EXIT,'&Exit','Exit Fprime')
        self.Bind(wx.EVT_MENU, self.OnExitMenu, id=wxID_EXIT)
        self.Bind(wx.EVT_MENU, self.OnNewMenu, id=wxID_NEW)
        self.Bind(wx.EVT_MENU, self.OnDeleteMenu, id=wxID_DELETE)
        
    def _init_utils(self):
        self.ABSORB = wx.Menu(title='')

        self.KALPHA = wx.Menu(title='')
        self.KALPHA.SetEvtHandlerEnabled(True)

        self.ABOUT = wx.Menu(title='')

        self.menuBar1 = wx.MenuBar()

        self._init_coll_ABSORB_Items(self.ABSORB)
        self._init_coll_KALPHA_Items(self.KALPHA)
        self._init_coll_ABOUT_Items(self.ABOUT)
        self._init_coll_menuBar1_Menus(self.menuBar1)

    def _init_ctrls(self, parent):
        wx.Frame.__init__(self, parent=parent,
              size=wx.Size(500, 400),style=wx.DEFAULT_FRAME_STYLE ^ wx.CLOSE_BOX, title='Absorb')              
        self._init_utils()
        self.SetMenuBar(self.menuBar1)
        self.DrawPanel()
        
    def SetSize(self):
        w,h = self.GetClientSize()
        self.panel.SetSize(wx.Size(w,h))

    def DrawPanel(self):
        self.panel = wx.Panel(self)

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.Results = wx.TextCtrl( parent=self.panel,
            style=wx.TE_MULTILINE|wx.TE_DONTWRAP )
        self.Results.SetEditable(False)
        mainSizer.Add(self.Results,1,wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND)
        mainSizer.Add((10,15),0)
        
        if self.Elems:
            lablSizer = wx.BoxSizer(wx.HORIZONTAL)
            lablSizer.Add((5,10),0)
            lablSizer.Add(wx.StaticText(parent=self.panel,label='Chemical Formula:'),0,
                wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_LEFT)
            mainSizer.Add(lablSizer,0)
            mainSizer.Add((5,5),0)
            nRow = len(self.Elems)/5
            compSizer = wx.FlexGridSizer(nRow+1,10,0,0)
            for Elem in self.Elems:
                compSizer.Add(wx.StaticText(parent=self.panel,label="  "+Elem[0].capitalize(),
                    size=wx.Size(30,20)),0,wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT)
                numElem = wx.TextCtrl(id=wxID_NUMELEM,parent=self.panel,name=Elem[0],
                    size=wx.Size(70,20),value="%.2f" % (Elem[2]),style=wx.TE_PROCESS_ENTER)
                compSizer.Add(numElem,0)
                numElem.Bind(wx.EVT_TEXT_ENTER, self.OnNumElem, id=wxID_NUMELEM)
            mainSizer.Add(compSizer,0)
            mainSizer.Add((10,15),0)           

        selSizer = wx.BoxSizer(wx.HORIZONTAL)
        selSizer.Add((5,10),0)
        selSizer.Add(wx.StaticText(parent=self.panel, label='Wavelength:'),0,
            wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        selSizer.Add((5,10),0)
        self.SpinText1 = wx.TextCtrl(id=wxID_SPINTEXT1, parent=self.panel, 
            size=wx.Size(100,20), value = "%.4f" % (self.Wave),style=wx.TE_PROCESS_ENTER )
        selSizer.Add(self.SpinText1,0)
        selSizer.Add((5,10),0)
        self.SpinText1.Bind(wx.EVT_TEXT_ENTER, self.OnSpinText1, id=wxID_SPINTEXT1)
        
        selSizer.Add(wx.StaticText(parent=self.panel, label='Energy:'),0,
            wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        selSizer.Add((5,10),0)
        self.SpinText2 = wx.TextCtrl(id=wxID_SPINTEXT2, parent=self.panel, 
            size=wx.Size(100,20), value = "%.4f" % (self.Energy),style=wx.TE_PROCESS_ENTER) 
        selSizer.Add(self.SpinText2,0)
        selSizer.Add((5,10),0)
        self.SpinText2.Bind(wx.EVT_TEXT_ENTER, self.OnSpinText2, id=wxID_SPINTEXT2)
        
        selSizer.Add(wx.StaticText(parent=self.panel, label='Plot scale:'),
            0,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        selSizer.Add((5,10),0)
        self.choice1 = wx.ComboBox(id=wxID_CHOICE1, parent=self.panel, value='Wavelength',
             choices=['Wavelength','Energy'],style=wx.CB_READONLY|wx.CB_DROPDOWN)
        selSizer.Add(self.choice1,0)
        selSizer.Add((10,10),0)
        self.choice1.Bind(wx.EVT_COMBOBOX, self.OnChoice1, id=wxID_CHOICE1)
        mainSizer.Add(selSizer,0)
        mainSizer.Add((10,10),0)
        
        slideSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SpinButton = wx.SpinButton(id=wxID_SPINBUTTON, parent=self.panel, 
              size=wx.Size(25,24), style=wx.SP_VERTICAL | wx.SP_ARROW_KEYS)
        slideSizer.Add(self.SpinButton,0,wx.ALIGN_RIGHT)
        self.SpinButton.SetRange(-1,1)
        self.SpinButton.SetValue(0)
        self.SpinButton.Bind(wx.EVT_SPIN, self.OnSpinButton, id=wxID_SPINBUTTON)

        self.slider1 = wx.Slider(id=wxID_SLIDER1, maxValue=int(1000.*self.Wmax),
            minValue=int(1000.*self.Wmin), parent=self.panel,style=wx.SL_HORIZONTAL,
            value=int(self.Wave*1000.), )
        slideSizer.Add(self.slider1,1,wx.EXPAND|wx.ALIGN_RIGHT)
        self.slider1.Bind(wx.EVT_SLIDER, self.OnSlider1, id=wxID_SLIDER1)
        mainSizer.Add(slideSizer,0,wx.EXPAND)
        mainSizer.Add((10,10),0)
        
        cellSizer = wx.BoxSizer(wx.HORIZONTAL)
        cellSizer.Add((5,10),0)
        cellSizer.Add(wx.StaticText(parent=self.panel, label='Volume:'),0,
            wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        cellSizer.Add((5,10),0)
        self.SpinText3 = wx.TextCtrl(id=wxID_SPINTEXT3, parent=self.panel, 
              size=wx.Size(100,20), value = "%.2f" % (self.Volume),style=wx.TE_PROCESS_ENTER )
        cellSizer.Add(self.SpinText3,0)
        cellSizer.Add((5,10),0)
        self.SpinText3.Bind(wx.EVT_TEXT_ENTER, self.OnSpinText3, id=wxID_SPINTEXT3)
        
        cellSizer.Add((5,10),0)
        cellSizer.Add(wx.StaticText(parent=self.panel, label='Z(vol):'),0,
            wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        cellSizer.Add((5,10),0)
        self.SpinText4 = wx.TextCtrl(id=wxID_SPINTEXT4, parent=self.panel, 
              size=wx.Size(50,20), value = "%d" % (self.Zcell),style=wx.TE_PROCESS_ENTER )
        cellSizer.Add(self.SpinText4,0)
        cellSizer.Add((5,10),0)
        self.SpinText4.Bind(wx.EVT_TEXT_ENTER, self.OnSpinText4, id=wxID_SPINTEXT4)
        
        cellSizer.Add((5,10),0)
        cellSizer.Add(wx.StaticText(parent=self.panel, label='Sample R:'),0,
            wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        cellSizer.Add((5,10),0)
        self.SpinText5 = wx.TextCtrl(id=wxID_SPINTEXT5, parent=self.panel, 
              size=wx.Size(50,20), value = "%.2f" % (self.Radius),style=wx.TE_PROCESS_ENTER )
        cellSizer.Add(self.SpinText5,0)
        cellSizer.Add((5,10),0)
        self.SpinText5.Bind(wx.EVT_TEXT_ENTER, self.OnSpinText5, id=wxID_SPINTEXT5)

        cellSizer.Add((5,10),0)
        cellSizer.Add(wx.StaticText(parent=self.panel, label='packing:'),0,
            wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)
        cellSizer.Add((5,10),0)
        self.SpinText6 = wx.TextCtrl(id=wxID_SPINTEXT6, parent=self.panel, 
              size=wx.Size(50,20), value = "%.2f" % (self.Pack),style=wx.TE_PROCESS_ENTER )
        cellSizer.Add(self.SpinText6,0)
        cellSizer.Add((5,10),0)
        self.SpinText6.Bind(wx.EVT_TEXT_ENTER, self.OnSpinText6, id=wxID_SPINTEXT6)

        mainSizer.Add(cellSizer,0)
        mainSizer.Add((10,10),0)
        self.panel.SetSizer(mainSizer)
        self.panel.Fit()
        self.panel.GetParent().SetSize()

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.parent = parent
        self.Lines = []
        self.Elems = []
        self.linePicked = None

    def OnExitMenu(self, event):
        self.parent.G2plotNB.Delete('Absorb')
        self.Close()
        self.Destroy()

    def OnNewMenu(self, event):
        ElList = []
        for Elem in self.Elems: ElList.append(Elem[0])
        PE = G2elemGUI.PickElements(self,ElList)
        if PE.ShowModal() == wx.ID_OK:
            Elem = PE.Elem
        PE.Destroy()
        if Elem:
            for El in Elem:
                ElemSym = El.strip().upper()
                if ElemSym not in ElList:
                    atomData = G2elem.GetAtomInfo(ElemSym.capitalize())
                    FormFactors = G2elem.GetFormFactorCoeff(ElemSym)
                    for FormFac in FormFactors:
                        FormSym = FormFac['Symbol'].strip()
                        if FormSym == ElemSym:
                            Z = FormFac['Z']                #At. No.
                            N = 1.                          #no atoms / formula unit
                            Orbs = G2elem.GetXsectionCoeff(ElemSym)
                            Elem = [ElemSym,Z,N,FormFac,Orbs,atomData]
                    self.Elems.append(Elem)
            self.Delete.Enable(True)
            self.panel.Destroy()
            self.DrawPanel()
            self.NewFPPlot = True
            self.SetWaveEnergy(self.Wave)
            
    def OnDeleteMenu(self, event):
        if len(self.Elems):
            ElList = []
            for Elem in self.Elems: ElList.append(Elem[0])
            S = []
            DE = G2elemGUI.DeleteElement(self,ElList)
            if DE.ShowModal() == wx.ID_OK:
                El = DE.GetDeleteElement().strip().upper()
                for Elem in self.Elems:
                    if Elem[0] != El:
                        S.append(Elem)
                self.Elems = S
                self.CalcFPPS()
                if not self.Elems:
                    self.Delete.Enable(False)
                self.panel.Destroy()
                self.DrawPanel()
                self.NewFPPlot = True
                self.SetWaveEnergy(self.Wave)
        
    def OnNumElem(self, event):
        for Elem in self.Elems:
            if event.GetEventObject().GetName() == Elem[0]:
                Elem[2] = float(event.GetEventObject().GetValue())
                event.GetEventObject().SetValue("%8.2f" % (Elem[2]))
                self.SetWaveEnergy(self.Wave)                
        
    def OnSpinText1(self, event):
        self.SetWaveEnergy(float(self.SpinText1.GetValue()))
        
    def OnSpinText2(self, event):
        self.SetWaveEnergy(self.Kev/(float(self.SpinText2.GetValue())))
        
    def OnSpinText3(self,event):
        self.Volume = max(10.,float(self.SpinText3.GetValue()))
        self.ifVol = True
        self.SetWaveEnergy(self.Wave)
        
    def OnSpinText4(self,event):
        self.Zcell = max(1,float(self.SpinText4.GetValue()))
        self.SetWaveEnergy(self.Wave)
        
    def OnSpinText5(self, event):
        self.Radius = max(0.01,float(self.SpinText5.GetValue()))
        self.SetWaveEnergy(self.Wave)
       
    def OnSpinText6(self, event):
        self.Pack = min(1.0,max(0.01,float(self.SpinText6.GetValue())))
        self.SetWaveEnergy(self.Wave)
       
    def OnSpinButton(self, event):
        move = self.SpinButton.GetValue()/10000.
        self.Wave = min(max(self.Wave+move,self.Wmin),self.Wmax)
        self.SpinButton.SetValue(0)
        self.SetWaveEnergy(self.Wave)

    def OnSlider1(self, event):
        if self.ifWave:
            Wave = float(self.slider1.GetValue())/1000.
        else:
            Wave = self.Kev/(float(self.slider1.GetValue())/1000.)
        self.SetWaveEnergy(Wave)
        
    def SetWaveEnergy(self,Wave):
        self.Wave = Wave
        self.Energy = self.Kev/self.Wave
        self.Energy = round(self.Energy,4)
        E = self.Energy
        DE = E*self.Eres                         #smear by defined source resolution
        self.SpinText1.SetValue("%.4f" % (self.Wave))
        self.SpinText2.SetValue("%.4f" % (self.Energy))
        self.SpinText1.Update()
        self.SpinText2.Update()
        if self.ifWave:
            self.slider1.SetValue(int(1000.*self.Wave))
        else:
            self.slider1.SetValue(int(1000.*self.Energy))
        Text = ''
        if not self.ifVol:
            self.Volume = 0
            for Elem in self.Elems:
                self.Volume += 10.*Elem[2]
        muT = 0
        Mass = 0
        Fo = 0
        Fop = 0
        for Elem in self.Elems:
            Mass += self.Zcell*Elem[2]*Elem[5]['Mass']
            r1 = G2elem.FPcalc(Elem[4],E+DE)
            r2 = G2elem.FPcalc(Elem[4],E-DE)
            Els = Elem[0]
            Els = Els.ljust(2).lower().capitalize()
            mu = 0
            Fo += Elem[2]*Elem[1]
            if Elem[1] > 78 and self.Energy+DE > self.Kev/0.16:
                mu = self.Zcell*Elem[2]*(r1[2]+r2[2])/2.0
                Text += "%s\t%s%8.2f  %s%6s  %s%6.3f  %s%10.2f %s\n" %    (
                    'Element= '+str(Els),"N = ",Elem[2]," f'=",'not valid',
                    ' f"=',(r1[1]+r2[1])/2.0,' '+Gkmu+'=',mu,'barns')
            elif Elem[1] > 94 and self.Energy-DE < self.Kev/2.67:
                mu = 0
                Text += "%s\t%s%8.2f  %s%6s  %s%6s  %s%10s%s\n" %    (
                    'Element= '+str(Els),"N = ",Elem[2]," f'=",'not valid',
                    ' f"=','not valid',' '+Gkmu+'=','not valid')
            else:
                mu = self.Zcell*Elem[2]*(r1[2]+r2[2])/2.0
                Fop += Elem[2]*(Elem[1]+(r1[0]+r2[0])/2.0)
                Text += "%s\t%s%8.2f  %s%6.3f  %s%6.3f  %s%10.2f %s\n" %    (
                    'Element= '+str(Els),"N = ",Elem[2]," f'=",(r1[0]+r2[0])/2.0,
                    ' f"=',(r1[1]+r2[1])/2.0,' '+Gkmu+'=',mu,'barns')
            muT += mu
        
        if self.Volume:
            Text += "%s %s%10.2f %s" % ("Total",' '+Gkmu+'=',self.Pack*muT/self.Volume,'cm'+Pwrm1+', ')
            Text += "%s%10.2f%s" % ('Total '+Gkmu+'R=',self.Radius*self.Pack*muT/(10.0*self.Volume),', ')
            Text += "%s%10.4f%s\n" % ('Transmission exp(-2'+Gkmu+'R)=', \
                100.0*math.exp(-2*self.Radius*self.Pack*muT/(10.0*self.Volume)),'%')
            self.Results.SetValue(Text)
            den = Mass/(0.602*self.Volume)                
            if self.ifVol:
                Text += '%s' % ('Theor. density=')
            else:  
                Text += '%s' % ('Est. density=')
            Text += '%6.3f %s%.3f %s\n' % (den,'g/cm'+Pwr3+', Powder density=',self.Pack*den,'g/cm'+Pwr3)
            Text += '%s%10.2f%s\n'%('X-ray small angle scattering contrast',(28.179*Fo/self.Volume)**2,'*10'+Pwr20+'/cm'+Pwr4)
            if Fop:
                Text += '%s%10.2f%s\n'%('Anomalous X-ray small angle scattering contrast',(28.179*Fop/self.Volume)**2,'*10'+Pwr20+'/cm'+Pwr4)
            self.Results.SetValue(Text)
        self.Results.Update()
        self.SpinText3.SetValue("%.2f" % (self.Volume))
        self.SpinText3.Update()
        self.SpinText4.SetValue("%d" % (self.Zcell))
        self.SpinText4.Update()
        self.SpinText5.SetValue("%.2f" % (self.Radius))
        self.SpinText5.Update()
        self.SpinText6.SetValue("%.2f" % (self.Pack))
        self.SpinText6.Update()
        if len(self.Elems):
            self.CalcFPPS()
            self.UpDateAbsPlot(Wave,rePlot=True)

    def CalcFPPS(self):
        """generate f" curves for selected elements
           does constant delta-lambda/lambda steps over defined range
        """
        FPPS = []
        if self.Elems:
            wx.BeginBusyCursor()
            Corr = self.Zcell*self.Radius*self.Pack/(10.0*self.Volume)
            try:
                muT = []
                for iE,Elem in enumerate(self.Elems):
                    Els = Elem[0]
                    Els = Els = Els.ljust(2).lower().capitalize()
                    Wmin = self.Wmin
                    Wmax = self.Wmax
                    lWmin = math.log(Wmin)
                    N = int(round(math.log(Wmax/Wmin)/self.Wres))    #number of constant delta-lam/lam steps
                    I = range(N+1)
                    Ws = []
                    for i in I: Ws.append(math.exp(i*self.Wres+lWmin))
                    mus = []
                    Es = []
                    for j,W in enumerate(Ws):
                        E = self.Kev/W
                        DE = E*self.Eres                         #smear by defined source resolution
                        res1 = G2elem.FPcalc(Elem[4],E+DE)
                        res2 = G2elem.FPcalc(Elem[4],E-DE)
                        muR = Corr*Elem[2]*(res1[2]+res2[2])/2.0
                        mus.append(muR)
                        if iE:
                            muT[j] += muR
                        else:
                            muT.append(muR)
                        Es.append(E)
                    if self.ifWave:
                        Fpps = (Els,Ws,mus)
                    else:
                        Fpps = (Els,Es,mus)
                    FPPS.append(Fpps)
                if self.ifWave:
                    Fpps = ('Total',Ws,muT)
                else:
                    Fpps = ('Total',Es,muT)
                FPPS.append(Fpps)
            finally:
                wx.EndBusyCursor()
        self.FPPS = FPPS

    def OnChoice1(self, event):
        if event.GetString() == "Wavelength":
            self.ifWave = True
            self.NewFPPlot = True
            self.Wave = round(self.Wave,4)
            self.slider1.SetRange(int(1000.*self.Wmin),int(1000.*self.Wmax))
            self.slider1.SetValue(int(1000.*self.Wave))
            self.SpinText1.SetValue("%6.4f" % (self.Wave))
            self.SpinText2.SetValue("%7.4f" % (self.Energy))
        else:
            self.ifWave = False
            self.NewFPPlot = True
            Emin = self.Kev/self.Wmax
            Emax = self.Kev/self.Wmin
            self.Energy = round(self.Energy,4)
            self.slider1.SetRange(int(1000.*Emin),int(1000.*Emax))
            self.slider1.SetValue(int(1000.*self.Energy))
            self.SpinText1.SetValue("%6.4f" % (self.Wave))
            self.SpinText2.SetValue("%7.4f" % (self.Energy))
        if len(self.Elems):
            self.CalcFPPS()
            self.UpDateAbsPlot(self.Wave,rePlot=False)
        
    def OnKeyPress(self,event):
        if event.key == 'g':
            mpl.rcParams['axes.grid'] = not mpl.rcParams['axes.grid']
            self.UpDateAbsPlot(self.Wave,rePlot=False)

    def UpDateAbsPlot(self,Wave,rePlot=True):
        """Plot mu vs wavelength 0.05-3.0A"""
        xylim = []
        try:
            if rePlot:
                asb = self.Page.figure.get_axes()[1]
                xylim = asb.get_xlim(),asb.get_ylim()
            newPlot = False
        except:
            new,plotNum,self.Page,self.fplot,lim = self.parent.G2plotNB.FindPlotTab('Absorb','mpl')
            self.Page.canvas.mpl_connect('pick_event', self.OnPick)
            self.Page.canvas.mpl_connect('button_release_event', self.OnRelease)
            self.Page.canvas.mpl_connect('motion_notify_event', self.OnMotion)
            self.Page.canvas.mpl_connect('key_press_event', self.OnKeyPress)
            newPlot = True
        ax = self.Page.figure.add_subplot(111,label='absorb')
        self.fplot.set_visible(False)
        self.Page.Choice = (' key press','g: toggle grid',)
        self.Page.keyPress = self.OnKeyPress    
        ax.clear()
        ax.set_title('X-Ray Absorption',x=0,ha='left')
        ax.set_ylabel(r"$\mu R$",fontsize=14)
        Ymin = 0.0
        Ymax = 0.0
        if self.FPPS: 
            for Fpps in self.FPPS:
                Ymin = min(Ymin,min(Fpps[2]))
                Ymax = max(Ymax,max(Fpps[2]))
                fppsP1 = np.array(Fpps[1])
                fppsP2 = np.array(Fpps[2])
                ax.plot(fppsP1,fppsP2,label=r'$\mu R$ '+Fpps[0])
        if self.ifWave: 
            ax.set_xlabel(r'$\mathsf{\lambda, \AA}$',fontsize=14)
            ax.axvline(x=Wave,picker=3,color='black')
        else:
            ax.set_xlabel(r'$\mathsf{E, keV}$',fontsize=14)
            ax.set_xscale('log')
            ax.axvline(x=self.Kev/Wave,picker=3,color='black')
        ax.axhline(y=1.0,color='b')
        ax.axhline(y=5.0,color='r')
        ax.set_ylim(Ymin,Ymax)
        if self.FPPS:
            ax.legend(loc='best')
        if newPlot:
            newPlot = False
            self.Page.canvas.draw()
        else:
            if rePlot:
                tb = self.Page.canvas.toolbar
                tb.push_current()
                ax.set_xlim(xylim[0])
                ax.set_ylim(xylim[1])
                xylim = []
                tb.push_current()
            self.Page.canvas.draw()
        
    def OnPick(self, event):
        self.linePicked = event.artist
        
    def OnMotion(self,event):
        xpos = event.xdata
        if xpos and xpos>0.1:
            ypos = event.ydata
            if self.ifWave:
                Wave = xpos
            else:
                Wave = self.Kev/xpos
            Wave = min(max(Wave,self.Wmin),self.Wmax)
            self.parent.G2plotNB.status.SetStatusText('Wavelength: %.4f, Energy: %.3f, %sR: %.3f'%(Wave,self.Kev/Wave,Gkmu,ypos),1)
        if self.linePicked:
            self.SetWaveEnergy(Wave)
                
    def OnRelease(self, event):
        if self.linePicked is None: return
        self.linePicked = None
        xpos = event.xdata
        if xpos:
            if self.ifWave:
                Wave = xpos
            else:
                Wave = self.Kev/xpos               
            self.SetWaveEnergy(Wave)
            
    def OnABOUTItems0Menu(self, event):
        ''' '''
        try:
            import wx.adv as wxadv  # AboutBox moved here in Phoenix
        except:
            wxadv = wx
        info = wxadv.AboutDialogInfo()
        info.Name = 'Absorb'
        info.Copyright = '''
Robert B. Von Dreele, 2009(C)
Argonne National Laboratory
This product includes software developed 
by the UChicago Argonne, LLC, as 
Operator of Argonne National Laboratory.        '''
        info.Description = '''
For calculating X-ray absorption factors to 250keV for cylindrical      
powder samples; based on Fortran program Fprime of Cromer & Liberman 
corrected for Kissel & Pratt energy term; Jensen term not included
        '''
        wxadv.AboutBox(info)

