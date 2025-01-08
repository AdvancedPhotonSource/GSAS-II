# -*- coding: utf-8 -*-
'''Routines for Periodic table wx.Frame follow.
'''
from __future__ import division, print_function
from . import GSASIIpath
import wx
import os
import wx.lib.colourselect as wscs
# for Python 3.10+, define a version of wxPoint that accepts float params
def wxPoint(x,y): return wx.Point(int(x),int(y))
class PickElement(wx.Dialog):
    '''Makes periodic table widget for picking element. Modes:
        oneOnly if True element symbols are provided, otherwise select valence
        ifNone if True show None button
        ifMag if True present magnetic scatters only
        ifOrbs if True present orbital form actors only
        multiple if True multiple elements can be selected
    '''
    Elem=None
    def _init_ctrls(self,prnt,ifMag=False,ifOrbs=False):
        wx.Dialog.__init__(self, id=-1, name='PickElement',
              parent=prnt, pos=wx.DefaultPosition,
              style=wx.DEFAULT_DIALOG_STYLE, title='Pick Element')
        from . import ElementTable as ET
        self.butWid = 90
        if 'nt' in os.name:
            self.butWid = 50
        self.SetClientSize(wx.Size(50+18*self.butWid, 250))
        self.Centre()
        i=0
        Elems = ET.ElTable
        if ifMag:
            Elems = ET.MagElTable
        if ifOrbs:
            Elems = ET.OrbsElTable
        for E in Elems:
            if E[1] < 0: continue
            if self.oneOnly:
                color=E[4]
            else:
                color=E[6]
            self.ElButton(name=E[0],
               pos=wxPoint(E[1]*self.butWid+25,E[2]*24+24),
               tip=E[3],color=color)
            i+=1
        if self.multiple:
            b = wx.Button(self,wx.ID_CLOSE,
                pos=wxPoint(16.5*self.butWid+25,7.75*24+24),label="Done")
            b.Bind(wx.EVT_BUTTON, self.OnClose)

    def __init__(self, parent,oneOnly=False,ifNone=False,ifMag=False,ifOrbs=False,multiple=False):
        self.oneOnly = oneOnly
        self.ifNone = ifNone
        self.multiple = multiple
        self._init_ctrls(parent,ifMag=ifMag,ifOrbs=ifOrbs)
        self.elementList = []

    def ElButton(self, name, pos, tip, color):
        'Creates an element button widget'
        self.color = color
        if not self.ifNone and name[0] == 'None':
            return
        if self.oneOnly:
            El = wscs.ColourSelect(label=name[0], parent=self,colour=color,
                pos=pos, size=wx.Size(self.butWid,23), style=wx.RAISED_BORDER)
            El.Bind(wx.EVT_BUTTON, self.OnElButton)
        else:
            butWid = self.butWid
            if name[0] == 'None':
                butWid *= 2
            # patch for wx 2.9+ where EVT_COMBOBOX happens only on a value change.
            i,j= wx.__version__.split('.')[0:2]
            if int(i)+int(j)/10. > 2.8:
                startname = name[0]+' '  # add an invisible space
            else:
                startname = name[0]
            # Not ideal because wx.CB_READONLY is better.
            El = wx.ComboBox(choices=name, parent=self, pos=pos, size=wx.Size(butWid,27),
                    style=wx.CB_DROPDOWN, value=startname)
            #El = wx.ComboBox(choices=name, parent=self, pos=pos, size=wx.Size(butWid,23),
            #        style=wx.CB_READONLY, value=name[0])
            if sum(color)/3 < 150 and color[1] < 150: # background is mostly dark, use white letters
                El.SetForegroundColour((255,255,255))
            else:
                El.SetForegroundColour((10,10,10))
            El.Bind(wx.EVT_COMBOBOX,self.OnElButton)

        El.SetBackgroundColour(color)
        if 'phoenix' in wx.version():
            El.SetToolTip(tip)
        else:
            El.SetToolTipString(tip)

    def OnElButton(self, event):
        if self.oneOnly:
            El = event.GetEventObject().GetLabel()
        else:
            El = event.GetEventObject().GetValue()
        self.Elem = El
        if self.multiple:
            if El in self.elementList:
                self.elementList.remove(El)
                event.GetEventObject().SetBackgroundColour(self.color) # Shows on Mac
            else:
                self.elementList.append(El)
                event.GetEventObject().SetBackgroundColour('black') # Shows on Mac
            event.GetEventObject().SetColour(
                wx.Colour(*[int(i/2) for i in event.GetEventObject().GetColour()]))
        else:
            self.EndModal(wx.ID_OK)

    def OnClose(self,event):
        self.EndModal(wx.ID_OK)

class PickElements(wx.Dialog):
    """Makes periodic table widget for picking elements - caller maintains element list"""
    Elem = []
    def _init_ctrls(self, prnt,list):
        wx.Dialog.__init__(self, id=-1, name='PickElements',
              parent=prnt, pos=wx.DefaultPosition, size=wx.Size(580, 360),
              style=wx.DEFAULT_DIALOG_STYLE, title='Pick Elements')
        panel = wx.Panel(self)

        REcolor = wx.Colour(128, 128, 255)
        Metcolor = wx.Colour(192, 192, 192)
        Noblecolor = wx.Colour(255, 128, 255)
        Alkcolor = wx.Colour(255, 255, 128)
        AlkEcolor = wx.Colour(255, 128, 0)
        SemMetcolor = wx.Colour(128, 255, 0)
        NonMetcolor = wx.Colour(0, 255, 255)
        White = wx.Colour(255, 255, 255)
        self.Elem = []
        for El in list:
            self.Elem.append(El.lower().capitalize())

        self.ElTable = [
            ("H",   0,0, "Hydrogen",    White,           0.0000),
            ("He", 17,0, "Helium",      Noblecolor,      0.0000),
            ("Li",  0,1, "Lithium",     Alkcolor,        0.0004),
            ("Be",  1,1, "Beryllium",   AlkEcolor,       0.0006),
            ("B",  12,1, "Boron",       NonMetcolor,     0.0012),
            ("C",  13,1, "Carbon",      NonMetcolor,     0.0018),
            ("N",  14,1, "Nitrogen",    NonMetcolor,     0.0030),
            ("O",  15,1, "Oxygen",      NonMetcolor,     0.0042),
            ("F",  16,1, "Fluorine",    NonMetcolor,     0.0054),
            ("Ne", 17,1, "Neon",        Noblecolor,      0.0066),
            ("Na",  0,2, "Sodium",      Alkcolor,        0.0084),
            ("Mg",  1,2, "Magnesium",   AlkEcolor,       0.0110),
            ("Al", 12,2, "Aluminum",    SemMetcolor,     0.0125),
            ("Si", 13,2, "Silicon",     NonMetcolor,     0.0158),
            ("P",  14,2, "Phosphorus",  NonMetcolor,     0.0180),
            ("S",  15,2, "Sulphur",     NonMetcolor,     0.0210),
            ("Cl", 16,2, "Chlorine",    NonMetcolor,     0.0250),
            ("Ar", 17,2, "Argon",       Noblecolor,      0.0285),
            ("K",   0,3, "Potassium",   Alkcolor,        0.0320),
            ("Ca",  1,3, "Calcium",     AlkEcolor,       0.0362),
            ("Sc",  2,3, "Scandium",    Metcolor,        0.0410),
            ("Ti",  3,3, "Titanium",    Metcolor,        0.0460),
            ("V",   4,3, "Vanadium",    Metcolor,        0.0510),
            ("Cr",  5,3, "Chromium",    Metcolor,        0.0560),
            ("Mn",  6,3, "Manganese",   Metcolor,        0.0616),
            ("Fe",  7,3, "Iron",        Metcolor,        0.0680),
            ("Co",  8,3, "Cobalt",      Metcolor,        0.0740),
            ("Ni",  9,3, "Nickel",      Metcolor,        0.0815),
            ("Cu", 10,3, "Copper",      Metcolor,        0.0878),
            ("Zn", 11,3, "Zinc",        Metcolor,        0.0960),
            ("Ga", 12,3, "Gallium",     SemMetcolor,      0.104),
            ("Ge", 13,3, "Germanium",   SemMetcolor,      0.114),
            ("As", 14,3, "Arsenic",     NonMetcolor,      0.120),
            ("Se", 15,3, "Selenium",    NonMetcolor,      0.132),
            ("Br", 16,3, "Bromine",     NonMetcolor,      0.141),
            ("Kr", 17,3, "Krypton",     Noblecolor,       0.150),
            ("Rb",  0,4, "Rubidium",    Alkcolor,         0.159),
            ("Sr",  1,4, "Strontium",   AlkEcolor,        0.171),
            ("Y",   2,4, "Yittrium",    Metcolor,         0.180),
            ("Zr",  3,4, "Zirconium",   Metcolor,         0.192),
            ("Nb",  4,4, "Niobium",     Metcolor,         0.204),
            ("Mo",  5,4, "Molybdenium", Metcolor,         0.216),
            ("Tc",  6,4, "Technetium",  Metcolor,         0.228),
            ("Ru",  7,4, "Ruthenium",   Metcolor,         0.246),
            ("Rh",  8,4, "Rhodium",     Metcolor,         0.258),
            ("Pd",  9,4, "Palladium",   Metcolor,         0.270),
            ("Ag", 10,4, "Silver",      Metcolor,         0.285),
            ("Cd", 11,4, "Cadmium",     Metcolor,         0.300),
            ("In", 12,4, "Indium",      SemMetcolor,      0.318),
            ("Sn", 13,4, "Tin",         SemMetcolor,      0.330),
            ("Sb", 14,4, "Antimony",    SemMetcolor,      0.348),
            ("Te", 15,4, "Tellurium",   NonMetcolor,      0.363),
            ("I",  16,4, "Iodine",      NonMetcolor,      0.384),
            ("Xe", 17,4, "Xenon",       Noblecolor,       0.396),
            ("Cs",  0,5, "Caesium",     Alkcolor,         0.414),
            ("Ba",  1,5, "Barium",      AlkEcolor,        0.438),
            ("La",  2,5, "Lanthanium",  Metcolor,         0.456),
            ("Ce",  3.5,6.5, "Cerium",      REcolor,      0.474),
            ("Pr",  4.5,6.5, "Praseodymium",REcolor,      0.492),
            ("Nd",  5.5,6.5, "Neodymium",   REcolor,      0.516),
            ("Pm",  6.5,6.5, "Promethium",  REcolor,      0.534),
            ("Sm",  7.5,6.5, "Samarium",    REcolor,      0.558),
            ("Eu",  8.5,6.5, "Europium",    REcolor,      0.582),
            ("Gd",  9.5,6.5, "Gadolinium",  REcolor,      0.610),
            ("Tb", 10.5,6.5, "Terbium",     REcolor,      0.624),
            ("Dy", 11.5,6.5, "Dysprosium",  REcolor,      0.648),
            ("Ho", 12.5,6.5, "Holmium",     REcolor,      0.672),
            ("Er", 13.5,6.5, "Erbium",      REcolor,      0.696),
            ("Tm", 14.5,6.5, "Thulium",     REcolor,      0.723),
            ("Yb", 15.5,6.5, "Ytterbium",   REcolor,      0.750),
            ("Lu", 16.5,6.5, "Lutetium",    REcolor,      0.780),
            ("Hf",  3,5, "Hafnium",     Metcolor,         0.804),
            ("Ta",  4,5, "Tantalum",    Metcolor,         0.834),
            ("W",   5,5, "Tungsten",    Metcolor,         0.864),
            ("Re",  6,5, "Rhenium",     Metcolor,         0.900),
            ("Os",  7,5, "Osmium",      Metcolor,         0.919),
            ("Ir",  8,5, "Iridium",     Metcolor,         0.948),
            ("Pt",  9,5, "Platinium",   Metcolor,         0.984),
            ("Au", 10,5, "Gold",        Metcolor,         1.014),
            ("Hg", 11,5, "Mercury",     Metcolor,         1.046),
            ("Tl", 12,5, "Thallium",    SemMetcolor,      1.080),
            ("Pb", 13,5, "Lead",        SemMetcolor,      1.116),
            ("Bi", 14,5, "Bismuth",     SemMetcolor,      1.149),
            ("Po", 15,5, "Polonium",    SemMetcolor,      1.189),
            ("At", 16,5, "Astatine",    NonMetcolor,      1.224),
            ("Rn", 17,5, "Radon",       Noblecolor,       1.260),
            ("Fr",  0,6, "Francium",    Alkcolor,         1.296),
            ("Ra",  1,6, "Radium",      AlkEcolor,        1.332),
            ("Ac",  2,6, "Actinium",    Metcolor,         1.374),
            ("Th",  3.5,7.5, "Thorium",     REcolor,      1.416),
            ("Pa",  4.5,7.5, "Protactinium",REcolor,      1.458),
            ("U",   5.5,7.5, "Uranium",     REcolor,      1.470),
            ("Np",  6.5,7.5, "Neptunium",   REcolor,      1.536),
            ("Pu",  7.5,7.5, "Plutonium",   REcolor,      1.584),
            ("Am",  8.5,7.5, "Americium",   REcolor,      1.626),
            ("Cm",  9.5,7.5, "Curium",      REcolor,      1.669),
            ("Bk", 10.5,7.5, "Berkelium",   REcolor,      1.716),
            ("Cf", 11.5,7.5, "Californium", REcolor,      1.764)]

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        elPanel = wx.Panel(panel)

        i=0
        for E in self.ElTable:
            PickElements.ElButton(self,parent=elPanel,name=E[0],
                pos=wxPoint(E[1]*30+20,E[2]*30+25),tip=E[3],color=E[4])
            i+=1
        mainSizer.Add(elPanel,0,wx.EXPAND)
        mainSizer.Add((10,10),0)

        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        OkBtn = wx.Button(panel,-1,"Ok")
        OkBtn.Bind(wx.EVT_BUTTON, self.OnOk)
        cancelBtn = wx.Button(panel,-1,"Cancel")
        cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)
        btnSizer.Add((20,20),1)
        btnSizer.Add(OkBtn)
        btnSizer.Add((20,20),1)
        btnSizer.Add(cancelBtn)
        btnSizer.Add((20,20),1)
        mainSizer.Add(btnSizer,0,wx.EXPAND|wx.BOTTOM, 10)
        panel.SetSizer(mainSizer)
        panel.Fit()

    def OnOk(self,event):
        if self.Elem:
            self.EndModal(wx.ID_OK)
        else:
            self.EndModal(wx.ID_CANCEL)

    def OnCancel(self,event):
        self.EndModal(wx.ID_CANCEL)

    def __init__(self, parent,list):
        self._init_ctrls(parent,list)

    def ElButton(self, parent, name, pos, tip, color):
        Black = wx.Colour(0,0,0)
        if name in self.Elem:
            color = Black
        El = wscs.ColourSelect(label=name, parent=parent,colour=color,
            pos=pos, size=wx.Size(32, 32), style=wx.RAISED_BORDER)
        El.SetBackgroundColour(color)
        El.SetLabel(name)
        if 'phoenix' in wx.version():
            El.SetToolTip(tip)
        else:
            El.SetToolTipString(tip)
        El.Bind(wx.EVT_BUTTON, self.OnElButton)

    def OnElButton(self, event):
        Black = wx.Colour(0,0,0)
        btn = event.GetEventObject()
        El = btn.GetLabel()
        if btn.GetColour() != Black:
            for Elem in self.ElTable:
                if El in Elem:
                    ElColor = Elem[4]
            if El in self.Elem:
                btn.SetColour(ElColor)
                self.Elem.remove(El)
            else:
                btn.SetColour(wx.Colour(255,0,0))
                self.Elem.append(El)

class DeleteElement(wx.Dialog):
    "Delete element from selected set widget"
    def _init_ctrls(self, parent,choice):
        l = len(choice)-1
        wx.Dialog.__init__(self, id=-1, name='Delete', parent=parent,
              pos=wx.DefaultPosition, size=wx.Size(max(128,64+l*24), 87),
              style=wx.DEFAULT_DIALOG_STYLE, title='Delete Element')
        self.Show(True)
        self.SetAutoLayout(True)
        self.SetHelpText('Select element to delete')
        self.SetWindowVariant(wx.WINDOW_VARIANT_SMALL)

        i = 0
        Elem = []
        for Elem in choice:
            self.ElButton(name=Elem,pos=wxPoint(16+i*24, 16))
            i+=1

    def __init__(self, parent,choice):
        DeleteElement.El = ' '
        self._init_ctrls(parent,choice)

    def ElButton(self, name, pos):
        'Needs a doc string'
        White = wx.Colour(255, 255, 255)
        El = wscs.ColourSelect(label=name, parent=self, colour = White,
            pos=pos, size=wx.Size(24, 23), style=wx.RAISED_BORDER)
        El.Bind(wx.EVT_BUTTON, self.OnDeleteButton)

    def OnDeleteButton(self, event):
        DeleteElement.El=event.GetEventObject().GetLabel()
        self.EndModal(wx.ID_OK)

    def GetDeleteElement(self):
        return DeleteElement.El

if __name__ == '__main__':
    app = wx.PySimpleApp()
    GSASIIpath.InvokeDebugOpts()
    G2frame = wx.Frame(None) # create a frame
    G2frame.Show(True)

    PE = PickElement(G2frame)
    if PE.ShowModal() == wx.ID_OK:
        print (PE.Elem)
    PE.Destroy()
