# -*- coding: utf-8 -*-
"""ElementGUI: class defn. for element GUIs
   Copyright: 2008, Robert B. Von Dreele & Brian H. Toby (Argonne National Laboratory)
"""
########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################
import wx
import os
import wx.lib.colourselect as wscs
class PickElement(wx.Dialog):
    "Makes periodic table widget for picking element - caller maintains element list"
    Elem=None
    def _init_ctrls(self, prnt):
        wx.Dialog.__init__(self, id=-1, name='PickElement',
              parent=prnt, pos=wx.DefaultPosition, 
              style=wx.DEFAULT_DIALOG_STYLE, title='Pick Element')
        import ElementTable as ET
        self.butWid = 55
        if 'nt' in os.name:
            self.butWid = 40
        self.SetClientSize(wx.Size(50+18*self.butWid, 250))
        
        i=0
        for E in ET.ElTable:
            if self.oneOnly:
                color=E[4]
            else:
                color=E[6]
            PickElement.ElButton(self,name=E[0],
               pos=wx.Point(E[1]*self.butWid+25,E[2]*24+24),tip=E[3],color=color)
            i+=1

    def __init__(self, parent,oneOnly=False,ifNone=False):
        self.oneOnly = oneOnly
        self.ifNone = ifNone
        self._init_ctrls(parent)
        
    def ElButton(self, name, pos, tip, color):
        Black = wx.Colour(0,0,0)
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
            El = wx.ComboBox(choices=name, parent=self, pos=pos, size=wx.Size(butWid,23),
                style=wx.CB_READONLY, value=name[0])
            El.Bind(wx.EVT_COMBOBOX,self.OnElButton)
        
        El.SetBackgroundColour(color)
        El.SetToolTipString(tip)

    def OnElButton(self, event):
        if self.oneOnly:
            El = event.GetEventObject().GetLabel()
        else:
            El = event.GetEventObject().GetValue()
        self.Elem = El
        self.EndModal(wx.ID_OK)        
        
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
            self.ElButton(id=-1,name=Elem,pos=wx.Point(16+i*24, 16))
            i+=1
              
    def __init__(self, parent,choice):
        DeleteElement.El = ' '
        self._init_ctrls(parent,choice)

    def ElButton(self, id, name, pos):
        White = wx.Colour(255, 255, 255)
        El = wscs.ColourSelect(label=name, parent=self, colour = White,
            pos=pos, size=wx.Size(24, 23), style=wx.RAISED_BORDER)
        El.Bind(wx.EVT_BUTTON, self.OnDeleteButton)
    
    def OnDeleteButton(self, event):
        DeleteElement.El=event.GetEventObject().GetLabel()
        self.EndModal(wx.ID_OK)
        
    def GetDeleteElement(self):
        return DeleteElement.El
        

